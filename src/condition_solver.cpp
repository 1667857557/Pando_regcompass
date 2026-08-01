#include "condition_sparse_input.h"
#include "condition_solver_internal.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

using Eigen::MatrixXd;
using Eigen::SparseMatrix;
using Eigen::VectorXd;

namespace {

struct SmoothResult {
    double value;
    MatrixXd gradient;
    VectorXd intercept;
};

std::vector<SparseMatrix<double>> as_sparse_list(const Rcpp::List& input) {
    std::vector<SparseMatrix<double>> answer;
    answer.reserve(input.size());
    for (R_xlen_t i = 0; i < input.size(); ++i) {
        SEXP value = input[i];
        if (Rf_isS4(value)) {
            answer.emplace_back(pando_condition_as_dgCMatrix(
                value, "X_list sparse element"
            ));
        } else {
            Rcpp::NumericMatrix dense(value);
            Eigen::Map<MatrixXd> mapped(
                dense.begin(), dense.nrow(), dense.ncol()
            );
            answer.emplace_back(mapped.sparseView());
        }
        answer.back().makeCompressed();
    }
    return answer;
}

std::vector<VectorXd> as_vector_list(const Rcpp::List& input) {
    std::vector<VectorXd> answer;
    answer.reserve(input.size());
    for (R_xlen_t i = 0; i < input.size(); ++i) {
        Rcpp::NumericVector value(input[i]);
        Eigen::Map<VectorXd> mapped(value.begin(), value.size());
        answer.emplace_back(mapped);
    }
    return answer;
}

SmoothResult profiled_smooth(
    const MatrixXd& B,
    const std::vector<SparseMatrix<double>>& X,
    const std::vector<VectorXd>& y,
    const VectorXd& loss_weights,
    double ridge
) {
    const int p = B.rows();
    const int tasks = B.cols();
    SmoothResult out{
        0.0,
        MatrixXd::Zero(p, tasks),
        VectorXd::Zero(tasks)
    };
    for (int task = 0; task < tasks; ++task) {
        VectorXd xb = X[task] * B.col(task);
        out.intercept[task] = (y[task] - xb).mean();
        VectorXd residual =
            y[task].array() - out.intercept[task] - xb.array();
        out.value +=
            0.5 * loss_weights[task] * residual.squaredNorm();
        out.gradient.col(task).noalias() =
            -loss_weights[task] * (X[task].transpose() * residual);
    }
    if (ridge > 0.0) {
        out.value += 0.5 * ridge * B.squaredNorm();
        out.gradient.noalias() += ridge * B;
    }
    return out;
}

double sparse_group_penalty(
    const MatrixXd& B,
    double lambda,
    double alpha,
    double condition_mix
) {
    const double strength = lambda * alpha;
    double group = 0.0;
    for (int row = 0; row < B.rows(); ++row) {
        group += B.row(row).norm();
    }
    const double element = B.cwiseAbs().sum();
    return strength * (
        (1.0 - condition_mix) * group + condition_mix * element
    );
}

MatrixXd sparse_group_prox(
    const MatrixXd& V,
    const Eigen::ArrayXXi& mask,
    double step,
    double lambda,
    double alpha,
    double condition_mix
) {
    MatrixXd U = MatrixXd::Zero(V.rows(), V.cols());
    const double element_threshold =
        step * lambda * alpha * condition_mix;
    const double group_threshold =
        step * lambda * alpha * (1.0 - condition_mix);
    for (int row = 0; row < V.rows(); ++row) {
        for (int task = 0; task < V.cols(); ++task) {
            if (!mask(row, task)) continue;
            const double value = V(row, task);
            const double magnitude = std::max(
                std::abs(value) - element_threshold,
                0.0
            );
            U(row, task) = std::copysign(magnitude, value);
        }
        const double norm = U.row(row).norm();
        if (norm > 0.0) {
            const double scale = std::max(
                1.0 - group_threshold / norm,
                0.0
            );
            U.row(row) *= scale;
        }
    }
    return U;
}

double centered_spectral_squared_norm(const SparseMatrix<double>& X) {
    if (X.rows() == 0 || X.cols() == 0) return 0.0;
    VectorXd direction(X.cols());
    for (int column = 0; column < X.cols(); ++column) {
        direction[column] = 1.0 + static_cast<double>(column % 7) / 7.0;
    }
    direction.normalize();
    double previous = -1.0;
    for (int iteration = 0; iteration < 30; ++iteration) {
        VectorXd projected = X * direction;
        projected.array() -= projected.mean();
        const double estimate = projected.squaredNorm();
        VectorXd next = X.transpose() * projected;
        const double next_norm = next.norm();
        if (!std::isfinite(estimate) || !std::isfinite(next_norm)) {
            Rcpp::stop("Centered spectral-norm iteration became non-finite.");
        }
        if (next_norm <= std::numeric_limits<double>::epsilon()) {
            return std::max(estimate, 0.0);
        }
        direction = next / next_norm;
        if (previous >= 0.0 &&
            std::abs(estimate - previous) <=
                1e-8 * std::max(1.0, estimate)) {
            break;
        }
        previous = estimate;
    }
    VectorXd projected = X * direction;
    projected.array() -= projected.mean();
    return std::max(projected.squaredNorm(), 0.0);
}

double initial_step_for_lambda(
    const VectorXd& centered_spectral_norm,
    const VectorXd& loss_weights,
    double ridge
) {
    double lipschitz = ridge;
    for (int task = 0; task < centered_spectral_norm.size(); ++task) {
        const double task_lipschitz = ridge +
            1.05 * loss_weights[task] * centered_spectral_norm[task];
        lipschitz = std::max(lipschitz, task_lipschitz);
    }
    return (!std::isfinite(lipschitz) || lipschitz <= 0.0) ?
        1.0 : 1.0 / lipschitz;
}

Rcpp::List fit_one_lambda(
    const std::vector<SparseMatrix<double>>& X,
    const std::vector<VectorXd>& y,
    const VectorXd& loss_weights,
    const Eigen::ArrayXXi& mask,
    double lambda,
    double alpha,
    double condition_mix,
    const MatrixXd& initial_B,
    double initial_step,
    double default_step,
    int max_iter,
    double tol_objective,
    double tol_coef,
    bool keep_history
) {
    const double ridge = lambda * (1.0 - alpha);
    MatrixXd B = initial_B;
    for (int row = 0; row < B.rows(); ++row) {
        for (int task = 0; task < B.cols(); ++task) {
            if (!mask(row, task)) B(row, task) = 0.0;
        }
    }
    MatrixXd Z = B;
    double acceleration = 1.0;
    double step = default_step;
    if (std::isfinite(initial_step) && initial_step > 0.0) {
        step = std::max(initial_step, default_step);
    }
    SmoothResult smooth_B =
        profiled_smooth(B, X, y, loss_weights, ridge);
    double objective_previous = smooth_B.value +
        sparse_group_penalty(B, lambda, alpha, condition_mix);
    std::vector<double> history;
    if (keep_history) history.push_back(objective_previous);
    bool converged = false;
    double coef_change = std::numeric_limits<double>::infinity();
    double objective_change = std::numeric_limits<double>::infinity();
    int iteration = 0;
    const double backtrack = 0.5;
    const double min_step = 1e-14;

    for (iteration = 1; iteration <= max_iter; ++iteration) {
        SmoothResult smooth_Z =
            Z.isApprox(B, 0.0) ? smooth_B :
            profiled_smooth(Z, X, y, loss_weights, ridge);
        MatrixXd candidate;
        SmoothResult smooth_candidate;
        while (true) {
            candidate = sparse_group_prox(
                Z - step * smooth_Z.gradient,
                mask,
                step,
                lambda,
                alpha,
                condition_mix
            );
            smooth_candidate = profiled_smooth(
                candidate, X, y, loss_weights, ridge
            );
            MatrixXd difference = candidate - Z;
            const double quadratic_bound = smooth_Z.value +
                (smooth_Z.gradient.array() * difference.array()).sum() +
                difference.squaredNorm() / (2.0 * step);
            if (smooth_candidate.value <= quadratic_bound + 1e-10) break;
            step *= backtrack;
            if (step < min_step) {
                Rcpp::stop("Backtracking line search reached min_step.");
            }
        }
        double objective_candidate = smooth_candidate.value +
            sparse_group_penalty(
                candidate, lambda, alpha, condition_mix
            );
        if (objective_candidate > objective_previous + 1e-10) {
            acceleration = 1.0;
            Z = B;
            smooth_Z = smooth_B;
            while (true) {
                candidate = sparse_group_prox(
                    Z - step * smooth_Z.gradient,
                    mask,
                    step,
                    lambda,
                    alpha,
                    condition_mix
                );
                smooth_candidate = profiled_smooth(
                    candidate, X, y, loss_weights, ridge
                );
                MatrixXd difference = candidate - Z;
                const double quadratic_bound = smooth_Z.value +
                    (smooth_Z.gradient.array() * difference.array()).sum() +
                    difference.squaredNorm() / (2.0 * step);
                if (smooth_candidate.value <= quadratic_bound + 1e-10) break;
                step *= backtrack;
                if (step < min_step) {
                    Rcpp::stop(
                        "Backtracking line search reached min_step after restart."
                    );
                }
            }
            objective_candidate = smooth_candidate.value +
                sparse_group_penalty(
                    candidate, lambda, alpha, condition_mix
                );
        }
        coef_change = (candidate - B).norm() /
            (B.norm() + std::numeric_limits<double>::epsilon());
        objective_change =
            std::abs(objective_candidate - objective_previous) /
            (std::abs(objective_previous) +
             std::numeric_limits<double>::epsilon());
        if (keep_history) history.push_back(objective_candidate);
        MatrixXd B_previous = B;
        B = candidate;
        smooth_B = smooth_candidate;
        objective_previous = objective_candidate;
        if (objective_change < tol_objective && coef_change < tol_coef) {
            converged = true;
            break;
        }
        const double acceleration_new =
            (1.0 + std::sqrt(
                1.0 + 4.0 * acceleration * acceleration
            )) / 2.0;
        MatrixXd Z_new = B +
            ((acceleration - 1.0) / acceleration_new) *
            (B - B_previous);
        const double restart =
            ((Z - B).array() * (B - B_previous).array()).sum();
        if (restart > 0.0) {
            acceleration = 1.0;
            Z = B;
        } else {
            acceleration = acceleration_new;
            Z = Z_new;
        }
    }
    if (iteration > max_iter) iteration = max_iter;
    SEXP history_output = R_NilValue;
    if (keep_history) history_output = Rcpp::wrap(history);
    return Rcpp::List::create(
        Rcpp::Named("beta") = B,
        Rcpp::Named("intercept") = smooth_B.intercept,
        Rcpp::Named("lambda") = lambda,
        Rcpp::Named("objective") = objective_previous,
        Rcpp::Named("objective_change") = objective_change,
        Rcpp::Named("coef_change") = coef_change,
        Rcpp::Named("iterations") = iteration,
        Rcpp::Named("converged") = converged,
        Rcpp::Named("step") = step,
        Rcpp::Named("history") = history_output
    );
}

} // namespace

Rcpp::List condition_fit_multitask_path_visit_internal(
    const std::vector<SparseMatrix<double>>& X,
    const std::vector<VectorXd>& y,
    const std::vector<double>& lambda,
    double alpha,
    double condition_mix,
    const std::string& condition_weight,
    const Eigen::ArrayXXi& coefficient_mask,
    int max_iter,
    double tol_objective,
    double tol_coef,
    bool keep_history,
    const ConditionFitVisitor& visitor,
    bool retain_fits
) {
    if (X.size() < 2 || X.size() != y.size()) {
        Rcpp::stop(
            "X and y must contain the same two or more conditions."
        );
    }
    if (!std::isfinite(alpha) || alpha < 0.0 || alpha > 1.0 ||
        !std::isfinite(condition_mix) ||
        condition_mix < 0.0 || condition_mix > 1.0) {
        Rcpp::stop(
            "alpha and condition_mix must be finite values in [0, 1]."
        );
    }
    if (max_iter < 1 || !std::isfinite(tol_objective) ||
        tol_objective <= 0.0 || !std::isfinite(tol_coef) ||
        tol_coef <= 0.0) {
        Rcpp::stop("Invalid solver convergence controls.");
    }
    const int tasks = static_cast<int>(X.size());
    const int p = X[0].cols();
    if (coefficient_mask.rows() != p ||
        coefficient_mask.cols() != tasks) {
        Rcpp::stop(
            "coefficient_mask dimensions do not match predictors and conditions."
        );
    }
    for (int row = 0; row < p; ++row) {
        int eligible = 0;
        for (int task = 0; task < tasks; ++task) {
            const int value = coefficient_mask(row, task);
            if (value != 0 && value != 1) {
                Rcpp::stop("coefficient_mask must contain only zero or one.");
            }
            eligible += value;
        }
        if (!eligible) {
            Rcpp::stop("Every predictor requires an eligible condition.");
        }
    }
    VectorXd loss_weights(tasks);
    VectorXd centered_spectral_norm(tasks);
    double total_n = 0.0;
    for (int task = 0; task < tasks; ++task) {
        if (X[task].cols() != p ||
            X[task].rows() != y[task].size() || X[task].rows() < 1) {
            Rcpp::stop("Condition matrices and responses are not aligned.");
        }
        if (!y[task].allFinite()) {
            Rcpp::stop("Condition responses contain non-finite values.");
        }
        total_n += X[task].rows();
        centered_spectral_norm[task] =
            centered_spectral_squared_norm(X[task]);
    }
    if (condition_weight == "equal") {
        for (int task = 0; task < tasks; ++task) {
            loss_weights[task] = 1.0 / X[task].rows();
        }
    } else if (condition_weight == "cell_count") {
        loss_weights.setConstant(1.0 / total_n);
    } else {
        Rcpp::stop("Unsupported condition_weight.");
    }
    MatrixXd B = MatrixXd::Zero(p, tasks);
    double step = std::numeric_limits<double>::quiet_NaN();
    Rcpp::List fits(retain_fits ? lambda.size() : 0);
    for (int index = 0; index < static_cast<int>(lambda.size()); ++index) {
        const double current_lambda = lambda[index];
        if (!std::isfinite(current_lambda) || current_lambda < 0.0) {
            Rcpp::stop("lambda must contain finite non-negative values.");
        }
        const double ridge = current_lambda * (1.0 - alpha);
        const double default_step = initial_step_for_lambda(
            centered_spectral_norm, loss_weights, ridge
        );
        Rcpp::List fit = fit_one_lambda(
            X,
            y,
            loss_weights,
            coefficient_mask,
            current_lambda,
            alpha,
            condition_mix,
            B,
            step,
            default_step,
            max_iter,
            tol_objective,
            tol_coef,
            keep_history
        );
        fit["backend"] = "cpp_eigen_sparse_matrix_free_fista";
        B = Rcpp::as<MatrixXd>(fit["beta"]);
        step = Rcpp::as<double>(fit["step"]);
        if (visitor) visitor(index, fit);
        if (retain_fits) fits[index] = fit;
    }
    return Rcpp::List::create(
        Rcpp::Named("lambda") = Rcpp::wrap(lambda),
        Rcpp::Named("fits") = fits,
        Rcpp::Named("backend") = "cpp_eigen_sparse_matrix_free_fista"
    );
}

Rcpp::List condition_fit_multitask_path_internal(
    const std::vector<SparseMatrix<double>>& X,
    const std::vector<VectorXd>& y,
    const std::vector<double>& lambda,
    double alpha,
    double condition_mix,
    const std::string& condition_weight,
    const Eigen::ArrayXXi& coefficient_mask,
    int max_iter,
    double tol_objective,
    double tol_coef,
    bool keep_history
) {
    return condition_fit_multitask_path_visit_internal(
        X,
        y,
        lambda,
        alpha,
        condition_mix,
        condition_weight,
        coefficient_mask,
        max_iter,
        tol_objective,
        tol_coef,
        keep_history,
        ConditionFitVisitor(),
        true
    );
}

// [[Rcpp::export]]
Rcpp::List condition_fit_multitask_path_cpp(
    Rcpp::List X_list,
    Rcpp::List y_list,
    Rcpp::NumericVector lambda,
    double alpha,
    double condition_mix,
    std::string condition_weight,
    Rcpp::LogicalMatrix coefficient_mask,
    int max_iter,
    double tol_objective,
    double tol_coef,
    bool keep_history
) {
    if (X_list.size() < 2 || X_list.size() != y_list.size()) {
        Rcpp::stop(
            "X_list and y_list must contain the same two or more conditions."
        );
    }
    if (!std::isfinite(alpha) || alpha < 0.0 || alpha > 1.0 ||
        !std::isfinite(condition_mix) ||
        condition_mix < 0.0 || condition_mix > 1.0) {
        Rcpp::stop(
            "alpha and condition_mix must be finite values in [0, 1]."
        );
    }
    if (max_iter < 1 || !std::isfinite(tol_objective) ||
        tol_objective <= 0.0 || !std::isfinite(tol_coef) ||
        tol_coef <= 0.0) {
        Rcpp::stop("Invalid solver convergence controls.");
    }
    std::vector<SparseMatrix<double>> X = as_sparse_list(X_list);
    std::vector<VectorXd> y = as_vector_list(y_list);
    const int tasks = X.size();
    const int p = X[0].cols();
    if (coefficient_mask.nrow() != p ||
        coefficient_mask.ncol() != tasks) {
        Rcpp::stop(
            "coefficient_mask dimensions do not match predictors and conditions."
        );
    }
    Eigen::ArrayXXi mask(p, tasks);
    for (int row = 0; row < p; ++row) {
        int eligible = 0;
        for (int task = 0; task < tasks; ++task) {
            if (Rcpp::LogicalVector::is_na(
                    coefficient_mask(row, task))) {
                Rcpp::stop("coefficient_mask cannot contain NA.");
            }
            mask(row, task) = coefficient_mask(row, task) ? 1 : 0;
            eligible += mask(row, task);
        }
        if (!eligible) {
            Rcpp::stop("Every predictor requires an eligible condition.");
        }
    }
    return condition_fit_multitask_path_internal(
        X,
        y,
        std::vector<double>(lambda.begin(), lambda.end()),
        alpha,
        condition_mix,
        condition_weight,
        mask,
        max_iter,
        tol_objective,
        tol_coef,
        keep_history
    );
}
