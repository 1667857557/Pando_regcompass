#include <RcppEigen.h>
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
        if (Rf_isS4(value) && Rf_inherits(value, "sparseMatrix")) {
            Eigen::MappedSparseMatrix<double> mapped =
                Rcpp::as<Eigen::MappedSparseMatrix<double>>(value);
            answer.emplace_back(mapped);
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

double centered_squared_norm(const SparseMatrix<double>& X) {
    if (X.rows() == 0 || X.cols() == 0) return 0.0;
    VectorXd ones = VectorXd::Ones(X.rows());
    VectorXd mean =
        (X.transpose() * ones) / static_cast<double>(X.rows());
    return std::max(
        X.squaredNorm() - X.rows() * mean.squaredNorm(),
        0.0
    );
}

double safe_step_for_lambda(
    const VectorXd& centered_norm,
    const VectorXd& loss_weights,
    double ridge
) {
    double upper = ridge;
    for (int task = 0; task < centered_norm.size(); ++task) {
        upper += loss_weights[task] * centered_norm[task];
    }
    return (!std::isfinite(upper) || upper <= 0.0) ?
        1.0 : 1.0 / upper;
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
    double safe_step,
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
    double step = safe_step;
    if (std::isfinite(initial_step) && initial_step > 0.0) {
        step = std::max(initial_step, safe_step);
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
    VectorXd loss_weights(tasks);
    VectorXd centered_norm(tasks);
    double total_n = 0.0;
    for (int task = 0; task < tasks; ++task) {
        if (X[task].cols() != p ||
            X[task].rows() != y[task].size() || X[task].rows() < 1) {
            Rcpp::stop(
                "Condition matrices and responses are not aligned."
            );
        }
        total_n += X[task].rows();
        centered_norm[task] = centered_squared_norm(X[task]);
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
    Rcpp::List fits(lambda.size());
    for (R_xlen_t index = 0; index < lambda.size(); ++index) {
        const double current_lambda = lambda[index];
        if (!std::isfinite(current_lambda) || current_lambda < 0.0) {
            Rcpp::stop(
                "lambda must contain finite non-negative values."
            );
        }
        const double ridge = current_lambda * (1.0 - alpha);
        const double safe_step = safe_step_for_lambda(
            centered_norm, loss_weights, ridge
        );
        Rcpp::List fit = fit_one_lambda(
            X,
            y,
            loss_weights,
            mask,
            current_lambda,
            alpha,
            condition_mix,
            B,
            step,
            safe_step,
            max_iter,
            tol_objective,
            tol_coef,
            keep_history
        );
        B = Rcpp::as<MatrixXd>(fit["beta"]);
        step = Rcpp::as<double>(fit["step"]);
        fits[index] = fit;
    }
    return Rcpp::List::create(
        Rcpp::Named("lambda") = lambda,
        Rcpp::Named("fits") = fits
    );
}
