#include <RcppEigen.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

using Eigen::LLT;
using Eigen::MatrixXd;
using Eigen::VectorXd;

namespace {

const double kJitter = std::sqrt(std::numeric_limits<double>::epsilon());

MatrixXd as_finite_matrix(
    SEXP value, int expected_rows, int expected_cols,
    const std::string& label
) {
    Rcpp::NumericMatrix input(value);
    if (input.nrow() != expected_rows || input.ncol() != expected_cols) {
        Rcpp::stop(label + " has incompatible dimensions.");
    }
    Eigen::Map<MatrixXd> mapped(
        input.begin(), input.nrow(), input.ncol()
    );
    MatrixXd answer(mapped);
    if (!answer.allFinite()) {
        Rcpp::stop(label + " contains non-finite values.");
    }
    return answer;
}

VectorXd as_finite_vector(
    SEXP value, int expected_size, const std::string& label
) {
    Rcpp::NumericVector input(value);
    if (input.size() != expected_size) {
        Rcpp::stop(label + " has incompatible length.");
    }
    Eigen::Map<VectorXd> mapped(input.begin(), input.size());
    VectorXd answer(mapped);
    if (!answer.allFinite()) {
        Rcpp::stop(label + " contains non-finite values.");
    }
    return answer;
}

std::vector<int> true_indices(
    const Rcpp::LogicalMatrix& mask, int column
) {
    std::vector<int> answer;
    answer.reserve(mask.nrow());
    for (int row = 0; row < mask.nrow(); ++row) {
        const int value = mask(row, column);
        if (value == NA_LOGICAL) {
            Rcpp::stop("estimability_mask cannot contain NA values.");
        }
        if (value) answer.push_back(row);
    }
    return answer;
}

MatrixXd subset_matrix(
    const MatrixXd& source,
    const std::vector<int>& rows,
    const std::vector<int>& columns
) {
    MatrixXd answer(rows.size(), columns.size());
    for (int column = 0; column < static_cast<int>(columns.size()); ++column) {
        for (int row = 0; row < static_cast<int>(rows.size()); ++row) {
            answer(row, column) = source(rows[row], columns[column]);
        }
    }
    return answer;
}

VectorXd subset_vector(
    const VectorXd& source, const std::vector<int>& index
) {
    VectorXd answer(index.size());
    for (int position = 0; position < static_cast<int>(index.size()); ++position) {
        answer[position] = source[index[position]];
    }
    return answer;
}

void add_submatrix(
    MatrixXd& destination,
    const std::vector<int>& index,
    const MatrixXd& value,
    double scale
) {
    for (int column = 0; column < static_cast<int>(index.size()); ++column) {
        for (int row = 0; row < static_cast<int>(index.size()); ++row) {
            destination(index[row], index[column]) += scale * value(row, column);
        }
    }
}

void add_subvector(
    VectorXd& destination,
    const std::vector<int>& index,
    const VectorXd& value,
    double scale
) {
    for (int position = 0; position < static_cast<int>(index.size()); ++position) {
        destination[index[position]] += scale * value[position];
    }
}

struct DenseSolve {
    MatrixXd lhs;
    MatrixXd solution;
};

DenseSolve solve_spd(
    MatrixXd lhs, const MatrixXd& rhs, const std::string& label
) {
    if (lhs.rows() != lhs.cols() || lhs.rows() != rhs.rows()) {
        Rcpp::stop(label + " system dimensions are incompatible.");
    }
    if (!lhs.allFinite() || !rhs.allFinite()) {
        Rcpp::stop(label + " system contains non-finite values.");
    }
    lhs = 0.5 * (lhs + lhs.transpose());
    lhs.diagonal().array() += kJitter;
    LLT<MatrixXd> factor(lhs);
    if (factor.info() != Eigen::Success) {
        Rcpp::stop(label + " Cholesky factorization failed.");
    }
    MatrixXd solution = factor.solve(rhs);
    if (factor.info() != Eigen::Success || !solution.allFinite()) {
        Rcpp::stop(label + " Cholesky solve failed.");
    }
    const double denominator =
        rhs.norm() + lhs.norm() * solution.norm() + 1.0;
    const double residual = (lhs * solution - rhs).norm() / denominator;
    if (!std::isfinite(residual) || residual > 1e-8) {
        Rcpp::stop(
            label + " failed the internal linear-system residual check."
        );
    }
    return DenseSolve{lhs, solution};
}

struct TaskCache {
    VectorXd x_mean;
    double y_mean;
    MatrixXd gram;
    VectorXd rhs;
};

struct RefitCache {
    std::vector<TaskCache> task;
    MatrixXd common_metric;
    VectorXd loss_weights;
    VectorXd average_weights;
};

RefitCache parse_cache(const Rcpp::List& cache, int p, int tasks) {
    if (!cache.containsElementNamed("task") ||
        !cache.containsElementNamed("common_metric") ||
        !cache.containsElementNamed("loss_weights") ||
        !cache.containsElementNamed("average_weights")) {
        Rcpp::stop("ConditionRefitCache is missing required fields.");
    }
    Rcpp::List task_input(cache["task"]);
    if (task_input.size() != tasks) {
        Rcpp::stop("ConditionRefitCache task count is not aligned.");
    }
    RefitCache answer{
        std::vector<TaskCache>(),
        as_finite_matrix(cache["common_metric"], p, p, "common_metric"),
        as_finite_vector(cache["loss_weights"], tasks, "loss_weights"),
        as_finite_vector(cache["average_weights"], tasks, "average_weights")
    };
    if ((answer.loss_weights.array() <= 0.0).any() ||
        (answer.average_weights.array() <= 0.0).any()) {
        Rcpp::stop("Refit weights must be finite and positive.");
    }
    answer.task.reserve(tasks);
    for (int task = 0; task < tasks; ++task) {
        Rcpp::List value(task_input[task]);
        if (!value.containsElementNamed("x_mean") ||
            !value.containsElementNamed("y_mean") ||
            !value.containsElementNamed("gram") ||
            !value.containsElementNamed("rhs")) {
            Rcpp::stop("ConditionRefitCache task is missing required fields.");
        }
        const double y_mean = Rcpp::as<double>(value["y_mean"]);
        if (!std::isfinite(y_mean)) {
            Rcpp::stop("ConditionRefitCache y_mean is non-finite.");
        }
        answer.task.push_back(TaskCache{
            as_finite_vector(value["x_mean"], p, "task x_mean"),
            y_mean,
            as_finite_matrix(value["gram"], p, p, "task gram"),
            as_finite_vector(value["rhs"], p, "task rhs")
        });
    }
    return answer;
}

struct TaskPlan {
    bool active;
    std::vector<int> active_index;
    std::vector<int> estimable_index;
    MatrixXd cross;
    VectorXd data_solution;
    MatrixXd cross_solution;
    MatrixXd lhs;
    VectorXd data_rhs;
};

Rcpp::List refit_one(
    const MatrixXd& beta_selection,
    const Rcpp::LogicalMatrix& estimability_mask,
    const std::vector<std::vector<int>>& estimable_index,
    const RefitCache& cache,
    const MatrixXd& shared_lhs,
    double ridge,
    double active_tol,
    double tolerance,
    int path_index
) {
    const int p = beta_selection.rows();
    const int tasks = beta_selection.cols();
    Rcpp::LogicalMatrix support_mask(p, tasks);
    std::vector<TaskPlan> plan(tasks);

    for (int task = 0; task < tasks; ++task) {
        std::vector<int> active;
        active.reserve(p);
        for (int predictor = 0; predictor < p; ++predictor) {
            const bool estimable = estimability_mask(predictor, task);
            const bool selected = estimable &&
                std::abs(beta_selection(predictor, task)) > active_tol;
            support_mask(predictor, task) = selected;
            if (selected) active.push_back(predictor);
        }
        plan[task].active = !active.empty();
        plan[task].active_index = active;
        plan[task].estimable_index = estimable_index[task];
        if (active.empty()) continue;

        const MatrixXd task_gram = subset_matrix(
            cache.task[task].gram, active, active
        );
        const MatrixXd common_active = subset_matrix(
            cache.common_metric, active, active
        );
        MatrixXd lhs = cache.loss_weights[task] * task_gram +
            ridge * cache.average_weights[task] * common_active;
        MatrixXd cross = subset_matrix(
            cache.common_metric, active, estimable_index[task]
        );
        VectorXd data_rhs = cache.loss_weights[task] *
            subset_vector(cache.task[task].rhs, active);
        MatrixXd rhs_block(active.size(), 1 + cross.cols());
        rhs_block.col(0) = data_rhs;
        if (cross.cols() > 0) {
            rhs_block.rightCols(cross.cols()) = cross;
        }
        DenseSolve solved = solve_spd(
            lhs,
            rhs_block,
            "condition refit task system at path index " +
                std::to_string(path_index + 1)
        );
        plan[task].cross = cross;
        plan[task].data_solution = solved.solution.col(0);
        plan[task].cross_solution = solved.solution.rightCols(cross.cols());
        plan[task].lhs = solved.lhs;
        plan[task].data_rhs = data_rhs;
    }

    MatrixXd schur = shared_lhs;
    VectorXd shared_rhs = VectorXd::Zero(p);
    for (int task = 0; task < tasks; ++task) {
        if (!plan[task].active) continue;
        const double weight = cache.average_weights[task];
        const VectorXd contribution =
            plan[task].cross.transpose() * plan[task].data_solution;
        add_subvector(
            shared_rhs,
            plan[task].estimable_index,
            contribution,
            weight
        );
        const MatrixXd correction =
            plan[task].cross.transpose() * plan[task].cross_solution;
        add_submatrix(
            schur,
            plan[task].estimable_index,
            correction,
            -ridge * weight * weight
        );
    }
    schur = 0.5 * (schur + schur.transpose());
    MatrixXd shared_rhs_matrix(shared_rhs.size(), 1);
    shared_rhs_matrix.col(0) = shared_rhs;
    DenseSolve shared_solve = solve_spd(
        schur,
        shared_rhs_matrix,
        "condition refit shared Schur system at path index " +
            std::to_string(path_index + 1)
    );
    VectorXd shared = shared_solve.solution.col(0);

    MatrixXd beta = MatrixXd::Zero(p, tasks);
    VectorXd intercept = VectorXd::Zero(tasks);
    VectorXd task_residual = VectorXd::Zero(tasks);
    for (int task = 0; task < tasks; ++task) {
        if (!plan[task].active) {
            intercept[task] = cache.task[task].y_mean;
            continue;
        }
        const VectorXd shared_estimable = subset_vector(
            shared, plan[task].estimable_index
        );
        const VectorXd coefficient = plan[task].data_solution +
            ridge * cache.average_weights[task] *
            (plan[task].cross_solution * shared_estimable);
        for (int position = 0;
             position < static_cast<int>(plan[task].active_index.size());
             ++position) {
            beta(plan[task].active_index[position], task) = coefficient[position];
        }
        const VectorXd x_mean_active = subset_vector(
            cache.task[task].x_mean, plan[task].active_index
        );
        intercept[task] = cache.task[task].y_mean -
            x_mean_active.dot(coefficient);
        const VectorXd lhs_value = plan[task].lhs * coefficient;
        const VectorXd rhs_value = plan[task].data_rhs +
            ridge * cache.average_weights[task] *
            (plan[task].cross * shared_estimable);
        task_residual[task] = (lhs_value - rhs_value).norm() /
            (rhs_value.norm() + 1.0);
    }

    VectorXd shared_residual = shared_lhs * shared;
    for (int task = 0; task < tasks; ++task) {
        if (estimable_index[task].empty()) continue;
        VectorXd product(estimable_index[task].size());
        for (int row = 0;
             row < static_cast<int>(estimable_index[task].size()); ++row) {
            product[row] = cache.common_metric.row(
                estimable_index[task][row]
            ).dot(beta.col(task));
        }
        add_subvector(
            shared_residual,
            estimable_index[task],
            product,
            -cache.average_weights[task]
        );
    }
    double coef_change = task_residual.maxCoeff();
    coef_change = std::max(
        coef_change,
        shared_residual.norm() / (shared_rhs.norm() + 1.0)
    );
    if (!std::isfinite(coef_change)) {
        Rcpp::stop(
            "Compiled condition refit produced a non-finite residual."
        );
    }
    const bool converged = coef_change < std::max(
        tolerance,
        100.0 * std::numeric_limits<double>::epsilon()
    );
    return Rcpp::List::create(
        Rcpp::Named("beta") = beta,
        Rcpp::Named("shared") = shared,
        Rcpp::Named("intercept") = intercept,
        Rcpp::Named("support_mask") = support_mask,
        Rcpp::Named("iterations") = 1,
        Rcpp::Named("coef_change") = coef_change,
        Rcpp::Named("converged") = converged,
        Rcpp::Named("backend") = "cpp_eigen_dense_llt_batched_path"
    );
}

} // namespace

// [[Rcpp::export]]
Rcpp::List condition_refit_path_cpp(
    Rcpp::List beta_path,
    Rcpp::LogicalMatrix estimability_mask,
    Rcpp::NumericVector ridge,
    double active_tol,
    Rcpp::List cache,
    double tolerance
) {
    if (beta_path.size() < 1 || beta_path.size() != ridge.size()) {
        Rcpp::stop(
            "beta_path and ridge must contain the same positive number of entries."
        );
    }
    if (!std::isfinite(active_tol) || active_tol < 0.0 ||
        !std::isfinite(tolerance) || tolerance <= 0.0) {
        Rcpp::stop("Invalid condition refit numerical controls.");
    }
    const int p = estimability_mask.nrow();
    const int tasks = estimability_mask.ncol();
    if (p < 1 || tasks < 2) {
        Rcpp::stop(
            "condition refit requires predictors and at least two conditions."
        );
    }
    std::vector<std::vector<int>> estimable_index(tasks);
    for (int task = 0; task < tasks; ++task) {
        estimable_index[task] = true_indices(estimability_mask, task);
    }
    for (int predictor = 0; predictor < p; ++predictor) {
        bool estimable = false;
        for (int task = 0; task < tasks; ++task) {
            estimable = estimable || estimability_mask(predictor, task);
        }
        if (!estimable) {
            Rcpp::stop(
                "Every refit predictor requires an estimable condition."
            );
        }
    }

    RefitCache parsed = parse_cache(cache, p, tasks);
    MatrixXd shared_lhs = MatrixXd::Zero(p, p);
    for (int task = 0; task < tasks; ++task) {
        if (estimable_index[task].empty()) continue;
        const MatrixXd task_metric = subset_matrix(
            parsed.common_metric,
            estimable_index[task],
            estimable_index[task]
        );
        add_submatrix(
            shared_lhs,
            estimable_index[task],
            task_metric,
            parsed.average_weights[task]
        );
    }

    Rcpp::List answer(beta_path.size());
    for (R_xlen_t index = 0; index < beta_path.size(); ++index) {
        const double current_ridge = ridge[index];
        if (!std::isfinite(current_ridge) || current_ridge <= 0.0) {
            Rcpp::stop("ridge must contain finite positive values.");
        }
        MatrixXd beta = as_finite_matrix(
            beta_path[index], p, tasks,
            "beta_path coefficient matrix"
        );
        answer[index] = refit_one(
            beta,
            estimability_mask,
            estimable_index,
            parsed,
            shared_lhs,
            current_ridge,
            active_tol,
            tolerance,
            static_cast<int>(index)
        );
    }
    return answer;
}
