#include "condition_matrix_free_refit.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

using Eigen::ArrayXXi;
using Eigen::LLT;
using Eigen::MatrixXd;
using Eigen::SparseMatrix;
using Eigen::Triplet;
using Eigen::VectorXd;

namespace {

constexpr double kMachineEps = std::numeric_limits<double>::epsilon();
const double kJitter = std::sqrt(kMachineEps);

struct SparseTaskStats {
    VectorXd x_mean;
    double y_mean;
    VectorXd rhs;
    VectorXd gram_diagonal;
    double loss_weight;
    double average_weight;
};

struct SparseRefitProblem {
    const std::vector<SparseMatrix<double>>* X;
    std::vector<SparseTaskStats> task;
    VectorXd common_diagonal;
    int predictors;
    int tasks;
    std::size_t nonzeros;
};

struct PCGResult {
    VectorXd solution;
    int iterations;
    double relative_residual;
    bool converged;
};

struct TaskPlan {
    int task_index;
    double ridge;
    bool active;
    bool dense_factor;
    std::vector<int> active_index;
    std::vector<int> estimable_index;
    MatrixXd lhs;
    LLT<MatrixXd> factor;
    VectorXd diagonal;
    VectorXd data_rhs;
    VectorXd data_solution;
    int data_iterations;
};

VectorXd subset_vector(
    const VectorXd& source,
    const std::vector<int>& index
) {
    VectorXd out(index.size());
    for (int position = 0; position < static_cast<int>(index.size()); ++position) {
        out[position] = source[index[position]];
    }
    return out;
}

void assign_subset(
    VectorXd& destination,
    const std::vector<int>& index,
    const VectorXd& value
) {
    if (index.size() != static_cast<std::size_t>(value.size())) {
        Rcpp::stop("Matrix-free refit subset assignment is not aligned.");
    }
    for (int position = 0; position < static_cast<int>(index.size()); ++position) {
        destination[index[position]] = value[position];
    }
}

void add_subset(
    VectorXd& destination,
    const std::vector<int>& index,
    const VectorXd& value,
    double scale
) {
    if (index.size() != static_cast<std::size_t>(value.size())) {
        Rcpp::stop("Matrix-free refit subset accumulation is not aligned.");
    }
    for (int position = 0; position < static_cast<int>(index.size()); ++position) {
        destination[index[position]] += scale * value[position];
    }
}

VectorXd apply_centered_gram(
    const SparseMatrix<double>& X,
    const VectorXd& mean,
    const VectorXd& value
) {
    if (X.cols() != value.size() || mean.size() != value.size()) {
        Rcpp::stop("Centered sparse Gram operator dimensions are invalid.");
    }
    VectorXd projected = X * value;
    VectorXd out = X.transpose() * projected;
    out.noalias() -= static_cast<double>(X.rows()) *
        mean * mean.dot(value);
    if (!out.allFinite()) {
        Rcpp::stop("Centered sparse Gram operator produced non-finite values.");
    }
    return out;
}

SparseRefitProblem make_problem(
    const std::vector<SparseMatrix<double>>& X,
    const std::vector<VectorXd>& y
) {
    if (X.size() < 2 || X.size() != y.size()) {
        Rcpp::stop("Matrix-free refit requires two or more aligned tasks.");
    }
    SparseRefitProblem out;
    out.X = &X;
    out.tasks = static_cast<int>(X.size());
    out.predictors = X[0].cols();
    out.nonzeros = 0;
    out.common_diagonal = VectorXd::Constant(out.predictors, kJitter);
    out.task.reserve(out.tasks);
    const double average = 1.0 / static_cast<double>(out.tasks);
    for (int task = 0; task < out.tasks; ++task) {
        if (X[task].cols() != out.predictors ||
            X[task].rows() != y[task].size() || X[task].rows() < 1 ||
            !y[task].allFinite()) {
            Rcpp::stop("Matrix-free refit task inputs are not aligned.");
        }
        out.nonzeros += static_cast<std::size_t>(X[task].nonZeros());
        const double n = static_cast<double>(X[task].rows());
        VectorXd sum = VectorXd::Zero(out.predictors);
        VectorXd square = VectorXd::Zero(out.predictors);
        for (int column = 0; column < X[task].outerSize(); ++column) {
            for (SparseMatrix<double>::InnerIterator it(X[task], column); it; ++it) {
                sum[column] += it.value();
                square[column] += it.value() * it.value();
            }
        }
        VectorXd mean = sum / n;
        VectorXd diagonal = square - n * mean.array().square().matrix();
        diagonal = diagonal.cwiseMax(0.0);
        const double y_mean = y[task].mean();
        VectorXd rhs = X[task].transpose() * y[task] -
            n * mean * y_mean;
        const double loss_weight = 1.0 / n;
        out.common_diagonal.noalias() += loss_weight * diagonal;
        out.task.push_back(SparseTaskStats{
            std::move(mean),
            y_mean,
            std::move(rhs),
            std::move(diagonal),
            loss_weight,
            average
        });
    }
    if (!out.common_diagonal.allFinite() ||
        (out.common_diagonal.array() < 0.0).any()) {
        Rcpp::stop("Matrix-free common-metric diagonal is invalid.");
    }
    return out;
}

long double matrix_free_base_bytes(
    const SparseRefitProblem& problem,
    const TargetEngineControl& control
) {
    const long double sparse_resident = 6.0L *
        static_cast<long double>(problem.nonzeros) *
        (sizeof(double) + 2.0L * sizeof(int));
    const long double coefficient_resident =
        static_cast<long double>(problem.predictors) *
        (static_cast<long double>(problem.tasks) *
             (14.0L + control.lambda_batch_size) + 20.0L) *
        sizeof(double);
    return sparse_resident + coefficient_resident;
}

long double dense_workspace_budget(
    const SparseRefitProblem& problem,
    const TargetEngineControl& control
) {
    const long double available = std::max(
        0.0L,
        static_cast<long double>(control.worker_budget_bytes) -
            matrix_free_base_bytes(problem, control)
    );
    return 0.5L * available;
}

VectorXd apply_common_metric(
    const SparseRefitProblem& problem,
    const VectorXd& value
) {
    if (value.size() != problem.predictors) {
        Rcpp::stop("Common-metric operator input is not aligned.");
    }
    VectorXd out = kJitter * value;
    for (int task = 0; task < problem.tasks; ++task) {
        out.noalias() += problem.task[task].loss_weight *
            apply_centered_gram(
                (*problem.X)[task], problem.task[task].x_mean, value
            );
    }
    if (!out.allFinite()) {
        Rcpp::stop("Common-metric operator produced non-finite values.");
    }
    return out;
}

SparseMatrix<double> subset_columns(
    const SparseMatrix<double>& source,
    const std::vector<int>& columns
) {
    std::vector<Triplet<double>> triplets;
    std::size_t reserve = 0;
    for (int column : columns) {
        if (column < 0 || column >= source.cols()) {
            Rcpp::stop("Matrix-free refit support index is out of bounds.");
        }
        reserve += static_cast<std::size_t>(
            source.outerIndexPtr()[column + 1] - source.outerIndexPtr()[column]
        );
    }
    triplets.reserve(reserve);
    for (int output = 0; output < static_cast<int>(columns.size()); ++output) {
        for (SparseMatrix<double>::InnerIterator it(source, columns[output]); it; ++it) {
            triplets.emplace_back(it.row(), output, it.value());
        }
    }
    SparseMatrix<double> out(source.rows(), columns.size());
    out.setFromTriplets(triplets.begin(), triplets.end());
    out.makeCompressed();
    return out;
}

MatrixXd centered_gram_subset(
    const SparseMatrix<double>& X,
    const VectorXd& mean,
    const std::vector<int>& index
) {
    SparseMatrix<double> selected = subset_columns(X, index);
    MatrixXd out = MatrixXd(selected.transpose() * selected);
    VectorXd selected_mean = subset_vector(mean, index);
    out.noalias() -= static_cast<double>(X.rows()) *
        selected_mean * selected_mean.transpose();
    out = (0.5 * (out + out.transpose())).eval();
    if (!out.allFinite()) {
        Rcpp::stop("Active-support Gram matrix contains non-finite values.");
    }
    return out;
}

MatrixXd common_metric_subset(
    const SparseRefitProblem& problem,
    const std::vector<int>& index
) {
    MatrixXd out = kJitter * MatrixXd::Identity(
        index.size(), index.size()
    );
    for (int task = 0; task < problem.tasks; ++task) {
        out.noalias() += problem.task[task].loss_weight *
            centered_gram_subset(
                (*problem.X)[task], problem.task[task].x_mean, index
            );
    }
    return out;
}

PCGResult pcg_solve_preconditioned(
    const std::function<VectorXd(const VectorXd&)>& apply,
    const std::function<VectorXd(const VectorXd&)>& precondition,
    const VectorXd& rhs,
    const VectorXd& initial,
    double tolerance,
    int max_iterations,
    const std::string& label
) {
    if (initial.size() != rhs.size() || !rhs.allFinite() ||
        !initial.allFinite()) {
        Rcpp::stop(label + " PCG inputs are invalid.");
    }
    VectorXd x = initial;
    VectorXd residual = rhs - apply(x);
    const double denominator = rhs.norm() + 1.0;
    double relative = residual.norm() / denominator;
    if (!std::isfinite(relative)) {
        Rcpp::stop(label + " PCG initial residual is non-finite.");
    }
    if (relative <= tolerance) {
        return PCGResult{std::move(x), 0, relative, true};
    }
    VectorXd preconditioned = precondition(residual);
    if (preconditioned.size() != residual.size() ||
        !preconditioned.allFinite()) {
        Rcpp::stop(label + " PCG preconditioner produced invalid values.");
    }
    VectorXd direction = preconditioned;
    double residual_product = residual.dot(preconditioned);
    if (!std::isfinite(residual_product) || residual_product <= 0.0) {
        Rcpp::stop(label + " PCG preconditioner is not positive.");
    }
    int iteration = 0;
    for (iteration = 1; iteration <= max_iterations; ++iteration) {
        VectorXd applied = apply(direction);
        const double curvature = direction.dot(applied);
        if (!std::isfinite(curvature) || curvature <= 0.0) {
            Rcpp::stop(label + " PCG encountered non-positive curvature.");
        }
        const double step = residual_product / curvature;
        x.noalias() += step * direction;
        residual.noalias() -= step * applied;
        relative = residual.norm() / denominator;
        if (!std::isfinite(relative) || !x.allFinite()) {
            Rcpp::stop(label + " PCG became non-finite.");
        }
        if (relative <= tolerance) {
            return PCGResult{std::move(x), iteration, relative, true};
        }
        preconditioned = precondition(residual);
        if (preconditioned.size() != residual.size() ||
            !preconditioned.allFinite()) {
            Rcpp::stop(
                label + " PCG preconditioner produced invalid values."
            );
        }
        const double next_product = residual.dot(preconditioned);
        if (!std::isfinite(next_product) || next_product <= 0.0) {
            Rcpp::stop(label + " PCG residual product is not positive.");
        }
        const double update = next_product / residual_product;
        direction = preconditioned + update * direction;
        residual_product = next_product;
    }
    return PCGResult{
        std::move(x), max_iterations, relative, false
    };
}

PCGResult pcg_solve(
    const std::function<VectorXd(const VectorXd&)>& apply,
    const VectorXd& diagonal,
    const VectorXd& rhs,
    const VectorXd& initial,
    double tolerance,
    int max_iterations,
    const std::string& label
) {
    if (diagonal.size() != rhs.size() || !diagonal.allFinite() ||
        (diagonal.array() <= 0.0).any()) {
        Rcpp::stop(label + " PCG diagonal preconditioner is invalid.");
    }
    return pcg_solve_preconditioned(
        apply,
        [&](const VectorXd& residual) {
            return residual.cwiseQuotient(diagonal);
        },
        rhs,
        initial,
        tolerance,
        max_iterations,
        label
    );
}

VectorXd apply_task_lhs(
    const SparseRefitProblem& problem,
    const TaskPlan& plan,
    const VectorXd& value
) {
    VectorXd full = VectorXd::Zero(problem.predictors);
    assign_subset(full, plan.active_index, value);
    VectorXd out = problem.task[plan.task_index].loss_weight *
        apply_centered_gram(
            (*problem.X)[plan.task_index],
            problem.task[plan.task_index].x_mean,
            full
        );
    out.noalias() += plan.ridge *
        problem.task[plan.task_index].average_weight *
        apply_common_metric(problem, full);
    out.noalias() += kJitter * full;
    return subset_vector(out, plan.active_index);
}

PCGResult solve_task(
    const SparseRefitProblem& problem,
    const TaskPlan& plan,
    const VectorXd& rhs,
    const TargetEngineControl& control,
    const std::string& label
) {
    if (plan.dense_factor) {
        VectorXd solution = plan.factor.solve(rhs);
        if (plan.factor.info() != Eigen::Success || !solution.allFinite()) {
            Rcpp::stop(label + " active-support LLT solve failed.");
        }
        const double relative = (plan.lhs * solution - rhs).norm() /
            (rhs.norm() + plan.lhs.norm() * solution.norm() + 1.0);
        return PCGResult{std::move(solution), 1, relative, relative <= 1e-8};
    }
    const double tolerance = std::min(
        1e-10, 0.1 * control.refit_pcg_tolerance
    );
    return pcg_solve(
        [&](const VectorXd& value) {
            return apply_task_lhs(problem, plan, value);
        },
        plan.diagonal,
        rhs,
        VectorXd::Zero(rhs.size()),
        tolerance,
        control.refit_pcg_max_iterations,
        label + " active-support"
    );
}

std::vector<int> true_indices(const ArrayXXi& mask, int column) {
    std::vector<int> out;
    out.reserve(mask.rows());
    for (int row = 0; row < mask.rows(); ++row) {
        if (mask(row, column)) out.push_back(row);
    }
    return out;
}

TaskPlan build_task_plan(
    const SparseRefitProblem& problem,
    const MatrixXd& beta_selection,
    const ArrayXXi& estimability_mask,
    int task,
    double ridge,
    double active_tol,
    const TargetEngineControl& control,
    int path_index,
    long double dense_budget,
    long double* dense_persistent,
    long double* dense_peak
) {
    TaskPlan out;
    out.task_index = task;
    out.ridge = ridge;
    out.active = false;
    out.dense_factor = false;
    out.estimable_index = true_indices(estimability_mask, task);
    out.data_iterations = 0;
    for (int predictor = 0; predictor < problem.predictors; ++predictor) {
        if (estimability_mask(predictor, task) &&
            std::abs(beta_selection(predictor, task)) > active_tol) {
            out.active_index.push_back(predictor);
        }
    }
    if (out.active_index.empty()) return out;
    out.active = true;
    const std::size_t support = out.active_index.size();
    const long double dense_construction_bytes = 6.0L * sizeof(double) *
        static_cast<long double>(support) * support;
    const long double dense_retained_bytes = 2.0L * sizeof(double) *
        static_cast<long double>(support) * support;
    out.dense_factor = *dense_persistent + dense_construction_bytes <=
        dense_budget;
    *dense_peak = std::max(
        *dense_peak,
        *dense_persistent +
            (out.dense_factor ? dense_construction_bytes : 0.0L)
    );
    const SparseTaskStats& stats = problem.task[task];
    out.data_rhs = stats.loss_weight *
        subset_vector(stats.rhs, out.active_index);
    out.diagonal = stats.loss_weight *
        subset_vector(stats.gram_diagonal, out.active_index) +
        ridge * stats.average_weight *
        subset_vector(problem.common_diagonal, out.active_index);
    out.diagonal.array() += kJitter;
    if (out.dense_factor) {
        out.lhs = stats.loss_weight * centered_gram_subset(
            (*problem.X)[task], stats.x_mean, out.active_index
        ) + ridge * stats.average_weight *
            common_metric_subset(problem, out.active_index);
        out.lhs = (0.5 * (out.lhs + out.lhs.transpose())).eval();
        out.lhs.diagonal().array() += kJitter;
        out.factor.compute(out.lhs);
        if (out.factor.info() != Eigen::Success) {
            Rcpp::stop(
                "Matrix-free refit active-support LLT failed at path index " +
                std::to_string(path_index + 1) + "."
            );
        }
        *dense_persistent += dense_retained_bytes;
    }
    PCGResult solved = solve_task(
        problem,
        out,
        out.data_rhs,
        control,
        "Matrix-free refit path index " + std::to_string(path_index + 1)
    );
    if (!solved.converged) {
        Rcpp::stop(
            "Matrix-free refit active-support solve did not converge at path index " +
            std::to_string(path_index + 1) + "."
        );
    }
    out.data_solution = std::move(solved.solution);
    out.data_iterations = solved.iterations;
    return out;
}

VectorXd apply_shared_lhs(
    const SparseRefitProblem& problem,
    const std::vector<TaskPlan>& plans,
    const VectorXd& value
) {
    VectorXd out = VectorXd::Zero(problem.predictors);
    for (const TaskPlan& plan : plans) {
        if (plan.estimable_index.empty()) continue;
        VectorXd masked = VectorXd::Zero(problem.predictors);
        assign_subset(
            masked,
            plan.estimable_index,
            subset_vector(value, plan.estimable_index)
        );
        VectorXd applied = apply_common_metric(problem, masked);
        add_subset(
            out,
            plan.estimable_index,
            subset_vector(applied, plan.estimable_index),
            problem.task[plan.task_index].average_weight
        );
    }
    return out;
}

VectorXd apply_schur(
    const SparseRefitProblem& problem,
    const std::vector<TaskPlan>& plans,
    const VectorXd& value,
    const TargetEngineControl& control,
    int path_index,
    int* maximum_inner_iterations
) {
    VectorXd out = apply_shared_lhs(problem, plans, value);
    for (const TaskPlan& plan : plans) {
        if (!plan.active) continue;
        VectorXd estimable = VectorXd::Zero(problem.predictors);
        assign_subset(
            estimable,
            plan.estimable_index,
            subset_vector(value, plan.estimable_index)
        );
        VectorXd common = apply_common_metric(problem, estimable);
        VectorXd cross_value = subset_vector(common, plan.active_index);
        PCGResult solved = solve_task(
            problem,
            plan,
            cross_value,
            control,
            "Matrix-free Schur path index " + std::to_string(path_index + 1)
        );
        if (!solved.converged) {
            Rcpp::stop(
                "Matrix-free active block did not converge during Schur application."
            );
        }
        *maximum_inner_iterations = std::max(
            *maximum_inner_iterations, solved.iterations
        );
        VectorXd active = VectorXd::Zero(problem.predictors);
        assign_subset(active, plan.active_index, solved.solution);
        VectorXd back = apply_common_metric(problem, active);
        const double average =
            problem.task[plan.task_index].average_weight;
        add_subset(
            out,
            plan.estimable_index,
            subset_vector(back, plan.estimable_index),
            -plan.ridge * average * average
        );
    }
    out.noalias() += kJitter * value;
    return out;
}

Rcpp::List refit_one(
    const SparseRefitProblem& problem,
    const MatrixXd& beta_selection,
    const ArrayXXi& estimability_mask,
    double ridge,
    double active_tol,
    const TargetEngineControl& control,
    int path_index,
    const VectorXd& initial_shared
) {
    std::vector<TaskPlan> plans;
    plans.reserve(problem.tasks);
    Rcpp::LogicalMatrix support_mask(problem.predictors, problem.tasks);
    int maximum_active = 0;
    const long double dense_budget = dense_workspace_budget(
        problem, control
    );
    long double dense_persistent = 0.0L;
    long double dense_peak = 0.0L;
    for (int task = 0; task < problem.tasks; ++task) {
        plans.emplace_back(build_task_plan(
            problem,
            beta_selection,
            estimability_mask,
            task,
            ridge,
            active_tol,
            control,
            path_index,
            dense_budget,
            &dense_persistent,
            &dense_peak
        ));
        maximum_active = std::max(
            maximum_active,
            static_cast<int>(plans.back().active_index.size())
        );
        for (int predictor : plans.back().active_index) {
            support_mask(predictor, task) = true;
        }
    }
    VectorXd shared_rhs = VectorXd::Zero(problem.predictors);
    for (const TaskPlan& plan : plans) {
        if (!plan.active) continue;
        VectorXd active = VectorXd::Zero(problem.predictors);
        assign_subset(active, plan.active_index, plan.data_solution);
        VectorXd contribution = apply_common_metric(problem, active);
        add_subset(
            shared_rhs,
            plan.estimable_index,
            subset_vector(contribution, plan.estimable_index),
            problem.task[plan.task_index].average_weight
        );
    }
    VectorXd shared_diagonal = VectorXd::Constant(
        problem.predictors, kJitter
    );
    for (const TaskPlan& plan : plans) {
        const double average =
            problem.task[plan.task_index].average_weight;
        for (int predictor : plan.estimable_index) {
            shared_diagonal[predictor] +=
                average * problem.common_diagonal[predictor];
        }
    }
    std::vector<int> active_union;
    for (const TaskPlan& plan : plans) {
        active_union.insert(
            active_union.end(),
            plan.active_index.begin(),
            plan.active_index.end()
        );
    }
    std::sort(active_union.begin(), active_union.end());
    active_union.erase(
        std::unique(active_union.begin(), active_union.end()),
        active_union.end()
    );
    int maximum_inner_iterations = 0;
    bool hybrid_preconditioner = false;
    LLT<MatrixXd> hybrid_factor;
    const std::size_t hybrid_active_limit = static_cast<std::size_t>(
        std::min(control.dense_max_p, 512)
    );
    if (control.preconditioner == "hybrid" && !active_union.empty() &&
        active_union.size() <= hybrid_active_limit) {
        const long double active_size =
            static_cast<long double>(active_union.size());
        const long double hybrid_construction_bytes = 8.0L *
            sizeof(double) * active_size * active_size;
        if (dense_persistent + hybrid_construction_bytes <= dense_budget) {
            dense_peak = std::max(
                dense_peak,
                dense_persistent + hybrid_construction_bytes
            );
            MatrixXd hybrid_block(
                active_union.size(), active_union.size()
            );
            hybrid_block.setZero();
            const MatrixXd common_union = common_metric_subset(
                problem, active_union
            );
            for (int task = 0; task < problem.tasks; ++task) {
                const double average =
                    problem.task[task].average_weight;
                for (int column = 0;
                     column < static_cast<int>(active_union.size());
                     ++column) {
                    if (!estimability_mask(active_union[column], task)) {
                        continue;
                    }
                    for (int row = 0;
                         row < static_cast<int>(active_union.size());
                         ++row) {
                        if (estimability_mask(active_union[row], task)) {
                            hybrid_block(row, column) +=
                                average * common_union(row, column);
                        }
                    }
                }
            }
            hybrid_block.diagonal().array() += kJitter;
            hybrid_block =
                (0.5 * (hybrid_block + hybrid_block.transpose())).eval();
            hybrid_factor.compute(hybrid_block);
            hybrid_preconditioner = hybrid_factor.info() == Eigen::Success;
        }
    }
    auto apply_preconditioner = [&](const VectorXd& residual) {
        VectorXd out = residual.cwiseQuotient(shared_diagonal);
        if (hybrid_preconditioner) {
            VectorXd solved = hybrid_factor.solve(
                subset_vector(residual, active_union)
            );
            if (hybrid_factor.info() != Eigen::Success ||
                !solved.allFinite()) {
                Rcpp::stop(
                    "Matrix-free hybrid preconditioner solve failed."
                );
            }
            assign_subset(out, active_union, solved);
        }
        return out;
    };
    PCGResult shared_solve = pcg_solve_preconditioned(
        [&](const VectorXd& value) {
            return apply_schur(
                problem,
                plans,
                value,
                control,
                path_index,
                &maximum_inner_iterations
            );
        },
        apply_preconditioner,
        shared_rhs,
        initial_shared,
        control.refit_pcg_tolerance,
        control.refit_pcg_max_iterations,
        "Matrix-free shared Schur path index " +
            std::to_string(path_index + 1)
    );
    if (!shared_solve.converged) {
        std::ostringstream message;
        message << "Matrix-free shared Schur PCG did not converge at path index "
                << path_index + 1 << " (p=" << problem.predictors
                << ", max_support=" << maximum_active
                << ", iterations=" << shared_solve.iterations
                << ", relative_residual="
                << shared_solve.relative_residual << ").";
        Rcpp::stop(message.str());
    }
    const VectorXd shared = shared_solve.solution;
    MatrixXd beta = MatrixXd::Zero(problem.predictors, problem.tasks);
    VectorXd intercept = VectorXd::Zero(problem.tasks);
    VectorXd task_residual = VectorXd::Zero(problem.tasks);
    int maximum_task_iterations = 0;
    bool any_iterative_active = false;
    for (const TaskPlan& plan : plans) {
        const int task = plan.task_index;
        if (!plan.active) {
            intercept[task] = problem.task[task].y_mean;
            continue;
        }
        VectorXd estimable = VectorXd::Zero(problem.predictors);
        assign_subset(
            estimable,
            plan.estimable_index,
            subset_vector(shared, plan.estimable_index)
        );
        VectorXd cross = subset_vector(
            apply_common_metric(problem, estimable), plan.active_index
        );
        PCGResult correction = solve_task(
            problem,
            plan,
            cross,
            control,
            "Matrix-free final active solve path index " +
                std::to_string(path_index + 1)
        );
        if (!correction.converged) {
            Rcpp::stop("Matrix-free final active-support solve did not converge.");
        }
        const double average = problem.task[task].average_weight;
        VectorXd coefficient = plan.data_solution +
            ridge * average * correction.solution;
        for (int position = 0;
             position < static_cast<int>(plan.active_index.size());
             ++position) {
            beta(plan.active_index[position], task) = coefficient[position];
        }
        intercept[task] = problem.task[task].y_mean -
            subset_vector(
                problem.task[task].x_mean, plan.active_index
            ).dot(coefficient);
        VectorXd expected = plan.data_rhs + ridge * average * cross;
        VectorXd observed = plan.dense_factor ?
            plan.lhs * coefficient :
            apply_task_lhs(problem, plan, coefficient);
        task_residual[task] = (observed - expected).norm() /
            (expected.norm() + 1.0);
        maximum_task_iterations = std::max(
            maximum_task_iterations,
            std::max(plan.data_iterations, correction.iterations)
        );
        any_iterative_active = any_iterative_active || !plan.dense_factor;
    }
    VectorXd shared_residual = apply_shared_lhs(problem, plans, shared);
    for (const TaskPlan& plan : plans) {
        if (plan.estimable_index.empty()) continue;
        VectorXd applied = apply_common_metric(
            problem, beta.col(plan.task_index)
        );
        add_subset(
            shared_residual,
            plan.estimable_index,
            subset_vector(applied, plan.estimable_index),
            -problem.task[plan.task_index].average_weight
        );
    }
    double coefficient_residual = task_residual.maxCoeff();
    coefficient_residual = std::max(
        coefficient_residual,
        shared_residual.norm() / (shared_rhs.norm() + 1.0)
    );
    if (!std::isfinite(coefficient_residual)) {
        Rcpp::stop("Matrix-free refit produced a non-finite residual.");
    }
    const bool converged = coefficient_residual < std::max(
        control.refit_pcg_tolerance,
        100.0 * kMachineEps
    );
    if (!converged) {
        std::ostringstream message;
        message << "Matrix-free refit failed coefficient residual verification "
                << "at path index " << path_index + 1
                << " (residual=" << coefficient_residual << ").";
        Rcpp::stop(message.str());
    }
    const long double estimated_peak =
        matrix_free_base_bytes(problem, control) + dense_peak;
    if (estimated_peak >
        static_cast<long double>(control.worker_budget_bytes)) {
        Rcpp::stop(
            "Matrix-free refit exceeded its internal worker-memory budget."
        );
    }
    return Rcpp::List::create(
        Rcpp::Named("beta") = beta,
        Rcpp::Named("shared") = shared,
        Rcpp::Named("intercept") = intercept,
        Rcpp::Named("support_mask") = support_mask,
        Rcpp::Named("iterations") = shared_solve.iterations,
        Rcpp::Named("coef_change") = coefficient_residual,
        Rcpp::Named("converged") = true,
        Rcpp::Named("backend") =
            "cpp_eigen_matrix_free_schur_pcg_exact_active_refit",
        Rcpp::Named("pcg_iterations") = shared_solve.iterations,
        Rcpp::Named("pcg_relative_residual") =
            shared_solve.relative_residual,
        Rcpp::Named("active_solve_max_iterations") =
            maximum_task_iterations,
        Rcpp::Named("active_solver") = any_iterative_active ?
            "matrix_free_pcg" : "dense_support_llt",
        Rcpp::Named("maximum_support_size") = maximum_active,
        Rcpp::Named("preconditioner") = hybrid_preconditioner ?
            "active_union_block_plus_diagonal" : "diagonal",
        Rcpp::Named("preconditioner_active_size") =
            hybrid_preconditioner ?
                static_cast<int>(active_union.size()) : 0,
        Rcpp::Named("dense_workspace_peak_bytes") =
            static_cast<double>(dense_peak),
        Rcpp::Named("estimated_peak_bytes") =
            static_cast<double>(estimated_peak),
        Rcpp::Named("budget_guard_passed") = true
    );
}

} // namespace

struct ConditionMatrixFreeRefitWorkspace::Impl {
    SparseRefitProblem problem;
    ArrayXXi estimability_mask;
    double active_tol;
    TargetEngineControl control;
    VectorXd warm_shared;

    Impl(
        const std::vector<SparseMatrix<double>>& X,
        const std::vector<VectorXd>& y,
        const ArrayXXi& mask,
        double tolerance,
        const TargetEngineControl& engine_control
    ) :
        problem(make_problem(X, y)),
        estimability_mask(mask),
        active_tol(tolerance),
        control(engine_control),
        warm_shared(VectorXd::Zero(problem.predictors)) {
        if (!std::isfinite(active_tol) || active_tol < 0.0) {
            Rcpp::stop("Invalid matrix-free active tolerance.");
        }
        if (estimability_mask.rows() != problem.predictors ||
            estimability_mask.cols() != problem.tasks) {
            Rcpp::stop("Matrix-free estimability mask is not aligned.");
        }
        for (int predictor = 0; predictor < problem.predictors; ++predictor) {
            if (!estimability_mask.row(predictor).any()) {
                Rcpp::stop(
                    "Every refit predictor requires an estimable condition."
                );
            }
        }
    }
};

ConditionMatrixFreeRefitWorkspace::ConditionMatrixFreeRefitWorkspace(
    const std::vector<SparseMatrix<double>>& X,
    const std::vector<VectorXd>& y,
    const ArrayXXi& estimability_mask,
    double active_tol,
    const TargetEngineControl& control
) : impl_(new Impl(X, y, estimability_mask, active_tol, control)) {}

ConditionMatrixFreeRefitWorkspace::~ConditionMatrixFreeRefitWorkspace() =
    default;

ConditionMatrixFreeRefitWorkspace::ConditionMatrixFreeRefitWorkspace(
    ConditionMatrixFreeRefitWorkspace&&
) noexcept = default;

ConditionMatrixFreeRefitWorkspace&
ConditionMatrixFreeRefitWorkspace::operator=(
    ConditionMatrixFreeRefitWorkspace&&
) noexcept = default;

Rcpp::List ConditionMatrixFreeRefitWorkspace::refit(
    const MatrixXd& beta,
    double ridge,
    int path_index
) {
    if (!std::isfinite(ridge) || ridge <= 0.0) {
        Rcpp::stop("Matrix-free refit ridge must be finite and positive.");
    }
    if (path_index < 0 || beta.rows() != impl_->problem.predictors ||
        beta.cols() != impl_->problem.tasks || !beta.allFinite()) {
        Rcpp::stop("Matrix-free beta path is not aligned or finite.");
    }
    Rcpp::List current = refit_one(
        impl_->problem,
        beta,
        impl_->estimability_mask,
        ridge,
        impl_->active_tol,
        impl_->control,
        path_index,
        impl_->warm_shared
    );
    impl_->warm_shared = Rcpp::as<VectorXd>(current["shared"]);
    return current;
}

Rcpp::List condition_refit_path_matrix_free(
    const std::vector<SparseMatrix<double>>& X,
    const std::vector<VectorXd>& y,
    const Rcpp::List& beta_path,
    const ArrayXXi& estimability_mask,
    const std::vector<double>& ridge,
    double active_tol,
    const TargetEngineControl& control
) {
    if (beta_path.size() < 1 ||
        beta_path.size() != static_cast<int>(ridge.size())) {
        Rcpp::stop("Invalid matrix-free refit path inputs.");
    }
    ConditionMatrixFreeRefitWorkspace workspace(
        X, y, estimability_mask, active_tol, control
    );
    Rcpp::List out(beta_path.size());
    for (int index = 0; index < beta_path.size(); ++index) {
        MatrixXd beta = Rcpp::as<MatrixXd>(beta_path[index]);
        out[index] = workspace.refit(beta, ridge[index], index);
    }
    return out;
}

double condition_validation_mse_sparse(
    const SparseMatrix<double>& X,
    const VectorXd& y,
    const VectorXd& beta,
    double intercept
) {
    if (X.rows() != y.size() || X.cols() != beta.size() ||
        y.size() < 1 || !y.allFinite() || !beta.allFinite() ||
        !std::isfinite(intercept)) {
        Rcpp::stop("Sparse validation inputs are not aligned or finite.");
    }
    VectorXd residual = y - X * beta;
    residual.array() -= intercept;
    const double mse = residual.squaredNorm() /
        static_cast<double>(y.size());
    if (!std::isfinite(mse)) {
        Rcpp::stop("Sparse residual validation MSE is non-finite.");
    }
    return mse;
}

Rcpp::List condition_matrix_free_refit_self_test() {
    std::vector<Triplet<double>> triplets{
        Triplet<double>(0, 0, 1.0),
        Triplet<double>(1, 1, 2.0),
        Triplet<double>(2, 0, 3.0),
        Triplet<double>(2, 2, 1.0),
        Triplet<double>(3, 1, 1.0)
    };
    SparseMatrix<double> X(4, 3);
    X.setFromTriplets(triplets.begin(), triplets.end());
    X.makeCompressed();
    VectorXd mean = VectorXd::Zero(3);
    for (int column = 0; column < X.outerSize(); ++column) {
        for (SparseMatrix<double>::InnerIterator it(X, column); it; ++it) {
            mean[column] += it.value();
        }
    }
    mean /= static_cast<double>(X.rows());
    VectorXd value(3);
    value << 0.5, -1.0, 2.0;
    MatrixXd dense(X);
    MatrixXd centered = dense.rowwise() - mean.transpose();
    VectorXd reference = centered.transpose() * centered * value;
    VectorXd operated = apply_centered_gram(X, mean, value);
    const double gram_error = (reference - operated).norm() /
        (reference.norm() + 1.0);

    MatrixXd dense_second(4, 3);
    dense_second <<
        0.0, 1.0, 0.5,
        2.0, 0.0, 0.0,
        0.0, 1.5, 1.0,
        1.0, 0.0, 2.0;
    SparseMatrix<double> X_second = dense_second.sparseView();
    X_second.makeCompressed();
    VectorXd y_first(4);
    y_first << 0.8, -0.1, 1.2, 0.3;
    VectorXd y_second(4);
    y_second << -0.2, 0.9, 0.4, 1.1;
    std::vector<SparseMatrix<double>> refit_X{X, X_second};
    std::vector<VectorXd> refit_y{y_first, y_second};
    SparseRefitProblem refit_problem = make_problem(refit_X, refit_y);
    MatrixXd dense_common = MatrixXd::Zero(3, 3);
    std::vector<MatrixXd> dense_gram(2);
    std::vector<VectorXd> dense_rhs(2);
    for (int task = 0; task < 2; ++task) {
        MatrixXd task_X(refit_X[task]);
        VectorXd x_mean = task_X.colwise().mean();
        VectorXd y_centered = refit_y[task].array() -
            refit_y[task].mean();
        MatrixXd centered_X = task_X.rowwise() - x_mean.transpose();
        dense_gram[task] = centered_X.transpose() * centered_X;
        dense_rhs[task] = centered_X.transpose() * y_centered;
        dense_common.noalias() += 0.25 * dense_gram[task];
    }
    dense_common.diagonal().array() += kJitter;
    VectorXd common_operated = apply_common_metric(refit_problem, value);
    VectorXd common_reference = dense_common * value;
    const double common_error =
        (common_operated - common_reference).norm() /
        (common_reference.norm() + 1.0);

    MatrixXd beta_selection = MatrixXd::Zero(3, 2);
    beta_selection(0, 0) = 0.8;
    beta_selection(1, 0) = -0.3;
    beta_selection(1, 1) = 0.4;
    beta_selection(2, 1) = 0.7;
    ArrayXXi estimability = ArrayXXi::Ones(3, 2);
    TargetEngineControl refit_control = default_target_engine_control();
    refit_control.refit_pcg_tolerance = 1e-11;
    refit_control.refit_pcg_max_iterations = 500;
    const double refit_ridge = 0.04;
    ConditionMatrixFreeRefitWorkspace refit_workspace(
        refit_X, refit_y, estimability, 1e-12, refit_control
    );
    Rcpp::List matrix_free_refit = refit_workspace.refit(
        beta_selection, refit_ridge, 0
    );
    const std::vector<std::vector<int>> active{{0, 1}, {1, 2}};
    const int shared_offset = 4;
    MatrixXd block_system = MatrixXd::Zero(7, 7);
    VectorXd block_rhs = VectorXd::Zero(7);
    int offset = 0;
    for (int task = 0; task < 2; ++task) {
        MatrixXd task_lhs = 0.25 * centered_gram_subset(
            refit_X[task], refit_problem.task[task].x_mean, active[task]
        ) + 0.5 * refit_ridge * common_metric_subset(
            refit_problem, active[task]
        );
        task_lhs.diagonal().array() += kJitter;
        MatrixXd cross(active[task].size(), 3);
        for (int row = 0; row < static_cast<int>(active[task].size()); ++row) {
            cross.row(row) = dense_common.row(active[task][row]);
        }
        block_system.block(
            offset, offset, active[task].size(), active[task].size()
        ) = task_lhs;
        block_system.block(
            offset, shared_offset, active[task].size(), 3
        ) = -0.5 * refit_ridge * cross;
        block_system.block(
            shared_offset, offset, 3, active[task].size()
        ) = -0.5 * cross.transpose();
        block_rhs.segment(offset, active[task].size()) = 0.25 *
            subset_vector(dense_rhs[task], active[task]);
        offset += active[task].size();
    }
    block_system.bottomRightCorner(3, 3) = dense_common;
    block_system.bottomRightCorner(3, 3).diagonal().array() += kJitter;
    VectorXd block_solution = block_system.fullPivLu().solve(block_rhs);
    MatrixXd dense_beta = MatrixXd::Zero(3, 2);
    offset = 0;
    for (int task = 0; task < 2; ++task) {
        for (int position = 0;
             position < static_cast<int>(active[task].size());
             ++position) {
            dense_beta(active[task][position], task) =
                block_solution[offset + position];
        }
        offset += active[task].size();
    }
    MatrixXd matrix_free_beta = Rcpp::as<MatrixXd>(
        matrix_free_refit["beta"]
    );
    VectorXd matrix_free_shared = Rcpp::as<VectorXd>(
        matrix_free_refit["shared"]
    );
    VectorXd dense_shared = block_solution.tail(3);
    const double schur_error = (
        (matrix_free_beta - dense_beta).norm() +
        (matrix_free_shared - dense_shared).norm()
    ) / (dense_beta.norm() + dense_shared.norm() + 1.0);

    MatrixXd system(2, 2);
    system << 4.0, 1.0, 1.0, 3.0;
    VectorXd rhs(2);
    rhs << 1.0, 2.0;
    PCGResult pcg = pcg_solve(
        [&](const VectorXd& x) { return system * x; },
        system.diagonal(),
        rhs,
        VectorXd::Zero(2),
        1e-12,
        20,
        "Native self-test"
    );
    VectorXd direct = system.llt().solve(rhs);
    const double pcg_error = (pcg.solution - direct).norm();
    VectorXd validation_beta(3);
    validation_beta << 0.2, -0.3, 0.4;
    VectorXd validation_y(4);
    validation_y << 0.5, -0.2, 1.1, 0.7;
    const double validation_intercept = -0.15;
    VectorXd validation_residual = validation_y -
        dense * validation_beta;
    validation_residual.array() -= validation_intercept;
    const double dense_validation_mse =
        validation_residual.squaredNorm() / validation_y.size();
    const double sparse_validation_mse = condition_validation_mse_sparse(
        X, validation_y, validation_beta, validation_intercept
    );
    const double validation_error = std::abs(
        dense_validation_mse - sparse_validation_mse
    );
    const bool hybrid_preconditioner =
        Rcpp::as<std::string>(matrix_free_refit["preconditioner"]) ==
            "active_union_block_plus_diagonal";
    const bool budget_guard_passed = Rcpp::as<bool>(
        matrix_free_refit["budget_guard_passed"]
    );
    TargetEngineControl bounded_control = refit_control;
    bounded_control.worker_budget_bytes = 2500;
    ConditionMatrixFreeRefitWorkspace bounded_workspace(
        refit_X, refit_y, estimability, 1e-12, bounded_control
    );
    Rcpp::List bounded_refit = bounded_workspace.refit(
        beta_selection, refit_ridge, 0
    );
    const bool bounded_diagonal_fallback =
        Rcpp::as<std::string>(bounded_refit["preconditioner"]) ==
            "diagonal" &&
        Rcpp::as<double>(bounded_refit["dense_workspace_peak_bytes"]) ==
            0.0 &&
        Rcpp::as<double>(bounded_refit["estimated_peak_bytes"]) <=
            static_cast<double>(bounded_control.worker_budget_bytes);
    const bool passed = gram_error < 1e-12 && common_error < 1e-12 &&
        schur_error < 1e-8 && pcg.converged &&
        pcg_error < 1e-10 && validation_error < 1e-12 &&
        hybrid_preconditioner && budget_guard_passed &&
        bounded_diagonal_fallback;
    if (!passed) {
        Rcpp::stop("Native matrix-free numerical self-test failed.");
    }
    return Rcpp::List::create(
        Rcpp::Named("passed") = passed,
        Rcpp::Named("centered_gram_relative_error") = gram_error,
        Rcpp::Named("common_metric_relative_error") = common_error,
        Rcpp::Named("schur_refit_relative_error") = schur_error,
        Rcpp::Named("pcg_absolute_error") = pcg_error,
        Rcpp::Named("pcg_relative_residual") = pcg.relative_residual,
        Rcpp::Named("validation_absolute_error") = validation_error,
        Rcpp::Named("hybrid_preconditioner") = hybrid_preconditioner,
        Rcpp::Named("bounded_diagonal_fallback") =
            bounded_diagonal_fallback,
        Rcpp::Named("budget_guard_passed") = budget_guard_passed,
        Rcpp::Named("estimated_peak_bytes") =
            Rcpp::as<double>(matrix_free_refit["estimated_peak_bytes"])
    );
}
