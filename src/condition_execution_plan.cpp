#include "condition_execution_plan.h"
#include <algorithm>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace {

std::size_t saturated_multiply(std::size_t left, std::size_t right) {
    if (left == 0 || right == 0) return 0;
    if (left > std::numeric_limits<std::size_t>::max() / right) {
        return std::numeric_limits<std::size_t>::max();
    }
    return left * right;
}

std::size_t saturated_add(std::size_t left, std::size_t right) {
    if (left > std::numeric_limits<std::size_t>::max() - right) {
        return std::numeric_limits<std::size_t>::max();
    }
    return left + right;
}

std::size_t square_bytes(int predictors, std::size_t copies) {
    const std::size_t p = static_cast<std::size_t>(predictors);
    return saturated_multiply(
        saturated_multiply(saturated_multiply(p, p), copies),
        sizeof(double)
    );
}

} // namespace

TargetEngineControl default_target_engine_control() {
    return TargetEngineControl{
        static_cast<std::size_t>(1024) * 1024 * 1024,
        2048,
        1,
        1e-8,
        2000,
        "hybrid",
        false
    };
}

TargetExecutionPlan plan_target_execution(
    int predictors,
    int tasks,
    std::size_t nonzeros,
    const TargetEngineControl& control
) {
    if (predictors < 1 || tasks < 2 || control.dense_max_p < 1 ||
        control.worker_budget_bytes < 1 || control.lambda_batch_size < 1 ||
        control.refit_pcg_tolerance <= 0.0 ||
        control.refit_pcg_max_iterations < 1) {
        throw std::invalid_argument("Invalid target execution-plan controls.");
    }
    const std::size_t k = static_cast<std::size_t>(tasks);
    const std::size_t p = static_cast<std::size_t>(predictors);
    const std::size_t sparse_values = saturated_multiply(
        nonzeros,
        sizeof(double) + 2 * sizeof(int)
    );
    const std::size_t sparse_resident = saturated_multiply(
        sparse_values, 6
    );
    const std::size_t dense_path = saturated_add(
        sparse_resident, square_bytes(predictors, 2 * k + 8)
    );
    const std::size_t dense_validation = saturated_add(
        sparse_resident, square_bytes(predictors, 2 * k + 10)
    );
    const std::size_t dense_refit = saturated_add(
        sparse_resident, square_bytes(predictors, 3 * k + 12)
    );
    const std::size_t coefficient_work = saturated_multiply(
        saturated_multiply(
            p,
            k * (14 + static_cast<std::size_t>(control.lambda_batch_size)) +
                20
        ),
        sizeof(double)
    );
    const std::size_t matrix_free = saturated_add(
        sparse_resident, coefficient_work
    );
    const std::size_t largest_dense = std::max(
        dense_path, std::max(dense_validation, dense_refit)
    );
    const bool dense_allowed = predictors <= control.dense_max_p &&
        largest_dense <= control.worker_budget_bytes;
    if (!dense_allowed && matrix_free > control.worker_budget_bytes) {
        std::ostringstream message;
        message << "No target-engine backend fits the per-worker memory budget "
                << "(p=" << predictors << ", tasks=" << tasks
                << ", nnz=" << nonzeros
                << ", matrix_free_estimated_bytes=" << matrix_free
                << ", worker_budget_bytes=" << control.worker_budget_bytes
                << ").";
        throw std::invalid_argument(message.str());
    }
    const long double gram_work = static_cast<long double>(tasks) *
        predictors * predictors;
    const long double sparse_work = 2.0L *
        static_cast<long double>(std::max<std::size_t>(nonzeros, 1));
    const bool dense_path_is_efficient = gram_work <= 4.0L * sparse_work;

    TargetExecutionPlan out;
    out.path = dense_allowed && dense_path_is_efficient ?
        PathBackend::DenseCenteredGram : PathBackend::SparseMatrixFree;
    out.validation = dense_allowed ?
        ValidationBackend::DenseSupportStats :
        ValidationBackend::SparseResidual;
    out.refit = dense_allowed ?
        RefitBackend::DenseDirectSchur :
        RefitBackend::MatrixFreeSchurPCG;
    out.predictors = predictors;
    out.tasks = tasks;
    out.lambda_batch_size = control.lambda_batch_size;
    out.nonzeros = nonzeros;
    out.dense_path_bytes = dense_path;
    out.dense_validation_bytes = dense_validation;
    out.dense_refit_bytes = dense_refit;
    out.matrix_free_bytes = matrix_free;
    out.worker_budget_bytes = control.worker_budget_bytes;
    out.full_predictor_square_allocated = dense_allowed;
    return out;
}

const char* path_backend_name(PathBackend backend) {
    return backend == PathBackend::DenseCenteredGram ?
        "dense_centered_gram" : "sparse_matrix_free";
}

const char* validation_backend_name(ValidationBackend backend) {
    return backend == ValidationBackend::DenseSupportStats ?
        "dense_support_statistics" : "sparse_residual";
}

const char* refit_backend_name(RefitBackend backend) {
    return backend == RefitBackend::DenseDirectSchur ?
        "dense_direct_schur" : "matrix_free_schur_pcg";
}
