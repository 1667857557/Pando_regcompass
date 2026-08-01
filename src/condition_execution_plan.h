#ifndef PANDO_CONDITION_EXECUTION_PLAN_H
#define PANDO_CONDITION_EXECUTION_PLAN_H

#include <cstddef>
#include <string>

enum class PathBackend {
    DenseCenteredGram,
    SparseMatrixFree
};

enum class ValidationBackend {
    DenseSupportStats,
    SparseResidual
};

enum class RefitBackend {
    DenseDirectSchur,
    MatrixFreeSchurPCG
};

struct TargetEngineControl {
    std::size_t worker_budget_bytes;
    int dense_max_p;
    int lambda_batch_size;
    double refit_pcg_tolerance;
    int refit_pcg_max_iterations;
    std::string preconditioner;
    bool compact_diagnostics;
};

struct TargetExecutionPlan {
    PathBackend path;
    ValidationBackend validation;
    RefitBackend refit;
    int predictors;
    int tasks;
    int lambda_batch_size;
    std::size_t nonzeros;
    std::size_t dense_path_bytes;
    std::size_t dense_validation_bytes;
    std::size_t dense_refit_bytes;
    std::size_t matrix_free_bytes;
    std::size_t worker_budget_bytes;
    bool full_predictor_square_allocated;
};

TargetEngineControl default_target_engine_control();

TargetExecutionPlan plan_target_execution(
    int predictors,
    int tasks,
    std::size_t nonzeros,
    const TargetEngineControl& control
);

const char* path_backend_name(PathBackend backend);
const char* validation_backend_name(ValidationBackend backend);
const char* refit_backend_name(RefitBackend backend);

#endif
