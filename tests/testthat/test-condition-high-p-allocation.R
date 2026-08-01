test_that("native high-p plan forbids every full predictor square", {
    expect_true(is.loaded(
        "_Pando_condition_native_self_test_cpp", PACKAGE = "Pando"
    ))
    result <- Pando:::.condition_native_self_test_cpp()
    expect_true(result$passed)
    expect_true(result$budget_guard_passed)
    expect_true(result$numerical$passed)
    expect_true(result$numerical$hybrid_preconditioner)
    expect_true(result$numerical$bounded_diagonal_fallback)
    expect_true(result$numerical$budget_guard_passed)
    expect_lt(result$numerical$validation_absolute_error, 1e-12)
    expect_lt(result$numerical$common_metric_relative_error, 1e-12)
    expect_lt(result$numerical$schur_refit_relative_error, 1e-8)
    expect_identical(
        result$execution_plan$path_backend, "sparse_matrix_free"
    )
    expect_identical(
        result$execution_plan$validation_backend, "sparse_residual"
    )
    expect_identical(
        result$execution_plan$refit_backend, "matrix_free_schur_pcg"
    )
    expect_false(
        result$execution_plan$full_predictor_square_allocated
    )
    expect_lte(
        result$execution_plan$matrix_free_bytes,
        result$execution_plan$worker_budget_bytes
    )
})

test_that("engine control rejects unknown or unsafe values", {
    expect_error(
        Pando:::.condition_normalize_engine_control(list(unknown = 1)),
        "Unknown engine_control"
    )
    expect_error(
        Pando:::.condition_normalize_engine_control(list(dense_max_p = 0)),
        "dense_max_p"
    )
    expect_error(
        Pando:::.condition_normalize_engine_control(list(dense_max_p = "2")),
        "dense_max_p"
    )
    expect_error(
        Pando:::.condition_normalize_engine_control(list(
            memory_budget_mb = "256"
        )),
        "memory_budget_mb"
    )
    control <- Pando:::.condition_normalize_engine_control(list(
        memory_budget_mb = 256,
        diagnostics_level = "compact",
        resume = TRUE
    ))
    expect_identical(control$dense_max_p, 2048L)
    expect_identical(control$diagnostics_level, "compact")
    expect_true(control$resume)
})
