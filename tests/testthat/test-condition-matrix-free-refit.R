test_that("forced matrix-free target engine matches dense reference", {
    set.seed(808)
    n <- 36L
    p <- 5L
    X <- list(
        A = Matrix::Matrix(matrix(rnorm(n * p), n, p), sparse = TRUE),
        B = Matrix::Matrix(matrix(rnorm(n * p), n, p), sparse = TRUE)
    )
    colnames(X$A) <- colnames(X$B) <- paste0("e", seq_len(p))
    y <- list(
        A = as.numeric(X$A[, 1] - 0.4 * X$A[, 3] + rnorm(n, sd = 0.1)),
        B = as.numeric(-0.8 * X$B[, 1] + 0.5 * X$B[, 2] + rnorm(n, sd = 0.1))
    )
    mask <- matrix(TRUE, p, 2L, dimnames = list(colnames(X$A), names(X)))
    lambda <- c(0.15, 0.05)
    seed <- 809L
    fit <- Pando:::.condition_fit_target_native_engine(
        X_raw_list = X,
        y_raw_list = y,
        coefficient_mask = mask,
        comparison_conditions = names(X),
        lambda = lambda,
        nlambda = length(lambda),
        lambda_min_ratio = NULL,
        alpha = 0.5,
        condition_mix = 0.5,
        active_tol = 1e-8,
        outer_nfolds = 2L,
        inner_nfolds = 2L,
        lambda_selection = "lambda.1se",
        seed = seed,
        max_iter = 2500L,
        tol_objective = 1e-8,
        tol_coef = 1e-7,
        engine_control = Pando:::.condition_normalize_engine_control(list(
            dense_max_p = 1L,
            memory_budget_mb = 128,
            refit_pcg_tol = 1e-9
        ))
    )
    reference <- .reference_complete_target_fit(
        X, y, mask, lambda,
        outer_nfolds = 2L,
        inner_nfolds = 2L,
        seed = seed,
        max_iter = 2500L,
        tol_objective = 1e-8,
        tol_coef = 1e-7
    )
    expect_identical(fit$execution$path_backend, "sparse_matrix_free")
    expect_identical(fit$execution$validation_backend, "sparse_residual")
    expect_identical(fit$execution$refit_backend, "matrix_free_schur_pcg")
    expect_false(fit$execution$full_predictor_square_allocated)
    expect_identical(
        fit$refit$common_metric,
        "pooled_weighted_predictor_gram_cpp_exact_schur"
    )
    expect_true(is.finite(fit$refit$pcg_iterations))
    expect_true(is.finite(fit$refit$pcg_relative_residual))
    expect_true(fit$refit$budget_guard_passed)
    expect_lte(
        fit$refit$estimated_peak_bytes,
        fit$execution$worker_budget_bytes
    )
    expect_identical(
        fit$refit$preconditioner,
        "active_union_block_plus_diagonal"
    )
    expect_gt(fit$refit$preconditioner_active_size, 0L)
    expect_identical(fit$full_cv$selected_index,
                     reference$full_cv$selected_index)
    expect_equal(fit$full_cv$cv_mean, reference$full_cv$cv_mean,
                 tolerance = 2e-6)
    expect_equal(fit$refit$beta, reference$refit$beta, tolerance = 2e-5)
    expect_equal(fit$refit$beta_shared, reference$refit$beta_shared,
                 tolerance = 2e-5)
    expect_equal(
        unlist(fit$cv$oof_prediction, use.names = FALSE),
        unlist(reference$cv$oof_prediction, use.names = FALSE),
        tolerance = 2e-5
    )
    expect_true(all(unlist(fit$cv$oof_assignment_count) == 1L))
})
