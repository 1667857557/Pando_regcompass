test_that("hybrid inner CV uses exact sufficient-statistic validation", {
    set.seed(20260801)
    n <- 40L
    p <- 5L
    X <- list(
        Control = Matrix::Matrix(matrix(rnorm(n * p), n, p), sparse = TRUE),
        Drug = Matrix::Matrix(matrix(rnorm(n * p), n, p), sparse = TRUE)
    )
    colnames(X$Control) <- colnames(X$Drug) <- paste0("e", seq_len(p))
    y <- list(
        Control = as.numeric(0.8 * X$Control[, 1] - 0.4 * X$Control[, 3] +
            rnorm(n, sd = 0.15)),
        Drug = as.numeric(-0.7 * X$Drug[, 1] + 0.5 * X$Drug[, 2] +
            rnorm(n, sd = 0.15))
    )
    mask <- matrix(TRUE, p, 2L,
                   dimnames = list(colnames(X$Control), names(X)))
    lambda <- c(0.2, 0.08, 0.03)
    seed <- 443L

    native <- Pando:::.condition_fit_target_native_engine(
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
        outer_nfolds = 3L,
        inner_nfolds = 2L,
        lambda_selection = "lambda.1se",
        seed = seed,
        max_iter = 2500L,
        tol_objective = 1e-8,
        tol_coef = 1e-7
    )
    reference <- .reference_complete_target_fit(
        X, y, mask, lambda,
        seed = seed,
        max_iter = 2500L,
        tol_objective = 1e-8,
        tol_coef = 1e-7
    )

    expect_identical(
        native$inner_cv_backend,
        "exact_refit_validation_dense_support_or_sparse_residual_v1"
    )
    expect_true(all(grepl(
        "centered_gram|sparse_matrix_free|intercept_only",
        native$full_cv$solver_backend
    )))
    expect_equal(native$full_cv$cv_mean, reference$full_cv$cv_mean,
                 tolerance = 8e-7)
    expect_equal(native$full_cv$cv_se, reference$full_cv$cv_se,
                 tolerance = 8e-7)
    expect_identical(native$full_cv$selected_index,
                     reference$full_cv$selected_index)
    expect_equal(native$cv$fold_selected_lambda,
                 reference$cv$fold_selected_lambda,
                 tolerance = 1e-10)
    expect_equal(
        unlist(native$cv$oof_prediction, use.names = FALSE),
        unlist(reference$cv$oof_prediction, use.names = FALSE),
        tolerance = 5e-6
    )
    expect_true(all(unlist(native$cv$oof_assignment_count) == 1L))
})

test_that("hybrid cost model retains sparse matrix-free path for sparse wide folds", {
    set.seed(97)
    n <- 36L
    p <- 96L
    X <- list(
        A = Matrix::rsparsematrix(n, p, density = 0.035),
        B = Matrix::rsparsematrix(n, p, density = 0.035)
    )
    colnames(X$A) <- colnames(X$B) <- paste0("e", seq_len(p))
    y <- list(A = rnorm(n), B = rnorm(n))
    mask <- matrix(TRUE, p, 2L,
                   dimnames = list(colnames(X$A), names(X)))
    fit <- Pando:::.condition_fit_target_native_engine(
        X_raw_list = X,
        y_raw_list = y,
        coefficient_mask = mask,
        comparison_conditions = names(X),
        lambda = c(0.2, 0.08),
        nlambda = 2L,
        lambda_min_ratio = NULL,
        alpha = 0.5,
        condition_mix = 0.5,
        active_tol = 1e-8,
        outer_nfolds = 2L,
        inner_nfolds = 2L,
        lambda_selection = "lambda.1se",
        seed = 98L,
        max_iter = 600L,
        tol_objective = 1e-7,
        tol_coef = 1e-6
    )
    expect_true(any(grepl(
        "sparse_matrix_free",
        fit$full_cv$solver_backend,
        fixed = TRUE
    )))
    expect_true(all(unlist(fit$cv$oof_assignment_count) == 1L))
})

test_that("native ABI metadata advertises the canonical hybrid engine", {
    description <- utils::packageDescription("Pando")
    expect_identical(description[["Config/Pando/NativeSparseABI"]], "6")
    expect_identical(
        description[["Config/Pando/ConditionInnerCVBackend"]],
        "exact-refit-validation-sparse-residual-v1"
    )
    expect_identical(
        description[["Config/Pando/ConditionTargetEngineBackend"]],
        "cpp-eigen-memory-bounded-hybrid-target-v1"
    )
    expect_identical(
        description[["Config/Pando/ConditionMemoryContract"]],
        "no-full-p2-on-high-p-path-v1"
    )
})
