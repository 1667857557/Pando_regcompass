test_that("the per-target engine is registered and is the canonical target path", {
    expect_true(is.loaded(
        "_Pando_condition_fit_target_engine_cpp", PACKAGE = "Pando"
    ))
    body_text <- paste(deparse(body(Pando:::.condition_fit_target)), collapse = "\n")
    expect_match(body_text, ".condition_fit_target_native_engine", fixed = TRUE)
    expect_false(grepl(
        ".condition_nested_crossfit_within_cell_type(",
        body_text,
        fixed = TRUE
    ))
    expect_false(grepl(
        ".condition_select_lambda_nested(",
        body_text,
        fixed = TRUE
    ))
    expect_false(grepl(
        ".condition_fit_multitask_path(",
        body_text,
        fixed = TRUE
    ))
    expect_false(grepl(
        ".condition_refit_shared_baseline(",
        body_text,
        fixed = TRUE
    ))
})

test_that("fused target CV, lambda selection and refit match the R architecture", {
    set.seed(20260801)
    n <- 36L
    p <- 4L
    X <- list(
        Control = Matrix::Matrix(matrix(rnorm(n * p), n, p), sparse = TRUE),
        Drug = Matrix::Matrix(matrix(rnorm(n * p), n, p), sparse = TRUE)
    )
    predictor <- paste0("e", seq_len(p))
    colnames(X$Control) <- colnames(X$Drug) <- predictor
    y <- list(
        Control = 0.9 * X$Control[, 1] - 0.35 * X$Control[, 3] +
            rnorm(n, sd = 0.12),
        Drug = -0.75 * X$Drug[, 1] + 0.45 * X$Drug[, 2] +
            rnorm(n, sd = 0.12)
    )
    y <- lapply(y, as.numeric)
    mask <- matrix(
        TRUE, p, 2L,
        dimnames = list(predictor, names(X))
    )
    lambda <- c(0.18, 0.07, 0.025)
    seed <- 137L

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
        max_iter = 2000L,
        tol_objective = 1e-8,
        tol_coef = 1e-7
    )
    reference <- .reference_complete_target_fit(
        X, y, mask, lambda,
        seed = seed,
        max_iter = 2000L,
        tol_objective = 1e-8,
        tol_coef = 1e-7
    )

    expect_identical(native$lambda_path, lambda)
    expect_equal(
        native$cv$fold_selected_lambda,
        reference$cv$fold_selected_lambda,
        tolerance = 1e-10
    )
    expect_equal(
        unlist(native$cv$oof_prediction, use.names = FALSE),
        unlist(reference$cv$oof_prediction, use.names = FALSE),
        tolerance = 3e-6
    )
    expect_equal(
        unlist(native$cv$projection_condition_full_oof, use.names = FALSE),
        unlist(reference$cv$projection_condition_full_oof, use.names = FALSE),
        tolerance = 3e-6
    )
    expect_equal(
        unlist(native$cv$projection_common_oof, use.names = FALSE),
        unlist(reference$cv$projection_common_oof, use.names = FALSE),
        tolerance = 3e-6
    )
    expect_identical(
        native$cv$oof_assignment_count,
        reference$cv$oof_assignment_count
    )
    expect_equal(native$full_cv$cv_mean, reference$full_cv$cv_mean,
                 tolerance = 3e-7)
    expect_equal(native$full_cv$cv_se, reference$full_cv$cv_se,
                 tolerance = 3e-7)
    expect_identical(
        native$full_cv$selected_index,
        reference$full_cv$selected_index
    )
    expect_equal(native$beta_selection, reference$selected$beta,
                 tolerance = 3e-6)
    expect_equal(native$refit$beta, reference$refit$beta,
                 tolerance = 3e-6)
    expect_equal(native$refit$beta_shared, reference$refit$beta_shared,
                 tolerance = 3e-6)
    expect_equal(native$refit$intercept, reference$refit$intercept,
                 tolerance = 3e-6)
    expect_identical(native$refit$support_mask, reference$refit$support_mask)
    expect_identical(
        native$refit$estimability_mask,
        reference$refit$estimability_mask
    )
    expect_equal(
        native$full_transform$predictor_center,
        reference$transform$predictor_center,
        tolerance = 1e-12
    )
    expect_equal(
        native$full_transform$predictor_scale,
        reference$transform$predictor_scale,
        tolerance = 1e-12
    )
    expect_identical(
        native$backend,
        "cpp_eigen_fused_target_nested_cv_hybrid_gram_refit_validation_stats"
    )
})

test_that("the fused engine preserves condition-specific signs and structural zeros", {
    set.seed(84)
    n <- 42L
    e1_control <- rnorm(n)
    e1_drug <- rnorm(n)
    X <- list(
        Control = Matrix::Matrix(cbind(
            e1 = e1_control,
            e2 = rnorm(n),
            condition_only = rnorm(n)
        ), sparse = TRUE),
        Drug = Matrix::Matrix(cbind(
            e1 = e1_drug,
            e2 = rnorm(n),
            condition_only = rep(0, n)
        ), sparse = TRUE)
    )
    y <- list(
        Control = 1.3 * e1_control + 0.8 * X$Control[, "condition_only"] +
            rnorm(n, sd = 0.06),
        Drug = -1.2 * e1_drug + rnorm(n, sd = 0.06)
    )
    y <- lapply(y, as.numeric)
    mask <- matrix(
        TRUE, 3L, 2L,
        dimnames = list(colnames(X$Control), names(X))
    )
    fit <- Pando:::.condition_fit_target_native_engine(
        X_raw_list = X,
        y_raw_list = y,
        coefficient_mask = mask,
        comparison_conditions = names(X),
        lambda = c(0.02, 0.005),
        nlambda = 2L,
        lambda_min_ratio = NULL,
        alpha = 0.5,
        condition_mix = 0.5,
        active_tol = 1e-8,
        outer_nfolds = 3L,
        inner_nfolds = 2L,
        lambda_selection = "lambda.min",
        seed = 91L,
        max_iter = 3000L,
        tol_objective = 1e-9,
        tol_coef = 1e-8
    )
    expect_lt(
        fit$refit$beta["e1", "Control"] *
            fit$refit$beta["e1", "Drug"],
        0
    )
    expect_true(fit$refit$estimability_mask["condition_only", "Control"])
    expect_false(fit$refit$estimability_mask["condition_only", "Drug"])
    expect_true(is.na(
        fit$refit$beta_condition["condition_only", "Drug"]
    ))
    expect_true(all(unlist(fit$cv$oof_assignment_count) == 1L))
    expect_true(all(is.finite(unlist(fit$cv$oof_prediction))))
})

test_that("automatic lambda construction remains numerically equivalent", {
    set.seed(501)
    X <- list(
        A = Matrix::Matrix(matrix(rnorm(240), 60, 4), sparse = TRUE),
        B = Matrix::Matrix(matrix(rnorm(240), 60, 4), sparse = TRUE)
    )
    colnames(X$A) <- colnames(X$B) <- paste0("e", 1:4)
    y <- list(A = rnorm(60), B = rnorm(60))
    mask <- matrix(TRUE, 4, 2, dimnames = list(colnames(X$A), names(X)))
    fold_stats <- Pando:::.condition_build_fold_statistics(X, y, mask)
    transform <- Pando:::.condition_build_balanced_transform(
        X, y, fold_statistics = fold_stats
    )
    scaled <- Pando:::.condition_apply_balanced_transform(X, y, transform)
    expected <- Pando:::.condition_make_lambda_path(
        scaled$X,
        scaled$y,
        alpha = 0.5,
        condition_mix = 0.5,
        condition_weight = "equal",
        coefficient_mask = mask,
        nlambda = 8L,
        lambda_min_ratio = 0.01
    )
    native <- Pando:::.condition_fit_target_native_engine(
        X_raw_list = X,
        y_raw_list = y,
        coefficient_mask = mask,
        comparison_conditions = names(X),
        lambda = NULL,
        nlambda = 8L,
        lambda_min_ratio = 0.01,
        alpha = 0.5,
        condition_mix = 0.5,
        active_tol = 1e-8,
        outer_nfolds = 3L,
        inner_nfolds = 2L,
        lambda_selection = "lambda.1se",
        seed = 73L,
        max_iter = 1000L,
        tol_objective = 1e-7,
        tol_coef = 1e-6
    )
    expect_equal(native$lambda_path, expected, tolerance = 1e-10)
})
