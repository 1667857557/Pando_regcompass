test_that("compact diagnostics do not alter fitted numerical fields", {
    set.seed(144)
    X <- list(
        A = Matrix::Matrix(matrix(rnorm(160), 40, 4), sparse = TRUE),
        B = Matrix::Matrix(matrix(rnorm(160), 40, 4), sparse = TRUE)
    )
    colnames(X$A) <- colnames(X$B) <- paste0("e", 1:4)
    y <- list(A = rnorm(40), B = rnorm(40))
    mask <- matrix(TRUE, 4, 2, dimnames = list(colnames(X$A), names(X)))
    run <- function(level) Pando:::.condition_fit_target_native_engine(
        X_raw_list = X, y_raw_list = y, coefficient_mask = mask,
        comparison_conditions = names(X), lambda = c(0.2, 0.08),
        nlambda = 2L, lambda_min_ratio = NULL, alpha = 0.5,
        condition_mix = 0.5, active_tol = 1e-8, outer_nfolds = 2L,
        inner_nfolds = 2L, lambda_selection = "lambda.1se", seed = 145L,
        max_iter = 1200L, tol_objective = 1e-7, tol_coef = 1e-6,
        engine_control = Pando:::.condition_normalize_engine_control(list(
            diagnostics_level = level
        ))
    )
    full <- run("full")
    compact <- run("compact")
    expect_equal(compact$refit$beta, full$refit$beta, tolerance = 1e-10)
    expect_equal(compact$cv$oof_prediction, full$cv$oof_prediction,
                 tolerance = 1e-10)
    expect_identical(compact$full_cv$selected_index,
                     full$full_cv$selected_index)
    expect_length(compact$cv$fold_transform, 0L)
    expect_length(compact$cv$fold_inner_cv, 0L)
    expect_length(compact$full_cv$fold_transform, 0L)
})
