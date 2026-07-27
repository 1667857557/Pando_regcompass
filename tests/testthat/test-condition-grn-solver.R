test_that('condition-wise screening retains cancelling condition effects', {
    x <- Matrix::Matrix(cbind(feature = c(1:8, 1:8)), sparse = TRUE)
    y <- Matrix::Matrix(c(c(1:8), rev(1:8)), ncol = 1, sparse = TRUE)
    condition <- factor(rep(c('Control', 'Drug'), each = 8))

    union_keep <- Pando:::.condition_screen_columns(
        x, y, condition, threshold = 0.8, candidate_screen = 'condition_union'
    )
    pooled_keep <- Pando:::.condition_screen_columns(
        x, y, condition, threshold = 0.8, candidate_screen = 'pooled'
    )

    expect_true(union_keep[[1]])
    expect_false(pooled_keep[[1]])
})

test_that('sparse-group proximal map supports group and element sparsity', {
    value <- matrix(c(3, 3, 0.2, 0.2), nrow = 2, byrow = TRUE)

    group_only <- Pando:::.condition_sparse_group_prox(
        value,
        step = 1,
        lambda = 1,
        alpha = 1,
        condition_mix = 0
    )
    element_only <- Pando:::.condition_sparse_group_prox(
        value,
        step = 1,
        lambda = 1,
        alpha = 1,
        condition_mix = 1
    )

    expect_true(all(group_only[2, ] == 0))
    expect_equal(element_only, matrix(c(2, 2, 0, 0), nrow = 2, byrow = TRUE))
})

test_that('multi-task solver recovers opposite condition signs', {
    set.seed(11)
    x_control <- Matrix::Matrix(matrix(stats::rnorm(80), ncol = 1), sparse = TRUE)
    x_drug <- Matrix::Matrix(matrix(stats::rnorm(80), ncol = 1), sparse = TRUE)
    y_control <- 2 * as.numeric(x_control[, 1]) + stats::rnorm(80, sd = 0.05)
    y_drug <- -2 * as.numeric(x_drug[, 1]) + stats::rnorm(80, sd = 0.05)

    fit <- Pando:::.condition_fit_multitask_lambda(
        X_list = list(Control = x_control, Drug = x_drug),
        y_list = list(Control = y_control, Drug = y_drug),
        lambda = 1e-4,
        alpha = 0.5,
        condition_mix = 0.5,
        max_iter = 2000,
        tol_objective = 1e-9,
        tol_coef = 1e-8
    )

    expect_gt(fit$beta[1, 1], 0)
    expect_lt(fit$beta[1, 2], 0)
    expect_true(is.finite(fit$objective))
})

test_that('large lambda produces an all-zero coefficient matrix', {
    set.seed(12)
    X <- list(
        A = Matrix::Matrix(matrix(stats::rnorm(80), ncol = 2), sparse = TRUE),
        B = Matrix::Matrix(matrix(stats::rnorm(80), ncol = 2), sparse = TRUE)
    )
    y <- list(A = stats::rnorm(40), B = stats::rnorm(40))

    fit <- Pando:::.condition_fit_multitask_lambda(
        X_list = X,
        y_list = y,
        lambda = 1e6,
        alpha = 1,
        condition_mix = 0.5,
        max_iter = 100
    )

    expect_equal(fit$beta, matrix(0, 2, 2))
})

test_that('condition-stratified folds contain every fold in every condition', {
    folds <- Pando:::.condition_make_folds(
        list(A = numeric(20), B = numeric(25)), nfolds = 5, seed = 99
    )
    expect_equal(sort(unique(folds$A)), 1:5)
    expect_equal(sort(unique(folds$B)), 1:5)
    expect_equal(length(folds$A), 20)
    expect_equal(length(folds$B), 25)
})

test_that('lambda path is decreasing and begins at a zero-model scale', {
    set.seed(13)
    X <- list(
        A = Matrix::Matrix(matrix(stats::rnorm(120), ncol = 3), sparse = TRUE),
        B = Matrix::Matrix(matrix(stats::rnorm(120), ncol = 3), sparse = TRUE)
    )
    y <- list(A = stats::rnorm(40), B = stats::rnorm(40))
    lambda <- Pando:::.condition_make_lambda_path(
        X, y, alpha = 0.5, condition_mix = 0.5, nlambda = 8
    )

    expect_length(lambda, 8)
    expect_true(all(diff(lambda) < 0))
    first_fit <- Pando:::.condition_fit_multitask_lambda(
        X, y, lambda = lambda[[1]], alpha = 0.5, condition_mix = 0.5
    )
    expect_true(all(abs(first_fit$beta) < 1e-8))
})
