test_that('shared screening retains opposite condition effects', {
    x <- Matrix::Matrix(cbind(feature = c(1:8, 1:8)), sparse = TRUE)
    y <- Matrix::Matrix(c(c(1:8), rev(1:8)), ncol = 1, sparse = TRUE)
    condition <- factor(rep(c('Control', 'Drug'), each = 8))

    score <- Pando:::.condition_within_association_score(x, y, condition)
    shared_keep <- Pando:::.condition_screen_columns(
        x, y, condition, threshold = 0.8,
        candidate_screen = 'pooled_within_condition'
    )

    expect_equal(score[['feature']], 1, tolerance = 1e-12)
    expect_true(shared_keep[[1]])
})

test_that('edge union never pairs TF and peak retained in different conditions', {
    edges <- data.frame(
        tf = 'TF1', region = 'peak1', target = 'GENE1',
        edge_id = 'TF1\001peak1\001GENE1'
    )
    peak_mask <- matrix(
        c(TRUE, FALSE), nrow = 1,
        dimnames = list('peak1', c('Control', 'Drug'))
    )
    tf_mask <- matrix(
        c(FALSE, TRUE), nrow = 1,
        dimnames = list('TF1', c('Control', 'Drug'))
    )

    mask <- Pando:::.condition_edge_mask(edges, peak_mask, tf_mask)
    expect_false(any(mask))
})

test_that('final TF by peak predictors use one pooled edge transform', {
    gene_data <- Matrix::Matrix(
        cbind(TF1 = c(1, 2, 4, 8)), sparse = TRUE
    )
    peak_data <- Matrix::Matrix(
        cbind(peak1 = c(1, 3, 2, 5)), sparse = TRUE
    )
    response <- Matrix::Matrix(cbind(GENE1 = c(2, 3, 5, 9)), sparse = TRUE)
    edges <- data.frame(
        tf = 'TF1', region = 'peak1', target = 'GENE1',
        edge_id = 'TF1\001peak1\001GENE1'
    )
    design <- Pando:::.condition_build_design(
        response, gene_data, peak_data, edges, scale = TRUE
    )
    raw_edge <- as.numeric(gene_data[, 'TF1'] * peak_data[, 'peak1'])

    expect_equal(design$predictor_center, mean(raw_edge))
    expect_equal(design$predictor_scale, stats::sd(raw_edge))
    expect_equal(as.numeric(Matrix::colMeans(design$X)), 0, tolerance = 1e-12)
    expect_equal(
        Pando:::.condition_column_variance(design$X), 1, tolerance = 1e-12
    )
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

test_that('shared-design independent coefficients are separable by condition', {
    set.seed(14)
    X <- list(
        Control = Matrix::Matrix(matrix(stats::rnorm(160), ncol = 4)),
        Drug = Matrix::Matrix(matrix(stats::rnorm(160), ncol = 4))
    )
    y <- list(Control = stats::rnorm(40), Drug = stats::rnorm(40))
    first <- Pando:::.condition_fit_multitask_lambda(
        X, y, lambda = 0.05, alpha = 0.5, condition_mix = 1,
        condition_weight = 'equal', max_iter = 5000,
        tol_objective = 1e-10, tol_coef = 1e-9
    )
    changed_y <- y
    changed_y$Drug <- changed_y$Drug + 10 * as.numeric(X$Drug[, 1])
    second <- Pando:::.condition_fit_multitask_lambda(
        X, changed_y, lambda = 0.05, alpha = 0.5, condition_mix = 1,
        condition_weight = 'equal', max_iter = 5000,
        tol_objective = 1e-10, tol_coef = 1e-9
    )

    expect_equal(
        first$beta[, 1L], second$beta[, 1L],
        tolerance = 1e-7
    )
    expect_false(isTRUE(all.equal(
        first$beta[, 2L], second$beta[, 2L], tolerance = 1e-3
    )))
})

test_that('eligibility mask fixes absent condition edges at zero', {
    set.seed(15)
    X <- list(
        Control = Matrix::Matrix(matrix(stats::rnorm(120), ncol = 3)),
        Drug = Matrix::Matrix(matrix(stats::rnorm(120), ncol = 3))
    )
    y <- list(Control = stats::rnorm(40), Drug = stats::rnorm(40))
    mask <- matrix(
        c(TRUE, TRUE, TRUE, FALSE, TRUE, FALSE),
        nrow = 3, ncol = 2,
        dimnames = list(NULL, c('Control', 'Drug'))
    )
    fit <- Pando:::.condition_fit_multitask_lambda(
        X, y, lambda = 0.01, alpha = 0.5, condition_mix = 1,
        coefficient_mask = mask
    )

    expect_equal(fit$beta[!mask], rep(0, sum(!mask)))
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
