test_that('pooled-within-condition screening removes pure condition shifts', {
    condition <- factor(rep(c('A', 'B'), each = 4L))
    within_x <- rep(c(-1, -1, 1, 1), 2L)
    within_y <- rep(c(-1, 1, -1, 1), 2L)
    x <- within_x + rep(c(0, 10), each = 4L)
    y <- within_y + rep(c(0, 10), each = 4L)
    x <- Matrix::Matrix(matrix(x, ncol = 1L), sparse = TRUE)
    colnames(x) <- 'feature'
    y <- Matrix::Matrix(matrix(y, ncol = 1L), sparse = TRUE)

    pooled <- abs(Pando:::.condition_named_cor(x, y)[['feature']])
    within <- abs(Pando:::.condition_within_named_cor(
        x, y, condition
    )[['feature']])
    mask <- Pando:::.condition_component_masks(
        x, y, condition, threshold = 0.1,
        candidate_screen = 'pooled_within_condition'
    )

    expect_gt(pooled, 0.9)
    expect_lt(within, 1e-10)
    expect_false(any(mask))
    expect_identical(colnames(mask), levels(condition))
})

test_that('support-constrained refit preserves condition-specific graphs', {
    set.seed(42)
    n <- 120L
    X_A <- Matrix::Matrix(
        cbind(rnorm(n), rnorm(n)), sparse = TRUE
    )
    X_B <- Matrix::Matrix(
        cbind(rnorm(n), rnorm(n)), sparse = TRUE
    )
    colnames(X_A) <- colnames(X_B) <- c('edge_1', 'edge_2')
    y_A <- 1.5 * as.numeric(X_A[, 1L]) + rnorm(n, sd = 0.05)
    y_B <- -1.25 * as.numeric(X_B[, 2L]) + rnorm(n, sd = 0.05)

    beta_selection <- matrix(
        c(1, 0, 0, -1),
        nrow = 2L,
        byrow = TRUE,
        dimnames = list(c('edge_1', 'edge_2'), c('A', 'B'))
    )
    eligibility <- matrix(
        TRUE,
        nrow = 2L,
        ncol = 2L,
        dimnames = dimnames(beta_selection)
    )

    fit <- Pando:::.condition_refit_shared_baseline(
        X_list = list(A = X_A, B = X_B),
        y_list = list(A = y_A, B = y_B),
        beta_selection = beta_selection,
        eligibility_mask = eligibility,
        ridge = 1e-3,
        condition_weight = 'equal'
    )

    expect_true(fit$converged)
    expect_gt(fit$beta['edge_1', 'A'], 0)
    expect_equal(fit$beta['edge_1', 'B'], 0)
    expect_equal(fit$beta['edge_2', 'A'], 0)
    expect_lt(fit$beta['edge_2', 'B'], 0)
    expect_false(fit$active_mask['edge_1', 'B'])
    expect_false(fit$active_mask['edge_2', 'A'])
    expect_equal(
        fit$beta_condition,
        sweep(fit$delta_condition, 1L, fit$beta_shared, '+'),
        tolerance = 1e-8
    )
})

test_that('refit separates unavailable edges from estimated zero', {
    set.seed(7)
    X <- Matrix::Matrix(matrix(rnorm(160), ncol = 2L), sparse = TRUE)
    colnames(X) <- c('edge_1', 'edge_2')
    beta_selection <- matrix(
        c(1, 0, 0, 0),
        nrow = 2L,
        byrow = TRUE,
        dimnames = list(c('edge_1', 'edge_2'), c('A', 'B'))
    )
    eligibility <- matrix(
        c(TRUE, TRUE, TRUE, FALSE),
        nrow = 2L,
        byrow = TRUE,
        dimnames = dimnames(beta_selection)
    )

    fit <- Pando:::.condition_refit_shared_baseline(
        X_list = list(A = X, B = X),
        y_list = list(A = as.numeric(X[, 1L]), B = as.numeric(X[, 1L])),
        beta_selection = beta_selection,
        eligibility_mask = eligibility,
        ridge = 1e-3,
        condition_weight = 'equal'
    )

    expect_equal(fit$beta_condition['edge_1', 'B'], 0)
    expect_true(is.na(fit$beta_condition['edge_2', 'B']))
    expect_false(fit$active_mask['edge_1', 'B'])
    expect_false(fit$estimability_mask['edge_2', 'B'])
})

test_that('condition GOF target fitter is patched with named condition lists', {
    body_text <- paste(
        deparse(body(Pando:::.condition_fit_target_v2_patched)),
        collapse = '\n'
    )
    expect_match(body_text, 'names(condition_gof)', fixed = TRUE)
})

test_that('new defaults select the condition-sparse engine', {
    defaults <- formals(Pando:::infer_condition_grn.GRNData)
    expect_identical(
        eval(defaults$method)[[1L]],
        'shared_baseline_condition_sparse'
    )
    expect_identical(
        eval(defaults$candidate_screen)[[1L]],
        'pooled_within_condition'
    )
    expect_equal(eval(defaults$condition_mix), 0.5)
})