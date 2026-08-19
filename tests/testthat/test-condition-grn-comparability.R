test_that("Scheme E reproduces the two-condition scalar closed form", {
    i_delta <- 4
    d_delta <- 1
    i_gamma <- matrix(2 * i_delta, 1, 1)
    d_gamma <- d_delta / sqrt(2)
    rhs <- as.numeric(i_gamma * d_gamma)
    block <- Pando:::.condition_scheme_e_block_inverse_roots(
        i_gamma, p = 1L, k_minus_one = 1L,
        edge_keep = TRUE, rank_tol = 1e-12
    )
    solved <- Pando:::.condition_scheme_e_fista(
        i_gamma, rhs, block,
        Pando:::.condition_scheme_e_control(list(
            solver_tol = 1e-12, solver_max_iter = 10000L
        ))
    )
    expected_delta <- sign(d_delta) *
        max(abs(d_delta) - 0.25 / sqrt(i_delta), 0)
    observed_delta <- sqrt(2) * solved$gamma[[1L]]
    expect_identical(solved$status, "ok")
    expect_equal(observed_delta, expected_delta, tolerance = 1e-9)
    expect_lte(solved$kkt_residual, 1e-9)
})

test_that("zero-information exact edges are constrained to exact sharing", {
    n <- 40L
    signal_a <- seq(-1, 1, length.out = n)
    signal_b <- seq(-0.8, 1.2, length.out = n)
    x <- list(
        A = cbind(signal = signal_a, no_information = rep(1, n)),
        B = cbind(signal = signal_b, no_information = rep(2, n))
    )
    y <- list(
        A = 1 + 1.2 * signal_a + sin(seq_len(n)) * 0.05,
        B = 1 + 1.2 * signal_b + cos(seq_len(n)) * 0.05
    )
    scaling <- Pando:::.condition_ridge_scaling(x, 1e-8)
    fit <- Pando:::.condition_scheme_e_fit(
        x, y, scaling, min_residual_df = 1L, inference = FALSE
    )
    expect_identical(fit$status, "ok")
    expect_false(fit$contrast_identifiable[[2L]])
    expect_true(fit$shared_by_boundary[[2L]])
    expect_equal(diff(range(fit$beta[, 2L])), 0, tolerance = 0)
})

test_that("cell abundance remains in raw condition information", {
    base_x <- cbind(edge = rep(c(-1, -0.5, 0, 0.5, 1), 20L))
    x <- list(
        Large = base_x,
        Small = base_x[seq_len(nrow(base_x) / 4L), , drop = FALSE]
    )
    y <- list(
        Large = 0.7 * x$Large[, 1] + sin(seq_len(nrow(x$Large))) * 0.2,
        Small = 0.7 * x$Small[, 1] + sin(seq_len(nrow(x$Small))) * 0.2
    )
    scaling <- Pando:::.condition_ridge_scaling(x, 1e-8)
    fit <- Pando:::.condition_scheme_e_fit(
        x, y, scaling, min_residual_df = 1L, inference = FALSE
    )
    expect_identical(fit$status, "ok")
    ratio <- fit$raw_information["Large", "edge"] /
        fit$raw_information["Small", "edge"]
    expect_equal(ratio, 4, tolerance = 1e-10)
    expect_equal(unname(fit$condition_weight), c(1, 1))
})

test_that("three-condition Scheme E is invariant to condition ordering", {
    set.seed(71)
    n <- c(A = 70L, B = 45L, C = 30L)
    make_x <- function(nn, shift) {
        e1 <- rnorm(nn, shift, 1)
        e2 <- 0.45 * e1 + rnorm(nn, sd = 0.8)
        cbind(e1 = e1, e2 = e2)
    }
    x <- Map(make_x, n, c(A = 0, B = 0.2, C = -0.1))
    y <- lapply(names(x), function(condition) {
        b <- switch(condition,
                    A = c(1.0, -0.4), B = c(1.1, -0.4), C = c(0.7, -0.1))
        as.numeric(x[[condition]] %*% b + rnorm(n[[condition]], sd = 0.3))
    })
    names(y) <- names(x)
    scaling <- Pando:::.condition_ridge_scaling(x, 1e-8)
    fit <- Pando:::.condition_scheme_e_fit(
        x, y, scaling, min_residual_df = 1L, inference = FALSE,
        control = list(solver_tol = 1e-10, solver_max_iter = 10000L)
    )
    order <- c("C", "A", "B")
    x2 <- x[order]
    y2 <- y[order]
    scaling2 <- Pando:::.condition_ridge_scaling(x2, 1e-8)
    fit2 <- Pando:::.condition_scheme_e_fit(
        x2, y2, scaling2, min_residual_df = 1L, inference = FALSE,
        control = list(solver_tol = 1e-10, solver_max_iter = 10000L)
    )
    expect_identical(fit$status, "ok")
    expect_identical(fit2$status, "ok")
    expect_equal(fit$beta, fit2$beta[rownames(fit$beta), , drop = FALSE],
                 tolerance = 1e-7)
})

test_that("correlated multi-edge Scheme E reports a converged KKT residual", {
    set.seed(101)
    n <- 80L
    xa <- cbind(e1 = rnorm(n), e2 = rnorm(n))
    xa[, 2] <- 0.8 * xa[, 1] + 0.6 * xa[, 2]
    xb <- cbind(e1 = rnorm(n), e2 = rnorm(n))
    xb[, 2] <- 0.8 * xb[, 1] + 0.6 * xb[, 2]
    x <- list(A = xa, B = xb)
    y <- list(
        A = as.numeric(xa %*% c(1.0, -0.4) + rnorm(n, sd = 0.25)),
        B = as.numeric(xb %*% c(1.25, -0.1) + rnorm(n, sd = 0.25))
    )
    scaling <- Pando:::.condition_ridge_scaling(x, 1e-8)
    fit <- Pando:::.condition_scheme_e_fit(
        x, y, scaling, min_residual_df = 1L, inference = TRUE,
        control = list(solver_tol = 1e-9, solver_max_iter = 10000L)
    )
    expect_identical(fit$status, "ok")
    expect_identical(fit$solver_status, "ok")
    expect_true(is.finite(fit$kkt_residual))
    expect_lte(fit$kkt_residual, 1e-7)
    expect_true(all(is.finite(fit$beta)))
})