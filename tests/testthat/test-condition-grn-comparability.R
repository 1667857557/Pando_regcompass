test_that("E-star reproduces the two-condition scalar closed form", {
    solved <- Pando:::.condition_E_star_fit(
        H = matrix(4, 1L, 1L), r = 1.2, information = 4,
        control = Pando:::.condition_E_star_control(list(
            solver_tol = 1e-12, solver_max_iter = 10000L
        ))
    )
    expected <- 0.3 - 0.25 / sqrt(4)
    expect_identical(solved$status, "ok")
    expect_equal(solved$delta[[1L]], expected, tolerance = 1e-9)
    expect_lte(solved$kkt_residual, 1e-9)
})

test_that("zero-information exact edges remain production boundary sharing", {
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
        x, y, scaling, min_residual_df = 1L,
        reference_condition = "A"
    )
    expect_identical(fit$status, "ok")
    expect_false(fit$contrast_identifiable[[2L]])
    expect_true(fit$shared_by_boundary[[2L]])
    expect_equal(diff(range(fit$beta[, 2L])), 0, tolerance = 0)
})

test_that("three-condition E-star is invariant to row order with fixed reference", {
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
        x, y, scaling, min_residual_df = 1L,
        control = list(solver_tol = 1e-10, solver_max_iter = 10000L),
        reference_condition = "A"
    )
    order <- c("C", "A", "B")
    scaling2 <- Pando:::.condition_ridge_scaling(x[order], 1e-8)
    fit2 <- Pando:::.condition_scheme_e_fit(
        x[order], y[order], scaling2, min_residual_df = 1L,
        control = list(solver_tol = 1e-10, solver_max_iter = 10000L),
        reference_condition = "A"
    )
    expect_identical(fit$status, "ok")
    expect_identical(fit2$status, "ok")
    expect_equal(fit$beta, fit2$beta[rownames(fit$beta), , drop = FALSE],
                 tolerance = 1e-7)
})

test_that("condition subgraphs use one common exact-edge topology", {
    fit <- list(
        schema_version = Pando:::.condition_common_dictionary_schema,
        condition_levels = c("A", "B"),
        coefficients = data.frame(
            edge_id = rep(c("g||t1||p1", "g||t2||p2"), each = 2L),
            condition = rep(c("A", "B"), times = 2L),
            significant = c(TRUE, TRUE, FALSE, FALSE),
            edge_supported = c(TRUE, TRUE, FALSE, FALSE),
            active_in_regcompass = c(TRUE, TRUE, FALSE, FALSE),
            penalty_effect = c(1.0, 0.8, 0.3, 0.2),
            stringsAsFactors = FALSE
        )
    )
    class(fit) <- c("ConditionGRNFit", "list")
    a <- condition_grn_subgraph(fit, "A", significant_only = TRUE)
    b <- condition_grn_subgraph(fit, "B", significant_only = TRUE)
    expect_identical(as.character(a$edge_id), "g||t1||p1")
    expect_identical(as.character(b$edge_id), "g||t1||p1")
    expect_equal(a$penalty_effect[[1L]], 1.0)
    expect_equal(b$penalty_effect[[1L]], 0.8)
})
