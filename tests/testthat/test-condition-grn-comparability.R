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
        x, y, scaling, min_residual_df = 1L, inference = FALSE,
        reference_condition = "A"
    )
    expect_identical(fit$status, "ok")
    expect_false(fit$contrast_identifiable[[2L]])
    expect_true(fit$shared_by_boundary[[2L]])
    expect_equal(diff(range(fit$beta[, 2L])), 0, tolerance = 0)
})

test_that("reference-unidentifiable coordinates keep a full boundary tree", {
    n <- 60L
    b <- seq(-1, 1, length.out = n)
    c <- seq(-0.8, 1.2, length.out = n)
    x <- list(
        A = cbind(edge = rep(1, n)),
        B = cbind(edge = b),
        C = cbind(edge = c)
    )
    y <- list(
        A = rep(0.3, n) + sin(seq_len(n)) * 0.02,
        B = 0.7 * b + sin(seq_len(n)) * 0.02,
        C = 1.1 * c + cos(seq_len(n)) * 0.02
    )
    scaling <- Pando:::.condition_ridge_scaling(x, 1e-8)
    fit <- Pando:::.condition_scheme_e_fit(
        x, y, scaling, min_residual_df = 1L, inference = TRUE,
        reference_condition = "A",
        control = list(solver_tol = 1e-10, solver_max_iter = 10000L)
    )
    expect_identical(fit$status, "ok")
    tree <- fit$contrast_tree
    expect_equal(nrow(tree), 2L)
    expect_equal(sum(tree$contrast_identifiable %in% TRUE), 1L)
    expect_true(any(tree$shared_by_boundary %in% TRUE))
    expect_equal(fit$dr_error, 0, tolerance = 1e-9)
    boundary <- tree[tree$shared_by_boundary %in% TRUE, , drop = FALSE]
    expect_equal(boundary$delta_standardized, 0, tolerance = 0)
    for (i in seq_len(nrow(boundary))) {
        a <- boundary$condition_a[[i]]
        bb <- boundary$condition_b[[i]]
        expect_equal(fit$beta[a, "edge"], fit$beta[bb, "edge"], tolerance = 1e-9)
    }
    expect_true(any(
        tree$contrast_identifiable %in% TRUE &
        !tree$condition_a %in% "A" & !tree$condition_b %in% "A"
    ))
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
        x, y, scaling, min_residual_df = 1L, inference = FALSE,
        reference_condition = "Large"
    )
    expect_identical(fit$status, "ok")
    ratio <- fit$raw_information["Large", "edge"] /
        fit$raw_information["Small", "edge"]
    expect_equal(ratio, 4, tolerance = 1e-10)
    expect_equal(unname(fit$condition_weight), c(1, 1))
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
        x, y, scaling, min_residual_df = 1L, inference = FALSE,
        control = list(solver_tol = 1e-10, solver_max_iter = 10000L),
        reference_condition = "A"
    )
    order <- c("C", "A", "B")
    x2 <- x[order]
    y2 <- y[order]
    scaling2 <- Pando:::.condition_ridge_scaling(x2, 1e-8)
    fit2 <- Pando:::.condition_scheme_e_fit(
        x2, y2, scaling2, min_residual_df = 1L, inference = FALSE,
        control = list(solver_tol = 1e-10, solver_max_iter = 10000L),
        reference_condition = "A"
    )
    expect_identical(fit$status, "ok")
    expect_identical(fit2$status, "ok")
    expect_equal(fit$beta, fit2$beta[rownames(fit$beta), , drop = FALSE],
                 tolerance = 1e-7)
})

test_that("correlated multi-edge E-star reports a converged KKT residual", {
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
        control = list(solver_tol = 1e-9, solver_max_iter = 10000L),
        reference_condition = "A"
    )
    expect_identical(fit$status, "ok")
    expect_identical(fit$solver_status, "ok")
    expect_true(is.finite(fit$kkt_residual))
    expect_lte(fit$kkt_residual, 1e-7)
    expect_true(all(is.finite(fit$beta)))
    expect_identical(fit$inference_schema,
                     "scheme_e_fusion_component_joint_refit_v1")
})
