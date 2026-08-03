test_that("separate fixed-dictionary fits equal a block-diagonal model", {
    set.seed(14)
    n1 <- 80L
    n2 <- 95L
    p <- 3L
    x1 <- matrix(rnorm(n1 * p), nrow = n1)
    x2 <- matrix(rnorm(n2 * p), nrow = n2)
    beta1 <- c(0.8, -1.1, 0.3)
    beta2 <- c(-0.4, 0.6, 1.2)
    y1 <- 0.7 + x1 %*% beta1 + rnorm(n1, sd = 0.2)
    y2 <- -0.5 + x2 %*% beta2 + rnorm(n2, sd = 0.2)
    terms <- paste0("edge_", seq_len(p))

    fit1 <- Pando:::.condition_fit_target_matrix(
        y1, x1, terms, rank_action = "error"
    )
    fit2 <- Pando:::.condition_fit_target_matrix(
        y2, x2, terms, rank_action = "error"
    )

    block <- as.matrix(Matrix::bdiag(cbind(1, x1), cbind(1, x2)))
    block_coef <- stats::qr.coef(stats::qr(block), c(y1, y2))
    expect_equal(fit1$coefs$estimate, block_coef[2:(p + 1L)],
                 tolerance = 1e-10)
    offset <- p + 1L
    expect_equal(fit2$coefs$estimate,
                 block_coef[(offset + 2L):(2L * offset)],
                 tolerance = 1e-10)
})

test_that("edge order changes neither aligned estimates nor directions", {
    set.seed(15)
    x <- matrix(rnorm(500), nrow = 100, ncol = 5)
    y <- 1 + x %*% c(0.2, -0.4, 0.7, 0, 1.1) + rnorm(100, sd = 0.1)
    terms <- paste0("edge_", 1:5)
    fit1 <- Pando:::.condition_fit_target_matrix(
        y, x, terms, rank_action = "error"
    )
    order2 <- c(5, 2, 4, 1, 3)
    fit2 <- Pando:::.condition_fit_target_matrix(
        y, x[, order2, drop = FALSE], terms[order2], rank_action = "error"
    )
    aligned <- fit2$coefs$estimate[match(terms, fit2$coefs$term)]
    expect_equal(fit1$coefs$estimate, aligned, tolerance = 1e-10)
    expect_identical(sign(fit1$coefs$estimate), sign(aligned))
})

test_that("union refit restores a direction reversed by omitted variables", {
    set.seed(16)
    n <- 1000L
    x1 <- rnorm(n)
    x2 <- 0.8 * x1 + sqrt(1 - 0.8^2) * rnorm(n)
    y <- x1 - 2 * x2 + rnorm(n, sd = 0.05)

    reduced <- stats::coef(stats::lm(y ~ x1))[["x1"]]
    full <- Pando:::.condition_fit_target_matrix(
        y, cbind(x1, x2), c("edge_x1", "edge_x2"),
        rank_action = "error"
    )
    expect_lt(reduced, 0)
    expect_gt(full$coefs$estimate[[1L]], 0)
    expect_lt(full$coefs$estimate[[2L]], 0)
})
