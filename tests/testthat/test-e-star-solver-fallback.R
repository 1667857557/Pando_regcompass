test_that("E-star falls back when accelerated updates exhaust max_iter", {
    fit <- Pando:::.condition_E_star_fit(
        H = matrix(4, 1L, 1L),
        r = 1.2,
        information = 4,
        control = Pando:::.condition_E_star_control(list(
            solver_max_iter = 1L
        ))
    )

    expect_identical(fit$status, "ok")
    expect_gt(fit$iterations, 1L)
    expect_equal(fit$delta, 0.175, tolerance = 1e-8)
    expect_lte(fit$kkt_residual, 1e-8)
})

test_that("E-star fallback converges for a near-collinear profile", {
    set.seed(177)
    x <- matrix(stats::rnorm(800), 40L, 20L)
    x[, 20L] <- x[, 1L] + 1e-7 * stats::rnorm(40L)
    fit <- Pando:::.condition_E_star_fit(
        H = crossprod(x),
        r = stats::rnorm(20L),
        information = rep(1, 20L),
        control = Pando:::.condition_E_star_control(list(
            solver_max_iter = 1L
        ))
    )

    expect_identical(fit$status, "ok")
    expect_lte(fit$kkt_residual, 1e-8)
    expect_true(all(is.finite(fit$delta)))
})

