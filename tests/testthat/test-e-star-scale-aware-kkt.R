test_that("E-star KKT residual is invariant to objective scaling", {
    delta <- c(0.2, 0)
    baseline_scale <- 1e6
    gradient <- c(-0.5 + 1e-9, 0.2) * baseline_scale
    weights <- c(2, 1) * baseline_scale
    active <- c(TRUE, TRUE)
    baseline <- Pando:::.condition_E_star_kkt(
        delta, gradient, weights, active, z = 0.25
    )
    objective_scale <- 1e12
    scaled <- Pando:::.condition_E_star_kkt(
        delta,
        gradient * objective_scale,
        weights * objective_scale,
        active,
        z = 0.25
    )

    expect_equal(scaled, baseline, tolerance = 1e-12)
})

test_that("E-star converges when a moderate system has a huge objective scale", {
    set.seed(177)
    x <- matrix(stats::rnorm(280), 40L, 7L)
    H <- crossprod(x)
    expected <- c(0.02, -0.01, 0.03, -0.015, 0.01, -0.025, 0.005)
    weights <- seq(0.5, 1.1, length.out = length(expected))
    r <- as.numeric(H %*% expected) +
        Pando:::.condition_E_star_z * weights * sign(expected)
    objective_scale <- 1e18

    fit <- Pando:::.condition_E_star_fit(
        H = H * objective_scale,
        r = r * objective_scale,
        information = (weights * objective_scale)^2,
        control = Pando:::.condition_E_star_control()
    )

    expect_identical(fit$status, "ok")
    expect_lte(fit$kkt_residual, 1e-8)
    expect_equal(fit$delta, expected, tolerance = 1e-8)
})

test_that("Q-orthogonality diagnostics are invariant to information scaling", {
    Q <- matrix(c(
        4, 1, 0, 0,
        1, 3, 0, 0,
        0, 0, 5, 2,
        0, 0, 2, 6
    ), 4L, 4L, byrow = TRUE)
    A <- matrix(c(1, 0, 1, 0, 0, 1, 0, 1), 4L, 2L)
    D <- matrix(c(-1, 0, 1, 0, 0, -1, 0, 1), 2L, 4L, byrow = TRUE)

    baseline <- Pando:::.condition_q_orthogonal_decomposition(
        Q, A, D, rank_tol = 1e-10
    )
    scaled <- Pando:::.condition_q_orthogonal_decomposition(
        Q * 1e18, A, D, rank_tol = 1e-10
    )

    expect_lte(baseline$orthogonality_error, 1e-12)
    expect_lte(scaled$orthogonality_error, 1e-12)
    expect_equal(
        scaled$orthogonality_error,
        baseline$orthogonality_error,
        tolerance = 1e-12
    )
    expect_equal(scaled$R, baseline$R, tolerance = 1e-12)
})

test_that("adaptive restart converges before the fixed E-star budget", {
    set.seed(177)
    x <- matrix(stats::rnorm(920), 40L, 23L)
    x[, 23L] <- x[, 1L] + 0.05 * stats::rnorm(40L)
    H <- crossprod(x)
    expected <- seq(-0.03, 0.03, length.out = 23L)
    expected[12L] <- 0.004
    weights <- seq(0.5, 1.5, length.out = 23L)
    r <- as.numeric(H %*% expected) +
        Pando:::.condition_E_star_z * weights * sign(expected)

    fit <- Pando:::.condition_E_star_fit(
        H = H, r = r, information = weights^2,
        control = Pando:::.condition_E_star_control(list(
            solver_max_iter = 5000L
        ))
    )

    expect_identical(fit$status, "ok")
    expect_lt(fit$iterations, 5000L)
    expect_lte(fit$kkt_residual, 1e-8)
    expect_equal(fit$delta, expected, tolerance = 1e-7)
})

