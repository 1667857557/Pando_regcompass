test_that("scaled Q geometry returns the dense shared MLE", {
    q1 <- matrix(c(8, 1.5, 1.5, 4), 2L, 2L)
    q2 <- matrix(c(3, -0.5, -0.5, 6), 2L, 2L)
    Q <- Pando:::.condition_blockdiag(list(q1, q2))
    A <- kronecker(matrix(1, 2L, 1L), diag(2L))
    D <- matrix(c(-1, 0, 1, 0, 0, -1, 0, 1), 2L, 4L,
                byrow = TRUE)
    h <- c(1.2, -0.7, 0.8, 1.1)

    geometry <- Pando:::.condition_q_orthogonal_decomposition(
        Q, A, D, rank_tol = 1e-12
    )
    obtained <- as.numeric(
        geometry$shared_inverse %*%
            (crossprod(A, h) / geometry$q_scale)
    )
    expected <- as.numeric(
        Pando:::.condition_symmetric_pinv(
            crossprod(A, Q %*% A), 1e-12
        )$inverse %*% crossprod(A, h)
    )

    expect_equal(obtained, expected, tolerance = 1e-12)
})

test_that("shared MLE and fitted beta are invariant to common information scaling", {
    Q <- matrix(c(
        4, 1, 0, 0,
        1, 3, 0, 0,
        0, 0, 5, 2,
        0, 0, 2, 6
    ), 4L, 4L, byrow = TRUE)
    A <- matrix(c(1, 0, 1, 0, 0, 1, 0, 1), 4L, 2L)
    D <- matrix(c(-1, 0, 1, 0, 0, -1, 0, 1), 2L, 4L,
                byrow = TRUE)
    h <- c(0.9, -0.3, 1.4, 0.6)
    delta <- c(0.12, -0.08)

    fit_at_scale <- function(multiplier) {
        geometry <- Pando:::.condition_q_orthogonal_decomposition(
            Q * multiplier, A, D, rank_tol = 1e-12
        )
        mu <- as.numeric(
            geometry$shared_inverse %*%
                (crossprod(A, h * multiplier) / geometry$q_scale)
        )
        list(mu = mu, beta = as.numeric(A %*% mu + geometry$R %*% delta),
             R = geometry$R)
    }
    baseline <- fit_at_scale(1)
    for (multiplier in c(1e6, 1e12, 1e18)) {
        scaled <- fit_at_scale(multiplier)
        expect_equal(scaled$mu, baseline$mu, tolerance = 1e-11)
        expect_equal(scaled$beta, baseline$beta, tolerance = 1e-11)
        expect_equal(scaled$R, baseline$R, tolerance = 1e-11)
    }
})

test_that("scaling Q alone changes the shared MLE inversely", {
    Q <- diag(c(2, 3, 5, 7))
    A <- matrix(c(1, 0, 1, 0, 0, 1, 0, 1), 4L, 2L)
    D <- matrix(c(-1, 0, 1, 0, 0, -1, 0, 1), 2L, 4L,
                byrow = TRUE)
    h <- c(1, 2, 3, 4)
    shared <- function(multiplier) {
        geometry <- Pando:::.condition_q_orthogonal_decomposition(
            Q * multiplier, A, D, rank_tol = 1e-12
        )
        as.numeric(
            geometry$shared_inverse %*%
                (crossprod(A, h) / geometry$q_scale)
        )
    }
    expect_equal(shared(1e6), shared(1) / 1e6, tolerance = 1e-12)
})
