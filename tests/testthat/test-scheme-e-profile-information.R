test_that("E-star scalar update matches the closed form at z=0.25", {
    fit <- Pando:::.condition_E_star_fit(
        H = matrix(4, 1L, 1L),
        r = 1.2,
        information = 4,
        control = Pando:::.condition_E_star_control()
    )
    expect_identical(fit$status, "ok")
    expect_equal(fit$delta, 0.175, tolerance = 1e-8)
    expect_lt(fit$kkt_residual, 1e-7)
})

test_that("Q-orthogonal decomposition preserves contrasts and shared orthogonality", {
    q1 <- matrix(c(7, 1.2, 1.2, 3.5), 2L, 2L)
    q2 <- matrix(c(2.5, -0.4, -0.4, 5), 2L, 2L)
    q3 <- matrix(c(4.2, 0.7, 0.7, 6.1), 2L, 2L)
    Q <- Pando:::.condition_blockdiag(list(q1, q2, q3))
    A <- kronecker(matrix(1, 3L, 1L), diag(2L))
    D <- matrix(0, 4L, 6L)
    D[1, c(1, 3)] <- c(-1, 1)
    D[2, c(1, 5)] <- c(-1, 1)
    D[3, c(2, 4)] <- c(-1, 1)
    D[4, c(2, 6)] <- c(-1, 1)
    geometry <- Pando:::.condition_q_orthogonal_decomposition(
        Q, A, D, rank_tol = 1e-12
    )
    expect_equal(D %*% geometry$R, diag(4L), tolerance = 1e-10)
    expect_equal(crossprod(A, Q %*% geometry$R),
                 matrix(0, 2L, 4L), tolerance = 1e-10)
    expect_lt(geometry$dr_error, 1e-10)
    expect_lt(geometry$orthogonality_error, 1e-10)
})

test_that("profile information uses the full correlated contrast Hessian", {
    H <- matrix(c(5, 2, 2, 4), 2L, 2L)
    info <- Pando:::.condition_profile_coordinate_information(H, 1e-12)
    covariance <- solve(H)
    expect_equal(info$information, 1 / diag(covariance), tolerance = 1e-12)
    expect_identical(info$estimable, c(TRUE, TRUE))
})

.make_joint_scalar_fixture <- function(information, beta = 0.7, sigma2 = 1) {
    k <- length(information)
    x <- lapply(information, function(one) {
        a <- sqrt(one / 2)
        matrix(c(-a, 0, a), ncol = 1L, dimnames = list(NULL, "edge"))
    })
    names(x) <- paste0("C", seq_len(k))
    residual_df <- 3L * k - (k + 1L)
    residual_scale <- sqrt(residual_df * sigma2 / (6 * k))
    residual <- residual_scale * c(1, -2, 1)
    y <- lapply(x, function(one) beta * one[, 1] + residual)
    names(y) <- names(x)
    list(x = x, y = y)
}

test_that("fully shared JSE pools large and small condition information", {
    fixture <- .make_joint_scalar_fixture(c(100, 10, 10))
    fusion <- list(component = matrix(1L, 3L, 1L))
    scaling <- list(center = c(edge = 0), scale = c(edge = 1))
    joint <- Pando:::.condition_joint_inference_design(
        fixture$x, fixture$y, scaling, keep = 1L,
        fusion = fusion, rank_tol = 1e-12
    )
    rows <- Pando:::.condition_joint_inference_fit(
        joint, scaling, keep = 1L, k = 3L, p_full = 1L,
        conditions = names(fixture$x)
    )
    expect_identical(joint$status, "ok")
    expect_equal(joint$sigma2, 1, tolerance = 1e-10)
    expect_equal(as.numeric(rows$se[, 1]),
                 rep(1 / sqrt(120), 3L), tolerance = 1e-10)
    expect_equal(rows$estimate[, 1], rep(0.7, 3L), tolerance = 1e-10)
})

test_that("partial fusion JSE pools only the selected component", {
    fixture <- .make_joint_scalar_fixture(c(100, 10, 10))
    residual_df <- 9L - 5L
    residual_scale <- sqrt(residual_df / 18)
    residual <- residual_scale * c(1, -2, 1)
    fixture$y <- lapply(seq_along(fixture$x), function(i) {
        beta <- if (i < 3L) 0.7 else -0.2
        beta * fixture$x[[i]][, 1] + residual
    })
    names(fixture$y) <- names(fixture$x)
    fusion <- list(component = matrix(c(1L, 1L, 2L), 3L, 1L))
    scaling <- list(center = c(edge = 0), scale = c(edge = 1))
    joint <- Pando:::.condition_joint_inference_design(
        fixture$x, fixture$y, scaling, keep = 1L,
        fusion = fusion, rank_tol = 1e-12
    )
    rows <- Pando:::.condition_joint_inference_fit(
        joint, scaling, keep = 1L, k = 3L, p_full = 1L,
        conditions = names(fixture$x)
    )
    expect_equal(joint$sigma2, 1, tolerance = 1e-10)
    expect_equal(rows$se[1, 1], 1 / sqrt(110), tolerance = 1e-10)
    expect_equal(rows$se[2, 1], 1 / sqrt(110), tolerance = 1e-10)
    expect_equal(rows$se[3, 1], 1 / sqrt(10), tolerance = 1e-10)
})
