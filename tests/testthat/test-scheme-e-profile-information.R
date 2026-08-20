test_that("E-star scalar update matches the closed form at z=0.25", {
    fit <- Pando:::.condition_E_star_fit(
        H = matrix(4, 1L, 1L), r = 1.2, information = 4,
        control = Pando:::.condition_E_star_control()
    )
    expect_identical(fit$status, "ok")
    expect_equal(fit$delta, 0.175, tolerance = 1e-8)
    expect_lt(fit$kkt_residual, 1e-7)
})

test_that("Q-orthogonal decomposition preserves contrasts", {
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
})

test_that("profile information uses the full correlated Hessian", {
    H <- matrix(c(5, 2, 2, 4), 2L, 2L)
    info <- Pando:::.condition_profile_coordinate_information(H, 1e-12)
    expect_equal(info$information, 1 / diag(solve(H)), tolerance = 1e-12)
    expect_identical(info$estimable, c(TRUE, TRUE))
})

test_that("E-star core no longer performs selected-fusion inference", {
    args <- names(formals(Pando:::.condition_scheme_e_fit))
    expect_false("inference" %in% args)
    ns <- asNamespace("Pando")
    expect_false(exists(".condition_joint_inference_design", ns, inherits = FALSE))
    expect_false(exists(".condition_joint_inference_fit", ns, inherits = FALSE))
    expect_false(exists(".condition_joint_contrast", ns, inherits = FALSE))
})
