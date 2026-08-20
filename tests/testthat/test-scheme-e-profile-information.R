test_that("two-condition Scheme-E profile information has the exact delta normalization", {
    qa <- matrix(c(7.0, 1.2, 1.2, 3.5), 2L, 2L)
    qb <- matrix(c(2.5, -0.4, -0.4, 5.0), 2L, 2L)
    basis <- Pando:::.condition_contrast_basis(2L)
    h_mu_mu <- qa + qb
    h_mu_gamma <- basis[1L, 1L] * qa + basis[2L, 1L] * qb
    h_gamma_gamma <- basis[1L, 1L]^2 * qa + basis[2L, 1L]^2 * qb
    profile_gamma <- h_gamma_gamma -
        crossprod(h_mu_gamma, solve(h_mu_mu, h_mu_gamma))
    i_delta <- solve(solve(qa) + solve(qb))

    expect_equal(profile_gamma / 2, i_delta, tolerance = 1e-12)

    raw_scale <- c(2.0, 0.5)
    fit_stub <- list(
        beta = matrix(0, nrow = 2L, ncol = 2L),
        profile_information = diag(profile_gamma) * raw_scale^2
    )
    exported <- Pando:::.condition_profile_information_export(fit_stub)
    expect_equal(exported, diag(i_delta) * raw_scale^2, tolerance = 1e-12)
})
