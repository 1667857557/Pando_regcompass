test_that("insufficient residual degrees of freedom are not estimable effects", {
    set.seed(13)
    predictor <- matrix(
        rnorm(12), nrow = 4,
        dimnames = list(NULL, paste0("x", 1:3))
    )
    response <- as.numeric(1 + predictor %*% c(1, -1, 0.5))
    fit <- Pando:::.condition_fit_target_matrix(
        response = response,
        predictor = predictor,
        terms = paste0("edge_", 1:3),
        rank_action = "mark",
        min_residual_df = 1L
    )
    expect_identical(fit$gof$fit_status[[1L]], "insufficient_df")
    expect_true(all(!fit$coefs$estimable))
    expect_true(all(is.na(fit$coefs$estimate)))
    expect_true(all(is.na(fit$coefs$pval)))
})
