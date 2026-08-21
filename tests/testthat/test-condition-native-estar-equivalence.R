test_that("native E-star solver matches the R reference", {
    set.seed(2755)
    for (dimension in c(1L, 7L, 23L)) {
        x <- matrix(stats::rnorm((dimension + 20L) * dimension),
                    dimension + 20L, dimension)
        H <- crossprod(x)
        expected <- seq(-0.03, 0.03, length.out = dimension)
        expected[abs(expected) < 1e-12] <- 0.01
        weights <- seq(0.5, 1.5, length.out = dimension)
        r <- as.numeric(H %*% expected) +
            Pando:::.condition_E_star_z * weights * sign(expected)
        information <- weights^2
        control <- Pando:::.condition_E_star_control()
        reference <- Pando:::.condition_E_star_fit_reference(
            H, r, information, control
        )
        native <- Pando:::.condition_E_star_fit(H, r, information, control)

        expect_identical(native$status, reference$status)
        expect_equal(native$delta, reference$delta, tolerance = 1e-8)
        expect_equal(native$objective, reference$objective, tolerance = 1e-9)
        expect_lte(native$kkt_residual, 1e-8)
    }
})
