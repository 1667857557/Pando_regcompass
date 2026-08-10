test_that("condition GRN reliability uses target-condition sqrt R-squared", {
    fit <- list(
        schema_version = "pando_condition_grn_common_dictionary_v1",
        fit = data.frame(
            target = c("G1", "G1", "G2", "G2", "G3", "G4"),
            condition = c("A", "B", "A", "B", "A", "A"),
            rsq = c(0.25, 0.81, 0.04, NA, 1.44, -0.2),
            fit_status = c("ok", "ok", "rank_deficient", "ok", "ok", "ok"),
            stringsAsFactors = FALSE
        ),
        coefficients = data.frame(
            target = c("G1", "G1", "G2", "G2", "G3", "G4"),
            condition = c("A", "B", "A", "B", "A", "A"),
            estimable = TRUE,
            significant = c(TRUE, TRUE, TRUE, FALSE, TRUE, TRUE),
            stringsAsFactors = FALSE
        )
    )
    class(fit) <- c("ConditionGRNFit", "list")

    reliability <- condition_grn_reliability(fit)

    expect_equal(reliability$reliability, c(0.5, 0.9, NA, NA, 1, 0))
    expect_equal(reliability$n_active_edges, c(1L, 1L, 1L, 0L, 1L, 1L))
})

test_that("condition GRN reliability can use estimable edges", {
    fit <- list(
        schema_version = "pando_condition_grn_common_dictionary_v1",
        fit = data.frame(
            target = "G1", condition = "A", rsq = 0.36,
            fit_status = "ok", stringsAsFactors = FALSE
        ),
        coefficients = data.frame(
            target = "G1", condition = "A", estimable = TRUE,
            significant = FALSE, stringsAsFactors = FALSE
        )
    )
    class(fit) <- c("ConditionGRNFit", "list")

    expect_true(is.na(condition_grn_reliability(fit)$reliability))
    expect_equal(
        condition_grn_reliability(fit, significant_only = FALSE)$reliability,
        0.6
    )
})
