test_that("condition ridge padj threshold is capped at 0.1", {
    expect_equal(Pando:::.condition_validate_padj_threshold(0.05), 0.05)
    expect_equal(Pando:::.condition_validate_padj_threshold(0.1), 0.1)
    expect_error(
        Pando:::.condition_validate_padj_threshold(0.10001),
        "\\(0, 0.1\\]"
    )
    expect_error(Pando:::.condition_validate_padj_threshold(0))
})

test_that("condition ridge penalty_effect follows BH significance gate", {
    fit <- list(
        padj_threshold = 0.1,
        coefficients = data.frame(
            estimate = c(0.5, 0.4, 0.3, NA_real_),
            estimable = c(TRUE, TRUE, TRUE, FALSE),
            padj = c(0.049, 0.08, 0.11, NA_real_),
            significant = c(FALSE, FALSE, FALSE, FALSE),
            penalty_effect = c(0, 0, 0, 0),
            stringsAsFactors = FALSE
        )
    )
    class(fit) <- c("ConditionGRNFit", "list")
    gated <- Pando:::.condition_apply_significance_gate(fit)

    expect_identical(gated$coefficients$significant,
                     c(TRUE, TRUE, FALSE, FALSE))
    expect_equal(gated$coefficients$penalty_effect,
                 c(0.5, 0.4, 0, 0))
    expect_identical(gated$projection_policy,
                     "padj_significant_ridge_effects")
})

test_that("default threshold remains 0.05 while 0.1 is allowed", {
    default <- formals(Pando:::.pando_infer_condition_grn_one)$padj_threshold
    expect_equal(eval(default), 0.05)

    body_text <- paste(
        deparse(body(Pando:::.pando_infer_condition_grn_one)),
        collapse = "\n"
    )
    expect_match(body_text, ".condition_validate_padj_threshold", fixed = TRUE)
})
