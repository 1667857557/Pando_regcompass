test_that("condition ridge inference requires BH and a valid threshold", {
    expect_identical(Pando:::.condition_validate_adjust_method("BH"), "BH")
    expect_identical(Pando:::.condition_validate_adjust_method("bh"), "BH")
    expect_error(
        Pando:::.condition_validate_adjust_method("holm"),
        'adjust_method = "BH"',
        fixed = TRUE
    )
    expect_equal(Pando:::.condition_validate_padj_threshold(0.05), 0.05)
    expect_equal(Pando:::.condition_validate_padj_threshold(0.5), 0.5)
    expect_error(Pando:::.condition_validate_padj_threshold(1))
    expect_error(Pando:::.condition_validate_padj_threshold(0))
})

.activity_fit_fixture <- function(padj_threshold = 0.05) {
    fit <- list(
        padj_threshold = padj_threshold,
        dictionary_support_table = data.frame(
            edge_id = c("E1", "E2", "E1"),
            condition = c("A", "A", "B"),
            peak_target_cor = c(0.2, 0.3, 0.25),
            tf_target_cor = c(0.4, 0.5, 0.35),
            stringsAsFactors = FALSE
        ),
        coefficients = data.frame(
            edge_id = rep(c("E1", "E2"), 2L),
            condition = rep(c("A", "B"), each = 2L),
            estimate = c(0.5, 0.4, 0.6, 0.7),
            estimable = TRUE,
            padj = c(0.01, 0.02, 0.01, 0.01),
            stringsAsFactors = FALSE
        )
    )
    class(fit) <- c("ConditionGRNFit", "list")
    fit
}

test_that("statistical support and local Pando support remain distinct", {
    fit <- Pando:::.condition_apply_activity_gate(.activity_fit_fixture())
    expect_identical(fit$coefficients$statistically_supported, rep(TRUE, 4L))
    expect_identical(fit$coefficients$local_support,
                     c(TRUE, TRUE, TRUE, FALSE))
    expect_identical(fit$coefficients$active,
                     c(TRUE, TRUE, TRUE, FALSE))
    expect_identical(fit$coefficients$significant, fit$coefficients$active)
    expect_equal(fit$coefficients$penalty_effect, c(0.5, 0.4, 0.6, 0))
    expect_identical(
        fit$projection_policy,
        "active_condition_pando_support_and_bh_ridge_effects"
    )
})

test_that("BH failure prevents activity even with local Pando support", {
    fixture <- .activity_fit_fixture()
    fixture$coefficients$padj[[2L]] <- 0.2
    fit <- Pando:::.condition_apply_activity_gate(fixture)
    expect_false(fit$coefficients$statistically_supported[[2L]])
    expect_true(fit$coefficients$local_support[[2L]])
    expect_false(fit$coefficients$active[[2L]])
    expect_equal(fit$coefficients$penalty_effect[[2L]], 0)
})

test_that("local support stores the condition-specific Pando correlations", {
    fit <- Pando:::.condition_apply_activity_gate(.activity_fit_fixture())
    row_a1 <- which(fit$coefficients$edge_id == "E1" &
                    fit$coefficients$condition == "A")
    row_b2 <- which(fit$coefficients$edge_id == "E2" &
                    fit$coefficients$condition == "B")
    expect_equal(fit$coefficients$peak_target_cor[[row_a1]], 0.2)
    expect_equal(fit$coefficients$tf_target_cor[[row_a1]], 0.4)
    expect_true(is.na(fit$coefficients$peak_target_cor[[row_b2]]))
    expect_true(is.na(fit$coefficients$tf_target_cor[[row_b2]]))
})

test_that("obsolete significant-union refit helpers are absent", {
    ns <- asNamespace("Pando")
    expect_false(exists(".condition_dictionary_screen", ns, inherits = FALSE))
    expect_false(exists(".condition_subset_dictionary", ns, inherits = FALSE))
    expect_false(exists(".condition_ridge_refit_contract", ns, inherits = FALSE))
    expect_false(exists(
        ".condition_ridge_refit_contract_one_pass", ns, inherits = FALSE
    ))
    expect_true(is.function(Pando:::.condition_ridge_fit_contract))
    expect_true(is.function(Pando:::.condition_ridge_fit_contract_one_pass))
})

test_that("default adjusted-P threshold remains 0.05", {
    default <- formals(Pando:::.pando_infer_condition_grn_multitask_ridge_one)$padj_threshold
    expect_equal(eval(default), 0.05)
})
