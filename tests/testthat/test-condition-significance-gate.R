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
            edge_id = c("E1", "E2", "E1", "E3"),
            source_type = c("global", "global", "condition", "condition"),
            condition = c(NA, NA, "A", "A"),
            peak_target_cor = c(0.22, 0.31, 0.2, 0.3),
            tf_target_cor = c(0.42, 0.51, 0.4, 0.5),
            stringsAsFactors = FALSE
        ),
        coefficients = data.frame(
            edge_id = rep(c("E1", "E2", "E3"), 2L),
            condition = rep(c("A", "B"), each = 3L),
            estimate = c(0.5, 0.4, 0.3, 0.6, 0.7, 0.8),
            estimable = TRUE,
            padj = c(0.01, 0.02, 0.03, 0.01, 0.01, 0.01),
            stringsAsFactors = FALSE
        )
    )
    class(fit) <- c("ConditionGRNFit", "list")
    fit
}

test_that("global and local Pando support remain provenance only", {
    fit <- Pando:::.condition_apply_activity_gate(.activity_fit_fixture())
    expect_identical(fit$coefficients$statistically_supported, rep(TRUE, 6L))
    expect_identical(
        fit$coefficients$global_support,
        c(TRUE, TRUE, FALSE, TRUE, TRUE, FALSE)
    )
    expect_identical(
        fit$coefficients$local_support,
        c(TRUE, FALSE, TRUE, FALSE, FALSE, FALSE)
    )
    expect_identical(fit$coefficients$dictionary_support, rep(TRUE, 6L))
    expect_identical(fit$coefficients$active, rep(TRUE, 6L))
    expect_identical(fit$coefficients$significant, fit$coefficients$active)
    expect_equal(
        fit$coefficients$penalty_effect,
        c(0.5, 0.4, 0.3, 0.6, 0.7, 0.8)
    )
    expect_identical(
        fit$projection_policy,
        "condition_bh_supported_common_dictionary_ridge_effects"
    )
})

test_that("local-only dictionary edge can be active in another condition", {
    fit <- Pando:::.condition_apply_activity_gate(.activity_fit_fixture())
    row <- which(fit$coefficients$edge_id == "E3" &
                 fit$coefficients$condition == "B")
    expect_true(fit$coefficients$statistically_supported[[row]])
    expect_false(fit$coefficients$global_support[[row]])
    expect_false(fit$coefficients$local_support[[row]])
    expect_true(fit$coefficients$dictionary_support[[row]])
    expect_true(fit$coefficients$active[[row]])
    expect_equal(fit$coefficients$penalty_effect[[row]], 0.8)
})

test_that("BH failure prevents activity regardless of candidate provenance", {
    fixture <- .activity_fit_fixture()
    fixture$coefficients$padj[[3L]] <- 0.2
    fit <- Pando:::.condition_apply_activity_gate(fixture)
    expect_false(fit$coefficients$statistically_supported[[3L]])
    expect_false(fit$coefficients$global_support[[3L]])
    expect_true(fit$coefficients$local_support[[3L]])
    expect_false(fit$coefficients$active[[3L]])
    expect_equal(fit$coefficients$penalty_effect[[3L]], 0)
})

test_that("local and global support store their own Pando correlations", {
    fit <- Pando:::.condition_apply_activity_gate(.activity_fit_fixture())
    row_a1 <- which(fit$coefficients$edge_id == "E1" &
                    fit$coefficients$condition == "A")
    row_b2 <- which(fit$coefficients$edge_id == "E2" &
                    fit$coefficients$condition == "B")
    row_b3 <- which(fit$coefficients$edge_id == "E3" &
                    fit$coefficients$condition == "B")
    expect_equal(fit$coefficients$peak_target_cor[[row_a1]], 0.2)
    expect_equal(fit$coefficients$tf_target_cor[[row_a1]], 0.4)
    expect_equal(fit$coefficients$global_peak_target_cor[[row_a1]], 0.22)
    expect_equal(fit$coefficients$global_tf_target_cor[[row_a1]], 0.42)
    expect_true(is.na(fit$coefficients$peak_target_cor[[row_b2]]))
    expect_true(is.na(fit$coefficients$tf_target_cor[[row_b2]]))
    expect_equal(fit$coefficients$global_peak_target_cor[[row_b2]], 0.31)
    expect_equal(fit$coefficients$global_tf_target_cor[[row_b2]], 0.51)
    expect_true(is.na(fit$coefficients$peak_target_cor[[row_b3]]))
    expect_true(is.na(fit$coefficients$global_peak_target_cor[[row_b3]]))
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
