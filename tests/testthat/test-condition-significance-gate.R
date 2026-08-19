test_that("condition diagnostic inference requires BH and a valid threshold", {
    expect_identical(Pando:::.condition_validate_adjust_method("BH"), "BH")
    expect_identical(Pando:::.condition_validate_adjust_method("bh"), "BH")
    expect_error(Pando:::.condition_validate_adjust_method("holm"),
                 'adjust_method = "BH"', fixed = TRUE)
    expect_equal(Pando:::.condition_validate_padj_threshold(0.05), 0.05)
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
            padj = c(0.01, 0.20, 0.03, 0.50, 0.01, 0.40),
            contrast_identifiable = TRUE,
            shared_by_boundary = FALSE,
            fused_by_penalty = FALSE,
            stringsAsFactors = FALSE
        )
    )
    class(fit) <- c("ConditionGRNFit", "list")
    fit
}

test_that("global and local Pando support remain provenance only", {
    fit <- Pando:::.condition_apply_activity_gate(.activity_fit_fixture())
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
    expect_identical(
        fit$coefficients$significant,
        fit$coefficients$statistically_supported
    )
    expect_equal(fit$coefficients$penalty_effect,
                 fit$coefficients$estimate)
    expect_identical(
        fit$projection_policy,
        "continuous_common_dictionary_scheme_e_effects"
    )
})

test_that("BH failure does not remove or zero a Scheme E dictionary edge", {
    fit <- Pando:::.condition_apply_activity_gate(.activity_fit_fixture())
    row <- which(fit$coefficients$edge_id == "E2" &
                 fit$coefficients$condition == "A")
    expect_false(fit$coefficients$statistically_supported[[row]])
    expect_true(fit$coefficients$active[[row]])
    expect_false(fit$coefficients$significant[[row]])
    expect_equal(fit$coefficients$penalty_effect[[row]],
                 fit$coefficients$estimate[[row]])
})

test_that("local-only dictionary edge remains present in every condition", {
    fit <- Pando:::.condition_apply_activity_gate(.activity_fit_fixture())
    row <- which(fit$coefficients$edge_id == "E3" &
                 fit$coefficients$condition == "B")
    expect_false(fit$coefficients$global_support[[row]])
    expect_false(fit$coefficients$local_support[[row]])
    expect_true(fit$coefficients$dictionary_support[[row]])
    expect_true(fit$coefficients$active[[row]])
    expect_equal(fit$coefficients$penalty_effect[[row]], 0.8)
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
    expect_equal(fit$coefficients$global_peak_target_cor[[row_b2]], 0.31)
    expect_true(is.na(fit$coefficients$peak_target_cor[[row_b3]]))
    expect_true(is.na(fit$coefficients$global_peak_target_cor[[row_b3]]))
})

test_that("default adjusted-P threshold remains diagnostic 0.05", {
    default <- formals(Pando:::.pando_infer_condition_grn_multitask_ridge_one)$padj_threshold
    expect_equal(eval(default), 0.05)
})