test_that("condition inference requires BH and a valid threshold", {
    expect_identical(Pando:::.condition_validate_adjust_method("BH"), "BH")
    expect_identical(Pando:::.condition_validate_adjust_method("bh"), "BH")
    expect_error(Pando:::.condition_validate_adjust_method("holm"),
                 'adjust_method = "BH"', fixed = TRUE)
    expect_equal(Pando:::.condition_validate_padj_threshold(0.05), 0.05)
    expect_error(Pando:::.condition_validate_padj_threshold(1))
    expect_error(Pando:::.condition_validate_padj_threshold(0))
})

.activity_fit_fixture <- function(padj_threshold = 0.05) {
    coefficient <- data.frame(
        edge_id = rep(c("E1", "E2", "E3"), 2L),
        target = "G",
        condition = rep(c("A", "B"), each = 3L),
        estimate = c(0.5, 0.4, 0.3, 0.6, 0.7, 0.8),
        penalty_effect = c(0.5, 0.4, 0.3, 0.6, 0.7, 0.8),
        pval = c(0.003, 0.1, 0.2, 0.3, 0.003, 0.2),
        padj = c(0.01, 0.20, 0.30, 0.50, 0.01, 0.40),
        inference_estimable = TRUE,
        condition_significant = c(TRUE, FALSE, FALSE, FALSE, TRUE, FALSE),
        fusion_component_id = "component1",
        shared_edge = FALSE,
        stringsAsFactors = FALSE
    )
    fit <- list(
        condition_levels = c("A", "B"),
        padj_threshold = padj_threshold,
        dictionary_support_table = data.frame(
            edge_id = c("E1", "E2", "E1", "E3"),
            source_type = c("global", "global", "condition", "condition"),
            condition = c(NA, NA, "A", "A"),
            peak_target_cor = c(0.22, 0.31, 0.2, 0.3),
            tf_target_cor = c(0.42, 0.51, 0.4, 0.5),
            stringsAsFactors = FALSE
        ),
        coefficients = coefficient,
        fit = data.frame(
            target = c("G", "G"), condition = c("A", "B"),
            fit_status = c("ok", "ok"), stringsAsFactors = FALSE
        )
    )
    class(fit) <- c("ConditionGRNFit", "list")
    fit
}

test_that("global and local candidate support remain provenance only", {
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
    expect_equal(fit$coefficients$penalty_effect,
                 fit$coefficients$estimate)
    expect_identical(
        fit$projection_policy,
        "any_condition_padj_exact_edge_union"
    )
})

test_that("one significant condition admits the exact edge in every condition", {
    fit <- Pando:::.condition_apply_activity_gate(.activity_fit_fixture())
    e1 <- fit$coefficients$edge_id == "E1"
    e2 <- fit$coefficients$edge_id == "E2"
    e3 <- fit$coefficients$edge_id == "E3"
    expect_true(all(fit$coefficients$active_in_regcompass[e1]))
    expect_true(all(fit$coefficients$active_in_regcompass[e2]))
    expect_false(any(fit$coefficients$active_in_regcompass[e3]))
    expect_identical(unique(fit$coefficients$supporting_conditions[e1]), "A")
    expect_identical(unique(fit$coefficients$supporting_conditions[e2]), "B")
    expect_identical(unique(fit$coefficients$supporting_conditions[e3]), "")
})

test_that("a nonsignificant condition keeps its continuous beta after union admission", {
    fit <- Pando:::.condition_apply_activity_gate(.activity_fit_fixture())
    row <- which(fit$coefficients$edge_id == "E1" &
                 fit$coefficients$condition == "B")
    expect_false(fit$coefficients$condition_significant[[row]])
    expect_true(fit$coefficients$active[[row]])
    expect_true(fit$coefficients$active_in_regcompass[[row]])
    expect_equal(fit$coefficients$penalty_effect[[row]], 0.6)
})

test_that("local-only dictionary edge remains fitted even when not handed off", {
    fit <- Pando:::.condition_apply_activity_gate(.activity_fit_fixture())
    row <- which(fit$coefficients$edge_id == "E3" &
                 fit$coefficients$condition == "B")
    expect_false(fit$coefficients$global_support[[row]])
    expect_false(fit$coefficients$local_support[[row]])
    expect_true(fit$coefficients$dictionary_support[[row]])
    expect_true(fit$coefficients$active[[row]])
    expect_false(fit$coefficients$active_in_regcompass[[row]])
    expect_equal(fit$coefficients$penalty_effect[[row]], 0.8)
})

test_that("local and global support retain their own correlations", {
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

test_that("default adjusted-P threshold remains 0.05", {
    default <- formals(Pando:::.pando_infer_condition_grn_multitask_ridge_one)$padj_threshold
    expect_equal(eval(default), 0.05)
})
