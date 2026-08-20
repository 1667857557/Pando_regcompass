test_that("public condition API exposes design controls but not z/CV controls", {
    args <- names(formals(Pando:::infer_condition_grn.GRNData))
    expect_true(all(c(
        "cell_type_col", "condition_col", "tf_cor", "peak_cor",
        "min_cells_per_condition", "rank_action", "min_residual_df",
        "parallel", "BPPARAM", "parallel_scope", "fallback_args"
    ) %in% args))
    expect_false(any(c(
        "ridge_control", "condition_ridge_control", "scheme_e_z", "z",
        "lambda_grid", "lambda_rule", "cv_folds"
    ) %in% args))
    expect_equal(eval(formals(Pando:::infer_condition_grn.GRNData)$tf_cor), 0.05)
    expect_equal(eval(formals(Pando:::infer_condition_grn.GRNData)$peak_cor), 0.05)
})

test_that("condition model is fixed E-star/JSE z=0.25", {
    expect_identical(
        Pando:::.condition_multitask_ridge_schema,
        "pando_condition_grn_Estar_jointse_v1"
    )
    expect_identical(
        Pando:::.condition_fit_engine,
        "condition_union_Estar_z025_jointse"
    )
    expect_equal(Pando:::.condition_E_star_z, 0.25)
    expect_identical(
        Pando:::.condition_E_star_penalty_family,
        "information_scaled_sparse_deviation"
    )
    control <- Pando:::.condition_E_star_control()
    expect_false(any(c(
        "scheme_e_z", "z", "lambda_grid", "lambda_rule", "cv_folds"
    ) %in% names(control)))
    expect_error(
        Pando:::.condition_E_star_control(list(scheme_e_z = 0.5)),
        "Unknown `condition_e_control` field", fixed = TRUE
    )
})

test_that("condition scaling is common and immune to condition mean shifts", {
    x <- list(
        A = cbind(edge = c(-1, 0, 1)),
        B = cbind(edge = c(99, 100, 101))
    )
    scaling <- Pando:::.condition_ridge_scaling(x, 1e-8)
    expected <- sqrt(mean(c(
        mean((x$A[, 1] - mean(x$A[, 1]))^2),
        mean((x$B[, 1] - mean(x$B[, 1]))^2)
    )))
    expect_equal(unname(scaling$scale[["edge"]]), expected)
    expect_identical(scaling$reference, "equal_condition_within_condition_rms")
    shifted <- x
    shifted$B[, 1] <- shifted$B[, 1] + 1e6
    expect_equal(
        Pando:::.condition_ridge_scaling(shifted, 1e-8)$scale,
        scaling$scale
    )
})

.contract_activity_fixture <- function() {
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
        fusion_component_id = rep(c("component1", "component1", "component1"), 2L),
        shared_edge = FALSE,
        stringsAsFactors = FALSE
    )
    fit <- list(
        condition_levels = c("A", "B"),
        padj_threshold = 0.05,
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

test_that("BH annotates conditions while RegCompass admission is exact-edge union", {
    gated <- Pando:::.condition_apply_activity_gate(.contract_activity_fixture())
    expect_identical(gated$coefficients$active, rep(TRUE, 6L))
    expect_identical(
        gated$coefficients$condition_significant,
        c(TRUE, FALSE, FALSE, FALSE, TRUE, FALSE)
    )
    expect_equal(gated$coefficients$penalty_effect,
                 gated$coefficients$estimate)
    expect_identical(
        gated$coefficients$active_in_regcompass,
        c(TRUE, TRUE, FALSE, TRUE, TRUE, FALSE)
    )
    expect_identical(
        gated$coefficients$supporting_conditions,
        c("A", "B", "", "A", "B", "")
    )
})

test_that("dictionary and handoff policies encode E-star/JSE union", {
    expect_identical(
        Pando:::.condition_fit_dictionary_policy,
        "global_and_condition_union_pando_correlation_supported_frozen_dictionary"
    )
    expect_identical(
        Pando:::.condition_significant_projection_policy,
        "any_condition_padj_exact_edge_union"
    )
})
