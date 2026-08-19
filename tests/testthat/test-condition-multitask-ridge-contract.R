test_that("public condition API remains canonical and explicit", {
    args <- names(formals(Pando:::infer_condition_grn.GRNData))
    expect_true(all(c(
        "cell_type_col", "condition_col", "tf_cor", "peak_cor",
        "min_cells_per_condition", "rank_action", "min_residual_df",
        "parallel", "BPPARAM", "parallel_scope", "fallback_args"
    ) %in% args))
    expect_false("ridge_control" %in% args)
    expect_equal(eval(formals(Pando:::infer_condition_grn.GRNData)$tf_cor), 0.05)
    expect_equal(eval(formals(Pando:::infer_condition_grn.GRNData)$peak_cor), 0.05)
})

test_that("condition model uses fixed Scheme E z=0.25", {
    expect_identical(
        Pando:::.condition_multitask_ridge_schema,
        "pando_condition_grn_sparse_deviation_v4"
    )
    expect_identical(
        Pando:::.condition_fit_engine,
        "condition_union_scheme_e_exact_edge_z025"
    )
    control <- Pando:::.condition_ridge_control()
    expect_equal(control$scheme_e_z, 0.25)
    expect_false(any(c("lambda_grid", "lambda_rule", "cv_folds") %in% names(control)))
    expect_error(
        Pando:::.condition_ridge_control(list(scheme_e_z = 0.5)),
        "Unknown `ridge_control` field", fixed = TRUE
    )
    expect_error(
        Pando:::.condition_ridge_control(list(fusion_ratio = 1)),
        "Unknown `ridge_control` field", fixed = TRUE
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

test_that("BH is diagnostic and cannot zero a Scheme E coefficient", {
    fit <- list(
        padj_threshold = 0.05,
        dictionary_support_table = data.frame(
            edge_id = c("G||TF1||P1", "G||TF2||P2"),
            source_type = c("global", "condition"),
            condition = c(NA, "A"),
            peak_target_cor = c(0.2, 0.3),
            tf_target_cor = c(0.4, 0.5), stringsAsFactors = FALSE
        ),
        coefficients = data.frame(
            edge_id = rep(c("G||TF1||P1", "G||TF2||P2"), 2L),
            condition = rep(c("A", "B"), each = 2L),
            estimate = c(0.5, 0.4, 0.6, 0.7),
            padj = c(0.01, 0.20, 0.50, 0.01),
            contrast_identifiable = TRUE,
            shared_by_boundary = FALSE,
            fused_by_penalty = FALSE,
            stringsAsFactors = FALSE
        )
    )
    class(fit) <- c("ConditionGRNFit", "list")
    gated <- Pando:::.condition_apply_activity_gate(fit)
    expect_identical(gated$coefficients$active, rep(TRUE, 4L))
    expect_identical(
        gated$coefficients$statistically_supported,
        c(TRUE, FALSE, FALSE, TRUE)
    )
    expect_identical(gated$coefficients$significant,
                     gated$coefficients$statistically_supported)
    expect_equal(gated$coefficients$penalty_effect,
                 gated$coefficients$estimate)
})

test_that("dictionary and projection policies encode continuous Scheme E", {
    expect_identical(
        Pando:::.condition_fit_dictionary_policy,
        "global_and_condition_union_pando_correlation_supported_frozen_dictionary"
    )
    expect_identical(
        Pando:::.condition_significant_projection_policy,
        "continuous_common_dictionary_scheme_e_effects"
    )
})