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
    expect_equal(eval(formals(Pando:::discover_grn_edges)$tf_cor), 0.05)
    expect_equal(eval(formals(Pando:::discover_grn_edges)$peak_cor), 0.05)
})

test_that("condition runtime keeps thresholds as adjustable inputs", {
    defaults <- formals(Pando:::.pando_infer_condition_grn_multitask_ridge_one)
    expect_equal(eval(defaults$tf_cor), 0.05)
    expect_equal(eval(defaults$peak_cor), 0.05)
    runtime_text <- paste(
        deparse(body(Pando:::.pando_infer_condition_grn_multitask_ridge_one)),
        collapse = "\n"
    )
    expect_match(runtime_text, "tf_cor = tf_cor", fixed = TRUE)
    expect_match(runtime_text, "peak_cor = peak_cor", fixed = TRUE)
    expect_match(runtime_text, "candidate_tf_cor = as.numeric(tf_cor)", fixed = TRUE)
    expect_match(runtime_text, "candidate_peak_cor = as.numeric(peak_cor)", fixed = TRUE)
})

test_that("condition runtime uses global-plus-condition union and one ridge fit", {
    runtime_text <- paste(
        deparse(body(Pando:::.pando_infer_condition_grn_multitask_ridge_one)),
        collapse = "\n"
    )
    expect_match(runtime_text, "global_edges <-", fixed = TRUE)
    expect_match(runtime_text, "condition_edges <- lapply", fixed = TRUE)
    expect_match(runtime_text, "global_edges = global_edges", fixed = TRUE)
    expect_match(runtime_text, ".condition_ridge_fit_contract", fixed = TRUE)
    expect_false(grepl(".condition_ridge_refit_contract", runtime_text, fixed = TRUE))

    fit_text <- paste(
        deparse(body(Pando:::.condition_ridge_fit_contract)), collapse = "\n"
    )
    expect_match(fit_text, ".condition_ridge_fit_contract_one_pass", fixed = TRUE)
    expect_false(grepl(".condition_dictionary_screen", fit_text, fixed = TRUE))
    expect_false(grepl(".condition_subset_dictionary", fit_text, fixed = TRUE))
    expect_false(grepl("ridge_preliminary", fit_text, fixed = TRUE))
    expect_false(grepl("ridge_final", fit_text, fixed = TRUE))
})

test_that("condition model uses the revised no-fusion schema", {
    expect_identical(
        Pando:::.condition_multitask_ridge_schema,
        "pando_condition_grn_multitask_ridge_v3"
    )
    expect_identical(
        Pando:::.condition_fit_engine,
        "condition_union_single_no_fusion_common_lambda_ridge"
    )
})

test_that("ridge controls contain no fusion parameter", {
    control <- Pando:::.condition_ridge_control()
    expect_true(all(control$lambda_grid > 0))
    expect_identical(control$lambda_rule, "1se")
    expect_false("fusion_ratio" %in% names(control))
    expect_equal(control$cv_folds, 5L)
    expect_equal(control$seed, 1L)
    expect_error(
        Pando:::.condition_ridge_control(list(fusion_ratio = 1)),
        "Unknown `ridge_control` field",
        fixed = TRUE
    )
})

test_that("condition ridge penalty is block-wise ordinary ridge", {
    k <- 3L
    p <- 4L
    expect_equal(Pando:::.condition_ridge_penalty(k, p), diag(k * p))
    expect_false("fusion_ratio" %in% names(formals(Pando:::.condition_ridge_fit)))
    expect_false("fusion_ratio" %in% names(formals(Pando:::.condition_ridge_penalty)))
})

test_that("scaling matches the condition-intercept loss geometry", {
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
    expect_identical(
        scaling$reference,
        "equal_condition_within_condition_rms"
    )
    shifted <- x
    shifted$B[, 1] <- shifted$B[, 1] + 1e6
    expect_equal(
        Pando:::.condition_ridge_scaling(shifted, 1e-8)$scale,
        scaling$scale
    )
})

test_that("candidate provenance never vetoes a BH-supported dictionary edge", {
    fit <- list(
        padj_threshold = 0.05,
        dictionary_support_table = data.frame(
            edge_id = c("G||TF1||P1", "G||TF2||P2", "G||TF3||P3"),
            source_type = c("global", "global", "condition"),
            condition = c(NA, NA, "A"),
            peak_target_cor = c(0.2, 0.3, 0.25),
            tf_target_cor = c(0.4, 0.5, 0.35),
            stringsAsFactors = FALSE
        ),
        coefficients = data.frame(
            edge_id = rep(c("G||TF1||P1", "G||TF2||P2", "G||TF3||P3"), 2L),
            condition = rep(c("A", "B"), each = 3L),
            estimate = c(0.5, 0.4, 0.3, 0.6, 0.7, 0.8),
            estimable = TRUE,
            padj = c(0.01, 0.02, 0.03, 0.01, 0.01, 0.01),
            stringsAsFactors = FALSE
        )
    )
    class(fit) <- c("ConditionGRNFit", "list")
    gated <- Pando:::.condition_apply_activity_gate(fit)
    expect_identical(gated$coefficients$statistically_supported,
                     rep(TRUE, 6L))
    expect_identical(
        gated$coefficients$global_support,
        c(TRUE, TRUE, FALSE, TRUE, TRUE, FALSE)
    )
    expect_identical(
        gated$coefficients$local_support,
        c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE)
    )
    expect_identical(gated$coefficients$active, rep(TRUE, 6L))
    expect_equal(gated$coefficients$penalty_effect,
                 c(0.5, 0.4, 0.3, 0.6, 0.7, 0.8))
})

test_that("fit dictionary policy is global-plus-condition and frozen", {
    expect_identical(
        Pando:::.condition_fit_dictionary_policy,
        "global_and_condition_union_pando_correlation_supported_frozen_dictionary"
    )
    expect_identical(
        Pando:::.condition_significant_projection_policy,
        "condition_bh_supported_common_dictionary_ridge_effects"
    )
})
