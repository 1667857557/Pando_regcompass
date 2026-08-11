test_that("public condition API remains canonical and explicit", {
    args <- names(formals(Pando:::infer_condition_grn.GRNData))
    expect_true(all(c(
        "cell_type_col", "condition_col", "tf_cor", "peak_cor",
        "min_cells_per_condition", "rank_action", "min_residual_df",
        "parallel", "BPPARAM", "parallel_scope", "fallback_args"
    ) %in% args))
    expect_false("ridge_control" %in% args)
})

test_that("multi-condition internal path goes directly from exact union to ridge", {
    body_text <- paste(
        deparse(body(Pando:::.pando_infer_condition_grn_one_without_padj_cap)),
        collapse = "\n"
    )
    expect_match(body_text, "union_grn_edges", fixed = TRUE)
    expect_match(body_text, ".condition_ridge_refit_contract", fixed = TRUE)
    expect_false(grepl(".condition_fit_dictionary_prepared", body_text,
                       fixed = TRUE))
    expect_false(grepl("stats::glm", body_text, fixed = TRUE))
})

test_that("multi-task model retains the external dictionary schema", {
    expect_identical(
        Pando:::.condition_common_dictionary_schema,
        "pando_condition_grn_common_dictionary_v1"
    )
    expect_identical(
        Pando:::.condition_multitask_ridge_schema,
        "pando_condition_grn_multitask_ridge_v2"
    )
})

test_that("ridge defaults are deterministic and strictly regularized", {
    control <- Pando:::.condition_ridge_control()
    expect_true(all(control$lambda_grid > 0))
    expect_identical(control$lambda_rule, "1se")
    expect_equal(control$fusion_ratio, 1)
    expect_equal(control$cv_folds, 5L)
    expect_equal(control$seed, 1L)
    expect_identical(
        Pando:::.condition_ridge_fallback_key,
        "condition_ridge_control"
    )
})

test_that("ridge penalty is positive definite", {
    k <- 3L
    p <- 4L
    penalty <- Pando:::.condition_ridge_penalty(k, p, fusion_ratio = 1)
    eigenvalues <- eigen(penalty, symmetric = TRUE, only.values = TRUE)$values
    expect_true(all(eigenvalues > 0))
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

test_that("condition projection reports the actual fitted engine", {
    body_text <- paste(
        deparse(body(Pando:::project_condition_grn_cells)),
        collapse = "\n"
    )
    expect_match(body_text, "fit$fit_engine", fixed = TRUE)
    expect_false(grepl("full_condition_fixed_dictionary_glm", body_text,
                       fixed = TRUE))
})

test_that("quantitative penalty is BH-gated after joint ridge fitting", {
    body_text <- paste(
        deparse(body(Pando:::.condition_apply_significance_gate)),
        collapse = "\n"
    )
    expect_match(body_text, "padj < threshold", fixed = TRUE)
    expect_match(body_text, "coefficient$significant", fixed = TRUE)
    expect_match(body_text, "coefficient$penalty_effect", fixed = TRUE)
    expect_identical(
        Pando:::.condition_significant_projection_policy,
        "padj_significant_ridge_effects"
    )
})
