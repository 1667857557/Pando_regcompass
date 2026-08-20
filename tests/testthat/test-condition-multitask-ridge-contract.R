test_that("public condition API exposes design controls but not z/CV controls", {
    args <- names(formals(Pando:::infer_condition_grn.GRNData))
    expect_true(all(c(
        "cell_type_col", "condition_col", "tf_cor", "peak_cor",
        "min_cells_per_condition", "rank_action", "min_residual_df",
        "reference_condition", "parallel", "BPPARAM", "parallel_scope",
        "fallback_args"
    ) %in% args))
    expect_false(any(c(
        "ridge_control", "condition_ridge_control", "scheme_e_z", "z",
        "lambda_grid", "lambda_rule", "cv_folds", "fusion_ratio"
    ) %in% args))
    expect_equal(eval(formals(Pando:::infer_condition_grn.GRNData)$tf_cor), 0.05)
    expect_equal(eval(formals(Pando:::infer_condition_grn.GRNData)$peak_cor), 0.05)
})

test_that("conditional production schema is fixed E-star z=0.25", {
    expect_identical(
        Pando:::.condition_multitask_ridge_schema,
        "pando_condition_grn_Estar_z025_inference_separated_v1"
    )
    expect_identical(
        Pando:::.condition_fit_engine,
        "condition_union_Estar_z025_inference_separated"
    )
    expect_identical(
        Pando:::.condition_inference_schema,
        "frozen_dictionary_condition_local_gaussian_lm_edge_omnibus_v1"
    )
    expect_equal(Pando:::.condition_E_star_z, 0.25)
    expect_identical(
        Pando:::.condition_E_star_penalty_family,
        "information_scaled_sparse_deviation"
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

test_that("no-fusion condition inference reproduces Gaussian lm coefficients", {
    set.seed(27)
    x <- cbind(e1 = rnorm(90), e2 = rnorm(90))
    x[, 2] <- 0.35 * x[, 1] + rnorm(90, sd = 0.9)
    y <- 1.1 + 0.8 * x[, 1] - 0.3 * x[, 2] + rnorm(90, sd = 0.4)
    scaling <- Pando:::.condition_ridge_scaling(list(A = x, B = x), 1e-8)
    fit <- Pando:::.condition_no_fusion_condition_fit(
        x, y, scaling, rank_tol = 1e-12, min_residual_df = 1L
    )
    expected <- summary(stats::lm(y ~ x))$coefficients[-1, , drop = FALSE]
    expect_identical(fit$status, "ok")
    expect_equal(fit$estimate, expected[, 1], tolerance = 1e-8)
    expect_equal(fit$se, expected[, 2], tolerance = 1e-8)
    expect_equal(fit$statistic, expected[, 3], tolerance = 1e-8)
    expect_equal(fit$pval, expected[, 4], tolerance = 1e-8)
})

test_that("dictionary and projection policies encode common edge topology", {
    expect_identical(
        Pando:::.condition_fit_dictionary_policy,
        "global_and_condition_union_pando_correlation_supported_frozen_dictionary"
    )
    expect_identical(
        Pando:::.condition_projection_policy,
        "exact_edge_whole_network_BH_common_topology"
    )
})
