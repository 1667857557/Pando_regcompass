test_that("public condition API remains canonical and explicit", {
    args <- names(formals(Pando:::infer_condition_grn.GRNData))
    expect_true(all(c(
        "cell_type_col", "condition_col", "tf_cor", "peak_cor",
        "min_cells_per_condition", "rank_action", "min_residual_df",
        "parallel", "BPPARAM", "parallel_scope", "fallback_args"
    ) %in% args))
    expect_false("ridge_control" %in% args)
})

test_that("multi-condition path retains exact condition union and one ridge fit", {
    body_text <- paste(
        deparse(body(Pando:::.pando_infer_condition_grn_one)),
        collapse = "\n"
    )
    expect_match(body_text, "union_grn_edges", fixed = TRUE)
    expect_match(body_text, ".condition_ridge_refit_contract", fixed = TRUE)
    expect_false(grepl("stats::glm", body_text, fixed = TRUE))

    refit_text <- paste(
        deparse(body(Pando:::.condition_ridge_refit_contract)),
        collapse = "\n"
    )
    expect_match(refit_text, ".condition_ridge_refit_contract_one_pass", fixed = TRUE)
    expect_false(grepl(".condition_dictionary_screen", refit_text, fixed = TRUE))
    expect_false(grepl(".condition_subset_dictionary", refit_text, fixed = TRUE))
    expect_false(grepl("ridge_preliminary", refit_text, fixed = TRUE))
    expect_false(grepl("ridge_final", refit_text, fixed = TRUE))
})

test_that("condition model advertises the revised schemas", {
    expect_identical(
        Pando:::.condition_common_dictionary_schema,
        "pando_condition_grn_common_dictionary_v2"
    )
    expect_identical(
        Pando:::.condition_multitask_ridge_schema,
        "pando_condition_grn_multitask_ridge_v3"
    )
    expect_identical(
        Pando:::.condition_fit_engine,
        "condition_union_single_no_fusion_common_lambda_ridge"
    )
})

test_that("ridge defaults are deterministic and no-fusion", {
    control <- Pando:::.condition_ridge_control()
    expect_true(all(control$lambda_grid > 0))
    expect_identical(control$lambda_rule, "1se")
    expect_equal(control$fusion_ratio, 0)
    expect_equal(control$cv_folds, 5L)
    expect_equal(control$seed, 1L)
    expect_error(
        Pando:::.condition_ridge_control(list(fusion_ratio = 1)),
        "fusion_ratio = 0",
        fixed = TRUE
    )
})

test_that("condition ridge penalty is ordinary ridge only", {
    k <- 3L
    p <- 4L
    penalty <- Pando:::.condition_ridge_penalty(k, p, fusion_ratio = 0)
    expect_equal(penalty, diag(k * p))
    expect_error(
        Pando:::.condition_ridge_penalty(k, p, fusion_ratio = 1),
        "does not permit fusion",
        fixed = TRUE
    )
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

test_that("active projection requires both BH and condition Pando support", {
    fit <- list(
        padj_threshold = 0.05,
        edge_dictionary = structure(
            data.frame(
                edge_id = c("G||TF1||P1", "G||TF2||P2"),
                stringsAsFactors = FALSE
            ),
            condition_support_table = data.frame(
                edge_id = c("G||TF1||P1", "G||TF2||P2", "G||TF1||P1"),
                condition = c("A", "A", "B"),
                peak_target_cor = c(0.2, 0.3, 0.25),
                tf_target_cor = c(0.4, 0.5, 0.35),
                stringsAsFactors = FALSE
            )
        ),
        coefficients = data.frame(
            edge_id = rep(c("G||TF1||P1", "G||TF2||P2"), 2L),
            condition = rep(c("A", "B"), each = 2L),
            estimate = c(0.5, 0.4, 0.6, 0.7),
            estimable = TRUE,
            padj = c(0.01, 0.02, 0.01, 0.01),
            stringsAsFactors = FALSE
        )
    )
    class(fit) <- c("ConditionGRNFit", "list")
    gated <- Pando:::.condition_annotate_local_pando_support(fit)
    expect_identical(gated$coefficients$statistically_supported,
                     rep(TRUE, 4L))
    expect_identical(gated$coefficients$local_support,
                     c(TRUE, TRUE, TRUE, FALSE))
    expect_identical(gated$coefficients$active,
                     c(TRUE, TRUE, TRUE, FALSE))
    expect_equal(gated$coefficients$penalty_effect,
                 c(0.5, 0.4, 0.6, 0))
})

test_that("fit dictionary policy is condition-union and frozen", {
    expect_identical(
        Pando:::.condition_fit_dictionary_policy,
        "condition_union_pando_correlation_supported_frozen_dictionary"
    )
    expect_identical(
        Pando:::.condition_significant_projection_policy,
        "active_condition_pando_support_and_bh_ridge_effects"
    )
})
