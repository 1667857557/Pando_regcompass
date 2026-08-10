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
        deparse(body(Pando:::.pando_infer_condition_grn_one)),
        collapse = "\n"
    )
    expect_match(body_text, "union_grn_edges", fixed = TRUE)
    expect_match(body_text, ".condition_ridge_refit_contract", fixed = TRUE)
    expect_false(grepl(".condition_fit_dictionary_prepared", body_text,
                       fixed = TRUE))
    expect_false(grepl("stats::glm", body_text, fixed = TRUE))
})

test_that("multi-task model retains the external common-dictionary schema", {
    expect_identical(
        Pando:::.condition_common_dictionary_schema,
        "pando_condition_grn_common_dictionary_v1"
    )
    expect_identical(
        Pando:::.condition_multitask_ridge_schema,
        "pando_condition_grn_multitask_ridge_v1"
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

test_that("runtime integration does not redefine the public S3 method", {
    runtime <- readLines(
        system.file("R", "condition_multitask_ridge_runtime.R", package = "Pando"),
        warn = FALSE
    )
    if (length(runtime)) {
        expect_false(any(grepl(
            "^infer_condition_grn\\.GRNData\\s*<-\\s*function",
            runtime
        )))
    } else {
        succeed()
    }
})
