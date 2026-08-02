test_that("exact edge union does not create Cartesian edges", {
    global <- data.frame(
        target = "G", tf = "TF1", region = "P1",
        atac_feature_id = "A1", peak_target_cor = 0.2,
        tf_target_cor = 0.3, stringsAsFactors = FALSE
    )
    condition <- list(
        A = data.frame(
            target = "G", tf = "TF2", region = "P2",
            atac_feature_id = "A2", peak_target_cor = 0.4,
            tf_target_cor = 0.5, stringsAsFactors = FALSE
        ),
        B = data.frame(
            target = "G", tf = "TF1", region = "P1",
            atac_feature_id = "A1", peak_target_cor = 0.1,
            tf_target_cor = 0.2, stringsAsFactors = FALSE
        )
    )
    dictionary <- union_grn_edges(global, condition)
    expect_identical(nrow(dictionary), 2L)
    expect_setequal(
        dictionary$edge_id,
        c("G||TF1||P1", "G||TF2||P2")
    )
    expect_false("G||TF1||P2" %in% dictionary$edge_id)
    expect_false("G||TF2||P1" %in% dictionary$edge_id)
})

test_that("fixed dictionary GLM permits opposite condition directions", {
    set.seed(11)
    n <- 300L
    x1 <- rnorm(n)
    x2 <- rnorm(n)
    X <- cbind(x1, x2)
    fit_a <- Pando:::.condition_fit_target_matrix(
        response = 1.5 * x1 - 0.7 * x2 + rnorm(n, sd = 0.05),
        predictor = X,
        terms = c("edge_1", "edge_2"),
        rank_action = "error"
    )
    fit_b <- Pando:::.condition_fit_target_matrix(
        response = -1.2 * x1 + 0.9 * x2 + rnorm(n, sd = 0.05),
        predictor = X,
        terms = c("edge_1", "edge_2"),
        rank_action = "error"
    )
    expect_gt(fit_a$coefs$estimate[[1L]], 0)
    expect_lt(fit_b$coefs$estimate[[1L]], 0)
    expect_lt(fit_a$coefs$estimate[[2L]], 0)
    expect_gt(fit_b$coefs$estimate[[2L]], 0)
})

test_that("zero variance and aliased edges are not fitted zeros", {
    set.seed(12)
    x <- rnorm(100)
    X <- cbind(variable = x, duplicate = x, closed = 0)
    fit <- Pando:::.condition_fit_target_matrix(
        response = 2 * x + rnorm(100, sd = 0.1),
        predictor = X,
        terms = c("edge_1", "edge_2", "edge_3"),
        rank_action = "mark"
    )
    expect_true(fit$coefs$zero_variance[[3L]])
    expect_false(fit$coefs$estimable[[3L]])
    expect_true(is.na(fit$coefs$estimate[[3L]]))
    expect_true(any(fit$coefs$aliased[1:2]))
    expect_true(any(is.na(fit$coefs$estimate[1:2])))
})

test_that("condition API contains only common dictionary controls", {
    formal_names <- names(formals(Pando:::infer_condition_grn.GRNData))
    expect_true(all(c(
        "tf_cor", "peak_cor", "adjust_method", "padj_threshold",
        "rank_action", "min_residual_df"
    ) %in% formal_names))
    expect_false(any(c(
        "candidate_screen", "condition_mix", "condition_weight",
        "nlambda", "lambda", "outer_nfolds", "inner_nfolds",
        "lambda_selection", "engine_control", "scale"
    ) %in% formal_names))
    description <- utils::packageDescription("Pando")
    expect_identical(
        description[["Config/Pando/ConditionGRNSchema"]],
        "pando_condition_grn_common_dictionary_v1"
    )
})
