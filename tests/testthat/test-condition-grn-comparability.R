test_that("condition-wise Pando candidates form the exact common union", {
    control <- data.frame(
        target = c("G", "G"),
        tf = c("TF1", "TF3"),
        region = c("P1", "P3"),
        atac_feature_id = c("A1", "A3"),
        peak_target_cor = c(0.4, 0.6),
        tf_target_cor = c(0.5, 0.4),
        stringsAsFactors = FALSE
    )
    drug <- data.frame(
        target = "G", tf = "TF1", region = "P1", atac_feature_id = "A1",
        peak_target_cor = 0.7, tf_target_cor = 0.6,
        stringsAsFactors = FALSE
    )
    dictionary <- union_grn_edges(
        global_edges = NULL,
        condition_edges = list(Control = control, Drug = drug)
    )
    support <- Pando:::.condition_candidate_support_table(
        list(Control = control, Drug = drug)
    )
    summary <- Pando:::.condition_candidate_support_summary(dictionary, support)

    expect_s3_class(dictionary, "PandoEdgeDictionary")
    expect_setequal(
        dictionary$edge_id,
        c("G||TF1||P1", "G||TF3||P3")
    )
    expect_identical(anyDuplicated(dictionary$edge_id), 0L)
    expect_equal(
        summary$n_support_conditions[summary$edge_id == "G||TF1||P1"], 2L
    )
    expect_equal(
        summary$n_support_conditions[summary$edge_id == "G||TF3||P3"], 1L
    )
})

.rank_deficient_multitask_fixture <- function(n = 30L) {
    x1_a <- seq(-1, 1, length.out = n)
    x2_a <- cos(seq(0, pi, length.out = n))
    x1_b <- seq(-0.8, 1.2, length.out = n)
    x2_b <- sin(seq(0, pi, length.out = n))
    xa <- cbind(edge1 = x1_a, edge2 = x2_a, edge3 = x1_a + x2_a)
    xb <- cbind(edge1 = x1_b, edge2 = x2_b, edge3 = x1_b + x2_b)
    x <- list(Control = xa, Drug = xb)
    y <- list(
        Control = 1 + 2 * x1_a - 0.8 * x2_a,
        Drug = 1.2 + 2.4 * x1_b - 0.4 * x2_b
    )
    list(x = x, y = y)
}

test_that("ridge fits a raw rank-deficient common dictionary", {
    fixture <- .rank_deficient_multitask_fixture()
    scaling <- Pando:::.condition_ridge_scaling(fixture$x, 1e-8)
    fit <- Pando:::.condition_ridge_fit(
        fixture$x, fixture$y, scaling,
        lambda = 0.1, min_residual_df = 1L, inference = TRUE
    )

    expect_identical(fit$status, "ok")
    expect_true(all(fit$raw_rank < sum(fit$informative)))
    expect_true(all(is.finite(fit$beta)))
    expect_true(all(is.finite(fit$beta_z)))
    expect_true(is.finite(fit$regularized_kappa))
    expect_true(is.matrix(fit$covariance_z))
})

test_that("no-fusion ridge retains condition-specific coefficient differences", {
    fixture <- .rank_deficient_multitask_fixture()
    scaling <- Pando:::.condition_ridge_scaling(fixture$x, 1e-8)
    fit <- Pando:::.condition_ridge_fit(
        fixture$x, fixture$y, scaling,
        lambda = 0.05, min_residual_df = 1L, inference = FALSE
    )
    expect_identical(fit$status, "ok")
    expect_gt(sum((fit$beta["Control", ] - fit$beta["Drug", ])^2), 0)
})

test_that("condition-specific zero variance is diagnostic rather than fatal", {
    n <- 30L
    control_x <- cbind(
        constant_here = rep(1, n),
        shared_signal = seq(-1, 1, length.out = n)
    )
    drug_x <- cbind(
        constant_here = seq(0.5, 1.5, length.out = n),
        shared_signal = seq(-0.8, 1.2, length.out = n)
    )
    x <- list(Control = control_x, Drug = drug_x)
    y <- list(
        Control = 1 + 1.5 * control_x[, "shared_signal"],
        Drug = 1 + 0.8 * drug_x[, "constant_here"] +
            1.5 * drug_x[, "shared_signal"]
    )
    scaling <- Pando:::.condition_ridge_scaling(x, 1e-8)
    fit <- Pando:::.condition_ridge_fit(
        x, y, scaling, lambda = 0.1,
        min_residual_df = 1L, inference = TRUE
    )

    expect_identical(fit$status, "ok")
    expect_true(fit$zero_variance["Control", "constant_here"])
    expect_false(fit$zero_variance["Drug", "constant_here"])
    expect_true(fit$informative[["constant_here"]])
    expect_true(all(is.finite(fit$beta[, "constant_here"])))
})

test_that("condition-stratified CV selects one lambda for every condition", {
    fixture <- .rank_deficient_multitask_fixture(n = 36L)
    cells <- list(
        Control = paste0("c", seq_len(nrow(fixture$x$Control))),
        Drug = paste0("d", seq_len(nrow(fixture$x$Drug)))
    )
    control <- Pando:::.condition_ridge_control(list(
        lambda_grid = c(0.01, 0.1, 1),
        lambda_rule = "1se",
        cv_folds = 4L,
        seed = 7L
    ))
    folds <- Pando:::.condition_ridge_folds(
        cells, control$cv_folds, control$seed
    )
    cv <- Pando:::.condition_ridge_cv(
        fixture$x, fixture$y, folds, control,
        min_residual_df = 1L
    )

    expect_true(cv$lambda %in% control$lambda_grid)
    expect_true(cv$lambda_min %in% control$lambda_grid)
    expect_identical(names(cv$rsq_oof), c("Control", "Drug"))
    expect_true(all(is.finite(cv$rsq_oof)))
    expect_equal(nrow(cv$curve), length(control$lambda_grid))
})

test_that("joint block covariance yields finite pairwise ridge contrasts", {
    fixture <- .rank_deficient_multitask_fixture()
    scaling <- Pando:::.condition_ridge_scaling(fixture$x, 1e-8)
    fit <- Pando:::.condition_ridge_fit(
        fixture$x, fixture$y, scaling,
        lambda = 0.1, min_residual_df = 1L, inference = TRUE
    )
    edges <- data.frame(
        tf = paste0("TF", 1:3), target = rep("G", 3),
        region = paste0("P", 1:3), edge_id = colnames(fixture$x$Control),
        atac_feature_id = paste0("A", 1:3), stringsAsFactors = FALSE
    )
    informative <- sweep(!fit$zero_variance, 2L, fit$informative, "&")
    contrast <- Pando:::.condition_ridge_contrasts(
        fit, edges, informative, scaling
    )
    expect_equal(nrow(contrast), 3L)
    expect_true(all(contrast$contrast_estimable))
    expect_true(all(is.finite(contrast$contrast_estimate)))
    expect_true(all(is.finite(contrast$contrast_se)))
})

test_that("rank auditing remains available under ridge", {
    fixture <- .rank_deficient_multitask_fixture()
    scaling <- Pando:::.condition_ridge_scaling(fixture$x, 1e-8)
    fit <- Pando:::.condition_ridge_fit(
        fixture$x, fixture$y, scaling,
        lambda = 0.1, min_residual_df = 1L, inference = FALSE
    )
    expect_identical(fit$status, "ok")
    expect_true(any(fit$raw_rank < sum(fit$informative)))
})
