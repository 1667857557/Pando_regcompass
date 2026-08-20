test_that("conditional edge inference requires BH and a valid threshold", {
    expect_identical(Pando:::.condition_validate_adjust_method("BH"), "BH")
    expect_identical(Pando:::.condition_validate_adjust_method("bh"), "BH")
    expect_error(Pando:::.condition_validate_adjust_method("holm"),
                 'adjust_method = "BH"', fixed = TRUE)
    expect_equal(Pando:::.condition_validate_padj_threshold(0.05), 0.05)
    expect_error(Pando:::.condition_validate_padj_threshold(1))
    expect_error(Pando:::.condition_validate_padj_threshold(0))
})

.edge_gate_fixture <- function() {
    edge <- data.frame(
        edge_id = c("G||TF1||P1", "G||TF2||P2", "G||TF3||P3"),
        target = "G", tf = c("TF1", "TF2", "TF3"),
        region = c("P1", "P2", "P3"),
        atac_feature_id = c("P1", "P2", "P3"),
        candidate_index = 1:3,
        stringsAsFactors = FALSE
    )
    coefficient <- data.frame(
        edge_id = rep(edge$edge_id, each = 2L),
        target = "G",
        condition = rep(c("A", "B"), 3L),
        estimate = c(0.7, 0.65, 0.6, 0.55, 0.10, 0.05),
        penalty_effect = c(0.7, 0.65, 0.6, 0.55, 0.10, 0.05),
        inference_estimate = c(0.7, NA, 0.6, 0.55, 0.10, 0.05),
        inference_se = c(0.10, NA, 0.10, 0.10, 0.20, 0.20),
        inference_variance = c(0.01, NA, 0.01, 0.01, 0.04, 0.04),
        inference_statistic = c(7, NA, 6, 5.5, 0.5, 0.25),
        condition_pval = c(0.001, NA, 1e-5, 2e-5, 0.62, 0.80),
        condition_inference_estimable =
            c(TRUE, FALSE, TRUE, TRUE, TRUE, TRUE),
        inference_estimable = c(TRUE, FALSE, TRUE, TRUE, TRUE, TRUE),
        stringsAsFactors = FALSE
    )
    fit <- list(
        edge_dictionary = edge,
        condition_levels = c("A", "B"),
        padj_threshold = 0.05,
        adjust_method = "BH",
        dictionary_support_table = data.frame(
            edge_id = c(edge$edge_id, edge$edge_id[c(1, 3)]),
            source_type = c(rep("global", 3), rep("condition", 2)),
            condition = c(NA, NA, NA, "A", "A"),
            peak_target_cor = c(0.2, 0.3, 0.4, 0.25, 0.45),
            tf_target_cor = c(0.3, 0.4, 0.5, 0.35, 0.55),
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

test_that("one-condition estimability degenerates to the exact t-test P value", {
    fit <- .edge_gate_fixture()
    edge <- Pando:::.condition_exact_edge_inference(fit, fit$coefficients)
    one <- edge[edge$edge_id == "G||TF1||P1", , drop = FALSE]
    expect_identical(one$edge_df, 1L)
    expect_identical(one$edge_inference_test, "single_condition_exact_t")
    expect_equal(one$edge_pval, 0.001, tolerance = 1e-15)
})

test_that("multi-condition edge uses an omnibus Wald chi-square", {
    fit <- .edge_gate_fixture()
    edge <- Pando:::.condition_exact_edge_inference(fit, fit$coefficients)
    two <- edge[edge$edge_id == "G||TF2||P2", , drop = FALSE]
    expected_statistic <- 0.6^2 / 0.01 + 0.55^2 / 0.01
    expect_identical(two$edge_df, 2L)
    expect_identical(two$edge_inference_test,
                     "independent_condition_wald_chisq")
    expect_equal(two$edge_statistic, expected_statistic, tolerance = 1e-12)
    expect_equal(two$edge_pval,
                 pchisq(expected_statistic, df = 2, lower.tail = FALSE),
                 tolerance = 1e-15)
})

test_that("BH is applied once across exact edges and topology is common", {
    gated <- Pando:::.condition_apply_activity_gate(.edge_gate_fixture())
    edge <- gated$edge_inference
    expect_identical(unique(edge$bh_scope),
                     "exact_edge_whole_cell_type_network_BH")
    expect_identical(unique(edge$bh_family_size), 3L)
    expect_true(edge$edge_supported[edge$edge_id == "G||TF1||P1"])
    expect_true(edge$edge_supported[edge$edge_id == "G||TF2||P2"])
    expect_false(edge$edge_supported[edge$edge_id == "G||TF3||P3"])
    for (id in edge$edge_id) {
        rows <- gated$coefficients$edge_id == id
        expect_length(unique(gated$coefficients$active_in_regcompass[rows]), 1L)
        expect_identical(
            unique(gated$coefficients$active_in_regcompass[rows]),
            edge$edge_supported[edge$edge_id == id]
        )
    }
})

test_that("candidate support is provenance and does not recreate local topology", {
    gated <- Pando:::.condition_apply_activity_gate(.edge_gate_fixture())
    row <- gated$coefficients$edge_id == "G||TF1||P1" &
        gated$coefficients$condition == "B"
    expect_false(gated$coefficients$local_support[row])
    expect_true(gated$coefficients$global_support[row])
    expect_true(gated$coefficients$active_in_regcompass[row])
    expect_equal(gated$coefficients$penalty_effect[row], 0.65)
})
