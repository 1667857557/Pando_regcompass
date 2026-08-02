test_that("fixed-dictionary networks retain accessors and exclude unavailable modules", {
    coefficients <- data.frame(
        target = c("G", "G"),
        tf = c("TF1", "TF2"),
        region = c("P1", "P2"),
        term = c("edge_0000001", "edge_0000002"),
        edge_id = c("G||TF1||P1", "G||TF2||P2"),
        atac_feature_id = c("A1", "A2"),
        estimate = c(1, NA_real_),
        std_err = c(0.1, NA_real_),
        statistic = c(10, NA_real_),
        pval = c(0.001, NA_real_),
        padj = c(0.001, NA_real_),
        significant = c(TRUE, FALSE),
        penalty_effect = c(1, 0),
        direction = c("positive", "undefined"),
        estimable = c(TRUE, FALSE),
        zero_variance = c(FALSE, FALSE),
        aliased = c(FALSE, TRUE),
        condition = "A",
        candidate_index = 1:2,
        source_global = TRUE,
        source_conditions = "A",
        n_sources = 2L,
        effect_definition = "fixed_dictionary_condition_glm_coefficient",
        inference_scope = "conditional_on_selected_edge_dictionary",
        stringsAsFactors = FALSE
    )
    fit <- data.frame(
        target = "G", condition = "A", rsq = 0.8,
        nvariables = 2L,
        nvariables_dictionary = 2L, nvariables_estimable = 1L,
        rank = 2L, residual_df = 20L, condition_number = 1,
        n_zero_variance = 0L, n_aliased = 1L,
        fit_status = "rank_deficient", intercept = 0.2,
        stringsAsFactors = FALSE
    )
    network <- methods::new(
        Class = "Network",
        features = "G",
        coefs = coefficients,
        fit = fit,
        params = list(method = "glm", fit_mode = "fixed_edge_dictionary")
    )
    expect_s4_class(network, "Network")
    tab_before <- stats::coef(network)
    expect_true(all(c("tf", "target", "region", "estimate", "padj") %in%
                    colnames(tab_before)))
    expect_true(is.na(tab_before$estimate[tab_before$tf == "TF2"]))
    expect_true(is.na(tab_before$padj[tab_before$tf == "TF2"]))
    expect_identical(Pando::gof(network)$nvariables[[1L]], 2L)

    network <- Pando::find_modules(
        network,
        p_thresh = 0.05,
        rsq_thresh = 0,
        nvar_thresh = 0,
        min_genes_per_module = 0,
        verbose = FALSE
    )
    module_meta <- Pando::NetworkModules(network)@meta
    expect_identical(unique(module_meta$tf), "TF1")
    expect_false(any(module_meta$tf == "TF2", na.rm = TRUE))
    tab_after <- stats::coef(network)
    expect_true(is.na(tab_after$estimate[tab_after$tf == "TF2"]))
    expect_true(is.na(tab_after$padj[tab_after$tf == "TF2"]))
})
