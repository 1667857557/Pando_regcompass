test_that("fixed-dictionary networks retain standard accessors and modules", {
    dictionary <- data.frame(
        edge_id = "G||TF1||P1", target = "G", tf = "TF1",
        region = "P1", atac_feature_id = "A1", candidate_index = 1L,
        stringsAsFactors = FALSE
    )
    coefficients <- data.frame(
        target = "G", tf = "TF1", region = "P1",
        edge_id = "G||TF1||P1", atac_feature_id = "A1",
        candidate_index = 1L, estimate = 1, std_err = 0.1,
        statistic = 10, pval = 0.001, padj = 0.001,
        direction = "positive", estimable = TRUE,
        zero_variance = FALSE, aliased = FALSE,
        stringsAsFactors = FALSE
    )
    fit <- data.frame(
        target = "G", condition = "A", rsq = 0.8,
        nvariables_dictionary = 1L, nvariables_estimable = 1L,
        rank = 2L, residual_df = 20L, condition_number = 1,
        n_zero_variance = 0L, n_aliased = 0L, fit_status = "ok",
        intercept = 0.2, stringsAsFactors = FALSE
    )
    network <- Pando:::.condition_make_network(
        coefficients = coefficients,
        fit = fit,
        dictionary = dictionary,
        condition_label = "A",
        params = list(method = "glm")
    )
    expect_s4_class(network, "Network")
    tab <- stats::coef(network)
    expect_true(all(c("tf", "target", "region", "estimate", "padj") %in%
                    colnames(tab)))
    expect_identical(tab$padj[tab$term == "(Intercept)"], 1)
    expect_identical(Pando::gof(network)$nvariables[[1L]], 1L)

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
    expect_false(any(module_meta$tf == "(Intercept)", na.rm = TRUE))
})
