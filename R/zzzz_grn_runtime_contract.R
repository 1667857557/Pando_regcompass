# Final runtime contract for the route-aware GRN wrappers.

.pando_infer_grn_routed_impl <- infer_grn.GRNData

infer_grn.GRNData <- function(
    object,
    genes = NULL,
    network_name = paste0(method, '_network'),
    peak_to_gene_method = c('Signac', 'GREAT'),
    upstream = 100000,
    downstream = 0,
    extend = 1000000,
    only_tss = FALSE,
    parallel = FALSE,
    tf_cor = 0.1,
    peak_cor = 0.,
    aggregate_rna_col = NULL,
    aggregate_peaks_col = NULL,
    method = c('glm', 'glmnet', 'cv.glmnet', 'brms', 'xgb',
               'bagging_ridge', 'bayesian_ridge'),
    alpha = 0.5,
    family = 'gaussian',
    interaction_term = ':',
    adjust_method = 'fdr',
    scale = FALSE,
    verbose = TRUE,
    ...
) {
    method <- match.arg(method)
    peak_to_gene_method <- match.arg(peak_to_gene_method)
    if (missing(network_name)) network_name <- paste0(method, "_network")
    do.call(.pando_infer_grn_routed_impl, c(list(
        object = object,
        genes = genes,
        network_name = network_name,
        peak_to_gene_method = peak_to_gene_method,
        upstream = upstream,
        downstream = downstream,
        extend = extend,
        only_tss = only_tss,
        parallel = parallel,
        tf_cor = tf_cor,
        peak_cor = peak_cor,
        aggregate_rna_col = aggregate_rna_col,
        aggregate_peaks_col = aggregate_peaks_col,
        method = method,
        alpha = alpha,
        family = family,
        interaction_term = interaction_term,
        adjust_method = adjust_method,
        scale = scale,
        verbose = verbose
    ), list(...)))
}

.pando_infer_condition_grn_routed_impl <- infer_condition_grn.GRNData

infer_condition_grn.GRNData <- function(
    object, cell_type_col = NULL, condition_col = NULL, cell_type = NULL,
    genes = NULL, network_name = "condition_grn",
    peak_to_gene_method = c("Signac", "GREAT"), upstream = 100000,
    downstream = 0, extend = 1000000, only_tss = FALSE,
    peak_to_gene_domains = NULL, rna_layer = "data", peak_layer = "data",
    peak_value_type = c("normalized", "probability", "other"),
    tf_cor = 0.1, peak_cor = 0,
    min_cells_per_condition = 50L,
    small_condition_action = c("error", "drop_condition", "skip_cell_type"),
    adjust_method = "BH", padj_threshold = 0.05,
    rank_action = c("mark", "error"), min_residual_df = 1L,
    parallel = FALSE, BPPARAM = NULL,
    parallel_scope = c("auto", "cell_type", "target"),
    overwrite = FALSE, fallback_args = list(), verbose = TRUE, ...) {
    .pando_validate_bpparam(BPPARAM)
    do.call(.pando_infer_condition_grn_routed_impl, c(list(
        object = object,
        cell_type_col = cell_type_col,
        condition_col = condition_col,
        cell_type = cell_type,
        genes = genes,
        network_name = network_name,
        peak_to_gene_method = peak_to_gene_method,
        upstream = upstream,
        downstream = downstream,
        extend = extend,
        only_tss = only_tss,
        peak_to_gene_domains = peak_to_gene_domains,
        rna_layer = rna_layer,
        peak_layer = peak_layer,
        peak_value_type = peak_value_type,
        tf_cor = tf_cor,
        peak_cor = peak_cor,
        min_cells_per_condition = min_cells_per_condition,
        small_condition_action = small_condition_action,
        adjust_method = adjust_method,
        padj_threshold = padj_threshold,
        rank_action = rank_action,
        min_residual_df = min_residual_df,
        parallel = parallel,
        BPPARAM = BPPARAM,
        parallel_scope = parallel_scope,
        overwrite = overwrite,
        fallback_args = fallback_args,
        verbose = verbose
    ), list(...)))
}

.pando_merge_cell_type_grn_results_impl <-
    .pando_merge_cell_type_grn_results

.pando_merge_cell_type_grn_results <- function(
    object, results, condition_col, cell_type_col) {
    answer <- .pando_merge_cell_type_grn_results_impl(
        object = object,
        results = results,
        condition_col = condition_col,
        cell_type_col = cell_type_col
    )
    condition_levels <- unique(unlist(lapply(results, function(value) {
        as.character(value@grn@params$condition_levels)
    }), use.names = FALSE))
    condition_levels <- condition_levels[
        !is.na(condition_levels) & nzchar(condition_levels)
    ]
    answer@grn@params$condition_levels <- condition_levels

    fallback_reasons <- unique(unlist(lapply(results, function(value) {
        as.character(value@grn@params$standard_fallback_reason)
    }), use.names = FALSE))
    fallback_reasons <- fallback_reasons[
        !is.na(fallback_reasons) & nzchar(fallback_reasons)
    ]
    answer@grn@params$standard_fallback_reason <- if (length(fallback_reasons)) {
        paste(fallback_reasons, collapse = ";")
    } else NULL

    added_networks <- setdiff(
        names(answer@grn@networks), names(object@grn@networks)
    )
    if (length(added_networks)) {
        answer@grn@active_network <- tail(added_networks, 1L)
    }
    answer
}
