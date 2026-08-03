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
