# Normalize condition-GRN output schema independently of execution mode.
#
# The common-dictionary implementation stores preprocessing provenance on the
# frozen edge dictionary. Every ConditionGRNFit must expose the corresponding
# verification flag directly because downstream consumers validate the fit
# contract without depending on data-frame attributes.

.pando_condition_runtime_impl <- infer_condition_grn.GRNData

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
    parallel_scope = c("auto", "condition_cell_type", "cell_type", "target"),
    overwrite = FALSE, fallback_args = list(), verbose = TRUE, ...) {
    answer <- .pando_condition_runtime_impl(
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
        verbose = verbose,
        ...
    )
    fits <- answer@grn@params$condition_grn_fits
    if (is.list(fits) && length(fits)) {
        fits <- lapply(fits, function(fit) {
            dictionary <- fit$edge_dictionary
            fit$dictionary_preprocessing_provenance_verified <- isTRUE(attr(
                dictionary, "preprocessing_provenance_verified", exact = TRUE
            ))
            fit
        })
        answer@grn@params$condition_grn_fits <- fits
    }
    answer
}
