# Preserve the explicit infer_condition_grn() API while adding ridge controls.

#' @rdname infer_condition_grn
#' @method infer_condition_grn GRNData
#' @export
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
    ridge_control = list(),
    parallel = FALSE, BPPARAM = NULL,
    parallel_scope = c("auto", "cell_type", "target"),
    overwrite = FALSE, fallback_args = list(), verbose = TRUE, ...) {
    dots <- list(...)
    if (length(dots)) {
        label <- names(dots)
        if (is.null(label)) label <- rep("<unnamed>", length(dots))
        label[!nzchar(label)] <- "<unnamed>"
        stop("Unused condition-GRN argument(s): ",
             paste(label, collapse = ", "), call. = FALSE)
    }
    if (!is.list(fallback_args)) {
        stop("`fallback_args` must be a list.", call. = FALSE)
    }
    fallback_args[[.condition_ridge_control_key]] <-
        .condition_ridge_control(ridge_control)
    .condition_legacy_infer_condition_grn_method(
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
    )
}
