# Finalize the common-dictionary fit contract returned to downstream packages.

.pando_complete_condition_fit_contract <- function(fit) {
    if (!inherits(fit, "ConditionGRNFit")) return(fit)
    field <- "dictionary_preprocessing_provenance_verified"
    if (!field %in% names(fit)) {
        fit[[field]] <- isTRUE(attr(
            fit$edge_dictionary,
            "preprocessing_provenance_verified",
            exact = TRUE
        ))
    }
    fit
}

.pando_complete_condition_fit_contracts <- function(fits) {
    if (inherits(fits, "ConditionGRNFit")) {
        return(.pando_complete_condition_fit_contract(fits))
    }
    if (!is.list(fits)) return(fits)
    lapply(fits, .pando_complete_condition_fit_contract)
}

.pando_complete_condition_fit_object <- function(object) {
    fits <- object@grn@params$condition_grn_fits
    if (is.list(fits) && length(fits)) {
        object@grn@params$condition_grn_fits <-
            .pando_complete_condition_fit_contracts(fits)
    }
    object
}

.pando_infer_condition_grn_complete_contract_impl <-
    infer_condition_grn.GRNData

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
    answer <- do.call(
        .pando_infer_condition_grn_complete_contract_impl,
        c(list(
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
        ), list(...))
    )
    .pando_complete_condition_fit_object(answer)
}

.pando_condition_grn_fit_complete_contract_impl <-
    condition_grn_fit.GRNData

condition_grn_fit.GRNData <- function(object, cell_type = NULL, ...) {
    answer <- .pando_condition_grn_fit_complete_contract_impl(
        object = object, cell_type = cell_type, ...
    )
    .pando_complete_condition_fit_contracts(answer)
}
