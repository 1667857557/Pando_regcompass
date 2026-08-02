# Projection from a frozen dictionary must not rerun candidate discovery.

.condition_prepare_projection_input <- function(object) {
    params <- Params(object)
    gene_data <- Matrix::t(LayerData(
        object, assay = params$rna_assay, layer = "data"
    ))
    peak_data_all <- Matrix::t(LayerData(
        object, assay = params$peak_assay, layer = "data"
    ))
    common_cells <- intersect(rownames(gene_data), rownames(peak_data_all))
    if (!length(common_cells)) {
        stop("RNA and ATAC assays do not share paired cells.", call. = FALSE)
    }
    gene_data <- gene_data[common_cells, , drop = FALSE]
    peak_data_all <- peak_data_all[common_cells, , drop = FALSE]
    regions <- NetworkRegions(object)
    if (!length(regions@peaks) || any(
        regions@peaks < 1L | regions@peaks > ncol(peak_data_all)
    )) {
        stop("Pando region-to-peak indices are invalid.", call. = FALSE)
    }
    peak_data <- peak_data_all[, regions@peaks, drop = FALSE]
    region_id <- rownames(regions@motifs@data)
    if (is.null(region_id) || length(region_id) != ncol(peak_data) ||
        anyNA(region_id) || any(!nzchar(region_id)) ||
        anyDuplicated(region_id)) {
        stop("Pando motif regions are not aligned to measured peaks.",
             call. = FALSE)
    }
    atac_feature_id <- colnames(peak_data)
    colnames(peak_data) <- region_id
    list(
        gene_data = gene_data,
        peak_data = peak_data,
        region_map = data.frame(
            region = region_id,
            atac_feature_id = atac_feature_id,
            stringsAsFactors = FALSE
        )
    )
}

#' @rdname project_condition_grn_cells
#' @export
project_condition_grn_cells <- function(
    object, fit, targets = NULL, significant_only = TRUE,
    return_edge_contributions = FALSE) {
    if (!inherits(object, "GRNData") || !inherits(fit, "ConditionGRNFit") ||
        !identical(fit$schema_version, .condition_common_dictionary_schema)) {
        stop("Common-dictionary object and fit are required.", call. = FALSE)
    }
    prepared <- .condition_prepare_projection_input(object)
    dictionary <- as.data.frame(fit$edge_dictionary, stringsAsFactors = FALSE)
    required <- c("edge_id", "target", "tf", "region", "atac_feature_id")
    if (!all(required %in% colnames(dictionary)) ||
        anyDuplicated(dictionary$edge_id) ||
        any(!dictionary$target %in% colnames(prepared$gene_data)) ||
        any(!dictionary$tf %in% colnames(prepared$gene_data)) ||
        any(!dictionary$region %in% colnames(prepared$peak_data))) {
        stop("The fitted edge dictionary is incompatible with the input object.",
             call. = FALSE)
    }
    mapped <- prepared$region_map$atac_feature_id[
        match(dictionary$region, prepared$region_map$region)
    ]
    if (anyNA(mapped) || !identical(
        as.character(mapped), as.character(dictionary$atac_feature_id)
    )) {
        stop("The fitted region-to-ATAC mapping changed after fitting.",
             call. = FALSE)
    }
    coefficient <- as.data.frame(fit$coefficients, stringsAsFactors = FALSE)
    if (!is.null(targets)) {
        targets <- intersect(as.character(targets), fit$target_genes)
        coefficient <- coefficient[
            coefficient$target %in% targets, , drop = FALSE
        ]
    }
    cells <- unique(unlist(fit$condition_cell_ids, use.names = FALSE))
    cells <- cells[cells %in% rownames(prepared$gene_data)]
    cell_condition <- rep(NA_character_, length(cells))
    names(cell_condition) <- cells
    for (condition in fit$condition_levels) {
        selected <- intersect(cells, fit$condition_cell_ids[[condition]])
        cell_condition[selected] <- condition
    }
    if (anyNA(cell_condition)) {
        stop("Fitted condition cells do not align to paired assay matrices.",
             call. = FALSE)
    }
    edge_contribution <- matrix(
        0, nrow = length(cells), ncol = nrow(coefficient),
        dimnames = list(
            cells,
            paste(coefficient$edge_id, coefficient$condition, sep = "@@")
        )
    )
    effect <- if (isTRUE(significant_only)) {
        coefficient$penalty_effect
    } else {
        coefficient$estimate
    }
    effect[!is.finite(effect)] <- 0
    for (j in seq_len(nrow(coefficient))) {
        selected <- cells[cell_condition == coefficient$condition[[j]]]
        if (!length(selected) || effect[[j]] == 0) next
        edge_contribution[selected, j] <-
            as.numeric(prepared$gene_data[selected, coefficient$tf[[j]]]) *
            as.numeric(prepared$peak_data[selected, coefficient$region[[j]]]) *
            effect[[j]]
    }
    target_names <- unique(as.character(coefficient$target))
    gene_score <- matrix(
        0, nrow = length(cells), ncol = length(target_names),
        dimnames = list(cells, target_names)
    )
    for (target in target_names) {
        index <- which(coefficient$target == target)
        gene_score[, target] <- rowSums(
            edge_contribution[, index, drop = FALSE]
        )
    }
    answer <- list(
        schema_version = "pando_condition_projection_common_dictionary_v1",
        gene_score = gene_score,
        edge_contribution = if (isTRUE(return_edge_contributions)) {
            edge_contribution
        } else NULL,
        cell_condition = cell_condition,
        cell_type = fit$cell_type,
        effect_column = if (isTRUE(significant_only)) {
            "penalty_effect"
        } else "estimate",
        projection_origin = "full_condition_fixed_dictionary_glm",
        coefficient_scale = fit$coefficient_scale,
        projection_policy = if (isTRUE(significant_only)) {
            "BH_adjusted_p_below_threshold_penalty_effect"
        } else "all_estimable_condition_coefficients"
    )
    class(answer) <- c("PandoConditionProjection", "list")
    answer
}
