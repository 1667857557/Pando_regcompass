# Projection from a frozen dictionary must not rerun candidate discovery.

.condition_prepare_projection_input <- function(object, fit) {
    params <- Params(object)
    gene_data <- Matrix::t(LayerData(
        object, assay = params$rna_assay, layer = fit$rna_layer
    ))
    peak_data_all <- Matrix::t(LayerData(
        object, assay = params$peak_assay, layer = fit$peak_layer
    ))
    common_cells <- intersect(rownames(gene_data), rownames(peak_data_all))
    if (!length(common_cells)) {
        stop("RNA and ATAC assays do not share paired cells.", call. = FALSE)
    }
    gene_data <- gene_data[common_cells, , drop = FALSE]
    peak_data_all <- peak_data_all[common_cells, , drop = FALSE]
    observed_fingerprint <- .condition_preprocessing_fingerprint(
        object = object, gene_data = gene_data, peak_data = peak_data_all,
        rna_layer = fit$rna_layer, peak_layer = fit$peak_layer,
        peak_value_type = fit$peak_value_type
    )
    if (!identical(observed_fingerprint, fit$preprocessing_fingerprint)) {
        stop("RNA/ATAC preprocessing identity changed after condition fitting.",
             call. = FALSE)
    }
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

.condition_resolve_projection_targets <- function(targets, fitted_targets) {
    fitted_targets <- unique(as.character(fitted_targets))
    if (!length(fitted_targets) || anyNA(fitted_targets) ||
        any(!nzchar(fitted_targets))) {
        stop("Fitted target names must be complete non-empty strings.",
             call. = FALSE)
    }
    fitted_key <- toupper(fitted_targets)
    if (anyDuplicated(fitted_key)) {
        duplicated_key <- unique(fitted_key[duplicated(fitted_key)])
        stop(
            "Fitted targets are ambiguous after case normalization: ",
            paste(utils::head(duplicated_key, 10L), collapse = ", "),
            call. = FALSE
        )
    }
    if (is.null(targets)) return(fitted_targets)
    requested <- unique(as.character(targets))
    requested <- requested[!is.na(requested) & nzchar(requested)]
    index <- match(toupper(requested), fitted_key)
    unique(fitted_targets[index[!is.na(index)]])
}

.condition_validate_projection_cells <- function(fit, available_cells) {
    levels <- as.character(fit$condition_levels)
    cell_lists <- fit$condition_cell_ids
    list_names <- names(cell_lists)
    if (!length(levels) || anyNA(levels) || any(!nzchar(levels)) ||
        anyDuplicated(levels) || !is.list(cell_lists) ||
        is.null(list_names) || anyNA(list_names) ||
        any(!nzchar(list_names)) || anyDuplicated(list_names) ||
        !all(levels %in% list_names)) {
        stop("Fitted condition cell lists are incomplete or ambiguously named.",
             call. = FALSE)
    }
    cells_by_condition <- cell_lists[levels]
    condition_sizes <- lengths(cells_by_condition)
    if (any(condition_sizes < 1L)) {
        stop("Every fitted condition must contain at least one paired cell.",
             call. = FALSE)
    }
    all_cells <- as.character(unlist(cells_by_condition, use.names = FALSE))
    if (!length(all_cells) || anyNA(all_cells) || any(!nzchar(all_cells))) {
        stop("Fitted condition cell IDs must be complete non-empty strings.",
             call. = FALSE)
    }
    if (anyDuplicated(all_cells)) {
        stop("A fitted cell is assigned to more than one condition.",
             call. = FALSE)
    }
    available_cells <- as.character(available_cells)
    if (anyNA(available_cells) || any(!nzchar(available_cells)) ||
        anyDuplicated(available_cells)) {
        stop("Projection assay cell IDs must be complete and unique.",
             call. = FALSE)
    }
    missing <- setdiff(all_cells, available_cells)
    if (length(missing)) {
        stop(
            "The projection object is missing ", length(missing),
            " fitted paired cell(s); first missing ID: ", missing[[1L]],
            call. = FALSE
        )
    }
    list(
        cells = all_cells,
        cells_by_condition = cells_by_condition,
        condition_levels = levels
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
    required_provenance <- c(
        "rna_layer", "peak_layer", "peak_value_type",
        "preprocessing_fingerprint"
    )
    if (!all(required_provenance %in% names(fit)) ||
        any(!nzchar(vapply(fit[required_provenance], as.character,
                           character(1))))) {
        stop("The fitted condition GRN lacks preprocessing provenance.",
             call. = FALSE)
    }
    prepared <- .condition_prepare_projection_input(object, fit)
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
    coefficient_required <- c(
        "edge_id", "target", "tf", "region", "condition",
        "estimate", "penalty_effect", "active_in_regcompass"
    )
    if (!all(coefficient_required %in% colnames(coefficient))) {
        stop("The condition fit lacks E-star/RegCompass projection fields.",
             call. = FALSE)
    }
    resolved_targets <- .condition_resolve_projection_targets(
        targets, fit$target_genes
    )
    if (!is.null(targets)) {
        coefficient <- coefficient[
            coefficient$target %in% resolved_targets, , drop = FALSE
        ]
    }
    cell_contract <- .condition_validate_projection_cells(
        fit, rownames(prepared$gene_data)
    )
    cells <- cell_contract$cells
    cell_condition <- rep(NA_character_, length(cells))
    names(cell_condition) <- cells
    for (condition in cell_contract$condition_levels) {
        selected <- cell_contract$cells_by_condition[[condition]]
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
    # `penalty_effect` is always the continuous production E-star coefficient.
    # With `significant_only = TRUE`, the same exact-edge whole-network BH
    # topology is applied to every condition before each condition's own
    # continuous production coefficient is used.
    effect <- if (isTRUE(significant_only)) {
        value <- as.numeric(coefficient$penalty_effect)
        value[!(coefficient$active_in_regcompass %in% TRUE)] <- 0
        value
    } else {
        as.numeric(coefficient$estimate)
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
    default_policy <- if (
        is.character(fit$projection_policy) &&
        length(fit$projection_policy) == 1L &&
        !is.na(fit$projection_policy) && nzchar(fit$projection_policy)
    ) fit$projection_policy else .condition_projection_policy
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
        projection_origin = if (
            is.character(fit$fit_engine) && length(fit$fit_engine) == 1L &&
            !is.na(fit$fit_engine) && nzchar(fit$fit_engine)
        ) fit$fit_engine else "condition_fit_engine_unspecified",
        coefficient_scale = fit$coefficient_scale,
        projection_policy = if (isTRUE(significant_only)) {
            default_policy
        } else "all_finite_condition_coefficients"
    )
    class(answer) <- c("PandoConditionProjection", "list")
    answer
}
