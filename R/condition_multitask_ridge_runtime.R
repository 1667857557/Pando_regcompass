# Multi-task ridge integration behind the canonical infer_condition_grn.GRNData
# definition in condition_grn.R. This file loads after condition_grn.R and
# changes only the internal multi-condition estimator.

.condition_ridge_fallback_key <- "condition_ridge_control"

.pando_infer_condition_grn_one <- function(
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
    parallel = FALSE, overwrite = FALSE, fallback_args = list(),
    verbose = TRUE, ...) {
    dots <- list(...)
    if (length(dots)) {
        label <- names(dots)
        if (is.null(label)) label <- rep("<unnamed>", length(dots))
        label[!nzchar(label)] <- "<unnamed>"
        stop("Unused condition-GRN argument(s): ",
             paste(label, collapse = ", "), call. = FALSE)
    }
    if (!inherits(object, "GRNData")) {
        stop("`object` must be a GRNData object.", call. = FALSE)
    }
    if (!is.list(fallback_args)) {
        stop("`fallback_args` must be a list.", call. = FALSE)
    }
    control <- fallback_args[[.condition_ridge_fallback_key]]
    fallback_args[[.condition_ridge_fallback_key]] <- NULL
    control <- .condition_ridge_control(control)
    peak_to_gene_method <- match.arg(peak_to_gene_method)
    small_condition_action <- match.arg(small_condition_action)
    rank_action <- match.arg(rank_action)
    peak_value_type <- match.arg(peak_value_type)

    metadata <- object@data@meta.data
    condition_levels <- .condition_resolve_levels(metadata, condition_col)
    if (length(condition_levels) < 2L) {
        standard <- .condition_standard_by_cell_type(
            object = object, metadata = metadata,
            cell_type_col = cell_type_col, cell_type = cell_type,
            genes = genes, network_name = network_name,
            peak_to_gene_method = peak_to_gene_method,
            upstream = upstream, downstream = downstream, extend = extend,
            only_tss = only_tss, parallel = parallel,
            tf_cor = tf_cor, peak_cor = peak_cor,
            fallback_args = fallback_args, verbose = verbose,
            overwrite = overwrite
        )
        object <- standard$object
        object@grn@params$analysis_mode <- "standard_grn"
        object@grn@params$condition_col <- condition_col
        object@grn@params$condition_levels <- condition_levels
        object@grn@params$condition_coefficients_calculated <- FALSE
        object@grn@params$standard_network_index <- standard$network_index
        object@grn@params$standard_fallback_reason <- if (is.null(condition_col)) {
            "condition_col_not_supplied"
        } else if (!condition_col %in% colnames(metadata)) {
            "condition_col_absent"
        } else {
            "fewer_than_two_condition_levels"
        }
        return(object)
    }

    if (is.null(cell_type_col)) {
        stop("`cell_type_col` is required when multiple conditions are present.",
             call. = FALSE)
    }
    .condition_validate_labels(metadata, c(cell_type_col, condition_col))
    if (!is.numeric(min_cells_per_condition) ||
        length(min_cells_per_condition) != 1L ||
        !is.finite(min_cells_per_condition) || min_cells_per_condition < 3L ||
        min_cells_per_condition != as.integer(min_cells_per_condition)) {
        stop("`min_cells_per_condition` must be an integer >= 3.",
             call. = FALSE)
    }
    if (!is.numeric(padj_threshold) || length(padj_threshold) != 1L ||
        !is.finite(padj_threshold) || padj_threshold <= 0 ||
        padj_threshold >= 1) {
        stop("`padj_threshold` must be one number in (0, 1).",
             call. = FALSE)
    }

    available_types <- unique(as.character(metadata[[cell_type_col]]))
    requested_types <- if (is.null(cell_type)) {
        available_types
    } else unique(as.character(cell_type))
    missing_types <- setdiff(requested_types, available_types)
    if (length(missing_types)) {
        stop("Requested cell type(s) were not found: ",
             paste(missing_types, collapse = ", "), call. = FALSE)
    }

    prepared <- .condition_prepare_common_input(
        object = object, genes = genes,
        peak_to_gene_method = peak_to_gene_method,
        upstream = upstream, downstream = downstream, extend = extend,
        only_tss = only_tss,
        peak_to_gene_domains = peak_to_gene_domains,
        rna_layer = rna_layer, peak_layer = peak_layer,
        peak_value_type = peak_value_type, verbose = verbose
    )

    fits <- list()
    network_index <- list()
    for (type_label in requested_types) {
        type_cells <- rownames(metadata)[
            as.character(metadata[[cell_type_col]]) == type_label
        ]
        type_conditions <- unique(as.character(
            metadata[type_cells, condition_col]
        ))
        counts <- vapply(type_conditions, function(condition) {
            sum(as.character(metadata[type_cells, condition_col]) == condition)
        }, integer(1))
        eligible <- type_conditions[counts >= as.integer(min_cells_per_condition)]
        small <- setdiff(type_conditions, eligible)
        if (length(small)) {
            detail <- paste0(small, "=", counts[small], collapse = ", ")
            if (identical(small_condition_action, "error")) {
                stop("Cell type `", type_label,
                     "` has undersized condition(s): ", detail, call. = FALSE)
            }
            if (identical(small_condition_action, "skip_cell_type")) {
                log_message("Skipping cell type ", type_label,
                            " because condition(s) are undersized: ", detail,
                            verbose = verbose)
                next
            }
            log_message("Dropping undersized condition(s) for cell type ",
                        type_label, ": ", detail, verbose = verbose)
        }
        if (length(eligible) < 2L) {
            if (identical(small_condition_action, "error")) {
                stop("Cell type `", type_label,
                     "` retains fewer than two eligible conditions.",
                     call. = FALSE)
            }
            next
        }

        cells_by_condition <- stats::setNames(lapply(eligible, function(condition) {
            type_cells[as.character(metadata[type_cells, condition_col]) == condition]
        }), eligible)
        global_cells <- unlist(cells_by_condition, use.names = FALSE)
        log_message("Discovering global candidates for cell type ", type_label,
                    verbose = verbose)
        global_edges <- .condition_discover_edges_prepared(
            prepared, global_cells, source_label = "global",
            source_type = "global", tf_cor = tf_cor, peak_cor = peak_cor,
            parallel = parallel, verbose = verbose
        )
        condition_edges <- lapply(eligible, function(condition) {
            log_message("Discovering candidates for ", type_label, " / ",
                        condition, verbose = verbose)
            .condition_discover_edges_prepared(
                prepared, cells_by_condition[[condition]],
                source_label = condition, source_type = "condition",
                tf_cor = tf_cor, peak_cor = peak_cor,
                parallel = parallel, verbose = verbose
            )
        })
        names(condition_edges) <- eligible
        dictionary <- union_grn_edges(global_edges, condition_edges)

        network_names <- stats::setNames(vapply(eligible, function(condition) {
            paste0(
                network_name, "__", .condition_safe_label(type_label),
                "__condition__", .condition_safe_label(condition)
            )
        }, character(1)), eligible)
        conflicts <- intersect(network_names, names(object@grn@networks))
        if (length(conflicts) && !isTRUE(overwrite)) {
            stop("Network `", conflicts[[1L]],
                 "` already exists; set overwrite=TRUE.", call. = FALSE)
        }

        skeleton <- list(
            schema_version = .condition_common_dictionary_schema,
            model_schema = .condition_multitask_ridge_schema,
            fit_engine = "two_stage_exact_edge_union_multitask_ridge",
            coefficient_scale = "raw_tf_atac_interaction_units",
            internal_predictor_scale = "equal_condition_within_condition_rms",
            inference_scope =
                "approximate_ridge_wald_diagnostic_conditional_on_dictionary_cv_lambda_and_fusion",
            cell_type = type_label,
            condition_levels = eligible,
            condition_col = condition_col,
            cell_type_col = cell_type_col,
            condition_cell_ids = cells_by_condition,
            edge_dictionary = dictionary,
            coefficients = NULL,
            contrasts = NULL,
            fit = NULL,
            network_names = network_names,
            padj_threshold = padj_threshold,
            adjust_method = adjust_method,
            scale = FALSE,
            interaction = ":",
            projection_effect_column = "penalty_effect",
            projection_policy = "continuous_estimable_ridge_effects",
            target_genes = unique(as.character(dictionary$target)),
            rna_assay = prepared$params$rna_assay,
            atac_assay = prepared$params$peak_assay,
            rna_layer = prepared$rna_layer,
            peak_layer = prepared$peak_layer,
            peak_value_type = prepared$peak_value_type,
            preprocessing_fingerprint = prepared$preprocessing_fingerprint,
            ridge_control = control
        )
        class(skeleton) <- c("ConditionGRNFit", "list")

        refitted <- .condition_ridge_refit_contract(
            object = object, fit = skeleton, prepared = prepared,
            control = control, rank_action = rank_action,
            min_residual_df = min_residual_df,
            parallel = parallel, verbose = verbose
        )
        object <- refitted$object
        fits[[type_label]] <- refitted$fit

        for (condition in eligible) {
            one <- refitted$fit$coefficients[
                refitted$fit$coefficients$condition == condition,
                , drop = FALSE
            ]
            network_index[[length(network_index) + 1L]] <- data.frame(
                cell_type = type_label,
                condition = condition,
                network_name = network_names[[condition]],
                n_cells = length(cells_by_condition[[condition]]),
                n_dictionary_edges = nrow(dictionary),
                n_projection_edges = sum(
                    one$estimable %in% TRUE & is.finite(one$estimate)
                ),
                n_significant_diagnostic_edges = sum(one$significant %in% TRUE),
                stringsAsFactors = FALSE
            )
        }
    }

    if (!length(fits)) {
        stop("No cell type retained at least two eligible conditions.",
             call. = FALSE)
    }
    object@grn@params$analysis_mode <- "condition_grn"
    object@grn@params$condition_col <- condition_col
    object@grn@params$condition_levels <- condition_levels
    object@grn@params$cell_type_col <- cell_type_col
    object@grn@params$condition_coefficients_calculated <- TRUE
    object@grn@params$condition_grn_schema <-
        .condition_common_dictionary_schema
    object@grn@params$condition_grn_model_schema <-
        .condition_multitask_ridge_schema
    object@grn@params$condition_grn_method <-
        "two_stage_exact_edge_union_multitask_ridge"
    object@grn@params$condition_projection_policy <-
        "continuous_estimable_ridge_effects"
    object@grn@params$condition_ridge_control <- control
    object@grn@params$condition_grn_fits <- fits
    object@grn@params$condition_network_index <- do.call(rbind, network_index)
    object
}
