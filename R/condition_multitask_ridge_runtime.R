# Global-plus-condition common-dictionary ridge integration behind the canonical
# infer_condition_grn.GRNData definition in condition_grn.R.

.condition_candidate_support_table <- function(
    global_edges = NULL, condition_edges) {
    required <- c(
        "target", "tf", "region", "atac_feature_id",
        "peak_target_cor", "tf_target_cor"
    )
    make_rows <- function(one, source_type, condition = NA_character_) {
        one <- as.data.frame(one, stringsAsFactors = FALSE)
        if (!nrow(one)) return(NULL)
        if (!all(required %in% colnames(one))) {
            stop("Pando candidate table is missing support columns.",
                 call. = FALSE)
        }
        out <- one[, required, drop = FALSE]
        out$edge_id <- paste(out$target, out$tf, out$region, sep = "||")
        out$source_type <- source_type
        out$condition <- condition
        out[, c(
            "edge_id", "target", "tf", "region", "atac_feature_id",
            "source_type", "condition", "peak_target_cor", "tf_target_cor"
        ), drop = FALSE]
    }
    rows <- list()
    if (!is.null(global_edges)) {
        rows[[length(rows) + 1L]] <- make_rows(global_edges, "global")
    }
    for (condition in names(condition_edges)) {
        rows[[length(rows) + 1L]] <- make_rows(
            condition_edges[[condition]], "condition", condition
        )
    }
    rows <- rows[!vapply(rows, is.null, logical(1))]
    if (!length(rows)) {
        return(data.frame(
            edge_id = character(), target = character(), tf = character(),
            region = character(), atac_feature_id = character(),
            source_type = character(), condition = character(),
            peak_target_cor = numeric(), tf_target_cor = numeric(),
            stringsAsFactors = FALSE
        ))
    }
    out <- do.call(rbind, rows)
    key <- paste(
        out$edge_id, out$source_type,
        ifelse(is.na(out$condition), "__global__", out$condition), sep = "\001"
    )
    out <- out[!duplicated(key), , drop = FALSE]
    out <- out[order(
        out$edge_id, out$source_type,
        ifelse(is.na(out$condition), "", out$condition)
    ), , drop = FALSE]
    rownames(out) <- NULL
    out
}

.condition_candidate_support_summary <- function(dictionary, support) {
    ids <- as.character(dictionary$edge_id)
    summary_rows <- lapply(ids, function(id) {
        one <- support[as.character(support$edge_id) == id, , drop = FALSE]
        global <- one$source_type == "global"
        condition <- one$source_type == "condition"
        supported_conditions <- sort(unique(
            as.character(one$condition[condition & !is.na(one$condition)])
        ))
        data.frame(
            edge_id = id,
            source_global = any(global),
            n_support_conditions = length(supported_conditions),
            support_conditions = paste(supported_conditions, collapse = ";"),
            n_sources = as.integer(any(global)) + length(supported_conditions),
            max_abs_peak_target_cor = if (nrow(one)) {
                max(abs(as.numeric(one$peak_target_cor)))
            } else NA_real_,
            max_abs_tf_target_cor = if (nrow(one)) {
                max(abs(as.numeric(one$tf_target_cor)))
            } else NA_real_,
            stringsAsFactors = FALSE
        )
    })
    summary <- do.call(rbind, summary_rows)
    rownames(summary) <- NULL
    summary
}

.pando_infer_condition_grn_multitask_ridge_one <- function(
    object, cell_type_col = NULL, condition_col = NULL, cell_type = NULL,
    genes = NULL, network_name = "condition_grn",
    peak_to_gene_method = c("Signac", "GREAT"), upstream = 100000,
    downstream = 0, extend = 1000000, only_tss = FALSE,
    peak_to_gene_domains = NULL, rna_layer = "data", peak_layer = "data",
    peak_value_type = c("normalized", "probability", "other"),
    tf_cor = 0.05, peak_cor = 0.05,
    min_cells_per_condition = 50L,
    small_condition_action = c("error", "drop_condition", "skip_cell_type"),
    adjust_method = "BH", padj_threshold = 0.05,
    rank_action = c("mark", "error"), min_residual_df = 1L,
    condition_ridge_control = list(),
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
    if (!is.list(condition_ridge_control)) {
        stop("`condition_ridge_control` must be a list.", call. = FALSE)
    }
    if (!is.list(fallback_args)) {
        stop("`fallback_args` must be a list.", call. = FALSE)
    }
    if ("condition_ridge_control" %in% names(fallback_args)) {
        stop(
            "`fallback_args$condition_ridge_control` has been removed; pass ",
            "`condition_ridge_control` directly.", call. = FALSE
        )
    }
    control <- .condition_ridge_control(condition_ridge_control)
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
    .condition_validate_adjust_method(adjust_method)
    padj_threshold <- .condition_validate_padj_threshold(padj_threshold)

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
        global_cells <- unique(unlist(cells_by_condition, use.names = FALSE))
        log_message(
            "Discovering pooled/global Pando candidates for cell type ",
            type_label, verbose = verbose
        )
        global_edges <- .condition_discover_edges_compact(
            prepared, global_cells,
            source_label = "global", source_type = "global",
            tf_cor = tf_cor, peak_cor = peak_cor,
            parallel = parallel, verbose = verbose
        )
        condition_edges <- lapply(eligible, function(condition) {
            log_message("Discovering Pando candidates for ", type_label, " / ",
                        condition, verbose = verbose)
            .condition_discover_edges_compact(
                prepared, cells_by_condition[[condition]],
                source_label = condition, source_type = "condition",
                tf_cor = tf_cor, peak_cor = peak_cor,
                parallel = parallel, verbose = verbose
            )
        })
        names(condition_edges) <- eligible

        dictionary <- union_grn_edges(
            global_edges = global_edges, condition_edges = condition_edges
        )
        if (anyDuplicated(as.character(dictionary$edge_id))) {
            stop("Common condition dictionary contains duplicated exact edges.",
                 call. = FALSE)
        }
        support_table <- .condition_candidate_support_table(
            global_edges = global_edges, condition_edges = condition_edges
        )
        support_summary <- .condition_candidate_support_summary(
            dictionary, support_table
        )
        if (!nrow(support_table) ||
            !setequal(as.character(dictionary$edge_id),
                      as.character(support_summary$edge_id)) ||
            any(support_summary$n_sources < 1L)) {
            stop("Global/condition Pando support audit is inconsistent.",
                 call. = FALSE)
        }
        index <- match(dictionary$edge_id, support_summary$edge_id)
        if (!identical(
            as.logical(dictionary$source_global),
            as.logical(support_summary$source_global[index])
        )) {
            stop("Dictionary global-support provenance is inconsistent.",
                 call. = FALSE)
        }

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
            fit_engine = .condition_fit_engine,
            coefficient_scale = "raw_tf_atac_interaction_units",
            internal_predictor_scale = "equal_condition_within_condition_rms",
            inference_scope =
                "approximate_ridge_wald_conditional_on_global_or_condition_pando_screened_dictionary_and_cv_lambda",
            cell_type = type_label,
            condition_levels = eligible,
            condition_col = condition_col,
            cell_type_col = cell_type_col,
            condition_cell_ids = cells_by_condition,
            edge_dictionary = dictionary,
            dictionary_support_table = support_table,
            dictionary_support_summary = support_summary,
            candidate_tf_cor = as.numeric(tf_cor),
            candidate_peak_cor = as.numeric(peak_cor),
            coefficients = NULL,
            contrasts = NULL,
            fit = NULL,
            network_names = network_names,
            padj_threshold = padj_threshold,
            adjust_method = "BH",
            scale = FALSE,
            interaction = ":",
            projection_effect_column = "penalty_effect",
            projection_policy = .condition_significant_projection_policy,
            fit_dictionary_policy = .condition_fit_dictionary_policy,
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

        fitted <- .condition_ridge_fit_contract(
            object = object, fit = skeleton, prepared = prepared,
            control = control, rank_action = rank_action,
            min_residual_df = min_residual_df,
            parallel = parallel, verbose = verbose
        )
        object <- fitted$object
        fits[[type_label]] <- fitted$fit

        for (condition in eligible) {
            one <- fitted$fit$coefficients[
                fitted$fit$coefficients$condition == condition,
                , drop = FALSE
            ]
            network_index[[length(network_index) + 1L]] <- data.frame(
                cell_type = type_label,
                condition = condition,
                network_name = network_names[[condition]],
                n_cells = length(cells_by_condition[[condition]]),
                n_dictionary_edges = nrow(fitted$fit$edge_dictionary),
                n_statistically_supported_edges =
                    sum(one$statistically_supported %in% TRUE),
                n_locally_supported_edges = sum(one$local_support %in% TRUE),
                n_global_supported_edges = sum(one$global_support %in% TRUE),
                n_projection_edges = sum(one$active %in% TRUE),
                n_active_edges = sum(one$active %in% TRUE),
                n_significant_edges = sum(one$active %in% TRUE),
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
        .condition_fit_dictionary_policy
    object@grn@params$condition_projection_policy <-
        .condition_significant_projection_policy
    object@grn@params$condition_ridge_control <- control
    object@grn@params$condition_grn_fits <- fits
    object@grn@params$condition_network_index <- do.call(rbind, network_index)
    object
}
