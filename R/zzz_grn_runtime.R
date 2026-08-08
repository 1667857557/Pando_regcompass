# Route standard-GRN arguments safely and parallelize condition-GRN work.

.pando_condition_only_standard_args <- c(
    "padj_threshold", "rank_action", "min_residual_df",
    "rna_layer", "peak_layer", "peak_value_type",
    "min_cells_per_condition", "small_condition_action",
    "condition_col", "cell_type_col", "cell_type",
    "BPPARAM", "parallel_scope"
)

.pando_sanitize_standard_infer_dots <- function(dots, warn = TRUE) {
    if (!is.list(dots)) stop("Standard Pando inference arguments must be a list.",
                             call. = FALSE)
    if (!length(dots)) {
        return(list(args = dots, disabled = character()))
    }
    if (is.null(names(dots))) names(dots) <- rep("", length(dots))
    unnamed <- which(!nzchar(names(dots)))
    if (length(unnamed)) {
        stop("Every additional standard Pando argument must be named.",
             call. = FALSE)
    }
    disabled <- intersect(names(dots), .pando_condition_only_standard_args)
    if (length(disabled)) {
        dots[disabled] <- NULL
        if (isTRUE(warn)) {
            message(
                "Standard Pando disabled condition-only argument(s): ",
                paste(disabled, collapse = ", ")
            )
        }
    }
    list(args = dots, disabled = disabled)
}

.pando_validate_bpparam <- function(BPPARAM) {
    if (identical(BPPARAM, FALSE) || is.null(BPPARAM)) return(invisible(TRUE))
    if (is.logical(BPPARAM)) {
        stop("`BPPARAM` must be NULL, FALSE, or a BiocParallelParam object.",
             call. = FALSE)
    }
    if (!requireNamespace("BiocParallel", quietly = TRUE) ||
        !methods::is(BPPARAM, "BiocParallelParam")) {
        stop("A valid BiocParallelParam is required for `BPPARAM`.",
             call. = FALSE)
    }
    invisible(TRUE)
}

.pando_default_bpparam <- function(workers) {
    workers <- max(1L, as.integer(workers[[1L]]))
    if (workers < 2L || !requireNamespace("BiocParallel", quietly = TRUE)) {
        return(NULL)
    }
    if (identical(.Platform$OS.type, "windows")) {
        BiocParallel::SnowParam(workers = workers, type = "SOCK")
    } else {
        BiocParallel::MulticoreParam(workers = workers)
    }
}

.pando_parallel_lapply <- function(X, FUN, parallel = FALSE, BPPARAM = NULL) {
    if (!is.function(FUN)) stop("`FUN` must be a function.", call. = FALSE)
    .pando_validate_bpparam(BPPARAM)
    if (!isTRUE(parallel) || length(X) <= 1L || identical(BPPARAM, FALSE)) {
        return(lapply(X, FUN))
    }
    if (is.null(BPPARAM)) {
        detected <- suppressWarnings(parallel::detectCores(logical = TRUE))
        if (!is.finite(detected) || detected < 2L) detected <- 2L
        BPPARAM <- .pando_default_bpparam(min(length(X), detected - 1L))
    }
    if (is.null(BPPARAM)) {
        warning("BiocParallel is unavailable; using serial execution.",
                call. = FALSE)
        return(lapply(X, FUN))
    }
    BiocParallel::bplapply(X, FUN, BPPARAM = BPPARAM)
}

.pando_infer_grn_grndata_impl <- infer_grn.GRNData

# Standard infer_grn accepts model-specific dots, but condition-GRN controls must
# never be forwarded into stats::glm() or another model backend.
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
    routed <- .pando_sanitize_standard_infer_dots(list(...), warn = verbose)
    args <- list(
        object = object, genes = genes, network_name = network_name,
        peak_to_gene_method = peak_to_gene_method,
        upstream = upstream, downstream = downstream, extend = extend,
        only_tss = only_tss, parallel = parallel,
        tf_cor = tf_cor, peak_cor = peak_cor,
        aggregate_rna_col = aggregate_rna_col,
        aggregate_peaks_col = aggregate_peaks_col,
        method = method, alpha = alpha, family = family,
        interaction_term = interaction_term,
        adjust_method = adjust_method, scale = scale, verbose = verbose
    )
    do.call(.pando_infer_grn_grndata_impl, c(args, routed$args))
}

.pando_infer_condition_grn_grndata_impl <- infer_condition_grn.GRNData

.pando_merge_cell_type_grn_results <- function(object, results, condition_col,
                                                cell_type_col) {
    modes <- character()
    condition_fits <- list()
    condition_index <- list()
    standard_index <- list()
    routing <- list()
    base_networks <- names(object@grn@networks)
    for (i in seq_along(results)) {
        value <- results[[i]]
        params <- value@grn@params
        mode <- if (is.null(params$analysis_mode)) "unknown" else as.character(params$analysis_mode)
        modes[[i]] <- mode
        for (network_name in setdiff(names(value@grn@networks), base_networks)) {
            if (network_name %in% names(object@grn@networks)) {
                stop("Duplicated network while merging parallel cell-type jobs: ",
                     network_name, call. = FALSE)
            }
            object@grn@networks[[network_name]] <- value@grn@networks[[network_name]]
        }
        fits <- params$condition_grn_fits
        if (is.list(fits) && length(fits)) {
            condition_fits <- c(condition_fits, fits)
        }
        if (is.data.frame(params$condition_network_index) &&
            nrow(params$condition_network_index)) {
            condition_index[[length(condition_index) + 1L]] <-
                params$condition_network_index
        }
        if (is.data.frame(params$standard_network_index) &&
            nrow(params$standard_network_index)) {
            standard_index[[length(standard_index) + 1L]] <-
                params$standard_network_index
        }
        type_value <- if (is.data.frame(params$condition_network_index) &&
                          nrow(params$condition_network_index)) {
            unique(as.character(params$condition_network_index$cell_type))
        } else if (is.data.frame(params$standard_network_index) &&
                   nrow(params$standard_network_index)) {
            unique(as.character(params$standard_network_index$cell_type))
        } else {
            NA_character_
        }
        if (!length(type_value)) type_value <- NA_character_
        routing[[length(routing) + 1L]] <- data.frame(
            cell_type = type_value[[1L]], analysis_mode = mode,
            stringsAsFactors = FALSE
        )
    }
    bind_rows <- function(values) {
        if (!length(values)) return(data.frame())
        out <- do.call(rbind, values)
        rownames(out) <- NULL
        out
    }
    object@grn@params$analysis_mode <- if (length(unique(modes)) == 1L) {
        unique(modes)
    } else {
        "mixed_grn"
    }
    object@grn@params$condition_col <- condition_col
    object@grn@params$cell_type_col <- cell_type_col
    object@grn@params$condition_coefficients_calculated <- length(condition_fits) > 0L
    object@grn@params$condition_grn_schema <- if (length(condition_fits)) {
        .condition_common_dictionary_schema
    } else NULL
    object@grn@params$condition_grn_fits <- condition_fits
    object@grn@params$condition_network_index <- bind_rows(condition_index)
    object@grn@params$standard_network_index <- bind_rows(standard_index)
    object@grn@params$cell_type_analysis_mode <- bind_rows(routing)
    object
}

.pando_condition_celltype_plan <- function(
    metadata, cell_type_col, condition_col, cell_type,
    min_cells_per_condition, small_condition_action, verbose) {
    .condition_validate_labels(metadata, c(cell_type_col, condition_col))
    available_types <- unique(as.character(metadata[[cell_type_col]]))
    requested_types <- if (is.null(cell_type)) {
        available_types
    } else unique(as.character(cell_type))
    missing_types <- setdiff(requested_types, available_types)
    if (length(missing_types)) {
        stop("Requested cell type(s) were not found: ",
             paste(missing_types, collapse = ", "), call. = FALSE)
    }

    plan <- list()
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
        plan[[type_label]] <- list(
            cell_type = type_label,
            eligible = eligible,
            cells_by_condition = cells_by_condition,
            global_cells = unlist(cells_by_condition, use.names = FALSE)
        )
    }
    if (!length(plan)) {
        stop("No cell type retained at least two eligible conditions.",
             call. = FALSE)
    }
    plan
}

.pando_condition_celltype_parallel_impl <- function(
    object, cell_type_col, condition_col, cell_type, genes, network_name,
    peak_to_gene_method, upstream, downstream, extend, only_tss,
    peak_to_gene_domains, rna_layer, peak_layer, peak_value_type,
    tf_cor, peak_cor, min_cells_per_condition, small_condition_action,
    adjust_method, padj_threshold, rank_action, min_residual_df,
    BPPARAM, overwrite, verbose) {
    metadata <- object@data@meta.data
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
    plan <- .pando_condition_celltype_plan(
        metadata = metadata,
        cell_type_col = cell_type_col,
        condition_col = condition_col,
        cell_type = cell_type,
        min_cells_per_condition = min_cells_per_condition,
        small_condition_action = small_condition_action,
        verbose = verbose
    )
    prepared <- .condition_prepare_common_input(
        object = object, genes = genes,
        peak_to_gene_method = peak_to_gene_method,
        upstream = upstream, downstream = downstream, extend = extend,
        only_tss = only_tss,
        peak_to_gene_domains = peak_to_gene_domains,
        rna_layer = rna_layer, peak_layer = peak_layer,
        peak_value_type = peak_value_type, verbose = verbose
    )

    discovery_tasks <- list()
    for (type_label in names(plan)) {
        one <- plan[[type_label]]
        discovery_tasks[[length(discovery_tasks) + 1L]] <- list(
            cell_type = type_label,
            condition = NA_character_,
            source_label = "global",
            source_type = "global",
            cells = one$global_cells
        )
        for (condition in one$eligible) {
            discovery_tasks[[length(discovery_tasks) + 1L]] <- list(
                cell_type = type_label,
                condition = condition,
                source_label = condition,
                source_type = "condition",
                cells = one$cells_by_condition[[condition]]
            )
        }
    }
    discover_one <- function(task) {
        edges <- .condition_discover_edges_prepared(
            prepared = prepared,
            cells = task$cells,
            source_label = task$source_label,
            source_type = task$source_type,
            tf_cor = tf_cor,
            peak_cor = peak_cor,
            parallel = FALSE,
            verbose = FALSE
        )
        list(
            cell_type = task$cell_type,
            condition = task$condition,
            source_type = task$source_type,
            edges = edges
        )
    }
    discovery_results <- .pando_parallel_lapply(
        discovery_tasks, discover_one, parallel = TRUE, BPPARAM = BPPARAM
    )

    dictionaries <- list()
    for (type_label in names(plan)) {
        global_result <- Filter(function(x) {
            identical(x$cell_type, type_label) && identical(x$source_type, "global")
        }, discovery_results)
        condition_result <- Filter(function(x) {
            identical(x$cell_type, type_label) && identical(x$source_type, "condition")
        }, discovery_results)
        if (length(global_result) != 1L ||
            length(condition_result) != length(plan[[type_label]]$eligible)) {
            stop("Condition candidate discovery returned an incomplete task set.",
                 call. = FALSE)
        }
        condition_edges <- lapply(condition_result, `[[`, "edges")
        names(condition_edges) <- vapply(
            condition_result, `[[`, character(1), "condition"
        )
        condition_edges <- condition_edges[plan[[type_label]]$eligible]
        dictionaries[[type_label]] <- union_grn_edges(
            global_edges = global_result[[1L]]$edges,
            condition_edges = condition_edges
        )
    }

    fit_tasks <- list()
    for (type_label in names(plan)) {
        dictionary <- dictionaries[[type_label]]
        for (condition in plan[[type_label]]$eligible) {
            one_name <- paste0(
                network_name, "__", .condition_safe_label(type_label),
                "__condition__", .condition_safe_label(condition)
            )
            if (one_name %in% names(object@grn@networks) && !isTRUE(overwrite)) {
                stop("Network `", one_name, "` already exists.", call. = FALSE)
            }
            fit_tasks[[length(fit_tasks) + 1L]] <- list(
                cell_type = type_label,
                condition = condition,
                cells = plan[[type_label]]$cells_by_condition[[condition]],
                network_name = one_name,
                dictionary = dictionary
            )
        }
    }
    fit_one <- function(task) {
        fitted <- .condition_fit_dictionary_prepared(
            object = object,
            prepared = prepared,
            edge_dictionary = task$dictionary,
            cells = task$cells,
            condition_label = task$condition,
            network_name = task$network_name,
            adjust_method = adjust_method,
            padj_threshold = padj_threshold,
            rank_action = rank_action,
            min_residual_df = min_residual_df,
            parallel = FALSE,
            verbose = FALSE,
            overwrite = overwrite
        )
        list(
            cell_type = task$cell_type,
            condition = task$condition,
            network_name = task$network_name,
            network = fitted$network,
            coefficients = fitted$coefficients,
            fit = fitted$fit
        )
    }
    fit_results <- .pando_parallel_lapply(
        fit_tasks, fit_one, parallel = TRUE, BPPARAM = BPPARAM
    )

    fits <- list()
    network_index <- list()
    for (result in fit_results) {
        if (result$network_name %in% names(object@grn@networks) &&
            !isTRUE(overwrite)) {
            stop("Duplicated condition network while merging parallel jobs: ",
                 result$network_name, call. = FALSE)
        }
        object@grn@networks[[result$network_name]] <- result$network
        object@grn@active_network <- result$network_name
        dictionary <- dictionaries[[result$cell_type]]
        network_index[[length(network_index) + 1L]] <- data.frame(
            cell_type = result$cell_type,
            condition = result$condition,
            network_name = result$network_name,
            n_cells = length(plan[[result$cell_type]]$cells_by_condition[[result$condition]]),
            n_dictionary_edges = nrow(dictionary),
            n_significant_edges = sum(result$coefficients$significant),
            stringsAsFactors = FALSE
        )
    }

    for (type_label in names(plan)) {
        eligible <- plan[[type_label]]$eligible
        one_result <- fit_results[vapply(fit_results, function(x) {
            identical(x$cell_type, type_label)
        }, logical(1))]
        names(one_result) <- vapply(one_result, `[[`, character(1), "condition")
        one_result <- one_result[eligible]
        coefficient_table <- do.call(rbind, lapply(one_result, `[[`, "coefficients"))
        fit_table <- do.call(rbind, lapply(one_result, `[[`, "fit"))
        rownames(coefficient_table) <- rownames(fit_table) <- NULL
        dictionary <- dictionaries[[type_label]]
        network_names <- stats::setNames(
            vapply(one_result, `[[`, character(1), "network_name"), eligible
        )
        fit_contract <- list(
            schema_version = .condition_common_dictionary_schema,
            fit_engine = "two_stage_exact_edge_union_fixed_dictionary_glm",
            coefficient_scale = "shared_preprocessed_input_units_unscaled",
            inference_scope = "conditional_on_selected_edge_dictionary",
            cell_type = type_label,
            condition_levels = eligible,
            condition_col = condition_col,
            cell_type_col = cell_type_col,
            condition_cell_ids = plan[[type_label]]$cells_by_condition,
            edge_dictionary = dictionary,
            coefficients = coefficient_table,
            fit = fit_table,
            network_names = network_names,
            padj_threshold = padj_threshold,
            adjust_method = adjust_method,
            scale = FALSE,
            interaction = ":",
            projection_effect_column = "penalty_effect",
            projection_policy = "padj_significant_effects_only",
            target_genes = unique(as.character(dictionary$target)),
            rna_assay = prepared$params$rna_assay,
            atac_assay = prepared$params$peak_assay,
            rna_layer = prepared$rna_layer,
            peak_layer = prepared$peak_layer,
            peak_value_type = prepared$peak_value_type,
            preprocessing_fingerprint = prepared$preprocessing_fingerprint,
            dictionary_preprocessing_provenance_verified = isTRUE(attr(
                dictionary, "preprocessing_provenance_verified", exact = TRUE
            ))
        )
        class(fit_contract) <- c("ConditionGRNFit", "list")
        fits[[type_label]] <- fit_contract
    }

    object@grn@params$analysis_mode <- "condition_grn"
    object@grn@params$condition_col <- condition_col
    object@grn@params$condition_levels <- .condition_resolve_levels(
        metadata, condition_col
    )
    object@grn@params$cell_type_col <- cell_type_col
    object@grn@params$condition_coefficients_calculated <- TRUE
    object@grn@params$condition_grn_schema <- .condition_common_dictionary_schema
    object@grn@params$condition_grn_fits <- fits
    object@grn@params$condition_network_index <- do.call(rbind, network_index)
    object@grn@params$parallel_plan <- list(
        scope = "condition_x_cell_type",
        candidate_discovery_tasks = length(discovery_tasks),
        fixed_dictionary_fit_tasks = length(fit_tasks),
        nested_target_parallel = FALSE,
        stage_barrier =
            "candidate_discovery_then_exact_union_then_fixed_dictionary_fit"
    )
    object
}

# Condition mode supports three non-nested parallel granularities. `auto`
# selects condition x cell type for multi-condition fits; `cell_type` retains
# the previous broad-cell-type scheduler; `target` delegates to the original
# target-level Pando path.
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
    dots <- list(...)
    if (length(dots)) {
        label <- names(dots)
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
    .pando_validate_bpparam(BPPARAM)
    parallel_scope <- match.arg(parallel_scope)
    peak_to_gene_method <- match.arg(peak_to_gene_method)
    peak_value_type <- match.arg(peak_value_type)
    small_condition_action <- match.arg(small_condition_action)
    rank_action <- match.arg(rank_action)
    metadata <- object@data@meta.data
    condition_levels <- .condition_resolve_levels(metadata, condition_col)

    if (!isTRUE(parallel) || length(condition_levels) < 2L ||
        identical(parallel_scope, "target")) {
        return(do.call(.pando_infer_condition_grn_grndata_impl, list(
            object = object, cell_type_col = cell_type_col,
            condition_col = condition_col, cell_type = cell_type,
            genes = genes, network_name = network_name,
            peak_to_gene_method = peak_to_gene_method,
            upstream = upstream, downstream = downstream, extend = extend,
            only_tss = only_tss, peak_to_gene_domains = peak_to_gene_domains,
            rna_layer = rna_layer, peak_layer = peak_layer,
            peak_value_type = peak_value_type,
            tf_cor = tf_cor, peak_cor = peak_cor,
            min_cells_per_condition = min_cells_per_condition,
            small_condition_action = small_condition_action,
            adjust_method = adjust_method, padj_threshold = padj_threshold,
            rank_action = rank_action, min_residual_df = min_residual_df,
            parallel = isTRUE(parallel) && identical(parallel_scope, "target"),
            overwrite = overwrite,
            fallback_args = fallback_args, verbose = verbose
        )))
    }

    if (is.null(cell_type_col) || !is.character(cell_type_col) ||
        length(cell_type_col) != 1L || is.na(cell_type_col) ||
        !cell_type_col %in% colnames(metadata)) {
        stop("`cell_type_col` is required when multiple conditions are present.",
             call. = FALSE)
    }

    if (identical(parallel_scope, "auto") ||
        identical(parallel_scope, "condition_cell_type")) {
        return(.pando_condition_celltype_parallel_impl(
            object = object, cell_type_col = cell_type_col,
            condition_col = condition_col, cell_type = cell_type,
            genes = genes, network_name = network_name,
            peak_to_gene_method = peak_to_gene_method,
            upstream = upstream, downstream = downstream, extend = extend,
            only_tss = only_tss,
            peak_to_gene_domains = peak_to_gene_domains,
            rna_layer = rna_layer, peak_layer = peak_layer,
            peak_value_type = peak_value_type,
            tf_cor = tf_cor, peak_cor = peak_cor,
            min_cells_per_condition = min_cells_per_condition,
            small_condition_action = small_condition_action,
            adjust_method = adjust_method, padj_threshold = padj_threshold,
            rank_action = rank_action, min_residual_df = min_residual_df,
            BPPARAM = BPPARAM, overwrite = overwrite, verbose = verbose
        ))
    }

    available <- unique(as.character(metadata[[cell_type_col]]))
    requested <- if (is.null(cell_type)) available else unique(as.character(cell_type))
    missing <- setdiff(requested, available)
    if (length(missing)) {
        stop("Requested cell type(s) were not found: ",
             paste(missing, collapse = ", "), call. = FALSE)
    }
    if (length(requested) <= 1L) {
        return(do.call(.pando_infer_condition_grn_grndata_impl, list(
            object = object, cell_type_col = cell_type_col,
            condition_col = condition_col, cell_type = cell_type,
            genes = genes, network_name = network_name,
            peak_to_gene_method = peak_to_gene_method,
            upstream = upstream, downstream = downstream, extend = extend,
            only_tss = only_tss, peak_to_gene_domains = peak_to_gene_domains,
            rna_layer = rna_layer, peak_layer = peak_layer,
            peak_value_type = peak_value_type,
            tf_cor = tf_cor, peak_cor = peak_cor,
            min_cells_per_condition = min_cells_per_condition,
            small_condition_action = small_condition_action,
            adjust_method = adjust_method, padj_threshold = padj_threshold,
            rank_action = rank_action, min_residual_df = min_residual_df,
            parallel = FALSE, overwrite = overwrite,
            fallback_args = fallback_args, verbose = verbose
        )))
    }
    run_one <- function(type_label) {
        cells <- rownames(metadata)[
            as.character(metadata[[cell_type_col]]) == type_label
        ]
        one <- object
        one@data <- subset(object@data, cells = cells)
        do.call(.pando_infer_condition_grn_grndata_impl, list(
            object = one, cell_type_col = cell_type_col,
            condition_col = condition_col, cell_type = type_label,
            genes = genes, network_name = network_name,
            peak_to_gene_method = peak_to_gene_method,
            upstream = upstream, downstream = downstream, extend = extend,
            only_tss = only_tss, peak_to_gene_domains = peak_to_gene_domains,
            rna_layer = rna_layer, peak_layer = peak_layer,
            peak_value_type = peak_value_type,
            tf_cor = tf_cor, peak_cor = peak_cor,
            min_cells_per_condition = min_cells_per_condition,
            small_condition_action = small_condition_action,
            adjust_method = adjust_method, padj_threshold = padj_threshold,
            rank_action = rank_action, min_residual_df = min_residual_df,
            parallel = FALSE, overwrite = overwrite,
            fallback_args = fallback_args, verbose = verbose
        ))
    }
    results <- .pando_parallel_lapply(
        requested, run_one, parallel = TRUE, BPPARAM = BPPARAM
    )
    object <- .pando_merge_cell_type_grn_results(
        object, results, condition_col = condition_col,
        cell_type_col = cell_type_col
    )
    object@grn@params$parallel_plan <- list(
        scope = "cell_type", n_jobs = length(requested),
        nested_parallel = FALSE, cell_types = requested
    )
    object
}
