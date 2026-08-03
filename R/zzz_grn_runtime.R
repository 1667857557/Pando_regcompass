# Route standard-GRN arguments safely and parallelize independent cell-type jobs.

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
        warning("BiocParallel is unavailable; using serial cell-type execution.",
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

# Parallelizes independent cell-type routes. Each worker runs its own Pando job
# serially, preventing nested target-level parallelism and oversubscription.
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
    dots <- list(...)
    if (length(dots)) {
        label <- names(dots)
        label[!nzchar(label)] <- "<unnamed>"
        stop("Unused condition-GRN argument(s): ",
             paste(label, collapse = ", "), call. = FALSE)
    }
    parallel_scope <- match.arg(parallel_scope)
    metadata <- object@data@meta.data
    can_split <- isTRUE(parallel) &&
        !identical(parallel_scope, "target") &&
        is.character(cell_type_col) && length(cell_type_col) == 1L &&
        !is.na(cell_type_col) && cell_type_col %in% colnames(metadata)
    available <- if (can_split) {
        unique(as.character(metadata[[cell_type_col]]))
    } else character()
    requested <- if (!can_split) {
        character()
    } else if (is.null(cell_type)) {
        available
    } else {
        unique(as.character(cell_type))
    }
    if (!can_split || length(requested) <= 1L) {
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
            parallel = parallel, overwrite = overwrite,
            fallback_args = fallback_args, verbose = verbose
        )))
    }
    missing <- setdiff(requested, available)
    if (length(missing)) {
        stop("Requested cell type(s) were not found: ",
             paste(missing, collapse = ", "), call. = FALSE)
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
