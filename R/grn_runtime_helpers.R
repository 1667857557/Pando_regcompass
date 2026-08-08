# Runtime helpers for standard and condition-GRN routing.

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

.pando_merge_cell_type_grn_results <- function(
    object, results, condition_col, cell_type_col) {
    modes <- character()
    condition_fits <- list()
    condition_index <- list()
    standard_index <- list()
    routing <- list()
    base_networks <- names(object@grn@networks)
    for (i in seq_along(results)) {
        value <- results[[i]]
        params <- value@grn@params
        mode <- if (is.null(params$analysis_mode)) {
  "unknown"
        } else {
  as.character(params$analysis_mode)
        }
        modes[[i]] <- mode
        for (network_name in setdiff(names(value@grn@networks), base_networks)) {
  if (network_name %in% names(object@grn@networks)) {
      stop("Duplicated network while merging parallel cell-type jobs: ",
           network_name, call. = FALSE)
  }
  object@grn@networks[[network_name]] <-
      value@grn@networks[[network_name]]
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
    object@grn@params$condition_coefficients_calculated <-
        length(condition_fits) > 0L
    object@grn@params$condition_grn_schema <- if (length(condition_fits)) {
        .condition_common_dictionary_schema
    } else NULL
    object@grn@params$condition_grn_fits <- condition_fits
    object@grn@params$condition_network_index <- bind_rows(condition_index)
    object@grn@params$standard_network_index <- bind_rows(standard_index)
    object@grn@params$cell_type_analysis_mode <- bind_rows(routing)

    condition_levels <- unique(unlist(lapply(results, function(value) {
        as.character(value@grn@params$condition_levels)
    }), use.names = FALSE))
    condition_levels <- condition_levels[
        !is.na(condition_levels) & nzchar(condition_levels)
    ]
    object@grn@params$condition_levels <- condition_levels

    fallback_reasons <- unique(unlist(lapply(results, function(value) {
        as.character(value@grn@params$standard_fallback_reason)
    }), use.names = FALSE))
    fallback_reasons <- fallback_reasons[
        !is.na(fallback_reasons) & nzchar(fallback_reasons)
    ]
    object@grn@params$standard_fallback_reason <- if (length(fallback_reasons)) {
        paste(fallback_reasons, collapse = ";")
    } else NULL

    added_networks <- setdiff(names(object@grn@networks), base_networks)
    if (length(added_networks)) {
        object@grn@active_network <- tail(added_networks, 1L)
    }
    object
}
