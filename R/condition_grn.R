# Condition-aware extension of the Pando TF-peak-target model.

#' Infer a GRN with automatic standard/condition-aware routing
#'
#' Uses the original [infer_grn()] workflow when `condition_col` is omitted,
#' absent from the metadata, or contains fewer than two non-missing levels.
#' With two or more conditions, fits condition-specific sub-GRNs on one shared
#' within-cell-type coordinate system. The fallback path does not construct a
#' `ConditionGRNFit` and does not calculate condition coefficients.
#'
#' @param object A GRNData object containing paired RNA and ATAC measurements.
#' @param cell_type_col Optional metadata column containing broad cell-type
#'   labels. It is required for condition-aware fitting. In standard fallback
#'   mode it is used only when `cell_type` is supplied.
#' @param condition_col Optional metadata column containing condition labels.
#' @param cell_type Optional cell-type label or labels to retain before fitting.
#' @param genes Target genes. Defaults to RNA variable features.
#' @param network_name Prefix used for generated network names.
#' @param peak_to_gene_method Peak-to-gene method, either Signac or GREAT.
#' @param upstream,downstream,extend,only_tss Regulatory-domain parameters.
#' @param peak_to_gene_domains Optional custom gene regulatory domains for the
#'   condition-aware path.
#' @param tf_cor,peak_cor Correlation screening thresholds.
#' @param candidate_screen Condition-aware candidate screening strategy.
#' @param alpha Elastic-net mixing parameter.
#' @param condition_mix Group-versus-condition sparsity mixing parameter.
#' @param comparison_conditions Conditions defining common-support projection.
#' @param condition_weight Must be `equal` in condition-aware mode.
#' @param nlambda,lambda,lambda_min_ratio Lambda-path controls.
#' @param outer_nfolds,inner_nfolds Nested cross-fitting fold counts.
#' @param lambda_selection Selection of lambda.1se or lambda.min.
#' @param min_cells_per_condition Minimum cells per condition.
#' @param small_condition_action Skip, drop, or error for small conditions.
#' @param scale Must be `TRUE` for condition-aware fitting. The standard
#'   fallback uses the original Pando default unless overridden in
#'   `fallback_args`.
#' @param active_tol Numerical activity threshold.
#' @param parallel Use the existing Pando foreach backend.
#' @param BPPARAM Optional BiocParallel parameter for condition-aware fitting.
#' @param overwrite Replace generated condition networks with matching names.
#' @param seed Random seed for condition-stratified folds.
#' @param max_iter,tol_objective,tol_coef Solver controls.
#' @param fallback_args Named arguments passed only to [infer_grn()] in standard
#'   fallback mode. Common managed fields cannot be overridden.
#' @param verbose Display progress messages.
#' @param ... Must be empty.
#' @return A GRNData object. `object@grn@params$analysis_mode` is either
#'   `"standard_grn"` or `"condition_grn"`.
#' @export
infer_condition_grn <- function(object, ...) {
    UseMethod(generic = 'infer_condition_grn', object = object)
}

#' @rdname infer_condition_grn
#' @method infer_condition_grn GRNData
#' @export
infer_condition_grn.GRNData <- function(
    object,
    cell_type_col = NULL,
    condition_col = NULL,
    cell_type = NULL,
    genes = NULL,
    network_name = 'condition_grn',
    peak_to_gene_method = c('Signac', 'GREAT'),
    upstream = 100000,
    downstream = 0,
    extend = 1000000,
    only_tss = FALSE,
    peak_to_gene_domains = NULL,
    tf_cor = 0.1,
    peak_cor = 0,
    candidate_screen = c('motif_domain', 'pooled_within_condition'),
    alpha = 0.5,
    condition_mix = 0.5,
    comparison_conditions = NULL,
    condition_weight = 'equal',
    nlambda = 50L,
    lambda = NULL,
    lambda_min_ratio = NULL,
    outer_nfolds = 5L,
    inner_nfolds = 5L,
    lambda_selection = c('lambda.1se', 'lambda.min'),
    min_cells_per_condition = 50L,
    small_condition_action = c('skip_cell_type', 'drop_condition', 'error'),
    scale = TRUE,
    active_tol = 1e-8,
    parallel = FALSE,
    BPPARAM = NULL,
    overwrite = FALSE,
    seed = 12345L,
    max_iter = 5000L,
    tol_objective = 1e-7,
    tol_coef = 1e-6,
    fallback_args = list(),
    verbose = TRUE,
    ...
) {
    dots <- list(...)
    if (length(dots)) {
        dot_names <- names(dots)
        if (is.null(dot_names) || anyNA(dot_names) || any(!nzchar(dot_names))) {
            dot_names <- paste0('..', seq_along(dots))
        }
        stop('Unused condition GRN argument(s): ', paste(dot_names, collapse = ', '), '.')
    }
    if (!inherits(object, 'GRNData')) stop('object must be a GRNData object.')
    if (!is.list(fallback_args)) stop('fallback_args must be a list.')

    peak_to_gene_method <- match.arg(peak_to_gene_method)
    metadata <- object@data@meta.data
    condition_levels <- .condition_resolve_levels(metadata, condition_col)
    analysis_mode <- if (length(condition_levels) >= 2L) {
        'condition_grn'
    } else {
        'standard_grn'
    }

    if (identical(analysis_mode, 'standard_grn')) {
        value <- .condition_run_standard_fallback(
            object = object,
            metadata = metadata,
            cell_type_col = cell_type_col,
            cell_type = cell_type,
            condition_col = condition_col,
            condition_levels = condition_levels,
            genes = genes,
            network_name = network_name,
            peak_to_gene_method = peak_to_gene_method,
            upstream = upstream,
            downstream = downstream,
            extend = extend,
            only_tss = only_tss,
            parallel = parallel,
            tf_cor = tf_cor,
            peak_cor = peak_cor,
            alpha = alpha,
            fallback_args = fallback_args,
            verbose = verbose
        )
        value@grn@params$analysis_mode <- 'standard_grn'
        value@grn@params$condition_col <- condition_col
        value@grn@params$condition_levels <- condition_levels
        value@grn@params$condition_coefficients_calculated <- FALSE
        value@grn@params$standard_fallback_reason <- if (is.null(condition_col)) {
            'condition_col_not_supplied'
        } else if (!condition_col %in% colnames(metadata)) {
            'condition_col_absent'
        } else {
            'fewer_than_two_condition_levels'
        }
        return(value)
    }

    candidate_screen <- match.arg(candidate_screen)
    lambda_selection <- match.arg(lambda_selection)
    small_condition_action <- match.arg(small_condition_action)
    if (!identical(condition_weight, 'equal')) {
        stop('condition_weight must be "equal" for condition comparability.')
    }
    if (!isTRUE(scale)) {
        stop('scale must be TRUE for condition-aware fitting.')
    }
    .condition_validate_public_args(
        object, cell_type_col, condition_col, scale, alpha, condition_mix,
        comparison_conditions, nlambda, lambda, outer_nfolds, inner_nfolds,
        min_cells_per_condition, active_tol, max_iter, tol_objective, tol_coef,
        tf_cor, peak_cor, lambda_min_ratio, parallel, overwrite, seed, verbose
    )
    prepared <- .condition_prepare_global_data(
        object = object,
        cell_type_col = cell_type_col,
        condition_col = condition_col,
        genes = genes,
        peak_to_gene_method = peak_to_gene_method,
        upstream = upstream,
        downstream = downstream,
        extend = extend,
        only_tss = only_tss,
        peak_to_gene_domains = peak_to_gene_domains,
        verbose = verbose
    )
    value <- .condition_fit_cell_type_models(
        object = object,
        prepared = prepared,
        cell_type_col = cell_type_col,
        cell_type = cell_type,
        condition_col = condition_col,
        network_name = network_name,
        peak_to_gene_method = peak_to_gene_method,
        upstream = upstream,
        downstream = downstream,
        extend = extend,
        only_tss = only_tss,
        tf_cor = tf_cor,
        peak_cor = peak_cor,
        candidate_screen = candidate_screen,
        alpha = alpha,
        condition_mix = condition_mix,
        reference_condition = NULL,
        comparison_conditions = comparison_conditions,
        condition_weight = condition_weight,
        nlambda = nlambda,
        lambda = lambda,
        lambda_min_ratio = lambda_min_ratio,
        outer_nfolds = outer_nfolds,
        inner_nfolds = inner_nfolds,
        lambda_selection = lambda_selection,
        min_cells_per_condition = min_cells_per_condition,
        small_condition_action = small_condition_action,
        active_tol = active_tol,
        parallel = parallel,
        BPPARAM = BPPARAM,
        overwrite = overwrite,
        seed = seed,
        max_iter = max_iter,
        tol_objective = tol_objective,
        tol_coef = tol_coef,
        verbose = verbose
    )
    value@grn@params$analysis_mode <- 'condition_grn'
    value@grn@params$condition_coefficients_calculated <- TRUE
    value
}

.condition_resolve_levels <- function(metadata, condition_col) {
    if (is.null(condition_col) || !is.character(condition_col) ||
        length(condition_col) != 1L || is.na(condition_col) ||
        !nzchar(trimws(condition_col)) || !condition_col %in% colnames(metadata)) {
        return(character())
    }
    values <- trimws(as.character(metadata[[condition_col]]))
    values <- values[!is.na(values) & nzchar(values)]
    unique(values)
}

.condition_run_standard_fallback <- function(
    object, metadata, cell_type_col, cell_type, condition_col, condition_levels,
    genes, network_name, peak_to_gene_method, upstream, downstream, extend,
    only_tss, parallel, tf_cor, peak_cor, alpha, fallback_args, verbose
) {
    if (!is.null(cell_type)) {
        if (is.null(cell_type_col) || !is.character(cell_type_col) ||
            length(cell_type_col) != 1L || !cell_type_col %in% colnames(metadata)) {
            stop('cell_type_col must identify a metadata column when cell_type is supplied.')
        }
        requested <- unique(trimws(as.character(cell_type)))
        available <- unique(as.character(metadata[[cell_type_col]]))
        missing <- setdiff(requested, available)
        if (length(missing)) {
            stop('Requested cell_type value(s) were not found: ', paste(missing, collapse = ', '), '.')
        }
        cells <- rownames(metadata)[as.character(metadata[[cell_type_col]]) %in% requested]
        object@data <- subset(object@data, cells = cells)
    }
    managed <- c(
        'object', 'genes', 'network_name', 'peak_to_gene_method', 'upstream',
        'downstream', 'extend', 'only_tss', 'parallel', 'tf_cor', 'peak_cor',
        'alpha', 'verbose'
    )
    conflict <- intersect(names(fallback_args), managed)
    if (length(conflict)) {
        stop('fallback_args cannot override managed fields: ', paste(conflict, collapse = ', '), '.')
    }
    defaults <- list(
        object = object,
        genes = genes,
        network_name = network_name,
        peak_to_gene_method = peak_to_gene_method,
        upstream = upstream,
        downstream = downstream,
        extend = extend,
        only_tss = only_tss,
        parallel = parallel,
        tf_cor = tf_cor,
        peak_cor = peak_cor,
        alpha = alpha,
        verbose = verbose
    )
    do.call(infer_grn, c(defaults, fallback_args))
}

.condition_validate_public_args <- function(
    object, cell_type_col, condition_col, scale, alpha, condition_mix,
    comparison_conditions, nlambda, lambda, outer_nfolds, inner_nfolds,
    min_cells_per_condition, active_tol, max_iter, tol_objective, tol_coef,
    tf_cor, peak_cor, lambda_min_ratio, parallel, overwrite, seed, verbose
) {
    if (!inherits(object, 'GRNData')) stop('object must be a GRNData object.')
    metadata <- object@data@meta.data
    .condition_validate_analysis_metadata(metadata, cell_type_col, condition_col)
    if (!isTRUE(scale)) stop('Condition GRN comparability requires scale = TRUE.')
    if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
        alpha < 0 || alpha > 1) stop('alpha must be between 0 and 1.')
    if (!is.numeric(condition_mix) || length(condition_mix) != 1L ||
        !is.finite(condition_mix) || condition_mix <= 0 || condition_mix > 1) {
        stop('condition_mix must be in (0, 1].')
    }
    if (!is.null(comparison_conditions)) {
        comparison_conditions <- as.character(comparison_conditions)
        if (length(comparison_conditions) < 2L || anyNA(comparison_conditions) ||
            any(!nzchar(trimws(comparison_conditions))) ||
            anyDuplicated(comparison_conditions)) {
            stop('comparison_conditions must contain at least two unique labels.')
        }
    }
    if (!is.numeric(nlambda) || length(nlambda) != 1L || !is.finite(nlambda) ||
        nlambda < 1L || nlambda != as.integer(nlambda)) stop('nlambda must be a positive integer.')
    if (is.null(lambda) && nlambda < 2L) stop('nlambda must be at least 2 when lambda is generated.')
    if (!is.null(lambda) && (!is.numeric(lambda) || !length(lambda) ||
        any(!is.finite(lambda)) || any(lambda < 0))) stop('lambda must contain finite non-negative values.')
    folds <- c(outer_nfolds, inner_nfolds)
    if (any(!is.finite(folds)) || any(folds < 2L) || any(folds != as.integer(folds))) {
        stop('outer_nfolds and inner_nfolds must be integers >= 2.')
    }
    if (!is.numeric(min_cells_per_condition) || length(min_cells_per_condition) != 1L ||
        !is.finite(min_cells_per_condition) || min_cells_per_condition < outer_nfolds ||
        min_cells_per_condition != as.integer(min_cells_per_condition)) {
        stop('min_cells_per_condition must be an integer >= outer_nfolds.')
    }
    correlation_thresholds <- c(tf_cor, peak_cor)
    if (any(!is.finite(correlation_thresholds)) ||
        any(correlation_thresholds < 0 | correlation_thresholds > 1)) {
        stop('tf_cor and peak_cor must be in [0, 1].')
    }
    if (!is.null(lambda_min_ratio) && (!is.numeric(lambda_min_ratio) ||
        length(lambda_min_ratio) != 1L || !is.finite(lambda_min_ratio) ||
        lambda_min_ratio <= 0 || lambda_min_ratio >= 1)) {
        stop('lambda_min_ratio must be NULL or in (0, 1).')
    }
    if (!is.numeric(active_tol) || length(active_tol) != 1L ||
        !is.finite(active_tol) || active_tol < 0) stop('active_tol must be non-negative.')
    if (!is.numeric(max_iter) || length(max_iter) != 1L || !is.finite(max_iter) ||
        max_iter < 1L || max_iter != as.integer(max_iter) ||
        !is.numeric(tol_objective) || length(tol_objective) != 1L ||
        !is.finite(tol_objective) || tol_objective <= 0 ||
        !is.numeric(tol_coef) || length(tol_coef) != 1L ||
        !is.finite(tol_coef) || tol_coef <= 0) {
        stop('Solver iteration and tolerance parameters are invalid.')
    }
    logical_controls <- list(parallel = parallel, overwrite = overwrite, verbose = verbose)
    if (any(vapply(logical_controls, function(x) !is.logical(x) || length(x) != 1L || is.na(x), logical(1)))) {
        stop('parallel, overwrite, and verbose must be TRUE or FALSE.')
    }
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
        seed != as.integer(seed)) stop('seed must be one finite integer.')
    invisible(TRUE)
}

.condition_validate_analysis_metadata <- function(metadata, cell_type_col, condition_col) {
    columns <- c(cell_type_col, condition_col)
    if (anyNA(columns) || any(!nzchar(trimws(columns))) ||
        length(unique(columns)) != 2L) {
        stop('cell_type_col and condition_col must be different non-empty column names.')
    }
    missing <- setdiff(columns, colnames(metadata))
    if (length(missing)) stop('Missing metadata column(s): ', paste(missing, collapse = ', '), '.')
    labels <- lapply(metadata[, columns, drop = FALSE], as.character)
    invalid <- vapply(labels, function(x) anyNA(x) || any(!nzchar(trimws(x))) || any(x != trimws(x)), logical(1))
    if (any(invalid)) {
        stop('Cell-type and condition labels must be complete and free of surrounding whitespace.')
    }
    if (length(unique(labels[[2L]])) < 2L) {
        stop('Condition-aware fitting requires at least two condition levels.')
    }
    invisible(TRUE)
}
