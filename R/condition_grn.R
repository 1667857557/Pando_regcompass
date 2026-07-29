# Condition-aware extension of the Pando TF-peak-target model.

#' Infer condition-specific sub-GRNs on one shared coordinate system
#'
#' Extends Pando's TF-expression by motif-bearing peak-accessibility model with
#' independent condition-aware workflows within broad cell types. Training,
#' validation and final refitting for a cell type use only cells from that type.
#' Within each type, OOF folds are stratified by condition.
#'
#' @param object A GRNData object containing paired RNA and ATAC measurements.
#' @param cell_type_col Metadata column containing cell-type labels.
#' @param cell_type Optional cell-type label or labels to fit. When omitted,
#' every observed cell type is fitted independently.
#' @param condition_col Metadata column containing condition labels.
#' @param genes Target genes. Defaults to RNA variable features.
#' @param network_name Prefix used for generated network names.
#' @param peak_to_gene_method Peak-to-gene method, either Signac or GREAT.
#' @param upstream,downstream,extend,only_tss Regulatory-domain parameters.
#' @param peak_to_gene_domains Optional custom gene regulatory domains.
#' @param tf_cor,peak_cor Correlation screening thresholds.
#' @param candidate_screen Candidate screening strategy. The default
#' `motif_domain` retains the response-independent shared motif/domain candidate
#' graph. `pooled_within_condition` is a response-dependent sensitivity screen.
#' @param alpha Elastic-net mixing parameter.
#' @param condition_mix Group-versus-condition sparsity mixing parameter.
#' @param reference_condition Condition used for explicit coefficient contrasts.
#' Defaults to the first condition level.
#' @param condition_weight Equal-condition or cell-count loss weighting.
#' @param nlambda Number of lambda values when lambda is not supplied.
#' @param lambda Optional fixed lambda value or decreasing lambda path.
#' @param lambda_min_ratio Smallest lambda relative to lambda maximum.
#' @param nfolds Number of condition-stratified cell folds within each cell
#' type.
#' @param lambda_selection Selection of lambda.1se or lambda.min.
#' @param min_cells_per_condition Minimum single cells per condition.
#' @param small_condition_action Skip the cell type, drop the small condition,
#' or stop.
#' @param scale Must be `TRUE`; all conditions use one pooled predictor and
#' response scale.
#' @param active_tol Numerical threshold used for network summaries.
#' @param parallel Use the existing Pando foreach mapping backend.
#' @param BPPARAM Optional BiocParallel parameter object.
#' @param overwrite Replace generated networks with matching names.
#' @param seed Random seed for condition-stratified folds.
#' @param max_iter,tol_objective,tol_coef Solver controls.
#' @param verbose Display progress messages.
#' @param ... Must be empty.
#' @return A GRNData object containing standard Pando Network objects and
#' versioned ConditionGRNFit contracts retrievable with condition_grn_fit().
#' @export
infer_condition_grn <- function(object, ...) {
    UseMethod(generic = 'infer_condition_grn', object = object)
}

#' @rdname infer_condition_grn
#' @method infer_condition_grn GRNData
#' @export
infer_condition_grn.GRNData <- function(
    object,
    cell_type_col,
    condition_col,
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
    reference_condition = NULL,
    condition_weight = c('equal', 'cell_count'),
    nlambda = 50L,
    lambda = NULL,
    lambda_min_ratio = NULL,
    nfolds = 5L,
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
    verbose = TRUE,
    ...
) {
    dots <- list(...)
    if (length(dots)) {
        dot_names <- names(dots)
        if (is.null(dot_names) || anyNA(dot_names) || any(!nzchar(dot_names))) {
            dot_names <- paste0('..', seq_along(dots))
        }
        stop(
            'Unused condition sub-GRN argument(s): ',
            paste(dot_names, collapse = ', '), '.'
        )
    }
    peak_to_gene_method <- match.arg(peak_to_gene_method)
    candidate_screen <- match.arg(candidate_screen)
    condition_weight <- match.arg(condition_weight)
    lambda_selection <- match.arg(lambda_selection)
    small_condition_action <- match.arg(small_condition_action)

    .condition_validate_public_args(
        object, cell_type_col, condition_col, scale, alpha, condition_mix,
        reference_condition, nlambda, lambda, nfolds, min_cells_per_condition,
        active_tol, max_iter, tol_objective, tol_coef, tf_cor, peak_cor,
        lambda_min_ratio, parallel, overwrite, seed, verbose
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
    .condition_fit_cell_type_models(
        object = object, prepared = prepared, cell_type_col = cell_type_col,
        cell_type = cell_type,
        condition_col = condition_col, network_name = network_name,
        peak_to_gene_method = peak_to_gene_method, upstream = upstream,
        downstream = downstream, extend = extend, only_tss = only_tss,
        tf_cor = tf_cor, peak_cor = peak_cor, candidate_screen = candidate_screen,
        alpha = alpha, condition_mix = condition_mix,
        reference_condition = reference_condition,
        condition_weight = condition_weight, nlambda = nlambda, lambda = lambda,
        lambda_min_ratio = lambda_min_ratio, nfolds = nfolds,
        lambda_selection = lambda_selection,
        min_cells_per_condition = min_cells_per_condition,
        small_condition_action = small_condition_action,
        active_tol = active_tol, parallel = parallel, BPPARAM = BPPARAM,
        overwrite = overwrite, seed = seed, max_iter = max_iter,
        tol_objective = tol_objective, tol_coef = tol_coef, verbose = verbose
    )
}

.condition_validate_public_args <- function(
    object,
    cell_type_col,
    condition_col,
    scale,
    alpha,
    condition_mix,
    reference_condition,
    nlambda,
    lambda,
    nfolds,
    min_cells_per_condition,
    active_tol,
    max_iter,
    tol_objective,
    tol_coef,
    tf_cor,
    peak_cor,
    lambda_min_ratio,
    parallel,
    overwrite,
    seed,
    verbose
) {
    if (!inherits(object, 'GRNData')) {
        stop('object must be a GRNData object.')
    }
    metadata <- object@data@meta.data
    .condition_validate_analysis_metadata(
        metadata, cell_type_col, condition_col
    )
    if (!isTRUE(scale)) {
        stop('Condition sub-GRN comparability requires scale = TRUE.')
    }
    if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
        alpha < 0 || alpha > 1) {
        stop('alpha must be between 0 and 1.')
    }
    if (!is.numeric(condition_mix) || length(condition_mix) != 1L ||
        !is.finite(condition_mix) || condition_mix <= 0 ||
        condition_mix > 1) {
        stop('condition_mix must be greater than zero and no larger than one.')
    }
    if (!is.null(reference_condition) &&
        (length(reference_condition) != 1L || is.na(reference_condition) ||
         !is.character(reference_condition) ||
         !nzchar(trimws(reference_condition)) ||
         reference_condition != trimws(reference_condition))) {
        stop(
            'reference_condition must be NULL or one exact non-empty ',
            'condition label without surrounding whitespace.'
        )
    }
    if (!is.numeric(nlambda) || length(nlambda) != 1L ||
        !is.finite(nlambda) ||
        abs(nlambda - round(nlambda)) > sqrt(.Machine$double.eps) ||
        nlambda > .Machine$integer.max ||
        nlambda < 1L) {
        stop('nlambda must be a positive integer.')
    }
    if (!is.null(lambda) &&
        (!is.numeric(lambda) || length(lambda) == 0L || any(!is.finite(lambda)) || any(lambda < 0))) {
        stop('lambda must contain finite non-negative values.')
    }
    if (is.null(lambda) && nlambda < 2L) {
        stop('nlambda must be at least 2 when a lambda path is generated.')
    }
    if (!is.numeric(nfolds) || length(nfolds) != 1L ||
        !is.finite(nfolds) ||
        abs(nfolds - round(nfolds)) > sqrt(.Machine$double.eps) ||
        nfolds > .Machine$integer.max ||
        nfolds < 2L) {
        stop('nfolds must be at least 2 for within-cell-type OOF.')
    }
    if (!is.numeric(min_cells_per_condition) ||
        length(min_cells_per_condition) != 1L ||
        !is.finite(min_cells_per_condition) ||
        abs(min_cells_per_condition - round(min_cells_per_condition)) >
            sqrt(.Machine$double.eps) ||
        min_cells_per_condition > .Machine$integer.max ||
        min_cells_per_condition < nfolds) {
        stop(
            'min_cells_per_condition must be an integer no smaller than ',
            'nfolds.'
        )
    }
    correlation_thresholds <- c(tf_cor = tf_cor, peak_cor = peak_cor)
    if (!is.numeric(correlation_thresholds) ||
        length(correlation_thresholds) != 2L ||
        any(!is.finite(correlation_thresholds)) ||
        any(correlation_thresholds < 0 | correlation_thresholds > 1)) {
        stop('tf_cor and peak_cor must be finite values between 0 and 1.')
    }
    if (!is.null(lambda_min_ratio) &&
        (!is.numeric(lambda_min_ratio) || length(lambda_min_ratio) != 1L ||
         !is.finite(lambda_min_ratio) || lambda_min_ratio <= 0 ||
         lambda_min_ratio >= 1)) {
        stop('lambda_min_ratio must be NULL or one finite value in (0, 1).')
    }
    if (!is.numeric(active_tol) || length(active_tol) != 1L ||
        active_tol < 0 || !is.finite(active_tol)) {
        stop('active_tol must be finite and non-negative.')
    }
    if (!is.numeric(max_iter) || length(max_iter) != 1L ||
        !is.finite(max_iter) ||
        abs(max_iter - round(max_iter)) > sqrt(.Machine$double.eps) ||
        max_iter > .Machine$integer.max ||
        max_iter < 1L ||
        !is.numeric(tol_objective) || length(tol_objective) != 1L ||
        !is.finite(tol_objective) || tol_objective <= 0 ||
        !is.numeric(tol_coef) || length(tol_coef) != 1L ||
        !is.finite(tol_coef) || tol_coef <= 0) {
        stop('Solver iteration and tolerance parameters are invalid.')
    }
    logical_controls <- list(
        parallel = parallel, overwrite = overwrite, verbose = verbose
    )
    invalid_logical <- vapply(logical_controls, function(value) {
        !is.logical(value) || length(value) != 1L || is.na(value)
    }, logical(1))
    if (any(invalid_logical)) {
        stop('parallel, overwrite, and verbose must be TRUE or FALSE.')
    }
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
        abs(seed - round(seed)) > sqrt(.Machine$double.eps) ||
        abs(seed) > .Machine$integer.max) {
        stop('seed must be one finite integer.')
    }
    invisible(TRUE)
}

.condition_validate_analysis_metadata <- function(
    metadata, cell_type_col, condition_col
) {
    column_arguments <- list(
        cell_type_col = cell_type_col,
        condition_col = condition_col
    )
    invalid_column_argument <- vapply(column_arguments, function(value) {
        !is.character(value) || length(value) != 1L || is.na(value) ||
            !nzchar(trimws(value))
    }, logical(1))
    if (any(invalid_column_argument)) {
        stop('cell_type_col and condition_col must be non-empty column names.')
    }
    if (identical(cell_type_col, condition_col)) {
        stop('cell_type_col and condition_col must be different columns.')
    }
    missing_columns <- setdiff(c(cell_type_col, condition_col), colnames(metadata))
    if (length(missing_columns) > 0L) {
        stop('Missing metadata column(s): ', paste(missing_columns, collapse = ', '), '.')
    }
    labels <- lapply(
        metadata[, c(cell_type_col, condition_col), drop = FALSE],
        as.character
    )
    invalid_labels <- vapply(labels, function(value) {
        anyNA(value) || any(!nzchar(trimws(value))) ||
            any(value != trimws(value))
    }, logical(1))
    if (any(invalid_labels)) {
        stop(
            'cell_type_col and condition_col labels must be complete, ',
            'non-empty, and free of surrounding whitespace.'
        )
    }
    invisible(TRUE)
}
