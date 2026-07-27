# Condition-aware extension of the Pando TF-peak-target model.

#' Infer cell-type Universal and condition-specific gene regulatory networks
#'
#' Extends Pando's TF-expression by motif-bearing peak-accessibility interaction
#' model to a sparse-group multi-task model. Cell types are fit independently
#' across conditions. All generated networks remain standard Pando Network
#' objects inside the returned GRNData object.
#'
#' @param object A GRNData object containing paired RNA and ATAC measurements.
#' @param cell_type_col Metadata column containing cell-type labels.
#' @param condition_col Metadata column containing condition labels.
#' @param genes Target genes. Defaults to RNA variable features.
#' @param network_name Prefix used for generated network names.
#' @param peak_to_gene_method Peak-to-gene method, either Signac or GREAT.
#' @param upstream,downstream,extend,only_tss Regulatory-domain parameters.
#' @param peak_to_gene_domains Optional custom gene regulatory domains.
#' @param tf_cor,peak_cor Correlation screening thresholds.
#' @param candidate_screen Candidate screening strategy.
#' @param aggregate_rna_col,aggregate_peaks_col Optional shared aggregation label.
#' @param method Must be multitask_glmnet.
#' @param family Must be gaussian.
#' @param interaction_term Must be colon, matching default Pando behavior.
#' @param alpha Elastic-net mixing parameter.
#' @param condition_mix Group-versus-condition sparsity mixing parameter.
#' @param condition_weight Equal-condition or cell-count loss weighting.
#' @param nlambda Number of lambda values when lambda is not supplied.
#' @param lambda Optional fixed lambda value or decreasing lambda path.
#' @param lambda_min_ratio Smallest lambda relative to lambda maximum.
#' @param nfolds Number of condition-stratified cell folds.
#' @param lambda_selection Selection of lambda.1se or lambda.min.
#' @param min_cells_per_condition Minimum cells or aggregate units per condition.
#' @param on_small_condition Action when a cell type contains a small condition.
#' @param scale Whether to reproduce Pando scaling before interaction creation.
#' @param active_tol Numerical threshold used for network summaries.
#' @param parallel Use the existing Pando foreach mapping backend.
#' @param BPPARAM Optional BiocParallel parameter object.
#' @param overwrite Replace generated networks with matching names.
#' @param seed Random seed for condition-stratified folds.
#' @param max_iter,tol_objective,tol_coef Solver controls.
#' @param verbose Display progress messages.
#' @param ... Reserved for future extensions.
#' @return A GRNData object containing standard Pando Network objects.
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
    candidate_screen = c('condition_union', 'pooled', 'motif_domain'),
    aggregate_rna_col = NULL,
    aggregate_peaks_col = NULL,
    method = 'multitask_glmnet',
    family = 'gaussian',
    interaction_term = ':',
    alpha = 0.5,
    condition_mix = 0.5,
    condition_weight = c('equal', 'cell_count'),
    nlambda = 50L,
    lambda = NULL,
    lambda_min_ratio = NULL,
    nfolds = 5L,
    lambda_selection = c('lambda.1se', 'lambda.min'),
    min_cells_per_condition = 50L,
    on_small_condition = c('skip_cell_type', 'drop_condition', 'error'),
    scale = FALSE,
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
    peak_to_gene_method <- match.arg(peak_to_gene_method)
    candidate_screen <- match.arg(candidate_screen)
    condition_weight <- match.arg(condition_weight)
    lambda_selection <- match.arg(lambda_selection)
    on_small_condition <- match.arg(on_small_condition)

    .condition_validate_public_args(
        object, cell_type_col, condition_col, aggregate_rna_col,
        aggregate_peaks_col, method, family, interaction_term, alpha,
        condition_mix, nlambda, lambda, nfolds, min_cells_per_condition,
        active_tol, max_iter, tol_objective, tol_coef
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
        aggregate_col = aggregate_rna_col,
        verbose = verbose
    )
    .condition_run_all_cell_types(
        object = object, prepared = prepared, cell_type_col = cell_type_col,
        condition_col = condition_col, network_name = network_name,
        peak_to_gene_method = peak_to_gene_method, upstream = upstream,
        downstream = downstream, extend = extend, only_tss = only_tss,
        tf_cor = tf_cor, peak_cor = peak_cor, candidate_screen = candidate_screen,
        alpha = alpha, condition_mix = condition_mix,
        condition_weight = condition_weight, nlambda = nlambda, lambda = lambda,
        lambda_min_ratio = lambda_min_ratio, nfolds = nfolds,
        lambda_selection = lambda_selection,
        min_cells_per_condition = min_cells_per_condition,
        on_small_condition = on_small_condition, scale = scale,
        active_tol = active_tol, parallel = parallel, BPPARAM = BPPARAM,
        overwrite = overwrite, seed = seed, max_iter = max_iter,
        tol_objective = tol_objective, tol_coef = tol_coef, verbose = verbose
    )
}

.condition_validate_public_args <- function(
    object,
    cell_type_col,
    condition_col,
    aggregate_rna_col,
    aggregate_peaks_col,
    method,
    family,
    interaction_term,
    alpha,
    condition_mix,
    nlambda,
    lambda,
    nfolds,
    min_cells_per_condition,
    active_tol,
    max_iter,
    tol_objective,
    tol_coef
) {
    if (!inherits(object, 'GRNData')) {
        stop('object must be a GRNData object.')
    }
    metadata <- object@data@meta.data
    missing_columns <- setdiff(c(cell_type_col, condition_col), colnames(metadata))
    if (length(missing_columns) > 0L) {
        stop('Missing metadata column(s): ', paste(missing_columns, collapse = ', '), '.')
    }
    if (anyNA(metadata[[cell_type_col]]) || anyNA(metadata[[condition_col]])) {
        stop('cell_type_col and condition_col cannot contain NA values.')
    }
    if (identical(cell_type_col, condition_col)) {
        stop('cell_type_col and condition_col must be different columns.')
    }
    if (xor(is.null(aggregate_rna_col), is.null(aggregate_peaks_col)) ||
        (!is.null(aggregate_rna_col) && !identical(aggregate_rna_col, aggregate_peaks_col))) {
        stop('aggregate_rna_col and aggregate_peaks_col must both be NULL or identical.')
    }
    if (!is.null(aggregate_rna_col) && !aggregate_rna_col %in% colnames(metadata)) {
        stop('Aggregation column was not found in metadata.')
    }
    if (!identical(method, 'multitask_glmnet')) {
        stop("condition-aware inference currently requires method = 'multitask_glmnet'.")
    }
    if (!identical(family, 'gaussian')) {
        stop("condition-aware inference currently requires family = 'gaussian'.")
    }
    if (!identical(interaction_term, ':')) {
        stop("condition-aware inference currently requires interaction_term = ':'.")
    }
    if (!is.numeric(alpha) || length(alpha) != 1L || alpha < 0 || alpha > 1) {
        stop('alpha must be between 0 and 1.')
    }
    if (!is.numeric(condition_mix) || length(condition_mix) != 1L ||
        condition_mix < 0 || condition_mix > 1) {
        stop('condition_mix must be between 0 and 1.')
    }
    if (!is.numeric(nlambda) || length(nlambda) != 1L || nlambda < 1L) {
        stop('nlambda must be a positive integer.')
    }
    if (!is.null(lambda) &&
        (!is.numeric(lambda) || length(lambda) == 0L || any(!is.finite(lambda)) || any(lambda < 0))) {
        stop('lambda must contain finite non-negative values.')
    }
    cv_required <- is.null(lambda) || length(unique(as.numeric(lambda))) > 1L
    if (cv_required && nlambda < 2L) {
        stop('nlambda must be at least 2 when a lambda path is generated.')
    }
    if (!is.numeric(nfolds) || length(nfolds) != 1L || nfolds < 1L) {
        stop('nfolds must be a positive integer.')
    }
    if (cv_required && nfolds < 2L) {
        stop('nfolds must be at least 2 when cross-validation is required.')
    }
    minimum_required <- if (cv_required) nfolds else 1L
    if (!is.numeric(min_cells_per_condition) || min_cells_per_condition < minimum_required) {
        stop('min_cells_per_condition is too small for the requested validation scheme.')
    }
    if (!is.numeric(active_tol) || active_tol < 0 || !is.finite(active_tol)) {
        stop('active_tol must be finite and non-negative.')
    }
    if (!is.numeric(max_iter) || max_iter < 1L || tol_objective <= 0 || tol_coef <= 0) {
        stop('Solver iteration and tolerance parameters are invalid.')
    }
    invisible(TRUE)
}
