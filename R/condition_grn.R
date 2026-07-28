# Condition-aware extension of the Pando TF-peak-target model.

#' Infer condition-specific sub-GRNs on one shared coordinate system
#'
#' Extends Pando's TF-expression by motif-bearing peak-accessibility model with
#' one condition-aware workflow. Each cell type uses one candidate
#' TF-peak-target supergraph and one pooled single-cell effect scale. Sparse
#' multi-task selection permits different active edges and nodes in each
#' condition, including sign reversals. A support-constrained refit supplies the
#' comparable absolute effects used by downstream projection.
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
#' @param candidate_screen Candidate screening strategy. The default removes
#' condition means before pooled screening. `motif_domain` disables correlation
#' screening but retains the shared motif/domain candidate graph.
#' @param method Condition sub-GRN engine. Retired condition method labels are
#' accepted as call-compatibility aliases and route to the same condition-sparse
#' common-scale engine.
#' @param alpha Elastic-net mixing parameter.
#' @param condition_mix Group-versus-condition sparsity mixing parameter.
#' @param reference_condition Condition used for explicit coefficient contrasts.
#' Defaults to the first condition level within each cell type.
#' @param condition_weight Equal-condition or cell-count loss weighting.
#' @param nlambda Number of lambda values when lambda is not supplied.
#' @param lambda Optional fixed lambda value or decreasing lambda path.
#' @param lambda_min_ratio Smallest lambda relative to lambda maximum.
#' @param nfolds Number of condition-stratified cell folds.
#' @param lambda_selection Selection of lambda.1se or lambda.min.
#' @param min_cells_per_condition Minimum single cells per condition.
#' @param on_small_condition Action when a cell type contains a small condition.
#' @param scale Must be `TRUE`; all conditions use one pooled predictor and
#' response scale.
#' @param active_tol Numerical threshold used for network summaries.
#' @param parallel Use the existing Pando foreach mapping backend.
#' @param BPPARAM Optional BiocParallel parameter object.
#' @param overwrite Replace generated networks with matching names.
#' @param seed Random seed for condition-stratified folds.
#' @param max_iter,tol_objective,tol_coef Solver controls.
#' @param verbose Display progress messages.
#' @param ... Must be empty; retired condition modules do not accept additional
#' arguments.
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
    candidate_screen = c('pooled_within_condition', 'motif_domain'),
    method = c(
        'shared_baseline_condition_sparse',
        'shared_design_independent',
        'multitask_glmnet'
    ),
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
    on_small_condition = c('skip_cell_type', 'drop_condition', 'error'),
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
    method <- match.arg(method)
    condition_weight <- match.arg(condition_weight)
    lambda_selection <- match.arg(lambda_selection)
    on_small_condition <- match.arg(on_small_condition)
    if (!identical(method, 'shared_baseline_condition_sparse')) {
        log_message(
            'Treating method = "', method,
            '" as an alias of shared_baseline_condition_sparse.',
            verbose = verbose
        )
    }

    .condition_validate_public_args(
        object, cell_type_col, condition_col, method, scale, alpha, condition_mix,
        reference_condition, nlambda, lambda, nfolds, min_cells_per_condition,
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
        verbose = verbose
    )
    .condition_run_all_cell_types(
        object = object, prepared = prepared, cell_type_col = cell_type_col,
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
        on_small_condition = on_small_condition,
        active_tol = active_tol, parallel = parallel, BPPARAM = BPPARAM,
        overwrite = overwrite, seed = seed, max_iter = max_iter,
        tol_objective = tol_objective, tol_coef = tol_coef, verbose = verbose
    )
}

.condition_validate_public_args <- function(
    object,
    cell_type_col,
    condition_col,
    method,
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
    if (!method %in% c(
        'shared_baseline_condition_sparse',
        'shared_design_independent',
        'multitask_glmnet'
    )) {
        stop('Unsupported condition sub-GRN method label.')
    }
    if (!isTRUE(scale)) {
        stop('Condition sub-GRN comparability requires scale = TRUE.')
    }
    if (!is.numeric(alpha) || length(alpha) != 1L || alpha < 0 || alpha > 1) {
        stop('alpha must be between 0 and 1.')
    }
    if (!is.numeric(condition_mix) || length(condition_mix) != 1L ||
        condition_mix <= 0 || condition_mix > 1) {
        stop('condition_mix must be greater than zero and no larger than one.')
    }
    if (!is.null(reference_condition) &&
        (length(reference_condition) != 1L || is.na(reference_condition) ||
         !nzchar(trimws(as.character(reference_condition))))) {
        stop('reference_condition must be NULL or one non-empty condition label.')
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
