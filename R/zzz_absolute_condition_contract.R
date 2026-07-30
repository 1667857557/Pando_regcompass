# Absolute condition-effect contract without a reference-condition contrast.

.condition_reference_fields <- c(
    'reference_condition', 'reference_beta', 'contrast', 'comparison_mask',
    'contrast_formula', 'pairwise_contrast_formula', 'comparison_mask_formula'
)

.condition_drop_reference_columns <- function(x) {
    if (!is.data.frame(x)) return(x)
    x[, setdiff(colnames(x), .condition_reference_fields), drop = FALSE]
}

.condition_absolute_fit_contract <- function(fit) {
    if (!is.list(fit)) return(fit)
    fit[intersect(names(fit), .condition_reference_fields)] <- NULL
    if (is.data.frame(fit$response_transform)) {
        fit$response_transform <- .condition_drop_reference_columns(
            fit$response_transform
        )
    }
    if (is.list(fit$direction_semantics)) {
        fit$direction_semantics$pairwise <- NULL
    }
    fit$contract_version <- 'condition_absolute_oof_v2'
    fit$coefficient_contract <- 'absolute_condition_effects_only'
    fit$condition_comparison_policy <- paste(
        'compare beta_condition on the shared equal-condition coordinate;',
        'no reference-condition coefficient or stored contrast is defined'
    )
    class(fit) <- unique(c('ConditionGRNFit', class(fit)))
    fit
}

.condition_absolute_object_contract <- function(object) {
    fits <- object@grn@params$condition_grn_fits
    if (is.list(fits) && length(fits)) {
        object@grn@params$condition_grn_fits <- lapply(
            fits, .condition_absolute_fit_contract
        )
    }
    index <- object@grn@params$condition_network_index
    if (is.data.frame(index)) {
        object@grn@params$condition_network_index <-
            .condition_drop_reference_columns(index)
    }
    object@grn@params[intersect(
        names(object@grn@params), .condition_reference_fields
    )] <- NULL
    object@grn@params$condition_coefficient_contract <-
        'absolute_condition_effects_only'
    object
}

.condition_combine_fit_contracts_with_reference <-
    .condition_combine_fit_contracts
.condition_combine_fit_contracts <- function(...) {
    .condition_absolute_fit_contract(
        .condition_combine_fit_contracts_with_reference(...)
    )
}

.condition_network_params_with_reference <- .condition_network_params
.condition_network_params <- function(...) {
    value <- .condition_network_params_with_reference(...)
    value[intersect(names(value), .condition_reference_fields)] <- NULL
    value$coefficient_contract <- 'absolute_condition_effects_only'
    value
}

.condition_index_row_with_reference <- .condition_index_row
.condition_index_row <- function(...) {
    .condition_drop_reference_columns(
        .condition_index_row_with_reference(...)
    )
}

.infer_condition_grn_with_reference <- infer_condition_grn.GRNData
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
    verbose = TRUE,
    ...
) {
    value <- .infer_condition_grn_with_reference(
        object = object,
        cell_type_col = cell_type_col,
        condition_col = condition_col,
        cell_type = cell_type,
        genes = genes,
        network_name = network_name,
        peak_to_gene_method = peak_to_gene_method,
        upstream = upstream,
        downstream = downstream,
        extend = extend,
        only_tss = only_tss,
        peak_to_gene_domains = peak_to_gene_domains,
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
        scale = scale,
        active_tol = active_tol,
        parallel = parallel,
        BPPARAM = BPPARAM,
        overwrite = overwrite,
        seed = seed,
        max_iter = max_iter,
        tol_objective = tol_objective,
        tol_coef = tol_coef,
        verbose = verbose,
        ...
    )
    .condition_absolute_object_contract(value)
}

.condition_grn_fit_with_reference <- condition_grn_fit.GRNData
condition_grn_fit.GRNData <- function(
    object, network_name = NULL, cell_type = NULL
) {
    value <- .condition_grn_fit_with_reference(
        object, network_name = network_name, cell_type = cell_type
    )
    if (inherits(value, 'ConditionGRNFit')) {
        return(.condition_absolute_fit_contract(value))
    }
    lapply(value, .condition_absolute_fit_contract)
}

if (exists('condition_grn_contrast', inherits = FALSE)) {
    rm(condition_grn_contrast)
}
