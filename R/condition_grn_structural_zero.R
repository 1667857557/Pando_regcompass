# Condition-full OOF projection with explicit projectable structural zeros.
#
# This file contains structural-zero helpers and the public primary-projection
# adapter only. Canonical lambda-path fitting, nested cross-fitting and refitting
# live in condition_grn_solver.R, condition_grn_crossfit.R and
# condition_grn_refit.R and must not be redefined here.

.condition_population_variance_base <- .condition_population_variance
.condition_column_variance_base <- .condition_column_variance
.condition_combine_fit_contracts_base <- .condition_combine_fit_contracts

.condition_exact_zero_columns <- function(x) {
    if (!ncol(x)) return(logical())
    as.numeric(Matrix::colSums(abs(x))) == 0
}

.condition_true_variance_mask <- function(X_list, coefficient_mask = NULL) {
    if (!is.list(X_list) || !length(X_list)) {
        stop('X_list must contain one predictor matrix per condition.')
    }
    p <- ncol(X_list[[1L]])
    k <- length(X_list)
    if (any(vapply(X_list, ncol, integer(1)) != p)) {
        stop('All condition matrices must share predictor columns.')
    }
    if (is.null(coefficient_mask)) {
        coefficient_mask <- matrix(
            TRUE, p, k,
            dimnames = list(colnames(X_list[[1L]]), names(X_list))
        )
    }
    coefficient_mask <- as.matrix(coefficient_mask)
    if (!is.logical(coefficient_mask) || anyNA(coefficient_mask) ||
        !identical(dim(coefficient_mask), c(p, k))) {
        stop('coefficient_mask must be a logical predictors-by-conditions matrix.')
    }
    variance_mask <- vapply(
        X_list,
        function(x) {
            value <- .condition_population_variance_base(x)
            is.finite(value) & value > .Machine$double.eps
        },
        logical(p)
    )
    if (is.null(dim(variance_mask))) {
        variance_mask <- matrix(variance_mask, nrow = p)
    }
    dimnames(variance_mask) <- dimnames(coefficient_mask)
    coefficient_mask & variance_mask
}

# Candidate-retention calls omit `center`; exact-zero columns therefore remain
# represented in the shared structural supergraph. Fold-statistics calls provide
# `center` explicitly and receive the true zero variance, so those columns are
# never coefficient-estimable.
.condition_population_variance <- function(x, center = NULL) {
    value <- .condition_population_variance_base(x, center = center)
    if (is.null(center)) {
        zero <- .condition_exact_zero_columns(x)
        value[zero] <- 4 * .Machine$double.eps
    }
    value
}

.condition_column_variance <- function(x) {
    value <- .condition_column_variance_base(x)
    zero <- .condition_exact_zero_columns(x)
    value[zero] <- 4 * .Machine$double.eps
    value
}

.condition_combine_fit_contracts <- function(...) {
    fit <- .condition_combine_fit_contracts_base(...)
    fit$contract_version <- 'condition_absolute_oof_v4'
    fit$coefficient_estimable_mask <- fit$estimability_mask
    fit$projectable_structural_zero_mask <-
        fit$topology_mask & !fit$estimability_mask
    fit$projection_support_mask <-
        fit$coefficient_estimable_mask |
        fit$projectable_structural_zero_mask
    fit$projection_used_for_penalty <- TRUE
    fit$primary_projection <- 'projection_condition_full_oof'
    fit$common_projection_role <- 'shared_estimable_component'
    fit$condition_unique_projection_role <-
        'projection_condition_full_oof - projection_common_oof'
    fit$nonestimable_projection_policy <- 'structural_zero_by_condition'
    fit$projection_contract$score <- paste(
        'outer-heldout condition-full sum(z_edge * beta_condition)',
        'with nonestimable edge contribution fixed at zero'
    )
    fit$projection_contract$primary_support_policy <- 'condition_full_oof'
    fit$projection_contract$common_support_role <-
        'shared_estimable_component'
    fit$projection_contract$condition_unique_role <-
        'condition_full_oof_minus_common_support_oof'
    fit$projection_contract$nonestimable <-
        'projectable_structural_zero_by_condition'
    fit$projection_contract$condition_full_role <- 'primary_penalty'
    fit
}

#' Project the primary condition-full OOF regulatory signal
#'
#' The primary score contains every edge estimable in the focal condition.
#' Jointly estimable edges form the common-support component. An edge that is
#' non-estimable in one or both conditions remains in the shared candidate
#' supergraph and contributes exactly zero in each non-estimable condition.
#'
#' @param object GRNData object containing the paired cells used for inference.
#' @param fit Optional canonical `pando_condition_grn_fit` object.
#' @param network_name,cell_type Optional fit filters when `fit` is omitted.
#' @param scale Standardized or raw coefficient scale.
#' @param targets Optional target subset.
#' @param nonestimable Use structural zeros or stop on unavailable target scores.
#' @param active_tol Activity threshold used in metadata.
#' @return A primary `ConditionGRNProjection`.
#' @export
project_condition_grn_primary_cells <- function(
    object,
    fit = NULL,
    network_name = NULL,
    cell_type = NULL,
    scale = c('std', 'raw'),
    targets = NULL,
    nonestimable = c('structural_zero', 'error'),
    active_tol = if (is.null(fit)) 1e-8 else fit$active_tol
) {
    scale <- match.arg(scale)
    nonestimable <- match.arg(nonestimable)
    projection <- project_condition_grn_cells(
        object = object,
        fit = fit,
        network_name = network_name,
        cell_type = cell_type,
        component = 'condition',
        scale = scale,
        output = 'gene_score',
        targets = targets,
        nonestimable = nonestimable,
        support_policy = 'condition_estimable',
        origin = 'oof',
        diagnostic_only = TRUE,
        active_tol = active_tol
    )
    projection$schema_version <-
        'pando_condition_grn_primary_projection_v1'
    projection$projection_used_for_penalty <- TRUE
    projection$projection_role <- 'primary_penalty'
    projection$score_comparability_class <-
        'condition_full_oof_on_shared_celltype_coordinate'
    projection$primary_support_policy <- 'condition_full_oof'
    projection$common_support_role <- 'shared_estimable_component'
    projection$condition_unique_role <-
        'condition_full_oof_minus_common_support_oof'
    projection$nonestimable_policy <- 'structural_zero_by_condition'
    projection$aggregation_contract$projection_role <- 'primary_penalty'
    projection$aggregation_contract$primary_support_policy <-
        'condition_full_oof'
    projection
}
