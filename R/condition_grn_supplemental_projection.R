# Explicit supplemental RegCompass handoff for condition-estimable OOF edges.

#' Project condition-specific supplemental OOF regulatory signal
#'
#' Uses each condition's own outer-fold estimability mask. Edges unique to one
#' condition can therefore contribute to a supplemental penalty route while the
#' pairwise/global common-support projection remains the primary comparison.
#'
#' @param object GRNData object containing the paired cells used for inference.
#' @param fit Optional canonical `pando_condition_grn_fit` object.
#' @param network_name,cell_type Optional fit filters when `fit` is omitted.
#' @param scale Standardized or raw coefficient scale.
#' @param output Return target scores or target scores plus edge contributions.
#' @param targets Optional target subset.
#' @param nonestimable Use structural zeros or stop on unavailable effects.
#' @param active_tol Activity threshold used in metadata.
#' @return A `ConditionGRNProjection` explicitly marked as supplemental rather
#'   than primary.
#' @export
project_condition_grn_supplemental_cells <- function(
    object,
    fit = NULL,
    network_name = NULL,
    cell_type = NULL,
    scale = c('std', 'raw'),
    output = c('gene_score', 'edge_contribution'),
    targets = NULL,
    nonestimable = c('structural_zero', 'error'),
    active_tol = if (is.null(fit)) 1e-8 else fit$active_tol
) {
    scale <- match.arg(scale)
    output <- match.arg(output)
    nonestimable <- match.arg(nonestimable)
    projection <- project_condition_grn_cells(
        object = object,
        fit = fit,
        network_name = network_name,
        cell_type = cell_type,
        component = 'condition',
        scale = scale,
        output = output,
        targets = targets,
        nonestimable = nonestimable,
        support_policy = 'condition_estimable',
        origin = 'oof',
        diagnostic_only = TRUE,
        active_tol = active_tol
    )
    projection$schema_version <-
        'pando_condition_grn_supplemental_projection_v1'
    projection$projection_used_for_penalty <- TRUE
    projection$full_fit_projection_used_for_penalty <- FALSE
    projection$projection_role <- 'supplemental_penalty'
    projection$score_comparability_class <-
        'condition_specific_supplement_on_shared_celltype_coordinate'
    projection$condition_unique_edges_allowed <- TRUE
    projection$common_support_primary <- TRUE
    projection$aggregation_contract$projection_role <-
        'supplemental_penalty'
    projection$aggregation_contract$condition_unique_edges_allowed <- TRUE
    projection$aggregation_contract$common_support_primary <- TRUE
    projection
}
