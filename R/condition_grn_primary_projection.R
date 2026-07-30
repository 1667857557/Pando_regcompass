# Canonical RegCompass handoff that retains condition-specific estimable edges.

#' Project the primary condition-specific OOF regulatory signal
#'
#' This is the canonical penalty handoff for condition-aware Pando fits. It uses
#' each condition's own outer-fold estimability mask, so an edge may contribute
#' in one condition and be a structural zero in another. Candidate identity,
#' coefficient units, fold-local transforms, targets, and cell-type scope remain
#' shared; common-edge support is not required.
#'
#' @param object GRNData object containing the paired cells used for inference.
#' @param fit Optional canonical `pando_condition_grn_fit` object.
#' @param network_name,cell_type Optional fit filters when `fit` is omitted.
#' @param scale Standardized or raw target units.
#' @param output Return target scores or target scores plus edge contributions.
#' @param targets Optional target subset. Matching is case-insensitive when a fit
#'   is supplied, while returned target names retain the fit's original case.
#' @param nonestimable Use structural zeros or stop on unavailable effects.
#' @param active_tol Activity threshold used in metadata.
#' @return A `ConditionGRNProjection` whose condition-specific OOF scores are
#'   explicitly eligible for the primary RegCompass penalty.
#' @export
project_condition_grn_primary_cells <- function(
    object,
    fit = NULL,
    network_name = NULL,
    cell_type = NULL,
    scale = c("std", "raw"),
    output = c("gene_score", "edge_contribution"),
    targets = NULL,
    nonestimable = c("structural_zero", "error"),
    active_tol = if (is.null(fit)) 1e-8 else fit$active_tol
) {
    scale <- match.arg(scale)
    output <- match.arg(output)
    nonestimable <- match.arg(nonestimable)
    if (!is.null(fit) && !is.null(targets)) {
        fit_targets <- unique(as.character(fit$edge_table$target))
        index <- match(tolower(as.character(targets)), tolower(fit_targets))
        if (anyNA(index)) {
            stop(
                "Target(s) were not found in the fit: ",
                paste(as.character(targets)[is.na(index)], collapse = ", "),
                "."
            )
        }
        targets <- fit_targets[index]
    }
    projection <- project_condition_grn_cells(
        object = object,
        fit = fit,
        network_name = network_name,
        cell_type = cell_type,
        component = "condition",
        scale = scale,
        output = output,
        targets = targets,
        nonestimable = nonestimable,
        support_policy = "condition_estimable",
        origin = "oof",
        diagnostic_only = TRUE,
        active_tol = active_tol
    )
    projection$schema_version <- "pando_condition_grn_primary_projection_v1"
    projection$projection_used_for_penalty <- TRUE
    projection$full_fit_projection_used_for_penalty <- FALSE
    projection$projection_role <- "primary_penalty"
    projection$score_comparability_class <-
        "condition_specific_effects_on_shared_celltype_coordinate"
    projection$condition_unique_edges_allowed <- TRUE
    projection$common_support_required <- FALSE
    projection$primary_support_policy <- "condition_estimable"
    projection$aggregation_contract$projection_role <- "primary_penalty"
    projection$aggregation_contract$condition_unique_edges_allowed <- TRUE
    projection$aggregation_contract$common_support_required <- FALSE
    projection
}
