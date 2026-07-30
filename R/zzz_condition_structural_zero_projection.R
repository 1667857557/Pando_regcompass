# Structural-zero semantics for non-estimable projection edges.

.project_condition_grn_cells_na <- project_condition_grn_cells
project_condition_grn_cells <- function(
    object,
    fit = NULL,
    network_name = NULL,
    cell_type = NULL,
    component = c('condition', 'shared', 'deviation'),
    scale = c('std', 'raw'),
    output = c('gene_score', 'edge_contribution'),
    targets = NULL,
    nonestimable = c('structural_zero', 'error'),
    support_policy = c(
        'pairwise_common', 'global_common', 'condition_estimable', 'strict'
    ),
    comparison_conditions = NULL,
    origin = c('oof', 'full_fit'),
    diagnostic_only = FALSE,
    active_tol = if (is.null(fit)) 1e-8 else fit$active_tol
) {
    component <- match.arg(component)
    scale <- match.arg(scale)
    output <- match.arg(output)
    nonestimable <- match.arg(nonestimable)
    support_policy <- match.arg(support_policy)
    origin <- match.arg(origin)
    value <- .project_condition_grn_cells_na(
        object = object,
        fit = fit,
        network_name = network_name,
        cell_type = cell_type,
        component = component,
        scale = scale,
        output = output,
        targets = targets,
        nonestimable = if (nonestimable == 'error') 'error' else 'propagate',
        support_policy = support_policy,
        comparison_conditions = comparison_conditions,
        origin = origin,
        diagnostic_only = diagnostic_only,
        active_tol = active_tol
    )
    if (nonestimable == 'error') return(value)

    score_zero <- !is.finite(value$gene_score)
    value$gene_score[score_zero] <- 0
    value$gene_direction <- sign(value$gene_score)
    value$gene_structural_zero_mask <- score_zero

    if (!is.null(value$edge_contribution)) {
        edge_zero <- !is.finite(value$edge_contribution)
        value$edge_contribution[edge_zero] <- 0
        value$edge_structural_zero_mask <- edge_zero
    } else {
        value$edge_structural_zero_mask <- NULL
    }
    if (is.data.frame(value$edge_effects)) {
        value$edge_effects$structural_zero <-
            !value$edge_effects$estimable |
            !is.finite(value$edge_effects$estimate)
        value$edge_effects$estimate[
            value$edge_effects$structural_zero
        ] <- 0
        value$edge_effects$direction <- sign(value$edge_effects$estimate)
    }
    if (is.data.frame(value$target_condition_status)) {
        value$target_condition_status$nonestimable_edge_policy <-
            'structural_zero'
        value$target_condition_status$structural_zero_edges <-
            value$target_condition_status$n_candidate_edges -
            value$target_condition_status$n_projection_support_edges
        value$target_condition_status$score_available <- TRUE
    }
    value$schema_version <- 'pando_condition_grn_projection_v4'
    value$nonestimable_policy <- 'structural_zero'
    value$structural_zero_definition <- paste(
        'an edge absent from the requested estimability/support set contributes',
        'exactly zero before target summation and remains eligible for downstream',
        'metacell aggregation'
    )
    value$aggregation_contract$propagate_NA <- FALSE
    value$aggregation_contract$nonestimable_edge <- 'structural_zero'
    value$aggregation_contract$structural_zero_enters_downstream <- TRUE
    value$aggregation_contract$operation <-
        'arithmetic_mean_by_target_including_structural_zeros'
    value
}

.aggregate_condition_grn_projection_na <- aggregate_condition_grn_projection
aggregate_condition_grn_projection <- function(
    projection, membership, group_col = 'metacell_id'
) {
    if (!inherits(projection, 'ConditionGRNProjection')) {
        stop('projection must inherit from ConditionGRNProjection.')
    }
    membership$cell_id <- as.character(membership$cell_id)
    membership[[group_col]] <- as.character(membership[[group_col]])
    cells <- rownames(projection$gene_score)
    membership <- membership[match(cells, membership$cell_id), , drop = FALSE]
    if (anyNA(membership$cell_id)) {
        stop('membership does not cover every projected single cell.')
    }
    structural_zero_fraction <- NULL
    if (!is.null(projection$gene_structural_zero_mask)) {
        groups <- unique(membership[[group_col]])
        structural_zero_fraction <- do.call(rbind, lapply(groups, function(group) {
            rows <- membership[[group_col]] == group
            colMeans(
                projection$gene_structural_zero_mask[rows, , drop = FALSE]
            )
        }))
        rownames(structural_zero_fraction) <- groups
    }
    projection$gene_score[!is.finite(projection$gene_score)] <- 0
    value <- .aggregate_condition_grn_projection_na(
        projection, membership, group_col = group_col
    )
    value$schema_version <- 'pando_condition_grn_group_projection_v2'
    value$gene_score[!is.finite(value$gene_score)] <- 0
    value$gene_direction <- sign(value$gene_score)
    value$gene_structural_zero_fraction <- structural_zero_fraction
    value$nonestimable_policy <- 'structural_zero'
    value$structural_zero_enters_downstream <- TRUE
    value$aggregation_order <- paste(
        'single_cell_TF_times_ATAC_then_transform_then_project_with',
        'nonestimable_edges_zero_then_mean'
    )
    value
}
