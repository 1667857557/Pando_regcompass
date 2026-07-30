# Fixed-transform single-cell projection for downstream metacell aggregation.

.condition_projection_coefficients <- function(fit, component, scale) {
    if (component == 'condition') {
        return(if (scale == 'std') fit$beta_condition_std else fit$beta_condition_raw)
    }
    if (component == 'deviation') {
        return(if (scale == 'std') fit$delta_condition_std else fit$delta_condition_raw)
    }
    shared <- if (scale == 'std') fit$beta_shared_std else fit$beta_shared_raw
    matrix(
        shared,
        nrow = length(shared),
        ncol = length(fit$condition_levels),
        dimnames = list(names(shared), fit$condition_levels)
    )
}

.condition_projection_predictor <- function(
    raw_predictor, center, predictor_scale, scale = c('std', 'raw')
) {
    scale <- match.arg(scale)
    predictor <- as.numeric(raw_predictor) - center
    if (scale == 'std') {
        if (!is.numeric(predictor_scale) || length(predictor_scale) != 1L ||
            !is.finite(predictor_scale) ||
            predictor_scale <= .Machine$double.eps) {
            stop('Stored predictor scale must be one finite positive value.')
        }
        predictor <- predictor / predictor_scale
    }
    predictor
}

.condition_projection_support <- function(
    fit, edge_index, support_policy, comparison_conditions
) {
    estimability <- fit$estimability_mask[edge_index, , drop = FALSE]
    projected_conditions <- fit$condition_levels
    support <- estimability
    if (support_policy == 'global_common') {
        common <- rowSums(estimability) == ncol(estimability)
        support[,] <- common
    } else if (support_policy == 'pairwise_common') {
        if (is.null(comparison_conditions)) {
            if (length(fit$condition_levels) != 2L) {
                stop('comparison_conditions must identify two fitted conditions.')
            }
            comparison_conditions <- fit$condition_levels
        }
        comparison_conditions <- unique(as.character(comparison_conditions))
        if (length(comparison_conditions) != 2L ||
            !all(comparison_conditions %in% fit$condition_levels)) {
            stop('comparison_conditions must identify exactly two fitted conditions.')
        }
        projected_conditions <- comparison_conditions
        common <- rowSums(
            estimability[, projected_conditions, drop = FALSE]
        ) == length(projected_conditions)
        support[,] <- FALSE
        support[, projected_conditions] <- common
    }
    list(
        support = support,
        estimability = estimability,
        projected_conditions = projected_conditions,
        comparison_conditions = comparison_conditions
    )
}

.condition_projection_edge_effects <- function(
    fit, coefficient, edge_index, component, active_tol
) {
    edge <- fit$edge_table[edge_index, , drop = FALSE]
    coefficient <- coefficient[edge_index, , drop = FALSE]
    activity <- .condition_projection_coefficients(
        fit, component = component, scale = 'std'
    )[edge_index, , drop = FALSE]
    do.call(rbind, lapply(fit$condition_levels, function(condition) {
        estimate <- coefficient[, condition]
        estimable <- if (component == 'shared') {
            rep(TRUE, nrow(edge))
        } else {
            fit$estimability_mask[edge_index, condition]
        }
        structural_zero <- !estimable | !is.finite(estimate)
        estimate[structural_zero] <- 0
        data.frame(
            edge,
            condition = condition,
            component = component,
            estimate = estimate,
            direction = sign(estimate),
            estimable = estimable,
            structural_zero = structural_zero,
            active = estimable & abs(activity[, condition]) > active_tol,
            stringsAsFactors = FALSE,
            check.names = FALSE
        )
    }))
}

.condition_projection_status <- function(
    fit, edges, edge_index, support, projected_conditions
) {
    do.call(rbind, lapply(unique(edges$target), function(target) {
        local <- which(edges$target == target)
        do.call(rbind, lapply(projected_conditions, function(condition) {
            candidate <- length(local)
            n_support <- sum(support[local, condition])
            n_estimable <- sum(fit$estimability_mask[edge_index[local], condition])
            data.frame(
                target = target,
                condition = condition,
                n_candidate_edges = candidate,
                n_estimable_edges = n_estimable,
                n_projection_support_edges = n_support,
                structural_zero_edges = candidate - n_support,
                estimable_fraction = if (candidate) n_estimable / candidate else NA_real_,
                score_completeness = if (candidate) n_support / candidate else NA_real_,
                score_available = TRUE,
                nonestimable_edge_policy = 'structural_zero',
                stringsAsFactors = FALSE
            )
        }))
    }))
}

#' Project a condition-specific GRN to paired single cells
#'
#' Non-estimable or unsupported edges contribute exactly zero before target
#' summation. Only stored outer-heldout common-support projections are eligible
#' for the primary downstream penalty.
#'
#' @param object GRNData object containing the paired cells used for inference.
#' @param fit Optional canonical `pando_condition_grn_fit` object.
#' @param network_name,cell_type Optional fit filters.
#' @param component Absolute condition, shared, or deviation coefficients.
#' @param scale Standardized or raw coefficient scale.
#' @param output Return target scores or target scores plus edge contributions.
#' @param targets Optional target subset.
#' @param nonestimable Use structural zeros or stop on unavailable effects.
#' @param support_policy Pairwise/global common support or diagnostic support.
#' @param comparison_conditions Pair used by pairwise common support.
#' @param origin Stored outer OOF projection or full-fit interpretation.
#' @param diagnostic_only Required for non-common support policies.
#' @param active_tol Activity threshold used in metadata.
#' @return A ConditionGRNProjection.
#' @export
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
    if (support_policy %in% c('condition_estimable', 'strict') &&
        !isTRUE(diagnostic_only)) {
        stop('condition_estimable and strict are diagnostic-only support policies.')
    }
    if (is.null(fit)) {
        fit <- condition_grn_fit(
            object, network_name = network_name, cell_type = cell_type
        )
    }
    .condition_require_fit(fit)
    if (!inherits(object, 'GRNData')) stop('object must be a GRNData object.')
    targets <- if (is.null(targets)) {
        unique(fit$edge_table$target)
    } else {
        unique(as.character(targets))
    }
    missing_targets <- setdiff(targets, fit$edge_table$target)
    if (length(missing_targets)) {
        stop('Target(s) were not found in the fit: ', paste(missing_targets, collapse = ', '), '.')
    }
    edges <- fit$edge_table[fit$edge_table$target %in% targets, , drop = FALSE]
    edge_index <- match(edges$edge_id, fit$edge_table$edge_id)
    support_info <- .condition_projection_support(
        fit, edge_index, support_policy, comparison_conditions
    )
    support <- support_info$support
    projected_conditions <- support_info$projected_conditions
    comparison_conditions <- support_info$comparison_conditions
    coefficient <- .condition_projection_coefficients(fit, component, scale)

    metadata <- object@data@meta.data
    cells <- fit$cell_ids
    if (is.null(cells) || !length(cells) ||
        !all(cells %in% rownames(metadata))) {
        stop('The object no longer contains every paired cell used by the fit.')
    }
    condition <- as.character(metadata[cells, fit$condition_col])
    keep <- condition %in% projected_conditions
    cells <- cells[keep]
    condition <- condition[keep]
    cell_metadata <- data.frame(
        cell_id = cells,
        cell_type = fit$cell_type,
        condition = condition,
        stringsAsFactors = FALSE,
        row.names = cells
    )
    status <- .condition_projection_status(
        fit, edges, edge_index, support, projected_conditions
    )

    if (origin == 'oof') {
        if (component != 'condition') {
            stop('OOF penalty projections are defined only for component = "condition".')
        }
        source <- if (support_policy == 'global_common') {
            fit$projection_global_common_oof
        } else if (support_policy == 'pairwise_common') {
            fit$projection_common_oof
        } else {
            fit$projection_condition_full_oof
        }
        if (is.null(source) || !all(cells %in% rownames(source)) ||
            !all(targets %in% colnames(source))) {
            stop('Stored OOF projection is incomplete for the requested cells or targets.')
        }
        score <- source[cells, targets, drop = FALSE]
        structural_zero <- !is.finite(score)
        if (nonestimable == 'error' && any(structural_zero)) {
            stop('Stored OOF projection contains unavailable target scores.')
        }
        score[structural_zero] <- 0
        if (scale == 'raw') {
            response_scale <- fit$response_transform$scale[
                match(targets, fit$response_transform$target)
            ]
            score <- sweep(score, 2L, response_scale, '*')
        }
        primary <- support_policy %in% c('pairwise_common', 'global_common') &&
            !isTRUE(diagnostic_only)
        return(structure(
            list(
                schema_version = 'pando_condition_grn_projection',
                network_name = fit$network_name,
                cell_type = fit$cell_type,
                component = component,
                scale = scale,
                support_policy = support_policy,
                comparison_conditions = comparison_conditions,
                projection_origin = 'outer_condition_stratified_cell_oof',
                projection_used_for_penalty = primary,
                full_fit_projection_used_for_penalty = FALSE,
                projection_role = if (primary) 'primary_penalty' else 'exploratory_only',
                score_comparability_class = if (primary) {
                    'primary_common_support_comparable'
                } else {
                    'exploratory_condition_full_not_strictly_comparable'
                },
                cell_metadata = cell_metadata,
                gene_score = score,
                gene_direction = sign(score),
                gene_structural_zero_mask = structural_zero,
                edge_contribution = NULL,
                edge_structural_zero_mask = NULL,
                edge_effects = .condition_projection_edge_effects(
                    fit, coefficient, edge_index, component, active_tol
                ),
                target_condition_status = status,
                nonestimable_policy = 'structural_zero',
                structural_zero_definition =
                    'unsupported or unavailable edge contributions equal zero before target summation',
                aggregation_contract = list(
                    group_within = c(fit$cell_type_col, fit$condition_col),
                    operation = 'arithmetic_mean_by_target_including_structural_zeros',
                    signed_scores = TRUE,
                    identical_target_columns = TRUE,
                    propagate_NA = FALSE,
                    nonestimable_edge = 'structural_zero',
                    structural_zero_enters_downstream = TRUE,
                    recompute_TF_peak_product = FALSE,
                    recompute_center_or_scale = FALSE,
                    refit_coefficients = FALSE,
                    projection_origin = 'outer_condition_stratified_cell_oof',
                    projection_role = if (primary) 'primary_penalty' else 'exploratory_only'
                )
            ),
            class = c('ConditionGRNProjection', 'list')
        ))
    }

    params <- Params(object)
    gene_data <- Matrix::t(LayerData(
        object, assay = params$rna_assay, layer = fit$assay_contract$rna_layer
    ))
    peak_data <- Matrix::t(LayerData(
        object, assay = params$peak_assay, layer = fit$assay_contract$peak_layer
    ))
    gene_data <- gene_data[cells, , drop = FALSE]
    peak_data <- peak_data[cells, , drop = FALSE]
    transform <- fit$predictor_transform[
        match(edges$edge_id, fit$predictor_transform$edge_id), , drop = FALSE
    ]
    if (length(setdiff(unique(edges$tf), colnames(gene_data))) ||
        length(setdiff(unique(edges$region), colnames(peak_data)))) {
        stop('The object no longer contains every TF and peak required by the fit.')
    }
    score <- matrix(
        0,
        nrow = length(cells),
        ncol = length(targets),
        dimnames = list(cells, targets)
    )
    edge_contribution <- if (output == 'edge_contribution') {
        matrix(
            0,
            nrow = length(cells),
            ncol = nrow(edges),
            dimnames = list(cells, edges$edge_id)
        )
    } else NULL
    edge_zero <- matrix(
        FALSE,
        nrow = length(cells),
        ncol = nrow(edges),
        dimnames = list(cells, edges$edge_id)
    )
    for (edge_pos in seq_len(nrow(edges))) {
        predictor <- .condition_projection_predictor(
            as.numeric(gene_data[, edges$tf[[edge_pos]]]) *
                as.numeric(peak_data[, edges$region[[edge_pos]]]),
            center = transform$center[[edge_pos]],
            predictor_scale = transform$scale[[edge_pos]],
            scale = scale
        )
        target_pos <- match(edges$target[[edge_pos]], targets)
        for (condition_name in projected_conditions) {
            cell_index <- condition == condition_name
            included <- support[edge_pos, condition_name]
            effect <- coefficient[edge_index[[edge_pos]], condition_name]
            unavailable <- !included || !is.finite(effect)
            if (unavailable) {
                edge_zero[cell_index, edge_pos] <- TRUE
                if (nonestimable == 'error') {
                    stop('A requested edge effect is unavailable for full-fit projection.')
                }
                next
            }
            contribution <- predictor[cell_index] * effect
            score[cell_index, target_pos] <- score[cell_index, target_pos] + contribution
            if (!is.null(edge_contribution)) {
                edge_contribution[cell_index, edge_pos] <- contribution
            }
        }
    }
    structure(
        list(
            schema_version = 'pando_condition_grn_projection',
            network_name = fit$network_name,
            cell_type = fit$cell_type,
            component = component,
            scale = scale,
            support_policy = support_policy,
            comparison_conditions = comparison_conditions,
            projection_origin = 'full_data_fit_interpretation_only',
            projection_used_for_penalty = FALSE,
            full_fit_projection_used_for_penalty = FALSE,
            projection_role = 'interpretation_only',
            cell_metadata = cell_metadata,
            gene_score = score,
            gene_direction = sign(score),
            gene_structural_zero_mask = matrix(
                FALSE, nrow(score), ncol(score), dimnames = dimnames(score)
            ),
            edge_contribution = edge_contribution,
            edge_structural_zero_mask = edge_zero,
            edge_effects = .condition_projection_edge_effects(
                fit, coefficient, edge_index, component, active_tol
            ),
            target_condition_status = status,
            nonestimable_policy = 'structural_zero',
            aggregation_contract = list(
                group_within = c(fit$cell_type_col, fit$condition_col),
                operation = 'arithmetic_mean_by_target_including_structural_zeros',
                signed_scores = TRUE,
                identical_target_columns = TRUE,
                propagate_NA = FALSE,
                nonestimable_edge = 'structural_zero',
                structural_zero_enters_downstream = TRUE,
                recompute_TF_peak_product = FALSE,
                recompute_center_or_scale = FALSE,
                refit_coefficients = FALSE
            )
        ),
        class = c('ConditionGRNProjection', 'list')
    )
}

#' Aggregate a cell-first condition GRN projection to metacells or groups
#'
#' @param projection A ConditionGRNProjection.
#' @param membership One row per cell with `cell_id` and `group_col`.
#' @param group_col Membership column identifying the output group.
#' @return A ConditionGRNGroupProjection.
#' @export
aggregate_condition_grn_projection <- function(
    projection, membership, group_col = 'metacell_id'
) {
    if (!inherits(projection, 'ConditionGRNProjection')) {
        stop('projection must inherit from ConditionGRNProjection.')
    }
    if (!is.data.frame(membership) ||
        !all(c('cell_id', group_col) %in% colnames(membership))) {
        stop('membership must contain cell_id and the requested group_col.')
    }
    membership$cell_id <- as.character(membership$cell_id)
    membership[[group_col]] <- as.character(membership[[group_col]])
    if (anyNA(membership$cell_id) || anyNA(membership[[group_col]]) ||
        any(!nzchar(trimws(membership$cell_id))) ||
        any(!nzchar(trimws(membership[[group_col]]))) ||
        anyDuplicated(membership$cell_id)) {
        stop('membership must map each cell exactly once to one non-empty group.')
    }
    cells <- rownames(projection$gene_score)
    membership <- membership[match(cells, membership$cell_id), , drop = FALSE]
    if (anyNA(membership$cell_id)) {
        stop('membership does not cover every projected single cell.')
    }
    groups <- unique(membership[[group_col]])
    score <- do.call(rbind, lapply(groups, function(group) {
        rows <- membership[[group_col]] == group
        colMeans(projection$gene_score[rows, , drop = FALSE])
    }))
    rownames(score) <- groups
    zero_fraction <- NULL
    if (!is.null(projection$gene_structural_zero_mask)) {
        zero_fraction <- do.call(rbind, lapply(groups, function(group) {
            rows <- membership[[group_col]] == group
            colMeans(projection$gene_structural_zero_mask[rows, , drop = FALSE])
        }))
        rownames(zero_fraction) <- groups
    }
    group_metadata <- do.call(rbind, lapply(groups, function(group) {
        rows <- membership[[group_col]] == group
        observed <- projection$cell_metadata[cells[rows], , drop = FALSE]
        values <- lapply(c('cell_type', 'condition'), function(field) {
            value <- unique(as.character(observed[[field]]))
            value <- value[!is.na(value) & nzchar(value)]
            if (length(value) != 1L) stop('Group ', group, ' mixes ', field, ' values.')
            value[[1L]]
        })
        data.frame(
            group_id = group,
            cell_type = values[[1L]],
            condition = values[[2L]],
            n_cells = sum(rows),
            stringsAsFactors = FALSE
        )
    }))
    rownames(group_metadata) <- groups
    structure(
        list(
            schema_version = 'pando_condition_grn_group_projection',
            source_projection = projection,
            projection_origin = projection$projection_origin,
            projection_used_for_penalty =
                isTRUE(projection$projection_used_for_penalty),
            full_fit_projection_used_for_penalty = FALSE,
            gene_score = score,
            gene_direction = sign(score),
            gene_structural_zero_fraction = zero_fraction,
            group_metadata = group_metadata,
            group_col = group_col,
            nonestimable_policy = 'structural_zero',
            structural_zero_enters_downstream = TRUE,
            aggregation_order = paste(
                'single_cell_TF_times_ATAC_then_transform_then_project_with',
                'nonestimable_edges_zero_then_mean'
            )
        ),
        class = c('ConditionGRNGroupProjection', 'list')
    )
}
