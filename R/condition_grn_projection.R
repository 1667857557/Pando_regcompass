# Fixed-transform single-cell projection for downstream metacell aggregation.

.condition_projection_coefficients <- function(fit, component, scale) {
    if (component == 'condition') {
        return(if (scale == 'std') {
            fit$beta_condition_std
        } else {
            fit$beta_condition_raw
        })
    }
    if (component == 'deviation') {
        return(if (scale == 'std') {
            fit$delta_condition_std
        } else {
            fit$delta_condition_raw
        })
    }
    shared <- if (scale == 'std') fit$beta_shared_std else fit$beta_shared_raw
    out <- matrix(
        shared,
        nrow = length(shared),
        ncol = length(fit$condition_levels),
        dimnames = list(names(shared), fit$condition_levels)
    )
    out
}

.condition_projection_edge_effects <- function(
    fit, coefficient, edge_index, component, active_tol
) {
    edge <- fit$edge_table[edge_index, , drop = FALSE]
    coefficient <- coefficient[edge_index, , drop = FALSE]
    activity_coefficient <- .condition_projection_coefficients(
        fit, component = component, scale = 'std'
    )[edge_index, , drop = FALSE]
    parts <- lapply(fit$condition_levels, function(condition) {
        estimate <- coefficient[, condition]
        estimable <- if (component == 'shared') {
            rep(TRUE, nrow(edge))
        } else {
            fit$estimability_mask[edge_index, condition]
        }
        activity_estimate <- activity_coefficient[, condition]
        data.frame(
            edge,
            condition = condition,
            component = component,
            estimate = estimate,
            direction = sign(estimate),
            estimable = estimable,
            active = estimable & abs(activity_estimate) > active_tol,
            stringsAsFactors = FALSE,
            check.names = FALSE
        )
    })
    do.call(rbind, parts)
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

#' Project a condition-specific GRN to paired single cells
#'
#' @param object The `GRNData` object containing the paired RNA and ATAC cells
#' used for inference.
#' @param fit Optional `pando_condition_grn_fit_v5` object. When omitted it is
#' retrieved from `object`.
#' @param network_name,cell_type Optional fit filters used when `fit` is
#' omitted.
#' @param component Project absolute condition effects, the shared effect, or
#' condition deviations.
#' @param scale Return scores in equal-condition standardized target units or
#' raw target units.
#' @param output Return target gene scores or also retain cell-by-edge
#' contributions.
#' @param targets Optional target gene subset.
#' @param nonestimable Propagate unavailable edge effects to an `NA` target score
#' or stop.
#' @param support_policy Edge support used for projection. Pairwise comparisons
#' use the intersection of estimable edges in two conditions; omnibus
#' comparisons use the global intersection.
#' @param origin Use outer held-out projections for downstream analysis or
#' full-fit coefficients for network interpretation only.
#' @param diagnostic_only Must be `TRUE` to request condition-specific or
#' strict support, which is not a primary cross-condition comparison.
#' @param comparison_conditions Two fitted conditions required by
#' `pairwise_common` when more than two conditions were fitted.
#' @param active_tol Minimum absolute standardized effect recorded as active in
#' the handoff metadata, independent of the returned score scale.
#' @return A `ConditionGRNProjection` containing signed cell-by-target scores,
#' condition-specific edge effects, estimability coverage and an explicit
#' metacell aggregation contract.
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
    nonestimable = c('propagate', 'error'),
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
        stop(
            'condition_estimable and strict are diagnostic-only support ',
            'policies; set diagnostic_only = TRUE explicitly.'
        )
    }
    if (is.null(fit)) {
        fit <- condition_grn_fit(
            object, network_name = network_name, cell_type = cell_type
        )
    }
    .condition_require_v5(fit)
    if (!inherits(object, 'GRNData')) {
        stop('object must be a GRNData object.')
    }
    if (!is.null(targets)) {
        targets <- unique(as.character(targets))
        missing_targets <- setdiff(targets, fit$edge_table$target)
        if (length(missing_targets)) {
            stop(
                'Target(s) were not found in the fit: ',
                paste(missing_targets, collapse = ', '), '.'
            )
        }
    } else {
        targets <- unique(fit$edge_table$target)
    }
    if (identical(origin, 'oof')) {
        if (!identical(component, 'condition')) {
            stop('OOF penalty projections are defined only for component = "condition".')
        }
        if (!identical(fit$projection_origin,
                       'outer_condition_stratified_cell_oof') ||
            !isTRUE(fit$projection_used_for_penalty) ||
            !identical(fit$full_fit_projection_used_for_penalty, FALSE)) {
            stop('The fit does not contain an eligible outer-heldout projection.')
        }
        if (support_policy == 'pairwise_common' &&
            is.null(comparison_conditions)) {
            comparison_conditions <- fit$comparison_conditions
        }
        full_projection <- project_condition_grn_cells(
            object = object,
            fit = fit,
            component = component,
            scale = scale,
            output = 'gene_score',
            targets = targets,
            nonestimable = nonestimable,
            support_policy = support_policy,
            comparison_conditions = comparison_conditions,
            origin = 'full_fit',
            diagnostic_only = diagnostic_only,
            active_tol = active_tol
        )
        source <- if (support_policy == 'global_common') {
            fit$projection_global_common_oof
        } else if (support_policy == 'pairwise_common') {
            requested <- unique(as.character(comparison_conditions))
            fitted <- unique(as.character(fit$comparison_conditions))
            if (!setequal(requested, fitted)) {
                stop(
                    'The requested pair differs from the pair used during ',
                    'outer cross-fitting; refit with those comparison_conditions.'
                )
            }
            fit$projection_common_oof
        } else {
            fit$projection_condition_full_oof
        }
        cells <- rownames(full_projection$gene_score)
        if (is.null(source) ||
            !all(cells %in% rownames(source)) ||
            !all(targets %in% colnames(source))) {
            stop('Stored OOF projection is incomplete for the requested cells or targets.')
        }
        score <- source[cells, targets, drop = FALSE]
        if (scale == 'raw') {
            response_scale <- fit$response_transform$scale[
                match(targets, fit$response_transform$target)
            ]
            score <- sweep(score, 2L, response_scale, '*')
        }
        full_projection$schema_version <-
            'pando_condition_grn_projection_v3'
        full_projection$gene_score <- score
        full_projection$gene_direction <- sign(score)
        full_projection$edge_contribution <- NULL
        full_projection$projection_origin <-
            'outer_condition_stratified_cell_oof'
        primary_common_support <- support_policy %in% c(
            'pairwise_common', 'global_common'
        )
        full_projection$projection_used_for_penalty <-
            primary_common_support && !isTRUE(diagnostic_only)
        full_projection$full_fit_projection_used_for_penalty <- FALSE
        full_projection$gene_score_available <- is.finite(score)
        full_projection$score_comparability_class <-
            if (primary_common_support) {
                'primary_common_support_comparable'
            } else {
                'exploratory_condition_full_not_strictly_comparable'
            }
        full_projection$projection_role <- if (
            isTRUE(full_projection$projection_used_for_penalty)
        ) {
            'primary_penalty'
        } else {
            'exploratory_only'
        }
        full_projection$aggregation_contract$projection_origin <-
            'outer_condition_stratified_cell_oof'
        full_projection$aggregation_contract$projection_role <-
            full_projection$projection_role
        return(full_projection)
    }

    metadata <- object@data@meta.data
    required_metadata <- c(fit$cell_type_col, fit$condition_col)
    if (!all(required_metadata %in% colnames(metadata))) {
        stop('The fitted cell-type or condition metadata column is missing.')
    }
    params <- Params(object)
    if (is.null(fit$assay_contract) ||
        !identical(params$rna_assay, fit$assay_contract$rna_assay) ||
        !identical(params$peak_assay, fit$assay_contract$peak_assay)) {
        stop('The object RNA or ATAC assay no longer matches the fitted contract.')
    }
    gene_data <- Matrix::t(LayerData(
        object,
        assay = params$rna_assay,
        layer = fit$assay_contract$rna_layer
    ))
    peak_data <- Matrix::t(LayerData(
        object,
        assay = params$peak_assay,
        layer = fit$assay_contract$peak_layer
    ))
    cells <- fit$cell_ids
    if (is.null(cells) || !length(cells) ||
        !all(cells %in% rownames(metadata)) ||
        !all(cells %in% rownames(gene_data)) ||
        !all(cells %in% rownames(peak_data))) {
        stop('The object no longer contains every paired single cell used by the fit.')
    }
    condition <- as.character(metadata[cells, fit$condition_col])
    if (any(
            as.character(metadata[cells, fit$cell_type_col]) !=
                fit$cell_type
        ) ||
        !identical(condition, unname(fit$cell_condition[cells]))) {
        stop('The fitted cell-type or condition labels have changed in the object.')
    }
    edges <- fit$edge_table[fit$edge_table$target %in% targets, , drop = FALSE]
    edge_index <- match(edges$edge_id, fit$edge_table$edge_id)
    projected_conditions <- fit$condition_levels
    if (support_policy == 'pairwise_common') {
        if (is.null(comparison_conditions)) {
            if (length(fit$condition_levels) != 2L) {
                stop(
                    'comparison_conditions must name two fitted conditions ',
                    'when pairwise_common is used with more than two conditions.'
                )
            }
            comparison_conditions <- fit$condition_levels
        }
        comparison_conditions <- unique(as.character(comparison_conditions))
        if (length(comparison_conditions) != 2L ||
            !all(comparison_conditions %in% fit$condition_levels)) {
            stop('comparison_conditions must identify exactly two fitted conditions.')
        }
        projected_conditions <- comparison_conditions
        keep_cells <- condition %in% projected_conditions
        cells <- cells[keep_cells]
        condition <- condition[keep_cells]
        gene_data <- gene_data[cells, , drop = FALSE]
        peak_data <- peak_data[cells, , drop = FALSE]
    }
    coefficient <- .condition_projection_coefficients(
        fit, component = component, scale = scale
    )[edge_index, , drop = FALSE]
    transform <- fit$predictor_transform[
        match(edges$edge_id, fit$predictor_transform$edge_id), , drop = FALSE
    ]
    missing_tfs <- setdiff(unique(edges$tf), colnames(gene_data))
    missing_peaks <- setdiff(unique(edges$region), colnames(peak_data))
    if (length(missing_tfs) || length(missing_peaks)) {
        stop('The object no longer contains every TF and peak required by the fit.')
    }

    gene_score <- matrix(
        0,
        nrow = length(cells),
        ncol = length(targets),
        dimnames = list(cells, targets)
    )
    edge_contribution <- if (output == 'edge_contribution') {
        matrix(
            NA_real_,
            nrow = length(cells),
            ncol = nrow(edges),
            dimnames = list(cells, edges$edge_id)
        )
    } else {
        NULL
    }
    target_status <- vector('list', length(targets))
    estimability_all <- fit$estimability_mask[edge_index, , drop = FALSE]
    common_global <- rowSums(estimability_all) == ncol(estimability_all)
    support_mask <- estimability_all
    if (support_policy == 'global_common') {
        support_mask[,] <- common_global
    } else if (support_policy == 'pairwise_common') {
        common_pair <- rowSums(
            estimability_all[, projected_conditions, drop = FALSE]
        ) == 2L
        support_mask[,] <- FALSE
        support_mask[, projected_conditions] <- common_pair
    }

    for (target_index in seq_along(targets)) {
        target <- targets[[target_index]]
        target_edges <- which(edges$target == target)
        missing_by_condition <- stats::setNames(
            rep(FALSE, length(projected_conditions)), projected_conditions
        )
        tf_activity <- matrix(NA_real_, nrow = length(cells),
                              ncol = length(target_edges))
        peak_activity <- tf_activity
        product_activity <- tf_activity
        for (local_edge in target_edges) {
            tf_value <- as.numeric(gene_data[cells, edges$tf[[local_edge]]])
            peak_value <- as.numeric(peak_data[cells, edges$region[[local_edge]]])
            raw_predictor <- tf_value * peak_value
            diagnostic_column <- match(local_edge, target_edges)
            tf_activity[, diagnostic_column] <- tf_value
            peak_activity[, diagnostic_column] <- peak_value
            product_activity[, diagnostic_column] <- raw_predictor
            predictor <- .condition_projection_predictor(
                raw_predictor = raw_predictor,
                center = transform$center[[local_edge]],
                predictor_scale = transform$scale[[local_edge]],
                scale = scale
            )
            for (condition_name in projected_conditions) {
                cell_index <- condition == condition_name
                effect <- coefficient[local_edge, condition_name]
                included <- support_mask[local_edge, condition_name]
                if (!included) {
                    if (support_policy == 'strict') {
                        missing_by_condition[[condition_name]] <- TRUE
                    }
                    next
                }
                if (is.na(effect)) {
                    missing_by_condition[[condition_name]] <- TRUE
                    next
                }
                contribution <- predictor[cell_index] * effect
                gene_score[cell_index, target_index] <-
                    gene_score[cell_index, target_index] + contribution
                if (!is.null(edge_contribution)) {
                    edge_contribution[cell_index, local_edge] <- contribution
                }
            }
        }
        estimability <- fit$estimability_mask[
            match(edges$edge_id[target_edges], fit$edge_table$edge_id),
            ,
            drop = FALSE
        ]
        active <- fit$active_mask[
            match(edges$edge_id[target_edges], fit$edge_table$edge_id),
            ,
            drop = FALSE
        ]
        target_status[[target_index]] <- do.call(rbind, lapply(
            projected_conditions,
            function(condition_name) {
                target_cell <- condition == condition_name
                target_support <- support_mask[target_edges, condition_name]
                product_variance <- apply(
                    product_activity[target_cell, , drop = FALSE],
                    2L, stats::var
                )
                n_common_global <- sum(
                    estimability_all[target_edges, condition_name] &
                    common_global[target_edges]
                )
                n_projection_support <- sum(target_support)
                n_common_policy <- if (support_policy %in% c(
                    'global_common', 'pairwise_common'
                )) {
                    n_projection_support
                } else {
                    n_common_global
                }
                data.frame(
                    target = target,
                    condition = condition_name,
                    n_candidate_edges = length(target_edges),
                    n_estimable_edges = sum(estimability[, condition_name]),
                    n_common_estimable_edges = n_common_policy,
                    n_global_common_estimable_edges = n_common_global,
                    n_projection_support_edges = n_projection_support,
                    n_active_edges = sum(active[, condition_name]),
                    estimable_fraction =
                        sum(estimability[, condition_name]) / length(target_edges),
                    common_estimable_fraction =
                        n_common_policy / length(target_edges),
                    global_common_estimable_fraction =
                        n_common_global / length(target_edges),
                    tf_detection_fraction = mean(
                        tf_activity[target_cell, , drop = FALSE] != 0
                    ),
                    peak_detection_fraction = mean(
                        peak_activity[target_cell, , drop = FALSE] != 0
                    ),
                    tf_peak_nonzero_fraction = mean(
                        product_activity[target_cell, , drop = FALSE] != 0
                    ),
                    tf_peak_variance = if (all(!is.finite(product_variance))) {
                        NA_real_
                    } else {
                        mean(product_variance[is.finite(product_variance)])
                    },
                    score_completeness =
                        sum(target_support) / length(target_edges),
                    fully_estimable = all(estimability[, condition_name]),
                    stringsAsFactors = FALSE
                )
            }
        ))
        if (any(missing_by_condition)) {
            if (nonestimable == 'error') {
                stop(
                    'Projection is not fully estimable for target ', target,
                    ' in condition(s): ',
                    paste(names(missing_by_condition)[missing_by_condition],
                          collapse = ', '), '.'
                )
            }
            for (condition_name in names(missing_by_condition)[missing_by_condition]) {
                gene_score[condition == condition_name, target_index] <- NA_real_
            }
        }
        for (condition_name in projected_conditions) {
            if (!any(support_mask[target_edges, condition_name])) {
                gene_score[condition == condition_name, target_index] <- NA_real_
            }
        }
    }

    cell_metadata <- data.frame(
        cell_id = cells,
        cell_type = fit$cell_type,
        condition = condition,
        stringsAsFactors = FALSE,
        row.names = cells
    )
    structure(
        list(
            schema_version = 'pando_condition_grn_projection_v2',
            network_name = fit$network_name,
            cell_type = fit$cell_type,
            component = component,
            scale = scale,
            support_policy = support_policy,
            comparison_conditions = comparison_conditions,
            projection_origin = 'full_data_fit_interpretation_only',
            projection_used_for_penalty = FALSE,
            full_fit_projection_used_for_penalty = FALSE,
            cell_metadata = cell_metadata,
            gene_score = gene_score,
            gene_direction = sign(gene_score),
            edge_contribution = edge_contribution,
            edge_effects = .condition_projection_edge_effects(
                fit, coefficient = .condition_projection_coefficients(
                    fit, component, scale
                ), edge_index = edge_index, component = component,
                active_tol = active_tol
            ),
            target_condition_status = do.call(rbind, target_status),
            aggregation_contract = list(
                group_within = c(fit$cell_type_col, fit$condition_col),
                operation = 'arithmetic_mean_by_target',
                signed_scores = TRUE,
                identical_target_columns = TRUE,
                propagate_NA = TRUE,
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
#' @param projection A `ConditionGRNProjection` returned by
#' `project_condition_grn_cells()`.
#' @param membership One row per cell with `cell_id` and `group_col`.
#' @param group_col Membership column identifying the output group.
#' @return A `ConditionGRNGroupProjection` whose scores are arithmetic means of
#' already-computed single-cell projections.
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
        stop(
            'membership must map each non-empty cell ID exactly once ',
            'to one non-empty group.'
        )
    }
    cells <- rownames(projection$gene_score)
    index <- match(cells, membership$cell_id)
    if (anyNA(index)) {
        stop('membership does not cover every projected single cell.')
    }
    membership <- membership[index, , drop = FALSE]
    groups <- unique(membership[[group_col]])
    score <- do.call(rbind, lapply(groups, function(group) {
        rows <- membership[[group_col]] == group
        colMeans(projection$gene_score[rows, , drop = FALSE], na.rm = FALSE)
    }))
    rownames(score) <- groups
    group_metadata <- do.call(rbind, lapply(groups, function(group) {
        rows <- membership[[group_col]] == group
        observed <- projection$cell_metadata[cells[rows], , drop = FALSE]
        invariant <- c('cell_type', 'condition')
        values <- lapply(invariant, function(field) {
            value <- unique(as.character(observed[[field]]))
            value <- value[!is.na(value) & nzchar(value)]
            if (length(value) > 1L) {
                stop('Group ', group, ' mixes ', field, ' values.')
            }
            if (length(value)) value else NA_character_
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
            schema_version = 'pando_condition_grn_group_projection_v1',
            source_projection = projection,
            projection_origin = projection$projection_origin,
            projection_used_for_penalty =
                isTRUE(projection$projection_used_for_penalty),
            full_fit_projection_used_for_penalty =
                isTRUE(projection$full_fit_projection_used_for_penalty),
            score_comparability_class =
                projection$score_comparability_class,
            group_col = group_col,
            group_metadata = group_metadata,
            gene_score = score,
            gene_direction = sign(score),
            aggregation_order =
                'single_cell_TF_times_ATAC_then_transform_then_project_then_mean'
        ),
        class = c('ConditionGRNGroupProjection', 'list')
    )
}
