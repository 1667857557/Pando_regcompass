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
#' @param fit Optional `pando_condition_grn_fit_v3` object. When omitted it is
#' retrieved from `object`.
#' @param network_name,cell_type Optional fit filters used when `fit` is omitted.
#' @param component Project absolute condition effects, the shared effect, or
#' condition deviations.
#' @param scale Return scores in pooled standardized target units or raw target
#' units.
#' @param output Return target gene scores or also retain cell-by-edge
#' contributions.
#' @param targets Optional target gene subset.
#' @param nonestimable Propagate unavailable edge effects to an `NA` target score
#' or stop.
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
    active_tol = if (is.null(fit)) 1e-8 else fit$active_tol
) {
    component <- match.arg(component)
    scale <- match.arg(scale)
    output <- match.arg(output)
    nonestimable <- match.arg(nonestimable)
    if (is.null(fit)) {
        fit <- condition_grn_fit(
            object, network_name = network_name, cell_type = cell_type
        )
    }
    .condition_require_v3(fit)
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
    if (any(as.character(metadata[cells, fit$cell_type_col]) != fit$cell_type) ||
        !identical(condition, unname(fit$cell_condition[cells]))) {
        stop('The fitted cell-type or condition labels have changed in the object.')
    }
    edges <- fit$edge_table[fit$edge_table$target %in% targets, , drop = FALSE]
    edge_index <- match(edges$edge_id, fit$edge_table$edge_id)
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

    for (target_index in seq_along(targets)) {
        target <- targets[[target_index]]
        target_edges <- which(edges$target == target)
        missing_by_condition <- stats::setNames(
            rep(FALSE, length(fit$condition_levels)), fit$condition_levels
        )
        for (local_edge in target_edges) {
            raw_predictor <- as.numeric(gene_data[cells, edges$tf[[local_edge]]]) *
                as.numeric(peak_data[cells, edges$region[[local_edge]]])
            predictor <- .condition_projection_predictor(
                raw_predictor = raw_predictor,
                center = transform$center[[local_edge]],
                predictor_scale = transform$scale[[local_edge]],
                scale = scale
            )
            for (condition_name in fit$condition_levels) {
                cell_index <- condition == condition_name
                effect <- coefficient[local_edge, condition_name]
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
            fit$condition_levels,
            function(condition_name) {
                data.frame(
                    target = target,
                    condition = condition_name,
                    n_candidate_edges = length(target_edges),
                    n_estimable_edges = sum(estimability[, condition_name]),
                    n_active_edges = sum(active[, condition_name]),
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
            schema_version = 'pando_condition_grn_projection_v1',
            network_name = fit$network_name,
            cell_type = fit$cell_type,
            component = component,
            scale = scale,
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
                group_within = c(
                    fit$cell_type_col, fit$condition_col, 'sample_or_donor'
                ),
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
