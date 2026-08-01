# Pando-compatible output helpers for condition-aware networks.

.PANDO_CONDITION_GRN_FIT_SCHEMA <- 'pando_condition_grn_fit'

.condition_format_coefs <- function(edges, estimate, corr) {
    data.frame(
        tf = edges$tf,
        target = edges$target,
        region = edges$region,
        term = edges$term,
        estimate = as.numeric(estimate),
        corr = as.numeric(corr),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
}

.condition_build_network <- function(features, coefs, fit, params) {
    coefs <- coefs[, c('tf', 'target', 'region', 'term', 'estimate', 'corr'), drop = FALSE]
    fit <- fit[, c('target', 'lambda', 'rsq', 'alpha', 'nvariables'), drop = FALSE]
    new(
        Class = 'Network',
        features = as.character(features),
        coefs = as.data.frame(coefs),
        fit = as.data.frame(fit),
        params = params
    )
}

.condition_network_params <- function(
    network_level,
    cell_type,
    condition,
    cell_type_col,
    condition_col,
    condition_levels,
    candidate_screen,
    fit_engine,
    condition_weight,
    alpha,
    condition_mix,
    reference_condition = NULL,
    comparison_conditions = NULL,
    fit_contract_key,
    lambda_selection,
    nlambda,
    outer_nfolds,
    inner_nfolds,
    oof_scheme,
    scale,
    active_tol,
    seed,
    upstream,
    downstream,
    extend,
    only_tss,
    peak_to_gene_method,
    tf_cor,
    peak_cor,
    engine_control = NULL
) {
    list(
        method = 'glmnet',
        family = 'gaussian',
        dist = c(upstream = upstream, downstream = downstream, extend = extend),
        only_tss = only_tss,
        interaction = ':',
        tf_cor = tf_cor,
        peak_cor = peak_cor,
        fit_engine = fit_engine,
        network_level = network_level,
        cell_type = cell_type,
        condition = condition,
        cell_type_col = cell_type_col,
        condition_col = condition_col,
        condition_levels = condition_levels,
        comparison_conditions = comparison_conditions,
        fit_contract_key = fit_contract_key,
        candidate_screen = candidate_screen,
        peak_to_gene_method = peak_to_gene_method,
        condition_weight = condition_weight,
        alpha = alpha,
        condition_mix = condition_mix,
        lambda_selection = lambda_selection,
        nlambda = nlambda,
        outer_nfolds = outer_nfolds,
        inner_nfolds = inner_nfolds,
        oof_scheme = oof_scheme,
        scale = scale,
        active_tol = active_tol,
        seed = seed,
        engine_control = engine_control,
        coefficient_contract = 'absolute_condition_effects_only'
    )
}

.condition_index_row <- function(
    network_id,
    network_name,
    cell_type,
    network_level,
    condition,
    n_cells,
    coefs,
    n_targets,
    active_tol,
    fit_engine,
    reference_condition = NULL,
    fit_contract_key
) {
    data.frame(
        network_id = network_id,
        network_name = network_name,
        cell_type = cell_type,
        network_level = network_level,
        condition = condition,
        fit_engine = fit_engine,
        fit_contract_key = fit_contract_key,
        n_cells = as.integer(n_cells),
        n_targets = as.integer(n_targets),
        n_candidate_edges = as.integer(nrow(coefs)),
        n_active_edges = as.integer(sum(abs(coefs$estimate) > active_tol, na.rm = TRUE)),
        coefficient_contract = 'absolute_condition_effects_only',
        stringsAsFactors = FALSE
    )
}

.condition_bind_contract_matrix <- function(contracts, field, edge_ids) {
    value <- do.call(rbind, lapply(contracts, `[[`, field))
    rownames(value) <- edge_ids
    value
}

.condition_bind_oof_projection <- function(contracts, field, cells, targets) {
    out <- matrix(
        NA_real_,
        nrow = length(cells),
        ncol = length(contracts),
        dimnames = list(cells, targets)
    )
    for (i in seq_along(contracts)) {
        value <- contracts[[i]][[field]]
        if (is.null(names(value)) || !setequal(names(value), cells)) {
            stop('OOF projection cells differ across target contracts.')
        }
        out[, i] <- as.numeric(value[cells])
    }
    out
}

.condition_combine_fit_contracts <- function(
    successful,
    network_name,
    cell_type,
    cell_type_col,
    condition_col,
    reference_condition = NULL,
    comparison_conditions = NULL,
    candidate_screen,
    scale,
    fit_engine
) {
    contracts <- lapply(successful, `[[`, 'fit_contract')
    edge_table <- do.call(rbind, lapply(contracts, `[[`, 'edge_table'))
    rownames(edge_table) <- NULL
    if (anyDuplicated(edge_table$edge_id)) {
        stop('ConditionGRNFit contains duplicated edge identifiers.')
    }
    edge_ids <- edge_table$edge_id
    beta <- .condition_bind_contract_matrix(contracts, 'beta', edge_ids)
    beta_selection <- .condition_bind_contract_matrix(
        contracts, 'beta_selection', edge_ids
    )
    beta_condition <- .condition_bind_contract_matrix(
        contracts, 'beta_condition', edge_ids
    )
    delta_condition <- .condition_bind_contract_matrix(
        contracts, 'delta_condition', edge_ids
    )
    structural_candidate_mask <- .condition_bind_contract_matrix(
        contracts, 'structural_candidate_mask', edge_ids
    )
    screening_mask <- .condition_bind_contract_matrix(
        contracts, 'screening_mask', edge_ids
    )
    estimability_mask <- .condition_bind_contract_matrix(
        contracts, 'estimability_mask', edge_ids
    )
    support_mask <- .condition_bind_contract_matrix(
        contracts, 'support_mask', edge_ids
    )
    active_mask <- .condition_bind_contract_matrix(
        contracts, 'active_mask', edge_ids
    )
    beta_shared <- unlist(lapply(contracts, `[[`, 'beta_shared'), use.names = FALSE)
    names(beta_shared) <- edge_ids
    predictor_center <- unlist(
        lapply(contracts, `[[`, 'predictor_center'), use.names = FALSE
    )
    predictor_scale <- unlist(
        lapply(contracts, `[[`, 'predictor_scale'), use.names = FALSE
    )
    names(predictor_center) <- names(predictor_scale) <- edge_ids
    response_transform <- data.frame(
        target = vapply(contracts, `[[`, character(1), 'target'),
        center = vapply(contracts, `[[`, numeric(1), 'response_center'),
        scale = vapply(contracts, `[[`, numeric(1), 'response_scale'),
        stringsAsFactors = FALSE
    )
    target_fit <- data.frame(
        target = response_transform$target,
        selected_lambda = vapply(contracts, `[[`, numeric(1), 'selected_lambda'),
        alpha = vapply(contracts, `[[`, numeric(1), 'alpha'),
        condition_mix = vapply(contracts, `[[`, numeric(1), 'condition_mix'),
        stringsAsFactors = FALSE
    )
    condition_rsq_train <- do.call(rbind, lapply(contracts, `[[`, 'condition_rsq_train'))
    condition_rsq_oof <- do.call(rbind, lapply(contracts, `[[`, 'condition_rsq_oof'))
    condition_rmse_oof <- do.call(rbind, lapply(contracts, `[[`, 'condition_rmse_oof'))
    rownames(condition_rsq_train) <- rownames(condition_rsq_oof) <-
        rownames(condition_rmse_oof) <- response_transform$target
    target_rsq_oof_pooled <- stats::setNames(
        vapply(contracts, `[[`, numeric(1), 'target_rsq_oof_pooled'),
        response_transform$target
    )
    predictive_oof_available <- stats::setNames(
        vapply(contracts, `[[`, logical(1), 'predictive_oof_available'),
        response_transform$target
    )
    oof_validation_level <- stats::setNames(
        vapply(contracts, `[[`, character(1), 'oof_validation_level'),
        response_transform$target
    )
    intercept <- do.call(rbind, lapply(contracts, `[[`, 'intercept'))
    rownames(intercept) <- response_transform$target
    response_scale_by_edge <- response_transform$scale[
        match(edge_table$target, response_transform$target)
    ]
    raw_factor <- response_scale_by_edge / predictor_scale
    beta_condition_raw <- beta_condition * raw_factor
    delta_condition_raw <- delta_condition * raw_factor
    beta_shared_raw <- beta_shared * raw_factor
    oof_cells <- names(contracts[[1L]]$projection_common_oof)
    if (is.null(oof_cells) || anyDuplicated(oof_cells)) {
        stop('OOF projection cells must be unique and named.')
    }
    projection_common_oof <- .condition_bind_oof_projection(
        contracts, 'projection_common_oof', oof_cells, response_transform$target
    )
    projection_condition_full_oof <- .condition_bind_oof_projection(
        contracts, 'projection_condition_full_oof', oof_cells,
        response_transform$target
    )
    projection_global_common_oof <- .condition_bind_oof_projection(
        contracts, 'projection_global_common_oof', oof_cells,
        response_transform$target
    )
    oof_assignment_count <- .condition_bind_oof_projection(
        contracts, 'oof_assignment_count', oof_cells, response_transform$target
    )
    if (any(oof_assignment_count != 1L)) {
        stop('Every target and cell must have exactly one outer-fold assignment.')
    }
    target_transform_contract <- data.frame(
        target = response_transform$target,
        transform_policy = vapply(contracts, `[[`, character(1), 'transform_policy'),
        predictor_center_hash = vapply(
            contracts, `[[`, character(1), 'predictor_center_hash'
        ),
        predictor_scale_hash = vapply(
            contracts, `[[`, character(1), 'predictor_scale_hash'
        ),
        training_fold_only = vapply(
            contracts, `[[`, logical(1), 'training_fold_only'
        ),
        stringsAsFactors = FALSE
    )
    has_refit_stability <- all(vapply(
        contracts, function(value) !is.null(value$refit_stability),
        logical(1)
    ))
    stability_fields <- if (has_refit_stability) {
        list(
            refit_stability = stats::setNames(
                lapply(contracts, `[[`, 'refit_stability'),
                response_transform$target
            ),
            refit_stability_edge = do.call(rbind, lapply(
                contracts, function(value) {
                    edge <- value$refit_stability$edge
                    if (is.null(edge) || !nrow(edge)) return(edge)
                    edge$target <- value$target
                    edge
                }
            ))
        )
    } else {
        list()
    }
    structure(
        c(list(
            schema_version = .PANDO_CONDITION_GRN_FIT_SCHEMA,
            schema_policy = 'single_unversioned_schema',
            contract_version = 'condition_absolute_oof_v3',
            coefficient_contract = 'absolute_condition_effects_only',
            network_name = network_name,
            cell_type = cell_type,
            cell_type_col = cell_type_col,
            condition_col = condition_col,
            condition_levels = colnames(beta),
            comparison_conditions = comparison_conditions,
            fit_engine = fit_engine,
            candidate_screen = candidate_screen,
            coefficient_scale =
                'equal_condition_within_variance_standardized_refit',
            edge_table = edge_table,
            beta = beta,
            beta_selection = beta_selection,
            beta_condition = beta_condition,
            beta_shared = beta_shared,
            delta_condition = delta_condition,
            beta_condition_std = beta_condition,
            beta_shared_std = beta_shared,
            delta_condition_std = delta_condition,
            beta_condition_raw = beta_condition_raw,
            beta_shared_raw = beta_shared_raw,
            delta_condition_raw = delta_condition_raw,
            topology_mask = structural_candidate_mask,
            structural_candidate_mask = structural_candidate_mask,
            screening_mask = screening_mask,
            eligibility_mask = estimability_mask,
            estimability_mask = estimability_mask,
            support_mask = support_mask,
            active_mask = active_mask,
            absolute_direction = sign(beta_condition),
            deviation_direction = sign(delta_condition),
            predictor_transform = data.frame(
                edge_id = edge_ids,
                center = predictor_center,
                scale = predictor_scale,
                stringsAsFactors = FALSE
            ),
            response_transform = response_transform,
            target_transform_contract = target_transform_contract,
            intercept = intercept,
            target_fit = target_fit,
            target_rsq = condition_rsq_train,
            condition_rsq_train = condition_rsq_train,
            condition_rsq_oof = condition_rsq_oof,
            condition_rmse_oof = condition_rmse_oof,
            target_rsq_oof_pooled = target_rsq_oof_pooled,
            target_rsq_oof_pooled_definition =
                '1 - pooled OOF SSE / pooled within-condition SST',
            cv_method = stats::setNames(
                vapply(contracts, `[[`, character(1), 'cv_method'),
                response_transform$target
            ),
            oof_model = stats::setNames(
                vapply(contracts, `[[`, character(1), 'oof_model'),
                response_transform$target
            ),
            predictive_oof_available = predictive_oof_available,
            oof_validation_level = oof_validation_level,
            projection_common_oof = projection_common_oof,
            projection_condition_full_oof = projection_condition_full_oof,
            projection_global_common_oof = projection_global_common_oof,
            projection_origin = 'outer_condition_stratified_cell_oof',
            projection_used_for_penalty = all(vapply(
                contracts, `[[`, logical(1), 'projection_used_for_penalty'
            )),
            full_fit_projection_used_for_penalty = FALSE,
            oof_cell_coverage = stats::setNames(
                vapply(contracts, `[[`, numeric(1), 'oof_cell_coverage'),
                response_transform$target
            ),
            oof_projection_available_fraction = stats::setNames(
                vapply(
                    contracts, `[[`, numeric(1),
                    'oof_projection_available_fraction'
                ),
                response_transform$target
            ),
            oof_assignment_count = oof_assignment_count,
            fold_transform_policy =
                'equal_condition_center_equal_condition_within_variance_v1',
            oof_fold = stats::setNames(
                lapply(contracts, `[[`, 'oof_fold'), response_transform$target
            ),
            cv_fold_transform = stats::setNames(
                lapply(contracts, `[[`, 'cv_fold_transform'),
                response_transform$target
            ),
            cv_effective_nfolds = stats::setNames(
                vapply(contracts, function(x) as.integer(x$cv_effective_nfolds), integer(1)),
                response_transform$target
            ),
            outer_nfolds = stats::setNames(
                vapply(contracts, `[[`, integer(1), 'outer_nfolds'),
                response_transform$target
            ),
            inner_nfolds = stats::setNames(
                vapply(contracts, `[[`, integer(1), 'inner_nfolds'),
                response_transform$target
            ),
            lambda_path = stats::setNames(
                lapply(contracts, `[[`, 'lambda_path'), response_transform$target
            ),
            cv = stats::setNames(
                lapply(contracts, function(x) list(mean = x$cv_mean, se = x$cv_se)),
                response_transform$target
            ),
            refit = stats::setNames(
                lapply(contracts, `[[`, 'refit'), response_transform$target
            ),
            execution = stats::setNames(
                lapply(contracts, `[[`, 'execution'),
                response_transform$target
            ),
            condition_weight = contracts[[1L]]$condition_weight,
            active_tol = contracts[[1L]]$active_tol,
            universal_summary = 'estimability-aware shared baseline',
            condition_subgraph_definition =
                'active_mask on the shared TF-peak-target candidate supergraph',
            direction_semantics = list(
                absolute = 'sign(beta_condition)',
                deviation = 'sign(delta_condition)'
            ),
            normalization_contract = list(
                input = 'paired_single_cells',
                predictor = 'RNA_TF * ATAC_peak',
                predictor_transform =
                    'equal-condition center and within-condition variance',
                response_transform =
                    'equal-condition center and within-condition variance',
                condition_weights = 'exactly 1/K',
                variance_denominator = 'population_for_duplication_invariance',
                condition_specific_renormalization = FALSE
            ),
            projection_contract = list(
                function_name = 'project_condition_grn_cells',
                score = 'outer-heldout sum(z_edge * beta_condition) by target',
                projection_origin = 'outer_condition_stratified_cell_oof',
                projection_used_for_penalty = all(predictive_oof_available),
                full_fit_projection_used_for_penalty = FALSE,
                signed = TRUE,
                support_policy = c(
                    'pairwise_common', 'global_common',
                    'condition_estimable_diagnostic', 'strict_diagnostic'
                ),
                primary_support_policy = 'pairwise_common or global_common only',
                condition_full_role = 'exploratory_only',
                nonestimable = 'structural_zero_before_target_summation',
                metacell_aggregation =
                    'mean of cell-first TF-times-ATAC projections',
                refit_after_aggregation = FALSE
            )
        ), stability_fields),
        class = c('ConditionGRNFit', 'list')
    )
}

#' Extract complete condition-aware GRN fit contracts
#'
#' @param object A GRNData object returned by [infer_condition_grn()].
#' @param network_name Optional network prefix.
#' @param cell_type Optional fitted cell-type label.
#' @return One canonical `pando_condition_grn_fit` object or a named list.
#' @export
condition_grn_fit <- function(object, network_name = NULL, cell_type = NULL) {
    UseMethod('condition_grn_fit')
}

#' @rdname condition_grn_fit
#' @method condition_grn_fit GRNData
#' @export
condition_grn_fit.GRNData <- function(
    object, network_name = NULL, cell_type = NULL
) {
    fits <- object@grn@params$condition_grn_fits
    if (is.null(fits) || !length(fits)) {
        mode <- object@grn@params$analysis_mode
        if (identical(mode, 'standard_grn')) {
            stop('This object used standard infer_grn(); no condition coefficients were calculated.')
        }
        stop('No ConditionGRNFit contracts were found in this GRNData object.')
    }
    keep <- vapply(fits, function(x) {
        (is.null(network_name) || identical(x$network_name, network_name)) &&
            (is.null(cell_type) || identical(x$cell_type, cell_type))
    }, logical(1))
    fits <- fits[keep]
    if (!length(fits)) stop('No ConditionGRNFit contract matched the requested filters.')
    invisible(lapply(fits, .condition_require_fit))
    if (length(fits) == 1L) fits[[1L]] else fits
}

.condition_require_fit <- function(fit) {
    if (!inherits(fit, 'ConditionGRNFit') ||
        !identical(fit$schema_version, .PANDO_CONDITION_GRN_FIT_SCHEMA)) {
        stop('fit must use the canonical pando_condition_grn_fit schema.')
    }
    invisible(TRUE)
}

#' Extract one active condition sub-GRN
#'
#' @param fit A canonical `pando_condition_grn_fit` object.
#' @param condition Fitted condition label.
#' @param scale Coefficient scale, standardized or raw assay units.
#' @param active_tol Minimum absolute standardized effect retained as active.
#' @return Active edges and induced TF, peak, and target node sets.
#' @export
condition_grn_subgraph <- function(
    fit, condition, scale = c('std', 'raw'), active_tol = fit$active_tol
) {
    .condition_require_fit(fit)
    scale <- match.arg(scale)
    condition <- as.character(condition)
    if (length(condition) != 1L || !condition %in% fit$condition_levels) {
        stop('condition must identify one fitted condition.')
    }
    beta <- if (scale == 'std') {
        fit$beta_condition_std[, condition]
    } else {
        fit$beta_condition_raw[, condition]
    }
    estimable <- fit$estimability_mask[, condition]
    active <- estimable & abs(fit$beta_condition_std[, condition]) > active_tol
    edges <- data.frame(
        fit$edge_table[active, , drop = FALSE],
        estimate = beta[active],
        absolute_direction = sign(beta[active]),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    structure(
        list(
            condition = condition,
            scale = scale,
            edges = edges,
            nodes = list(
                tf = unique(edges$tf),
                peak = unique(edges$region),
                target = unique(edges$target)
            )
        ),
        class = c('ConditionSubGRN', 'list')
    )
}

.condition_map <- function(x, fun, parallel, BPPARAM, verbose) {
    if (!is.null(BPPARAM)) {
        if (!requireNamespace('BiocParallel', quietly = TRUE)) {
            stop('BiocParallel is required when BPPARAM is supplied.')
        }
        result <- BiocParallel::bplapply(x, fun, BPPARAM = BPPARAM)
        names(result) <- names(x)
        return(result)
    }
    map_par(x, fun, parallel = parallel, verbose = verbose)
}

.condition_levels <- function(x) {
    if (is.factor(x)) return(levels(droplevels(x)))
    unique(as.character(x))
}

.condition_safe_id <- function(x) {
    stringr::str_replace_all(as.character(x), '[^A-Za-z0-9_.-]', '_')
}

.condition_model_name <- function(x) {
    stringr::str_replace_all(as.character(x), '-', '_')
}

.condition_seed_for <- function(label, seed) {
    offset <- sum(utf8ToInt(enc2utf8(as.character(label))))
    as.integer((as.double(seed) + offset) %% .Machine$integer.max)
}
