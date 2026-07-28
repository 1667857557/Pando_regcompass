# Pando-compatible output helpers for condition-aware networks.

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
    reference_condition,
    fit_contract_key,
    lambda_selection,
    nlambda,
    nfolds,
    scale,
    active_tol,
    seed,
    upstream,
    downstream,
    extend,
    only_tss,
    peak_to_gene_method,
    tf_cor,
    peak_cor
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
        reference_condition = reference_condition,
        fit_contract_key = fit_contract_key,
        candidate_screen = candidate_screen,
        peak_to_gene_method = peak_to_gene_method,
        condition_weight = condition_weight,
        alpha = alpha,
        condition_mix = condition_mix,
        lambda_selection = lambda_selection,
        nlambda = nlambda,
        nfolds = nfolds,
        scale = scale,
        active_tol = active_tol,
        seed = seed
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
    reference_condition,
    fit_contract_key
) {
    data.frame(
        network_id = network_id,
        network_name = network_name,
        cell_type = cell_type,
        network_level = network_level,
        condition = condition,
        fit_engine = fit_engine,
        reference_condition = reference_condition,
        fit_contract_key = fit_contract_key,
        n_cells = as.integer(n_cells),
        n_targets = as.integer(n_targets),
        n_candidate_edges = as.integer(nrow(coefs)),
        n_active_edges = as.integer(sum(
            abs(coefs$estimate) > active_tol, na.rm = TRUE
        )),
        stringsAsFactors = FALSE
    )
}

.condition_combine_fit_contracts <- function(
    successful,
    network_name,
    cell_type,
    cell_type_col,
    condition_col,
    reference_condition,
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
    bind_matrix <- function(field) {
        value <- do.call(rbind, lapply(contracts, `[[`, field))
        rownames(value) <- edge_table$edge_id
        value
    }
    beta <- bind_matrix('beta')
    beta_selection <- bind_matrix('beta_selection')
    beta_condition <- bind_matrix('beta_condition')
    delta_condition <- bind_matrix('delta_condition')
    contrast <- bind_matrix('contrast')
    estimability_mask <- bind_matrix('estimability_mask')
    support_mask <- bind_matrix('support_mask')
    active_mask <- bind_matrix('active_mask')
    beta_shared <- unlist(
        lapply(contracts, `[[`, 'beta_shared'), use.names = FALSE
    )
    names(beta_shared) <- edge_table$edge_id
    predictor_center <- unlist(
        lapply(contracts, `[[`, 'predictor_center'), use.names = FALSE
    )
    predictor_scale <- unlist(
        lapply(contracts, `[[`, 'predictor_scale'), use.names = FALSE
    )
    names(predictor_center) <- names(predictor_scale) <- edge_table$edge_id
    response_transform <- data.frame(
        target = vapply(contracts, `[[`, character(1), 'target'),
        center = vapply(contracts, `[[`, numeric(1), 'response_center'),
        scale = vapply(contracts, `[[`, numeric(1), 'response_scale'),
        reference_condition = reference_condition,
        stringsAsFactors = FALSE
    )
    target_fit <- data.frame(
        target = vapply(contracts, `[[`, character(1), 'target'),
        selected_lambda = vapply(contracts, `[[`, numeric(1), 'selected_lambda'),
        alpha = vapply(contracts, `[[`, numeric(1), 'alpha'),
        condition_mix = vapply(contracts, `[[`, numeric(1), 'condition_mix'),
        stringsAsFactors = FALSE
    )
    target_rsq <- do.call(rbind, lapply(contracts, `[[`, 'condition_rsq'))
    rownames(target_rsq) <- response_transform$target
    intercept <- do.call(rbind, lapply(contracts, `[[`, 'intercept'))
    rownames(intercept) <- response_transform$target
    response_scale_by_edge <- response_transform$scale[
        match(edge_table$target, response_transform$target)
    ]
    raw_factor <- response_scale_by_edge / predictor_scale
    beta_condition_raw <- beta_condition * raw_factor
    delta_condition_raw <- delta_condition * raw_factor
    beta_shared_raw <- beta_shared * raw_factor
    topology_mask <- matrix(
        TRUE, nrow(edge_table), ncol(beta_condition),
        dimnames = dimnames(beta_condition)
    )
    comparison_mask <- .condition_reference_comparison_mask(
        estimability_mask, reference_condition
    )
    structure(
        list(
            schema_version = 'pando_condition_grn_fit_v3',
            network_name = network_name,
            cell_type = cell_type,
            cell_type_col = cell_type_col,
            condition_col = condition_col,
            condition_levels = colnames(beta),
            reference_condition = reference_condition,
            fit_engine = fit_engine,
            candidate_screen = candidate_screen,
            coefficient_scale = 'pooled_cell_type_standardized_refit',
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
            contrast = contrast,
            topology_mask = topology_mask,
            eligibility_mask = estimability_mask,
            estimability_mask = estimability_mask,
            support_mask = support_mask,
            active_mask = active_mask,
            absolute_direction = sign(beta_condition),
            deviation_direction = sign(delta_condition),
            comparison_mask = comparison_mask,
            predictor_transform = data.frame(
                edge_id = edge_table$edge_id,
                center = predictor_center,
                scale = predictor_scale,
                stringsAsFactors = FALSE
            ),
            response_transform = response_transform,
            intercept = intercept,
            target_fit = target_fit,
            target_rsq = target_rsq,
            lambda_path = stats::setNames(
                lapply(contracts, `[[`, 'lambda_path'),
                response_transform$target
            ),
            cv = stats::setNames(
                lapply(contracts, function(x) {
                    list(mean = x$cv_mean, se = x$cv_se)
                }),
                response_transform$target
            ),
            refit = stats::setNames(
                lapply(contracts, `[[`, 'refit'),
                response_transform$target
            ),
            condition_weight = contracts[[1L]]$condition_weight,
            active_tol = contracts[[1L]]$active_tol,
            universal_summary = 'estimability-aware shared baseline',
            contrast_formula = 'beta_condition - beta_reference',
            pairwise_contrast_formula =
                'beta_condition[, condition_2] - beta_condition[, condition_1]',
            comparison_mask_formula =
                'estimable_in_condition_1 AND estimable_in_condition_2',
            condition_subgraph_definition =
                'active_mask on the shared TF-peak-target candidate supergraph',
            direction_semantics = list(
                absolute = 'sign(beta_condition)',
                deviation = 'sign(delta_condition)',
                pairwise =
                    'sign(beta_condition[, condition_2] - beta_condition[, condition_1])'
            ),
            normalization_contract = list(
                input = 'paired_single_cells',
                predictor = 'RNA_TF * ATAC_peak',
                predictor_transform = 'one pooled center and scale per edge',
                response_transform = 'one pooled center and scale per target',
                condition_specific_renormalization = FALSE
            ),
            projection_contract = list(
                function_name = 'project_condition_grn_cells',
                score = 'sum(z_edge * beta_condition) by target and observed condition',
                signed = TRUE,
                nonestimable = 'propagate_NA',
                metacell_aggregation = 'column mean of single-cell scores',
                refit_after_aggregation = FALSE
            )
        ),
        class = c('ConditionGRNFit', 'list')
    )
}

.condition_reference_comparison_mask <- function(
    estimability_mask, reference_condition
) {
    estimability_mask <- as.matrix(estimability_mask)
    if (!is.logical(estimability_mask) || anyNA(estimability_mask) ||
        is.null(rownames(estimability_mask)) ||
        is.null(colnames(estimability_mask))) {
        stop('estimability_mask must be a named logical matrix without NA values.')
    }
    if (!is.character(reference_condition) ||
        length(reference_condition) != 1L ||
        !reference_condition %in% colnames(estimability_mask)) {
        stop('reference_condition must identify one estimability-mask column.')
    }
    reference_estimable <- estimability_mask[, reference_condition]
    out <- estimability_mask & matrix(
        reference_estimable,
        nrow = nrow(estimability_mask),
        ncol = ncol(estimability_mask)
    )
    dimnames(out) <- dimnames(estimability_mask)
    out
}

#' Extract complete condition-aware GRN fit contracts
#'
#' @param object A GRNData object returned by infer_condition_grn().
#' @param network_name Optional network prefix.
#' @param cell_type Optional cell-type label.
#' @return One ConditionGRNFit object when exactly one fit matches, otherwise
#' a named list of matching ConditionGRNFit objects. Each fit contains the
#' exact edge dictionary, condition coefficient matrix, reference contrasts,
#' edge-by-condition estimability and activity masks, pooled transformations,
#' the fitted cell/assay contract, target-specific lambda selection,
#' condition-level fit quality, and reproducibility metadata.
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
        stop('No ConditionGRNFit contracts were found in this GRNData object.')
    }
    keep <- vapply(fits, function(x) {
        (is.null(network_name) || identical(x$network_name, network_name)) &&
            (is.null(cell_type) || identical(x$cell_type, cell_type))
    }, logical(1))
    fits <- fits[keep]
    if (!length(fits)) {
        stop('No ConditionGRNFit contract matched the requested filters.')
    }
    if (length(fits) == 1L) fits[[1L]] else fits
}

.condition_require_v3 <- function(fit) {
    if (!inherits(fit, 'ConditionGRNFit') ||
        !identical(fit$schema_version, 'pando_condition_grn_fit_v3')) {
        stop('fit must be a pando_condition_grn_fit_v3 ConditionGRNFit object.')
    }
    invisible(TRUE)
}

#' Compare two condition-specific sub-GRNs
#'
#' @param fit A `pando_condition_grn_fit_v3` object.
#' @param condition_1 Baseline condition label.
#' @param condition_2 Comparison condition label.
#' @param scale Coefficient scale, standardized or raw assay units.
#' @param active_tol Minimum absolute standardized effect used for activity and
#' sign-switch calls, independent of the coefficient display scale.
#' @return A data frame on the shared edge dictionary with absolute effects,
#' pairwise differences, directions, activity, estimability and sign switches.
#' @export
condition_grn_contrast <- function(
    fit,
    condition_1,
    condition_2,
    scale = c('std', 'raw'),
    active_tol = fit$active_tol
) {
    .condition_require_v3(fit)
    scale <- match.arg(scale)
    conditions <- c(as.character(condition_1), as.character(condition_2))
    if (length(unique(conditions)) != 2L ||
        !all(conditions %in% fit$condition_levels)) {
        stop('condition_1 and condition_2 must identify two fitted conditions.')
    }
    beta <- if (scale == 'std') {
        fit$beta_condition_std
    } else {
        fit$beta_condition_raw
    }
    beta_1 <- beta[, conditions[[1L]]]
    beta_2 <- beta[, conditions[[2L]]]
    estimable_1 <- fit$estimability_mask[, conditions[[1L]]]
    estimable_2 <- fit$estimability_mask[, conditions[[2L]]]
    comparable <- estimable_1 & estimable_2
    delta_beta <- beta_2 - beta_1
    delta_beta[!comparable] <- NA_real_
    beta_activity_1 <- fit$beta_condition_std[, conditions[[1L]]]
    beta_activity_2 <- fit$beta_condition_std[, conditions[[2L]]]
    active_1 <- estimable_1 & abs(beta_activity_1) > active_tol
    active_2 <- estimable_2 & abs(beta_activity_2) > active_tol
    sign_switch <- comparable & active_1 & active_2 &
        beta_1 * beta_2 < 0
    data.frame(
        fit$edge_table,
        condition_1 = conditions[[1L]],
        condition_2 = conditions[[2L]],
        beta_condition_1 = beta_1,
        beta_condition_2 = beta_2,
        delta_beta = delta_beta,
        absolute_direction_1 = sign(beta_1),
        absolute_direction_2 = sign(beta_2),
        pairwise_direction = sign(delta_beta),
        estimable_1 = estimable_1,
        estimable_2 = estimable_2,
        comparable = comparable,
        active_1 = active_1,
        active_2 = active_2,
        sign_switch = sign_switch,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
}

#' Extract one active condition sub-GRN
#'
#' @param fit A `pando_condition_grn_fit_v3` object.
#' @param condition Fitted condition label.
#' @param scale Coefficient scale, standardized or raw assay units.
#' @param active_tol Minimum absolute standardized effect retained as active,
#' independent of the coefficient display scale.
#' @return A list containing the active edge table and its induced TF, peak and
#' target node sets.
#' @export
condition_grn_subgraph <- function(
    fit,
    condition,
    scale = c('std', 'raw'),
    active_tol = fit$active_tol
) {
    .condition_require_v3(fit)
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
    active <- estimable &
        abs(fit$beta_condition_std[, condition]) > active_tol
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
    if (is.factor(x)) {
        return(levels(droplevels(x)))
    }
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
