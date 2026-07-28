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
        n_active_edges = as.integer(sum(abs(coefs$estimate) > active_tol)),
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
    contrast <- bind_matrix('contrast')
    eligibility_mask <- bind_matrix('eligibility_mask')
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
    structure(
        list(
            schema_version = 'pando_condition_grn_fit_v2',
            network_name = network_name,
            cell_type = cell_type,
            cell_type_col = cell_type_col,
            condition_col = condition_col,
            condition_levels = colnames(beta),
            reference_condition = reference_condition,
            fit_engine = fit_engine,
            candidate_screen = candidate_screen,
            coefficient_scale = if (isTRUE(scale)) {
                'pooled_cell_type_edge_and_target_standardized'
            } else {
                'input_assay_scale'
            },
            edge_table = edge_table,
            beta = beta,
            contrast = contrast,
            eligibility_mask = eligibility_mask,
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
            universal_summary = 'equal-condition coefficient mean for Network compatibility only',
            contrast_formula = 'beta_condition - beta_reference'
        ),
        class = c('ConditionGRNFit', 'list')
    )
}

#' Extract complete condition-aware GRN fit contracts
#'
#' @param object A GRNData object returned by infer_condition_grn().
#' @param network_name Optional network prefix.
#' @param cell_type Optional cell-type label.
#' @return One ConditionGRNFit object when exactly one fit matches, otherwise
#' a named list of matching ConditionGRNFit objects. Each fit contains the
#' exact edge dictionary, condition coefficient matrix, reference contrasts,
#' edge-by-condition eligibility mask, pooled transformations, target-specific
#' lambda selection, condition-level fit quality, and reproducibility metadata.
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
