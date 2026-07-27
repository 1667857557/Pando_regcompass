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
    condition_weight,
    alpha,
    condition_mix,
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
        fit_engine = 'multitask_sparse_group_glmnet',
        network_level = network_level,
        cell_type = cell_type,
        condition = condition,
        cell_type_col = cell_type_col,
        condition_col = condition_col,
        condition_levels = condition_levels,
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
    active_tol
) {
    data.frame(
        network_id = network_id,
        network_name = network_name,
        cell_type = cell_type,
        network_level = network_level,
        condition = condition,
        fit_engine = 'multitask_sparse_group_glmnet',
        n_cells = as.integer(n_cells),
        n_targets = as.integer(n_targets),
        n_candidate_edges = as.integer(nrow(coefs)),
        n_active_edges = as.integer(sum(abs(coefs$estimate) > active_tol)),
        stringsAsFactors = FALSE
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
