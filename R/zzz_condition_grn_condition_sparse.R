# Shared-baseline, condition-sparse extension for condition-aware Pando.
#
# This file intentionally loads after the v2 implementation. It keeps the
# existing Network/output surface while adding:
# - pooled-within-condition candidate screening,
# - sparse-group support selection with condition-specific supports,
# - support-constrained hierarchical-ridge refitting on the pooled scale,
# - explicit shared, condition, and deviation coefficient layers.

.condition_infer_condition_grn_grndata_v2 <- infer_condition_grn.GRNData
.condition_validate_public_args_v2 <- .condition_validate_public_args
.condition_component_masks_v2 <- .condition_component_masks
.condition_run_cell_type_v2 <- .condition_run_cell_type
.condition_combine_fit_contracts_v2 <- .condition_combine_fit_contracts

.condition_replace_call <- function(node, predicate, replacement) {
    if (predicate(node)) {
        return(replacement)
    }
    if (is.call(node) || is.pairlist(node) || is.expression(node)) {
        for (index in seq_along(node)) {
            node[[index]] <- .condition_replace_call(
                node[[index]], predicate, replacement
            )
        }
    }
    node
}

.condition_patch_condition_gof_names <- function(fun) {
    predicate <- function(node) {
        is.call(node) &&
            identical(node[[1L]], as.name('<-')) &&
            identical(node[[2L]], as.name('condition_gof')) &&
            is.call(node[[3L]]) &&
            identical(node[[3L]][[1L]], as.name('vector'))
    }
    replacement <- quote({
        condition_gof <- vector('list', length(condition_levels))
        names(condition_gof) <- condition_levels
    })
    body(fun) <- .condition_replace_call(body(fun), predicate, replacement)
    fun
}

.condition_fit_target_v2_patched <-
    .condition_patch_condition_gof_names(.condition_fit_target)

.condition_within_named_cor <- function(x, y, condition) {
    if (ncol(x) == 0L) {
        return(stats::setNames(numeric(), character()))
    }
    condition <- droplevels(condition)
    y <- as.numeric(y[, 1L])
    cross_total <- numeric(ncol(x))
    ssx_total <- numeric(ncol(x))
    ssy_total <- 0

    for (level in levels(condition)) {
        keep <- condition == level
        n_level <- sum(keep)
        if (n_level <= 1L) {
            next
        }
        x_level <- x[keep, , drop = FALSE]
        y_level <- y[keep]
        sum_x <- as.numeric(Matrix::colSums(x_level))
        sum_y <- sum(y_level)
        cross_total <- cross_total +
            as.numeric(crossprod(x_level, y_level)) -
            sum_x * sum_y / n_level
        ssx_total <- ssx_total +
            as.numeric(Matrix::colSums(x_level * x_level)) -
            sum_x * sum_x / n_level
        ssy_total <- ssy_total +
            sum(y_level * y_level) -
            sum_y * sum_y / n_level
    }

    denominator <- sqrt(pmax(ssx_total, 0) * max(ssy_total, 0))
    correlation <- cross_total / denominator
    correlation[!is.finite(correlation)] <- 0
    stats::setNames(correlation, colnames(x))
}

.condition_component_masks <- function(
    x, y, condition, threshold, candidate_screen
) {
    if (!identical(candidate_screen, 'pooled_within_condition')) {
        return(.condition_component_masks_v2(
            x, y, condition, threshold, candidate_screen
        ))
    }
    score <- abs(.condition_within_named_cor(x, y, condition))
    keep <- is.finite(score) & score > threshold
    matrix(
        keep,
        nrow = length(keep),
        ncol = nlevels(condition),
        dimnames = list(names(score), levels(condition))
    )
}

.condition_average_weights <- function(
    X_list, condition_weight = c('equal', 'cell_count')
) {
    condition_weight <- match.arg(condition_weight)
    n_condition <- vapply(X_list, nrow, integer(1))
    if (condition_weight == 'equal') {
        return(rep(1 / length(X_list), length(X_list)))
    }
    n_condition / sum(n_condition)
}

.condition_weighted_row_mean <- function(B, eligibility, weights) {
    out <- numeric(nrow(B))
    for (edge in seq_len(nrow(B))) {
        use <- eligibility[edge, ]
        if (!any(use)) {
            out[[edge]] <- 0
            next
        }
        edge_weights <- weights[use]
        edge_weights <- edge_weights / sum(edge_weights)
        out[[edge]] <- sum(edge_weights * B[edge, use])
    }
    out
}

.condition_ridge_solve <- function(gram, rhs) {
    if (!length(rhs)) {
        return(numeric())
    }
    jitter <- max(1e-10, sqrt(.Machine$double.eps))
    diag(gram) <- diag(gram) + jitter
    answer <- tryCatch(
        as.numeric(solve(gram, rhs)),
        error = function(error) {
            as.numeric(qr.solve(as.matrix(gram), as.numeric(rhs)))
        }
    )
    if (any(!is.finite(answer))) {
        stop('Hierarchical refit produced non-finite coefficients.')
    }
    answer
}

.condition_refit_shared_baseline <- function(
    X_list,
    y_list,
    beta_selection,
    eligibility_mask,
    ridge,
    condition_weight = c('equal', 'cell_count'),
    max_iter = 200L,
    tol = 1e-8
) {
    condition_weight <- match.arg(condition_weight)
    beta_selection <- as.matrix(beta_selection)
    eligibility_mask <- as.matrix(eligibility_mask)
    if (!identical(dim(beta_selection), dim(eligibility_mask))) {
        stop('beta_selection and eligibility_mask must have identical dimensions.')
    }
    if (!is.logical(eligibility_mask) || anyNA(eligibility_mask)) {
        stop('eligibility_mask must be a logical matrix without NA values.')
    }
    if (!is.numeric(ridge) || length(ridge) != 1L ||
        !is.finite(ridge) || ridge <= 0) {
        stop('ridge must be one finite positive value.')
    }

    p <- nrow(beta_selection)
    n_tasks <- ncol(beta_selection)
    active_mask <- eligibility_mask &
        abs(beta_selection) > sqrt(.Machine$double.eps)
    beta <- matrix(
        0,
        nrow = p,
        ncol = n_tasks,
        dimnames = dimnames(beta_selection)
    )
    beta[active_mask] <- beta_selection[active_mask]
    average_weights <- .condition_average_weights(X_list, condition_weight)
    loss_weights <- .condition_loss_weights(X_list, condition_weight)
    shared <- .condition_weighted_row_mean(
        beta, eligibility_mask, average_weights
    )
    intercept <- numeric(n_tasks)
    converged <- FALSE
    coef_change <- Inf

    for (iteration in seq_len(max_iter)) {
        beta_previous <- beta
        shared_previous <- shared

        for (task in seq_len(n_tasks)) {
            active <- active_mask[, task]
            beta[, task] <- 0
            y_task <- y_list[[task]]
            if (!any(active)) {
                intercept[[task]] <- mean(y_task)
                next
            }

            X_task <- X_list[[task]][, active, drop = FALSE]
            x_mean <- as.numeric(Matrix::colMeans(X_task))
            y_mean <- mean(y_task)
            X_centered <- sweep(as.matrix(X_task), 2L, x_mean, '-')
            y_centered <- y_task - y_mean
            gram <- loss_weights[[task]] * crossprod(X_centered) +
                ridge * diag(sum(active))
            rhs <- loss_weights[[task]] *
                as.numeric(crossprod(X_centered, y_centered)) +
                ridge * shared[active]
            coefficient <- .condition_ridge_solve(gram, rhs)
            beta[active, task] <- coefficient
            intercept[[task]] <- y_mean - sum(x_mean * coefficient)
        }

        shared <- .condition_weighted_row_mean(
            beta, eligibility_mask, average_weights
        )
        numerator <- sqrt(sum((beta - beta_previous)^2)) +
            sqrt(sum((shared - shared_previous)^2))
        denominator <- sqrt(sum(beta_previous^2)) +
            sqrt(sum(shared_previous^2)) +
            .Machine$double.eps
        coef_change <- numerator / denominator
        if (coef_change < tol) {
            converged <- TRUE
            break
        }
    }

    beta_condition <- beta
    beta_condition[!eligibility_mask] <- NA_real_
    deviation <- sweep(beta_condition, 1L, shared, '-')
    list(
        beta = beta,
        beta_condition = beta_condition,
        beta_shared = shared,
        delta_condition = deviation,
        active_mask = active_mask,
        estimability_mask = eligibility_mask,
        intercept = intercept,
        ridge = ridge,
        iterations = iteration,
        coef_change = coef_change,
        converged = converged
    )
}

.condition_rebuild_target_design <- function(
    result, gene_data, peak_data, condition
) {
    contract <- result$fit_contract
    edges <- contract$edge_table
    tf_matrix <- gene_data[
        , match(edges$tf, colnames(gene_data)), drop = FALSE
    ]
    peak_matrix <- peak_data[
        , match(edges$region, colnames(peak_data)), drop = FALSE
    ]
    X <- tf_matrix * peak_matrix
    center <- as.numeric(contract$predictor_center[edges$edge_id])
    scale <- as.numeric(contract$predictor_scale[edges$edge_id])
    if (any(!is.finite(scale)) || any(scale <= 0)) {
        stop('Stored predictor scales are invalid for hierarchical refitting.')
    }
    X <- sweep(as.matrix(X), 2L, center, '-')
    X <- sweep(X, 2L, scale, '/')
    X <- Matrix::Matrix(X, sparse = FALSE)
    colnames(X) <- edges$term

    y <- as.numeric(gene_data[, contract$target])
    y <- (y - contract$response_center) / contract$response_scale
    condition_levels <- levels(condition)
    X_list <- lapply(condition_levels, function(level) {
        X[condition == level, , drop = FALSE]
    })
    y_list <- lapply(condition_levels, function(level) {
        y[condition == level]
    })
    names(X_list) <- names(y_list) <- condition_levels
    list(X = X, y = y, X_list = X_list, y_list = y_list)
}

.condition_apply_shared_baseline_refit <- function(
    result,
    gene_data,
    peak_data,
    condition,
    alpha,
    condition_weight
) {
    contract <- result$fit_contract
    design <- .condition_rebuild_target_design(
        result, gene_data, peak_data, condition
    )
    ridge <- max(
        contract$selected_lambda * (1 - alpha),
        1e-6
    )
    refit <- .condition_refit_shared_baseline(
        X_list = design$X_list,
        y_list = design$y_list,
        beta_selection = contract$beta,
        eligibility_mask = contract$eligibility_mask,
        ridge = ridge,
        condition_weight = condition_weight
    )

    condition_levels <- colnames(refit$beta)
    beta_selection <- contract$beta
    reference_condition <- contract$reference_condition
    reference_beta <- refit$beta[, reference_condition]
    contrast <- sweep(refit$beta, 1L, reference_beta, '-')

    result$universal_coefs$estimate <- refit$beta_shared
    for (level in condition_levels) {
        result$condition_coefs[[level]]$estimate <- refit$beta[, level]
        prediction <- refit$intercept[[level]] +
            as.numeric(design$X_list[[level]] %*% refit$beta[, level])
        result$condition_gof[[level]]$rsq <-
            .condition_rsq(design$y_list[[level]], prediction)
    }

    universal_rsq <- vapply(condition_levels, function(level) {
        universal_intercept <- mean(
            design$y_list[[level]] -
                as.numeric(design$X_list[[level]] %*% refit$beta_shared)
        )
        prediction <- universal_intercept +
            as.numeric(design$X_list[[level]] %*% refit$beta_shared)
        .condition_rsq(design$y_list[[level]], prediction)
    }, numeric(1))
    result$universal_gof$rsq <- if (any(is.finite(universal_rsq))) {
        mean(universal_rsq, na.rm = TRUE)
    } else {
        NA_real_
    }

    contract$beta_selection <- beta_selection
    contract$beta <- refit$beta
    contract$beta_condition <- refit$beta_condition
    contract$beta_shared <- refit$beta_shared
    contract$delta_condition <- refit$delta_condition
    contract$active_mask <- refit$active_mask
    contract$estimability_mask <- refit$estimability_mask
    contract$reference_beta <- reference_beta
    contract$contrast <- contrast
    contract$intercept <- stats::setNames(
        refit$intercept, condition_levels
    )
    contract$condition_rsq <- vapply(
        result$condition_gof, function(x) x$rsq[[1L]], numeric(1)
    )
    contract$fit_engine <-
        'shared_support_condition_sparse_hierarchical_refit'
    contract$refit <- list(
        method = 'support_constrained_hierarchical_ridge',
        ridge = refit$ridge,
        converged = refit$converged,
        iterations = refit$iterations,
        coef_change = refit$coef_change,
        support_source = 'sparse_group_selection',
        inactive_semantics = 'estimable_zero',
        unavailable_semantics = 'NA_in_beta_condition_and_masks'
    )
    result$fit_contract <- contract
    result$diagnostics$converged <-
        result$diagnostics$converged & refit$converged
    result$diagnostics$coef_change <- refit$coef_change
    result
}

.condition_fit_target <- function(..., method) {
    arguments <- list(...)
    arguments$method <- method
    result <- do.call(.condition_fit_target_v2_patched, arguments)
    if (!identical(method, 'shared_baseline_condition_sparse')) {
        return(result)
    }
    .condition_apply_shared_baseline_refit(
        result = result,
        gene_data = arguments$gene_data,
        peak_data = arguments$peak_data,
        condition = arguments$condition,
        alpha = arguments$alpha,
        condition_weight = arguments$condition_weight
    )
}

.condition_validate_public_args <- function(
    object,
    cell_type_col,
    condition_col,
    aggregate_rna_col,
    aggregate_peaks_col,
    method,
    family,
    interaction_term,
    alpha,
    condition_mix,
    reference_condition,
    condition_weight,
    nlambda,
    lambda,
    nfolds,
    min_cells_per_condition,
    active_tol,
    max_iter,
    tol_objective,
    tol_coef
) {
    validation_method <- if (
        identical(method, 'shared_baseline_condition_sparse')
    ) {
        'multitask_glmnet'
    } else {
        method
    }
    .condition_validate_public_args_v2(
        object = object,
        cell_type_col = cell_type_col,
        condition_col = condition_col,
        aggregate_rna_col = aggregate_rna_col,
        aggregate_peaks_col = aggregate_peaks_col,
        method = validation_method,
        family = family,
        interaction_term = interaction_term,
        alpha = alpha,
        condition_mix = condition_mix,
        reference_condition = reference_condition,
        condition_weight = condition_weight,
        nlambda = nlambda,
        lambda = lambda,
        nfolds = nfolds,
        min_cells_per_condition = min_cells_per_condition,
        active_tol = active_tol,
        max_iter = max_iter,
        tol_objective = tol_objective,
        tol_coef = tol_coef
    )
    if (identical(method, 'shared_baseline_condition_sparse') &&
        condition_mix <= 0) {
        stop(
            'shared_baseline_condition_sparse requires condition_mix > 0 ',
            'to permit condition-specific edge support.'
        )
    }
    invisible(TRUE)
}

#' @rdname infer_condition_grn
#' @method infer_condition_grn GRNData
#' @export
infer_condition_grn.GRNData <- function(
    object,
    cell_type_col,
    condition_col,
    genes = NULL,
    network_name = 'condition_grn',
    peak_to_gene_method = c('Signac', 'GREAT'),
    upstream = 100000,
    downstream = 0,
    extend = 1000000,
    only_tss = FALSE,
    peak_to_gene_domains = NULL,
    tf_cor = 0.1,
    peak_cor = 0,
    candidate_screen = c(
        'pooled_within_condition', 'pooled', 'motif_domain',
        'condition_union'
    ),
    aggregate_rna_col = NULL,
    aggregate_peaks_col = NULL,
    method = c(
        'shared_baseline_condition_sparse',
        'shared_design_independent',
        'multitask_glmnet'
    ),
    family = 'gaussian',
    interaction_term = ':',
    alpha = 0.5,
    condition_mix = 0.5,
    reference_condition = NULL,
    condition_weight = c('equal', 'cell_count'),
    nlambda = 50L,
    lambda = NULL,
    lambda_min_ratio = NULL,
    nfolds = 5L,
    lambda_selection = c('lambda.1se', 'lambda.min'),
    min_cells_per_condition = 50L,
    on_small_condition = c('skip_cell_type', 'drop_condition', 'error'),
    scale = TRUE,
    active_tol = 1e-8,
    parallel = FALSE,
    BPPARAM = NULL,
    overwrite = FALSE,
    seed = 12345L,
    max_iter = 5000L,
    tol_objective = 1e-7,
    tol_coef = 1e-6,
    verbose = TRUE,
    ...
) {
    peak_to_gene_method <- match.arg(peak_to_gene_method)
    candidate_screen <- match.arg(candidate_screen)
    method <- match.arg(method)
    condition_weight <- match.arg(condition_weight)
    lambda_selection <- match.arg(lambda_selection)
    on_small_condition <- match.arg(on_small_condition)

    .condition_validate_public_args(
        object, cell_type_col, condition_col, aggregate_rna_col,
        aggregate_peaks_col, method, family, interaction_term, alpha,
        condition_mix, reference_condition, condition_weight, nlambda,
        lambda, nfolds, min_cells_per_condition, active_tol, max_iter,
        tol_objective, tol_coef
    )
    prepared <- .condition_prepare_global_data(
        object = object,
        cell_type_col = cell_type_col,
        condition_col = condition_col,
        genes = genes,
        peak_to_gene_method = peak_to_gene_method,
        upstream = upstream,
        downstream = downstream,
        extend = extend,
        only_tss = only_tss,
        peak_to_gene_domains = peak_to_gene_domains,
        aggregate_col = aggregate_rna_col,
        verbose = verbose
    )
    .condition_run_all_cell_types(
        object = object,
        prepared = prepared,
        cell_type_col = cell_type_col,
        condition_col = condition_col,
        network_name = network_name,
        peak_to_gene_method = peak_to_gene_method,
        upstream = upstream,
        downstream = downstream,
        extend = extend,
        only_tss = only_tss,
        tf_cor = tf_cor,
        peak_cor = peak_cor,
        candidate_screen = candidate_screen,
        method = method,
        alpha = alpha,
        condition_mix = condition_mix,
        reference_condition = reference_condition,
        condition_weight = condition_weight,
        nlambda = nlambda,
        lambda = lambda,
        lambda_min_ratio = lambda_min_ratio,
        nfolds = nfolds,
        lambda_selection = lambda_selection,
        min_cells_per_condition = min_cells_per_condition,
        on_small_condition = on_small_condition,
        scale = scale,
        active_tol = active_tol,
        parallel = parallel,
        BPPARAM = BPPARAM,
        overwrite = overwrite,
        seed = seed,
        max_iter = max_iter,
        tol_objective = tol_objective,
        tol_coef = tol_coef,
        verbose = verbose
    )
}

.condition_combine_fit_contracts <- function(...) {
    arguments <- list(...)
    fit <- do.call(.condition_combine_fit_contracts_v2, arguments)
    contracts <- lapply(arguments$successful, `[[`, 'fit_contract')
    if (!all(vapply(
        contracts,
        function(contract) !is.null(contract$beta_shared),
        logical(1)
    ))) {
        return(fit)
    }

    edge_ids <- fit$edge_table$edge_id
    bind_matrix <- function(field) {
        value <- do.call(rbind, lapply(contracts, `[[`, field))
        rownames(value) <- edge_ids
        value
    }
    beta_shared <- unlist(
        lapply(contracts, `[[`, 'beta_shared'),
        use.names = FALSE
    )
    names(beta_shared) <- edge_ids

    fit$schema_version <- 'pando_condition_grn_fit_v3'
    fit$fit_engine <-
        'shared_support_condition_sparse_hierarchical_refit'
    fit$beta_selection <- bind_matrix('beta_selection')
    fit$beta_condition <- bind_matrix('beta_condition')
    fit$beta_shared <- beta_shared
    fit$delta_condition <- bind_matrix('delta_condition')
    fit$active_mask <- bind_matrix('active_mask')
    fit$estimability_mask <- bind_matrix('estimability_mask')
    fit$coefficient_scale <-
        'pooled_cell_type_standardized_support_constrained_refit'
    fit$universal_summary <-
        'explicit estimability-weighted shared baseline from refitted condition effects'
    fit$condition_subgraph_definition <-
        'active_mask on the shared TF-peak-target candidate supergraph'
    fit$direction_semantics <- list(
        absolute = 'sign(beta_condition)',
        deviation = 'sign(delta_condition)',
        pairwise = 'sign(beta_condition_2 - beta_condition_1)'
    )
    fit$refit <- stats::setNames(
        lapply(contracts, `[[`, 'refit'),
        vapply(contracts, `[[`, character(1), 'target')
    )
    fit
}

.condition_run_cell_type <- function(..., method) {
    result <- .condition_run_cell_type_v2(..., method = method)
    if (!identical(method, 'shared_baseline_condition_sparse') ||
        is.null(result$index)) {
        return(result)
    }
    fit_engine <- 'shared_support_condition_sparse_hierarchical_refit'
    result$index$fit_engine <- fit_engine
    for (network_id in result$index$network_id) {
        result$object@grn@networks[[network_id]]@params$fit_engine <-
            fit_engine
    }
    contract_keys <- unique(result$index$fit_contract_key)
    for (key in contract_keys) {
        result$object@grn@params$condition_grn_fits[[key]]$fit_engine <-
            fit_engine
    }
    result
}
