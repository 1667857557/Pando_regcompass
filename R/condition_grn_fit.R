# Target-level fitting for condition-aware Pando.

.condition_target_skip <- function(message) {
    stop(structure(
        list(message = message, call = NULL),
        class = c('condition_target_skip', 'error', 'condition')
    ))
}

.condition_fit_targets <- function(
    features,
    gene_data,
    peak_data,
    condition,
    peaks2gene,
    peaks2motif,
    motif2tf,
    candidate_screen,
    tf_cor,
    peak_cor,
    scale,
    alpha,
    condition_mix,
    active_tol,
    reference_condition,
    comparison_conditions,
    condition_weight,
    nlambda,
    lambda,
    lambda_min_ratio,
    outer_nfolds,
    inner_nfolds,
    lambda_selection,
    seed,
    max_iter,
    tol_objective,
    tol_coef,
    parallel,
    BPPARAM,
    verbose
) {
    names(features) <- features
    condition_index <- split(seq_along(condition), condition, drop = TRUE)
    fit_one <- function(gene) {
        tryCatch(
            .condition_fit_target(
                target = gene,
                gene_data = gene_data,
                peak_data = peak_data,
                condition = condition,
                condition_index = condition_index,
                peaks2gene = peaks2gene,
                peaks2motif = peaks2motif,
                motif2tf = motif2tf,
                candidate_screen = candidate_screen,
                tf_cor = tf_cor,
                peak_cor = peak_cor,
                scale = scale,
                alpha = alpha,
                condition_mix = condition_mix,
                active_tol = active_tol,
                reference_condition = reference_condition,
                comparison_conditions = comparison_conditions,
                condition_weight = condition_weight,
                nlambda = nlambda,
                lambda = lambda,
                lambda_min_ratio = lambda_min_ratio,
                outer_nfolds = outer_nfolds,
                inner_nfolds = inner_nfolds,
                lambda_selection = lambda_selection,
                seed = .condition_seed_for(gene, seed),
                max_iter = max_iter,
                tol_objective = tol_objective,
                tol_coef = tol_coef
            ),
            condition_target_skip = function(error) {
                list(
                    error = conditionMessage(error),
                    diagnostics = data.frame(
                        target = gene,
                        stage = 'target_skip',
                        converged = FALSE,
                        iterations = NA_integer_,
                        objective = NA_real_,
                        coef_change = NA_real_,
                        selected_lambda = NA_real_,
                        cv_mean = NA_real_,
                        cv_se = NA_real_,
                        error_message = conditionMessage(error),
                        stringsAsFactors = FALSE
                    )
                )
            }
        )
    }
    results <- .condition_map(
        features, fit_one, parallel = parallel, BPPARAM = BPPARAM, verbose = verbose
    )
    success <- vapply(results, function(x) is.null(x$error), logical(1))
    fits <- results[success]
    diagnostics <- do.call(rbind, lapply(results, function(x) x$diagnostics))
    rownames(diagnostics) <- NULL
    list(fits = fits, diagnostics = diagnostics)
}

.condition_fit_target <- function(
    target,
    gene_data,
    peak_data,
    condition,
    peaks2gene,
    peaks2motif,
    motif2tf,
    candidate_screen,
    tf_cor,
    peak_cor,
    scale,
    alpha,
    condition_mix,
    active_tol,
    reference_condition,
    comparison_conditions,
    condition_weight,
    nlambda,
    lambda,
    lambda_min_ratio,
    outer_nfolds,
    inner_nfolds,
    lambda_selection,
    seed,
    max_iter,
    tol_objective,
    tol_coef,
    condition_index = NULL
) {
    if (!target %in% rownames(peaks2gene)) {
        .condition_target_skip(
            'Target was not found in the peak-to-gene domain matrix.'
        )
    }
    gene_peak_flag <- as.logical(peaks2gene[target, ])
    candidate_peaks <- colnames(peaks2gene)[gene_peak_flag]
    if (length(candidate_peaks) == 0L) {
        .condition_target_skip('No candidate regulatory peaks were found.')
    }
    response_raw <- gene_data[, target, drop = FALSE]
    condition_levels <- levels(condition)
    if (is.null(condition_index)) {
        condition_index <- split(seq_along(condition), condition, drop = TRUE)
    }
    peak_raw <- peak_data[, candidate_peaks, drop = FALSE]
    peak_mask <- .condition_component_masks(
        peak_raw, response_raw, condition, peak_cor, candidate_screen
    )
    peak_motifs <- peaks2motif[candidate_peaks, , drop = FALSE]
    peak_tfs <- lapply(candidate_peaks, function(peak) {
        motif_flag <- as.logical(peak_motifs[peak, ])
        if (!any(motif_flag)) {
            return(character())
        }
        tf_flag <- as.logical(
            sparseMatrixStats::colMaxs(
                motif2tf[motif_flag, , drop = FALSE]
            ) > 0
        )
        setdiff(colnames(motif2tf)[tf_flag], target)
    })
    names(peak_tfs) <- candidate_peaks
    gene_tfs <- unique(unlist(peak_tfs, use.names = FALSE))
    if (length(gene_tfs) == 0L) {
        .condition_target_skip(
            'No motif-mapped TFs were found for candidate peaks.'
        )
    }
    tf_raw <- gene_data[, gene_tfs, drop = FALSE]
    tf_mask <- .condition_component_masks(
        tf_raw, response_raw, condition, tf_cor, candidate_screen
    )
    edge_parts <- lapply(names(peak_tfs), function(peak) {
        tfs <- peak_tfs[[peak]]
        if (length(tfs) == 0L) {
            return(NULL)
        }
        data.frame(
            tf = tfs,
            target = target,
            region = peak,
            stringsAsFactors = FALSE
        )
    })
    edge_parts <- edge_parts[!vapply(edge_parts, is.null, logical(1))]
    if (length(edge_parts) == 0L) {
        .condition_target_skip(
            'No TF-peak-target edges remained after screening.'
        )
    }
    edges <- unique(do.call(rbind, edge_parts))
    rownames(edges) <- NULL
    edges$edge_id <- paste(edges$tf, edges$region, edges$target, sep = '\001')
    edge_mask <- .condition_edge_mask(edges, peak_mask, tf_mask)
    edges <- edges[rowSums(edge_mask) > 0L, , drop = FALSE]
    edge_mask <- edge_mask[edges$edge_id, , drop = FALSE]
    if (nrow(edges) == 0L) {
        .condition_target_skip(
            'No TF-peak-target edge passed shared candidate screening.'
        )
    }
    tf_variance <- .condition_column_variance(
        gene_data[, unique(edges$tf), drop = FALSE]
    )
    names(tf_variance) <- unique(edges$tf)
    peak_variance <- .condition_column_variance(
        peak_data[, unique(edges$region), drop = FALSE]
    )
    names(peak_variance) <- unique(edges$region)
    valid_tfs <- names(tf_variance)[
        is.finite(tf_variance) & tf_variance > .Machine$double.eps
    ]
    valid_peaks <- names(peak_variance)[
        is.finite(peak_variance) & peak_variance > .Machine$double.eps
    ]
    edges <- edges[
        edges$tf %in% valid_tfs & edges$region %in% valid_peaks,
        , drop = FALSE
    ]
    edge_mask <- edge_mask[edges$edge_id, , drop = FALSE]
    if (nrow(edges) == 0L) {
        .condition_target_skip(
            'No edges remained after TF and peak variance checks.'
        )
    }
    prepared_design <- .condition_build_design(
        response_raw = response_raw,
        gene_data = gene_data,
        peak_data = peak_data,
        edges = edges,
        condition = condition,
        scale = scale,
        condition_index = condition_index
    )
    edges <- prepared_design$edges
    X <- prepared_design$X
    y <- prepared_design$y
    X_raw <- prepared_design$X_raw
    y_raw <- prepared_design$y_raw
    edge_mask <- edge_mask[edges$edge_id, , drop = FALSE]
    screening_mask <- edge_mask
    if (ncol(X) == 0L) {
        .condition_target_skip(
            'No non-constant interaction predictors remained.'
        )
    }
    for (level in condition_levels) {
        index <- condition_index[[level]]
        variance <- .condition_population_variance(
            X[index, , drop = FALSE]
        )
        edge_mask[, level] <- edge_mask[, level] &
            is.finite(variance) & variance > .Machine$double.eps
    }
    keep_edge <- rowSums(edge_mask) > 0L
    X <- X[, keep_edge, drop = FALSE]
    X_raw <- X_raw[, keep_edge, drop = FALSE]
    edges <- edges[keep_edge, , drop = FALSE]
    edge_mask <- edge_mask[keep_edge, , drop = FALSE]
    screening_mask <- screening_mask[keep_edge, , drop = FALSE]
    prepared_design$predictor_center <-
        prepared_design$predictor_center[keep_edge]
    prepared_design$predictor_scale <-
        prepared_design$predictor_scale[keep_edge]
    if (!ncol(X)) {
        .condition_target_skip(
            'No edge has a non-constant predictor in an eligible condition.'
        )
    }
    if (is.null(reference_condition)) {
        reference_condition <- condition_levels[[1L]]
    }
    reference_condition <- as.character(reference_condition)
    if (!reference_condition %in% condition_levels) {
        stop(
            'reference_condition was not found for target ', target, ': ',
            reference_condition, '.'
        )
    }
    comparison_conditions <- if (is.null(comparison_conditions)) {
        condition_levels
    } else {
        unique(as.character(comparison_conditions))
    }
    if (length(comparison_conditions) < 2L ||
        !all(comparison_conditions %in% condition_levels)) {
        stop(
            'comparison_conditions were not found for target ', target, '.'
        )
    }
    X_list <- lapply(condition_levels, function(level) {
        X[condition_index[[level]], , drop = FALSE]
    })
    y_list <- lapply(condition_levels, function(level) {
        y[condition_index[[level]]]
    })
    X_raw_list <- lapply(condition_levels, function(level) {
        X_raw[condition_index[[level]], , drop = FALSE]
    })
    y_raw_list <- lapply(condition_levels, function(level) {
        y_raw[condition_index[[level]]]
    })
    cell_id_list <- lapply(condition_levels, function(level) {
        rownames(gene_data)[condition_index[[level]]]
    })
    names(X_list) <- names(y_list) <- condition_levels
    names(X_raw_list) <- names(y_raw_list) <- condition_levels
    names(cell_id_list) <- condition_levels
    engine <- .condition_fit_target_native_engine(
        X_raw_list = X_raw_list,
        y_raw_list = y_raw_list,
        coefficient_mask = edge_mask,
        comparison_conditions = comparison_conditions,
        lambda = lambda,
        nlambda = nlambda,
        lambda_min_ratio = lambda_min_ratio,
        alpha = alpha,
        condition_mix = condition_mix,
        active_tol = active_tol,
        outer_nfolds = outer_nfolds,
        inner_nfolds = inner_nfolds,
        lambda_selection = lambda_selection,
        seed = seed,
        max_iter = max_iter,
        tol_objective = tol_objective,
        tol_coef = tol_coef
    )
    lambda_path <- engine$lambda_path
    cv <- engine$cv
    full_cv <- engine$full_cv
    selected <- engine$selected
    beta_selection <- engine$beta_selection
    refit <- engine$refit
    native_transform <- engine$full_transform
    transform_checks <- list(
        predictor_center = all.equal(
            unname(native_transform$predictor_center),
            unname(prepared_design$predictor_center),
            tolerance = 1e-12,
            check.attributes = FALSE
        ),
        predictor_scale = all.equal(
            unname(native_transform$predictor_scale),
            unname(prepared_design$predictor_scale),
            tolerance = 1e-12,
            check.attributes = FALSE
        ),
        response_center = all.equal(
            native_transform$response_center,
            prepared_design$response_center,
            tolerance = 1e-12,
            check.attributes = FALSE
        ),
        response_scale = all.equal(
            native_transform$response_scale,
            prepared_design$response_scale,
            tolerance = 1e-12,
            check.attributes = FALSE
        )
    )
    if (!all(vapply(transform_checks, isTRUE, logical(1)))) {
        stop(
            'Compiled target engine did not preserve the canonical full-data ',
            'condition-balanced transform.',
            call. = FALSE
        )
    }
    beta <- refit$beta
    rownames(beta) <- edges$edge_id
    pooled_tf_corr <- .condition_named_cor(
        gene_data[, unique(edges$tf), drop = FALSE], response_raw
    )
    condition_tf_corr <- lapply(condition_levels, function(level) {
        index <- condition_index[[level]]
        .condition_named_cor(
            gene_data[index, unique(edges$tf), drop = FALSE],
            response_raw[index, , drop = FALSE]
        )
    })
    names(condition_tf_corr) <- condition_levels
    universal_beta <- refit$beta_shared
    universal_coefs <- .condition_format_coefs(
        edges, universal_beta, pooled_tf_corr[edges$tf]
    )
    condition_coefs <- lapply(condition_levels, function(level) {
        .condition_format_coefs(
            edges,
            refit$beta_condition[, level],
            condition_tf_corr[[level]][edges$tf]
        )
    })
    names(condition_coefs) <- condition_levels
    condition_gof <- stats::setNames(
        vector('list', length(condition_levels)), condition_levels
    )
    universal_rsq <- numeric(length(condition_levels))
    for (task in seq_along(condition_levels)) {
        level <- condition_levels[[task]]
        condition_prediction <- refit$intercept[[task]] + as.numeric(
            X_list[[task]] %*% beta[, task]
        )
        condition_rsq <- .condition_rsq(y_list[[task]], condition_prediction)
        condition_gof[[level]] <- data.frame(
            target = target,
            lambda = selected$lambda,
            rsq = condition_rsq,
            alpha = alpha,
            nvariables = nrow(edges),
            stringsAsFactors = FALSE
        )
        universal_linear <- as.numeric(X_list[[task]] %*% universal_beta)
        universal_intercept <- mean(y_list[[task]] - universal_linear)
        universal_prediction <- universal_intercept + universal_linear
        universal_rsq[[task]] <- .condition_rsq(
            y_list[[task]], universal_prediction
        )
    }
    universal_gof <- data.frame(
        target = target,
        lambda = selected$lambda,
        rsq = mean(universal_rsq, na.rm = TRUE),
        alpha = alpha,
        nvariables = nrow(edges),
        stringsAsFactors = FALSE
    )
    if (!any(is.finite(universal_rsq))) {
        universal_gof$rsq <- NA_real_
    }
    selected_cv_mean <- if (all(is.na(full_cv$cv_mean))) {
        NA_real_
    } else {
        full_cv$cv_mean[[full_cv$selected_index]]
    }
    selected_cv_se <- if (all(is.na(full_cv$cv_se))) {
        NA_real_
    } else {
        full_cv$cv_se[[full_cv$selected_index]]
    }
    diagnostics <- data.frame(
        target = target,
        stage = 'complete',
        converged = selected$converged && refit$converged,
        iterations = selected$iterations + refit$iterations,
        objective = selected$objective,
        coef_change = max(selected$coef_change, refit$coef_change),
        selected_lambda = selected$lambda,
        cv_mean = selected_cv_mean,
        cv_se = selected_cv_se,
        error_message = NA_character_,
        stringsAsFactors = FALSE
    )
    reference_beta <- refit$beta_condition[, reference_condition]
    contrast <- sweep(refit$beta_condition, 1L, reference_beta, '-')
    fit_engine <- 'condition_sparse_within_cell_type_oof_refit'
    fit_contract <- list(
        target = target,
        edge_table = edges[, c(
            'edge_id', 'tf', 'target', 'region', 'term'
        ), drop = FALSE],
        beta = beta,
        beta_selection = beta_selection,
        beta_condition = refit$beta_condition,
        beta_shared = refit$beta_shared,
        delta_condition = refit$delta_condition,
        reference_beta = reference_beta,
        contrast = contrast,
        eligibility_mask = edge_mask,
        estimability_mask = refit$estimability_mask,
        support_mask = refit$support_mask,
        active_mask = refit$active_mask,
        structural_candidate_mask = matrix(
            TRUE, nrow(edge_mask), ncol(edge_mask), dimnames = dimnames(edge_mask)
        ),
        screening_mask = screening_mask,
        predictor_center = native_transform$predictor_center,
        predictor_scale = native_transform$predictor_scale,
        response_center = native_transform$response_center,
        response_scale = native_transform$response_scale,
        transform_policy = native_transform$transform_policy,
        condition_transform_weights = native_transform$condition_weights,
        predictor_center_hash = native_transform$predictor_center_hash,
        predictor_scale_hash = native_transform$predictor_scale_hash,
        training_fold_only = TRUE,
        intercept = refit$intercept,
        condition_rsq = vapply(
            condition_gof, function(x) x$rsq[[1L]], numeric(1)
        ),
        condition_rsq_train = vapply(
            condition_gof, function(x) x$rsq[[1L]], numeric(1)
        ),
        condition_rsq_oof = engine$condition_rsq_oof,
        condition_rmse_oof = engine$condition_rmse_oof,
        target_rsq_oof_pooled = engine$target_rsq_oof_pooled,
        oof_fold = if (is.null(cv$oof_fold)) NULL else cv$oof_fold,
        cv_fold_transform = if (is.null(cv$fold_transform)) NULL else cv$fold_transform,
        outer_nfolds = cv$outer_nfolds,
        inner_nfolds = cv$inner_nfolds,
        cv_effective_nfolds = cv$outer_nfolds,
        cv_method = cv$cv_method,
        oof_model = cv$oof_model,
        projection_condition_full_oof = unlist(
            unname(cv$projection_condition_full_oof), use.names = TRUE
        ),
        projection_common_oof = unlist(
            unname(cv$projection_common_oof), use.names = TRUE
        ),
        projection_global_common_oof = unlist(
            unname(cv$projection_global_common_oof), use.names = TRUE
        ),
        projection_origin = cv$projection_origin,
        projection_used_for_penalty =
            cv$projection_used_for_penalty &&
            identical(candidate_screen, 'motif_domain'),
        full_fit_projection_used_for_penalty =
            cv$full_fit_projection_used_for_penalty,
        fold_transform_policy = cv$fold_transform_policy,
        oof_cell_coverage = mean(unlist(cv$oof_cell_coverage)),
        oof_projection_available_fraction =
            mean(unlist(cv$oof_projection_available_fraction)),
        oof_assignment_count = unlist(
            unname(cv$oof_assignment_count), use.names = TRUE
        ),
        outer_selected_lambda = cv$fold_selected_lambda,
        comparison_conditions = comparison_conditions,
        predictive_oof_available = identical(candidate_screen, 'motif_domain'),
        oof_validation_level = 'outer_condition_stratified_heldout_cells',
        selected_lambda = selected$lambda,
        lambda_path = lambda_path,
        cv_mean = full_cv$cv_mean,
        cv_se = full_cv$cv_se,
        alpha = alpha,
        condition_mix = condition_mix,
        condition_weight = condition_weight,
        active_tol = active_tol,
        fit_engine = fit_engine,
        reference_condition = reference_condition,
        refit = list(
            method = 'support_constrained_common_metric_ridge_direct_schur',
            numerical_backend = engine$backend,
            ridge = refit$ridge,
            common_metric = refit$common_metric,
            converged = refit$converged,
            iterations = refit$iterations,
            coef_change = refit$coef_change,
            support_source = 'condition_sparse_selection',
            inactive_semantics = 'estimable_zero',
            unavailable_semantics = 'NA'
        )
    )
    list(
        universal_coefs = universal_coefs,
        condition_coefs = condition_coefs,
        universal_gof = universal_gof,
        condition_gof = condition_gof,
        diagnostics = diagnostics,
        fit_contract = fit_contract
    )
}

.condition_build_design <- function(
    response_raw, gene_data, peak_data, edges, condition, scale,
    condition_index = NULL
) {
    y <- as.numeric(response_raw[, 1L])
    tf_names <- unique(edges$tf)
    peak_names <- unique(edges$region)
    tf_matrix <- methods::as(
        gene_data[, tf_names, drop = FALSE], 'dgCMatrix'
    )
    peak_matrix <- methods::as(
        peak_data[, peak_names, drop = FALSE], 'dgCMatrix'
    )
    tf_index <- match(edges$tf, colnames(tf_matrix))
    peak_index <- match(edges$region, colnames(peak_matrix))
    if (anyNA(tf_index) || anyNA(peak_index)) {
        stop('TF and peak edge indices are not aligned to the design matrices.')
    }
    X <- .condition_product_matrix_cpp(
        tf_matrix,
        peak_matrix,
        as.integer(tf_index),
        as.integer(peak_index)
    )
    dimnames(X) <- list(
        rownames(tf_matrix), colnames(tf_matrix)[tf_index]
    )
    edge_variance <- .condition_population_variance(X)
    keep <- is.finite(edge_variance) & edge_variance > .Machine$double.eps
    X <- X[, keep, drop = FALSE]
    edges <- edges[keep, , drop = FALSE]
    X_raw <- X
    y_raw <- y
    response_center <- 0
    response_scale <- 1
    predictor_center <- rep(0, ncol(X))
    predictor_scale <- rep(1, ncol(X))
    transform <- NULL
    if (scale) {
        levels_condition <- levels(condition)
        if (is.null(condition_index)) {
            condition_index <- split(seq_along(condition), condition, drop = TRUE)
        }
        X_list <- lapply(levels_condition, function(level) {
            X_raw[condition_index[[level]], , drop = FALSE]
        })
        y_list <- lapply(levels_condition, function(level) {
            y_raw[condition_index[[level]]]
        })
        names(X_list) <- names(y_list) <- levels_condition
        fold_stats <- .condition_build_fold_statistics(X_list, y_list)
        transform <- tryCatch(
            .condition_build_balanced_transform(
                X_list, y_list, fold_statistics = fold_stats
            ),
            error = function(error) {
                .condition_target_skip(conditionMessage(error))
            }
        )
        response_center <- transform$response_center
        response_scale <- transform$response_scale
        predictor_center <- transform$predictor_center
        predictor_scale <- transform$predictor_scale
        scaled <- .condition_apply_balanced_transform(
            list(all = X_raw), list(all = y_raw), transform
        )
        X <- scaled$X[[1L]]
        y <- scaled$y[[1L]]
    } else {
        response_variance <- stats::var(y)
        if (!is.finite(response_variance) ||
            response_variance <= .Machine$double.eps) {
            .condition_target_skip('Target response has zero variance.')
        }
    }
    edges$term <- paste(
        .condition_model_name(edges$region),
        .condition_model_name(edges$tf),
        sep = ':'
    )
    colnames(X) <- edges$term
    names(predictor_center) <- names(predictor_scale) <- edges$edge_id
    list(
        X = X,
        y = y,
        X_raw = X_raw,
        y_raw = y_raw,
        edges = edges,
        predictor_center = predictor_center,
        predictor_scale = predictor_scale,
        response_center = response_center,
        response_scale = response_scale,
        transform = transform
    )
}

.condition_scale_vector <- function(x) {
    standard_deviation <- stats::sd(x)
    if (!is.finite(standard_deviation) ||
        standard_deviation <= .Machine$double.eps) {
        return(NULL)
    }
    (x - mean(x)) / standard_deviation
}

.condition_scale_matrix <- function(x) {
    dense <- as.matrix(x)
    means <- colMeans(dense)
    standard_deviations <- apply(dense, 2, stats::sd)
    valid <- is.finite(standard_deviations) &
        standard_deviations > .Machine$double.eps
    dense <- dense[, valid, drop = FALSE]
    means <- means[valid]
    standard_deviations <- standard_deviations[valid]
    if (ncol(dense) == 0L) {
        return(Matrix::Matrix(matrix(numeric(), nrow(x), 0), sparse = FALSE))
    }
    dense <- sweep(dense, 2, means, '-')
    dense <- sweep(dense, 2, standard_deviations, '/')
    Matrix::Matrix(dense, sparse = FALSE)
}

.condition_screen_columns <- function(x, y, condition, threshold, candidate_screen) {
    rowSums(.condition_component_masks(
        x, y, condition, threshold, candidate_screen
    )) > 0L
}

.condition_edge_mask <- function(edges, peak_mask, tf_mask) {
    condition_levels <- colnames(peak_mask)
    if (!identical(condition_levels, colnames(tf_mask))) {
        stop('TF and peak screening masks must use identical condition levels.')
    }
    edge_mask <- vapply(condition_levels, function(level) {
        peak_mask[edges$region, level] & tf_mask[edges$tf, level]
    }, logical(nrow(edges)))
    if (is.null(dim(edge_mask))) {
        edge_mask <- matrix(edge_mask, nrow = nrow(edges))
    }
    dimnames(edge_mask) <- list(edges$edge_id, condition_levels)
    edge_mask
}

.condition_component_masks <- function(
    x, y, condition, threshold, candidate_screen
) {
    levels_condition <- levels(condition)
    if (candidate_screen == 'motif_domain') {
        return(matrix(
            TRUE, nrow = ncol(x), ncol = length(levels_condition),
            dimnames = list(colnames(x), levels_condition)
        ))
    }
    if (candidate_screen != 'pooled_within_condition') {
        stop('Unsupported candidate screening strategy.')
    }
    score <- .condition_within_association_score(x, y, condition)
    keep <- is.finite(score) & score > threshold
    matrix(
        keep,
        nrow = length(keep),
        ncol = length(levels_condition),
        dimnames = list(names(score), levels_condition)
    )
}

.condition_within_association_score <- function(x, y, condition) {
    if (ncol(x) == 0L) {
        return(stats::setNames(numeric(), character()))
    }
    condition <- droplevels(condition)
    squared <- numeric(ncol(x))
    n_contributing <- 0L
    for (level in levels(condition)) {
        keep <- condition == level
        if (sum(keep) <= 1L) {
            next
        }
        correlation <- .condition_named_cor(
            x[keep, , drop = FALSE], y[keep, , drop = FALSE]
        )
        squared <- squared + correlation^2
        n_contributing <- n_contributing + 1L
    }
    score <- if (n_contributing > 0L) {
        sqrt(squared / n_contributing)
    } else {
        rep(0, ncol(x))
    }
    score[!is.finite(score)] <- 0
    stats::setNames(score, colnames(x))
}

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
        ssy_total <- ssy_total + sum(y_level * y_level) -
            sum_y * sum_y / n_level
    }
    denominator <- sqrt(pmax(ssx_total, 0) * max(ssy_total, 0))
    correlation <- cross_total / denominator
    correlation[!is.finite(correlation)] <- 0
    stats::setNames(correlation, colnames(x))
}

.condition_named_cor <- function(x, y) {
    if (ncol(x) == 0L) {
        return(setNames(numeric(), character()))
    }
    correlation <- as.numeric(sparse_cor(x, y)[, 1L])
    correlation[!is.finite(correlation)] <- 0
    names(correlation) <- colnames(x)
    correlation
}
