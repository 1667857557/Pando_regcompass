# Target-level fitting for condition-aware Pando.

.condition_target_skip <- function(message) {
    stop(structure(
        list(message = message, call = NULL),
        class = c('condition_target_skip', 'error', 'condition')
    ))
}

.condition_fit_cell_type <- function(
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
    condition_weight,
    nlambda,
    lambda,
    lambda_min_ratio,
    nfolds,
    cv_block,
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
    fit_one <- function(gene) {
        tryCatch(
            .condition_fit_target(
                target = gene,
                gene_data = gene_data,
                peak_data = peak_data,
                condition = condition,
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
                condition_weight = condition_weight,
                nlambda = nlambda,
                lambda = lambda,
                lambda_min_ratio = lambda_min_ratio,
                nfolds = nfolds,
                cv_block = cv_block,
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
    condition_weight,
    nlambda,
    lambda,
    lambda_min_ratio,
    nfolds,
    cv_block,
    lambda_selection,
    seed,
    max_iter,
    tol_objective,
    tol_coef
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
            sparseMatrixStats::colMaxs(motif2tf[motif_flag, , drop = FALSE]) > 0
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
        scale = scale
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
        variance <- .condition_column_variance(
            X[condition == level, , drop = FALSE]
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
    X_list <- lapply(condition_levels, function(level) {
        X[condition == level, , drop = FALSE]
    })
    y_list <- lapply(condition_levels, function(level) y[condition == level])
    X_raw_list <- lapply(condition_levels, function(level) {
        X_raw[condition == level, , drop = FALSE]
    })
    y_raw_list <- lapply(condition_levels, function(level) {
        y_raw[condition == level]
    })
    cell_id_list <- lapply(condition_levels, function(level) {
        rownames(gene_data)[condition == level]
    })
    names(X_list) <- names(y_list) <- condition_levels
    names(X_raw_list) <- names(y_raw_list) <- condition_levels
    names(cell_id_list) <- condition_levels
    block_list <- if (is.null(cv_block)) {
        NULL
    } else {
        if (length(cv_block) != length(condition) || anyNA(cv_block)) {
            stop('cv_block must align to condition without missing values.')
        }
        lapply(condition_levels, function(level) {
            as.character(cv_block[condition == level])
        })
    }
    if (!is.null(block_list)) names(block_list) <- condition_levels
    sample_block_status <- .condition_sample_block_status(block_list)
    sample_blocked_oof_available <- sample_block_status$available
    cv_block_list_used <- if (sample_blocked_oof_available) {
        block_list
    } else {
        NULL
    }

    lambda_path <- if (is.null(lambda)) {
        .condition_make_lambda_path(
            X_list = X_list,
            y_list = y_list,
            alpha = alpha,
            condition_mix = condition_mix,
            condition_weight = condition_weight,
            coefficient_mask = edge_mask,
            nlambda = nlambda,
            lambda_min_ratio = lambda_min_ratio
        )
    } else {
        sort(unique(as.numeric(lambda)), decreasing = TRUE)
    }

    if (length(lambda_path) == 1L && is.null(cv_block_list_used)) {
        cv <- list(
            lambda = lambda_path,
            cv_mean = NA_real_,
            cv_se = NA_real_,
            selected_index = 1L,
            selected_lambda = lambda_path[[1L]],
            lambda_min = lambda_path[[1L]],
            lambda_1se = lambda_path[[1L]]
        )
    } else {
        cv <- .condition_cv_multitask_path(
            X_list = X_raw_list,
            y_list = y_raw_list,
            lambda = lambda_path,
            alpha = alpha,
            condition_mix = condition_mix,
            condition_weight = condition_weight,
            coefficient_mask = edge_mask,
            nfolds = nfolds,
            block_list = cv_block_list_used,
            standardize = scale,
            active_tol = active_tol,
            lambda_selection = lambda_selection,
            seed = seed,
            max_iter = max_iter,
            tol_objective = tol_objective,
            tol_coef = tol_coef
        )
    }
    if (!is.null(block_list) && !sample_blocked_oof_available) {
        cv$selection_oof_prediction <- cv$oof_prediction
        cv$selection_oof_fold <- cv$oof_fold
        cv$oof_prediction <- NULL
        cv$oof_fold <- NULL
        cv$block_to_fold <- NULL
        cv$fold_transform <- NULL
        cv$effective_nfolds <- NA_integer_
    }
    cv$cv_method <- if (sample_blocked_oof_available) {
        'biological_sample_blocked'
    } else if (!is.null(block_list)) {
        'cell_level_lambda_selection_no_sample_blocked_oof'
    } else {
        'cell_level'
    }
    if (!is.null(cv$oof_fold)) {
        cv$oof_fold <- lapply(seq_along(cv$oof_fold), function(task) {
            stats::setNames(
                as.integer(cv$oof_fold[[task]]), cell_id_list[[task]]
            )
        })
        names(cv$oof_fold) <- condition_levels
    }
    if (!is.null(cv$oof_prediction)) {
        cv$oof_prediction <- lapply(
            seq_along(cv$oof_prediction), function(task) {
                stats::setNames(
                    as.numeric(cv$oof_prediction[[task]]),
                    cell_id_list[[task]]
                )
            }
        )
        names(cv$oof_prediction) <- condition_levels
    }

    full_path <- .condition_fit_multitask_path(
        X_list = X_list,
        y_list = y_list,
        lambda = lambda_path,
        alpha = alpha,
        condition_mix = condition_mix,
        condition_weight = condition_weight,
        coefficient_mask = edge_mask,
        max_iter = max_iter,
        tol_objective = tol_objective,
        tol_coef = tol_coef
    )
    selected_index <- match(cv$selected_lambda, full_path$lambda)
    selected <- full_path$fits[[selected_index]]
    beta_selection <- selected$beta
    colnames(beta_selection) <- condition_levels
    rownames(beta_selection) <- edges$edge_id
    refit <- .condition_refit_shared_baseline(
        X_list = X_list,
        y_list = y_list,
        beta_selection = beta_selection,
        estimability_mask = edge_mask,
        ridge = max(selected$lambda * (1 - alpha), 1e-6),
        active_tol = active_tol,
        condition_weight = condition_weight
    )
    beta <- refit$beta
    rownames(beta) <- edges$edge_id

    pooled_tf_corr <- .condition_named_cor(
        gene_data[, unique(edges$tf), drop = FALSE], response_raw
    )
    condition_tf_corr <- lapply(condition_levels, function(level) {
        .condition_named_cor(
            gene_data[condition == level, unique(edges$tf), drop = FALSE],
            response_raw[condition == level, , drop = FALSE]
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
        universal_intercept <- mean(
            y_list[[task]] - as.numeric(X_list[[task]] %*% universal_beta)
        )
        universal_prediction <- universal_intercept + as.numeric(
            X_list[[task]] %*% universal_beta
        )
        universal_rsq[[task]] <- .condition_rsq(y_list[[task]], universal_prediction)
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

    selected_cv_mean <- if (all(is.na(cv$cv_mean))) {
        NA_real_
    } else {
        cv$cv_mean[[cv$selected_index]]
    }
    selected_cv_se <- if (all(is.na(cv$cv_se))) {
        NA_real_
    } else {
        cv$cv_se[[cv$selected_index]]
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
    fit_engine <- 'condition_sparse_common_scale_refit'
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
        predictor_center = stats::setNames(
            prepared_design$predictor_center, edges$edge_id
        ),
        predictor_scale = stats::setNames(
            prepared_design$predictor_scale, edges$edge_id
        ),
        response_center = prepared_design$response_center,
        response_scale = prepared_design$response_scale,
        intercept = refit$intercept,
        condition_rsq = vapply(
            condition_gof, function(x) x$rsq[[1L]], numeric(1)
        ),
        condition_rsq_train = vapply(
            condition_gof, function(x) x$rsq[[1L]], numeric(1)
        ),
        condition_rsq_oof = if (is.null(cv$oof_prediction)) {
            stats::setNames(rep(NA_real_, length(condition_levels)), condition_levels)
        } else {
            stats::setNames(vapply(seq_along(condition_levels), function(task) {
                .condition_rsq(
                    y_raw_list[[task]], cv$oof_prediction[[task]]
                )
            }, numeric(1)), condition_levels)
        },
        condition_rmse_oof = if (is.null(cv$oof_prediction)) {
            stats::setNames(rep(NA_real_, length(condition_levels)), condition_levels)
        } else {
            stats::setNames(vapply(seq_along(condition_levels), function(task) {
                sqrt(mean(
                    (y_raw_list[[task]] - cv$oof_prediction[[task]])^2
                ))
            }, numeric(1)), condition_levels)
        },
        target_rsq_oof_pooled = if (is.null(cv$oof_prediction)) {
            NA_real_
        } else {
            .condition_pooled_task_rsq(
                y_raw_list, cv$oof_prediction
            )
        },
        oof_fold = if (is.null(cv$oof_fold)) NULL else cv$oof_fold,
        cv_block_to_fold = if (is.null(cv$block_to_fold)) {
            NULL
        } else {
            cv$block_to_fold
        },
        cv_fold_transform = if (is.null(cv$fold_transform)) {
            NULL
        } else {
            cv$fold_transform
        },
        cv_effective_nfolds = if (is.null(cv$effective_nfolds)) {
            NA_integer_
        } else {
            cv$effective_nfolds
        },
        cv_method = cv$cv_method,
        oof_model = if (is.null(cv$oof_model)) {
            'fixed_lambda_without_cross_validation'
        } else {
            cv$oof_model
        },
        sample_blocked_oof_available = sample_blocked_oof_available,
        sample_block_status = sample_block_status,
        selected_lambda = selected$lambda,
        lambda_path = lambda_path,
        cv_mean = cv$cv_mean,
        cv_se = cv$cv_se,
        alpha = alpha,
        condition_mix = condition_mix,
        condition_weight = condition_weight,
        active_tol = active_tol,
        fit_engine = fit_engine,
        reference_condition = reference_condition,
        refit = list(
            method = 'support_constrained_common_metric_ridge',
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

.condition_build_design <- function(response_raw, gene_data, peak_data, edges, scale) {
    y <- as.numeric(response_raw[, 1L])
    tf_names <- unique(edges$tf)
    peak_names <- unique(edges$region)
    tf_matrix <- gene_data[, tf_names, drop = FALSE]
    peak_matrix <- peak_data[, peak_names, drop = FALSE]

    tf_edge <- tf_matrix[, match(edges$tf, colnames(tf_matrix)), drop = FALSE]
    peak_edge <- peak_matrix[, match(edges$region, colnames(peak_matrix)), drop = FALSE]
    X <- tf_edge * peak_edge
    edge_variance <- .condition_column_variance(X)
    keep <- is.finite(edge_variance) & edge_variance > .Machine$double.eps
    X <- X[, keep, drop = FALSE]
    edges <- edges[keep, , drop = FALSE]
    X_raw <- X
    y_raw <- y
    response_center <- 0
    response_scale <- 1
    predictor_center <- rep(0, ncol(X))
    predictor_scale <- rep(1, ncol(X))
    if (scale) {
        response_center <- mean(y)
        response_scale <- stats::sd(y)
        if (!is.finite(response_scale) ||
            response_scale <= .Machine$double.eps) {
            .condition_target_skip('Target response has zero variance.')
        }
        y <- (y - response_center) / response_scale
        predictor_center <- as.numeric(Matrix::colMeans(X))
        predictor_scale <- sqrt(.condition_column_variance(X))
        X <- sweep(as.matrix(X), 2L, predictor_center, '-')
        X <- sweep(X, 2L, predictor_scale, '/')
        X <- Matrix::Matrix(X, sparse = FALSE)
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
        response_scale = response_scale
    )
}

.condition_scale_vector <- function(x) {
    standard_deviation <- stats::sd(x)
    if (!is.finite(standard_deviation) || standard_deviation <= .Machine$double.eps) {
        return(NULL)
    }
    (x - mean(x)) / standard_deviation
}

.condition_scale_matrix <- function(x) {
    dense <- as.matrix(x)
    means <- colMeans(dense)
    standard_deviations <- apply(dense, 2, stats::sd)
    valid <- is.finite(standard_deviations) & standard_deviations > .Machine$double.eps
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
        out <- matrix(
            TRUE, nrow = ncol(x), ncol = length(levels_condition),
            dimnames = list(colnames(x), levels_condition)
        )
        return(out)
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
    total_weight <- 0
    for (level in levels(condition)) {
        keep <- condition == level
        weight <- max(sum(keep) - 1L, 0L)
        if (weight == 0L) {
            next
        }
        correlation <- .condition_named_cor(
            x[keep, , drop = FALSE], y[keep, , drop = FALSE]
        )
        squared <- squared + weight * correlation^2
        total_weight <- total_weight + weight
    }
    score <- if (total_weight > 0) {
        sqrt(squared / total_weight)
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
