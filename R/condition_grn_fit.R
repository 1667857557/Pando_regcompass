# Target-level fitting for condition-aware Pando.

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
    condition_weight,
    nlambda,
    lambda,
    lambda_min_ratio,
    nfolds,
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
                condition_weight = condition_weight,
                nlambda = nlambda,
                lambda = lambda,
                lambda_min_ratio = lambda_min_ratio,
                nfolds = nfolds,
                lambda_selection = lambda_selection,
                seed = .condition_seed_for(gene, seed),
                max_iter = max_iter,
                tol_objective = tol_objective,
                tol_coef = tol_coef
            ),
            error = function(error) {
                list(
                    error = conditionMessage(error),
                    diagnostics = data.frame(
                        target = gene,
                        stage = 'target_fit',
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
    condition_weight,
    nlambda,
    lambda,
    lambda_min_ratio,
    nfolds,
    lambda_selection,
    seed,
    max_iter,
    tol_objective,
    tol_coef
) {
    if (!target %in% rownames(peaks2gene)) {
        stop('Target was not found in the peak-to-gene domain matrix.')
    }
    gene_peak_flag <- as.logical(peaks2gene[target, ])
    candidate_peaks <- colnames(peaks2gene)[gene_peak_flag]
    if (length(candidate_peaks) == 0L) {
        stop('No candidate regulatory peaks were found.')
    }

    response_raw <- gene_data[, target, drop = FALSE]
    peak_raw <- peak_data[, candidate_peaks, drop = FALSE]
    peak_keep <- .condition_screen_columns(
        peak_raw, response_raw, condition, peak_cor, candidate_screen
    )
    candidate_peaks <- candidate_peaks[peak_keep]
    if (length(candidate_peaks) == 0L) {
        stop('No peaks passed candidate screening.')
    }

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
        stop('No motif-mapped TFs were found for candidate peaks.')
    }

    tf_raw <- gene_data[, gene_tfs, drop = FALSE]
    tf_keep <- .condition_screen_columns(
        tf_raw, response_raw, condition, tf_cor, candidate_screen
    )
    retained_tfs <- gene_tfs[tf_keep]
    if (length(retained_tfs) == 0L) {
        stop('No TFs passed candidate screening.')
    }

    edge_parts <- lapply(names(peak_tfs), function(peak) {
        tfs <- intersect(peak_tfs[[peak]], retained_tfs)
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
        stop('No TF-peak-target edges remained after screening.')
    }
    edges <- unique(do.call(rbind, edge_parts))
    rownames(edges) <- NULL

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
    if (nrow(edges) == 0L) {
        stop('No edges remained after TF and peak variance checks.')
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
    if (ncol(X) == 0L) {
        stop('No non-constant interaction predictors remained.')
    }

    condition_levels <- levels(condition)
    X_list <- lapply(condition_levels, function(level) {
        X[condition == level, , drop = FALSE]
    })
    y_list <- lapply(condition_levels, function(level) y[condition == level])
    names(X_list) <- names(y_list) <- condition_levels

    lambda_path <- if (is.null(lambda)) {
        .condition_make_lambda_path(
            X_list = X_list,
            y_list = y_list,
            alpha = alpha,
            condition_mix = condition_mix,
            condition_weight = condition_weight,
            nlambda = nlambda,
            lambda_min_ratio = lambda_min_ratio
        )
    } else {
        sort(unique(as.numeric(lambda)), decreasing = TRUE)
    }

    if (length(lambda_path) == 1L) {
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
            X_list = X_list,
            y_list = y_list,
            lambda = lambda_path,
            alpha = alpha,
            condition_mix = condition_mix,
            condition_weight = condition_weight,
            nfolds = nfolds,
            lambda_selection = lambda_selection,
            seed = seed,
            max_iter = max_iter,
            tol_objective = tol_objective,
            tol_coef = tol_coef
        )
    }

    full_path <- .condition_fit_multitask_path(
        X_list = X_list,
        y_list = y_list,
        lambda = lambda_path,
        alpha = alpha,
        condition_mix = condition_mix,
        condition_weight = condition_weight,
        max_iter = max_iter,
        tol_objective = tol_objective,
        tol_coef = tol_coef
    )
    selected_index <- match(cv$selected_lambda, full_path$lambda)
    selected <- full_path$fits[[selected_index]]
    beta <- selected$beta
    colnames(beta) <- condition_levels
    rownames(beta) <- edges$term

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

    universal_beta <- rowMeans(beta)
    universal_coefs <- .condition_format_coefs(
        edges, universal_beta, pooled_tf_corr[edges$tf]
    )
    condition_coefs <- lapply(condition_levels, function(level) {
        .condition_format_coefs(
            edges, beta[, level], condition_tf_corr[[level]][edges$tf]
        )
    })
    names(condition_coefs) <- condition_levels

    condition_gof <- vector('list', length(condition_levels))
    universal_rsq <- numeric(length(condition_levels))
    for (task in seq_along(condition_levels)) {
        level <- condition_levels[[task]]
        condition_prediction <- selected$intercept[[task]] + as.numeric(
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
        converged = selected$converged,
        iterations = selected$iterations,
        objective = selected$objective,
        coef_change = selected$coef_change,
        selected_lambda = selected$lambda,
        cv_mean = selected_cv_mean,
        cv_se = selected_cv_se,
        error_message = NA_character_,
        stringsAsFactors = FALSE
    )

    list(
        universal_coefs = universal_coefs,
        condition_coefs = condition_coefs,
        universal_gof = universal_gof,
        condition_gof = condition_gof,
        diagnostics = diagnostics
    )
}

.condition_build_design <- function(response_raw, gene_data, peak_data, edges, scale) {
    y <- as.numeric(response_raw[, 1L])
    tf_names <- unique(edges$tf)
    peak_names <- unique(edges$region)
    tf_matrix <- gene_data[, tf_names, drop = FALSE]
    peak_matrix <- peak_data[, peak_names, drop = FALSE]

    if (scale) {
        scaled_y <- .condition_scale_vector(y)
        if (is.null(scaled_y)) {
            stop('Target response has zero variance.')
        }
        y <- scaled_y
        tf_matrix <- .condition_scale_matrix(tf_matrix)
        peak_matrix <- .condition_scale_matrix(peak_matrix)
        edges <- edges[
            edges$tf %in% colnames(tf_matrix) & edges$region %in% colnames(peak_matrix),
            , drop = FALSE
        ]
    } else {
        response_variance <- stats::var(y)
        if (!is.finite(response_variance) || response_variance <= .Machine$double.eps) {
            stop('Target response has zero variance.')
        }
    }
    if (nrow(edges) == 0L) {
        stop('No edges remained after variance checks.')
    }

    tf_edge <- tf_matrix[, match(edges$tf, colnames(tf_matrix)), drop = FALSE]
    peak_edge <- peak_matrix[, match(edges$region, colnames(peak_matrix)), drop = FALSE]
    X <- tf_edge * peak_edge
    edge_variance <- .condition_column_variance(X)
    keep <- is.finite(edge_variance) & edge_variance > .Machine$double.eps
    X <- X[, keep, drop = FALSE]
    edges <- edges[keep, , drop = FALSE]
    edges$term <- paste(
        .condition_model_name(edges$region),
        .condition_model_name(edges$tf),
        sep = ':'
    )
    colnames(X) <- edges$term
    list(X = X, y = y, edges = edges)
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
    if (candidate_screen == 'motif_domain') {
        return(rep(TRUE, ncol(x)))
    }
    if (candidate_screen == 'pooled') {
        score <- abs(.condition_named_cor(x, y))
        return(is.finite(score) & score > threshold)
    }
    scores <- lapply(levels(condition), function(level) {
        abs(.condition_named_cor(
            x[condition == level, , drop = FALSE],
            y[condition == level, , drop = FALSE]
        ))
    })
    score <- Reduce(pmax, scores)
    is.finite(score) & score > threshold
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
