# Support-constrained refitting for condition-specific sub-GRNs.

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

.condition_ridge_solve <- function(gram, rhs) {
    if (!length(rhs)) {
        return(numeric())
    }
    diag(gram) <- diag(gram) + sqrt(.Machine$double.eps)
    answer <- tryCatch(
        as.numeric(solve(gram, rhs)),
        error = function(error) {
            as.numeric(qr.solve(as.matrix(gram), as.numeric(rhs)))
        }
    )
    if (any(!is.finite(answer))) {
        stop('Common-scale refit produced non-finite coefficients.')
    }
    answer
}

.condition_common_metric <- function(X_list, loss_weights) {
    p <- ncol(X_list[[1L]])
    metric <- matrix(0, p, p)
    for (task in seq_along(X_list)) {
        X_task <- X_list[[task]]
        centered <- sweep(
            as.matrix(X_task), 2L, as.numeric(Matrix::colMeans(X_task)), '-'
        )
        metric <- metric + loss_weights[[task]] * crossprod(centered)
    }
    diag(metric) <- diag(metric) + sqrt(.Machine$double.eps)
    metric
}

.condition_update_shared <- function(
    beta, estimability_mask, weights, common_metric
) {
    p <- nrow(beta)
    lhs <- matrix(0, p, p)
    rhs <- numeric(p)
    for (task in seq_len(ncol(beta))) {
        estimable <- estimability_mask[, task]
        if (!any(estimable)) {
            next
        }
        task_metric <- common_metric[
            estimable, estimable, drop = FALSE
        ]
        lhs[estimable, estimable] <-
            lhs[estimable, estimable, drop = FALSE] +
            weights[[task]] * task_metric
        rhs[estimable] <- rhs[estimable] +
            weights[[task]] * as.numeric(
                task_metric %*% beta[estimable, task]
            )
    }
    .condition_ridge_solve(lhs, rhs)
}

.condition_refit_shared_baseline <- function(
    X_list,
    y_list,
    beta_selection,
    estimability_mask,
    ridge,
    active_tol = 1e-8,
    condition_weight = c('equal', 'cell_count'),
    max_iter = 200L,
    tol = 1e-8
) {
    condition_weight <- match.arg(condition_weight)
    beta_selection <- as.matrix(beta_selection)
    estimability_mask <- as.matrix(estimability_mask)
    if (!identical(dim(beta_selection), dim(estimability_mask))) {
        stop('beta_selection and estimability_mask must have identical dimensions.')
    }
    if (!is.logical(estimability_mask) || anyNA(estimability_mask)) {
        stop('estimability_mask must be a logical matrix without NA values.')
    }
    if (!is.numeric(ridge) || length(ridge) != 1L ||
        !is.finite(ridge) || ridge <= 0) {
        stop('ridge must be one finite positive value.')
    }
    if (!is.numeric(active_tol) || length(active_tol) != 1L ||
        !is.finite(active_tol) || active_tol < 0) {
        stop('active_tol must be finite and non-negative.')
    }

    p <- nrow(beta_selection)
    n_tasks <- ncol(beta_selection)
    support_mask <- estimability_mask & abs(beta_selection) > active_tol
    beta <- matrix(
        0, p, n_tasks, dimnames = dimnames(beta_selection)
    )
    beta[support_mask] <- beta_selection[support_mask]
    average_weights <- .condition_average_weights(X_list, condition_weight)
    loss_weights <- .condition_loss_weights(X_list, condition_weight)
    common_metric <- .condition_common_metric(X_list, loss_weights)
    shared <- .condition_update_shared(
        beta, estimability_mask, average_weights, common_metric
    )
    intercept <- numeric(n_tasks)
    converged <- FALSE
    coef_change <- Inf

    for (iteration in seq_len(max_iter)) {
        beta_previous <- beta
        shared_previous <- shared
        for (task in seq_len(n_tasks)) {
            active <- support_mask[, task]
            beta[, task] <- 0
            y_task <- y_list[[task]]
            if (!any(active)) {
                intercept[[task]] <- mean(y_task)
                next
            }
            estimable <- estimability_mask[, task]
            X_task <- X_list[[task]][, active, drop = FALSE]
            x_mean <- as.numeric(Matrix::colMeans(X_task))
            y_mean <- mean(y_task)
            X_centered <- sweep(as.matrix(X_task), 2L, x_mean, '-')
            y_centered <- y_task - y_mean
            gram <- loss_weights[[task]] * crossprod(X_centered) +
                ridge * average_weights[[task]] *
                common_metric[active, active, drop = FALSE]
            rhs <- loss_weights[[task]] *
                as.numeric(crossprod(X_centered, y_centered)) +
                ridge * average_weights[[task]] * as.numeric(
                    common_metric[active, estimable, drop = FALSE] %*%
                        shared[estimable]
                )
            coefficient <- .condition_ridge_solve(gram, rhs)
            beta[active, task] <- coefficient
            intercept[[task]] <- y_mean - sum(x_mean * coefficient)
        }
        shared <- .condition_update_shared(
            beta, estimability_mask, average_weights, common_metric
        )
        numerator <- sqrt(sum((beta - beta_previous)^2)) +
            sqrt(sum((shared - shared_previous)^2))
        denominator <- sqrt(sum(beta_previous^2)) +
            sqrt(sum(shared_previous^2)) + .Machine$double.eps
        coef_change <- numerator / denominator
        if (coef_change < tol) {
            converged <- TRUE
            break
        }
    }

    beta_condition <- beta
    beta_condition[!estimability_mask] <- NA_real_
    delta_condition <- sweep(beta_condition, 1L, shared, '-')
    active_mask <- estimability_mask & abs(beta) > active_tol
    list(
        beta = beta,
        beta_condition = beta_condition,
        beta_shared = stats::setNames(shared, rownames(beta)),
        delta_condition = delta_condition,
        support_mask = support_mask,
        active_mask = active_mask,
        estimability_mask = estimability_mask,
        intercept = stats::setNames(intercept, colnames(beta)),
        ridge = ridge,
        common_metric = 'pooled_weighted_predictor_gram',
        iterations = iteration,
        coef_change = coef_change,
        converged = converged
    )
}

.condition_refit_stability_diagnostics <- function(
    X_list,
    y_list,
    beta_selection,
    estimability_mask,
    ridge,
    comparison_conditions,
    selection_lambda,
    alpha,
    condition_mix,
    active_tol = 1e-8,
    condition_weight = 'equal',
    bootstrap_replicates = 20L,
    seed = 12345L,
    max_iter = 5000L,
    tol_objective = 1e-7,
    tol_coef = 1e-6
) {
    conditions <- names(X_list)
    comparison_conditions <- unique(as.character(comparison_conditions))
    if (length(comparison_conditions) != 2L ||
        !all(comparison_conditions %in% conditions)) {
        return(list(
            edge = data.frame(),
            status = 'requires_exactly_two_comparison_conditions'
        ))
    }
    pair_index <- match(comparison_conditions, conditions)
    ridge_grid <- ridge * c(0.1, 1, 10)
    ridge_fits <- lapply(ridge_grid, function(value) {
        .condition_refit_shared_baseline(
            X_list = X_list,
            y_list = y_list,
            beta_selection = beta_selection,
            estimability_mask = estimability_mask,
            ridge = value,
            active_tol = active_tol,
            condition_weight = condition_weight
        )
    })
    unshrunk <- .condition_refit_shared_baseline(
        X_list = X_list,
        y_list = y_list,
        beta_selection = beta_selection,
        estimability_mask = estimability_mask,
        ridge = max(ridge * 1e-6, sqrt(.Machine$double.eps)),
        active_tol = active_tol,
        condition_weight = condition_weight
    )
    contrast_matrix <- vapply(ridge_fits, function(fit) {
        fit$beta_condition[, pair_index[[2L]]] -
            fit$beta_condition[, pair_index[[1L]]]
    }, numeric(nrow(beta_selection)))
    if (is.null(dim(contrast_matrix))) {
        contrast_matrix <- matrix(
            contrast_matrix, nrow = nrow(beta_selection)
        )
    }
    main_contrast <- contrast_matrix[, 2L]
    unshrunk_contrast <-
        unshrunk$beta_condition[, pair_index[[2L]]] -
        unshrunk$beta_condition[, pair_index[[1L]]]
    contrast_retention <- rep(NA_real_, length(main_contrast))
    retained <- is.finite(unshrunk_contrast) &
        abs(unshrunk_contrast) > active_tol
    contrast_retention[retained] <-
        main_contrast[retained] / unshrunk_contrast[retained]
    sign_stable <- apply(contrast_matrix, 1L, function(value) {
        value <- sign(value[is.finite(value) & abs(value) > active_tol])
        length(value) > 0L && length(unique(value)) == 1L
    })
    sensitivity_range <- apply(
        contrast_matrix, 1L, function(value) {
            value <- value[is.finite(value)]
            if (length(value)) diff(range(value)) else NA_real_
        }
    )
    variance_by_condition <- vapply(
        X_list, .condition_population_variance,
        numeric(nrow(beta_selection))
    )
    variance_ratio <- apply(
        variance_by_condition, 1L, function(value) {
            value <- value[is.finite(value) & value > .Machine$double.eps]
            if (length(value) >= 2L) max(value) / min(value) else NA_real_
        }
    )
    bootstrap_replicates <- as.integer(bootstrap_replicates)
    selected_count <- matrix(
        0, nrow(beta_selection), ncol(beta_selection),
        dimnames = dimnames(beta_selection)
    )
    successful <- 0L
    set.seed(seed)
    for (iteration in seq_len(bootstrap_replicates)) {
        index <- lapply(X_list, function(x) {
            sample.int(nrow(x), nrow(x), replace = TRUE)
        })
        boot_X <- lapply(seq_along(X_list), function(task) {
            X_list[[task]][index[[task]], , drop = FALSE]
        })
        boot_y <- lapply(seq_along(y_list), function(task) {
            y_list[[task]][index[[task]]]
        })
        names(boot_X) <- names(X_list)
        names(boot_y) <- names(y_list)
        boot <- tryCatch({
            boot_selection <- .condition_fit_multitask_path(
                X_list = boot_X,
                y_list = boot_y,
                lambda = selection_lambda,
                alpha = alpha,
                condition_mix = condition_mix,
                condition_weight = condition_weight,
                coefficient_mask = estimability_mask,
                max_iter = max_iter,
                tol_objective = tol_objective,
                tol_coef = tol_coef
            )$fits[[1L]]
            .condition_refit_shared_baseline(
                X_list = boot_X,
                y_list = boot_y,
                beta_selection = boot_selection$beta,
                estimability_mask = estimability_mask,
                ridge = ridge,
                active_tol = active_tol,
                condition_weight = condition_weight
            )
        },
            error = function(error) NULL
        )
        if (is.null(boot)) next
        selected_count <- selected_count +
            (is.finite(boot$beta_condition) &
             abs(boot$beta_condition) > active_tol)
        successful <- successful + 1L
    }
    bootstrap_frequency <- if (successful) {
        selected_count / successful
    } else {
        selected_count * NA_real_
    }
    pair_frequency <- apply(
        bootstrap_frequency[, pair_index, drop = FALSE],
        1L, min, na.rm = TRUE
    )
    pair_frequency[!is.finite(pair_frequency)] <- NA_real_
    list(
        edge = data.frame(
            edge_id = rownames(beta_selection),
            condition_a = comparison_conditions[[1L]],
            condition_b = comparison_conditions[[2L]],
            contrast_retention = contrast_retention,
            contrast_sign_stable = sign_stable,
            predictor_variance_ratio = variance_ratio,
            ridge_sensitivity_range = sensitivity_range,
            bootstrap_selection_frequency = pair_frequency,
            stringsAsFactors = FALSE
        ),
        ridge_grid = ridge_grid,
        bootstrap_replicates_requested = bootstrap_replicates,
        bootstrap_replicates_successful = successful,
        bootstrap_definition =
            paste(
                'within-condition cell bootstrap with sparse support',
                'reselection at the selected full-data lambda'
            ),
        unshrunk_definition =
            'same selected support with near-zero shared-baseline ridge',
        status = if (successful == bootstrap_replicates) {
            'complete'
        } else {
            'partial_bootstrap'
        }
    )
}
