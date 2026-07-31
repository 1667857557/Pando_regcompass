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
    gram <- as.matrix(gram)
    diag(gram) <- diag(gram) + sqrt(.Machine$double.eps)
    answer <- tryCatch(
        solve(gram, rhs),
        error = function(error) qr.solve(gram, rhs)
    )
    answer <- if (is.matrix(rhs)) as.matrix(answer) else as.numeric(answer)
    if (any(!is.finite(answer))) {
        stop('Common-scale refit produced non-finite coefficients.')
    }
    answer
}

.condition_factorize_system <- function(lhs) {
    lhs <- as.matrix(lhs)
    if (!nrow(lhs)) {
        return(list(method = 'empty', factor = lhs, lhs = lhs))
    }
    diag(lhs) <- diag(lhs) + sqrt(.Machine$double.eps)
    factor <- tryCatch(chol(lhs), error = function(error) NULL)
    if (!is.null(factor)) {
        return(list(method = 'chol', factor = factor, lhs = lhs))
    }
    list(method = 'solve', factor = NULL, lhs = lhs)
}

.condition_factor_solve <- function(factorization, rhs) {
    if (!length(rhs)) {
        return(if (is.matrix(rhs)) rhs else numeric())
    }
    if (identical(factorization$method, 'chol')) {
        answer <- backsolve(
            factorization$factor,
            forwardsolve(t(factorization$factor), rhs)
        )
    } else {
        answer <- tryCatch(
            solve(factorization$lhs, rhs),
            error = function(error) qr.solve(factorization$lhs, rhs)
        )
    }
    answer <- if (is.matrix(rhs)) as.matrix(answer) else as.numeric(answer)
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
        n <- nrow(X_task)
        x_mean <- as.numeric(Matrix::colMeans(X_task))
        metric <- metric + loss_weights[[task]] * (
            as.matrix(Matrix::crossprod(X_task)) - n * tcrossprod(x_mean)
        )
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
        task_metric <- common_metric[estimable, estimable, drop = FALSE]
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

.condition_make_refit_cache <- function(
    X_list, y_list, condition_weight = c('equal', 'cell_count')
) {
    condition_weight <- match.arg(condition_weight)
    if (!length(X_list) || length(X_list) != length(y_list)) {
        stop('Refit cache requires aligned non-empty condition lists.')
    }
    p <- ncol(X_list[[1L]])
    loss_weights <- .condition_loss_weights(X_list, condition_weight)
    average_weights <- .condition_average_weights(X_list, condition_weight)
    task <- lapply(seq_along(X_list), function(index) {
        X <- X_list[[index]]
        y <- as.numeric(y_list[[index]])
        if (ncol(X) != p || nrow(X) != length(y)) {
            stop('Refit cache inputs are not aligned.')
        }
        n <- nrow(X)
        x_mean <- as.numeric(Matrix::colMeans(X))
        y_mean <- mean(y)
        gram <- as.matrix(Matrix::crossprod(X)) - n * tcrossprod(x_mean)
        rhs <- as.numeric(Matrix::crossprod(X, y)) - n * x_mean * y_mean
        gram[abs(gram) < .Machine$double.eps] <- 0
        list(
            x_mean = x_mean,
            y_mean = y_mean,
            gram = gram,
            rhs = rhs
        )
    })
    common_metric <- matrix(0, p, p)
    for (index in seq_along(task)) {
        common_metric <- common_metric +
            loss_weights[[index]] * task[[index]]$gram
    }
    diag(common_metric) <- diag(common_metric) + sqrt(.Machine$double.eps)
    structure(
        list(
            task = task,
            common_metric = common_metric,
            loss_weights = loss_weights,
            average_weights = average_weights,
            condition_weight = condition_weight,
            predictor_names = colnames(X_list[[1L]])
        ),
        class = c('ConditionRefitCache', 'list')
    )
}

.condition_subset_refit_cache <- function(cache, keep) {
    if (is.null(cache)) {
        return(NULL)
    }
    keep <- as.logical(keep)
    if (length(keep) == 0L || all(keep)) {
        return(cache)
    }
    cache$task <- lapply(cache$task, function(value) {
        list(
            x_mean = value$x_mean[keep],
            y_mean = value$y_mean,
            gram = value$gram[keep, keep, drop = FALSE],
            rhs = value$rhs[keep]
        )
    })
    cache$common_metric <- cache$common_metric[keep, keep, drop = FALSE]
    cache$predictor_names <- cache$predictor_names[keep]
    cache
}

.condition_refit_core_alternating <- function(
    beta_selection,
    estimability_mask,
    ridge,
    active_tol,
    cache,
    max_iter = 200L,
    tol = 1e-8
) {
    p <- nrow(beta_selection)
    n_tasks <- ncol(beta_selection)
    support_mask <- estimability_mask & abs(beta_selection) > active_tol
    beta <- matrix(0, p, n_tasks, dimnames = dimnames(beta_selection))
    beta[support_mask] <- beta_selection[support_mask]
    shared <- .condition_update_shared(
        beta, estimability_mask, cache$average_weights, cache$common_metric
    )
    intercept <- numeric(n_tasks)
    converged <- FALSE
    coef_change <- Inf
    iteration <- 0L
    for (iteration in seq_len(max_iter)) {
        beta_previous <- beta
        shared_previous <- shared
        for (task in seq_len(n_tasks)) {
            active <- support_mask[, task]
            beta[, task] <- 0
            task_cache <- cache$task[[task]]
            if (!any(active)) {
                intercept[[task]] <- task_cache$y_mean
                next
            }
            estimable <- estimability_mask[, task]
            lhs <- cache$loss_weights[[task]] *
                task_cache$gram[active, active, drop = FALSE] +
                ridge * cache$average_weights[[task]] *
                cache$common_metric[active, active, drop = FALSE]
            rhs <- cache$loss_weights[[task]] * task_cache$rhs[active] +
                ridge * cache$average_weights[[task]] * as.numeric(
                    cache$common_metric[active, estimable, drop = FALSE] %*%
                        shared[estimable]
                )
            coefficient <- .condition_ridge_solve(lhs, rhs)
            beta[active, task] <- coefficient
            intercept[[task]] <- task_cache$y_mean -
                sum(task_cache$x_mean[active] * coefficient)
        }
        shared <- .condition_update_shared(
            beta, estimability_mask, cache$average_weights,
            cache$common_metric
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
    list(
        beta = beta,
        shared = shared,
        intercept = intercept,
        support_mask = support_mask,
        iterations = iteration,
        coef_change = coef_change,
        converged = converged
    )
}

.condition_refit_core_direct <- function(
    beta_selection,
    estimability_mask,
    ridge,
    active_tol,
    cache
) {
    p <- nrow(beta_selection)
    n_tasks <- ncol(beta_selection)
    support_mask <- estimability_mask & abs(beta_selection) > active_tol
    shared_lhs <- matrix(0, p, p)
    task_plan <- vector('list', n_tasks)

    for (task in seq_len(n_tasks)) {
        estimable <- estimability_mask[, task]
        if (any(estimable)) {
            task_metric <- cache$common_metric[
                estimable, estimable, drop = FALSE
            ]
            shared_lhs[estimable, estimable] <-
                shared_lhs[estimable, estimable, drop = FALSE] +
                cache$average_weights[[task]] * task_metric
        }
        active <- support_mask[, task]
        if (!any(active)) {
            next
        }
        lhs <- cache$loss_weights[[task]] *
            cache$task[[task]]$gram[active, active, drop = FALSE] +
            ridge * cache$average_weights[[task]] *
            cache$common_metric[active, active, drop = FALSE]
        cross <- cache$common_metric[active, estimable, drop = FALSE]
        data_rhs <- cache$loss_weights[[task]] *
            cache$task[[task]]$rhs[active]
        factor <- .condition_factorize_system(lhs)
        task_plan[[task]] <- list(
            active = active,
            estimable = estimable,
            cross = cross,
            data_solution = .condition_factor_solve(factor, data_rhs),
            cross_solution = .condition_factor_solve(factor, cross),
            factor = factor,
            data_rhs = data_rhs
        )
    }

    schur <- shared_lhs
    shared_rhs <- numeric(p)
    for (task in seq_len(n_tasks)) {
        plan <- task_plan[[task]]
        if (is.null(plan)) {
            next
        }
        estimable <- plan$estimable
        weight <- cache$average_weights[[task]]
        shared_rhs[estimable] <- shared_rhs[estimable] +
            weight * as.numeric(crossprod(plan$cross, plan$data_solution))
        schur[estimable, estimable] <-
            schur[estimable, estimable, drop = FALSE] -
            ridge * weight^2 *
            crossprod(plan$cross, plan$cross_solution)
    }
    schur <- 0.5 * (schur + t(schur))
    shared <- .condition_factor_solve(
        .condition_factorize_system(schur), shared_rhs
    )

    beta <- matrix(0, p, n_tasks, dimnames = dimnames(beta_selection))
    intercept <- numeric(n_tasks)
    task_residual <- numeric(n_tasks)
    for (task in seq_len(n_tasks)) {
        plan <- task_plan[[task]]
        if (is.null(plan)) {
            intercept[[task]] <- cache$task[[task]]$y_mean
            next
        }
        coefficient <- plan$data_solution +
            ridge * cache$average_weights[[task]] * as.numeric(
                plan$cross_solution %*% shared[plan$estimable]
            )
        beta[plan$active, task] <- coefficient
        intercept[[task]] <- cache$task[[task]]$y_mean -
            sum(cache$task[[task]]$x_mean[plan$active] * coefficient)
        lhs_value <- plan$factor$lhs %*% coefficient
        rhs_value <- plan$data_rhs +
            ridge * cache$average_weights[[task]] * as.numeric(
                plan$cross %*% shared[plan$estimable]
            )
        task_residual[[task]] <- sqrt(sum((lhs_value - rhs_value)^2)) /
            (sqrt(sum(rhs_value^2)) + 1)
    }
    shared_residual <- shared_lhs %*% shared
    for (task in seq_len(n_tasks)) {
        estimable <- estimability_mask[, task]
        if (!any(estimable)) {
            next
        }
        shared_residual[estimable] <- shared_residual[estimable] -
            cache$average_weights[[task]] * as.numeric(
                cache$common_metric[estimable, , drop = FALSE] %*%
                    beta[, task]
            )
    }
    coef_change <- max(
        c(task_residual, sqrt(sum(shared_residual^2)) /
            (sqrt(sum(shared_rhs^2)) + 1)),
        na.rm = TRUE
    )
    list(
        beta = beta,
        shared = shared,
        intercept = intercept,
        support_mask = support_mask,
        iterations = 1L,
        coef_change = coef_change,
        converged = is.finite(coef_change) && coef_change < 1e-6
    )
}

.condition_finish_refit <- function(
    core, estimability_mask, ridge, predictor_names, conditions,
    common_metric_label
) {
    beta <- core$beta
    dimnames(beta) <- list(predictor_names, conditions)
    shared <- stats::setNames(as.numeric(core$shared), predictor_names)
    beta_condition <- beta
    beta_condition[!estimability_mask] <- NA_real_
    delta_condition <- sweep(beta_condition, 1L, shared, '-')
    active_mask <- estimability_mask & abs(beta) > 1e-8
    list(
        beta = beta,
        beta_condition = beta_condition,
        beta_shared = shared,
        delta_condition = delta_condition,
        support_mask = core$support_mask,
        active_mask = active_mask,
        estimability_mask = estimability_mask,
        intercept = stats::setNames(core$intercept, conditions),
        ridge = ridge,
        common_metric = common_metric_label,
        iterations = core$iterations,
        coef_change = core$coef_change,
        converged = core$converged
    )
}

.condition_refit_shared_baseline_reference <- function(
    X_list,
    y_list,
    beta_selection,
    estimability_mask,
    ridge,
    active_tol = 1e-8,
    condition_weight = c('equal', 'cell_count'),
    max_iter = 200L,
    tol = 1e-8,
    cache = NULL,
    estimability_verified = FALSE
) {
    condition_weight <- match.arg(condition_weight)
    beta_selection <- as.matrix(beta_selection)
    actual <- if (isTRUE(estimability_verified)) {
        as.matrix(estimability_mask)
    } else {
        .condition_true_variance_mask(X_list, estimability_mask)
    }
    keep <- rowSums(actual) > 0L
    predictor_names <- rownames(beta_selection)
    conditions <- colnames(beta_selection)
    if (!any(keep)) {
        return(.condition_empty_refit(
            y_list, predictor_names, conditions, actual, ridge
        ))
    }
    cache <- .condition_subset_refit_cache(cache, keep)
    if (is.null(cache)) {
        cache <- .condition_make_refit_cache(
            lapply(X_list, function(x) x[, keep, drop = FALSE]),
            y_list,
            condition_weight
        )
    }
    core <- .condition_refit_core_alternating(
        beta_selection[keep, , drop = FALSE],
        actual[keep, , drop = FALSE],
        ridge,
        active_tol,
        cache,
        max_iter,
        tol
    )
    fitted <- .condition_finish_refit(
        core,
        actual[keep, , drop = FALSE],
        ridge,
        predictor_names[keep],
        conditions,
        'pooled_weighted_predictor_gram_alternating_reference'
    )
    .condition_expand_refit(fitted, keep, predictor_names, conditions, actual)
}

.condition_empty_refit <- function(
    y_list, predictor_names, conditions, estimability_mask, ridge
) {
    beta <- matrix(
        0, length(predictor_names), length(conditions),
        dimnames = list(predictor_names, conditions)
    )
    unavailable <- beta
    unavailable[,] <- NA_real_
    list(
        beta = beta,
        beta_condition = unavailable,
        beta_shared = stats::setNames(rep(0, length(predictor_names)), predictor_names),
        delta_condition = unavailable,
        support_mask = matrix(FALSE, nrow(beta), ncol(beta), dimnames = dimnames(beta)),
        active_mask = matrix(FALSE, nrow(beta), ncol(beta), dimnames = dimnames(beta)),
        estimability_mask = estimability_mask,
        intercept = stats::setNames(vapply(y_list, mean, numeric(1)), conditions),
        ridge = ridge,
        common_metric = 'pooled_weighted_predictor_gram_direct_schur',
        iterations = 0L,
        coef_change = 0,
        converged = TRUE
    )
}

.condition_expand_refit <- function(
    fitted, keep, predictor_names, conditions, estimability_mask
) {
    beta <- matrix(
        0, length(predictor_names), length(conditions),
        dimnames = list(predictor_names, conditions)
    )
    beta_condition <- delta_condition <- beta
    beta_condition[,] <- delta_condition[,] <- NA_real_
    support_mask <- active_mask <- matrix(
        FALSE, length(predictor_names), length(conditions),
        dimnames = dimnames(beta)
    )
    beta_shared <- stats::setNames(rep(0, length(predictor_names)), predictor_names)
    beta[keep, ] <- fitted$beta
    beta_condition[keep, ] <- fitted$beta_condition
    delta_condition[keep, ] <- fitted$delta_condition
    support_mask[keep, ] <- fitted$support_mask
    active_mask[keep, ] <- fitted$active_mask
    beta_shared[keep] <- fitted$beta_shared
    fitted$beta <- beta
    fitted$beta_condition <- beta_condition
    fitted$beta_shared <- beta_shared
    fitted$delta_condition <- delta_condition
    fitted$support_mask <- support_mask
    fitted$active_mask <- active_mask
    fitted$estimability_mask <- estimability_mask
    fitted
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
    tol = 1e-8,
    cache = NULL,
    estimability_verified = FALSE,
    solver = getOption('Pando.condition_refit_solver', 'direct')
) {
    condition_weight <- match.arg(condition_weight)
    solver <- match.arg(solver, c('direct', 'alternating'))
    beta_selection <- as.matrix(beta_selection)
    estimability_mask <- as.matrix(estimability_mask)
    if (!identical(dim(beta_selection), dim(estimability_mask))) {
        stop('beta_selection and estimability_mask must have identical dimensions.')
    }
    if (!is.logical(estimability_mask) || anyNA(estimability_mask)) {
        stop('estimability_mask must be logical without NA values.')
    }
    if (!is.numeric(ridge) || length(ridge) != 1L ||
        !is.finite(ridge) || ridge <= 0) {
        stop('ridge must be one finite positive value.')
    }
    actual <- if (isTRUE(estimability_verified)) {
        estimability_mask
    } else {
        .condition_true_variance_mask(X_list, estimability_mask)
    }
    predictor_names <- rownames(beta_selection)
    conditions <- colnames(beta_selection)
    keep <- rowSums(actual) > 0L
    if (!any(keep)) {
        return(.condition_empty_refit(
            y_list, predictor_names, conditions, actual, ridge
        ))
    }
    subset_X <- if (all(keep)) X_list else {
        lapply(X_list, function(x) x[, keep, drop = FALSE])
    }
    subset_cache <- .condition_subset_refit_cache(cache, keep)
    if (is.null(subset_cache)) {
        subset_cache <- .condition_make_refit_cache(
            subset_X, y_list, condition_weight
        )
    }
    core <- if (solver == 'direct') {
        .condition_refit_core_direct(
            beta_selection[keep, , drop = FALSE],
            actual[keep, , drop = FALSE],
            ridge,
            active_tol,
            subset_cache
        )
    } else {
        .condition_refit_core_alternating(
            beta_selection[keep, , drop = FALSE],
            actual[keep, , drop = FALSE],
            ridge,
            active_tol,
            subset_cache,
            max_iter,
            tol
        )
    }
    fitted <- .condition_finish_refit(
        core,
        actual[keep, , drop = FALSE],
        ridge,
        predictor_names[keep],
        conditions,
        if (solver == 'direct') {
            'pooled_weighted_predictor_gram_direct_schur'
        } else {
            'pooled_weighted_predictor_gram_alternating'
        }
    )
    .condition_expand_refit(fitted, keep, predictor_names, conditions, actual)
}
