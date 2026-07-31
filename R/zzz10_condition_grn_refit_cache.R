# Cached sufficient statistics for support-constrained hierarchical refits.

.condition_refit_shared_baseline_reference <-
    .condition_refit_shared_baseline_unmodified
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
    if (is.null(cache)) return(NULL)
    keep <- as.logical(keep)
    task <- lapply(cache$task, function(value) {
        list(
            x_mean = value$x_mean[keep],
            y_mean = value$y_mean,
            gram = value$gram[keep, keep, drop = FALSE],
            rhs = value$rhs[keep]
        )
    })
    cache$task <- task
    cache$common_metric <- cache$common_metric[keep, keep, drop = FALSE]
    cache$predictor_names <- cache$predictor_names[keep]
    cache
}

.condition_refit_core_cached <- function(
    X_list,
    y_list,
    beta_selection,
    estimability_mask,
    ridge,
    active_tol = 1e-8,
    condition_weight = c('equal', 'cell_count'),
    max_iter = 200L,
    tol = 1e-8,
    cache = NULL
) {
    condition_weight <- match.arg(condition_weight)
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
    if (is.null(cache)) {
        cache <- .condition_make_refit_cache(
            X_list, y_list, condition_weight = condition_weight
        )
    }
    p <- nrow(beta_selection)
    n_tasks <- ncol(beta_selection)
    if (!identical(dim(cache$common_metric), c(p, p)) ||
        length(cache$task) != n_tasks) {
        stop('Refit cache dimensions do not match the model.')
    }
    support_mask <- estimability_mask & abs(beta_selection) > active_tol
    beta <- matrix(0, p, n_tasks, dimnames = dimnames(beta_selection))
    beta[support_mask] <- beta_selection[support_mask]
    shared <- .condition_update_shared(
        beta,
        estimability_mask,
        cache$average_weights,
        cache$common_metric
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
            gram <- cache$loss_weights[[task]] *
                task_cache$gram[active, active, drop = FALSE] +
                ridge * cache$average_weights[[task]] *
                cache$common_metric[active, active, drop = FALSE]
            rhs <- cache$loss_weights[[task]] * task_cache$rhs[active] +
                ridge * cache$average_weights[[task]] * as.numeric(
                    cache$common_metric[active, estimable, drop = FALSE] %*%
                        shared[estimable]
                )
            coefficient <- .condition_ridge_solve(gram, rhs)
            beta[active, task] <- coefficient
            intercept[[task]] <- task_cache$y_mean -
                sum(task_cache$x_mean[active] * coefficient)
        }
        shared <- .condition_update_shared(
            beta,
            estimability_mask,
            cache$average_weights,
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
        common_metric = 'pooled_weighted_predictor_gram_cached',
        iterations = iteration,
        coef_change = coef_change,
        converged = converged
    )
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
    cache = NULL
) {
    condition_weight <- match.arg(condition_weight)
    beta_selection <- as.matrix(beta_selection)
    actual <- .condition_true_variance_mask(X_list, estimability_mask)
    predictor_names <- rownames(beta_selection)
    conditions <- colnames(beta_selection)
    keep <- rowSums(actual) > 0L
    if (!any(keep)) {
        beta <- matrix(
            0, length(predictor_names), length(conditions),
            dimnames = list(predictor_names, conditions)
        )
        unavailable <- beta
        unavailable[,] <- NA_real_
        return(list(
            beta = beta,
            beta_condition = unavailable,
            beta_shared = stats::setNames(rep(0, length(predictor_names)), predictor_names),
            delta_condition = unavailable,
            support_mask = matrix(FALSE, length(predictor_names), length(conditions), dimnames = dimnames(beta)),
            active_mask = matrix(FALSE, length(predictor_names), length(conditions), dimnames = dimnames(beta)),
            estimability_mask = actual,
            intercept = stats::setNames(vapply(y_list, mean, numeric(1)), conditions),
            ridge = ridge,
            common_metric = 'pooled_weighted_predictor_gram_cached',
            iterations = 0L,
            coef_change = 0,
            converged = TRUE
        ))
    }
    subset_cache <- .condition_subset_refit_cache(cache, keep)
    fitted <- .condition_refit_core_cached(
        X_list = lapply(X_list, function(x) x[, keep, drop = FALSE]),
        y_list = y_list,
        beta_selection = beta_selection[keep, , drop = FALSE],
        estimability_mask = actual[keep, , drop = FALSE],
        ridge = ridge,
        active_tol = active_tol,
        condition_weight = condition_weight,
        max_iter = max_iter,
        tol = tol,
        cache = subset_cache
    )
    beta <- matrix(0, length(predictor_names), length(conditions),
                   dimnames = list(predictor_names, conditions))
    beta_condition <- delta_condition <- beta
    beta_condition[,] <- delta_condition[,] <- NA_real_
    support_mask <- active_mask <- matrix(
        FALSE, length(predictor_names), length(conditions), dimnames = dimnames(beta)
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
    fitted$estimability_mask <- actual
    fitted
}
