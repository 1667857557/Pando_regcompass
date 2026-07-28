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
