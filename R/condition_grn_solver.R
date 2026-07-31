# Internal numerical routines for condition-aware Pando.

.condition_sqnorm <- function(x) {
    if (inherits(x, 'sparseMatrix')) {
        return(sum(x@x * x@x))
    }
    sum(x * x)
}

.condition_centered_sqnorm <- function(x) {
    n <- nrow(x)
    if (!n || !ncol(x)) {
        return(0)
    }
    raw <- .condition_sqnorm(x)
    mu <- as.numeric(Matrix::colMeans(x))
    max(raw - n * sum(mu * mu), 0)
}

.condition_column_variance <- function(x) {
    n <- nrow(x)
    if (n <= 1L) {
        return(rep(0, ncol(x)))
    }
    mu <- as.numeric(Matrix::colMeans(x))
    ss <- as.numeric(Matrix::colSums(x * x)) - n * mu * mu
    pmax(ss / (n - 1), 0)
}

.condition_rsq <- function(y, prediction) {
    denominator <- sum((y - mean(y))^2)
    if (!is.finite(denominator) || denominator <= .Machine$double.eps) {
        return(NA_real_)
    }
    1 - sum((y - prediction)^2) / denominator
}

.condition_pooled_task_rsq <- function(y_list, prediction_list) {
    if (!is.list(y_list) || !is.list(prediction_list) ||
        length(y_list) != length(prediction_list) || !length(y_list)) {
        stop('y_list and prediction_list must contain matching tasks.')
    }
    residual_ss <- 0
    within_task_ss <- 0
    for (task in seq_along(y_list)) {
        y <- as.numeric(y_list[[task]])
        prediction <- as.numeric(prediction_list[[task]])
        if (length(y) != length(prediction) ||
            any(!is.finite(y)) || any(!is.finite(prediction))) {
            stop('Pooled task R-squared requires aligned finite values.')
        }
        residual_ss <- residual_ss + sum((y - prediction)^2)
        within_task_ss <- within_task_ss + sum((y - mean(y))^2)
    }
    if (!is.finite(within_task_ss) ||
        within_task_ss <= .Machine$double.eps) {
        return(NA_real_)
    }
    1 - residual_ss / within_task_ss
}

.condition_task_validation_loss <- function(
    task_mse, n_validation, condition_weight = c('equal', 'cell_count')
) {
    condition_weight <- match.arg(condition_weight)
    task_mse <- as.numeric(task_mse)
    n_validation <- as.numeric(n_validation)
    if (!length(task_mse) || length(task_mse) != length(n_validation) ||
        any(!is.finite(task_mse)) || any(task_mse < 0) ||
        any(!is.finite(n_validation)) || any(n_validation <= 0)) {
        stop(
            'Validation losses require finite task MSE and positive task sizes.'
        )
    }
    if (condition_weight == 'equal') {
        return(mean(task_mse))
    }
    stats::weighted.mean(task_mse, n_validation)
}

.condition_loss_weights <- function(
    X_list, condition_weight = c('equal', 'cell_count')
) {
    condition_weight <- match.arg(condition_weight)
    n_condition <- vapply(X_list, nrow, integer(1))
    if (any(n_condition <= 0L)) {
        stop('Every condition must contain at least one observation.')
    }
    if (condition_weight == 'equal') {
        return(1 / n_condition)
    }
    rep(1 / sum(n_condition), length(n_condition))
}

.condition_profiled_smooth <- function(B, X_list, y_list, loss_weights, ridge) {
    p <- nrow(B)
    n_tasks <- ncol(B)
    gradient <- matrix(0, nrow = p, ncol = n_tasks)
    intercept <- numeric(n_tasks)
    value <- 0
    for (task in seq_len(n_tasks)) {
        xb <- as.numeric(X_list[[task]] %*% B[, task])
        intercept[[task]] <- mean(y_list[[task]] - xb)
        residual <- y_list[[task]] - intercept[[task]] - xb
        value <- value + 0.5 * loss_weights[[task]] * sum(residual * residual)
        gradient[, task] <- -loss_weights[[task]] * as.numeric(
            crossprod(X_list[[task]], residual)
        )
    }
    if (ridge > 0) {
        value <- value + 0.5 * ridge * sum(B * B)
        gradient <- gradient + ridge * B
    }
    list(value = value, gradient = gradient, intercept = intercept)
}

.condition_sparse_group_penalty <- function(
    B,
    lambda,
    alpha,
    condition_mix,
    group_factor = NULL,
    element_factor = NULL,
    coefficient_mask = NULL
) {
    if (is.null(group_factor)) {
        group_factor <- rep(1, nrow(B))
    }
    if (is.null(element_factor)) {
        element_factor <- matrix(1, nrow(B), ncol(B))
    }
    sparse_strength <- lambda * alpha
    group_penalty <- sparse_strength * (1 - condition_mix) *
        sum(group_factor * sqrt(rowSums(B * B)))
    element_penalty <- sparse_strength * condition_mix *
        sum(element_factor * abs(B))
    group_penalty + element_penalty
}

.condition_sparse_group_prox <- function(
    V,
    step,
    lambda,
    alpha,
    condition_mix,
    group_factor = NULL,
    element_factor = NULL,
    coefficient_mask = NULL
) {
    if (is.null(group_factor)) {
        group_factor <- rep(1, nrow(V))
    }
    if (is.null(element_factor)) {
        element_factor <- matrix(1, nrow(V), ncol(V))
    }
    if (!is.null(coefficient_mask)) {
        V[!coefficient_mask] <- 0
    }
    element_threshold <- step * lambda * alpha * condition_mix * element_factor
    U <- sign(V) * pmax(abs(V) - element_threshold, 0)
    group_threshold <- step * lambda * alpha * (1 - condition_mix) * group_factor
    row_norm <- sqrt(rowSums(U * U))
    row_scale <- rep(0, length(row_norm))
    nonzero <- row_norm > 0
    row_scale[nonzero] <- pmax(
        1 - group_threshold[nonzero] / row_norm[nonzero], 0
    )
    U <- U * row_scale
    if (!is.null(coefficient_mask)) {
        U[!coefficient_mask] <- 0
    }
    U
}

.condition_initial_step <- function(X_list, loss_weights, ridge) {
    upper_bound <- ridge
    for (task in seq_along(X_list)) {
        upper_bound <- upper_bound +
            loss_weights[[task]] * .condition_centered_sqnorm(X_list[[task]])
    }
    if (!is.finite(upper_bound) || upper_bound <= 0) {
        return(1)
    }
    1 / upper_bound
}

.condition_fit_multitask_lambda <- function(
    X_list,
    y_list,
    lambda,
    alpha = 0.5,
    condition_mix = 0.5,
    condition_weight = c('equal', 'cell_count'),
    initial_B = NULL,
    initial_step = NULL,
    group_factor = NULL,
    element_factor = NULL,
    coefficient_mask = NULL,
    max_iter = 5000L,
    tol_objective = 1e-7,
    tol_coef = 1e-6,
    backtrack = 0.5,
    min_step = 1e-14,
    keep_history = FALSE
) {
    condition_weight <- match.arg(condition_weight)
    if (length(X_list) < 2L || length(X_list) != length(y_list)) {
        stop('X_list and y_list must contain the same two or more conditions.')
    }
    if (!is.numeric(lambda) || length(lambda) != 1L ||
        !is.finite(lambda) || lambda < 0) {
        stop('lambda must be one finite non-negative number.')
    }
    if (!is.numeric(alpha) || length(alpha) != 1L ||
        alpha < 0 || alpha > 1) {
        stop('alpha must be between 0 and 1.')
    }
    if (!is.numeric(condition_mix) || length(condition_mix) != 1L ||
        condition_mix < 0 || condition_mix > 1) {
        stop('condition_mix must be between 0 and 1.')
    }
    p <- ncol(X_list[[1L]])
    n_tasks <- length(X_list)
    if (p == 0L) {
        stop('At least one predictor is required.')
    }
    for (task in seq_len(n_tasks)) {
        if (ncol(X_list[[task]]) != p ||
            nrow(X_list[[task]]) != length(y_list[[task]])) {
            stop(
                'All task matrices must share predictor columns and ',
                'match their response lengths.'
            )
        }
    }
    if (is.null(group_factor)) {
        group_factor <- rep(1, p)
    }
    if (is.null(element_factor)) {
        element_factor <- matrix(1, p, n_tasks)
    }
    if (is.null(coefficient_mask)) {
        coefficient_mask <- matrix(
            TRUE, p, n_tasks,
            dimnames = list(colnames(X_list[[1L]]), names(X_list))
        )
    }
    coefficient_mask <- as.matrix(coefficient_mask)
    if (!identical(dim(coefficient_mask), c(p, n_tasks)) ||
        !is.logical(coefficient_mask) || anyNA(coefficient_mask) ||
        any(rowSums(coefficient_mask) == 0L)) {
        stop(
            'coefficient_mask must be a logical predictors-by-conditions ',
            'matrix with one eligible condition per predictor.'
        )
    }
    loss_weights <- .condition_loss_weights(X_list, condition_weight)
    ridge <- lambda * (1 - alpha)
    B <- if (is.null(initial_B)) matrix(0, p, n_tasks) else as.matrix(initial_B)
    if (!identical(dim(B), c(p, n_tasks))) {
        stop('initial_B must have dimensions predictors by conditions.')
    }
    B[!coefficient_mask] <- 0
    Z <- B
    acceleration <- 1
    safe_step <- .condition_initial_step(X_list, loss_weights, ridge)
    step <- if (is.null(initial_step)) {
        safe_step
    } else {
        max(as.numeric(initial_step), safe_step)
    }
    smooth_B <- .condition_profiled_smooth(
        B, X_list, y_list, loss_weights, ridge
    )
    objective_previous <- smooth_B$value +
        .condition_sparse_group_penalty(
            B, lambda, alpha, condition_mix, group_factor, element_factor
        )
    history <- if (keep_history) objective_previous else NULL
    converged <- FALSE
    coef_change <- Inf
    objective_change <- Inf
    iteration <- 0L
    for (iteration in seq_len(max_iter)) {
        smooth_Z <- if (identical(Z, B)) {
            smooth_B
        } else {
            .condition_profiled_smooth(Z, X_list, y_list, loss_weights, ridge)
        }
        repeat {
            candidate <- .condition_sparse_group_prox(
                Z - step * smooth_Z$gradient,
                step,
                lambda,
                alpha,
                condition_mix,
                group_factor,
                element_factor,
                coefficient_mask
            )
            smooth_candidate <- .condition_profiled_smooth(
                candidate, X_list, y_list, loss_weights, ridge
            )
            difference <- candidate - Z
            quadratic_bound <- smooth_Z$value +
                sum(smooth_Z$gradient * difference) +
                sum(difference * difference) / (2 * step)
            if (smooth_candidate$value <= quadratic_bound + 1e-10) {
                break
            }
            step <- step * backtrack
            if (step < min_step) {
                stop('Backtracking line search reached min_step.')
            }
        }
        objective_candidate <- smooth_candidate$value +
            .condition_sparse_group_penalty(
                candidate, lambda, alpha, condition_mix,
                group_factor, element_factor
            )
        if (objective_candidate > objective_previous + 1e-10) {
            acceleration <- 1
            Z <- B
            smooth_Z <- smooth_B
            repeat {
                candidate <- .condition_sparse_group_prox(
                    Z - step * smooth_Z$gradient,
                    step,
                    lambda,
                    alpha,
                    condition_mix,
                    group_factor,
                    element_factor,
                    coefficient_mask
                )
                smooth_candidate <- .condition_profiled_smooth(
                    candidate, X_list, y_list, loss_weights, ridge
                )
                difference <- candidate - Z
                quadratic_bound <- smooth_Z$value +
                    sum(smooth_Z$gradient * difference) +
                    sum(difference * difference) / (2 * step)
                if (smooth_candidate$value <= quadratic_bound + 1e-10) {
                    break
                }
                step <- step * backtrack
                if (step < min_step) {
                    stop(
                        'Backtracking line search reached min_step after restart.'
                    )
                }
            }
            objective_candidate <- smooth_candidate$value +
                .condition_sparse_group_penalty(
                    candidate, lambda, alpha, condition_mix,
                    group_factor, element_factor
                )
        }
        coef_change <- sqrt(sum((candidate - B)^2)) /
            (sqrt(sum(B^2)) + .Machine$double.eps)
        objective_change <- abs(objective_candidate - objective_previous) /
            (abs(objective_previous) + .Machine$double.eps)
        if (keep_history) {
            history <- c(history, objective_candidate)
        }
        B_previous <- B
        B <- candidate
        smooth_B <- smooth_candidate
        objective_previous <- objective_candidate
        if (objective_change < tol_objective && coef_change < tol_coef) {
            converged <- TRUE
            break
        }
        acceleration_new <- (1 + sqrt(1 + 4 * acceleration^2)) / 2
        Z_new <- B + ((acceleration - 1) / acceleration_new) *
            (B - B_previous)
        restart <- sum((Z - B) * (B - B_previous)) > 0
        if (restart) {
            acceleration <- 1
            Z <- B
        } else {
            acceleration <- acceleration_new
            Z <- Z_new
        }
    }
    list(
        beta = B,
        intercept = smooth_B$intercept,
        lambda = lambda,
        objective = objective_previous,
        objective_change = objective_change,
        coef_change = coef_change,
        iterations = iteration,
        converged = converged,
        step = step,
        history = history
    )
}

.condition_lambda_max <- function(
    X_list,
    y_list,
    alpha,
    condition_mix,
    condition_weight = c('equal', 'cell_count'),
    coefficient_mask = NULL
) {
    condition_weight <- match.arg(condition_weight)
    p <- ncol(X_list[[1L]])
    n_tasks <- length(X_list)
    zero <- matrix(0, p, n_tasks)
    loss_weights <- .condition_loss_weights(X_list, condition_weight)
    gradient <- .condition_profiled_smooth(
        zero, X_list, y_list, loss_weights, ridge = 0
    )$gradient
    if (!is.null(coefficient_mask)) {
        gradient[!coefficient_mask] <- 0
    }
    if (alpha <= .Machine$double.eps) {
        value <- max(abs(gradient))
    } else if (condition_mix <= .Machine$double.eps) {
        value <- max(sqrt(rowSums(gradient * gradient))) / alpha
    } else {
        value <- max(abs(gradient)) / (alpha * condition_mix)
    }
    if (!is.finite(value) || value <= 0) {
        value <- 1
    }
    value * 1.001
}

.condition_make_lambda_path <- function(
    X_list,
    y_list,
    alpha = 0.5,
    condition_mix = 0.5,
    condition_weight = c('equal', 'cell_count'),
    coefficient_mask = NULL,
    nlambda = 50L,
    lambda_min_ratio = NULL
) {
    condition_weight <- match.arg(condition_weight)
    if (nlambda < 2L) {
        stop('nlambda must be at least 2 when lambda is not supplied.')
    }
    if (is.null(lambda_min_ratio)) {
        n_cells <- sum(vapply(X_list, nrow, integer(1)))
        lambda_min_ratio <- if (n_cells < ncol(X_list[[1L]])) 0.01 else 1e-4
    }
    if (!is.finite(lambda_min_ratio) || lambda_min_ratio <= 0 ||
        lambda_min_ratio >= 1) {
        stop('lambda_min_ratio must be between 0 and 1.')
    }
    lambda_max <- .condition_lambda_max(
        X_list, y_list, alpha, condition_mix, condition_weight,
        coefficient_mask = coefficient_mask
    )
    exp(seq(
        log(lambda_max),
        log(lambda_max * lambda_min_ratio),
        length.out = nlambda
    ))
}

.condition_fit_multitask_path_reference <- function(
    X_list,
    y_list,
    lambda,
    alpha = 0.5,
    condition_mix = 0.5,
    condition_weight = c('equal', 'cell_count'),
    coefficient_mask = NULL,
    max_iter = 5000L,
    tol_objective = 1e-7,
    tol_coef = 1e-6,
    keep_history = FALSE
) {
    condition_weight <- match.arg(condition_weight)
    lambda <- sort(unique(as.numeric(lambda)), decreasing = TRUE)
    if (length(lambda) == 0L || any(!is.finite(lambda)) || any(lambda < 0)) {
        stop('lambda must contain finite non-negative values.')
    }
    fits <- vector('list', length(lambda))
    initial_B <- NULL
    initial_step <- NULL
    for (index in seq_along(lambda)) {
        fits[[index]] <- .condition_fit_multitask_lambda(
            X_list = X_list,
            y_list = y_list,
            lambda = lambda[[index]],
            alpha = alpha,
            condition_mix = condition_mix,
            condition_weight = condition_weight,
            coefficient_mask = coefficient_mask,
            initial_B = initial_B,
            initial_step = initial_step,
            max_iter = max_iter,
            tol_objective = tol_objective,
            tol_coef = tol_coef,
            keep_history = keep_history
        )
        initial_B <- fits[[index]]$beta
        initial_step <- fits[[index]]$step
    }
    list(lambda = lambda, fits = fits)
}

.condition_native_solver_available <- function() {
    is.loaded('_Pando_condition_fit_multitask_path_cpp', PACKAGE = 'Pando')
}

.condition_fit_multitask_path_core <- function(
    X_list,
    y_list,
    lambda,
    alpha,
    condition_mix,
    condition_weight,
    coefficient_mask,
    max_iter,
    tol_objective,
    tol_coef,
    keep_history,
    backend = getOption('Pando.condition_solver', 'auto')
) {
    backend <- match.arg(backend, c('auto', 'cpp', 'R'))
    use_cpp <- backend != 'R' && .condition_native_solver_available()
    if (backend == 'cpp' && !use_cpp) {
        stop('The compiled condition-GRN solver is not loaded.')
    }
    if (use_cpp) {
        answer <- .condition_fit_multitask_path_cpp(
            X_list,
            y_list,
            lambda,
            alpha,
            condition_mix,
            condition_weight,
            coefficient_mask,
            as.integer(max_iter),
            tol_objective,
            tol_coef,
            keep_history
        )
        predictor_names <- colnames(X_list[[1L]])
        conditions <- names(X_list)
        answer$fits <- lapply(answer$fits, function(fit) {
            dimnames(fit$beta) <- list(predictor_names, conditions)
            names(fit$intercept) <- conditions
            fit$backend <- 'cpp_eigen_sparse_fista'
            fit
        })
        return(answer)
    }
    answer <- .condition_fit_multitask_path_reference(
        X_list,
        y_list,
        lambda,
        alpha,
        condition_mix,
        condition_weight,
        coefficient_mask,
        max_iter,
        tol_objective,
        tol_coef,
        keep_history
    )
    predictor_names <- colnames(X_list[[1L]])
    conditions <- names(X_list)
    answer$fits <- lapply(answer$fits, function(fit) {
        dimnames(fit$beta) <- list(predictor_names, conditions)
        names(fit$intercept) <- conditions
        fit$backend <- 'R_reference_fista'
        fit
    })
    answer
}

.condition_fit_multitask_path <- function(
    X_list,
    y_list,
    lambda,
    alpha = 0.5,
    condition_mix = 0.5,
    condition_weight = c('equal', 'cell_count'),
    coefficient_mask = NULL,
    max_iter = 5000L,
    tol_objective = 1e-7,
    tol_coef = 1e-6,
    keep_history = FALSE,
    backend = getOption('Pando.condition_solver', 'auto'),
    verified_estimability_mask = NULL
) {
    condition_weight <- match.arg(condition_weight)
    if (is.null(names(X_list))) names(X_list) <- names(y_list)
    if (is.null(names(y_list))) names(y_list) <- names(X_list)
    conditions <- names(X_list)
    predictor_names <- colnames(X_list[[1L]])
    lambda <- sort(unique(as.numeric(lambda)), decreasing = TRUE)
    actual <- if (is.null(verified_estimability_mask)) {
        .condition_true_variance_mask(X_list, coefficient_mask)
    } else {
        verified <- as.matrix(verified_estimability_mask)
        if (!is.logical(verified) || anyNA(verified) ||
            !identical(dim(verified), c(length(predictor_names), length(conditions)))) {
            stop('verified_estimability_mask is not aligned with the path inputs.')
        }
        verified
    }
    keep <- rowSums(actual) > 0L
    if (!any(keep)) {
        loss_weights <- .condition_loss_weights(X_list, condition_weight)
        fits <- lapply(lambda, function(value) {
            intercept <- vapply(y_list, mean, numeric(1))
            objective <- sum(vapply(seq_along(y_list), function(task) {
                residual <- as.numeric(y_list[[task]]) - intercept[[task]]
                0.5 * loss_weights[[task]] * sum(residual * residual)
            }, numeric(1)))
            list(
                beta = matrix(
                    0, length(predictor_names), length(conditions),
                    dimnames = list(predictor_names, conditions)
                ),
                intercept = intercept,
                lambda = value,
                objective = objective,
                objective_change = 0,
                coef_change = 0,
                iterations = 0L,
                converged = TRUE,
                step = 1,
                history = if (keep_history) objective else NULL,
                backend = 'intercept_only'
            )
        })
        return(list(lambda = lambda, fits = fits))
    }
    X_core <- if (all(keep)) X_list else {
        lapply(X_list, function(x) x[, keep, drop = FALSE])
    }
    answer <- .condition_fit_multitask_path_core(
        X_core,
        y_list,
        lambda,
        alpha,
        condition_mix,
        condition_weight,
        actual[keep, , drop = FALSE],
        max_iter,
        tol_objective,
        tol_coef,
        keep_history,
        backend
    )
    if (!all(keep)) {
        answer$fits <- lapply(
            answer$fits,
            .condition_expand_path_fit,
            keep = keep,
            predictor_names = predictor_names,
            conditions = conditions
        )
    }
    answer
}

.condition_make_within_cell_type_folds <- function(
    y_list, nfolds = 5L, seed = 12345L
) {
    if (nfolds < 2L) {
        stop('nfolds must be at least 2.')
    }
    minimum_n <- min(vapply(y_list, length, integer(1)))
    if (minimum_n < nfolds) {
        stop(
            'Every condition in the current cell type must contain at least ',
            'nfolds cells.'
        )
    }
    set.seed(seed)
    folds <- lapply(y_list, function(y) {
        sample(rep(seq_len(nfolds), length.out = length(y)))
    })
    list(folds = folds, effective_nfolds = as.integer(nfolds))
}

.condition_crossfit_within_cell_type <- function(
    X_list,
    y_list,
    lambda,
    alpha = 0.5,
    condition_mix = 0.5,
    condition_weight = c('equal', 'cell_count'),
    coefficient_mask = NULL,
    nfolds = 5L,
    standardize = FALSE,
    active_tol = 1e-8,
    lambda_selection = c('lambda.1se', 'lambda.min'),
    seed = 12345L,
    max_iter = 5000L,
    tol_objective = 1e-7,
    tol_coef = 1e-6
) {
    if (!isTRUE(standardize)) {
        stop('Comparable nested cross-fitting requires standardize = TRUE.')
    }
    .condition_nested_crossfit_within_cell_type(
        X_list = X_list,
        y_list = y_list,
        lambda = lambda,
        alpha = alpha,
        condition_mix = condition_mix,
        condition_weight = match.arg(condition_weight),
        coefficient_mask = coefficient_mask,
        outer_nfolds = nfolds,
        inner_nfolds = nfolds,
        active_tol = active_tol,
        lambda_selection = match.arg(lambda_selection),
        seed = seed,
        max_iter = max_iter,
        tol_objective = tol_objective,
        tol_coef = tol_coef
    )
}
