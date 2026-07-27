# Internal numerical routines for condition-aware Pando.

.condition_sqnorm <- function(x) {
    if (inherits(x, 'sparseMatrix')) {
        return(sum(x@x * x@x))
    }
    sum(x * x)
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

.condition_loss_weights <- function(X_list, condition_weight = c('equal', 'cell_count')) {
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
    element_factor = NULL
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
    element_factor = NULL
) {
    if (is.null(group_factor)) {
        group_factor <- rep(1, nrow(V))
    }
    if (is.null(element_factor)) {
        element_factor <- matrix(1, nrow(V), ncol(V))
    }

    element_threshold <- step * lambda * alpha * condition_mix * element_factor
    U <- sign(V) * pmax(abs(V) - element_threshold, 0)

    group_threshold <- step * lambda * alpha * (1 - condition_mix) * group_factor
    row_norm <- sqrt(rowSums(U * U))
    row_scale <- rep(0, length(row_norm))
    nonzero <- row_norm > 0
    row_scale[nonzero] <- pmax(1 - group_threshold[nonzero] / row_norm[nonzero], 0)
    U * row_scale
}

.condition_objective <- function(
    B,
    X_list,
    y_list,
    loss_weights,
    lambda,
    alpha,
    condition_mix,
    group_factor = NULL,
    element_factor = NULL
) {
    smooth <- .condition_profiled_smooth(
        B = B,
        X_list = X_list,
        y_list = y_list,
        loss_weights = loss_weights,
        ridge = lambda * (1 - alpha)
    )
    smooth$value + .condition_sparse_group_penalty(
        B = B,
        lambda = lambda,
        alpha = alpha,
        condition_mix = condition_mix,
        group_factor = group_factor,
        element_factor = element_factor
    )
}

.condition_initial_step <- function(X_list, loss_weights, ridge) {
    upper_bound <- ridge
    for (task in seq_along(X_list)) {
        upper_bound <- upper_bound + loss_weights[[task]] * .condition_sqnorm(X_list[[task]])
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
    if (!is.numeric(lambda) || length(lambda) != 1L || !is.finite(lambda) || lambda < 0) {
        stop('lambda must be one finite non-negative number.')
    }
    if (!is.numeric(alpha) || length(alpha) != 1L || alpha < 0 || alpha > 1) {
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
        if (ncol(X_list[[task]]) != p || nrow(X_list[[task]]) != length(y_list[[task]])) {
            stop('All task matrices must share predictor columns and match their response lengths.')
        }
    }

    if (is.null(group_factor)) {
        group_factor <- rep(1, p)
    }
    if (is.null(element_factor)) {
        element_factor <- matrix(1, p, n_tasks)
    }
    if (!identical(dim(element_factor), c(p, n_tasks))) {
        stop('element_factor must have dimensions predictors by conditions.')
    }

    loss_weights <- .condition_loss_weights(X_list, condition_weight)
    ridge <- lambda * (1 - alpha)
    B <- if (is.null(initial_B)) matrix(0, p, n_tasks) else as.matrix(initial_B)
    if (!identical(dim(B), c(p, n_tasks))) {
        stop('initial_B must have dimensions predictors by conditions.')
    }
    Z <- B
    acceleration <- 1
    step <- if (is.null(initial_step)) {
        .condition_initial_step(X_list, loss_weights, ridge)
    } else {
        initial_step
    }
    objective_previous <- .condition_objective(
        B, X_list, y_list, loss_weights, lambda, alpha, condition_mix,
        group_factor, element_factor
    )
    history <- if (keep_history) objective_previous else NULL
    converged <- FALSE
    coef_change <- Inf
    objective_change <- Inf

    for (iteration in seq_len(max_iter)) {
        smooth_Z <- .condition_profiled_smooth(
            Z, X_list, y_list, loss_weights, ridge
        )

        repeat {
            candidate <- .condition_sparse_group_prox(
                V = Z - step * smooth_Z$gradient,
                step = step,
                lambda = lambda,
                alpha = alpha,
                condition_mix = condition_mix,
                group_factor = group_factor,
                element_factor = element_factor
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
            smooth_Z <- .condition_profiled_smooth(
                Z, X_list, y_list, loss_weights, ridge
            )
            repeat {
                candidate <- .condition_sparse_group_prox(
                    V = Z - step * smooth_Z$gradient,
                    step = step,
                    lambda = lambda,
                    alpha = alpha,
                    condition_mix = condition_mix,
                    group_factor = group_factor,
                    element_factor = element_factor
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
                    stop('Backtracking line search reached min_step after restart.')
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
        if (objective_change < tol_objective && coef_change < tol_coef) {
            converged <- TRUE
            objective_previous <- objective_candidate
            break
        }

        acceleration_new <- (1 + sqrt(1 + 4 * acceleration^2)) / 2
        Z_new <- B + ((acceleration - 1) / acceleration_new) * (B - B_previous)
        if (sum((B - B_previous) * (Z_new - B)) > 0) {
            acceleration_new <- 1
            Z_new <- B
        }
        acceleration <- acceleration_new
        Z <- Z_new
        objective_previous <- objective_candidate
    }

    final_smooth <- .condition_profiled_smooth(
        B, X_list, y_list, loss_weights, ridge
    )
    list(
        beta = B,
        intercept = final_smooth$intercept,
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
    condition_weight = c('equal', 'cell_count')
) {
    condition_weight <- match.arg(condition_weight)
    p <- ncol(X_list[[1L]])
    n_tasks <- length(X_list)
    zero <- matrix(0, p, n_tasks)
    loss_weights <- .condition_loss_weights(X_list, condition_weight)
    gradient <- .condition_profiled_smooth(
        zero, X_list, y_list, loss_weights, ridge = 0
    )$gradient

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
    if (!is.finite(lambda_min_ratio) || lambda_min_ratio <= 0 || lambda_min_ratio >= 1) {
        stop('lambda_min_ratio must be between 0 and 1.')
    }
    lambda_max <- .condition_lambda_max(
        X_list, y_list, alpha, condition_mix, condition_weight
    )
    exp(seq(log(lambda_max), log(lambda_max * lambda_min_ratio), length.out = nlambda))
}

.condition_fit_multitask_path <- function(
    X_list,
    y_list,
    lambda,
    alpha = 0.5,
    condition_mix = 0.5,
    condition_weight = c('equal', 'cell_count'),
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

.condition_make_folds <- function(y_list, nfolds = 5L, seed = 12345L) {
    minimum_n <- min(vapply(y_list, length, integer(1)))
    if (nfolds < 2L) {
        stop('nfolds must be at least 2.')
    }
    if (minimum_n < nfolds) {
        stop('Every condition must contain at least nfolds observations.')
    }
    set.seed(seed)
    lapply(y_list, function(y) {
        sample(rep(seq_len(nfolds), length.out = length(y)))
    })
}

.condition_cv_multitask_path <- function(
    X_list,
    y_list,
    lambda,
    alpha = 0.5,
    condition_mix = 0.5,
    condition_weight = c('equal', 'cell_count'),
    nfolds = 5L,
    lambda_selection = c('lambda.1se', 'lambda.min'),
    seed = 12345L,
    max_iter = 5000L,
    tol_objective = 1e-7,
    tol_coef = 1e-6
) {
    condition_weight <- match.arg(condition_weight)
    lambda_selection <- match.arg(lambda_selection)
    lambda <- sort(unique(as.numeric(lambda)), decreasing = TRUE)
    folds <- .condition_make_folds(y_list, nfolds, seed)
    fold_loss <- matrix(NA_real_, nrow = nfolds, ncol = length(lambda))

    for (fold in seq_len(nfolds)) {
        X_train <- vector('list', length(X_list))
        y_train <- vector('list', length(y_list))
        X_test <- vector('list', length(X_list))
        y_test <- vector('list', length(y_list))
        for (task in seq_along(X_list)) {
            test_index <- folds[[task]] == fold
            X_train[[task]] <- X_list[[task]][!test_index, , drop = FALSE]
            y_train[[task]] <- y_list[[task]][!test_index]
            X_test[[task]] <- X_list[[task]][test_index, , drop = FALSE]
            y_test[[task]] <- y_list[[task]][test_index]
        }

        path <- .condition_fit_multitask_path(
            X_list = X_train,
            y_list = y_train,
            lambda = lambda,
            alpha = alpha,
            condition_mix = condition_mix,
            condition_weight = condition_weight,
            max_iter = max_iter,
            tol_objective = tol_objective,
            tol_coef = tol_coef
        )

        for (lambda_index in seq_along(lambda)) {
            fit <- path$fits[[lambda_index]]
            task_mse <- numeric(length(X_list))
            for (task in seq_along(X_list)) {
                prediction <- fit$intercept[[task]] + as.numeric(
                    X_test[[task]] %*% fit$beta[, task]
                )
                task_mse[[task]] <- mean((y_test[[task]] - prediction)^2)
            }
            fold_loss[fold, lambda_index] <- mean(task_mse)
        }
    }

    cv_mean <- colMeans(fold_loss)
    cv_se <- apply(fold_loss, 2, stats::sd) / sqrt(nfolds)
    minimum_index <- which.min(cv_mean)
    if (lambda_selection == 'lambda.min') {
        selected_index <- minimum_index
    } else {
        threshold <- cv_mean[[minimum_index]] + cv_se[[minimum_index]]
        eligible <- which(cv_mean <= threshold)
        selected_index <- eligible[[1L]]
    }

    list(
        lambda = lambda,
        fold_loss = fold_loss,
        cv_mean = cv_mean,
        cv_se = cv_se,
        selected_index = selected_index,
        selected_lambda = lambda[[selected_index]],
        lambda_min = lambda[[minimum_index]],
        lambda_1se = lambda[[which(cv_mean <= cv_mean[[minimum_index]] + cv_se[[minimum_index]])[[1L]]]]
    )
}
