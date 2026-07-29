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

.condition_sample_block_status <- function(block_list) {
    if (is.null(block_list)) {
        return(list(
            available = FALSE,
            n_blocks = NULL,
            reason = 'cv_block_not_supplied'
        ))
    }
    if (!is.list(block_list) || !length(block_list) ||
        is.null(names(block_list)) || any(names(block_list) == '')) {
        stop('block_list must be a named non-empty list.')
    }
    n_blocks <- vapply(block_list, function(x) {
        x <- trimws(as.character(x))
        if (!length(x) || anyNA(x) || any(x == '')) {
            stop('Every condition must have complete non-empty CV blocks.')
        }
        length(unique(x))
    }, integer(1))
    available <- all(n_blocks >= 2L)
    list(
        available = available,
        n_blocks = n_blocks,
        reason = if (available) {
            'at_least_two_biological_samples_per_condition'
        } else {
            'one_or_more_conditions_have_fewer_than_two_biological_samples'
        }
    )
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
    row_scale[nonzero] <- pmax(1 - group_threshold[nonzero] / row_norm[nonzero], 0)
    U <- U * row_scale
    if (!is.null(coefficient_mask)) {
        U[!coefficient_mask] <- 0
    }
    U
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
    if (is.null(coefficient_mask)) {
        coefficient_mask <- matrix(TRUE, p, n_tasks)
    }
    coefficient_mask <- as.matrix(coefficient_mask)
    if (!identical(dim(coefficient_mask), c(p, n_tasks)) ||
        !is.logical(coefficient_mask) || anyNA(coefficient_mask) ||
        any(rowSums(coefficient_mask) == 0L)) {
        stop('coefficient_mask must be a logical predictors-by-conditions matrix with one eligible condition per predictor.')
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
                element_factor = element_factor,
                coefficient_mask = coefficient_mask
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
                    element_factor = element_factor,
                    coefficient_mask = coefficient_mask
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
        restart <- sum((Z - B) * (B - B_previous)) > 0
        if (restart) {
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
    if (!is.finite(lambda_min_ratio) || lambda_min_ratio <= 0 || lambda_min_ratio >= 1) {
        stop('lambda_min_ratio must be between 0 and 1.')
    }
    lambda_max <- .condition_lambda_max(
        X_list, y_list, alpha, condition_mix, condition_weight,
        coefficient_mask = coefficient_mask
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

.condition_make_folds <- function(
    y_list, nfolds = 5L, seed = 12345L, block_list = NULL
) {
    minimum_n <- min(vapply(y_list, length, integer(1)))
    if (nfolds < 2L) {
        stop('nfolds must be at least 2.')
    }
    set.seed(seed)
    if (is.null(block_list)) {
        if (minimum_n < nfolds) {
            stop('Every condition must contain at least nfolds observations.')
        }
        folds <- lapply(y_list, function(y) {
            sample(rep(seq_len(nfolds), length.out = length(y)))
        })
        return(list(
            folds = folds,
            effective_nfolds = as.integer(nfolds),
            block_to_fold = NULL
        ))
    }
    if (!is.list(block_list) || length(block_list) != length(y_list) ||
        any(vapply(seq_along(y_list), function(i) {
            length(block_list[[i]]) != length(y_list[[i]]) ||
                anyNA(block_list[[i]]) ||
                any(!nzchar(trimws(as.character(block_list[[i]]))))
        }, logical(1)))) {
        stop('block_list must provide one complete block label per observation.')
    }
    block_list <- lapply(block_list, function(x) trimws(as.character(x)))
    minimum_blocks <- min(vapply(block_list, function(x) {
        length(unique(x))
    }, integer(1)))
    if (minimum_blocks < 2L) {
        stop('Every condition must contain at least two biological CV blocks.')
    }
    effective_nfolds <- min(as.integer(nfolds), minimum_blocks)
    all_blocks <- unique(unlist(block_list, use.names = FALSE))
    membership <- vapply(all_blocks, function(block) {
        sum(vapply(block_list, function(x) block %in% x, logical(1)))
    }, integer(1))
    size <- vapply(all_blocks, function(block) {
        sum(vapply(block_list, function(x) sum(x == block), integer(1)))
    }, integer(1))
    order_blocks <- order(-membership, -size, runif(length(all_blocks)))
    assignment <- stats::setNames(integer(length(all_blocks)), all_blocks)
    condition_load <- matrix(
        0L, nrow = length(block_list), ncol = effective_nfolds
    )
    total_load <- integer(effective_nfolds)
    for (block in all_blocks[order_blocks]) {
        tasks <- which(vapply(block_list, function(x) block %in% x, logical(1)))
        score <- colSums(condition_load[tasks, , drop = FALSE])
        candidates <- which(score == min(score))
        candidates <- candidates[total_load[candidates] ==
            min(total_load[candidates])]
        fold <- sample(candidates, 1L)
        assignment[[block]] <- fold
        for (task in tasks) {
            condition_load[task, fold] <- condition_load[task, fold] +
                sum(block_list[[task]] == block)
        }
        total_load[[fold]] <- total_load[[fold]] + size[[block]]
    }
    folds <- lapply(block_list, function(x) unname(assignment[x]))
    missing_fold <- vapply(folds, function(x) {
        length(unique(x)) != effective_nfolds
    }, logical(1))
    if (any(missing_fold)) {
        stop(
            'Unable to assign biological blocks so every condition is ',
            'represented in every validation fold.'
        )
    }
    list(
        folds = folds,
        effective_nfolds = effective_nfolds,
        block_to_fold = data.frame(
            block = names(assignment),
            fold = as.integer(assignment),
            stringsAsFactors = FALSE
        )
    )
}

.condition_cv_multitask_path <- function(
    X_list,
    y_list,
    lambda,
    alpha = 0.5,
    condition_mix = 0.5,
    condition_weight = c('equal', 'cell_count'),
    coefficient_mask = NULL,
    nfolds = 5L,
    block_list = NULL,
    standardize = FALSE,
    active_tol = 1e-8,
    lambda_selection = c('lambda.1se', 'lambda.min'),
    seed = 12345L,
    max_iter = 5000L,
    tol_objective = 1e-7,
    tol_coef = 1e-6
) {
    condition_weight <- match.arg(condition_weight)
    lambda_selection <- match.arg(lambda_selection)
    lambda <- sort(unique(as.numeric(lambda)), decreasing = TRUE)
    fold_info <- .condition_make_folds(
        y_list, nfolds, seed, block_list = block_list
    )
    folds <- fold_info$folds
    nfolds <- fold_info$effective_nfolds
    fold_loss <- matrix(NA_real_, nrow = nfolds, ncol = length(lambda))
    oof_path <- lapply(y_list, function(y) {
        matrix(NA_real_, nrow = length(y), ncol = length(lambda))
    })
    fold_transform <- vector('list', nfolds)

    for (fold in seq_len(nfolds)) {
        X_train <- vector('list', length(X_list))
        y_train <- vector('list', length(y_list))
        X_test <- vector('list', length(X_list))
        y_test <- vector('list', length(y_list))
        test_index_list <- vector('list', length(X_list))
        for (task in seq_along(X_list)) {
            test_index <- folds[[task]] == fold
            test_index_list[[task]] <- test_index
            X_train[[task]] <- X_list[[task]][!test_index, , drop = FALSE]
            y_train[[task]] <- y_list[[task]][!test_index]
            X_test[[task]] <- X_list[[task]][test_index, , drop = FALSE]
            y_test[[task]] <- y_list[[task]][test_index]
        }
        if (isTRUE(standardize)) {
            pooled_X_train <- do.call(rbind, X_train)
            predictor_center <- as.numeric(Matrix::colMeans(pooled_X_train))
            predictor_scale <- sqrt(
                .condition_column_variance(pooled_X_train)
            )
            predictor_scale[
                !is.finite(predictor_scale) |
                    predictor_scale <= .Machine$double.eps
            ] <- 1
            pooled_y_train <- unlist(y_train, use.names = FALSE)
            response_center <- mean(pooled_y_train)
            response_scale <- stats::sd(pooled_y_train)
            if (!is.finite(response_scale) ||
                response_scale <= .Machine$double.eps) {
                stop(
                    'A training fold has zero target variance; reduce nfolds ',
                    'or revise the biological CV blocks.'
                )
            }
            X_train <- lapply(X_train, function(x) {
                value <- sweep(
                    as.matrix(x), 2L, predictor_center, '-'
                )
                Matrix::Matrix(
                    sweep(value, 2L, predictor_scale, '/'), sparse = FALSE
                )
            })
            X_test <- lapply(X_test, function(x) {
                value <- sweep(
                    as.matrix(x), 2L, predictor_center, '-'
                )
                Matrix::Matrix(
                    sweep(value, 2L, predictor_scale, '/'), sparse = FALSE
                )
            })
            y_train <- lapply(y_train, function(y) {
                (y - response_center) / response_scale
            })
            y_test_model <- lapply(y_test, function(y) {
                (y - response_center) / response_scale
            })
        } else {
            predictor_center <- rep(0, ncol(X_list[[1L]]))
            predictor_scale <- rep(1, ncol(X_list[[1L]]))
            response_center <- 0
            response_scale <- 1
            y_test_model <- y_test
        }
        fold_transform[[fold]] <- list(
            predictor_center = predictor_center,
            predictor_scale = predictor_scale,
            response_center = response_center,
            response_scale = response_scale,
            training_only = TRUE
        )

        path <- .condition_fit_multitask_path(
            X_list = X_train,
            y_list = y_train,
            lambda = lambda,
            alpha = alpha,
            condition_mix = condition_mix,
            condition_weight = condition_weight,
            coefficient_mask = coefficient_mask,
            max_iter = max_iter,
            tol_objective = tol_objective,
            tol_coef = tol_coef
        )

        for (lambda_index in seq_along(lambda)) {
            selection_fit <- path$fits[[lambda_index]]
            refit <- .condition_refit_shared_baseline(
                X_list = X_train,
                y_list = y_train,
                beta_selection = selection_fit$beta,
                estimability_mask = coefficient_mask,
                ridge = max(
                    selection_fit$lambda * (1 - alpha), 1e-6
                ),
                active_tol = active_tol,
                condition_weight = condition_weight
            )
            task_mse <- numeric(length(X_list))
            for (task in seq_along(X_list)) {
                prediction <- refit$intercept[[task]] + as.numeric(
                    X_test[[task]] %*% refit$beta[, task]
                )
                task_mse[[task]] <- mean(
                    (y_test_model[[task]] - prediction)^2
                )
                oof_path[[task]][
                    test_index_list[[task]], lambda_index
                ] <-
                    prediction * response_scale + response_center
            }
            fold_loss[fold, lambda_index] <-
                .condition_task_validation_loss(
                    task_mse = task_mse,
                    n_validation = vapply(y_test, length, integer(1)),
                    condition_weight = condition_weight
                )
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
        lambda_1se = lambda[[which(cv_mean <= cv_mean[[minimum_index]] + cv_se[[minimum_index]])[[1L]]]],
        oof_prediction = lapply(
            oof_path, function(x) as.numeric(x[, selected_index])
        ),
        oof_fold = folds,
        effective_nfolds = fold_info$effective_nfolds,
        block_to_fold = fold_info$block_to_fold,
        fold_transform = fold_transform,
        oof_model = 'condition_sparse_selection_plus_common_metric_refit'
    )
}
