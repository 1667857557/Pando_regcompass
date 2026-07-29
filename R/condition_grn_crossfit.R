# Nested cross-fitting and condition-balanced transformations.

.condition_population_variance <- function(x) {
    if (!ncol(x)) {
        return(numeric())
    }
    center <- as.numeric(Matrix::colMeans(x))
    second_moment <- as.numeric(Matrix::colMeans(x * x))
    pmax(second_moment - center * center, 0)
}

.condition_stable_hash <- function(x) {
    path <- tempfile('pando-condition-hash-')
    on.exit(unlink(path), add = TRUE)
    writeBin(serialize(x, NULL, version = 2L), path)
    unname(tools::md5sum(path))
}

.condition_build_balanced_transform <- function(
    X_list, y_list, condition_weights = NULL
) {
    if (!is.list(X_list) || !length(X_list) ||
        !is.list(y_list) || length(X_list) != length(y_list)) {
        stop('X_list and y_list must be non-empty aligned condition lists.')
    }
    p <- ncol(X_list[[1L]])
    if (any(vapply(X_list, ncol, integer(1)) != p) ||
        any(vapply(seq_along(X_list), function(i) {
            nrow(X_list[[i]]) != length(y_list[[i]]) ||
                nrow(X_list[[i]]) < 2L
        }, logical(1)))) {
        stop('Each condition requires aligned predictors and at least two responses.')
    }
    k <- length(X_list)
    if (is.null(condition_weights)) {
        condition_weights <- rep(1 / k, k)
    }
    condition_weights <- as.numeric(condition_weights)
    if (length(condition_weights) != k ||
        any(!is.finite(condition_weights)) ||
        any(condition_weights <= 0) ||
        abs(sum(condition_weights) - 1) > sqrt(.Machine$double.eps)) {
        stop('condition_weights must be positive and sum to one.')
    }
    equal_weights <- rep(1 / k, k)
    if (!isTRUE(all.equal(
            condition_weights, equal_weights, tolerance = 0
        ))) {
        stop(
            'Comparable condition transforms require exactly equal condition weights.'
        )
    }
    predictor_means <- vapply(
        X_list, function(x) as.numeric(Matrix::colMeans(x)), numeric(p)
    )
    predictor_variances <- vapply(
        X_list, .condition_population_variance, numeric(p)
    )
    if (is.null(dim(predictor_means))) {
        predictor_means <- matrix(predictor_means, nrow = p)
        predictor_variances <- matrix(predictor_variances, nrow = p)
    }
    predictor_center <- as.numeric(predictor_means %*% condition_weights)
    predictor_scale <- sqrt(as.numeric(
        predictor_variances %*% condition_weights
    ))
    response_means <- vapply(y_list, mean, numeric(1))
    response_variances <- vapply(y_list, function(y) {
        mean((as.numeric(y) - mean(y))^2)
    }, numeric(1))
    response_center <- sum(condition_weights * response_means)
    response_scale <- sqrt(sum(condition_weights * response_variances))
    predictor_estimable <- is.finite(predictor_scale) &
        predictor_scale > .Machine$double.eps
    if (!is.finite(response_scale) ||
        response_scale <= .Machine$double.eps) {
        stop('The equal-condition within-condition response variance is zero.')
    }
    safe_scale <- predictor_scale
    safe_scale[!predictor_estimable] <- 1
    names(predictor_center) <- names(predictor_scale) <-
        colnames(X_list[[1L]])
    structure(
        list(
            predictor_center = predictor_center,
            predictor_scale = safe_scale,
            predictor_estimable = predictor_estimable,
            response_center = response_center,
            response_scale = response_scale,
            condition_weights = stats::setNames(
                condition_weights,
                if (is.null(names(X_list))) {
                    as.character(seq_len(k))
                } else {
                    names(X_list)
                }
            ),
            center_policy = 'equal_condition_mean',
            scale_policy = 'equal_condition_within_population_variance',
            transform_policy =
                'equal_condition_center_equal_condition_within_variance_v1',
            training_fold_only = TRUE,
            predictor_center_hash = .condition_stable_hash(
                unname(predictor_center)
            ),
            predictor_scale_hash = .condition_stable_hash(
                unname(safe_scale)
            )
        ),
        class = c('ConditionBalancedTransform', 'list')
    )
}

.condition_apply_balanced_transform <- function(
    X_list, y_list = NULL, transform
) {
    X_scaled <- lapply(X_list, function(x) {
        value <- sweep(
            as.matrix(x), 2L, transform$predictor_center, '-'
        )
        Matrix::Matrix(
            sweep(value, 2L, transform$predictor_scale, '/'),
            sparse = FALSE
        )
    })
    names(X_scaled) <- names(X_list)
    y_scaled <- if (is.null(y_list)) {
        NULL
    } else {
        lapply(y_list, function(y) {
            (as.numeric(y) - transform$response_center) /
                transform$response_scale
        })
    }
    if (!is.null(y_scaled)) names(y_scaled) <- names(y_list)
    list(X = X_scaled, y = y_scaled)
}

.condition_training_estimability <- function(X_list, coefficient_mask) {
    p <- ncol(X_list[[1L]])
    k <- length(X_list)
    if (is.null(coefficient_mask)) {
        coefficient_mask <- matrix(
            TRUE, p, k,
            dimnames = list(colnames(X_list[[1L]]), names(X_list))
        )
    }
    coefficient_mask <- as.matrix(coefficient_mask)
    if (!is.logical(coefficient_mask) ||
        !identical(dim(coefficient_mask), c(p, k)) ||
        anyNA(coefficient_mask)) {
        stop('coefficient_mask must be a logical predictors-by-conditions matrix.')
    }
    variance_mask <- vapply(
        X_list,
        function(x) {
            value <- .condition_population_variance(x)
            is.finite(value) & value > .Machine$double.eps
        },
        logical(p)
    )
    if (is.null(dim(variance_mask))) {
        variance_mask <- matrix(variance_mask, nrow = p)
    }
    dimnames(variance_mask) <- dimnames(coefficient_mask)
    coefficient_mask & variance_mask
}

.condition_select_lambda_nested <- function(
    X_list, y_list, lambda, alpha, condition_mix, condition_weight,
    coefficient_mask, nfolds, active_tol, lambda_selection, seed,
    max_iter, tol_objective, tol_coef
) {
    if (!identical(condition_weight, 'equal')) {
        stop('Nested comparable selection requires condition_weight = "equal".')
    }
    effective_nfolds <- min(
        as.integer(nfolds),
        min(vapply(y_list, length, integer(1)))
    )
    if (effective_nfolds < 2L) {
        stop('Inner cross-validation requires at least two cells per condition.')
    }
    fold_info <- .condition_make_within_cell_type_folds(
        y_list, nfolds = effective_nfolds, seed = seed
    )
    fold_loss <- matrix(
        NA_real_, nrow = effective_nfolds, ncol = length(lambda)
    )
    fold_transform <- vector('list', effective_nfolds)
    for (fold in seq_len(effective_nfolds)) {
        train <- lapply(seq_along(X_list), function(task) {
            fold_info$folds[[task]] != fold
        })
        X_train_raw <- lapply(seq_along(X_list), function(task) {
            X_list[[task]][train[[task]], , drop = FALSE]
        })
        y_train_raw <- lapply(seq_along(y_list), function(task) {
            y_list[[task]][train[[task]]]
        })
        X_valid_raw <- lapply(seq_along(X_list), function(task) {
            X_list[[task]][!train[[task]], , drop = FALSE]
        })
        y_valid_raw <- lapply(seq_along(y_list), function(task) {
            y_list[[task]][!train[[task]]]
        })
        names(X_train_raw) <- names(X_valid_raw) <- names(X_list)
        names(y_train_raw) <- names(y_valid_raw) <- names(y_list)
        transform <- .condition_build_balanced_transform(
            X_train_raw, y_train_raw
        )
        fold_transform[[fold]] <- transform
        train_scaled <- .condition_apply_balanced_transform(
            X_train_raw, y_train_raw, transform
        )
        valid_scaled <- .condition_apply_balanced_transform(
            X_valid_raw, y_valid_raw, transform
        )
        fold_mask <- .condition_training_estimability(
            X_train_raw, coefficient_mask
        )
        keep <- rowSums(fold_mask) > 0L & transform$predictor_estimable
        if (!any(keep)) {
            next
        }
        path <- .condition_fit_multitask_path(
            X_list = lapply(train_scaled$X, function(x) {
                x[, keep, drop = FALSE]
            }),
            y_list = train_scaled$y,
            lambda = lambda,
            alpha = alpha,
            condition_mix = condition_mix,
            condition_weight = 'equal',
            coefficient_mask = fold_mask[keep, , drop = FALSE],
            max_iter = max_iter,
            tol_objective = tol_objective,
            tol_coef = tol_coef
        )
        for (lambda_index in seq_along(lambda)) {
            selection_fit <- path$fits[[lambda_index]]
            refit <- .condition_refit_shared_baseline(
                X_list = lapply(train_scaled$X, function(x) {
                    x[, keep, drop = FALSE]
                }),
                y_list = train_scaled$y,
                beta_selection = selection_fit$beta,
                estimability_mask = fold_mask[keep, , drop = FALSE],
                ridge = max(selection_fit$lambda * (1 - alpha), 1e-6),
                active_tol = active_tol,
                condition_weight = 'equal'
            )
            mse <- vapply(seq_along(X_list), function(task) {
                prediction <- refit$intercept[[task]] + as.numeric(
                    valid_scaled$X[[task]][, keep, drop = FALSE] %*%
                        refit$beta[, task]
                )
                mean((valid_scaled$y[[task]] - prediction)^2)
            }, numeric(1))
            fold_loss[fold, lambda_index] <- mean(mse)
        }
    }
    cv_mean <- colMeans(fold_loss, na.rm = TRUE)
    cv_se <- apply(fold_loss, 2L, function(x) {
        x <- x[is.finite(x)]
        if (length(x) < 2L) 0 else stats::sd(x) / sqrt(length(x))
    })
    if (!any(is.finite(cv_mean))) {
        stop('Inner cross-validation produced no finite validation loss.')
    }
    minimum_index <- which.min(cv_mean)
    selected_index <- if (lambda_selection == 'lambda.min') {
        minimum_index
    } else {
        eligible <- which(
            cv_mean <= cv_mean[[minimum_index]] + cv_se[[minimum_index]]
        )
        eligible[[1L]]
    }
    list(
        selected_index = selected_index,
        selected_lambda = lambda[[selected_index]],
        lambda_min = lambda[[minimum_index]],
        lambda_1se = lambda[[which(
            cv_mean <= cv_mean[[minimum_index]] + cv_se[[minimum_index]]
        )[[1L]]]],
        cv_mean = cv_mean,
        cv_se = cv_se,
        fold_loss = fold_loss,
        fold_transform = fold_transform,
        effective_nfolds = effective_nfolds
    )
}

.condition_nested_crossfit_within_cell_type <- function(
    X_list,
    y_list,
    lambda,
    lambda_auto = FALSE,
    nlambda = length(lambda),
    lambda_min_ratio = NULL,
    alpha = 0.5,
    condition_mix = 0.5,
    condition_weight = 'equal',
    coefficient_mask = NULL,
    outer_nfolds = 5L,
    inner_nfolds = 5L,
    active_tol = 1e-8,
    lambda_selection = c('lambda.1se', 'lambda.min'),
    comparison_conditions = NULL,
    seed = 12345L,
    max_iter = 5000L,
    tol_objective = 1e-7,
    tol_coef = 1e-6
) {
    lambda_selection <- match.arg(lambda_selection)
    if (!identical(condition_weight, 'equal')) {
        stop('Nested comparable cross-fitting requires equal condition weights.')
    }
    if (is.null(names(X_list))) names(X_list) <- names(y_list)
    if (is.null(names(y_list))) names(y_list) <- names(X_list)
    conditions <- names(X_list)
    if (is.null(comparison_conditions)) {
        comparison_conditions <- if (length(conditions) == 2L) {
            conditions
        } else {
            conditions
        }
    }
    comparison_conditions <- unique(as.character(comparison_conditions))
    if (length(comparison_conditions) < 2L ||
        !all(comparison_conditions %in% conditions)) {
        stop('comparison_conditions must contain at least two fitted conditions.')
    }
    p <- ncol(X_list[[1L]])
    if (is.null(coefficient_mask)) {
        coefficient_mask <- matrix(
            TRUE, p, length(X_list),
            dimnames = list(colnames(X_list[[1L]]), conditions)
        )
    }
    coefficient_mask <- as.matrix(coefficient_mask)
    if (!is.logical(coefficient_mask) ||
        !identical(dim(coefficient_mask), c(p, length(X_list))) ||
        anyNA(coefficient_mask)) {
        stop('coefficient_mask must be a logical predictors-by-conditions matrix.')
    }
    lambda <- sort(unique(as.numeric(lambda)), decreasing = TRUE)
    outer <- .condition_make_within_cell_type_folds(
        y_list, nfolds = outer_nfolds, seed = seed
    )
    prediction <- lapply(y_list, function(y) rep(NA_real_, length(y)))
    projection_full <- lapply(y_list, function(y) rep(NA_real_, length(y)))
    projection_common <- lapply(y_list, function(y) rep(NA_real_, length(y)))
    projection_global <- lapply(y_list, function(y) rep(NA_real_, length(y)))
    assignment_count <- lapply(y_list, function(y) integer(length(y)))
    fold_transform <- vector('list', outer$effective_nfolds)
    fold_selected_lambda <- rep(NA_real_, outer$effective_nfolds)
    fold_inner_cv <- vector('list', outer$effective_nfolds)
    fold_support <- vector('list', outer$effective_nfolds)

    for (fold in seq_len(outer$effective_nfolds)) {
        test <- lapply(outer$folds, `==`, fold)
        X_train_raw <- lapply(seq_along(X_list), function(task) {
            X_list[[task]][!test[[task]], , drop = FALSE]
        })
        y_train_raw <- lapply(seq_along(y_list), function(task) {
            y_list[[task]][!test[[task]]]
        })
        X_test_raw <- lapply(seq_along(X_list), function(task) {
            X_list[[task]][test[[task]], , drop = FALSE]
        })
        y_test_raw <- lapply(seq_along(y_list), function(task) {
            y_list[[task]][test[[task]]]
        })
        names(X_train_raw) <- names(X_test_raw) <- conditions
        names(y_train_raw) <- names(y_test_raw) <- conditions
        transform <- .condition_build_balanced_transform(
            X_train_raw, y_train_raw
        )
        train_scaled <- .condition_apply_balanced_transform(
            X_train_raw, y_train_raw, transform
        )
        test_scaled <- .condition_apply_balanced_transform(
            X_test_raw, y_test_raw, transform
        )
        fold_mask <- .condition_training_estimability(
            X_train_raw, coefficient_mask
        )
        keep <- rowSums(fold_mask) > 0L & transform$predictor_estimable
        fold_transform[[fold]] <- transform
        if (!any(keep)) {
            fold_support[[fold]] <- list(
                estimability_mask = fold_mask,
                comparison_conditions = comparison_conditions,
                pairwise_or_requested_common = stats::setNames(
                    logical(), character()
                ),
                global_common = stats::setNames(logical(), character()),
                projection_unavailable_reason =
                    'no_outer_training_estimable_predictor'
            )
            for (task in seq_along(X_list)) {
                intercept_scaled <- mean(train_scaled$y[[task]])
                prediction[[task]][test[[task]]] <-
                    intercept_scaled * transform$response_scale +
                    transform$response_center
                assignment_count[[task]][test[[task]]] <-
                    assignment_count[[task]][test[[task]]] + 1L
            }
            next
        }
        fold_lambda <- if (isTRUE(lambda_auto)) {
            .condition_make_lambda_path(
                X_list = lapply(train_scaled$X, function(x) {
                    x[, keep, drop = FALSE]
                }),
                y_list = train_scaled$y,
                alpha = alpha,
                condition_mix = condition_mix,
                condition_weight = 'equal',
                coefficient_mask = fold_mask[keep, , drop = FALSE],
                nlambda = nlambda,
                lambda_min_ratio = lambda_min_ratio
            )
        } else {
            lambda
        }
        inner <- .condition_select_lambda_nested(
            lapply(X_train_raw, function(x) {
                x[, keep, drop = FALSE]
            }),
            y_train_raw,
            fold_lambda,
            alpha,
            condition_mix,
            condition_weight,
            fold_mask[keep, , drop = FALSE],
            inner_nfolds,
            active_tol,
            lambda_selection,
            .condition_seed_for(paste0('inner-', fold), seed),
            max_iter,
            tol_objective,
            tol_coef
        )
        fold_selected_lambda[[fold]] <- inner$selected_lambda
        fold_inner_cv[[fold]] <- inner
        selected <- .condition_fit_multitask_path(
            X_list = lapply(train_scaled$X, function(x) {
                x[, keep, drop = FALSE]
            }),
            y_list = train_scaled$y,
            lambda = inner$selected_lambda,
            alpha = alpha,
            condition_mix = condition_mix,
            condition_weight = 'equal',
            coefficient_mask = fold_mask[keep, , drop = FALSE],
            max_iter = max_iter,
            tol_objective = tol_objective,
            tol_coef = tol_coef
        )$fits[[1L]]
        refit <- .condition_refit_shared_baseline(
            X_list = lapply(train_scaled$X, function(x) {
                x[, keep, drop = FALSE]
            }),
            y_list = train_scaled$y,
            beta_selection = selected$beta,
            estimability_mask = fold_mask[keep, , drop = FALSE],
            ridge = max(selected$lambda * (1 - alpha), 1e-6),
            active_tol = active_tol,
            condition_weight = 'equal'
        )
        common_pair <- rowSums(
            fold_mask[keep, comparison_conditions, drop = FALSE]
        ) == length(comparison_conditions)
        common_global <- rowSums(
            fold_mask[keep, , drop = FALSE]
        ) == length(conditions)
        fold_support[[fold]] <- list(
            estimability_mask = fold_mask,
            comparison_conditions = comparison_conditions,
            pairwise_or_requested_common = stats::setNames(
                common_pair, rownames(fold_mask)[keep]
            ),
            global_common = stats::setNames(
                common_global, rownames(fold_mask)[keep]
            )
        )
        for (task in seq_along(X_list)) {
            X_test <- test_scaled$X[[task]][, keep, drop = FALSE]
            beta <- refit$beta[, task]
            estimable <- refit$estimability_mask[, task]
            full_score <- if (any(estimable)) {
                as.numeric(X_test[, estimable, drop = FALSE] %*% beta[estimable])
            } else {
                rep(NA_real_, nrow(X_test))
            }
            common_score <- if (any(common_pair)) {
                as.numeric(X_test[, common_pair, drop = FALSE] %*%
                    beta[common_pair])
            } else {
                rep(NA_real_, nrow(X_test))
            }
            global_score <- if (any(common_global)) {
                as.numeric(X_test[, common_global, drop = FALSE] %*%
                    beta[common_global])
            } else {
                rep(NA_real_, nrow(X_test))
            }
            raw_prediction <- (
                refit$intercept[[task]] + full_score
            ) * transform$response_scale + transform$response_center
            prediction[[task]][test[[task]]] <- raw_prediction
            projection_full[[task]][test[[task]]] <- full_score
            projection_common[[task]][test[[task]]] <- common_score
            projection_global[[task]][test[[task]]] <- global_score
            assignment_count[[task]][test[[task]]] <-
                assignment_count[[task]][test[[task]]] + 1L
        }
    }
    if (any(vapply(assignment_count, function(x) any(x != 1L), logical(1)))) {
        stop('Every cell must receive exactly one outer-fold projection.')
    }
    projection_available_fraction <- vapply(projection_common, function(x) {
        mean(is.finite(x))
    }, numeric(1))
    cell_coverage <- vapply(assignment_count, function(x) {
        mean(x == 1L)
    }, numeric(1))
    list(
        oof_prediction = prediction,
        projection_condition_full_oof = projection_full,
        projection_common_oof = projection_common,
        projection_global_common_oof = projection_global,
        oof_fold = outer$folds,
        outer_nfolds = outer$effective_nfolds,
        inner_nfolds = as.integer(inner_nfolds),
        fold_transform = fold_transform,
        fold_selected_lambda = fold_selected_lambda,
        fold_inner_cv = fold_inner_cv,
        fold_support = fold_support,
        oof_assignment_count = assignment_count,
        oof_cell_coverage = cell_coverage,
        oof_projection_available_fraction =
            projection_available_fraction,
        comparison_conditions = comparison_conditions,
        projection_origin = 'outer_condition_stratified_cell_oof',
        projection_used_for_penalty = TRUE,
        full_fit_projection_used_for_penalty = FALSE,
        fold_transform_policy =
            'equal_condition_center_equal_condition_within_variance_v1',
        oof_model =
            'nested_selection_shared_baseline_refit_heldout_projection'
    )
}
