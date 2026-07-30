# Condition-full OOF projection with explicit projectable structural zeros.
# Exact-zero predictors are retained in the shared candidate supergraph, are
# excluded from coefficient estimation, and contribute deterministically zero.

.condition_population_variance_unmodified <- .condition_population_variance
.condition_column_variance_unmodified <- .condition_column_variance
.condition_make_lambda_path_unmodified <- .condition_make_lambda_path
.condition_fit_multitask_path_unmodified <- .condition_fit_multitask_path
.condition_select_lambda_nested_unmodified <- .condition_select_lambda_nested
.condition_refit_shared_baseline_unmodified <-
    .condition_refit_shared_baseline
.condition_combine_fit_contracts_unmodified <-
    .condition_combine_fit_contracts

.condition_exact_zero_columns <- function(x) {
    if (!ncol(x)) return(logical())
    as.numeric(Matrix::colSums(abs(x))) == 0
}

.condition_true_variance_mask <- function(X_list, coefficient_mask = NULL) {
    p <- ncol(X_list[[1L]])
    k <- length(X_list)
    if (is.null(coefficient_mask)) {
        coefficient_mask <- matrix(
            TRUE, p, k,
            dimnames = list(colnames(X_list[[1L]]), names(X_list))
        )
    }
    coefficient_mask <- as.matrix(coefficient_mask)
    variance_mask <- vapply(
        X_list,
        function(x) {
            value <- .condition_population_variance_unmodified(x)
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

# Keep exact-zero candidates in the shared supergraph. A small sentinel is used
# only by candidate-retention checks; coefficient estimability always uses the
# unmodified variance functions above.
.condition_population_variance <- function(x) {
    value <- .condition_population_variance_unmodified(x)
    zero <- .condition_exact_zero_columns(x)
    value[zero] <- 4 * .Machine$double.eps
    value
}

.condition_column_variance <- function(x) {
    value <- .condition_column_variance_unmodified(x)
    zero <- .condition_exact_zero_columns(x)
    value[zero] <- 4 * .Machine$double.eps
    value
}

.condition_training_estimability <- function(X_list, coefficient_mask) {
    .condition_true_variance_mask(X_list, coefficient_mask)
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
    actual <- .condition_true_variance_mask(X_list, coefficient_mask)
    keep <- rowSums(actual) > 0L
    if (!any(keep)) {
        if (is.null(lambda_min_ratio)) lambda_min_ratio <- 1e-4
        return(exp(seq(
            log(1), log(lambda_min_ratio), length.out = as.integer(nlambda)
        )))
    }
    .condition_make_lambda_path_unmodified(
        X_list = lapply(X_list, function(x) x[, keep, drop = FALSE]),
        y_list = y_list,
        alpha = alpha,
        condition_mix = condition_mix,
        condition_weight = condition_weight,
        coefficient_mask = actual[keep, , drop = FALSE],
        nlambda = nlambda,
        lambda_min_ratio = lambda_min_ratio
    )
}

.condition_expand_path_fit <- function(fit, keep, predictor_names, conditions) {
    beta <- matrix(
        0, length(keep), length(conditions),
        dimnames = list(predictor_names, conditions)
    )
    beta[keep, ] <- fit$beta
    fit$beta <- beta
    fit
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
    if (is.null(names(X_list))) names(X_list) <- names(y_list)
    if (is.null(names(y_list))) names(y_list) <- names(X_list)
    conditions <- names(X_list)
    predictor_names <- colnames(X_list[[1L]])
    lambda <- sort(unique(as.numeric(lambda)), decreasing = TRUE)
    actual <- .condition_true_variance_mask(X_list, coefficient_mask)
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
                history = if (keep_history) objective else NULL
            )
        })
        return(list(lambda = lambda, fits = fits))
    }
    answer <- .condition_fit_multitask_path_unmodified(
        X_list = lapply(X_list, function(x) x[, keep, drop = FALSE]),
        y_list = y_list,
        lambda = lambda,
        alpha = alpha,
        condition_mix = condition_mix,
        condition_weight = condition_weight,
        coefficient_mask = actual[keep, , drop = FALSE],
        max_iter = max_iter,
        tol_objective = tol_objective,
        tol_coef = tol_coef,
        keep_history = keep_history
    )
    answer$fits <- lapply(
        answer$fits,
        .condition_expand_path_fit,
        keep = keep,
        predictor_names = predictor_names,
        conditions = conditions
    )
    answer
}

.condition_select_lambda_nested <- function(
    X_list, y_list, lambda, alpha, condition_mix, condition_weight,
    coefficient_mask, nfolds, active_tol, lambda_selection, seed,
    max_iter, tol_objective, tol_coef
) {
    actual <- .condition_true_variance_mask(X_list, coefficient_mask)
    if (any(rowSums(actual) > 0L)) {
        return(.condition_select_lambda_nested_unmodified(
            X_list, y_list, lambda, alpha, condition_mix,
            condition_weight, coefficient_mask, nfolds, active_tol,
            lambda_selection, seed, max_iter, tol_objective, tol_coef
        ))
    }
    lambda <- sort(unique(as.numeric(lambda)), decreasing = TRUE)
    effective_nfolds <- min(
        as.integer(nfolds), min(vapply(y_list, length, integer(1)))
    )
    list(
        selected_index = 1L,
        selected_lambda = lambda[[1L]],
        lambda_min = lambda[[1L]],
        lambda_1se = lambda[[1L]],
        cv_mean = rep(NA_real_, length(lambda)),
        cv_se = rep(NA_real_, length(lambda)),
        fold_loss = matrix(
            NA_real_, nrow = effective_nfolds, ncol = length(lambda)
        ),
        fold_transform = vector('list', effective_nfolds),
        effective_nfolds = effective_nfolds,
        selection_reason = 'intercept_only_all_predictors_structural_zero'
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
    tol = 1e-8
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
            beta_shared = stats::setNames(
                rep(0, length(predictor_names)), predictor_names
            ),
            delta_condition = unavailable,
            support_mask = matrix(
                FALSE, length(predictor_names), length(conditions),
                dimnames = dimnames(beta)
            ),
            active_mask = matrix(
                FALSE, length(predictor_names), length(conditions),
                dimnames = dimnames(beta)
            ),
            estimability_mask = actual,
            intercept = stats::setNames(
                vapply(y_list, mean, numeric(1)), conditions
            ),
            ridge = ridge,
            common_metric = 'pooled_weighted_predictor_gram',
            iterations = 0L,
            coef_change = 0,
            converged = TRUE
        ))
    }
    fitted <- .condition_refit_shared_baseline_unmodified(
        X_list = lapply(X_list, function(x) x[, keep, drop = FALSE]),
        y_list = y_list,
        beta_selection = beta_selection[keep, , drop = FALSE],
        estimability_mask = actual[keep, , drop = FALSE],
        ridge = ridge,
        active_tol = active_tol,
        condition_weight = condition_weight,
        max_iter = max_iter,
        tol = tol
    )
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
    beta_shared <- stats::setNames(
        rep(0, length(predictor_names)), predictor_names
    )
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
    if (is.null(comparison_conditions)) comparison_conditions <- conditions
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
        names(X_train_raw) <- names(X_test_raw) <- conditions
        names(y_train_raw) <- conditions
        transform <- .condition_build_balanced_transform(
            X_train_raw, y_train_raw
        )
        train_scaled <- .condition_apply_balanced_transform(
            X_train_raw, y_train_raw, transform
        )
        test_scaled <- .condition_apply_balanced_transform(
            X_test_raw, transform = transform
        )
        fold_estimable <- .condition_training_estimability(
            X_train_raw, coefficient_mask
        )
        fold_structural_zero <- !fold_estimable
        dimnames(fold_structural_zero) <- dimnames(fold_estimable)
        keep <- rowSums(fold_estimable) > 0L & transform$predictor_estimable
        common_pair_full <- rowSums(
            fold_estimable[, comparison_conditions, drop = FALSE]
        ) == length(comparison_conditions)
        common_global_full <- rowSums(fold_estimable) == length(conditions)
        fold_transform[[fold]] <- transform
        fold_support[[fold]] <- list(
            coefficient_estimable_mask = fold_estimable,
            projectable_structural_zero_mask = fold_structural_zero,
            projection_support_mask = fold_estimable | fold_structural_zero,
            comparison_conditions = comparison_conditions,
            pairwise_or_requested_common = stats::setNames(
                common_pair_full, rownames(fold_estimable)
            ),
            global_common = stats::setNames(
                common_global_full, rownames(fold_estimable)
            )
        )
        if (!any(keep)) {
            fold_support[[fold]]$projection_status <-
                'intercept_only_all_predictors_structural_zero'
            for (task in seq_along(X_list)) {
                intercept_scaled <- mean(train_scaled$y[[task]])
                prediction[[task]][test[[task]]] <-
                    intercept_scaled * transform$response_scale +
                    transform$response_center
                projection_full[[task]][test[[task]]] <- 0
                projection_common[[task]][test[[task]]] <- 0
                projection_global[[task]][test[[task]]] <- 0
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
                coefficient_mask = fold_estimable[keep, , drop = FALSE],
                nlambda = nlambda,
                lambda_min_ratio = lambda_min_ratio
            )
        } else {
            lambda
        }
        inner <- .condition_select_lambda_nested(
            lapply(X_train_raw, function(x) x[, keep, drop = FALSE]),
            y_train_raw,
            fold_lambda,
            alpha,
            condition_mix,
            condition_weight,
            fold_estimable[keep, , drop = FALSE],
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
            coefficient_mask = fold_estimable[keep, , drop = FALSE],
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
            estimability_mask = fold_estimable[keep, , drop = FALSE],
            ridge = max(selected$lambda * (1 - alpha), 1e-6),
            active_tol = active_tol,
            condition_weight = 'equal'
        )
        common_pair <- common_pair_full[keep]
        common_global <- common_global_full[keep]
        for (task in seq_along(X_list)) {
            X_test <- test_scaled$X[[task]][, keep, drop = FALSE]
            beta <- refit$beta[, task]
            estimable <- refit$estimability_mask[, task]
            full_score <- rep(0, nrow(X_test))
            if (any(estimable)) {
                full_score <- as.numeric(
                    X_test[, estimable, drop = FALSE] %*% beta[estimable]
                )
            }
            common_score <- rep(0, nrow(X_test))
            common_estimable <- common_pair & estimable
            if (any(common_estimable)) {
                common_score <- as.numeric(
                    X_test[, common_estimable, drop = FALSE] %*%
                    beta[common_estimable]
                )
            }
            global_score <- rep(0, nrow(X_test))
            global_estimable <- common_global & estimable
            if (any(global_estimable)) {
                global_score <- as.numeric(
                    X_test[, global_estimable, drop = FALSE] %*%
                    beta[global_estimable]
                )
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
    projection_available_fraction <- vapply(projection_full, function(x) {
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
        oof_projection_available_fraction = projection_available_fraction,
        comparison_conditions = comparison_conditions,
        projection_origin = 'outer_condition_stratified_cell_oof',
        primary_projection = 'condition_full_oof',
        common_projection_role = 'shared_estimable_component',
        condition_unique_projection_role =
            'condition_full_oof_minus_common_support_oof',
        nonestimable_projection_policy = 'structural_zero_by_condition',
        projection_used_for_penalty = TRUE,
        full_fit_projection_used_for_penalty = FALSE,
        fold_transform_policy =
            'equal_condition_center_equal_condition_within_variance_v1',
        oof_model =
            'nested_selection_shared_baseline_refit_heldout_condition_full_projection'
    )
}

.condition_combine_fit_contracts <- function(...) {
    fit <- .condition_combine_fit_contracts_unmodified(...)
    fit$contract_version <- 'condition_absolute_oof_v4'
    fit$coefficient_estimable_mask <- fit$estimability_mask
    fit$projectable_structural_zero_mask <-
        fit$topology_mask & !fit$estimability_mask
    fit$projection_support_mask <-
        fit$coefficient_estimable_mask |
        fit$projectable_structural_zero_mask
    fit$projection_used_for_penalty <- TRUE
    fit$primary_projection <- 'projection_condition_full_oof'
    fit$common_projection_role <- 'shared_estimable_component'
    fit$condition_unique_projection_role <-
        'projection_condition_full_oof - projection_common_oof'
    fit$nonestimable_projection_policy <- 'structural_zero_by_condition'
    fit$projection_contract$score <- paste(
        'outer-heldout condition-full sum(z_edge * beta_condition)',
        'with nonestimable edge contribution fixed at zero'
    )
    fit$projection_contract$primary_support_policy <- 'condition_full_oof'
    fit$projection_contract$common_support_role <-
        'shared_estimable_component'
    fit$projection_contract$condition_unique_role <-
        'condition_full_oof_minus_common_support_oof'
    fit$projection_contract$nonestimable <-
        'projectable_structural_zero_by_condition'
    fit$projection_contract$condition_full_role <- 'primary_penalty'
    fit
}

#' Project the primary condition-full OOF regulatory signal
#'
#' The primary score contains every edge estimable in the focal condition.
#' Jointly estimable edges form the common-support component. An edge that is
#' non-estimable in one or both conditions remains in the shared candidate
#' supergraph and contributes exactly zero in each non-estimable condition.
#'
#' @param object GRNData object containing the paired cells used for inference.
#' @param fit Optional canonical `pando_condition_grn_fit` object.
#' @param network_name,cell_type Optional fit filters when `fit` is omitted.
#' @param scale Standardized or raw coefficient scale.
#' @param targets Optional target subset.
#' @param nonestimable Use structural zeros or stop on unavailable target scores.
#' @param active_tol Activity threshold used in metadata.
#' @return A primary `ConditionGRNProjection`.
#' @export
project_condition_grn_primary_cells <- function(
    object,
    fit = NULL,
    network_name = NULL,
    cell_type = NULL,
    scale = c('std', 'raw'),
    targets = NULL,
    nonestimable = c('structural_zero', 'error'),
    active_tol = if (is.null(fit)) 1e-8 else fit$active_tol
) {
    scale <- match.arg(scale)
    nonestimable <- match.arg(nonestimable)
    projection <- project_condition_grn_cells(
        object = object,
        fit = fit,
        network_name = network_name,
        cell_type = cell_type,
        component = 'condition',
        scale = scale,
        output = 'gene_score',
        targets = targets,
        nonestimable = nonestimable,
        support_policy = 'condition_estimable',
        origin = 'oof',
        diagnostic_only = TRUE,
        active_tol = active_tol
    )
    projection$schema_version <-
        'pando_condition_grn_primary_projection_v1'
    projection$projection_used_for_penalty <- TRUE
    projection$projection_role <- 'primary_penalty'
    projection$score_comparability_class <-
        'condition_full_oof_on_shared_celltype_coordinate'
    projection$primary_support_policy <- 'condition_full_oof'
    projection$common_support_role <- 'shared_estimable_component'
    projection$condition_unique_role <-
        'condition_full_oof_minus_common_support_oof'
    projection$nonestimable_policy <- 'structural_zero_by_condition'
    projection$aggregation_contract$projection_role <- 'primary_penalty'
    projection$aggregation_contract$primary_support_policy <-
        'condition_full_oof'
    projection
}
