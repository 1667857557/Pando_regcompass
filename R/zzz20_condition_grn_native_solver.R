# Compiled sparse-group FISTA backend with an R reference fallback.

.condition_fit_multitask_path_reference <-
    .condition_fit_multitask_path_unmodified
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
            X_list = X_list,
            y_list = y_list,
            lambda = lambda,
            alpha = alpha,
            condition_mix = condition_mix,
            condition_weight = condition_weight,
            coefficient_mask = coefficient_mask,
            max_iter = as.integer(max_iter),
            tol_objective = tol_objective,
            tol_coef = tol_coef,
            keep_history = keep_history
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
        X_list = X_list,
        y_list = y_list,
        lambda = lambda,
        alpha = alpha,
        condition_mix = condition_mix,
        condition_weight = condition_weight,
        coefficient_mask = coefficient_mask,
        max_iter = max_iter,
        tol_objective = tol_objective,
        tol_coef = tol_coef,
        keep_history = keep_history
    )
    answer$fits <- lapply(answer$fits, function(fit) {
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
    backend = getOption('Pando.condition_solver', 'auto')
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
                beta = matrix(0, length(predictor_names), length(conditions),
                              dimnames = list(predictor_names, conditions)),
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
    answer <- .condition_fit_multitask_path_core(
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
        keep_history = keep_history,
        backend = backend
    )
    answer$fits <- lapply(answer$fits, .condition_expand_path_fit,
                          keep = keep,
                          predictor_names = predictor_names,
                          conditions = conditions)
    answer
}
