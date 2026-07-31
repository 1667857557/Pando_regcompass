# Sparse equal-condition transformation without cell-by-edge densification.

.condition_scale_columns_sparse <- function(x, scale) {
    scale <- as.numeric(scale)
    if (length(scale) != ncol(x) || any(!is.finite(scale)) || any(scale <= 0)) {
        stop('Predictor scale must be finite, positive and aligned to columns.')
    }
    if (inherits(x, 'sparseMatrix')) {
        value <- x %*% Matrix::Diagonal(x = 1 / scale)
        dimnames(value) <- dimnames(x)
        return(methods::as(value, 'dgCMatrix'))
    }
    value <- sweep(as.matrix(x), 2L, scale, '/')
    dimnames(value) <- dimnames(x)
    value
}

.condition_apply_balanced_transform <- function(
    X_list, y_list = NULL, transform
) {
    X_scaled <- lapply(X_list, function(x) {
        .condition_scale_columns_sparse(x, transform$predictor_scale)
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
    list(
        X = X_scaled,
        y = y_scaled,
        predictor_center_implementation =
            'implicit_in_condition_intercept_and_projection_shift'
    )
}

.condition_centered_sqnorm <- function(x) {
    n <- nrow(x)
    if (!n || !ncol(x)) return(0)
    raw <- if (inherits(x, 'sparseMatrix')) {
        sum(x@x * x@x)
    } else {
        sum(x * x)
    }
    mu <- as.numeric(Matrix::colMeans(x))
    max(raw - n * sum(mu * mu), 0)
}

.condition_initial_step <- function(X_list, loss_weights, ridge) {
    upper_bound <- ridge
    for (task in seq_along(X_list)) {
        upper_bound <- upper_bound +
            loss_weights[[task]] * .condition_centered_sqnorm(X_list[[task]])
    }
    if (!is.finite(upper_bound) || upper_bound <= 0) return(1)
    1 / upper_bound
}
