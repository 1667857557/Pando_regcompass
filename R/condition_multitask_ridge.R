# Common-dictionary conditional GRN solver entry points.
#
# Conditional fits use Scheme E exact-edge sparse deviations with z = 0.25.
# The public condition-GRN API is unchanged; `condition_ridge_control` remains
# the compatibility name for numerical solver controls only.

.condition_multitask_ridge_schema <- "pando_condition_grn_sparse_deviation_v4"

.condition_ridge_control <- function(control = list()) {
    .condition_scheme_e_control(control)
}

.condition_ridge_predictors <- function(prepared, edges, cells_by_condition) {
    out <- lapply(cells_by_condition, function(cells) {
        x <- vapply(seq_len(nrow(edges)), function(j) {
            as.numeric(prepared$gene_data[cells, edges$tf[[j]]]) *
                as.numeric(prepared$peak_data[cells, edges$region[[j]]])
        }, numeric(length(cells)))
        if (is.null(dim(x))) x <- matrix(x, ncol = 1L)
        colnames(x) <- edges$edge_id
        x
    })
    names(out) <- names(cells_by_condition)
    out
}

.condition_ridge_scaling <- function(x, floor) {
    columns <- colnames(x[[1L]])
    if (any(vapply(x, function(one) {
        !is.matrix(one) || !identical(colnames(one), columns) ||
            any(!is.finite(one))
    }, logical(1)))) {
        stop("Conditions do not share one finite ordered predictor dictionary.",
             call. = FALSE)
    }

    # Condition-specific intercepts remove between-condition mean shifts. Every
    # exact edge then receives one common scale computed from the equal-condition
    # RMS of its within-condition variation. Cell number remains in X'X and is
    # deliberately not cancelled by a condition-size weight.
    pooled <- do.call(rbind, x)
    center <- colMeans(pooled)
    within_variance <- vapply(seq_len(ncol(x[[1L]])), function(j) {
        mean(vapply(x, function(one) {
            value <- as.numeric(one[, j])
            mean((value - mean(value))^2)
        }, numeric(1)))
    }, numeric(1))
    names(within_variance) <- columns
    scale <- sqrt(within_variance)
    informative <- is.finite(scale) & scale > floor
    scale[!informative] <- 1
    list(
        center = center,
        scale = scale,
        informative = informative,
        floor = floor,
        reference = "equal_condition_within_condition_rms"
    )
}

.condition_ridge_kappa <- function(x) {
    if (!is.matrix(x)) x <- as.matrix(x)
    if (!nrow(x) || !ncol(x) || ncol(x) > 2000L) return(NA_real_)
    tryCatch(as.numeric(kappa(x, exact = FALSE)),
             error = function(e) NA_real_)
}

.condition_ridge_fit <- function(
    x, y, scaling, min_residual_df = 1L, inference = TRUE,
    control = list()) {
    .condition_scheme_e_fit(
        x = x,
        y = y,
        scaling = scaling,
        min_residual_df = min_residual_df,
        inference = inference,
        control = control
    )
}
