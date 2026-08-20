# Common-dictionary conditional GRN entry points.
#
# The existing condition-GRN path is implemented directly as E-star/JSE with
# the fixed production threshold z = 0.25. Predictor construction and scaling
# stay on the exact TF-peak-target dictionary used by every condition.

.condition_multitask_ridge_schema <- "pando_condition_grn_Estar_jointse_v1"

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
    control = list(), reference_condition = NULL) {
    .condition_scheme_e_fit(
        x = x,
        y = y,
        scaling = scaling,
        min_residual_df = min_residual_df,
        inference = inference,
        control = control,
        reference_condition = reference_condition
    )
}
