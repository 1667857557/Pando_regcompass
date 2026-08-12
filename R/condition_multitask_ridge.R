# Common-lambda no-fusion ridge solver for condition GRNs.
#
# Every condition uses the same ordered TF-by-ATAC predictor dictionary, the
# same predictor scaling convention and the same target-specific CV lambda. The
# coefficient blocks are otherwise independent: no penalty couples beta values
# across conditions.

.condition_multitask_ridge_schema <- "pando_condition_grn_multitask_ridge_v3"

.condition_ridge_control <- function(control = list()) {
    if (is.null(control)) control <- list()
    if (!is.list(control)) {
        stop("`ridge_control` must be a list.", call. = FALSE)
    }
    defaults <- list(
        lambda_grid = 10^seq(-3, 2, length.out = 9L),
        lambda_rule = "1se",
        cv_folds = 5L,
        seed = 1L,
        scale_floor = 1e-8
    )
    unknown <- setdiff(names(control), names(defaults))
    if (length(unknown)) {
        stop("Unknown `ridge_control` field(s): ",
             paste(unknown, collapse = ", "), call. = FALSE)
    }
    out <- utils::modifyList(defaults, control)
    out$lambda_grid <- sort(unique(as.numeric(out$lambda_grid)))
    if (!length(out$lambda_grid) || any(!is.finite(out$lambda_grid)) ||
        any(out$lambda_grid <= 0)) {
        stop("`ridge_control$lambda_grid` must contain positive finite values.",
             call. = FALSE)
    }
    out$lambda_rule <- match.arg(out$lambda_rule, c("1se", "min"))
    if (!is.numeric(out$cv_folds) || length(out$cv_folds) != 1L ||
        !is.finite(out$cv_folds) || out$cv_folds < 2L ||
        out$cv_folds != as.integer(out$cv_folds)) {
        stop("`ridge_control$cv_folds` must be an integer >= 2.",
             call. = FALSE)
    }
    if (!is.numeric(out$seed) || length(out$seed) != 1L ||
        !is.finite(out$seed) || out$seed != as.integer(out$seed)) {
        stop("`ridge_control$seed` must be a finite integer.", call. = FALSE)
    }
    if (!is.numeric(out$scale_floor) || length(out$scale_floor) != 1L ||
        !is.finite(out$scale_floor) || out$scale_floor <= 0) {
        stop("`ridge_control$scale_floor` must be positive and finite.",
             call. = FALSE)
    }
    out$cv_folds <- as.integer(out$cv_folds)
    out$seed <- as.integer(out$seed)
    out
}

.condition_ridge_folds <- function(cells_by_condition, nfolds, seed) {
    nfolds <- min(as.integer(nfolds), min(lengths(cells_by_condition)))
    if (nfolds < 2L) stop("At least two CV folds are required.", call. = FALSE)
    old <- if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) {
        get(".Random.seed", .GlobalEnv)
    } else NULL
    on.exit({
        if (is.null(old)) {
            if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) {
                rm(".Random.seed", envir = .GlobalEnv)
            }
        } else assign(".Random.seed", old, envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(seed)
    out <- lapply(cells_by_condition, function(cells) {
        order <- sample.int(length(cells))
        fold <- rep(seq_len(nfolds), length.out = length(cells))
        answer <- integer(length(cells))
        answer[order] <- fold
        answer
    })
    names(out) <- names(cells_by_condition)
    attr(out, "nfolds") <- nfolds
    out
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

    # Condition-specific intercepts remove between-condition mean shifts from the
    # slope problem. Scale every shared predictor by the equal-condition RMS of
    # its within-condition SD so the ridge metric matches the equal-condition
    # loss rather than pooled sample abundance or between-condition mean shifts.
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

.condition_ridge_penalty <- function(k, p) {
    diag(k * p)
}

.condition_ridge_kappa <- function(x) {
    if (!is.matrix(x)) x <- as.matrix(x)
    if (!nrow(x) || !ncol(x) || ncol(x) > 2000L) return(NA_real_)
    tryCatch(as.numeric(kappa(x, exact = FALSE)),
             error = function(e) NA_real_)
}

.condition_ridge_fit <- function(
    x, y, scaling, lambda, min_residual_df = 1L, inference = TRUE) {
    conditions <- names(x)
    k <- length(conditions)
    p_full <- ncol(x[[1L]])
    n <- lengths(y)
    n_total <- sum(n)
    if (!identical(names(y), conditions) ||
        any(vapply(seq_along(x), function(i) {
            nrow(x[[i]]) != length(y[[i]]) || any(!is.finite(y[[i]]))
        }, logical(1)))) {
        return(list(status = "nonfinite_or_misaligned_input"))
    }
    if (n_total - k < min_residual_df) {
        return(list(status = "insufficient_df"))
    }
    weight <- n_total / (k * n)
    names(weight) <- conditions
    informative <- as.logical(scaling$informative)

    beta <- matrix(0, k, p_full,
                   dimnames = list(conditions, colnames(x[[1L]])))
    beta_z <- beta
    se <- matrix(NA_real_, k, p_full, dimnames = dimnames(beta))
    statistic <- pval <- se
    zero_variance_values <- vapply(seq_along(x), function(i) {
        value <- apply(x[[i]], 2L, stats::var)
        !is.finite(value) | value <= scaling$floor^2
    }, logical(p_full))
    zero_variance <- matrix(
        as.logical(zero_variance_values),
        nrow = k,
        ncol = p_full,
        byrow = TRUE,
        dimnames = list(conditions, colnames(beta))
    )

    if (!any(informative)) {
        intercept <- vapply(y, mean, numeric(1))
        prediction <- Map(function(value, a) rep(a, length(value)),
                          y, intercept)
        rsq <- vapply(y, function(value) {
            if (stats::var(value) > 0) 0 else NA_real_
        }, numeric(1))
        return(list(
            status = "ok", beta = beta, beta_z = beta_z, se = se,
            statistic = statistic, pval = pval, intercept = intercept,
            prediction = prediction, rsq = rsq,
            raw_rank = setNames(rep(0L, k), conditions),
            effective_df = setNames(rep(0, k), conditions),
            residual_df = n_total - k,
            raw_kappa = setNames(rep(NA_real_, k), conditions),
            regularized_kappa = NA_real_, zero_variance = zero_variance,
            informative = informative, informative_index = integer(),
            covariance_z = NULL, weight = weight
        ))
    }

    keep <- which(informative)
    p <- length(keep)
    gram <- matrix(0, k * p, k * p)
    rhs <- numeric(k * p)
    xc <- yc <- vector("list", k)
    ybar <- numeric(k)
    raw_rank <- integer(k)
    raw_kappa <- numeric(k)
    blocks <- vector("list", k)

    for (i in seq_len(k)) {
        z <- sweep(x[[i]][, keep, drop = FALSE], 2L,
                   scaling$center[keep], "-")
        z <- sweep(z, 2L, scaling$scale[keep], "/")
        xc[[i]] <- sweep(z, 2L, colMeans(z), "-")
        ybar[[i]] <- mean(y[[i]])
        yc[[i]] <- as.numeric(y[[i]]) - ybar[[i]]
        block <- ((i - 1L) * p + 1L):(i * p)
        blocks[[i]] <- block
        gram[block, block] <- weight[[i]] * crossprod(xc[[i]])
        rhs[block] <- weight[[i]] *
            as.numeric(crossprod(xc[[i]], yc[[i]]))
        raw_rank[[i]] <- qr(xc[[i]])$rank
        raw_kappa[[i]] <- .condition_ridge_kappa(cbind(1, z))
    }
    names(ybar) <- names(raw_rank) <- names(raw_kappa) <- conditions

    system <- gram + n_total * lambda * .condition_ridge_penalty(k, p)
    r <- tryCatch(chol(system), error = function(e) NULL)
    if (is.null(r)) return(list(status = "failed"))
    theta <- backsolve(r, forwardsolve(t(r), rhs))
    inv <- chol2inv(r)
    df <- vapply(seq_len(k), function(i) {
        b <- blocks[[i]]
        sum(inv[b, b, drop = FALSE] * gram[b, b, drop = FALSE])
    }, numeric(1))
    names(df) <- conditions
    residual_df <- n_total - k - sum(df)
    if (!is.finite(residual_df) || residual_df < min_residual_df) {
        return(list(status = "insufficient_df"))
    }

    beta_z_small <- matrix(theta, nrow = k, byrow = TRUE)
    beta_small <- sweep(beta_z_small, 2L, scaling$scale[keep], "/")
    beta_z[, keep] <- beta_z_small
    beta[, keep] <- beta_small

    intercept <- prediction <- vector("list", k)
    sse <- numeric(k)
    rsq <- numeric(k)
    for (i in seq_len(k)) {
        intercept[[i]] <- mean(y[[i]]) -
            sum(colMeans(x[[i]][, keep, drop = FALSE]) * beta_small[i, ])
        prediction[[i]] <- as.numeric(
            intercept[[i]] + x[[i]][, keep, drop = FALSE] %*% beta_small[i, ]
        )
        residual <- as.numeric(y[[i]]) - prediction[[i]]
        sse[[i]] <- sum(residual^2)
        tss <- sum((as.numeric(y[[i]]) - mean(y[[i]]))^2)
        rsq[[i]] <- if (tss > 0) 1 - sse[[i]] / tss else NA_real_
    }
    intercept <- unlist(intercept, use.names = FALSE)
    names(intercept) <- names(prediction) <- names(sse) <- names(rsq) <- conditions

    covariance <- NULL
    if (isTRUE(inference)) {
        meat <- matrix(0, k * p, k * p)
        for (i in seq_len(k)) {
            b <- blocks[[i]]
            residual <- as.numeric(y[[i]]) - prediction[[i]]
            xr <- sweep(xc[[i]], 1L, residual, "*")
            meat[b, b] <- weight[[i]]^2 * crossprod(xr)
        }
        covariance <- n_total / max(1, residual_df) *
            (inv %*% meat %*% inv)
        se_z <- matrix(sqrt(pmax(0, diag(covariance))),
                       nrow = k, byrow = TRUE)
        se_small <- sweep(se_z, 2L, scaling$scale[keep], "/")
        se[, keep] <- se_small
        zstat <- beta_small / se_small
        zstat[!is.finite(zstat)] <- NA_real_
        statistic[, keep] <- zstat
        pval[, keep] <- 2 * stats::pnorm(-abs(zstat))
    }

    list(
        status = "ok", beta = beta, beta_z = beta_z, se = se,
        statistic = statistic, pval = pval, intercept = intercept,
        prediction = prediction, rsq = rsq, raw_rank = raw_rank,
        effective_df = df, residual_df = residual_df, raw_kappa = raw_kappa,
        regularized_kappa = .condition_ridge_kappa(system),
        zero_variance = zero_variance, informative = informative,
        informative_index = keep, covariance_z = covariance, weight = weight
    )
}

.condition_ridge_cv <- function(
    x, y, folds, control, min_residual_df = 1L) {
    nfolds <- attr(folds, "nfolds")
    grid <- control$lambda_grid
    loss <- matrix(NA_real_, nfolds, length(grid))
    for (fold in seq_len(nfolds)) {
        train_x <- train_y <- valid_x <- valid_y <- vector("list", length(x))
        names(train_x) <- names(train_y) <- names(valid_x) <-
            names(valid_y) <- names(x)
        for (i in seq_along(x)) {
            train <- folds[[i]] != fold
            valid <- !train
            train_x[[i]] <- x[[i]][train, , drop = FALSE]
            train_y[[i]] <- y[[i]][train]
            valid_x[[i]] <- x[[i]][valid, , drop = FALSE]
            valid_y[[i]] <- y[[i]][valid]
        }
        scaling <- .condition_ridge_scaling(train_x, control$scale_floor)
        for (j in seq_along(grid)) {
            fit <- .condition_ridge_fit(
                train_x, train_y, scaling, grid[[j]],
                min_residual_df, inference = FALSE
            )
            if (!identical(fit$status, "ok")) next
            mse <- vapply(seq_along(valid_x), function(i) {
                pred <- fit$intercept[[i]] +
                    valid_x[[i]] %*% fit$beta[i, ]
                mean((as.numeric(valid_y[[i]]) - as.numeric(pred))^2)
            }, numeric(1))
            if (all(is.finite(mse))) loss[fold, j] <- mean(mse)
        }
    }
    mean_loss <- colMeans(loss, na.rm = TRUE)
    n_valid <- colSums(is.finite(loss))
    mean_loss[n_valid == 0L | !is.finite(mean_loss)] <- Inf
    se_loss <- vapply(seq_along(grid), function(j) {
        value <- loss[, j]
        value <- value[is.finite(value)]
        if (length(value) < 2L) NA_real_ else stats::sd(value) / sqrt(length(value))
    }, numeric(1))
    if (!any(is.finite(mean_loss))) {
        stop("No ridge lambda produced a valid CV fit.", call. = FALSE)
    }
    i_min <- which.min(mean_loss)
    i_pick <- i_min
    if (identical(control$lambda_rule, "1se")) {
        limit <- mean_loss[[i_min]] +
            ifelse(is.finite(se_loss[[i_min]]), se_loss[[i_min]], 0)
        eligible <- which(is.finite(mean_loss) & mean_loss <= limit)
        i_pick <- eligible[[which.max(grid[eligible])]]
    }
    lambda <- grid[[i_pick]]

    oof <- lapply(y, function(value) rep(NA_real_, length(value)))
    names(oof) <- names(y)
    for (fold in seq_len(nfolds)) {
        train_x <- train_y <- vector("list", length(x))
        names(train_x) <- names(train_y) <- names(x)
        for (i in seq_along(x)) {
            train <- folds[[i]] != fold
            train_x[[i]] <- x[[i]][train, , drop = FALSE]
            train_y[[i]] <- y[[i]][train]
        }
        scaling <- .condition_ridge_scaling(train_x, control$scale_floor)
        fit <- .condition_ridge_fit(
            train_x, train_y, scaling, lambda,
            min_residual_df, inference = FALSE
        )
        if (!identical(fit$status, "ok")) next
        for (i in seq_along(x)) {
            index <- which(folds[[i]] == fold)
            oof[[i]][index] <- fit$intercept[[i]] +
                x[[i]][index, , drop = FALSE] %*% fit$beta[i, ]
        }
    }
    rsq_oof <- vapply(seq_along(y), function(i) {
        value <- as.numeric(y[[i]])
        pred <- as.numeric(oof[[i]])
        valid <- is.finite(value) & is.finite(pred)
        if (sum(valid) < 2L) return(NA_real_)
        tss <- sum((value[valid] - mean(value[valid]))^2)
        if (tss <= 0) NA_real_ else
            1 - sum((value[valid] - pred[valid])^2) / tss
    }, numeric(1))
    names(rsq_oof) <- names(y)
    list(
        lambda = lambda, lambda_min = grid[[i_min]],
        cv_mse = mean_loss[[i_pick]], cv_se = se_loss[[i_pick]],
        rsq_oof = rsq_oof,
        curve = data.frame(lambda = grid, mean_mse = mean_loss,
                           se_mse = se_loss, n_folds = n_valid)
    )
}
