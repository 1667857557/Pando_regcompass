# Canonical single-task ridge helpers for standard Pando.
#
# Conditional GRNs now use Scheme E. Standard method='ridge' remains an
# independent single-task ridge path so the conditional redesign cannot change
# ordinary Pando GRNs.

.pando_standard_ridge_family_ok <- function(family) {
    if (is.character(family) && length(family) == 1L) {
        return(tolower(family) == "gaussian")
    }
    identical(tryCatch(family$family, error = function(e) NULL), "gaussian") &&
        identical(tryCatch(family$link, error = function(e) NULL), "identity")
}

.pando_standard_ridge_control <- function(control = list()) {
    if (is.null(control)) control <- list()
    if (!is.list(control)) stop("`ridge_control` must be a list.", call. = FALSE)
    defaults <- list(
        lambda_grid = 10^seq(-3, 2, length.out = 9L),
        lambda_rule = "1se", cv_folds = 5L, seed = 1L, scale_floor = 1e-8
    )
    unknown <- setdiff(names(control), names(defaults))
    if (length(unknown)) {
        stop("Unknown standard `ridge_control` field(s): ",
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
        stop("`ridge_control$cv_folds` must be an integer >= 2.", call. = FALSE)
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

.pando_standard_ridge_folds <- function(n, nfolds, seed) {
    nfolds <- min(as.integer(nfolds), as.integer(n))
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
    order <- sample.int(n)
    fold <- rep(seq_len(nfolds), length.out = n)
    answer <- integer(n)
    answer[order] <- fold
    attr(answer, "nfolds") <- nfolds
    answer
}

.pando_standard_ridge_numeric <- function(
    x, y, lambda, scale_floor = 1e-8, min_residual_df = 1L,
    inference = TRUE) {
    x <- as.matrix(x)
    y <- as.numeric(y)
    if (!nrow(x) || nrow(x) != length(y) || any(!is.finite(x)) ||
        any(!is.finite(y))) return(list(status = "nonfinite_or_misaligned_input"))
    center <- colMeans(x)
    xc_raw <- sweep(x, 2L, center, "-")
    scale <- sqrt(colMeans(xc_raw^2))
    informative <- is.finite(scale) & scale > scale_floor
    scale[!informative] <- 1
    beta <- beta_z <- rep(0, ncol(x))
    se <- statistic <- pval <- rep(NA_real_, ncol(x))
    zero_variance <- !informative
    intercept_only <- function() {
        prediction <- rep(mean(y), length(y))
        tss <- sum((y - mean(y))^2)
        list(
            status = "ok", beta = beta, beta_z = beta_z, se = se,
            statistic = statistic, pval = pval, intercept = mean(y),
            prediction = prediction, rsq = if (tss > 0) 0 else NA_real_,
            raw_rank = 0L, effective_df = 0, residual_df = length(y) - 1,
            raw_kappa = NA_real_, regularized_kappa = NA_real_,
            zero_variance = zero_variance, informative = informative,
            scale = scale, center = center
        )
    }
    if (!any(informative)) return(intercept_only())
    keep <- which(informative)
    z <- sweep(x[, keep, drop = FALSE], 2L, center[keep], "-")
    z <- sweep(z, 2L, scale[keep], "/")
    z <- sweep(z, 2L, colMeans(z), "-")
    yc <- y - mean(y)
    gram <- crossprod(z)
    system <- gram + length(y) * lambda * diag(length(keep))
    chol_system <- tryCatch(chol(system), error = function(e) NULL)
    if (is.null(chol_system)) return(list(status = "failed"))
    inverse <- chol2inv(chol_system)
    theta <- backsolve(
        chol_system,
        forwardsolve(t(chol_system), as.numeric(crossprod(z, yc)))
    )
    effective_df <- sum(inverse * gram)
    residual_df <- length(y) - 1 - effective_df
    if (!is.finite(residual_df) || residual_df < min_residual_df) {
        return(list(status = "insufficient_df"))
    }
    beta_z[keep] <- theta
    beta[keep] <- theta / scale[keep]
    intercept <- mean(y) - sum(colMeans(x[, keep, drop = FALSE]) * beta[keep])
    prediction <- as.numeric(intercept + x[, keep, drop = FALSE] %*% beta[keep])
    residual <- y - prediction
    tss <- sum((y - mean(y))^2)
    rsq <- if (tss > 0) 1 - sum(residual^2) / tss else NA_real_
    if (isTRUE(inference)) {
        xr <- sweep(z, 1L, residual, "*")
        covariance <- length(y) / max(1, residual_df) *
            (inverse %*% crossprod(xr) %*% inverse)
        se_z <- sqrt(pmax(0, diag(covariance)))
        se[keep] <- se_z / scale[keep]
        statistic[keep] <- beta[keep] / se[keep]
        statistic[!is.finite(statistic)] <- NA_real_
        pval[keep] <- 2 * stats::pnorm(-abs(statistic[keep]))
    }
    list(
        status = "ok", beta = beta, beta_z = beta_z, se = se,
        statistic = statistic, pval = pval, intercept = intercept,
        prediction = prediction, rsq = rsq,
        raw_rank = qr(z)$rank, effective_df = effective_df,
        residual_df = residual_df,
        raw_kappa = .condition_ridge_kappa(cbind(1, z)),
        regularized_kappa = .condition_ridge_kappa(system),
        zero_variance = zero_variance, informative = informative,
        scale = scale, center = center
    )
}

.pando_standard_ridge_cv <- function(
    x, y, folds, control, min_residual_df = 1L) {
    nfolds <- attr(folds, "nfolds")
    loss <- matrix(NA_real_, nfolds, length(control$lambda_grid))
    for (fold in seq_len(nfolds)) {
        train <- folds != fold
        valid <- !train
        for (j in seq_along(control$lambda_grid)) {
            fit <- .pando_standard_ridge_numeric(
                x[train, , drop = FALSE], y[train], control$lambda_grid[[j]],
                scale_floor = control$scale_floor,
                min_residual_df = min_residual_df, inference = FALSE
            )
            if (!identical(fit$status, "ok")) next
            pred <- fit$intercept + x[valid, , drop = FALSE] %*% fit$beta
            residual <- y[valid] - as.numeric(pred)
            if (length(residual) && all(is.finite(residual))) {
                loss[fold, j] <- mean(residual^2)
            }
        }
    }
    mean_loss <- colMeans(loss, na.rm = TRUE)
    n_valid <- colSums(is.finite(loss))
    mean_loss[n_valid == 0L | !is.finite(mean_loss)] <- Inf
    se_loss <- vapply(seq_along(control$lambda_grid), function(j) {
        value <- loss[, j]
        value <- value[is.finite(value)]
        if (length(value) < 2L) NA_real_ else stats::sd(value) / sqrt(length(value))
    }, numeric(1))
    if (!any(is.finite(mean_loss))) {
        stop("No standard ridge lambda produced a valid CV fit.", call. = FALSE)
    }
    i_min <- which.min(mean_loss)
    i_pick <- i_min
    if (identical(control$lambda_rule, "1se")) {
        limit <- mean_loss[[i_min]] +
            ifelse(is.finite(se_loss[[i_min]]), se_loss[[i_min]], 0)
        eligible <- which(is.finite(mean_loss) & mean_loss <= limit)
        i_pick <- eligible[[which.max(control$lambda_grid[eligible])]]
    }
    list(
        lambda = control$lambda_grid[[i_pick]],
        lambda_min = control$lambda_grid[[i_min]],
        cv_mse = mean_loss[[i_pick]], cv_se = se_loss[[i_pick]],
        curve = data.frame(
            lambda = control$lambda_grid, mean_mse = mean_loss,
            se_mse = se_loss, n_folds = n_valid
        )
    )
}

.pando_standard_ridge_target <- function(
    prepared, edges, cells, folds, control, min_residual_df, rank_action) {
    x <- .condition_ridge_predictors(
        prepared, edges, list(standard = cells)
    )[[1L]]
    y <- as.numeric(prepared$gene_data[cells, edges$target[[1L]]])
    cv <- .pando_standard_ridge_cv(x, y, folds, control, min_residual_df)
    fit <- .pando_standard_ridge_numeric(
        x, y, cv$lambda, scale_floor = control$scale_floor,
        min_residual_df = min_residual_df, inference = TRUE
    )
    if (!identical(fit$status, "ok")) {
        stop("Final standard ridge failed for target `", edges$target[[1L]],
             "`: ", fit$status, ".", call. = FALSE)
    }
    rank_deficient <- fit$raw_rank < sum(fit$informative)
    if (identical(rank_action, "error") && rank_deficient) {
        stop("Raw standard ridge design is rank deficient for target `",
             edges$target[[1L]], "`.", call. = FALSE)
    }
    coefficient <- data.frame(
        tf = edges$tf, target = edges$target, region = edges$region,
        term = sprintf("edge_%07d", edges$candidate_index),
        edge_id = edges$edge_id, atac_feature_id = edges$atac_feature_id,
        estimate = fit$beta, estimate_standardized = fit$beta_z,
        std_err = fit$se, statistic = fit$statistic, pval = fit$pval,
        estimable = fit$informative & is.finite(fit$beta),
        zero_variance = fit$zero_variance, aliased = FALSE,
        candidate_index = edges$candidate_index,
        source_global = edges$source_global,
        source_conditions = edges$source_conditions,
        n_sources = edges$n_sources, stringsAsFactors = FALSE
    )
    fit_table <- data.frame(
        target = edges$target[[1L]], rsq = fit$rsq,
        rank = as.integer(1L + fit$raw_rank), raw_rank = as.integer(fit$raw_rank),
        residual_df = as.integer(max(0, floor(length(cells) - 1 - fit$effective_df))),
        effective_df = fit$effective_df,
        residual_df_effective_joint = fit$residual_df,
        condition_number = fit$raw_kappa,
        condition_number_regularized = fit$regularized_kappa,
        fit_status = "ok", intercept = fit$intercept,
        lambda = cv$lambda, lambda_min = cv$lambda_min,
        lambda_rule = control$lambda_rule, cv_mse = cv$cv_mse, cv_se = cv$cv_se,
        design_rank_deficient = rank_deficient,
        nvariables = nrow(edges), nvariables_dictionary = nrow(edges),
        nvariables_estimable = sum(coefficient$estimable),
        n_zero_variance = sum(fit$zero_variance), n_aliased = 0L,
        predictor_scale_reference = "within_task_rms",
        stringsAsFactors = FALSE
    )
    list(coefficients = coefficient, fit = fit_table, cv = cv)
}

.pando_standard_ridge_target_worker <- function(payload) {
    .pando_standard_ridge_target(
        prepared = payload$prepared, edges = payload$edges,
        cells = payload$cells[[1L]], folds = payload$folds[[1L]],
        control = payload$control,
        min_residual_df = payload$min_residual_df,
        rank_action = payload$rank_action
    )
}

.pando_standard_ridge_fit <- function(
    object, genes, network_name,
    peak_to_gene_method, upstream, downstream, extend, only_tss,
    tf_cor, peak_cor, adjust_method, padj_threshold,
    rank_action, min_residual_df, ridge_control,
    parallel, verbose) {
    control <- .pando_standard_ridge_control(ridge_control)
    rank_action <- match.arg(rank_action, c("mark", "error"))
    if (!is.numeric(min_residual_df) || length(min_residual_df) != 1L ||
        !is.finite(min_residual_df) || min_residual_df < 1L) {
        stop("`min_residual_df` must be one positive number.", call. = FALSE)
    }
    if (!is.numeric(padj_threshold) || length(padj_threshold) != 1L ||
        !is.finite(padj_threshold) || padj_threshold <= 0 ||
        padj_threshold >= 1) {
        stop("`padj_threshold` must be one finite number in (0, 1).",
             call. = FALSE)
    }
    prepared <- .condition_prepare_common_input(
        object = object, genes = genes,
        peak_to_gene_method = peak_to_gene_method,
        upstream = upstream, downstream = downstream, extend = extend,
        only_tss = only_tss, rna_layer = "data", peak_layer = "data",
        peak_value_type = "normalized", verbose = verbose
    )
    cells <- rownames(prepared$gene_data)
    if (length(cells) < 3L) {
        stop("Standard ridge requires at least three paired cells.", call. = FALSE)
    }
    candidates <- .condition_discover_edges_compact(
        prepared = prepared, cells = cells, source_label = "standard",
        source_type = "global", tf_cor = tf_cor, peak_cor = peak_cor,
        parallel = parallel, verbose = verbose
    )
    if (!is.data.frame(candidates) || !nrow(candidates)) {
        stop("Standard ridge candidate discovery produced no TF-peak-target edges.",
             call. = FALSE)
    }
    dictionary <- union_grn_edges(
        global_edges = NULL, condition_edges = list(standard = candidates)
    )
    folds <- .pando_standard_ridge_folds(
        length(cells), control$cv_folds, control$seed
    )
    target_names <- unique(as.character(dictionary$target))
    names(target_names) <- target_names
    result <- .pando_target_payload_map(
        keys = target_names,
        build_payload = function(target) {
            .pando_ridge_target_payload(
                prepared = prepared, edge_dictionary = dictionary,
                target = target, cells = list(standard = cells),
                folds = list(standard = folds), control = control,
                min_residual_df = min_residual_df, rank_action = rank_action
            )
        },
        worker_name = ".pando_standard_ridge_target_worker",
        parallel = parallel, verbose = verbose,
        phase = "ridge_standard", label = network_name
    )
    coefficient <- do.call(rbind, lapply(result, `[[`, "coefficients"))
    fit_table <- do.call(rbind, lapply(result, `[[`, "fit"))
    rownames(coefficient) <- rownames(fit_table) <- NULL
    coefficient$padj <- NA_real_
    valid <- which(coefficient$estimable %in% TRUE & is.finite(coefficient$pval))
    if (length(valid)) {
        coefficient$padj[valid] <- stats::p.adjust(
            coefficient$pval[valid], method = adjust_method
        )
    }
    coefficient$statistically_supported <-
        coefficient$estimable %in% TRUE & is.finite(coefficient$padj) &
        coefficient$padj < padj_threshold
    coefficient$significant <- coefficient$statistically_supported
    coefficient$active <- coefficient$statistically_supported
    coefficient$penalty_effect <- ifelse(
        coefficient$active, coefficient$estimate, 0
    )
    network <- methods::new(
        Class = "Network",
        features = unique(as.character(coefficient$target)),
        coefs = coefficient, fit = fit_table,
        params = list(
            method = "ridge", family = "gaussian_identity",
            fit_mode = "single_task_ridge",
            edge_dictionary = dictionary, scale = FALSE,
            internal_scale_reference = "within_task_rms",
            exported_coefficient_scale = "raw_tf_atac_interaction_units",
            interaction = ":", adjust_method = adjust_method,
            padj_threshold = padj_threshold,
            projection_policy = "standard_ridge_bh_diagnostic",
            ridge_solver = "single_task_cv_ridge",
            ridge_control = control
        )
    )
    object@grn@networks[[network_name]] <- network
    object@grn@active_network <- network_name
    object@grn@params$analysis_mode <- "standard_grn"
    object@grn@params$standard_fit_method <- "ridge"
    object@grn@params$standard_ridge_schema <-
        "pando_standard_grn_single_task_ridge_v1"
    object@grn@params$standard_ridge_contract <- list(
        schema = "pando_standard_grn_single_task_ridge_v1",
        solver = "single_task_cv_ridge",
        candidate_edges = nrow(dictionary),
        fitted_targets = length(unique(as.character(fit_table$target))),
        coefficient_scale = "raw_tf_atac_interaction_units",
        predictor_scale_reference = "within_task_rms",
        adjust_method = as.character(adjust_method),
        padj_threshold = as.numeric(padj_threshold), rank_action = rank_action,
        min_residual_df = min_residual_df, ridge_control = control
    )
    object
}
