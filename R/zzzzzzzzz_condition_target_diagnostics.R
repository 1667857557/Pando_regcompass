# Stage-aware target diagnostics and single-edge ridge shape fix.
#
# This file loads after zzzzzzzz_target_payload_parallel.R. It does not alter
# candidate filters, ridge penalties, CV, BH screening, or projection equations.

.pando_compact_discover_edges_progress_impl <- .condition_discover_edges_prepared
.pando_compact_ridge_one_pass_progress_impl <- .condition_ridge_refit_contract_one_pass

.pando_diag_or <- function(value, fallback) {
    if (is.null(value) || !length(value)) fallback else value
}

.pando_target_progress_context <- function(default = "target_work") {
    context <- getOption("Pando.target_progress_context", NULL)
    if (!is.list(context)) context <- list()
    stage <- as.character(.pando_diag_or(context$stage, default))
    if (length(stage) != 1L || is.na(stage) || !nzchar(stage)) stage <- default
    label <- as.character(.pando_diag_or(context$label, ""))
    if (length(label) != 1L || is.na(label)) label <- ""
    list(stage = stage, label = label)
}

.pando_target_payload_map <- function(
    keys, build_payload, worker, parallel = FALSE, verbose = TRUE) {
    if (!length(keys)) return(list())
    if (!is.function(build_payload) || !is.function(worker)) {
        stop("Target payload mapping requires builder and worker functions.",
             call. = FALSE)
    }
    key_names <- names(keys)
    if (is.null(key_names) || any(!nzchar(key_names))) {
        key_names <- as.character(keys)
    }
    context <- .pando_target_progress_context()
    phase <- context$stage
    phase_label <- if (nzchar(context$label)) {
        paste0(phase, ":", context$label)
    } else phase
    out <- vector("list", length(keys))
    names(out) <- key_names
    batch_size <- min(length(keys), .pando_target_worker_limit(parallel))
    starts <- seq.int(1L, length(keys), by = batch_size)
    n_batches <- length(starts)
    if (isTRUE(verbose)) {
        message(
            "Pando target phase=", phase_label,
            " | started | targets=", length(keys),
            ";batch_size=", batch_size,
            ";batches=", n_batches,
            ";parallel=", isTRUE(parallel)
        )
    }
    for (batch_index in seq_along(starts)) {
        start <- starts[[batch_index]]
        index <- seq.int(start, min(length(keys), start + batch_size - 1L))
        payloads <- lapply(index, function(i) build_payload(keys[[i]]))
        tasks <- lapply(seq_along(payloads), function(j) {
            list(key = key_names[index[[j]]], payload = payloads[[j]])
        })
        names(tasks) <- key_names[index]
        run_one <- function(task) {
            tryCatch(
                worker(task$payload),
                error = function(error) {
                    stop(
                        "Pando target task failed [phase=", phase_label,
                        "; target=", task$key, "]: ",
                        conditionMessage(error),
                        call. = FALSE
                    )
                }
            )
        }
        chunk <- if (isTRUE(parallel)) {
            map_par(tasks, run_one, parallel = TRUE, verbose = FALSE)
        } else {
            lapply(tasks, run_one)
        }
        out[index] <- chunk
        if (isTRUE(verbose)) {
            message(
                "Pando target phase=", phase_label,
                " | completed=", max(index), "/", length(keys),
                " (", sprintf("%.1f", 100 * max(index) / length(keys)), "%)",
                ";batch=", batch_index, "/", n_batches
            )
        }
        payloads <- NULL
        tasks <- NULL
        chunk <- NULL
        invisible(gc(verbose = FALSE, full = TRUE))
    }
    out
}

.pando_discovery_target_worker <- function(payload) {
    on.exit(invisible(gc(verbose = FALSE, full = TRUE)), add = TRUE)
    if (isTRUE(payload$skip)) return(NULL)
    .pando_full_discover_edges_prepared(
        prepared = payload$prepared,
        cells = payload$cells,
        source_label = payload$source_label,
        source_type = payload$source_type,
        tf_cor = payload$tf_cor,
        peak_cor = payload$peak_cor,
        parallel = FALSE,
        verbose = FALSE
    )
}

.condition_discover_edges_prepared <- function(
    prepared, cells, source_label, source_type, tf_cor, peak_cor,
    parallel = FALSE, verbose = TRUE) {
    old_context <- getOption("Pando.target_progress_context", NULL)
    stage <- if (identical(as.character(source_type), "global")) {
        "candidate_global"
    } else {
        "candidate_condition"
    }
    options(Pando.target_progress_context = list(
        stage = stage,
        label = as.character(source_label)
    ))
    on.exit(
        options(Pando.target_progress_context = old_context),
        add = TRUE
    )
    .pando_compact_discover_edges_progress_impl(
        prepared = prepared,
        cells = cells,
        source_label = source_label,
        source_type = source_type,
        tf_cor = tf_cor,
        peak_cor = peak_cor,
        parallel = parallel,
        verbose = verbose
    )
}

.condition_ridge_fit <- function(
    x, y, scaling, lambda, fusion_ratio, min_residual_df = 1L,
    inference = TRUE) {
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

    penalty <- .condition_ridge_penalty(k, p, fusion_ratio)
    system <- gram + n_total * lambda * penalty
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

.condition_ridge_refit_contract_one_pass <- function(
    object, fit, prepared, control, rank_action = "mark",
    min_residual_df = 1L, parallel = FALSE, verbose = TRUE) {
    old_context <- getOption("Pando.target_progress_context", NULL)
    if (!is.list(old_context)) {
        stage <- if (length(fit$condition_levels) == 1L) {
            "ridge_standard"
        } else {
            "ridge_single_pass"
        }
        options(Pando.target_progress_context = list(
            stage = stage,
            label = as.character(.pando_diag_or(fit$cell_type, ""))
        ))
        on.exit(
            options(Pando.target_progress_context = old_context),
            add = TRUE
        )
    }
    .pando_compact_ridge_one_pass_progress_impl(
        object = object,
        fit = fit,
        prepared = prepared,
        control = control,
        rank_action = rank_action,
        min_residual_df = min_residual_df,
        parallel = parallel,
        verbose = verbose
    )
}

.condition_ridge_refit_contract <- function(
    object, fit, prepared, control, rank_action = "mark",
    min_residual_df = 1L, parallel = FALSE, verbose = TRUE) {
    fit$adjust_method <- .condition_validate_adjust_method(fit$adjust_method)
    fit$padj_threshold <- .condition_validate_padj_threshold(fit$padj_threshold)
    candidate_dictionary <- fit$edge_dictionary
    old_context <- getOption("Pando.target_progress_context", NULL)
    on.exit(
        options(Pando.target_progress_context = old_context),
        add = TRUE
    )

    options(Pando.target_progress_context = list(
        stage = "ridge_preliminary",
        label = as.character(.pando_diag_or(fit$cell_type, ""))
    ))
    preliminary <- .condition_ridge_refit_contract_one_pass(
        object = object, fit = fit, prepared = prepared, control = control,
        rank_action = rank_action, min_residual_df = min_residual_df,
        parallel = parallel, verbose = verbose
    )
    screen <- .condition_dictionary_screen(preliminary$fit)
    if (isTRUE(verbose)) {
        message(
            "Pando condition phase=bh_dictionary_screen",
            " | cell_type=", as.character(.pando_diag_or(fit$cell_type, "")),
            ";candidate_edges=", nrow(candidate_dictionary),
            ";supported_edges=", length(screen$keep_edge_ids),
            ";threshold=", format(screen$threshold, trim = TRUE)
        )
    }
    if (!length(screen$keep_edge_ids)) {
        stop(
            "No condition-GRN edge passes BH padj < ",
            format(screen$threshold, trim = TRUE),
            " in any condition; no statistically supported common fit ",
            "dictionary can be constructed.", call. = FALSE
        )
    }

    final_dictionary <- .condition_subset_dictionary(
        candidate_dictionary, screen$keep_edge_ids, screen$summary
    )
    if (nrow(final_dictionary) == nrow(candidate_dictionary) &&
        identical(sort(as.character(final_dictionary$edge_id)),
                  sort(as.character(candidate_dictionary$edge_id)))) {
        final <- preliminary
        final$fit$edge_dictionary <- final_dictionary
        if (isTRUE(verbose)) {
            message(
                "Pando condition phase=ridge_final",
                " | cell_type=", as.character(.pando_diag_or(fit$cell_type, "")),
                ";action=reuse_preliminary",
                ";fit_edges=", nrow(final_dictionary)
            )
        }
    } else {
        final_skeleton <- fit
        final_skeleton$edge_dictionary <- final_dictionary
        final_skeleton$target_genes <-
            unique(as.character(final_dictionary$target))
        final_skeleton$coefficients <- NULL
        final_skeleton$contrasts <- NULL
        final_skeleton$fit <- NULL
        options(Pando.target_progress_context = list(
            stage = "ridge_final",
            label = as.character(.pando_diag_or(fit$cell_type, ""))
        ))
        final <- .condition_ridge_refit_contract_one_pass(
            object = object, fit = final_skeleton, prepared = prepared,
            control = control, rank_action = rank_action,
            min_residual_df = min_residual_df,
            parallel = parallel, verbose = verbose
        )
    }

    final$fit <- .condition_apply_significance_gate(final$fit)
    final$fit$candidate_edge_count <- nrow(candidate_dictionary)
    final$fit$fit_dictionary_edge_count <- nrow(final_dictionary)
    final$fit$fit_dictionary_policy <- .condition_fit_dictionary_policy
    final$fit$dictionary_screening_threshold <- screen$threshold
    final$fit$dictionary_screening_summary <- screen$summary
    final$fit$edge_dictionary <- final_dictionary
    final$object <- .condition_update_network_significance(final$object, final$fit)
    if (isTRUE(verbose)) {
        message(
            "Pando condition phase=condition_fit_complete",
            " | cell_type=", as.character(.pando_diag_or(fit$cell_type, "")),
            ";candidate_edges=", nrow(candidate_dictionary),
            ";fit_edges=", nrow(final_dictionary),
            ";targets=", length(unique(as.character(final_dictionary$target)))
        )
    }
    final
}
