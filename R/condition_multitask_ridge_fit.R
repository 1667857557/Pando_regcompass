# One-pass common-dictionary E-star production plus separate inference.

.condition_fit_engine <-
    "condition_union_Estar_z025_inference_separated"
.condition_inference_schema <-
    "frozen_dictionary_condition_local_gaussian_lm_edge_omnibus_v1"

.condition_no_fusion_condition_fit <- function(
    x, y, scaling, rank_tol = 1e-10, min_residual_df = 1L) {
    p_full <- ncol(x)
    estimate <- se <- statistic <- pval <- variance <- rep(NA_real_, p_full)
    estimable <- rep(FALSE, p_full)
    keep <- which(as.logical(scaling$informative))
    if (!length(keep)) {
        return(list(
            status = "no_informative_predictor", estimate = estimate,
            se = se, statistic = statistic, pval = pval,
            variance = variance, estimable = estimable,
            sigma2 = NA_real_, residual_df = length(y) - 1L, rank = 1L
        ))
    }
    z <- sweep(x[, keep, drop = FALSE], 2L, scaling$center[keep], "-")
    z <- sweep(z, 2L, scaling$scale[keep], "/")
    design <- cbind(`(Intercept)` = 1, z)
    decomposition <- .condition_symmetric_pinv(
        crossprod(design), rank_tol = rank_tol
    )
    theta <- as.numeric(
        decomposition$inverse %*% crossprod(design, as.numeric(y))
    )
    residual <- as.numeric(y) - as.numeric(design %*% theta)
    residual_df <- nrow(design) - decomposition$rank
    if (!is.finite(residual_df) || residual_df < min_residual_df) {
        return(list(
            status = "insufficient_df", estimate = estimate,
            se = se, statistic = statistic, pval = pval,
            variance = variance, estimable = estimable,
            sigma2 = NA_real_, residual_df = as.integer(residual_df),
            rank = as.integer(decomposition$rank)
        ))
    }
    y_scale <- stats::var(as.numeric(y))
    variance_floor <- .Machine$double.eps * max(1, y_scale, na.rm = TRUE)
    sigma2 <- max(sum(residual^2) / residual_df, variance_floor)
    covariance <- sigma2 * decomposition$inverse
    projector <- decomposition$projector
    for (local_edge in seq_along(keep)) {
        full_edge <- keep[[local_edge]]
        column <- 1L + local_edge
        unit <- numeric(ncol(design)); unit[[column]] <- 1
        projection_error <- max(abs(unit - projector[column, ]))
        ok <- is.finite(projection_error) &&
            projection_error <= max(1e-8, 20 * rank_tol)
        raw_scale <- scaling$scale[[full_edge]]
        value <- theta[[column]] / raw_scale
        var_value <- covariance[column, column] / raw_scale^2
        if (!ok || !is.finite(var_value) || var_value <= 0) next
        std <- sqrt(var_value)
        stat <- value / std
        estimate[[full_edge]] <- value
        se[[full_edge]] <- std
        statistic[[full_edge]] <- stat
        pval[[full_edge]] <- 2 * stats::pt(
            -abs(stat), df = residual_df
        )
        variance[[full_edge]] <- var_value
        estimable[[full_edge]] <- TRUE
    }
    list(
        status = "ok", estimate = estimate, se = se,
        statistic = statistic, pval = pval, variance = variance,
        estimable = estimable, sigma2 = sigma2,
        residual_df = as.integer(residual_df),
        rank = as.integer(decomposition$rank)
    )
}

.condition_no_fusion_inference <- function(
    x, y, scaling, rank_tol = 1e-10, min_residual_df = 1L) {
    conditions <- names(x)
    p <- ncol(x[[1L]])
    estimate <- se <- statistic <- pval <- variance <- matrix(
        NA_real_, length(conditions), p,
        dimnames = list(conditions, colnames(x[[1L]]))
    )
    estimable <- matrix(FALSE, length(conditions), p, dimnames = dimnames(estimate))
    sigma2 <- residual_df <- rank <- rep(NA_real_, length(conditions))
    status <- rep(NA_character_, length(conditions))
    names(sigma2) <- names(residual_df) <- names(rank) <- names(status) <- conditions
    for (i in seq_along(conditions)) {
        one <- .condition_no_fusion_condition_fit(
            x = x[[i]], y = y[[i]], scaling = scaling,
            rank_tol = rank_tol, min_residual_df = min_residual_df
        )
        estimate[i, ] <- one$estimate
        se[i, ] <- one$se
        statistic[i, ] <- one$statistic
        pval[i, ] <- one$pval
        variance[i, ] <- one$variance
        estimable[i, ] <- one$estimable
        sigma2[[i]] <- one$sigma2
        residual_df[[i]] <- one$residual_df
        rank[[i]] <- one$rank
        status[[i]] <- one$status
    }
    list(
        schema = .condition_inference_schema,
        estimate = estimate, se = se, statistic = statistic,
        pval = pval, variance = variance, estimable = estimable,
        sigma2 = sigma2, residual_df = residual_df,
        rank = rank, status = status
    )
}

.condition_ridge_contrasts <- function(fit, edges) {
    conditions <- rownames(fit$beta)
    if (length(conditions) < 2L) return(data.frame())
    pairs <- utils::combn(seq_along(conditions), 2L, simplify = FALSE)
    out <- list()
    for (pair in pairs) {
        a <- pair[[1L]]
        b <- pair[[2L]]
        for (edge in seq_len(nrow(edges))) {
            estimate_a <- as.numeric(fit$beta[a, edge])
            estimate_b <- as.numeric(fit$beta[b, edge])
            same_component <- identical(
                as.character(fit$fusion_component_id[a, edge]),
                as.character(fit$fusion_component_id[b, edge])
            )
            out[[length(out) + 1L]] <- data.frame(
                tf = edges$tf[[edge]],
                target = edges$target[[edge]],
                region = edges$region[[edge]],
                edge_id = edges$edge_id[[edge]],
                atac_feature_id = edges$atac_feature_id[[edge]],
                condition_a = conditions[[a]],
                condition_b = conditions[[b]],
                contrast = paste0(conditions[[b]], "-", conditions[[a]]),
                estimate_a = estimate_a,
                estimate_b = estimate_b,
                contrast_estimate = estimate_b - estimate_a,
                contrast_estimable = is.finite(estimate_a) && is.finite(estimate_b),
                contrast_identifiable =
                    fit$contrast_identifiable[[edge]] %in% TRUE,
                contrast_status = if (same_component) {
                    "same_Estar_component"
                } else "ok",
                shared_by_boundary =
                    fit$shared_by_boundary[[edge]] %in% TRUE,
                fused_by_penalty =
                    fit$fused_by_penalty[[edge]] %in% TRUE,
                shared_edge = fit$shared_edge[[edge]] %in% TRUE,
                profile_information_delta =
                    as.numeric(fit$profile_information[[edge]]),
                penalty_family = fit$penalty_family,
                penalty_value = fit$penalty_value,
                reference_condition = fit$reference_condition,
                solver_status = fit$solver_status,
                objective = fit$objective,
                kkt_residual = fit$kkt_residual,
                iterations = fit$iterations,
                stringsAsFactors = FALSE
            )
        }
    }
    do.call(rbind, out)
}

.condition_ridge_target <- function(
    prepared, edges, cells_by_condition, folds = NULL, control,
    min_residual_df, rank_action, reference_condition = NULL) {
    x <- .condition_ridge_predictors(prepared, edges, cells_by_condition)
    y <- lapply(cells_by_condition, function(cells) {
        as.numeric(prepared$gene_data[cells, edges$target[[1L]]])
    })
    names(y) <- names(cells_by_condition)
    scaling <- .condition_ridge_scaling(x, control$scale_floor)
    fit <- .condition_ridge_fit(
        x = x, y = y, scaling = scaling,
        min_residual_df = min_residual_df, control = control,
        reference_condition = reference_condition
    )
    if (!identical(fit$status, "ok")) {
        detail <- if (!is.null(fit$solver_status)) {
            paste0(
                " (solver=", fit$solver_status,
                ", kkt=", format(fit$kkt_residual, digits = 4), ")"
            )
        } else ""
        stop(
            "Final condition E-star fit failed for target `",
            edges$target[[1L]], "`: ", fit$status, detail, ".",
            call. = FALSE
        )
    }
    raw_rank_deficient <- fit$raw_rank < sum(fit$informative)
    if (identical(rank_action, "error") && any(raw_rank_deficient)) {
        stop(
            "Raw common-dictionary design is rank deficient for target `",
            edges$target[[1L]],
            "`; rank_action='error' requested failure.",
            call. = FALSE
        )
    }
    inference <- .condition_no_fusion_inference(
        x = x, y = y, scaling = scaling,
        rank_tol = control$rank_tol,
        min_residual_df = min_residual_df
    )

    shared <- fit$shared_z / scaling$scale
    deviation <- sweep(fit$deviation_z, 2L, scaling$scale, "/")
    condition_informative <- sweep(
        !fit$zero_variance, 2L, fit$informative, "&"
    )
    terms <- sprintf("edge_%07d", edges$candidate_index)
    coefficients <- fits <- vector("list", length(cells_by_condition))
    names(coefficients) <- names(fits) <- names(cells_by_condition)
    for (i in seq_along(cells_by_condition)) {
        condition <- names(cells_by_condition)[[i]]
        estimate <- as.numeric(fit$beta[i, ])
        inference_estimate <- as.numeric(inference$estimate[i, ])
        inference_se <- as.numeric(inference$se[i, ])
        inference_statistic <- as.numeric(inference$statistic[i, ])
        condition_pval <- as.numeric(inference$pval[i, ])
        inference_estimable <- as.logical(inference$estimable[i, ])
        coefficients[[i]] <- data.frame(
            tf = edges$tf,
            target = edges$target,
            region = edges$region,
            term = terms,
            edge_id = edges$edge_id,
            atac_feature_id = edges$atac_feature_id,
            condition = condition,
            z = fit$penalty_value,
            estimate = estimate,
            penalty_effect = estimate,
            estimate_standardized = as.numeric(fit$beta_z[i, ]),
            beta_shared = as.numeric(shared),
            shared_estimate = as.numeric(shared),
            condition_deviation = as.numeric(deviation[i, ]),
            delta_beta = as.numeric(deviation[i, ]),
            contrast_identifiable = as.logical(fit$contrast_identifiable),
            shared_by_boundary = as.logical(fit$shared_by_boundary),
            boundary_condition = as.logical(fit$boundary_condition[i, ]),
            fused_by_penalty = as.logical(fit$fused_by_penalty),
            fusion_component_id = as.character(fit$fusion_component_id[i, ]),
            shared_edge = as.logical(fit$shared_edge),
            raw_information_condition = as.numeric(fit$raw_information[i, ]),
            profile_information = as.numeric(fit$profile_information),
            profile_information_delta = as.numeric(fit$profile_information),
            condition_number = fit$raw_kappa[[i]],
            sigma2_common = fit$sigma2,
            inference_schema = .condition_inference_schema,
            inference_estimate = inference_estimate,
            inference_se = inference_se,
            inference_variance = as.numeric(inference$variance[i, ]),
            inference_statistic = inference_statistic,
            inference_estimable = inference_estimable,
            condition_inference_estimable = inference_estimable,
            condition_pval = condition_pval,
            condition_inference_residual_df = inference$residual_df[[i]],
            condition_inference_sigma2 = inference$sigma2[[i]],
            condition_inference_status = inference$status[[i]],
            std_err = inference_se,
            statistic = inference_statistic,
            pval = condition_pval,
            padj = NA_real_,
            estimable = is.finite(estimate),
            zero_variance = as.logical(fit$zero_variance[i, ]),
            condition_informative = as.logical(condition_informative[i, ]),
            aliased = !inference_estimable,
            reference_condition = fit$reference_condition,
            penalty_family = fit$penalty_family,
            penalty_value = fit$penalty_value,
            solver_status = fit$solver_status,
            objective = fit$objective,
            kkt_residual = fit$kkt_residual,
            iterations = fit$iterations,
            candidate_index = edges$candidate_index,
            source_global = edges$source_global,
            source_conditions = edges$source_conditions,
            n_sources = edges$n_sources,
            stringsAsFactors = FALSE
        )
        local_residual <- length(cells_by_condition[[i]]) -
            1 - fit$raw_rank[[i]]
        fits[[i]] <- data.frame(
            target = edges$target[[1L]],
            condition = condition,
            z = fit$penalty_value,
            rsq = fit$rsq[[i]],
            rsq_definition = "scheme_e_z025_full_data_R2_diagnostic",
            rank = as.integer(1L + fit$raw_rank[[i]]),
            raw_rank = as.integer(fit$raw_rank[[i]]),
            residual_df = as.integer(max(0, floor(local_residual))),
            residual_df_effective_joint = fit$residual_df,
            condition_number = fit$raw_kappa[[i]],
            fit_status = "ok",
            intercept = fit$intercept[[i]],
            sigma2_common = fit$sigma2,
            inference_sigma2 = inference$sigma2[[i]],
            inference_residual_df = inference$residual_df[[i]],
            inference_rank = inference$rank[[i]],
            inference_status = inference$status[[i]],
            reference_condition = fit$reference_condition,
            deviation_z = fit$penalty_value,
            penalty_family = fit$penalty_family,
            solver_status = fit$solver_status,
            objective = fit$objective,
            kkt_residual = fit$kkt_residual,
            iterations = fit$iterations,
            nvariables = nrow(edges),
            nvariables_dictionary = nrow(edges),
            nvariables_contrast_identifiable =
                sum(fit$contrast_identifiable %in% TRUE),
            n_shared_by_boundary = sum(fit$shared_by_boundary %in% TRUE),
            n_fused_by_penalty = sum(fit$fused_by_penalty %in% TRUE),
            n_shared_edges = sum(fit$shared_edge %in% TRUE),
            n_zero_variance = sum(fit$zero_variance[i, ]),
            condition_weight = 1,
            predictor_scale_reference = scaling$reference,
            inference_schema = .condition_inference_schema,
            orthogonality_error = fit$orthogonality_error,
            dr_error = fit$dr_error,
            stringsAsFactors = FALSE
        )
    }
    list(
        coefficients = do.call(rbind, coefficients),
        contrasts = .condition_ridge_contrasts(fit, edges),
        fit = do.call(rbind, fits),
        scaling = scaling,
        contrast_tree = fit$contrast_tree,
        solver = list(
            status = fit$solver_status,
            objective = fit$objective,
            kkt_residual = fit$kkt_residual,
            iterations = fit$iterations,
            penalty_family = fit$penalty_family,
            penalty_value = fit$penalty_value,
            orthogonality_error = fit$orthogonality_error,
            dr_error = fit$dr_error
        )
    )
}

.condition_ridge_fit_contract_one_pass <- function(
    object, fit, prepared, control, rank_action = "mark",
    min_residual_df = 1L, parallel = FALSE, verbose = TRUE,
    progress_phase = NULL, progress_label = NULL,
    checkpoint_dir = NULL, resume = TRUE) {
    .condition_validate_dictionary(fit$edge_dictionary, prepared)
    cells <- fit$condition_cell_ids[fit$condition_levels]
    if (length(cells) < 2L) {
        stop("Conditional E-star requires at least two conditions.",
             call. = FALSE)
    }
    if (any(lengths(cells) < 3L)) {
        stop("Every condition E-star fit needs at least three cells.",
             call. = FALSE)
    }
    targets <- unique(as.character(fit$edge_dictionary$target))
    names(targets) <- targets
    if (is.null(progress_phase)) progress_phase <- "condition_Estar_z025"
    if (is.null(progress_label)) progress_label <- fit$cell_type %||% ""
    checkpoint_fingerprint <- if (is.null(checkpoint_dir)) NULL else {
        .condition_hash_object(list(
            schema = "pando_Estar_target_checkpoint_input_v3_qscale_block_pair",
            preprocessing_fingerprint = prepared$preprocessing_fingerprint,
            edge_id = as.character(fit$edge_dictionary$edge_id),
            condition_cells = cells, control = control,
            rank_action = rank_action,
            min_residual_df = as.numeric(min_residual_df),
            reference_condition = fit$reference_condition
        ))
    }

    result <- .pando_target_payload_map(
        keys = targets,
        build_payload = function(target) {
            .pando_ridge_target_payload(
                prepared = prepared,
                edge_dictionary = fit$edge_dictionary,
                target = target,
                cells = cells,
                control = control,
                min_residual_df = min_residual_df,
                rank_action = rank_action,
                reference_condition = fit$reference_condition
            )
        },
        worker_name = ".pando_ridge_target_worker",
        parallel = parallel,
        verbose = verbose,
        phase = progress_phase,
        label = progress_label,
        checkpoint_dir = checkpoint_dir,
        resume = resume,
        checkpoint_fingerprint = checkpoint_fingerprint
    )

    coefficient <- do.call(rbind, lapply(result, `[[`, "coefficients"))
    contrast <- do.call(rbind, lapply(result, `[[`, "contrasts"))
    fit_table <- do.call(rbind, lapply(result, `[[`, "fit"))
    rownames(coefficient) <- rownames(fit_table) <- NULL
    if (is.data.frame(contrast) && nrow(contrast)) rownames(contrast) <- NULL

    coefficient$penalty_effect <- coefficient$estimate
    coefficient$direction <- ifelse(
        coefficient$estimate > 0, "positive",
        ifelse(coefficient$estimate < 0, "negative", "zero")
    )
    coefficient$effect_definition <-
        "E_star_z025_continuous_condition_coefficient_raw_tf_atac_units"
    coefficient$inference_scope <-
        "frozen_dictionary_no_fusion_condition_local_gaussian_lm"

    fit$model_schema <- .condition_multitask_ridge_schema
    fit$fit_engine <- .condition_fit_engine
    fit$coefficient_scale <- "raw_tf_atac_interaction_units"
    fit$internal_predictor_scale <- "equal_condition_within_condition_rms"
    fit$inference_schema <- .condition_inference_schema
    fit$inference_scope <-
        "separate_no_fusion_condition_local_lm_then_exact_edge_omnibus"
    fit$coefficients <- coefficient
    fit$contrasts <- contrast
    fit$edge_inference <- NULL
    fit$fit <- fit_table
    fit$scale <- FALSE
    fit$projection_effect_column <- "penalty_effect"
    fit$projection_policy <- "exact_edge_whole_network_BH_common_topology"
    fit$condition_e_control <- control
    fit$deviation_penalty <- list(
        family = .condition_E_star_penalty_family,
        z = .condition_E_star_z
    )
    fit$target_solver <- lapply(result, `[[`, "solver")
    fit$target_scaling <- lapply(result, `[[`, "scaling")
    fit$target_contrast_tree <- lapply(result, `[[`, "contrast_tree")
    fit$rsq_definition <- "scheme_e_z025_full_data_R2_diagnostic"
    class(fit) <- c("ConditionGRNFit", "list")

    for (condition in fit$condition_levels) {
        network_name <- fit$network_names[[condition]]
        coefs_one <- coefficient[
            coefficient$condition == condition, , drop = FALSE
        ]
        fit_one <- fit_table[
            fit_table$condition == condition, , drop = FALSE
        ]
        network <- methods::new(
            Class = "Network",
            features = unique(as.character(coefs_one$target)),
            coefs = coefs_one,
            fit = fit_one,
            params = list(
                method = "condition_Estar_z025",
                family = "gaussian_identity",
                fit_mode = "frozen_common_dictionary_Estar_z025",
                condition = condition,
                reference_condition = fit$reference_condition,
                edge_dictionary = fit$edge_dictionary,
                scale = FALSE,
                internal_scale_reference =
                    "equal_condition_within_condition_rms",
                exported_coefficient_scale =
                    "raw_tf_atac_interaction_units",
                interaction = ":",
                rna_layer = prepared$rna_layer,
                peak_layer = prepared$peak_layer,
                peak_value_type = prepared$peak_value_type,
                preprocessing_fingerprint = prepared$preprocessing_fingerprint,
                adjust_method = fit$adjust_method,
                padj_threshold = fit$padj_threshold,
                bh_scope = "exact_edge_whole_cell_type_network_BH",
                projection_policy = fit$projection_policy,
                deviation_penalty = fit$deviation_penalty,
                inference_schema = fit$inference_schema,
                condition_e_control = control
            )
        )
        object@grn@networks[[network_name]] <- network
        object@grn@active_network <- network_name
    }

    result <- NULL
    invisible(gc(verbose = FALSE, full = TRUE))
    list(object = object, fit = fit)
}
