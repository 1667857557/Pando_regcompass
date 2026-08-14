# Single-pass common-dictionary ridge fit for condition GRNs.
#
# Candidate discovery is completed before this file is called. Every condition
# is fitted on the same frozen TF-peak-target dictionary with one target-specific
# CV lambda, shared predictor scaling and no cross-condition coefficient fusion.

.condition_fit_engine <- "condition_union_single_no_fusion_common_lambda_ridge"

.condition_ridge_contrasts <- function(
    fit, edges, condition_informative, scaling) {
    conditions <- rownames(fit$beta)
    if (length(conditions) < 2L) return(data.frame())
    pairs <- utils::combn(seq_along(conditions), 2L, simplify = FALSE)
    keep <- fit$informative_index
    p <- length(keep)
    out <- lapply(pairs, function(pair) {
        a <- pair[[1L]]
        b <- pair[[2L]]
        delta <- as.numeric(fit$beta[a, ] - fit$beta[b, ])
        estimable <- as.logical(condition_informative[a, ]) &
            as.logical(condition_informative[b, ]) & is.finite(delta)
        se <- rep(NA_real_, nrow(edges))
        if (!is.null(fit$covariance_z) && p > 0L) {
            for (j in seq_along(keep)) {
                edge_index <- keep[[j]]
                ia <- (a - 1L) * p + j
                ib <- (b - 1L) * p + j
                variance_z <- fit$covariance_z[ia, ia] +
                    fit$covariance_z[ib, ib] -
                    2 * fit$covariance_z[ia, ib]
                if (is.finite(variance_z) && variance_z >= 0) {
                    se[[edge_index]] <- sqrt(variance_z) /
                        scaling$scale[[edge_index]]
                }
            }
        }
        statistic <- delta / se
        statistic[!is.finite(statistic) | !estimable] <- NA_real_
        pval <- 2 * stats::pnorm(-abs(statistic))
        data.frame(
            tf = edges$tf,
            target = edges$target,
            region = edges$region,
            edge_id = edges$edge_id,
            atac_feature_id = edges$atac_feature_id,
            condition_a = conditions[[a]],
            condition_b = conditions[[b]],
            contrast = paste0(conditions[[a]], "-", conditions[[b]]),
            estimate_a = as.numeric(fit$beta[a, ]),
            estimate_b = as.numeric(fit$beta[b, ]),
            contrast_estimate = delta,
            contrast_se = se,
            contrast_statistic = statistic,
            contrast_pval = pval,
            contrast_estimable = estimable,
            stringsAsFactors = FALSE
        )
    })
    do.call(rbind, out)
}

.condition_ridge_target <- function(
    prepared, edges, cells_by_condition, folds, control,
    min_residual_df, rank_action) {
    x <- .condition_ridge_predictors(prepared, edges, cells_by_condition)
    y <- lapply(cells_by_condition, function(cells) {
        as.numeric(prepared$gene_data[cells, edges$target[[1L]]])
    })
    names(y) <- names(cells_by_condition)
    cv <- .condition_ridge_cv(x, y, folds, control, min_residual_df)
    scaling <- .condition_ridge_scaling(x, control$scale_floor)
    fit <- .condition_ridge_fit(
        x, y, scaling, cv$lambda, min_residual_df, inference = TRUE
    )
    if (!identical(fit$status, "ok")) {
        stop("Final condition ridge failed for target `",
             edges$target[[1L]], "`: ", fit$status, ".", call. = FALSE)
    }
    raw_rank_deficient <- fit$raw_rank < sum(fit$informative)
    if (identical(rank_action, "error") && any(raw_rank_deficient)) {
        stop("Raw common-dictionary design is rank deficient for target `",
             edges$target[[1L]],
             "`; ridge succeeded but rank_action='error' requested failure.",
             call. = FALSE)
    }

    shared <- colMeans(fit$beta)
    deviation <- sweep(fit$beta, 2L, shared, "-")
    condition_informative <- sweep(
        !fit$zero_variance, 2L, fit$informative, "&"
    )
    terms <- sprintf("edge_%07d", edges$candidate_index)
    coefficients <- fits <- vector("list", length(cells_by_condition))
    names(coefficients) <- names(fits) <- names(cells_by_condition)
    for (i in seq_along(cells_by_condition)) {
        condition <- names(cells_by_condition)[[i]]
        estimate <- as.numeric(fit$beta[i, ])
        informative_here <- as.logical(condition_informative[i, ])
        std_err <- as.numeric(fit$se[i, ])
        statistic <- as.numeric(fit$statistic[i, ])
        pval <- as.numeric(fit$pval[i, ])
        std_err[!informative_here] <- NA_real_
        statistic[!informative_here] <- NA_real_
        pval[!informative_here] <- NA_real_
        coefficients[[i]] <- data.frame(
            tf = edges$tf,
            target = edges$target,
            region = edges$region,
            term = terms,
            edge_id = edges$edge_id,
            atac_feature_id = edges$atac_feature_id,
            estimate = estimate,
            estimate_standardized = as.numeric(fit$beta_z[i, ]),
            shared_estimate = as.numeric(shared),
            condition_deviation = as.numeric(deviation[i, ]),
            std_err = std_err,
            statistic = statistic,
            pval = pval,
            estimable = informative_here & is.finite(estimate),
            zero_variance = as.logical(fit$zero_variance[i, ]),
            condition_informative = informative_here,
            aliased = FALSE,
            condition = condition,
            candidate_index = edges$candidate_index,
            source_global = edges$source_global,
            source_conditions = edges$source_conditions,
            n_sources = edges$n_sources,
            stringsAsFactors = FALSE
        )
        local_residual <- length(cells_by_condition[[i]]) - 1 -
            fit$effective_df[[i]]
        fits[[i]] <- data.frame(
            target = edges$target[[1L]],
            condition = condition,
            rsq = fit$rsq[[i]],
            rsq_oof = cv$rsq_oof[[i]],
            rsq_in_sample = fit$rsq[[i]],
            rank = as.integer(1L + fit$raw_rank[[i]]),
            raw_rank = as.integer(fit$raw_rank[[i]]),
            residual_df = as.integer(max(0, floor(local_residual))),
            effective_df = fit$effective_df[[i]],
            residual_df_effective_joint = fit$residual_df,
            condition_number = fit$raw_kappa[[i]],
            condition_number_regularized = fit$regularized_kappa,
            fit_status = "ok",
            intercept = fit$intercept[[i]],
            lambda = cv$lambda,
            lambda_min = cv$lambda_min,
            lambda_rule = control$lambda_rule,
            cv_mse = cv$cv_mse,
            cv_se = cv$cv_se,
            design_rank_deficient = raw_rank_deficient[[i]],
            nvariables = nrow(edges),
            nvariables_dictionary = nrow(edges),
            nvariables_estimable = sum(informative_here),
            n_zero_variance = sum(fit$zero_variance[i, ]),
            n_aliased = 0L,
            condition_weight = fit$weight[[i]],
            predictor_scale_reference = scaling$reference,
            stringsAsFactors = FALSE
        )
    }
    list(
        coefficients = do.call(rbind, coefficients),
        contrasts = .condition_ridge_contrasts(
            fit, edges, condition_informative, scaling
        ),
        fit = do.call(rbind, fits),
        cv = cv,
        scaling = scaling
    )
}

.condition_ridge_fit_contract_one_pass <- function(
    object, fit, prepared, control, rank_action = "mark",
    min_residual_df = 1L, parallel = FALSE, verbose = TRUE,
    progress_phase = NULL, progress_label = NULL) {
    .condition_validate_dictionary(fit$edge_dictionary, prepared)
    cells <- fit$condition_cell_ids[fit$condition_levels]
    if (any(lengths(cells) < 3L)) {
        stop("Every condition ridge fit needs at least three cells.",
             call. = FALSE)
    }
    folds <- .condition_ridge_folds(cells, control$cv_folds, control$seed)
    targets <- unique(as.character(fit$edge_dictionary$target))
    names(targets) <- targets

    if (is.null(progress_phase)) {
        progress_phase <- if (length(fit$condition_levels) == 1L) {
            "ridge_standard"
        } else {
            "ridge_condition"
        }
    }
    if (is.null(progress_label)) progress_label <- fit$cell_type %||% ""

    result <- .pando_target_payload_map(
        keys = targets,
        build_payload = function(target) {
            .pando_ridge_target_payload(
                prepared = prepared,
                edge_dictionary = fit$edge_dictionary,
                target = target,
                cells = cells,
                folds = folds,
                control = control,
                min_residual_df = min_residual_df,
                rank_action = rank_action
            )
        },
        worker_name = ".pando_ridge_target_worker",
        parallel = parallel,
        verbose = verbose,
        phase = progress_phase,
        label = progress_label
    )

    coefficient <- do.call(rbind, lapply(result, `[[`, "coefficients"))
    contrast <- do.call(rbind, lapply(result, `[[`, "contrasts"))
    fit_table <- do.call(rbind, lapply(result, `[[`, "fit"))
    rownames(coefficient) <- rownames(fit_table) <- NULL
    if (is.data.frame(contrast) && nrow(contrast)) rownames(contrast) <- NULL

    coefficient$padj <- NA_real_
    for (condition in fit$condition_levels) {
        index <- which(coefficient$condition == condition)
        valid <- index[coefficient$estimable[index] %in% TRUE &
                       is.finite(coefficient$pval[index])]
        if (length(valid)) {
            coefficient$padj[valid] <- stats::p.adjust(
                coefficient$pval[valid], method = fit$adjust_method
            )
        }
    }
    coefficient$statistically_supported <- coefficient$estimable &
        is.finite(coefficient$padj) &
        coefficient$padj < fit$padj_threshold
    coefficient$significant <- coefficient$statistically_supported
    coefficient$penalty_effect <- ifelse(
        coefficient$statistically_supported & is.finite(coefficient$estimate),
        coefficient$estimate, 0
    )
    coefficient$direction <- ifelse(
        !coefficient$estimable, "undefined",
        ifelse(coefficient$estimate > 0, "positive",
               ifelse(coefficient$estimate < 0, "negative", "zero"))
    )
    coefficient$effect_definition <-
        "no_fusion_common_lambda_ridge_condition_coefficient_raw_tf_atac_units"
    coefficient$inference_scope <-
        "approximate_ridge_wald_conditional_on_frozen_dictionary_and_cv_lambda"

    if (is.data.frame(contrast) && nrow(contrast)) {
        contrast$contrast_padj <- NA_real_
        pair_key <- paste(contrast$condition_a, contrast$condition_b, sep = "\001")
        for (key in unique(pair_key)) {
            index <- which(pair_key == key)
            valid <- index[contrast$contrast_estimable[index] %in% TRUE &
                           is.finite(contrast$contrast_pval[index])]
            if (length(valid)) {
                contrast$contrast_padj[valid] <- stats::p.adjust(
                    contrast$contrast_pval[valid], method = fit$adjust_method
                )
            }
        }
        contrast$contrast_significant <- contrast$contrast_estimable &
            is.finite(contrast$contrast_padj) &
            contrast$contrast_padj < fit$padj_threshold
        contrast$inference_scope <-
            "approximate_no_fusion_ridge_wald_contrast_diagnostic"
    }

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
                method = "condition_ridge",
                family = "gaussian_identity",
                fit_mode = "frozen_common_dictionary_no_fusion",
                condition = condition,
                edge_dictionary = fit$edge_dictionary,
                scale = FALSE,
                internal_scale_reference =
                    "equal_condition_within_condition_rms",
                exported_coefficient_scale = "raw_tf_atac_interaction_units",
                interaction = ":",
                rna_layer = prepared$rna_layer,
                peak_layer = prepared$peak_layer,
                peak_value_type = prepared$peak_value_type,
                preprocessing_fingerprint = prepared$preprocessing_fingerprint,
                adjust_method = fit$adjust_method,
                padj_threshold = fit$padj_threshold,
                projection_policy =
                    "condition_bh_supported_common_dictionary_ridge_effects",
                ridge_control = control
            )
        )
        object@grn@networks[[network_name]] <- network
        object@grn@active_network <- network_name
    }

    fit$model_schema <- .condition_multitask_ridge_schema
    fit$fit_engine <- .condition_fit_engine
    fit$coefficient_scale <- "raw_tf_atac_interaction_units"
    fit$internal_predictor_scale <- "equal_condition_within_condition_rms"
    fit$inference_scope <-
        "approximate_ridge_wald_conditional_on_frozen_dictionary_and_cv_lambda"
    fit$coefficients <- coefficient
    fit$contrasts <- contrast
    fit$fit <- fit_table
    fit$scale <- FALSE
    fit$projection_effect_column <- "penalty_effect"
    fit$projection_policy <-
        "condition_bh_supported_common_dictionary_ridge_effects"
    fit$ridge_control <- control
    fit$target_cv <- lapply(result, function(one) {
        one$cv[c("lambda", "lambda_min", "cv_mse", "cv_se",
                 "rsq_oof", "curve")]
    })
    fit$target_scaling <- lapply(result, `[[`, "scaling")
    fit$target_parallel_memory_policy <- list(
        payload = "target_specific_rna_atac_edges",
        batching = "worker_sized",
        worker_dispatch = "namespace_level_function_name",
        worker_gc = TRUE,
        master_batch_gc = TRUE
    )
    fit$rsq_definition <- "selected_lambda_full_data_R2"
    fit$rsq_oof_role <- "cross_validated_prediction_diagnostic_only"
    class(fit) <- c("ConditionGRNFit", "list")
    result <- NULL
    invisible(gc(verbose = FALSE, full = TRUE))
    list(object = object, fit = fit)
}
