# Joint multi-task ridge refit of one exact-union ConditionGRNFit skeleton.

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
        x, y, scaling, cv$lambda, control$fusion_ratio,
        min_residual_df, inference = TRUE
    )
    if (!identical(fit$status, "ok")) {
        stop("Final multi-task ridge failed for target `",
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
    terms <- sprintf("edge_%07d", edges$candidate_index)
    coefficients <- fits <- vector("list", length(cells_by_condition))
    names(coefficients) <- names(fits) <- names(cells_by_condition)
    for (i in seq_along(cells_by_condition)) {
        condition <- names(cells_by_condition)[[i]]
        estimate <- as.numeric(fit$beta[i, ])
        condition_informative <- fit$informative &
            !as.logical(fit$zero_variance[i, ])
        std_err <- as.numeric(fit$se[i, ])
        statistic <- as.numeric(fit$statistic[i, ])
        pval <- as.numeric(fit$pval[i, ])
        std_err[!condition_informative] <- NA_real_
        statistic[!condition_informative] <- NA_real_
        pval[!condition_informative] <- NA_real_
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
            estimable = condition_informative & is.finite(estimate),
            zero_variance = as.logical(fit$zero_variance[i, ]),
            condition_informative = condition_informative,
            borrowed_by_fusion = !condition_informative &
                fit$informative & is.finite(estimate) & estimate != 0,
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
            rsq = cv$rsq_oof[[i]],
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
            fusion_ratio = control$fusion_ratio,
            cv_mse = cv$cv_mse,
            cv_se = cv$cv_se,
            design_rank_deficient = raw_rank_deficient[[i]],
            nvariables_dictionary = nrow(edges),
            nvariables_estimable = sum(condition_informative),
            n_zero_variance = sum(fit$zero_variance[i, ]),
            n_aliased = 0L,
            condition_weight = fit$weight[[i]],
            stringsAsFactors = FALSE
        )
    }
    list(
        coefficients = do.call(rbind, coefficients),
        fit = do.call(rbind, fits),
        cv = cv,
        scaling = scaling
    )
}

.condition_ridge_refit_contract <- function(
    object, fit, prepared, control, rank_action = "mark",
    min_residual_df = 1L, parallel = FALSE, verbose = TRUE) {
    .condition_validate_dictionary(fit$edge_dictionary, prepared)
    cells <- fit$condition_cell_ids[fit$condition_levels]
    if (any(lengths(cells) < 3L)) {
        stop("Every multi-task ridge condition needs at least three cells.",
             call. = FALSE)
    }
    folds <- .condition_ridge_folds(cells, control$cv_folds, control$seed)
    targets <- unique(as.character(fit$edge_dictionary$target))
    names(targets) <- targets
    result <- map_par(targets, function(target) {
        edges <- fit$edge_dictionary[
            fit$edge_dictionary$target == target, , drop = FALSE
        ]
        .condition_ridge_target(
            prepared, edges, cells, folds, control,
            min_residual_df, rank_action
        )
    }, verbose = verbose, parallel = parallel)

    coefficient <- do.call(rbind, lapply(result, `[[`, "coefficients"))
    fit_table <- do.call(rbind, lapply(result, `[[`, "fit"))
    rownames(coefficient) <- rownames(fit_table) <- NULL

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
    coefficient$significant <- coefficient$estimable &
        is.finite(coefficient$padj) &
        coefficient$padj < fit$padj_threshold
    coefficient$penalty_effect <- ifelse(
        coefficient$significant, coefficient$estimate, 0
    )
    coefficient$direction <- ifelse(
        !coefficient$estimable, "undefined",
        ifelse(coefficient$estimate > 0, "positive",
               ifelse(coefficient$estimate < 0, "negative", "zero"))
    )
    coefficient$effect_definition <-
        "multitask_ridge_condition_coefficient_raw_tf_atac_units"
    coefficient$inference_scope <-
        "ridge_wald_conditional_on_dictionary_cv_lambda_and_fusion"

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
                method = "multitask_ridge",
                family = "gaussian_identity",
                fit_mode = "fixed_edge_dictionary_joint_conditions",
                condition = condition,
                edge_dictionary = fit$edge_dictionary,
                scale = FALSE,
                internal_scale_reference = "pooled_conditions_common_zscore",
                exported_coefficient_scale = "raw_tf_atac_interaction_units",
                interaction = ":",
                rna_layer = prepared$rna_layer,
                peak_layer = prepared$peak_layer,
                peak_value_type = prepared$peak_value_type,
                preprocessing_fingerprint = prepared$preprocessing_fingerprint,
                adjust_method = fit$adjust_method,
                padj_threshold = fit$padj_threshold,
                ridge_control = control
            )
        )
        object@grn@networks[[network_name]] <- network
        object@grn@active_network <- network_name
    }

    fit$model_schema <- .condition_multitask_ridge_schema
    fit$fit_engine <- "two_stage_exact_edge_union_multitask_ridge"
    fit$coefficient_scale <- "raw_tf_atac_interaction_units"
    fit$internal_predictor_scale <- "pooled_conditions_common_zscore"
    fit$inference_scope <-
        "ridge_wald_conditional_on_dictionary_cv_lambda_and_fusion"
    fit$coefficients <- coefficient
    fit$fit <- fit_table
    fit$scale <- FALSE
    fit$projection_effect_column <- "penalty_effect"
    fit$projection_policy <- "padj_significant_effects_only"
    fit$ridge_control <- control
    fit$target_cv <- lapply(result, function(one) {
        one$cv[c("lambda", "lambda_min", "cv_mse", "cv_se",
                 "rsq_oof", "curve")]
    })
    fit$target_scaling <- lapply(result, `[[`, "scaling")
    class(fit) <- c("ConditionGRNFit", "list")
    list(object = object, fit = fit)
}
