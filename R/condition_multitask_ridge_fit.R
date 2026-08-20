# One-pass common-dictionary E-star/JSE fit for conditional GRNs.

.condition_fit_engine <- "condition_union_Estar_z025_jointse"

.condition_ridge_contrasts <- function(
    fit, edges, scaling) {
    conditions <- rownames(fit$beta)
    if (length(conditions) < 2L) return(data.frame())
    pairs <- utils::combn(seq_along(conditions), 2L, simplify = FALSE)
    keep <- fit$informative_index
    out <- list()
    for (pair in pairs) {
        a <- pair[[1L]]
        b <- pair[[2L]]
        for (edge in seq_len(nrow(edges))) {
            production_delta <- as.numeric(
                fit$beta[b, edge] - fit$beta[a, edge]
            )
            local_edge <- match(edge, keep)
            joint <- if (!is.na(local_edge)) {
                .condition_joint_contrast(
                    fit$joint_inference, scaling, keep,
                    condition_a = a, condition_b = b,
                    local_edge = local_edge
                )
            } else {
                list(
                    estimate = NA_real_, se = NA_real_,
                    statistic = NA_real_, pval = NA_real_,
                    estimable = FALSE, same_component = FALSE,
                    boundary_shared = FALSE
                )
            }
            contrast_status <- if (isTRUE(joint$same_component) &&
                                   isTRUE(joint$boundary_shared)) {
                "shared_by_boundary"
            } else if (isTRUE(joint$same_component)) {
                "fused_by_E"
            } else if (!isTRUE(joint$estimable)) {
                "not_estimable"
            } else {
                "ok"
            }
            out[[length(out) + 1L]] <- data.frame(
                tf = edges$tf[[edge]],
                target = edges$target[[edge]],
                region = edges$region[[edge]],
                edge_id = edges$edge_id[[edge]],
                atac_feature_id = edges$atac_feature_id[[edge]],
                condition_a = conditions[[a]],
                condition_b = conditions[[b]],
                contrast = paste0(conditions[[b]], "-", conditions[[a]]),
                estimate_a = as.numeric(fit$beta[a, edge]),
                estimate_b = as.numeric(fit$beta[b, edge]),
                contrast_estimate = production_delta,
                inference_contrast_estimate = joint$estimate,
                contrast_se = joint$se,
                contrast_statistic = joint$statistic,
                contrast_pval = joint$pval,
                contrast_estimable = isTRUE(joint$estimable) &&
                    !isTRUE(joint$same_component),
                contrast_identifiable =
                    fit$contrast_identifiable[[edge]] %in% TRUE,
                contrast_status = contrast_status,
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
                inference_schema = fit$inference_schema,
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
        x = x,
        y = y,
        scaling = scaling,
        min_residual_df = min_residual_df,
        inference = TRUE,
        control = control,
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
            "`; E-star marks non-identifiable contrasts instead of drifting, ",
            "but rank_action='error' requested failure.",
            call. = FALSE
        )
    }

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
        inference_estimate <- as.numeric(fit$inference_estimate[i, ])
        inference_se <- as.numeric(fit$inference_se[i, ])
        inference_statistic <- as.numeric(fit$inference_statistic[i, ])
        pval <- as.numeric(fit$inference_pval[i, ])
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
            contrast_coordinate = NA_character_,
            contrast_identifiable =
                as.logical(fit$contrast_identifiable),
            shared_by_boundary =
                as.logical(fit$shared_by_boundary),
            boundary_condition =
                as.logical(fit$boundary_condition[i, ]),
            fused_by_penalty =
                as.logical(fit$fused_by_penalty),
            fusion_component_id =
                as.character(fit$fusion_component_id[i, ]),
            shared_edge = as.logical(fit$shared_edge),
            raw_information_condition =
                as.numeric(fit$raw_information[i, ]),
            profile_information =
                as.numeric(fit$profile_information),
            profile_information_delta =
                as.numeric(fit$profile_information),
            condition_number = fit$raw_kappa[[i]],
            sigma2_common = fit$sigma2,
            inference_schema = fit$inference_schema,
            inference_component_id =
                as.character(fit$inference_component_id[i, ]),
            inference_hypothesis_id =
                as.character(fit$inference_hypothesis_id[i, ]),
            inference_estimate = inference_estimate,
            inference_se = inference_se,
            inference_statistic = inference_statistic,
            inference_estimable =
                as.logical(fit$inference_estimable[i, ]),
            std_err = inference_se,
            statistic = inference_statistic,
            pval = pval,
            estimable = is.finite(estimate),
            zero_variance = as.logical(fit$zero_variance[i, ]),
            condition_informative =
                as.logical(condition_informative[i, ]),
            aliased = !as.logical(fit$contrast_identifiable),
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
            inference_sigma2 = fit$inference_sigma2,
            inference_residual_df = fit$inference_residual_df,
            inference_rank = fit$inference_rank,
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
            n_shared_by_boundary =
                sum(fit$shared_by_boundary %in% TRUE),
            n_fused_by_penalty =
                sum(fit$fused_by_penalty %in% TRUE),
            n_shared_edges = sum(fit$shared_edge %in% TRUE),
            n_zero_variance = sum(fit$zero_variance[i, ]),
            condition_weight = 1,
            predictor_scale_reference = scaling$reference,
            inference_schema = fit$inference_schema,
            orthogonality_error = fit$orthogonality_error,
            dr_error = fit$dr_error,
            stringsAsFactors = FALSE
        )
    }
    list(
        coefficients = do.call(rbind, coefficients),
        contrasts = .condition_ridge_contrasts(fit, edges, scaling),
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

.condition_bh_by_condition_target <- function(
    coefficient, method = "BH") {
    coefficient$padj <- NA_real_
    key <- paste(
        as.character(coefficient$condition),
        as.character(coefficient$target),
        sep = "\001"
    )
    for (family in unique(key)) {
        index <- which(key == family)
        valid <- index[
            coefficient$inference_estimable[index] %in% TRUE &
            is.finite(as.numeric(coefficient$pval[index]))
        ]
        if (length(valid)) {
            coefficient$padj[valid] <- stats::p.adjust(
                as.numeric(coefficient$pval[valid]), method = method
            )
        }
    }
    coefficient$bh_scope <- "condition_target_BH"
    family_size <- integer(nrow(coefficient))
    for (family in unique(key)) {
        index <- which(key == family)
        family_size[index] <- sum(
            coefficient$inference_estimable[index] %in% TRUE &
            is.finite(as.numeric(coefficient$pval[index]))
        )
    }
    coefficient$bh_family_size <- family_size
    coefficient
}

.condition_ridge_fit_contract_one_pass <- function(
    object, fit, prepared, control, rank_action = "mark",
    min_residual_df = 1L, parallel = FALSE, verbose = TRUE,
    progress_phase = NULL, progress_label = NULL) {
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
    if (is.null(progress_phase)) progress_phase <- "condition_Estar_JSE"
    if (is.null(progress_label)) progress_label <- fit$cell_type %||% ""

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
        label = progress_label
    )

    coefficient <- do.call(rbind, lapply(result, `[[`, "coefficients"))
    contrast <- do.call(rbind, lapply(result, `[[`, "contrasts"))
    fit_table <- do.call(rbind, lapply(result, `[[`, "fit"))
    rownames(coefficient) <- rownames(fit_table) <- NULL
    if (is.data.frame(contrast) && nrow(contrast)) rownames(contrast) <- NULL

    coefficient <- .condition_bh_by_condition_target(
        coefficient, method = fit$adjust_method
    )
    coefficient$condition_significant <-
        coefficient$inference_estimable %in% TRUE &
        is.finite(coefficient$padj) &
        coefficient$padj < fit$padj_threshold
    coefficient$statistically_supported <- coefficient$condition_significant
    coefficient$significant <- coefficient$condition_significant
    coefficient$penalty_effect <- coefficient$estimate
    coefficient$direction <- ifelse(
        coefficient$estimate > 0, "positive",
        ifelse(coefficient$estimate < 0, "negative", "zero")
    )
    coefficient$effect_definition <-
        "E_star_z025_continuous_condition_coefficient_raw_tf_atac_units"
    coefficient$inference_scope <-
        "fusion_component_joint_refit_selected_structure_wald"

    if (is.data.frame(contrast) && nrow(contrast)) {
        contrast$contrast_padj <- NA_real_
        key <- paste(
            contrast$condition_a, contrast$condition_b,
            contrast$target, sep = "\001"
        )
        for (family in unique(key)) {
            index <- which(key == family)
            valid <- index[
                contrast$contrast_estimable[index] %in% TRUE &
                is.finite(contrast$contrast_pval[index])
            ]
            if (length(valid)) {
                contrast$contrast_padj[valid] <- stats::p.adjust(
                    contrast$contrast_pval[valid],
                    method = fit$adjust_method
                )
            }
        }
        contrast$contrast_significant <-
            contrast$contrast_estimable %in% TRUE &
            is.finite(contrast$contrast_padj) &
            contrast$contrast_padj < fit$padj_threshold
        contrast$inference_scope <-
            "fusion_component_joint_refit_pairwise_contrast"
    }

    fit$model_schema <- .condition_multitask_ridge_schema
    fit$fit_engine <- .condition_fit_engine
    fit$coefficient_scale <- "raw_tf_atac_interaction_units"
    fit$internal_predictor_scale <- "equal_condition_within_condition_rms"
    fit$inference_schema <- .condition_E_star_inference_schema
    fit$inference_scope <-
        "E_star_z025_primary_fusion_component_joint_refit"
    fit$coefficients <- coefficient
    fit$contrasts <- contrast
    fit$fit <- fit_table
    fit$scale <- FALSE
    fit$projection_effect_column <- "penalty_effect"
    fit$projection_policy <- "any_condition_padj_exact_edge_union"
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
                method = "condition_Estar_JSE",
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
                preprocessing_fingerprint =
                    prepared$preprocessing_fingerprint,
                adjust_method = fit$adjust_method,
                padj_threshold = fit$padj_threshold,
                bh_scope = "condition_target_BH",
                projection_policy =
                    "any_condition_padj_exact_edge_union",
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
