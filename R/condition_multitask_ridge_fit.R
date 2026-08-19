# Single-pass common-dictionary Scheme-E fit for condition GRNs.
#
# Candidate discovery is completed before this file is called. Every condition
# is fitted on the same frozen TF-peak-target dictionary and common predictor
# scale. Scheme E (z=0.25) couples condition coefficients through an exact-edge
# sparse deviation penalty and returns continuous coefficients for projection.

.condition_fit_engine <- "condition_union_scheme_e_exact_edge_z025"

.condition_profile_information_export <- function(fit) {
    information <- as.numeric(fit$profile_information)
    if (nrow(fit$beta) == 2L) information <- information / 2
    information
}

.condition_ridge_contrasts <- function(
    fit, edges, condition_informative, scaling) {
    conditions <- rownames(fit$beta)
    if (length(conditions) < 2L) return(data.frame())
    profile_information <- .condition_profile_information_export(fit)
    pairs <- utils::combn(seq_along(conditions), 2L, simplify = FALSE)
    out <- lapply(pairs, function(pair) {
        a <- pair[[1L]]
        b <- pair[[2L]]
        delta <- as.numeric(fit$beta[a, ] - fit$beta[b, ])
        identifiable <- as.logical(fit$contrast_identifiable) & is.finite(delta)
        se <- sqrt(as.numeric(fit$se[a, ])^2 + as.numeric(fit$se[b, ])^2)
        se[!is.finite(se) | se <= 0] <- NA_real_
        statistic <- delta / se
        statistic[!is.finite(statistic) | !identifiable] <- NA_real_
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
            contrast_estimable = identifiable,
            contrast_identifiable = identifiable,
            shared_by_boundary = as.logical(fit$shared_by_boundary),
            fused_by_penalty = as.logical(fit$fused_by_penalty),
            profile_information_delta = profile_information,
            penalty_family = fit$penalty_family,
            penalty_value = fit$penalty_value,
            solver_status = fit$solver_status,
            kkt_residual = fit$kkt_residual,
            iterations = fit$iterations,
            stringsAsFactors = FALSE
        )
    })
    do.call(rbind, out)
}

.condition_ridge_target <- function(
    prepared, edges, cells_by_condition, folds = NULL, control,
    min_residual_df, rank_action) {
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
        control = control
    )
    if (!identical(fit$status, "ok")) {
        detail <- if (!is.null(fit$solver_status)) {
            paste0(" (solver=", fit$solver_status,
                   ", kkt=", format(fit$kkt_residual, digits = 4), ")")
        } else ""
        stop("Final condition Scheme-E fit failed for target `",
             edges$target[[1L]], "`: ", fit$status, detail, ".",
             call. = FALSE)
    }
    raw_rank_deficient <- fit$raw_rank < sum(fit$informative)
    if (identical(rank_action, "error") && any(raw_rank_deficient)) {
        stop("Raw common-dictionary design is rank deficient for target `",
             edges$target[[1L]],
             "`; Scheme E retained only identifiable condition contrasts but ",
             "rank_action='error' requested failure.", call. = FALSE)
    }

    shared <- fit$shared_z / scaling$scale
    deviation <- sweep(fit$deviation_z, 2L, scaling$scale, "/")
    profile_information <- .condition_profile_information_export(fit)
    profile_information_definition <- if (length(cells_by_condition) == 2L) {
        "pairwise_delta_profile_information"
    } else {
        "minimum_orthonormal_condition_contrast_information"
    }
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
            beta_shared = as.numeric(shared),
            condition_deviation = as.numeric(deviation[i, ]),
            delta_beta = as.numeric(deviation[i, ]),
            std_err = std_err,
            statistic = statistic,
            pval = pval,
            estimable = is.finite(estimate),
            zero_variance = as.logical(fit$zero_variance[i, ]),
            condition_informative = informative_here,
            contrast_identifiable = as.logical(fit$contrast_identifiable),
            shared_by_boundary = as.logical(fit$shared_by_boundary),
            fused_by_penalty = as.logical(fit$fused_by_penalty),
            raw_information_condition = as.numeric(fit$raw_information[i, ]),
            profile_information_delta = profile_information,
            profile_information_definition = profile_information_definition,
            penalty_family = fit$penalty_family,
            penalty_value = fit$penalty_value,
            solver_status = fit$solver_status,
            kkt_residual = fit$kkt_residual,
            iterations = fit$iterations,
            aliased = !as.logical(fit$contrast_identifiable),
            condition = condition,
            candidate_index = edges$candidate_index,
            source_global = edges$source_global,
            source_conditions = edges$source_conditions,
            n_sources = edges$n_sources,
            stringsAsFactors = FALSE
        )
        local_residual <- length(cells_by_condition[[i]]) - 1 - fit$raw_rank[[i]]
        fits[[i]] <- data.frame(
            target = edges$target[[1L]],
            condition = condition,
            rsq = fit$rsq[[i]],
            rank = as.integer(1L + fit$raw_rank[[i]]),
            raw_rank = as.integer(fit$raw_rank[[i]]),
            residual_df = as.integer(max(0, floor(local_residual))),
            residual_df_effective_joint = fit$residual_df,
            condition_number = fit$raw_kappa[[i]],
            fit_status = "ok",
            intercept = fit$intercept[[i]],
            sigma2_common = fit$sigma2,
            deviation_z = fit$penalty_value,
            penalty_family = fit$penalty_family,
            solver_status = fit$solver_status,
            kkt_residual = fit$kkt_residual,
            iterations = fit$iterations,
            nvariables = nrow(edges),
            nvariables_dictionary = nrow(edges),
            nvariables_contrast_identifiable =
                sum(fit$contrast_identifiable %in% TRUE),
            n_shared_by_boundary = sum(fit$shared_by_boundary %in% TRUE),
            n_fused_by_penalty = sum(fit$fused_by_penalty %in% TRUE),
            n_zero_variance = sum(fit$zero_variance[i, ]),
            condition_weight = 1,
            predictor_scale_reference = scaling$reference,
            profile_information_definition = profile_information_definition,
            stringsAsFactors = FALSE
        )
    }
    list(
        coefficients = do.call(rbind, coefficients),
        contrasts = .condition_ridge_contrasts(
            fit, edges, condition_informative, scaling
        ),
        fit = do.call(rbind, fits),
        scaling = scaling,
        solver = list(
            status = fit$solver_status,
            kkt_residual = fit$kkt_residual,
            iterations = fit$iterations,
            penalty_family = fit$penalty_family,
            penalty_value = fit$penalty_value
        )
    )
}

.condition_ridge_fit_contract_one_pass <- function(
    object, fit, prepared, control, rank_action = "mark",
    min_residual_df = 1L, parallel = FALSE, verbose = TRUE,
    progress_phase = NULL, progress_label = NULL) {
    .condition_validate_dictionary(fit$edge_dictionary, prepared)
    cells <- fit$condition_cell_ids[fit$condition_levels]
    if (length(cells) < 2L) {
        stop("Conditional Scheme E requires at least two conditions.", call. = FALSE)
    }
    if (any(lengths(cells) < 3L)) {
        stop("Every condition Scheme-E fit needs at least three cells.",
             call. = FALSE)
    }
    targets <- unique(as.character(fit$edge_dictionary$target))
    names(targets) <- targets
    if (is.null(progress_phase)) progress_phase <- "scheme_e_condition"
    if (is.null(progress_label)) progress_label <- fit$cell_type %||% ""

    result <- .pando_target_payload_map(
        keys = targets,
        build_payload = function(target) {
            .pando_ridge_target_payload(
                prepared = prepared,
                edge_dictionary = fit$edge_dictionary,
                target = target,
                cells = cells,
                folds = NULL,
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
        valid <- index[is.finite(coefficient$pval[index])]
        if (length(valid)) {
            coefficient$padj[valid] <- stats::p.adjust(
                coefficient$pval[valid], method = fit$adjust_method
            )
        }
    }
    coefficient$statistically_supported <-
        is.finite(coefficient$padj) & coefficient$padj < fit$padj_threshold
    coefficient$significant <- coefficient$statistically_supported
    coefficient$penalty_effect <- coefficient$estimate
    coefficient$direction <- ifelse(
        coefficient$estimate > 0, "positive",
        ifelse(coefficient$estimate < 0, "negative", "zero")
    )
    coefficient$effect_definition <-
        "scheme_e_z025_continuous_condition_coefficient_raw_tf_atac_units"
    coefficient$inference_scope <-
        "diagnostic_wald_only;scheme_e_coefficients_define_projection"

    if (is.data.frame(contrast) && nrow(contrast)) {
        contrast$contrast_padj <- NA_real_
        pair_key <- paste(contrast$condition_a, contrast$condition_b, sep = "\001")
        for (key in unique(pair_key)) {
            index <- which(pair_key == key)
            valid <- index[
                contrast$contrast_identifiable[index] %in% TRUE &
                is.finite(contrast$contrast_pval[index])
            ]
            if (length(valid)) {
                contrast$contrast_padj[valid] <- stats::p.adjust(
                    contrast$contrast_pval[valid], method = fit$adjust_method
                )
            }
        }
        contrast$contrast_significant <-
            contrast$contrast_identifiable &
            is.finite(contrast$contrast_padj) &
            contrast$contrast_padj < fit$padj_threshold
        contrast$inference_scope <-
            "diagnostic_wald_contrast_on_scheme_e_continuous_coefficients"
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
                method = "condition_scheme_e",
                family = "gaussian_identity",
                fit_mode = "frozen_common_dictionary_scheme_e_z025",
                condition = condition,
                edge_dictionary = fit$edge_dictionary,
                scale = FALSE,
                internal_scale_reference = "equal_condition_within_condition_rms",
                exported_coefficient_scale = "raw_tf_atac_interaction_units",
                interaction = ":",
                rna_layer = prepared$rna_layer,
                peak_layer = prepared$peak_layer,
                peak_value_type = prepared$peak_value_type,
                preprocessing_fingerprint = prepared$preprocessing_fingerprint,
                adjust_method = fit$adjust_method,
                padj_threshold = fit$padj_threshold,
                projection_policy = "continuous_common_dictionary_scheme_e_effects",
                deviation_penalty = list(
                    family = .condition_scheme_e_penalty_family,
                    z = .condition_scheme_e_z
                ),
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
    fit$inference_scope <- "scheme_e_z025_primary;BH_and_R2_are_diagnostics_only"
    fit$coefficients <- coefficient
    fit$contrasts <- contrast
    fit$fit <- fit_table
    fit$scale <- FALSE
    fit$projection_effect_column <- "penalty_effect"
    fit$projection_policy <- "continuous_common_dictionary_scheme_e_effects"
    fit$ridge_control <- control
    fit$deviation_penalty <- list(
        family = .condition_scheme_e_penalty_family,
        z = .condition_scheme_e_z
    )
    fit$target_solver <- lapply(result, `[[`, "solver")
    fit$target_scaling <- lapply(result, `[[`, "scaling")
    fit$target_parallel_memory_policy <- list(
        payload = "target_specific_rna_atac_edges",
        batching = "worker_sized",
        worker_dispatch = "namespace_level_function_name",
        worker_gc = TRUE,
        master_batch_gc = TRUE
    )
    fit$rsq_definition <- "scheme_e_z025_full_data_R2_diagnostic"
    class(fit) <- c("ConditionGRNFit", "list")
    result <- NULL
    invisible(gc(verbose = FALSE, full = TRUE))
    list(object = object, fit = fit)
}
