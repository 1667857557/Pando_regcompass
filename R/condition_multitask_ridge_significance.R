# Statistical support and RegCompass handoff for conditional E-star GRNs.
#
# Candidate correlation support defines the frozen exact-edge dictionary.
# E-star z=0.25 supplies continuous condition-specific production coefficients.
# Formal inference is fit separately without fusion on that frozen dictionary.
# Each exact edge receives one omnibus P value across its estimable conditions;
# BH is then applied once across the complete exact-edge network for the cell
# type. The resulting edge topology is common to every retained condition.

.condition_projection_policy <-
    "exact_edge_whole_network_BH_common_topology"
.condition_fit_dictionary_policy <-
    "global_and_condition_union_pando_correlation_supported_frozen_dictionary"

.condition_validate_adjust_method <- function(adjust_method) {
    value <- toupper(as.character(adjust_method))
    if (length(value) != 1L || is.na(value) || !identical(value, "BH")) {
        stop(
            "Conditional exact-edge inference requires `adjust_method = \"BH\"`.",
            call. = FALSE
        )
    }
    "BH"
}

.condition_validate_padj_threshold <- function(padj_threshold) {
    if (!is.numeric(padj_threshold) || length(padj_threshold) != 1L ||
        !is.finite(padj_threshold) || padj_threshold <= 0 ||
        padj_threshold >= 1) {
        stop(
            "`padj_threshold` must be one finite number in (0, 1).",
            call. = FALSE
        )
    }
    as.numeric(padj_threshold)
}

.condition_exact_edge_inference <- function(fit, coefficient) {
    edge <- as.data.frame(fit$edge_dictionary, stringsAsFactors = FALSE)
    required <- c(
        "edge_id", "target", "tf", "region", "atac_feature_id",
        "candidate_index"
    )
    if (!all(required %in% colnames(edge)) || anyDuplicated(edge$edge_id)) {
        stop("Frozen exact-edge dictionary is incomplete.", call. = FALSE)
    }
    rows <- vector("list", nrow(edge))
    for (i in seq_len(nrow(edge))) {
        edge_id <- as.character(edge$edge_id[[i]])
        index <- which(as.character(coefficient$edge_id) == edge_id)
        if (length(index) != length(fit$condition_levels) ||
            anyDuplicated(as.character(coefficient$condition[index]))) {
            stop("Every exact edge must occur once in every fitted condition.",
                 call. = FALSE)
        }
        valid <- index[
            coefficient$condition_inference_estimable[index] %in% TRUE &
            is.finite(as.numeric(coefficient$inference_estimate[index])) &
            is.finite(as.numeric(coefficient$inference_variance[index])) &
            as.numeric(coefficient$inference_variance[index]) > 0 &
            is.finite(as.numeric(coefficient$condition_pval[index]))
        ]
        m <- length(valid)
        statistic <- pval <- NA_real_
        test <- "not_estimable"
        if (m == 1L) {
            one <- valid[[1L]]
            statistic <- as.numeric(coefficient$inference_statistic[[one]])^2
            pval <- as.numeric(coefficient$condition_pval[[one]])
            test <- "single_condition_exact_t"
        } else if (m > 1L) {
            beta <- as.numeric(coefficient$inference_estimate[valid])
            variance <- as.numeric(coefficient$inference_variance[valid])
            statistic <- sum(beta^2 / variance)
            pval <- stats::pchisq(statistic, df = m, lower.tail = FALSE)
            test <- "independent_condition_wald_chisq"
        }
        production_valid <- all(
            as.character(coefficient$fit_status[index]) == "ok" &
            is.finite(as.numeric(coefficient$penalty_effect[index]))
        )
        inference_conditions <- fit$condition_levels[
            fit$condition_levels %in% as.character(coefficient$condition[valid])
        ]
        rows[[i]] <- data.frame(
            edge_id = edge_id,
            target = as.character(edge$target[[i]]),
            tf = as.character(edge$tf[[i]]),
            region = as.character(edge$region[[i]]),
            atac_feature_id = as.character(edge$atac_feature_id[[i]]),
            candidate_index = as.integer(edge$candidate_index[[i]]),
            edge_df = as.integer(m),
            edge_statistic = statistic,
            edge_pval = pval,
            edge_padj = NA_real_,
            edge_inference_estimable = m > 0L,
            edge_inference_test = test,
            inference_conditions = paste(inference_conditions, collapse = ";"),
            all_conditions_fit_valid = production_valid,
            edge_supported = FALSE,
            stringsAsFactors = FALSE
        )
    }
    out <- do.call(rbind, rows)
    valid <- which(
        out$edge_inference_estimable %in% TRUE & is.finite(out$edge_pval)
    )
    if (length(valid)) {
        out$edge_padj[valid] <- stats::p.adjust(
            out$edge_pval[valid], method = fit$adjust_method
        )
    }
    threshold <- .condition_validate_padj_threshold(fit$padj_threshold)
    out$edge_supported <-
        out$all_conditions_fit_valid %in% TRUE &
        out$edge_inference_estimable %in% TRUE &
        is.finite(out$edge_padj) & out$edge_padj < threshold
    out$bh_scope <- "exact_edge_whole_cell_type_network_BH"
    out$bh_family_size <- length(valid)
    out$padj_threshold <- threshold
    out
}

.condition_apply_activity_gate <- function(fit) {
    if (!inherits(fit, "ConditionGRNFit")) {
        stop("A ConditionGRNFit is required for condition activity annotation.",
             call. = FALSE)
    }
    .condition_validate_padj_threshold(fit$padj_threshold)
    coefficient <- as.data.frame(fit$coefficients, stringsAsFactors = FALSE)
    support <- fit$dictionary_support_table
    fit_table <- as.data.frame(fit$fit, stringsAsFactors = FALSE)
    required_coefficient <- c(
        "edge_id", "target", "condition", "estimate", "penalty_effect",
        "inference_estimate", "inference_se", "inference_variance",
        "inference_statistic", "condition_pval",
        "condition_inference_estimable"
    )
    required_support <- c(
        "edge_id", "source_type", "condition",
        "peak_target_cor", "tf_target_cor"
    )
    required_fit <- c("target", "condition", "fit_status")
    if (!all(required_coefficient %in% colnames(coefficient)) ||
        !is.data.frame(support) ||
        !all(required_support %in% colnames(support)) ||
        !all(required_fit %in% colnames(fit_table))) {
        stop("Conditional E-star/inference metadata are incomplete.",
             call. = FALSE)
    }

    global <- support[support$source_type == "global", , drop = FALSE]
    local <- support[support$source_type == "condition", , drop = FALSE]
    if (anyDuplicated(as.character(global$edge_id))) {
        stop("Global Pando support rows must be unique by exact edge.",
             call. = FALSE)
    }
    local_key <- paste(
        as.character(local$edge_id), as.character(local$condition), sep = "\001"
    )
    if (anyNA(local_key) || anyDuplicated(local_key)) {
        stop("Condition-local Pando support rows must be unique.",
             call. = FALSE)
    }
    coefficient_key <- paste(
        as.character(coefficient$edge_id),
        as.character(coefficient$condition), sep = "\001"
    )
    local_index <- match(coefficient_key, local_key)
    global_index <- match(
        as.character(coefficient$edge_id), as.character(global$edge_id)
    )
    local_support <- !is.na(local_index)
    global_support <- !is.na(global_index)
    coefficient$peak_target_cor <- NA_real_
    coefficient$tf_target_cor <- NA_real_
    coefficient$peak_target_cor[local_support] <-
        as.numeric(local$peak_target_cor[local_index[local_support]])
    coefficient$tf_target_cor[local_support] <-
        as.numeric(local$tf_target_cor[local_index[local_support]])
    coefficient$global_peak_target_cor <- NA_real_
    coefficient$global_tf_target_cor <- NA_real_
    coefficient$global_peak_target_cor[global_support] <-
        as.numeric(global$peak_target_cor[global_index[global_support]])
    coefficient$global_tf_target_cor[global_support] <-
        as.numeric(global$tf_target_cor[global_index[global_support]])
    coefficient$local_support <- local_support
    coefficient$global_support <- global_support
    coefficient$dictionary_support <- TRUE
    coefficient$current_scope_correlation_support <-
        local_support | global_support
    coefficient$peak_cor_pass <- local_support
    coefficient$tf_cor_pass <- local_support

    fit_key <- paste(
        as.character(fit_table$target),
        as.character(fit_table$condition), sep = "\001"
    )
    if (anyDuplicated(fit_key)) {
        stop("Target-condition fit diagnostics must be unique.",
             call. = FALSE)
    }
    coefficient_fit_key <- paste(
        as.character(coefficient$target),
        as.character(coefficient$condition), sep = "\001"
    )
    fit_index <- match(coefficient_fit_key, fit_key)
    if (anyNA(fit_index)) {
        stop("Coefficient rows cannot be aligned to target-condition fits.",
             call. = FALSE)
    }
    coefficient$fit_status <- as.character(fit_table$fit_status[fit_index])

    estimate <- suppressWarnings(as.numeric(coefficient$estimate))
    coefficient$pando_estimation_active <- is.finite(estimate)
    coefficient$active <- coefficient$pando_estimation_active
    coefficient$penalty_effect <- ifelse(
        coefficient$pando_estimation_active, estimate, NA_real_
    )

    edge_inference <- .condition_exact_edge_inference(fit, coefficient)
    edge_index <- match(
        as.character(coefficient$edge_id), as.character(edge_inference$edge_id)
    )
    if (anyNA(edge_index)) {
        stop("Coefficient rows cannot be aligned to exact-edge inference.",
             call. = FALSE)
    }
    coefficient$edge_df <- edge_inference$edge_df[edge_index]
    coefficient$edge_statistic <- edge_inference$edge_statistic[edge_index]
    coefficient$edge_pval <- edge_inference$edge_pval[edge_index]
    coefficient$edge_padj <- edge_inference$edge_padj[edge_index]
    coefficient$edge_inference_estimable <-
        edge_inference$edge_inference_estimable[edge_index]
    coefficient$edge_inference_test <-
        edge_inference$edge_inference_test[edge_index]
    coefficient$inference_conditions <-
        edge_inference$inference_conditions[edge_index]
    coefficient$all_conditions_fit_valid <-
        edge_inference$all_conditions_fit_valid[edge_index]
    coefficient$edge_supported <- edge_inference$edge_supported[edge_index]
    coefficient$active_in_regcompass <- coefficient$edge_supported
    coefficient$statistically_supported <- coefficient$edge_supported
    coefficient$significant <- coefficient$edge_supported
    coefficient$pval <- coefficient$edge_pval
    coefficient$padj <- coefficient$edge_padj
    coefficient$bh_scope <- edge_inference$bh_scope[edge_index]
    coefficient$bh_family_size <- edge_inference$bh_family_size[edge_index]

    fit$coefficients <- coefficient
    fit$edge_inference <- edge_inference
    fit$projection_effect_column <- "penalty_effect"
    fit$projection_policy <- .condition_projection_policy
    fit$dictionary_support_role <-
        "global_plus_condition_candidate_provenance"
    fit$local_support_role <-
        "condition_specific_candidate_provenance_only"
    fit$global_support_role <-
        "pooled_candidate_provenance_only"
    fit$statistical_support_role <-
        "no_fusion_condition_local_gaussian_lm_then_exact_edge_whole_network_BH"
    fit$activity_role <-
        "Pando_keeps_continuous_Estar_estimates;RegCompass_uses_common_exact_edge_topology"
    fit$regcompass_edge_gate <-
        "all_conditions_fit_valid AND exact_edge_whole_network_BH_padj_lt_threshold"
    fit
}

.condition_ridge_fit_contract <- function(
    object, fit, prepared, control, rank_action = "mark",
    min_residual_df = 1L, parallel = FALSE, verbose = TRUE,
    checkpoint_dir = NULL, resume = TRUE) {
    fit$adjust_method <- .condition_validate_adjust_method(fit$adjust_method)
    fit$padj_threshold <- .condition_validate_padj_threshold(fit$padj_threshold)
    control <- .condition_E_star_control(control)
    dictionary <- fit$edge_dictionary
    progress_label <- fit$cell_type %||% ""

    if (isTRUE(verbose)) {
        message(
            "Pando condition phase=E_star_z025_start",
            " | cell_type=", as.character(progress_label),
            ";fit_edges=", nrow(dictionary),
            ";targets=", length(unique(as.character(dictionary$target))),
            ";reference=", fit$reference_condition
        )
    }
    final <- .condition_ridge_fit_contract_one_pass(
        object = object, fit = fit, prepared = prepared, control = control,
        rank_action = rank_action, min_residual_df = min_residual_df,
        parallel = parallel, verbose = verbose,
        progress_phase = "E_star_z025_and_no_fusion_inference",
        progress_label = progress_label,
        checkpoint_dir = checkpoint_dir, resume = resume
    )
    final$fit$fit_dictionary_policy <- .condition_fit_dictionary_policy
    final$fit$candidate_edge_count <- nrow(dictionary)
    final$fit$fit_dictionary_edge_count <- nrow(dictionary)
    final$fit$edge_dictionary <- dictionary
    final$fit$dictionary_support_table <- fit$dictionary_support_table
    final$fit$dictionary_support_summary <- fit$dictionary_support_summary
    final$fit$candidate_tf_cor <- fit$candidate_tf_cor
    final$fit$candidate_peak_cor <- fit$candidate_peak_cor
    final$fit$inference_scope <-
        "no_fusion_condition_local_lm_exact_edge_omnibus_whole_network_BH"
    final$fit <- .condition_apply_activity_gate(final$fit)
    final$object <- .condition_update_network_significance(
        final$object, final$fit
    )
    if (isTRUE(verbose)) {
        message(
            "Pando condition phase=condition_fit_complete",
            " | cell_type=", as.character(progress_label),
            ";fit_edges=", nrow(dictionary),
            ";regcompass_edges=",
            sum(final$fit$edge_inference$edge_supported %in% TRUE),
            ";inference_edges=",
            sum(final$fit$edge_inference$edge_inference_estimable %in% TRUE),
            ";targets=", length(unique(as.character(dictionary$target)))
        )
    }
    final
}
