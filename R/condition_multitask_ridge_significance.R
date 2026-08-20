# Significance and RegCompass handoff contract for conditional E-star/JSE GRNs.
#
# Candidate correlation support defines the frozen exact-edge dictionary.
# Production coefficients remain continuous on that dictionary. Joint-refit
# Wald p-values are adjusted within condition x target families. An exact edge
# is admitted to RegCompass iff at least one condition has padj < threshold and
# every fitted condition has a valid continuous production coefficient.

.condition_significant_projection_policy <-
    "any_condition_padj_exact_edge_union"
.condition_fit_dictionary_policy <-
    "global_and_condition_union_pando_correlation_supported_frozen_dictionary"

.condition_validate_adjust_method <- function(adjust_method) {
    value <- toupper(as.character(adjust_method))
    if (length(value) != 1L || is.na(value) || !identical(value, "BH")) {
        stop(
            "Conditional E-star inference requires `adjust_method = \"BH\"`.",
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

.condition_apply_activity_gate <- function(fit) {
    if (!inherits(fit, "ConditionGRNFit")) {
        stop("A ConditionGRNFit is required for condition activity annotation.",
             call. = FALSE)
    }
    threshold <- .condition_validate_padj_threshold(fit$padj_threshold)
    coefficient <- as.data.frame(fit$coefficients, stringsAsFactors = FALSE)
    support <- fit$dictionary_support_table
    fit_table <- as.data.frame(fit$fit, stringsAsFactors = FALSE)
    required_coefficient <- c(
        "edge_id", "target", "condition", "estimate", "penalty_effect",
        "pval", "padj", "inference_estimable", "condition_significant",
        "fusion_component_id", "shared_edge"
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
        stop("Conditional E-star/JSE metadata are incomplete.",
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
    padj <- suppressWarnings(as.numeric(coefficient$padj))
    expected_condition_significant <-
        coefficient$inference_estimable %in% TRUE &
        is.finite(padj) & padj < threshold
    if (!identical(
        as.logical(coefficient$condition_significant),
        expected_condition_significant
    )) {
        stop("Condition significance does not match condition-target BH.",
             call. = FALSE)
    }
    coefficient$statistically_supported <- expected_condition_significant
    coefficient$significant <- expected_condition_significant
    coefficient$pando_estimation_active <- is.finite(estimate)
    coefficient$active <- coefficient$pando_estimation_active
    coefficient$penalty_effect <- ifelse(
        coefficient$pando_estimation_active, estimate, NA_real_
    )

    coefficient$edge_union_supported <- FALSE
    coefficient$supporting_conditions <- ""
    coefficient$n_supporting_conditions <- 0L
    coefficient$all_conditions_fit_valid <- FALSE
    coefficient$active_in_regcompass <- FALSE
    edge_ids <- unique(as.character(coefficient$edge_id))
    for (edge_id in edge_ids) {
        index <- which(as.character(coefficient$edge_id) == edge_id)
        if (length(index) != length(fit$condition_levels) ||
            anyDuplicated(as.character(coefficient$condition[index]))) {
            stop("Every exact edge must occur once in every fitted condition.",
                 call. = FALSE)
        }
        valid_edge <- all(
            as.character(coefficient$fit_status[index]) == "ok" &
            is.finite(as.numeric(coefficient$penalty_effect[index]))
        )
        supporting <- as.character(coefficient$condition[index][
            coefficient$condition_significant[index] %in% TRUE
        ])
        supporting <- fit$condition_levels[
            fit$condition_levels %in% supporting
        ]
        union_support <- valid_edge && length(supporting) > 0L
        coefficient$all_conditions_fit_valid[index] <- valid_edge
        coefficient$supporting_conditions[index] <-
            paste(supporting, collapse = ";")
        coefficient$n_supporting_conditions[index] <- length(supporting)
        coefficient$edge_union_supported[index] <- union_support
        coefficient$active_in_regcompass[index] <- union_support
    }

    fit$coefficients <- coefficient
    fit$projection_effect_column <- "penalty_effect"
    fit$projection_policy <- .condition_significant_projection_policy
    fit$dictionary_support_role <-
        "global_plus_condition_candidate_provenance"
    fit$local_support_role <-
        "condition_specific_candidate_provenance_only"
    fit$global_support_role <-
        "pooled_candidate_provenance_only"
    fit$statistical_support_role <-
        "joint_refit_Wald_BH_within_condition_target"
    fit$activity_role <-
        "Pando_keeps_continuous_estimates;RegCompass_uses_exact_edge_any_condition_BH_union"
    fit$regcompass_edge_gate <-
        "all_conditions_fit_valid AND any_condition_padj_lt_threshold"
    fit
}

.condition_ridge_fit_contract <- function(
    object, fit, prepared, control, rank_action = "mark",
    min_residual_df = 1L, parallel = FALSE, verbose = TRUE) {
    fit$adjust_method <- .condition_validate_adjust_method(fit$adjust_method)
    fit$padj_threshold <- .condition_validate_padj_threshold(fit$padj_threshold)
    control <- .condition_E_star_control(control)
    dictionary <- fit$edge_dictionary
    progress_label <- fit$cell_type %||% ""

    if (isTRUE(verbose)) {
        message(
            "Pando condition phase=E_star_z025_JSE_start",
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
        progress_phase = "E_star_z025_JSE",
        progress_label = progress_label
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
        "E_star_z025_primary_fusion_component_joint_refit_condition_target_BH"
    final$fit <- .condition_apply_activity_gate(final$fit)
    final$object <- .condition_update_network_significance(
        final$object, final$fit
    )
    if (isTRUE(verbose)) {
        message(
            "Pando condition phase=condition_fit_complete",
            " | cell_type=", as.character(progress_label),
            ";fit_edges=", nrow(dictionary),
            ";regcompass_union_edges=",
            length(unique(final$fit$coefficients$edge_id[
                final$fit$coefficients$edge_union_supported %in% TRUE
            ])),
            ";condition_significant_rows=",
            sum(final$fit$coefficients$condition_significant %in% TRUE),
            ";targets=", length(unique(as.character(dictionary$target)))
        )
    }
    final
}
