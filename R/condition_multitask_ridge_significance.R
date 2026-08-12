# Activity and inference contract for global-plus-condition common-dictionary ridge.
#
# The frozen dictionary is the exact union of TF-peak-target edges passing the
# Pando peak-target and TF-target correlation gates either in the pooled/global
# cell set or in at least one condition. Ridge is fitted once. An edge is active
# in condition c only when beta_c is BH-supported and the edge has either global
# Pando support or condition-c local Pando support. Global support therefore
# protects candidate admission against small-condition correlation instability,
# but never makes a condition edge active without condition-specific beta evidence.

.condition_significant_projection_policy <-
    "active_global_or_local_pando_support_and_condition_bh_ridge_effects"
.condition_fit_dictionary_policy <-
    "global_and_condition_union_pando_correlation_supported_frozen_dictionary"

.condition_validate_adjust_method <- function(adjust_method) {
    value <- toupper(as.character(adjust_method))
    if (length(value) != 1L || is.na(value) || !identical(value, "BH")) {
        stop(
            "Condition ridge inference requires `adjust_method = \"BH\"`.",
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
        stop("A ConditionGRNFit is required for condition activity gating.",
             call. = FALSE)
    }
    threshold <- .condition_validate_padj_threshold(fit$padj_threshold)
    coefficient <- as.data.frame(fit$coefficients, stringsAsFactors = FALSE)
    support <- fit$dictionary_support_table
    required_coefficient <- c(
        "edge_id", "condition", "estimate", "estimable", "padj"
    )
    required_support <- c(
        "edge_id", "source_type", "condition",
        "peak_target_cor", "tf_target_cor"
    )
    if (!all(required_coefficient %in% colnames(coefficient)) ||
        !is.data.frame(support) ||
        !all(required_support %in% colnames(support))) {
        stop("Global/condition Pando support metadata are incomplete.",
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
    global_index <- match(as.character(coefficient$edge_id),
                          as.character(global$edge_id))
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
    coefficient$dictionary_support <- local_support | global_support
    coefficient$peak_cor_pass <- local_support
    coefficient$tf_cor_pass <- local_support

    estimate <- suppressWarnings(as.numeric(coefficient$estimate))
    padj <- suppressWarnings(as.numeric(coefficient$padj))
    statistically_supported <- coefficient$estimable %in% TRUE &
        is.finite(estimate) & is.finite(padj) & padj < threshold
    active <- statistically_supported & coefficient$dictionary_support
    coefficient$statistically_supported <- statistically_supported
    coefficient$active <- active
    # Compatibility alias for existing Network consumers. The explicit
    # statistically_supported column is the condition-wise ridge-BH result.
    coefficient$significant <- active
    coefficient$penalty_effect <- ifelse(active, estimate, 0)
    fit$coefficients <- coefficient
    fit$projection_effect_column <- "penalty_effect"
    fit$projection_policy <- .condition_significant_projection_policy
    fit$dictionary_support_role <-
        "pooled_global_or_condition_local_pando_correlation_candidate_support"
    fit$local_support_role <-
        "condition_specific_pando_correlation_support_annotation"
    fit$global_support_role <-
        "pooled_all_eligible_conditions_pando_correlation_support_annotation"
    fit$statistical_support_role <-
        "condition_wise_BH_adjusted_approximate_ridge_wald"
    fit
}

.condition_ridge_fit_contract <- function(
    object, fit, prepared, control, rank_action = "mark",
    min_residual_df = 1L, parallel = FALSE, verbose = TRUE) {
    fit$adjust_method <- .condition_validate_adjust_method(fit$adjust_method)
    fit$padj_threshold <- .condition_validate_padj_threshold(fit$padj_threshold)
    control <- .condition_ridge_control(control)
    dictionary <- fit$edge_dictionary
    progress_label <- fit$cell_type %||% ""

    if (isTRUE(verbose)) {
        message(
            "Pando condition phase=single_no_fusion_ridge_start",
            " | cell_type=", as.character(progress_label),
            ";fit_edges=", nrow(dictionary),
            ";targets=", length(unique(as.character(dictionary$target)))
        )
    }
    final <- .condition_ridge_fit_contract_one_pass(
        object = object, fit = fit, prepared = prepared, control = control,
        rank_action = rank_action, min_residual_df = min_residual_df,
        parallel = parallel, verbose = verbose,
        progress_phase = "ridge_single",
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
        "approximate_ridge_wald_conditional_on_global_or_condition_pando_screened_dictionary_and_cv_lambda"
    final$fit <- .condition_apply_activity_gate(final$fit)
    final$object <- .condition_update_network_significance(
        final$object, final$fit
    )
    if (isTRUE(verbose)) {
        message(
            "Pando condition phase=condition_fit_complete",
            " | cell_type=", as.character(progress_label),
            ";fit_edges=", nrow(dictionary),
            ";active_edges=", sum(final$fit$coefficients$active %in% TRUE),
            ";targets=", length(unique(as.character(dictionary$target)))
        )
    }
    final
}
