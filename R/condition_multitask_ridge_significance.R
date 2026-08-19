# Activity and diagnostic-inference contract for common-dictionary Scheme E.
#
# Candidate correlation support determines admission to the frozen exact-edge
# dictionary and remains provenance only. Scheme-E coefficients are continuous
# on that fixed dictionary. Condition-wise BH statistics are retained as model
# diagnostics; they never change network membership and never overwrite a fitted
# coefficient with zero.

.condition_significant_projection_policy <-
    "continuous_common_dictionary_scheme_e_effects"
.condition_fit_dictionary_policy <-
    "global_and_condition_union_pando_correlation_supported_frozen_dictionary"

.condition_validate_adjust_method <- function(adjust_method) {
    value <- toupper(as.character(adjust_method))
    if (length(value) != 1L || is.na(value) || !identical(value, "BH")) {
        stop(
            "Condition diagnostic inference requires `adjust_method = \"BH\"`.",
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
    required_coefficient <- c(
        "edge_id", "condition", "estimate", "padj",
        "contrast_identifiable", "shared_by_boundary", "fused_by_penalty"
    )
    required_support <- c(
        "edge_id", "source_type", "condition",
        "peak_target_cor", "tf_target_cor"
    )
    if (!all(required_coefficient %in% colnames(coefficient)) ||
        !is.data.frame(support) ||
        !all(required_support %in% colnames(support))) {
        stop("Scheme-E condition-GRN support metadata are incomplete.",
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
    coefficient$dictionary_support <- TRUE
    coefficient$current_scope_correlation_support <-
        local_support | global_support
    coefficient$peak_cor_pass <- local_support
    coefficient$tf_cor_pass <- local_support

    estimate <- suppressWarnings(as.numeric(coefficient$estimate))
    padj <- suppressWarnings(as.numeric(coefficient$padj))
    statistically_supported <- is.finite(estimate) & is.finite(padj) &
        padj < threshold
    # Fixed topology: every finite coefficient on the frozen dictionary remains
    # active for projection. BH is diagnostic only.
    active <- is.finite(estimate)
    coefficient$statistically_supported <- statistically_supported
    coefficient$active <- active
    coefficient$significant <- statistically_supported
    coefficient$penalty_effect <- ifelse(active, estimate, NA_real_)
    fit$coefficients <- coefficient
    fit$projection_effect_column <- "penalty_effect"
    fit$projection_policy <- .condition_significant_projection_policy
    fit$dictionary_support_role <-
        "membership_in_global_plus_condition_correlation_screened_common_dictionary"
    fit$local_support_role <-
        "condition_specific_pando_correlation_provenance_only"
    fit$global_support_role <-
        "pooled_all_eligible_conditions_pando_correlation_provenance_only"
    fit$statistical_support_role <-
        "condition_wise_BH_adjusted_wald_diagnostic_only"
    fit$activity_role <-
        "all_finite_scheme_e_coefficients_on_fixed_common_dictionary"
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
            "Pando condition phase=scheme_e_z025_start",
            " | cell_type=", as.character(progress_label),
            ";fit_edges=", nrow(dictionary),
            ";targets=", length(unique(as.character(dictionary$target)))
        )
    }
    final <- .condition_ridge_fit_contract_one_pass(
        object = object, fit = fit, prepared = prepared, control = control,
        rank_action = rank_action, min_residual_df = min_residual_df,
        parallel = parallel, verbose = verbose,
        progress_phase = "scheme_e_z025",
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
        "scheme_e_z025_primary;BH_and_target_R2_are_diagnostics_only"
    final$fit <- .condition_apply_activity_gate(final$fit)
    final$object <- .condition_update_network_significance(
        final$object, final$fit
    )
    if (isTRUE(verbose)) {
        message(
            "Pando condition phase=condition_fit_complete",
            " | cell_type=", as.character(progress_label),
            ";fit_edges=", nrow(dictionary),
            ";projection_edges=", sum(final$fit$coefficients$active %in% TRUE),
            ";bh_supported_diagnostic=",
            sum(final$fit$coefficients$statistically_supported %in% TRUE),
            ";targets=", length(unique(as.character(dictionary$target)))
        )
    }
    final
}
