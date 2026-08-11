# Significance-screened fit dictionary for multi-task condition GRNs.
#
# Candidate discovery deliberately remains broad: global and condition-specific
# Pando candidates are exact-unioned so a condition-specific edge cannot be lost
# before statistical estimation. A preliminary joint multi-task ridge is then
# used as the statistical screen. The actual fit dictionary is the union of
# edges with BH-adjusted ridge-Wald P below the configured threshold in at least
# one condition. The joint ridge is finally refit on that shared screened
# dictionary. Thus every condition is estimated on the same statistically
# supported edge set, while the downstream condition-specific projection still
# requires significance in that condition.

.condition_max_padj_threshold <- 0.1
.condition_significant_projection_policy <- "padj_significant_ridge_effects"
.condition_fit_dictionary_policy <-
    "preliminary_joint_ridge_bh_significant_union_then_joint_refit"

.condition_validate_adjust_method <- function(adjust_method) {
    value <- toupper(as.character(adjust_method))
    if (length(value) != 1L || is.na(value) || !identical(value, "BH")) {
        stop(
            "Multi-condition ridge significance screening requires ",
            "`adjust_method = \"BH\"`.", call. = FALSE
        )
    }
    "BH"
}

.condition_validate_padj_threshold <- function(padj_threshold) {
    if (!is.numeric(padj_threshold) || length(padj_threshold) != 1L ||
        !is.finite(padj_threshold) || padj_threshold <= 0 ||
        padj_threshold > .condition_max_padj_threshold) {
        stop(
            "`padj_threshold` must be one finite number in (0, 0.1].",
            call. = FALSE
        )
    }
    as.numeric(padj_threshold)
}

.condition_apply_significance_gate <- function(fit) {
    if (!inherits(fit, "ConditionGRNFit")) {
        stop("A ConditionGRNFit is required for significance gating.",
             call. = FALSE)
    }
    threshold <- .condition_validate_padj_threshold(fit$padj_threshold)
    coefficient <- as.data.frame(fit$coefficients, stringsAsFactors = FALSE)
    required <- c("estimate", "estimable", "padj")
    if (!all(required %in% colnames(coefficient))) {
        stop(
            "Condition ridge coefficients require estimate, estimable and padj ",
            "before significance gating.", call. = FALSE
        )
    }
    estimate <- suppressWarnings(as.numeric(coefficient$estimate))
    padj <- suppressWarnings(as.numeric(coefficient$padj))
    significant <- coefficient$estimable %in% TRUE &
        is.finite(estimate) & is.finite(padj) & padj < threshold
    coefficient$significant <- significant
    coefficient$penalty_effect <- ifelse(significant, estimate, 0)
    fit$coefficients <- coefficient
    fit$projection_effect_column <- "penalty_effect"
    fit$projection_policy <- .condition_significant_projection_policy
    fit
}

.condition_dictionary_screen <- function(fit) {
    .condition_validate_adjust_method(fit$adjust_method)
    threshold <- .condition_validate_padj_threshold(fit$padj_threshold)
    coefficient <- as.data.frame(fit$coefficients, stringsAsFactors = FALSE)
    required <- c("edge_id", "condition", "estimate", "estimable", "pval", "padj")
    if (!all(required %in% colnames(coefficient))) {
        stop("Preliminary ridge coefficients are incomplete for dictionary screening.",
             call. = FALSE)
    }
    estimate <- suppressWarnings(as.numeric(coefficient$estimate))
    padj <- suppressWarnings(as.numeric(coefficient$padj))
    supported <- coefficient$estimable %in% TRUE &
        is.finite(estimate) & is.finite(padj) & padj < threshold

    ids <- unique(as.character(coefficient$edge_id))
    summary <- lapply(ids, function(id) {
        index <- which(as.character(coefficient$edge_id) == id)
        sig_index <- index[supported[index]]
        finite_padj <- padj[index][is.finite(padj[index])]
        data.frame(
            edge_id = id,
            screening_min_padj = if (length(finite_padj)) min(finite_padj) else NA_real_,
            screening_n_significant_conditions = length(sig_index),
            screening_significant_conditions = if (length(sig_index)) {
                paste(as.character(coefficient$condition[sig_index]), collapse = ";")
            } else "",
            stringsAsFactors = FALSE
        )
    })
    summary <- do.call(rbind, summary)
    rownames(summary) <- NULL
    list(
        keep_edge_ids = summary$edge_id[
            summary$screening_n_significant_conditions > 0L
        ],
        summary = summary,
        threshold = threshold
    )
}

.condition_subset_dictionary <- function(dictionary, keep_edge_ids, screen_summary) {
    keep <- as.character(dictionary$edge_id) %in% as.character(keep_edge_ids)
    out <- dictionary[keep, , drop = FALSE]
    if (!nrow(out)) return(out)

    screen_index <- match(as.character(out$edge_id), screen_summary$edge_id)
    out$screening_min_padj <- screen_summary$screening_min_padj[screen_index]
    out$screening_n_significant_conditions <-
        screen_summary$screening_n_significant_conditions[screen_index]
    out$screening_significant_conditions <-
        screen_summary$screening_significant_conditions[screen_index]

    source_attributes <- attributes(dictionary)
    structural <- c("names", "row.names", "class")
    for (name in setdiff(names(source_attributes), structural)) {
        attr(out, name) <- source_attributes[[name]]
    }
    class(out) <- class(dictionary)
    out
}

# Preserve the numerical one-pass estimator from
# condition_multitask_ridge_refit.R. The canonical wrapper runs a preliminary
# screen and, when required, a second fit on the significant-union dictionary.
.condition_ridge_refit_contract_one_pass <- .condition_ridge_refit_contract

.condition_ridge_refit_contract <- function(
    object, fit, prepared, control, rank_action = "mark",
    min_residual_df = 1L, parallel = FALSE, verbose = TRUE) {
    fit$adjust_method <- .condition_validate_adjust_method(fit$adjust_method)
    fit$padj_threshold <- .condition_validate_padj_threshold(fit$padj_threshold)
    candidate_dictionary <- fit$edge_dictionary

    preliminary <- .condition_ridge_refit_contract_one_pass(
        object = object, fit = fit, prepared = prepared, control = control,
        rank_action = rank_action, min_residual_df = min_residual_df,
        parallel = parallel, verbose = verbose
    )
    screen <- .condition_dictionary_screen(preliminary$fit)
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
    } else {
        final_skeleton <- fit
        final_skeleton$edge_dictionary <- final_dictionary
        final_skeleton$target_genes <- unique(as.character(final_dictionary$target))
        final_skeleton$coefficients <- NULL
        final_skeleton$contrasts <- NULL
        final_skeleton$fit <- NULL
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
    final
}
