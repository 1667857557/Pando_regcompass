# Screening-supported fit dictionary for multi-task condition GRNs.
#
# Candidate discovery remains broad: global and condition-specific Pando
# candidates are exact-unioned before statistical estimation. A preliminary
# joint multi-task ridge supplies the only edge-level ridge-Wald P values and BH
# adjustment. Edges supported in at least one condition define the final shared
# dictionary. That dictionary is then refit jointly with inference disabled.
# Final coefficients are therefore comparable quantitative effects on one common
# model layer; the preliminary BH result remains the condition-specific support
# mask and is not recomputed after selection.

.condition_significant_projection_policy <-
    "screen_bh_supported_refit_ridge_effects"
.condition_fit_dictionary_policy <-
    "preliminary_joint_ridge_bh_supported_union_then_effect_only_joint_refit"

.condition_validate_adjust_method <- function(adjust_method) {
    value <- toupper(as.character(adjust_method))
    if (length(value) != 1L || is.na(value) || !identical(value, "BH")) {
        stop(
            "Multi-condition ridge dictionary screening requires ",
            "`adjust_method = \"BH\"`.", call. = FALSE
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

.condition_dictionary_screen <- function(fit) {
    .condition_validate_adjust_method(fit$adjust_method)
    threshold <- .condition_validate_padj_threshold(fit$padj_threshold)
    coefficient <- as.data.frame(fit$coefficients, stringsAsFactors = FALSE)
    required <- c(
        "edge_id", "condition", "estimate", "std_err", "statistic",
        "estimable", "pval", "padj"
    )
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

    audit <- data.frame(
        edge_id = as.character(coefficient$edge_id),
        condition = as.character(coefficient$condition),
        screen_estimate = estimate,
        screen_std_err = suppressWarnings(as.numeric(coefficient$std_err)),
        screen_statistic = suppressWarnings(as.numeric(coefficient$statistic)),
        screen_pval = suppressWarnings(as.numeric(coefficient$pval)),
        screen_padj = padj,
        screen_estimable = coefficient$estimable %in% TRUE,
        screen_significant = supported,
        stringsAsFactors = FALSE
    )
    if (anyDuplicated(paste(audit$edge_id, audit$condition, sep = "\001"))) {
        stop("Preliminary ridge screening produced duplicated edge-condition rows.",
             call. = FALSE)
    }

    list(
        keep_edge_ids = summary$edge_id[
            summary$screening_n_significant_conditions > 0L
        ],
        summary = summary,
        coefficients = audit,
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

.condition_apply_screen_support <- function(fit, screen_coefficients) {
    if (!inherits(fit, "ConditionGRNFit")) {
        stop("A ConditionGRNFit is required for screening-support projection.",
             call. = FALSE)
    }
    coefficient <- as.data.frame(fit$coefficients, stringsAsFactors = FALSE)
    required_final <- c(
        "edge_id", "condition", "estimate", "estimable", "pval", "padj"
    )
    required_screen <- c(
        "edge_id", "condition", "screen_estimate", "screen_std_err",
        "screen_statistic", "screen_pval", "screen_padj",
        "screen_estimable", "screen_significant"
    )
    if (!all(required_final %in% colnames(coefficient)) ||
        !is.data.frame(screen_coefficients) ||
        !all(required_screen %in% colnames(screen_coefficients))) {
        stop("Condition ridge screening and final-refit coefficients are incomplete.",
             call. = FALSE)
    }
    if (any(is.finite(suppressWarnings(as.numeric(coefficient$pval)))) ||
        any(is.finite(suppressWarnings(as.numeric(coefficient$padj))))) {
        stop(
            "Final screened-dictionary ridge refit must not contain second-round ",
            "P values or adjusted P values.", call. = FALSE
        )
    }

    final_key <- paste(
        as.character(coefficient$edge_id),
        as.character(coefficient$condition), sep = "\001"
    )
    screen_key <- paste(
        as.character(screen_coefficients$edge_id),
        as.character(screen_coefficients$condition), sep = "\001"
    )
    if (anyDuplicated(final_key) || anyDuplicated(screen_key)) {
        stop("Edge-condition keys must be unique across screening and refit.",
             call. = FALSE)
    }
    index <- match(final_key, screen_key)
    if (anyNA(index)) {
        stop("A final common-dictionary edge lacks preliminary screening evidence.",
             call. = FALSE)
    }

    fields <- setdiff(required_screen, c("edge_id", "condition"))
    for (field in fields) {
        coefficient[[field]] <- screen_coefficients[[field]][index]
    }
    estimate <- suppressWarnings(as.numeric(coefficient$estimate))
    supported <- coefficient$screen_significant %in% TRUE
    active <- supported & coefficient$estimable %in% TRUE & is.finite(estimate)
    coefficient$significant <- active
    coefficient$penalty_effect <- ifelse(active, estimate, 0)
    coefficient$inference_scope <-
        "post_screen_joint_ridge_effect_estimation_only_no_final_inference"

    fit$coefficients <- coefficient
    fit$projection_effect_column <- "penalty_effect"
    fit$projection_policy <- .condition_significant_projection_policy
    fit$inference_scope <-
        "post_screen_joint_ridge_effect_estimation_only_no_final_inference"
    fit$inference_performed <- FALSE
    fit
}

.condition_ridge_refit_contract <- function(
    object, fit, prepared, control, rank_action = "mark",
    min_residual_df = 1L, parallel = FALSE, verbose = TRUE) {
    fit$adjust_method <- .condition_validate_adjust_method(fit$adjust_method)
    fit$padj_threshold <- .condition_validate_padj_threshold(fit$padj_threshold)
    candidate_dictionary <- fit$edge_dictionary
    progress_label <- fit$cell_type %||% ""

    preliminary <- .condition_ridge_refit_contract_one_pass(
        object = object, fit = fit, prepared = prepared, control = control,
        rank_action = rank_action, min_residual_df = min_residual_df,
        parallel = parallel, verbose = verbose,
        progress_phase = "ridge_preliminary_screen",
        progress_label = progress_label,
        inference = TRUE
    )
    screen <- .condition_dictionary_screen(preliminary$fit)
    if (isTRUE(verbose)) {
        message(
            "Pando condition phase=bh_dictionary_screen",
            " | cell_type=", as.character(progress_label),
            ";candidate_edges=", nrow(candidate_dictionary),
            ";supported_edges=", length(screen$keep_edge_ids),
            ";threshold=", format(screen$threshold, trim = TRUE)
        )
    }
    if (!length(screen$keep_edge_ids)) {
        stop(
            "No condition-GRN edge passes preliminary BH padj < ",
            format(screen$threshold, trim = TRUE),
            " in any condition; no statistically supported common fit ",
            "dictionary can be constructed.", call. = FALSE
        )
    }

    final_dictionary <- .condition_subset_dictionary(
        candidate_dictionary, screen$keep_edge_ids, screen$summary
    )
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
        parallel = parallel, verbose = verbose,
        progress_phase = "ridge_final_effect_refit",
        progress_label = progress_label,
        inference = FALSE
    )

    final$fit <- .condition_apply_screen_support(
        final$fit, screen$coefficients
    )
    final$fit$candidate_edge_count <- nrow(candidate_dictionary)
    final$fit$fit_dictionary_edge_count <- nrow(final_dictionary)
    final$fit$fit_dictionary_policy <- .condition_fit_dictionary_policy
    final$fit$dictionary_screening_threshold <- screen$threshold
    final$fit$dictionary_screening_summary <- screen$summary
    final$fit$screening_inference_scope <- preliminary$fit$inference_scope
    final$fit$screening_adjust_method <- preliminary$fit$adjust_method
    final$fit$screening_padj_threshold <- screen$threshold
    final$fit$edge_dictionary <- final_dictionary
    final$object <- .condition_update_network_significance(final$object, final$fit)
    if (isTRUE(verbose)) {
        message(
            "Pando condition phase=condition_fit_complete",
            " | cell_type=", as.character(progress_label),
            ";candidate_edges=", nrow(candidate_dictionary),
            ";fit_edges=", nrow(final_dictionary),
            ";targets=", length(unique(as.character(final_dictionary$target))),
            ";final_inference=disabled"
        )
    }
    preliminary <- NULL
    invisible(gc(verbose = FALSE, full = TRUE))
    final
}
