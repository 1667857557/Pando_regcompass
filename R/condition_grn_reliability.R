#' Compute target-level reliability for condition GRN fits
#'
#' Converts the target-by-condition goodness-of-fit statistic stored by a
#' common-dictionary condition GRN into the amplitude-scale reliability used by
#' downstream regulatory projections. Reliability is defined as
#' `sqrt(pmin(1, pmax(0, rsq)))` for `fit_status == "ok"` target-condition fits
#' that retain at least one active edge. Targets without an active edge, with a
#' non-ok fit status, or with non-finite R-squared remain unavailable (`NA`).
#'
#' @param fit A `ConditionGRNFit` returned by [condition_grn_fit()].
#' @param significant_only If `TRUE`, an active target requires at least one
#'   coefficient with `significant == TRUE`. If `FALSE`, at least one estimable
#'   coefficient is sufficient.
#' @return A data frame with one row per target-condition fit and columns
#'   `target`, `condition`, `rsq`, `fit_status`, `n_active_edges`, and
#'   `reliability`.
#' @export
condition_grn_reliability <- function(fit, significant_only = TRUE) {
    if (!inherits(fit, "ConditionGRNFit") ||
        !identical(fit$schema_version, "pando_condition_grn_common_dictionary_v1")) {
        stop("`fit` must be a common-dictionary ConditionGRNFit.", call. = FALSE)
    }
    if (!is.logical(significant_only) || length(significant_only) != 1L ||
        is.na(significant_only)) {
        stop("`significant_only` must be TRUE or FALSE.", call. = FALSE)
    }

    fit_table <- as.data.frame(fit$fit, stringsAsFactors = FALSE)
    coefficient <- as.data.frame(fit$coefficients, stringsAsFactors = FALSE)
    required_fit <- c("target", "condition", "rsq", "fit_status")
    required_coefficient <- c("target", "condition", "estimable", "significant")
    if (!all(required_fit %in% colnames(fit_table)) || !nrow(fit_table) ||
        !all(required_coefficient %in% colnames(coefficient))) {
        stop("Condition GRN reliability requires aligned fit and coefficient tables.",
             call. = FALSE)
    }

    fit_target <- toupper(trimws(as.character(fit_table$target)))
    fit_condition <- as.character(fit_table$condition)
    fit_key <- paste(fit_target, fit_condition, sep = "\001")
    if (anyNA(fit_key) || any(!nzchar(fit_target)) ||
        any(!nzchar(fit_condition)) || anyDuplicated(fit_key)) {
        stop("Target-condition fit diagnostics must be complete and unique.",
             call. = FALSE)
    }

    coefficient_target <- toupper(trimws(as.character(coefficient$target)))
    coefficient_condition <- as.character(coefficient$condition)
    coefficient_key <- paste(
        coefficient_target, coefficient_condition, sep = "\001"
    )
    if (anyNA(coefficient_key) || any(!nzchar(coefficient_target)) ||
        any(!nzchar(coefficient_condition))) {
        stop("Condition coefficients contain incomplete target-condition labels.",
             call. = FALSE)
    }
    coefficient_index <- match(coefficient_key, fit_key)
    if (anyNA(coefficient_index)) {
        stop("Condition coefficients cannot be aligned to target fit diagnostics.",
             call. = FALSE)
    }

    active_edge <- if (isTRUE(significant_only)) {
        coefficient$significant %in% TRUE
    } else {
        coefficient$estimable %in% TRUE
    }
    n_active_edges <- tabulate(
        coefficient_index[active_edge], nbins = nrow(fit_table)
    )

    rsq <- suppressWarnings(as.numeric(fit_table$rsq))
    fit_status <- trimws(as.character(fit_table$fit_status))
    if (anyNA(fit_status) || any(!nzchar(fit_status))) {
        stop("Condition fit_status values must be complete.", call. = FALSE)
    }
    reliability <- rep(NA_real_, nrow(fit_table))
    eligible <- fit_status == "ok" & n_active_edges > 0L & is.finite(rsq)
    reliability[eligible] <- sqrt(pmin(1, pmax(0, rsq[eligible])))

    data.frame(
        target = as.character(fit_table$target),
        condition = fit_condition,
        rsq = rsq,
        fit_status = fit_status,
        n_active_edges = as.integer(n_active_edges),
        reliability = reliability,
        stringsAsFactors = FALSE
    )
}
