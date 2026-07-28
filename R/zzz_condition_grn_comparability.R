# Comparison support for condition-aware GRN contrasts.

.condition_reference_comparison_mask <- function(
    eligibility_mask, reference_condition) {
    eligibility_mask <- as.matrix(eligibility_mask)
    if (!is.logical(eligibility_mask) || anyNA(eligibility_mask) ||
        is.null(rownames(eligibility_mask)) ||
        is.null(colnames(eligibility_mask))) {
        stop('eligibility_mask must be a named logical matrix without NA values.')
    }
    if (!is.character(reference_condition) ||
        length(reference_condition) != 1L ||
        !reference_condition %in% colnames(eligibility_mask)) {
        stop('reference_condition must identify one eligibility-mask column.')
    }
    reference_eligible <- eligibility_mask[, reference_condition]
    comparison_mask <- eligibility_mask & matrix(
        reference_eligible,
        nrow = nrow(eligibility_mask),
        ncol = ncol(eligibility_mask)
    )
    dimnames(comparison_mask) <- dimnames(eligibility_mask)
    comparison_mask
}

.condition_combine_fit_contracts_without_comparison_mask <-
    .condition_combine_fit_contracts

.condition_combine_fit_contracts <- function(...) {
    fit <- .condition_combine_fit_contracts_without_comparison_mask(...)
    fit$comparison_mask <- .condition_reference_comparison_mask(
        fit$eligibility_mask, fit$reference_condition
    )
    fit$comparison_mask_formula <-
        'eligible_in_condition AND eligible_in_reference'
    fit
}
