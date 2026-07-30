# Canonical unversioned schema for condition-aware GRN fits.

.PANDO_CONDITION_GRN_FIT_SCHEMA <- 'pando_condition_grn_fit'

.condition_absolute_fit_contract_canonical <- .condition_absolute_fit_contract
.condition_absolute_fit_contract <- function(fit) {
    fit <- .condition_absolute_fit_contract_canonical(fit)
    if (!is.list(fit)) return(fit)
    fit$schema_version <- .PANDO_CONDITION_GRN_FIT_SCHEMA
    fit$schema_policy <- 'single_unversioned_schema'
    fit$legacy_schema_aliases <- NULL
    fit
}

.condition_require_fit <- function(fit) {
    if (!inherits(fit, 'ConditionGRNFit') ||
        !identical(fit$schema_version, .PANDO_CONDITION_GRN_FIT_SCHEMA)) {
        stop(
            'fit must be a pando_condition_grn_fit ConditionGRNFit object; ',
            'version-suffixed schemas are not supported.'
        )
    }
    invisible(TRUE)
}

# Internal compatibility for original source paths that still call the old
# validator symbol. It validates only the canonical unversioned schema and does
# not accept pando_condition_grn_fit_v4, pando_condition_grn_fit_v5, or aliases.
.condition_require_v5 <- function(fit) {
    .condition_require_fit(fit)
}
