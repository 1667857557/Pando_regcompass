# Canonical unversioned schema for condition-aware GRN fits.

.PANDO_CONDITION_GRN_FIT_SCHEMA <- 'pando_condition_grn_fit'

if (exists('.condition_absolute_fit_contract', inherits = FALSE)) {
    if (!exists(
            '.condition_absolute_fit_contract_pre_schema',
            inherits = FALSE
        )) {
        .condition_absolute_fit_contract_pre_schema <-
            .condition_absolute_fit_contract
    }
    .condition_absolute_fit_contract <- function(fit) {
        fit <- .condition_absolute_fit_contract_pre_schema(fit)
        if (!is.list(fit)) return(fit)
        fit$schema_version <- .PANDO_CONDITION_GRN_FIT_SCHEMA
        fit$schema_policy <- 'single_unversioned_schema'
        fit$legacy_schema_aliases <- NULL
        fit
    }
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

.condition_replace_validator_symbol <- function(x) {
    if (is.symbol(x) && identical(as.character(x), '.condition_require_v5')) {
        return(as.name('.condition_require_fit'))
    }
    if (is.call(x) || is.pairlist(x) || is.expression(x)) {
        for (i in seq_along(x)) {
            x[[i]] <- .condition_replace_validator_symbol(x[[i]])
        }
    }
    x
}

.condition_rewritten_functions <- character()
for (.condition_function_name in c(
    '.project_condition_grn_cells_na',
    'project_condition_grn_cells',
    'condition_grn_subgraph'
)) {
    if (exists(.condition_function_name, inherits = FALSE)) {
        .condition_function <- get(
            .condition_function_name,
            inherits = FALSE
        )
        if (is.function(.condition_function)) {
            body(.condition_function) <- .condition_replace_validator_symbol(
                body(.condition_function)
            )
            assign(
                .condition_function_name,
                .condition_function,
                inherits = FALSE
            )
            .condition_rewritten_functions <- c(
                .condition_rewritten_functions,
                .condition_function_name
            )
        }
    }
}

# Retired version-specific and contrast helpers are absent from the usable API.
if (exists('.condition_require_v5', inherits = FALSE)) {
    rm(.condition_require_v5)
}
if (exists('condition_grn_contrast', inherits = FALSE)) {
    rm(condition_grn_contrast)
}
rm(.condition_function_name)
if (exists('.condition_function', inherits = FALSE)) {
    rm(.condition_function)
}
