# Complete condition-fit provenance at the exported extraction boundary.

.pando_complete_condition_fit_contract <- function(fit) {
    if (!inherits(fit, "ConditionGRNFit")) return(fit)
    field <- "dictionary_preprocessing_provenance_verified"
    if (!field %in% names(fit)) {
        fit[[field]] <- isTRUE(attr(
            fit$edge_dictionary,
            "preprocessing_provenance_verified",
            exact = TRUE
        ))
    }
    fit
}

.pando_complete_condition_fit_contracts <- function(fits) {
    if (inherits(fits, "ConditionGRNFit")) {
        return(.pando_complete_condition_fit_contract(fits))
    }
    if (!is.list(fits)) return(fits)
    lapply(fits, .pando_complete_condition_fit_contract)
}

.pando_condition_grn_fit_complete_contract_method <- function(
    object, cell_type = NULL, ...) {
    answer <- .pando_condition_grn_fit_base_impl(
        object = object, cell_type = cell_type, ...
    )
    .pando_complete_condition_fit_contracts(answer)
}

.onLoad <- function(libname, pkgname) {
    namespace <- asNamespace(pkgname)
    base_method <- get(
        "condition_grn_fit.GRNData", envir = namespace, inherits = FALSE
    )
    assign(
        ".pando_condition_grn_fit_base_impl", base_method,
        envir = namespace
    )
    registerS3method(
        "condition_grn_fit", "GRNData",
        get(
            ".pando_condition_grn_fit_complete_contract_method",
            envir = namespace, inherits = FALSE
        ),
        envir = namespace
    )
}
