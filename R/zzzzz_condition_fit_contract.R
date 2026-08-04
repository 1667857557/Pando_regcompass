# Complete condition-fit provenance on fitted objects and extraction.

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

.pando_complete_condition_fit_object <- function(object) {
    if (!inherits(object, "GRNData")) return(object)
    fits <- object@grn@params$condition_grn_fits
    if (is.list(fits) && length(fits)) {
        object@grn@params$condition_grn_fits <-
            .pando_complete_condition_fit_contracts(fits)
    }
    object
}

.pando_infer_condition_grn_complete_contract_method <- function(
    object, ...) {
    answer <- .pando_infer_condition_grn_base_impl(object = object, ...)
    .pando_complete_condition_fit_object(answer)
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
    infer_base <- get(
        "infer_condition_grn.GRNData", envir = namespace, inherits = FALSE
    )
    fit_base <- get(
        "condition_grn_fit.GRNData", envir = namespace, inherits = FALSE
    )
    assign(
        ".pando_infer_condition_grn_base_impl", infer_base,
        envir = namespace
    )
    assign(
        ".pando_condition_grn_fit_base_impl", fit_base,
        envir = namespace
    )
    registerS3method(
        "infer_condition_grn", "GRNData",
        get(
            ".pando_infer_condition_grn_complete_contract_method",
            envir = namespace, inherits = FALSE
        ),
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
