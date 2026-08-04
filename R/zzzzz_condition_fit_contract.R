# Complete condition-fit provenance and strict downstream aggregation.

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

#' Aggregate a condition projection without dropping fitted cells
#'
#' Requires every projected paired cell to have exactly one complete membership
#' row. Extra membership rows from other cell types are allowed.
#'
#' @param projection A `PandoConditionProjection`.
#' @param membership Data frame containing `cell_id` and `group_col`.
#' @param group_col Membership grouping column.
#' @return Aggregated target-by-group regulatory scores.
#' @export
aggregate_condition_grn_projection_strict <- function(
    projection, membership, group_col = "metacell_id") {
    if (!inherits(projection, "PandoConditionProjection") ||
        !is.matrix(projection$gene_score) ||
        is.null(rownames(projection$gene_score)) ||
        anyNA(rownames(projection$gene_score)) ||
        any(!nzchar(rownames(projection$gene_score))) ||
        anyDuplicated(rownames(projection$gene_score)) ||
        !is.data.frame(membership) ||
        !all(c("cell_id", group_col) %in% colnames(membership))) {
        stop("Projection and membership inputs are invalid.", call. = FALSE)
    }
    cell_id <- as.character(membership$cell_id)
    group_id <- as.character(membership[[group_col]])
    if (anyNA(cell_id) || any(!nzchar(trimws(cell_id))) ||
        anyDuplicated(cell_id) || anyNA(group_id) ||
        any(!nzchar(trimws(group_id)))) {
        stop(
            "Membership cell and group IDs must be complete, with one row per cell.",
            call. = FALSE
        )
    }
    projected_cells <- rownames(projection$gene_score)
    missing <- setdiff(projected_cells, cell_id)
    if (length(missing)) {
        stop(
            "Membership is missing ", length(missing),
            " projected paired cell(s); first missing ID: ", missing[[1L]],
            call. = FALSE
        )
    }
    selected <- membership[match(projected_cells, cell_id), , drop = FALSE]
    groups <- unique(as.character(selected[[group_col]]))
    if (!ncol(projection$gene_score)) {
        return(list(
            gene_score = matrix(
                numeric(), nrow = 0L, ncol = length(groups),
                dimnames = list(character(), groups)
            ),
            group_col = group_col,
            source_projection = projection,
            aggregation =
                "arithmetic_mean_of_all_projected_paired_cell_regulatory_scores"
        ))
    }
    answer <- aggregate_condition_grn_projection(
        projection = projection,
        membership = selected,
        group_col = group_col
    )
    answer$aggregation <-
        "arithmetic_mean_of_all_projected_paired_cell_regulatory_scores"
    answer
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
