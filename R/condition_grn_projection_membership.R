# Complete-cell membership validation for condition-GRN projection aggregation.

#' Validate condition-GRN projection membership
#'
#' Requires every projected paired cell to have exactly one complete membership
#' row. Membership rows for cells outside the projection are allowed and are not
#' returned.
#'
#' @param projection A `PandoConditionProjection`.
#' @param membership Data frame containing `cell_id` and `group_col`.
#' @param group_col Membership grouping column, typically `metacell_id`.
#' @return `membership` restricted to all projected cells and ordered exactly as
#'   `rownames(projection$gene_score)`.
#' @export
validate_condition_grn_projection_membership <- function(
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
    if (!identical(as.character(selected$cell_id), projected_cells)) {
        stop("Projection membership could not be aligned exactly to projected cells.",
             call. = FALSE)
    }
    selected
}
