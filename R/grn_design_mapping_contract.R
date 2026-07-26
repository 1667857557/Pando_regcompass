# Supporting-region mapping validation loaded after grn_design_contract.R.

.pando_validate_grn_design_core <- validate_grn_design

validate_grn_design <- function(design) {
    .pando_validate_grn_design_core(design)
    edges <- design$candidate_edges
    if (!is.data.frame(edges) || !nrow(edges)) return(invisible(TRUE))

    region_map <- design$region_map
    lookup <- stats::setNames(
        as.character(region_map$atac_feature_id),
        as.character(region_map$region)
    )
    support <- .pando_split_supporting_regions(edges$supporting_regions)
    consistent <- vapply(seq_len(nrow(edges)), function(i) {
        mapped <- unname(lookup[support[[i]]])
        length(mapped) == length(support[[i]]) &&
            !anyNA(mapped) &&
            all(mapped == as.character(edges$atac_feature_id[[i]]))
    }, logical(1))
    if (!all(consistent)) {
        stop(
            paste(
                'Every supporting regulatory region must map to the exact',
                'measured ATAC feature recorded for its candidate edge.'
            ),
            call. = FALSE
        )
    }
    invisible(TRUE)
}
