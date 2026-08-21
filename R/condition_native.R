# Native numerical kernels for the conditional common-dictionary path.

.condition_native_predictors_scaling <- function(
    prepared, edges, cells_by_condition, floor) {
    tfs <- unique(as.character(edges$tf))
    regions <- unique(as.character(edges$region))
    tf_index <- match(as.character(edges$tf), tfs)
    peak_index <- match(as.character(edges$region), regions)
    if (anyNA(tf_index) || anyNA(peak_index)) {
        stop("Native predictor indices are incomplete.", call. = FALSE)
    }
    gene <- lapply(cells_by_condition, function(cells) {
        as.matrix(prepared$gene_data[cells, tfs, drop = FALSE])
    })
    peak <- lapply(cells_by_condition, function(cells) {
        as.matrix(prepared$peak_data[cells, regions, drop = FALSE])
    })
    native <- condition_cpp_predictors_scaling(
        gene, peak, as.integer(tf_index - 1L), as.integer(peak_index - 1L),
        as.numeric(floor)
    )
    x <- native$x
    names(x) <- names(cells_by_condition)
    x <- lapply(x, function(one) {
        colnames(one) <- as.character(edges$edge_id)
        one
    })
    columns <- as.character(edges$edge_id)
    center <- stats::setNames(as.numeric(native$center), columns)
    scale <- stats::setNames(as.numeric(native$scale), columns)
    informative <- stats::setNames(as.logical(native$informative), columns)
    list(
        x = x,
        scaling = list(
            center = center, scale = scale, informative = informative,
            floor = floor,
            reference = "equal_condition_within_condition_rms"
        )
    )
}
