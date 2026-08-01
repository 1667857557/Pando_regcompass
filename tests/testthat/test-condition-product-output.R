test_that("canonical design preserves raw arithmetic and output fields", {
    cells <- paste0("c", seq_len(10L))
    gene <- Matrix::Matrix(cbind(
        TF1 = c(0, 1, 2, 0, 3, 1, 0, 2, 1, 4),
        TF2 = c(1, 0, 1, 2, 0, 3, 1, 0, 2, 1),
        G = seq_len(10L)
    ), sparse = TRUE, dimnames = list(cells, c("TF1", "TF2", "G")))
    peak <- Matrix::Matrix(cbind(
        P1 = c(1, 0, 1, 0, 2, 0, 1, 0, 2, 1),
        P2 = c(0, 1, 0, 2, 0, 1, 0, 2, 1, 0)
    ), sparse = TRUE, dimnames = list(cells, c("P1", "P2")))
    edges <- data.frame(
        tf = c("TF1", "TF2", "TF1"), target = "G",
        region = c("P1", "P2", "P2"), stringsAsFactors = FALSE
    )
    edges$edge_id <- paste(
        edges$tf, edges$region, edges$target, sep = "\001"
    )
    reference <- gene[, c(1L, 2L, 1L), drop = FALSE] *
        peak[, c(1L, 2L, 2L), drop = FALSE]
    variance <- Pando:::.condition_population_variance(reference)
    reference <- reference[, variance > .Machine$double.eps, drop = FALSE]
    result <- Pando:::.condition_build_design(
        gene[, "G", drop = FALSE], gene, peak, edges,
        factor(rep(c("A", "B"), each = 5L)), scale = FALSE
    )
    colnames(reference) <- result$edges$edge_id
    expect_identical(as.matrix(result$X_raw), as.matrix(reference))
    expect_identical(dimnames(result$X_raw), dimnames(reference))
    expect_identical(result$y_raw, as.numeric(gene[, "G"]))
    expect_identical(names(result), c("X_raw", "y_raw", "edges"))
})
