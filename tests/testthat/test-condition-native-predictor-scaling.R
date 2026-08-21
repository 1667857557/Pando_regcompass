test_that("native predictor and scaling match the R reference", {
    set.seed(177)
    cells <- paste0("cell", seq_len(17L))
    gene <- matrix(stats::rnorm(17L * 4L), 17L, 4L,
                   dimnames = list(cells, paste0("TF", 1:4)))
    peak <- matrix(stats::rnorm(17L * 3L), 17L, 3L,
                   dimnames = list(cells, paste0("peak", 1:3)))
    prepared <- list(
        gene_data = Matrix::Matrix(gene, sparse = TRUE),
        peak_data = Matrix::Matrix(peak, sparse = TRUE)
    )
    edges <- data.frame(
        tf = c("TF1", "TF2", "TF1", "TF4"),
        region = c("peak1", "peak2", "peak3", "peak1"),
        edge_id = paste0("edge", 1:4), stringsAsFactors = FALSE
    )
    condition_cells <- list(A = cells[1:7], B = cells[8:17])
    reference_x <- Pando:::.condition_ridge_predictors(
        prepared, edges, condition_cells
    )
    reference_scaling <- Pando:::.condition_ridge_scaling(reference_x, 1e-8)
    native <- Pando:::.condition_native_predictors_scaling(
        prepared, edges, condition_cells, 1e-8
    )

    expect_equal(native$x, reference_x, tolerance = 1e-12)
    expect_equal(native$scaling$center, reference_scaling$center,
                 tolerance = 1e-12)
    expect_equal(native$scaling$scale, reference_scaling$scale,
                 tolerance = 1e-12)
    expect_identical(native$scaling$informative,
                     reference_scaling$informative)
    expect_identical(native$scaling$reference,
                     reference_scaling$reference)
})
