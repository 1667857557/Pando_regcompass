context('Seurat 4 compatibility tests')

library(testthat)
library(Pando)

testthat::test_that("Seurat v4 assay accessor is compatible", {
    testthat::skip_if_not_installed("Seurat")
    testthat::skip_if_not_installed("SeuratObject")

    counts <- Matrix::Matrix(
        matrix(
            c(
                1, 0, 2,
                3, 0, 4,
                0, 5, 1,
                2, 1, 0
            ),
            nrow = 3
        ),
        sparse = TRUE
    )

    rownames(counts) <- c("G1", "G2", "G3")
    colnames(counts) <- paste0("Cell", seq_len(ncol(counts)))

    object <- Seurat::CreateSeuratObject(counts = counts)
    object <- Seurat::NormalizeData(object, verbose = FALSE)
    object$pool <- c("A", "A", "B", "B")

    mat <- Pando:::.get_assay_data_compat(
        object,
        assay = "RNA",
        layer = "data"
    )

    testthat::expect_equal(dim(mat), dim(counts))
    testthat::expect_equal(rownames(mat), rownames(counts))
    testthat::expect_equal(colnames(mat), colnames(counts))

    object <- aggregate_assay(
        object,
        group_name = "pool",
        assay = "RNA",
        slot = "data"
    )

    summary_mat <- object[["RNA"]]@misc$summary$pool

    testthat::expect_equal(rownames(summary_mat), c("A", "B"))
    testthat::expect_equal(ncol(summary_mat), nrow(counts))
})
