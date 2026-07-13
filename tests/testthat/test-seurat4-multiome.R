context('Seurat 4 multiome compatibility tests')

library(testthat)
library(Pando)

testthat::test_that("GRNData accepts a Seurat v4 multiome object", {
    testthat::skip_if_not_installed("Seurat")
    testthat::skip_if_not_installed("Signac")
    testthat::skip_if_not_installed("GenomicRanges")
    testthat::skip_if_not_installed("IRanges")

    rna_counts <- Matrix::Matrix(
        matrix(
            c(
                2, 0, 1, 3,
                0, 2, 4, 1,
                1, 1, 0, 2
            ),
            nrow = 3
        ),
        sparse = TRUE
    )

    rownames(rna_counts) <- c("G1", "G2", "TF1")
    colnames(rna_counts) <- paste0("Cell", seq_len(ncol(rna_counts)))

    peak_counts <- Matrix::Matrix(
        matrix(
            c(
                1, 0, 2, 1,
                0, 1, 1, 3
            ),
            nrow = 2
        ),
        sparse = TRUE
    )

    rownames(peak_counts) <- c(
        "chr1-100-200",
        "chr1-500-600"
    )
    colnames(peak_counts) <- colnames(rna_counts)

    annotation <- GenomicRanges::GRanges(
        seqnames = c("chr1", "chr1", "chr1"),
        ranges = IRanges::IRanges(
            start = c(50, 450, 800),
            end = c(300, 700, 1000)
        ),
        strand = c("+", "+", "+"),
        type = c("gene", "gene", "gene"),
        gene_name = c("G1", "G2", "TF1"),
        gene_id = c("G1", "G2", "TF1")
    )

    object <- Seurat::CreateSeuratObject(counts = rna_counts)
    object <- Seurat::NormalizeData(object, verbose = FALSE)
    object <- Seurat::FindVariableFeatures(
        object,
        nfeatures = 3,
        verbose = FALSE
    )

    object[["peaks"]] <- Signac::CreateChromatinAssay(
        counts = peak_counts,
        sep = c("-", "-"),
        annotation = annotation
    )

    grn_object <- initiate_grn(
        object,
        peak_assay = "peaks",
        rna_assay = "RNA",
        exclude_exons = FALSE
    )

    testthat::expect_s4_class(grn_object, "GRNData")
    testthat::expect_s4_class(grn_object@data, "Seurat")

    rna_data <- LayerData(
        grn_object,
        assay = "RNA",
        layer = "data"
    )

    testthat::expect_equal(dim(rna_data), dim(rna_counts))
    testthat::expect_equal(
        Params(grn_object)$peak_assay,
        "peaks"
    )
    testthat::expect_equal(
        Params(grn_object)$rna_assay,
        "RNA"
    )
})
