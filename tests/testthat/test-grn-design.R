test_that("candidate edge builder preserves TF-peak-target triplets", {
    peaks2gene <- Matrix::Matrix(
        matrix(c(1, 0, 1, 1), nrow = 2, byrow = TRUE), sparse = TRUE,
        dimnames = list(c("G1", "G2"), c("R1", "R2"))
    )
    peaks2motif <- Matrix::Matrix(
        matrix(c(1, 1, 0, 1), nrow = 2, byrow = TRUE), sparse = TRUE,
        dimnames = list(c("R1", "R2"), c("M1", "M2"))
    )
    motif2tf <- Matrix::Matrix(
        matrix(c(1, 0, 0, 1), nrow = 2, byrow = TRUE), sparse = TRUE,
        dimnames = list(c("M1", "M2"), c("TF1", "TF2"))
    )
    edges <- Pando:::.pando_candidate_edge_table(
        peaks2gene = peaks2gene,
        peaks2motif = peaks2motif,
        motif2tf = motif2tf,
        region_to_peak = c(R1 = "chr1-1-10", R2 = "chr1-20-30")
    )

    expect_equal(
        edges$edge_id,
        c("TF1::R1::G1", "TF2::R1::G1", "TF2::R2::G2")
    )
    expect_equal(
        edges$atac_feature_id,
        c("chr1-1-10", "chr1-1-10", "chr1-20-30")
    )
    expect_false(anyDuplicated(edges$edge_id))
})

test_that("GRN design validation enforces the matrix contract", {
    edges <- data.frame(
        edge_id = "TF1::R1::G1",
        tf = "TF1",
        region = "R1",
        target = "G1",
        tf_feature_id = "TF1",
        atac_feature_id = "P1",
        target_feature_id = "G1",
        stringsAsFactors = FALSE
    )
    design <- structure(list(
        schema_version = "pando_grn_design_v1",
        candidate_edges = edges,
        region_map = data.frame(region = "R1", atac_feature_id = "P1"),
        target_diagnostics = data.frame(target = "G1"),
        feature_contract = list(
            cell_ids = c("C1", "C2"),
            rna_feature_ids = c("TF1", "G1"),
            atac_feature_ids = "P1"
        ),
        params = list(),
        design_id = "pando_grn_design_v1_test"
    ), class = c("PandoGRNDesign", "list"))

    expect_invisible(validate_grn_design(design))
    design$candidate_edges$atac_feature_id <- "missing"
    expect_error(validate_grn_design(design), "matrix contract")
})
