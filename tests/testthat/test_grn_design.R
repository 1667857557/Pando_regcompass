test_that("shared GRN design validation enforces the edge contract", {
  design <- structure(list(
    schema_version = "pando_grn_design_v1",
    candidate_edges = data.frame(
      edge_id = "TF1::chr1-1-10::GENE1",
      tf = "TF1",
      region = "chr1-1-10",
      target = "GENE1",
      atac_feature_id = "chr1-1-10",
      tf_feature_id = "TF1",
      target_feature_id = "GENE1",
      stringsAsFactors = FALSE
    ),
    feature_contract = list(
      cell_ids = c("c1", "c2"),
      rna_feature_ids = c("TF1", "GENE1"),
      atac_feature_ids = "chr1-1-10",
      rna_assay = "RNA",
      atac_assay = "ATAC"
    )
  ), class = c("PandoGRNDesign", "list"))

  expect_true(validate_grn_design(design))
  design$candidate_edges$edge_id <- ""
  expect_error(validate_grn_design(design), "unique and non-empty")
})

test_that("prepare_grn_design is a public S3 generic", {
  expect_true(is.function(prepare_grn_design))
  expect_true(is.function(prepare_grn_design.GRNData))
})
