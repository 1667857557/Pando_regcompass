make_test_grn_design <- function() {
  edges <- data.frame(
    edge_id = "TF1::peak1::GENE1",
    candidate_index = 1L,
    tf = "TF1",
    region = "chr1-1-10",
    target = "GENE1",
    atac_feature_id = "peak1",
    tf_feature_id = "TF1",
    target_feature_id = "GENE1",
    motif_supported = TRUE,
    peak_to_gene_supported = TRUE,
    supporting_regions = "chr1-1-10;chr1-2-9",
    n_supporting_regions = 2L,
    tf_detection = 0.5,
    peak_detection = 0.4,
    target_detection = 0.6,
    stringsAsFactors = FALSE
  )
  design <- structure(list(
    schema_version = "pando_grn_design_v1",
    candidate_edges = edges,
    region_map = data.frame(
      region = c("chr1-1-10", "chr1-2-9"),
      atac_feature_id = c("peak1", "peak1"),
      stringsAsFactors = FALSE
    ),
    target_diagnostics = data.frame(
      target = "GENE1",
      target_detection = 0.6,
      n_candidate_edges = 1L,
      n_candidate_tfs = 1L,
      n_candidate_regions = 1L,
      stringsAsFactors = FALSE
    ),
    feature_contract = list(
      cell_ids = c("c1", "c2"),
      rna_feature_ids = c("TF1", "GENE1"),
      atac_feature_ids = "peak1",
      rna_assay = "RNA",
      atac_assay = "ATAC",
      rna_layer = "data",
      atac_layer = "data"
    ),
    params = list(candidate_policy = "structural_shared_before_model_fitting"),
    design_fingerprint = "legacy"
  ), class = c("PandoGRNDesign", "list"))
  .pando_refresh_grn_design_contract(design)
}

test_that("version-2 shared GRN design validates its exact feature contract", {
  design <- make_test_grn_design()

  expect_identical(design$schema_version, "pando_grn_design_v2")
  expect_match(design$design_fingerprint, "^md5:")
  expect_true(validate_grn_design(design))
  expect_identical(design$target_diagnostics$n_candidate_regions, 2L)
  expect_identical(design$target_diagnostics$n_candidate_atac_features, 1L)
})

test_that("design validation rejects edge, feature, mapping and fingerprint drift", {
  design <- make_test_grn_design()

  broken <- design
  broken$candidate_edges$edge_id <- "wrong"
  expect_error(validate_grn_design(broken), "uniquely encode")

  broken <- design
  broken$candidate_edges$atac_feature_id <- "missing_peak"
  expect_error(validate_grn_design(broken), "feature contract")

  broken <- design
  broken$region_map$atac_feature_id[[1L]] <- "other_peak"
  expect_error(validate_grn_design(broken), "Representative regulatory regions")

  broken <- design
  broken$params$candidate_policy <- "changed"
  expect_error(validate_grn_design(broken), "fingerprint")
})

test_that("prepare_grn_design is a public S3 generic", {
  expect_true(is.function(prepare_grn_design))
  expect_true(is.function(prepare_grn_design.GRNData))
  expect_true(is.function(validate_grn_design))
})
