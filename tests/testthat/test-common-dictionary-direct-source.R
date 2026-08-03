test_that("common-dictionary functions have one direct definition", {
  files <- list.files("../../R", pattern = "[.]R$", full.names = TRUE)
  text <- vapply(files, paste, collapse = "\n", FUN.VALUE = character(1))
  names(text) <- files
  pattern <- "(?m)^([A-Za-z.][A-Za-z0-9._]*)\\s*<-\\s*function\\b"
  definitions <- unlist(lapply(text, function(value) {
    matches <- gregexpr(pattern, value, perl = TRUE)
    raw <- regmatches(value, matches)[[1L]]
    if (!length(raw) || identical(raw, character(0))) return(character())
    sub("\\s*<-\\s*function.*$", "", raw)
  }), use.names = FALSE)
  expect_false(anyDuplicated(definitions) > 0L)
  expect_false(any(basename(files) %in% c(
    "zz_common_dictionary_validation.R",
    "zz_common_dictionary_projection.R",
    "zzz_fit_schema_contract.R"
  )))
})

test_that("fixed-dictionary model controls are explicit and strict", {
  expect_match(
    paste(deparse(formals(Pando::fit_grn_from_edges)$method), collapse = ""),
    "glm", fixed = TRUE
  )
  expect_identical(formals(Pando::fit_grn_from_edges)$interaction_term, ":")
  expect_identical(formals(Pando::fit_grn_from_edges)$scale, FALSE)
  expect_true(all(c("rna_layer", "peak_layer", "peak_value_type") %in%
                  names(formals(Pando::fit_grn_from_edges))))
})

test_that("exact edge union enforces one preprocessing reference", {
  candidate <- data.frame(
    target = "G", tf = "TF", region = "chr1-1-2",
    atac_feature_id = "chr1-1-2", peak_target_cor = 0.2,
    tf_target_cor = 0.3, source_label = "x",
    source_type = "condition", edge_id = "G||TF||chr1-1-2",
    stringsAsFactors = FALSE
  )
  class(candidate) <- c("PandoEdgeDictionary", "data.frame")
  stamp <- function(value, fingerprint) {
    attr(value, "rna_layer") <- "data"
    attr(value, "peak_layer") <- "data"
    attr(value, "peak_value_type") <- "normalized"
    attr(value, "preprocessing_fingerprint") <- fingerprint
    attr(value, "dictionary_input_schema") <-
      "pando_candidate_input_provenance_v1"
    value
  }
  global <- stamp(candidate, "same")
  condition <- stamp(candidate, "same")
  union <- Pando::union_grn_edges(
    global, list(control = condition)
  )
  expect_true(isTRUE(attr(
    union, "preprocessing_provenance_verified", exact = TRUE
  )))
  expect_identical(
    attr(union, "preprocessing_fingerprint", exact = TRUE), "same"
  )

  mismatched <- stamp(candidate, "different")
  expect_error(
    Pando::union_grn_edges(global, list(control = mismatched)),
    "different RNA/ATAC"
  )
  unverified <- candidate
  attributes(unverified)[c(
    "rna_layer", "peak_layer", "peak_value_type",
    "preprocessing_fingerprint", "dictionary_input_schema"
  )] <- NULL
  expect_error(
    Pando::union_grn_edges(global, list(control = unverified)),
    "mix verified and unverified"
  )
})

test_that("condition fit extraction has no compatibility arguments", {
  method <- getS3method("condition_grn_fit", "GRNData")
  expect_false("network_name" %in% names(formals(method)))
  expect_error(
    method(methods::new("GRNData"), network_name = "legacy"),
    "Unused condition_grn_fit argument"
  )
})

test_that("retired condition APIs are absent", {
  namespace <- asNamespace("Pando")
  retired <- c(
    "condition_grn_contrast",
    "project_condition_grn_primary_cells",
    "project_condition_grn_to_cells"
  )
  expect_false(any(vapply(retired, exists, logical(1), envir = namespace,
                          inherits = FALSE)))
})

