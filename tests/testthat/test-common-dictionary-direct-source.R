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
