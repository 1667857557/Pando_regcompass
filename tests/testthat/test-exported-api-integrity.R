test_that("every exported symbol is defined in the Pando namespace", {
  ns <- asNamespace("Pando")
  exports <- getNamespaceExports("Pando")
  missing <- exports[!vapply(exports, function(name) {
    exists(name, envir = ns, inherits = FALSE)
  }, logical(1))]
  expect_identical(missing, character())
})

test_that("key GRN S3 implementations remain defined after override removal", {
  ns <- asNamespace("Pando")
  methods <- c(
    "infer_grn.GRNData",
    "infer_condition_grn.GRNData",
    "condition_grn_fit.GRNData",
    "find_modules.Network",
    "find_modules.GRNData"
  )
  missing <- methods[!vapply(methods, function(name) {
    exists(name, envir = ns, inherits = FALSE) &&
      is.function(get(name, envir = ns, inherits = FALSE))
  }, logical(1))]
  expect_identical(missing, character())
})

test_that("condition projection public API remains available", {
  expect_true(is.function(project_condition_grn_cells))
  expect_true(is.function(aggregate_condition_grn_projection))
  expect_true(is.function(validate_condition_grn_projection_membership))
})
