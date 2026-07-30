test_that("condition levels select the intended analysis mode", {
  meta <- data.frame(
    condition = c("A", "A", "B"),
    cell_type = c("T", "T", "T"),
    row.names = paste0("cell", 1:3)
  )
  expect_identical(.condition_resolve_levels(meta, NULL), character())
  expect_identical(.condition_resolve_levels(meta, "missing"), character())
  expect_identical(.condition_resolve_levels(meta[1:2, , drop = FALSE], "condition"), "A")
  expect_setequal(.condition_resolve_levels(meta, "condition"), c("A", "B"))
})

test_that("public condition inference exposes explicit fallback arguments", {
  arguments <- names(formals(infer_condition_grn.GRNData))
  expect_true(all(c("condition_col", "cell_type_col", "fallback_args") %in% arguments))
  expect_null(formals(infer_condition_grn.GRNData)$condition_col)
  expect_null(formals(infer_condition_grn.GRNData)$cell_type_col)
})

test_that("canonical functions are direct definitions", {
  expect_true(exists(".condition_require_fit", inherits = FALSE))
  expect_false(exists(".condition_require_v5", inherits = FALSE))
  expect_false(exists(".infer_condition_grn_with_reference", inherits = FALSE))
  expect_false(exists(".project_condition_grn_cells_na", inherits = FALSE))
  expect_false(exists("condition_grn_contrast", inherits = FALSE))
})
