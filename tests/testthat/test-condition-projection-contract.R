test_that("projection targets resolve across gene-symbol case", {
  expect_identical(
    .condition_resolve_projection_targets(
      c("tp53", "slc22a17"), c("TP53", "Slc22a17")
    ),
    c("TP53", "Slc22a17")
  )
})

test_that("case-ambiguous fitted targets are rejected", {
  expect_error(
    .condition_resolve_projection_targets("gene", c("GENE", "gene")),
    "ambiguous after case normalization"
  )
})

test_that("fitted cells are complete and condition-disjoint", {
  fit <- list(
    condition_levels = c("control", "treated"),
    condition_cell_ids = list(
      control = c("cell1", "cell2"),
      treated = c("cell3", "cell4")
    )
  )
  contract <- .condition_validate_projection_cells(
    fit, paste0("cell", 1:4)
  )
  expect_identical(contract$cells, paste0("cell", 1:4))

  overlapping <- fit
  overlapping$condition_cell_ids$treated[[1L]] <- "cell2"
  expect_error(
    .condition_validate_projection_cells(overlapping, paste0("cell", 1:4)),
    "more than one condition"
  )

  expect_error(
    .condition_validate_projection_cells(fit, paste0("cell", 1:3)),
    "missing 1 fitted paired cell"
  )
})
