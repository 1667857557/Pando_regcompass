test_that("projection membership requires every projected cell", {
  projection <- structure(list(
    gene_score = matrix(
      c(1, 3, 2, 4), nrow = 2, byrow = FALSE,
      dimnames = list(c("cell1", "cell2"), c("GENE1", "GENE2"))
    )
  ), class = c("PandoConditionProjection", "list"))

  incomplete <- data.frame(
    cell_id = "cell1", metacell_id = "mc1",
    stringsAsFactors = FALSE
  )
  expect_error(
    validate_condition_grn_projection_membership(projection, incomplete),
    "missing 1 projected paired cell"
  )

  membership <- data.frame(
    cell_id = c("other_cell", "cell2", "cell1"),
    metacell_id = c("other_mc", "mc1", "mc1"),
    stringsAsFactors = FALSE
  )
  selected <- validate_condition_grn_projection_membership(
    projection, membership
  )
  expect_identical(selected$cell_id, c("cell1", "cell2"))
  expect_identical(selected$metacell_id, c("mc1", "mc1"))

  aggregated <- aggregate_condition_grn_projection(projection, selected)
  expect_identical(colnames(aggregated$gene_score), "mc1")
  expect_equal(as.numeric(aggregated$gene_score[, "mc1"]), c(2, 3))
})

test_that("projection membership rejects duplicate membership rows", {
  projection <- structure(list(
    gene_score = matrix(
      1, nrow = 1, dimnames = list("cell1", "GENE1")
    )
  ), class = c("PandoConditionProjection", "list"))
  membership <- data.frame(
    cell_id = c("cell1", "cell1"),
    metacell_id = c("mc1", "mc1"),
    stringsAsFactors = FALSE
  )
  expect_error(
    validate_condition_grn_projection_membership(projection, membership),
    "one row per cell"
  )
})

test_that("projection membership rejects incomplete identifiers", {
  projection <- structure(list(
    gene_score = matrix(
      1, nrow = 1, dimnames = list("cell1", "GENE1")
    )
  ), class = c("PandoConditionProjection", "list"))
  membership <- data.frame(
    cell_id = "cell1", metacell_id = "",
    stringsAsFactors = FALSE
  )
  expect_error(
    validate_condition_grn_projection_membership(projection, membership),
    "complete"
  )
})
