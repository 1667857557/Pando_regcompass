test_that("standard infer dots remove condition-only controls", {
  routed <- .pando_sanitize_standard_infer_dots(
    list(
      padj_threshold = 0.05,
      rank_action = "mark",
      min_residual_df = 1L,
      maxit = 100L
    ),
    warn = FALSE
  )

  expect_setequal(
    routed$disabled,
    c("padj_threshold", "rank_action", "min_residual_df")
  )
  expect_equal(routed$args$maxit, 100L)
  expect_false(any(routed$disabled %in% names(routed$args)))
})

test_that("standard infer dots require names", {
  expect_error(
    .pando_sanitize_standard_infer_dots(list(0.05), warn = FALSE),
    "must be named"
  )
})

test_that("serial parallel helper preserves order", {
  observed <- .pando_parallel_lapply(
    1:4,
    function(x) x * 2L,
    parallel = TRUE,
    BPPARAM = FALSE
  )
  expect_identical(observed, as.list(c(2L, 4L, 6L, 8L)))
})

test_that("condition method exposes cell-type parallel controls", {
  args <- names(formals(infer_condition_grn.GRNData))
  expect_true(all(c("parallel", "BPPARAM", "parallel_scope") %in% args))
})

test_that("logical TRUE is not a BPPARAM", {
  expect_error(
    .pando_validate_bpparam(TRUE),
    "BiocParallelParam"
  )
})
