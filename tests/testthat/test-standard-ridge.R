test_that("standard ridge is the K=1 condition ridge penalty", {
  p <- 4L
  expect_equal(
    .condition_ridge_penalty(1L, p, fusion_ratio = 100),
    diag(p)
  )
})

test_that("standard infer_grn keeps glm default and exposes ridge controls", {
  args <- formals(infer_grn.GRNData)
  expect_true(all(c(
    "BPPARAM", "ridge_control", "rank_action", "min_residual_df",
    "padj_threshold"
  ) %in% names(args)))
  methods <- eval(args$method)
  expect_identical(methods[[1L]], "glm")
  expect_true("ridge" %in% methods)
})

test_that("standard ridge accepts Gaussian identity only", {
  expect_true(.pando_standard_ridge_family_ok("gaussian"))
  expect_true(.pando_standard_ridge_family_ok(stats::gaussian()))
  expect_false(.pando_standard_ridge_family_ok(stats::poisson()))
})
