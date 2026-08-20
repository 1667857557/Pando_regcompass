test_that("standard ridge has an independent CV ridge solver", {
  control <- Pando:::.pando_standard_ridge_control()
  expect_true(all(control$lambda_grid > 0))
  expect_identical(control$lambda_rule, "1se")
  expect_false("scheme_e_z" %in% names(control))
  expect_true(is.function(Pando:::.pando_standard_ridge_numeric))
  expect_true(is.function(Pando:::.pando_standard_ridge_cv))
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

test_that("standard ridge routing remains separate from Scheme E", {
  body_text <- paste(deparse(body(infer_grn.GRNData)), collapse = "\n")
  expect_match(body_text, "identical(method, \"ridge\")", fixed = TRUE)
  expect_match(body_text, ".pando_standard_ridge_fit", fixed = TRUE)
  standard_text <- paste(
    deparse(body(Pando:::.pando_standard_ridge_fit)), collapse = "\n"
  )
  expect_match(standard_text, ".pando_standard_ridge_control", fixed = TRUE)
  expect_match(standard_text, ".pando_standard_ridge_target_worker", fixed = TRUE)
  expect_false(grepl("condition_scheme_e", standard_text, fixed = TRUE))
})