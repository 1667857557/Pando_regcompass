test_that("condition target map uses supplied BiocParallel pool", {
  skip_if_not_installed("BiocParallel")
  old <- getOption("Pando.condition_target_BPPARAM", NULL)
  on.exit(options(Pando.condition_target_BPPARAM = old), add = TRUE)
  param <- BiocParallel::SerialParam()
  options(Pando.condition_target_BPPARAM = param)
  observed <- map_par(1:4, function(x) x * 3L, parallel = TRUE, verbose = FALSE)
  expect_identical(observed, as.list(c(3L, 6L, 9L, 12L)))
})

test_that("condition target parallel public controls remain stable", {
  args <- names(formals(infer_condition_grn.GRNData))
  expect_true(all(c("parallel", "BPPARAM", "parallel_scope") %in% args))
  expect_identical(formals(infer_condition_grn.GRNData)$parallel, FALSE)
})
