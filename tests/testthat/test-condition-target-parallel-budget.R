test_that("condition target payload mapper uses supplied BiocParallel pool", {
  skip_if_not_installed("BiocParallel")
  old <- getOption("Pando.condition_target_BPPARAM", NULL)
  on.exit(options(Pando.condition_target_BPPARAM = old), add = TRUE)
  param <- BiocParallel::SerialParam()
  options(Pando.condition_target_BPPARAM = param)

  keys <- stats::setNames(c("a", "b"), c("a", "b"))
  observed <- .pando_target_payload_map(
    keys = keys,
    build_payload = function(key) list(skip = TRUE, target = key),
    worker_name = ".pando_discovery_target_worker",
    parallel = TRUE,
    verbose = FALSE,
    phase = "unit"
  )
  expect_identical(names(observed), names(keys))
  expect_identical(length(observed), 2L)
})

test_that("ordinary map_par is no longer overwritten by condition target routing", {
  body_text <- paste(deparse(body(map_par)), collapse = "\n")
  expect_false(grepl("Pando.condition_target_BPPARAM", body_text, fixed = TRUE))
  expect_false(exists(
    ".pando_condition_parallel_map_impl",
    envir = asNamespace("Pando"), inherits = FALSE
  ))
})

test_that("condition target parallel public controls remain stable", {
  args <- names(formals(infer_condition_grn.GRNData))
  expect_true(all(c("parallel", "BPPARAM", "parallel_scope") %in% args))
  expect_identical(formals(infer_condition_grn.GRNData)$parallel, FALSE)
})

test_that("condition public method is the single canonical definition", {
  expect_false(exists(
    ".pando_condition_parallel_method_impl",
    envir = asNamespace("Pando"), inherits = FALSE
  ))
  expect_false(exists(
    ".pando_condition_target_parallel_method",
    envir = asNamespace("Pando"), inherits = FALSE
  ))
  body_text <- paste(deparse(body(infer_condition_grn.GRNData)), collapse = "\n")
  expect_match(body_text, ".pando_infer_condition_grn_multitask_ridge_one", fixed = TRUE)
  expect_match(body_text, "Pando.condition_target_BPPARAM", fixed = TRUE)
})
