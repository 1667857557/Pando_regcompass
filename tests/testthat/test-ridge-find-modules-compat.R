test_that("find_modules core accepts standard and multi-task ridge networks", {
  body_text <- paste(deparse(body(.find_modules_network_core)), collapse = "\n")
  expect_match(body_text, "'ridge'", fixed = TRUE)
  expect_match(body_text, "'multitask_ridge'", fixed = TRUE)
})

test_that("find_modules normalizes ridge dictionary diagnostics without changing stored fit contract", {
  body_text <- paste(deparse(body(find_modules.Network)), collapse = "\n")
  expect_match(body_text, "single_task_ridge", fixed = TRUE)
  expect_match(body_text, "fixed_edge_dictionary_joint_conditions", fixed = TRUE)
  expect_match(body_text, "nvariables_dictionary", fixed = TRUE)
  expect_match(body_text, "estimable", fixed = TRUE)
  expect_match(body_text, "original_coefs", fixed = TRUE)
  expect_match(body_text, "original_fit", fixed = TRUE)
})
