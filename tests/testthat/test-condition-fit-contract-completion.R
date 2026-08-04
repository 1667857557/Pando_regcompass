test_that("missing dictionary provenance is materialized on condition fits", {
  dictionary <- data.frame(edge_id = "G||TF||chr1:1-2")
  attr(dictionary, "preprocessing_provenance_verified") <- TRUE
  fit <- structure(
    list(edge_dictionary = dictionary),
    class = c("ConditionGRNFit", "list")
  )

  completed <- .pando_complete_condition_fit_contract(fit)

  expect_true(completed$dictionary_preprocessing_provenance_verified)
})

test_that("explicit provenance values are not overwritten", {
  dictionary <- data.frame(edge_id = "G||TF||chr1:1-2")
  attr(dictionary, "preprocessing_provenance_verified") <- TRUE
  fit <- structure(
    list(
      edge_dictionary = dictionary,
      dictionary_preprocessing_provenance_verified = FALSE
    ),
    class = c("ConditionGRNFit", "list")
  )

  completed <- .pando_complete_condition_fit_contract(fit)

  expect_false(completed$dictionary_preprocessing_provenance_verified)
})

test_that("named condition fit lists retain their names", {
  dictionary <- data.frame(edge_id = "G||TF||chr1:1-2")
  attr(dictionary, "preprocessing_provenance_verified") <- TRUE
  fit <- structure(
    list(edge_dictionary = dictionary),
    class = c("ConditionGRNFit", "list")
  )

  completed <- .pando_complete_condition_fit_contracts(
    list(epithelial_like = fit)
  )

  expect_identical(names(completed), "epithelial_like")
  expect_true(
    completed$epithelial_like$dictionary_preprocessing_provenance_verified
  )
})
