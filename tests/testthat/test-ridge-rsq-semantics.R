test_that("ridge fit exposes final full-data R2 as rsq", {
    body_text <- paste(
        deparse(body(Pando:::.condition_ridge_target)), collapse = "\n"
    )
    expect_match(body_text, "rsq = fit$rsq[[i]]", fixed = TRUE)
    expect_match(body_text, "rsq_oof = cv$rsq_oof[[i]]", fixed = TRUE)
    expect_match(body_text, "rsq_in_sample = fit$rsq[[i]]", fixed = TRUE)
    expect_false(grepl("rsq = cv$rsq_oof[[i]]", body_text, fixed = TRUE))
})
