test_that("condition-estimable OOF projection is supplemental", {
    body_text <- paste(
        deparse(body(project_condition_grn_supplemental_cells)), collapse = "\n"
    )
    expect_match(body_text, "support_policy = 'condition_estimable'", fixed = TRUE)
    expect_match(body_text, "diagnostic_only = TRUE", fixed = TRUE)
    expect_match(body_text, "projection_used_for_penalty <- TRUE", fixed = TRUE)
    expect_match(body_text, "projection_role <- 'supplemental_penalty'", fixed = TRUE)
    expect_match(body_text, "common_support_primary <- TRUE", fixed = TRUE)
})

test_that("supplemental projection is exported and documented", {
    expect_true(is.function(project_condition_grn_supplemental_cells))
    root <- Sys.getenv("GITHUB_WORKSPACE", unset = "")
    candidates <- unique(c(
        if (nzchar(root)) root else character(), ".", "..", file.path("..", "..")
    ))
    candidates <- candidates[file.exists(file.path(candidates, "NAMESPACE"))]
    if (!length(candidates)) skip("Package source is unavailable.")
    root <- normalizePath(candidates[[1L]], mustWork = TRUE)
    namespace <- paste(
        readLines(file.path(root, "NAMESPACE"), warn = FALSE), collapse = "\n"
    )
    expect_match(
        namespace,
        "export(project_condition_grn_supplemental_cells)",
        fixed = TRUE
    )
    expect_true(file.exists(file.path(
        root, "man", "project_condition_grn_supplemental_cells.Rd"
    )))
})
