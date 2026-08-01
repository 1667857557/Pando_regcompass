test_that("checkpoint names and diagnostic upgrades are deterministic", {
    path <- Pando:::.condition_checkpoint_path(
        tempdir(), "TF/target", "0123456789abcdef-extra"
    )
    expect_identical(
        basename(path), "TF_target__0123456789abcdef.rds"
    )

    old <- data.frame(
        target = "G1", stage = "complete", stringsAsFactors = FALSE
    )
    new <- data.frame(
        target = "G2", stage = "complete",
        path_backend = "sparse_matrix_free",
        stringsAsFactors = FALSE
    )
    combined <- Pando:::.condition_bind_schema_rows(old, new)
    expect_identical(
        names(combined), c("target", "stage", "path_backend")
    )
    expect_true(is.na(combined$path_backend[[1L]]))
    expect_identical(combined$path_backend[[2L]], "sparse_matrix_free")
})

test_that("engine controls reject duplicate field names", {
    duplicate <- structure(list(128, 256), names = c(
        "memory_budget_mb", "memory_budget_mb"
    ))
    expect_error(
        Pando:::.condition_normalize_engine_control(duplicate),
        "unique non-empty"
    )
})
