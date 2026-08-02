test_that("native runtime preflight is lightweight", {
    description <- utils::packageDescription("Pando")
    expect_identical(
        description[["Config/Pando/RuntimePreflight"]],
        "lightweight-metadata-symbol-check-v1"
    )

    preflight <- Pando:::.condition_native_self_test_cpp()
    expect_true(preflight$passed)
    expect_identical(
        preflight$mode,
        "lightweight_metadata_symbol_check_v1"
    )
    expect_false(preflight$numerical_test_run)
    expect_identical(
        preflight$memory_contract,
        "no-full-p2-on-high-p-path-v1"
    )
    expect_identical(length(preflight$native_symbols), 4L)
    expect_null(preflight$numerical)
})
