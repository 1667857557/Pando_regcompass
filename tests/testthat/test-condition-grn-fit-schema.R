test_that("condition GRN fits use one canonical unversioned schema", {
    fit <- structure(list(
        schema_version = "pando_condition_grn_fit_v4",
        beta_condition_std = matrix(
            c(1, -1), nrow = 1,
            dimnames = list("edge_1", c("A", "B"))
        ),
        response_transform = data.frame(
            target = "G1", center = 0, scale = 1
        )
    ), class = c("ConditionGRNFit", "list"))

    canonical <- Pando:::.condition_absolute_fit_contract(fit)
    expect_identical(
        canonical$schema_version,
        "pando_condition_grn_fit"
    )
    expect_identical(canonical$schema_policy, "single_unversioned_schema")
    expect_silent(Pando:::.condition_require_fit(canonical))
})

test_that("version-suffixed fit schemas are not usable", {
    for (schema in c(
        "pando_condition_grn_fit_v4",
        "pando_condition_grn_fit_v5"
    )) {
        fit <- structure(
            list(schema_version = schema),
            class = c("ConditionGRNFit", "list")
        )
        expect_error(
            Pando:::.condition_require_fit(fit),
            "version-suffixed schemas are not supported",
            fixed = TRUE
        )
    }
})

test_that("only the stable extractor is exported", {
    exports <- getNamespaceExports("Pando")
    expect_true("condition_grn_fit" %in% exports)
    expect_false(any(grepl(
        "^pando_condition_grn_fit_v[0-9]+$|^condition_grn_fit_v[0-9]+$",
        exports
    )))
})
