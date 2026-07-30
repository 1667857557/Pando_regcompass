test_that("condition API has no reference-condition argument or contrast export", {
    expect_false("reference_condition" %in%
        names(formals(Pando:::infer_condition_grn.GRNData)))
    expect_false("condition_grn_contrast" %in% getNamespaceExports("Pando"))
    expect_false(exists(
        "condition_grn_contrast", envir = asNamespace("Pando"),
        inherits = FALSE
    ))
})

test_that("stored condition contract retains absolute effects only", {
    fit <- structure(list(
        schema_version = "pando_condition_grn_fit_v5",
        contract_version = "legacy",
        beta_condition_std = matrix(
            c(1, -1), nrow = 1,
            dimnames = list("edge_1", c("A", "B"))
        ),
        reference_condition = "A",
        reference_beta = 1,
        contrast = matrix(c(0, -2), nrow = 1),
        comparison_mask = matrix(TRUE, nrow = 1, ncol = 2),
        response_transform = data.frame(
            target = "G1", center = 0, scale = 1,
            reference_condition = "A"
        ),
        direction_semantics = list(
            absolute = "sign(beta_condition)",
            pairwise = "legacy"
        )
    ), class = c("ConditionGRNFit", "list"))

    out <- Pando:::.condition_absolute_fit_contract(fit)
    expect_identical(out$beta_condition_std, fit$beta_condition_std)
    expect_identical(out$contract_version, "condition_absolute_oof_v2")
    expect_identical(
        out$coefficient_contract, "absolute_condition_effects_only"
    )
    expect_false(any(c(
        "reference_condition", "reference_beta", "contrast",
        "comparison_mask"
    ) %in% names(out)))
    expect_false("reference_condition" %in%
        colnames(out$response_transform))
    expect_null(out$direction_semantics$pairwise)
})
