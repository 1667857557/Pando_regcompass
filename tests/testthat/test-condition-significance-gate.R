test_that("condition ridge screening policy requires BH and any valid padj threshold", {
    expect_identical(Pando:::.condition_validate_adjust_method("BH"), "BH")
    expect_identical(Pando:::.condition_validate_adjust_method("bh"), "BH")
    expect_error(
        Pando:::.condition_validate_adjust_method("holm"),
        'adjust_method = "BH"',
        fixed = TRUE
    )
    expect_equal(Pando:::.condition_validate_padj_threshold(0.05), 0.05)
    expect_equal(Pando:::.condition_validate_padj_threshold(0.1), 0.1)
    expect_equal(Pando:::.condition_validate_padj_threshold(0.5), 0.5)
    expect_equal(Pando:::.condition_validate_padj_threshold(0.999), 0.999)
    expect_error(
        Pando:::.condition_validate_padj_threshold(1),
        "\\(0, 1\\)"
    )
    expect_error(Pando:::.condition_validate_padj_threshold(0))
})

test_that("preliminary joint ridge builds the fit dictionary by significant union", {
    fit <- list(
        adjust_method = "BH",
        padj_threshold = 0.05,
        coefficients = data.frame(
            edge_id = rep(c("E1", "E2", "E3"), 2L),
            condition = rep(c("A", "B"), each = 3L),
            estimate = c(1, 0.5, 0.2, 0.8, 0.3, 0.1),
            std_err = 0.1,
            statistic = c(10, 5, 2, 8, 3, 1),
            estimable = TRUE,
            pval = c(0.001, 0.2, 0.8, 0.01, 0.7, 0.9),
            padj = c(0.003, 0.3, 0.8, 0.03, 0.9, 0.9),
            stringsAsFactors = FALSE
        )
    )
    class(fit) <- c("ConditionGRNFit", "list")
    screen <- Pando:::.condition_dictionary_screen(fit)

    expect_setequal(screen$keep_edge_ids, "E1")
    e1 <- screen$summary[screen$summary$edge_id == "E1", , drop = FALSE]
    expect_equal(e1$screening_n_significant_conditions, 2L)
    expect_equal(e1$screening_min_padj, 0.003)
    expect_identical(e1$screening_significant_conditions, "A;B")
    expect_true(all(c(
        "screen_estimate", "screen_std_err", "screen_statistic",
        "screen_pval", "screen_padj", "screen_estimable",
        "screen_significant"
    ) %in% colnames(screen$coefficients)))
    expect_identical(screen$coefficients$screen_significant,
                     c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE))
})

test_that("final refit uses screening support without second-round P values", {
    fit <- list(
        coefficients = data.frame(
            edge_id = c("E1", "E1", "E2", "E2"),
            condition = c("A", "B", "A", "B"),
            estimate = c(0.5, 0.4, 0.3, NA_real_),
            estimable = c(TRUE, TRUE, TRUE, FALSE),
            pval = NA_real_,
            padj = NA_real_,
            significant = FALSE,
            penalty_effect = 0,
            stringsAsFactors = FALSE
        )
    )
    class(fit) <- c("ConditionGRNFit", "list")
    screen <- data.frame(
        edge_id = c("E1", "E1", "E2", "E2"),
        condition = c("A", "B", "A", "B"),
        screen_estimate = c(0.45, 0.35, 0.25, 0.1),
        screen_std_err = 0.1,
        screen_statistic = c(4.5, 3.5, 2.5, 1),
        screen_pval = c(0.001, 0.01, 0.02, 0.4),
        screen_padj = c(0.002, 0.02, 0.04, 0.4),
        screen_estimable = TRUE,
        screen_significant = c(TRUE, TRUE, TRUE, FALSE),
        stringsAsFactors = FALSE
    )
    gated <- Pando:::.condition_apply_screen_support(fit, screen)

    expect_identical(gated$coefficients$significant,
                     c(TRUE, TRUE, TRUE, FALSE))
    expect_equal(gated$coefficients$penalty_effect,
                 c(0.5, 0.4, 0.3, 0))
    expect_true(all(is.na(gated$coefficients$pval)))
    expect_true(all(is.na(gated$coefficients$padj)))
    expect_identical(
        gated$projection_policy,
        "screen_bh_supported_refit_ridge_effects"
    )
    expect_false(gated$inference_performed)
})

test_that("final refit rejects accidental second-round inference", {
    fit <- list(
        coefficients = data.frame(
            edge_id = "E1", condition = "A", estimate = 0.5,
            estimable = TRUE, pval = 0.01, padj = 0.01,
            stringsAsFactors = FALSE
        )
    )
    class(fit) <- c("ConditionGRNFit", "list")
    screen <- data.frame(
        edge_id = "E1", condition = "A", screen_estimate = 0.4,
        screen_std_err = 0.1, screen_statistic = 4,
        screen_pval = 0.001, screen_padj = 0.001,
        screen_estimable = TRUE, screen_significant = TRUE,
        stringsAsFactors = FALSE
    )
    expect_error(
        Pando:::.condition_apply_screen_support(fit, screen),
        "must not contain second-round"
    )
})

test_that("screened dictionary preserves provenance and screening audit", {
    dictionary <- data.frame(
        edge_id = c("E1", "E2"),
        target = c("G", "G"),
        tf = c("TF1", "TF2"),
        region = c("P1", "P2"),
        candidate_index = 1:2,
        stringsAsFactors = FALSE
    )
    class(dictionary) <- c("PandoEdgeDictionary", "data.frame")
    attr(dictionary, "preprocessing_provenance_verified") <- TRUE
    summary <- data.frame(
        edge_id = c("E1", "E2"),
        screening_min_padj = c(0.01, 0.4),
        screening_n_significant_conditions = c(1L, 0L),
        screening_significant_conditions = c("A", ""),
        stringsAsFactors = FALSE
    )
    out <- Pando:::.condition_subset_dictionary(
        dictionary, "E1", summary
    )

    expect_s3_class(out, "PandoEdgeDictionary")
    expect_identical(out$edge_id, "E1")
    expect_equal(out$screening_min_padj, 0.01)
    expect_equal(out$screening_n_significant_conditions, 1L)
    expect_true(isTRUE(attr(out, "preprocessing_provenance_verified")))
})

test_that("default screening threshold remains 0.05 and is configurable", {
    default <- formals(Pando:::.pando_infer_condition_grn_one)$padj_threshold
    expect_equal(eval(default), 0.05)
})
