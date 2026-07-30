test_that("exact-zero predictors are retained but not coefficient-estimable", {
    x <- list(
        A = matrix(0, 12, 1, dimnames = list(NULL, "peak_closed")),
        B = matrix(seq_len(12), 12, 1,
                   dimnames = list(NULL, "peak_closed"))
    )
    mask <- matrix(
        TRUE, 1, 2,
        dimnames = list("peak_closed", c("A", "B"))
    )
    estimable <- Pando:::.condition_training_estimability(x, mask)
    expect_false(estimable["peak_closed", "A"])
    expect_true(estimable["peak_closed", "B"])
    expect_gt(
        Pando:::.condition_population_variance(x$A)[[1L]],
        .Machine$double.eps
    )
})

test_that("condition-full OOF uses zero on the non-estimable side", {
    set.seed(9)
    common_a <- rnorm(30)
    common_b <- rnorm(30)
    unique_b <- rnorm(30)
    x <- list(
        A = cbind(common = common_a, unique = 0, closed = 0),
        B = cbind(common = common_b, unique = unique_b, closed = 0)
    )
    y <- list(
        A = 0.8 * common_a,
        B = 0.8 * common_b + 2 * unique_b
    )
    mask <- matrix(
        TRUE, 3, 2,
        dimnames = list(c("common", "unique", "closed"), c("A", "B"))
    )
    fit <- Pando:::.condition_nested_crossfit_within_cell_type(
        X_list = x,
        y_list = y,
        lambda = c(1e-4, 1e-6),
        coefficient_mask = mask,
        outer_nfolds = 3L,
        inner_nfolds = 2L,
        comparison_conditions = c("A", "B"),
        seed = 27L,
        max_iter = 3000L
    )
    full_a <- fit$projection_condition_full_oof$A
    common_a_score <- fit$projection_common_oof$A
    full_b <- fit$projection_condition_full_oof$B
    common_b_score <- fit$projection_common_oof$B
    expect_true(all(is.finite(c(full_a, common_a_score, full_b, common_b_score))))
    expect_equal(full_a, common_a_score, tolerance = 1e-8)
    expect_gt(max(abs(full_b - common_b_score)), 1e-4)
    expect_identical(fit$primary_projection, "condition_full_oof")
    expect_identical(
        fit$nonestimable_projection_policy,
        "structural_zero_by_condition"
    )
    expect_true(all(vapply(fit$fold_support, function(one) {
        all(one$projectable_structural_zero_mask["closed", ])
    }, logical(1))))
})

test_that("a fully closed predictor produces finite zero OOF projections", {
    x <- list(
        A = matrix(0, 12, 1, dimnames = list(NULL, "closed")),
        B = matrix(0, 12, 1, dimnames = list(NULL, "closed"))
    )
    y <- list(A = seq_len(12), B = seq_len(12) + 1)
    mask <- matrix(
        TRUE, 1, 2,
        dimnames = list("closed", c("A", "B"))
    )
    fit <- Pando:::.condition_nested_crossfit_within_cell_type(
        X_list = x,
        y_list = y,
        lambda = c(0.1, 0.01),
        coefficient_mask = mask,
        outer_nfolds = 3L,
        inner_nfolds = 2L,
        comparison_conditions = c("A", "B"),
        seed = 31L
    )
    expect_equal(unlist(fit$projection_condition_full_oof), rep(0, 24))
    expect_equal(unlist(fit$projection_common_oof), rep(0, 24))
    expect_equal(unlist(fit$projection_global_common_oof), rep(0, 24))
    expect_true(all(unlist(fit$oof_assignment_count) == 1L))
    expect_true(all(vapply(fit$fold_support, function(one) {
        all(one$projection_support_mask["closed", ])
    }, logical(1))))
})

test_that("refit keeps structural-zero coefficients unavailable", {
    x <- list(
        A = matrix(0, 12, 1, dimnames = list(NULL, "edge")),
        B = matrix(seq_len(12), 12, 1, dimnames = list(NULL, "edge"))
    )
    y <- list(A = seq_len(12), B = seq_len(12))
    selection <- matrix(
        c(0, 0.5), 1, 2,
        dimnames = list("edge", c("A", "B"))
    )
    fit <- Pando:::.condition_refit_shared_baseline(
        X_list = x,
        y_list = y,
        beta_selection = selection,
        estimability_mask = matrix(
            TRUE, 1, 2,
            dimnames = list("edge", c("A", "B"))
        ),
        ridge = 0.1
    )
    expect_false(fit$estimability_mask["edge", "A"])
    expect_true(fit$estimability_mask["edge", "B"])
    expect_true(is.na(fit$beta_condition["edge", "A"]))
    expect_true(is.finite(fit$beta_condition["edge", "B"]))
})
