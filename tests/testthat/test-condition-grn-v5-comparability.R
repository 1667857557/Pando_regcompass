test_that("equal-condition transforms are invariant to whole-condition duplication", {
    x_a <- matrix(
        c(1, 2, 3, 4, 2, 4, 6, 8),
        ncol = 2,
        dimnames = list(NULL, c("x1", "x2"))
    )
    x_b <- matrix(
        c(10, 12, 14, 16, 3, 5, 7, 9),
        ncol = 2,
        dimnames = list(NULL, c("x1", "x2"))
    )
    y_a <- c(1, 2, 4, 5)
    y_b <- c(11, 12, 14, 15)
    original <- Pando:::.condition_build_balanced_transform(
        list(A = x_a, B = x_b),
        list(A = y_a, B = y_b)
    )
    duplicated <- Pando:::.condition_build_balanced_transform(
        list(A = x_a[rep(seq_len(nrow(x_a)), 2L), ], B = x_b),
        list(A = rep(y_a, 2L), B = y_b)
    )
    expect_equal(
        original$predictor_center,
        duplicated$predictor_center,
        tolerance = 1e-12
    )
    expect_equal(
        original$predictor_scale,
        duplicated$predictor_scale,
        tolerance = 1e-12
    )
    expect_equal(
        original$response_center,
        duplicated$response_center,
        tolerance = 1e-12
    )
    expect_equal(
        original$response_scale,
        duplicated$response_scale,
        tolerance = 1e-12
    )
    expect_identical(
        original$condition_weights,
        c(A = 0.5, B = 0.5)
    )
    original_scaled <- Pando:::.condition_apply_balanced_transform(
        list(A = x_a, B = x_b), list(A = y_a, B = y_b), original
    )
    duplicated_scaled <- Pando:::.condition_apply_balanced_transform(
        list(A = x_a[rep(seq_len(nrow(x_a)), 2L), ], B = x_b),
        list(A = rep(y_a, 2L), B = y_b),
        duplicated
    )
    expect_equal(
        original_scaled$X$A,
        duplicated_scaled$X$A[seq_len(nrow(x_a)), , drop = FALSE],
        tolerance = 1e-12
    )
    expect_equal(
        original_scaled$y$A,
        duplicated_scaled$y$A[seq_along(y_a)],
        tolerance = 1e-12
    )
    mask <- matrix(
        TRUE, 2, 2, dimnames = list(c("x1", "x2"), c("A", "B"))
    )
    original_fit <- Pando:::.condition_fit_multitask_path(
        original_scaled$X, original_scaled$y, lambda = 0.1,
        alpha = 0.5, condition_mix = 0.5, condition_weight = "equal",
        coefficient_mask = mask, max_iter = 2000L
    )$fits[[1L]]
    duplicated_fit <- Pando:::.condition_fit_multitask_path(
        duplicated_scaled$X, duplicated_scaled$y, lambda = 0.1,
        alpha = 0.5, condition_mix = 0.5, condition_weight = "equal",
        coefficient_mask = mask, max_iter = 2000L
    )$fits[[1L]]
    expect_equal(
        original_fit$beta, duplicated_fit$beta, tolerance = 1e-8
    )
    expect_identical(
        abs(original_fit$beta) > 1e-8,
        abs(duplicated_fit$beta) > 1e-8
    )
})

test_that("nested cross-fitting assigns every cell exactly once", {
    set.seed(4)
    x <- list(
        A = matrix(rnorm(48), nrow = 24, ncol = 2),
        B = matrix(rnorm(48), nrow = 24, ncol = 2)
    )
    colnames(x$A) <- colnames(x$B) <- c("e1", "e2")
    y <- list(
        A = 0.7 * x$A[, 1] + rnorm(24, sd = 0.2),
        B = -0.4 * x$B[, 1] + rnorm(24, sd = 0.2)
    )
    mask <- matrix(
        TRUE, 2, 2,
        dimnames = list(c("e1", "e2"), c("A", "B"))
    )
    fit <- Pando:::.condition_nested_crossfit_within_cell_type(
        X_list = x,
        y_list = y,
        lambda = c(0.2, 0.05),
        coefficient_mask = mask,
        outer_nfolds = 3L,
        inner_nfolds = 2L,
        comparison_conditions = c("A", "B"),
        seed = 19L,
        max_iter = 300L
    )
    expect_true(all(unlist(fit$oof_assignment_count) == 1L))
    expect_true(all(is.finite(unlist(fit$oof_prediction))))
    expect_identical(
        fit$projection_origin,
        "outer_condition_stratified_cell_oof"
    )
    expect_identical(fit$primary_projection, "condition_full_oof")
    expect_true(fit$projection_used_for_penalty)
    expect_false(fit$full_fit_projection_used_for_penalty)
    expect_length(fit$fold_transform, 3L)
    expect_true(all(vapply(
        fit$fold_transform,
        function(x) isTRUE(x$training_fold_only),
        logical(1)
    )))
})

test_that("OOF assignment uses structural zeros when no edge is estimable", {
    x <- list(
        A = matrix(0, nrow = 12, ncol = 1,
                   dimnames = list(NULL, "closed_edge")),
        B = matrix(0, nrow = 12, ncol = 1,
                   dimnames = list(NULL, "closed_edge"))
    )
    y <- list(A = seq_len(12), B = seq_len(12) + 2)
    mask <- matrix(
        TRUE, 1, 2,
        dimnames = list("closed_edge", c("A", "B"))
    )
    fit <- Pando:::.condition_nested_crossfit_within_cell_type(
        X_list = x,
        y_list = y,
        lambda = 0.1,
        coefficient_mask = mask,
        outer_nfolds = 3L,
        inner_nfolds = 2L,
        comparison_conditions = c("A", "B"),
        seed = 23L
    )
    expect_true(all(unlist(fit$oof_assignment_count) == 1L))
    expect_true(all(is.finite(unlist(fit$oof_prediction))))
    expect_equal(unlist(fit$projection_condition_full_oof), rep(0, 24))
    expect_equal(unlist(fit$projection_common_oof), rep(0, 24))
    expect_true(all(vapply(
        fit$fold_support,
        function(x) identical(
            x$projection_status,
            "intercept_only_all_predictors_structural_zero"
        ),
        logical(1)
    )))
})

test_that("a NULL coefficient mask becomes a named logical matrix", {
    x <- list(
        A = matrix(rnorm(20), ncol = 2),
        B = matrix(rnorm(20), ncol = 2)
    )
    colnames(x$A) <- colnames(x$B) <- c("e1", "e2")
    mask <- Pando:::.condition_training_estimability(x, NULL)
    expect_type(mask, "logical")
    expect_identical(dimnames(mask), list(c("e1", "e2"), c("A", "B")))
    expect_false(anyNA(mask))
})
