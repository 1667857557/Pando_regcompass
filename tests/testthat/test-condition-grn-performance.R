test_that("implicit predictor centering preserves sparse matrices and scores", {
    x <- Matrix::Matrix(
        matrix(c(0, 1, 0, 2, 3, 0, 4, 0), nrow = 4),
        sparse = TRUE
    )
    colnames(x) <- c("e1", "e2")
    y <- c(1, 2, 3, 5)
    transform <- Pando:::.condition_build_balanced_transform(
        list(A = x, B = x + 1),
        list(A = y, B = y + 2)
    )
    scaled <- Pando:::.condition_apply_balanced_transform(
        list(A = x), list(A = y), transform
    )
    expect_s4_class(scaled$X$A, "dgCMatrix")
    beta <- c(0.4, -0.2)
    explicit <- sweep(as.matrix(x), 2L, transform$predictor_center, "-")
    explicit <- sweep(explicit, 2L, transform$predictor_scale, "/")
    implicit <- as.numeric(scaled$X$A %*% beta) -
        sum((transform$predictor_center / transform$predictor_scale) * beta)
    expect_equal(implicit, as.numeric(explicit %*% beta), tolerance = 1e-12)
})

test_that("fold statistics preserve the original transform arithmetic", {
    set.seed(4)
    x <- list(
        A = Matrix::Matrix(matrix(rnorm(120), 30, 4), sparse = TRUE),
        B = Matrix::Matrix(matrix(rnorm(120), 30, 4), sparse = TRUE)
    )
    colnames(x$A) <- colnames(x$B) <- paste0("e", seq_len(4))
    y <- list(A = rnorm(30), B = rnorm(30))
    stats <- Pando:::.condition_build_fold_statistics(x, y)
    cached <- Pando:::.condition_build_balanced_transform(
        x, y, fold_statistics = stats
    )
    uncached <- Pando:::.condition_build_balanced_transform(x, y)
    expect_identical(cached$predictor_center, uncached$predictor_center)
    expect_identical(cached$predictor_scale, uncached$predictor_scale)
    expect_identical(cached$response_center, uncached$response_center)
    expect_identical(cached$response_scale, uncached$response_scale)
})

test_that("verified estimability avoids recomputation without changing a path", {
    set.seed(9)
    x <- list(
        A = Matrix::Matrix(matrix(rnorm(160), 40, 4), sparse = TRUE),
        B = Matrix::Matrix(matrix(rnorm(160), 40, 4), sparse = TRUE)
    )
    colnames(x$A) <- colnames(x$B) <- paste0("e", seq_len(4))
    y <- list(A = rnorm(40), B = rnorm(40))
    mask <- matrix(TRUE, 4, 2, dimnames = list(colnames(x$A), c("A", "B")))
    verified <- Pando:::.condition_true_variance_mask(x, mask)
    lambda <- c(0.3, 0.12, 0.05)
    ordinary <- Pando:::.condition_fit_multitask_path(
        x, y, lambda, coefficient_mask = mask, backend = "R",
        max_iter = 2000L
    )
    cached <- Pando:::.condition_fit_multitask_path(
        x, y, lambda, coefficient_mask = mask, backend = "R",
        verified_estimability_mask = verified, max_iter = 2000L
    )
    for (index in seq_along(lambda)) {
        expect_identical(cached$fits[[index]]$beta, ordinary$fits[[index]]$beta)
        expect_identical(
            cached$fits[[index]]$intercept,
            ordinary$fits[[index]]$intercept
        )
    }
})

test_that("a lambda prefix preserves the selected warm-start fit", {
    set.seed(12)
    x <- list(
        A = Matrix::Matrix(matrix(rnorm(240), 60, 4), sparse = TRUE),
        B = Matrix::Matrix(matrix(rnorm(240), 60, 4), sparse = TRUE)
    )
    colnames(x$A) <- colnames(x$B) <- paste0("e", seq_len(4))
    y <- list(A = rnorm(60), B = rnorm(60))
    mask <- matrix(TRUE, 4, 2, dimnames = list(colnames(x$A), c("A", "B")))
    lambda <- c(0.4, 0.2, 0.1, 0.04)
    full <- Pando:::.condition_fit_multitask_path(
        x, y, lambda, coefficient_mask = mask, backend = "R",
        max_iter = 2000L
    )
    prefix <- Pando:::.condition_fit_multitask_path(
        x, y, lambda[1:3], coefficient_mask = mask, backend = "R",
        max_iter = 2000L
    )
    expect_identical(prefix$fits[[3]]$beta, full$fits[[3]]$beta)
    expect_identical(prefix$fits[[3]]$intercept, full$fits[[3]]$intercept)
})

test_that("direct Schur refit matches the alternating reference objective", {
    set.seed(17)
    x <- list(
        A = matrix(rnorm(180), nrow = 45, ncol = 4),
        B = matrix(rnorm(180), nrow = 45, ncol = 4)
    )
    colnames(x$A) <- colnames(x$B) <- paste0("e", seq_len(4))
    y <- list(
        A = 0.8 * x$A[, 1] - 0.3 * x$A[, 3] + rnorm(45, sd = 0.1),
        B = -0.4 * x$B[, 1] + 0.5 * x$B[, 2] + rnorm(45, sd = 0.1)
    )
    beta <- matrix(
        c(0.8, -0.4, 0, 0.5, -0.3, 0, 0, 0),
        nrow = 4,
        dimnames = list(colnames(x$A), c("A", "B"))
    )
    mask <- matrix(TRUE, 4, 2, dimnames = dimnames(beta))
    cache <- Pando:::.condition_make_refit_cache(x, y, "equal")
    reference <- Pando:::.condition_refit_shared_baseline_reference(
        x, y, beta, mask, ridge = 0.05, active_tol = 1e-8,
        condition_weight = "equal", cache = cache
    )
    direct <- Pando:::.condition_refit_shared_baseline(
        x, y, beta, mask, ridge = 0.05, active_tol = 1e-8,
        condition_weight = "equal", cache = cache, solver = "direct"
    )
    expect_equal(direct$beta, reference$beta, tolerance = 1e-6)
    expect_equal(direct$beta_shared, reference$beta_shared, tolerance = 1e-6)
    expect_equal(direct$intercept, reference$intercept, tolerance = 1e-6)
    expect_identical(direct$support_mask, reference$support_mask)
    expect_identical(direct$estimability_mask, reference$estimability_mask)
})

test_that("full projection reuses the prediction matrix product exactly", {
    set.seed(19)
    x <- Matrix::Matrix(matrix(rnorm(100), 25, 4), sparse = TRUE)
    beta <- c(0.2, -0.3, 0.5, 0)
    center <- c(0.1, -0.2, 0.3, 0.4)
    scale <- c(1.1, 0.8, 1.4, 0.9)
    mask <- c(TRUE, TRUE, TRUE, FALSE)
    linear <- as.numeric(x[, mask, drop = FALSE] %*% beta[mask])
    reused <- linear - sum((center[mask] / scale[mask]) * beta[mask])
    reference <- Pando:::.condition_projection_score(
        x, beta, center, scale, mask
    )
    expect_identical(reused, reference)
})

test_that("compiled and R sparse-group paths are numerically equivalent", {
    skip_if(!Pando:::.condition_native_solver_available())
    set.seed(21)
    x <- list(
        A = Matrix::Matrix(matrix(rnorm(180), 45, 4), sparse = TRUE),
        B = Matrix::Matrix(matrix(rnorm(180), 45, 4), sparse = TRUE)
    )
    colnames(x$A) <- colnames(x$B) <- paste0("e", seq_len(4))
    y <- list(
        A = 0.7 * x$A[, 1] - 0.2 * x$A[, 4] + rnorm(45, sd = 0.15),
        B = -0.5 * x$B[, 1] + 0.3 * x$B[, 2] + rnorm(45, sd = 0.15)
    )
    y <- lapply(y, as.numeric)
    mask <- matrix(
        TRUE, 4, 2,
        dimnames = list(colnames(x$A), c("A", "B"))
    )
    lambda <- c(0.2, 0.08)
    fit_r <- Pando:::.condition_fit_multitask_path(
        x, y, lambda, alpha = 0.5, condition_mix = 0.5,
        condition_weight = "equal", coefficient_mask = mask,
        max_iter = 2000L, backend = "R"
    )
    fit_cpp <- Pando:::.condition_fit_multitask_path(
        x, y, lambda, alpha = 0.5, condition_mix = 0.5,
        condition_weight = "equal", coefficient_mask = mask,
        max_iter = 2000L, backend = "cpp"
    )
    for (index in seq_along(lambda)) {
        expect_equal(
            fit_cpp$fits[[index]]$beta,
            fit_r$fits[[index]]$beta,
            tolerance = 1e-6
        )
        expect_equal(
            fit_cpp$fits[[index]]$intercept,
            fit_r$fits[[index]]$intercept,
            tolerance = 1e-6
        )
        expect_identical(
            abs(fit_cpp$fits[[index]]$beta) > 1e-8,
            abs(fit_r$fits[[index]]$beta) > 1e-8
        )
    }
})

test_that("condition GRN runtime has no zzz function overrides", {
    expect_length(list.files("R", pattern = "^zzz.*condition_grn.*[.]R$"), 0L)
    expect_false(exists(
        ".condition_refit_stability_diagnostics",
        envir = asNamespace("Pando"),
        inherits = FALSE
    ))
})
