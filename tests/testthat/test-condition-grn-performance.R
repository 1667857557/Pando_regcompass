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

test_that("cached hierarchical refit matches the reference implementation", {
    set.seed(17)
    x <- list(
        A = matrix(rnorm(120), nrow = 30, ncol = 4),
        B = matrix(rnorm(120), nrow = 30, ncol = 4)
    )
    colnames(x$A) <- colnames(x$B) <- paste0("e", seq_len(4))
    y <- list(
        A = 0.8 * x$A[, 1] - 0.3 * x$A[, 3] + rnorm(30, sd = 0.1),
        B = -0.4 * x$B[, 1] + 0.5 * x$B[, 2] + rnorm(30, sd = 0.1)
    )
    beta <- matrix(
        c(0.8, -0.4, 0, 0.5, -0.3, 0, 0, 0),
        nrow = 4,
        dimnames = list(colnames(x$A), c("A", "B"))
    )
    mask <- matrix(TRUE, 4, 2, dimnames = dimnames(beta))
    reference <- Pando:::.condition_refit_shared_baseline_reference(
        x, y, beta, mask, ridge = 0.05, active_tol = 1e-8,
        condition_weight = "equal"
    )
    cache <- Pando:::.condition_make_refit_cache(x, y, "equal")
    cached <- Pando:::.condition_refit_shared_baseline(
        x, y, beta, mask, ridge = 0.05, active_tol = 1e-8,
        condition_weight = "equal", cache = cache
    )
    expect_equal(cached$beta, reference$beta, tolerance = 1e-8)
    expect_equal(cached$beta_shared, reference$beta_shared, tolerance = 1e-8)
    expect_equal(cached$intercept, reference$intercept, tolerance = 1e-8)
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

test_that("stability and sensitivity refits are absent from the runtime path", {
    expect_null(Pando:::.condition_refit_stability_diagnostics())
})
