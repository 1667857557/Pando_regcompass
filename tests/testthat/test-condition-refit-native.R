test_that("compiled refit path is registered and fail-fast", {
    expect_true(exists(
        ".condition_refit_path_cpp",
        envir = asNamespace("Pando"),
        inherits = FALSE
    ))
    expect_true(is.loaded(
        "_Pando_condition_refit_path_cpp", PACKAGE = "Pando"
    ))
    expect_true(Pando:::.condition_native_refit_available())

    x <- list(
        A = matrix(rnorm(80), nrow = 20, ncol = 4),
        B = matrix(rnorm(80), nrow = 20, ncol = 4)
    )
    colnames(x$A) <- colnames(x$B) <- paste0("e", seq_len(4))
    y <- list(A = rnorm(20), B = rnorm(20))
    beta <- matrix(
        c(0.5, 0, -0.3, 0, -0.4, 0.2, 0, 0),
        nrow = 4,
        dimnames = list(colnames(x$A), names(x))
    )
    mask <- matrix(TRUE, 4, 2, dimnames = dimnames(beta))
    cache <- Pando:::.condition_make_refit_cache(x, y, "equal")

    expect_error(
        Pando:::.condition_refit_shared_baseline(
            x, y, beta, mask, ridge = 0.05,
            condition_weight = "equal", cache = cache,
            solver = "alternating"
        ),
        "requires the compiled C\\+\\+ direct refit kernel"
    )
    expect_error(
        Pando:::.condition_refit_shared_baseline_path(
            x, y, list(beta), mask, ridge = 0,
            condition_weight = "equal", cache = cache
        ),
        "finite positive"
    )
})

test_that("batched C++ refit matches the R numerical oracle", {
    set.seed(20260801)
    n <- 80L
    p <- 6L
    x <- list(
        Control = matrix(rnorm(n * p), nrow = n, ncol = p),
        Treatment = matrix(rnorm(n * p), nrow = n, ncol = p)
    )
    predictor_names <- paste0("edge", seq_len(p))
    colnames(x$Control) <- colnames(x$Treatment) <- predictor_names
    y <- list(
        Control = 0.9 * x$Control[, 1] - 0.35 * x$Control[, 3] +
            rnorm(n, sd = 0.08),
        Treatment = -0.65 * x$Treatment[, 1] + 0.5 * x$Treatment[, 2] +
            rnorm(n, sd = 0.08)
    )
    beta_path <- list(
        matrix(
            c(
                0.9, 0, -0.35, 0, 0, 0,
                -0.65, 0.5, 0, 0, 0, 0
            ),
            nrow = p,
            dimnames = list(predictor_names, names(x))
        ),
        matrix(
            c(
                0.7, 0, -0.2, 0.1, 0, 0,
                -0.45, 0.35, 0, 0, 0.08, 0
            ),
            nrow = p,
            dimnames = list(predictor_names, names(x))
        ),
        matrix(
            c(
                0.4, 0, 0, 0, 0, 0,
                -0.3, 0.2, 0, 0, 0, 0
            ),
            nrow = p,
            dimnames = list(predictor_names, names(x))
        )
    )
    mask <- matrix(
        TRUE, p, length(x),
        dimnames = list(predictor_names, names(x))
    )
    mask[6, "Treatment"] <- FALSE
    ridge <- c(0.08, 0.04, 0.015)
    cache <- Pando:::.condition_make_refit_cache(x, y, "equal")

    native <- Pando:::.condition_refit_shared_baseline_path(
        X_list = x,
        y_list = y,
        beta_path = beta_path,
        estimability_mask = mask,
        ridge = ridge,
        active_tol = 1e-8,
        condition_weight = "equal",
        cache = cache,
        estimability_verified = TRUE
    )
    reference <- lapply(seq_along(beta_path), function(index) {
        Pando:::.condition_refit_shared_baseline_reference(
            X_list = x,
            y_list = y,
            beta_selection = beta_path[[index]],
            estimability_mask = mask,
            ridge = ridge[[index]],
            active_tol = 1e-8,
            condition_weight = "equal",
            max_iter = 2000L,
            tol = 1e-10,
            cache = cache,
            estimability_verified = TRUE
        )
    })

    expect_length(native, length(beta_path))
    for (index in seq_along(native)) {
        expect_identical(native[[index]]$support_mask,
                         reference[[index]]$support_mask)
        expect_identical(native[[index]]$estimability_mask,
                         reference[[index]]$estimability_mask)
        expect_equal(native[[index]]$beta,
                     reference[[index]]$beta, tolerance = 2e-6)
        expect_equal(native[[index]]$beta_shared,
                     reference[[index]]$beta_shared, tolerance = 2e-6)
        expect_equal(native[[index]]$intercept,
                     reference[[index]]$intercept, tolerance = 2e-6)
        expect_true(native[[index]]$converged)
        expect_identical(
            native[[index]]$backend,
            "cpp_eigen_dense_llt_batched_path"
        )
    }
    expect_lt(native[[1]]$beta[1, "Control"] *
              native[[1]]$beta[1, "Treatment"], 0)
    expect_true(is.na(native[[1]]$beta_condition[6, "Treatment"]))
})

test_that("single-refit API delegates to the compiled path without changing output", {
    set.seed(42)
    x <- list(
        A = matrix(rnorm(240), nrow = 60, ncol = 4),
        B = matrix(rnorm(240), nrow = 60, ncol = 4)
    )
    colnames(x$A) <- colnames(x$B) <- paste0("e", seq_len(4))
    y <- list(A = rnorm(60), B = rnorm(60))
    beta <- matrix(
        c(0.6, 0, -0.2, 0, -0.4, 0.3, 0, 0),
        nrow = 4,
        dimnames = list(colnames(x$A), names(x))
    )
    mask <- matrix(TRUE, 4, 2, dimnames = dimnames(beta))
    cache <- Pando:::.condition_make_refit_cache(x, y, "equal")
    single <- Pando:::.condition_refit_shared_baseline(
        x, y, beta, mask, ridge = 0.03,
        condition_weight = "equal", cache = cache,
        estimability_verified = TRUE
    )
    batched <- Pando:::.condition_refit_shared_baseline_path(
        x, y, list(beta), mask, ridge = 0.03,
        condition_weight = "equal", cache = cache,
        estimability_verified = TRUE
    )[[1L]]
    expect_identical(single, batched)
})
