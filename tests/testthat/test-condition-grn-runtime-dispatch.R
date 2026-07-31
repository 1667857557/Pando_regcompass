test_that("structural-zero helpers do not override canonical numerical paths", {
    source <- readLines(
        file.path("R", "condition_grn_structural_zero.R"),
        warn = FALSE
    )
    forbidden <- c(
        "[.]condition_make_lambda_path[[:space:]]*<-[[:space:]]*function",
        "[.]condition_fit_multitask_path[[:space:]]*<-[[:space:]]*function",
        "[.]condition_select_lambda_nested[[:space:]]*<-[[:space:]]*function",
        "[.]condition_refit_shared_baseline[[:space:]]*<-[[:space:]]*function",
        "[.]condition_nested_crossfit_within_cell_type[[:space:]]*<-[[:space:]]*function"
    )
    expect_false(any(vapply(
        forbidden,
        function(pattern) any(grepl(pattern, source)),
        logical(1)
    )))

    expect_true(all(
        c("backend", "verified_estimability_mask") %in%
            names(formals(Pando:::.condition_fit_multitask_path))
    ))
    expect_true(all(
        c("cache", "estimability_verified", "solver") %in%
            names(formals(Pando:::.condition_refit_shared_baseline))
    ))
})

test_that("candidate retention does not make exact-zero columns estimable", {
    zero <- Matrix::Matrix(matrix(0, nrow = 12, ncol = 1), sparse = TRUE)
    signal_a <- Matrix::Matrix(matrix(seq_len(12), ncol = 1), sparse = TRUE)
    signal_b <- Matrix::Matrix(matrix(rev(seq_len(12)), ncol = 1), sparse = TRUE)
    x <- list(
        A = cbind(zero, signal_a),
        B = cbind(zero, signal_b)
    )
    colnames(x$A) <- colnames(x$B) <- c("zero", "signal")
    y <- list(A = seq_len(12), B = rev(seq_len(12)))
    mask <- matrix(
        TRUE, 2, 2,
        dimnames = list(c("zero", "signal"), c("A", "B"))
    )

    center <- as.numeric(Matrix::colMeans(x$A))
    fold_variance <- Pando:::.condition_population_variance(
        x$A, center = center
    )
    candidate_variance <- Pando:::.condition_population_variance(x$A)
    expect_equal(fold_variance[[1L]], 0)
    expect_gt(candidate_variance[[1L]], .Machine$double.eps)

    statistics <- Pando:::.condition_build_fold_statistics(x, y, mask)
    expect_false(any(statistics$estimability_mask["zero", ]))
    expect_true(all(statistics$estimability_mask["signal", ]))

    verified <- Pando:::.condition_true_variance_mask(x, mask)
    expect_identical(verified, statistics$estimability_mask)
})

test_that("optimized solver and refit arguments survive installed dispatch", {
    set.seed(29)
    x <- list(
        A = Matrix::Matrix(matrix(rnorm(120), 30, 4), sparse = TRUE),
        B = Matrix::Matrix(matrix(rnorm(120), 30, 4), sparse = TRUE)
    )
    colnames(x$A) <- colnames(x$B) <- paste0("e", seq_len(4))
    y <- list(A = rnorm(30), B = rnorm(30))
    mask <- matrix(
        TRUE, 4, 2,
        dimnames = list(colnames(x$A), c("A", "B"))
    )
    verified <- Pando:::.condition_true_variance_mask(x, mask)
    path <- Pando:::.condition_fit_multitask_path(
        x,
        y,
        lambda = c(0.2, 0.08),
        coefficient_mask = mask,
        verified_estimability_mask = verified,
        backend = "R",
        max_iter = 2000L
    )
    expect_identical(path$fits[[1L]]$backend, "R_reference_fista")

    cache <- Pando:::.condition_make_refit_cache(x, y, "equal")
    refit <- Pando:::.condition_refit_shared_baseline(
        x,
        y,
        beta_selection = path$fits[[2L]]$beta,
        estimability_mask = verified,
        ridge = 0.01,
        condition_weight = "equal",
        cache = cache,
        estimability_verified = TRUE,
        solver = "direct"
    )
    expect_identical(
        refit$common_metric,
        "pooled_weighted_predictor_gram_direct_schur"
    )
    expect_true(is.logical(refit$estimability_mask))
})
