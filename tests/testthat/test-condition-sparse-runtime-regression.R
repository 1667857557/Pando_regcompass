test_that("installed sparse condition kernels accept dgCMatrix inputs", {
    skip_if_not(exists(
        ".condition_product_matrix_cpp",
        envir = asNamespace("Pando"),
        inherits = FALSE
    ))
    cells <- paste0("c", seq_len(20L))
    base <- seq_len(20L) - mean(seq_len(20L))
    X <- Matrix::Matrix(
        cbind(x1 = base, x2 = 2 * base, x3 = -base),
        sparse = TRUE,
        dimnames = list(cells, c("x1", "x2", "x3"))
    )
    X <- methods::as(X, "dgCMatrix")
    product <- Pando:::.condition_product_matrix_cpp(
        X, X, as.integer(1:3), as.integer(1:3)
    )
    expect_s4_class(product, "dgCMatrix")
    expect_identical(
        as.matrix(product),
        as.matrix(X[, 1:3, drop = FALSE] * X[, 1:3, drop = FALSE])
    )

    X_list <- list(A = X, B = X)
    y_list <- list(A = base, B = rev(base))
    mask <- matrix(
        TRUE, nrow = ncol(X), ncol = length(X_list),
        dimnames = list(colnames(X), names(X_list))
    )
    lambda <- 1
    alpha <- 0.5
    fit <- Pando:::.condition_fit_multitask_path_cpp(
        X_list = X_list,
        y_list = y_list,
        lambda = lambda,
        alpha = alpha,
        condition_mix = 0.5,
        condition_weight = "equal",
        coefficient_mask = mask,
        max_iter = 1L,
        tol_objective = 1e-7,
        tol_coef = 1e-6,
        keep_history = FALSE
    )
    expect_true(is.list(fit$fits[[1L]]))
    expect_true(all(is.finite(fit$fits[[1L]]$beta)))
    centered_frobenius <- sum(X@x * X@x)
    old_frobenius_step <- 1 / (
        lambda * (1 - alpha) +
            2 * centered_frobenius / nrow(X)
    )
    expect_gt(fit$fits[[1L]]$step, old_frobenius_step)
})

test_that("sparse condition kernels reject unsupported S4 matrices clearly", {
    symmetric <- methods::as(
        Matrix::Diagonal(3), "dsCMatrix"
    )
    expect_error(
        Pando:::.condition_product_matrix_cpp(
            symmetric, symmetric, as.integer(1:3), as.integer(1:3)
        ),
        "dgCMatrix"
    )
})
