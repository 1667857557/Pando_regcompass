test_that("installed sparse condition kernels use strict C++ execution", {
    expect_true(exists(
        ".condition_product_matrix_cpp",
        envir = asNamespace("Pando"),
        inherits = FALSE
    ))
    expect_false(exists(
        ".condition_product_matrix_r",
        envir = asNamespace("Pando"),
        inherits = FALSE
    ))
    expect_true(is.loaded(
        "_Pando_condition_product_matrix_cpp", PACKAGE = "Pando"
    ))
    expect_true(is.loaded(
        "_Pando_condition_fit_multitask_path_cpp", PACKAGE = "Pando"
    ))
    cells <- paste0("c", seq_len(20L))
    base <- seq_len(20L) - mean(seq_len(20L))
    X <- methods::as(Matrix::Matrix(
        cbind(x1 = base, x2 = 2 * base, x3 = -base),
        sparse = TRUE,
        dimnames = list(cells, c("x1", "x2", "x3"))
    ), "dgCMatrix")
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
    fit <- Pando:::.condition_fit_multitask_path_core(
        X_list = X_list,
        y_list = y_list,
        lambda = 1,
        alpha = 0.5,
        condition_mix = 0.5,
        condition_weight = "equal",
        coefficient_mask = mask,
        max_iter = 1L,
        tol_objective = 1e-7,
        tol_coef = 1e-6,
        keep_history = FALSE,
        backend = "cpp"
    )
    expect_identical(fit$fits[[1L]]$backend, "cpp_eigen_sparse_fista")
    expect_true(all(is.finite(fit$fits[[1L]]$beta)))
    expect_error(
        Pando:::.condition_fit_multitask_path_core(
            X_list, y_list, 1, 0.5, 0.5, "equal", mask,
            1L, 1e-7, 1e-6, FALSE, backend = "R"
        ),
        "requires the compiled C\\+\\+ solver"
    )
})

test_that("strict native kernels reject unsupported or invalid sparse layouts", {
    symmetric <- methods::as(Matrix::Diagonal(3), "dsCMatrix")
    expect_error(
        Pando:::.condition_product_matrix_cpp(
            symmetric, symmetric, as.integer(1:3), as.integer(1:3)
        ),
        "dgCMatrix"
    )

    invalid <- methods::as(Matrix::sparseMatrix(
        i = c(1L, 3L),
        j = c(1L, 1L),
        x = c(1, 2),
        dims = c(3L, 1L)
    ), "dgCMatrix")
    invalid@i <- rev(invalid@i)
    expect_error(
        Pando:::.condition_product_matrix_cpp(
            invalid, invalid, 1L, 1L
        ),
        "non-increasing row indices"
    )
})

test_that("strict native product indices and arithmetic fail fast", {
    X <- methods::as(Matrix::Diagonal(2), "dgCMatrix")
    expect_error(
        Pando:::.condition_product_matrix_cpp(X, X, NA_integer_, 1L),
        "must not contain missing values"
    )
    huge <- methods::as(Matrix::sparseMatrix(
        i = 1L, j = 1L, x = .Machine$double.xmax, dims = c(1L, 1L)
    ), "dgCMatrix")
    expect_error(
        Pando:::.condition_product_matrix_cpp(huge, huge, 1L, 1L),
        "produced a non-finite value"
    )
})
