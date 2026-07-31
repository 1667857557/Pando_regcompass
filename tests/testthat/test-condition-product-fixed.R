test_that("compiled sparse product matches the original expression", {
    cells <- paste0("c", seq_len(8L))
    left <- Matrix::sparseMatrix(
        i = c(1, 2, 4, 7, 1, 3, 6, 8),
        j = c(1, 1, 1, 1, 2, 2, 2, 2),
        x = c(1, -2, 0.5, 3, 4, -1, 2, 5),
        dims = c(8, 3),
        dimnames = list(cells, c("TF1", "TF2", "TF0"))
    )
    right <- Matrix::sparseMatrix(
        i = c(1, 3, 4, 7, 2, 3, 6, 8),
        j = c(1, 1, 1, 1, 2, 2, 2, 2),
        x = c(2, 1, -2, 0.5, 3, -4, 1.5, 2),
        dims = c(8, 3),
        dimnames = list(cells, c("P1", "P2", "P0"))
    )
    li <- c(1L, 2L, 1L, 2L, 1L, 3L)
    ri <- c(1L, 1L, 2L, 2L, 1L, 3L)
    reference <- left[, li, drop = FALSE] * right[, ri, drop = FALSE]
    compiled <- Pando:::.condition_product_matrix_cpp(left, right, li, ri)
    dimnames(compiled) <- dimnames(reference)
    expect_s4_class(compiled, "dgCMatrix")
    expect_identical(dimnames(compiled), dimnames(reference))
    expect_identical(as.matrix(compiled), as.matrix(reference))
    expect_identical(
        which(compiled != 0, arr.ind = TRUE),
        which(reference != 0, arr.ind = TRUE)
    )
})
