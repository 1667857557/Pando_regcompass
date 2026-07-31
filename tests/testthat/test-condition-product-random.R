test_that("compiled sparse products match 64 randomized designs", {
    vl <- c(-3, -1, 0.5, 2, 4)
    vr <- c(-2, -0.5, 1, 3, 5)
    for (seed in seq_len(64L)) {
        set.seed(seed)
        nr <- sample(9:41, 1L)
        nl <- sample(2:9, 1L)
        nc <- sample(2:11, 1L)
        ne <- sample(1:35, 1L)
        left <- Matrix::rsparsematrix(
            nr, nl, density = runif(1L, 0.05, 0.6),
            rand.x = function(n) sample(vl, n, replace = TRUE)
        )
        right <- Matrix::rsparsematrix(
            nr, nc, density = runif(1L, 0.03, 0.5),
            rand.x = function(n) sample(vr, n, replace = TRUE)
        )
        rownames(left) <- rownames(right) <- paste0("c", seq_len(nr))
        colnames(left) <- paste0("L", seq_len(nl))
        colnames(right) <- paste0("R", seq_len(nc))
        li <- sample(seq_len(nl), ne, replace = TRUE)
        ri <- sample(seq_len(nc), ne, replace = TRUE)
        reference <- left[, li, drop = FALSE] * right[, ri, drop = FALSE]
        compiled <- Pando:::.condition_product_matrix_cpp(
            left, right, as.integer(li), as.integer(ri)
        )
        dimnames(compiled) <- dimnames(reference)
        expect_identical(
            as.matrix(compiled), as.matrix(reference), info = seed
        )
        expect_identical(
            which(compiled != 0, arr.ind = TRUE),
            which(reference != 0, arr.ind = TRUE), info = seed
        )
    }
})
