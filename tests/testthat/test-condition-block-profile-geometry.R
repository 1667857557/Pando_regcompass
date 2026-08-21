test_that("block profile geometry matches the dense reference", {
    set.seed(177)
    for (k in 2:4) {
        p <- 6L
        q <- lapply(seq_len(k), function(index) {
            x <- matrix(stats::rnorm((p + index) * p), p + index, p)
            crossprod(x)
        })
        Q <- Pando:::.condition_blockdiag(q)
        A <- kronecker(matrix(1, k, 1L), diag(p))
        D <- matrix(0, (k - 1L) * p, k * p)
        cursor <- 0L
        for (edge in seq_len(p)) {
            for (condition in 2:k) {
                cursor <- cursor + 1L
                D[cursor, edge] <- -1
                D[cursor, (condition - 1L) * p + edge] <- 1
            }
        }
        h_blocks <- lapply(seq_len(k), function(index) {
            stats::rnorm(p)
        })
        h <- unlist(h_blocks, use.names = FALSE)
        dense <- Pando:::.condition_q_orthogonal_decomposition(
            Q, A, D, 1e-10
        )
        blocked <- Pando:::.condition_q_orthogonal_decomposition_blocks(
            q, D, 1e-10
        )

        expect_equal(blocked$q_scale, dense$q_scale, tolerance = 0)
        expect_equal(blocked$shared_inverse, dense$shared_inverse,
                     tolerance = 1e-11)
        expect_equal(blocked$B0, dense$B0, tolerance = 1e-12)
        expect_equal(blocked$R, dense$R, tolerance = 1e-10)
        expect_equal(blocked$dr_error, dense$dr_error, tolerance = 1e-10)
        expect_equal(blocked$orthogonality_error,
                     dense$orthogonality_error, tolerance = 1e-10)

        H_dense <- crossprod(dense$R, Q %*% dense$R)
        H_block <- Reduce(`+`, Map(function(R_block, q_block) {
            crossprod(R_block, q_block %*% R_block)
        }, blocked$R_blocks, q))
        r_dense <- as.numeric(crossprod(dense$R, h))
        r_block <- as.numeric(Reduce(`+`, Map(function(R_block, h_block) {
            crossprod(R_block, h_block)
        }, blocked$R_blocks, h_blocks)))
        mu_dense <- as.numeric(
            dense$shared_inverse %*%
                (crossprod(A, h) / dense$q_scale)
        )
        mu_block <- as.numeric(
            blocked$shared_inverse %*%
                (Reduce(`+`, h_blocks) / blocked$q_scale)
        )
        delta <- stats::rnorm(nrow(D))

        expect_equal(H_block, H_dense, tolerance = 1e-9)
        expect_equal(r_block, r_dense, tolerance = 1e-10)
        expect_equal(mu_block, mu_dense, tolerance = 1e-11)
        expect_equal(
            as.numeric(rep(mu_block, times = k) + blocked$R %*% delta),
            as.numeric(A %*% mu_dense + dense$R %*% delta),
            tolerance = 1e-9
        )
    }
})

test_that("block profile geometry preserves rank-deficient shared systems", {
    q <- list(diag(c(4, 0, 2)), diag(c(3, 0, 5)), diag(c(6, 0, 1)))
    k <- length(q); p <- nrow(q[[1L]])
    Q <- Pando:::.condition_blockdiag(q)
    A <- kronecker(matrix(1, k, 1L), diag(p))
    D <- matrix(0, (k - 1L) * p, k * p)
    cursor <- 0L
    for (edge in seq_len(p)) for (condition in 2:k) {
        cursor <- cursor + 1L
        D[cursor, edge] <- -1
        D[cursor, (condition - 1L) * p + edge] <- 1
    }
    dense <- Pando:::.condition_q_orthogonal_decomposition(Q, A, D, 1e-10)
    blocked <- Pando:::.condition_q_orthogonal_decomposition_blocks(
        q, D, 1e-10
    )
    expect_equal(blocked$shared_inverse, dense$shared_inverse,
                 tolerance = 1e-12)
    expect_equal(blocked$R, dense$R, tolerance = 1e-10)
})
