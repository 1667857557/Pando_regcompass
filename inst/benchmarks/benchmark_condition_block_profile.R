library(Pando)

benchmark_profile <- function(k, p, iterations = 3L) {
    set.seed(2755 + k + p)
    q <- lapply(seq_len(k), function(index) {
        x <- matrix(stats::rnorm((p + index) * p), p + index, p)
        crossprod(x)
    })
    h_blocks <- lapply(seq_len(k), function(index) stats::rnorm(p))
    h <- unlist(h_blocks, use.names = FALSE)
    D <- matrix(0, (k - 1L) * p, k * p)
    cursor <- 0L
    for (edge in seq_len(p)) for (condition in 2:k) {
        cursor <- cursor + 1L
        D[cursor, edge] <- -1
        D[cursor, (condition - 1L) * p + edge] <- 1
    }
    elapsed <- function(fun) {
        min(replicate(iterations, system.time(fun())[["elapsed"]]))
    }
    dense_fun <- function() {
        Q <- Pando:::.condition_blockdiag(q)
        A <- kronecker(matrix(1, k, 1L), diag(p))
        geometry <- Pando:::.condition_q_orthogonal_decomposition(
            Q, A, D, 1e-10
        )
        list(H = crossprod(geometry$R, Q %*% geometry$R),
             r = crossprod(geometry$R, h))
    }
    block_fun <- function() {
        geometry <- Pando:::.condition_q_orthogonal_decomposition_blocks(
            q, D, 1e-10
        )
        list(
            H = Reduce(`+`, Map(function(R_block, q_block) {
                crossprod(R_block, q_block %*% R_block)
            }, geometry$R_blocks, q)),
            r = Reduce(`+`, Map(function(R_block, h_block) {
                crossprod(R_block, h_block)
            }, geometry$R_blocks, h_blocks))
        )
    }
    dense <- elapsed(dense_fun)
    blocked <- elapsed(block_fun)
    Q <- Pando:::.condition_blockdiag(q)
    data.frame(
        conditions = k, edges = p, dense_seconds = dense,
        block_seconds = blocked, speedup = dense / blocked,
        dense_Q_mib = as.numeric(object.size(Q)) / 1024^2,
        block_Q_mib = sum(vapply(
            q, function(block) as.numeric(object.size(block)), numeric(1)
        )) / 1024^2
    )
}

result <- do.call(rbind, list(
    benchmark_profile(2L, 20L),
    benchmark_profile(3L, 75L),
    benchmark_profile(4L, 150L, iterations = 1L)
))
print(result, row.names = FALSE)
