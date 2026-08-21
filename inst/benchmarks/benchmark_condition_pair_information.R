library(Pando)

benchmark_one <- function(k, p, iterations = 3L) {
    set.seed(177 + k + p)
    blocks <- lapply(seq_len(k), function(index) {
        x <- matrix(stats::rnorm((p + index) * p), p + index, p)
        crossprod(x)
    })
    Q <- Pando:::.condition_blockdiag(blocks)
    conditions <- paste0("condition", seq_len(k))
    elapsed <- function(fun) {
        min(replicate(iterations, system.time(fun())[["elapsed"]]))
    }
    dense <- elapsed(function() Pando:::.condition_pair_information(
        Q, p, conditions, 1e-10
    ))
    blocked <- elapsed(function() Pando:::.condition_pair_information_blocks(
        blocks, conditions, 1e-10
    ))
    data.frame(
        conditions = k, edges = p, dense_seconds = dense,
        block_seconds = blocked, speedup = dense / blocked,
        dense_matrix_mib = as.numeric(object.size(Q)) / 1024^2,
        block_matrix_mib = sum(vapply(
            blocks, function(block) as.numeric(object.size(block)), numeric(1)
        )) / 1024^2
    )
}

result <- do.call(rbind, list(
    benchmark_one(2L, 20L),
    benchmark_one(3L, 100L),
    benchmark_one(4L, 300L, iterations = 1L)
))
print(result, row.names = FALSE)
