test_that("block pair information matches the dense reference", {
    set.seed(177)
    for (k in 2:4) {
        p <- 7L
        blocks <- lapply(seq_len(k), function(index) {
            x <- matrix(stats::rnorm((p + index) * p), p + index, p)
            crossprod(x)
        })
        conditions <- paste0("condition", seq_len(k))
        dense <- Pando:::.condition_pair_information(
            Pando:::.condition_blockdiag(blocks), p, conditions, 1e-10
        )
        blocked <- Pando:::.condition_pair_information_blocks(
            blocks, conditions, 1e-10
        )
        expect_identical(blocked[, setdiff(names(blocked),
                                           "pair_information")],
                         dense[, setdiff(names(dense), "pair_information")])
        expect_equal(blocked$pair_information, dense$pair_information,
                     tolerance = 1e-10)
    }
})

test_that("block pair information preserves dense rank-boundary semantics", {
    rank_tol <- 1e-10
    k <- 3L
    p <- 4L
    full_dimension <- k * p
    threshold <- rank_tol * full_dimension
    eigenvalues <- list(
        c(4, 2, 1, 0.5 * threshold),
        c(5, 3, 1, 0.99 * threshold),
        c(6, 2, 1, 1.01 * threshold)
    )
    blocks <- lapply(eigenvalues, diag)
    conditions <- paste0("condition", seq_len(k))
    dense <- Pando:::.condition_pair_information(
        Pando:::.condition_blockdiag(blocks), p, conditions, rank_tol
    )
    blocked <- Pando:::.condition_pair_information_blocks(
        blocks, conditions, rank_tol
    )
    expect_identical(blocked$pair_estimable, dense$pair_estimable)
    expect_equal(blocked$pair_information, dense$pair_information,
                 tolerance = 1e-12)
})

test_that("block pair information preserves the deterministic contrast tree", {
    set.seed(2755)
    k <- 4L
    p <- 9L
    blocks <- lapply(seq_len(k), function(index) {
        x <- matrix(stats::rnorm((p + index) * p), p + index, p)
        crossprod(x)
    })
    Q <- Pando:::.condition_blockdiag(blocks)
    conditions <- paste0("condition", seq_len(k))
    dense <- Pando:::.condition_identifiable_contrast_tree(
        Q, p, conditions, conditions[[1L]], rank_tol = 1e-10
    )
    blocked <- Pando:::.condition_identifiable_contrast_tree(
        Q, p, conditions, conditions[[1L]], rank_tol = 1e-10,
        q_blocks = blocks
    )

    expect_identical(blocked$D, dense$D)
    metadata_discrete <- setdiff(
        names(dense$metadata), "pair_information"
    )
    expect_identical(blocked$metadata[, metadata_discrete],
                     dense$metadata[, metadata_discrete])
    expect_equal(blocked$metadata$pair_information,
                 dense$metadata$pair_information, tolerance = 1e-10)
    expect_equal(blocked$pair_information$pair_information,
                 dense$pair_information$pair_information,
                 tolerance = 1e-10)
})
