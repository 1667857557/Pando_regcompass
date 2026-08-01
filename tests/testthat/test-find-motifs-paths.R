test_that("find_motifs defaults to presence/absence motif matching", {
    method <- Pando:::find_motifs.GRNData

    expect_identical(formals(method)$exact_positions, FALSE)

    body_text <- paste(deparse(body(method)), collapse = "\n")
    expect_match(body_text, "if \\(exact_positions\\)")
    expect_match(body_text, "Signac::CreateMotifMatrix")
    expect_match(body_text, "score = FALSE", fixed = TRUE)
    expect_match(body_text, "use.counts = FALSE", fixed = TRUE)
    expect_match(body_text, "Signac::CreateMotifObject")
    expect_match(body_text, "positions = NULL", fixed = TRUE)

    expect_error(
        method(NULL, pfm = NULL, genome = NULL, exact_positions = NA),
        "must be either TRUE or FALSE",
        fixed = TRUE
    )
})

test_that("find_motifs retains the exact-position Signac path", {
    body_text <- paste(
        deparse(body(Pando:::find_motifs.GRNData)),
        collapse = "\n"
    )

    expect_match(body_text, "Signac::AddMotifs")
    expect_match(body_text, "verbose = verbose", fixed = TRUE)
})
