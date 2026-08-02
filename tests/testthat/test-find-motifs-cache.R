test_that("find_motifs exposes persistent cache controls", {
    method <- Pando:::find_motifs.GRNData
    expect_identical(formals(method)$reuse_cache, TRUE)
    expect_match(
        deparse(formals(method)$cache_dir),
        "Pando.motif_cache_dir",
        fixed = TRUE
    )
    body_text <- paste(deparse(body(method)), collapse = "\n")
    expect_match(body_text, ".pando_motif_cache_key", fixed = TRUE)
    expect_match(body_text, ".pando_read_motif_cache", fixed = TRUE)
    expect_match(body_text, ".pando_write_motif_cache", fixed = TRUE)
})

test_that("motif cache digests are deterministic and content-sensitive", {
    first <- Pando:::.pando_motif_cache_digest(list(a = 1, b = "x"))
    second <- Pando:::.pando_motif_cache_digest(list(a = 1, b = "x"))
    changed <- Pando:::.pando_motif_cache_digest(list(a = 2, b = "x"))
    expect_identical(first, second)
    expect_false(identical(first, changed))
})

test_that("motif genome signatures accept S4 genome-like objects", {
    class_name <- "PandoMotifGenomeSignatureFixture"
    if (!methods::isClass(class_name)) {
        methods::setClass(
            class_name,
            slots = c(
                pkgname = "character",
                organism = "character",
                provider = "character"
            )
        )
    }
    genome <- methods::new(
        class_name,
        pkgname = "fixture.genome",
        organism = "Homo sapiens",
        provider = "test"
    )

    signature <- Pando:::.pando_motif_genome_signature(genome)

    expect_identical(signature$fields$pkgname, "fixture.genome")
    expect_identical(signature$fields$organism, "Homo sapiens")
    expect_identical(signature$fields$provider, "test")
})

test_that("invalid motif cache controls fail before motif scanning", {
    method <- Pando:::find_motifs.GRNData
    expect_error(
        method(NULL, pfm = NULL, genome = NULL, reuse_cache = NA),
        "`reuse_cache` must be either TRUE or FALSE.",
        fixed = TRUE
    )
    expect_error(
        method(NULL, pfm = NULL, genome = NULL, cache_dir = ""),
        "`cache_dir` must be NULL or one non-empty path.",
        fixed = TRUE
    )
})
