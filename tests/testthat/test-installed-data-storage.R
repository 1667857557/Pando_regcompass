test_that("bundled datasets use standalone installed files", {
    description <- utils::packageDescription("Pando")
    expect_identical(tolower(as.character(description$LazyData)), "false")
    expect_identical(
        as.character(description[["Config/Pando/DataStorage"]]),
        "standalone-RData"
    )

    data_dir <- system.file("data", package = "Pando")
    expect_true(nzchar(data_dir))
    expect_true(dir.exists(data_dir))
    expected <- c(
        "EnsDb.Hsapiens.v93.annot.UCSC.hg38.RData",
        "SCREEN.ccRE.UCSC.hg38.RData",
        "motif2tf.RData",
        "motifs.RData",
        "phastConsElements20Mammals.UCSC.hg38.RData"
    )
    expect_true(all(file.exists(file.path(data_dir, expected))))
    expect_false(file.exists(file.path(data_dir, "Rdata.rdb")))
    expect_false(file.exists(file.path(data_dir, "Rdata.rdx")))

    data_environment <- new.env(parent = emptyenv())
    expect_silent(utils::data(
        list = "motif2tf",
        package = "Pando",
        envir = data_environment
    ))
    expect_true(exists(
        "motif2tf", envir = data_environment, inherits = FALSE
    ))
    motif2tf <- get(
        "motif2tf", envir = data_environment, inherits = FALSE
    )
    expect_true(is.data.frame(motif2tf))
    expect_true(ncol(motif2tf) >= 2L)
})
