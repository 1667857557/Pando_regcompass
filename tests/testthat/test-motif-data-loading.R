test_that("bundled motif-to-TF data load from the Pando package", {
    motif2tf <- Pando:::.load_default_motif2tf()

    expect_s3_class(motif2tf, "data.frame")
    expect_gt(nrow(motif2tf), 0L)
    expect_gte(ncol(motif2tf), 2L)
    expect_true(all(!is.na(motif2tf[[1L]])))
    expect_true(all(!is.na(motif2tf[[2L]])))
})
