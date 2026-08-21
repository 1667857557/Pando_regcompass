test_that("target checkpoints skip completed payload construction on resume", {
    directory <- tempfile("pando_checkpoint_")
    dir.create(directory)
    on.exit(unlink(directory, recursive = TRUE), add = TRUE)
    keys <- stats::setNames(c("A", "B"), c("A", "B"))
    builds <- 0L
    builder <- function(key) {
        builds <<- builds + 1L
        list(skip = TRUE, target = key)
    }

    first <- Pando:::.pando_target_payload_map(
        keys = keys,
        build_payload = builder,
        worker_name = ".pando_discovery_target_worker",
        parallel = FALSE,
        verbose = FALSE,
        checkpoint_dir = directory,
        checkpoint_fingerprint = "fixture-v1"
    )
    expect_identical(builds, 2L)
    expect_length(list.files(directory, recursive = TRUE), 2L)

    builds <- 0L
    second <- Pando:::.pando_target_payload_map(
        keys = keys,
        build_payload = builder,
        worker_name = ".pando_discovery_target_worker",
        parallel = FALSE,
        verbose = FALSE,
        checkpoint_dir = directory,
        checkpoint_fingerprint = "fixture-v1"
    )
    expect_identical(builds, 0L)
    expect_identical(second, first)
})

test_that("target checkpoints are isolated by exact input fingerprint", {
    directory <- tempfile("pando_checkpoint_")
    dir.create(directory)
    on.exit(unlink(directory, recursive = TRUE), add = TRUE)
    keys <- stats::setNames("A", "A")
    builds <- 0L
    builder <- function(key) {
        builds <<- builds + 1L
        list(skip = TRUE, target = key)
    }

    for (fingerprint in c("input-one", "input-two")) {
        Pando:::.pando_target_payload_map(
            keys = keys,
            build_payload = builder,
            worker_name = ".pando_discovery_target_worker",
            parallel = FALSE,
            verbose = FALSE,
            checkpoint_dir = directory,
            checkpoint_fingerprint = fingerprint
        )
    }
    expect_identical(builds, 2L)
    expect_length(list.dirs(directory, recursive = FALSE), 2L)
})

test_that("condition GRN exposes checkpoint and resume controls", {
    method <- getS3method("infer_condition_grn", "GRNData")
    expect_identical(formals(method)$checkpoint_dir, NULL)
    expect_identical(formals(method)$resume, TRUE)
})
