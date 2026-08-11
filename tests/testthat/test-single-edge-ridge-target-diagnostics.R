test_that("single-edge multi-task ridge keeps condition-by-edge dimensions", {
    x <- list(
        A = matrix(seq_len(8), ncol = 1L,
                   dimnames = list(paste0("A", seq_len(8)), "edge1")),
        B = matrix(seq_len(9) + 0.5, ncol = 1L,
                   dimnames = list(paste0("B", seq_len(9)), "edge1"))
    )
    y <- list(
        A = 1 + 0.3 * as.numeric(x$A[, 1L]),
        B = 2 + 0.5 * as.numeric(x$B[, 1L])
    )
    scaling <- .condition_ridge_scaling(x, 1e-8)
    fit <- .condition_ridge_fit(
        x = x, y = y, scaling = scaling,
        lambda = 0.1, fusion_ratio = 1,
        min_residual_df = 1L, inference = TRUE
    )
    expect_identical(fit$status, "ok")
    expect_equal(dim(fit$zero_variance), c(2L, 1L))
    expect_identical(rownames(fit$zero_variance), c("A", "B"))
    expect_identical(colnames(fit$zero_variance), "edge1")
    expect_equal(dim(fit$beta), c(2L, 1L))
    expect_equal(dim(fit$se), c(2L, 1L))
})

test_that("target payload errors identify phase and target", {
    old <- getOption("Pando.target_progress_context", NULL)
    options(Pando.target_progress_context = list(
        stage = "unit_phase", label = "celltype"
    ))
    on.exit(options(Pando.target_progress_context = old), add = TRUE)
    keys <- stats::setNames(c("good", "bad"), c("good", "bad"))
    expect_error(
        .pando_target_payload_map(
            keys = keys,
            build_payload = function(key) list(key = key),
            worker = function(payload) {
                if (identical(payload$key, "bad")) stop("boom")
                payload$key
            },
            parallel = FALSE,
            verbose = FALSE
        ),
        "phase=unit_phase:celltype; target=bad"
    )
})

test_that("target payload progress reports semantic phase and completion", {
    old <- getOption("Pando.target_progress_context", NULL)
    options(Pando.target_progress_context = list(
        stage = "unit_phase", label = "celltype"
    ))
    on.exit(options(Pando.target_progress_context = old), add = TRUE)
    keys <- stats::setNames(c("a", "b"), c("a", "b"))
    expect_message(
        .pando_target_payload_map(
            keys = keys,
            build_payload = function(key) list(key = key),
            worker = function(payload) payload$key,
            parallel = FALSE,
            verbose = TRUE
        ),
        "Pando target phase=unit_phase:celltype"
    )
})
