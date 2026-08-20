test_that("single-edge E-star/JSE keeps condition-by-edge dimensions", {
    set.seed(12)
    x <- list(
        A = matrix(seq_len(20), ncol = 1L,
                   dimnames = list(paste0("A", seq_len(20)), "edge1")),
        B = matrix(seq_len(24) + 0.5, ncol = 1L,
                   dimnames = list(paste0("B", seq_len(24)), "edge1"))
    )
    y <- list(
        A = 1 + 0.3 * as.numeric(x$A[, 1L]) + rnorm(20, sd = 0.1),
        B = 2 + 0.5 * as.numeric(x$B[, 1L]) + rnorm(24, sd = 0.1)
    )
    scaling <- .condition_ridge_scaling(x, 1e-8)
    fit <- .condition_scheme_e_fit(
        x = x, y = y, scaling = scaling,
        min_residual_df = 1L, inference = TRUE,
        reference_condition = "A"
    )
    expect_identical(fit$status, "ok")
    expect_equal(dim(fit$zero_variance), c(2L, 1L))
    expect_identical(rownames(fit$zero_variance), c("A", "B"))
    expect_identical(colnames(fit$zero_variance), "edge1")
    expect_equal(dim(fit$beta), c(2L, 1L))
    expect_equal(dim(fit$inference_se), c(2L, 1L))
    expect_length(fit$contrast_identifiable, 1L)
    expect_equal(fit$penalty_value, 0.25)
    expect_identical(fit$inference_schema,
                     "scheme_e_fusion_component_joint_refit_v1")
})

test_that("target payload errors identify phase and target without local worker closures", {
    keys <- stats::setNames("bad", "bad")
    expect_error(
        .pando_target_payload_map(
            keys = keys,
            build_payload = function(key) list(),
            worker_name = ".pando_ridge_target_worker",
            parallel = FALSE,
            verbose = FALSE,
            phase = "unit_phase",
            label = "celltype"
        ),
        "phase=unit_phase:celltype; target=bad"
    )
})

test_that("target payload progress reports semantic phase and completion", {
    keys <- stats::setNames(c("a", "b"), c("a", "b"))
    expect_message(
        .pando_target_payload_map(
            keys = keys,
            build_payload = function(key) list(skip = TRUE, target = key),
            worker_name = ".pando_discovery_target_worker",
            parallel = FALSE,
            verbose = TRUE,
            phase = "unit_phase",
            label = "celltype"
        ),
        "Pando target phase=unit_phase:celltype"
    )
})

test_that("canonical conditional target calls E-star directly", {
    body_text <- paste(deparse(body(.condition_ridge_target)), collapse = "\n")
    expect_match(body_text, ".condition_ridge_fit", fixed = TRUE)
    expect_match(body_text, "penalty_family", fixed = TRUE)
    expect_false(grepl("lambda =", body_text, fixed = TRUE))
    expect_false(any(c("fusion_ratio", "scheme_e_z", "z") %in%
                     names(formals(.condition_scheme_e_fit))))
})
