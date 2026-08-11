test_that("candidate target payload keeps only target-relevant multiome data", {
    cells <- paste0("c", 1:4)
    gene_data <- Matrix::Matrix(matrix(
        seq_len(16), nrow = 4,
        dimnames = list(cells, c("TARGET1", "TARGET2", "TF1", "OTHER"))
    ), sparse = TRUE)
    peak_data <- Matrix::Matrix(matrix(
        seq_len(12), nrow = 4,
        dimnames = list(cells, c("r1", "r2", "r3"))
    ), sparse = TRUE)
    peaks2gene <- Matrix::Matrix(
        matrix(c(1, 0, 0, 0, 1, 0), nrow = 2, byrow = TRUE,
               dimnames = list(c("TARGET1", "TARGET2"), c("r1", "r2", "r3"))),
        sparse = TRUE
    )
    peaks2motif <- Matrix::Matrix(
        matrix(c(1, 0, 0, 1, 0, 1), nrow = 3, byrow = TRUE,
               dimnames = list(c("r1", "r2", "r3"), c("m1", "m2"))),
        sparse = TRUE
    )
    motif2tf <- Matrix::Matrix(
        matrix(c(1, 0, 0, 1), nrow = 2, byrow = TRUE,
               dimnames = list(c("m1", "m2"), c("TF1", "OTHER"))),
        sparse = TRUE
    )
    prepared <- list(
        gene_data = gene_data,
        peak_data = peak_data,
        features = c("TARGET1", "TARGET2"),
        peaks2gene = peaks2gene,
        peaks2motif = peaks2motif,
        motif2tf = motif2tf,
        region_map = data.frame(
            region = c("r1", "r2", "r3"),
            atac_feature_id = c("p1", "p2", "p3"),
            stringsAsFactors = FALSE
        )
    )

    payload <- Pando:::.pando_discovery_target_payload(
        prepared = prepared, cells = cells, target = "TARGET1",
        source_label = "standard", source_type = "global",
        tf_cor = 0.1, peak_cor = 0.05
    )

    expect_false(payload$skip)
    expect_identical(payload$prepared$features, "TARGET1")
    expect_setequal(colnames(payload$prepared$gene_data), c("TARGET1", "TF1"))
    expect_identical(colnames(payload$prepared$peak_data), "r1")
    expect_identical(rownames(payload$prepared$peaks2motif), "r1")
    expect_identical(colnames(payload$prepared$motif2tf), "TF1")
    expect_false("TARGET2" %in% colnames(payload$prepared$gene_data))
    expect_false("r2" %in% colnames(payload$prepared$peak_data))
})

test_that("ridge target payload excludes other targets and regions", {
    cells <- paste0("c", 1:6)
    prepared <- list(
        gene_data = Matrix::Matrix(matrix(
            seq_len(24), nrow = 6,
            dimnames = list(cells, c("TARGET1", "TARGET2", "TF1", "TF2"))
        ), sparse = TRUE),
        peak_data = Matrix::Matrix(matrix(
            seq_len(18), nrow = 6,
            dimnames = list(cells, c("r1", "r2", "r3"))
        ), sparse = TRUE)
    )
    dictionary <- data.frame(
        target = c("TARGET1", "TARGET2"),
        tf = c("TF1", "TF2"),
        region = c("r1", "r2"),
        edge_id = c("e1", "e2"),
        candidate_index = 1:2,
        stringsAsFactors = FALSE
    )
    groups <- list(A = cells[1:3], B = cells[4:6])
    folds <- list(A = c(1L, 2L, 1L), B = c(1L, 2L, 1L))
    attr(folds, "nfolds") <- 2L

    payload <- Pando:::.pando_ridge_target_payload(
        prepared = prepared,
        edge_dictionary = dictionary,
        target = "TARGET1",
        cells = groups,
        folds = folds,
        control = list(),
        min_residual_df = 1L,
        rank_action = "mark"
    )

    expect_identical(payload$edges$target, "TARGET1")
    expect_setequal(colnames(payload$prepared$gene_data), c("TARGET1", "TF1"))
    expect_identical(colnames(payload$prepared$peak_data), "r1")
    expect_false("TARGET2" %in% colnames(payload$prepared$gene_data))
    expect_false("r2" %in% colnames(payload$prepared$peak_data))
})

test_that("parallel target dispatcher is namespace-level and does not capture mapper frame", {
    expect_identical(
        environment(Pando:::.pando_target_execute_task),
        asNamespace("Pando")
    )
    body_text <- paste(
        deparse(body(Pando:::.pando_target_payload_map)), collapse = "\n"
    )
    expect_match(body_text, ".pando_target_execute_task", fixed = TRUE)
    expect_false(grepl("run_one <-", body_text, fixed = TRUE))
    expect_false(grepl("worker(task$payload)", body_text, fixed = TRUE))
})

test_that("target payload mapper preserves target names without one-pass alias overrides", {
    keys <- stats::setNames(c("a", "b", "c"), c("a", "b", "c"))
    result <- Pando:::.pando_target_payload_map(
        keys = keys,
        build_payload = function(key) list(skip = TRUE, target = key),
        worker_name = ".pando_discovery_target_worker",
        parallel = FALSE,
        verbose = FALSE,
        phase = "unit"
    )
    expect_identical(names(result), names(keys))
    expect_true(is.function(Pando:::.condition_ridge_refit_contract_one_pass))
    expect_false(exists(
        ".condition_ridge_refit_contract_compact",
        envir = asNamespace("Pando"), inherits = FALSE
    ))
})
