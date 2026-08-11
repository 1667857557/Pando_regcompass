# Canonical target-level execution helpers for condition and standard ridge GRNs.
#
# These helpers do not replace or wrap existing public methods. They are called
# directly by the condition candidate-discovery and ridge-refit implementations.
# Worker tasks contain only one target-specific compact payload plus scalar
# diagnostic context; namespace-level worker functions avoid serializing the
# mapper frame, accumulated outputs, or full prepared multiome matrices.

.pando_condition_target_bpparam <- function() {
    value <- getOption("Pando.condition_target_BPPARAM", NULL)
    if (identical(value, FALSE) || is.null(value)) return(NULL)
    .pando_validate_bpparam(value)
    value
}

.pando_target_worker_limit <- function(parallel = FALSE) {
    if (!isTRUE(parallel)) return(1L)
    BPPARAM <- tryCatch(.pando_condition_target_bpparam(), error = function(e) NULL)
    if (!is.null(BPPARAM)) {
        value <- tryCatch(BiocParallel::bpnworkers(BPPARAM), error = function(e) 1L)
        return(max(1L, as.integer(value)))
    }
    value <- tryCatch(foreach::getDoParWorkers(), error = function(e) 1L)
    max(1L, as.integer(value))
}

.pando_target_phase_label <- function(phase, label = NULL) {
    phase <- as.character(phase)
    if (length(phase) != 1L || is.na(phase) || !nzchar(phase)) {
        phase <- "target_work"
    }
    label <- as.character(label %||% "")
    if (length(label) != 1L || is.na(label)) label <- ""
    if (nzchar(label)) paste0(phase, ":", label) else phase
}

.pando_target_execute_task <- function(task) {
    on.exit(invisible(gc(verbose = FALSE, full = TRUE)), add = TRUE)
    if (!is.list(task) || !is.character(task$worker_name) ||
        length(task$worker_name) != 1L || !nzchar(task$worker_name) ||
        !is.character(task$key) || length(task$key) != 1L ||
        !is.character(task$phase_label) || length(task$phase_label) != 1L) {
        stop("Invalid compact Pando target task.", call. = FALSE)
    }
    worker <- get(task$worker_name, envir = asNamespace("Pando"), inherits = FALSE)
    if (!is.function(worker)) {
        stop("Unknown Pando target worker: ", task$worker_name, call. = FALSE)
    }
    tryCatch(
        worker(task$payload),
        error = function(error) {
            stop(
                "Pando target task failed [phase=", task$phase_label,
                "; target=", task$key, "]: ", conditionMessage(error),
                call. = FALSE
            )
        }
    )
}

.pando_target_payload_map <- function(
    keys, build_payload, worker_name, parallel = FALSE, verbose = TRUE,
    phase = "target_work", label = NULL) {
    if (!length(keys)) return(list())
    if (!is.function(build_payload)) {
        stop("Target payload mapping requires a payload builder.", call. = FALSE)
    }
    if (!is.character(worker_name) || length(worker_name) != 1L ||
        !nzchar(worker_name)) {
        stop("`worker_name` must name one namespace-level target worker.",
             call. = FALSE)
    }
    if (!exists(worker_name, envir = asNamespace("Pando"), inherits = FALSE) ||
        !is.function(get(worker_name, envir = asNamespace("Pando"),
                         inherits = FALSE))) {
        stop("Unknown namespace-level target worker: ", worker_name,
             call. = FALSE)
    }

    key_names <- names(keys)
    if (is.null(key_names) || any(!nzchar(key_names))) {
        key_names <- as.character(keys)
    }
    phase_label <- .pando_target_phase_label(phase, label)
    out <- vector("list", length(keys))
    names(out) <- key_names
    batch_size <- min(length(keys), .pando_target_worker_limit(parallel))
    starts <- seq.int(1L, length(keys), by = batch_size)
    n_batches <- length(starts)

    if (isTRUE(verbose)) {
        message(
            "Pando target phase=", phase_label,
            " | started | targets=", length(keys),
            ";batch_size=", batch_size,
            ";batches=", n_batches,
            ";parallel=", isTRUE(parallel)
        )
    }

    for (batch_index in seq_along(starts)) {
        start <- starts[[batch_index]]
        index <- seq.int(start, min(length(keys), start + batch_size - 1L))
        payloads <- lapply(index, function(i) build_payload(keys[[i]]))
        tasks <- lapply(seq_along(payloads), function(j) {
            list(
                key = key_names[index[[j]]],
                phase_label = phase_label,
                worker_name = worker_name,
                payload = payloads[[j]]
            )
        })
        names(tasks) <- key_names[index]

        chunk <- if (isTRUE(parallel)) {
            BPPARAM <- .pando_condition_target_bpparam()
            if (!is.null(BPPARAM)) {
                BiocParallel::bplapply(
                    tasks, .pando_target_execute_task, BPPARAM = BPPARAM
                )
            } else {
                map_par(
                    tasks, .pando_target_execute_task,
                    parallel = TRUE, verbose = FALSE
                )
            }
        } else {
            lapply(tasks, .pando_target_execute_task)
        }
        names(chunk) <- names(tasks)
        out[index] <- chunk

        if (isTRUE(verbose)) {
            message(
                "Pando target phase=", phase_label,
                " | completed=", max(index), "/", length(keys),
                " (", sprintf("%.1f", 100 * max(index) / length(keys)), "%)",
                ";batch=", batch_index, "/", n_batches
            )
        }

        payloads <- NULL
        tasks <- NULL
        chunk <- NULL
        invisible(gc(verbose = FALSE, full = TRUE))
    }
    out
}

.condition_discover_one_target_prepared <- function(
    prepared, cells, target, source_label, source_type, tf_cor, peak_cor) {
    if (!target %in% rownames(prepared$peaks2gene)) return(NULL)
    gene_peak_mask <- as.logical(prepared$peaks2gene[target, , drop = TRUE])
    if (!any(gene_peak_mask)) return(NULL)

    target_x <- prepared$gene_data[cells, target, drop = FALSE]
    peak_x <- prepared$peak_data[cells, gene_peak_mask, drop = FALSE]
    peak_cor_matrix <- methods::as(
        sparse_cor(peak_x, target_x), "generalMatrix"
    )
    peak_cor_matrix[is.na(peak_cor_matrix)] <- 0
    selected_regions <- rownames(peak_cor_matrix)[
        abs(peak_cor_matrix[, 1L]) > peak_cor
    ]
    if (!length(selected_regions)) return(NULL)

    peak_motifs <- prepared$peaks2motif[selected_regions, , drop = FALSE]
    region_tfs <- lapply(selected_regions, function(region) {
        motif_present <- as.logical(peak_motifs[region, , drop = TRUE])
        if (!any(motif_present)) return(character())
        present <- sparseMatrixStats::colMaxs(
            prepared$motif2tf[motif_present, , drop = FALSE]
        )
        setdiff(colnames(prepared$motif2tf)[as.logical(present)], target)
    })
    names(region_tfs) <- selected_regions
    candidate_tfs <- unique(unlist(region_tfs, use.names = FALSE))
    if (!length(candidate_tfs)) return(NULL)

    tf_x <- prepared$gene_data[cells, candidate_tfs, drop = FALSE]
    tf_cor_matrix <- methods::as(
        sparse_cor(tf_x, target_x), "generalMatrix"
    )
    tf_cor_matrix[is.na(tf_cor_matrix)] <- 0
    selected_tfs <- rownames(tf_cor_matrix)[
        abs(tf_cor_matrix[, 1L]) > tf_cor
    ]
    if (!length(selected_tfs)) return(NULL)

    edge_rows <- lapply(selected_regions, function(region) {
        tf <- intersect(region_tfs[[region]], selected_tfs)
        if (!length(tf)) return(NULL)
        data.frame(
            target = target,
            tf = tf,
            region = region,
            atac_feature_id = prepared$region_map$atac_feature_id[
                match(region, prepared$region_map$region)
            ],
            peak_target_cor = as.numeric(peak_cor_matrix[region, 1L]),
            tf_target_cor = as.numeric(tf_cor_matrix[tf, 1L]),
            source_label = source_label,
            source_type = source_type,
            stringsAsFactors = FALSE
        )
    })
    edge_rows <- edge_rows[!vapply(edge_rows, is.null, logical(1))]
    if (!length(edge_rows)) return(NULL)
    do.call(rbind, edge_rows)
}

.pando_discovery_target_payload <- function(
    prepared, cells, target, source_label, source_type, tf_cor, peak_cor) {
    target <- as.character(target)
    skip <- !target %in% rownames(prepared$peaks2gene) ||
        !target %in% colnames(prepared$gene_data)
    if (skip) return(list(skip = TRUE, target = target))

    gene_peak_mask <- as.logical(prepared$peaks2gene[target, , drop = TRUE])
    if (!any(gene_peak_mask)) return(list(skip = TRUE, target = target))
    regions <- colnames(prepared$peaks2gene)[gene_peak_mask]
    peak_motifs <- prepared$peaks2motif[regions, , drop = FALSE]
    motif_keep <- as.logical(sparseMatrixStats::colMaxs(peak_motifs * 1))
    motif_ids <- colnames(peak_motifs)[motif_keep]
    if (!length(motif_ids)) return(list(skip = TRUE, target = target))

    tf_supported <- sparseMatrixStats::colMaxs(
        prepared$motif2tf[motif_ids, , drop = FALSE]
    )
    tfs <- colnames(prepared$motif2tf)[as.logical(tf_supported)]
    tfs <- setdiff(intersect(tfs, colnames(prepared$gene_data)), target)
    if (!length(tfs)) return(list(skip = TRUE, target = target))

    region_map <- prepared$region_map[
        prepared$region_map$region %in% regions, , drop = FALSE
    ]
    compact <- list(
        gene_data = prepared$gene_data[
            cells, unique(c(target, tfs)), drop = FALSE
        ],
        peak_data = prepared$peak_data[cells, regions, drop = FALSE],
        features = target,
        peaks2gene = prepared$peaks2gene[target, regions, drop = FALSE],
        peaks2motif = prepared$peaks2motif[regions, motif_ids, drop = FALSE],
        motif2tf = prepared$motif2tf[motif_ids, tfs, drop = FALSE],
        region_map = region_map
    )
    list(
        skip = FALSE,
        prepared = compact,
        cells = rownames(compact$gene_data),
        target = target,
        source_label = source_label,
        source_type = source_type,
        tf_cor = tf_cor,
        peak_cor = peak_cor
    )
}

.pando_discovery_target_worker <- function(payload) {
    if (isTRUE(payload$skip)) return(NULL)
    .condition_discover_one_target_prepared(
        prepared = payload$prepared,
        cells = payload$cells,
        target = payload$target,
        source_label = payload$source_label,
        source_type = payload$source_type,
        tf_cor = payload$tf_cor,
        peak_cor = payload$peak_cor
    )
}

.condition_discover_edges_compact <- function(
    prepared, cells, source_label, source_type, tf_cor, peak_cor,
    parallel = FALSE, verbose = TRUE) {
    cells <- intersect(as.character(cells), rownames(prepared$gene_data))
    if (length(cells) < 3L) {
        stop("Candidate discovery requires at least three paired cells.",
             call. = FALSE)
    }
    if (!is.numeric(tf_cor) || length(tf_cor) != 1L || !is.finite(tf_cor) ||
        tf_cor < 0 || tf_cor > 1 || !is.numeric(peak_cor) ||
        length(peak_cor) != 1L || !is.finite(peak_cor) ||
        peak_cor < 0 || peak_cor > 1) {
        stop("`tf_cor` and `peak_cor` must be finite values in [0, 1].",
             call. = FALSE)
    }

    features <- prepared$features
    names(features) <- features
    phase <- if (identical(as.character(source_type), "condition")) {
        "candidate_condition"
    } else {
        "candidate_global"
    }
    rows <- .pando_target_payload_map(
        keys = features,
        build_payload = function(target) {
            .pando_discovery_target_payload(
                prepared = prepared,
                cells = cells,
                target = target,
                source_label = source_label,
                source_type = source_type,
                tf_cor = tf_cor,
                peak_cor = peak_cor
            )
        },
        worker_name = ".pando_discovery_target_worker",
        parallel = parallel,
        verbose = verbose,
        phase = phase,
        label = source_label
    )
    rows <- rows[!vapply(rows, is.null, logical(1))]
    if (!length(rows)) {
        out <- data.frame(
            target = character(), tf = character(), region = character(),
            atac_feature_id = character(), peak_target_cor = numeric(),
            tf_target_cor = numeric(), source_label = character(),
            source_type = character(), edge_id = character(),
            stringsAsFactors = FALSE
        )
    } else {
        out <- unique(do.call(rbind, lapply(rows, as.data.frame)))
        out$edge_id <- paste(out$target, out$tf, out$region, sep = "||")
        out <- out[order(out$target, out$tf, out$region), , drop = FALSE]
        rownames(out) <- NULL
    }
    class(out) <- c("PandoEdgeDictionary", "data.frame")
    attr(out, "source_label") <- source_label
    attr(out, "source_type") <- source_type
    attr(out, "rna_layer") <- prepared$rna_layer
    attr(out, "peak_layer") <- prepared$peak_layer
    attr(out, "peak_value_type") <- prepared$peak_value_type
    attr(out, "preprocessing_fingerprint") <- prepared$preprocessing_fingerprint
    attr(out, "dictionary_input_schema") <-
        "pando_candidate_input_provenance_v1"
    out
}

.pando_ridge_target_payload <- function(
    prepared, edge_dictionary, target, cells, folds, control,
    min_residual_df, rank_action) {
    edges <- edge_dictionary[
        edge_dictionary$target == target, , drop = FALSE
    ]
    all_cells <- unique(unlist(cells, use.names = FALSE))
    compact <- list(
        gene_data = prepared$gene_data[
            all_cells, unique(c(target, as.character(edges$tf))), drop = FALSE
        ],
        peak_data = prepared$peak_data[
            all_cells, unique(as.character(edges$region)), drop = FALSE
        ]
    )
    list(
        prepared = compact,
        edges = edges,
        cells = cells,
        folds = folds,
        control = control,
        min_residual_df = min_residual_df,
        rank_action = rank_action
    )
}

.pando_ridge_target_worker <- function(payload) {
    .condition_ridge_target(
        prepared = payload$prepared,
        edges = payload$edges,
        cells_by_condition = payload$cells,
        folds = payload$folds,
        control = payload$control,
        min_residual_df = payload$min_residual_df,
        rank_action = payload$rank_action
    )
}
