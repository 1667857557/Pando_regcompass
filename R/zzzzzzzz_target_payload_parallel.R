# Memory-bounded target-level parallelism for candidate discovery and ridge refits.
#
# Target workers receive only the RNA, ATAC, motif and dictionary slices required
# by that target. Payloads are dispatched in worker-sized batches so the master
# does not materialize one compact copy for every target at once. Worker-local
# temporaries and completed master-side batches are released immediately.

.pando_full_discover_edges_prepared <- .condition_discover_edges_prepared

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

.pando_target_payload_map <- function(
    keys, build_payload, worker, parallel = FALSE, verbose = TRUE) {
    if (!length(keys)) return(list())
    if (!is.function(build_payload) || !is.function(worker)) {
        stop("Target payload mapping requires builder and worker functions.",
             call. = FALSE)
    }
    key_names <- names(keys)
    if (is.null(key_names) || any(!nzchar(key_names))) {
        key_names <- as.character(keys)
    }
    out <- vector("list", length(keys))
    names(out) <- key_names
    batch_size <- min(length(keys), .pando_target_worker_limit(parallel))
    starts <- seq.int(1L, length(keys), by = batch_size)
    for (start in starts) {
        index <- seq.int(start, min(length(keys), start + batch_size - 1L))
        payloads <- lapply(index, function(i) build_payload(keys[[i]]))
        names(payloads) <- key_names[index]
        chunk <- if (isTRUE(parallel)) {
            map_par(payloads, worker, parallel = TRUE, verbose = FALSE)
        } else {
            lapply(payloads, worker)
        }
        out[index] <- chunk
        payloads <- NULL
        chunk <- NULL
        invisible(gc(verbose = FALSE, full = TRUE))
    }
    out
}

.pando_discovery_target_payload <- function(
    prepared, cells, target, source_label, source_type,
    tf_cor, peak_cor, verbose) {
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
        peaks2gene = prepared$peaks2gene[
            target, regions, drop = FALSE
        ],
        peaks2motif = prepared$peaks2motif[
            regions, motif_ids, drop = FALSE
        ],
        motif2tf = prepared$motif2tf[
            motif_ids, tfs, drop = FALSE
        ],
        region_map = region_map,
        rna_layer = prepared$rna_layer,
        peak_layer = prepared$peak_layer,
        peak_value_type = prepared$peak_value_type,
        preprocessing_fingerprint = prepared$preprocessing_fingerprint
    )
    list(
        skip = FALSE,
        prepared = compact,
        cells = rownames(compact$gene_data),
        source_label = source_label,
        source_type = source_type,
        tf_cor = tf_cor,
        peak_cor = peak_cor,
        verbose = verbose
    )
}

.pando_discovery_target_worker <- function(payload) {
    on.exit(invisible(gc(verbose = FALSE, full = TRUE)), add = TRUE)
    if (isTRUE(payload$skip)) return(NULL)
    .pando_full_discover_edges_prepared(
        prepared = payload$prepared,
        cells = payload$cells,
        source_label = payload$source_label,
        source_type = payload$source_type,
        tf_cor = payload$tf_cor,
        peak_cor = payload$peak_cor,
        parallel = FALSE,
        verbose = payload$verbose
    )
}

.condition_discover_edges_prepared <- function(
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
    rows <- .pando_target_payload_map(
        keys = features,
        build_payload = function(target) {
            .pando_discovery_target_payload(
                prepared = prepared, cells = cells, target = target,
                source_label = source_label, source_type = source_type,
                tf_cor = tf_cor, peak_cor = peak_cor, verbose = verbose
            )
        },
        worker = .pando_discovery_target_worker,
        parallel = parallel,
        verbose = verbose
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
        rownames(out) <- NULL
        out <- out[order(out$target, out$tf, out$region), , drop = FALSE]
    }
    class(out) <- c("PandoEdgeDictionary", "data.frame")
    attr(out, "source_label") <- source_label
    attr(out, "source_type") <- source_type
    attr(out, "rna_layer") <- prepared$rna_layer
    attr(out, "peak_layer") <- prepared$peak_layer
    attr(out, "peak_value_type") <- prepared$peak_value_type
    attr(out, "preprocessing_fingerprint") <-
        prepared$preprocessing_fingerprint
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
    on.exit(invisible(gc(verbose = FALSE, full = TRUE)), add = TRUE)
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

.condition_ridge_refit_contract_compact <- function(
    object, fit, prepared, control, rank_action = "mark",
    min_residual_df = 1L, parallel = FALSE, verbose = TRUE) {
    .condition_validate_dictionary(fit$edge_dictionary, prepared)
    cells <- fit$condition_cell_ids[fit$condition_levels]
    if (any(lengths(cells) < 3L)) {
        stop("Every multi-task ridge condition needs at least three cells.",
             call. = FALSE)
    }
    folds <- .condition_ridge_folds(cells, control$cv_folds, control$seed)
    targets <- unique(as.character(fit$edge_dictionary$target))
    names(targets) <- targets
    result <- .pando_target_payload_map(
        keys = targets,
        build_payload = function(target) {
            .pando_ridge_target_payload(
                prepared = prepared,
                edge_dictionary = fit$edge_dictionary,
                target = target,
                cells = cells,
                folds = folds,
                control = control,
                min_residual_df = min_residual_df,
                rank_action = rank_action
            )
        },
        worker = .pando_ridge_target_worker,
        parallel = parallel,
        verbose = verbose
    )

    coefficient <- do.call(rbind, lapply(result, `[[`, "coefficients"))
    contrast <- do.call(rbind, lapply(result, `[[`, "contrasts"))
    fit_table <- do.call(rbind, lapply(result, `[[`, "fit"))
    rownames(coefficient) <- rownames(fit_table) <- NULL
    if (is.data.frame(contrast) && nrow(contrast)) rownames(contrast) <- NULL

    coefficient$padj <- NA_real_
    for (condition in fit$condition_levels) {
        index <- which(coefficient$condition == condition)
        valid <- index[coefficient$estimable[index] %in% TRUE &
                       is.finite(coefficient$pval[index])]
        if (length(valid)) {
            coefficient$padj[valid] <- stats::p.adjust(
                coefficient$pval[valid], method = fit$adjust_method
            )
        }
    }
    coefficient$significant <- coefficient$estimable &
        is.finite(coefficient$padj) &
        coefficient$padj < fit$padj_threshold
    coefficient$penalty_effect <- ifelse(
        coefficient$estimable & is.finite(coefficient$estimate),
        coefficient$estimate,
        0
    )
    coefficient$direction <- ifelse(
        !coefficient$estimable, "undefined",
        ifelse(coefficient$estimate > 0, "positive",
               ifelse(coefficient$estimate < 0, "negative", "zero"))
    )
    coefficient$effect_definition <-
        "multitask_ridge_condition_coefficient_raw_tf_atac_units"
    coefficient$inference_scope <-
        "approximate_ridge_wald_diagnostic_conditional_on_dictionary_cv_lambda_and_fusion"

    if (is.data.frame(contrast) && nrow(contrast)) {
        contrast$contrast_padj <- NA_real_
        pair_key <- paste(contrast$condition_a, contrast$condition_b, sep = "\001")
        for (key in unique(pair_key)) {
            index <- which(pair_key == key)
            valid <- index[contrast$contrast_estimable[index] %in% TRUE &
                           is.finite(contrast$contrast_pval[index])]
            if (length(valid)) {
                contrast$contrast_padj[valid] <- stats::p.adjust(
                    contrast$contrast_pval[valid], method = fit$adjust_method
                )
            }
        }
        contrast$contrast_significant <- contrast$contrast_estimable &
            is.finite(contrast$contrast_padj) &
            contrast$contrast_padj < fit$padj_threshold
        contrast$inference_scope <-
            "approximate_joint_ridge_wald_contrast_diagnostic"
    }

    for (condition in fit$condition_levels) {
        network_name <- fit$network_names[[condition]]
        coefs_one <- coefficient[
            coefficient$condition == condition, , drop = FALSE
        ]
        fit_one <- fit_table[
            fit_table$condition == condition, , drop = FALSE
        ]
        network <- methods::new(
            Class = "Network",
            features = unique(as.character(coefs_one$target)),
            coefs = coefs_one,
            fit = fit_one,
            params = list(
                method = "multitask_ridge",
                family = "gaussian_identity",
                fit_mode = "fixed_edge_dictionary_joint_conditions",
                condition = condition,
                edge_dictionary = fit$edge_dictionary,
                scale = FALSE,
                internal_scale_reference =
                    "equal_condition_within_condition_rms",
                exported_coefficient_scale = "raw_tf_atac_interaction_units",
                interaction = ":",
                rna_layer = prepared$rna_layer,
                peak_layer = prepared$peak_layer,
                peak_value_type = prepared$peak_value_type,
                preprocessing_fingerprint = prepared$preprocessing_fingerprint,
                adjust_method = fit$adjust_method,
                padj_threshold = fit$padj_threshold,
                projection_policy = "continuous_estimable_ridge_effects",
                ridge_control = control
            )
        )
        object@grn@networks[[network_name]] <- network
        object@grn@active_network <- network_name
    }

    fit$model_schema <- .condition_multitask_ridge_schema
    fit$fit_engine <- "two_stage_exact_edge_union_multitask_ridge"
    fit$coefficient_scale <- "raw_tf_atac_interaction_units"
    fit$internal_predictor_scale <- "equal_condition_within_condition_rms"
    fit$inference_scope <-
        "approximate_ridge_wald_diagnostic_conditional_on_dictionary_cv_lambda_and_fusion"
    fit$coefficients <- coefficient
    fit$contrasts <- contrast
    fit$fit <- fit_table
    fit$scale <- FALSE
    fit$projection_effect_column <- "penalty_effect"
    fit$projection_policy <- "continuous_estimable_ridge_effects"
    fit$ridge_control <- control
    fit$target_cv <- lapply(result, function(one) {
        one$cv[c("lambda", "lambda_min", "cv_mse", "cv_se",
                 "rsq_oof", "curve")]
    })
    fit$target_scaling <- lapply(result, `[[`, "scaling")
    fit$target_parallel_memory_policy <- list(
        payload = "target_specific_rna_atac_edges",
        batching = "worker_sized",
        worker_gc = TRUE,
        master_batch_gc = TRUE
    )
    class(fit) <- c("ConditionGRNFit", "list")
    result <- NULL
    invisible(gc(verbose = FALSE, full = TRUE))
    list(object = object, fit = fit)
}

# The significant-union wrapper resolves this symbol at runtime. Replacing the
# one-pass estimator therefore applies the compact payload policy to both the
# preliminary screen and final refit, and standard ridge reuses it directly.
.condition_ridge_refit_contract_one_pass <- .condition_ridge_refit_contract_compact
