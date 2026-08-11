# Canonical single-task ridge helpers for standard Pando.
#
# Standard ridge is the K = 1 specialization of the condition multi-task ridge
# solver. Public method routing remains in the original infer_grn.GRNData
# definition in grn.R; this file contains only numerical/helper implementation.

.pando_standard_ridge_family_ok <- function(family) {
    if (is.character(family) && length(family) == 1L) {
        return(tolower(family) == "gaussian")
    }
    identical(tryCatch(family$family, error = function(e) NULL), "gaussian") &&
        identical(tryCatch(family$link, error = function(e) NULL), "identity")
}

.pando_standard_ridge_fit <- function(
    object, genes, network_name,
    peak_to_gene_method, upstream, downstream, extend, only_tss,
    tf_cor, peak_cor, adjust_method, padj_threshold,
    rank_action, min_residual_df, ridge_control,
    parallel, verbose) {
    control <- .condition_ridge_control(ridge_control)
    rank_action <- match.arg(rank_action, c("mark", "error"))
    if (!is.numeric(min_residual_df) || length(min_residual_df) != 1L ||
        !is.finite(min_residual_df) || min_residual_df < 1L) {
        stop("`min_residual_df` must be one positive number.", call. = FALSE)
    }
    if (!is.numeric(padj_threshold) || length(padj_threshold) != 1L ||
        !is.finite(padj_threshold) || padj_threshold <= 0 ||
        padj_threshold >= 1) {
        stop("`padj_threshold` must be one finite number in (0, 1).",
             call. = FALSE)
    }

    prepared <- .condition_prepare_common_input(
        object = object,
        genes = genes,
        peak_to_gene_method = peak_to_gene_method,
        upstream = upstream,
        downstream = downstream,
        extend = extend,
        only_tss = only_tss,
        rna_layer = "data",
        peak_layer = "data",
        peak_value_type = "normalized",
        verbose = verbose
    )
    cells <- rownames(prepared$gene_data)
    if (length(cells) < 3L) {
        stop("Standard ridge requires at least three paired cells.",
             call. = FALSE)
    }

    candidates <- .condition_discover_edges_compact(
        prepared = prepared,
        cells = cells,
        source_label = "standard",
        source_type = "global",
        tf_cor = tf_cor,
        peak_cor = peak_cor,
        parallel = parallel,
        verbose = verbose
    )
    if (!is.data.frame(candidates) || !nrow(candidates)) {
        stop("Standard ridge candidate discovery produced no TF-peak-target edges.",
             call. = FALSE)
    }
    dictionary <- union_grn_edges(
        global_edges = NULL,
        condition_edges = list(standard = candidates)
    )
    network_names <- stats::setNames(network_name, "standard")
    skeleton <- list(
        schema_version = .condition_common_dictionary_schema,
        model_schema = .condition_multitask_ridge_schema,
        fit_engine = "single_task_ridge_same_condition_solver",
        coefficient_scale = "raw_tf_atac_interaction_units",
        internal_predictor_scale = "equal_condition_within_condition_rms",
        inference_scope =
            "approximate_ridge_wald_diagnostic_conditional_on_candidates_and_cv_lambda",
        cell_type = NA_character_,
        condition_levels = "standard",
        condition_col = NULL,
        cell_type_col = NULL,
        condition_cell_ids = list(standard = cells),
        edge_dictionary = dictionary,
        coefficients = NULL,
        contrasts = NULL,
        fit = NULL,
        network_names = network_names,
        padj_threshold = as.numeric(padj_threshold),
        adjust_method = as.character(adjust_method),
        scale = FALSE,
        interaction = ":",
        projection_effect_column = "penalty_effect",
        projection_policy = "standard_ridge_bh_diagnostic",
        target_genes = unique(as.character(dictionary$target)),
        rna_assay = prepared$params$rna_assay,
        atac_assay = prepared$params$peak_assay,
        rna_layer = prepared$rna_layer,
        peak_layer = prepared$peak_layer,
        peak_value_type = prepared$peak_value_type,
        preprocessing_fingerprint = prepared$preprocessing_fingerprint,
        ridge_control = control
    )
    class(skeleton) <- c("ConditionGRNFit", "list")

    refitted <- .condition_ridge_refit_contract_one_pass(
        object = object,
        fit = skeleton,
        prepared = prepared,
        control = control,
        rank_action = rank_action,
        min_residual_df = min_residual_df,
        parallel = parallel,
        verbose = verbose,
        progress_phase = "ridge_standard",
        progress_label = network_name
    )
    answer <- refitted$object
    fit <- refitted$fit

    network <- answer@grn@networks[[network_name]]
    if (!is.null(network)) {
        network@params$method <- "ridge"
        network@params$fit_mode <- "single_task_ridge"
        network@params$condition <- NULL
        network@params$projection_policy <- "standard_ridge_bh_diagnostic"
        network@params$ridge_solver <- "condition_ridge_k1"
        network@params$ridge_control <- control
        answer@grn@networks[[network_name]] <- network
        answer@grn@active_network <- network_name
    }
    answer@grn@params$analysis_mode <- "standard_grn"
    answer@grn@params$standard_fit_method <- "ridge"
    answer@grn@params$standard_ridge_schema <-
        "pando_standard_grn_single_task_ridge_v1"
    answer@grn@params$standard_ridge_contract <- list(
        schema = "pando_standard_grn_single_task_ridge_v1",
        solver = "condition_ridge_k1",
        candidate_edges = nrow(dictionary),
        fitted_targets = length(unique(as.character(fit$fit$target))),
        coefficient_scale = "raw_tf_atac_interaction_units",
        predictor_scale_reference = "within_task_rms",
        adjust_method = as.character(adjust_method),
        padj_threshold = as.numeric(padj_threshold),
        rank_action = rank_action,
        min_residual_df = min_residual_df,
        ridge_control = control
    )
    answer
}
