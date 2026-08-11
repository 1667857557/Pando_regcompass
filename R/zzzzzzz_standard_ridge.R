# Optional single-task ridge for standard Pando.
#
# This is the K = 1 special case of the condition-GRN ridge solver. The fusion
# Laplacian is exactly zero for K = 1, so the numerical path, predictor scaling,
# CV lambda selection, effective degrees of freedom and ridge-Wald diagnostics
# are identical to the condition implementation without introducing a second
# ridge backend. Original Gaussian GLM remains the default standard-Pando path.

.pando_standard_ridge_infer_impl <- infer_grn.GRNData

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

    candidates <- .condition_discover_edges_prepared(
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
        verbose = verbose
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

infer_grn.GRNData <- function(
    object,
    genes = NULL,
    network_name = paste0(method, "_network"),
    peak_to_gene_method = c("Signac", "GREAT"),
    upstream = 100000,
    downstream = 0,
    extend = 1000000,
    only_tss = FALSE,
    parallel = FALSE,
    tf_cor = 0.1,
    peak_cor = 0,
    aggregate_rna_col = NULL,
    aggregate_peaks_col = NULL,
    method = c("glm", "ridge", "glmnet", "cv.glmnet", "brms", "xgb",
               "bagging_ridge", "bayesian_ridge"),
    alpha = 0.5,
    family = "gaussian",
    interaction_term = ":",
    adjust_method = "fdr",
    scale = FALSE,
    verbose = TRUE,
    BPPARAM = NULL,
    ridge_control = list(),
    rank_action = c("mark", "error"),
    min_residual_df = 1L,
    padj_threshold = 0.05,
    ...) {
    method <- match.arg(method)
    peak_to_gene_method <- match.arg(peak_to_gene_method)
    if (!identical(method, "ridge")) {
        ridge_requested <- length(ridge_control) ||
            !identical(match.arg(rank_action), "mark") ||
            !identical(as.integer(min_residual_df), 1L) ||
            !isTRUE(all.equal(as.numeric(padj_threshold), 0.05)) ||
            (!is.null(BPPARAM) && !identical(BPPARAM, FALSE))
        if (ridge_requested) {
            stop(
                "`BPPARAM`, `ridge_control`, `rank_action`, `min_residual_df` ",
                "and `padj_threshold` are standard-ridge controls; set ",
                "`method = \"ridge\"` or remove them.", call. = FALSE
            )
        }
        return(.pando_standard_ridge_infer_impl(
            object = object,
            genes = genes,
            network_name = network_name,
            peak_to_gene_method = peak_to_gene_method,
            upstream = upstream,
            downstream = downstream,
            extend = extend,
            only_tss = only_tss,
            parallel = parallel,
            tf_cor = tf_cor,
            peak_cor = peak_cor,
            aggregate_rna_col = aggregate_rna_col,
            aggregate_peaks_col = aggregate_peaks_col,
            method = method,
            alpha = alpha,
            family = family,
            interaction_term = interaction_term,
            adjust_method = adjust_method,
            scale = scale,
            verbose = verbose,
            ...
        ))
    }

    if (!is.null(aggregate_rna_col) || !is.null(aggregate_peaks_col)) {
        stop(
            "Standard `method = \"ridge\"` uses paired single-cell RNA and ATAC; ",
            "aggregate_rna_col/aggregate_peaks_col are not supported.",
            call. = FALSE
        )
    }
    if (!.pando_standard_ridge_family_ok(family) ||
        !identical(interaction_term, ":") || !identical(scale, FALSE)) {
        stop(
            "Standard ridge requires Gaussian identity, interaction_term=':', ",
            "and scale=FALSE.", call. = FALSE
        )
    }
    .pando_validate_bpparam(BPPARAM)
    old_target_param <- getOption("Pando.condition_target_BPPARAM", NULL)
    started_here <- FALSE
    if (isTRUE(parallel) && !is.null(BPPARAM) &&
        !identical(BPPARAM, FALSE)) {
        if (!isTRUE(BiocParallel::bpisup(BPPARAM))) {
            BPPARAM <- BiocParallel::bpstart(BPPARAM)
            started_here <- TRUE
        }
        options(Pando.condition_target_BPPARAM = BPPARAM)
    }
    on.exit({
        options(Pando.condition_target_BPPARAM = old_target_param)
        if (started_here) try(BiocParallel::bpstop(BPPARAM), silent = TRUE)
        invisible(gc(verbose = FALSE, full = TRUE))
    }, add = TRUE)

    .pando_standard_ridge_fit(
        object = object,
        genes = genes,
        network_name = network_name,
        peak_to_gene_method = peak_to_gene_method,
        upstream = upstream,
        downstream = downstream,
        extend = extend,
        only_tss = only_tss,
        tf_cor = tf_cor,
        peak_cor = peak_cor,
        adjust_method = adjust_method,
        padj_threshold = padj_threshold,
        rank_action = rank_action,
        min_residual_df = min_residual_df,
        ridge_control = ridge_control,
        parallel = parallel,
        verbose = verbose
    )
}
