# Independent broad-cell-type orchestration for condition-aware Pando.

.condition_bind_schema_rows <- function(old, new) {
    if (is.null(old) || !nrow(old)) return(new)
    if (is.null(new) || !nrow(new)) return(old)
    columns <- union(names(old), names(new))
    add_missing <- function(value) {
        missing <- setdiff(columns, names(value))
        for (name in missing) value[[name]] <- NA
        value[, columns, drop = FALSE]
    }
    rbind(add_missing(old), add_missing(new))
}

.condition_resolve_cell_types <- function(
    metadata, cell_type_col, cell_type = NULL
) {
    available_cell_types <- .condition_levels(metadata[[cell_type_col]])
    cell_types <- if (is.null(cell_type)) available_cell_types else
        trimws(as.character(cell_type))
    if (!length(cell_types) || anyNA(cell_types) || any(!nzchar(cell_types))) {
        stop('cell_type must be NULL or contain non-empty labels.')
    }
    if (anyDuplicated(cell_types)) {
        stop('cell_type labels must not be duplicated.')
    }
    missing_cell_types <- setdiff(cell_types, available_cell_types)
    if (length(missing_cell_types)) {
        stop(
            'Requested cell_type value(s) were not found: ',
            paste(missing_cell_types, collapse = ', '), '.'
        )
    }
    safe_cell_types <- vapply(
        cell_types, .condition_safe_id, character(1)
    )
    if (anyDuplicated(safe_cell_types)) {
        stop('cell_type labels are not unique after name sanitization.')
    }
    cell_types
}

.condition_fit_cell_type_models <- function(
    object,
    prepared,
    cell_type_col,
    cell_type,
    condition_col,
    network_name,
    peak_to_gene_method,
    upstream,
    downstream,
    extend,
    only_tss,
    tf_cor,
    peak_cor,
    candidate_screen,
    alpha,
    condition_mix,
    reference_condition,
    comparison_conditions,
    condition_weight,
    nlambda,
    lambda,
    lambda_min_ratio,
    outer_nfolds,
    inner_nfolds,
    lambda_selection,
    min_cells_per_condition,
    small_condition_action,
    active_tol,
    parallel,
    BPPARAM,
    overwrite,
    seed,
    max_iter,
    tol_objective,
    tol_coef,
    engine_control,
    verbose
) {
    cell_types <- .condition_resolve_cell_types(
        prepared$metadata, cell_type_col, cell_type
    )
    safe_cell_types <- vapply(
        cell_types, .condition_safe_id, character(1)
    )

    generated_index <- list()
    generated_diagnostics <- list()
    first_shared_id <- NULL
    for (cell_type_index in seq_along(cell_types)) {
        result <- .condition_fit_one_cell_type(
            object = object,
            prepared = prepared,
            cell_type = cell_types[[cell_type_index]],
            safe_cell_type = safe_cell_types[[cell_type_index]],
            cell_type_col = cell_type_col,
            condition_col = condition_col,
            network_name = network_name,
            peak_to_gene_method = peak_to_gene_method,
            upstream = upstream,
            downstream = downstream,
            extend = extend,
            only_tss = only_tss,
            tf_cor = tf_cor,
            peak_cor = peak_cor,
            candidate_screen = candidate_screen,
            alpha = alpha,
            condition_mix = condition_mix,
            reference_condition = reference_condition,
            comparison_conditions = comparison_conditions,
            condition_weight = condition_weight,
            nlambda = nlambda,
            lambda = lambda,
            lambda_min_ratio = lambda_min_ratio,
            outer_nfolds = outer_nfolds,
            inner_nfolds = inner_nfolds,
            lambda_selection = lambda_selection,
            min_cells_per_condition = min_cells_per_condition,
            small_condition_action = small_condition_action,
            active_tol = active_tol,
            parallel = parallel,
            BPPARAM = BPPARAM,
            overwrite = overwrite,
            seed = seed,
            max_iter = max_iter,
            tol_objective = tol_objective,
            tol_coef = tol_coef,
            engine_control = engine_control,
            verbose = verbose
        )
        object <- result$object
        if (!is.null(result$index)) {
            generated_index[[length(generated_index) + 1L]] <- result$index
        }
        if (!is.null(result$diagnostics)) {
            generated_diagnostics[[length(generated_diagnostics) + 1L]] <-
                result$diagnostics
        }
        if (is.null(first_shared_id) && !is.null(result$shared_id)) {
            first_shared_id <- result$shared_id
        }
    }

    if (length(generated_index)) {
        new_index <- do.call(rbind, generated_index)
        old_index <- object@grn@params$condition_network_index
        if (!is.null(old_index) && nrow(old_index)) {
            new_key <- paste(
                new_index$network_name, new_index$cell_type, sep = '\001'
            )
            old_key <- paste(
                old_index$network_name, old_index$cell_type, sep = '\001'
            )
            old_index <- old_index[
                !old_key %in% new_key, , drop = FALSE
            ]
            new_index <- rbind(old_index, new_index)
        }
        rownames(new_index) <- NULL
        object@grn@params$condition_network_index <- new_index
    }
    if (length(generated_diagnostics)) {
        diagnostics <- do.call(rbind, generated_diagnostics)
        old_diagnostics <- object@grn@params$condition_fit_diagnostics
        if (!is.null(old_diagnostics) && nrow(old_diagnostics)) {
            new_key <- paste(
                diagnostics$network_name,
                diagnostics$cell_type,
                sep = '\001'
            )
            old_key <- paste(
                old_diagnostics$network_name,
                old_diagnostics$cell_type,
                sep = '\001'
            )
            old_diagnostics <- old_diagnostics[
                !old_key %in% new_key, , drop = FALSE
            ]
            diagnostics <- .condition_bind_schema_rows(
                old_diagnostics, diagnostics
            )
        }
        rownames(diagnostics) <- NULL
        object@grn@params$condition_fit_diagnostics <- diagnostics
    }
    if (!is.null(first_shared_id)) {
        object@grn@active_network <- first_shared_id
    }
    object
}

.condition_fit_one_cell_type <- function(
    object,
    prepared,
    cell_type,
    safe_cell_type,
    cell_type_col,
    condition_col,
    network_name,
    peak_to_gene_method,
    upstream,
    downstream,
    extend,
    only_tss,
    tf_cor,
    peak_cor,
    candidate_screen,
    alpha,
    condition_mix,
    reference_condition,
    comparison_conditions,
    condition_weight,
    nlambda,
    lambda,
    lambda_min_ratio,
    outer_nfolds,
    inner_nfolds,
    lambda_selection,
    min_cells_per_condition,
    small_condition_action,
    active_tol,
    parallel,
    BPPARAM,
    overwrite,
    seed,
    max_iter,
    tol_objective,
    tol_coef,
    engine_control,
    verbose
) {
    cell_rows <- which(
        as.character(prepared$metadata[[cell_type_col]]) == cell_type
    )
    metadata <- prepared$metadata[cell_rows, , drop = FALSE]
    condition_levels <- .condition_levels(metadata[[condition_col]])
    condition_counts <- table(factor(
        as.character(metadata[[condition_col]]), levels = condition_levels
    ))
    small_conditions <- names(
        condition_counts[condition_counts < min_cells_per_condition]
    )
    if (length(small_conditions)) {
        message_text <- paste0(
            'Cell type ', cell_type,
            ' has condition(s) below min_cells_per_condition: ',
            paste(small_conditions, collapse = ', '), '.'
        )
        if (small_condition_action == 'error') {
            stop(message_text)
        }
        if (small_condition_action == 'skip_cell_type') {
            return(list(
                object = object,
                index = NULL,
                diagnostics = .condition_cell_type_skip_diagnostic(
                    network_name, cell_type, message_text
                ),
                shared_id = NULL
            ))
        }
        keep_conditions <- setdiff(condition_levels, small_conditions)
        keep <- as.character(metadata[[condition_col]]) %in%
            keep_conditions
        cell_rows <- cell_rows[keep]
        metadata <- metadata[keep, , drop = FALSE]
        condition_levels <- keep_conditions
    }
    if (length(condition_levels) < 2L) {
        message_text <- paste0(
            'Cell type ', cell_type,
            ' has fewer than two eligible conditions.'
        )
        if (small_condition_action == 'error') {
            stop(message_text)
        }
        return(list(
            object = object,
            index = NULL,
            diagnostics = .condition_cell_type_skip_diagnostic(
                network_name, cell_type, message_text
            ),
            shared_id = NULL
        ))
    }

    safe_conditions <- vapply(
        condition_levels, .condition_safe_id, character(1)
    )
    if (anyDuplicated(safe_conditions)) {
        stop(
            'Condition labels are not unique after name sanitization in ',
            'cell type ', cell_type, '.'
        )
    }
    reference_condition_cell <- if (is.null(reference_condition)) {
        condition_levels[[1L]]
    } else {
        as.character(reference_condition)
    }
    if (!reference_condition_cell %in% condition_levels) {
        stop(
            'reference_condition was not found in cell type ',
            cell_type, '.'
        )
    }
    comparison_conditions_cell <- if (is.null(comparison_conditions)) {
        condition_levels
    } else {
        as.character(comparison_conditions)
    }
    if (length(comparison_conditions_cell) < 2L ||
        !all(comparison_conditions_cell %in% condition_levels)) {
        stop(
            'comparison_conditions were not all found in cell type ',
            cell_type, '.'
        )
    }

    shared_id <- paste(
        network_name, safe_cell_type, 'shared', sep = '__'
    )
    condition_ids <- paste(
        network_name,
        safe_cell_type,
        'condition',
        safe_conditions,
        sep = '__'
    )
    conflicts <- intersect(
        c(shared_id, condition_ids), names(object@grn@networks)
    )
    if (length(conflicts) && !overwrite) {
        stop(
            'The following networks already exist: ',
            paste(conflicts, collapse = ', '),
            '. Set overwrite = TRUE to replace them.'
        )
    }

    log_message(
        'Fitting cell type ', cell_type, ' using only its own cells across ',
        length(condition_levels), ' conditions',
        verbose = verbose
    )
    condition_factor <- factor(
        as.character(metadata[[condition_col]]), levels = condition_levels
    )
    cell_engine_control <- engine_control
    if (!is.null(cell_engine_control$checkpoint_dir)) {
        cell_engine_control$checkpoint_dir <- file.path(
            cell_engine_control$checkpoint_dir, safe_cell_type
        )
    }
    target_results <- .condition_fit_targets(
        features = prepared$features,
        gene_data = prepared$gene_data[cell_rows, , drop = FALSE],
        peak_data = prepared$peak_data[cell_rows, , drop = FALSE],
        condition = condition_factor,
        peaks2gene = prepared$peaks2gene,
        peaks2motif = prepared$peaks2motif,
        motif2tf = prepared$motif2tf,
        candidate_screen = candidate_screen,
        tf_cor = tf_cor,
        peak_cor = peak_cor,
        scale = TRUE,
        alpha = alpha,
        condition_mix = condition_mix,
        active_tol = active_tol,
        reference_condition = reference_condition_cell,
        comparison_conditions = comparison_conditions_cell,
        condition_weight = condition_weight,
        nlambda = nlambda,
        lambda = lambda,
        lambda_min_ratio = lambda_min_ratio,
        outer_nfolds = outer_nfolds,
        inner_nfolds = inner_nfolds,
        lambda_selection = lambda_selection,
        seed = .condition_seed_for(cell_type, seed),
        max_iter = max_iter,
        tol_objective = tol_objective,
        tol_coef = tol_coef,
        engine_control = cell_engine_control,
        cell_type = cell_type,
        parallel = parallel,
        BPPARAM = BPPARAM,
        verbose = verbose
    )
    diagnostics <- target_results$diagnostics
    diagnostics$cell_type <- cell_type
    diagnostics$network_name <- network_name
    diagnostics <- diagnostics[, c(
        'network_name', 'cell_type', 'target', 'stage', 'converged',
        'iterations', 'objective', 'coef_change', 'selected_lambda',
        'cv_mean', 'cv_se', 'error_message', 'predictors', 'nonzeros',
        'path_backend', 'validation_backend', 'refit_backend',
        'pcg_iterations', 'pcg_residual', 'estimated_peak_bytes'
    ), drop = FALSE]
    successful <- target_results$fits
    if (!length(successful)) {
        return(list(
            object = object,
            index = NULL,
            diagnostics = diagnostics,
            shared_id = NULL
        ))
    }

    fit_engine <- 'condition_sparse_within_cell_type_oof_refit'
    fit_contract_key <- paste(
        network_name, safe_cell_type, sep = '__'
    )
    fit_contract <- .condition_combine_fit_contracts(
        successful = successful,
        network_name = network_name,
        cell_type = cell_type,
        cell_type_col = cell_type_col,
        condition_col = condition_col,
        reference_condition = reference_condition_cell,
        comparison_conditions = comparison_conditions_cell,
        candidate_screen = candidate_screen,
        scale = TRUE,
        fit_engine = fit_engine
    )
    fit_contract$cell_ids <- rownames(metadata)
    fit_contract$cell_condition <- stats::setNames(
        as.character(metadata[[condition_col]]), rownames(metadata)
    )
    fit_contract$predictive_oof_complete <-
        all(fit_contract$predictive_oof_available)
    fit_contract$cell_provenance <- data.frame(
        cell_id = rownames(metadata),
        condition = as.character(metadata[[condition_col]]),
        cell_type = cell_type,
        stringsAsFactors = FALSE
    )
    object_params <- Params(object)
    fit_contract$assay_contract <- list(
        rna_assay = object_params$rna_assay,
        peak_assay = object_params$peak_assay,
        rna_layer = 'data',
        peak_layer = 'data'
    )
    fit_contract$projection_contract$cell_scope <-
        'exact paired cells of one broad cell type used for fitting'
    fit_contracts <- object@grn@params$condition_grn_fits
    if (is.null(fit_contracts)) {
        fit_contracts <- list()
    }
    fit_contracts[[fit_contract_key]] <- fit_contract
    object@grn@params$condition_grn_fits <- fit_contracts

    common_params <- list(
        cell_type = cell_type,
        cell_type_col = cell_type_col,
        condition_col = condition_col,
        condition_levels = condition_levels,
        candidate_screen = candidate_screen,
        fit_engine = fit_engine,
        condition_weight = condition_weight,
        alpha = alpha,
        condition_mix = condition_mix,
        reference_condition = reference_condition_cell,
        comparison_conditions = comparison_conditions_cell,
        fit_contract_key = fit_contract_key,
        lambda_selection = lambda_selection,
        nlambda = nlambda,
        outer_nfolds = outer_nfolds,
        inner_nfolds = inner_nfolds,
        oof_scheme = 'nested_outer_condition_stratified_cell_oof',
        scale = TRUE,
        active_tol = active_tol,
        seed = seed,
        upstream = upstream,
        downstream = downstream,
        extend = extend,
        only_tss = only_tss,
        peak_to_gene_method = peak_to_gene_method,
        tf_cor = tf_cor,
        peak_cor = peak_cor,
        engine_control = engine_control
    )
    successful_features <- names(successful)
    shared_coefs <- do.call(
        rbind, lapply(successful, function(x) x$universal_coefs)
    )
    shared_gof <- do.call(
        rbind, lapply(successful, function(x) x$universal_gof)
    )
    shared_params <- do.call(
        .condition_network_params,
        c(
            list(network_level = 'shared', condition = NA_character_),
            common_params
        )
    )
    object@grn@networks[[shared_id]] <- .condition_build_network(
        successful_features, shared_coefs, shared_gof, shared_params
    )
    index_rows <- list(.condition_index_row(
        network_id = shared_id,
        network_name = network_name,
        cell_type = cell_type,
        network_level = 'shared',
        condition = NA_character_,
        n_cells = nrow(metadata),
        coefs = shared_coefs,
        n_targets = length(successful_features),
        active_tol = active_tol,
        fit_engine = fit_engine,
        reference_condition = reference_condition_cell,
        fit_contract_key = fit_contract_key
    ))

    for (condition_index in seq_along(condition_levels)) {
        condition_name <- condition_levels[[condition_index]]
        condition_id <- condition_ids[[condition_index]]
        condition_coefs <- do.call(rbind, lapply(
            successful, function(x) x$condition_coefs[[condition_name]]
        ))
        condition_gof <- do.call(rbind, lapply(
            successful, function(x) x$condition_gof[[condition_name]]
        ))
        condition_params <- do.call(
            .condition_network_params,
            c(
                list(
                    network_level = 'condition',
                    condition = condition_name
                ),
                common_params
            )
        )
        object@grn@networks[[condition_id]] <- .condition_build_network(
            successful_features,
            condition_coefs,
            condition_gof,
            condition_params
        )
        index_rows[[length(index_rows) + 1L]] <- .condition_index_row(
            network_id = condition_id,
            network_name = network_name,
            cell_type = cell_type,
            network_level = 'condition',
            condition = condition_name,
            n_cells = sum(condition_factor == condition_name),
            coefs = condition_coefs,
            n_targets = length(successful_features),
            active_tol = active_tol,
            fit_engine = fit_engine,
            reference_condition = reference_condition_cell,
            fit_contract_key = fit_contract_key
        )
    }
    list(
        object = object,
        index = do.call(rbind, index_rows),
        diagnostics = diagnostics,
        shared_id = shared_id
    )
}

.condition_cell_type_skip_diagnostic <- function(
    network_name, cell_type, message_text
) {
    data.frame(
        network_name = network_name,
        cell_type = cell_type,
        target = NA_character_,
        stage = 'cell_type_validation',
        converged = FALSE,
        iterations = NA_integer_,
        objective = NA_real_,
        coef_change = NA_real_,
        selected_lambda = NA_real_,
        cv_mean = NA_real_,
        cv_se = NA_real_,
        error_message = message_text,
        predictors = NA_integer_,
        nonzeros = NA_real_,
        path_backend = NA_character_,
        validation_backend = NA_character_,
        refit_backend = NA_character_,
        pcg_iterations = NA_integer_,
        pcg_residual = NA_real_,
        estimated_peak_bytes = NA_real_,
        stringsAsFactors = FALSE
    )
}
