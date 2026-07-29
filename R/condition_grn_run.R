# Cell-type orchestration for condition-aware Pando.

.condition_run_all_cell_types <- function(
    object,
    prepared,
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
    condition_weight,
    nlambda,
    lambda,
    lambda_min_ratio,
    nfolds,
    cv_block_col,
    lambda_selection,
    min_cells_per_condition,
    on_small_condition,
    active_tol,
    parallel,
    BPPARAM,
    overwrite,
    seed,
    max_iter,
    tol_objective,
    tol_coef,
    verbose
) {
    cell_types <- .condition_levels(prepared$metadata[[cell_type_col]])
    safe_cell_types <- vapply(cell_types, .condition_safe_id, character(1))
    if (anyDuplicated(safe_cell_types)) {
        stop('cell_type labels are not unique after network-name sanitization.')
    }

    generated_index <- list()
    generated_diagnostics <- list()
    first_universal <- NULL

    for (cell_type_index in seq_along(cell_types)) {
        result <- .condition_run_cell_type(
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
            condition_weight = condition_weight,
            nlambda = nlambda,
            lambda = lambda,
            lambda_min_ratio = lambda_min_ratio,
            nfolds = nfolds,
            cv_block_col = cv_block_col,
            lambda_selection = lambda_selection,
            min_cells_per_condition = min_cells_per_condition,
            on_small_condition = on_small_condition,
            active_tol = active_tol,
            parallel = parallel,
            BPPARAM = BPPARAM,
            overwrite = overwrite,
            seed = seed,
            max_iter = max_iter,
            tol_objective = tol_objective,
            tol_coef = tol_coef,
            verbose = verbose
        )
        object <- result$object
        if (!is.null(result$index)) {
            generated_index[[length(generated_index) + 1L]] <- result$index
        }
        if (!is.null(result$diagnostics)) {
            generated_diagnostics[[length(generated_diagnostics) + 1L]] <- result$diagnostics
        }
        if (is.null(first_universal) && !is.null(result$universal_id)) {
            first_universal <- result$universal_id
        }
    }

    if (length(generated_index) > 0L) {
        new_index <- do.call(rbind, generated_index)
        old_index <- object@grn@params$condition_network_index
        if (!is.null(old_index) && nrow(old_index) > 0L) {
            old_index <- old_index[!old_index$network_id %in% new_index$network_id, , drop = FALSE]
            new_index <- rbind(old_index, new_index)
        }
        rownames(new_index) <- NULL
        object@grn@params$condition_network_index <- new_index
    }

    if (length(generated_diagnostics) > 0L) {
        diagnostics <- do.call(rbind, generated_diagnostics)
        old_diagnostics <- object@grn@params$condition_fit_diagnostics
        if (!is.null(old_diagnostics) && nrow(old_diagnostics) > 0L &&
            all(c('network_name', 'cell_type') %in% colnames(old_diagnostics))) {
            replaced_key <- paste(diagnostics$network_name, diagnostics$cell_type, sep = '\r')
            old_key <- paste(old_diagnostics$network_name, old_diagnostics$cell_type, sep = '\r')
            old_diagnostics <- old_diagnostics[!old_key %in% replaced_key, , drop = FALSE]
            diagnostics <- rbind(old_diagnostics, diagnostics)
        }
        rownames(diagnostics) <- NULL
        object@grn@params$condition_fit_diagnostics <- diagnostics
    }

    if (!is.null(first_universal)) {
        object@grn@active_network <- first_universal
    }
    object
}

.condition_run_cell_type <- function(
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
    condition_weight,
    nlambda,
    lambda,
    lambda_min_ratio,
    nfolds,
    cv_block_col,
    lambda_selection,
    min_cells_per_condition,
    on_small_condition,
    active_tol,
    parallel,
    BPPARAM,
    overwrite,
    seed,
    max_iter,
    tol_objective,
    tol_coef,
    verbose
) {
    cell_rows <- which(as.character(prepared$metadata[[cell_type_col]]) == cell_type)
    metadata_cell_type <- prepared$metadata[cell_rows, , drop = FALSE]
    condition_levels <- .condition_levels(metadata_cell_type[[condition_col]])
    condition_counts <- table(factor(
        as.character(metadata_cell_type[[condition_col]]), levels = condition_levels
    ))
    small_conditions <- names(condition_counts)[condition_counts < min_cells_per_condition]

    if (length(small_conditions) > 0L) {
        message_text <- paste0(
            'Cell type ', cell_type, ' has condition(s) below min_cells_per_condition: ',
            paste(small_conditions, collapse = ', '), '.'
        )
        if (on_small_condition == 'error') {
            stop(message_text)
        }
        if (on_small_condition == 'skip_cell_type') {
            log_message('Skipping ', message_text, verbose = verbose)
            return(list(
                object = object,
                index = NULL,
                diagnostics = .condition_skip_diagnostic(network_name, cell_type, message_text),
                universal_id = NULL
            ))
        }
        keep_conditions <- setdiff(condition_levels, small_conditions)
        keep_rows <- as.character(metadata_cell_type[[condition_col]]) %in% keep_conditions
        cell_rows <- cell_rows[keep_rows]
        metadata_cell_type <- prepared$metadata[cell_rows, , drop = FALSE]
        condition_levels <- keep_conditions
    }

    if (length(condition_levels) < 2L) {
        message_text <- paste0('Cell type ', cell_type, ' has fewer than two eligible conditions.')
        if (on_small_condition == 'error') {
            stop(message_text)
        }
        log_message('Skipping ', message_text, verbose = verbose)
        return(list(
            object = object,
            index = NULL,
            diagnostics = .condition_skip_diagnostic(network_name, cell_type, message_text),
            universal_id = NULL
        ))
    }

    safe_conditions <- vapply(condition_levels, .condition_safe_id, character(1))
    if (anyDuplicated(safe_conditions)) {
        stop('Condition labels are not unique after network-name sanitization in cell type ', cell_type, '.')
    }
    reference_condition_cell <- if (is.null(reference_condition)) {
        condition_levels[[1L]]
    } else {
        as.character(reference_condition)
    }
    if (!reference_condition_cell %in% condition_levels) {
        stop(
            'Reference condition ', reference_condition_cell,
            ' was not found in cell type ', cell_type, '.'
        )
    }
    universal_id <- paste(network_name, safe_cell_type, 'shared', sep = '__')
    condition_ids <- paste(network_name, safe_cell_type, 'condition', safe_conditions, sep = '__')
    intended_ids <- c(universal_id, condition_ids)
    conflicts <- intersect(intended_ids, names(object@grn@networks))
    if (length(conflicts) > 0L && !overwrite) {
        stop(
            'The following networks already exist: ', paste(conflicts, collapse = ', '),
            '. Set overwrite = TRUE to replace them.'
        )
    }

    log_message(
        'Fitting condition-aware GRN for cell type ', cell_type,
        ' across ', length(condition_levels), ' conditions', verbose = verbose
    )
    cv_block_status <- if (is.null(cv_block_col)) {
        NULL
    } else {
        block_list <- split(
            trimws(as.character(metadata_cell_type[[cv_block_col]])),
            factor(
                as.character(metadata_cell_type[[condition_col]]),
                levels = condition_levels
            )
        )
        status <- .condition_sample_block_status(block_list)
        data.frame(
            condition = names(status$n_blocks),
            n_cells = as.integer(table(factor(
                as.character(metadata_cell_type[[condition_col]]),
                levels = condition_levels
            ))),
            n_biological_samples = as.integer(status$n_blocks),
            sample_blocked_oof_available = status$n_blocks >= 2L,
            stringsAsFactors = FALSE
        )
    }
    if (!is.null(cv_block_status) &&
        any(!cv_block_status$sample_blocked_oof_available)) {
        warning(
            'Cell type ', cell_type, ' has fewer than two biological samples ',
            'in condition(s): ',
            paste(
                cv_block_status$condition[
                    !cv_block_status$sample_blocked_oof_available
                ],
                collapse = ', '
            ),
            '. Cell-level folds will select lambda, but sample-blocked OOF ',
            'performance is unavailable and must not be interpreted as ',
            'replicate-level validation.',
            call. = FALSE
        )
    }
    fit_cell_type <- .condition_fit_cell_type(
        features = prepared$features,
        gene_data = prepared$gene_data[cell_rows, , drop = FALSE],
        peak_data = prepared$peak_data[cell_rows, , drop = FALSE],
        condition = factor(
            as.character(metadata_cell_type[[condition_col]]), levels = condition_levels
        ),
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
        condition_weight = condition_weight,
        nlambda = nlambda,
        lambda = lambda,
        lambda_min_ratio = lambda_min_ratio,
        nfolds = nfolds,
        cv_block = if (is.null(cv_block_col)) {
            NULL
        } else {
            stats::setNames(
                trimws(as.character(metadata_cell_type[[cv_block_col]])),
                rownames(metadata_cell_type)
            )
        },
        lambda_selection = lambda_selection,
        seed = .condition_seed_for(cell_type, seed),
        max_iter = max_iter,
        tol_objective = tol_objective,
        tol_coef = tol_coef,
        parallel = parallel,
        BPPARAM = BPPARAM,
        verbose = verbose
    )

    diagnostics <- fit_cell_type$diagnostics
    diagnostics$cell_type <- cell_type
    diagnostics$network_name <- network_name
    diagnostics <- diagnostics[, c(
        'network_name', 'cell_type', 'target', 'stage', 'converged',
        'iterations', 'objective', 'coef_change', 'selected_lambda',
        'cv_mean', 'cv_se', 'error_message'
    ), drop = FALSE]
    successful <- fit_cell_type$fits
    if (length(successful) == 0L) {
        log_message('No target genes were successfully fit for ', cell_type, '.', verbose = verbose)
        return(list(object = object, index = NULL, diagnostics = diagnostics, universal_id = NULL))
    }

    successful_features <- names(successful)
    fit_engine <- 'condition_sparse_common_scale_refit'
    fit_contract_key <- paste(network_name, safe_cell_type, sep = '__')
    fit_contract <- .condition_combine_fit_contracts(
        successful = successful,
        network_name = network_name,
        cell_type = cell_type,
        cell_type_col = cell_type_col,
        condition_col = condition_col,
        reference_condition = reference_condition_cell,
        candidate_screen = candidate_screen,
        scale = TRUE,
        fit_engine = fit_engine
    )
    fit_contract$cell_ids <- rownames(metadata_cell_type)
    fit_contract$cell_condition <- stats::setNames(
        as.character(metadata_cell_type[[condition_col]]),
        rownames(metadata_cell_type)
    )
    fit_contract$cv_block_col <- cv_block_col
    fit_contract$cv_block_status <- cv_block_status
    fit_contract$sample_blocked_oof_complete_for_cell_type <-
        !is.null(cv_block_status) &&
        all(cv_block_status$sample_blocked_oof_available)
    fit_contract$cell_block <- if (is.null(cv_block_col)) {
        NULL
    } else {
        stats::setNames(
            trimws(as.character(metadata_cell_type[[cv_block_col]])),
            rownames(metadata_cell_type)
        )
    }
    fit_contract$cell_provenance <- data.frame(
        cell_id = rownames(metadata_cell_type),
        cell_type = cell_type,
        condition = as.character(metadata_cell_type[[condition_col]]),
        cv_block = if (is.null(cv_block_col)) {
            NA_character_
        } else {
            trimws(as.character(metadata_cell_type[[cv_block_col]]))
        },
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
        'exact paired cell_ids used for fitting'
    fit_contracts <- object@grn@params$condition_grn_fits
    if (is.null(fit_contracts)) fit_contracts <- list()
    fit_contracts[[fit_contract_key]] <- fit_contract
    object@grn@params$condition_grn_fits <- fit_contracts
    universal_coefs <- do.call(rbind, lapply(successful, function(x) x$universal_coefs))
    universal_gof <- do.call(rbind, lapply(successful, function(x) x$universal_gof))
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
        fit_contract_key = fit_contract_key,
        lambda_selection = lambda_selection,
        nlambda = nlambda,
        nfolds = nfolds,
        cv_block_col = cv_block_col,
        scale = TRUE,
        active_tol = active_tol,
        seed = seed,
        upstream = upstream,
        downstream = downstream,
        extend = extend,
        only_tss = only_tss,
        peak_to_gene_method = peak_to_gene_method,
        tf_cor = tf_cor,
        peak_cor = peak_cor
    )

    universal_params <- do.call(
        .condition_network_params,
        c(list(network_level = 'shared', condition = NA_character_), common_params)
    )
    object@grn@networks[[universal_id]] <- .condition_build_network(
        successful_features, universal_coefs, universal_gof, universal_params
    )
    index_rows <- list(.condition_index_row(
        network_id = universal_id,
        network_name = network_name,
        cell_type = cell_type,
        network_level = 'shared',
        condition = NA_character_,
        n_cells = nrow(metadata_cell_type),
        coefs = universal_coefs,
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
            c(list(network_level = 'condition', condition = condition_name), common_params)
        )
        object@grn@networks[[condition_id]] <- .condition_build_network(
            successful_features, condition_coefs, condition_gof, condition_params
        )
        index_rows[[length(index_rows) + 1L]] <- .condition_index_row(
            network_id = condition_id,
            network_name = network_name,
            cell_type = cell_type,
            network_level = 'condition',
            condition = condition_name,
            n_cells = sum(as.character(metadata_cell_type[[condition_col]]) == condition_name),
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
        universal_id = universal_id
    )
}

.condition_skip_diagnostic <- function(network_name, cell_type, message_text) {
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
        stringsAsFactors = FALSE
    )
}
