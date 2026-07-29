test_that('OOF folds are condition-stratified within one cell type', {
    X <- list(
        A = Matrix::Matrix(matrix(c(0, 1, 10, 11), ncol = 1L)),
        B = Matrix::Matrix(matrix(c(2, 3, 12, 13), ncol = 1L))
    )
    y <- list(A = c(0, 1, 10, 11), B = c(2, 3, 12, 13))
    fit <- Pando:::.condition_crossfit_within_cell_type(
        X_list = X,
        y_list = y,
        lambda = 0.1,
        coefficient_mask = matrix(TRUE, 1L, 2L),
        nfolds = 2L,
        standardize = TRUE,
        seed = 11L
    )

    for (task in seq_along(X)) {
        expect_equal(sort(unique(fit$oof_fold[[task]])), 1:2)
    }
    for (fold in seq_len(fit$effective_nfolds)) {
        training <- unlist(lapply(seq_along(X), function(task) {
            as.numeric(X[[task]][fit$oof_fold[[task]] != fold, 1L])
        }))
        expect_equal(
            fit$fold_transform[[fold]]$predictor_center,
            mean(training)
        )
        expect_true(fit$fold_transform[[fold]]$training_only)
    }
    expect_identical(
        fit$oof_model,
        'condition_sparse_selection_plus_common_metric_refit'
    )
})

test_that('OOF assignment aligns independently for unequal condition sizes', {
    X <- list(
        A = Matrix::Matrix(matrix(seq_len(4), ncol = 1L)),
        B = Matrix::Matrix(matrix(seq_len(6), ncol = 1L))
    )
    y <- list(A = seq_len(4), B = seq_len(6))
    fit <- Pando:::.condition_crossfit_within_cell_type(
        X_list = X,
        y_list = y,
        lambda = 0.1,
        coefficient_mask = matrix(TRUE, 1L, 2L),
        nfolds = 2L,
        standardize = TRUE,
        seed = 17L
    )

    expect_length(fit$oof_prediction$A, 4L)
    expect_length(fit$oof_prediction$B, 6L)
    expect_true(all(is.finite(unlist(fit$oof_prediction))))
})

test_that('pooled OOF R-squared excludes between-condition mean shifts', {
    y <- list(A = c(0, 1), B = c(100, 101))
    intercept_only <- list(A = c(0.5, 0.5), B = c(100.5, 100.5))

    expect_equal(
        Pando:::.condition_pooled_task_rsq(y, intercept_only),
        0
    )
    expect_gt(
        Pando:::.condition_rsq(
            unlist(y, use.names = FALSE),
            unlist(intercept_only, use.names = FALSE)
        ),
        0.99
    )
})

test_that('validation loss matches the requested condition weighting', {
    task_mse <- c(A = 1, B = 9)
    n_validation <- c(A = 90, B = 10)
    expect_equal(
        Pando:::.condition_task_validation_loss(
            task_mse, n_validation, condition_weight = 'equal'
        ),
        5
    )
    expect_equal(
        Pando:::.condition_task_validation_loss(
            task_mse, n_validation, condition_weight = 'cell_count'
        ),
        1.8
    )
})

test_that('OOF uses cells only and has no external blocking argument', {
    fit_text <- paste(
        deparse(body(Pando:::.condition_fit_target)), collapse = '\n'
    )
    expect_match(
        fit_text,
        'within_cell_type_condition_stratified_cell_oof',
        fixed = TRUE
    )
    expect_match(fit_text, 'predictive_oof_available', fixed = TRUE)
    expect_false(any(
        c('cell_type_block', 'cv_block') %in%
            names(formals(Pando:::.condition_fit_target))
    ))
})

test_that('cell_type selects exact independent analysis scopes', {
    metadata <- data.frame(
        cell_type = c('T', 'T', 'B', 'Myeloid'),
        condition = c('Control', 'Drug', 'Control', 'Drug'),
        stringsAsFactors = FALSE
    )

    expect_identical(
        Pando:::.condition_resolve_cell_types(
            metadata, 'cell_type', cell_type = 'T'
        ),
        'T'
    )
    expect_identical(
        Pando:::.condition_resolve_cell_types(
            metadata, 'cell_type', cell_type = c('B', 'T')
        ),
        c('B', 'T')
    )
    expect_identical(
        Pando:::.condition_resolve_cell_types(
            metadata, 'cell_type', cell_type = NULL
        ),
        c('T', 'B', 'Myeloid')
    )
    expect_error(
        Pando:::.condition_resolve_cell_types(
            metadata, 'cell_type', cell_type = 'Unknown'
        ),
        'not found'
    )
})

test_that('analysis metadata rejects ambiguous biological labels', {
    valid <- data.frame(
        cell_type = c('T', 'B'),
        condition = c('Control', 'Drug'),
        stringsAsFactors = FALSE
    )
    expect_silent(Pando:::.condition_validate_analysis_metadata(
        valid, 'cell_type', 'condition'
    ))

    whitespace <- valid
    whitespace$cell_type[[1L]] <- ' T'
    expect_error(
        Pando:::.condition_validate_analysis_metadata(
            whitespace, 'cell_type', 'condition'
        ),
        'surrounding whitespace'
    )
    blank <- valid
    blank$condition[[2L]] <- ''
    expect_error(
        Pando:::.condition_validate_analysis_metadata(
            blank, 'cell_type', 'condition'
        ),
        'non-empty'
    )
})

test_that('cell-first aggregation retains interaction covariance', {
    tf <- c(cell_1 = 0, cell_2 = 2)
    atac <- c(cell_1 = 0, cell_2 = 2)
    projection <- structure(
        list(
            gene_score = matrix(
                tf * atac,
                ncol = 1L,
                dimnames = list(names(tf), 'GENE')
            ),
            cell_metadata = data.frame(
                cell_id = names(tf),
                cell_type = 'T',
                condition = 'Control',
                row.names = names(tf)
            )
        ),
        class = c('ConditionGRNProjection', 'list')
    )
    membership <- data.frame(
        cell_id = names(tf),
        metacell_id = 'SC_1',
        stringsAsFactors = FALSE
    )

    aggregated <- aggregate_condition_grn_projection(
        projection, membership
    )

    expect_equal(aggregated$gene_score['SC_1', 'GENE'], mean(tf * atac))
    expect_false(
        isTRUE(all.equal(
            aggregated$gene_score['SC_1', 'GENE'], mean(tf) * mean(atac)
        ))
    )
    expect_equal(aggregated$group_metadata['SC_1', 'n_cells'], 2L)
})
