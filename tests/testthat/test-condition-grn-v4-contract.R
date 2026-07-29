test_that('biological blocks never cross validation folds', {
    X <- list(
        A = Matrix::Matrix(matrix(c(0, 1, 10, 11), ncol = 1L)),
        B = Matrix::Matrix(matrix(c(2, 3, 12, 13), ncol = 1L))
    )
    y <- list(A = c(0, 1, 10, 11), B = c(2, 3, 12, 13))
    blocks <- list(
        A = c('donor_1', 'donor_1', 'donor_2', 'donor_2'),
        B = c('donor_1', 'donor_1', 'donor_2', 'donor_2')
    )
    fit <- Pando:::.condition_cv_multitask_path(
        X_list = X,
        y_list = y,
        lambda = 0.1,
        coefficient_mask = matrix(TRUE, 1L, 2L),
        nfolds = 2L,
        block_list = blocks,
        standardize = TRUE,
        seed = 11L
    )

    for (task in seq_along(blocks)) {
        observed <- split(fit$oof_fold[[task]], blocks[[task]])
        expect_true(all(vapply(observed, function(x) {
            length(unique(x)) == 1L
        }, logical(1))))
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
    fit <- Pando:::.condition_cv_multitask_path(
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

test_that('single-sample conditions disable replicate-level OOF only', {
    status <- Pando:::.condition_sample_block_status(list(
        Control = rep('sample_1', 20L),
        Drug = rep(c('sample_2', 'sample_3'), each = 10L)
    ))

    expect_false(status$available)
    expect_identical(
        status$reason,
        'one_or_more_conditions_have_fewer_than_two_biological_samples'
    )
    expect_equal(status$n_blocks, c(Control = 1L, Drug = 2L))
})

test_that('cell-first aggregation permits mixed donors and retains covariance', {
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
                cv_block = c('donor_1', 'donor_2'),
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
    expect_equal(aggregated$group_metadata['SC_1', 'n_cv_blocks'], 2L)
})
