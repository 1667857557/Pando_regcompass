test_that('condition-aware API is exported', {
    expect_true(is.function(infer_condition_grn))
})

test_that('condition coefficient tables preserve Pando columns', {
    edges <- data.frame(
        tf = c('TF1', 'TF2'),
        target = c('GENE', 'GENE'),
        region = c('chr1-1-10', 'chr1-20-30'),
        term = c('chr1_1_10:TF1', 'chr1_20_30:TF2'),
        stringsAsFactors = FALSE
    )
    coefs <- Pando:::.condition_format_coefs(
        edges, estimate = c(0.5, 0), corr = c(0.2, -0.1)
    )

    expect_identical(
        colnames(coefs),
        c('tf', 'target', 'region', 'term', 'estimate', 'corr')
    )
    expect_equal(coefs$estimate, c(0.5, 0))
})

test_that('generated networks remain standard Pando Network objects', {
    coefs <- data.frame(
        tf = 'TF1',
        target = 'GENE',
        region = 'chr1-1-10',
        term = 'chr1_1_10:TF1',
        estimate = 0.5,
        corr = 0.2,
        stringsAsFactors = FALSE
    )
    fit <- data.frame(
        target = 'GENE',
        lambda = 0.1,
        rsq = 0.5,
        alpha = 0.5,
        nvariables = 1L,
        stringsAsFactors = FALSE
    )
    params <- list(method = 'glmnet', fit_engine = 'multitask_sparse_group_glmnet')
    network <- Pando:::.condition_build_network('GENE', coefs, fit, params)

    expect_s4_class(network, 'Network')
    expect_identical(colnames(coef(network)), colnames(coefs))
    expect_identical(colnames(gof(network)), colnames(fit))
    expect_identical(NetworkParams(network)$method, 'glmnet')
})

test_that('Universal coefficients are equal-condition means', {
    beta <- matrix(c(1, 3, -2, 2), nrow = 2, byrow = TRUE)
    expect_equal(rowMeans(beta), c(2, 0))
})

test_that('network names are sanitized deterministically', {
    expect_identical(Pando:::.condition_safe_id('Tumor cell/A'), 'Tumor_cell_A')
    expect_identical(Pando:::.condition_model_name('chr1-10-20'), 'chr1_10_20')
})
