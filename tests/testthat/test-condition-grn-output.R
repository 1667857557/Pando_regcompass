test_that('condition-aware API is exported', {
    expect_true(is.function(infer_condition_grn))
    expect_true(is.function(condition_grn_fit))
})

test_that('shared-design independent inference is the public default', {
    defaults <- formals(Pando:::infer_condition_grn.GRNData)
    expect_identical(
        eval(defaults$method),
        c('shared_design_independent', 'multitask_glmnet')
    )
    expect_identical(defaults$condition_mix, 1)
    expect_null(eval(defaults$reference_condition))
    expect_identical(eval(defaults$condition_weight), c('equal', 'cell_count'))
    expect_true(defaults$scale)
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

test_that('Universal coefficients remain compatibility summaries only', {
    beta <- matrix(c(1, 3, -2, 2), nrow = 2, byrow = TRUE)
    expect_equal(rowMeans(beta), c(2, 0))
})

test_that('fit contracts retain reference contrasts and pooled transforms', {
    make_contract <- function(target, edge_id, beta) {
        list(
            fit_contract = list(
                target = target,
                edge_table = data.frame(
                    edge_id = edge_id, tf = 'TF1', target = target,
                    region = 'peak1', term = 'peak1:TF1'
                ),
                beta = beta,
                contrast = sweep(beta, 1L, beta[, 'Control'], '-'),
                eligibility_mask = matrix(
                    TRUE, nrow = 1, ncol = 2,
                    dimnames = list(edge_id, c('Control', 'Drug'))
                ),
                predictor_center = stats::setNames(2, edge_id),
                predictor_scale = stats::setNames(3, edge_id),
                response_center = 4,
                response_scale = 5,
                intercept = c(Control = 0, Drug = 0),
                condition_rsq = c(Control = 0.4, Drug = 0.5),
                selected_lambda = 0.1,
                lambda_path = c(1, 0.1),
                cv_mean = c(2, 1),
                cv_se = c(0.2, 0.1),
                alpha = 0.5,
                condition_mix = 1
            )
        )
    }
    beta <- matrix(
        c(0.25, 0.75), nrow = 1,
        dimnames = list(NULL, c('Control', 'Drug'))
    )
    fit <- Pando:::.condition_combine_fit_contracts(
        successful = list(make_contract('GENE1', 'edge1', beta)),
        network_name = 'condition_grn',
        cell_type = 'Tumor',
        cell_type_col = 'cell_type',
        condition_col = 'condition',
        reference_condition = 'Control',
        candidate_screen = 'condition_union',
        scale = TRUE,
        fit_engine = 'shared_design_independent_elastic_net'
    )

    expect_s3_class(fit, 'ConditionGRNFit')
    expect_identical(fit$schema_version, 'pando_condition_grn_fit_v2')
    expect_equal(fit$contrast[, 'Control'], 0)
    expect_equal(fit$contrast[, 'Drug'], 0.5)
    expect_equal(fit$predictor_transform$center, 2)
    expect_equal(fit$predictor_transform$scale, 3)
})

test_that('network names are sanitized deterministically', {
    expect_identical(Pando:::.condition_safe_id('Tumor cell/A'), 'Tumor_cell_A')
    expect_identical(Pando:::.condition_model_name('chr1-10-20'), 'chr1_10_20')
})
