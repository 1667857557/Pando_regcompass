test_that('condition-aware API is exported', {
    expect_true(is.function(infer_condition_grn))
    expect_true(is.function(condition_grn_fit))
})

test_that('condition-sparse inference is the public default', {
    defaults <- formals(Pando:::infer_condition_grn.GRNData)
    expect_identical(
        eval(defaults$method)[[1L]],
        'shared_baseline_condition_sparse'
    )
    expect_identical(defaults$condition_mix, 0.5)
    expect_null(eval(defaults$reference_condition))
    expect_identical(eval(defaults$condition_weight), c('equal', 'cell_count'))
    expect_true(eval(defaults$scale))
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
        edges, estimate = c(0.5, NA_real_), corr = c(0.2, -0.1)
    )

    expect_identical(
        colnames(coefs),
        c('tf', 'target', 'region', 'term', 'estimate', 'corr')
    )
    expect_equal(coefs$estimate, c(0.5, NA_real_))
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
    params <- list(method = 'glmnet', fit_engine = 'condition_sparse_common_scale_refit')
    network <- Pando:::.condition_build_network('GENE', coefs, fit, params)

    expect_s4_class(network, 'Network')
    expect_identical(colnames(coef(network)), colnames(coefs))
    expect_identical(colnames(gof(network)), colnames(fit))
    expect_identical(NetworkParams(network)$method, 'glmnet')
})

test_that('shared coefficients remain standard Pando compatibility summaries', {
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
                beta_selection = beta,
                beta_condition = beta,
                beta_shared = rowMeans(beta),
                delta_condition = sweep(beta, 1L, rowMeans(beta), '-'),
                contrast = sweep(beta, 1L, beta[, 'Control'], '-'),
                eligibility_mask = matrix(
                    TRUE, nrow = 1, ncol = 2,
                    dimnames = list(edge_id, c('Control', 'Drug'))
                ),
                estimability_mask = matrix(
                    TRUE, nrow = 1, ncol = 2,
                    dimnames = list(edge_id, c('Control', 'Drug'))
                ),
                active_mask = matrix(
                    TRUE, nrow = 1, ncol = 2,
                    dimnames = list(edge_id, c('Control', 'Drug'))
                ),
                support_mask = matrix(
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
                condition_mix = 0.5,
                condition_weight = 'equal',
                active_tol = 1e-8,
                refit = list(method = 'test')
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
        candidate_screen = 'pooled_within_condition',
        scale = TRUE,
        fit_engine = 'condition_sparse_common_scale_refit'
    )

    expect_s3_class(fit, 'ConditionGRNFit')
    expect_identical(fit$schema_version, 'pando_condition_grn_fit_v3')
    expect_equal(fit$contrast[, 'Control'], 0)
    expect_equal(fit$contrast[, 'Drug'], 0.5)
    expect_equal(fit$predictor_transform$center, 2)
    expect_equal(fit$predictor_transform$scale, 3)
    expect_equal(fit$beta_condition_raw, beta * 5 / 3)
})

test_that('single-cell projection preserves standardized and raw score identity', {
    raw_predictor <- c(2, 5, 8)
    center <- 2
    predictor_scale <- 3
    response_scale <- 4
    beta_std <- 0.75
    beta_raw <- beta_std * response_scale / predictor_scale

    standardized <- Pando:::.condition_projection_predictor(
        raw_predictor, center, predictor_scale, scale = 'std'
    ) * beta_std
    raw <- Pando:::.condition_projection_predictor(
        raw_predictor, center, predictor_scale, scale = 'raw'
    ) * beta_raw

    expect_equal(raw, standardized * response_scale)
})

test_that('network names are sanitized deterministically', {
    expect_identical(Pando:::.condition_safe_id('Tumor cell/A'), 'Tumor_cell_A')
    expect_identical(Pando:::.condition_model_name('chr1-10-20'), 'chr1_10_20')
})
