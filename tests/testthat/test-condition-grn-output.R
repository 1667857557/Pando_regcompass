test_that('condition-aware API is exported', {
    expect_true(is.function(infer_condition_grn))
    expect_true(is.function(condition_grn_fit))
})

test_that('public inference exposes only the canonical functional contract', {
    defaults <- formals(Pando:::infer_condition_grn.GRNData)
    expect_false('method' %in% names(defaults))
    expect_null(eval(defaults$cell_type))
    expect_identical(defaults$condition_mix, 0.5)
    expect_null(eval(defaults$reference_condition))
    expect_identical(eval(defaults$condition_weight), 'equal')
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
    params <- list(
        method = 'glmnet',
        fit_engine = 'condition_sparse_within_cell_type_oof_refit'
    )
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
                structural_candidate_mask = matrix(
                    TRUE, nrow = 1, ncol = 2,
                    dimnames = list(edge_id, c('Control', 'Drug'))
                ),
                screening_mask = matrix(
                    TRUE, nrow = 1, ncol = 2,
                    dimnames = list(edge_id, c('Control', 'Drug'))
                ),
                predictor_center = stats::setNames(2, edge_id),
                predictor_scale = stats::setNames(3, edge_id),
                response_center = 4,
                response_scale = 5,
                transform_policy =
                    'equal_condition_center_equal_condition_within_variance_v1',
                predictor_center_hash = 'center-hash',
                predictor_scale_hash = 'scale-hash',
                training_fold_only = TRUE,
                intercept = c(Control = 0, Drug = 0),
                condition_rsq = c(Control = 0.4, Drug = 0.5),
                condition_rsq_train = c(Control = 0.4, Drug = 0.5),
                condition_rsq_oof = c(Control = 0.2, Drug = 0.3),
                condition_rmse_oof = c(Control = 1, Drug = 1.2),
                target_rsq_oof_pooled = 0.25,
                cv_method =
                    'nested_outer_condition_stratified_cell_oof',
                oof_model =
                    'nested_selection_shared_baseline_refit_heldout_projection',
                predictive_oof_available = TRUE,
                oof_validation_level =
                    'outer_condition_stratified_heldout_cells',
                projection_common_oof = c(
                    c1 = 0.1, c2 = 0.2, c3 = 0.3, c4 = 0.4
                ),
                projection_condition_full_oof = c(
                    c1 = 0.1, c2 = 0.2, c3 = 0.3, c4 = 0.4
                ),
                projection_global_common_oof = c(
                    c1 = 0.1, c2 = 0.2, c3 = 0.3, c4 = 0.4
                ),
                projection_origin =
                    'outer_condition_stratified_cell_oof',
                projection_used_for_penalty = TRUE,
                full_fit_projection_used_for_penalty = FALSE,
                fold_transform_policy =
                    'equal_condition_center_equal_condition_within_variance_v1',
                oof_cell_coverage = 1,
                oof_projection_available_fraction = 1,
                oof_assignment_count = c(
                    c1 = 1L, c2 = 1L, c3 = 1L, c4 = 1L
                ),
                oof_fold = list(
                    Control = c(c1 = 1L, c2 = 2L),
                    Drug = c(c3 = 1L, c4 = 2L)
                ),
                cv_fold_transform = list(),
                cv_effective_nfolds = 2L,
                outer_nfolds = 2L,
                inner_nfolds = 2L,
                selected_lambda = 0.1,
                lambda_path = c(1, 0.1),
                cv_mean = c(2, 1),
                cv_se = c(0.2, 0.1),
                alpha = 0.5,
                condition_mix = 0.5,
                condition_weight = 'equal',
                active_tol = 1e-8,
                refit = list(method = 'test'),
                refit_stability = list(edge = data.frame(), status = 'test')
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
        cell_type = 'T',
        cell_type_col = 'cell_type',
        condition_col = 'condition',
        reference_condition = 'Control',
        comparison_conditions = c('Control', 'Drug'),
        candidate_screen = 'motif_domain',
        scale = TRUE,
        fit_engine = 'condition_sparse_within_cell_type_oof_refit'
    )

    expect_s3_class(fit, 'ConditionGRNFit')
    expect_identical(fit$schema_version, 'pando_condition_grn_fit_v5')
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
