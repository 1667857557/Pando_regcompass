.reference_complete_target_fit <- function(
    X_raw,
    y_raw,
    mask,
    lambda,
    alpha = 0.5,
    condition_mix = 0.5,
    active_tol = 1e-8,
    outer_nfolds = 3L,
    inner_nfolds = 2L,
    lambda_selection = "lambda.1se",
    seed = 31L,
    max_iter = 2000L,
    tol_objective = 1e-8,
    tol_coef = 1e-7
) {
    cv <- Pando:::.condition_nested_crossfit_within_cell_type(
        X_list = X_raw,
        y_list = y_raw,
        lambda = lambda,
        lambda_auto = FALSE,
        alpha = alpha,
        condition_mix = condition_mix,
        condition_weight = "equal",
        coefficient_mask = mask,
        outer_nfolds = outer_nfolds,
        inner_nfolds = inner_nfolds,
        active_tol = active_tol,
        lambda_selection = lambda_selection,
        comparison_conditions = names(X_raw),
        seed = seed,
        max_iter = max_iter,
        tol_objective = tol_objective,
        tol_coef = tol_coef
    )
    full_cv <- Pando:::.condition_select_lambda_nested(
        X_list = X_raw,
        y_list = y_raw,
        lambda = lambda,
        alpha = alpha,
        condition_mix = condition_mix,
        condition_weight = "equal",
        coefficient_mask = mask,
        nfolds = inner_nfolds,
        active_tol = active_tol,
        lambda_selection = lambda_selection,
        seed = Pando:::.condition_seed_for("full-inner", seed),
        max_iter = max_iter,
        tol_objective = tol_objective,
        tol_coef = tol_coef
    )
    fold_stats <- Pando:::.condition_build_fold_statistics(
        X_raw, y_raw, mask
    )
    transform <- Pando:::.condition_build_balanced_transform(
        X_raw, y_raw, fold_statistics = fold_stats
    )
    scaled <- Pando:::.condition_apply_balanced_transform(
        X_raw, y_raw, transform
    )
    estimability <- Pando:::.condition_true_variance_mask(
        scaled$X, mask
    )
    required <- lambda[seq_len(full_cv$selected_index)]
    path <- Pando:::.condition_fit_multitask_path(
        X_list = scaled$X,
        y_list = scaled$y,
        lambda = required,
        alpha = alpha,
        condition_mix = condition_mix,
        condition_weight = "equal",
        coefficient_mask = mask,
        max_iter = max_iter,
        tol_objective = tol_objective,
        tol_coef = tol_coef,
        verified_estimability_mask = estimability
    )
    selected <- path$fits[[length(path$fits)]]
    refit <- Pando:::.condition_refit_shared_baseline(
        X_list = scaled$X,
        y_list = scaled$y,
        beta_selection = selected$beta,
        estimability_mask = estimability,
        ridge = max(selected$lambda * (1 - alpha), 1e-6),
        active_tol = active_tol,
        condition_weight = "equal",
        cache = Pando:::.condition_make_refit_cache(
            scaled$X, scaled$y, "equal"
        ),
        estimability_verified = TRUE
    )
    list(
        cv = cv,
        full_cv = full_cv,
        transform = transform,
        selected = selected,
        refit = refit
    )
}
