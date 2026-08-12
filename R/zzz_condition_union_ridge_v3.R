# Canonical condition-GRN v3 contract.
#
# The fitted dictionary is the union of condition-wise Pando candidates. An edge
# is admitted only when peak-target and TF-target correlation gates are both
# satisfied in at least one condition. The frozen union dictionary is fitted
# exactly once with a common-lambda ridge model and no cross-condition fusion.

.condition_common_dictionary_schema <- "pando_condition_grn_common_dictionary_v2"
.condition_multitask_ridge_schema <- "pando_condition_grn_multitask_ridge_v3"
.condition_fit_dictionary_policy <-
    "condition_union_pando_correlation_supported_frozen_dictionary"
.condition_significant_projection_policy <-
    "active_condition_pando_support_and_bh_ridge_effects"
.condition_fit_engine <- "condition_union_single_no_fusion_common_lambda_ridge"

.condition_ridge_control <- function(control = list()) {
    if (is.null(control)) control <- list()
    if (!is.list(control)) {
        stop("`ridge_control` must be a list.", call. = FALSE)
    }
    defaults <- list(
        lambda_grid = 10^seq(-3, 2, length.out = 9L),
        lambda_rule = "1se",
        fusion_ratio = 0,
        cv_folds = 5L,
        seed = 1L,
        scale_floor = 1e-8
    )
    unknown <- setdiff(names(control), names(defaults))
    if (length(unknown)) {
        stop("Unknown `ridge_control` field(s): ",
             paste(unknown, collapse = ", "), call. = FALSE)
    }
    out <- utils::modifyList(defaults, control)
    out$lambda_grid <- sort(unique(as.numeric(out$lambda_grid)))
    if (!length(out$lambda_grid) || any(!is.finite(out$lambda_grid)) ||
        any(out$lambda_grid <= 0)) {
        stop("`ridge_control$lambda_grid` must contain positive finite values.",
             call. = FALSE)
    }
    out$lambda_rule <- match.arg(out$lambda_rule, c("1se", "min"))
    if (!is.numeric(out$fusion_ratio) || length(out$fusion_ratio) != 1L ||
        !is.finite(out$fusion_ratio) || out$fusion_ratio != 0) {
        stop(
            "Canonical condition GRNs require `ridge_control$fusion_ratio = 0` ",
            "so condition coefficients are not shrunk toward one another.",
            call. = FALSE
        )
    }
    if (!is.numeric(out$cv_folds) || length(out$cv_folds) != 1L ||
        !is.finite(out$cv_folds) || out$cv_folds < 2L ||
        out$cv_folds != as.integer(out$cv_folds)) {
        stop("`ridge_control$cv_folds` must be an integer >= 2.",
             call. = FALSE)
    }
    if (!is.numeric(out$seed) || length(out$seed) != 1L ||
        !is.finite(out$seed) || out$seed != as.integer(out$seed)) {
        stop("`ridge_control$seed` must be a finite integer.", call. = FALSE)
    }
    if (!is.numeric(out$scale_floor) || length(out$scale_floor) != 1L ||
        !is.finite(out$scale_floor) || out$scale_floor <= 0) {
        stop("`ridge_control$scale_floor` must be positive and finite.",
             call. = FALSE)
    }
    out$cv_folds <- as.integer(out$cv_folds)
    out$seed <- as.integer(out$seed)
    out$fusion_ratio <- 0
    out
}

.condition_ridge_penalty <- function(k, p, fusion_ratio = 0) {
    if (!is.numeric(fusion_ratio) || length(fusion_ratio) != 1L ||
        !is.finite(fusion_ratio) || fusion_ratio != 0) {
        stop("Condition-GRN ridge penalty does not permit fusion.", call. = FALSE)
    }
    diag(k * p)
}

union_grn_edges <- function(global_edges = NULL, condition_edges) {
    if (!is.list(condition_edges) || !length(condition_edges) ||
        is.null(names(condition_edges)) || any(!nzchar(names(condition_edges)))) {
        stop("`condition_edges` must be a non-empty named list.",
             call. = FALSE)
    }
    required <- c(
        "target", "tf", "region", "atac_feature_id",
        "peak_target_cor", "tf_target_cor"
    )
    if (any(vapply(condition_edges, function(x) {
        !is.data.frame(x) || !all(required %in% colnames(x))
    }, logical(1)))) {
        stop(
            "Every condition candidate table requires target, tf, region, ATAC ",
            "IDs and both Pando correlation statistics.", call. = FALSE
        )
    }

    provenance_fields <- c(
        "rna_layer", "peak_layer", "peak_value_type",
        "preprocessing_fingerprint", "dictionary_input_schema"
    )
    provenance <- lapply(condition_edges, function(candidate) {
        stats::setNames(lapply(provenance_fields, function(field) {
            attr(candidate, field, exact = TRUE)
        }), provenance_fields)
    })
    complete_provenance <- vapply(provenance, function(value) {
        all(vapply(value, function(item) {
            is.character(item) && length(item) == 1L &&
                !is.na(item) && nzchar(item)
        }, logical(1)))
    }, logical(1))
    if (any(complete_provenance) && !all(complete_provenance)) {
        stop("Condition candidates mix verified and unverified preprocessing provenance.",
             call. = FALSE)
    }
    if (all(complete_provenance)) {
        reference <- provenance[[1L]]
        same_reference <- vapply(provenance, function(value) {
            identical(value, reference)
        }, logical(1))
        if (!all(same_reference)) {
            stop(
                "Condition candidates use different RNA/ATAC layers, value ",
                "semantics, or preprocessing fingerprints.", call. = FALSE
            )
        }
    }

    all_rows <- do.call(rbind, lapply(seq_along(condition_edges), function(i) {
        x <- as.data.frame(condition_edges[[i]], stringsAsFactors = FALSE)
        if (!nrow(x)) return(NULL)
        x$.union_condition <- names(condition_edges)[[i]]
        x
    }))
    if (is.null(all_rows) || !nrow(all_rows)) {
        stop(
            "No TF-peak-target edge passes the Pando peak_cor and tf_cor gates ",
            "in any condition.", call. = FALSE
        )
    }
    if (any(!is.finite(as.numeric(all_rows$peak_target_cor))) ||
        any(!is.finite(as.numeric(all_rows$tf_target_cor)))) {
        stop("Condition candidate correlation statistics must be finite.",
             call. = FALSE)
    }

    all_rows$edge_id <- paste(
        all_rows$target, all_rows$tf, all_rows$region, sep = "||"
    )
    support_table <- unique(all_rows[, c(
        "edge_id", "target", "tf", "region", "atac_feature_id",
        ".union_condition", "peak_target_cor", "tf_target_cor"
    ), drop = FALSE])
    colnames(support_table)[
        colnames(support_table) == ".union_condition"
    ] <- "condition"
    support_table <- support_table[order(
        support_table$edge_id, support_table$condition
    ), , drop = FALSE]
    rownames(support_table) <- NULL

    groups <- split(seq_len(nrow(all_rows)), all_rows$edge_id)
    dictionary <- do.call(rbind, lapply(groups, function(index) {
        one <- all_rows[index, , drop = FALSE]
        out <- one[1L, c("target", "tf", "region", "atac_feature_id"),
                   drop = FALSE]
        conditions <- sort(unique(as.character(one$.union_condition)))
        out$source_global <- FALSE
        out$source_conditions <- paste(conditions, collapse = ";")
        out$n_sources <- length(conditions)
        out$n_support_conditions <- length(conditions)
        out$support_conditions <- paste(conditions, collapse = ";")
        out$max_abs_peak_target_cor <- max(abs(as.numeric(one$peak_target_cor)))
        out$max_abs_tf_target_cor <- max(abs(as.numeric(one$tf_target_cor)))
        out
    }))
    dictionary$edge_id <- paste(
        dictionary$target, dictionary$tf, dictionary$region, sep = "||"
    )
    dictionary <- dictionary[order(
        dictionary$target, dictionary$tf, dictionary$region
    ), , drop = FALSE]
    dictionary$candidate_index <- seq_len(nrow(dictionary))
    rownames(dictionary) <- NULL
    if (anyDuplicated(dictionary$edge_id)) {
        stop("Condition-union dictionary produced duplicated triples.",
             call. = FALSE)
    }
    class(dictionary) <- c("PandoEdgeDictionary", "data.frame")
    attr(dictionary, "preprocessing_provenance_verified") <-
        all(complete_provenance)
    if (all(complete_provenance)) {
        for (field in provenance_fields) {
            attr(dictionary, field) <- provenance[[1L]][[field]]
        }
    }
    attr(dictionary, "condition_support_table") <- support_table
    attr(dictionary, "dictionary_policy") <- .condition_fit_dictionary_policy
    attr(dictionary, "global_candidate_input_ignored") <- !is.null(global_edges)
    dictionary
}

.condition_annotate_local_pando_support <- function(fit) {
    if (!inherits(fit, "ConditionGRNFit")) {
        stop("A ConditionGRNFit is required for Pando support annotation.",
             call. = FALSE)
    }
    coefficient <- as.data.frame(fit$coefficients, stringsAsFactors = FALSE)
    support <- attr(
        fit$edge_dictionary, "condition_support_table", exact = TRUE
    )
    required_coef <- c("edge_id", "condition", "estimate", "estimable", "padj")
    required_support <- c(
        "edge_id", "condition", "peak_target_cor", "tf_target_cor"
    )
    if (!all(required_coef %in% colnames(coefficient)) ||
        !is.data.frame(support) ||
        !all(required_support %in% colnames(support))) {
        stop("Condition-union Pando support metadata are incomplete.",
             call. = FALSE)
    }
    support_key <- paste(support$edge_id, support$condition, sep = "\001")
    if (anyDuplicated(support_key)) {
        stop("Condition-union support contains duplicated edge-condition rows.",
             call. = FALSE)
    }
    coefficient_key <- paste(
        coefficient$edge_id, coefficient$condition, sep = "\001"
    )
    index <- match(coefficient_key, support_key)
    local_support <- !is.na(index)
    coefficient$peak_target_cor <- NA_real_
    coefficient$tf_target_cor <- NA_real_
    coefficient$peak_target_cor[local_support] <-
        as.numeric(support$peak_target_cor[index[local_support]])
    coefficient$tf_target_cor[local_support] <-
        as.numeric(support$tf_target_cor[index[local_support]])
    coefficient$peak_cor_pass <- local_support
    coefficient$tf_cor_pass <- local_support
    coefficient$local_support <- local_support

    threshold <- .condition_validate_padj_threshold(fit$padj_threshold)
    estimate <- suppressWarnings(as.numeric(coefficient$estimate))
    padj <- suppressWarnings(as.numeric(coefficient$padj))
    statistically_supported <- coefficient$estimable %in% TRUE &
        is.finite(estimate) & is.finite(padj) & padj < threshold
    active <- statistically_supported & local_support
    coefficient$statistically_supported <- statistically_supported
    coefficient$active <- active
    # Compatibility alias used by legacy Pando/RegCompass accessors. `padj`
    # remains the condition-wise BH-adjusted ridge-Wald quantity; `significant`
    # here means eligible in the active condition GRN.
    coefficient$significant <- active
    coefficient$penalty_effect <- ifelse(active, estimate, 0)
    fit$coefficients <- coefficient
    fit$projection_effect_column <- "penalty_effect"
    fit$projection_policy <- .condition_significant_projection_policy
    fit$local_support_role <-
        "condition_specific_pando_peak_and_tf_correlation_activity_gate"
    fit$statistical_support_role <-
        "condition_wise_BH_adjusted_approximate_ridge_wald"
    fit
}

.condition_update_network_significance <- function(object, fit) {
    for (condition in fit$condition_levels) {
        network_name <- fit$network_names[[condition]]
        network <- object@grn@networks[[network_name]]
        if (is.null(network)) next
        one <- fit$coefficients[
            as.character(fit$coefficients$condition) == condition,
            , drop = FALSE
        ]
        network@coefs <- one
        network@params$projection_policy <- fit$projection_policy
        network@params$fit_engine <- fit$fit_engine
        network@params$fit_dictionary_policy <- fit$fit_dictionary_policy
        network@params$inference_scope <- fit$inference_scope
        network@params$local_support_role <- fit$local_support_role
        object@grn@networks[[network_name]] <- network
    }
    object
}

.condition_apply_significance_gate <- function(fit) {
    .condition_annotate_local_pando_support(fit)
}

.condition_ridge_refit_contract <- function(
    object, fit, prepared, control, rank_action = "mark",
    min_residual_df = 1L, parallel = FALSE, verbose = TRUE) {
    fit$adjust_method <- .condition_validate_adjust_method(fit$adjust_method)
    fit$padj_threshold <- .condition_validate_padj_threshold(fit$padj_threshold)
    control <- .condition_ridge_control(control)
    dictionary <- fit$edge_dictionary
    progress_label <- fit$cell_type %||% ""

    if (isTRUE(verbose)) {
        message(
            "Pando condition phase=single_no_fusion_ridge_start",
            " | cell_type=", as.character(progress_label),
            ";fit_edges=", nrow(dictionary),
            ";targets=", length(unique(as.character(dictionary$target)))
        )
    }
    final <- .condition_ridge_refit_contract_one_pass(
        object = object, fit = fit, prepared = prepared, control = control,
        rank_action = rank_action, min_residual_df = min_residual_df,
        parallel = parallel, verbose = verbose,
        progress_phase = "ridge_single",
        progress_label = progress_label
    )
    final$fit <- .condition_annotate_local_pando_support(final$fit)
    final$fit$schema_version <- .condition_common_dictionary_schema
    final$fit$model_schema <- .condition_multitask_ridge_schema
    final$fit$fit_engine <- .condition_fit_engine
    final$fit$fit_dictionary_policy <- .condition_fit_dictionary_policy
    final$fit$candidate_edge_count <- nrow(dictionary)
    final$fit$fit_dictionary_edge_count <- nrow(dictionary)
    final$fit$edge_dictionary <- dictionary
    final$fit$dictionary_support_summary <- dictionary[, intersect(c(
        "edge_id", "support_conditions", "n_support_conditions",
        "max_abs_peak_target_cor", "max_abs_tf_target_cor"
    ), colnames(dictionary)), drop = FALSE]
    final$fit$dictionary_support_table <- attr(
        dictionary, "condition_support_table", exact = TRUE
    )
    final$fit$inference_scope <-
        "approximate_ridge_wald_conditional_on_condition_union_pando_screened_dictionary_and_cv_lambda"
    final$fit$ridge_control <- control
    final$fit$projection_effect_column <- "penalty_effect"
    final$fit$projection_policy <- .condition_significant_projection_policy
    final$object <- .condition_update_network_significance(
        final$object, final$fit
    )
    if (isTRUE(verbose)) {
        message(
            "Pando condition phase=condition_fit_complete",
            " | cell_type=", as.character(progress_label),
            ";fit_edges=", nrow(dictionary),
            ";active_edges=",
            sum(final$fit$coefficients$active %in% TRUE),
            ";targets=", length(unique(as.character(dictionary$target)))
        )
    }
    final
}
