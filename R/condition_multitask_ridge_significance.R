# Significance policy for multi-task ridge condition GRNs.
#
# The exact dictionary and the joint ridge fit always retain every candidate
# edge. Statistical significance is applied only after fitting, when deciding
# which condition-specific coefficients may contribute to downstream
# projections. This keeps coefficient estimation comparable across conditions
# while requiring explicit statistical support for a realized edge effect.

.condition_max_padj_threshold <- 0.1
.condition_significant_projection_policy <- "padj_significant_ridge_effects"

.condition_validate_padj_threshold <- function(padj_threshold) {
    if (!is.numeric(padj_threshold) || length(padj_threshold) != 1L ||
        !is.finite(padj_threshold) || padj_threshold <= 0 ||
        padj_threshold > .condition_max_padj_threshold) {
        stop(
            "`padj_threshold` must be one finite number in (0, 0.1].",
            call. = FALSE
        )
    }
    as.numeric(padj_threshold)
}

.condition_apply_significance_gate <- function(fit) {
    if (!inherits(fit, "ConditionGRNFit")) {
        stop("A ConditionGRNFit is required for significance gating.",
             call. = FALSE)
    }
    threshold <- .condition_validate_padj_threshold(fit$padj_threshold)
    coefficient <- as.data.frame(fit$coefficients, stringsAsFactors = FALSE)
    required <- c("estimate", "estimable", "padj")
    if (!all(required %in% colnames(coefficient))) {
        stop(
            "Condition ridge coefficients require estimate, estimable and padj ",
            "before significance gating.", call. = FALSE
        )
    }
    estimate <- suppressWarnings(as.numeric(coefficient$estimate))
    padj <- suppressWarnings(as.numeric(coefficient$padj))
    significant <- coefficient$estimable %in% TRUE &
        is.finite(estimate) & is.finite(padj) & padj < threshold
    coefficient$significant <- significant
    coefficient$penalty_effect <- ifelse(significant, estimate, 0)
    fit$coefficients <- coefficient
    fit$projection_effect_column <- "penalty_effect"
    fit$projection_policy <- .condition_significant_projection_policy
    fit
}

# Keep the numerical estimator in condition_multitask_ridge_refit.R unchanged.
# This wrapper applies the inferential gate only after the joint fit and BH
# adjustment have been completed.
.condition_ridge_refit_contract_without_significance_gate <-
    .condition_ridge_refit_contract

.condition_ridge_refit_contract <- function(...) {
    answer <- .condition_ridge_refit_contract_without_significance_gate(...)
    fit <- .condition_apply_significance_gate(answer$fit)
    object <- answer$object

    for (condition in fit$condition_levels) {
        network_name <- fit$network_names[[condition]]
        network <- object@grn@networks[[network_name]]
        if (is.null(network)) next
        coefs_one <- fit$coefficients[
            as.character(fit$coefficients$condition) == condition,
            , drop = FALSE
        ]
        methods::slot(network, "coefs") <- coefs_one
        params <- methods::slot(network, "params")
        params$padj_threshold <- fit$padj_threshold
        params$projection_policy <- .condition_significant_projection_policy
        methods::slot(network, "params") <- params
        object@grn@networks[[network_name]] <- network
    }

    answer$object <- object
    answer$fit <- fit
    answer
}

# The pre-existing runtime accepts any threshold in (0, 1). The condition
# workflow intentionally caps the supported FDR threshold at 0.1: 0.05 remains
# the default, while 0.1 is the most permissive supported inclusion rule.
.pando_infer_condition_grn_one_without_padj_cap <-
    .pando_infer_condition_grn_one

.pando_infer_condition_grn_one <- function(..., padj_threshold = 0.05) {
    threshold <- .condition_validate_padj_threshold(padj_threshold)
    object <- .pando_infer_condition_grn_one_without_padj_cap(
        ..., padj_threshold = threshold
    )

    params <- object@grn@params
    if (!identical(params$analysis_mode, "condition_grn")) return(object)

    params$condition_projection_policy <-
        .condition_significant_projection_policy
    fits <- params$condition_grn_fits
    index <- params$condition_network_index
    if (is.list(fits) && length(fits) && is.data.frame(index) && nrow(index) &&
        all(c("cell_type", "condition") %in% colnames(index))) {
        n_projection <- integer(nrow(index))
        for (i in seq_len(nrow(index))) {
            type <- as.character(index$cell_type[[i]])
            condition <- as.character(index$condition[[i]])
            fit <- fits[[type]]
            if (is.null(fit)) next
            coefficient <- as.data.frame(fit$coefficients, stringsAsFactors = FALSE)
            n_projection[[i]] <- sum(
                as.character(coefficient$condition) == condition &
                    coefficient$significant %in% TRUE,
                na.rm = TRUE
            )
        }
        index$n_projection_edges <- n_projection
        index$n_significant_edges <- n_projection
        params$condition_network_index <- index
    }
    object@grn@params <- params
    object
}
