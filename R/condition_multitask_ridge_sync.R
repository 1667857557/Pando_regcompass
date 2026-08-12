# Keep each stored Pando Network synchronized with the screened ConditionGRNFit.
# This file loads after condition_multitask_ridge_significance.R.

.condition_update_network_significance <- function(object, fit) {
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
        params$edge_dictionary <- fit$edge_dictionary
        params$padj_threshold <- fit$padj_threshold
        params$projection_policy <- .condition_significant_projection_policy
        params$fit_dictionary_policy <- .condition_fit_dictionary_policy
        params$candidate_edge_count <- fit$candidate_edge_count
        params$fit_dictionary_edge_count <- fit$fit_dictionary_edge_count
        params$screening_inference_scope <- fit$screening_inference_scope
        params$screening_adjust_method <- fit$screening_adjust_method
        params$screening_padj_threshold <- fit$screening_padj_threshold
        params$inference_performed <- FALSE
        params$inference_scope <- fit$inference_scope
        methods::slot(network, "params") <- params
        object@grn@networks[[network_name]] <- network
    }
    object
}
