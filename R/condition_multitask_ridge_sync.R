# Keep each stored Pando Network synchronized with the final active
# ConditionGRNFit after the single no-fusion ridge fit.

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
        params$projection_policy <- fit$projection_policy
        params$fit_dictionary_policy <- fit$fit_dictionary_policy
        params$fit_engine <- fit$fit_engine
        params$candidate_edge_count <- fit$candidate_edge_count
        params$fit_dictionary_edge_count <- fit$fit_dictionary_edge_count
        params$candidate_tf_cor <- fit$candidate_tf_cor
        params$candidate_peak_cor <- fit$candidate_peak_cor
        params$local_support_role <- fit$local_support_role
        params$statistical_support_role <- fit$statistical_support_role
        params$inference_scope <- fit$inference_scope
        methods::slot(network, "params") <- params
        object@grn@networks[[network_name]] <- network
    }
    object
}
