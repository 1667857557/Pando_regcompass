# Synchronize stored Pando Network objects with the final E-star/JSE fit.

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
        params$model_schema <- fit$model_schema
        params$inference_schema <- fit$inference_schema
        params$reference_condition <- fit$reference_condition
        params$candidate_edge_count <- fit$candidate_edge_count
        params$fit_dictionary_edge_count <- fit$fit_dictionary_edge_count
        params$candidate_tf_cor <- fit$candidate_tf_cor
        params$candidate_peak_cor <- fit$candidate_peak_cor
        params$dictionary_support_role <- fit$dictionary_support_role
        params$local_support_role <- fit$local_support_role
        params$global_support_role <- fit$global_support_role
        params$statistical_support_role <- fit$statistical_support_role
        params$regcompass_edge_gate <- fit$regcompass_edge_gate
        params$inference_scope <- fit$inference_scope
        methods::slot(network, "params") <- params
        object@grn@networks[[network_name]] <- network
    }
    object
}
