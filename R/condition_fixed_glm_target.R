# Compact target helpers for the public fixed-dictionary Gaussian GLM path.
#
# This keeps fit_grn_from_edges(..., parallel = TRUE) from serializing the full
# prepared RNA/ATAC matrices to every worker.

.pando_fixed_glm_target_payload <- function(
    prepared, edge_dictionary, target, cells, rank_action, min_residual_df) {
    edges <- edge_dictionary[
        edge_dictionary$target == target, , drop = FALSE
    ]
    compact <- list(
        gene_data = prepared$gene_data[
            cells, unique(c(target, as.character(edges$tf))), drop = FALSE
        ],
        peak_data = prepared$peak_data[
            cells, unique(as.character(edges$region)), drop = FALSE
        ]
    )
    list(
        prepared = compact,
        edges = edges,
        cells = cells,
        target = target,
        rank_action = rank_action,
        min_residual_df = min_residual_df
    )
}

.pando_fixed_glm_target_worker <- function(payload) {
    edges <- payload$edges
    cells <- payload$cells
    target <- payload$target
    prepared <- payload$prepared
    terms <- sprintf("edge_%07d", edges$candidate_index)
    predictor <- vapply(seq_len(nrow(edges)), function(j) {
        as.numeric(prepared$gene_data[cells, edges$tf[[j]]]) *
            as.numeric(prepared$peak_data[cells, edges$region[[j]]])
    }, numeric(length(cells)))
    if (is.null(dim(predictor))) {
        predictor <- matrix(predictor, ncol = 1L)
    }
    fit <- .condition_fit_target_matrix(
        response = prepared$gene_data[cells, target],
        predictor = predictor,
        terms = terms,
        rank_action = payload$rank_action,
        min_residual_df = payload$min_residual_df
    )
    fit$coefs$target <- target
    fit$coefs$tf <- edges$tf
    fit$coefs$region <- edges$region
    fit$coefs$atac_feature_id <- edges$atac_feature_id
    fit$coefs$edge_id <- edges$edge_id
    fit$coefs$candidate_index <- edges$candidate_index
    fit$coefs$source_global <- edges$source_global
    fit$coefs$source_conditions <- edges$source_conditions
    fit$coefs$n_sources <- edges$n_sources
    fit$gof$target <- target
    fit$gof$nvariables_dictionary <- nrow(edges)
    fit$gof$nvariables_estimable <- sum(fit$coefs$estimable)
    fit$gof$n_zero_variance <- sum(fit$coefs$zero_variance)
    fit$gof$n_aliased <- sum(fit$coefs$aliased)
    fit
}
