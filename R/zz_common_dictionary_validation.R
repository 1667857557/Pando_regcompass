# Final validation for externally supplied common edge dictionaries.
# Fixed-dictionary fitting must preserve Pando's biological domain and motif
# support in addition to exact feature and region identity.

.condition_validate_dictionary <- function(dictionary, prepared) {
    required <- c(
        "edge_id", "target", "tf", "region", "atac_feature_id",
        "candidate_index"
    )
    if (!is.data.frame(dictionary) || !nrow(dictionary) ||
        !all(required %in% colnames(dictionary)) ||
        anyNA(dictionary$edge_id) || any(!nzchar(dictionary$edge_id)) ||
        anyDuplicated(dictionary$edge_id)) {
        stop("The common edge dictionary is invalid.", call. = FALSE)
    }
    expected <- paste(
        dictionary$target, dictionary$tf, dictionary$region, sep = "||"
    )
    if (!identical(as.character(dictionary$edge_id), expected)) {
        stop("Dictionary edge IDs do not match exact target-TF-region triples.",
             call. = FALSE)
    }
    candidate_index <- suppressWarnings(as.integer(dictionary$candidate_index))
    if (anyNA(candidate_index) || any(candidate_index < 1L) ||
        any(candidate_index != dictionary$candidate_index) ||
        anyDuplicated(candidate_index)) {
        stop("Dictionary candidate indices must be unique positive integers.",
             call. = FALSE)
    }
    if (anyNA(dictionary$atac_feature_id) ||
        any(!nzchar(as.character(dictionary$atac_feature_id))) ||
        any(!dictionary$target %in% colnames(prepared$gene_data)) ||
        any(!dictionary$tf %in% colnames(prepared$gene_data)) ||
        any(!dictionary$region %in% colnames(prepared$peak_data))) {
        stop("Dictionary targets, TFs, regions or ATAC IDs are absent.",
             call. = FALSE)
    }
    mapped <- prepared$region_map$atac_feature_id[
        match(dictionary$region, prepared$region_map$region)
    ]
    if (anyNA(mapped) || !identical(
        as.character(mapped), as.character(dictionary$atac_feature_id)
    )) {
        stop("Dictionary region-to-ATAC mappings differ from the fitted object.",
             call. = FALSE)
    }
    domain_supported <- vapply(seq_len(nrow(dictionary)), function(i) {
        value <- as.numeric(prepared$peaks2gene[
            dictionary$target[[i]], dictionary$region[[i]], drop = TRUE
        ])
        length(value) == 1L && is.finite(value) && value != 0
    }, logical(1))
    if (any(!domain_supported)) {
        bad <- dictionary$edge_id[!domain_supported]
        stop(
            "Dictionary contains target-region pairs outside the Pando domain: ",
            paste(utils::head(bad, 10L), collapse = ", "),
            call. = FALSE
        )
    }
    motif_supported <- vapply(seq_len(nrow(dictionary)), function(i) {
        motif_row <- prepared$peaks2motif[
            dictionary$region[[i]], , drop = FALSE
        ]
        motif_index <- which(as.numeric(motif_row) != 0)
        if (!length(motif_index)) return(FALSE)
        tf_support <- prepared$motif2tf[
            motif_index, dictionary$tf[[i]], drop = FALSE
        ]
        any(as.numeric(tf_support) != 0)
    }, logical(1))
    if (any(!motif_supported)) {
        bad <- dictionary$edge_id[!motif_supported]
        stop(
            "Dictionary contains peak-TF pairs without motif support: ",
            paste(utils::head(bad, 10L), collapse = ", "),
            call. = FALSE
        )
    }
    invisible(TRUE)
}

.condition_fit_target_matrix_common_dictionary_base <-
    .condition_fit_target_matrix

.condition_fit_target_matrix <- function(
    response, predictor, terms, rank_action = c("mark", "error"),
    min_residual_df = 1L) {
    rank_action <- match.arg(rank_action)
    result <- .condition_fit_target_matrix_common_dictionary_base(
        response = response,
        predictor = predictor,
        terms = terms,
        rank_action = rank_action,
        min_residual_df = min_residual_df
    )
    residual_df <- as.numeric(result$gof$residual_df[[1L]])
    if (is.finite(residual_df) && residual_df < as.integer(min_residual_df)) {
        result$coefs$estimate[] <- NA_real_
        result$coefs$std_err[] <- NA_real_
        result$coefs$statistic[] <- NA_real_
        result$coefs$pval[] <- NA_real_
        result$coefs$estimable[] <- FALSE
        result$gof$fit_status <- "insufficient_df"
        result$gof$intercept <- NA_real_
    }
    result
}

.condition_make_network_common_dictionary_base <- .condition_make_network

.condition_make_network <- function(
    coefficients, fit, dictionary, condition_label, params) {
    network <- .condition_make_network_common_dictionary_base(
        coefficients = coefficients,
        fit = fit,
        dictionary = dictionary,
        condition_label = condition_label,
        params = params
    )
    intercept <- as.character(network@coefs$term) == "(Intercept)"
    network@coefs$padj[intercept] <- 1
    network
}
