# Common-dictionary condition-specific extension of Pando.

.condition_common_dictionary_schema <- "pando_condition_grn_common_dictionary_v1"

.condition_hash_object <- function(value) {
    file <- tempfile("pando_condition_fingerprint_", fileext = ".rds")
    on.exit(unlink(file), add = TRUE)
    saveRDS(value, file, version = 2L, compress = FALSE)
    unname(tools::md5sum(file)[[1L]])
}

.condition_preprocessing_fingerprint <- function(
    object, gene_data, peak_data, rna_layer, peak_layer,
    peak_value_type = "normalized") {
    params <- Params(object)
    command_provenance <- tryCatch({
        commands <- object@data@commands
        lapply(commands, function(command) {
            list(
                name = tryCatch(as.character(command@name),
                                error = function(error) NA_character_),
                assay = tryCatch(as.character(command@assay.used),
                                 error = function(error) NA_character_),
                call = tryCatch(as.character(command@call.string),
                                error = function(error) NA_character_)
            )
        })
    }, error = function(error) list())
    .condition_hash_object(list(
        schema = "pando_common_dictionary_preprocessing_v1",
        rna_assay = params$rna_assay,
        peak_assay = params$peak_assay,
        rna_layer = rna_layer,
        peak_layer = peak_layer,
        peak_value_type = peak_value_type,
        gene_matrix_class = class(gene_data),
        peak_matrix_class = class(peak_data),
        gene_dim = dim(gene_data),
        peak_dim = dim(peak_data),
        gene_cells = rownames(gene_data),
        peak_cells = rownames(peak_data),
        gene_features = colnames(gene_data),
        peak_features = colnames(peak_data),
        preprocessing_commands = command_provenance
    ))
}

.condition_safe_label <- function(x) {
    value <- gsub("[^[:alnum:]_.-]+", "_", as.character(x))
    value[!nzchar(value)] <- "unnamed"
    value
}

.condition_resolve_levels <- function(metadata, condition_col) {
    if (is.null(condition_col) || !is.character(condition_col) ||
        length(condition_col) != 1L || is.na(condition_col) ||
        !nzchar(trimws(condition_col)) || !condition_col %in% colnames(metadata)) {
        return(character())
    }
    value <- trimws(as.character(metadata[[condition_col]]))
    value <- value[!is.na(value) & nzchar(value)]
    unique(value)
}

.condition_validate_labels <- function(metadata, columns) {
    missing <- setdiff(columns, colnames(metadata))
    if (length(missing)) {
        stop("Missing metadata column(s): ", paste(missing, collapse = ", "),
             call. = FALSE)
    }
    for (column in columns) {
        value <- as.character(metadata[[column]])
        if (anyNA(value) || any(!nzchar(trimws(value))) || any(value != trimws(value))) {
            stop("Metadata column `", column,
                 "` must contain complete labels without surrounding whitespace.",
                 call. = FALSE)
        }
    }
    invisible(TRUE)
}

.condition_prepare_common_input <- function(
    object, genes = NULL, peak_to_gene_method = c("Signac", "GREAT"),
    upstream = 100000, downstream = 0, extend = 1000000,
    only_tss = FALSE, peak_to_gene_domains = NULL,
    rna_layer = "data", peak_layer = "data",
    peak_value_type = c("normalized", "probability", "other"),
    verbose = TRUE) {
    if (!inherits(object, "GRNData")) {
        stop("`object` must be a GRNData object.", call. = FALSE)
    }
    peak_to_gene_method <- match.arg(peak_to_gene_method)
    peak_value_type <- match.arg(peak_value_type)
    valid_layer <- function(value) {
        is.character(value) && length(value) == 1L && !is.na(value) &&
            nzchar(trimws(value))
    }
    if (!valid_layer(rna_layer) || !valid_layer(peak_layer)) {
        stop("RNA and ATAC layer names must be non-empty strings.",
             call. = FALSE)
    }
    params <- Params(object)
    motif2tf <- NetworkTFs(object)
    if (is.null(motif2tf)) {
        stop("Motif matches are missing; run `find_motifs()` first.",
             call. = FALSE)
    }
    gene_annot <- Signac::Annotation(GetAssay(object, params$peak_assay))
    if (is.null(gene_annot)) {
        stop("The peak assay requires a gene annotation.", call. = FALSE)
    }
    if (is.null(genes)) {
        genes <- VariableFeatures(object, assay = params$rna_assay)
    }
    genes <- unique(as.character(genes))
    if (!length(genes)) {
        stop("No target genes were supplied or found as variable features.",
             call. = FALSE)
    }

    gene_data <- Matrix::t(LayerData(
        object, assay = params$rna_assay, layer = rna_layer
    ))
    peak_data_all <- Matrix::t(LayerData(
        object, assay = params$peak_assay, layer = peak_layer
    ))
    common_cells <- intersect(rownames(gene_data), rownames(peak_data_all))
    if (!length(common_cells)) {
        stop("RNA and ATAC assays do not share paired cells.", call. = FALSE)
    }
    gene_data <- gene_data[common_cells, , drop = FALSE]
    peak_data_all <- peak_data_all[common_cells, , drop = FALSE]
    if (identical(peak_value_type, "probability")) {
        observed <- if (inherits(peak_data_all, "sparseMatrix")) {
            peak_data_all@x
        } else {
            as.numeric(peak_data_all)
        }
        if (any(!is.finite(observed)) || any(observed < 0 | observed > 1)) {
            stop("Probability-valued ATAC layers must be finite and in [0, 1].",
                 call. = FALSE)
        }
    }
    preprocessing_fingerprint <- .condition_preprocessing_fingerprint(
        object = object, gene_data = gene_data, peak_data = peak_data_all,
        rna_layer = rna_layer, peak_layer = peak_layer,
        peak_value_type = peak_value_type
    )

    features <- intersect(gene_annot$gene_name, genes)
    features <- intersect(features, colnames(gene_data))
    features <- unique(as.character(features))
    if (!length(features)) {
        stop("No requested target overlaps the RNA assay and annotation.",
             call. = FALSE)
    }
    gene_annot <- gene_annot[gene_annot$gene_name %in% features, ]

    regions <- NetworkRegions(object)
    if (!length(regions@peaks) || is.null(regions@motifs) ||
        !nrow(regions@motifs@data)) {
        stop("Pando regions lack measured-peak or motif mappings.", call. = FALSE)
    }
    if (any(regions@peaks < 1L | regions@peaks > ncol(peak_data_all))) {
        stop("Pando region-to-peak indices are out of range.", call. = FALSE)
    }
    peak_data <- peak_data_all[, regions@peaks, drop = FALSE]
    atac_feature_id <- colnames(peak_data)
    peaks2motif <- regions@motifs@data
    region_id <- rownames(peaks2motif)
    if (is.null(region_id) || length(region_id) != ncol(peak_data) ||
        anyNA(region_id) || any(!nzchar(region_id)) || anyDuplicated(region_id)) {
        stop("Motif regions require unique IDs aligned to measured peaks.",
             call. = FALSE)
    }
    colnames(peak_data) <- region_id

    log_message("Selecting candidate regulatory regions near targets",
                verbose = verbose)
    if (is.null(peak_to_gene_domains)) {
        peaks_near_gene <- find_peaks_near_genes(
            peaks = regions@ranges,
            method = peak_to_gene_method,
            genes = gene_annot,
            upstream = upstream,
            downstream = downstream,
            extend = extend,
            only_tss = only_tss
        )
    } else {
        peaks_near_gene <- find_peaks_near_genes(
            peaks = regions@ranges,
            method = "Signac",
            genes = peak_to_gene_domains,
            upstream = 0,
            downstream = 0,
            extend = extend,
            only_tss = FALSE
        )
    }
    peaks2gene <- aggregate_matrix(
        t(peaks_near_gene), groups = colnames(peaks_near_gene), fun = "sum"
    )

    common_regions <- Reduce(intersect, list(
        colnames(peaks2gene), rownames(peaks2motif), colnames(peak_data)
    ))
    common_regions <- colnames(peaks2gene)[
        colnames(peaks2gene) %in% common_regions
    ]
    if (!length(common_regions) || anyDuplicated(common_regions)) {
        stop("Peak-to-gene, motif and ATAC matrices do not share unique regions.",
             call. = FALSE)
    }
    peaks2gene <- peaks2gene[, common_regions, drop = FALSE]
    peaks2motif <- peaks2motif[common_regions, , drop = FALSE]
    peak_data <- peak_data[, common_regions, drop = FALSE]

    tfs <- intersect(colnames(motif2tf), colnames(gene_data))
    motif2tf <- motif2tf[, tfs, drop = FALSE]
    motif_ids <- intersect(colnames(peaks2motif), rownames(motif2tf))
    motif_ids <- colnames(peaks2motif)[colnames(peaks2motif) %in% motif_ids]
    if (!length(motif_ids) || anyDuplicated(motif_ids)) {
        stop("Motif matches and motif-to-TF mappings do not share unique IDs.",
             call. = FALSE)
    }
    peaks2motif <- peaks2motif[, motif_ids, drop = FALSE]
    motif2tf <- motif2tf[motif_ids, , drop = FALSE]

    peaks_at_gene <- as.logical(sparseMatrixStats::colMaxs(peaks2gene))
    peaks_with_motif <- as.logical(sparseMatrixStats::rowMaxs(peaks2motif * 1))
    keep <- peaks_at_gene & peaks_with_motif
    peaks2gene <- peaks2gene[, keep, drop = FALSE]
    peaks2motif <- peaks2motif[keep, , drop = FALSE]
    peak_data <- peak_data[, keep, drop = FALSE]

    region_lookup <- stats::setNames(atac_feature_id, region_id)
    region_map <- data.frame(
        region = colnames(peak_data),
        atac_feature_id = unname(region_lookup[colnames(peak_data)]),
        stringsAsFactors = FALSE
    )
    if (anyNA(region_map$atac_feature_id)) {
        stop("A regulatory region lacks a measured ATAC feature mapping.",
             call. = FALSE)
    }
    list(
        gene_data = gene_data,
        peak_data = peak_data,
        peak_data_all = peak_data_all,
        features = features,
        peaks2gene = peaks2gene,
        peaks2motif = peaks2motif,
        motif2tf = motif2tf,
        region_map = region_map,
        params = params,
        rna_layer = rna_layer,
        peak_layer = peak_layer,
        peak_value_type = peak_value_type,
        preprocessing_fingerprint = preprocessing_fingerprint,
        peak_to_gene_method = peak_to_gene_method,
        upstream = upstream,
        downstream = downstream,
        extend = extend,
        only_tss = only_tss
    )
}

.condition_discover_edges_prepared <- function(
    prepared, cells, source_label, source_type, tf_cor, peak_cor,
    parallel = FALSE, verbose = TRUE) {
    .condition_discover_edges_compact(
        prepared = prepared,
        cells = cells,
        source_label = source_label,
        source_type = source_type,
        tf_cor = tf_cor,
        peak_cor = peak_cor,
        parallel = parallel,
        verbose = verbose
    )
}

#' Discover Pando TF-peak-target candidate edges
#'
#' Runs Pando domain, motif, peak-target correlation and TF-target correlation
#' candidate steps without using the resulting coefficients as the final
#' condition effect.
#'
#' @param object A `GRNData` object after `find_motifs()`.
#' @param genes Target genes.
#' @param cells Paired cells used for candidate discovery.
#' @param source_label Label stored with the candidate source.
#' @param source_type Either `"global"` or `"condition"`.
#' @param tf_cor,peak_cor Absolute correlation thresholds; defaults are 0.05
#'   for both and user-supplied values in [0, 1] are preserved.
#' @param ... Regulatory-domain and execution arguments.
#' @return A `PandoEdgeDictionary` data frame.
#' @export
discover_grn_edges <- function(
    object, genes = NULL, cells = NULL, source_label = "global",
    source_type = c("global", "condition"), tf_cor = 0.05, peak_cor = 0.05,
    peak_to_gene_method = c("Signac", "GREAT"), upstream = 100000,
    downstream = 0, extend = 1000000, only_tss = FALSE,
    peak_to_gene_domains = NULL, rna_layer = "data", peak_layer = "data",
    peak_value_type = c("normalized", "probability", "other"),
    parallel = FALSE, verbose = TRUE) {
    source_type <- match.arg(source_type)
    peak_value_type <- match.arg(peak_value_type)
    prepared <- .condition_prepare_common_input(
        object = object, genes = genes,
        peak_to_gene_method = peak_to_gene_method,
        upstream = upstream, downstream = downstream, extend = extend,
        only_tss = only_tss,
        peak_to_gene_domains = peak_to_gene_domains,
        rna_layer = rna_layer, peak_layer = peak_layer,
        peak_value_type = peak_value_type, verbose = verbose
    )
    if (is.null(cells)) cells <- rownames(prepared$gene_data)
    .condition_discover_edges_prepared(
        prepared = prepared, cells = cells, source_label = source_label,
        source_type = source_type, tf_cor = tf_cor, peak_cor = peak_cor,
        parallel = parallel, verbose = verbose
    )
}

#' Union exact Pando TF-peak-target candidate triples
#'
#' @param global_edges Optional global candidate table.
#' @param condition_edges Named list of condition candidate tables.
#' @return A common `PandoEdgeDictionary`.
#' @export
union_grn_edges <- function(global_edges = NULL, condition_edges) {
    if (!is.list(condition_edges) || !length(condition_edges) ||
        is.null(names(condition_edges)) || any(!nzchar(names(condition_edges)))) {
        stop("`condition_edges` must be a non-empty named list.",
             call. = FALSE)
    }
    inputs <- condition_edges
    if (!is.null(global_edges)) inputs <- c(list(global = global_edges), inputs)
    required <- c("target", "tf", "region", "atac_feature_id")
    if (any(vapply(inputs, function(x) {
        !is.data.frame(x) || !all(required %in% colnames(x))
    }, logical(1)))) {
        stop("Every candidate table requires target, tf, region and ATAC IDs.",
             call. = FALSE)
    }
    provenance_fields <- c(
        "rna_layer", "peak_layer", "peak_value_type",
        "preprocessing_fingerprint", "dictionary_input_schema"
    )
    provenance <- lapply(inputs, function(candidate) {
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
        stop(
            "Candidate tables mix verified and unverified preprocessing ",
            "provenance.", call. = FALSE
        )
    }
    if (all(complete_provenance)) {
        reference <- provenance[[1L]]
        same_reference <- vapply(provenance, function(value) {
            identical(value, reference)
        }, logical(1))
        if (!all(same_reference)) {
            stop(
                "Global and condition candidates use different RNA/ATAC ",
                "layers, value semantics, or preprocessing fingerprints.",
                call. = FALSE
            )
        }
    }

    all_rows <- do.call(rbind, lapply(seq_along(inputs), function(i) {
        x <- as.data.frame(inputs[[i]], stringsAsFactors = FALSE)
        if (!nrow(x)) return(NULL)
        x$.union_source <- names(inputs)[[i]]
        x
    }))
    if (is.null(all_rows) || !nrow(all_rows)) {
        stop("Candidate discovery returned no TF-peak-target edge.",
             call. = FALSE)
    }
    key <- paste(all_rows$target, all_rows$tf, all_rows$region, sep = "\001")
    groups <- split(seq_len(nrow(all_rows)), key)
    dictionary <- do.call(rbind, lapply(groups, function(index) {
        one <- all_rows[index, , drop = FALSE]
        out <- one[1L, c("target", "tf", "region", "atac_feature_id"),
                   drop = FALSE]
        source <- unique(as.character(one$.union_source))
        condition_source <- sort(setdiff(source, "global"))
        out$source_global <- "global" %in% source
        out$source_conditions <- paste(condition_source, collapse = ";")
        out$n_sources <- length(source)
        safe_max_abs <- function(value) {
            value <- abs(as.numeric(value))
            value <- value[is.finite(value)]
            if (length(value)) max(value) else NA_real_
        }
        out$max_abs_peak_target_cor <- if ("peak_target_cor" %in% colnames(one)) {
            safe_max_abs(one$peak_target_cor)
        } else NA_real_
        out$max_abs_tf_target_cor <- if ("tf_target_cor" %in% colnames(one)) {
            safe_max_abs(one$tf_target_cor)
        } else NA_real_
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
        stop("Exact edge union produced duplicated triples.", call. = FALSE)
    }
    class(dictionary) <- c("PandoEdgeDictionary", "data.frame")
    attr(dictionary, "preprocessing_provenance_verified") <-
        all(complete_provenance)
    if (all(complete_provenance)) {
        for (field in provenance_fields) {
            attr(dictionary, field) <- provenance[[1L]][[field]]
        }
    }
    dictionary
}

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
    dictionary_fingerprint <- attr(
        dictionary, "preprocessing_fingerprint", exact = TRUE
    )
    dictionary_rna_layer <- attr(dictionary, "rna_layer", exact = TRUE)
    dictionary_peak_layer <- attr(dictionary, "peak_layer", exact = TRUE)
    dictionary_peak_value_type <- attr(
        dictionary, "peak_value_type", exact = TRUE
    )
    provenance_values <- list(
        dictionary_fingerprint, dictionary_rna_layer,
        dictionary_peak_layer, dictionary_peak_value_type
    )
    has_provenance <- vapply(provenance_values, function(value) {
        is.character(value) && length(value) == 1L &&
            !is.na(value) && nzchar(value)
    }, logical(1))
    if (any(has_provenance) && !all(has_provenance)) {
        stop("The dictionary contains incomplete preprocessing provenance.",
             call. = FALSE)
    }
    if (all(has_provenance) &&
        (!identical(dictionary_fingerprint,
                    prepared$preprocessing_fingerprint) ||
         !identical(dictionary_rna_layer, prepared$rna_layer) ||
         !identical(dictionary_peak_layer, prepared$peak_layer) ||
         !identical(dictionary_peak_value_type,
                    prepared$peak_value_type))) {
        stop(
            "The frozen dictionary and fixed fit use different RNA/ATAC ",
            "preprocessing references.", call. = FALSE
        )
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

.condition_standard_by_cell_type <- function(
    object, metadata, cell_type_col, cell_type, genes, network_name,
    peak_to_gene_method, upstream, downstream, extend, only_tss,
    parallel, tf_cor, peak_cor, fallback_args, verbose, overwrite) {
    if (is.null(cell_type_col)) {
        args <- list(
            object = object, genes = genes, network_name = network_name,
            peak_to_gene_method = peak_to_gene_method,
            upstream = upstream, downstream = downstream, extend = extend,
            only_tss = only_tss, parallel = parallel,
            tf_cor = tf_cor, peak_cor = peak_cor, scale = FALSE,
            interaction_term = ":", method = "glm", verbose = verbose
        )
        conflict <- intersect(names(fallback_args), names(args))
        if (length(conflict)) {
            stop("`fallback_args` cannot override: ",
                 paste(conflict, collapse = ", "), call. = FALSE)
        }
        return(list(
            object = do.call(infer_grn, c(args, fallback_args)),
            network_index = data.frame(
                cell_type = NA_character_, network_name = network_name,
                stringsAsFactors = FALSE
            )
        ))
    }
    .condition_validate_labels(metadata, cell_type_col)
    available <- unique(as.character(metadata[[cell_type_col]]))
    requested <- if (is.null(cell_type)) available else unique(as.character(cell_type))
    missing <- setdiff(requested, available)
    if (length(missing)) {
        stop("Requested cell type(s) were not found: ",
             paste(missing, collapse = ", "), call. = FALSE)
    }
    index <- list()
    for (label in requested) {
        cells <- rownames(metadata)[as.character(metadata[[cell_type_col]]) == label]
        one <- object
        one@data <- subset(object@data, cells = cells)
        one_name <- paste0(
            network_name, "__", .condition_safe_label(label), "__standard"
        )
        if (one_name %in% names(object@grn@networks) && !isTRUE(overwrite)) {
            stop("Network `", one_name, "` already exists.", call. = FALSE)
        }
        args <- list(
            object = one, genes = genes, network_name = one_name,
            peak_to_gene_method = peak_to_gene_method,
            upstream = upstream, downstream = downstream, extend = extend,
            only_tss = only_tss, parallel = parallel,
            tf_cor = tf_cor, peak_cor = peak_cor, scale = FALSE,
            interaction_term = ":", method = "glm", verbose = verbose
        )
        conflict <- intersect(names(fallback_args), names(args))
        if (length(conflict)) {
            stop("`fallback_args` cannot override: ",
                 paste(conflict, collapse = ", "), call. = FALSE)
        }
        one <- do.call(infer_grn, c(args, fallback_args))
        object@grn@networks[[one_name]] <- GetNetwork(one, one_name)
        object@grn@active_network <- one_name
        index[[length(index) + 1L]] <- data.frame(
            cell_type = label, network_name = one_name,
            n_cells = length(cells), stringsAsFactors = FALSE
        )
    }
    list(object = object, network_index = do.call(rbind, index))
}

#' Infer comparable condition-specific Pando GRNs
#'
#' With two or more conditions in a cell type, Pando candidate discovery is run
#' on all eligible-condition cells pooled together and independently in every
#' condition. Exact TF-peak-target triples passing both configured correlation
#' gates in either scope are deduplicated into one frozen common dictionary.
#' Production coefficients are estimated jointly by E-star at the single fixed
#' deviation threshold z=0.25. E-star fusion/boundary metadata are retained only
#' as properties of the production estimator. Formal significance is computed
#' separately, without fusion, from condition-local Gaussian linear models on the
#' same frozen target dictionaries. Each exact edge receives one omnibus P value
#' across its estimable conditions, followed by one BH correction across all
#' exact edges in the broad-cell-type network. The resulting edge topology is
#' common to all retained conditions; condition specificity remains in each
#' condition's continuous E-star coefficient. With fewer than two usable
#' conditions, standard Pando is used instead.
#'
#' @param object A `GRNData` object after motif matching.
#' @param cell_type_col Broad cell-type metadata column.
#' @param condition_col Condition metadata column; `NULL`, absent or one level
#'   selects direct per-cell-type Pando.
#' @param cell_type Optional cell-type subset.
#' @param genes Target genes.
#' @param network_name Prefix for generated networks.
#' @param tf_cor,peak_cor Candidate-discovery thresholds. Both default to 0.05;
#'   user-supplied values are used unchanged for pooled/global and condition
#'   discovery after ordinary range validation.
#' @param min_cells_per_condition Minimum cells retained per condition.
#' @param small_condition_action Error, drop the condition, or skip the cell type.
#' @param adjust_method,padj_threshold Exact-edge inference requires BH. A raw
#'   edge P value is built from the estimable condition-local no-fusion tests,
#'   and BH is applied once across all estimable exact edges in the cell type.
#'   Edge support requires strict `edge_padj < padj_threshold`.
#' @param rank_action,min_residual_df Identifiable-subspace and residual-degree-
#'   of-freedom controls for production and no-fusion inference.
#' @param reference_condition Optional predefined reference condition for the
#'   K-condition E-star contrast-tree geometry. It must be present and retained
#'   in every fitted cell type. When `NULL`, the first retained condition is used.
#'   This is a production-model coordinate, not an inference tuning parameter.
#' @param BPPARAM Optional BiocParallel parameter.
#' @param parallel_scope Automatic, cell-type, or target-level parallel scope.
#' @param fallback_args Arguments used only by standard Pando fallback. The
#'   conditional route does not accept ridge-CV/lambda, alternative-z, or
#'   fusion-ratio controls.
#' @param ... Must be empty.
#' @return A `GRNData` object with a frozen common dictionary, E-star z=0.25
#'   production coefficients, separate no-fusion inference, and one common
#'   exact-edge topology for all retained conditions.
#' @export
infer_condition_grn <- function(object, ...) {
    UseMethod(generic = "infer_condition_grn", object = object)
}

#' @rdname infer_condition_grn
#' @method infer_condition_grn GRNData
#' @export
infer_condition_grn.GRNData <- function(
    object, cell_type_col = NULL, condition_col = NULL, cell_type = NULL,
    genes = NULL, network_name = "condition_grn",
    peak_to_gene_method = c("Signac", "GREAT"), upstream = 100000,
    downstream = 0, extend = 1000000, only_tss = FALSE,
    peak_to_gene_domains = NULL, rna_layer = "data", peak_layer = "data",
    peak_value_type = c("normalized", "probability", "other"),
    tf_cor = 0.05, peak_cor = 0.05,
    min_cells_per_condition = 50L,
    small_condition_action = c("error", "drop_condition", "skip_cell_type"),
    adjust_method = "BH", padj_threshold = 0.05,
    rank_action = c("mark", "error"), min_residual_df = 1L,
    reference_condition = NULL,
    parallel = FALSE, BPPARAM = NULL,
    parallel_scope = c("auto", "cell_type", "target"),
    overwrite = FALSE, fallback_args = list(), verbose = TRUE, ...) {
    dots <- list(...)
    if (length(dots)) {
        label <- names(dots)
        if (is.null(label)) label <- rep("<unnamed>", length(dots))
        label[!nzchar(label)] <- "<unnamed>"
        stop("Unused condition-GRN argument(s): ",
             paste(label, collapse = ", "), call. = FALSE)
    }
    .pando_validate_bpparam(BPPARAM)
    parallel_scope <- match.arg(parallel_scope)

    use_pool <- isTRUE(parallel) && !identical(BPPARAM, FALSE) &&
        !is.null(BPPARAM)
    started_here <- FALSE
    old_target_param <- getOption("Pando.condition_target_BPPARAM", NULL)
    on.exit(
        options(Pando.condition_target_BPPARAM = old_target_param),
        add = TRUE
    )
    if (use_pool) {
        if (!isTRUE(BiocParallel::bpisup(BPPARAM))) {
            BPPARAM <- BiocParallel::bpstart(BPPARAM)
            started_here <- TRUE
        }
        options(Pando.condition_target_BPPARAM = BPPARAM)
        if (started_here) {
            on.exit({
                try(BiocParallel::bpstop(BPPARAM), silent = TRUE)
                invisible(gc(verbose = FALSE, full = TRUE))
            }, add = TRUE)
        }
    }

    metadata <- object@data@meta.data
    can_split <- isTRUE(parallel) &&
        !identical(parallel_scope, "target") &&
        is.character(cell_type_col) && length(cell_type_col) == 1L &&
        !is.na(cell_type_col) && cell_type_col %in% colnames(metadata)
    available <- if (can_split) {
        unique(as.character(metadata[[cell_type_col]]))
    } else character()
    requested <- if (!can_split) {
        character()
    } else if (is.null(cell_type)) {
        available
    } else {
        unique(as.character(cell_type))
    }

    run_one <- function(input_object, type_label = cell_type,
                        inner_parallel = parallel) {
        .pando_infer_condition_grn_multitask_ridge_one(
            object = input_object,
            cell_type_col = cell_type_col,
            condition_col = condition_col,
            cell_type = type_label,
            genes = genes,
            network_name = network_name,
            peak_to_gene_method = peak_to_gene_method,
            upstream = upstream,
            downstream = downstream,
            extend = extend,
            only_tss = only_tss,
            peak_to_gene_domains = peak_to_gene_domains,
            rna_layer = rna_layer,
            peak_layer = peak_layer,
            peak_value_type = peak_value_type,
            tf_cor = tf_cor,
            peak_cor = peak_cor,
            min_cells_per_condition = min_cells_per_condition,
            small_condition_action = small_condition_action,
            adjust_method = adjust_method,
            padj_threshold = padj_threshold,
            rank_action = rank_action,
            min_residual_df = min_residual_df,
            reference_condition = reference_condition,
            parallel = inner_parallel,
            overwrite = overwrite,
            fallback_args = fallback_args,
            verbose = verbose
        )
    }

    if (!can_split || length(requested) <= 1L) {
        answer <- run_one(object)
        if (use_pool && inherits(answer, "GRNData")) {
            answer@grn@params$condition_target_parallel_plan <- list(
                scope = "target",
                workers = as.integer(BiocParallel::bpnworkers(BPPARAM)),
                backend = class(BPPARAM)[[1L]],
                pool_reused_across_target_stages = TRUE,
                nested_pool_created_by_pando = started_here,
                pool_release_policy = "stop_on_method_exit_if_started_here"
            )
        }
        return(answer)
    }

    missing <- setdiff(requested, available)
    if (length(missing)) {
        stop("Requested cell type(s) were not found: ",
             paste(missing, collapse = ", "), call. = FALSE)
    }
    run_type <- function(type_label) {
        cells <- rownames(metadata)[
            as.character(metadata[[cell_type_col]]) == type_label
        ]
        one <- object
        one@data <- subset(object@data, cells = cells)
        run_one(one, type_label = type_label, inner_parallel = FALSE)
    }
    results <- .pando_parallel_lapply(
        requested, run_type, parallel = TRUE, BPPARAM = BPPARAM
    )
    answer <- .pando_merge_cell_type_grn_results(
        object, results, condition_col = condition_col,
        cell_type_col = cell_type_col
    )
    answer@grn@params$parallel_plan <- list(
        scope = "cell_type", n_jobs = length(requested),
        nested_parallel = FALSE, cell_types = requested
    )
    answer
}

#' Extract common-dictionary condition GRN fits
#'
#' @param object A fitted `GRNData` object.
#' @param cell_type Optional fitted cell type.
#' @return One `ConditionGRNFit` or a named list.
#' @export
condition_grn_fit <- function(object, ...) {
    UseMethod(generic = "condition_grn_fit", object = object)
}

#' @rdname condition_grn_fit
#' @method condition_grn_fit GRNData
#' @export
condition_grn_fit.GRNData <- function(
    object, cell_type = NULL, ...) {
    dots <- list(...)
    if (length(dots)) {
        label <- names(dots)
        if (is.null(label)) label <- rep("<unnamed>", length(dots))
        label[!nzchar(label)] <- "<unnamed>"
        stop(
            "Unused condition_grn_fit argument(s): ",
            paste(label, collapse = ", "), call. = FALSE
        )
    }
    fits <- object@grn@params$condition_grn_fits
    if (!is.list(fits) || !length(fits)) {
        stop("No common-dictionary condition GRN fit is stored.",
             call. = FALSE)
    }
    if (is.null(cell_type)) return(fits)
    cell_type <- as.character(cell_type)
    missing <- setdiff(cell_type, names(fits))
    if (length(missing)) {
        stop("Condition GRN fit was not found for: ",
             paste(missing, collapse = ", "), call. = FALSE
    }
    answer <- fits[cell_type]
    if (length(answer) == 1L) answer[[1L]] else answer
}

#' Extract one condition subgraph
#'
#' @param fit A `ConditionGRNFit`.
#' @param condition Fitted condition label.
#' @param significant_only If `TRUE`, return only exact edges supported by the
#'   whole-network edge-level BH topology, retaining the selected condition's own
#'   continuous `penalty_effect`. Because topology is common, every condition
#'   returns the same supported exact-edge IDs.
#' @return Edge table for the selected condition.
#' @export
condition_grn_subgraph <- function(fit, condition, significant_only = TRUE) {
    if (!inherits(fit, "ConditionGRNFit") ||
        !identical(fit$schema_version, .condition_common_dictionary_schema)) {
        stop("`fit` is not a common-dictionary ConditionGRNFit.",
             call. = FALSE)
    }
    condition <- as.character(condition)[[1L]]
    if (!condition %in% fit$condition_levels) {
        stop("Unknown fitted condition: ", condition, call. = FALSE)
    }
    answer <- fit$coefficients[
        as.character(fit$coefficients$condition) == condition, , drop = FALSE
    ]
    if (isTRUE(significant_only)) {
        if (!"active_in_regcompass" %in% colnames(answer)) {
            stop("The condition fit lacks common exact-edge topology flags.",
                 call. = FALSE)
        }
        answer <- answer[
            answer$active_in_regcompass %in% TRUE, , drop = FALSE
        ]
    }
    answer
}

#' Project common-dictionary condition effects on paired cells
#'
#' Reconstructs TF RNA multiplied by peak ATAC on the original unscaled input and
#' applies the continuous condition-specific E-star production coefficient.
#' With `significant_only = TRUE`, projection first applies the common exact-edge
#' topology defined by whole-network edge-level BH and then retains every
#' condition's own continuous `penalty_effect`.
#'
#' @param object Fitted `GRNData` object.
#' @param fit A `ConditionGRNFit`.
#' @param targets Optional target subset.
#' @param significant_only Use the common exact-edge BH topology before applying
#'   continuous `penalty_effect`; if `FALSE`, project all finite estimates.
#' @param return_edge_contributions Return the cell-by-edge matrix.
#' @return A `PandoConditionProjection` list.
#' @export

#' Aggregate paired-cell condition projection
#'
#' @param projection A `PandoConditionProjection`.
#' @param membership Data frame with `cell_id` and a grouping column.
#' @param group_col Grouping column, typically `metacell_id`.
#' @return Aggregated gene scores and provenance.
#' @export
aggregate_condition_grn_projection <- function(
    projection, membership, group_col = "metacell_id") {
    if (!inherits(projection, "PandoConditionProjection") ||
        !is.data.frame(membership) ||
        !all(c("cell_id", group_col) %in% colnames(membership))) {
        stop("Projection and membership inputs are invalid.", call. = FALSE)
    }
    membership$cell_id <- as.character(membership$cell_id)
    membership[[group_col]] <- as.character(membership[[group_col]])
    if (anyDuplicated(membership$cell_id)) {
        stop("Every cell must map to exactly one aggregation group.",
             call. = FALSE)
    }
    cells <- intersect(rownames(projection$gene_score), membership$cell_id)
    if (!length(cells)) {
        stop("Projection and membership share no cells.", call. = FALSE)
    }
    membership <- membership[match(cells, membership$cell_id), , drop = FALSE]
    groups <- unique(membership[[group_col]])
    gene_score <- vapply(groups, function(group) {
        rows <- membership[[group_col]] == group
        colMeans(projection$gene_score[cells[rows], , drop = FALSE])
    }, numeric(ncol(projection$gene_score)))
    if (is.null(dim(gene_score))) {
        gene_score <- matrix(gene_score, ncol = 1L)
    }
    rownames(gene_score) <- colnames(projection$gene_score)
    colnames(gene_score) <- groups
    list(
        gene_score = gene_score,
        group_col = group_col,
        source_projection = projection,
        aggregation = "arithmetic_mean_of_paired_cell_regulatory_scores"
    )
}
