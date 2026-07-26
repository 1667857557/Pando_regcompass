# Explicit matrix alignment for the shared-design builder. This file is loaded
# after grn_design.R and before the versioned contract wrappers.

.pando_align_design_matrices <- function(
    peaks2gene, peaks2motif, peak_data, motif2tf) {
    region_sets <- list(
        colnames(peaks2gene),
        rownames(peaks2motif),
        colnames(peak_data)
    )
    if (any(vapply(region_sets, is.null, logical(1)))) {
        stop('Candidate matrices require explicit regulatory-region IDs.',
             call. = FALSE)
    }
    common_regions <- Reduce(intersect, lapply(region_sets, as.character))
    if (!length(common_regions)) {
        stop('Peak-to-gene, motif, and ATAC matrices share no regulatory regions.',
             call. = FALSE)
    }
    common_regions <- as.character(colnames(peaks2gene))[
        as.character(colnames(peaks2gene)) %in% common_regions
    ]
    if (anyDuplicated(common_regions)) {
        stop('Regulatory-region IDs must be unique before candidate alignment.',
             call. = FALSE)
    }
    peaks2gene <- peaks2gene[, common_regions, drop = FALSE]
    peaks2motif <- peaks2motif[common_regions, , drop = FALSE]
    peak_data <- peak_data[, common_regions, drop = FALSE]

    motif_ids <- intersect(
        as.character(colnames(peaks2motif)),
        as.character(rownames(motif2tf))
    )
    if (!length(motif_ids)) {
        stop('Motif matches and motif-to-TF mappings share no motif IDs.',
             call. = FALSE)
    }
    motif_ids <- as.character(colnames(peaks2motif))[
        as.character(colnames(peaks2motif)) %in% motif_ids
    ]
    if (anyDuplicated(motif_ids)) {
        stop('Motif IDs must be unique before candidate alignment.',
             call. = FALSE)
    }
    peaks2motif <- peaks2motif[, motif_ids, drop = FALSE]
    motif2tf <- motif2tf[motif_ids, , drop = FALSE]

    list(
        peaks2gene = peaks2gene,
        peaks2motif = peaks2motif,
        peak_data = peak_data,
        motif2tf = motif2tf
    )
}

#' @rdname prepare_grn_design
#' @method prepare_grn_design GRNData
#' @export
prepare_grn_design.GRNData <- function(
    object,
    genes = NULL,
    peak_to_gene_method = c('Signac', 'GREAT'),
    upstream = 100000,
    downstream = 0,
    extend = 1000000,
    only_tss = FALSE,
    peak_to_gene_domains = NULL,
    min_tf_detection = 0,
    min_peak_detection = 0,
    min_target_detection = 0,
    max_edges_per_target = Inf,
    verbose = TRUE,
    ...
) {
    peak_to_gene_method <- match.arg(peak_to_gene_method)
    validate_fraction <- function(x, name) {
        if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
            x < 0 || x > 1) {
            stop('`', name, '` must be one finite number in [0, 1].',
                 call. = FALSE)
        }
        as.numeric(x)
    }
    min_tf_detection <- validate_fraction(
        min_tf_detection, 'min_tf_detection'
    )
    min_peak_detection <- validate_fraction(
        min_peak_detection, 'min_peak_detection'
    )
    min_target_detection <- validate_fraction(
        min_target_detection, 'min_target_detection'
    )
    if (!is.numeric(max_edges_per_target) ||
        length(max_edges_per_target) != 1L ||
        is.na(max_edges_per_target) || max_edges_per_target <= 0 ||
        (is.finite(max_edges_per_target) &&
         abs(max_edges_per_target - round(max_edges_per_target)) >
             sqrt(.Machine$double.eps))) {
        stop(
            '`max_edges_per_target` must be one positive integer or `Inf`.',
            call. = FALSE
        )
    }
    if (is.finite(max_edges_per_target)) {
        max_edges_per_target <- as.integer(max_edges_per_target)
    }

    params <- Params(object)
    motif2tf <- NetworkTFs(object)
    if (is.null(motif2tf)) {
        stop('Motif matches have not been found. Run `find_motifs()` first.',
             call. = FALSE)
    }
    gene_annot <- Signac::Annotation(GetAssay(object, params$peak_assay))
    if (is.null(gene_annot)) {
        stop('Please provide a gene annotation for the ChromatinAssay.',
             call. = FALSE)
    }
    if (is.null(genes)) {
        genes <- VariableFeatures(object, assay = params$rna_assay)
        if (is.null(genes) || !length(genes)) {
            stop('Please provide target genes or run `FindVariableFeatures()`.',
                 call. = FALSE)
        }
    }
    genes <- unique(as.character(genes))

    gene_data <- Matrix::t(LayerData(
        object, assay = params$rna_assay, layer = 'data'
    ))
    peak_data_all <- Matrix::t(LayerData(
        object, assay = params$peak_assay, layer = 'data'
    ))
    if (!identical(rownames(gene_data), rownames(peak_data_all))) {
        common_cells <- intersect(
            rownames(gene_data), rownames(peak_data_all)
        )
        if (!length(common_cells)) {
            stop('RNA and ATAC assays do not share cells.', call. = FALSE)
        }
        gene_data <- gene_data[common_cells, , drop = FALSE]
        peak_data_all <- peak_data_all[common_cells, , drop = FALSE]
    }

    features <- intersect(gene_annot$gene_name, genes)
    features <- intersect(features, colnames(gene_data))
    features <- unique(as.character(features))
    gene_annot <- gene_annot[gene_annot$gene_name %in% features, ]
    if (!length(features)) {
        stop('No requested targets overlap RNA features and gene annotations.',
             call. = FALSE)
    }

    regions <- NetworkRegions(object)
    if (!length(regions@peaks) || !nrow(regions@motifs@data)) {
        stop('Pando regions do not contain peak mappings and motif matches.',
             call. = FALSE)
    }
    if (any(regions@peaks < 1 | regions@peaks > ncol(peak_data_all))) {
        stop('Pando region-to-ATAC peak indices are out of range.',
             call. = FALSE)
    }
    peak_data <- peak_data_all[, regions@peaks, drop = FALSE]
    atac_feature_id <- colnames(peak_data)
    peaks2motif <- regions@motifs@data
    region_id <- rownames(peaks2motif)
    if (is.null(region_id) || length(region_id) != ncol(peak_data) ||
        anyNA(region_id) || any(!nzchar(region_id)) ||
        anyDuplicated(region_id)) {
        stop('Pando motif regions must have unique, complete IDs aligned to peaks.',
             call. = FALSE)
    }
    colnames(peak_data) <- region_id

    log_message(
        'Selecting structural candidate regulatory regions',
        verbose = verbose
    )
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
            method = 'Signac',
            genes = peak_to_gene_domains,
            upstream = 0,
            downstream = 0,
            extend = extend,
            only_tss = FALSE
        )
    }
    peaks2gene <- aggregate_matrix(
        t(peaks_near_gene),
        groups = colnames(peaks_near_gene),
        fun = 'sum'
    )

    tfs <- intersect(colnames(motif2tf), colnames(gene_data))
    motif2tf <- motif2tf[, tfs, drop = FALSE]
    aligned <- .pando_align_design_matrices(
        peaks2gene = peaks2gene,
        peaks2motif = peaks2motif,
        peak_data = peak_data,
        motif2tf = motif2tf
    )
    peaks2gene <- aligned$peaks2gene
    peaks2motif <- aligned$peaks2motif
    peak_data <- aligned$peak_data
    motif2tf <- aligned$motif2tf

    peaks_at_gene <- as.logical(
        sparseMatrixStats::colMaxs(peaks2gene)
    )
    peaks_with_motif <- as.logical(
        sparseMatrixStats::rowMaxs(peaks2motif * 1)
    )
    peaks_use <- peaks_at_gene & peaks_with_motif
    peaks2gene <- peaks2gene[, peaks_use, drop = FALSE]
    peaks2motif <- peaks2motif[peaks_use, , drop = FALSE]
    peak_data <- peak_data[, peaks_use, drop = FALSE]

    region_lookup <- stats::setNames(atac_feature_id, region_id)
    region_map <- data.frame(
        region = colnames(peak_data),
        atac_feature_id = unname(region_lookup[colnames(peak_data)]),
        stringsAsFactors = FALSE
    )
    if (anyNA(region_map$atac_feature_id)) {
        stop('Aligned regulatory regions lack measured ATAC feature IDs.',
             call. = FALSE)
    }
    peak_detection <- Matrix::colMeans(peak_data > 0)
    names(peak_detection) <- colnames(peak_data)
    target_detection <- Matrix::colMeans(
        gene_data[, features, drop = FALSE] > 0
    )
    names(target_detection) <- features
    tf_detection <- Matrix::colMeans(
        gene_data[, tfs, drop = FALSE] > 0
    )
    names(tf_detection) <- tfs

    build_target <- function(g) {
        if (!g %in% rownames(peaks2gene)) return(data.frame())
        if (!is.finite(target_detection[[g]]) ||
            target_detection[[g]] < min_target_detection) {
            return(data.frame())
        }
        gene_peaks <- as.logical(peaks2gene[g, ])
        selected_regions <- colnames(peaks2gene)[gene_peaks]
        selected_regions <- selected_regions[
            is.finite(peak_detection[selected_regions]) &
                peak_detection[selected_regions] >= min_peak_detection
        ]
        if (!length(selected_regions)) return(data.frame())
        rows <- lapply(selected_regions, function(region) {
            motif_present <- as.logical(
                peaks2motif[region, , drop = TRUE]
            )
            if (!any(motif_present)) return(NULL)
            peak_tfs <- sparseMatrixStats::colMaxs(
                motif2tf[motif_present, , drop = FALSE]
            )
            peak_tfs <- colnames(motif2tf)[as.logical(peak_tfs)]
            peak_tfs <- setdiff(peak_tfs, g)
            peak_tfs <- peak_tfs[
                is.finite(tf_detection[peak_tfs]) &
                    tf_detection[peak_tfs] >= min_tf_detection
            ]
            if (!length(peak_tfs)) return(NULL)
            data.frame(
                tf = peak_tfs,
                target = g,
                region = region,
                atac_feature_id = region_map$atac_feature_id[
                    match(region, region_map$region)
                ],
                tf_feature_id = peak_tfs,
                target_feature_id = g,
                motif_supported = TRUE,
                peak_to_gene_supported = TRUE,
                tf_detection = as.numeric(tf_detection[peak_tfs]),
                peak_detection = as.numeric(peak_detection[[region]]),
                target_detection = as.numeric(target_detection[[g]]),
                stringsAsFactors = FALSE
            )
        })
        rows <- rows[!vapply(rows, is.null, logical(1))]
        if (!length(rows)) return(data.frame())
        out <- unique(do.call(rbind, rows))
        out[order(out$region, out$tf), , drop = FALSE]
    }

    candidate_rows <- lapply(features, build_target)
    candidate_rows <- candidate_rows[
        vapply(candidate_rows, nrow, integer(1)) > 0L
    ]
    if (length(candidate_rows)) {
        candidate_edges <- unique(do.call(rbind, candidate_rows))
        candidate_edges <- candidate_edges[order(
            candidate_edges$target,
            candidate_edges$tf,
            candidate_edges$atac_feature_id,
            candidate_edges$region
        ), , drop = FALSE]
        predictor_key <- paste(
            candidate_edges$tf,
            candidate_edges$atac_feature_id,
            candidate_edges$target,
            sep = '\001'
        )
        predictor_rows <- split(seq_len(nrow(candidate_edges)), predictor_key)
        candidate_edges <- do.call(rbind, lapply(
            predictor_rows, function(index) {
                one <- candidate_edges[index, , drop = FALSE]
                out <- one[1L, , drop = FALSE]
                support <- sort(unique(as.character(one$region)))
                out$region <- support[[1L]]
                out$supporting_regions <- paste(support, collapse = ';')
                out$n_supporting_regions <- length(support)
                out$tf_detection <- max(one$tf_detection, na.rm = TRUE)
                out$peak_detection <- max(one$peak_detection, na.rm = TRUE)
                out$target_detection <- max(one$target_detection, na.rm = TRUE)
                out
            }
        ))
        rownames(candidate_edges) <- NULL
        if (is.finite(max_edges_per_target)) {
            target_rows <- split(
                seq_len(nrow(candidate_edges)), candidate_edges$target
            )
            candidate_edges <- do.call(rbind, lapply(
                target_rows, function(index) {
                    index <- index[seq_len(min(
                        length(index), max_edges_per_target
                    ))]
                    candidate_edges[index, , drop = FALSE]
                }
            ))
            rownames(candidate_edges) <- NULL
        }
    } else {
        candidate_edges <- data.frame(
            tf = character(), target = character(), region = character(),
            atac_feature_id = character(), tf_feature_id = character(),
            target_feature_id = character(), motif_supported = logical(),
            peak_to_gene_supported = logical(), supporting_regions = character(),
            n_supporting_regions = integer(), tf_detection = numeric(),
            peak_detection = numeric(), target_detection = numeric(),
            stringsAsFactors = FALSE
        )
    }
    if (nrow(candidate_edges)) {
        candidate_edges$edge_id <- paste(
            candidate_edges$tf,
            candidate_edges$atac_feature_id,
            candidate_edges$target,
            sep = '::'
        )
        candidate_edges$candidate_index <- seq_len(nrow(candidate_edges))
        candidate_edges <- candidate_edges[, c(
            'edge_id', 'candidate_index', 'tf', 'region', 'target',
            'atac_feature_id', 'tf_feature_id', 'target_feature_id',
            'motif_supported', 'peak_to_gene_supported', 'supporting_regions',
            'n_supporting_regions', 'tf_detection', 'peak_detection',
            'target_detection'
        ), drop = FALSE]
    }

    target_diagnostics <- data.frame(
        target = features,
        target_detection = as.numeric(target_detection[features]),
        n_candidate_edges = vapply(features, function(g) {
            sum(candidate_edges$target == g)
        }, integer(1)),
        n_candidate_tfs = vapply(features, function(g) {
            length(unique(candidate_edges$tf[candidate_edges$target == g]))
        }, integer(1)),
        n_candidate_regions = vapply(features, function(g) {
            length(unique(candidate_edges$region[candidate_edges$target == g]))
        }, integer(1)),
        stringsAsFactors = FALSE
    )
    design_key <- if (nrow(candidate_edges)) {
        paste(candidate_edges$edge_id, collapse = '\n')
    } else {
        'empty'
    }
    answer <- list(
        schema_version = 'pando_grn_design_v1',
        candidate_edges = candidate_edges,
        region_map = region_map,
        target_diagnostics = target_diagnostics,
        feature_contract = list(
            cell_ids = rownames(gene_data),
            rna_feature_ids = colnames(gene_data),
            atac_feature_ids = colnames(peak_data_all),
            rna_assay = params$rna_assay,
            atac_assay = params$peak_assay,
            rna_layer = 'data',
            atac_layer = 'data'
        ),
        params = list(
            peak_to_gene_method = peak_to_gene_method,
            upstream = upstream,
            downstream = downstream,
            extend = extend,
            only_tss = only_tss,
            custom_peak_to_gene_domains = !is.null(peak_to_gene_domains),
            candidate_policy = 'structural_shared_before_model_fitting',
            min_tf_detection = min_tf_detection,
            min_peak_detection = min_peak_detection,
            min_target_detection = min_target_detection,
            max_edges_per_target = max_edges_per_target
        ),
        design_fingerprint = paste0(
            length(features), ':', nrow(candidate_edges), ':',
            sum(utf8ToInt(design_key)) %% 2147483647
        )
    )
    class(answer) <- c('PandoGRNDesign', 'list')
    answer
}
