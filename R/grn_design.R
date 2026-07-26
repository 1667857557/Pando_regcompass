#' Prepare a shared Pando TF-peak-target candidate design
#'
#' Builds the structural TF-peak-target edge universe used before model fitting.
#' The returned design is intentionally condition agnostic: callers may fit one
#' model across several experimental conditions while retaining one identical
#' candidate-edge dictionary. No coefficient, p-value, or pooled-expression
#' significance filter is applied by this function.
#'
#' @param object A `GRNData` object after [find_motifs()] has been run.
#' @param genes Target genes. Defaults to variable RNA features.
#' @param peak_to_gene_method Peak-to-gene domain method, `"Signac"` or `"GREAT"`.
#' @param upstream,downstream,extend,only_tss Regulatory-domain arguments passed
#'   to [find_peaks_near_genes()].
#' @param peak_to_gene_domains Optional custom regulatory domains. When supplied,
#'   these replace the distance-based domains.
#' @param min_tf_detection,min_peak_detection,min_target_detection Minimum
#'   fractions of cells with positive normalized signal. These filters are
#'   evaluated across the complete input object, not separately by condition.
#' @param max_edges_per_target Optional hard structural edge cap per target.
#'   The default, `Inf`, keeps every structurally valid edge. Finite caps retain
#'   edges in deterministic TF/region order and are intended only as a safety
#'   valve; data-dependent screening should be performed by the downstream
#'   multitask model using one shared rule across conditions.
#' @param verbose Display progress messages.
#' @param ... Reserved for future design builders.
#'
#' @return A `PandoGRNDesign` list containing `candidate_edges`, `region_map`,
#'   `target_diagnostics`, `feature_contract`, and design parameters.
#' @export
prepare_grn_design <- function(object, ...) {
    UseMethod(generic = 'prepare_grn_design', object = object)
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
        if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 0 || x > 1) {
            stop('`', name, '` must be one finite number in [0, 1].', call. = FALSE)
        }
        as.numeric(x)
    }
    min_tf_detection <- validate_fraction(min_tf_detection, 'min_tf_detection')
    min_peak_detection <- validate_fraction(min_peak_detection, 'min_peak_detection')
    min_target_detection <- validate_fraction(min_target_detection, 'min_target_detection')
    if (!is.numeric(max_edges_per_target) || length(max_edges_per_target) != 1L ||
        is.na(max_edges_per_target) || max_edges_per_target <= 0) {
        stop('`max_edges_per_target` must be one positive number or `Inf`.', call. = FALSE)
    }

    params <- Params(object)
    motif2tf <- NetworkTFs(object)
    if (is.null(motif2tf)) {
        stop('Motif matches have not been found. Run `find_motifs()` first.', call. = FALSE)
    }
    gene_annot <- Signac::Annotation(GetAssay(object, params$peak_assay))
    if (is.null(gene_annot)) {
        stop('Please provide a gene annotation for the ChromatinAssay.', call. = FALSE)
    }
    if (is.null(genes)) {
        genes <- VariableFeatures(object, assay = params$rna_assay)
        if (is.null(genes) || !length(genes)) {
            stop('Please provide target genes or run `FindVariableFeatures()`.', call. = FALSE)
        }
    }
    genes <- unique(as.character(genes))

    gene_data <- Matrix::t(LayerData(
        object,
        assay = params$rna_assay,
        layer = 'data'
    ))
    peak_data_all <- Matrix::t(LayerData(
        object,
        assay = params$peak_assay,
        layer = 'data'
    ))
    if (!identical(rownames(gene_data), rownames(peak_data_all))) {
        common_cells <- intersect(rownames(gene_data), rownames(peak_data_all))
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
        stop('No requested targets overlap RNA features and gene annotations.', call. = FALSE)
    }

    regions <- NetworkRegions(object)
    if (!length(regions@peaks) || !nrow(regions@motifs@data)) {
        stop('Pando regions do not contain peak mappings and motif matches.', call. = FALSE)
    }
    if (any(regions@peaks < 1 | regions@peaks > ncol(peak_data_all))) {
        stop('Pando region-to-ATAC peak indices are out of range.', call. = FALSE)
    }
    peak_data <- peak_data_all[, regions@peaks, drop = FALSE]
    atac_feature_id <- colnames(peak_data)
    region_id <- rownames(regions@motifs@data)
    if (length(region_id) != ncol(peak_data)) {
        stop('Pando motif regions and mapped ATAC features have different lengths.', call. = FALSE)
    }
    colnames(peak_data) <- region_id
    peaks2motif <- regions@motifs@data

    log_message('Selecting structural candidate regulatory regions', verbose = verbose)
    if (is.null(peak_to_gene_domains)) {
        peaks_near_gene <- find_peaks_near_genes(
            peaks = regions@ranges,
            method = peak_to_gene_method,
            genes = gene_annot,
            upstream = upstream,
            downstream = downstream,
            only_tss = only_tss
        )
    } else {
        peaks_near_gene <- find_peaks_near_genes(
            peaks = regions@ranges,
            method = 'Signac',
            genes = peak_to_gene_domains,
            upstream = 0,
            downstream = 0,
            only_tss = FALSE
        )
    }
    peaks2gene <- aggregate_matrix(
        t(peaks_near_gene),
        groups = colnames(peaks_near_gene),
        fun = 'sum'
    )
    peaks_at_gene <- as.logical(sparseMatrixStats::colMaxs(peaks2gene))
    peaks_with_motif <- as.logical(sparseMatrixStats::rowMaxs(peaks2motif * 1))
    peaks_use <- peaks_at_gene & peaks_with_motif
    peaks2gene <- peaks2gene[, peaks_use, drop = FALSE]
    peaks2motif <- peaks2motif[peaks_use, , drop = FALSE]
    peak_data <- peak_data[, peaks_use, drop = FALSE]

    region_map <- data.frame(
        region = region_id,
        atac_feature_id = atac_feature_id,
        stringsAsFactors = FALSE
    )
    region_map <- region_map[region_map$region %in% colnames(peak_data), , drop = FALSE]
    region_map <- region_map[match(colnames(peak_data), region_map$region), , drop = FALSE]
    peak_detection <- Matrix::colMeans(peak_data > 0)
    names(peak_detection) <- colnames(peak_data)
    target_detection <- Matrix::colMeans(gene_data[, features, drop = FALSE] > 0)
    names(target_detection) <- features

    tfs <- intersect(colnames(motif2tf), colnames(gene_data))
    motif2tf <- motif2tf[, tfs, drop = FALSE]
    tf_detection <- Matrix::colMeans(gene_data[, tfs, drop = FALSE] > 0)
    names(tf_detection) <- tfs

    build_target <- function(g) {
        if (!g %in% rownames(peaks2gene)) return(data.frame())
        if (!is.finite(target_detection[[g]]) || target_detection[[g]] < min_target_detection) {
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
            motif_present <- as.logical(peaks2motif[region, , drop = TRUE])
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
        out <- out[order(out$region, out$tf), , drop = FALSE]
        out
    }

    candidate_rows <- lapply(features, build_target)
    candidate_rows <- candidate_rows[vapply(candidate_rows, nrow, integer(1)) > 0L]
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
            sep = "\001"
        )
        predictor_rows <- split(seq_len(nrow(candidate_edges)), predictor_key)
        candidate_edges <- do.call(rbind, lapply(predictor_rows, function(index) {
            one <- candidate_edges[index, , drop = FALSE]
            out <- one[1L, , drop = FALSE]
            regions <- sort(unique(as.character(one$region)))
            out$region <- regions[[1L]]
            out$supporting_regions <- paste(regions, collapse = ";")
            out$n_supporting_regions <- length(regions)
            out$tf_detection <- max(one$tf_detection, na.rm = TRUE)
            out$peak_detection <- max(one$peak_detection, na.rm = TRUE)
            out$target_detection <- max(one$target_detection, na.rm = TRUE)
            out
        }))
        rownames(candidate_edges) <- NULL
        if (is.finite(max_edges_per_target)) {
            target_rows <- split(
                seq_len(nrow(candidate_edges)), candidate_edges$target
            )
            candidate_edges <- do.call(rbind, lapply(target_rows, function(index) {
                index <- index[seq_len(min(
                    length(index), as.integer(max_edges_per_target)
                ))]
                candidate_edges[index, , drop = FALSE]
            }))
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
    design_fingerprint <- paste0(
        length(features), ':', nrow(candidate_edges), ':',
        sum(utf8ToInt(design_key)) %% 2147483647
    )
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
        design_fingerprint = design_fingerprint
    )
    class(answer) <- c('PandoGRNDesign', 'list')
    answer
}

#' Validate a Pando shared candidate design
#'
#' @param design A `PandoGRNDesign` object.
#' @return `TRUE` invisibly; otherwise an error is raised.
#' @export
validate_grn_design <- function(design) {
    if (!inherits(design, 'PandoGRNDesign') || !is.list(design)) {
        stop('`design` must be a PandoGRNDesign object.', call. = FALSE)
    }
    if (!identical(design$schema_version, 'pando_grn_design_v1')) {
        stop('Unsupported Pando GRN design schema.', call. = FALSE)
    }
    edges <- design$candidate_edges
    required <- c(
        'edge_id', 'tf', 'region', 'target', 'atac_feature_id',
        'tf_feature_id', 'target_feature_id'
    )
    if (!is.data.frame(edges) || !all(required %in% colnames(edges))) {
        stop('Pando GRN design candidate edges are incomplete.', call. = FALSE)
    }
    if (nrow(edges) && (
        anyNA(edges$edge_id) || any(!nzchar(edges$edge_id)) ||
        anyDuplicated(edges$edge_id)
    )) {
        stop('Pando GRN design edge IDs must be unique and non-empty.', call. = FALSE)
    }
    contract <- design$feature_contract
    if (!is.list(contract) || !all(c(
        'cell_ids', 'rna_feature_ids', 'atac_feature_ids',
        'rna_assay', 'atac_assay'
    ) %in% names(contract))) {
        stop('Pando GRN design feature contract is incomplete.', call. = FALSE)
    }
    invisible(TRUE)
}
