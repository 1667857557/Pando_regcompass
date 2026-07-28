# Data preparation and target-level fitting for condition-aware Pando.

.condition_prepare_global_data <- function(
    object,
    cell_type_col,
    condition_col,
    genes,
    peak_to_gene_method,
    upstream,
    downstream,
    extend,
    only_tss,
    peak_to_gene_domains,
    verbose
) {
    params <- Params(object)
    motif2tf <- NetworkTFs(object)
    regions <- NetworkRegions(object)
    if (is.null(motif2tf) || length(motif2tf) == 0L) {
        stop('Motif matches have not been found. Please run find_motifs() first.')
    }
    if (is.null(regions) || is.null(regions@motifs)) {
        stop('Network regions and motif matches are required.')
    }
    gene_annotation <- Signac::Annotation(GetAssay(object, params$peak_assay))
    if (is.null(gene_annotation) || is.null(gene_annotation$gene_name)) {
        stop('The ChromatinAssay must contain gene annotation with gene_name.')
    }

    model_data <- .condition_get_model_data(
        object, params
    )
    gene_data <- model_data$gene_data
    peak_data <- model_data$peak_data
    metadata <- model_data$metadata

    if (is.null(genes)) {
        genes <- VariableFeatures(object, assay = params$rna_assay)
    }
    if (length(genes) == 0L) {
        stop('Please provide genes or run FindVariableFeatures() on the RNA assay.')
    }
    features <- Reduce(intersect, list(
        unique(as.character(gene_annotation$gene_name)),
        as.character(genes),
        colnames(gene_data)
    ))
    if (length(features) == 0L) {
        stop('No target genes overlap annotation, requested genes, and RNA data.')
    }
    gene_annotation <- gene_annotation[gene_annotation$gene_name %in% features]

    peaks2motif <- regions@motifs@data
    if (length(regions@peaks) != nrow(peaks2motif)) {
        stop('Network region peak indices and motif rows have different lengths.')
    }
    peak_data <- peak_data[, regions@peaks, drop = FALSE]
    colnames(peak_data) <- rownames(peaks2motif)

    log_message('Selecting candidate regulatory regions near genes', verbose = verbose)
    if (is.null(peak_to_gene_domains)) {
        peaks_near_gene <- find_peaks_near_genes(
            peaks = regions@ranges,
            genes = gene_annotation,
            method = peak_to_gene_method,
            upstream = upstream,
            downstream = downstream,
            extend = extend,
            only_tss = only_tss,
            verbose = verbose
        )
    } else {
        peaks_near_gene <- find_peaks_near_genes(
            peaks = regions@ranges,
            genes = peak_to_gene_domains,
            method = 'Signac',
            upstream = 0,
            downstream = 0,
            only_tss = FALSE,
            verbose = verbose
        )
    }
    peaks2gene <- aggregate_matrix(
        Matrix::t(peaks_near_gene),
        groups = colnames(peaks_near_gene),
        fun = 'sum'
    )

    common_peaks <- Reduce(intersect, list(
        colnames(peaks2gene), rownames(peaks2motif), colnames(peak_data)
    ))
    peaks2gene <- peaks2gene[, common_peaks, drop = FALSE]
    peaks2motif <- peaks2motif[common_peaks, , drop = FALSE]
    peak_data <- peak_data[, common_peaks, drop = FALSE]
    peaks_at_gene <- as.logical(sparseMatrixStats::colMaxs(peaks2gene) > 0)
    peaks_with_motif <- as.logical(sparseMatrixStats::rowMaxs(peaks2motif * 1) > 0)
    peaks_use <- peaks_at_gene & peaks_with_motif
    peaks2gene <- peaks2gene[, peaks_use, drop = FALSE]
    peaks2motif <- peaks2motif[peaks_use, , drop = FALSE]
    peak_data <- peak_data[, peaks_use, drop = FALSE]

    motif_ids <- intersect(colnames(peaks2motif), rownames(motif2tf))
    if (length(motif_ids) == 0L) {
        stop('Motif match columns do not overlap the motif-to-TF mapping.')
    }
    peaks2motif <- peaks2motif[, motif_ids, drop = FALSE]
    motif2tf <- motif2tf[motif_ids, , drop = FALSE]
    available_tfs <- intersect(colnames(motif2tf), colnames(gene_data))
    if (length(available_tfs) == 0L) {
        stop('No motif-mapped transcription factors were found in the RNA assay.')
    }
    motif2tf <- motif2tf[, available_tfs, drop = FALSE]

    list(
        gene_data = gene_data,
        peak_data = peak_data,
        metadata = metadata,
        features = features,
        peaks2gene = peaks2gene,
        peaks2motif = peaks2motif,
        motif2tf = motif2tf
    )
}

.condition_get_model_data <- function(
    object,
    params
) {
    metadata <- object@data@meta.data
    gene_data <- Matrix::t(LayerData(
        object, assay = params$rna_assay, layer = 'data'
    ))
    peak_data <- Matrix::t(LayerData(
        object, assay = params$peak_assay, layer = 'data'
    ))
    gene_cells <- rownames(gene_data)
    peak_cells <- rownames(peak_data)
    if (is.null(gene_cells) || is.null(peak_cells) ||
        anyDuplicated(gene_cells) || anyDuplicated(peak_cells) ||
        !setequal(gene_cells, peak_cells)) {
        stop('RNA and ATAC data must contain the same paired single cells.')
    }
    if (is.null(rownames(metadata)) || anyDuplicated(rownames(metadata)) ||
        !all(gene_cells %in% rownames(metadata))) {
        stop('Metadata must contain one unique row for every paired single cell.')
    }
    common_cells <- gene_cells
    list(
        gene_data = gene_data[common_cells, , drop = FALSE],
        peak_data = peak_data[common_cells, , drop = FALSE],
        metadata = metadata[common_cells, , drop = FALSE]
    )
}
