# Shared pre-fit TF-peak-gene design construction.

.pando_design_hash <- function(x) {
    text <- paste(as.character(x), collapse = "\n")
    values <- utf8ToInt(enc2utf8(text))
    hash <- 0
    modulus <- 2147483647
    if (length(values)) {
        for (value in values) {
            hash <- (hash * 131 + value) %% modulus
        }
    }
    sprintf("%08x", as.integer(hash))
}

.pando_detection_fraction <- function(x) {
    x <- as.numeric(x)
    if (!length(x)) return(0)
    sum(is.finite(x) & x > 0) / length(x)
}

.pando_group_max_abs_cor <- function(x, y, groups) {
    x <- as.matrix(x)
    y <- as.numeric(y)
    groups <- as.character(groups)
    answer <- rep(0, ncol(x))
    if (!ncol(x)) return(answer)
    for (group in unique(groups)) {
        index <- which(groups == group & is.finite(y))
        if (length(index) < 3L) next
        y_group <- y[index]
        if (!is.finite(stats::sd(y_group)) || stats::sd(y_group) <= 0) next
        x_group <- x[index, , drop = FALSE]
        value <- suppressWarnings(stats::cor(x_group, y_group, use = "pairwise.complete.obs"))
        value <- as.numeric(value)
        value[!is.finite(value)] <- 0
        answer <- pmax(answer, abs(value))
    }
    answer
}

.pando_candidate_edge_table <- function(
    peaks2gene, peaks2motif, motif2tf, region_to_peak,
    peak_to_gene_source = "Signac") {
    peaks2gene <- as.matrix(peaks2gene)
    peaks2motif <- as.matrix(peaks2motif)
    motif2tf <- as.matrix(motif2tf)
    if (is.null(rownames(peaks2gene)) || is.null(colnames(peaks2gene)) ||
        is.null(rownames(peaks2motif)) || is.null(colnames(peaks2motif)) ||
        is.null(rownames(motif2tf)) || is.null(colnames(motif2tf))) {
        stop("Candidate design matrices require complete dimnames.", call. = FALSE)
    }
    regions <- Reduce(intersect, list(
        colnames(peaks2gene), rownames(peaks2motif), names(region_to_peak)
    ))
    motifs <- intersect(colnames(peaks2motif), rownames(motif2tf))
    if (!length(regions) || !length(motifs)) return(data.frame())
    peaks2gene <- peaks2gene[, regions, drop = FALSE]
    peaks2motif <- peaks2motif[regions, motifs, drop = FALSE]
    motif2tf <- motif2tf[motifs, , drop = FALSE]

    output <- list()
    output_index <- 0L
    for (target in rownames(peaks2gene)) {
        target_regions <- regions[as.numeric(peaks2gene[target, ]) > 0]
        for (region in target_regions) {
            region_motifs <- motifs[as.numeric(peaks2motif[region, ]) > 0]
            if (!length(region_motifs)) next
            tf_present <- colSums(motif2tf[region_motifs, , drop = FALSE] != 0) > 0
            tfs <- setdiff(colnames(motif2tf)[tf_present], target)
            for (tf in tfs) {
                matched <- region_motifs[
                    as.numeric(motif2tf[region_motifs, tf, drop = TRUE]) != 0
                ]
                output_index <- output_index + 1L
                output[[output_index]] <- data.frame(
                    edge_id = paste(tf, region, target, sep = "::"),
                    tf = tf,
                    region = region,
                    target = target,
                    tf_feature_id = tf,
                    atac_feature_id = unname(region_to_peak[[region]]),
                    target_feature_id = target,
                    motif_id = paste(sort(unique(matched)), collapse = ";"),
                    peak_to_gene_source = peak_to_gene_source,
                    stringsAsFactors = FALSE
                )
            }
        }
    }
    if (!length(output)) return(data.frame())
    answer <- unique(do.call(rbind, output))
    answer <- answer[order(answer$target, answer$region, answer$tf), , drop = FALSE]
    rownames(answer) <- NULL
    answer
}

#' Prepare a shared pre-fit TF-peak-gene design
#'
#' Builds the structural candidate universe used before coefficient estimation.
#' The returned edge table is suitable for multi-task models in which every
#' condition must use the same TF-peak-target columns. By default no target-RNA
#' correlation filter is applied, so condition-specific sign reversals cannot
#' be removed by pooled cancellation before fitting.
#'
#' @param object A `GRNData` object after [find_motifs()].
#' @param genes Target genes. Defaults to RNA variable features.
#' @param peak_to_gene_method Peak-to-gene method passed to
#'   [find_peaks_near_genes()].
#' @param upstream,downstream,extend,only_tss Regulatory-domain parameters.
#' @param peak_to_gene_domains Optional supplied regulatory domains.
#' @param screen_method `"structural"` or
#'   `"union_within_group_correlation"`. The latter retains an edge when its
#'   peak-target and TF-target correlations pass thresholds in at least one
#'   group, while still returning one shared edge universe.
#' @param screen_group_col Metadata column defining screening groups.
#' @param tf_cor,peak_cor Absolute within-group screening thresholds.
#' @param min_tf_detection,min_peak_detection,min_target_detection Minimum
#'   fractions of cells with positive normalized signal.
#' @param max_edges_per_target Optional deterministic cap per target.
#' @param verbose Display progress messages.
#'
#' @return A `PandoGRNDesign` list containing `candidate_edges`, `region_map`,
#'   `target_diagnostics`, a matrix/feature contract, parameters and a stable
#'   design identifier.
#' @rdname prepare_grn_design
#' @export
#' @method prepare_grn_design GRNData
prepare_grn_design.GRNData <- function(
    object,
    genes = NULL,
    peak_to_gene_method = c("Signac", "GREAT"),
    upstream = 100000,
    downstream = 0,
    extend = 1000000,
    only_tss = FALSE,
    peak_to_gene_domains = NULL,
    screen_method = c("structural", "union_within_group_correlation"),
    screen_group_col = NULL,
    tf_cor = 0,
    peak_cor = 0,
    min_tf_detection = 0.01,
    min_peak_detection = 0.01,
    min_target_detection = 0.01,
    max_edges_per_target = Inf,
    verbose = TRUE) {
    peak_to_gene_method <- match.arg(peak_to_gene_method)
    screen_method <- match.arg(screen_method)
    numeric_probability <- function(value, name) {
        if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
            value < 0 || value > 1) {
            stop("`", name, "` must be one finite number in [0, 1].",
                 call. = FALSE)
        }
    }
    numeric_probability(min_tf_detection, "min_tf_detection")
    numeric_probability(min_peak_detection, "min_peak_detection")
    numeric_probability(min_target_detection, "min_target_detection")
    if (!is.numeric(tf_cor) || length(tf_cor) != 1L || !is.finite(tf_cor) ||
        tf_cor < 0 || tf_cor > 1 || !is.numeric(peak_cor) ||
        length(peak_cor) != 1L || !is.finite(peak_cor) || peak_cor < 0 ||
        peak_cor > 1) {
        stop("`tf_cor` and `peak_cor` must be finite numbers in [0, 1].",
             call. = FALSE)
    }
    if (!is.numeric(max_edges_per_target) || length(max_edges_per_target) != 1L ||
        is.na(max_edges_per_target) || max_edges_per_target <= 0) {
        stop("`max_edges_per_target` must be positive or Inf.", call. = FALSE)
    }

    params <- Params(object)
    motif2tf <- NetworkTFs(object)
    if (is.null(motif2tf)) {
        stop("Motif matches have not been found. Run `find_motifs()` first.",
             call. = FALSE)
    }
    gene_annot <- Signac::Annotation(GetAssay(object, params$peak_assay))
    if (is.null(gene_annot)) {
        stop("Please provide a gene annotation for the ChromatinAssay.",
             call. = FALSE)
    }
    if (is.null(genes)) {
        genes <- VariableFeatures(object, assay = params$rna_assay)
    }
    genes <- unique(as.character(genes))
    if (!length(genes)) stop("No target genes were supplied.", call. = FALSE)

    gene_data <- Matrix::t(LayerData(
        object, assay = params$rna_assay, layer = "data"
    ))
    peak_source <- Matrix::t(LayerData(
        object, assay = params$peak_assay, layer = "data"
    ))
    cells <- rownames(gene_data)
    if (!identical(cells, rownames(peak_source))) {
        stop("RNA and peak matrices must contain cells in identical order.",
             call. = FALSE)
    }
    features <- Reduce(intersect, list(
        as.character(gene_annot$gene_name), genes, colnames(gene_data)
    ))
    gene_annot <- gene_annot[gene_annot$gene_name %in% features, ]
    if (!length(features)) {
        stop("No requested genes overlap RNA features and the annotation.",
             call. = FALSE)
    }

    regions <- NetworkRegions(object)
    if (!length(regions@peaks)) {
        stop("Pando contains no regulatory regions.", call. = FALSE)
    }
    region_ids <- rownames(regions@motifs@data)
    if (length(region_ids) != length(regions@peaks)) {
        stop("Pando region and peak mappings have different lengths.",
             call. = FALSE)
    }
    source_peak_ids <- colnames(peak_source)[regions@peaks]
    if (anyNA(source_peak_ids) || any(!nzchar(source_peak_ids))) {
        stop("Pando regions could not be mapped to source ATAC features.",
             call. = FALSE)
    }
    region_to_peak <- stats::setNames(source_peak_ids, region_ids)
    peak_data <- peak_source[, regions@peaks, drop = FALSE]
    colnames(peak_data) <- region_ids
    peaks2motif <- regions@motifs@data

    log_message("Selecting shared candidate regulatory regions", verbose = verbose)
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
            method = "Signac",
            genes = peak_to_gene_domains,
            upstream = 0,
            downstream = 0,
            only_tss = FALSE
        )
    }
    peaks2gene <- aggregate_matrix(
        t(peaks_near_gene), groups = colnames(peaks_near_gene), fun = "sum"
    )
    peaks_at_gene <- as.logical(sparseMatrixStats::colMaxs(peaks2gene))
    peaks_with_motif <- as.logical(
        sparseMatrixStats::rowMaxs(peaks2motif * 1)
    )
    peaks_use <- peaks_at_gene & peaks_with_motif
    peaks2gene <- peaks2gene[, peaks_use, drop = FALSE]
    peaks2motif <- peaks2motif[peaks_use, , drop = FALSE]
    peak_data <- peak_data[, peaks_use, drop = FALSE]
    region_to_peak <- region_to_peak[colnames(peak_data)]

    expressed_tfs <- intersect(colnames(motif2tf), colnames(gene_data))
    motif2tf <- motif2tf[, expressed_tfs, drop = FALSE]
    candidates <- .pando_candidate_edge_table(
        peaks2gene = peaks2gene,
        peaks2motif = peaks2motif,
        motif2tf = motif2tf,
        region_to_peak = region_to_peak,
        peak_to_gene_source = if (is.null(peak_to_gene_domains)) {
            peak_to_gene_method
        } else {
            "provided_domains"
        }
    )
    if (!nrow(candidates)) {
        stop("No structural TF-peak-target candidates were found.",
             call. = FALSE)
    }

    tf_detection <- vapply(
        unique(candidates$tf),
        function(tf) .pando_detection_fraction(gene_data[, tf]),
        numeric(1)
    )
    target_detection <- vapply(
        unique(candidates$target),
        function(target) .pando_detection_fraction(gene_data[, target]),
        numeric(1)
    )
    peak_detection <- vapply(
        unique(candidates$region),
        function(region) .pando_detection_fraction(peak_data[, region]),
        numeric(1)
    )
    candidates$tf_detection <- unname(tf_detection[candidates$tf])
    candidates$target_detection <- unname(
        target_detection[candidates$target]
    )
    candidates$peak_detection <- unname(peak_detection[candidates$region])
    candidates <- candidates[
        candidates$tf_detection >= min_tf_detection &
        candidates$target_detection >= min_target_detection &
        candidates$peak_detection >= min_peak_detection,
        , drop = FALSE
    ]
    if (!nrow(candidates)) {
        stop("No candidate edges passed the detection filters.", call. = FALSE)
    }

    candidates$peak_target_screen_score <- NA_real_
    candidates$tf_target_screen_score <- NA_real_
    if (identical(screen_method, "union_within_group_correlation")) {
        if (is.null(screen_group_col) || length(screen_group_col) != 1L ||
            !screen_group_col %in% colnames(object@data@meta.data)) {
            stop("A valid `screen_group_col` is required for grouped screening.",
                 call. = FALSE)
        }
        groups <- object@data@meta.data[cells, screen_group_col, drop = TRUE]
        if (anyNA(groups) || any(!nzchar(trimws(as.character(groups))))) {
            stop("Grouped screening metadata must be complete.", call. = FALSE)
        }
        for (target in unique(candidates$target)) {
            index <- which(candidates$target == target)
            regions_target <- unique(candidates$region[index])
            tfs_target <- unique(candidates$tf[index])
            peak_score <- .pando_group_max_abs_cor(
                peak_data[, regions_target, drop = FALSE],
                gene_data[, target], groups
            )
            tf_score <- .pando_group_max_abs_cor(
                gene_data[, tfs_target, drop = FALSE],
                gene_data[, target], groups
            )
            names(peak_score) <- regions_target
            names(tf_score) <- tfs_target
            candidates$peak_target_screen_score[index] <-
                unname(peak_score[candidates$region[index]])
            candidates$tf_target_screen_score[index] <-
                unname(tf_score[candidates$tf[index]])
        }
        candidates <- candidates[
            candidates$peak_target_screen_score > peak_cor &
            candidates$tf_target_screen_score > tf_cor,
            , drop = FALSE
        ]
        if (!nrow(candidates)) {
            stop("No candidate edges passed grouped correlation screening.",
                 call. = FALSE)
        }
    }

    if (is.finite(max_edges_per_target)) {
        rows <- split(seq_len(nrow(candidates)), candidates$target)
        keep <- unlist(lapply(rows, function(index) {
            score <- pmax(
                candidates$peak_target_screen_score[index],
                candidates$tf_target_screen_score[index],
                na.rm = TRUE
            )
            score[!is.finite(score)] <-
                candidates$peak_detection[index] *
                candidates$tf_detection[index]
            order(-score, candidates$edge_id[index])[
                seq_len(min(length(index), as.integer(max_edges_per_target)))
            ] |> function(local) index[local]
        }), use.names = FALSE)
        candidates <- candidates[sort(keep), , drop = FALSE]
    }
    candidates <- candidates[
        order(candidates$target, candidates$region, candidates$tf),
        , drop = FALSE
    ]
    rownames(candidates) <- NULL

    target_rows <- split(seq_len(nrow(candidates)), candidates$target)
    target_diagnostics <- do.call(rbind, lapply(target_rows, function(index) {
        one <- candidates[index, , drop = FALSE]
        data.frame(
            target = one$target[[1L]],
            n_candidate_edges = nrow(one),
            n_candidate_regions = length(unique(one$region)),
            n_candidate_tfs = length(unique(one$tf)),
            target_detection = one$target_detection[[1L]],
            stringsAsFactors = FALSE
        )
    }))
    rownames(target_diagnostics) <- NULL
    region_map <- unique(candidates[, c("region", "atac_feature_id"), drop = FALSE])
    design_id <- paste0(
        "pando_grn_design_v1_",
        .pando_design_hash(candidates$edge_id)
    )
    answer <- list(
        schema_version = "pando_grn_design_v1",
        candidate_edges = candidates,
        region_map = region_map,
        target_diagnostics = target_diagnostics,
        feature_contract = list(
            cell_ids = cells,
            rna_feature_ids = colnames(gene_data),
            atac_feature_ids = colnames(peak_source),
            rna_assay = params$rna_assay,
            atac_assay = params$peak_assay,
            rna_layer = "data",
            atac_layer = "data"
        ),
        params = list(
            peak_to_gene_method = peak_to_gene_method,
            upstream = upstream,
            downstream = downstream,
            extend = extend,
            only_tss = only_tss,
            screen_method = screen_method,
            screen_group_col = screen_group_col,
            tf_cor = tf_cor,
            peak_cor = peak_cor,
            min_tf_detection = min_tf_detection,
            min_peak_detection = min_peak_detection,
            min_target_detection = min_target_detection,
            max_edges_per_target = max_edges_per_target
        ),
        design_id = design_id
    )
    class(answer) <- c("PandoGRNDesign", "list")
    answer
}

#' Validate a Pando shared GRN design
#'
#' @param design A `PandoGRNDesign` object.
#' @return The design invisibly.
#' @export
validate_grn_design <- function(design) {
    if (!inherits(design, "PandoGRNDesign") || !is.list(design)) {
        stop("`design` must inherit from PandoGRNDesign.", call. = FALSE)
    }
    required <- c(
        "schema_version", "candidate_edges", "region_map",
        "target_diagnostics", "feature_contract", "params", "design_id"
    )
    missing <- setdiff(required, names(design))
    if (length(missing)) {
        stop("GRN design is missing fields: ", paste(missing, collapse = ", "),
             call. = FALSE)
    }
    edges <- design$candidate_edges
    edge_columns <- c(
        "edge_id", "tf", "region", "target", "tf_feature_id",
        "atac_feature_id", "target_feature_id"
    )
    if (!is.data.frame(edges) || !nrow(edges) ||
        !all(edge_columns %in% colnames(edges))) {
        stop("GRN design candidate edge table is incomplete.", call. = FALSE)
    }
    if (anyNA(edges$edge_id) || any(!nzchar(edges$edge_id)) ||
        anyDuplicated(edges$edge_id)) {
        stop("GRN design edge IDs must be unique and non-empty.", call. = FALSE)
    }
    contract <- design$feature_contract
    if (!is.list(contract) ||
        !all(c("cell_ids", "rna_feature_ids", "atac_feature_ids") %in%
             names(contract))) {
        stop("GRN design feature contract is incomplete.", call. = FALSE)
    }
    if (any(!edges$tf_feature_id %in% contract$rna_feature_ids) ||
        any(!edges$target_feature_id %in% contract$rna_feature_ids) ||
        any(!edges$atac_feature_id %in% contract$atac_feature_ids)) {
        stop("GRN design edge features violate the matrix contract.",
             call. = FALSE)
    }
    invisible(design)
}
