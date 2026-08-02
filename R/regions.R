#' @importFrom Signac StringToGRanges
#' @importFrom Seurat GetAssay VariableFeatures
#' @importFrom S4Vectors subjectHits queryHits
#' @importFrom IRanges findOverlaps pintersect
NULL


#' Initiate the \code{RegulatoryNetwork} object.
#'
#' @param regions Candidate regions to consider for binding site inference.
#' If \code{NULL}, all peaks regions are considered.
#' @param peak_assay A character vector indicating the name of the chromatin
#' accessibility assay in the \code{Seurat} object.
#' @param rna_assay A character vector indicating the name of the gene expression
#' assay in the \code{Seurat} object.
#' @param exclude_exons Logical. Whether to consider exons for binding site inference.
#'
#' @return A GRNData object containing a RegulatoryNetwork object.
#'
#' @rdname initiate_grn
#' @export
#' @method initiate_grn Seurat
initiate_grn.Seurat <- function(
    object,
    regions = NULL,
    peak_assay = 'peaks',
    rna_assay = 'RNA',
    exclude_exons = TRUE
){
    gene_annot <- Signac::Annotation(object[[peak_assay]])
    # Through error if NULL
    if (is.null(gene_annot)){
        stop('Please provide a gene annotation for the ChromatinAssay.')
    }
    peak_ranges <- StringToGRanges(
        rownames(Seurat::GetAssay(object, assay=peak_assay))
    )

    # Find candidate ranges by intersecting the supplied regions with peaks
    # Per default take all peaks as candidate regions
    if (!is.null(regions)){
        cand_olaps <- findOverlaps(regions, peak_ranges)
        cand_ranges <- pintersect(
            peak_ranges[subjectHits(cand_olaps)],
            regions[queryHits(cand_olaps)]
        )
    } else {
        cand_ranges <- peak_ranges
    }

    # Exclude exons because they are usually conserved
    if (exclude_exons){
        exon_ranges <- gene_annot[gene_annot$type=='exon', ]
        names(exon_ranges@ranges) <- NULL
        exon_ranges <- IRanges::intersect(exon_ranges, exon_ranges)
        exon_ranges <- GenomicRanges::GRanges(
            seqnames = exon_ranges@seqnames,
            ranges = exon_ranges@ranges
        )
        cand_ranges <- GenomicRanges::subtract(
            cand_ranges, exon_ranges, ignore.strand=TRUE
        ) %>% unlist()
    }

    # Match candidate ranges to peaks
    peak_overlaps <- findOverlaps(cand_ranges, peak_ranges)
    peak_matches <- subjectHits(peak_overlaps)

    regions_obj <- new(
        Class = 'Regions',
        ranges = cand_ranges,
        peaks = peak_matches,
        motifs = NULL
    )

    params <- list(
        peak_assay = peak_assay,
        rna_assay = rna_assay,
        exclude_exons = exclude_exons
    )

    grn_obj <- new(
        Class = 'RegulatoryNetwork',
        regions = regions_obj,
        params = params
    )

    object <- new(
        Class = 'GRNData',
        grn = grn_obj,
        data = object
    )

    return(object)
}

#' Initiate the \code{RegulatoryNetwork} object.
#'
#' @param object An object.
#'
#' @rdname initiate_grn
#' @export
#' @method initiate_grn GRNData
initiate_grn.GRNData <- function(
        object,
        regions = NULL,
        peak_assay = 'peaks',
        rna_assay = 'RNA',
        exclude_exons = TRUE
){
    return(initiate_grn(
        object = object@data,
        regions = regions,
        peak_assay = peak_assay,
        rna_assay = rna_assay,
        exclude_exons = exclude_exons
    ))
}

.load_default_motif2tf <- function(){
    data_env <- new.env(parent = emptyenv())
    utils::data(
        "motif2tf",
        package = "Pando",
        envir = data_env
    )
    if (!exists("motif2tf", envir = data_env, inherits = FALSE)){
        stop(
            "The bundled Pando `motif2tf` dataset could not be loaded.",
            call. = FALSE
        )
    }
    motif2tf <- get("motif2tf", envir = data_env, inherits = FALSE)
    if (!is.data.frame(motif2tf) || ncol(motif2tf) < 2L){
        stop(
            "The bundled Pando `motif2tf` dataset must be a data frame with at least two columns.",
            call. = FALSE
        )
    }
    motif2tf
}

.pando_motif_cache_digest <- function(value){
    path <- tempfile("pando-motif-signature-", fileext = ".rds")
    on.exit(unlink(path, force = TRUE), add = TRUE)
    saveRDS(value, path, version = 2L, compress = FALSE)
    unname(tools::md5sum(path)[[1L]])
}

.pando_motif_genome_signature <- function(genome){
    slots <- if (base::isS4(genome)) methods::slotNames(genome) else character()
    fields <- intersect(
        slots,
        c("pkgname", "organism", "common_name", "provider", "provider_version", "release_name")
    )
    values <- lapply(fields, function(field){
        value <- methods::slot(genome, field)
        if (length(value) > 20L) value <- value[seq_len(20L)]
        value
    })
    names(values) <- fields
    list(
        class = class(genome),
        class_package = attr(class(genome), "package"),
        fields = values
    )
}

.pando_motif_cache_key <- function(cand_ranges, pfm, genome, exact_positions){
    .pando_motif_cache_digest(list(
        schema_version = "pando_motif_cache_v1",
        candidate_ranges = Signac::GRangesToString(cand_ranges),
        pfm = pfm,
        genome = .pando_motif_genome_signature(genome),
        exact_positions = isTRUE(exact_positions)
    ))
}

.pando_read_motif_cache <- function(path, candidate_ranges, exact_positions){
    if (!file.exists(path)) return(NULL)
    payload <- tryCatch(readRDS(path), error = function(error) NULL)
    if (!is.list(payload) ||
        !identical(payload$schema_version, "pando_motif_cache_v1") ||
        !identical(payload$candidate_ranges, candidate_ranges) ||
        !identical(payload$exact_positions, isTRUE(exact_positions)) ||
        is.null(payload$motif_object)){
        return(NULL)
    }
    payload$motif_object
}

.pando_write_motif_cache <- function(path, motif_object, candidate_ranges, exact_positions){
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    temporary <- tempfile(
        paste0(basename(path), "."),
        tmpdir = dirname(path),
        fileext = ".tmp"
    )
    on.exit(unlink(temporary, force = TRUE), add = TRUE)
    saveRDS(
        list(
            schema_version = "pando_motif_cache_v1",
            candidate_ranges = candidate_ranges,
            exact_positions = isTRUE(exact_positions),
            motif_object = motif_object
        ),
        temporary,
        version = 2L
    )
    if (!file.exists(path)){
        if (!file.rename(temporary, path) && !file.exists(path)){
            stop("Could not publish the Pando motif cache atomically.", call. = FALSE)
        }
    }
    invisible(path)
}

#' Scan for motifs in candidate regions.
#'
#' @importFrom dplyr distinct mutate select
#'
#' @param pfm A \code{PFMatrixList} object with position weight matrices.
#' @param genome A \code{BSgenome} object with the genome of interest.
#' @param motif_tfs A data frame matching motifs with TFs. The first column is assumed
#' to be the name of the motif, the second the name of the TF.
#' @param exact_positions Whether to retain the exact genomic positions of motif
#' matches. The default, \code{FALSE}, stores only the sparse motif
#' presence/absence matrix used to construct TF-peak-target candidates. Set to
#' \code{TRUE} when exact match positions are needed, for example for TF
#' footprinting.
#' @param cache_dir Optional directory for persistent motif-match caching. Cache
#' entries are keyed by the ordered candidate ranges, complete PFM collection,
#' genome identity and \code{exact_positions}. Identical region layouts therefore
#' reuse one motif scan across cell types and later runs.
#' @param reuse_cache Logical. Read a valid existing cache entry before scanning.
#' A newly computed entry is still written when this is \code{FALSE}.
#' @param verbose Display messages.
#'
#' @return A GRNData object with updated motif info.
#'
#' @rdname find_motifs
#' @export
#' @method find_motifs GRNData
find_motifs.GRNData <- function(
    object,
    pfm,
    genome,
    motif_tfs = NULL,
    verbose = TRUE,
    exact_positions = FALSE,
    cache_dir = getOption("Pando.motif_cache_dir", NULL),
    reuse_cache = TRUE
){
    if (!is.logical(exact_positions) || length(exact_positions) != 1L ||
        is.na(exact_positions)){
        stop("`exact_positions` must be either TRUE or FALSE.", call. = FALSE)
    }
    if (!is.logical(reuse_cache) || length(reuse_cache) != 1L ||
        is.na(reuse_cache)){
        stop("`reuse_cache` must be either TRUE or FALSE.", call. = FALSE)
    }
    if (!is.null(cache_dir) &&
        (!is.character(cache_dir) || length(cache_dir) != 1L ||
         is.na(cache_dir) || !nzchar(trimws(cache_dir)))){
        stop("`cache_dir` must be NULL or one non-empty path.", call. = FALSE)
    }
    params <- Params(object)

    # Add TF info for motifs
    log_message('Adding TF info', verbose=verbose)
    if (!is.null(motif_tfs)){
        motif2tf <- motif_tfs
    } else {
        motif2tf <- .load_default_motif2tf()
    }
    if (!is.data.frame(motif2tf) || ncol(motif2tf) < 2L){
        stop(
            "`motif_tfs` must be a data frame with motif names in the first column and TF names in the second column.",
            call. = FALSE
        )
    }

    # Spread dataframe to sparse matrix
    motif2tf <- motif2tf %>% select('motif'=1,'tf'=2) %>%
        distinct() %>% mutate(val=1) %>%
        tidyr::pivot_wider(names_from = 'tf', values_from=val, values_fill=0) %>%
        tibble::column_to_rownames('motif') %>%
        as.matrix() %>% Matrix::Matrix(sparse=TRUE)
    tfs_use <- intersect(
        rownames(GetAssay(object, params$rna_assay)),
        colnames(motif2tf)
    )
    if (length(tfs_use)==0){
        stop('None of the provided TFs were found in the dataset. Consider providing a custom motif-to-TF map as `motif_tfs`')
    }
    object@grn@regions@tfs <- motif2tf[, tfs_use, drop = FALSE]

    # Find motif matches with Signac/motifmatchr. Exact positions are optional;
    # Pando's candidate construction only requires the presence/absence matrix.
    cand_ranges <- object@grn@regions@ranges
    candidate_ranges <- Signac::GRangesToString(cand_ranges)
    cache_key <- NULL
    cache_file <- NULL
    motif_object <- NULL
    cache_hit <- FALSE
    if (!is.null(cache_dir)){
        cache_key <- .pando_motif_cache_key(
            cand_ranges = cand_ranges,
            pfm = pfm,
            genome = genome,
            exact_positions = exact_positions
        )
        cache_file <- file.path(
            cache_dir,
            paste0("motif_matches__", cache_key, ".rds")
        )
        if (isTRUE(reuse_cache)){
            motif_object <- .pando_read_motif_cache(
                cache_file,
                candidate_ranges = candidate_ranges,
                exact_positions = exact_positions
            )
            cache_hit <- !is.null(motif_object)
        }
    }

    if (is.null(motif_object)){
        if (exact_positions){
            motif_object <- Signac::AddMotifs(
                object = cand_ranges,
                genome = genome,
                pfm = pfm,
                verbose = verbose
            )
        } else {
            motif_matrix <- Signac::CreateMotifMatrix(
                features = cand_ranges,
                pwm = pfm,
                genome = genome,
                score = FALSE,
                use.counts = FALSE
            )
            motif_object <- Signac::CreateMotifObject(
                data = motif_matrix,
                pwm = pfm,
                positions = NULL
            )
        }
        if (!is.null(cache_file)){
            .pando_write_motif_cache(
                cache_file,
                motif_object = motif_object,
                candidate_ranges = candidate_ranges,
                exact_positions = exact_positions
            )
        }
    } else {
        log_message('Reusing cached motif matches', verbose=verbose)
    }
    object@grn@regions@motifs <- motif_object
    object@grn@params$motif_cache <- list(
        enabled = !is.null(cache_dir),
        hit = cache_hit,
        cache_key = cache_key,
        cache_file = cache_file,
        schema_version = "pando_motif_cache_v1"
    )

    return(object)
}
