# Hardened shared-design contract loaded after grn_design.R.

.pando_grn_design_base <- prepare_grn_design.GRNData

.pando_design_md5 <- function(x) {
    path <- tempfile(fileext = '.rds')
    on.exit(unlink(path), add = TRUE)
    saveRDS(x, path, version = 2, compress = FALSE)
    unname(as.character(tools::md5sum(path)))
}

.pando_split_supporting_regions <- function(x) {
    x <- as.character(x)
    out <- strsplit(x, ';', fixed = TRUE)
    lapply(out, function(value) {
        value <- trimws(value)
        sort(unique(value[!is.na(value) & nzchar(value)]))
    })
}

.pando_design_payload <- function(design) {
    edges <- design$candidate_edges
    if (is.data.frame(edges) && nrow(edges)) {
        order_columns <- intersect(c(
            'candidate_index', 'edge_id', 'tf', 'region', 'target',
            'atac_feature_id', 'tf_feature_id', 'target_feature_id',
            'supporting_regions', 'n_supporting_regions',
            'motif_supported', 'peak_to_gene_supported',
            'tf_detection', 'peak_detection', 'target_detection'
        ), colnames(edges))
        edges <- edges[, order_columns, drop = FALSE]
        rownames(edges) <- NULL
    }
    region_map <- design$region_map
    if (is.data.frame(region_map) && nrow(region_map)) {
        region_map <- region_map[order(
            as.character(region_map$region),
            as.character(region_map$atac_feature_id)
        ), , drop = FALSE]
        rownames(region_map) <- NULL
    }
    diagnostics <- design$target_diagnostics
    if (is.data.frame(diagnostics) && nrow(diagnostics)) {
        diagnostics <- diagnostics[order(as.character(diagnostics$target)), , drop = FALSE]
        rownames(diagnostics) <- NULL
    }
    contract <- design$feature_contract
    if (is.list(contract)) {
        contract$cell_ids <- as.character(contract$cell_ids)
        contract$rna_feature_ids <- as.character(contract$rna_feature_ids)
        contract$atac_feature_ids <- as.character(contract$atac_feature_ids)
    }
    list(
        schema_version = as.character(design$schema_version),
        candidate_edges = edges,
        region_map = region_map,
        target_diagnostics = diagnostics,
        feature_contract = contract,
        params = design$params
    )
}

.pando_refresh_grn_design_contract <- function(design) {
    edges <- design$candidate_edges
    targets <- as.character(design$target_diagnostics$target)
    if (is.data.frame(edges) && nrow(edges)) {
        support <- .pando_split_supporting_regions(edges$supporting_regions)
        all_regions <- function(index) {
            values <- unique(unlist(support[index], use.names = FALSE))
            length(values[!is.na(values) & nzchar(values)])
        }
        diagnostics <- data.frame(
            target = targets,
            target_detection = as.numeric(
                design$target_diagnostics$target_detection[
                    match(targets, design$target_diagnostics$target)
                ]
            ),
            n_candidate_edges = vapply(targets, function(g) {
                sum(as.character(edges$target) == g)
            }, integer(1)),
            n_candidate_tfs = vapply(targets, function(g) {
                index <- which(as.character(edges$target) == g)
                length(unique(as.character(edges$tf[index])))
            }, integer(1)),
            n_candidate_atac_features = vapply(targets, function(g) {
                index <- which(as.character(edges$target) == g)
                length(unique(as.character(edges$atac_feature_id[index])))
            }, integer(1)),
            n_candidate_regions = vapply(targets, function(g) {
                all_regions(which(as.character(edges$target) == g))
            }, integer(1)),
            stringsAsFactors = FALSE
        )
        design$target_diagnostics <- diagnostics
    } else {
        diagnostics <- design$target_diagnostics
        diagnostics$n_candidate_atac_features <- integer(nrow(diagnostics))
        if (!'n_candidate_regions' %in% colnames(diagnostics)) {
            diagnostics$n_candidate_regions <- integer(nrow(diagnostics))
        }
        design$target_diagnostics <- diagnostics
    }
    design$schema_version <- 'pando_grn_design_v2'
    design$design_fingerprint <- paste0(
        'md5:', .pando_design_md5(.pando_design_payload(design))
    )
    class(design) <- c('PandoGRNDesign', 'list')
    design
}

#' @rdname prepare_grn_design
#' @method prepare_grn_design GRNData
#' @export
prepare_grn_design.GRNData <- function(object, ...) {
    design <- .pando_grn_design_base(object = object, ...)
    .pando_refresh_grn_design_contract(design)
}

#' Validate a Pando shared candidate design
#'
#' Validation checks the structural edge identity, feature contract,
#' region-to-measured-peak mapping, candidate ordering, and the deterministic
#' design fingerprint. Version-1 objects remain readable but do not carry the
#' strengthened version-2 fingerprint contract.
#'
#' @param design A `PandoGRNDesign` object.
#' @return `TRUE` invisibly; otherwise an error is raised.
#' @export
validate_grn_design <- function(design) {
    if (!inherits(design, 'PandoGRNDesign') || !is.list(design)) {
        stop('`design` must be a PandoGRNDesign object.', call. = FALSE)
    }
    schema <- as.character(design$schema_version)
    if (!schema %in% c('pando_grn_design_v1', 'pando_grn_design_v2')) {
        stop('Unsupported Pando GRN design schema.', call. = FALSE)
    }
    required_objects <- c(
        'candidate_edges', 'region_map', 'target_diagnostics',
        'feature_contract', 'params', 'design_fingerprint'
    )
    if (!all(required_objects %in% names(design))) {
        stop('Pando GRN design object is incomplete.', call. = FALSE)
    }
    edges <- design$candidate_edges
    required_edges <- c(
        'edge_id', 'candidate_index', 'tf', 'region', 'target',
        'atac_feature_id', 'tf_feature_id', 'target_feature_id',
        'supporting_regions', 'n_supporting_regions'
    )
    if (!is.data.frame(edges) || !all(required_edges %in% colnames(edges))) {
        stop('Pando GRN design candidate edges are incomplete.', call. = FALSE)
    }
    character_columns <- setdiff(required_edges, c(
        'candidate_index', 'n_supporting_regions'
    ))
    if (nrow(edges) && any(vapply(
        edges[, character_columns, drop = FALSE],
        function(x) anyNA(x) || any(!nzchar(trimws(as.character(x)))),
        logical(1)
    ))) {
        stop('Pando GRN design edge identifiers must be complete.', call. = FALSE)
    }
    if (nrow(edges)) {
        expected_edge_id <- paste(
            as.character(edges$tf),
            as.character(edges$atac_feature_id),
            as.character(edges$target),
            sep = '::'
        )
        if (!identical(as.character(edges$edge_id), expected_edge_id) ||
            anyDuplicated(edges$edge_id)) {
            stop(
                'Pando GRN design edge IDs must uniquely encode TF, ATAC feature, and target.',
                call. = FALSE
            )
        }
        index <- suppressWarnings(as.integer(edges$candidate_index))
        if (anyNA(index) || !identical(index, seq_len(nrow(edges)))) {
            stop('Pando GRN design candidate indices must be sequential.', call. = FALSE)
        }
        support <- .pando_split_supporting_regions(edges$supporting_regions)
        support_n <- vapply(support, length, integer(1))
        if (!identical(as.integer(edges$n_supporting_regions), support_n) ||
            any(!mapply(function(region, values) region %in% values,
                        as.character(edges$region), support))) {
            stop('Pando supporting-region provenance is inconsistent.', call. = FALSE)
        }
    }
    contract <- design$feature_contract
    required_contract <- c(
        'cell_ids', 'rna_feature_ids', 'atac_feature_ids',
        'rna_assay', 'atac_assay', 'rna_layer', 'atac_layer'
    )
    if (!is.list(contract) || !all(required_contract %in% names(contract))) {
        stop('Pando GRN design feature contract is incomplete.', call. = FALSE)
    }
    if (anyDuplicated(contract$cell_ids) ||
        anyDuplicated(contract$rna_feature_ids) ||
        anyDuplicated(contract$atac_feature_ids)) {
        stop('Pando GRN design feature-contract identifiers must be unique.', call. = FALSE)
    }
    if (nrow(edges)) {
        missing_rna <- setdiff(unique(c(
            as.character(edges$tf_feature_id),
            as.character(edges$target_feature_id)
        )), as.character(contract$rna_feature_ids))
        missing_atac <- setdiff(
            unique(as.character(edges$atac_feature_id)),
            as.character(contract$atac_feature_ids)
        )
        if (length(missing_rna) || length(missing_atac)) {
            stop('Candidate edges reference features absent from the feature contract.',
                 call. = FALSE)
        }
    }
    region_map <- design$region_map
    if (!is.data.frame(region_map) ||
        !all(c('region', 'atac_feature_id') %in% colnames(region_map)) ||
        anyNA(region_map$region) || anyNA(region_map$atac_feature_id) ||
        any(!nzchar(as.character(region_map$region))) ||
        any(!nzchar(as.character(region_map$atac_feature_id))) ||
        anyDuplicated(as.character(region_map$region))) {
        stop('Pando GRN design region map is incomplete or ambiguous.', call. = FALSE)
    }
    if (nrow(edges)) {
        mapped <- as.character(region_map$atac_feature_id[
            match(as.character(edges$region), as.character(region_map$region))
        ])
        if (anyNA(mapped) || !identical(mapped, as.character(edges$atac_feature_id))) {
            stop('Representative regulatory regions do not map to their recorded ATAC features.',
                 call. = FALSE)
        }
        support <- unique(unlist(
            .pando_split_supporting_regions(edges$supporting_regions),
            use.names = FALSE
        ))
        if (length(setdiff(support, as.character(region_map$region)))) {
            stop('Supporting regulatory regions are absent from the region map.',
                 call. = FALSE)
        }
    }
    if (identical(schema, 'pando_grn_design_v2')) {
        expected <- paste0(
            'md5:', .pando_design_md5(.pando_design_payload(design))
        )
        if (!identical(as.character(design$design_fingerprint), expected)) {
            stop('Pando GRN design fingerprint does not match its contents.',
                 call. = FALSE)
        }
    }
    invisible(TRUE)
}
