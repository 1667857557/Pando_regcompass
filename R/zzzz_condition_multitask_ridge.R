# Load last so the public API keeps the existing exact-union discovery path
# while replacing only its multi-condition estimator.

.condition_legacy_infer_condition_grn_one <- .pando_infer_condition_grn_one
.condition_legacy_infer_condition_grn_method <- infer_condition_grn.GRNData
.condition_ridge_control_key <- ".pando_multitask_ridge_control"

.pando_infer_condition_grn_one <- function(..., fallback_args = list()) {
    args <- list(...)
    requested_rank_action <- if ("rank_action" %in% names(args)) {
        match.arg(args$rank_action, c("mark", "error"))
    } else {
        "mark"
    }
    control <- fallback_args[[.condition_ridge_control_key]]
    fallback_args[[.condition_ridge_control_key]] <- NULL
    control <- .condition_ridge_control(control)

    legacy_args <- args
    legacy_args$rank_action <- "mark"
    legacy_args$fallback_args <- fallback_args
    object <- do.call(
        .condition_legacy_infer_condition_grn_one,
        legacy_args
    )
    if (!identical(object@grn@params$analysis_mode, "condition_grn")) {
        return(object)
    }

    arg <- function(name, default) {
        if (name %in% names(args)) args[[name]] else default
    }
    fits <- object@grn@params$condition_grn_fits
    if (!is.list(fits) || !length(fits)) {
        stop("Legacy exact-union condition discovery returned no fit contract.",
             call. = FALSE)
    }
    prepared <- .condition_prepare_common_input(
        object = object,
        genes = unique(unlist(lapply(fits, `[[`, "target_genes"),
                              use.names = FALSE)),
        peak_to_gene_method = arg("peak_to_gene_method", "Signac"),
        upstream = arg("upstream", 100000),
        downstream = arg("downstream", 0),
        extend = arg("extend", 1000000),
        only_tss = arg("only_tss", FALSE),
        peak_to_gene_domains = arg("peak_to_gene_domains", NULL),
        rna_layer = arg("rna_layer", "data"),
        peak_layer = arg("peak_layer", "data"),
        peak_value_type = arg("peak_value_type", "normalized"),
        verbose = arg("verbose", TRUE)
    )

    for (cell_type in names(fits)) {
        refitted <- .condition_ridge_refit_contract(
            object = object,
            fit = fits[[cell_type]],
            prepared = prepared,
            control = control,
            rank_action = requested_rank_action,
            min_residual_df = arg("min_residual_df", 1L),
            parallel = arg("parallel", FALSE),
            verbose = arg("verbose", TRUE)
        )
        object <- refitted$object
        fits[[cell_type]] <- refitted$fit
    }

    object@grn@params$condition_grn_fits <- fits
    object@grn@params$condition_grn_model_schema <-
        .condition_multitask_ridge_schema
    object@grn@params$condition_grn_method <-
        "two_stage_exact_edge_union_multitask_ridge"
    object@grn@params$condition_ridge_control <- control

    index <- object@grn@params$condition_network_index
    if (is.data.frame(index) && nrow(index)) {
        for (i in seq_len(nrow(index))) {
            fit <- fits[[as.character(index$cell_type[[i]])]]
            if (is.null(fit)) next
            condition <- as.character(index$condition[[i]])
            one <- fit$coefficients[fit$coefficients$condition == condition,
                                    , drop = FALSE]
            index$n_significant_edges[[i]] <- sum(one$significant %in% TRUE)
        }
        object@grn@params$condition_network_index <- index
    }
    object
}

#' @rdname infer_condition_grn
#' @method infer_condition_grn GRNData
#' @export
infer_condition_grn.GRNData <- function(object, ..., ridge_control = list()) {
    dots <- list(...)
    fallback_args <- if ("fallback_args" %in% names(dots)) {
        dots$fallback_args
    } else {
        list()
    }
    if (!is.list(fallback_args)) {
        stop("`fallback_args` must be a list.", call. = FALSE)
    }
    fallback_args[[.condition_ridge_control_key]] <-
        .condition_ridge_control(ridge_control)
    dots$fallback_args <- fallback_args
    do.call(
        .condition_legacy_infer_condition_grn_method,
        c(list(object = object), dots)
    )
}
