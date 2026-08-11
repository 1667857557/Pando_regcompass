# Bounded target-level BiocParallel routing for condition GRNs.
#
# The canonical condition-GRN direct definition remains in R/condition_grn.R.
# This adapter installs an extended method through an alias so the source tree
# keeps exactly one direct S3 method definition while target-scope BPPARAM can be
# reused across candidate discovery and multi-task ridge fitting.

.pando_condition_parallel_method_impl <- infer_condition_grn.GRNData
.pando_condition_parallel_map_impl <- map_par

.pando_condition_target_bpparam <- function() {
    value <- getOption("Pando.condition_target_BPPARAM", NULL)
    if (identical(value, FALSE) || is.null(value)) return(NULL)
    .pando_validate_bpparam(value)
    value
}

map_par <- function(x, fun, parallel = FALSE, verbose = TRUE) {
    if (!isTRUE(parallel)) {
        return(.pando_condition_parallel_map_impl(
            x = x, fun = fun, parallel = FALSE, verbose = verbose
        ))
    }
    BPPARAM <- .pando_condition_target_bpparam()
    if (is.null(BPPARAM)) {
        return(.pando_condition_parallel_map_impl(
            x = x, fun = fun, parallel = TRUE, verbose = verbose
        ))
    }
    out <- BiocParallel::bplapply(x, fun, BPPARAM = BPPARAM)
    names(out) <- names(x)
    out
}

.pando_condition_target_parallel_method <- function(
    object, cell_type_col = NULL, condition_col = NULL, cell_type = NULL,
    genes = NULL, network_name = "condition_grn",
    peak_to_gene_method = c("Signac", "GREAT"), upstream = 100000,
    downstream = 0, extend = 1000000, only_tss = FALSE,
    peak_to_gene_domains = NULL, rna_layer = "data", peak_layer = "data",
    peak_value_type = c("normalized", "probability", "other"),
    tf_cor = 0.1, peak_cor = 0,
    min_cells_per_condition = 50L,
    small_condition_action = c("error", "drop_condition", "skip_cell_type"),
    adjust_method = "BH", padj_threshold = 0.05,
    rank_action = c("mark", "error"), min_residual_df = 1L,
    parallel = FALSE, BPPARAM = NULL,
    parallel_scope = c("auto", "cell_type", "target"),
    overwrite = FALSE, fallback_args = list(), verbose = TRUE, ...) {
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

    dots <- list(...)
    args <- list(
        object = object,
        cell_type_col = cell_type_col,
        condition_col = condition_col,
        cell_type = cell_type,
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
        parallel = parallel,
        BPPARAM = BPPARAM,
        parallel_scope = parallel_scope,
        overwrite = overwrite,
        fallback_args = fallback_args,
        verbose = verbose
    )
    answer <- do.call(.pando_condition_parallel_method_impl, c(args, dots))
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
    answer
}

infer_condition_grn.GRNData <- .pando_condition_target_parallel_method
