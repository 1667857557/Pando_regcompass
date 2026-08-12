![Build](https://github.com/quadbiolab/Pando/workflows/build/badge.svg?branch=main)

# Pando

This fork retains the original `infer_grn()` workflow and provides one canonical condition-comparable extension that preserves Pando's regulatory-domain, motif, peak-target-correlation, TF-target-correlation, and TF×ATAC predictor logic while using ridge regularization for a shared condition dictionary.

## Install

```r
remotes::install_github("1667857557/Pando_regcompass")
```

## Multiple conditions: pooled/global + condition common dictionary

```r
library(BiocParallel)

upstream_bp <- if (.Platform$OS.type == "windows") {
  SnowParam(workers = 6L, type = "SOCK")
} else {
  MulticoreParam(workers = 6L)
}

condition_object <- infer_condition_grn(
  grn_object,
  cell_type_col = "cell_type",
  condition_col = "condition",
  genes = metabolic_genes,
  tf_cor = 0.05,
  peak_cor = 0.05,
  min_cells_per_condition = 300L,
  padj_threshold = 0.05,
  parallel = TRUE,
  BPPARAM = upstream_bp,
  parallel_scope = "cell_type"
)

fit <- condition_grn_fit(
  condition_object,
  cell_type = "T_cell"
)
```

`tf_cor` and `peak_cor` both default to `0.05`. These are adjustable candidate-discovery thresholds, not fixed constants. The same supplied values are used for the pooled/global calculation and for every condition-specific calculation.

With multiple requested broad cell types, `parallel_scope = "cell_type"` runs independent cell-type jobs through `BPPARAM`. `parallel_scope = "target"` uses target-level parallel work. The default `"auto"` chooses the available non-nested route.

For every broad cell type, the canonical condition workflow performs these steps:

1. restrict the structural candidate space by Pando regulatory domains, measured regulatory regions, and TF motif support;
2. compute Pando `peak_cor` and `tf_cor` candidate support using all eligible-condition cells pooled together;
3. compute the same candidate support independently within every condition;
4. form `unique(D_global ∪ D_condition1 ∪ ...)` on the exact `(target, TF, region)` triple;
5. freeze that deduplicated target-specific dictionary;
6. build the original Pando interaction predictor `RNA(TF) × ATAC(peak)` for every dictionary edge in every condition;
7. use one common predictor-scaling convention and one target-specific cross-validated ridge `lambda`, but fit independent condition coefficient blocks with **no cross-condition fusion penalty**;
8. BH-adjust the approximate ridge-Wald P values within condition;
9. mark an edge active in condition `c` when that condition's coefficient is estimable and BH-supported;
10. retain `global_support`, `local_support`, and the corresponding correlations as candidate provenance only;
11. expose `penalty_effect = estimate` for active edges and zero otherwise, while retaining the complete coefficient table for diagnostics and direct condition contrasts.

The common dictionary therefore provides comparable coefficient coordinates without forcing condition coefficients to be similar. Ridge stabilizes the penalized system under severe collinearity and raw rank deficiency, but highly correlated predictors can still make individual biological attribution uncertain.

Pooled/global and condition-specific screening are used only to decide which exact edges are allowed into the common model space. This is important for small conditions: an edge can enter through pooled/global support or through another condition and still be estimated in every condition. If condition `B` has weak or noisy marginal `tf_cor`/`peak_cor` but its own ridge coefficient is estimable and BH-supported, the edge is allowed to be active in `B`. Thus the correlation screen is not applied a second time after ridge fitting.

## Candidate dictionary inspection

Candidate discovery and exact-union helpers can be used to inspect dictionary provenance directly:

```r
edges_global <- discover_grn_edges(
  grn_object,
  genes = metabolic_genes,
  cells = cells_of_one_cell_type,
  source_label = "global",
  source_type = "global",
  tf_cor = 0.05,
  peak_cor = 0.05
)

edges_condition <- lapply(condition_levels, function(level) {
  discover_grn_edges(
    grn_object,
    genes = metabolic_genes,
    cells = cells_by_condition[[level]],
    source_label = level,
    source_type = "condition",
    tf_cor = 0.05,
    peak_cor = 0.05
  )
})
names(edges_condition) <- condition_levels

edge_dictionary <- union_grn_edges(edges_global, edges_condition)
stopifnot(!anyDuplicated(edge_dictionary$edge_id))
```

The canonical coefficient-fitting API for multiple conditions is `infer_condition_grn()`. Candidate tables are provenance/inspection objects and are not substituted for the final condition-specific ridge coefficients.

## No condition or one condition

When `condition_col` is `NULL`, absent, or has fewer than two observed levels, `infer_condition_grn()` uses the standard Pando route rather than creating a multi-condition fit contract.

Condition-GRN controls are not forwarded blindly into unrelated standard model backends. This prevents condition-only arguments from producing repeated per-target model failures.

## RegCompass handoff

```r
projection <- project_condition_grn_cells(
  condition_object,
  fit = fit,
  significant_only = TRUE
)

metacell_projection <- aggregate_condition_grn_projection(
  projection,
  membership,
  group_col = "metacell_id"
)
```

Pando's projection reconstructs the paired-cell predictor `TF RNA × peak ATAC` from the fitted assay layers and applies `penalty_effect`. RegCompass Layer 1 intentionally uses a different downstream metacell estimand, `beta × mean(TF) × mean(ATAC)`, using the same active condition coefficients. Do not confuse that product of means with `beta × mean(TF × ATAC)`.

The reported coefficient and contrast P values are approximate and conditional on data-dependent candidate discovery and CV-selected ridge regularization; the workflow does not claim exact selective inference.

## Citation

Fleck, J.S., Jansen, S.M.J., Wollny, D. et al. Inferring and perturbing cell fate regulomes in human brain organoids. *Nature* 621, 365–372 (2023).
