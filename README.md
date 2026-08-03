![Build](https://github.com/quadbiolab/Pando/workflows/build/badge.svg?branch=main)

# Pando

This fork retains the original `infer_grn()` workflow and provides one canonical
condition-specific extension based on the original Gaussian identity Pando model.

## Install

```r
remotes::install_github("1667857557/Pando_regcompass")
```

## Multiple conditions: two-stage common dictionary

```r
condition_object <- infer_condition_grn(
  grn_object,
  cell_type_col = "cell_type",
  condition_col = "condition",
  genes = metabolic_genes,
  tf_cor = 0.1,
  peak_cor = 0,
  min_cells_per_condition = 300L,
  adjust_method = "BH",
  padj_threshold = 0.05,
  rank_action = "mark",
  parallel = TRUE
)

fit <- condition_grn_fit(
  condition_object,
  cell_type = "T_cell"
)
```

For every broad cell type, Pando performs exactly these steps:

1. discover candidate edges on all eligible cells of that type;
2. discover candidate edges separately in every condition;
3. union exact `(target, TF, region)` triples;
4. freeze that target-specific edge dictionary;
5. refit every condition with the same predictor columns using
   `target ~ TF:peak`, Gaussian identity GLM and `scale = FALSE`;
6. calculate network-wide adjusted P values and expose `penalty_effect`, which is
   the fitted coefficient only when `padj < padj_threshold`.

The unfiltered condition coefficient remains in `estimate`. Non-significant
coefficients are not rewritten; only `penalty_effect` is zeroed for downstream
penalty construction. Zero-variance and aliased coefficients remain `NA` with
explicit diagnostics.

Candidate discovery uses the original Pando peak-to-gene domain, motif mapping,
peak-target correlation and TF-target correlation logic. Candidate coefficients
from this first stage are never used as final condition effects. Global
coefficients are not used to rescale condition coefficients.

## Explicit APIs

The automatic workflow is composed from three public functions:

```r
edges_global <- discover_grn_edges(
  grn_object,
  genes = metabolic_genes,
  cells = cells_of_one_cell_type,
  source_label = "global",
  source_type = "global"
)

edges_condition <- lapply(condition_levels, function(level) {
  discover_grn_edges(
    grn_object,
    genes = metabolic_genes,
    cells = cells_by_condition[[level]],
    source_label = level,
    source_type = "condition"
  )
})
names(edges_condition) <- condition_levels

edge_dictionary <- union_grn_edges(edges_global, edges_condition)

condition_object <- fit_grn_from_edges(
  grn_object,
  edge_dictionary = edge_dictionary,
  cells = cells_by_condition[["Control"]],
  condition_label = "Control",
  adjust_method = "BH",
  padj_threshold = 0.05
)
```

## No condition or one condition

When `condition_col` is `NULL`, absent, or has fewer than two observed levels,
`infer_condition_grn()` directly runs the original Pando Gaussian interaction
GRN independently for every requested broad cell type. It creates no condition
coefficient or condition fit contract.

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

Projection reconstructs the paired-cell predictor `TF RNA × peak ATAC` from the
same globally preprocessed assay layers used for fitting and applies
`penalty_effect`. Projection occurs before exact metacell aggregation. Do not
renormalize each condition, regenerate a different edge dictionary, or refit
coefficients after aggregation.

The reported P values are conditional on the selected edge dictionary; the
workflow does not claim selective-inference correction for candidate discovery.

## Citation

Fleck, J.S., Jansen, S.M.J., Wollny, D. et al. Inferring and perturbing cell
fate regulomes in human brain organoids. *Nature* 621, 365–372 (2023).
