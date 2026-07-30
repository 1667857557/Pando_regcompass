![Build](https://github.com/quadbiolab/Pando/workflows/build/badge.svg?branch=main)

# Pando

This fork retains the original `infer_grn()` workflow and adds
`infer_condition_grn()` for condition-specific TF–peak–target effects within each
broad cell type.

## Install

```r
remotes::install_github("1667857557/Pando_regcompass")
```

## Condition-aware fit

```r
condition_object <- infer_condition_grn(
  grn_object,
  cell_type_col = "cell_type",
  condition_col = "condition",
  genes = metabolic_genes,
  candidate_screen = "motif_domain",
  condition_weight = "equal",
  outer_nfolds = 5L,
  inner_nfolds = 5L,
  scale = TRUE
)

fit <- condition_grn_fit(
  condition_object,
  network_name = "condition_grn",
  cell_type = "T_cell"
)
```

Each broad cell type is fitted separately. Conditions share the candidate
supergraph, coefficient coordinate, target columns and fold-local transforms.
They may have different active edges, exact zeros and opposite directions.

## Primary RegCompass handoff

```r
projection <- project_condition_grn_primary_cells(
  condition_object,
  fit = fit,
  scale = "std",
  targets = metabolic_genes
)
```

The primary route is `condition-full OOF`:

- an edge estimable in the focal condition contributes its outer-heldout
  condition coefficient;
- jointly estimable edges form the common-support component;
- an edge non-estimable in one or both conditions contributes exactly zero in
  each non-estimable condition;
- an exact-zero predictor, including a peak closed in every input cell, remains
  in the shared candidate supergraph as a projectable structural-zero edge and
  is not assigned a fitted coefficient.

The common-support decomposition remains available through
`project_condition_grn_cells(..., support_policy = "pairwise_common")` or
`"global_common"`. It is not the primary penalty.

Project paired cells before exact metacell aggregation:

```r
metacell_projection <- aggregate_condition_grn_projection(
  projection,
  membership,
  group_col = "metacell_id"
)
```

Do not recompute TF×ATAC from metacell averages, renormalize separately by
condition, or refit after aggregation.

All equations are maintained in one location:
[RegCompass mathematical model](https://github.com/1667857557/Regcompass/blob/Main/docs/mathematical-model.md).

## Citation

Fleck, J.S., Jansen, S.M.J., Wollny, D. et al. Inferring and perturbing cell
fate regulomes in human brain organoids. *Nature* 621, 365–372 (2023).
