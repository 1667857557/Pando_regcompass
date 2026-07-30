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

Conditions share the candidate dictionary, coefficient coordinate, target
columns and fold-local transforms. They may still have different active edges,
different node sets, exact zeros and opposite coefficient directions.

## RegCompass handoff

The primary cross-condition route remains common support:

```r
primary <- project_condition_grn_cells(
  condition_object,
  fit = fit,
  component = "condition",
  scale = "std",
  support_policy = "pairwise_common",
  comparison_conditions = c("Control", "Drug"),
  origin = "oof"
)
```

For more than two conditions, use `support_policy = "global_common"`.

Condition-unique estimable edges are retained in a separate supplemental route:

```r
supplemental <- project_condition_grn_supplemental_cells(
  condition_object,
  fit = fit,
  scale = "std",
  targets = metabolic_genes
)
```

The supplemental route is eligible for a secondary RegCompass penalty and its
increment relative to the common-support penalty. It does not replace the
common-support primary ranking or primary condition statistics.

Both routes are projected on paired single cells before exact metacell
aggregation:

```r
primary_metacell <- aggregate_condition_grn_projection(
  primary, membership, group_col = "metacell_id"
)
```

Do not recompute TF×ATAC from metacell averages, renormalize separately by
condition, or refit coefficients after aggregation.

All mathematical definitions are maintained in one location:
[RegCompass mathematical model](https://github.com/1667857557/Regcompass/blob/Main/docs/mathematical-model.md).

## Citation

Fleck, J.S., Jansen, S.M.J., Wollny, D. et al. Inferring and perturbing cell
fate regulomes in human brain organoids. *Nature* 621, 365–372 (2023).
