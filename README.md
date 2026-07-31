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

## Numerical core

The canonical condition workflow contains only the main nested-CV and final
support-constrained refit. Bootstrap stability, ridge-grid sensitivity and other
sensitivity refits are not run or stored.

The TF-by-ATAC product, sparse-group lambda path and support-constrained refit are
all compiled C++17/RcppEigen kernels. Missing native registration, unsupported
matrix layouts, invalid dimensions, failed factorization, non-finite arithmetic
or failed residual verification stop the analysis immediately. There is no
runtime R fallback.

Predictor scaling remains sparse; equal-condition centering is handled
algebraically by the condition intercept and the projection shift. Fold-level
centered Gram matrices and response cross-products are cached. Each refit is
sent to the native path-capable kernel with the cached sufficient statistics,
preserving the original double-precision equations and output schema while
removing R-level Schur construction and matrix slicing.

The alternating R refit remains available only as the internal numerical oracle
`Pando:::.condition_refit_shared_baseline_reference()` for package regression
tests. It is not a selectable analysis backend.

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

- an edge estimable in the focal condition contributes its outerheldout
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
[RegCompass Tutorial 3: mathematical model](https://github.com/1667857557/Regcompass/blob/Main/docs/tutorial-03-mathematical-model.md).

## Citation

Fleck, J.S., Jansen, S.M.J., Wollny, D. et al. Inferring and perturbing cell
fate regulomes in human brain organoids. *Nature* 621, 365–372 (2023).
