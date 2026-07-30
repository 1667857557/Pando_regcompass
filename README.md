![Build](https://github.com/quadbiolab/Pando/workflows/build/badge.svg?branch=main)

# Pando

Pando uses paired single-cell RNA and chromatin-accessibility measurements to
infer TF–peak–target regulatory relationships.

This RegCompass fork keeps two routes:

1. the original `infer_grn()` workflow;
2. `infer_condition_grn()` for condition-specific sub-GRNs fitted independently
   within each broad cell type.

## Installation

```r
remotes::install_github("1667857557/Pando_regcompass")
```

## Original workflow

```r
condition_object <- initiate_grn(seurat_object)
condition_object <- find_motifs(
  condition_object,
  pfm = motifs,
  genome = BSgenome.Hsapiens.UCSC.hg38
)
condition_object <- infer_grn(condition_object)
```

## Condition-aware workflow

```r
condition_object <- infer_condition_grn(
  condition_object,
  cell_type_col = "cell_type",
  condition_col = "condition",
  genes = metabolic_genes,
  candidate_screen = "motif_domain",
  condition_mix = 0.5,
  condition_weight = "equal",
  outer_nfolds = 5L,
  inner_nfolds = 5L,
  scale = TRUE
)
```

Each broad cell type is fitted separately. Conditions share the candidate edge
dictionary, coefficient coordinate, fold-local transforms and target columns,
but may have different active edges, different node sets, exact zeros and
opposite coefficient directions.

The canonical fit extractor is:

```r
fit <- condition_grn_fit(
  condition_object,
  network_name = "condition_grn",
  cell_type = "T_cell"
)
```

`fit$schema_version` is `pando_condition_grn_fit`. Unavailable coefficients
remain `NA`; estimable inactive coefficients are numeric zero.

## Primary RegCompass handoff

```r
projection <- project_condition_grn_primary_cells(
  condition_object,
  fit = fit,
  scale = "std",
  targets = metabolic_genes,
  nonestimable = "structural_zero"
)
```

The primary projection uses each condition's own outer-fold estimability mask.
A condition-unique edge therefore contributes to that condition and is a
structural zero where it is not estimable. Common-edge support is not required.
The shared candidate dictionary, shared units, fold-local transforms and exact
single-cell membership preserve the coordinate system used for downstream
condition comparison.

The projection is computed on paired single cells before aggregation:

```r
metacell_projection <- aggregate_condition_grn_projection(
  projection,
  membership,
  group_col = "metacell_id"
)
```

Do not recompute TF×ATAC from metacell averages, renormalize separately by
condition, or refit coefficients after aggregation.

The lower-level `project_condition_grn_cells()` remains available for diagnostic
support policies and full-fit interpretation. RegCompass uses only
`project_condition_grn_primary_cells()` for its main penalty.

All equations are maintained in the single
[RegCompass mathematical tutorial](https://github.com/1667857557/Regcompass/blob/Main/docs/mathematical-model.md).

## Citation

Fleck, J.S., Jansen, S.M.J., Wollny, D. et al. Inferring and perturbing cell
fate regulomes in human brain organoids. *Nature* 621, 365–372 (2023).
