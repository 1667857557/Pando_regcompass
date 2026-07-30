![Build](https://github.com/quadbiolab/Pando/workflows/build/badge.svg?branch=main)

# Pando <img src="man/figures/logo.png" align="right" width="180"/>

Pando uses paired single-cell RNA and chromatin-accessibility measurements to
infer gene regulatory networks through TF-expression by motif-bearing-region
accessibility predictors.

This RegCompass fork exposes two analysis routes:

1. the unchanged original Pando `infer_grn()` workflow;
2. `infer_condition_grn()` for condition-specific sub-GRNs on one shared
   single-cell coordinate system.

## Installation

```r
remotes::install_github("1667857557/Pando_regcompass")
```

## Original Pando workflow

```r
library(Pando)
library(Seurat)
library(BSgenome.Hsapiens.UCSC.hg38)

data(motifs)
seurat_object <- Seurat::FindVariableFeatures(
  seurat_object,
  assay = "RNA"
)
grn_object <- initiate_grn(seurat_object)
grn_object <- find_motifs(
  grn_object,
  pfm = motifs,
  genome = BSgenome.Hsapiens.UCSC.hg38
)
grn_object <- infer_grn(grn_object)
coef(grn_object)
```

## Condition-specific sub-GRNs

```r
condition_object <- infer_condition_grn(
  grn_object,
  cell_type_col = "cell_type",
  condition_col = "condition",
  cell_type = c("T_cell", "B_cell"),
  genes = metabolic_genes,
  candidate_screen = "motif_domain",
  comparison_conditions = c("Control", "Drug"),
  condition_mix = 0.5,
  condition_weight = "equal",
  outer_nfolds = 5L,
  inner_nfolds = 5L,
  scale = TRUE
)
```

Each selected broad cell type is fitted independently. Within one cell type,
the condition model uses:

- the same motif/domain TF–peak–target candidate dictionary;
- an outer-training transform whose center is the mean of condition means;
- an outer-training scale equal to the square root of the mean
  condition-specific population variances;
- condition-specific sparse support;
- nested condition-stratified cell OOF for preprocessing, lambda selection and
  reliability;
- a support-constrained common-metric refit.

Different conditions may therefore have different active TF, peak, target and
edge sets, different effect magnitudes, exact zeros, or opposite directions
without changing coefficient units.

### Standard Pando-compatible output

Each result is still a standard Pando `Network`:

```r
Params(condition_object)$condition_network_index

drug_network <- GetNetwork(
  condition_object,
  network = "condition_grn__T_cell__condition__Drug"
)
coef(drug_network)
gof(drug_network)
NetworkParams(drug_network)
```

Condition coefficient tables retain the original columns
`tf`, `target`, `region`, `term`, `estimate` and `corr`. A shared network is
named, for example, `condition_grn__T_cell__shared`. In a condition network, an
unavailable edge has `estimate = NA`; an estimable inactive edge has
`estimate = 0`.

### Lossless absolute-condition contract

```r
fit <- condition_grn_fit(
  condition_object,
  network_name = "condition_grn",
  cell_type = "T_cell"
)

fit$beta_condition_std
fit$beta_shared_std
fit$delta_condition_std
fit$active_mask
fit$estimability_mask
fit$absolute_direction

condition_grn_subgraph(fit, "Drug")
```

`beta_condition` is the absolute regulatory effect on the common
within-cell-type coordinate. The public and stored contract contains no
reference-condition coefficient, reference contrast, comparison mask, or
contrast helper. Conditions are compared directly by their absolute
coefficients or their condition-specific OOF projections. Unavailable model
coefficients remain `NA` in the fit contract; estimable inactive coefficients
are exact zero.

### RegCompass handoff

```r
projection <- project_condition_grn_cells(
  condition_object,
  fit = fit,
  component = "condition",
  scale = "std",
  targets = metabolic_genes,
  nonestimable = "structural_zero",
  support_policy = "pairwise_common",
  comparison_conditions = c("Control", "Drug"),
  origin = "oof"
)
```

Pando reconstructs each interaction on the paired single cells using the stored
outer-fold training transform and coefficient, then returns signed held-out
cell-by-target regulatory scores. Every fitted cell is assigned exactly one
outer-fold prediction. `origin = "full_fit"` is interpretation-only and is
explicitly ineligible for penalty construction.

At the projection-contribution layer, an edge that is not estimable under the
requested support policy contributes exactly zero. It is retained as an
auditable structural zero in `edge_structural_zero_mask`; target scores remain
finite and structural zeros enter metacell means and downstream RegCompass
analysis. RegCompass must aggregate the cell-first scores within condition ×
broad cell type, without recomputing TF×ATAC from metacell averages,
renormalizing by condition, or refitting coefficients.

## Citation

For the Pando method, cite:

Fleck, J.S., Jansen, S.M.J., Wollny, D. et al. Inferring and perturbing cell
fate regulomes in human brain organoids. Nature 621, 365–372 (2023).
https://doi.org/10.1038/s41586-022-05279-8
