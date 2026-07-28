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
  genes = metabolic_genes,
  candidate_screen = "pooled_within_condition",
  condition_mix = 0.5,
  condition_weight = "equal",
  reference_condition = "Control",
  scale = TRUE
)
```

Every condition of one cell type uses:

- the same motif/domain TF–peak–target candidate dictionary;
- the same pooled `RNA_TF × ATAC_peak` center and scale;
- the same pooled target scale;
- condition-specific sparse support;
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
  network = "condition_grn__Tumor_cell__condition__Drug"
)
coef(drug_network)
gof(drug_network)
NetworkParams(drug_network)
```

Condition coefficient tables retain the original columns
`tf`, `target`, `region`, `term`, `estimate` and `corr`. The shared network is
named `condition_grn__<cell_type>__shared`. In a condition network, an
unavailable edge has `estimate = NA`; an estimable inactive edge has
`estimate = 0`.

### Lossless condition contract

```r
fit <- condition_grn_fit(
  condition_object,
  network_name = "condition_grn",
  cell_type = "Tumor_cell"
)

fit$beta_condition_std
fit$beta_shared_std
fit$delta_condition_std
fit$active_mask
fit$estimability_mask
fit$absolute_direction

condition_grn_subgraph(fit, "Drug")
condition_grn_contrast(fit, "Control", "Drug")
```

`beta_condition` is the absolute regulatory effect. `delta_condition` is the
deviation from the shared effect and must not be interpreted as the absolute
direction. Unavailable effects are `NA`; estimable inactive effects are exact
zero.

### RegCompass handoff

```r
projection <- project_condition_grn_cells(
  condition_object,
  fit = fit,
  component = "condition",
  scale = "std",
  targets = metabolic_genes,
  nonestimable = "propagate"
)
```

Pando reconstructs each interaction on the paired single cells using the stored
pooled transform, then returns signed cell-by-target regulatory scores.
RegCompass should aggregate these scores within
cell type × condition × sample/donor after projection. It must not recompute
TF×ATAC from metacell averages, renormalize by condition, refit coefficients or
replace unavailable values with zero.

## Citation

For the Pando method, cite:

Fleck, J.S., Jansen, S.M.J., Wollny, D. et al. Inferring and perturbing cell
fate regulomes in human brain organoids. Nature 621, 365–372 (2023).
https://doi.org/10.1038/s41586-022-05279-8
