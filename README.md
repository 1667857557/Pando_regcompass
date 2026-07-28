![Build](https://github.com/quadbiolab/Pando/workflows/build/badge.svg?branch=main)

# Pando <img src="man/figures/logo.png" align="right" width="180"/>

Pando uses paired single-cell RNA and chromatin-accessibility measurements to infer gene regulatory networks through TF-expression by motif-bearing-region accessibility predictors.

This RegCompass fork preserves the original Pando workflow and adds two public
interfaces for comparable condition models: a pre-fit structural-design
contract and a lossless post-fit `ConditionGRNFit` contract.

## Installation

```r
remotes::install_github("1667857557/Pando_regcompass")
```

## Standard Pando workflow

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

## Shared structural GRN design API

`prepare_grn_design()` exposes the TF–peak–target candidate universe before coefficient fitting. This supports downstream joint or multitask models in which every condition must use the same edge dictionary.

```r
design <- prepare_grn_design(
  grn_object,
  genes = metabolic_genes,
  peak_to_gene_method = "Signac",
  min_tf_detection = 0.01,
  min_peak_detection = 0.01,
  min_target_detection = 0.01,
  max_edges_per_target = Inf
)

validate_grn_design(design)
design$candidate_edges
design$target_diagnostics
design$feature_contract
design$design_fingerprint
```

The design API:

- uses Pando peak-to-gene domains and motif-to-TF mappings;
- contains no fitted coefficient, p-value or pooled-significance requirement;
- records the exact measured ATAC feature used for each candidate;
- deduplicates predictors by `(TF, ATAC feature, target)`;
- preserves all regulatory regions supporting a deduplicated predictor in `supporting_regions`;
- distinguishes the number of predictor columns, measured ATAC features and supporting regulatory regions per target;
- validates edge IDs, feature IDs, region-to-peak mappings and candidate order;
- records an `md5:` fingerprint of the canonical candidate, feature and parameter contract;
- leaves the original `infer_grn()` interface unchanged.

The principal candidate table fields are:

```text
edge_id
candidate_index
tf
region
supporting_regions
n_supporting_regions
target
atac_feature_id
tf_feature_id
target_feature_id
tf_detection
peak_detection
target_detection
```

The structural API does not interpret condition metadata. Use
`infer_condition_grn()` when Pando should construct and fit the comparable
condition models itself.

## Condition-comparable GRNs

The default condition-aware engine estimates each condition independently
inside a cell type, while sharing the quantities required for direct
comparison:

- one exact TF–peak–target edge dictionary;
- an edge-by-condition eligibility mask;
- pooled scaling of the final `TF expression × peak accessibility` predictor;
- one target-specific lambda path and one selected lambda shared by conditions.

```r
grn_object <- infer_condition_grn(
  grn_object,
  cell_type_col = "cell_type",
  condition_col = "condition",
  genes = metabolic_genes,
  method = "shared_design_independent",
  candidate_screen = "condition_union",
  condition_mix = 1,
  condition_weight = "equal",
  reference_condition = "Control",
  scale = TRUE
)

fit <- condition_grn_fit(
  grn_object,
  network_name = "condition_grn",
  cell_type = "Tumor_cell"
)
fit$edge_table
fit$beta
fit$contrast
fit$eligibility_mask
fit$predictor_transform
fit$target_rsq
```

`fit$contrast` is exactly
`beta_condition - beta_reference`. Candidate union is performed at edge level:
the TF and peak must both pass inside the same condition, so Pando never creates
an edge by combining a TF retained only in one condition with a peak retained
only in another. Standard Pando `Network` objects are still written for
compatibility. Their Universal coefficient is an equal-condition summary for
visualization, not the condition-effect baseline.

## Citation

For the Pando method, cite:

Fleck, J.S., Jansen, S.M.J., Wollny, D. et al. Inferring and perturbing cell fate regulomes in human brain organoids. Nature 621, 365–372 (2023). https://doi.org/10.1038/s41586-022-05279-8
