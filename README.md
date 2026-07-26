![Build](https://github.com/quadbiolab/Pando/workflows/build/badge.svg?branch=main)

# Pando <img src="man/figures/logo.png" align="right" width="180"/>

Pando uses paired single-cell RNA and chromatin-accessibility measurements to infer gene regulatory networks through TF-expression by motif-bearing-region accessibility predictors.

This RegCompass fork preserves the original Pando workflow and adds a public pre-fit structural-design contract in version 1.1.2.

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

`prepare_grn_design()` exposes the TF–peak–target candidate universe before coefficient fitting. This is intended for downstream joint or multitask models that require every condition to use the same edge dictionary.

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
- returns one deterministic edge ID and candidate order;
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

No condition metadata is interpreted by Pando. The caller is responsible for using the same `PandoGRNDesign` in every task or condition model.

## Citation

For the Pando method, cite:

Fleck, J.S., Jansen, S.M.J., Wollny, D. et al. Inferring and perturbing cell fate regulomes in human brain organoids. Nature 621, 365–372 (2023). https://doi.org/10.1038/s41586-022-05279-8
