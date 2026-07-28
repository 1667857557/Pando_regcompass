# Condition-specific sub-GRNs

Pando exposes two analysis routes:

- `infer_grn()` is the unchanged original Pando workflow.
- `infer_condition_grn()` fits condition-specific sub-GRNs on one shared
  single-cell coordinate system.

The condition route writes one standard Pando `Network` per condition plus a
shared network. `coef()`, `gof()` and `NetworkParams()` therefore retain their
original interfaces. A versioned `ConditionGRNFit` stores additional masks,
directions, transforms and projection metadata without changing the standard
network tables.

## Shared dictionary, different active graphs

All conditions of a cell type use the same motif/domain TF–peak–target
dictionary. Candidate screening uses a weighted RMS of within-condition
correlations. Pure condition mean shifts do not generate candidates, while
opposite within-condition effects cannot cancel each other during screening.

Sparse selection and a support-constrained common-metric refit permit:

- different active TF, peak and target nodes;
- different active edges;
- exact estimated zeros;
- different effect magnitudes;
- sign reversals;
- explicit `NA` for unavailable effects.

The v3 contract exposes standardized and raw layers:

```text
beta_condition = beta_shared + delta_condition
absolute_direction = sign(beta_condition)
deviation_direction = sign(delta_condition)
```

Use `condition_grn_contrast()` for reference-free pairwise effects and
`condition_grn_subgraph()` for the active edge and node sets of one condition.

## RegCompass handoff

`project_condition_grn_cells()` reconstructs the fitted interaction from the
original paired single cells and the stored pooled transform:

```text
z_edge = (RNA_TF * ATAC_peak - pooled_center) / pooled_scale
gene_score = sum(z_edge * beta_condition)
```

The result has identical target columns for every condition, signed scores,
target/condition estimability status and an explicit aggregation contract.
Projection is restricted to the exact paired cells, assays and condition labels
stored by the fit.
RegCompass should average these single-cell target scores within
cell type × condition × sample/donor to form metacell evidence.

RegCompass must not rebuild the interaction from aggregated TF and ATAC values,
renormalize a condition, refit coefficients, use selection coefficients, or
replace unavailable values with zero.
