# Shared-baseline condition-specific sub-GRNs

## Scope

This implementation fits condition-aware TF–peak–target models exclusively on
paired RNA+ATAC **single cells**. It does not construct or consume metacells.
Metacells remain a downstream RegCompass responsibility after Pando has fixed
the candidate graph, active support, effect sizes, directions and transforms.

## Model architecture

Each cell type uses one shared candidate TF–peak–target supergraph. Candidate
edges are defined by the peak-to-gene domain and TF motif mapping. The default
`pooled_within_condition` screen removes condition-specific means before
calculating the screening correlation, so a pure condition shift cannot create
an apparent regulatory edge.

For condition `c`, the absolute edge effect is represented as

```text
beta_condition[e, c] = beta_shared[e] + delta_condition[e, c]
```

The implementation uses two stages.

### 1. Sparse-group support selection

The existing multitask sparse-group solver selects an edge-by-condition
coefficient matrix on the pooled single-cell predictor scale:

```text
loss(B)
+ lambda * alpha * (1 - condition_mix) * sum_e ||B[e, ]||_2
+ lambda * alpha * condition_mix * sum_e,c |B[e, c]|
+ lambda * (1 - alpha) / 2 * ||B||_F^2
```

The row-group term favours shared edge support. The elementwise term permits a
candidate edge to be active in only a subset of conditions. Coefficients remain
unconstrained in sign, so the same TF–peak–target edge may be positive in one
condition and negative in another.

### 2. Support-constrained hierarchical-ridge refit

Regularized selection coefficients are not exported as the primary biological
effect. For every target, Pando fixes the selected edge-by-condition support and
refits condition coefficients on the common pooled predictor and response scale:

```text
sum_c loss_c(beta_c)
+ rho / 2 * sum_c ||beta_c - beta_shared||_2^2
```

Inactive but estimable coefficients remain exactly zero. Structurally or
numerically unavailable coefficients are represented as `NA` in
`beta_condition` and by `estimability_mask == FALSE`. The compatibility `beta`
matrix keeps unavailable entries at zero and must be interpreted together with
the masks.

The shared coefficient is the estimability-aware weighted mean of the refitted
condition coefficients. Deviations are then defined as

```text
delta_condition[e, c] = beta_condition[e, c] - beta_shared[e]
```

## Different nodes and sign reversals

The candidate supergraph is shared, but each condition sub-GRN is induced by its
own active edge mask:

```text
active_mask[e, c] = TRUE
```

Consequently, conditions may have different active TF nodes, peaks, targets and
edges. The model also permits:

- different magnitudes with the same sign;
- one condition active and another estimated zero;
- positive-to-negative or negative-to-positive sign reversals;
- an edge estimable in one condition but unavailable in another.

## Output contract

The new engine writes `pando_condition_grn_fit_v3` with the following principal
fields:

- `edge_table`: shared candidate TF–peak–target dictionary;
- `beta_selection`: sparse-group selection coefficients;
- `beta`: refitted compatibility matrix, unavailable entries stored as zero;
- `beta_condition`: refitted absolute effects, unavailable entries stored as
  `NA`;
- `beta_shared`: shared baseline effects;
- `delta_condition`: condition deviations from the shared baseline;
- `eligibility_mask` and `estimability_mask`;
- `active_mask`;
- `contrast` and `comparison_mask`;
- `predictor_transform` and `response_transform`;
- target-level selection and refit diagnostics.

Direction semantics are explicit:

```text
absolute direction = sign(beta_condition)
deviation direction = sign(delta_condition)
pairwise direction = sign(beta_condition[, c2] - beta_condition[, c1])
```

A negative deviation is not automatically a negative absolute regulatory effect.
It can represent a positive edge that is weaker than the shared baseline.

## Recommended call

```r
object <- infer_condition_grn(
  object,
  cell_type_col = "cell_type",
  condition_col = "condition",
  genes = metabolic_genes,
  candidate_screen = "pooled_within_condition",
  method = "shared_baseline_condition_sparse",
  alpha = 0.5,
  condition_mix = 0.5,
  condition_weight = "equal",
  scale = TRUE,
  reference_condition = "Control"
)

fit <- condition_grn_fit(object, cell_type = "Tumor_cell")
fit$beta_shared
fit$beta_condition
fit$delta_condition
fit$active_mask
fit$estimability_mask
```

## RegCompass handoff

RegCompass should consume the refitted coefficient layers and masks, not
`beta_selection`. Metacells must be constructed only after the Pando fit is
complete. RegCompass must not reinterpret an unavailable coefficient as a
biological zero and must not refit the Pando coefficients on metacells.

The exact single-cell regulatory projection remains a separate downstream
interface task. This pull request establishes the coefficient and mask contract
required for that projection without introducing metacells into Pando.