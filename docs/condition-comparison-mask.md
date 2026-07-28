# Condition-comparison eligibility mask

`ConditionGRNFit` stores an edge-by-condition `eligibility_mask`. A coefficient fixed to zero because an edge is not estimable in one condition is not biologically equivalent to an estimated zero coefficient.

For a reference-condition contrast, Pando now also records

```text
comparison_mask[e, c] = eligibility_mask[e, c] &&
                        eligibility_mask[e, reference_condition]
```

Only entries with `comparison_mask[e, c] == TRUE` support an interpretable coefficient contrast

```text
beta[e, c] - beta[e, reference_condition]
```

The numeric `contrast` matrix is retained for backward compatibility. Downstream consumers must use `comparison_mask` to distinguish an estimated contrast from a structurally unavailable comparison.

This change does not alter coefficient estimation, the shared edge dictionary, pooled predictor scaling, or target-specific shared lambda selection. It only makes comparison support explicit and lossless.
