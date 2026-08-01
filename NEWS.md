# Pando 1.6.0

- Moves the complete numerical workflow for each condition-aware target into one
  registered C++17/RcppEigen engine call. Candidate-edge construction and exact
  R fold labels remain outside the engine; outer/inner nested CV, fold-local
  transforms, lambda construction and warm-start paths, fixed-support refits,
  held-out prediction/projection, final full-data selection/refit, and OOF
  validation execute natively without returning through the R interpreter.
- Preserves the existing condition-stratified folds, equal-condition transform,
  lambda grid and selection rule, sparse-group objective, support threshold,
  ridge definition, direct-Schur equations, intercept recovery, estimability,
  structural-zero semantics, condition-specific signs, and output schema.
- Reuses one fold-local centered Gram/RHS cache across the complete inner lambda
  refit path. Invalid folds, dimensions, registration, arithmetic, Cholesky
  systems, residual checks, or OOF assignment stop immediately; no R numerical
  fallback is selected by the canonical target path.
- Retains the previous R nested-CV functions and alternating refit only as
  explicit numerical regression oracles. `.condition_fit_target()` contains one
  native target-engine call and no R-level CV, lambda-path, refit, or validation
  loop.
- Adds native-versus-reference tests for fold-selected lambda, CV loss, OOF
  prediction and projections, final sparse coefficients, direct-Schur refit,
  opposite condition directions, condition-specific structural zeros, and
  automatic lambda construction.
- Bumps the native condition ABI to 4 and publishes the fused target-engine
  backend in package metadata for strict downstream validation.

# Pando 1.5.2

- Replaces the canonical R direct-Schur support-constrained refit with a registered C++17/RcppEigen double-precision kernel; there is no runtime R fallback.
- Moves the complete fixed-support direct-Schur refit into a native path-capable kernel and reuses the fold-level Gram/RHS cache already constructed by the nested-CV workflow.
- Preserves the existing centered sufficient statistics, support rule, equal-condition weights, ridge definition, shared-baseline Schur equations, intercept recovery, structural-zero semantics and output schema.
- Stops immediately on missing native registration, malformed caches, non-finite inputs, failed Cholesky systems or failed residual verification.
- Retains the alternating R implementation only as an explicitly named numerical test oracle and adds direct native-versus-reference equivalence coverage.
- Bumps the native condition ABI to 3 so downstream packages can require the compiled refit contract.

# Pando 1.5.1

- Reads condition-model sparse inputs directly from canonical `dgCMatrix` slots instead of coercing S4 objects through `Rcpp::NumericMatrix` or `Rcpp::as<Eigen::MappedSparseMatrix<double>>`.
- Uses only the registered C++ TF-by-ATAC product kernel and C++ multitask solver. Missing symbols, unsupported sparse layouts, and native execution errors stop the analysis immediately.
- Initializes native FISTA from centered per-condition spectral-norm estimates rather than the summed Frobenius-norm upper bound.
- Adds installed-package CI for R 4.5.1, Matrix 1.7.4, native `.Call` registration, strict sparse type validation, the product kernel, and the multitask solver.
- Bumps the package version so installations from before and after the native ABI change are distinguishable.
