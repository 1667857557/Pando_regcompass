# Pando 2.0.0

- Replaces the retired condition-aware nested-CV, sparse-group, native solver,
  structural-zero and OOF projection engines with one two-stage
  common-dictionary workflow.
- Discovers exact TF-peak-target candidates in the complete cell type and in
  each condition, freezes their exact union, and fits the same Gaussian
  identity `TF:peak` design in every condition with `scale = FALSE`.
- Records RNA/ATAC layers, ATAC value semantics and a preprocessing fingerprint
  in candidate tables, the frozen dictionary, fitted networks and projections.
  Candidate unions and fixed fits reject mixed preprocessing references.
- Uses BH-adjusted condition coefficients for `penalty_effect`; unavailable or
  aliased coefficients remain unavailable rather than becoming fitted zeros.
- Removes obsolete condition APIs and compatibility arguments. The supported
  condition API is `discover_grn_edges()`, `union_grn_edges()`,
  `fit_grn_from_edges()`, `infer_condition_grn()`, `condition_grn_fit()`,
  `condition_grn_subgraph()`, `project_condition_grn_cells()` and
  `aggregate_condition_grn_projection()`.

# Pando 1.6.3

- Plans every target before any dense Gram, validation, Schur, or LLT
  allocation. High-dimensional targets use sparse matrix-free FISTA, exact
  residual validation, and a full-coordinate matrix-free Schur PCG refit; the
  active-support solve also becomes matrix-free when its dense factor would
  exceed the cumulative target-level worker budget. The default hybrid
  preconditioner uses an active-union block only when its retained and
  construction workspaces fit that same budget, otherwise using the exact same
  full-coordinate equations with a diagonal preconditioner.
- Calls the sparse solver directly inside C++ and removes the former
  Eigen-to-`dgCMatrix`-to-Eigen conversion. Condition-balanced transforms and
  training goodness-of-fit are now computed once by the native target engine.
- Adds strict `engine_control`, compact diagnostics, bounded target batches,
  resumable target checkpoints, allocation diagnostics, and a registered
  native numerical self-test. Failed iterative solves stop with numerical
  context and never return an unconverged approximation.
- Preserves candidate edges, double precision, equal-condition scaling,
  sparse-group objective, nested fold plans, per-lambda exact refits, lambda
  selection, shared-baseline equations, single-cell OOF projections, and the
  downstream condition-GRN coefficient contract. NativeSparseABI is now 6.

# Pando 1.6.1

- Adds an `exact_positions = FALSE` option to `find_motifs.GRNData()`. The
  default uses Signac's sparse motif presence/absence matrix without retaining
  genomic match positions, while `exact_positions = TRUE` preserves the
  original `AddMotifs()` path for footprinting workflows.

- Replaces repeated inner-fold sparse matrix prediction with exact validation
  sufficient statistics. For each condition and lambda, validation SSE is
  computed from `n`, column sums, `X'X`, `X'y`, response sum and response square
  sum, giving the same intercept-aware MSE without a per-cell lambda loop.
- Reuses the centered training Gram/RHS cache for both sparse-group selection and
  fixed-support direct-Schur refitting. A deterministic cost model chooses a
  centered-Gram FISTA path or the existing sparse matrix-free FISTA path.
- Keeps the public `infer_condition_grn()` API, fold plans, training-fold-only
  transforms, equal-condition weights, lambda rules, estimability masks,
  structural zeros, OOF projections and `ConditionGRNFit` schema unchanged.
- Raises the strict native condition ABI to 5 and stops on incompatible backend
  metadata, non-finite sufficient statistics, invalid SSE, failed factorization
  or incomplete OOF assignment. No R runtime fallback is introduced.

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
