# Pando 1.5.2

- Replaces the canonical R direct-Schur support-constrained refit with a registered C++17/RcppEigen double-precision kernel; there is no runtime R fallback.
- Batches every inner-fold lambda refit in one native call, reusing the fold-level Gram/RHS cache and the estimability-dependent shared system across the complete lambda path.
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
