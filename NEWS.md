# Pando 1.5.1

- Reads condition-model sparse inputs directly from canonical `dgCMatrix` slots instead of coercing S4 objects through `Rcpp::NumericMatrix` or `Rcpp::as<Eigen::MappedSparseMatrix<double>>`.
- Adds an exact sparse R fallback for TF-by-ATAC products. Use `options(Pando.condition_product = "auto")` (default), `"cpp"`, or `"R"`.
- Retains the existing R reference fallback for the multitask solver.
- Initializes FISTA from centered per-condition spectral-norm estimates rather than the summed Frobenius-norm upper bound.
- Adds installed-package CI for R 4.5.1, Matrix 1.7.4, native `.Call` registration, the sparse product kernel, and the multitask solver.
- Bumps the package version so installations from before and after the native ABI change are distinguishable.
