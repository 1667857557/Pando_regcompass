# Pando 1.5.1

- Reads condition-model sparse inputs directly from canonical `dgCMatrix` slots instead of coercing S4 objects through `Rcpp::NumericMatrix` or `Rcpp::as<Eigen::MappedSparseMatrix<double>>`.
- Uses only the registered C++ TF-by-ATAC product kernel and C++ multitask solver. Missing symbols, unsupported sparse layouts, and native execution errors stop the analysis immediately.
- Initializes native FISTA from centered per-condition spectral-norm estimates rather than the summed Frobenius-norm upper bound.
- Adds installed-package CI for R 4.5.1, Matrix 1.7.4, native `.Call` registration, strict sparse type validation, the product kernel, and the multitask solver.
- Bumps the package version so installations from before and after the native ABI change are distinguishable.
