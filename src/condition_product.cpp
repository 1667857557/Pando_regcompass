#include "condition_sparse_input.h"
#include <vector>

using Eigen::SparseMatrix;
using Eigen::Triplet;

// [[Rcpp::export]]
SEXP condition_product_matrix_cpp(
    SEXP left_matrix, SEXP right_matrix,
    Rcpp::IntegerVector left_index, Rcpp::IntegerVector right_index
) {
    SparseMatrix<double> left = pando_condition_as_dgCMatrix(
        left_matrix, "left_matrix"
    );
    SparseMatrix<double> right = pando_condition_as_dgCMatrix(
        right_matrix, "right_matrix"
    );
    const int n = left_index.size();
    if (left.rows() != right.rows() || right_index.size() != n) {
        Rcpp::stop("Input dimensions are not aligned.");
    }
    std::vector<Triplet<double>> entries;
    entries.reserve(static_cast<std::size_t>(right.nonZeros()));
    for (int e = 0; e < n; ++e) {
        const int a = left_index[e] - 1;
        const int b = right_index[e] - 1;
        if (a < 0 || a >= left.cols() || b < 0 || b >= right.cols()) {
            Rcpp::stop("Column index is out of range.");
        }
        SparseMatrix<double>::InnerIterator i(left, a);
        SparseMatrix<double>::InnerIterator j(right, b);
        while (i && j) {
            if (i.row() < j.row()) ++i;
            else if (j.row() < i.row()) ++j;
            else {
                const double value = i.value() * j.value();
                if (value != 0.0) entries.emplace_back(i.row(), e, value);
                ++i;
                ++j;
            }
        }
        if ((e & 255) == 0) Rcpp::checkUserInterrupt();
    }
    SparseMatrix<double> out(left.rows(), n);
    out.setFromTriplets(entries.begin(), entries.end());
    out.makeCompressed();
    return Rcpp::wrap(out);
}
