#include <RcppEigen.h>
#include "condition_product_core.h"

using Eigen::SparseMatrix;

// [[Rcpp::export]]
SEXP condition_product_matrix_cpp(
    SEXP left_matrix,
    SEXP right_matrix,
    Rcpp::IntegerVector left_index,
    Rcpp::IntegerVector right_index
) {
    Eigen::MappedSparseMatrix<double> left_mapped =
        Rcpp::as<Eigen::MappedSparseMatrix<double>>(left_matrix);
    Eigen::MappedSparseMatrix<double> right_mapped =
        Rcpp::as<Eigen::MappedSparseMatrix<double>>(right_matrix);
    SparseMatrix<double> left(left_mapped);
    SparseMatrix<double> right(right_mapped);
    left.makeCompressed();
    right.makeCompressed();

    if (left.rows() != right.rows()) {
        Rcpp::stop("Input matrices must have identical row counts.");
    }
    if (left_index.size() != right_index.size()) {
        Rcpp::stop("Column index vectors must have equal lengths.");
    }
    SparseMatrix<double> answer = condition_product_core(
        left, right, left_index, right_index
    );
    return Rcpp::wrap(answer);
}
