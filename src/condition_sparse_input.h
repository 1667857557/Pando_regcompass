#ifndef PANDO_CONDITION_SPARSE_INPUT_H
#define PANDO_CONDITION_SPARSE_INPUT_H

#include <RcppEigen.h>
#include <cmath>
#include <limits>
#include <string>

inline Eigen::SparseMatrix<double> pando_condition_as_dgCMatrix(
    SEXP value, const std::string& label
) {
    if (!Rf_isS4(value) || !Rf_inherits(value, "dgCMatrix")) {
        Rcpp::stop(label + " must be a canonical Matrix::dgCMatrix.");
    }
    Rcpp::S4 matrix(value);
    Rcpp::IntegerVector dimensions = matrix.slot("Dim");
    Rcpp::IntegerVector column_pointer = matrix.slot("p");
    Rcpp::IntegerVector row_index = matrix.slot("i");
    Rcpp::NumericVector values = matrix.slot("x");
    if (dimensions.size() != 2 || dimensions[0] < 0 || dimensions[1] < 0) {
        Rcpp::stop(label + " has invalid dimensions.");
    }
    const int rows = dimensions[0];
    const int columns = dimensions[1];
    if (values.size() > std::numeric_limits<int>::max()) {
        Rcpp::stop(label + " contains too many nonzero entries for Eigen.");
    }
    const int nonzero = static_cast<int>(values.size());
    if (column_pointer.size() != columns + 1 ||
        row_index.size() != values.size() ||
        column_pointer[0] != 0 || column_pointer[columns] != nonzero) {
        Rcpp::stop(label + " has inconsistent dgCMatrix slots.");
    }
    for (int column = 0; column < columns; ++column) {
        if (column_pointer[column] > column_pointer[column + 1]) {
            Rcpp::stop(label + " has a non-monotone column pointer.");
        }
    }
    for (int index = 0; index < nonzero; ++index) {
        if (row_index[index] < 0 || row_index[index] >= rows ||
            !std::isfinite(values[index])) {
            Rcpp::stop(label + " contains an invalid sparse entry.");
        }
    }
    Eigen::MappedSparseMatrix<double> mapped(
        rows,
        columns,
        nonzero,
        column_pointer.begin(),
        row_index.begin(),
        values.begin()
    );
    Eigen::SparseMatrix<double> answer(mapped);
    answer.makeCompressed();
    return answer;
}

#endif
