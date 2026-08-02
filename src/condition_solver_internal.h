#ifndef PANDO_CONDITION_SOLVER_INTERNAL_H
#define PANDO_CONDITION_SOLVER_INTERNAL_H

#include <RcppEigen.h>
#include <functional>
#include <string>
#include <vector>

using ConditionFitVisitor = std::function<void(
    int, const Rcpp::List&
)>;

Rcpp::List condition_fit_multitask_path_internal(
    const std::vector<Eigen::SparseMatrix<double>>& X,
    const std::vector<Eigen::VectorXd>& y,
    const std::vector<double>& lambda,
    double alpha,
    double condition_mix,
    const std::string& condition_weight,
    const Eigen::ArrayXXi& coefficient_mask,
    int max_iter,
    double tol_objective,
    double tol_coef,
    bool keep_history
);

Rcpp::List condition_fit_multitask_path_visit_internal(
    const std::vector<Eigen::SparseMatrix<double>>& X,
    const std::vector<Eigen::VectorXd>& y,
    const std::vector<double>& lambda,
    double alpha,
    double condition_mix,
    const std::string& condition_weight,
    const Eigen::ArrayXXi& coefficient_mask,
    int max_iter,
    double tol_objective,
    double tol_coef,
    bool keep_history,
    const ConditionFitVisitor& visitor,
    bool retain_fits
);

#endif
