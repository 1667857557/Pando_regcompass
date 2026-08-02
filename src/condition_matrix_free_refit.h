#ifndef PANDO_CONDITION_MATRIX_FREE_REFIT_H
#define PANDO_CONDITION_MATRIX_FREE_REFIT_H

#include "condition_execution_plan.h"
#include <RcppEigen.h>
#include <memory>
#include <vector>

class ConditionMatrixFreeRefitWorkspace {
public:
    ConditionMatrixFreeRefitWorkspace(
        const std::vector<Eigen::SparseMatrix<double>>& X,
        const std::vector<Eigen::VectorXd>& y,
        const Eigen::ArrayXXi& estimability_mask,
        double active_tol,
        const TargetEngineControl& control
    );
    ~ConditionMatrixFreeRefitWorkspace();
    ConditionMatrixFreeRefitWorkspace(
        const ConditionMatrixFreeRefitWorkspace&
    ) = delete;
    ConditionMatrixFreeRefitWorkspace& operator=(
        const ConditionMatrixFreeRefitWorkspace&
    ) = delete;
    ConditionMatrixFreeRefitWorkspace(
        ConditionMatrixFreeRefitWorkspace&&
    ) noexcept;
    ConditionMatrixFreeRefitWorkspace& operator=(
        ConditionMatrixFreeRefitWorkspace&&
    ) noexcept;

    Rcpp::List refit(
        const Eigen::MatrixXd& beta,
        double ridge,
        int path_index
    );

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

Rcpp::List condition_refit_path_matrix_free(
    const std::vector<Eigen::SparseMatrix<double>>& X,
    const std::vector<Eigen::VectorXd>& y,
    const Rcpp::List& beta_path,
    const Eigen::ArrayXXi& estimability_mask,
    const std::vector<double>& ridge,
    double active_tol,
    const TargetEngineControl& control
);

double condition_validation_mse_sparse(
    const Eigen::SparseMatrix<double>& X,
    const Eigen::VectorXd& y,
    const Eigen::VectorXd& beta,
    double intercept
);

Rcpp::List condition_matrix_free_refit_self_test();

#endif
