#include "condition_sparse_input.h"
#include <RcppEigen.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

using Eigen::ArrayXXi;
using Eigen::MatrixXd;
using Eigen::SparseMatrix;
using Eigen::Triplet;
using Eigen::VectorXd;

using ColSparse = SparseMatrix<double, Eigen::ColMajor, int>;
using RowSparse = SparseMatrix<double, Eigen::RowMajor, int>;

// Native kernels implemented in the package's existing translation units.
Rcpp::List condition_fit_multitask_path_cpp(
    Rcpp::List X_list,
    Rcpp::List y_list,
    Rcpp::NumericVector lambda,
    double alpha,
    double condition_mix,
    std::string condition_weight,
    Rcpp::LogicalMatrix coefficient_mask,
    int max_iter,
    double tol_objective,
    double tol_coef,
    bool keep_history
);

Rcpp::List condition_refit_path_cpp(
    Rcpp::List beta_path,
    Rcpp::LogicalMatrix estimability_mask,
    Rcpp::NumericVector ridge,
    double active_tol,
    Rcpp::List cache,
    double tolerance
);

namespace {

constexpr double kEps = std::numeric_limits<double>::epsilon();
const double kMetricJitter = std::sqrt(kEps);

struct RawTargetData {
    std::vector<RowSparse> X;
    std::vector<VectorXd> y;
    int predictors;
    int tasks;
};

struct Transform {
    std::vector<VectorXd> x_mean;
    std::vector<VectorXd> x_variance;
    VectorXd y_mean;
    VectorXd y_variance;
    VectorXd predictor_center;
    VectorXd predictor_scale;
    Eigen::ArrayXi predictor_estimable;
    double response_center;
    double response_scale;
    ArrayXXi estimability;
};

struct ScaledData {
    std::vector<ColSparse> X;
    std::vector<VectorXd> y;
};

struct PathBundle {
    bool intercept_only;
    std::vector<int> core_relative;
    ArrayXXi estimability_core;
    ScaledData training;
    Rcpp::List solver;
    Rcpp::List refits;
};

struct LambdaSelection {
    int selected_index;
    int minimum_index;
    int one_se_index;
    std::vector<double> cv_mean;
    std::vector<double> cv_se;
    std::vector<std::string> fold_backend;
    MatrixXd fold_loss;
    Rcpp::List fold_transform;
    int effective_folds;
    std::string reason;
};

std::vector<int> sequence(int size) {
    std::vector<int> out(size);
    std::iota(out.begin(), out.end(), 0);
    return out;
}

std::vector<int> rows_from_fold(
    const Rcpp::IntegerVector& fold,
    int value,
    bool equal
) {
    std::vector<int> rows;
    rows.reserve(fold.size());
    for (int i = 0; i < fold.size(); ++i) {
        if (Rcpp::IntegerVector::is_na(fold[i])) {
            Rcpp::stop("Fold assignments cannot contain NA.");
        }
        const bool match = fold[i] == value;
        if ((equal && match) || (!equal && !match)) rows.push_back(i);
    }
    return rows;
}

std::vector<int> map_positions(
    const std::vector<int>& base,
    const Rcpp::IntegerVector& fold,
    int value,
    bool equal
) {
    if (static_cast<int>(base.size()) != fold.size()) {
        Rcpp::stop("Nested fold assignments are not aligned to training rows.");
    }
    std::vector<int> rows;
    rows.reserve(base.size());
    for (int i = 0; i < fold.size(); ++i) {
        if (Rcpp::IntegerVector::is_na(fold[i])) {
            Rcpp::stop("Fold assignments cannot contain NA.");
        }
        const bool match = fold[i] == value;
        if ((equal && match) || (!equal && !match)) rows.push_back(base[i]);
    }
    return rows;
}

RawTargetData parse_raw_target(
    const Rcpp::List& X_list,
    const Rcpp::List& y_list
) {
    if (X_list.size() < 2 || X_list.size() != y_list.size()) {
        Rcpp::stop(
            "Target engine requires aligned predictor and response lists for "
            "at least two conditions."
        );
    }
    RawTargetData out;
    out.tasks = X_list.size();
    out.predictors = -1;
    out.X.reserve(out.tasks);
    out.y.reserve(out.tasks);
    for (int task = 0; task < out.tasks; ++task) {
        SEXP value = X_list[task];
        if (!Rf_isS4(value)) {
            Rcpp::stop(
                "Target engine requires canonical dgCMatrix predictors."
            );
        }
        ColSparse column = pando_condition_as_dgCMatrix(
            value, "target engine predictor matrix"
        );
        column.makeCompressed();
        RowSparse row = column;
        row.makeCompressed();
        Rcpp::NumericVector y_input(y_list[task]);
        Eigen::Map<VectorXd> y_map(y_input.begin(), y_input.size());
        VectorXd y(y_map);
        if (!y.allFinite()) {
            Rcpp::stop("Target engine response contains non-finite values.");
        }
        if (row.rows() != y.size() || row.rows() < 2) {
            Rcpp::stop(
                "Each target condition requires aligned predictors and at "
                "least two responses."
            );
        }
        if (out.predictors < 0) out.predictors = row.cols();
        if (row.cols() != out.predictors) {
            Rcpp::stop("Target condition matrices do not share predictors.");
        }
        out.X.emplace_back(std::move(row));
        out.y.emplace_back(std::move(y));
    }
    if (out.predictors < 1) {
        Rcpp::stop("Target engine requires at least one predictor.");
    }
    return out;
}

ArrayXXi parse_mask(
    const Rcpp::LogicalMatrix& input,
    int predictors,
    int tasks
) {
    if (input.nrow() != predictors || input.ncol() != tasks) {
        Rcpp::stop("coefficient_mask is not aligned to target predictors.");
    }
    ArrayXXi mask(predictors, tasks);
    for (int row = 0; row < predictors; ++row) {
        bool any = false;
        for (int task = 0; task < tasks; ++task) {
            const int value = input(row, task);
            if (value == NA_LOGICAL) {
                Rcpp::stop("coefficient_mask cannot contain NA.");
            }
            mask(row, task) = value ? 1 : 0;
            any = any || value;
        }
        if (!any) {
            Rcpp::stop(
                "Every target predictor requires at least one eligible condition."
            );
        }
    }
    return mask;
}

std::vector<std::vector<int>> all_rows(const RawTargetData& data) {
    std::vector<std::vector<int>> rows(data.tasks);
    for (int task = 0; task < data.tasks; ++task) {
        rows[task] = sequence(data.X[task].rows());
    }
    return rows;
}

ArrayXXi raw_estimability_mask(
    const RawTargetData& data,
    const std::vector<std::vector<int>>& rows,
    const std::vector<int>& columns,
    const ArrayXXi& mask
) {
    const int p = columns.size();
    if (mask.rows() != p || mask.cols() != data.tasks ||
        static_cast<int>(rows.size()) != data.tasks) {
        Rcpp::stop("Raw estimability inputs are not aligned.");
    }
    std::vector<int> column_map(data.predictors, -1);
    for (int index = 0; index < p; ++index) {
        column_map[columns[index]] = index;
    }
    ArrayXXi out = ArrayXXi::Zero(p, data.tasks);
    for (int task = 0; task < data.tasks; ++task) {
        const int n = rows[task].size();
        if (n < 1) Rcpp::stop("Raw estimability requires condition cells.");
        VectorXd sum = VectorXd::Zero(p);
        VectorXd square = VectorXd::Zero(p);
        for (int source_row : rows[task]) {
            for (RowSparse::InnerIterator it(data.X[task], source_row); it; ++it) {
                const int mapped = column_map[it.col()];
                if (mapped < 0) continue;
                sum[mapped] += it.value();
                square[mapped] += it.value() * it.value();
            }
        }
        VectorXd mean = sum / static_cast<double>(n);
        VectorXd variance = square / static_cast<double>(n) -
            mean.array().square().matrix();
        variance = variance.cwiseMax(0.0);
        for (int predictor = 0; predictor < p; ++predictor) {
            out(predictor, task) = mask(predictor, task) &&
                std::isfinite(variance[predictor]) &&
                variance[predictor] > kEps;
        }
    }
    return out;
}

Transform compute_transform(
    const RawTargetData& data,
    const std::vector<std::vector<int>>& rows,
    const std::vector<int>& columns,
    const ArrayXXi& mask
) {
    const int p = columns.size();
    const int tasks = data.tasks;
    if (mask.rows() != p || mask.cols() != tasks ||
        static_cast<int>(rows.size()) != tasks) {
        Rcpp::stop("Transform inputs are not aligned.");
    }
    Transform out;
    out.x_mean.resize(tasks);
    out.x_variance.resize(tasks);
    out.y_mean = VectorXd::Zero(tasks);
    out.y_variance = VectorXd::Zero(tasks);
    out.predictor_center = VectorXd::Zero(p);
    out.predictor_scale = VectorXd::Zero(p);
    out.predictor_estimable = Eigen::ArrayXi::Zero(p);
    out.estimability = ArrayXXi::Zero(p, tasks);
    out.response_center = 0.0;
    out.response_scale = 0.0;

    std::vector<int> column_map(data.predictors, -1);
    for (int index = 0; index < p; ++index) {
        const int column = columns[index];
        if (column < 0 || column >= data.predictors ||
            column_map[column] >= 0) {
            Rcpp::stop("Transform predictor map is invalid.");
        }
        column_map[column] = index;
    }

    for (int task = 0; task < tasks; ++task) {
        const int n = rows[task].size();
        if (n < 2) {
            Rcpp::stop(
                "Every transform training condition requires at least two cells."
            );
        }
        VectorXd sum = VectorXd::Zero(p);
        VectorXd square = VectorXd::Zero(p);
        double y_sum = 0.0;
        for (int source_row : rows[task]) {
            if (source_row < 0 || source_row >= data.X[task].rows()) {
                Rcpp::stop("Transform row index is out of bounds.");
            }
            const double response = data.y[task][source_row];
            y_sum += response;
            for (RowSparse::InnerIterator it(data.X[task], source_row); it; ++it) {
                const int mapped = column_map[it.col()];
                if (mapped < 0) continue;
                const double value = it.value();
                sum[mapped] += value;
                square[mapped] += value * value;
            }
        }
        VectorXd mean = sum / static_cast<double>(n);
        VectorXd variance = square / static_cast<double>(n) -
            mean.array().square().matrix();
        variance = variance.cwiseMax(0.0);
        const double y_mean = y_sum / static_cast<double>(n);
        double y_centered_square = 0.0;
        for (int source_row : rows[task]) {
            const double difference = data.y[task][source_row] - y_mean;
            y_centered_square += difference * difference;
        }
        const double y_variance = y_centered_square /
            static_cast<double>(n);
        out.x_mean[task] = mean;
        out.x_variance[task] = variance;
        out.y_mean[task] = y_mean;
        out.y_variance[task] = y_variance;
        out.predictor_center.noalias() +=
            mean / static_cast<double>(tasks);
        out.predictor_scale.noalias() +=
            variance / static_cast<double>(tasks);
        out.response_center += y_mean / static_cast<double>(tasks);
        out.response_scale += y_variance / static_cast<double>(tasks);
        for (int predictor = 0; predictor < p; ++predictor) {
            out.estimability(predictor, task) =
                mask(predictor, task) &&
                std::isfinite(variance[predictor]) &&
                variance[predictor] > kEps;
        }
    }
    out.predictor_scale = out.predictor_scale.array().sqrt().matrix();
    out.response_scale = std::sqrt(out.response_scale);
    if (!std::isfinite(out.response_scale) || out.response_scale <= kEps) {
        Rcpp::stop(
            "The equal-condition within-condition response variance is zero."
        );
    }
    for (int predictor = 0; predictor < p; ++predictor) {
        const bool estimable = std::isfinite(out.predictor_scale[predictor]) &&
            out.predictor_scale[predictor] > kEps;
        out.predictor_estimable[predictor] = estimable ? 1 : 0;
        if (!estimable) out.predictor_scale[predictor] = 1.0;
    }
    return out;
}

std::vector<int> transform_keep(const Transform& transform) {
    std::vector<int> keep;
    keep.reserve(transform.estimability.rows());
    for (int predictor = 0; predictor < transform.estimability.rows(); ++predictor) {
        bool any = false;
        for (int task = 0; task < transform.estimability.cols(); ++task) {
            any = any || transform.estimability(predictor, task);
        }
        if (any && transform.predictor_estimable[predictor]) {
            keep.push_back(predictor);
        }
    }
    return keep;
}

std::vector<int> map_columns(
    const std::vector<int>& base,
    const std::vector<int>& relative
) {
    std::vector<int> out;
    out.reserve(relative.size());
    for (int index : relative) {
        if (index < 0 || index >= static_cast<int>(base.size())) {
            Rcpp::stop("Relative predictor index is out of bounds.");
        }
        out.push_back(base[index]);
    }
    return out;
}

ArrayXXi subset_mask_rows(
    const ArrayXXi& mask,
    const std::vector<int>& rows
) {
    ArrayXXi out(rows.size(), mask.cols());
    for (int row = 0; row < static_cast<int>(rows.size()); ++row) {
        out.row(row) = mask.row(rows[row]);
    }
    return out;
}

VectorXd subset_vector(const VectorXd& input, const std::vector<int>& index) {
    VectorXd out(index.size());
    for (int i = 0; i < static_cast<int>(index.size()); ++i) {
        out[i] = input[index[i]];
    }
    return out;
}

ColSparse build_scaled_matrix(
    const RowSparse& source,
    const std::vector<int>& rows,
    const std::vector<int>& global_columns,
    const VectorXd& scale
) {
    if (global_columns.size() != static_cast<std::size_t>(scale.size())) {
        Rcpp::stop("Scaled matrix columns and scales are not aligned.");
    }
    std::vector<int> column_map(source.cols(), -1);
    for (int column = 0; column < static_cast<int>(global_columns.size()); ++column) {
        const int global = global_columns[column];
        if (global < 0 || global >= source.cols()) {
            Rcpp::stop("Scaled matrix column is out of bounds.");
        }
        column_map[global] = column;
    }
    std::vector<Triplet<double>> triplets;
    triplets.reserve(std::min<std::size_t>(
        source.nonZeros(),
        static_cast<std::size_t>(rows.size()) * global_columns.size()
    ));
    for (int output_row = 0; output_row < static_cast<int>(rows.size()); ++output_row) {
        const int source_row = rows[output_row];
        if (source_row < 0 || source_row >= source.rows()) {
            Rcpp::stop("Scaled matrix row is out of bounds.");
        }
        for (RowSparse::InnerIterator it(source, source_row); it; ++it) {
            const int output_column = column_map[it.col()];
            if (output_column < 0) continue;
            const double value = it.value() / scale[output_column];
            if (value != 0.0) {
                triplets.emplace_back(output_row, output_column, value);
            }
        }
    }
    ColSparse out(rows.size(), global_columns.size());
    out.setFromTriplets(triplets.begin(), triplets.end());
    out.makeCompressed();
    return out;
}

VectorXd build_scaled_response(
    const VectorXd& source,
    const std::vector<int>& rows,
    double center,
    double scale
) {
    VectorXd out(rows.size());
    for (int i = 0; i < static_cast<int>(rows.size()); ++i) {
        out[i] = (source[rows[i]] - center) / scale;
    }
    return out;
}

ScaledData build_scaled_data(
    const RawTargetData& data,
    const std::vector<std::vector<int>>& rows,
    const std::vector<int>& global_columns,
    const VectorXd& scale,
    double response_center,
    double response_scale
) {
    ScaledData out;
    out.X.reserve(data.tasks);
    out.y.reserve(data.tasks);
    for (int task = 0; task < data.tasks; ++task) {
        out.X.emplace_back(build_scaled_matrix(
            data.X[task], rows[task], global_columns, scale
        ));
        out.y.emplace_back(build_scaled_response(
            data.y[task], rows[task], response_center, response_scale
        ));
    }
    return out;
}

VectorXd population_variance(const ColSparse& matrix) {
    VectorXd sum = VectorXd::Zero(matrix.cols());
    VectorXd square = VectorXd::Zero(matrix.cols());
    for (int column = 0; column < matrix.outerSize(); ++column) {
        for (ColSparse::InnerIterator it(matrix, column); it; ++it) {
            sum[column] += it.value();
            square[column] += it.value() * it.value();
        }
    }
    const double n = static_cast<double>(matrix.rows());
    VectorXd mean = sum / n;
    return (square / n - mean.array().square().matrix()).cwiseMax(0.0);
}

ArrayXXi true_variance_mask(
    const std::vector<ColSparse>& X,
    const ArrayXXi& mask
) {
    ArrayXXi out = mask;
    for (int task = 0; task < static_cast<int>(X.size()); ++task) {
        VectorXd variance = population_variance(X[task]);
        for (int predictor = 0; predictor < variance.size(); ++predictor) {
            out(predictor, task) = out(predictor, task) &&
                std::isfinite(variance[predictor]) &&
                variance[predictor] > kEps;
        }
    }
    return out;
}

std::vector<int> any_estimable(const ArrayXXi& mask) {
    std::vector<int> keep;
    keep.reserve(mask.rows());
    for (int predictor = 0; predictor < mask.rows(); ++predictor) {
        bool any = false;
        for (int task = 0; task < mask.cols(); ++task) {
            any = any || mask(predictor, task);
        }
        if (any) keep.push_back(predictor);
    }
    return keep;
}

ColSparse subset_columns(
    const ColSparse& source,
    const std::vector<int>& columns
) {
    std::vector<Triplet<double>> triplets;
    triplets.reserve(source.nonZeros());
    for (int output_column = 0;
         output_column < static_cast<int>(columns.size()); ++output_column) {
        const int source_column = columns[output_column];
        for (ColSparse::InnerIterator it(source, source_column); it; ++it) {
            triplets.emplace_back(it.row(), output_column, it.value());
        }
    }
    ColSparse out(source.rows(), columns.size());
    out.setFromTriplets(triplets.begin(), triplets.end());
    out.makeCompressed();
    return out;
}

Rcpp::S4 to_dgCMatrix(const ColSparse& input) {
    ColSparse matrix = input;
    matrix.makeCompressed();
    Rcpp::S4 out("dgCMatrix");
    if (matrix.nonZeros() > 0) {
        out.slot("i") = Rcpp::IntegerVector(
            matrix.innerIndexPtr(),
            matrix.innerIndexPtr() + matrix.nonZeros()
        );
        out.slot("x") = Rcpp::NumericVector(
            matrix.valuePtr(),
            matrix.valuePtr() + matrix.nonZeros()
        );
    } else {
        out.slot("i") = Rcpp::IntegerVector(0);
        out.slot("x") = Rcpp::NumericVector(0);
    }
    out.slot("p") = Rcpp::IntegerVector(
        matrix.outerIndexPtr(),
        matrix.outerIndexPtr() + matrix.outerSize() + 1
    );
    out.slot("Dim") = Rcpp::IntegerVector::create(
        matrix.rows(), matrix.cols()
    );
    out.slot("Dimnames") = Rcpp::List::create(R_NilValue, R_NilValue);
    return out;
}

Rcpp::List to_matrix_list(const std::vector<ColSparse>& matrices) {
    Rcpp::List out(matrices.size());
    for (int task = 0; task < static_cast<int>(matrices.size()); ++task) {
        out[task] = to_dgCMatrix(matrices[task]);
    }
    return out;
}

Rcpp::List to_vector_list(const std::vector<VectorXd>& vectors) {
    Rcpp::List out(vectors.size());
    for (int task = 0; task < static_cast<int>(vectors.size()); ++task) {
        out[task] = Rcpp::wrap(vectors[task]);
    }
    return out;
}

Rcpp::LogicalMatrix to_logical_matrix(const ArrayXXi& mask) {
    Rcpp::LogicalMatrix out(mask.rows(), mask.cols());
    for (int row = 0; row < mask.rows(); ++row) {
        for (int column = 0; column < mask.cols(); ++column) {
            out(row, column) = mask(row, column) != 0;
        }
    }
    return out;
}

Rcpp::NumericVector to_numeric(const std::vector<double>& values) {
    return Rcpp::NumericVector(values.begin(), values.end());
}

std::vector<double> sorted_lambda(Rcpp::Nullable<Rcpp::NumericVector> input) {
    if (input.isNull()) return {};
    Rcpp::NumericVector value(input.get());
    std::vector<double> lambda(value.begin(), value.end());
    if (lambda.empty()) return {};
    for (double current : lambda) {
        if (!std::isfinite(current) || current < 0.0) {
            Rcpp::stop("lambda must contain finite non-negative values.");
        }
    }
    std::sort(lambda.begin(), lambda.end(), std::greater<double>());
    lambda.erase(std::unique(lambda.begin(), lambda.end()), lambda.end());
    return lambda;
}

double lambda_max(
    const ScaledData& data,
    const ArrayXXi& mask,
    double alpha,
    double condition_mix
) {
    const int p = data.X[0].cols();
    const int tasks = data.X.size();
    MatrixXd gradient = MatrixXd::Zero(p, tasks);
    for (int task = 0; task < tasks; ++task) {
        const double y_mean = data.y[task].mean();
        VectorXd residual = data.y[task].array() - y_mean;
        gradient.col(task).noalias() =
            -(data.X[task].transpose() * residual) /
            static_cast<double>(data.X[task].rows());
        for (int predictor = 0; predictor < p; ++predictor) {
            if (!mask(predictor, task)) gradient(predictor, task) = 0.0;
        }
    }
    double value = 0.0;
    if (alpha <= kEps) {
        value = gradient.cwiseAbs().maxCoeff();
    } else if (condition_mix <= kEps) {
        for (int predictor = 0; predictor < p; ++predictor) {
            value = std::max(value, gradient.row(predictor).norm() / alpha);
        }
    } else {
        value = gradient.cwiseAbs().maxCoeff() / (alpha * condition_mix);
    }
    if (!std::isfinite(value) || value <= 0.0) value = 1.0;
    return value * 1.001;
}

std::vector<double> make_lambda_path(
    const ScaledData& data,
    const ArrayXXi& mask,
    double alpha,
    double condition_mix,
    int nlambda,
    double lambda_min_ratio
) {
    if (nlambda < 2) {
        Rcpp::stop("nlambda must be at least two when lambda is automatic.");
    }
    if (!std::isfinite(lambda_min_ratio)) {
        int total_cells = 0;
        for (const auto& matrix : data.X) total_cells += matrix.rows();
        lambda_min_ratio = total_cells < data.X[0].cols() ? 0.01 : 1e-4;
    }
    if (lambda_min_ratio <= 0.0 || lambda_min_ratio >= 1.0) {
        Rcpp::stop("lambda_min_ratio must be between zero and one.");
    }
    const double maximum = lambda_max(
        data, mask, alpha, condition_mix
    );
    std::vector<double> out(nlambda);
    const double start = std::log(maximum);
    const double end = std::log(maximum * lambda_min_ratio);
    for (int index = 0; index < nlambda; ++index) {
        const double fraction = static_cast<double>(index) /
            static_cast<double>(nlambda - 1);
        out[index] = std::exp(start + fraction * (end - start));
    }
    return out;
}

Rcpp::List make_refit_cache(const ScaledData& data) {
    const int tasks = data.X.size();
    const int p = data.X[0].cols();
    Rcpp::List task_list(tasks);
    MatrixXd common = MatrixXd::Zero(p, p);
    Rcpp::NumericVector loss(tasks);
    Rcpp::NumericVector average(tasks, 1.0 / static_cast<double>(tasks));
    for (int task = 0; task < tasks; ++task) {
        const double n = static_cast<double>(data.X[task].rows());
        VectorXd x_sum = VectorXd::Zero(p);
        for (int column = 0; column < data.X[task].outerSize(); ++column) {
            for (ColSparse::InnerIterator it(data.X[task], column); it; ++it) {
                x_sum[column] += it.value();
            }
        }
        VectorXd x_mean = x_sum / n;
        const double y_mean = data.y[task].mean();
        MatrixXd gram = MatrixXd(data.X[task].transpose() * data.X[task]);
        gram.noalias() -= n * (x_mean * x_mean.transpose());
        for (int row = 0; row < p; ++row) {
            for (int column = 0; column < p; ++column) {
                if (std::abs(gram(row, column)) < kEps) {
                    gram(row, column) = 0.0;
                }
            }
        }
        VectorXd rhs = data.X[task].transpose() * data.y[task] -
            n * x_mean * y_mean;
        loss[task] = 1.0 / n;
        common.noalias() += loss[task] * gram;
        task_list[task] = Rcpp::List::create(
            Rcpp::Named("x_mean") = x_mean,
            Rcpp::Named("y_mean") = y_mean,
            Rcpp::Named("gram") = gram,
            Rcpp::Named("rhs") = rhs
        );
    }
    common.diagonal().array() += kMetricJitter;
    return Rcpp::List::create(
        Rcpp::Named("task") = task_list,
        Rcpp::Named("common_metric") = common,
        Rcpp::Named("loss_weights") = loss,
        Rcpp::Named("average_weights") = average
    );
}


struct GramTaskProblem {
    MatrixXd gram;
    VectorXd rhs;
    VectorXd x_mean;
    double y_mean;
    double response_ss;
};

struct GramProblem {
    std::vector<GramTaskProblem> task;
    VectorXd loss_weights;
    int predictors;
    int tasks;
};

struct GramSmoothResult {
    double value;
    MatrixXd gradient;
    VectorXd intercept;
};

struct ValidationTaskStats {
    int n;
    VectorXd sum_x;
    MatrixXd xtx;
    VectorXd xty;
    double sum_y;
    double sum_y2;
};

double gram_sparse_group_penalty(
    const MatrixXd& B,
    double lambda,
    double alpha,
    double condition_mix
) {
    const double strength = lambda * alpha;
    double group = 0.0;
    for (int row = 0; row < B.rows(); ++row) {
        group += B.row(row).norm();
    }
    return strength * (
        (1.0 - condition_mix) * group +
        condition_mix * B.cwiseAbs().sum()
    );
}

MatrixXd gram_sparse_group_prox(
    const MatrixXd& value,
    const ArrayXXi& mask,
    double step,
    double lambda,
    double alpha,
    double condition_mix
) {
    MatrixXd out = MatrixXd::Zero(value.rows(), value.cols());
    const double element_threshold =
        step * lambda * alpha * condition_mix;
    const double group_threshold =
        step * lambda * alpha * (1.0 - condition_mix);
    for (int row = 0; row < value.rows(); ++row) {
        for (int task = 0; task < value.cols(); ++task) {
            if (!mask(row, task)) continue;
            const double current = value(row, task);
            const double magnitude = std::max(
                std::abs(current) - element_threshold, 0.0
            );
            out(row, task) = std::copysign(magnitude, current);
        }
        const double norm = out.row(row).norm();
        if (norm > 0.0) {
            out.row(row) *= std::max(
                1.0 - group_threshold / norm, 0.0
            );
        }
    }
    return out;
}

GramProblem parse_gram_problem(
    const Rcpp::List& cache,
    const std::vector<VectorXd>& response
) {
    if (!cache.containsElementNamed("task") ||
        !cache.containsElementNamed("loss_weights")) {
        Rcpp::stop("Gram solver received an incomplete sufficient-statistic cache.");
    }
    Rcpp::List task_input(cache["task"]);
    Rcpp::NumericVector weight_input(cache["loss_weights"]);
    if (task_input.size() < 2 ||
        task_input.size() != weight_input.size() ||
        task_input.size() != static_cast<int>(response.size())) {
        Rcpp::stop("Gram solver cache is not condition-aligned.");
    }
    GramProblem out;
    out.tasks = task_input.size();
    out.predictors = -1;
    out.loss_weights = VectorXd(out.tasks);
    out.task.reserve(out.tasks);
    for (int task = 0; task < out.tasks; ++task) {
        Rcpp::List current(task_input[task]);
        GramTaskProblem value;
        value.gram = Rcpp::as<MatrixXd>(current["gram"]);
        value.rhs = Rcpp::as<VectorXd>(current["rhs"]);
        value.x_mean = Rcpp::as<VectorXd>(current["x_mean"]);
        value.y_mean = Rcpp::as<double>(current["y_mean"]);
        if (out.predictors < 0) out.predictors = value.gram.rows();
        if (value.gram.rows() != out.predictors ||
            value.gram.cols() != out.predictors ||
            value.rhs.size() != out.predictors ||
            value.x_mean.size() != out.predictors ||
            response[task].size() < 1) {
            Rcpp::stop("Gram solver cache dimensions are invalid.");
        }
        value.gram = (
            0.5 * (value.gram + value.gram.transpose())
        ).eval();
        value.response_ss = (
            response[task].array() - value.y_mean
        ).square().sum();
        out.loss_weights[task] = weight_input[task];
        if (!value.gram.allFinite() || !value.rhs.allFinite() ||
            !value.x_mean.allFinite() || !std::isfinite(value.y_mean) ||
            !std::isfinite(value.response_ss) ||
            !std::isfinite(out.loss_weights[task]) ||
            out.loss_weights[task] <= 0.0) {
            Rcpp::stop("Gram solver cache contains non-finite values.");
        }
        out.task.emplace_back(std::move(value));
    }
    return out;
}

GramSmoothResult gram_profiled_smooth(
    const MatrixXd& B,
    const GramProblem& problem,
    double ridge
) {
    GramSmoothResult out{
        0.0,
        MatrixXd::Zero(problem.predictors, problem.tasks),
        VectorXd::Zero(problem.tasks)
    };
    for (int task = 0; task < problem.tasks; ++task) {
        const GramTaskProblem& value = problem.task[task];
        VectorXd gram_beta = value.gram * B.col(task);
        const double quadratic =
            B.col(task).dot(gram_beta) -
            2.0 * value.rhs.dot(B.col(task)) +
            value.response_ss;
        out.value += 0.5 * problem.loss_weights[task] * quadratic;
        out.gradient.col(task).noalias() =
            problem.loss_weights[task] * (gram_beta - value.rhs);
        out.intercept[task] =
            value.y_mean - value.x_mean.dot(B.col(task));
    }
    if (ridge > 0.0) {
        out.value += 0.5 * ridge * B.squaredNorm();
        out.gradient.noalias() += ridge * B;
    }
    if (!std::isfinite(out.value) || !out.gradient.allFinite() ||
        !out.intercept.allFinite()) {
        Rcpp::stop("Gram solver produced non-finite smooth quantities.");
    }
    return out;
}

double gram_spectral_squared_norm(const MatrixXd& gram) {
    if (gram.rows() == 0) return 0.0;
    VectorXd direction(gram.cols());
    for (int column = 0; column < gram.cols(); ++column) {
        direction[column] = 1.0 + static_cast<double>(column % 7) / 7.0;
    }
    direction.normalize();
    double previous = -1.0;
    for (int iteration = 0; iteration < 30; ++iteration) {
        VectorXd next = gram * direction;
        const double estimate = direction.dot(next);
        const double next_norm = next.norm();
        if (!std::isfinite(estimate) || !std::isfinite(next_norm)) {
            Rcpp::stop("Gram spectral-norm iteration became non-finite.");
        }
        if (next_norm <= kEps) return std::max(estimate, 0.0);
        direction = next / next_norm;
        if (previous >= 0.0 &&
            std::abs(estimate - previous) <=
                1e-8 * std::max(1.0, estimate)) {
            break;
        }
        previous = estimate;
    }
    return std::max(direction.dot(gram * direction), 0.0);
}

double gram_initial_step(
    const VectorXd& spectral,
    const VectorXd& weights,
    double ridge
) {
    double lipschitz = ridge;
    for (int task = 0; task < spectral.size(); ++task) {
        lipschitz = std::max(
            lipschitz,
            ridge + 1.05 * weights[task] * spectral[task]
        );
    }
    return (!std::isfinite(lipschitz) || lipschitz <= 0.0) ?
        1.0 : 1.0 / lipschitz;
}

Rcpp::List fit_one_lambda_gram(
    const GramProblem& problem,
    const ArrayXXi& mask,
    double lambda,
    double alpha,
    double condition_mix,
    const MatrixXd& initial_B,
    double initial_step,
    double default_step,
    int max_iter,
    double tol_objective,
    double tol_coef
) {
    const double ridge = lambda * (1.0 - alpha);
    MatrixXd B = initial_B;
    for (int row = 0; row < B.rows(); ++row) {
        for (int task = 0; task < B.cols(); ++task) {
            if (!mask(row, task)) B(row, task) = 0.0;
        }
    }
    MatrixXd Z = B;
    double acceleration = 1.0;
    double step = default_step;
    if (std::isfinite(initial_step) && initial_step > 0.0) {
        step = std::max(initial_step, default_step);
    }
    GramSmoothResult smooth_B = gram_profiled_smooth(B, problem, ridge);
    double objective_previous = smooth_B.value + gram_sparse_group_penalty(
        B, lambda, alpha, condition_mix
    );
    bool converged = false;
    double coef_change = std::numeric_limits<double>::infinity();
    double objective_change = std::numeric_limits<double>::infinity();
    int iteration = 0;
    const double backtrack = 0.5;
    const double min_step = 1e-14;
    for (iteration = 1; iteration <= max_iter; ++iteration) {
        GramSmoothResult smooth_Z = Z.isApprox(B, 0.0) ?
            smooth_B : gram_profiled_smooth(Z, problem, ridge);
        MatrixXd candidate;
        GramSmoothResult smooth_candidate;
        while (true) {
            candidate = gram_sparse_group_prox(
                Z - step * smooth_Z.gradient,
                mask,
                step,
                lambda,
                alpha,
                condition_mix
            );
            smooth_candidate = gram_profiled_smooth(
                candidate, problem, ridge
            );
            MatrixXd difference = candidate - Z;
            const double bound = smooth_Z.value +
                (smooth_Z.gradient.array() * difference.array()).sum() +
                difference.squaredNorm() / (2.0 * step);
            if (smooth_candidate.value <= bound + 1e-10) break;
            step *= backtrack;
            if (step < min_step) {
                Rcpp::stop("Gram backtracking line search reached min_step.");
            }
        }
        double objective_candidate = smooth_candidate.value +
            gram_sparse_group_penalty(
                candidate, lambda, alpha, condition_mix
            );
        if (objective_candidate > objective_previous + 1e-10) {
            acceleration = 1.0;
            Z = B;
            smooth_Z = smooth_B;
            while (true) {
                candidate = gram_sparse_group_prox(
                    Z - step * smooth_Z.gradient,
                    mask,
                    step,
                    lambda,
                    alpha,
                    condition_mix
                );
                smooth_candidate = gram_profiled_smooth(
                    candidate, problem, ridge
                );
                MatrixXd difference = candidate - Z;
                const double bound = smooth_Z.value +
                    (smooth_Z.gradient.array() * difference.array()).sum() +
                    difference.squaredNorm() / (2.0 * step);
                if (smooth_candidate.value <= bound + 1e-10) break;
                step *= backtrack;
                if (step < min_step) {
                    Rcpp::stop(
                        "Gram backtracking reached min_step after restart."
                    );
                }
            }
            objective_candidate = smooth_candidate.value +
                gram_sparse_group_penalty(
                    candidate, lambda, alpha, condition_mix
                );
        }
        coef_change = (candidate - B).norm() /
            (B.norm() + kEps);
        objective_change = std::abs(
            objective_candidate - objective_previous
        ) / (std::abs(objective_previous) + kEps);
        MatrixXd B_previous = B;
        B = candidate;
        smooth_B = smooth_candidate;
        objective_previous = objective_candidate;
        if (objective_change < tol_objective && coef_change < tol_coef) {
            converged = true;
            break;
        }
        const double acceleration_new = (
            1.0 + std::sqrt(1.0 + 4.0 * acceleration * acceleration)
        ) / 2.0;
        MatrixXd Z_new = B +
            ((acceleration - 1.0) / acceleration_new) *
            (B - B_previous);
        const double restart =
            ((Z - B).array() * (B - B_previous).array()).sum();
        if (restart > 0.0) {
            acceleration = 1.0;
            Z = B;
        } else {
            acceleration = acceleration_new;
            Z = Z_new;
        }
    }
    if (iteration > max_iter) iteration = max_iter;
    return Rcpp::List::create(
        Rcpp::Named("beta") = B,
        Rcpp::Named("intercept") = smooth_B.intercept,
        Rcpp::Named("lambda") = lambda,
        Rcpp::Named("objective") = objective_previous,
        Rcpp::Named("objective_change") = objective_change,
        Rcpp::Named("coef_change") = coef_change,
        Rcpp::Named("iterations") = iteration,
        Rcpp::Named("converged") = converged,
        Rcpp::Named("step") = step,
        Rcpp::Named("history") = R_NilValue,
        Rcpp::Named("backend") = "cpp_eigen_centered_gram_fista"
    );
}

Rcpp::List fit_multitask_path_gram(
    const ScaledData& training,
    const Rcpp::List& cache,
    const ArrayXXi& mask,
    const std::vector<double>& lambda,
    double alpha,
    double condition_mix,
    int max_iter,
    double tol_objective,
    double tol_coef
) {
    GramProblem problem = parse_gram_problem(cache, training.y);
    if (mask.rows() != problem.predictors || mask.cols() != problem.tasks) {
        Rcpp::stop("Gram solver coefficient mask is not aligned.");
    }
    VectorXd spectral(problem.tasks);
    for (int task = 0; task < problem.tasks; ++task) {
        spectral[task] = gram_spectral_squared_norm(
            problem.task[task].gram
        );
    }
    MatrixXd B = MatrixXd::Zero(problem.predictors, problem.tasks);
    double step = std::numeric_limits<double>::quiet_NaN();
    Rcpp::List fits(lambda.size());
    for (int index = 0; index < static_cast<int>(lambda.size()); ++index) {
        const double current = lambda[index];
        if (!std::isfinite(current) || current < 0.0) {
            Rcpp::stop("Gram solver lambda is invalid.");
        }
        const double default_step = gram_initial_step(
            spectral,
            problem.loss_weights,
            current * (1.0 - alpha)
        );
        Rcpp::List fit = fit_one_lambda_gram(
            problem,
            mask,
            current,
            alpha,
            condition_mix,
            B,
            step,
            default_step,
            max_iter,
            tol_objective,
            tol_coef
        );
        B = Rcpp::as<MatrixXd>(fit["beta"]);
        step = Rcpp::as<double>(fit["step"]);
        fits[index] = fit;
    }
    return Rcpp::List::create(
        Rcpp::Named("lambda") = to_numeric(lambda),
        Rcpp::Named("fits") = fits,
        Rcpp::Named("backend") = "cpp_eigen_centered_gram_fista"
    );
}

bool use_gram_solver(const ScaledData& training) {
    const double p = static_cast<double>(training.X[0].cols());
    const double tasks = static_cast<double>(training.X.size());
    double nonzero = 0.0;
    for (const auto& matrix : training.X) {
        nonzero += static_cast<double>(matrix.nonZeros());
    }
    const double gram_work = tasks * p * p;
    const double sparse_work = 2.0 * std::max(nonzero, 1.0);
    const double gram_bytes = tasks * p * p * sizeof(double);
    return p <= 2048.0 &&
        gram_bytes <= 512.0 * 1024.0 * 1024.0 &&
        gram_work <= 4.0 * sparse_work;
}

std::vector<ValidationTaskStats> make_validation_stats(
    const ScaledData& validation
) {
    std::vector<ValidationTaskStats> out;
    out.reserve(validation.X.size());
    for (int task = 0; task < static_cast<int>(validation.X.size()); ++task) {
        const ColSparse& X = validation.X[task];
        const VectorXd& y = validation.y[task];
        if (X.rows() != y.size() || X.rows() < 1) {
            Rcpp::stop("Validation sufficient statistics are not aligned.");
        }
        ValidationTaskStats value;
        value.n = X.rows();
        value.sum_x = VectorXd::Zero(X.cols());
        for (int column = 0; column < X.outerSize(); ++column) {
            for (ColSparse::InnerIterator it(X, column); it; ++it) {
                value.sum_x[column] += it.value();
            }
        }
        value.xtx = MatrixXd(X.transpose() * X);
        value.xty = X.transpose() * y;
        value.sum_y = y.sum();
        value.sum_y2 = y.squaredNorm();
        if (!value.sum_x.allFinite() || !value.xtx.allFinite() ||
            !value.xty.allFinite() || !std::isfinite(value.sum_y) ||
            !std::isfinite(value.sum_y2)) {
            Rcpp::stop("Validation sufficient statistics are non-finite.");
        }
        out.emplace_back(std::move(value));
    }
    return out;
}

double validation_mse_from_stats(
    const ValidationTaskStats& stats,
    const VectorXd& beta,
    double intercept
) {
    if (beta.size() != stats.sum_x.size() || !beta.allFinite() ||
        !std::isfinite(intercept)) {
        Rcpp::stop("Validation coefficients are not aligned or finite.");
    }
    const double intercept_response =
        -2.0 * intercept * stats.sum_y;
    const double coefficient_response =
        -2.0 * beta.dot(stats.xty);
    const double intercept_square =
        static_cast<double>(stats.n) * intercept * intercept;
    const double intercept_coefficient =
        2.0 * intercept * beta.dot(stats.sum_x);
    const double coefficient_square = beta.dot(stats.xtx * beta);
    double sse = stats.sum_y2 + intercept_response +
        coefficient_response + intercept_square +
        intercept_coefficient + coefficient_square;
    const double scale = std::max(
        1.0,
        std::abs(stats.sum_y2) +
        std::abs(intercept_response) +
        std::abs(coefficient_response) +
        std::abs(intercept_square) +
        std::abs(intercept_coefficient) +
        std::abs(coefficient_square)
    );
    if (sse < 0.0 && sse >= -1e-10 * scale) sse = 0.0;
    if (!std::isfinite(sse) || sse < 0.0) {
        Rcpp::stop("Validation sufficient-statistic SSE is invalid.");
    }
    return sse / static_cast<double>(stats.n);
}

PathBundle fit_and_refit_path(
    const RawTargetData& data,
    const std::vector<std::vector<int>>& training_rows,
    const std::vector<int>& base_columns,
    const Transform& transform,
    const std::vector<int>& initial_keep,
    const ArrayXXi& initial_mask,
    const std::vector<double>& lambda,
    double alpha,
    double condition_mix,
    double active_tol,
    int max_iter,
    double tol_objective,
    double tol_coef
) {
    PathBundle out;
    out.intercept_only = false;
    std::vector<int> kept_columns = map_columns(base_columns, initial_keep);
    VectorXd kept_scale = subset_vector(transform.predictor_scale, initial_keep);
    ScaledData initial = build_scaled_data(
        data,
        training_rows,
        kept_columns,
        kept_scale,
        transform.response_center,
        transform.response_scale
    );
    ArrayXXi actual = true_variance_mask(initial.X, initial_mask);
    out.core_relative = any_estimable(actual);
    if (out.core_relative.empty()) {
        out.intercept_only = true;
        out.estimability_core = ArrayXXi(0, data.tasks);
        out.training = std::move(initial);
        return out;
    }
    out.estimability_core = subset_mask_rows(actual, out.core_relative);
    out.training.y = initial.y;
    out.training.X.reserve(data.tasks);
    for (int task = 0; task < data.tasks; ++task) {
        out.training.X.emplace_back(subset_columns(
            initial.X[task], out.core_relative
        ));
    }
    Rcpp::List cache = make_refit_cache(out.training);
    Rcpp::List solver;
    if (use_gram_solver(out.training)) {
        solver = fit_multitask_path_gram(
            out.training,
            cache,
            out.estimability_core,
            lambda,
            alpha,
            condition_mix,
            max_iter,
            tol_objective,
            tol_coef
        );
    } else {
        solver = condition_fit_multitask_path_cpp(
            to_matrix_list(out.training.X),
            to_vector_list(out.training.y),
            to_numeric(lambda),
            alpha,
            condition_mix,
            "equal",
            to_logical_matrix(out.estimability_core),
            max_iter,
            tol_objective,
            tol_coef,
            false
        );
        solver["backend"] = "cpp_eigen_sparse_matrix_free_fista";
    }
    Rcpp::List fits = solver["fits"];
    Rcpp::List beta_path(fits.size());
    Rcpp::NumericVector ridge(fits.size());
    for (int index = 0; index < fits.size(); ++index) {
        Rcpp::List fit = fits[index];
        beta_path[index] = fit["beta"];
        const double current_lambda = Rcpp::as<double>(fit["lambda"]);
        ridge[index] = std::max(current_lambda * (1.0 - alpha), 1e-6);
    }
    Rcpp::List refits = condition_refit_path_cpp(
        beta_path,
        to_logical_matrix(out.estimability_core),
        ridge,
        active_tol,
        cache,
        1e-8
    );
    for (int index = 0; index < refits.size(); ++index) {
        Rcpp::List refit = refits[index];
        if (!Rcpp::as<bool>(refit["converged"])) {
            Rcpp::stop(
                "Compiled target engine refit failed residual verification."
            );
        }
    }
    out.solver = solver;
    out.refits = refits;
    return out;
}

Rcpp::List transform_to_list(const Transform& transform) {
    Rcpp::List task(transform.x_mean.size());
    for (int index = 0; index < static_cast<int>(transform.x_mean.size()); ++index) {
        task[index] = Rcpp::List::create(
            Rcpp::Named("x_mean") = transform.x_mean[index],
            Rcpp::Named("x_variance") = transform.x_variance[index],
            Rcpp::Named("y_mean") = transform.y_mean[index],
            Rcpp::Named("y_variance") = transform.y_variance[index]
        );
    }
    Rcpp::LogicalVector predictor_estimable(
        transform.predictor_estimable.size()
    );
    for (int i = 0; i < predictor_estimable.size(); ++i) {
        predictor_estimable[i] = transform.predictor_estimable[i] != 0;
    }
    Rcpp::NumericVector condition_weights(
        transform.x_mean.size(),
        1.0 / static_cast<double>(transform.x_mean.size())
    );
    return Rcpp::List::create(
        Rcpp::Named("task") = task,
        Rcpp::Named("predictor_center") = transform.predictor_center,
        Rcpp::Named("predictor_scale") = transform.predictor_scale,
        Rcpp::Named("predictor_estimable") = predictor_estimable,
        Rcpp::Named("response_center") = transform.response_center,
        Rcpp::Named("response_scale") = transform.response_scale,
        Rcpp::Named("condition_weights") = condition_weights,
        Rcpp::Named("center_policy") = "equal_condition_mean",
        Rcpp::Named("scale_policy") =
            "equal_condition_within_population_variance",
        Rcpp::Named("transform_policy") =
            "equal_condition_center_equal_condition_within_variance_v1",
        Rcpp::Named("training_fold_only") = true,
        Rcpp::Named("predictor_center_implementation") =
            "implicit_in_condition_intercept_and_projection_shift"
    );
}

std::vector<std::vector<int>> rows_for_inner_fold(
    const std::vector<std::vector<int>>& base_rows,
    const Rcpp::List& folds,
    int fold,
    bool validation
) {
    if (folds.size() != static_cast<int>(base_rows.size())) {
        Rcpp::stop("Inner fold plan task count is not aligned.");
    }
    std::vector<std::vector<int>> out(base_rows.size());
    for (int task = 0; task < static_cast<int>(base_rows.size()); ++task) {
        out[task] = map_positions(
            base_rows[task],
            Rcpp::IntegerVector(folds[task]),
            fold,
            validation
        );
    }
    return out;
}

int validate_fold_plan(
    const Rcpp::List& folds,
    const std::vector<std::vector<int>>& base_rows,
    const std::string& label
) {
    if (folds.size() != static_cast<int>(base_rows.size())) {
        Rcpp::stop(label + " task count is not aligned.");
    }
    int maximum = 0;
    for (int task = 0; task < folds.size(); ++task) {
        Rcpp::IntegerVector value(folds[task]);
        if (value.size() != static_cast<int>(base_rows[task].size())) {
            Rcpp::stop(label + " is not aligned to condition rows.");
        }
        for (int current : value) {
            if (current == NA_INTEGER || current < 1) {
                Rcpp::stop(label + " contains an invalid fold label.");
            }
            maximum = std::max(maximum, current);
        }
    }
    if (maximum < 2) Rcpp::stop(label + " requires at least two folds.");
    for (int fold = 1; fold <= maximum; ++fold) {
        for (int task = 0; task < folds.size(); ++task) {
            Rcpp::IntegerVector value(folds[task]);
            if (std::find(value.begin(), value.end(), fold) == value.end()) {
                Rcpp::stop(label + " omits a fold in one condition.");
            }
        }
    }
    return maximum;
}

LambdaSelection select_lambda_nested(
    const RawTargetData& data,
    const std::vector<std::vector<int>>& base_rows,
    const std::vector<int>& base_columns,
    const ArrayXXi& coefficient_mask,
    const Rcpp::List& fold_plan,
    const std::vector<double>& lambda,
    double alpha,
    double condition_mix,
    double active_tol,
    const std::string& lambda_selection,
    int max_iter,
    double tol_objective,
    double tol_coef
) {
    LambdaSelection out;
    out.effective_folds = validate_fold_plan(
        fold_plan, base_rows, "Inner fold plan"
    );
    out.fold_loss = MatrixXd::Constant(
        out.effective_folds,
        lambda.size(),
        std::numeric_limits<double>::quiet_NaN()
    );
    out.fold_transform = Rcpp::List(out.effective_folds);
    out.fold_backend.assign(out.effective_folds, "not_run");
    out.reason = "";

    ArrayXXi actual = raw_estimability_mask(
        data, base_rows, base_columns, coefficient_mask
    );
    bool any_actual = false;
    for (int row = 0; row < actual.rows(); ++row) {
        for (int task = 0; task < actual.cols(); ++task) {
            any_actual = any_actual || actual(row, task);
        }
    }
    if (!any_actual) {
        out.selected_index = 0;
        out.minimum_index = 0;
        out.one_se_index = 0;
        out.cv_mean.assign(
            lambda.size(), std::numeric_limits<double>::quiet_NaN()
        );
        out.cv_se = out.cv_mean;
        std::fill(
            out.fold_backend.begin(), out.fold_backend.end(),
            "intercept_only_all_predictors_structural_zero"
        );
        out.reason = "intercept_only_all_predictors_structural_zero";
        return out;
    }

    for (int fold = 1; fold <= out.effective_folds; ++fold) {
        std::vector<std::vector<int>> train_rows = rows_for_inner_fold(
            base_rows, fold_plan, fold, false
        );
        std::vector<std::vector<int>> valid_rows = rows_for_inner_fold(
            base_rows, fold_plan, fold, true
        );
        Transform transform = compute_transform(
            data, train_rows, base_columns, coefficient_mask
        );
        out.fold_transform[fold - 1] = transform_to_list(transform);
        std::vector<int> keep = transform_keep(transform);
        if (keep.empty()) {
            out.fold_backend[fold - 1] = "intercept_only";
            continue;
        }
        ArrayXXi raw_mask = subset_mask_rows(transform.estimability, keep);
        PathBundle path = fit_and_refit_path(
            data,
            train_rows,
            base_columns,
            transform,
            keep,
            raw_mask,
            lambda,
            alpha,
            condition_mix,
            active_tol,
            max_iter,
            tol_objective,
            tol_coef
        );
        out.fold_backend[fold - 1] = path.intercept_only ?
            "intercept_only" :
            Rcpp::as<std::string>(path.solver["backend"]);
        std::vector<int> kept_global = map_columns(base_columns, keep);
        VectorXd kept_scale = subset_vector(transform.predictor_scale, keep);
        ScaledData validation = build_scaled_data(
            data,
            valid_rows,
            kept_global,
            kept_scale,
            transform.response_center,
            transform.response_scale
        );
        if (path.intercept_only) {
            double mean_mse = 0.0;
            for (int task = 0; task < data.tasks; ++task) {
                const double intercept = path.training.y[task].mean();
                mean_mse += (
                    validation.y[task].array() - intercept
                ).square().mean() / static_cast<double>(data.tasks);
            }
            for (int lambda_index = 0;
                 lambda_index < static_cast<int>(lambda.size()); ++lambda_index) {
                out.fold_loss(fold - 1, lambda_index) = mean_mse;
            }
            continue;
        }
        ScaledData validation_core;
        validation_core.y = validation.y;
        validation_core.X.reserve(data.tasks);
        for (int task = 0; task < data.tasks; ++task) {
            validation_core.X.emplace_back(subset_columns(
                validation.X[task], path.core_relative
            ));
        }
        std::vector<ValidationTaskStats> validation_stats =
            make_validation_stats(validation_core);
        for (int lambda_index = 0;
             lambda_index < static_cast<int>(lambda.size()); ++lambda_index) {
            Rcpp::List refit = path.refits[lambda_index];
            MatrixXd beta = Rcpp::as<MatrixXd>(refit["beta"]);
            VectorXd intercept = Rcpp::as<VectorXd>(refit["intercept"]);
            double mean_mse = 0.0;
            for (int task = 0; task < data.tasks; ++task) {
                mean_mse += validation_mse_from_stats(
                    validation_stats[task],
                    beta.col(task),
                    intercept[task]
                ) / static_cast<double>(data.tasks);
            }
            if (!std::isfinite(mean_mse)) {
                Rcpp::stop("Inner validation loss is non-finite.");
            }
            out.fold_loss(fold - 1, lambda_index) = mean_mse;
        }
    }

    out.cv_mean.resize(lambda.size());
    out.cv_se.resize(lambda.size());
    bool any_finite = false;
    for (int column = 0; column < static_cast<int>(lambda.size()); ++column) {
        std::vector<double> finite;
        finite.reserve(out.effective_folds);
        for (int row = 0; row < out.effective_folds; ++row) {
            const double value = out.fold_loss(row, column);
            if (std::isfinite(value)) finite.push_back(value);
        }
        if (finite.empty()) {
            out.cv_mean[column] = std::numeric_limits<double>::quiet_NaN();
            out.cv_se[column] = std::numeric_limits<double>::quiet_NaN();
            continue;
        }
        any_finite = true;
        const double mean = std::accumulate(
            finite.begin(), finite.end(), 0.0
        ) / static_cast<double>(finite.size());
        out.cv_mean[column] = mean;
        if (finite.size() < 2) {
            out.cv_se[column] = 0.0;
        } else {
            double square = 0.0;
            for (double value : finite) {
                square += (value - mean) * (value - mean);
            }
            const double sd = std::sqrt(
                square / static_cast<double>(finite.size() - 1)
            );
            out.cv_se[column] = sd / std::sqrt(
                static_cast<double>(finite.size())
            );
        }
    }
    if (!any_finite) {
        Rcpp::stop("Inner cross-validation produced no finite validation loss.");
    }
    out.minimum_index = -1;
    double minimum = std::numeric_limits<double>::infinity();
    for (int index = 0; index < static_cast<int>(lambda.size()); ++index) {
        if (std::isfinite(out.cv_mean[index]) && out.cv_mean[index] < minimum) {
            minimum = out.cv_mean[index];
            out.minimum_index = index;
        }
    }
    if (out.minimum_index < 0) {
        Rcpp::stop("Inner cross-validation did not identify a minimum.");
    }
    const double threshold = out.cv_mean[out.minimum_index] +
        out.cv_se[out.minimum_index];
    out.one_se_index = -1;
    for (int index = 0; index < static_cast<int>(lambda.size()); ++index) {
        if (std::isfinite(out.cv_mean[index]) &&
            out.cv_mean[index] <= threshold) {
            out.one_se_index = index;
            break;
        }
    }
    if (out.one_se_index < 0) {
        Rcpp::stop("Inner cross-validation did not identify lambda.1se.");
    }
    out.selected_index = lambda_selection == "lambda.min" ?
        out.minimum_index : out.one_se_index;
    return out;
}

Rcpp::List selection_to_list(
    const LambdaSelection& selection,
    const std::vector<double>& lambda
) {
    Rcpp::List out = Rcpp::List::create(
        Rcpp::Named("selected_index") = selection.selected_index + 1,
        Rcpp::Named("selected_lambda") = lambda[selection.selected_index],
        Rcpp::Named("lambda_min") = lambda[selection.minimum_index],
        Rcpp::Named("lambda_1se") = lambda[selection.one_se_index],
        Rcpp::Named("cv_mean") = to_numeric(selection.cv_mean),
        Rcpp::Named("cv_se") = to_numeric(selection.cv_se),
        Rcpp::Named("fold_loss") = selection.fold_loss,
        Rcpp::Named("fold_transform") = selection.fold_transform,
        Rcpp::Named("solver_backend") = Rcpp::wrap(selection.fold_backend),
        Rcpp::Named("effective_nfolds") = selection.effective_folds
    );
    if (!selection.reason.empty()) {
        out["selection_reason"] = selection.reason;
    }
    return out;
}

MatrixXd expand_matrix(
    const MatrixXd& core,
    int full_rows,
    int tasks,
    const std::vector<int>& initial_keep,
    const std::vector<int>& core_relative
) {
    MatrixXd out = MatrixXd::Zero(full_rows, tasks);
    for (int row = 0; row < static_cast<int>(core_relative.size()); ++row) {
        out.row(initial_keep[core_relative[row]]) = core.row(row);
    }
    return out;
}

VectorXd expand_vector(
    const VectorXd& core,
    int full_rows,
    const std::vector<int>& initial_keep,
    const std::vector<int>& core_relative
) {
    VectorXd out = VectorXd::Zero(full_rows);
    for (int row = 0; row < static_cast<int>(core_relative.size()); ++row) {
        out[initial_keep[core_relative[row]]] = core[row];
    }
    return out;
}

ArrayXXi expand_mask(
    const ArrayXXi& core,
    int full_rows,
    int tasks,
    const std::vector<int>& initial_keep,
    const std::vector<int>& core_relative
) {
    ArrayXXi out = ArrayXXi::Zero(full_rows, tasks);
    for (int row = 0; row < static_cast<int>(core_relative.size()); ++row) {
        out.row(initial_keep[core_relative[row]]) = core.row(row);
    }
    return out;
}

Rcpp::List intercept_only_fit(
    const ScaledData& training,
    double lambda,
    int predictors,
    int tasks
) {
    VectorXd intercept(tasks);
    double objective = 0.0;
    for (int task = 0; task < tasks; ++task) {
        intercept[task] = training.y[task].mean();
        VectorXd residual = training.y[task].array() - intercept[task];
        objective += 0.5 * residual.squaredNorm() /
            static_cast<double>(training.y[task].size());
    }
    return Rcpp::List::create(
        Rcpp::Named("beta") = MatrixXd::Zero(predictors, tasks),
        Rcpp::Named("intercept") = intercept,
        Rcpp::Named("lambda") = lambda,
        Rcpp::Named("objective") = objective,
        Rcpp::Named("objective_change") = 0.0,
        Rcpp::Named("coef_change") = 0.0,
        Rcpp::Named("iterations") = 0,
        Rcpp::Named("converged") = true,
        Rcpp::Named("step") = 1.0,
        Rcpp::Named("backend") = "intercept_only"
    );
}

Rcpp::List finish_full_refit(
    const PathBundle& path,
    const Rcpp::List& selected_refit,
    const Rcpp::List& selected_fit,
    const Transform& transform,
    const std::vector<int>& initial_keep,
    int full_predictors,
    int tasks,
    double active_tol,
    double ridge
) {
    MatrixXd beta = MatrixXd::Zero(full_predictors, tasks);
    VectorXd shared = VectorXd::Zero(full_predictors);
    ArrayXXi support = ArrayXXi::Zero(full_predictors, tasks);
    VectorXd intercept(tasks);
    double coef_change = 0.0;
    int iterations = 0;
    bool converged = true;
    std::string backend = "intercept_only";
    if (path.intercept_only) {
        for (int task = 0; task < tasks; ++task) {
            intercept[task] = path.training.y[task].mean();
        }
    } else {
        MatrixXd beta_core = Rcpp::as<MatrixXd>(selected_refit["beta"]);
        VectorXd shared_core = Rcpp::as<VectorXd>(selected_refit["shared"]);
        Rcpp::LogicalMatrix support_input(selected_refit["support_mask"]);
        ArrayXXi support_core(support_input.nrow(), support_input.ncol());
        for (int row = 0; row < support_input.nrow(); ++row) {
            for (int task = 0; task < support_input.ncol(); ++task) {
                support_core(row, task) = support_input(row, task) ? 1 : 0;
            }
        }
        beta = expand_matrix(
            beta_core, full_predictors, tasks,
            initial_keep, path.core_relative
        );
        shared = expand_vector(
            shared_core, full_predictors,
            initial_keep, path.core_relative
        );
        support = expand_mask(
            support_core, full_predictors, tasks,
            initial_keep, path.core_relative
        );
        intercept = Rcpp::as<VectorXd>(selected_refit["intercept"]);
        coef_change = Rcpp::as<double>(selected_refit["coef_change"]);
        iterations = Rcpp::as<int>(selected_refit["iterations"]);
        converged = Rcpp::as<bool>(selected_refit["converged"]);
        backend = Rcpp::as<std::string>(selected_refit["backend"]);
    }
    ArrayXXi estimability = ArrayXXi::Zero(full_predictors, tasks);
    if (!path.intercept_only) {
        estimability = expand_mask(
            path.estimability_core,
            full_predictors,
            tasks,
            initial_keep,
            path.core_relative
        );
    }
    MatrixXd beta_condition = beta;
    MatrixXd delta = MatrixXd::Zero(full_predictors, tasks);
    ArrayXXi active = ArrayXXi::Zero(full_predictors, tasks);
    for (int row = 0; row < full_predictors; ++row) {
        for (int task = 0; task < tasks; ++task) {
            if (!estimability(row, task)) {
                beta_condition(row, task) = NA_REAL;
                delta(row, task) = NA_REAL;
            } else {
                delta(row, task) = beta(row, task) - shared[row];
                active(row, task) = std::abs(beta(row, task)) > active_tol;
            }
        }
    }
    return Rcpp::List::create(
        Rcpp::Named("beta") = beta,
        Rcpp::Named("beta_condition") = beta_condition,
        Rcpp::Named("beta_shared") = shared,
        Rcpp::Named("delta_condition") = delta,
        Rcpp::Named("support_mask") = to_logical_matrix(support),
        Rcpp::Named("active_mask") = to_logical_matrix(active),
        Rcpp::Named("estimability_mask") = to_logical_matrix(estimability),
        Rcpp::Named("intercept") = intercept,
        Rcpp::Named("ridge") = ridge,
        Rcpp::Named("common_metric") =
            "pooled_weighted_predictor_gram_cpp_direct_schur",
        Rcpp::Named("iterations") = iterations,
        Rcpp::Named("coef_change") = coef_change,
        Rcpp::Named("converged") = converged,
        Rcpp::Named("backend") = backend
    );
}

double rsq(const VectorXd& y, const VectorXd& prediction) {
    const double mean = y.mean();
    const double denominator = (y.array() - mean).square().sum();
    if (!std::isfinite(denominator) || denominator <= kEps) {
        return NA_REAL;
    }
    return 1.0 - (y - prediction).squaredNorm() / denominator;
}

Rcpp::List run_outer_cv(
    const RawTargetData& data,
    const ArrayXXi& full_mask,
    const Rcpp::List& outer_folds,
    const Rcpp::List& outer_inner,
    const std::vector<int>& comparison,
    const std::vector<double>& global_lambda,
    bool lambda_auto,
    int nlambda,
    double lambda_min_ratio,
    double alpha,
    double condition_mix,
    double active_tol,
    const std::string& lambda_selection,
    int max_iter,
    double tol_objective,
    double tol_coef
) {
    std::vector<std::vector<int>> complete_rows = all_rows(data);
    const int outer_count = validate_fold_plan(
        outer_folds, complete_rows, "Outer fold plan"
    );
    if (outer_inner.size() != outer_count) {
        Rcpp::stop("Outer-inner fold plan count is not aligned.");
    }
    std::vector<int> full_columns = sequence(data.predictors);
    std::vector<VectorXd> prediction(data.tasks);
    std::vector<VectorXd> projection_full(data.tasks);
    std::vector<VectorXd> projection_common(data.tasks);
    std::vector<VectorXd> projection_global(data.tasks);
    std::vector<Eigen::VectorXi> assignment(data.tasks);
    for (int task = 0; task < data.tasks; ++task) {
        prediction[task] = VectorXd::Constant(
            data.y[task].size(), std::numeric_limits<double>::quiet_NaN()
        );
        projection_full[task] = prediction[task];
        projection_common[task] = prediction[task];
        projection_global[task] = prediction[task];
        assignment[task] = Eigen::VectorXi::Zero(data.y[task].size());
    }
    Rcpp::List fold_transform(outer_count);
    Rcpp::NumericVector fold_selected_lambda(
        outer_count, NA_REAL
    );
    Rcpp::List fold_inner_cv(outer_count);
    Rcpp::List fold_support(outer_count);

    for (int fold = 1; fold <= outer_count; ++fold) {
        std::vector<std::vector<int>> train_rows(data.tasks);
        std::vector<std::vector<int>> test_rows(data.tasks);
        for (int task = 0; task < data.tasks; ++task) {
            Rcpp::IntegerVector values(outer_folds[task]);
            train_rows[task] = rows_from_fold(values, fold, false);
            test_rows[task] = rows_from_fold(values, fold, true);
        }
        Transform transform = compute_transform(
            data, train_rows, full_columns, full_mask
        );
        fold_transform[fold - 1] = transform_to_list(transform);
        std::vector<int> keep = transform_keep(transform);
        ArrayXXi structural = (transform.estimability == 0).cast<int>();
        Rcpp::LogicalVector common_pair(data.predictors);
        Rcpp::LogicalVector common_all(data.predictors);
        for (int predictor = 0; predictor < data.predictors; ++predictor) {
            bool pair = true;
            for (int task : comparison) {
                pair = pair && transform.estimability(predictor, task);
            }
            bool global = true;
            for (int task = 0; task < data.tasks; ++task) {
                global = global && transform.estimability(predictor, task);
            }
            common_pair[predictor] = pair;
            common_all[predictor] = global;
        }
        ArrayXXi projection_support = ArrayXXi::Ones(
            data.predictors, data.tasks
        );
        Rcpp::List support = Rcpp::List::create(
            Rcpp::Named("coefficient_estimable_mask") =
                to_logical_matrix(transform.estimability),
            Rcpp::Named("projectable_structural_zero_mask") =
                to_logical_matrix(structural),
            Rcpp::Named("projection_support_mask") =
                to_logical_matrix(projection_support),
            Rcpp::Named("pairwise_or_requested_common") = common_pair,
            Rcpp::Named("global_common") = common_all
        );
        fold_support[fold - 1] = support;
        if (keep.empty()) {
            support["projection_status"] =
                "intercept_only_all_predictors_structural_zero";
            fold_support[fold - 1] = support;
            for (int task = 0; task < data.tasks; ++task) {
                const double intercept_scaled =
                    (transform.y_mean[task] - transform.response_center) /
                    transform.response_scale;
                const double raw_prediction =
                    intercept_scaled * transform.response_scale +
                    transform.response_center;
                for (int row : test_rows[task]) {
                    prediction[task][row] = raw_prediction;
                    projection_full[task][row] = 0.0;
                    projection_common[task][row] = 0.0;
                    projection_global[task][row] = 0.0;
                    assignment[task][row] += 1;
                }
            }
            continue;
        }
        std::vector<int> outer_columns = map_columns(full_columns, keep);
        ArrayXXi raw_mask = subset_mask_rows(transform.estimability, keep);
        VectorXd outer_scale = subset_vector(transform.predictor_scale, keep);
        ScaledData outer_scaled = build_scaled_data(
            data,
            train_rows,
            outer_columns,
            outer_scale,
            transform.response_center,
            transform.response_scale
        );
        std::vector<double> fold_lambda = lambda_auto ?
            make_lambda_path(
                outer_scaled,
                raw_mask,
                alpha,
                condition_mix,
                nlambda,
                lambda_min_ratio
            ) : global_lambda;
        LambdaSelection inner = select_lambda_nested(
            data,
            train_rows,
            outer_columns,
            raw_mask,
            Rcpp::List(outer_inner[fold - 1]),
            fold_lambda,
            alpha,
            condition_mix,
            active_tol,
            lambda_selection,
            max_iter,
            tol_objective,
            tol_coef
        );
        fold_selected_lambda[fold - 1] =
            fold_lambda[inner.selected_index];
        Rcpp::List inner_output = selection_to_list(inner, fold_lambda);
        Rcpp::IntegerVector inner_predictor_index(outer_columns.size());
        for (int index = 0; index < static_cast<int>(outer_columns.size()); ++index) {
            inner_predictor_index[index] = outer_columns[index] + 1;
        }
        inner_output["predictor_index"] = inner_predictor_index;
        fold_inner_cv[fold - 1] = inner_output;
        std::vector<double> selected_lambda{
            fold_lambda[inner.selected_index]
        };
        PathBundle path = fit_and_refit_path(
            data,
            train_rows,
            full_columns,
            transform,
            keep,
            raw_mask,
            selected_lambda,
            alpha,
            condition_mix,
            active_tol,
            max_iter,
            tol_objective,
            tol_coef
        );
        ScaledData test_scaled = build_scaled_data(
            data,
            test_rows,
            outer_columns,
            outer_scale,
            transform.response_center,
            transform.response_scale
        );
        if (path.intercept_only) {
            for (int task = 0; task < data.tasks; ++task) {
                const double intercept = path.training.y[task].mean();
                const double raw_prediction =
                    intercept * transform.response_scale +
                    transform.response_center;
                for (int index = 0; index < static_cast<int>(test_rows[task].size()); ++index) {
                    const int row = test_rows[task][index];
                    prediction[task][row] = raw_prediction;
                    projection_full[task][row] = 0.0;
                    projection_common[task][row] = 0.0;
                    projection_global[task][row] = 0.0;
                    assignment[task][row] += 1;
                }
            }
            continue;
        }
        Rcpp::List refit = path.refits[0];
        MatrixXd beta = Rcpp::as<MatrixXd>(refit["beta"]);
        VectorXd intercept = Rcpp::as<VectorXd>(refit["intercept"]);
        std::vector<ColSparse> test_core;
        test_core.reserve(data.tasks);
        for (int task = 0; task < data.tasks; ++task) {
            test_core.emplace_back(subset_columns(
                test_scaled.X[task], path.core_relative
            ));
        }
        std::vector<int> core_outer_relative;
        core_outer_relative.reserve(path.core_relative.size());
        for (int core : path.core_relative) {
            core_outer_relative.push_back(keep[core]);
        }
        VectorXd core_center(path.core_relative.size());
        VectorXd core_scale(path.core_relative.size());
        for (int index = 0; index < static_cast<int>(path.core_relative.size()); ++index) {
            core_center[index] = transform.predictor_center[
                core_outer_relative[index]
            ];
            core_scale[index] = transform.predictor_scale[
                core_outer_relative[index]
            ];
        }
        for (int task = 0; task < data.tasks; ++task) {
            VectorXd linear = test_core[task] * beta.col(task);
            VectorXd full_score = linear;
            double full_shift = 0.0;
            double common_shift = 0.0;
            double global_shift = 0.0;
            VectorXd common_beta = VectorXd::Zero(beta.rows());
            VectorXd global_beta = VectorXd::Zero(beta.rows());
            for (int predictor = 0; predictor < beta.rows(); ++predictor) {
                if (!path.estimability_core(predictor, task)) continue;
                const double coefficient = beta(predictor, task);
                const double shift =
                    core_center[predictor] / core_scale[predictor] * coefficient;
                full_shift += shift;
                const int full_predictor = core_outer_relative[predictor];
                if (common_pair[full_predictor]) {
                    common_beta[predictor] = coefficient;
                    common_shift += shift;
                }
                if (common_all[full_predictor]) {
                    global_beta[predictor] = coefficient;
                    global_shift += shift;
                }
            }
            VectorXd common_score = test_core[task] * common_beta;
            VectorXd global_score = test_core[task] * global_beta;
            full_score.array() -= full_shift;
            common_score.array() -= common_shift;
            global_score.array() -= global_shift;
            VectorXd raw_prediction = (
                linear.array() + intercept[task]
            ) * transform.response_scale + transform.response_center;
            for (int index = 0; index < static_cast<int>(test_rows[task].size()); ++index) {
                const int row = test_rows[task][index];
                prediction[task][row] = raw_prediction[index];
                projection_full[task][row] = full_score[index];
                projection_common[task][row] = common_score[index];
                projection_global[task][row] = global_score[index];
                assignment[task][row] += 1;
            }
        }
    }

    Rcpp::List prediction_out(data.tasks);
    Rcpp::List full_out(data.tasks);
    Rcpp::List common_out(data.tasks);
    Rcpp::List global_out(data.tasks);
    Rcpp::List assignment_out(data.tasks);
    Rcpp::List fold_out(data.tasks);
    Rcpp::NumericVector coverage(data.tasks);
    Rcpp::NumericVector available(data.tasks);
    for (int task = 0; task < data.tasks; ++task) {
        for (int row = 0; row < assignment[task].size(); ++row) {
            if (assignment[task][row] != 1) {
                Rcpp::stop(
                    "Every cell must receive exactly one outer-fold projection."
                );
            }
        }
        prediction_out[task] = prediction[task];
        full_out[task] = projection_full[task];
        common_out[task] = projection_common[task];
        global_out[task] = projection_global[task];
        assignment_out[task] = assignment[task];
        fold_out[task] = outer_folds[task];
        coverage[task] = 1.0;
        int finite = 0;
        for (int row = 0; row < projection_full[task].size(); ++row) {
            finite += std::isfinite(projection_full[task][row]);
        }
        available[task] = static_cast<double>(finite) /
            static_cast<double>(projection_full[task].size());
    }
    return Rcpp::List::create(
        Rcpp::Named("oof_prediction") = prediction_out,
        Rcpp::Named("projection_condition_full_oof") = full_out,
        Rcpp::Named("projection_common_oof") = common_out,
        Rcpp::Named("projection_global_common_oof") = global_out,
        Rcpp::Named("oof_fold") = fold_out,
        Rcpp::Named("outer_nfolds") = outer_count,
        Rcpp::Named("fold_transform") = fold_transform,
        Rcpp::Named("fold_selected_lambda") = fold_selected_lambda,
        Rcpp::Named("fold_inner_cv") = fold_inner_cv,
        Rcpp::Named("fold_support") = fold_support,
        Rcpp::Named("oof_assignment_count") = assignment_out,
        Rcpp::Named("oof_cell_coverage") = coverage,
        Rcpp::Named("oof_projection_available_fraction") = available,
        Rcpp::Named("projection_origin") =
            "outer_condition_stratified_cell_oof",
        Rcpp::Named("primary_projection") = "condition_full_oof",
        Rcpp::Named("common_projection_role") =
            "shared_estimable_component",
        Rcpp::Named("condition_unique_projection_role") =
            "condition_full_oof_minus_common_support_oof",
        Rcpp::Named("nonestimable_projection_policy") =
            "structural_zero_by_condition",
        Rcpp::Named("projection_used_for_penalty") = true,
        Rcpp::Named("full_fit_projection_used_for_penalty") = false,
        Rcpp::Named("fold_transform_policy") =
            "equal_condition_center_equal_condition_within_variance_v1",
        Rcpp::Named("predictor_center_implementation") =
            "implicit_in_condition_intercept_and_projection_shift",
        Rcpp::Named("oof_model") =
            "nested_selection_cached_refit_heldout_condition_full_projection",
        Rcpp::Named("engine_backend") =
            "cpp_eigen_fused_target_nested_cv_hybrid_gram_sufficient_statistics"
    );
}

} // namespace

// [[Rcpp::export]]
Rcpp::List condition_fit_target_engine_cpp(
    Rcpp::List X_list,
    Rcpp::List y_list,
    Rcpp::LogicalMatrix coefficient_mask,
    Rcpp::IntegerVector comparison_index,
    Rcpp::Nullable<Rcpp::NumericVector> lambda,
    bool lambda_auto,
    int nlambda,
    double lambda_min_ratio,
    double alpha,
    double condition_mix,
    double active_tol,
    Rcpp::List fold_plan,
    std::string lambda_selection,
    int max_iter,
    double tol_objective,
    double tol_coef
) {
    if (!std::isfinite(alpha) || alpha < 0.0 || alpha > 1.0 ||
        !std::isfinite(condition_mix) || condition_mix < 0.0 ||
        condition_mix > 1.0 || !std::isfinite(active_tol) ||
        active_tol < 0.0 || max_iter < 1 ||
        !std::isfinite(tol_objective) || tol_objective <= 0.0 ||
        !std::isfinite(tol_coef) || tol_coef <= 0.0) {
        Rcpp::stop("Invalid fused target-engine controls.");
    }
    if (lambda_selection != "lambda.1se" &&
        lambda_selection != "lambda.min") {
        Rcpp::stop("Unsupported lambda selection rule.");
    }
    RawTargetData data = parse_raw_target(X_list, y_list);
    ArrayXXi mask = parse_mask(
        coefficient_mask, data.predictors, data.tasks
    );
    std::vector<int> comparison;
    comparison.reserve(comparison_index.size());
    for (int value : comparison_index) {
        if (value == NA_INTEGER || value < 0 || value >= data.tasks) {
            Rcpp::stop("comparison_index is out of bounds.");
        }
        if (std::find(comparison.begin(), comparison.end(), value) ==
            comparison.end()) {
            comparison.push_back(value);
        }
    }
    if (comparison.size() < 2) {
        Rcpp::stop("At least two comparison conditions are required.");
    }
    if (!fold_plan.containsElementNamed("outer") ||
        !fold_plan.containsElementNamed("outer_inner") ||
        !fold_plan.containsElementNamed("full_inner")) {
        Rcpp::stop("Target fold plan is incomplete.");
    }
    std::vector<std::vector<int>> complete_rows = all_rows(data);
    std::vector<int> full_columns = sequence(data.predictors);
    Transform full_transform = compute_transform(
        data, complete_rows, full_columns, mask
    );
    std::vector<int> full_keep = transform_keep(full_transform);
    if (full_keep.empty()) {
        Rcpp::stop(
            "No target predictor remains estimable in the full condition design."
        );
    }
    std::vector<int> full_kept_columns = map_columns(full_columns, full_keep);
    VectorXd full_kept_scale = subset_vector(
        full_transform.predictor_scale, full_keep
    );
    ScaledData full_scaled = build_scaled_data(
        data,
        complete_rows,
        full_kept_columns,
        full_kept_scale,
        full_transform.response_center,
        full_transform.response_scale
    );
    ArrayXXi full_raw_mask = subset_mask_rows(
        full_transform.estimability, full_keep
    );
    std::vector<double> lambda_path = sorted_lambda(lambda);
    if (lambda_auto) {
        lambda_path = make_lambda_path(
            full_scaled,
            full_raw_mask,
            alpha,
            condition_mix,
            nlambda,
            lambda_min_ratio
        );
    } else if (lambda_path.empty()) {
        Rcpp::stop("A non-empty lambda vector is required when lambda_auto is false.");
    }

    Rcpp::List cv = run_outer_cv(
        data,
        mask,
        Rcpp::List(fold_plan["outer"]),
        Rcpp::List(fold_plan["outer_inner"]),
        comparison,
        lambda_path,
        lambda_auto,
        nlambda,
        lambda_min_ratio,
        alpha,
        condition_mix,
        active_tol,
        lambda_selection,
        max_iter,
        tol_objective,
        tol_coef
    );

    LambdaSelection full_cv = select_lambda_nested(
        data,
        complete_rows,
        full_columns,
        mask,
        Rcpp::List(fold_plan["full_inner"]),
        lambda_path,
        alpha,
        condition_mix,
        active_tol,
        lambda_selection,
        max_iter,
        tol_objective,
        tol_coef
    );
    const int selected_index = full_cv.selected_index;
    std::vector<double> required_lambda(
        lambda_path.begin(), lambda_path.begin() + selected_index + 1
    );
    PathBundle final_path = fit_and_refit_path(
        data,
        complete_rows,
        full_columns,
        full_transform,
        full_keep,
        full_raw_mask,
        required_lambda,
        alpha,
        condition_mix,
        active_tol,
        max_iter,
        tol_objective,
        tol_coef
    );
    Rcpp::List selected_fit;
    Rcpp::List selected_refit;
    MatrixXd beta_selection = MatrixXd::Zero(
        data.predictors, data.tasks
    );
    if (final_path.intercept_only) {
        selected_fit = intercept_only_fit(
            final_path.training,
            lambda_path[selected_index],
            data.predictors,
            data.tasks
        );
        selected_refit = Rcpp::List::create();
    } else {
        Rcpp::List fits = final_path.solver["fits"];
        selected_fit = Rcpp::List(fits[fits.size() - 1]);
        selected_refit = Rcpp::List(
            final_path.refits[final_path.refits.size() - 1]
        );
        MatrixXd beta_core = Rcpp::as<MatrixXd>(selected_fit["beta"]);
        beta_selection = expand_matrix(
            beta_core,
            data.predictors,
            data.tasks,
            full_keep,
            final_path.core_relative
        );
        selected_fit["beta"] = beta_selection;
        selected_fit["backend"] = "cpp_eigen_fused_target_engine";
    }
    const double selected_ridge = std::max(
        lambda_path[selected_index] * (1.0 - alpha), 1e-6
    );
    Rcpp::List refit = finish_full_refit(
        final_path,
        selected_refit,
        selected_fit,
        full_transform,
        full_keep,
        data.predictors,
        data.tasks,
        active_tol,
        selected_ridge
    );

    Rcpp::List prediction = cv["oof_prediction"];
    Rcpp::NumericVector condition_rsq_oof(data.tasks);
    Rcpp::NumericVector condition_rmse_oof(data.tasks);
    double residual_ss = 0.0;
    double within_ss = 0.0;
    for (int task = 0; task < data.tasks; ++task) {
        VectorXd fitted = Rcpp::as<VectorXd>(prediction[task]);
        condition_rsq_oof[task] = rsq(data.y[task], fitted);
        condition_rmse_oof[task] = std::sqrt(
            (data.y[task] - fitted).squaredNorm() /
            static_cast<double>(data.y[task].size())
        );
        residual_ss += (data.y[task] - fitted).squaredNorm();
        within_ss += (
            data.y[task].array() - data.y[task].mean()
        ).square().sum();
    }
    const double pooled_rsq = within_ss <= kEps ? NA_REAL :
        1.0 - residual_ss / within_ss;

    Rcpp::List full_cv_output = selection_to_list(full_cv, lambda_path);
    Rcpp::IntegerVector full_predictor_index(data.predictors);
    for (int predictor = 0; predictor < data.predictors; ++predictor) {
        full_predictor_index[predictor] = predictor + 1;
    }
    full_cv_output["predictor_index"] = full_predictor_index;

    return Rcpp::List::create(
        Rcpp::Named("lambda_path") = to_numeric(lambda_path),
        Rcpp::Named("cv") = cv,
        Rcpp::Named("full_cv") = full_cv_output,
        Rcpp::Named("selected") = selected_fit,
        Rcpp::Named("beta_selection") = beta_selection,
        Rcpp::Named("refit") = refit,
        Rcpp::Named("full_transform") = transform_to_list(full_transform),
        Rcpp::Named("condition_rsq_oof") = condition_rsq_oof,
        Rcpp::Named("condition_rmse_oof") = condition_rmse_oof,
        Rcpp::Named("target_rsq_oof_pooled") = pooled_rsq,
        Rcpp::Named("inner_cv_backend") =
            "hybrid_centered_gram_or_sparse_fista_with_validation_sufficient_statistics",
        Rcpp::Named("backend") =
            "cpp_eigen_fused_target_nested_cv_hybrid_gram_refit_validation_stats"
    );
}
