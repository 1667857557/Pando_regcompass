#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

using namespace Rcpp;

// [[Rcpp::export]]
List condition_cpp_predictors_scaling(
        List gene_conditions, List peak_conditions,
        IntegerVector tf_index, IntegerVector peak_index,
        double scale_floor) {
    const int k = gene_conditions.size();
    const int p = tf_index.size();
    if (k < 1 || peak_conditions.size() != k || p < 1 ||
        peak_index.size() != p || !R_finite(scale_floor) || scale_floor <= 0) {
        stop("Invalid native predictor/scaling input.");
    }
    List predictors(k);
    NumericVector center(p, 0.0), within_variance(p, 0.0);
    double total_n = 0.0;
    std::vector<int> sizes(k);

    for (int condition = 0; condition < k; ++condition) {
        NumericMatrix gene = gene_conditions[condition];
        NumericMatrix peak = peak_conditions[condition];
        const int n = gene.nrow();
        if (n < 1 || peak.nrow() != n) {
            stop("Native condition matrices are misaligned.");
        }
        sizes[condition] = n;
        total_n += n;
        NumericMatrix x(n, p);
        for (int edge = 0; edge < p; ++edge) {
            const int tf = tf_index[edge];
            const int region = peak_index[edge];
            if (tf < 0 || tf >= gene.ncol() ||
                region < 0 || region >= peak.ncol()) {
                stop("Native predictor column index is out of range.");
            }
            double sum = 0.0;
            for (int row = 0; row < n; ++row) {
                const double value = gene(row, tf) * peak(row, region);
                if (!R_finite(value)) stop("Native predictor is non-finite.");
                x(row, edge) = value;
                sum += value;
            }
            center[edge] += sum;
        }
        predictors[condition] = x;
    }
    for (int edge = 0; edge < p; ++edge) center[edge] /= total_n;

    for (int condition = 0; condition < k; ++condition) {
        NumericMatrix x = predictors[condition];
        const int n = sizes[condition];
        for (int edge = 0; edge < p; ++edge) {
            double mean = 0.0;
            for (int row = 0; row < n; ++row) mean += x(row, edge);
            mean /= static_cast<double>(n);
            double sumsq = 0.0;
            for (int row = 0; row < n; ++row) {
                const double difference = x(row, edge) - mean;
                sumsq += difference * difference;
            }
            within_variance[edge] += sumsq / static_cast<double>(n);
        }
    }

    NumericVector scale(p);
    LogicalVector informative(p);
    for (int edge = 0; edge < p; ++edge) {
        within_variance[edge] /= static_cast<double>(k);
        const double value = std::sqrt(std::max(0.0, within_variance[edge]));
        informative[edge] = R_finite(value) && value > scale_floor;
        scale[edge] = informative[edge] ? value : 1.0;
    }
    return List::create(
        _["x"] = predictors, _["center"] = center, _["scale"] = scale,
        _["informative"] = informative,
        _["within_variance"] = within_variance
    );
}

// [[Rcpp::export]]
List condition_cpp_estar_solver(
        NumericMatrix H, NumericVector r, NumericVector information,
        double z, double solver_tol, int solver_max_iter,
        double solver_scale, double lipschitz) {
    const int dimension = r.size();
    if (H.nrow() != dimension || H.ncol() != dimension ||
        information.size() != dimension || !R_finite(z) || z <= 0 ||
        !R_finite(solver_tol) || solver_tol <= 0 || solver_max_iter < 1 ||
        !R_finite(solver_scale) || solver_scale <= 0 ||
        !R_finite(lipschitz) || lipschitz <= 0) {
        stop("Invalid native E-star solver input.");
    }
    std::vector<int> active;
    for (int index = 0; index < dimension; ++index) {
        if (information[index] > 0) active.push_back(index);
    }
    NumericVector delta(dimension, 0.0);
    if (active.empty()) {
        return List::create(
            _["delta"] = delta, _["status"] = "ok", _["iterations"] = 0,
            _["kkt_residual"] = 0.0, _["objective"] = 0.0
        );
    }
    const int n = active.size();
    std::vector<double> current(n, 0.0), accelerated(n, 0.0);
    std::vector<double> next(n), next_accelerated(n), gradient(n), weights(n);
    for (int i = 0; i < n; ++i) {
        weights[i] = std::sqrt(information[active[i]]) / solver_scale;
    }
    double momentum = 1.0;
    double kkt = std::numeric_limits<double>::infinity();
    double objective = std::numeric_limits<double>::infinity();
    bool converged = false;
    int iteration = 0;
    for (iteration = 1; iteration <= solver_max_iter; ++iteration) {
        for (int i = 0; i < n; ++i) {
            double value = -r[active[i]] / solver_scale;
            for (int j = 0; j < n; ++j) {
                value += (H(active[i], active[j]) / solver_scale) * accelerated[j];
            }
            gradient[i] = value;
            const double trial = accelerated[i] - value / lipschitz;
            const double threshold = z * weights[i] / lipschitz;
            next[i] = std::copysign(std::max(std::abs(trial) - threshold, 0.0),
                                    trial);
        }
        const double next_momentum =
            (1.0 + std::sqrt(1.0 + 4.0 * momentum * momentum)) / 2.0;
        double restart_inner = 0.0;
        for (int i = 0; i < n; ++i) {
            next_accelerated[i] = next[i] +
                ((momentum - 1.0) / next_momentum) * (next[i] - current[i]);
            restart_inner += (accelerated[i] - next[i]) *
                (next[i] - current[i]);
        }
        double updated_momentum = next_momentum;
        if (restart_inner > 0) {
            updated_momentum = 1.0;
            next_accelerated = next;
        }
        kkt = 0.0;
        double step_squared = 0.0, current_squared = 0.0;
        for (int i = 0; i < n; ++i) {
            double value = -r[active[i]] / solver_scale;
            for (int j = 0; j < n; ++j) {
                value += (H(active[i], active[j]) / solver_scale) * next[j];
            }
            const double scale = std::max(
                1.0, std::max(std::abs(value), z * weights[i])
            );
            double residual;
            if (std::abs(next[i]) > 1e-12) {
                residual = std::abs(value + z * weights[i] *
                    (next[i] > 0 ? 1.0 : -1.0)) / scale;
            } else {
                residual = std::max(0.0, std::abs(value) - z * weights[i]) /
                    scale;
            }
            kkt = std::max(kkt, residual);
            const double difference = next[i] - current[i];
            step_squared += difference * difference;
            current_squared += current[i] * current[i];
        }
        current = next;
        accelerated = next_accelerated;
        momentum = updated_momentum;
        if (std::sqrt(step_squared) <= solver_tol *
                (1.0 + std::sqrt(current_squared)) &&
            kkt <= std::max(solver_tol, 1e-10)) {
            converged = true;
            break;
        }
    }
    for (int i = 0; i < n; ++i) delta[active[i]] = current[i];
    objective = 0.0;
    for (int i = 0; i < dimension; ++i) {
        objective -= r[i] * delta[i];
        objective += z * std::sqrt(information[i]) * std::abs(delta[i]);
        for (int j = 0; j < dimension; ++j) {
            objective += 0.5 * delta[i] * H(i, j) * delta[j];
        }
    }
    return List::create(
        _["delta"] = delta,
        _["status"] = converged ? "ok" : "max_iter",
        _["iterations"] = std::min(iteration, solver_max_iter),
        _["kkt_residual"] = kkt, _["objective"] = objective
    );
}
