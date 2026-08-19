# Native conditional-GRN Scheme E: exact-edge sparse condition deviations.
#
# The fitted dictionary is fixed before this solver is called. Predictors are
# centered within condition and share one cross-condition scale. The Gaussian
# information is therefore allowed to grow with the number of paired cells in
# each condition; no n_total/(K*n_c) reweighting is applied. A single pooled
# residual variance puts every condition on the same information scale.
#
# For K conditions, beta_c = mu + C_c gamma where C is an orthonormal basis of
# the centered condition subspace. Scheme E is extended without a condition
# graph by applying one group penalty per exact edge to its full (K-1)-vector of
# condition contrasts. For K=2 this is exactly
#   sign(d) * (abs(d) - z / sqrt(I_delta))_+
# with z = 0.25.

.condition_scheme_e_z <- 0.25
.condition_scheme_e_penalty_family <- "exact_edge_sparse_deviation"

.condition_scheme_e_control <- function(control = list()) {
    if (is.null(control)) control <- list()
    if (!is.list(control)) {
        stop("`ridge_control` must be a list.", call. = FALSE)
    }
    defaults <- list(
        scale_floor = 1e-8,
        rank_tol = 1e-10,
        solver_tol = 1e-8,
        solver_max_iter = 5000L,
        fusion_tol = 1e-8
    )
    # Accept the previous ridge-CV fields only as inert compatibility inputs.
    # They are not used to select or tune Scheme E.
    legacy <- c("lambda_grid", "lambda_rule", "cv_folds", "seed")
    unknown <- setdiff(names(control), c(names(defaults), legacy))
    if (length(unknown)) {
        stop("Unknown `ridge_control` field(s): ",
             paste(unknown, collapse = ", "), call. = FALSE)
    }
    control[intersect(names(control), legacy)] <- NULL
    out <- utils::modifyList(defaults, control)
    scalar_positive <- c("scale_floor", "rank_tol", "solver_tol", "fusion_tol")
    for (field in scalar_positive) {
        value <- suppressWarnings(as.numeric(out[[field]]))
        if (length(value) != 1L || !is.finite(value) || value <= 0) {
            stop("`ridge_control$", field,
                 "` must be one positive finite number.", call. = FALSE)
        }
        out[[field]] <- value
    }
    value <- suppressWarnings(as.numeric(out$solver_max_iter))
    if (length(value) != 1L || !is.finite(value) || value < 1L ||
        value != as.integer(value)) {
        stop("`ridge_control$solver_max_iter` must be a positive integer.",
             call. = FALSE)
    }
    out$solver_max_iter <- as.integer(value)
    out$scheme_e_z <- .condition_scheme_e_z
    out
}

.condition_symmetric_eigen <- function(x) {
    x <- (as.matrix(x) + t(as.matrix(x))) / 2
    eigen(x, symmetric = TRUE)
}

.condition_symmetric_pinv <- function(x, rank_tol = 1e-10) {
    x <- (as.matrix(x) + t(as.matrix(x))) / 2
    if (!nrow(x)) {
        return(list(
            inverse = x, rank = 0L, values = numeric(), vectors = x,
            keep = logical(), tolerance = 0
        ))
    }
    eig <- .condition_symmetric_eigen(x)
    scale <- max(1, max(abs(eig$values)))
    tolerance <- rank_tol * max(dim(x)) * scale
    keep <- eig$values > tolerance
    inverse <- if (any(keep)) {
        eig$vectors[, keep, drop = FALSE] %*%
            (t(eig$vectors[, keep, drop = FALSE]) /
                 eig$values[keep])
    } else matrix(0, nrow(x), ncol(x))
    list(
        inverse = inverse,
        rank = sum(keep),
        values = eig$values,
        vectors = eig$vectors,
        keep = keep,
        tolerance = tolerance
    )
}

.condition_contrast_basis <- function(k) {
    k <- as.integer(k)
    if (k < 2L) return(matrix(numeric(), nrow = k, ncol = 0L))
    basis <- stats::contr.helmert(k)
    sweep(basis, 2L, sqrt(colSums(basis^2)), "/")
}

.condition_scheme_e_identifiable_edges <- function(info, p, k_minus_one,
                                                     rank_tol = 1e-10) {
    if (!nrow(info) || p < 1L || k_minus_one < 1L) {
        return(list(
            edge = rep(FALSE, p), projector = matrix(0, nrow(info), ncol(info)),
            rank = 0L, tolerance = 0
        ))
    }
    decomposition <- .condition_symmetric_pinv(info, rank_tol)
    projector <- if (any(decomposition$keep)) {
        v <- decomposition$vectors[, decomposition$keep, drop = FALSE]
        v %*% t(v)
    } else matrix(0, nrow(info), ncol(info))
    edge <- vapply(seq_len(p), function(e) {
        index <- (seq_len(k_minus_one) - 1L) * p + e
        block <- projector[index, index, drop = FALSE]
        max(abs(block - diag(k_minus_one))) <=
            max(1e-7, 10 * sqrt(.Machine$double.eps))
    }, logical(1))
    list(
        edge = edge,
        projector = projector,
        rank = decomposition$rank,
        tolerance = decomposition$tolerance
    )
}

.condition_scheme_e_block_inverse_roots <- function(info, p, k_minus_one,
                                                     edge_keep, rank_tol) {
    retained_edges <- which(edge_keep)
    if (!length(retained_edges)) {
        return(list(
            coordinate = integer(), inverse_root = matrix(0, 0, 0),
            groups = list(), information = list()
        ))
    }
    coordinate <- unlist(lapply(retained_edges, function(e) {
        (seq_len(k_minus_one) - 1L) * p + e
    }), use.names = FALSE)
    sub <- info[coordinate, coordinate, drop = FALSE]
    q <- k_minus_one
    inverse_root <- matrix(0, nrow(sub), ncol(sub))
    groups <- vector("list", length(retained_edges))
    information <- vector("list", length(retained_edges))
    for (g in seq_along(retained_edges)) {
        index <- ((g - 1L) * q + 1L):(g * q)
        groups[[g]] <- index
        block <- (sub[index, index, drop = FALSE] +
                      t(sub[index, index, drop = FALSE])) / 2
        eig <- eigen(block, symmetric = TRUE)
        tolerance <- rank_tol * max(1, q) * max(1, max(abs(eig$values)))
        if (any(eig$values <= tolerance)) {
            stop("Identifiable Scheme-E edge has singular profile information.",
                 call. = FALSE)
        }
        inverse_root[index, index] <- eig$vectors %*%
            diag(1 / sqrt(eig$values), nrow = q) %*% t(eig$vectors)
        information[[g]] <- block
    }
    names(groups) <- names(information) <- as.character(retained_edges)
    list(
        coordinate = coordinate,
        inverse_root = inverse_root,
        groups = groups,
        information = information,
        retained_edges = retained_edges,
        sub_information = sub
    )
}

.condition_scheme_e_group_kkt <- function(eta, gradient, groups, z) {
    if (!length(groups)) return(0)
    max(vapply(groups, function(index) {
        value <- eta[index]
        size <- sqrt(sum(value^2))
        if (size > 1e-12) {
            sqrt(sum((gradient[index] + z * value / size)^2))
        } else {
            max(0, sqrt(sum(gradient[index]^2)) - z)
        }
    }, numeric(1)))
}

.condition_scheme_e_fista <- function(info, rhs, block, control) {
    if (!length(block$coordinate)) {
        return(list(
            gamma = numeric(nrow(info)), status = "ok", iterations = 0L,
            kkt_residual = 0, objective = 0, fused = logical()
        ))
    }
    coordinate <- block$coordinate
    sub_info <- block$sub_information
    sub_rhs <- rhs[coordinate]
    winv <- block$inverse_root
    hessian <- crossprod(winv, sub_info %*% winv)
    hessian <- (hessian + t(hessian)) / 2
    transformed_rhs <- as.numeric(crossprod(winv, sub_rhs))
    lipschitz <- max(eigen(hessian, symmetric = TRUE, only.values = TRUE)$values)
    if (!is.finite(lipschitz) || lipschitz <= 0) {
        stop("Scheme-E transformed profile Hessian is not positive definite.",
             call. = FALSE)
    }

    eta <- numeric(length(coordinate))
    accelerated <- eta
    momentum <- 1
    z <- .condition_scheme_e_z
    converged <- FALSE
    kkt <- Inf
    objective <- Inf
    for (iteration in seq_len(control$solver_max_iter)) {
        gradient <- as.numeric(hessian %*% accelerated - transformed_rhs)
        trial <- accelerated - gradient / lipschitz
        next_eta <- numeric(length(trial))
        for (index in block$groups) {
            size <- sqrt(sum(trial[index]^2))
            threshold <- z / lipschitz
            if (is.finite(size) && size > threshold) {
                next_eta[index] <- (1 - threshold / size) * trial[index]
            }
        }
        next_momentum <- (1 + sqrt(1 + 4 * momentum^2)) / 2
        next_accelerated <- next_eta +
            ((momentum - 1) / next_momentum) * (next_eta - eta)
        current_gradient <- as.numeric(hessian %*% next_eta - transformed_rhs)
        kkt <- .condition_scheme_e_group_kkt(
            next_eta, current_gradient, block$groups, z
        )
        objective <- 0.5 * drop(crossprod(next_eta, hessian %*% next_eta)) -
            drop(crossprod(transformed_rhs, next_eta)) +
            z * sum(vapply(block$groups, function(index) {
                sqrt(sum(next_eta[index]^2))
            }, numeric(1)))
        step <- sqrt(sum((next_eta - eta)^2))
        scale <- 1 + sqrt(sum(eta^2))
        eta <- next_eta
        accelerated <- next_accelerated
        momentum <- next_momentum
        if (step <= control$solver_tol * scale &&
            kkt <= max(control$solver_tol, 1e-10)) {
            converged <- TRUE
            break
        }
    }
    gamma <- numeric(nrow(info))
    gamma[coordinate] <- as.numeric(winv %*% eta)
    fused <- vapply(block$groups, function(index) {
        sqrt(sum(gamma[coordinate[index]]^2)) <= control$fusion_tol
    }, logical(1))
    list(
        gamma = gamma,
        status = if (converged) "ok" else "max_iter",
        iterations = as.integer(iteration),
        kkt_residual = as.numeric(kkt),
        objective = as.numeric(objective),
        fused = fused
    )
}

.condition_scheme_e_fit <- function(x, y, scaling, min_residual_df = 1L,
                                    inference = TRUE, control = list()) {
    control <- .condition_scheme_e_control(control)
    conditions <- names(x)
    k <- length(conditions)
    p_full <- ncol(x[[1L]])
    if (k < 2L || !identical(names(y), conditions) ||
        any(vapply(seq_along(x), function(i) {
            !is.matrix(x[[i]]) || nrow(x[[i]]) != length(y[[i]]) ||
                any(!is.finite(x[[i]])) || any(!is.finite(y[[i]]))
        }, logical(1)))) {
        return(list(status = "nonfinite_or_misaligned_input"))
    }

    informative <- as.logical(scaling$informative)
    beta <- matrix(0, k, p_full,
                   dimnames = list(conditions, colnames(x[[1L]])))
    beta_z <- beta
    shared_z <- numeric(p_full)
    deviation_z <- beta
    se <- statistic <- pval <- matrix(
        NA_real_, k, p_full, dimnames = dimnames(beta)
    )
    raw_information <- matrix(0, k, p_full, dimnames = dimnames(beta))
    zero_variance_values <- vapply(seq_along(x), function(i) {
        value <- apply(x[[i]], 2L, stats::var)
        !is.finite(value) | value <= scaling$floor^2
    }, logical(p_full))
    zero_variance <- matrix(
        as.logical(zero_variance_values), nrow = k, ncol = p_full,
        byrow = TRUE, dimnames = dimnames(beta)
    )

    if (!any(informative)) {
        intercept <- vapply(y, mean, numeric(1))
        prediction <- Map(function(value, a) rep(a, length(value)),
                          y, intercept)
        rsq <- vapply(y, function(value) {
            if (stats::var(value) > 0) 0 else NA_real_
        }, numeric(1))
        return(list(
            status = "ok", beta = beta, beta_z = beta_z,
            shared_z = shared_z, deviation_z = deviation_z,
            se = se, statistic = statistic, pval = pval,
            intercept = intercept, prediction = prediction, rsq = rsq,
            raw_rank = stats::setNames(rep(0L, k), conditions),
            residual_df = sum(lengths(y)) - k,
            raw_kappa = stats::setNames(rep(NA_real_, k), conditions),
            zero_variance = zero_variance, informative = informative,
            informative_index = integer(), sigma2 = NA_real_,
            raw_information = raw_information,
            profile_information = rep(0, p_full),
            contrast_identifiable = rep(FALSE, p_full),
            shared_by_boundary = rep(TRUE, p_full),
            fused_by_penalty = rep(FALSE, p_full),
            solver_status = "ok", kkt_residual = 0, iterations = 0L,
            penalty_family = .condition_scheme_e_penalty_family,
            penalty_value = .condition_scheme_e_z,
            condition_weight = stats::setNames(rep(1, k), conditions),
            diagnostic_covariance_z = vector("list", k)
        ))
    }

    keep <- which(informative)
    p <- length(keep)
    xc <- yc <- vector("list", k)
    raw_rank <- integer(k)
    raw_kappa <- numeric(k)
    gram_raw <- rhs_raw <- vector("list", k)
    sse_raw <- numeric(k)
    diagnostic_covariance_z <- vector("list", k)
    for (i in seq_len(k)) {
        zmat <- sweep(x[[i]][, keep, drop = FALSE], 2L,
                      scaling$center[keep], "-")
        zmat <- sweep(zmat, 2L, scaling$scale[keep], "/")
        xc[[i]] <- sweep(zmat, 2L, colMeans(zmat), "-")
        yc[[i]] <- as.numeric(y[[i]]) - mean(y[[i]])
        gram_raw[[i]] <- crossprod(xc[[i]])
        rhs_raw[[i]] <- as.numeric(crossprod(xc[[i]], yc[[i]]))
        decomposition <- .condition_symmetric_pinv(
            gram_raw[[i]], control$rank_tol
        )
        raw_rank[[i]] <- decomposition$rank
        raw_kappa[[i]] <- .condition_ridge_kappa(cbind(1, zmat))
        ols <- as.numeric(decomposition$inverse %*% rhs_raw[[i]])
        residual <- yc[[i]] - as.numeric(xc[[i]] %*% ols)
        sse_raw[[i]] <- sum(residual^2)
        diagnostic_covariance_z[[i]] <- decomposition$inverse
    }
    names(raw_rank) <- names(raw_kappa) <- conditions
    names(diagnostic_covariance_z) <- conditions
    residual_df <- sum(lengths(y) - 1L - raw_rank)
    if (!is.finite(residual_df) || residual_df < min_residual_df) {
        return(list(status = "insufficient_df"))
    }
    pooled_y_scale <- stats::var(unlist(y, use.names = FALSE))
    variance_floor <- .Machine$double.eps * max(1, pooled_y_scale, na.rm = TRUE)
    sigma2 <- max(sum(sse_raw) / residual_df, variance_floor)

    q <- lapply(gram_raw, function(one) one / sigma2)
    r <- lapply(rhs_raw, function(one) one / sigma2)
    for (i in seq_len(k)) {
        raw_information[i, keep] <- diag(q[[i]]) * scaling$scale[keep]^2
        diagnostic_covariance_z[[i]] <-
            diagnostic_covariance_z[[i]] * sigma2
    }

    contrast_basis <- .condition_contrast_basis(k)
    km1 <- k - 1L
    m <- km1 * p
    h_mu <- Reduce(`+`, r)
    h_mg <- matrix(0, p, m)
    h_gg <- matrix(0, m, m)
    rhs_gamma <- numeric(m)
    h_mu_mu <- Reduce(`+`, q)
    for (condition_index in seq_len(k)) {
        for (j in seq_len(km1)) {
            j_index <- ((j - 1L) * p + 1L):(j * p)
            h_mg[, j_index] <- h_mg[, j_index] +
                contrast_basis[condition_index, j] * q[[condition_index]]
            rhs_gamma[j_index] <- rhs_gamma[j_index] +
                contrast_basis[condition_index, j] * r[[condition_index]]
            for (ell in seq_len(km1)) {
                ell_index <- ((ell - 1L) * p + 1L):(ell * p)
                h_gg[j_index, ell_index] <- h_gg[j_index, ell_index] +
                    contrast_basis[condition_index, j] *
                    contrast_basis[condition_index, ell] * q[[condition_index]]
            }
        }
    }
    shared_inverse <- .condition_symmetric_pinv(
        h_mu_mu, control$rank_tol
    )$inverse
    profile_information <- h_gg -
        crossprod(h_mg, shared_inverse %*% h_mg)
    profile_information <-
        (profile_information + t(profile_information)) / 2
    profile_rhs <- rhs_gamma -
        as.numeric(crossprod(h_mg, shared_inverse %*% h_mu))

    identifiability <- .condition_scheme_e_identifiable_edges(
        profile_information, p, km1, control$rank_tol
    )
    blocks <- .condition_scheme_e_block_inverse_roots(
        profile_information, p, km1, identifiability$edge,
        control$rank_tol
    )
    solved <- .condition_scheme_e_fista(
        profile_information, profile_rhs, blocks, control
    )
    if (!identical(solved$status, "ok")) {
        return(list(
            status = "solver_failed", solver_status = solved$status,
            kkt_residual = solved$kkt_residual,
            iterations = solved$iterations
        ))
    }
    gamma <- solved$gamma
    mu <- as.numeric(shared_inverse %*% (h_mu - h_mg %*% gamma))
    beta_z_small <- matrix(0, k, p)
    for (i in seq_len(k)) {
        beta_z_small[i, ] <- mu
        for (j in seq_len(km1)) {
            index <- ((j - 1L) * p + 1L):(j * p)
            beta_z_small[i, ] <- beta_z_small[i, ] +
                contrast_basis[i, j] * gamma[index]
        }
    }
    beta_small <- sweep(beta_z_small, 2L, scaling$scale[keep], "/")
    beta_z[, keep] <- beta_z_small
    beta[, keep] <- beta_small
    shared_z[keep] <- mu
    deviation_z[, keep] <- sweep(beta_z_small, 2L, mu, "-")

    edge_profile_information <- rep(0, p_full)
    for (e in seq_len(p)) {
        index <- (seq_len(km1) - 1L) * p + e
        block <- profile_information[index, index, drop = FALSE]
        eig <- eigen((block + t(block)) / 2, symmetric = TRUE,
                     only.values = TRUE)$values
        # Convert from standardized beta coordinates to raw coefficient units.
        edge_profile_information[keep[[e]]] <-
            max(0, min(eig)) * scaling$scale[keep[[e]]]^2
    }
    contrast_identifiable <- rep(FALSE, p_full)
    contrast_identifiable[keep] <- identifiability$edge
    shared_by_boundary <- !contrast_identifiable
    fused_by_penalty <- rep(FALSE, p_full)
    if (length(blocks$retained_edges)) {
        fused_by_penalty[keep[blocks$retained_edges]] <- solved$fused
    }

    intercept <- numeric(k)
    prediction <- vector("list", k)
    rsq <- numeric(k)
    for (i in seq_len(k)) {
        intercept[[i]] <- mean(y[[i]]) -
            sum(colMeans(x[[i]][, keep, drop = FALSE]) * beta_small[i, ])
        prediction[[i]] <- as.numeric(
            intercept[[i]] + x[[i]][, keep, drop = FALSE] %*% beta_small[i, ]
        )
        residual <- as.numeric(y[[i]]) - prediction[[i]]
        tss <- sum((as.numeric(y[[i]]) - mean(y[[i]]))^2)
        rsq[[i]] <- if (tss > 0) 1 - sum(residual^2) / tss else NA_real_
    }
    names(intercept) <- names(prediction) <- names(rsq) <- conditions

    if (isTRUE(inference)) {
        for (i in seq_len(k)) {
            covariance <- diagnostic_covariance_z[[i]]
            decomposition <- .condition_symmetric_pinv(
                gram_raw[[i]], control$rank_tol
            )
            projector <- if (any(decomposition$keep)) {
                v <- decomposition$vectors[, decomposition$keep, drop = FALSE]
                v %*% t(v)
            } else matrix(0, p, p)
            estimable_here <- vapply(seq_len(p), function(e) {
                max(abs(projector[e, e, drop = FALSE] - 1)) <= 1e-7
            }, logical(1))
            se_z <- sqrt(pmax(0, diag(covariance)))
            se_raw <- se_z / scaling$scale[keep]
            se[i, keep[estimable_here]] <- se_raw[estimable_here]
            zstat <- beta_small[i, estimable_here] / se_raw[estimable_here]
            zstat[!is.finite(zstat)] <- NA_real_
            statistic[i, keep[estimable_here]] <- zstat
            pval[i, keep[estimable_here]] <- 2 * stats::pnorm(-abs(zstat))
        }
    }

    list(
        status = "ok", beta = beta, beta_z = beta_z,
        shared_z = shared_z, deviation_z = deviation_z,
        se = se, statistic = statistic, pval = pval,
        intercept = intercept, prediction = prediction, rsq = rsq,
        raw_rank = raw_rank, residual_df = residual_df,
        raw_kappa = raw_kappa, zero_variance = zero_variance,
        informative = informative, informative_index = keep,
        sigma2 = sigma2, raw_information = raw_information,
        profile_information = edge_profile_information,
        contrast_identifiable = contrast_identifiable,
        shared_by_boundary = shared_by_boundary,
        fused_by_penalty = fused_by_penalty,
        solver_status = solved$status,
        kkt_residual = solved$kkt_residual,
        iterations = solved$iterations,
        objective = solved$objective,
        penalty_family = .condition_scheme_e_penalty_family,
        penalty_value = .condition_scheme_e_z,
        condition_weight = stats::setNames(rep(1, k), conditions),
        diagnostic_covariance_z = diagnostic_covariance_z,
        contrast_basis = contrast_basis,
        rank_profile = identifiability$rank
    )
}
