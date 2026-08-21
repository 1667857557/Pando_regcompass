# Conditional-GRN E-star production core for one frozen exact-edge dictionary.
#
# Production uses one fixed deviation threshold z = 0.25. The shared component
# is the joint MLE; only identifiable condition contrasts are sparsified.
# A full K-1 contrast tree is retained for every exact edge. Tree coordinates
# that are not identifiable remain in the geometry with information zero and are
# fixed to delta = 0 as a deterministic boundary convention. E-star equality is
# retained as production-estimator metadata only; formal edge inference is run
# separately on the frozen dictionary without using the selected fusion graph.

.condition_E_star_z <- 0.25
.condition_E_star_penalty_family <- "information_scaled_sparse_deviation"

.condition_E_star_control <- function(control = list()) {
    if (is.null(control)) control <- list()
    if (!is.list(control)) {
        stop("`condition_e_control` must be a list.", call. = FALSE)
    }
    defaults <- list(
        scale_floor = 1e-8,
        rank_tol = 1e-10,
        solver_tol = 1e-8,
        solver_max_iter = 5000L,
        fusion_tol = 1e-8
    )
    unknown <- setdiff(names(control), names(defaults))
    if (length(unknown)) {
        stop("Unknown `condition_e_control` field(s): ",
             paste(unknown, collapse = ", "), call. = FALSE)
    }
    out <- utils::modifyList(defaults, control)
    for (field in c("scale_floor", "rank_tol", "solver_tol", "fusion_tol")) {
        value <- suppressWarnings(as.numeric(out[[field]]))
        if (length(value) != 1L || !is.finite(value) || value <= 0) {
            stop("`condition_e_control$", field,
                 "` must be one positive finite number.", call. = FALSE)
        }
        out[[field]] <- value
    }
    value <- suppressWarnings(as.numeric(out$solver_max_iter))
    if (length(value) != 1L || !is.finite(value) || value < 1L ||
        value != as.integer(value)) {
        stop("`condition_e_control$solver_max_iter` must be a positive integer.",
             call. = FALSE)
    }
    out$solver_max_iter <- as.integer(value)
    out
}

.condition_symmetric_pinv <- function(x, rank_tol = 1e-10) {
    x <- (as.matrix(x) + t(as.matrix(x))) / 2
    if (!nrow(x)) {
        return(list(
            inverse = x, rank = 0L, values = numeric(), vectors = x,
            keep = logical(), projector = x, tolerance = 0
        ))
    }
    eig <- eigen(x, symmetric = TRUE)
    scale <- max(1, max(abs(eig$values)))
    tolerance <- rank_tol * max(dim(x)) * scale
    keep <- eig$values > tolerance
    inverse <- if (any(keep)) {
        v <- eig$vectors[, keep, drop = FALSE]
        v %*% (t(v) / eig$values[keep])
    } else matrix(0, nrow(x), ncol(x))
    projector <- if (any(keep)) {
        v <- eig$vectors[, keep, drop = FALSE]
        v %*% t(v)
    } else matrix(0, nrow(x), ncol(x))
    list(
        inverse = inverse, rank = sum(keep), values = eig$values,
        vectors = eig$vectors, keep = keep, projector = projector,
        tolerance = tolerance
    )
}

.condition_blockdiag <- function(blocks) {
    if (!length(blocks)) return(matrix(numeric(), 0L, 0L))
    size <- vapply(blocks, nrow, integer(1))
    out <- matrix(0, sum(size), sum(size))
    offset <- 0L
    for (i in seq_along(blocks)) {
        index <- offset + seq_len(size[[i]])
        out[index, index] <- blocks[[i]]
        offset <- offset + size[[i]]
    }
    out
}

.condition_pair_information <- function(Q, p, conditions, rank_tol) {
    k <- length(conditions)
    decomposition <- .condition_symmetric_pinv(Q, rank_tol)
    projector <- decomposition$projector
    inverse <- decomposition$inverse
    pairs <- utils::combn(seq_len(k), 2L, simplify = FALSE)
    rows <- vector("list", p * length(pairs))
    cursor <- 0L
    for (edge in seq_len(p)) {
        for (pair in pairs) {
            a <- pair[[1L]]
            b <- pair[[2L]]
            l <- numeric(k * p)
            l[(a - 1L) * p + edge] <- -1
            l[(b - 1L) * p + edge] <- 1
            residual <- as.numeric(l - l %*% projector)
            estimable <- max(abs(residual)) <= max(
                1e-8, 20 * rank_tol * (1 + max(abs(l)))
            )
            variance <- if (estimable) {
                as.numeric(l %*% inverse %*% l)
            } else NA_real_
            information <- if (estimable && is.finite(variance) &&
                               variance > 0) 1 / variance else 0
            cursor <- cursor + 1L
            rows[[cursor]] <- data.frame(
                edge_index = edge,
                condition_a_index = a,
                condition_b_index = b,
                condition_a = conditions[[a]],
                condition_b = conditions[[b]],
                pair_information = information,
                pair_estimable = information > 0,
                stringsAsFactors = FALSE
            )
        }
    }
    do.call(rbind, rows)
}

.condition_maximum_spanning_tree <- function(vertices, edges, reference_index) {
    vertices <- sort(unique(as.integer(vertices)))
    if (length(vertices) < 2L) return(edges[FALSE, , drop = FALSE])
    edges <- edges[
        edges$condition_a_index %in% vertices &
        edges$condition_b_index %in% vertices, , drop = FALSE
    ]
    if (!nrow(edges)) return(edges)
    reference_edge <- edges$condition_a_index == reference_index |
        edges$condition_b_index == reference_index
    edges <- edges[order(
        !(edges$pair_estimable %in% TRUE),
        !reference_edge,
        -as.numeric(edges$pair_information),
        edges$condition_a_index, edges$condition_b_index
    ), , drop = FALSE]
    parent <- seq_len(max(vertices))
    find_root <- function(x) {
        while (parent[[x]] != x) x <- parent[[x]]
        x
    }
    selected <- logical(nrow(edges))
    for (i in seq_len(nrow(edges))) {
        a <- edges$condition_a_index[[i]]
        b <- edges$condition_b_index[[i]]
        ra <- find_root(a)
        rb <- find_root(b)
        if (ra == rb) next
        selected[[i]] <- TRUE
        parent[[rb]] <- ra
        if (sum(selected) == length(vertices) - 1L) break
    }
    answer <- edges[selected, , drop = FALSE]
    if (nrow(answer) != length(vertices) - 1L) {
        stop("Unable to construct a full condition contrast tree.",
             call. = FALSE)
    }
    answer
}

.condition_identifiable_contrast_tree <- function(
    Q, p, conditions, reference_condition, rank_tol = 1e-10) {
    k <- length(conditions)
    reference_index <- match(reference_condition, conditions)
    if (is.na(reference_index)) {
        stop("`reference_condition` is not one of the fitted conditions.",
             call. = FALSE)
    }
    pair_information <- .condition_pair_information(
        Q, p, conditions, rank_tol
    )
    selected_rows <- vector("list", p)
    boundary <- vector("list", p)
    vertices_used <- vector("list", p)
    for (edge in seq_len(p)) {
        candidate <- pair_information[
            pair_information$edge_index == edge, , drop = FALSE
        ]
        tree <- .condition_maximum_spanning_tree(
            seq_len(k), candidate, reference_index
        )
        tree$contrast_coordinate <- paste0(
            tree$condition_b, "-", tree$condition_a
        )
        tree$tree_edge_index <- seq_len(nrow(tree))
        selected_rows[[edge]] <- tree
        identifiable_vertices <- sort(unique(c(
            tree$condition_a_index[tree$pair_estimable %in% TRUE],
            tree$condition_b_index[tree$pair_estimable %in% TRUE]
        )))
        vertices_used[[edge]] <- identifiable_vertices
        boundary[[edge]] <- sort(unique(c(
            tree$condition_a_index[!(tree$pair_estimable %in% TRUE)],
            tree$condition_b_index[!(tree$pair_estimable %in% TRUE)]
        )))
    }
    selected <- do.call(rbind, selected_rows)
    rownames(selected) <- NULL
    selected$contrast_index <- seq_len(nrow(selected))
    D <- matrix(0, nrow(selected), k * p)
    for (i in seq_len(nrow(selected))) {
        edge <- selected$edge_index[[i]]
        a <- selected$condition_a_index[[i]]
        b <- selected$condition_b_index[[i]]
        D[i, (a - 1L) * p + edge] <- -1
        D[i, (b - 1L) * p + edge] <- 1
    }
    list(
        D = D,
        metadata = selected,
        pair_information = pair_information,
        boundary_conditions = boundary,
        vertices_used = vertices_used,
        reference_condition = reference_condition,
        reference_index = reference_index
    )
}

.condition_q_orthogonal_decomposition <- function(Q, A, D, rank_tol) {
    shared <- .condition_symmetric_pinv(crossprod(A, Q %*% A), rank_tol)
    if (!nrow(D)) {
        return(list(
            R = matrix(0, nrow(Q), 0L),
            B0 = matrix(0, nrow(Q), 0L),
            shared_inverse = shared$inverse,
            dr_error = 0,
            orthogonality_error = 0
        ))
    }
    dd <- tcrossprod(D)
    dd_inverse <- .condition_symmetric_pinv(dd, rank_tol)$inverse
    B0 <- t(D) %*% dd_inverse
    R <- B0 - A %*% shared$inverse %*% crossprod(A, Q %*% B0)
    dr_error <- max(abs(D %*% R - diag(nrow(D))))
    orthogonality_error <- max(abs(crossprod(A, Q %*% R)))
    list(
        R = R,
        B0 = B0,
        shared_inverse = shared$inverse,
        dr_error = dr_error,
        orthogonality_error = orthogonality_error
    )
}

.condition_profile_coordinate_information <- function(H, rank_tol) {
    if (!nrow(H)) {
        return(list(
            information = numeric(), estimable = logical(),
            inverse = H, projector = H
        ))
    }
    decomposition <- .condition_symmetric_pinv(H, rank_tol)
    inverse <- decomposition$inverse
    projector <- decomposition$projector
    information <- numeric(nrow(H))
    estimable <- logical(nrow(H))
    for (j in seq_len(nrow(H))) {
        unit <- numeric(nrow(H)); unit[[j]] <- 1
        residual <- unit - projector[j, ]
        ok <- max(abs(residual)) <= max(1e-8, 20 * rank_tol)
        variance <- inverse[j, j]
        if (ok && is.finite(variance) && variance > 0) {
            estimable[[j]] <- TRUE
            information[[j]] <- 1 / variance
        }
    }
    list(
        information = information,
        estimable = estimable,
        inverse = inverse,
        projector = projector
    )
}

.condition_E_star_kkt <- function(delta, gradient, weights, active, z) {
    if (!length(delta) || !any(active)) return(0)
    value <- numeric(length(delta))
    nonzero <- active & abs(delta) > 1e-12
    zero <- active & !nonzero
    value[nonzero] <- abs(
        gradient[nonzero] + z * weights[nonzero] * sign(delta[nonzero])
    )
    value[zero] <- pmax(0, abs(gradient[zero]) - z * weights[zero])
    max(value[active])
}

.condition_E_star_coordinate_descent <- function(
    H, r, weights, initial, z, tol, max_iter) {
    diagonal <- diag(H)
    diagonal_floor <- .Machine$double.eps *
        max(1, max(abs(diagonal), na.rm = TRUE))
    usable <- is.finite(diagonal) & diagonal > diagonal_floor
    current <- as.numeric(initial)
    current[!is.finite(current) | !usable] <- 0
    gradient <- as.numeric(H %*% current - r)
    kkt <- Inf
    converged <- FALSE

    for (iteration in seq_len(max_iter)) {
        for (j in which(usable)) {
            old <- current[[j]]
            partial <- diagonal[[j]] * old - gradient[[j]]
            updated <- sign(partial) *
                max(abs(partial) - z * weights[[j]], 0) / diagonal[[j]]
            difference <- updated - old
            if (difference != 0) {
                current[[j]] <- updated
                gradient <- gradient + H[, j] * difference
            }
        }
        kkt <- .condition_E_star_kkt(
            current, gradient, weights,
            rep(TRUE, length(current)), z
        )
        if (is.finite(kkt) && kkt <= tol) {
            converged <- TRUE
            break
        }
    }

    objective <- 0.5 * drop(crossprod(current, H %*% current)) -
        drop(crossprod(r, current)) + z * sum(weights * abs(current))
    list(
        delta = current,
        status = if (converged) "ok" else "max_iter",
        iterations = as.integer(iteration),
        kkt_residual = as.numeric(kkt),
        objective = as.numeric(objective)
    )
}

.condition_E_star_fit <- function(
    H, r, information, control = .condition_E_star_control()) {
    control <- .condition_E_star_control(control)
    if (!length(r)) {
        return(list(
            delta = numeric(), status = "ok", iterations = 0L,
            kkt_residual = 0, objective = 0
        ))
    }
    H <- (as.matrix(H) + t(as.matrix(H))) / 2
    r <- as.numeric(r)
    information <- as.numeric(information)
    if (!identical(dim(H), c(length(r), length(r))) ||
        length(information) != length(r) || any(!is.finite(H)) ||
        any(!is.finite(r)) || any(!is.finite(information)) ||
        any(information < 0)) {
        stop("Invalid E-star profile system.", call. = FALSE)
    }
    active <- information > 0
    delta <- numeric(length(r))
    if (!any(active)) {
        return(list(
            delta = delta, status = "ok", iterations = 0L,
            kkt_residual = 0, objective = 0
        ))
    }
    index <- which(active)
    Hs <- H[index, index, drop = FALSE]
    rs <- r[index]
    weights <- sqrt(information[index])
    lipschitz <- max(eigen(
        (Hs + t(Hs)) / 2, symmetric = TRUE, only.values = TRUE
    )$values)
    if (!is.finite(lipschitz) || lipschitz <= 0) {
        return(list(
            delta = delta, status = "invalid_profile_information",
            iterations = 0L, kkt_residual = Inf, objective = Inf
        ))
    }
    current <- numeric(length(index))
    accelerated <- current
    momentum <- 1
    converged <- FALSE
    kkt <- Inf
    objective <- Inf
    for (iteration in seq_len(control$solver_max_iter)) {
        gradient <- as.numeric(Hs %*% accelerated - rs)
        trial <- accelerated - gradient / lipschitz
        threshold <- .condition_E_star_z * weights / lipschitz
        next_value <- sign(trial) * pmax(abs(trial) - threshold, 0)
        next_momentum <- (1 + sqrt(1 + 4 * momentum^2)) / 2
        next_accelerated <- next_value +
            ((momentum - 1) / next_momentum) * (next_value - current)
        current_gradient <- as.numeric(Hs %*% next_value - rs)
        kkt <- .condition_E_star_kkt(
            next_value, current_gradient, weights,
            rep(TRUE, length(index)), .condition_E_star_z
        )
        full_delta <- numeric(length(r)); full_delta[index] <- next_value
        objective <- 0.5 * drop(crossprod(full_delta, H %*% full_delta)) -
            drop(crossprod(r, full_delta)) +
            .condition_E_star_z * sum(sqrt(information) * abs(full_delta))
        step <- sqrt(sum((next_value - current)^2))
        scale <- 1 + sqrt(sum(current^2))
        current <- next_value
        accelerated <- next_accelerated
        momentum <- next_momentum
        if (step <= control$solver_tol * scale &&
            kkt <= max(control$solver_tol, 1e-10)) {
            converged <- TRUE
            break
        }
    }
    if (!converged) {
        fallback_max_iter <- as.integer(min(
            .Machine$integer.max,
            max(50000, 10 * as.numeric(control$solver_max_iter))
        ))
        fallback <- .condition_E_star_coordinate_descent(
            H = Hs, r = rs, weights = weights, initial = current,
            z = .condition_E_star_z,
            tol = max(control$solver_tol, 1e-10),
            max_iter = fallback_max_iter
        )
        current <- fallback$delta
        kkt <- fallback$kkt_residual
        objective <- fallback$objective
        iteration <- as.integer(iteration + fallback$iterations)
        converged <- identical(fallback$status, "ok")
    }
    delta[index] <- current
    list(
        delta = delta,
        status = if (converged) "ok" else "max_iter",
        iterations = as.integer(iteration),
        kkt_residual = as.numeric(kkt),
        objective = as.numeric(objective)
    )
}

.condition_union_find_components <- function(k, links) {
    parent <- seq_len(k)
    root <- function(x) {
        while (parent[[x]] != x) x <- parent[[x]]
        x
    }
    if (length(links)) {
        for (link in links) {
            a <- as.integer(link[[1L]])
            b <- as.integer(link[[2L]])
            ra <- root(a); rb <- root(b)
            if (ra != rb) parent[[rb]] <- ra
        }
    }
    roots <- vapply(seq_len(k), root, integer(1))
    match(roots, unique(roots))
}

.condition_fusion_components <- function(tree, delta, p, conditions, control) {
    k <- length(conditions)
    component <- matrix(NA_integer_, k, p)
    boundary_condition <- matrix(FALSE, k, p)
    fused_by_penalty <- logical(p)
    shared_by_boundary <- logical(p)
    tree_metadata <- tree$metadata
    tree_metadata$delta <- as.numeric(delta)
    if (!"contrast_identifiable" %in% colnames(tree_metadata)) {
        tree_metadata$contrast_identifiable <- FALSE
    }
    tree_metadata$shared_by_boundary <-
        !(tree_metadata$contrast_identifiable %in% TRUE)
    tree_metadata$fused_by_penalty <-
        tree_metadata$contrast_identifiable %in% TRUE &
        abs(tree_metadata$delta) <= control$fusion_tol
    tree_metadata$equality_selected <-
        tree_metadata$shared_by_boundary | tree_metadata$fused_by_penalty
    for (edge in seq_len(p)) {
        links <- list()
        index <- which(tree_metadata$edge_index == edge)
        for (j in index) {
            if (!tree_metadata$equality_selected[[j]]) next
            a <- tree_metadata$condition_a_index[[j]]
            b <- tree_metadata$condition_b_index[[j]]
            links[[length(links) + 1L]] <- c(a, b)
            if (tree_metadata$shared_by_boundary[[j]]) {
                shared_by_boundary[[edge]] <- TRUE
                boundary_condition[c(a, b), edge] <- TRUE
            } else if (tree_metadata$fused_by_penalty[[j]]) {
                fused_by_penalty[[edge]] <- TRUE
            }
        }
        component[, edge] <- .condition_union_find_components(k, links)
    }
    shared_edge <- vapply(seq_len(p), function(edge) {
        length(unique(component[, edge])) == 1L
    }, logical(1))
    list(
        component = component,
        boundary_condition = boundary_condition,
        fused_by_penalty = fused_by_penalty,
        shared_by_boundary = shared_by_boundary,
        shared_edge = shared_edge,
        contrast_tree = tree_metadata
    )
}

.condition_scheme_e_fit <- function(
    x, y, scaling, min_residual_df = 1L, control = list(),
    reference_condition = NULL) {
    control <- .condition_E_star_control(control)
    conditions <- names(x)
    k <- length(conditions)
    p_full <- ncol(x[[1L]])
    if (k < 2L || !identical(names(y), conditions) ||
        any(vapply(seq_along(x), function(i) {
            !is.matrix(x[[i]]) || nrow(x[[i]]) != length(y[[i]]) ||
                !identical(colnames(x[[i]]), colnames(x[[1L]])) ||
                any(!is.finite(x[[i]])) || any(!is.finite(y[[i]]))
        }, logical(1)))) {
        return(list(status = "nonfinite_or_misaligned_input"))
    }
    if (is.null(reference_condition)) reference_condition <- conditions[[1L]]
    reference_condition <- as.character(reference_condition)
    if (length(reference_condition) != 1L ||
        !reference_condition %in% conditions) {
        return(list(status = "invalid_reference_condition"))
    }

    informative <- as.logical(scaling$informative)
    beta <- beta_z <- deviation_z <- matrix(
        0, k, p_full, dimnames = list(conditions, colnames(x[[1L]]))
    )
    shared_z <- numeric(p_full); names(shared_z) <- colnames(beta)
    raw_information <- matrix(0, k, p_full, dimnames = dimnames(beta))
    zero_variance <- matrix(FALSE, k, p_full, dimnames = dimnames(beta))
    for (i in seq_len(k)) {
        value <- apply(x[[i]], 2L, stats::var)
        zero_variance[i, ] <- !is.finite(value) | value <= scaling$floor^2
    }

    if (!any(informative)) {
        intercept <- vapply(y, mean, numeric(1))
        prediction <- Map(function(value, a) rep(a, length(value)), y, intercept)
        rsq <- vapply(y, function(value) {
            if (stats::var(value) > 0) 0 else NA_real_
        }, numeric(1))
        return(list(
            status = "ok", beta = beta, beta_z = beta_z,
            shared_z = shared_z, deviation_z = deviation_z,
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
            shared_edge = rep(TRUE, p_full),
            fusion_component_id = matrix(
                "component1", k, p_full, dimnames = dimnames(beta)
            ),
            boundary_condition = matrix(
                TRUE, k, p_full, dimnames = dimnames(beta)
            ),
            contrast_tree = data.frame(), solver_status = "ok",
            kkt_residual = 0, iterations = 0L, objective = 0,
            penalty_family = .condition_E_star_penalty_family,
            penalty_value = .condition_E_star_z,
            reference_condition = reference_condition,
            condition_weight = stats::setNames(rep(1, k), conditions),
            orthogonality_error = 0, dr_error = 0
        ))
    }

    keep <- which(informative)
    p <- length(keep)
    gram_raw <- rhs_raw <- xc <- yc <- vector("list", k)
    raw_rank <- integer(k); raw_kappa <- numeric(k); sse_raw <- numeric(k)
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
    }
    names(raw_rank) <- names(raw_kappa) <- conditions
    residual_df <- sum(lengths(y) - 1L - raw_rank)
    if (!is.finite(residual_df) || residual_df < min_residual_df) {
        return(list(status = "insufficient_df"))
    }
    pooled_y_scale <- stats::var(unlist(y, use.names = FALSE))
    variance_floor <- .Machine$double.eps * max(1, pooled_y_scale, na.rm = TRUE)
    sigma2 <- max(sum(sse_raw) / residual_df, variance_floor)
    q <- lapply(gram_raw, function(one) one / sigma2)
    h_condition <- lapply(rhs_raw, function(one) one / sigma2)
    for (i in seq_len(k)) {
        raw_information[i, keep] <- diag(q[[i]]) * scaling$scale[keep]^2
    }

    Q <- .condition_blockdiag(q)
    h <- unlist(h_condition, use.names = FALSE)
    A <- kronecker(matrix(1, nrow = k, ncol = 1), diag(p))
    tree <- .condition_identifiable_contrast_tree(
        Q = Q, p = p, conditions = conditions,
        reference_condition = reference_condition, rank_tol = control$rank_tol
    )
    geometry <- .condition_q_orthogonal_decomposition(
        Q, A, tree$D, control$rank_tol
    )
    shared_score <- as.numeric(crossprod(A, h))
    mu <- as.numeric(geometry$shared_inverse %*% shared_score)
    H <- crossprod(geometry$R, Q %*% geometry$R)
    H <- (H + t(H)) / 2
    r <- as.numeric(crossprod(geometry$R, h))
    coordinate_information <- .condition_profile_coordinate_information(
        H, control$rank_tol
    )
    solved <- .condition_E_star_fit(
        H, r, coordinate_information$information, control
    )
    if (!identical(solved$status, "ok")) {
        return(list(
            status = "solver_failed", solver_status = solved$status,
            kkt_residual = solved$kkt_residual, iterations = solved$iterations
        ))
    }
    delta <- solved$delta
    beta_vector <- as.numeric(A %*% mu + geometry$R %*% delta)
    tree$metadata$profile_information <- coordinate_information$information
    tree$metadata$contrast_identifiable <- coordinate_information$estimable
    tree$metadata$shared_by_boundary <- !coordinate_information$estimable
    tree$metadata$delta_standardized <- delta
    tree$metadata$delta <- delta /
        scaling$scale[keep[tree$metadata$edge_index]]

    beta_z_small <- matrix(beta_vector, nrow = k, byrow = TRUE)
    beta_small <- sweep(beta_z_small, 2L, scaling$scale[keep], "/")
    beta_z[, keep] <- beta_z_small
    beta[, keep] <- beta_small
    shared_z[keep] <- mu
    deviation_z[, keep] <- sweep(beta_z_small, 2L, mu, "-")

    fusion <- .condition_fusion_components(tree, delta, p, conditions, control)
    fusion_component_id <- matrix(
        NA_character_, k, p_full, dimnames = dimnames(beta)
    )
    boundary_condition <- matrix(TRUE, k, p_full, dimnames = dimnames(beta))
    for (local_edge in seq_len(p)) {
        full_edge <- keep[[local_edge]]
        fusion_component_id[, full_edge] <- paste0(
            "component", fusion$component[, local_edge]
        )
        boundary_condition[, full_edge] <-
            fusion$boundary_condition[, local_edge]
    }
    profile_information <- rep(0, p_full)
    contrast_identifiable <- rep(FALSE, p_full)
    for (local_edge in seq_len(p)) {
        index <- fusion$contrast_tree$edge_index == local_edge
        value <- fusion$contrast_tree$profile_information[index]
        identifiable_here <-
            fusion$contrast_tree$contrast_identifiable[index] %in% TRUE
        if (any(identifiable_here)) {
            full_edge <- keep[[local_edge]]
            raw_info <- value[identifiable_here] * scaling$scale[[full_edge]]^2
            raw_info <- raw_info[is.finite(raw_info) & raw_info > 0]
            if (length(raw_info)) {
                profile_information[[full_edge]] <- min(raw_info)
            }
            contrast_identifiable[[full_edge]] <- TRUE
        }
    }
    shared_by_boundary <- rep(TRUE, p_full)
    fused_by_penalty <- rep(FALSE, p_full)
    shared_edge <- rep(TRUE, p_full)
    shared_by_boundary[keep] <- fusion$shared_by_boundary
    fused_by_penalty[keep] <- fusion$fused_by_penalty
    shared_edge[keep] <- fusion$shared_edge

    intercept <- numeric(k); prediction <- vector("list", k); rsq <- numeric(k)
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

    list(
        status = "ok", beta = beta, beta_z = beta_z,
        shared_z = shared_z, deviation_z = deviation_z,
        intercept = intercept, prediction = prediction, rsq = rsq,
        raw_rank = raw_rank, residual_df = residual_df,
        raw_kappa = raw_kappa, zero_variance = zero_variance,
        informative = informative, informative_index = keep,
        sigma2 = sigma2, raw_information = raw_information,
        profile_information = profile_information,
        contrast_identifiable = contrast_identifiable,
        shared_by_boundary = shared_by_boundary,
        fused_by_penalty = fused_by_penalty,
        shared_edge = shared_edge,
        fusion_component_id = fusion_component_id,
        boundary_condition = boundary_condition,
        contrast_tree = fusion$contrast_tree,
        solver_status = solved$status, kkt_residual = solved$kkt_residual,
        iterations = solved$iterations, objective = solved$objective,
        penalty_family = .condition_E_star_penalty_family,
        penalty_value = .condition_E_star_z,
        reference_condition = reference_condition,
        condition_weight = stats::setNames(rep(1, k), conditions),
        orthogonality_error = geometry$orthogonality_error,
        dr_error = geometry$dr_error
    )
}

