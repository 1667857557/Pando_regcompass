# Conditional-GRN E-star/JSE core for one frozen exact-edge dictionary.
#
# Production uses one fixed deviation threshold z = 0.25. The shared component
# is the joint MLE; only identifiable condition contrasts are sparsified.
# Equality selected by E-star is then used to build a reduced joint model for
# covariance, Wald inference and downstream condition-by-target BH adjustment.

.condition_E_star_z <- 0.25
.condition_E_star_penalty_family <- "information_scaled_sparse_deviation"
.condition_E_star_inference_schema <- "scheme_e_fusion_component_joint_refit_v1"

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

.condition_graph_component <- function(vertices, edges) {
    vertices <- as.integer(vertices)
    if (!length(vertices)) return(list())
    adjacency <- stats::setNames(vector("list", length(vertices)), vertices)
    if (nrow(edges)) {
        for (i in seq_len(nrow(edges))) {
            a <- as.character(edges$condition_a_index[[i]])
            b <- as.character(edges$condition_b_index[[i]])
            adjacency[[a]] <- unique(c(adjacency[[a]], as.integer(b)))
            adjacency[[b]] <- unique(c(adjacency[[b]], as.integer(a)))
        }
    }
    seen <- integer()
    components <- list()
    for (root in vertices) {
        if (root %in% seen) next
        queue <- root
        component <- integer()
        while (length(queue)) {
            node <- queue[[1L]]
            queue <- queue[-1L]
            if (node %in% component) next
            component <- c(component, node)
            next_node <- adjacency[[as.character(node)]]
            queue <- c(queue, setdiff(next_node, component))
        }
        seen <- c(seen, component)
        components[[length(components) + 1L]] <- sort(unique(component))
    }
    components
}

.condition_maximum_spanning_tree <- function(vertices, edges, reference_index) {
    vertices <- sort(unique(as.integer(vertices)))
    if (length(vertices) < 2L) return(edges[FALSE, , drop = FALSE])
    edges <- edges[
        edges$condition_a_index %in% vertices &
        edges$condition_b_index %in% vertices &
        edges$pair_estimable %in% TRUE, , drop = FALSE
    ]
    if (!nrow(edges)) return(edges)
    reference_edge <- edges$condition_a_index == reference_index |
        edges$condition_b_index == reference_index
    edges <- edges[order(
        !reference_edge, -as.numeric(edges$pair_information),
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
    edges[selected, , drop = FALSE]
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
    selected_rows <- list()
    boundary <- vector("list", p)
    vertices_used <- vector("list", p)
    for (edge in seq_len(p)) {
        candidate <- pair_information[
            pair_information$edge_index == edge &
            pair_information$pair_estimable %in% TRUE, , drop = FALSE
        ]
        components <- .condition_graph_component(seq_len(k), candidate)
        informative_components <- components[lengths(components) >= 2L]
        if (!length(informative_components)) {
            boundary[[edge]] <- seq_len(k)
            vertices_used[[edge]] <- integer()
            next
        }
        contains_reference <- vapply(
            informative_components,
            function(x) reference_index %in% x,
            logical(1)
        )
        if (any(contains_reference)) {
            component <- informative_components[[which(contains_reference)[[1L]]]]
        } else {
            sizes <- lengths(informative_components)
            component <- informative_components[[which.max(sizes)]]
        }
        tree <- .condition_maximum_spanning_tree(
            component, candidate, reference_index
        )
        if (nrow(tree)) {
            tree$contrast_coordinate <- paste0(
                tree$condition_b, "-", tree$condition_a
            )
            selected_rows[[length(selected_rows) + 1L]] <- tree
            vertices_used[[edge]] <- sort(unique(c(
                tree$condition_a_index, tree$condition_b_index
            )))
            boundary[[edge]] <- setdiff(seq_len(k), vertices_used[[edge]])
        } else {
            vertices_used[[edge]] <- integer()
            boundary[[edge]] <- seq_len(k)
        }
    }
    selected <- if (length(selected_rows)) {
        do.call(rbind, selected_rows)
    } else {
        pair_information[FALSE, , drop = FALSE]
    }
    if (nrow(selected)) {
        rownames(selected) <- NULL
        D <- matrix(0, nrow(selected), k * p)
        for (i in seq_len(nrow(selected))) {
            edge <- selected$edge_index[[i]]
            a <- selected$condition_a_index[[i]]
            b <- selected$condition_b_index[[i]]
            D[i, (a - 1L) * p + edge] <- -1
            D[i, (b - 1L) * p + edge] <- 1
        }
    } else {
        D <- matrix(0, 0L, k * p)
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
        unit <- numeric(nrow(H))
        unit[[j]] <- 1
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
        gradient[nonzero] +
        z * weights[nonzero] * sign(delta[nonzero])
    )
    value[zero] <- pmax(
        0, abs(gradient[zero]) - z * weights[zero]
    )
    max(value[active])
}

.condition_E_star_fit <- function(H, r, information, control) {
    if (!length(r)) {
        return(list(
            delta = numeric(), status = "ok", iterations = 0L,
            kkt_residual = 0, objective = 0
        ))
    }
    active <- is.finite(information) & information > 0
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
        objective <- 0.5 * drop(crossprod(
            next_value, Hs %*% next_value
        )) - drop(crossprod(rs, next_value)) +
            .condition_E_star_z * sum(weights * abs(next_value))
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
            ra <- root(a)
            rb <- root(b)
            if (ra != rb) parent[[rb]] <- ra
        }
    }
    roots <- vapply(seq_len(k), root, integer(1))
    unique_root <- unique(roots)
    match(roots, unique_root)
}

.condition_fusion_components <- function(
    tree, delta, p, conditions, control) {
    k <- length(conditions)
    component <- matrix(NA_integer_, k, p)
    boundary_condition <- matrix(FALSE, k, p)
    fused_by_penalty <- logical(p)
    shared_by_boundary <- logical(p)
    tree_metadata <- tree$metadata
    if (nrow(tree_metadata)) {
        tree_metadata$delta <- delta
        if (!"contrast_identifiable" %in% colnames(tree_metadata)) {
            tree_metadata$contrast_identifiable <- TRUE
        }
        tree_metadata$fused_by_penalty <-
            tree_metadata$contrast_identifiable %in% TRUE &
            abs(delta) <= control$fusion_tol
    } else {
        tree_metadata$delta <- numeric()
        tree_metadata$contrast_identifiable <- logical()
        tree_metadata$fused_by_penalty <- logical()
    }
    for (edge in seq_len(p)) {
        links <- list()
        index <- which(tree_metadata$edge_index == edge)
        if (length(index)) {
            for (j in index) {
                a <- tree_metadata$condition_a_index[[j]]
                b <- tree_metadata$condition_b_index[[j]]
                if (!(tree_metadata$contrast_identifiable[[j]] %in% TRUE)) {
                    links[[length(links) + 1L]] <- c(a, b)
                    shared_by_boundary[[edge]] <- TRUE
                    boundary_condition[c(a, b), edge] <- TRUE
                } else if (tree_metadata$fused_by_penalty[[j]]) {
                    links[[length(links) + 1L]] <- c(a, b)
                    fused_by_penalty[[edge]] <- TRUE
                }
            }
        }
        boundary <- tree$boundary_conditions[[edge]]
        if (length(boundary)) {
            shared_by_boundary[[edge]] <- TRUE
            boundary_condition[boundary, edge] <- TRUE
            vertices <- tree$vertices_used[[edge]]
            anchor <- if (tree$reference_index %in% vertices) {
                tree$reference_index
            } else if (length(vertices)) {
                vertices[[1L]]
            } else {
                tree$reference_index
            }
            for (node in setdiff(boundary, anchor)) {
                links[[length(links) + 1L]] <- c(anchor, node)
            }
            if (!length(vertices)) {
                for (node in setdiff(seq_len(k), anchor)) {
                    links[[length(links) + 1L]] <- c(anchor, node)
                }
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

.condition_joint_inference_design <- function(
    x, y, scaling, keep, fusion, rank_tol) {
    conditions <- names(x)
    k <- length(conditions)
    p <- length(keep)
    n <- lengths(y)
    n_total <- sum(n)
    component_key <- character()
    component_index <- matrix(NA_integer_, k, p)
    for (edge in seq_len(p)) {
        ids <- fusion$component[, edge]
        for (condition in seq_len(k)) {
            key <- paste(edge, ids[[condition]], sep = "::")
            position <- match(key, component_key)
            if (is.na(position)) {
                component_key <- c(component_key, key)
                position <- length(component_key)
            }
            component_index[condition, edge] <- position
        }
    }
    z <- matrix(0, n_total, k + length(component_key))
    response <- numeric(n_total)
    row_start <- 0L
    standardized <- vector("list", k)
    for (condition in seq_len(k)) {
        rows <- row_start + seq_len(n[[condition]])
        z[rows, condition] <- 1
        one <- sweep(
            x[[condition]][, keep, drop = FALSE],
            2L, scaling$center[keep], "-"
        )
        one <- sweep(one, 2L, scaling$scale[keep], "/")
        standardized[[condition]] <- one
        for (edge in seq_len(p)) {
            column <- k + component_index[condition, edge]
            z[rows, column] <- z[rows, column] + one[, edge]
        }
        response[rows] <- as.numeric(y[[condition]])
        row_start <- row_start + n[[condition]]
    }
    decomposition <- .condition_symmetric_pinv(crossprod(z), rank_tol)
    theta <- as.numeric(decomposition$inverse %*% crossprod(z, response))
    residual <- response - as.numeric(z %*% theta)
    residual_df <- n_total - decomposition$rank
    if (residual_df <= 0) {
        return(list(
            status = "insufficient_joint_inference_df",
            component_index = component_index,
            component_key = component_key
        ))
    }
    sigma2 <- sum(residual^2) / residual_df
    covariance <- sigma2 * decomposition$inverse
    list(
        status = "ok",
        design = z,
        theta = theta,
        covariance = covariance,
        projector = decomposition$projector,
        rank = decomposition$rank,
        residual_df = residual_df,
        sigma2 = sigma2,
        component_index = component_index,
        component_key = component_key,
        standardized = standardized
    )
}

.condition_joint_inference_fit <- function(
    joint, scaling, keep, k, p_full, conditions) {
    estimate <- se <- statistic <- pval <- matrix(
        NA_real_, k, p_full,
        dimnames = list(conditions, names(scaling$scale))
    )
    estimable <- matrix(
        FALSE, k, p_full,
        dimnames = dimnames(estimate)
    )
    component_id <- matrix(
        NA_character_, k, p_full,
        dimnames = dimnames(estimate)
    )
    hypothesis_id <- matrix(
        NA_character_, k, p_full,
        dimnames = dimnames(estimate)
    )
    if (!identical(joint$status, "ok")) {
        return(list(
            estimate = estimate, se = se, statistic = statistic,
            pval = pval, estimable = estimable,
            component_id = component_id, hypothesis_id = hypothesis_id,
            sigma2 = NA_real_, residual_df = NA_integer_, rank = NA_integer_
        ))
    }
    n_condition <- k
    p <- length(keep)
    for (condition in seq_len(k)) {
        for (local_edge in seq_len(p)) {
            full_edge <- keep[[local_edge]]
            component <- joint$component_index[condition, local_edge]
            column <- n_condition + component
            unit <- numeric(length(joint$theta))
            unit[[column]] <- 1
            residual <- unit - joint$projector[column, ]
            ok <- max(abs(residual)) <= 1e-7
            component_id[condition, full_edge] <- paste0(
                "edge", full_edge, "_component", component
            )
            hypothesis_id[condition, full_edge] <- paste0(
                "H0_", component_id[condition, full_edge], "_eq_0"
            )
            if (!ok) next
            value_z <- joint$theta[[column]]
            var_z <- joint$covariance[column, column]
            if (!is.finite(var_z) || var_z < 0) next
            raw_scale <- scaling$scale[[full_edge]]
            value <- value_z / raw_scale
            std <- sqrt(max(0, var_z)) / raw_scale
            estimate[condition, full_edge] <- value
            se[condition, full_edge] <- std
            estimable[condition, full_edge] <- is.finite(std) && std > 0
            if (estimable[condition, full_edge]) {
                statistic[condition, full_edge] <- value / std
                pval[condition, full_edge] <- 2 * stats::pnorm(
                    -abs(statistic[condition, full_edge])
                )
            }
        }
    }
    list(
        estimate = estimate, se = se, statistic = statistic,
        pval = pval, estimable = estimable,
        component_id = component_id, hypothesis_id = hypothesis_id,
        sigma2 = joint$sigma2, residual_df = joint$residual_df,
        rank = joint$rank
    )
}

.condition_joint_contrast <- function(
    joint, scaling, keep, condition_a, condition_b, local_edge) {
    if (!identical(joint$status, "ok")) {
        return(list(
            estimate = NA_real_, se = NA_real_, statistic = NA_real_,
            pval = NA_real_, estimable = FALSE, same_component = FALSE
        ))
    }
    component_a <- joint$component_index[condition_a, local_edge]
    component_b <- joint$component_index[condition_b, local_edge]
    if (component_a == component_b) {
        return(list(
            estimate = 0, se = NA_real_, statistic = NA_real_,
            pval = NA_real_, estimable = TRUE, same_component = TRUE
        ))
    }
    k <- nrow(joint$component_index)
    column_a <- k + component_a
    column_b <- k + component_b
    l <- numeric(length(joint$theta))
    l[[column_a]] <- -1
    l[[column_b]] <- 1
    residual <- as.numeric(l - l %*% joint$projector)
    ok <- max(abs(residual)) <= 1e-7
    if (!ok) {
        return(list(
            estimate = NA_real_, se = NA_real_, statistic = NA_real_,
            pval = NA_real_, estimable = FALSE, same_component = FALSE
        ))
    }
    raw_scale <- scaling$scale[[keep[[local_edge]]]]
    estimate <- as.numeric(l %*% joint$theta) / raw_scale
    variance <- as.numeric(l %*% joint$covariance %*% l) / raw_scale^2
    if (!is.finite(variance) || variance <= 0) {
        return(list(
            estimate = estimate, se = NA_real_, statistic = NA_real_,
            pval = NA_real_, estimable = FALSE, same_component = FALSE
        ))
    }
    se <- sqrt(variance)
    statistic <- estimate / se
    list(
        estimate = estimate, se = se, statistic = statistic,
        pval = 2 * stats::pnorm(-abs(statistic)),
        estimable = TRUE, same_component = FALSE
    )
}

.condition_scheme_e_fit <- function(
    x, y, scaling, min_residual_df = 1L, inference = TRUE,
    control = list(), reference_condition = NULL) {
    control <- .condition_E_star_control(control)
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
    if (is.null(reference_condition)) reference_condition <- conditions[[1L]]
    reference_condition <- as.character(reference_condition)
    if (length(reference_condition) != 1L ||
        !reference_condition %in% conditions) {
        return(list(status = "invalid_reference_condition"))
    }

    informative <- as.logical(scaling$informative)
    beta <- beta_z <- deviation_z <- matrix(
        0, k, p_full,
        dimnames = list(conditions, colnames(x[[1L]]))
    )
    shared_z <- numeric(p_full)
    raw_information <- matrix(0, k, p_full, dimnames = dimnames(beta))
    zero_variance <- matrix(
        FALSE, k, p_full, dimnames = dimnames(beta)
    )
    for (i in seq_len(k)) {
        value <- apply(x[[i]], 2L, stats::var)
        zero_variance[i, ] <- !is.finite(value) |
            value <= scaling$floor^2
    }

    if (!any(informative)) {
        intercept <- vapply(y, mean, numeric(1))
        prediction <- Map(
            function(value, a) rep(a, length(value)), y, intercept
        )
        rsq <- vapply(y, function(value) {
            if (stats::var(value) > 0) 0 else NA_real_
        }, numeric(1))
        empty_inference <- matrix(
            NA_real_, k, p_full, dimnames = dimnames(beta)
        )
        empty_bool <- matrix(
            FALSE, k, p_full, dimnames = dimnames(beta)
        )
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
            contrast_tree = data.frame(),
            solver_status = "ok", kkt_residual = 0,
            iterations = 0L, objective = 0,
            penalty_family = .condition_E_star_penalty_family,
            penalty_value = .condition_E_star_z,
            reference_condition = reference_condition,
            condition_weight = stats::setNames(rep(1, k), conditions),
            inference_schema = .condition_E_star_inference_schema,
            inference_estimate = empty_inference,
            inference_se = empty_inference,
            inference_statistic = empty_inference,
            inference_pval = empty_inference,
            inference_estimable = empty_bool,
            inference_component_id = matrix(
                NA_character_, k, p_full, dimnames = dimnames(beta)
            ),
            inference_hypothesis_id = matrix(
                NA_character_, k, p_full, dimnames = dimnames(beta)
            ),
            inference_sigma2 = NA_real_,
            inference_residual_df = NA_integer_,
            inference_rank = NA_integer_,
            orthogonality_error = 0, dr_error = 0,
            joint_inference = NULL
        ))
    }

    keep <- which(informative)
    p <- length(keep)
    gram_raw <- rhs_raw <- vector("list", k)
    xc <- yc <- vector("list", k)
    raw_rank <- integer(k)
    raw_kappa <- numeric(k)
    sse_raw <- numeric(k)
    for (i in seq_len(k)) {
        zmat <- sweep(
            x[[i]][, keep, drop = FALSE],
            2L, scaling$center[keep], "-"
        )
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
    variance_floor <- .Machine$double.eps * max(
        1, pooled_y_scale, na.rm = TRUE
    )
    sigma2 <- max(sum(sse_raw) / residual_df, variance_floor)
    q <- lapply(gram_raw, function(one) one / sigma2)
    h_condition <- lapply(rhs_raw, function(one) one / sigma2)
    for (i in seq_len(k)) {
        raw_information[i, keep] <- diag(q[[i]]) *
            scaling$scale[keep]^2
    }

    Q <- .condition_blockdiag(q)
    h <- unlist(h_condition, use.names = FALSE)
    A <- kronecker(matrix(1, nrow = k, ncol = 1), diag(p))
    tree <- .condition_identifiable_contrast_tree(
        Q = Q, p = p, conditions = conditions,
        reference_condition = reference_condition,
        rank_tol = control$rank_tol
    )
    geometry <- .condition_q_orthogonal_decomposition(
        Q, A, tree$D, control$rank_tol
    )
    shared_score <- as.numeric(crossprod(A, h))
    mu <- as.numeric(geometry$shared_inverse %*% shared_score)
    if (ncol(geometry$R)) {
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
                status = "solver_failed",
                solver_status = solved$status,
                kkt_residual = solved$kkt_residual,
                iterations = solved$iterations
            ))
        }
        delta <- solved$delta
        beta_vector <- as.numeric(A %*% mu + geometry$R %*% delta)
        tree$metadata$profile_information <-
            coordinate_information$information
        tree$metadata$contrast_identifiable <-
            coordinate_information$estimable
    } else {
        H <- matrix(0, 0L, 0L)
        delta <- numeric()
        solved <- list(
            status = "ok", iterations = 0L,
            kkt_residual = 0, objective = 0
        )
        beta_vector <- as.numeric(A %*% mu)
    }

    beta_z_small <- matrix(beta_vector, nrow = k, byrow = TRUE)
    beta_small <- sweep(
        beta_z_small, 2L, scaling$scale[keep], "/"
    )
    beta_z[, keep] <- beta_z_small
    beta[, keep] <- beta_small
    shared_z[keep] <- mu
    deviation_z[, keep] <- sweep(beta_z_small, 2L, mu, "-")

    fusion <- .condition_fusion_components(
        tree, delta, p, conditions, control
    )
    fusion_component_id <- matrix(
        NA_character_, k, p_full, dimnames = dimnames(beta)
    )
    boundary_condition <- matrix(
        TRUE, k, p_full, dimnames = dimnames(beta)
    )
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
    if (nrow(fusion$contrast_tree)) {
        for (local_edge in seq_len(p)) {
            index <- fusion$contrast_tree$edge_index == local_edge &
                fusion$contrast_tree$contrast_identifiable %in% TRUE
            value <- fusion$contrast_tree$profile_information[index]
            value <- value[is.finite(value) & value > 0]
            if (length(value)) {
                full_edge <- keep[[local_edge]]
                profile_information[[full_edge]] <- min(value)
                contrast_identifiable[[full_edge]] <- TRUE
            }
        }
    }
    shared_by_boundary <- rep(TRUE, p_full)
    fused_by_penalty <- rep(FALSE, p_full)
    shared_edge <- rep(TRUE, p_full)
    shared_by_boundary[keep] <- fusion$shared_by_boundary
    fused_by_penalty[keep] <- fusion$fused_by_penalty
    shared_edge[keep] <- fusion$shared_edge

    joint <- if (isTRUE(inference)) {
        .condition_joint_inference_design(
            x, y, scaling, keep, fusion, control$rank_tol
        )
    } else {
        list(status = "not_requested")
    }
    joint_rows <- .condition_joint_inference_fit(
        joint, scaling, keep, k, p_full, conditions
    )

    intercept <- numeric(k)
    prediction <- vector("list", k)
    rsq <- numeric(k)
    for (i in seq_len(k)) {
        intercept[[i]] <- mean(y[[i]]) -
            sum(colMeans(x[[i]][, keep, drop = FALSE]) * beta_small[i, ])
        prediction[[i]] <- as.numeric(
            intercept[[i]] +
            x[[i]][, keep, drop = FALSE] %*% beta_small[i, ]
        )
        residual <- as.numeric(y[[i]]) - prediction[[i]]
        tss <- sum((as.numeric(y[[i]]) - mean(y[[i]]))^2)
        rsq[[i]] <- if (tss > 0) {
            1 - sum(residual^2) / tss
        } else NA_real_
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
        solver_status = solved$status,
        kkt_residual = solved$kkt_residual,
        iterations = solved$iterations,
        objective = solved$objective,
        penalty_family = .condition_E_star_penalty_family,
        penalty_value = .condition_E_star_z,
        reference_condition = reference_condition,
        condition_weight = stats::setNames(rep(1, k), conditions),
        inference_schema = .condition_E_star_inference_schema,
        inference_estimate = joint_rows$estimate,
        inference_se = joint_rows$se,
        inference_statistic = joint_rows$statistic,
        inference_pval = joint_rows$pval,
        inference_estimable = joint_rows$estimable,
        inference_component_id = joint_rows$component_id,
        inference_hypothesis_id = joint_rows$hypothesis_id,
        inference_sigma2 = joint_rows$sigma2,
        inference_residual_df = joint_rows$residual_df,
        inference_rank = joint_rows$rank,
        joint_inference = joint,
        orthogonality_error = geometry$orthogonality_error,
        dr_error = geometry$dr_error
    )
}
