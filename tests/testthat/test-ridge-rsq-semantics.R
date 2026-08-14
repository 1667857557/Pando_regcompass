test_that("canonical ridge rsq is selected-lambda full-data R2", {
  set.seed(11)
  n <- 30L
  cells <- list(
    A = paste0("a", seq_len(n)),
    B = paste0("b", seq_len(n))
  )
  all_cells <- unlist(cells, use.names = FALSE)
  tf <- rnorm(2L * n)
  peak <- runif(2L * n, 0.2, 2)
  interaction <- tf * peak
  target <- c(
    1 + 1.8 * interaction[seq_len(n)] + rnorm(n, sd = 0.2),
    1 - 0.9 * interaction[n + seq_len(n)] + rnorm(n, sd = 0.2)
  )
  prepared <- list(
    gene_data = cbind(TF1 = tf, G = target),
    peak_data = cbind(P1 = peak)
  )
  rownames(prepared$gene_data) <- all_cells
  rownames(prepared$peak_data) <- all_cells
  edges <- data.frame(
    tf = "TF1", target = "G", region = "P1",
    edge_id = "G||TF1||P1", atac_feature_id = "P1",
    candidate_index = 1L, source_global = TRUE,
    source_conditions = "A;B", n_sources = 3L,
    stringsAsFactors = FALSE
  )
  control <- Pando:::.condition_ridge_control(list(
    lambda_grid = c(0.01, 0.1, 1),
    lambda_rule = "min", cv_folds = 5L, seed = 5L
  ))
  folds <- Pando:::.condition_ridge_folds(cells, 5L, 5L)
  result <- Pando:::.condition_ridge_target(
    prepared = prepared,
    edges = edges,
    cells_by_condition = cells,
    folds = folds,
    control = control,
    min_residual_df = 1L,
    rank_action = "mark"
  )
  expect_equal(result$fit$rsq, result$fit$rsq_in_sample, tolerance = 1e-12)
  expect_equal(result$fit$rsq_oof, unname(result$cv$rsq_oof), tolerance = 1e-12)
  expect_true(all(is.finite(result$fit$rsq)))
  expect_true(all(result$fit$rsq <= 1 + 1e-12))
})

test_that("OOF R2 is recovered from the lambda-selection pass", {
  set.seed(23)
  n <- 25L
  x <- list(
    A = cbind(e1 = rnorm(n), e2 = rnorm(n)),
    B = cbind(e1 = rnorm(n), e2 = rnorm(n))
  )
  y <- list(
    A = 0.8 * x$A[, "e1"] - 0.2 * x$A[, "e2"] + rnorm(n, sd = 0.3),
    B = -0.5 * x$B[, "e1"] + 0.6 * x$B[, "e2"] + rnorm(n, sd = 0.3)
  )
  cells <- list(A = seq_len(n), B = seq_len(n))
  control <- Pando:::.condition_ridge_control(list(
    lambda_grid = c(0.01, 0.1, 1),
    lambda_rule = "1se", cv_folds = 5L, seed = 7L
  ))
  folds <- Pando:::.condition_ridge_folds(cells, 5L, 7L)
  cv <- Pando:::.condition_ridge_cv(
    x, y, folds, control, min_residual_df = 1L
  )

  # Reconstruct the historical second pass only inside the test. It must be
  # numerically identical to the sufficient-statistic result returned above.
  legacy_oof <- lapply(y, function(value) rep(NA_real_, length(value)))
  for (fold in seq_len(attr(folds, "nfolds"))) {
    train_x <- train_y <- vector("list", length(x))
    names(train_x) <- names(train_y) <- names(x)
    for (i in seq_along(x)) {
      train <- folds[[i]] != fold
      train_x[[i]] <- x[[i]][train, , drop = FALSE]
      train_y[[i]] <- y[[i]][train]
    }
    scaling <- Pando:::.condition_ridge_scaling(
      train_x, control$scale_floor
    )
    fit <- Pando:::.condition_ridge_fit(
      train_x, train_y, scaling, cv$lambda,
      min_residual_df = 1L, inference = FALSE
    )
    expect_identical(fit$status, "ok")
    for (i in seq_along(x)) {
      index <- which(folds[[i]] == fold)
      legacy_oof[[i]][index] <- fit$intercept[[i]] +
        x[[i]][index, , drop = FALSE] %*% fit$beta[i, ]
    }
  }
  legacy_rsq <- vapply(seq_along(y), function(i) {
    value <- as.numeric(y[[i]])
    pred <- as.numeric(legacy_oof[[i]])
    1 - sum((value - pred)^2) / sum((value - mean(value))^2)
  }, numeric(1))
  names(legacy_rsq) <- names(y)
  expect_equal(cv$rsq_oof, legacy_rsq, tolerance = 1e-12)

  body_text <- paste(deparse(body(Pando:::.condition_ridge_cv)), collapse = "\n")
  expect_match(body_text, "oof_sse", fixed = TRUE)
  expect_false(grepl("oof <- lapply", body_text, fixed = TRUE))
  expect_equal(
    lengths(regmatches(
      body_text,
      gregexpr("for (fold in seq_len(nfolds))", body_text, fixed = TRUE)
    )),
    1L
  )
})
