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

  x <- Pando:::.condition_ridge_predictors(prepared, edges, cells)
  y <- lapply(cells, function(one) as.numeric(prepared$gene_data[one, "G"]))
  scaling <- Pando:::.condition_ridge_scaling(x, control$scale_floor)
  full_fit <- Pando:::.condition_ridge_fit(
    x, y, scaling, result$cv$lambda,
    min_residual_df = 1L, inference = TRUE
  )

  expect_equal(result$fit$rsq, unname(full_fit$rsq), tolerance = 1e-12)
  expect_true(all(is.finite(result$fit$rsq)))
  expect_true(all(result$fit$rsq <= 1 + 1e-12))
  expect_false(any(c("rsq_oof", "rsq_in_sample") %in% names(result$fit)))
})

test_that("ridge CV only selects lambda and exports CV loss diagnostics", {
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

  expect_named(cv, c("lambda", "lambda_min", "cv_mse", "cv_se", "curve"))
  expect_true(cv$lambda %in% control$lambda_grid)
  expect_true(cv$lambda_min %in% control$lambda_grid)
  expect_true(is.finite(cv$cv_mse))
  expect_equal(nrow(cv$curve), length(control$lambda_grid))

  body_text <- paste(deparse(body(Pando:::.condition_ridge_cv)), collapse = "\n")
  expect_false(grepl("oof", body_text, ignore.case = TRUE))
  expect_equal(
    lengths(regmatches(
      body_text,
      gregexpr("for (fold in seq_len(nfolds))", body_text, fixed = TRUE)
    )),
    1L
  )
})
