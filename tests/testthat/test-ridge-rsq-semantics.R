test_that("conditional R2 is full-data Scheme E diagnostic", {
  set.seed(11)
  n <- 40L
  cells <- list(A = paste0("a", seq_len(n)), B = paste0("b", seq_len(n)))
  all_cells <- unlist(cells, use.names = FALSE)
  tf <- rnorm(2L * n)
  peak <- runif(2L * n, 0.2, 2)
  interaction <- tf * peak
  target <- c(
    1 + 1.8 * interaction[seq_len(n)] + rnorm(n, sd = 0.2),
    1 - 0.9 * interaction[n + seq_len(n)] + rnorm(n, sd = 0.2)
  )
  prepared <- list(
    gene_data = cbind(TF1 = tf, G = target), peak_data = cbind(P1 = peak)
  )
  rownames(prepared$gene_data) <- all_cells
  rownames(prepared$peak_data) <- all_cells
  edges <- data.frame(
    tf = "TF1", target = "G", region = "P1", edge_id = "G||TF1||P1",
    atac_feature_id = "P1", candidate_index = 1L, source_global = TRUE,
    source_conditions = "A;B", n_sources = 3L, stringsAsFactors = FALSE
  )
  control <- Pando:::.condition_ridge_control()
  result <- Pando:::.condition_ridge_target(
    prepared = prepared, edges = edges, cells_by_condition = cells,
    folds = NULL, control = control, min_residual_df = 1L,
    rank_action = "mark"
  )
  expect_true(all(is.finite(result$fit$rsq)))
  expect_true(all(result$fit$rsq <= 1 + 1e-12))
  expect_true(all(result$fit$penalty_family == "exact_edge_sparse_deviation"))
  expect_true(all(result$fit$deviation_z == 0.25))
  expect_false(any(c("lambda", "lambda_min", "cv_mse", "rsq_oof") %in%
                   names(result$fit)))
})

test_that("standard ridge alone retains lambda CV", {
  set.seed(23)
  n <- 35L
  x <- cbind(e1 = rnorm(n), e2 = rnorm(n))
  y <- 0.8 * x[, "e1"] - 0.2 * x[, "e2"] + rnorm(n, sd = 0.3)
  control <- Pando:::.pando_standard_ridge_control(list(
    lambda_grid = c(0.01, 0.1, 1), lambda_rule = "1se",
    cv_folds = 5L, seed = 7L
  ))
  folds <- Pando:::.pando_standard_ridge_folds(n, 5L, 7L)
  cv <- Pando:::.pando_standard_ridge_cv(
    x, y, folds, control, min_residual_df = 1L
  )
  expect_named(cv, c("lambda", "lambda_min", "cv_mse", "cv_se", "curve"))
  expect_true(cv$lambda %in% control$lambda_grid)
  expect_true(is.finite(cv$cv_mse))
})