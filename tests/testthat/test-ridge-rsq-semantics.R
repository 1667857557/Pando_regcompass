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
