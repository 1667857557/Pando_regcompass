test_that("native wrappers resolve registered symbols explicitly", {
    wrappers <- c(
        ".condition_product_matrix_cpp",
        ".condition_fit_multitask_path_cpp",
        ".condition_refit_path_cpp",
        ".condition_fit_target_engine_cpp"
    )
    namespace <- asNamespace("Pando")
    helper <- get(".pando_registered_symbol", namespace, inherits = FALSE)
    symbols <- c(
        "_Pando_condition_product_matrix_cpp",
        "_Pando_condition_fit_multitask_path_cpp",
        "_Pando_condition_refit_path_cpp",
        "_Pando_condition_fit_target_engine_cpp"
    )

    for (wrapper in wrappers) {
        value <- get(wrapper, namespace, inherits = FALSE)
        expect_match(
            paste(deparse(body(value)), collapse = "\n"),
            ".pando_registered_call",
            fixed = TRUE,
            info = wrapper
        )
    }
    for (symbol in symbols) {
        info <- helper(symbol)
        expect_true(is.list(info) && !is.null(info$address), info = symbol)
    }
})

test_that("registered symbols resolve in a clean PSOCK worker", {
    skip_on_cran()
    cluster <- parallel::makePSOCKcluster(1L)
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    symbols <- c(
        "_Pando_condition_product_matrix_cpp",
        "_Pando_condition_fit_multitask_path_cpp",
        "_Pando_condition_refit_path_cpp",
        "_Pando_condition_fit_target_engine_cpp"
    )
    result <- parallel::clusterCall(
        cluster,
        function(symbols, library_paths) {
            .libPaths(library_paths)
            loadNamespace("Pando")
            helper <- getFromNamespace(".pando_registered_symbol", "Pando")
            vapply(symbols, function(symbol) {
                info <- helper(symbol)
                is.list(info) && !is.null(info$address)
            }, logical(1))
        },
        symbols,
        .libPaths()
    )
    expect_true(all(unlist(result, use.names = FALSE)))
})
