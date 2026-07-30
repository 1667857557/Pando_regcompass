test_that("projection API defaults to structural-zero semantics", {
    formals_projection <- formals(Pando::project_condition_grn_cells)
    expect_identical(
        eval(formals_projection$nonestimable),
        c("structural_zero", "error")
    )
})

test_that("structural zeros enter metacell aggregation", {
    projection <- structure(list(
        gene_score = matrix(
            c(1, 0, NA, 2), nrow = 2, byrow = TRUE,
            dimnames = list(c("c1", "c2"), c("g1", "g2"))
        ),
        gene_structural_zero_mask = matrix(
            c(FALSE, TRUE, TRUE, FALSE), nrow = 2, byrow = TRUE,
            dimnames = list(c("c1", "c2"), c("g1", "g2"))
        ),
        cell_metadata = data.frame(
            cell_id = c("c1", "c2"),
            cell_type = c("T", "T"),
            condition = c("A", "A"),
            row.names = c("c1", "c2"),
            stringsAsFactors = FALSE
        ),
        projection_origin = "outer_condition_stratified_cell_oof",
        projection_used_for_penalty = TRUE,
        full_fit_projection_used_for_penalty = FALSE,
        score_comparability_class = "primary_common_support_comparable"
    ), class = c("ConditionGRNProjection", "list"))
    membership <- data.frame(
        cell_id = c("c1", "c2"),
        metacell_id = c("m1", "m1"),
        stringsAsFactors = FALSE
    )
    out <- aggregate_condition_grn_projection(projection, membership)
    expect_true(all(is.finite(out$gene_score)))
    expect_equal(unname(out$gene_score["m1", ]), c(0.5, 1))
    expect_equal(
        unname(out$gene_structural_zero_fraction["m1", ]),
        c(0.5, 0.5)
    )
    expect_identical(out$nonestimable_policy, "structural_zero")
    expect_true(out$structural_zero_enters_downstream)
})
