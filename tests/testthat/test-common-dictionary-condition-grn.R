test_that("exact edge union does not create Cartesian edges", {
    global <- data.frame(
        target = "G", tf = "TF1", region = "P1",
        atac_feature_id = "A1", peak_target_cor = 0.2,
        tf_target_cor = 0.3, stringsAsFactors = FALSE
    )
    condition <- list(
        A = data.frame(
            target = "G", tf = "TF2", region = "P2",
            atac_feature_id = "A2", peak_target_cor = 0.4,
            tf_target_cor = 0.5, stringsAsFactors = FALSE
        ),
        B = data.frame(
            target = "G", tf = "TF1", region = "P1",
            atac_feature_id = "A1", peak_target_cor = 0.1,
            tf_target_cor = 0.2, stringsAsFactors = FALSE
        )
    )
    dictionary <- union_grn_edges(global, condition)
    expect_identical(nrow(dictionary), 2L)
    expect_setequal(dictionary$edge_id, c("G||TF1||P1", "G||TF2||P2"))
    expect_identical(anyDuplicated(dictionary$edge_id), 0L)
    expect_false("G||TF1||P2" %in% dictionary$edge_id)
    expect_false("G||TF2||P1" %in% dictionary$edge_id)
})

test_that("exact edge union records pooled/global and condition provenance", {
    global <- data.frame(
        target = c("G", "G"), tf = c("TF1", "TF2"),
        region = c("P1", "P2"), atac_feature_id = c("A1", "A2"),
        peak_target_cor = c(0.2, 0.3), tf_target_cor = c(0.3, 0.4),
        stringsAsFactors = FALSE
    )
    condition <- list(
        A = data.frame(
            target = c("G", "G"), tf = c("TF1", "TF3"),
            region = c("P1", "P3"), atac_feature_id = c("A1", "A3"),
            peak_target_cor = c(0.25, 0.35), tf_target_cor = c(0.35, 0.45),
            stringsAsFactors = FALSE
        )
    )
    dictionary <- union_grn_edges(global, condition)
    expect_setequal(dictionary$edge_id,
                    c("G||TF1||P1", "G||TF2||P2", "G||TF3||P3"))
    expect_true(dictionary$source_global[dictionary$edge_id == "G||TF1||P1"])
    expect_true(dictionary$source_global[dictionary$edge_id == "G||TF2||P2"])
    expect_false(dictionary$source_global[dictionary$edge_id == "G||TF3||P3"])
    expect_identical(
        dictionary$source_conditions[dictionary$edge_id == "G||TF1||P1"], "A"
    )
})

test_that("external dictionaries preserve domain and motif support", {
    prepared <- list(
        gene_data = matrix(
            seq_len(9), nrow = 3,
            dimnames = list(paste0("c", 1:3), c("G", "TF1", "TF2"))
        ),
        peak_data = matrix(
            seq_len(6), nrow = 3,
            dimnames = list(paste0("c", 1:3), c("P1", "P2"))
        ),
        region_map = data.frame(
            region = c("P1", "P2"), atac_feature_id = c("A1", "A2"),
            stringsAsFactors = FALSE
        ),
        peaks2gene = Matrix::Matrix(
            matrix(c(1, 0), nrow = 1,
                   dimnames = list("G", c("P1", "P2"))), sparse = TRUE
        ),
        peaks2motif = Matrix::Matrix(
            matrix(c(1, 0), nrow = 2,
                   dimnames = list(c("P1", "P2"), "M1")), sparse = TRUE
        ),
        motif2tf = Matrix::Matrix(
            matrix(c(1, 0), nrow = 1,
                   dimnames = list("M1", c("TF1", "TF2"))), sparse = TRUE
        )
    )
    valid <- data.frame(
        edge_id = "G||TF1||P1", target = "G", tf = "TF1",
        region = "P1", atac_feature_id = "A1", candidate_index = 1L,
        stringsAsFactors = FALSE
    )
    expect_invisible(Pando:::.condition_validate_dictionary(valid, prepared))
    outside_domain <- valid
    outside_domain$region <- "P2"
    outside_domain$atac_feature_id <- "A2"
    outside_domain$edge_id <- "G||TF1||P2"
    expect_error(Pando:::.condition_validate_dictionary(outside_domain, prepared),
                 "outside the Pando domain")
    unsupported_tf <- valid
    unsupported_tf$tf <- "TF2"
    unsupported_tf$edge_id <- "G||TF2||P1"
    expect_error(Pando:::.condition_validate_dictionary(unsupported_tf, prepared),
                 "without motif support")
})

test_that("condition API exposes one fixed E-star/JSE production design", {
    formal_names <- names(formals(Pando:::infer_condition_grn.GRNData))
    expect_true(all(c(
        "tf_cor", "peak_cor", "adjust_method", "padj_threshold",
        "rank_action", "min_residual_df"
    ) %in% formal_names))
    expect_equal(eval(formals(Pando:::infer_condition_grn.GRNData)$tf_cor), 0.05)
    expect_equal(eval(formals(Pando:::infer_condition_grn.GRNData)$peak_cor), 0.05)
    expect_false(any(c(
        "candidate_screen", "condition_mix", "condition_weight",
        "nlambda", "lambda", "outer_nfolds", "inner_nfolds",
        "lambda_selection", "engine_control", "scale", "fusion_ratio",
        "scheme_e_z", "z", "condition_ridge_control"
    ) %in% formal_names))
    description <- utils::packageDescription("Pando")
    expect_identical(
        description[["Config/Pando/ConditionGRNSchema"]],
        "pando_condition_grn_common_dictionary_v1"
    )
    expect_identical(
        description[["Config/Pando/ConditionGRNModelSchema"]],
        "pando_condition_grn_Estar_jointse_v1"
    )
    expect_identical(
        description[["Config/Pando/ConditionGRNMethod"]],
        "global-condition-union-Estar-z025-JSE"
    )
    expect_identical(
        description[["Config/Pando/ConditionProjectionPolicy"]],
        "any-condition-padj-exact-edge-union"
    )
})
