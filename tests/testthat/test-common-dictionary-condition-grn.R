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
    expect_false("G||TF1||P2" %in% dictionary$edge_id)
    expect_false("G||TF2||P1" %in% dictionary$edge_id)
})

test_that("candidate union records provenance but not significance", {
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
    expect_false(dictionary$source_global[dictionary$edge_id == "G||TF3||P3"])
    expect_false(any(c("pval", "padj", "significant") %in%
                     colnames(dictionary)))
})

test_that("package metadata advertises separated inference", {
    description <- utils::packageDescription("Pando")
    expect_identical(
        description[["Config/Pando/ConditionGRNSchema"]],
        "pando_condition_grn_common_dictionary_v1"
    )
    expect_identical(
        description[["Config/Pando/ConditionGRNModelSchema"]],
        "pando_condition_grn_Estar_z025_inference_separated_v1"
    )
    expect_identical(
        description[["Config/Pando/ConditionGRNMethod"]],
        "global-condition-union-Estar-z025-separated-inference"
    )
    expect_identical(
        description[["Config/Pando/ConditionProjectionPolicy"]],
        "exact-edge-whole-network-BH-common-topology"
    )
    expect_identical(
        description[["Config/Pando/ConditionInference"]],
        "no-fusion-condition-local-lm-edge-omnibus-whole-network-BH"
    )
})
