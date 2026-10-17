# 通路富集置换校准补充表（返修：审稿人2第三段） ####
# 读取6.2.16保存的结果表，整理为投稿用补充表（英文表头），不做新的统计分析

rm(list = ls())
gc()

# 1. 加载配置 ####
source("Code/config.R")

# 2. 结果目录与输入文件 ####
results_dir_name <- "6.2.18通路富集置换校准补充表"
script_results_path <- file.path(results_path, results_dir_name)

if (!dir.exists(script_results_path)) {
    dir.create(script_results_path, recursive = TRUE, showWarnings = FALSE)
}

input_xlsx_file <- file.path(results_path, "6.2.16通路富集置换校准分析", "通路富集置换校准结果.xlsx")
output_xlsx_file <- file.path(script_results_path, "通路富集置换校准补充表.xlsx")

if (!file.exists(input_xlsx_file)) {
    stop("输入结果文件不存在: ", input_xlsx_file)
}

pattern_levels <- c("Ocular", "Hematologic", "Neurologic", "Systemic-like", "Visceral-like")
database_labels <- c(
    "Hallmark" = "MSigDB Hallmark",
    "GO_BP" = "GO Biological Process",
    "KEGG" = "KEGG",
    "Reactome" = "Reactome",
    "C7_IMMUNESIGDB" = "MSigDB C7 (ImmuneSigDB)"
)

message("========== 通路富集置换校准补充表整理开始 ==========")


# 3. 数据读取 ####
calibration_df <- readxl::read_xlsx(input_xlsx_file, sheet = "通路校准结果")


# 4. 补充表：Figure 7通路的原ORA结果与置换校准结果 ####
calibration_table <- calibration_df %>%
    dplyr::mutate(
        pattern = factor(pattern, levels = pattern_levels),
        hit_direction = factor(hit_direction, levels = c("positive", "negative"))
    ) %>%
    dplyr::arrange(pattern, hit_direction, enrichment_pvalue) %>%
    dplyr::transmute(
        Pattern = as.character(pattern),
        Direction = ifelse(hit_direction == "positive", "Positive Δρ", "Negative Δρ"),
        Database = unname(database_labels[database]),
        Pathway = term_id,
        `Candidate genes in pathway` = n_hit_in_term,
        `Pathway genes in background` = n_term_in_universe,
        `Fold enrichment` = round(fold_enrichment, 2),
        `ORA P` = signif(enrichment_pvalue, 3),
        `ORA BH-adjusted P` = signif(enrichment_padj_BH, 3),
        `Permutations with ORA P ≤ observed` = n_permutations_as_extreme,
        `Calibrated P` = signif(calibrated_p, 3),
        `Calibrated BH-adjusted P` = signif(calibrated_p_BH, 3)
    )


# 5. 结果保存 ####
writexl::write_xlsx(
    list("通路富集置换校准" = calibration_table),
    path = output_xlsx_file
)

message("已生成补充表: ", output_xlsx_file)
message("========== 通路富集置换校准补充表整理结束 ==========")
