# ICI类别分层补充表（返修：审稿人1第12条） ####
# 读取6.2.7与6.2.9保存的结果表，整理为投稿用补充表（英文表头），不做新的统计分析

rm(list = ls())
gc()

# 1. 加载配置 ####
source("Code/config.R")

# 2. 结果目录与输入文件 ####
results_dir_name <- "6.2.11ICI类别分层补充表"
script_results_path <- file.path(results_path, results_dir_name)

if (!dir.exists(script_results_path)) {
    dir.create(script_results_path, recursive = TRUE, showWarnings = FALSE)
}

stratification_xlsx_file <- file.path(results_path, "6.2.7ICI类别分层五模式稳定性分析", "ICI类别分层五模式稳定性分析结果.xlsx")
sparsity_xlsx_file <- file.path(results_path, "6.2.9ICI类别分层数据稀疏性分析", "ICI类别分层数据稀疏性分析结果.xlsx")
output_xlsx_file <- file.path(script_results_path, "ICI类别分层补充表.xlsx")

for (current_file in c(stratification_xlsx_file, sparsity_xlsx_file)) {
    if (!file.exists(current_file)) {
        stop("输入结果文件不存在: ", current_file)
    }
}

dataset_labels <- c(
    "All ICI (main analysis)" = "All ICI (main analysis)",
    "anti-PD-1" = "Anti-PD-1 monotherapy",
    "anti-PD-L1" = "Anti-PD-L1 monotherapy",
    "anti-CTLA-4+anti-PD-1" = "Anti-CTLA-4 + anti-PD-1"
)
pattern_levels <- c("Ocular", "Hematologic", "Neurologic", "Systemic-like", "Visceral-like")

message("========== ICI类别分层补充表整理开始 ==========")


# 3. 数据读取 ####
dilution_df <- readxl::read_xlsx(sparsity_xlsx_file, sheet = "稀释汇总_全部癌种")
precision_df <- readxl::read_xlsx(sparsity_xlsx_file, sheet = "精度汇总_纳入癌种")
cells_df <- readxl::read_xlsx(sparsity_xlsx_file, sheet = "单元明细")
stratum_summary_df <- readxl::read_xlsx(stratification_xlsx_file, sheet = "分层投影汇总")
denovo_summary_df <- readxl::read_xlsx(stratification_xlsx_file, sheet = "抗PD1重新NMF汇总")
denovo_match_df <- readxl::read_xlsx(stratification_xlsx_file, sheet = "抗PD1新模式与最终模式匹配")


# 4. A表：各分层的数据稀疏性与ROR精度 ####
zero_event_df <- cells_df %>%
    dplyr::filter(n_events == 0) %>%
    dplyr::group_by(dataset) %>%
    dplyr::summarise(zero_event_ror_above_1 = paste0(sum(ror > 1), "/", dplyr::n()), .groups = "drop")

sparsity_table <- dilution_df %>%
    dplyr::select(dataset, n_cancer_types_ge_threshold, median_events_per_cell, prop_cells_0_events, prop_cells_lt5_events, median_ci_ratio) %>%
    dplyr::left_join(zero_event_df, by = "dataset") %>%
    dplyr::left_join(
        precision_df %>% dplyr::select(dataset, median_abs_log2_diff_vs_main, spearman_vs_main),
        by = "dataset"
    ) %>%
    dplyr::mutate(dataset = factor(dataset, levels = names(dataset_labels))) %>%
    dplyr::arrange(dataset) %>%
    dplyr::transmute(
        `ICI class` = unname(dataset_labels[as.character(dataset)]),
        `Cancer types with >=100 reports` = n_cancer_types_ge_threshold,
        `Median events per cancer-SOC cell` = median_events_per_cell,
        `Cells with 0 events (%)` = round(100 * prop_cells_0_events, 1),
        `Cells with <5 events (%)` = round(100 * prop_cells_lt5_events, 1),
        `Zero-event cells with ROR > 1` = zero_event_ror_above_1,
        `Median upper/lower 95% CI ratio of ROR` = round(median_ci_ratio, 2),
        `Median |log2 ROR difference| vs main analysis` = ifelse(dataset == "All ICI (main analysis)", NA, round(median_abs_log2_diff_vs_main, 2)),
        `Spearman correlation of log2 ROR with main analysis` = ifelse(dataset == "All ICI (main analysis)", NA, round(spearman_vs_main, 2))
    )


# 5. B表：各分层固定W投影结果 ####
projection_table <- stratum_summary_df %>%
    dplyr::mutate(
        ici_regimen = factor(ici_regimen, levels = names(dataset_labels)),
        sparse_handling = factor(sparse_handling, levels = c("sparse_to_1", "raw"))
    ) %>%
    dplyr::arrange(sparse_handling, dplyr::desc(min_cancer_n), ici_regimen) %>%
    dplyr::transmute(
        `Sparse-cell handling` = ifelse(sparse_handling == "sparse_to_1", "RORs of cells with <3 events set to 1 (primary)", "Original RORs"),
        `ICI class` = unname(dataset_labels[as.character(ici_regimen)]),
        `Inclusion threshold (reports)` = paste0("≥", min_cancer_n),
        `Cancer types included` = n_cancers,
        `Cancer types retaining final group` = paste0(n_cancers_retained, "/", n_cancers),
        `ARI vs final groups` = round(ari_vs_final, 3),
        `Variance explained by final basis patterns` = round(explained_variance_fixed_w, 3),
        `Cancer types reassigned (final -> projected group)` = ifelse(is.na(reassigned_cancers) | reassigned_cancers == "", "None", reassigned_cancers)
    )


# 6. C表：抗PD-1单药层重新NMF ####
denovo_pattern_wide <- denovo_match_df %>%
    dplyr::mutate(final_pattern = factor(final_pattern, levels = pattern_levels)) %>%
    dplyr::arrange(final_pattern) %>%
    dplyr::mutate(final_pattern = paste0("Cosine similarity: ", final_pattern)) %>%
    dplyr::select(sparse_handling, final_pattern, cosine_similarity) %>%
    dplyr::mutate(cosine_similarity = round(cosine_similarity, 3)) %>%
    tidyr::pivot_wider(names_from = final_pattern, values_from = cosine_similarity)

denovo_table <- denovo_summary_df %>%
    dplyr::left_join(denovo_pattern_wide, by = "sparse_handling") %>%
    dplyr::mutate(sparse_handling = factor(sparse_handling, levels = c("sparse_to_1", "raw"))) %>%
    dplyr::arrange(sparse_handling) %>%
    dplyr::transmute(
        `Sparse-cell handling` = ifelse(sparse_handling == "sparse_to_1", "RORs of cells with <3 events set to 1 (primary)", "Original RORs"),
        `Cancer types included` = n_cancers,
        `Cosine similarity: Ocular` = `Cosine similarity: Ocular`,
        `Cosine similarity: Hematologic` = `Cosine similarity: Hematologic`,
        `Cosine similarity: Neurologic` = `Cosine similarity: Neurologic`,
        `Cosine similarity: Systemic-like` = `Cosine similarity: Systemic-like`,
        `Cosine similarity: Visceral-like` = `Cosine similarity: Visceral-like`,
        `Cancer types retaining final group` = paste0(n_cancers_retained, "/", n_cancers),
        `ARI vs final groups` = round(ari_vs_final, 3),
        `Cancer types reassigned` = ifelse(is.na(reassigned_cancers) | reassigned_cancers == "", "None", reassigned_cancers)
    )


# 7. 结果保存 ####
writexl::write_xlsx(
    list(
        "各分层数据稀疏性" = sparsity_table,
        "各分层固定模式投影" = projection_table,
        "抗PD1重新NMF" = denovo_table
    ),
    path = output_xlsx_file
)

message("已生成补充表: ", output_xlsx_file)
message("========== ICI类别分层补充表整理结束 ==========")
