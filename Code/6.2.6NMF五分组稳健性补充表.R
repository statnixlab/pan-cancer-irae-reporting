# NMF五分组稳健性补充表（返修：审稿人2第一段） ####
# 读取6.2.4保存的结果表，整理为投稿用补充表（英文表头），不做新的统计分析
# 条件名称与补充图S15、S16保持一致

rm(list = ls())
gc()

# 1. 加载配置 ####
source("Code/config.R")

# 2. 结果目录与输入文件 ####
results_dir_name <- "6.2.6NMF五分组稳健性补充表"
script_results_path <- file.path(results_path, results_dir_name)

if (!dir.exists(script_results_path)) {
    dir.create(script_results_path, recursive = TRUE, showWarnings = FALSE)
}

input_xlsx_file <- file.path(results_path, "6.2.4NMF五分组稳健性分析", "NMF五分组稳健性分析结果.xlsx")
output_xlsx_file <- file.path(script_results_path, "NMF五分组稳健性补充表.xlsx")

if (!file.exists(input_xlsx_file)) {
    stop("6.2.4结果文件不存在，请先运行 Code/6.2.4NMF五分组稳健性分析.R: ", input_xlsx_file)
}

reference_group_levels <- c("Ocular", "Hematologic", "Neurologic", "Systemic-like", "Visceral-like")
analysis_type_levels <- c(
    "Random seed",
    "Profile-preserving transformation",
    "Inclusion threshold",
    "Leave one cancer out",
    "Excluding uncertain cancers",
    "Profile-altering transformation"
)

message("========== NMF五分组稳健性补充表整理开始 ==========")


# 3. 数据读取 ####
cancer_stability_df <- readxl::read_xlsx(input_xlsx_file, sheet = "癌种归属稳定性")
cancer_stability_by_type_df <- readxl::read_xlsx(input_xlsx_file, sheet = "癌种归属稳定性_按检验类别")
run_summary_df <- readxl::read_xlsx(input_xlsx_file, sheet = "稳健性运行汇总")
exclusion_summary_df <- readxl::read_xlsx(input_xlsx_file, sheet = "剔除不确定癌种汇总")


# 4. A表：各癌种的归属稳定性 ####
stability_by_type_wide <- cancer_stability_by_type_df %>%
    dplyr::mutate(
        retained = paste0(n_same_as_final, "/", n_runs),
        condition_type = factor(condition_type, levels = analysis_type_levels)
    ) %>%
    dplyr::arrange(condition_type) %>%
    dplyr::select(cancer_type, condition_type, retained) %>%
    tidyr::pivot_wider(names_from = condition_type, values_from = retained)

cancer_table <- cancer_stability_df %>%
    dplyr::select(cancer_type, reference_group, final_silhouette, dominant_h_share) %>%
    dplyr::left_join(stability_by_type_wide, by = "cancer_type") %>%
    dplyr::mutate(reference_group = factor(reference_group, levels = reference_group_levels)) %>%
    dplyr::arrange(reference_group, cancer_type) %>%
    dplyr::transmute(
        `Cancer type` = cancer_type,
        `Final group` = as.character(reference_group),
        `Silhouette width (final model)` = round(final_silhouette, 3),
        `Dominant-pattern share` = round(dominant_h_share, 3),
        `Random seed` = `Random seed`,
        `Profile-preserving transformation` = `Profile-preserving transformation`,
        `Inclusion threshold` = `Inclusion threshold`,
        `Leave one cancer out` = `Leave one cancer out`,
        `Profile-altering transformation` = `Profile-altering transformation`
    )


# 5. B表：各次稳健性分析汇总 ####
analysis_table <- dplyr::bind_rows(
    run_summary_df %>%
        dplyr::select(condition_id, condition_type, condition_label, n_cancers, ari_vs_final, n_cancers_reassigned, reassigned_cancers, status),
    exclusion_summary_df %>%
        dplyr::mutate(condition_type = "Excluding uncertain cancers") %>%
        dplyr::select(condition_id, condition_type, condition_label, n_cancers, ari_vs_final, n_cancers_reassigned, reassigned_cancers, status, excluded_cancers)
) %>%
    dplyr::filter(status == "completed") %>%
    dplyr::mutate(
        analysis_label = dplyr::case_when(
            condition_id == "transform_log2_minmax" ~ "log2 ROR, min-max",
            condition_id == "transform_max_scaled" ~ "ROR / maximum",
            condition_id == "transform_zscore_shifted" ~ "z-score",
            condition_id == "transform_log2_channel" ~ "log2 ROR, two-channel",
            condition_id == "transform_log2_channel_sigonly" ~ "log2 ROR, sig. only",
            condition_id == "transform_rank" ~ "Within-cancer rank",
            condition_id == "exclude_low_silhouette" ~ paste0("Excluding 5 lowest-silhouette cancer types (", excluded_cancers, ")"),
            condition_type == "Inclusion threshold" ~ sub("^>= ", "≥", condition_label),
            TRUE ~ condition_label
        ),
        condition_type = factor(condition_type, levels = analysis_type_levels),
        sort_key = suppressWarnings(as.numeric(gsub("\\D", "", condition_id)))
    ) %>%
    dplyr::arrange(condition_type, sort_key, condition_id) %>%
    dplyr::transmute(
        `Type of analysis` = as.character(condition_type),
        `Analysis` = analysis_label,
        `Cancer types included` = n_cancers,
        `ARI vs final groups` = round(ari_vs_final, 3),
        `Cancer types reassigned (n)` = n_cancers_reassigned,
        `Cancer types reassigned` = ifelse(is.na(reassigned_cancers) | reassigned_cancers == "", "None", reassigned_cancers)
    )


# 6. 结果保存 ####
writexl::write_xlsx(
    list(
        "各癌种归属稳定性" = cancer_table,
        "各次稳健性分析汇总" = analysis_table
    ),
    path = output_xlsx_file
)

message("已生成补充表: ", output_xlsx_file)
message("========== NMF五分组稳健性补充表整理结束 ==========")
