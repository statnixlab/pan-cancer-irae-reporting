# NMF秩选择补充表（返修：审稿人1第1条、审稿人2秩选择） ####
# 读取6.2.1保存的结果表，整理为投稿用补充表（英文表头），不做新的统计分析

rm(list = ls())
gc()

# 1. 加载配置 ####
source("Code/config.R")

# 2. 结果目录与输入文件 ####
results_dir_name <- "6.2.3NMF秩选择补充表"
script_results_path <- file.path(results_path, results_dir_name)

if (!dir.exists(script_results_path)) {
    dir.create(script_results_path, recursive = TRUE, showWarnings = FALSE)
}

input_xlsx_file <- file.path(results_path, "6.2.1NMF秩选择敏感性分析", "NMF秩选择敏感性分析结果.xlsx")
output_xlsx_file <- file.path(script_results_path, "NMF秩选择补充表.xlsx")

if (!file.exists(input_xlsx_file)) {
    stop("6.2.1结果文件不存在，请先运行 Code/6.2.1NMF秩选择敏感性分析.R: ", input_xlsx_file)
}

message("========== NMF秩选择补充表整理开始 ==========")


# 3. 数据读取 ####
quality_df <- readxl::read_xlsx(input_xlsx_file, sheet = "各秩聚类质量汇总")
rss_drop_df <- readxl::read_xlsx(input_xlsx_file, sheet = "RSS下降_实际与零分布")
ari_df <- readxl::read_xlsx(input_xlsx_file, sheet = "分组一致性ARI") %>%
    dplyr::filter(comparison == "vs final rank-5 groups (500 runs)")
w_long_df <- readxl::read_xlsx(input_xlsx_file, sheet = "各秩W矩阵")


# 4. 秩选择指标汇总表 ####
rank_metrics_table <- quality_df %>%
    dplyr::left_join(rss_drop_df, by = "rank") %>%
    dplyr::left_join(ari_df %>% dplyr::select(rank, ari), by = "rank") %>%
    dplyr::arrange(rank) %>%
    dplyr::transmute(
        `Rank` = rank,
        `Number of clusters` = n_clusters,
        `Minimum cluster size` = min_cluster_size,
        `Singleton clusters` = n_singleton_clusters,
        `Cophenetic correlation` = round(cophenetic, 3),
        `Dispersion` = round(dispersion, 3),
        `Mean silhouette` = round(mean_silhouette, 3),
        `Minimum silhouette` = round(min_silhouette, 3),
        `Cancer types with silhouette < 0.5` = n_cancers_silhouette_below_threshold,
        `RSS` = round(rss, 3),
        `Explained variance` = round(evar, 3),
        `Decrease in RSS from previous rank (observed)` = round(rss_drop_observed, 3),
        `Decrease in RSS from previous rank (permutation mean)` = round(rss_drop_null_mean, 3),
        `Decrease in RSS from previous rank (permutation minimum)` = round(rss_drop_null_min, 3),
        `Decrease in RSS from previous rank (permutation maximum)` = round(rss_drop_null_max, 3),
        `Excess decrease in RSS (observed - permutation mean)` = round(observed_minus_null_mean, 3),
        `ARI vs final rank-5 groups` = round(ari, 3)
    )

## 4.1 转置：行为指标，列为各秩 ####
rank_metrics_table <- rank_metrics_table %>%
    tidyr::pivot_longer(-Rank, names_to = "Metric", values_to = "value") %>%
    dplyr::mutate(
        Metric = factor(Metric, levels = setdiff(names(rank_metrics_table), "Rank")),
        Rank = paste0("Rank ", Rank)
    ) %>%
    tidyr::pivot_wider(names_from = Rank, values_from = value) %>%
    dplyr::arrange(Metric) %>%
    dplyr::mutate(Metric = as.character(Metric))


# 5. 各秩W矩阵宽表 ####
soc_column_order <- sort(unique(w_long_df$soc_abbrev))

w_matrix_table <- w_long_df %>%
    dplyr::mutate(w_value = round(w_value, 3)) %>%
    tidyr::pivot_wider(
        id_cols = c(rank, pattern, pattern_label),
        names_from = soc_abbrev,
        values_from = w_value
    ) %>%
    dplyr::mutate(pattern_order = as.integer(sub("^P", "", pattern))) %>%
    dplyr::arrange(rank, pattern_order) %>%
    dplyr::select(
        `Rank` = rank,
        `Pattern` = pattern,
        `Dominant SOCs` = pattern_label,
        dplyr::all_of(soc_column_order)
    )


# 6. 结果保存 ####
writexl::write_xlsx(
    list(
        "秩选择指标汇总" = rank_metrics_table,
        "各秩W矩阵" = w_matrix_table
    ),
    path = output_xlsx_file
)

message("已生成补充表: ", output_xlsx_file)
message("========== NMF秩选择补充表整理结束 ==========")
