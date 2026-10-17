# 致死结局与发病时间补充表（返修：审稿人2安全性结局） ####
# 读取6.2.12保存的结果表，整理为投稿用补充表（英文表头），不做新的统计分析

rm(list = ls())
gc()

# 1. 加载配置 ####
source("Code/config.R")

# 2. 结果目录与输入文件 ####
results_dir_name <- "6.2.14致死结局与发病时间补充表"
script_results_path <- file.path(results_path, results_dir_name)

if (!dir.exists(script_results_path)) {
    dir.create(script_results_path, recursive = TRUE, showWarnings = FALSE)
}

input_xlsx_file <- file.path(results_path, "6.2.12致死结局与发病时间不确定性分析", "致死结局与发病时间不确定性分析结果.xlsx")
output_xlsx_file <- file.path(script_results_path, "致死结局与发病时间补充表.xlsx")

if (!file.exists(input_xlsx_file)) {
    stop("输入结果文件不存在: ", input_xlsx_file)
}

pattern_levels <- c("Visceral-like", "Neurologic", "Ocular", "Systemic-like", "Hematologic")

# 6.2.12中的定义说明 → 补充表中的分析名称与分母说明 ####
definition_labels <- tibble::tribble(
    ~definition_prefix, ~analysis, ~fatal_definition, ~denominator,
    "FAERS primary", "Primary", "Report-level death (outcome code DE, any cause)", "All irAE reports, including reports without outcome codes",
    "JADER primary", "Primary", "Outcome of death for an irAE reaction", "All irAE reports, including reports with unknown irAE outcomes",
    "FAERS sensitivity", "Sensitivity: known outcomes only", "Report-level death (outcome code DE, any cause)", "Reports with at least one outcome code",
    "JADER sensitivity", "Sensitivity: known outcomes only", "Outcome of death for an irAE reaction", "Reports with a known outcome for at least one irAE reaction"
)

jader_aligned_label <- tibble::tibble(
    analysis = "Sensitivity: definition aligned with FAERS",
    fatal_definition = "Outcome of death for any reaction in the report",
    denominator = "Reports with a known outcome for at least one reaction"
)

level_order <- c("Whole database", "ICI-exposed solid-tumor reports", "irAE reports in the five groups")

message("========== 致死结局与发病时间补充表整理开始 ==========")


# 3. 数据读取 ####
group_fatality_df <- readxl::read_xlsx(input_xlsx_file, sheet = "各组致死比例")
database_level_df <- readxl::read_xlsx(input_xlsx_file, sheet = "数据库层面死亡比例")
tto_df <- readxl::read_xlsx(input_xlsx_file, sheet = "FAERS发病时间有效例数")


# 4. 表1：各癌种群致死比例（主分析与敏感性分析） ####
# 主分析的n_unknown已因二分类而归零；缺失数量应取相同结局定义的已知结局分析，不能把纳入分母误写成无缺失。
missing_outcome_df <- group_fatality_df %>%
    dplyr::filter(grepl("sensitivity:", definition), !grepl("aligned with FAERS", definition)) %>%
    dplyr::select(dataset, pattern_group, n_original_unknown = n_unknown, pct_original_unknown = pct_unknown)

stopifnot(nrow(missing_outcome_df) == 10L,
          !anyDuplicated(missing_outcome_df[c("dataset", "pattern_group")]))

fatality_table <- group_fatality_df %>%
    dplyr::mutate(definition_prefix = sub(":.*$", "", definition)) %>%
    dplyr::left_join(definition_labels, by = "definition_prefix") %>%
    dplyr::left_join(missing_outcome_df, by = c("dataset", "pattern_group")) %>%
    dplyr::mutate(
        n_unknown = ifelse(grepl("primary:", definition), n_original_unknown, n_unknown),
        pct_unknown = ifelse(grepl("primary:", definition), pct_original_unknown, pct_unknown),
        analysis = ifelse(grepl("aligned with FAERS", definition), jader_aligned_label$analysis, analysis),
        fatal_definition = ifelse(grepl("aligned with FAERS", definition), jader_aligned_label$fatal_definition, fatal_definition),
        denominator = ifelse(grepl("aligned with FAERS", definition), jader_aligned_label$denominator, denominator),
        analysis = factor(analysis, levels = c("Primary", "Sensitivity: known outcomes only", "Sensitivity: definition aligned with FAERS")),
        dataset = factor(dataset, levels = c("FAERS", "JADER")),
        pattern_group = factor(pattern_group, levels = pattern_levels)
    ) %>%
    dplyr::arrange(analysis, dataset, pattern_group) %>%
    dplyr::transmute(
        Analysis = as.character(analysis),
        Database = as.character(dataset),
        `Fatal outcome` = fatal_definition,
        Denominator = denominator,
        Group = as.character(pattern_group),
        `irAE reports` = n_reports,
        `Reports with unknown outcome, n (%)` = paste0(n_unknown, " (", sprintf("%.1f", pct_unknown), ")"),
        `Reports in denominator` = n_known,
        `Fatal reports` = n_fatal,
        `Fatal outcome reporting, % (95% CI)` = paste0(sprintf("%.1f", pct_fatal_known), " (", sprintf("%.1f", ci_low_pct), "–", sprintf("%.1f", ci_high_pct), ")")
    )


# 5. 表2：数据库层面死亡比例（全库、ICI实体瘤报告、irAE队列） ####
database_table <- database_level_df %>%
    dplyr::mutate(
        level = factor(level, levels = level_order),
        dataset = factor(dataset, levels = c("FAERS", "JADER"))
    ) %>%
    dplyr::arrange(level, dataset) %>%
    dplyr::transmute(
        `Report set` = as.character(level),
        Database = as.character(dataset),
        `Fatal outcome` = ifelse(dataset == "FAERS", "Report-level death (outcome code DE, any cause)", "Outcome of death for any reaction in the report"),
        Reports = n_reports,
        `Reports with unknown outcome, n (%)` = paste0(n_unknown, " (", sprintf("%.1f", pct_unknown), ")"),
        `Fatal reports` = n_fatal,
        `Reports with known outcome` = n_known,
        `Fatal outcome reporting among all reports, %` = round(pct_fatal_all_reports, 1),
        `Fatal outcome reporting among reports with known outcome, % (95% CI)` = paste0(sprintf("%.1f", pct_fatal_known), " (", sprintf("%.1f", ci_low_pct), "–", sprintf("%.1f", ci_high_pct), ")")
    )


# 6. 表3：FAERS发病时间的有效例数 ####
tto_table <- tto_df %>%
    dplyr::mutate(pattern_group = factor(pattern_group, levels = c(pattern_levels, "All five groups"))) %>%
    dplyr::arrange(pattern_group) %>%
    dplyr::transmute(
        Group = as.character(pattern_group),
        `irAE reports` = n_irae_reports,
        `Reports with missing or invalid onset intervals` = n_missing_or_invalid_dates,
        `Reports with valid onset intervals, n (%)` = paste0(n_valid_dates, " (", sprintf("%.1f", pct_valid_dates), ")"),
        `Excluded above the 99th percentile` = n_above_99th_percentile,
        `Reports analysed` = n_analysed,
        `Time to onset, median (IQR), days` = paste0(median_days, " (", q1_days, "–", q3_days, ")")
    )


# 7. 结果保存 ####
writexl::write_xlsx(
    list(
        "各组致死比例" = fatality_table,
        "数据库层面死亡比例" = database_table,
        "FAERS发病时间有效例数" = tto_table
    ),
    path = output_xlsx_file
)

message("已生成补充表: ", output_xlsx_file)
message("========== 致死结局与发病时间补充表整理结束 ==========")
