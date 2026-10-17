# ICI类别分层数据稀疏性分析（返修：审稿人1第12条） ####
# 目的：定量说明按ICI类别分层后，各癌种-SOC单元的事件数被稀释、ROR估计精度下降，
#      从而说明分层模式分析在小分层中不可靠。
#      1) 比较主分析与各分层中每个癌种-SOC单元的事件数、0事件单元比例、事件数<5单元比例；
#      2) 比较95%置信区间宽度（上限/下限之比）；
#      3) 比较分层ROR与主分析ROR的一致性，并按事件数分组；
#      4) 统计0事件但ROR>1的单元。
# 只读取3.1与3.3已保存的结果，不重新计算ROR。

rm(list = ls())
gc()

# 1. 加载配置 ####
source("Code/config.R")
source("Code/mapping_rules.R")

# 2. 结果目录与输入文件 ####
results_dir_name <- "6.2.9ICI类别分层数据稀疏性分析"
script_results_path <- file.path(results_path, results_dir_name)

if (!dir.exists(script_results_path)) {
    dir.create(script_results_path, recursive = TRUE, showWarnings = FALSE)
}

main_detail_file <- file.path(results_path, "3.1primary_organ不成比例分析", "FAERS_primary_organ_SOC_M0.xlsx")
regimen_detail_file <- file.path(results_path, "3.3同一ICI_regimen下primary_organ不成比例分析", "FAERS_within_ICI_regimen_primary_organ_SOC_M0.xlsx")
main_wh_matrix_file <- file.path(results_path, "6.1癌种NMF分群分析_rank5", "FAERS_NMF_WH_matrix.xlsx")
output_xlsx_file <- file.path(script_results_path, "ICI类别分层数据稀疏性分析结果.xlsx")

stratum_min_cancer_n <- 100L
event_bins <- c(-Inf, 0, 4, 19, Inf)
event_bin_labels <- c("0", "1-4", "5-19", ">=20")
dataset_levels <- c("All ICI (main analysis)", "anti-PD-1", "anti-PD-L1", "anti-CTLA-4+anti-PD-1")

message("========== ICI类别分层数据稀疏性分析开始 ==========")


# 3. 辅助函数 ####

## 函数：build_params_dataframe ####
# 功能：将参数列表转换为可写入Excel的数据框####
# 输入：params_list，命名列表####
# 输出：包含参数和值两列的数据框####
build_params_dataframe <- function(params_list) {
    data.frame(
        参数 = names(params_list),
        值 = vapply(params_list, function(param_value) {
            paste(as.character(param_value), collapse = ", ")
        }, character(1)),
        stringsAsFactors = FALSE
    )
}


# 4. 读取单元层面的事件数与ROR ####
final_cancers <- readxl::read_xlsx(main_wh_matrix_file, sheet = "Pattern group assignment") %>%
    dplyr::pull(cancer_type) %>%
    as.character()

main_cells_df <- readxl::read_xlsx(main_detail_file, sheet = "详细结果") %>%
    dplyr::filter(cancer_type %in% final_cancers) %>%
    dplyr::transmute(
        dataset = "All ICI (main analysis)",
        cancer_type = as.character(cancer_type),
        soc_abbrev = vapply(gsub("_", " ", irae), function(x) {
            current_abbrev <- get_soc_abbrev(x)
            if (is.na(current_abbrev) || current_abbrev == "") x else current_abbrev
        }, character(1)),
        n_reports = as.numeric(n_cancer),
        n_events = as.numeric(a),
        ror = as.numeric(ROR),
        ror_lower = as.numeric(ROR_lower),
        ror_upper = as.numeric(ROR_upper)
    )

stratum_cells_df <- readxl::read_xlsx(regimen_detail_file, sheet = "详细结果") %>%
    dplyr::filter(primary_organ %in% final_cancers) %>%
    dplyr::transmute(
        dataset = as.character(ici_regimen_raw),
        cancer_type = as.character(primary_organ),
        soc_abbrev = as.character(irae_abbrev),
        n_reports = as.numeric(organ_n),
        n_events = as.numeric(a),
        ror = as.numeric(ROR),
        ror_lower = as.numeric(ROR_lower),
        ror_upper = as.numeric(ROR_upper)
    )

cells_df <- dplyr::bind_rows(main_cells_df, stratum_cells_df) %>%
    dplyr::mutate(
        dataset = factor(dataset, levels = dataset_levels),
        included_in_pattern_analysis = n_reports >= stratum_min_cancer_n,
        event_bin = cut(n_events, breaks = event_bins, labels = event_bin_labels),
        ci_ratio = ror_upper / ror_lower,
        zero_event_ror_above_1 = n_events == 0 & ror > 1
    ) %>%
    dplyr::left_join(
        main_cells_df %>% dplyr::select(cancer_type, soc_abbrev, main_ror = ror, main_n_events = n_events),
        by = c("cancer_type", "soc_abbrev")
    ) %>%
    dplyr::mutate(abs_log2_diff_vs_main = abs(log2(ror) - log2(main_ror)))


# 5. 汇总：各数据集的单元事件数与ROR精度 ####

## 5.1 纳入模式分析的癌种（报告数≥阈值）####
precision_summary_df <- cells_df %>%
    dplyr::filter(included_in_pattern_analysis) %>%
    dplyr::group_by(dataset) %>%
    dplyr::summarise(
        n_cancer_types = dplyr::n_distinct(cancer_type),
        n_cells = dplyr::n(),
        median_reports_per_cancer = stats::median(n_reports[!duplicated(cancer_type)]),
        median_events_per_cell = stats::median(n_events),
        prop_cells_0_events = mean(n_events == 0),
        prop_cells_lt5_events = mean(n_events < 5),
        n_zero_event_ror_above_1 = sum(zero_event_ror_above_1),
        median_ci_ratio = stats::median(ci_ratio),
        q3_ci_ratio = stats::quantile(ci_ratio, 0.75, names = FALSE),
        median_abs_log2_diff_vs_main = stats::median(abs_log2_diff_vs_main),
        spearman_vs_main = suppressWarnings(stats::cor(log2(ror), log2(main_ror), method = "spearman")),
        .groups = "drop"
    ) %>%
    dplyr::mutate(scope = paste0("Cancer types with >=", stratum_min_cancer_n, " reports"), .before = 1)

## 5.2 全部25个癌种（显示稀释程度）####
dilution_summary_df <- cells_df %>%
    dplyr::group_by(dataset) %>%
    dplyr::summarise(
        n_cancer_types = dplyr::n_distinct(cancer_type),
        n_cancer_types_ge_threshold = dplyr::n_distinct(cancer_type[included_in_pattern_analysis]),
        n_cells = dplyr::n(),
        median_events_per_cell = stats::median(n_events),
        prop_cells_0_events = mean(n_events == 0),
        prop_cells_lt5_events = mean(n_events < 5),
        n_zero_event_ror_above_1 = sum(zero_event_ror_above_1),
        median_ci_ratio = stats::median(ci_ratio),
        .groups = "drop"
    ) %>%
    dplyr::mutate(scope = "All 25 cancer types", .before = 1)

## 5.3 分层ROR与主分析ROR的偏离程度，按分层内事件数分组 ####
deviation_by_events_df <- cells_df %>%
    dplyr::filter(dataset != "All ICI (main analysis)") %>%
    dplyr::group_by(dataset, event_bin) %>%
    dplyr::summarise(
        n_cells = dplyr::n(),
        median_abs_log2_diff_vs_main = stats::median(abs_log2_diff_vs_main),
        median_ci_ratio = stats::median(ci_ratio),
        prop_ror_above_1 = mean(ror > 1),
        .groups = "drop"
    )

message("纳入模式分析的癌种（报告数>=", stratum_min_cancer_n, "）:")
print(as.data.frame(precision_summary_df), digits = 3)
message("全部25个癌种:")
print(as.data.frame(dilution_summary_df), digits = 3)
message("分层ROR与主分析ROR的偏离（按事件数）:")
print(as.data.frame(deviation_by_events_df), digits = 3)


# 6. 结果保存 ####
writexl::write_xlsx(
    list(
        "精度汇总_纳入癌种" = precision_summary_df,
        "稀释汇总_全部癌种" = dilution_summary_df,
        "与主分析偏离_按事件数" = deviation_by_events_df,
        "单元明细" = cells_df %>% dplyr::mutate(dataset = as.character(dataset), event_bin = as.character(event_bin)),
        "分析参数" = build_params_dataframe(list(
            main_source = main_detail_file,
            stratum_source = regimen_detail_file,
            stratum_min_cancer_n = stratum_min_cancer_n,
            event_bins = paste(event_bin_labels, collapse = ", "),
            ci_ratio = "95%置信区间上限/下限",
            deviation = "|log2(分层ROR) - log2(主分析ROR)|"
        ))
    ),
    path = output_xlsx_file
)

message("已生成结果文件: ", output_xlsx_file)
message("========== ICI类别分层数据稀疏性分析结束 ==========")
