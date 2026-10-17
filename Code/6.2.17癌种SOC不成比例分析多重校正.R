# 癌种-SOC不成比例分析的多重校正（返修：审稿人2第二段，多重检验与效应量） ####
# 目的：原3.1分析未做多重校正（p_adjust_method = "none"），热图与补充表中的显著性按名义P < 0.05判定。
#      本脚本读取3.1已保存的FAERS、JADER主分析（M0）结果，在每个数据库的全部癌种-SOC检验内做BH校正，
#      输出加入BH校正P值的Table S4、S5，并汇总校正前后的显著关联数及未通过校正的关联。
# 不重新计算ROR，不改变NMF输入（NMF使用ROR数值，不依赖显著性）。

rm(list = ls())
gc()

# 1. 加载配置 ####
source("Code/config.R")
source("Code/mapping_rules.R")


# 2. 结果目录与输入文件 ####
results_dir_name <- "6.2.17癌种SOC不成比例分析多重校正"
script_results_path <- file.path(results_path, results_dir_name)

if (!dir.exists(script_results_path)) {
    dir.create(script_results_path, recursive = TRUE, showWarnings = FALSE)
}

input_files <- c(
    FAERS = file.path(results_path, "3.1primary_organ不成比例分析", "FAERS_primary_organ_SOC_M0.xlsx"),
    JADER = file.path(results_path, "3.1primary_organ不成比例分析", "JADER_primary_organ_SOC_M0.xlsx")
)
output_xlsx_file <- file.path(script_results_path, "癌种SOC不成比例分析多重校正结果.xlsx")

for (current_file in input_files) {
    if (!file.exists(current_file)) {
        stop("输入结果文件不存在: ", current_file)
    }
}

significance_threshold <- 0.05

message("========== 癌种SOC不成比例分析多重校正开始 ==========")


# 3. 数据读取与BH校正 ####
detail_list <- lapply(names(input_files), function(dataset_name) {
    readxl::read_xlsx(input_files[[dataset_name]], sheet = "详细结果") %>%
        dplyr::mutate(
            dataset = dataset_name,
            irae_abbrev = vapply(gsub("_", " ", irae), function(x) {
                current_abbrev <- get_soc_abbrev(x)
                if (is.na(current_abbrev) || current_abbrev == "") x else current_abbrev
            }, character(1)),
            p_value_BH = stats::p.adjust(p_value, method = "BH"),
            nominal_significant = p_value < significance_threshold,
            BH_significant = p_value_BH < significance_threshold,
            direction = ifelse(ROR > 1, "Higher reporting", "Lower reporting")
        )
})
names(detail_list) <- names(input_files)

# 原补充表中的显著性与名义P < 0.05一致性核查 ####
for (dataset_name in names(detail_list)) {
    if (!identical(detail_list[[dataset_name]]$signal_detected, detail_list[[dataset_name]]$nominal_significant)) {
        stop(dataset_name, "：原signal_detected与名义P < 0.05不一致，请核查")
    }
}


# 4. 统计汇总 ####

## 4.1 校正前后显著关联数 ####
summary_df <- dplyr::bind_rows(lapply(detail_list, function(df) {
    df %>%
        dplyr::summarise(
            dataset = dplyr::first(dataset),
            n_tests = dplyr::n(),
            n_nominal_significant = sum(nominal_significant),
            n_nominal_higher = sum(nominal_significant & ROR > 1),
            n_nominal_lower = sum(nominal_significant & ROR < 1),
            n_BH_significant = sum(BH_significant),
            n_BH_higher = sum(BH_significant & ROR > 1),
            n_BH_lower = sum(BH_significant & ROR < 1),
            n_lost_after_BH = sum(nominal_significant & !BH_significant),
            n_lost_higher = sum(nominal_significant & !BH_significant & ROR > 1),
            n_lost_lower = sum(nominal_significant & !BH_significant & ROR < 1),
            largest_nominal_p_retained = max(p_value[BH_significant])
        )
}))

## 4.2 未通过BH校正的名义显著关联 ####
lost_df <- dplyr::bind_rows(detail_list) %>%
    dplyr::filter(nominal_significant, !BH_significant) %>%
    dplyr::transmute(dataset, cancer_type, irae_abbrev, direction, a, ROR, ROR_lower, ROR_upper, p_value, p_value_BH) %>%
    dplyr::arrange(dataset, p_value)

message("校正前后显著关联数:")
print(as.data.frame(summary_df))


# 5. 补充表：在原Table S4、S5的列后加入BH校正P值 ####

## 函数：build_supplementary_table ####
# 功能：按原Table S4、S5的列名与顺序整理详细结果，并在最后加入BH校正P值 ####
# 输入：df，某一数据库加入BH校正后的详细结果 ####
# 输出：补充表数据框 ####
build_supplementary_table <- function(df) {
    df %>%
        dplyr::transmute(
            `Cancer types` = cancer_type,
            irAEs = irae_abbrev,
            `a (target cancer with target irAE)` = a,
            `b (target cancer without target irAE)` = b,
            `c (other cancers with target irAE)` = c,
            `d (other cancers without target irAE)` = d,
            ROR = ROR,
            ROR_lower = ROR_lower,
            ROR_upper = ROR_upper,
            p_value = p_value,
            p_value_BH = p_value_BH
        )
}

table_s4_df <- build_supplementary_table(detail_list$FAERS)
table_s5_df <- build_supplementary_table(detail_list$JADER)


# 6. 结果保存 ####
writexl::write_xlsx(
    list(
        "FAERS补充表_含BH" = table_s4_df,
        "JADER补充表_含BH" = table_s5_df,
        "校正前后显著关联数" = summary_df,
        "未通过BH的名义显著关联" = lost_df,
        "分析参数" = tibble::tibble(
            参数 = c("输入", "校正方法", "校正范围", "显著性阈值", "说明"),
            值 = c(
                "3.1primary_organ不成比例分析的FAERS、JADER M0（未调整）详细结果",
                "Benjamini–Hochberg",
                "每个数据库内全部癌种-SOC检验（各308个），双侧P值",
                "名义P < 0.05；BH校正P < 0.05",
                "原分析显著性（signal_detected）即名义P < 0.05，未做多重校正；本脚本不改变ROR及NMF输入"
            )
        )
    ),
    path = output_xlsx_file
)

message("已生成结果文件: ", output_xlsx_file)
message("========== 癌种SOC不成比例分析多重校正结束 ==========")
