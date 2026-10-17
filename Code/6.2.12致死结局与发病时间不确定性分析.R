# 致死结局与发病时间的不确定性及缺失分析（返修：审稿人2安全性结局） ####
# 目的：
#      1) 各癌种群致死比例给出分母、95%置信区间和结局缺失情况；主分析两库分母统一为全部irAE报告（结局缺失计为非致死），
#         原稿6.1.2中FAERS排除了无结局编码的报告、JADER将“不明”计为非致死，两库分母不一致；
#      2) 敏感性分析：两库只用结局已知的报告，以及JADER与FAERS定义对齐的报告层面致死；
#      3) 在全库、ICI实体瘤报告和irAE队列三个层面比较FAERS与JADER的死亡比例，说明两库差异来自数据库结构与结局定义；
#      4) 发病时间给出各组有效日期例数与排除例数。
# FAERS的OUTC表为报告层面结局（DE表示患者死亡，不限死因），JADER的转归为反应层面结局。
# 只读取已清洗的数据与6.1.1分群结果，不修改6.1.2的结果。

rm(list = ls())
gc()

# 1. 加载配置 ####
source("Code/config.R")
source("Code/mapping_rules.R")


# 2. 结果目录与参数 ####
results_dir_name <- "6.2.12致死结局与发病时间不确定性分析"
script_results_path <- file.path(results_path, results_dir_name)

if (!dir.exists(script_results_path)) {
    dir.create(script_results_path, recursive = TRUE, showWarnings = FALSE)
}

pattern_assignment_file <- file.path(results_path, "6.1癌种NMF分群分析_rank5", "FAERS_NMF_WH_matrix.xlsx")
faers_demo_file <- file.path(output_path, "faers_demo.RData")
faers_outc_file <- file.path(output_path, "faers_outc.RData")
faers_solid_file <- file.path(output_path, "faers_solid.RData")
jader_demo_file <- file.path(output_path, "jader_demo.RData")
jader_reac_file <- file.path(output_path, "jader_reac.RData")
jader_solid_file <- file.path(output_path, "jader_solid.RData")
irae_dict_file <- file.path(data_path, "irAE列表.xlsx")

output_xlsx_file <- file.path(script_results_path, "致死结局与发病时间不确定性分析结果.xlsx")

pattern_display_order <- c("Ocular", "Hematologic", "Neurologic", "Systemic-like", "Visceral-like")

# JADER转归：已知转归按严重程度排序，“不明/未知/空值”视为结局未知 ####
jader_known_outcomes <- c("死亡", "後遺症あり", "未回復", "軽快", "回復")
jader_outcome_labels <- c(
    "死亡" = "Death",
    "後遺症あり" = "Sequelae",
    "未回復" = "Not recovered",
    "軽快" = "Improving",
    "回復" = "Recovered",
    "Unknown" = "Unknown"
)

faers_outcome_labels <- c(
    "DE" = "Death",
    "LT" = "Life-threatening",
    "HO" = "Hospitalization",
    "DS" = "Disability",
    "CA" = "Congenital anomaly",
    "RI" = "Required intervention",
    "OT" = "Other serious"
)

level_labels <- c(
    "whole" = "Whole database",
    "ici_solid" = "ICI-exposed solid-tumor reports",
    "irae_cohort" = "irAE reports in the five groups"
)

message("========== 致死结局与发病时间不确定性分析开始 ==========")
message("结果目录: ", script_results_path)


# 3. 辅助函数 ####

## 函数：load_named_rdata_object ####
# 功能：从RData文件中读取指定对象 ####
# 输入：input_path，RData文件路径；object_name，对象名称 ####
# 输出：返回读取出的对象 ####
load_named_rdata_object <- function(input_path, object_name) {
    input_env <- new.env(parent = emptyenv())
    load(input_path, envir = input_env)

    if (!exists(object_name, envir = input_env, inherits = FALSE)) {
        stop("RData中未找到对象: ", object_name)
    }

    input_env[[object_name]]
}

## 函数：standardize_medra_term ####
# 功能：标准化MedDRA术语文本（与6.1.2一致），便于与irAE词典匹配 ####
# 输入：x，字符型术语向量 ####
# 输出：大写的标准化字符向量 ####
standardize_medra_term <- function(x) {
    x <- as.character(x)
    x <- stringr::str_replace_all(x, "_", " ")
    x <- stringr::str_replace_all(x, "-", " ")
    x <- stringr::str_replace_all(x, "\\s+", " ")
    x <- trimws(x)
    toupper(x)
}

## 函数：standardize_jader_outcome_known ####
# 功能：将JADER转归中的“不明/未知/空值”统一记为Unknown，已知转归保留原值 ####
# 输入：x，字符型转归向量 ####
# 输出：字符向量，取值为已知转归或"Unknown" ####
standardize_jader_outcome_known <- function(x) {
    x <- trimws(as.character(x))
    ifelse(x %in% jader_known_outcomes, x, "Unknown")
}

## 函数：calc_wilson_ci ####
# 功能：计算比例的Wilson 95%置信区间 ####
# 输入：x，事件数；n，分母 ####
# 输出：含下限与上限（百分比）的数据框 ####
calc_wilson_ci <- function(x, n) {
    z <- stats::qnorm(0.975)
    p <- ifelse(n > 0, x / n, NA_real_)
    denominator <- 1 + z^2 / n
    centre <- (p + z^2 / (2 * n)) / denominator
    half_width <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denominator

    tibble::tibble(
        ci_low_pct = ifelse(n > 0, (centre - half_width) * 100, NA_real_),
        ci_high_pct = ifelse(n > 0, (centre + half_width) * 100, NA_real_)
    )
}

## 函数：load_pattern_group_assignment ####
# 功能：读取6.1.1导出的癌种群分配表 ####
# 输入：input_file，分群Excel路径 ####
# 输出：primary_organ与pattern_group对照表 ####
load_pattern_group_assignment <- function(input_file) {
    assignment_df <- readxl::read_xlsx(input_file, sheet = "Pattern group assignment")
    names(assignment_df)[1:2] <- c("primary_organ", "pattern_group")

    assignment_df %>%
        dplyr::transmute(primary_organ = as.character(primary_organ), pattern_group = as.character(pattern_group)) %>%
        dplyr::filter(!is.na(primary_organ), primary_organ != "") %>%
        dplyr::distinct(primary_organ, .keep_all = TRUE)
}

## 函数：summarise_group_fatality ####
# 功能：按癌种群汇总报告层面的致死状态（Fatal/Non-fatal/Unknown），给出分母、缺失比例、致死比例及Wilson 95%CI ####
# 输入：status_df，含pattern_group与fatality_status的报告层面数据框；dataset_name，数据库；definition_name，致死定义说明 ####
# 输出：每组一行的汇总数据框 ####
summarise_group_fatality <- function(status_df, dataset_name, definition_name) {
    summary_df <- tibble::tibble(pattern_group = pattern_display_order) %>%
        dplyr::left_join(
            status_df %>%
                dplyr::group_by(pattern_group) %>%
                dplyr::summarise(
                    n_reports = dplyr::n(),
                    n_fatal = sum(fatality_status == "Fatal"),
                    n_non_fatal = sum(fatality_status == "Non-fatal"),
                    n_unknown = sum(fatality_status == "Unknown"),
                    .groups = "drop"
                ),
            by = "pattern_group"
        ) %>%
        dplyr::mutate(
            n_known = n_fatal + n_non_fatal,
            pct_unknown = n_unknown / n_reports * 100,
            pct_fatal_known = n_fatal / n_known * 100,
            pct_fatal_all_reports = n_fatal / n_reports * 100
        )

    dplyr::bind_cols(summary_df, calc_wilson_ci(summary_df$n_fatal, summary_df$n_known)) %>%
        dplyr::mutate(
            dataset = dataset_name,
            definition = definition_name,
            rank_fatal_known = rank(-pct_fatal_known, ties.method = "min"),
            .before = 1
        ) %>%
        dplyr::select(
            dataset, definition, pattern_group, n_reports, n_known, n_unknown, pct_unknown,
            n_fatal, pct_fatal_known, ci_low_pct, ci_high_pct, pct_fatal_all_reports, rank_fatal_known
        )
}

## 函数：test_group_fatality ####
# 功能：在结局已知的报告中，用卡方检验比较各癌种群的致死比例 ####
# 输入：status_df，报告层面致死状态数据框；dataset_name，数据库；definition_name，致死定义说明 ####
# 输出：一行的检验结果数据框 ####
test_group_fatality <- function(status_df, dataset_name, definition_name) {
    known_df <- status_df %>% dplyr::filter(fatality_status != "Unknown")
    test_result <- stats::chisq.test(table(known_df$pattern_group, known_df$fatality_status))

    tibble::tibble(
        dataset = dataset_name,
        definition = definition_name,
        n_known = nrow(known_df),
        chi_squared = unname(test_result$statistic),
        df = unname(test_result$parameter),
        p_value = test_result$p.value
    )
}

## 函数：test_pairwise_fatality ####
# 功能：在结局已知的报告中，用pairwise.prop.test（Bonferroni校正，与6.1.2一致）两两比较各癌种群的致死比例 ####
# 输入：status_df，报告层面致死状态数据框；dataset_name，数据库；definition_name，致死定义说明 ####
# 输出：每对比较一行的数据框 ####
test_pairwise_fatality <- function(status_df, dataset_name, definition_name) {
    known_df <- status_df %>%
        dplyr::filter(fatality_status != "Unknown") %>%
        dplyr::mutate(pattern_group = factor(pattern_group, levels = pattern_display_order))
    group_counts <- known_df %>%
        dplyr::group_by(pattern_group) %>%
        dplyr::summarise(n_fatal = sum(fatality_status == "Fatal"), n_known = dplyr::n(), .groups = "drop")
    pairwise_result <- stats::pairwise.prop.test(
        x = stats::setNames(group_counts$n_fatal, group_counts$pattern_group),
        n = stats::setNames(group_counts$n_known, group_counts$pattern_group),
        p.adjust.method = "bonferroni"
    )

    as.data.frame(as.table(pairwise_result$p.value), stringsAsFactors = FALSE) %>%
        dplyr::filter(!is.na(Freq)) %>%
        dplyr::transmute(
            dataset = dataset_name,
            definition = definition_name,
            group1 = Var1,
            group2 = Var2,
            p_value_bonferroni = Freq,
            significance = dplyr::case_when(
                Freq < 0.001 ~ "***",
                Freq < 0.01 ~ "**",
                Freq < 0.05 ~ "*",
                TRUE ~ "ns"
            )
        )
}

## 函数：summarise_faers_level ####
# 功能：汇总FAERS某一层面报告的死亡比例（报告层面DE）及结局编码缺失情况 ####
# 输入：record_ids，该层面的primaryid；outc_df，FAERS OUTC表；level_name，层面名称 ####
# 输出：一行的汇总数据框 ####
summarise_faers_level <- function(record_ids, outc_df, level_name) {
    record_ids <- unique(as.character(record_ids))
    level_outc <- outc_df %>% dplyr::filter(primaryid %in% record_ids)
    ids_any_outc <- unique(level_outc$primaryid)
    ids_death <- unique(level_outc$primaryid[level_outc$outc_cod == "DE"])

    n_reports <- length(record_ids)
    n_known <- length(ids_any_outc)
    n_fatal <- length(ids_death)
    ci_df <- calc_wilson_ci(n_fatal, n_known)

    tibble::tibble(
        level = level_name,
        dataset = "FAERS",
        definition = "Report-level death (OUTC = DE, any cause)",
        n_reports = n_reports,
        n_known = n_known,
        n_unknown = n_reports - n_known,
        pct_unknown = (n_reports - n_known) / n_reports * 100,
        n_fatal = n_fatal,
        pct_fatal_known = n_fatal / n_known * 100,
        ci_low_pct = ci_df$ci_low_pct,
        ci_high_pct = ci_df$ci_high_pct,
        pct_fatal_all_reports = n_fatal / n_reports * 100
    )
}

## 函数：summarise_jader_level ####
# 功能：汇总JADER某一层面报告的死亡比例（任一反应转归为死亡即计为死亡），以及反应层面的转归构成 ####
# 输入：record_ids，该层面的primaryid；reac_df，含outcome_std的JADER反应表；level_name，层面名称 ####
# 输出：一行的汇总数据框 ####
summarise_jader_level <- function(record_ids, reac_df, level_name) {
    record_ids <- unique(as.character(record_ids))
    level_reac <- reac_df %>% dplyr::filter(primaryid %in% record_ids)

    report_status <- level_reac %>%
        dplyr::group_by(primaryid) %>%
        dplyr::summarise(
            any_death = any(outcome_std == "死亡"),
            any_known = any(outcome_std != "Unknown"),
            .groups = "drop"
        )

    n_reports <- length(record_ids)
    n_known <- sum(report_status$any_known)
    n_fatal <- sum(report_status$any_death)
    ci_df <- calc_wilson_ci(n_fatal, n_known)

    tibble::tibble(
        level = level_name,
        dataset = "JADER",
        definition = "Report-level death (any reaction with outcome of death)",
        n_reports = n_reports,
        n_known = n_known,
        n_unknown = n_reports - n_known,
        pct_unknown = (n_reports - n_known) / n_reports * 100,
        n_fatal = n_fatal,
        pct_fatal_known = n_fatal / n_known * 100,
        ci_low_pct = ci_df$ci_low_pct,
        ci_high_pct = ci_df$ci_high_pct,
        pct_fatal_all_reports = n_fatal / n_reports * 100,
        n_reactions = nrow(level_reac),
        pct_reactions_death = mean(level_reac$outcome_std == "死亡") * 100,
        pct_reactions_unknown = mean(level_reac$outcome_std == "Unknown") * 100
    )
}

## 函数：tabulate_faers_outcome_codes ####
# 功能：统计FAERS某一层面报告中各结局编码的报告数及比例（一份报告可有多个编码），并统计无任何编码的报告 ####
# 输入：record_ids，该层面的primaryid；outc_df，FAERS OUTC表；level_name，层面名称 ####
# 输出：每个结局编码一行的数据框 ####
tabulate_faers_outcome_codes <- function(record_ids, outc_df, level_name) {
    record_ids <- unique(as.character(record_ids))
    level_outc <- outc_df %>%
        dplyr::filter(primaryid %in% record_ids, outc_cod %in% names(faers_outcome_labels)) %>%
        dplyr::distinct(primaryid, outc_cod)

    code_counts <- tibble::tibble(outc_cod = names(faers_outcome_labels)) %>%
        dplyr::left_join(dplyr::count(level_outc, outc_cod, name = "n_reports_with_code"), by = "outc_cod") %>%
        dplyr::mutate(n_reports_with_code = tidyr::replace_na(n_reports_with_code, 0L))

    dplyr::bind_rows(
        code_counts,
        tibble::tibble(outc_cod = "None", n_reports_with_code = length(setdiff(record_ids, level_outc$primaryid)))
    ) %>%
        dplyr::mutate(
            level = level_name,
            outcome = dplyr::coalesce(unname(faers_outcome_labels[outc_cod]), "No outcome code (non-serious or not reported)"),
            n_reports = length(record_ids),
            pct_reports = n_reports_with_code / n_reports * 100,
            .before = 1
        ) %>%
        dplyr::select(level, outc_cod, outcome, n_reports_with_code, n_reports, pct_reports)
}

## 函数：tabulate_jader_worst_outcome ####
# 功能：统计JADER某一层面报告的最严重转归（死亡 > 後遺症あり > 未回復 > 軽快 > 回復 > 不明）分布 ####
# 输入：record_ids，该层面的primaryid；reac_df，含outcome_std的JADER反应表；level_name，层面名称 ####
# 输出：每种转归一行的数据框 ####
tabulate_jader_worst_outcome <- function(record_ids, reac_df, level_name) {
    record_ids <- unique(as.character(record_ids))
    outcome_order <- c(jader_known_outcomes, "Unknown")

    worst_df <- reac_df %>%
        dplyr::filter(primaryid %in% record_ids) %>%
        dplyr::mutate(outcome_rank = match(outcome_std, outcome_order)) %>%
        dplyr::group_by(primaryid) %>%
        dplyr::summarise(worst_outcome = outcome_order[min(outcome_rank)], .groups = "drop")

    tibble::tibble(worst_outcome = outcome_order) %>%
        dplyr::left_join(dplyr::count(worst_df, worst_outcome, name = "n_reports_with_outcome"), by = "worst_outcome") %>%
        dplyr::mutate(
            n_reports_with_outcome = tidyr::replace_na(n_reports_with_outcome, 0L),
            level = level_name,
            outcome = unname(jader_outcome_labels[worst_outcome]),
            n_reports = length(record_ids),
            pct_reports = n_reports_with_outcome / n_reports * 100
        ) %>%
        dplyr::select(level, worst_outcome, outcome, n_reports_with_outcome, n_reports, pct_reports)
}


# 4. 数据读取 ####
pattern_assignment_df <- load_pattern_group_assignment(pattern_assignment_file)

irae_pt_keys <- readxl::read_xlsx(irae_dict_file) %>%
    dplyr::transmute(pt_key = standardize_medra_term(PT)) %>%
    dplyr::filter(!is.na(pt_key), pt_key != "") %>%
    dplyr::distinct() %>%
    dplyr::pull(pt_key)

faers_demo_ids <- unique(as.character(load_named_rdata_object(faers_demo_file, "faers_demo")$primaryid))
gc()

faers_outc <- load_named_rdata_object(faers_outc_file, "faers_outc") %>%
    dplyr::transmute(primaryid = as.character(primaryid), outc_cod = trimws(as.character(outc_cod))) %>%
    dplyr::filter(!is.na(outc_cod), outc_cod != "")

faers_solid <- load_named_rdata_object(faers_solid_file, "faers_solid") %>%
    dplyr::mutate(primaryid = as.character(primaryid), primary_organ = as.character(primary_organ))

jader_demo <- load_named_rdata_object(jader_demo_file, "jader_demo") %>%
    dplyr::mutate(primaryid = as.character(primaryid)) %>%
    dplyr::filter(report_year >= 2015, report_year <= 2025)

jader_reac <- load_named_rdata_object(jader_reac_file, "jader_reac") %>%
    dplyr::mutate(
        primaryid = as.character(primaryid),
        outcome_std = standardize_jader_outcome_known(outcome),
        pt_key = dplyr::coalesce(standardize_medra_term(pt_std), standardize_medra_term(pt_term_en)),
        is_irae_reaction = pt_key %in% irae_pt_keys
    )

jader_solid <- load_named_rdata_object(jader_solid_file, "jader_solid") %>%
    dplyr::mutate(primaryid = as.character(primaryid), primary_organ = as.character(primary_organ))


# 5. 队列构建（与6.1.2一致：ICI实体瘤报告 → 映射到五组 → 有irAE） ####
faers_ici <- faers_solid %>% dplyr::filter(has_ici == TRUE)
faers_irae_cohort <- faers_ici %>%
    dplyr::inner_join(pattern_assignment_df, by = "primary_organ") %>%
    dplyr::filter(has_irae == TRUE)

jader_ici <- jader_solid %>% dplyr::filter(has_ici == TRUE)
jader_irae_cohort <- jader_ici %>%
    dplyr::inner_join(pattern_assignment_df, by = "primary_organ") %>%
    dplyr::filter(has_irae == TRUE)

message("FAERS: 全库 ", length(faers_demo_ids), "；ICI实体瘤 ", nrow(faers_ici), "；irAE队列 ", nrow(faers_irae_cohort))
message("JADER: 全库 ", nrow(jader_demo), "；ICI实体瘤 ", nrow(jader_ici), "；irAE队列 ", nrow(jader_irae_cohort))


# 6. 统计分析 ####

## 6.1 各癌种群报告层面致死状态 ####

### FAERS：OUTC中有DE为致死；有任一结局编码但无DE为非致死；无任何编码为结局未知（与6.1.2一致） ####
faers_outc_cohort <- faers_outc %>% dplyr::filter(primaryid %in% faers_irae_cohort$primaryid)

faers_status_known <- faers_irae_cohort %>%
    dplyr::select(primaryid, pattern_group) %>%
    dplyr::mutate(
        fatality_status = dplyr::case_when(
            primaryid %in% faers_outc_cohort$primaryid[faers_outc_cohort$outc_cod == "DE"] ~ "Fatal",
            primaryid %in% faers_outc_cohort$primaryid ~ "Non-fatal",
            TRUE ~ "Unknown"
        )
    )

# 主分析（返修）：分母为全部irAE报告，无结局编码的报告计为非致死，与JADER原稿做法（“不明”计为非致死）一致 ####
faers_status_all_reports <- faers_status_known %>%
    dplyr::mutate(fatality_status = ifelse(fatality_status == "Unknown", "Non-fatal", fatality_status))

### JADER：三种定义 ####
jader_reac_cohort <- jader_reac %>% dplyr::filter(primaryid %in% jader_irae_cohort$primaryid)

jader_irae_reac_status <- jader_reac_cohort %>%
    dplyr::filter(is_irae_reaction) %>%
    dplyr::group_by(primaryid) %>%
    dplyr::summarise(
        irae_death = any(outcome_std == "死亡"),
        irae_any_known = any(outcome_std != "Unknown"),
        irae_any_recorded = any(!is.na(outcome) & trimws(outcome) != "" & trimws(outcome) != "未知"),
        .groups = "drop"
    )

jader_any_reac_status <- jader_reac_cohort %>%
    dplyr::group_by(primaryid) %>%
    dplyr::summarise(
        any_death = any(outcome_std == "死亡"),
        any_known = any(outcome_std != "Unknown"),
        .groups = "drop"
    )

jader_base_df <- jader_irae_cohort %>%
    dplyr::select(primaryid, pattern_group) %>%
    dplyr::left_join(jader_irae_reac_status, by = "primaryid") %>%
    dplyr::left_join(jader_any_reac_status, by = "primaryid") %>%
    dplyr::mutate(dplyr::across(c(irae_death, irae_any_known, irae_any_recorded, any_death, any_known), ~ tidyr::replace_na(.x, FALSE)))

# 定义1（原稿6.1.2）：irAE反应转归为死亡即致死；“不明”计为非致死 ####
jader_status_original <- jader_base_df %>%
    dplyr::mutate(fatality_status = dplyr::case_when(
        irae_death ~ "Fatal",
        irae_any_recorded ~ "Non-fatal",
        TRUE ~ "Unknown"
    ))

# 定义2（修正）：irAE反应转归为死亡即致死；所有irAE反应转归均为“不明”时记为结局未知 ####
jader_status_corrected <- jader_base_df %>%
    dplyr::mutate(fatality_status = dplyr::case_when(
        irae_death ~ "Fatal",
        irae_any_known ~ "Non-fatal",
        TRUE ~ "Unknown"
    ))

# 定义3（与FAERS对齐）：报告中任一反应转归为死亡即致死；所有反应转归均为“不明”时记为结局未知 ####
jader_status_report_level <- jader_base_df %>%
    dplyr::mutate(fatality_status = dplyr::case_when(
        any_death ~ "Fatal",
        any_known ~ "Non-fatal",
        TRUE ~ "Unknown"
    ))

definition_faers_primary <- "FAERS primary: report-level death (OUTC = DE, any cause); denominator all irAE reports, reports without outcome codes counted as non-fatal"
definition_faers_known <- "FAERS sensitivity: report-level death (OUTC = DE, any cause); denominator reports with any outcome code"
definition_jader_primary <- "JADER primary: death as outcome of an irAE reaction; denominator all irAE reports, unknown outcomes counted as non-fatal"
definition_jader_known <- "JADER sensitivity: death as outcome of an irAE reaction; denominator reports with a known irAE outcome"
definition_jader_report_level <- "JADER sensitivity: death as outcome of any reaction in the report (aligned with FAERS); denominator reports with a known outcome"

fatality_definitions <- list(
    list(df = faers_status_all_reports, dataset = "FAERS", definition = definition_faers_primary),
    list(df = faers_status_known, dataset = "FAERS", definition = definition_faers_known),
    list(df = jader_status_original, dataset = "JADER", definition = definition_jader_primary),
    list(df = jader_status_corrected, dataset = "JADER", definition = definition_jader_known),
    list(df = jader_status_report_level, dataset = "JADER", definition = definition_jader_report_level)
)

group_fatality_df <- dplyr::bind_rows(lapply(fatality_definitions, function(x) summarise_group_fatality(x$df, x$dataset, x$definition)))
group_fatality_test_df <- dplyr::bind_rows(lapply(fatality_definitions, function(x) test_group_fatality(x$df, x$dataset, x$definition)))
group_fatality_pairwise_df <- dplyr::bind_rows(lapply(fatality_definitions, function(x) test_pairwise_fatality(x$df, x$dataset, x$definition)))

## 6.2 两库各组致死比例排序的一致性（FAERS定义与JADER定义配对比较） ####
definition_pairs <- tibble::tibble(
    comparison = c("Primary analysis (all irAE reports)", "Known outcomes only", "Known outcomes only, JADER aligned with FAERS"),
    faers_definition = c(definition_faers_primary, definition_faers_known, definition_faers_known),
    jader_definition = c(definition_jader_primary, definition_jader_known, definition_jader_report_level)
)

rank_concordance_df <- dplyr::bind_rows(lapply(seq_len(nrow(definition_pairs)), function(i) {
    faers_df <- group_fatality_df %>% dplyr::filter(definition == definition_pairs$faers_definition[i])
    jader_df <- group_fatality_df %>% dplyr::filter(definition == definition_pairs$jader_definition[i])
    paired_df <- dplyr::inner_join(
        faers_df %>% dplyr::select(pattern_group, faers_pct = pct_fatal_known, faers_rank = rank_fatal_known),
        jader_df %>% dplyr::select(pattern_group, jader_pct = pct_fatal_known, jader_rank = rank_fatal_known),
        by = "pattern_group"
    )

    tibble::tibble(
        comparison = definition_pairs$comparison[i],
        faers_order = paste(paired_df$pattern_group[order(paired_df$faers_rank)], collapse = " > "),
        jader_order = paste(paired_df$pattern_group[order(paired_df$jader_rank)], collapse = " > "),
        identical_order = identical(paired_df$pattern_group[order(paired_df$faers_rank)], paired_df$pattern_group[order(paired_df$jader_rank)]),
        spearman_rho = stats::cor(paired_df$faers_pct, paired_df$jader_pct, method = "spearman"),
        ratio_faers_to_jader_min = min(paired_df$faers_pct / paired_df$jader_pct),
        ratio_faers_to_jader_max = max(paired_df$faers_pct / paired_df$jader_pct)
    )
}))

## 6.3 数据库层面死亡比例：全库、ICI实体瘤报告、irAE队列 ####
database_level_df <- dplyr::bind_rows(
    summarise_faers_level(faers_demo_ids, faers_outc, level_labels[["whole"]]),
    summarise_faers_level(faers_ici$primaryid, faers_outc, level_labels[["ici_solid"]]),
    summarise_faers_level(faers_irae_cohort$primaryid, faers_outc, level_labels[["irae_cohort"]]),
    summarise_jader_level(jader_demo$primaryid, jader_reac, level_labels[["whole"]]),
    summarise_jader_level(jader_ici$primaryid, jader_reac, level_labels[["ici_solid"]]),
    summarise_jader_level(jader_irae_cohort$primaryid, jader_reac, level_labels[["irae_cohort"]])
) %>%
    dplyr::mutate(level = factor(level, levels = unname(level_labels))) %>%
    dplyr::arrange(level, dataset) %>%
    dplyr::mutate(level = as.character(level))

database_ratio_df <- database_level_df %>%
    dplyr::select(level, dataset, pct_fatal_known, pct_fatal_all_reports) %>%
    tidyr::pivot_wider(names_from = dataset, values_from = c(pct_fatal_known, pct_fatal_all_reports)) %>%
    dplyr::mutate(
        ratio_faers_to_jader_known = pct_fatal_known_FAERS / pct_fatal_known_JADER,
        ratio_faers_to_jader_all_reports = pct_fatal_all_reports_FAERS / pct_fatal_all_reports_JADER
    )

## 6.4 结局构成：FAERS结局编码与JADER最严重转归 ####
faers_outcome_codes_df <- dplyr::bind_rows(
    tabulate_faers_outcome_codes(faers_demo_ids, faers_outc, level_labels[["whole"]]),
    tabulate_faers_outcome_codes(faers_ici$primaryid, faers_outc, level_labels[["ici_solid"]]),
    tabulate_faers_outcome_codes(faers_irae_cohort$primaryid, faers_outc, level_labels[["irae_cohort"]])
)

jader_worst_outcome_df <- dplyr::bind_rows(
    tabulate_jader_worst_outcome(jader_demo$primaryid, jader_reac, level_labels[["whole"]]),
    tabulate_jader_worst_outcome(jader_ici$primaryid, jader_reac, level_labels[["ici_solid"]]),
    tabulate_jader_worst_outcome(jader_irae_cohort$primaryid, jader_reac, level_labels[["irae_cohort"]])
)

## 6.5 FAERS发病时间：各组有效日期例数与排除例数（截尾规则与6.1.2一致） ####
faers_tto_df <- faers_irae_cohort %>%
    dplyr::transmute(
        primaryid,
        pattern_group,
        ici_tto = suppressWarnings(as.numeric(as.character(ici_tto)))
    ) %>%
    dplyr::mutate(valid_tto = !is.na(ici_tto) & ici_tto >= 0)

tto_cutoff_99 <- as.numeric(stats::quantile(faers_tto_df$ici_tto[faers_tto_df$valid_tto], probs = 0.99, names = FALSE))

faers_tto_df <- faers_tto_df %>%
    dplyr::mutate(analysed = valid_tto & ici_tto <= tto_cutoff_99)

tto_by_group_df <- faers_tto_df %>%
    dplyr::mutate(pattern_group = factor(pattern_group, levels = pattern_display_order)) %>%
    dplyr::group_by(pattern_group) %>%
    dplyr::summarise(
        n_irae_reports = dplyr::n(),
        n_missing_or_invalid_dates = sum(!valid_tto),
        n_valid_dates = sum(valid_tto),
        pct_valid_dates = mean(valid_tto) * 100,
        n_above_99th_percentile = sum(valid_tto & !analysed),
        n_analysed = sum(analysed),
        median_days = stats::median(ici_tto[analysed]),
        q1_days = stats::quantile(ici_tto[analysed], 0.25, names = FALSE),
        q3_days = stats::quantile(ici_tto[analysed], 0.75, names = FALSE),
        .groups = "drop"
    ) %>%
    dplyr::mutate(pattern_group = as.character(pattern_group))

tto_by_group_df <- dplyr::bind_rows(
    tto_by_group_df,
    faers_tto_df %>%
        dplyr::summarise(
            pattern_group = "All five groups",
            n_irae_reports = dplyr::n(),
            n_missing_or_invalid_dates = sum(!valid_tto),
            n_valid_dates = sum(valid_tto),
            pct_valid_dates = mean(valid_tto) * 100,
            n_above_99th_percentile = sum(valid_tto & !analysed),
            n_analysed = sum(analysed),
            median_days = stats::median(ici_tto[analysed]),
            q1_days = stats::quantile(ici_tto[analysed], 0.25, names = FALSE),
            q3_days = stats::quantile(ici_tto[analysed], 0.75, names = FALSE)
        )
) %>%
    dplyr::mutate(tto_cutoff_99th_percentile_days = tto_cutoff_99)

# 发病时间日期缺失是否在各组间不同 ####
tto_missing_test <- stats::chisq.test(table(faers_tto_df$pattern_group, faers_tto_df$valid_tto))
tto_missing_test_df <- tibble::tibble(
    comparison = "Proportion of irAE reports with valid onset dates across groups (FAERS)",
    chi_squared = unname(tto_missing_test$statistic),
    df = unname(tto_missing_test$parameter),
    p_value = tto_missing_test$p.value
)

message("各组致死比例:")
print(as.data.frame(group_fatality_df), digits = 3)
message("两库排序一致性:")
print(as.data.frame(rank_concordance_df), digits = 3)
message("数据库层面死亡比例:")
print(as.data.frame(database_level_df), digits = 3)
message("FAERS发病时间有效例数:")
print(as.data.frame(tto_by_group_df), digits = 3)


# 7. 结果保存 ####
writexl::write_xlsx(
    list(
        "各组致死比例" = group_fatality_df,
        "各组致死比例总体检验" = group_fatality_test_df,
        "各组致死比例两两比较" = group_fatality_pairwise_df,
        "两库各组排序一致性" = rank_concordance_df,
        "数据库层面死亡比例" = database_level_df,
        "数据库层面FAERS与JADER之比" = database_ratio_df,
        "FAERS结局编码分布" = faers_outcome_codes_df,
        "JADER最严重转归分布" = jader_worst_outcome_df,
        "FAERS发病时间有效例数" = tto_by_group_df,
        "发病时间日期缺失组间检验" = tto_missing_test_df,
        "分析参数" = tibble::tibble(
            参数 = c(
                "FAERS致死定义", "FAERS结局未知", "JADER转归未知", "主分析（返修）", "FAERS敏感性分析", "JADER敏感性分析1", "JADER敏感性分析2（对齐FAERS）",
                "置信区间与检验", "全库范围", "ICI实体瘤报告", "irAE队列", "发病时间"
            ),
            值 = c(
                "OUTC表中outc_cod = DE（报告层面，不限死因）",
                "该报告在OUTC表中无任何结局编码",
                "转归为“不明”“未知”或空值",
                "两库分母均为全部irAE报告，结局缺失计为非致死；JADER为irAE反应转归为死亡（与原稿6.1.2一致），FAERS为报告层面DE",
                "FAERS分母为有任一结局编码的报告（原稿6.1.2做法）",
                "JADER irAE反应转归为死亡；全部irAE反应转归未知时记为结局未知，不计入分母",
                "JADER报告中任一反应转归为死亡即致死；全部反应转归未知时记为结局未知，不计入分母",
                "Wilson 95%置信区间；总体卡方检验；两两比较为pairwise.prop.test，Bonferroni校正（与6.1.2一致）；表中“n_known”为各定义下的分母",
                "FAERS：去重后2015–2025年全部报告（faers_demo）；JADER：2015–2025年全部报告（jader_demo）",
                "faers_solid / jader_solid中has_ici为TRUE的报告",
                "ICI实体瘤报告中映射到五组且has_irae为TRUE的报告（与6.1.2一致）",
                "FAERS ici_tto为数值且≥0为有效日期；超过有效值第99百分位者排除（与6.1.2一致）；JADER无可用发病时间"
            )
        )
    ),
    path = output_xlsx_file
)

message("已生成结果文件: ", output_xlsx_file)
message("========== 致死结局与发病时间不确定性分析结束 ==========")
