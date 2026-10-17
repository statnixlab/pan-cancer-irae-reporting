# 癌种群 irAE 表型分析（FAERS + JADER）####

rm(list = ls())
gc()

# 1. 加载配置与依赖 ####
source("Code/config.R")
source("Code/mapping_rules.R")

if (!requireNamespace("dunn.test", quietly = TRUE)) {
    install.packages("dunn.test", quiet = TRUE)
}

if (!requireNamespace("MASS", quietly = TRUE)) {
    install.packages("MASS", quiet = TRUE)
}


# 2. 结果目录与参数 ####
script_name <- "6.1.2癌种群irAE表型分析"
script_results_path <- file.path(results_path, script_name)

if (!dir.exists(script_results_path)) {
    dir.create(script_results_path, recursive = TRUE, showWarnings = FALSE)
}

pattern_assignment_file <- file.path(
    results_path,
    "6.1癌种NMF分群分析_rank5",
    "FAERS_NMF_WH_matrix.xlsx"
)

faers_solid_file <- file.path(output_path, "faers_solid.RData")
faers_outc_file <- file.path(output_path, "faers_outc.RData")
jader_solid_file <- file.path(output_path, "jader_solid.RData")
jader_reac_file <- file.path(output_path, "jader_reac.RData")
irae_dict_file <- file.path(data_path, "irAE列表.xlsx")

faers_output_file <- file.path(script_results_path, "FAERS_癌种群irAE表型汇总.xlsx")
jader_output_file <- file.path(script_results_path, "JADER_癌种群irAE表型汇总.xlsx")

pattern_display_order <- c(
    "Ocular",
    "Hematologic",
    "Neurologic",
    "Systemic-like",
    "Visceral-like"
)

message("========== 癌种群 irAE 表型分析开始 ==========")
message("结果目录: ", script_results_path)
message("分群表: ", pattern_assignment_file)


# 3. 辅助函数 ####

## 函数：load_named_rdata_object ####
## 功能：从RData文件中读取指定对象 ####
## 输入：input_path，RData文件路径；object_name，对象名称 ####
## 输出：返回读取出的对象 ####
load_named_rdata_object <- function(input_path, object_name) {
    input_env <- new.env(parent = emptyenv())

    if (!file.exists(input_path)) {
        stop("输入文件不存在: ", input_path)
    }

    load(input_path, envir = input_env)

    if (!exists(object_name, envir = input_env, inherits = FALSE)) {
        stop("RData中未找到对象: ", object_name)
    }

    input_env[[object_name]]
}

## 函数：standardize_medra_term ####
## 功能：标准化MedDRA术语文本，便于不同数据源之间匹配 ####
## 输入：x，字符型术语向量 ####
## 输出：返回标准化后的大写字符向量 ####
standardize_medra_term <- function(x) {
    x <- as.character(x)
    x <- stringr::str_replace_all(x, "_", " ")
    x <- stringr::str_replace_all(x, "-", " ")
    x <- stringr::str_replace_all(x, "\\s+", " ")
    x <- trimws(x)
    toupper(x)
}

## 函数：standardize_jader_outcome ####
## 功能：标准化JADER反应级结局文本，将“不明/未知/空值”统一为缺失 ####
## 输入：x，字符型结局向量 ####
## 输出：返回标准化后的字符向量 ####
standardize_jader_outcome <- function(x) {
    x <- as.character(x)
    x <- trimws(x)
    x[x %in% c("未知")] <- NA_character_
    x
}

## 函数：get_soc_abbrev_safe ####
## 功能：将SOC全称转换为缩写，无法映射时返回原始SOC名称 ####
## 输入：soc_full，字符型SOC全称 ####
## 输出：返回SOC缩写或原始全称 ####
get_soc_abbrev_safe <- function(soc_full) {
    current_abbrev <- get_soc_abbrev(soc_full)

    if (is.na(current_abbrev) || current_abbrev == "") {
        return(as.character(soc_full))
    }

    as.character(current_abbrev)
}

## 函数：build_significance_label ####
## 功能：根据P值生成显著性标签 ####
## 输入：p_value，数值型P值 ####
## 输出：返回显著性标签字符值 ####
build_significance_label <- function(p_value) {
    dplyr::case_when(
        is.na(p_value) ~ NA_character_,
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ "ns"
    )
}

## 函数：load_pattern_group_assignment ####
## 功能：读取6.1.1导出的癌种群分配表并校验分群顺序 ####
## 输入：input_file，分群Excel路径；pattern_order，固定的pattern顺序 ####
## 输出：返回标准化后的primary_organ与pattern_group对照表 ####
load_pattern_group_assignment <- function(input_file, pattern_order) {
    if (!file.exists(input_file)) {
        stop("分群结果文件不存在: ", input_file)
    }

    assignment_df <- readxl::read_xlsx(
        input_file,
        sheet = "Pattern group assignment"
    )

    if (ncol(assignment_df) < 2) {
        stop("Pattern group assignment工作表字段不足")
    }

    names(assignment_df)[1:2] <- c("primary_organ", "pattern_group")

    assignment_df <- assignment_df %>%
        dplyr::transmute(
            primary_organ = as.character(primary_organ),
            pattern_group = as.character(pattern_group)
        ) %>%
        dplyr::filter(!is.na(primary_organ), primary_organ != "") %>%
        dplyr::distinct(primary_organ, .keep_all = TRUE)

    if (anyDuplicated(assignment_df$primary_organ) > 0) {
        stop("分群表中存在重复的primary_organ")
    }

    if (!setequal(unique(assignment_df$pattern_group), pattern_order)) {
        stop("分群表中的pattern_group与预设pattern_display_order不一致")
    }

    assignment_df %>%
        dplyr::mutate(
            pattern_group = factor(pattern_group, levels = pattern_order)
        ) %>%
        dplyr::arrange(pattern_group, primary_organ) %>%
        dplyr::mutate(pattern_group = as.character(pattern_group))
}

## 函数：load_soc_order_from_w_matrix ####
## 功能：读取6.1.1的W矩阵行顺序，作为SOC谱固定展示顺序 ####
## 输入：input_file，分群Excel路径 ####
## 输出：返回SOC缩写顺序字符向量 ####
load_soc_order_from_w_matrix <- function(input_file) {
    w_matrix_df <- readxl::read_xlsx(input_file, sheet = "W_matrix")

    if (!"soc_abbrev" %in% names(w_matrix_df)) {
        stop("W_matrix工作表缺少soc_abbrev列")
    }

    soc_order <- as.character(w_matrix_df$soc_abbrev)
    soc_order <- soc_order[!is.na(soc_order) & soc_order != ""]
    unique(soc_order)
}

## 函数：extract_ordered_soc_metadata ####
## 功能：从分析集提取SOC列信息并对齐到NMF输出的SOC顺序 ####
## 输入：data，分析数据框；soc_order_vec，SOC缩写顺序向量 ####
## 输出：返回包含SOC变量名、全称、缩写和顺序的数据框 ####
extract_ordered_soc_metadata <- function(data, soc_order_vec) {
    soc_vars <- grep("^SOC_", names(data), value = TRUE)

    if (length(soc_vars) == 0) {
        stop("数据中未检测到SOC列")
    }

    soc_metadata <- tibble::tibble(
        soc_var = soc_vars,
        soc_full = gsub("_", " ", sub("^SOC_", "", soc_vars))
    ) %>%
        dplyr::mutate(
            soc_abbrev = vapply(soc_full, get_soc_abbrev_safe, character(1)),
            soc_order = match(soc_abbrev, soc_order_vec)
        )

    missing_order_index <- which(is.na(soc_metadata$soc_order))

    if (length(missing_order_index) > 0) {
        max_order <- length(soc_order_vec)
        soc_metadata$soc_order[missing_order_index] <- max_order + seq_along(missing_order_index)
    }

    soc_metadata %>%
        dplyr::arrange(soc_order, soc_abbrev, soc_full)
}

## 函数：prepare_mapped_ici_dataset ####
## 功能：筛选ICI病例、映射pattern group并生成队列流程表 ####
## 输入：data，原始分析集；id_col，主键列名；dataset_name，数据库名称；assignment_df，分群对照表；pattern_order，固定pattern顺序 ####
## 输出：返回包含mapped队列、irAE子集和队列流程表的列表 ####
prepare_mapped_ici_dataset <- function(data, id_col, dataset_name, assignment_df, pattern_order) {
    data_step0 <- data %>%
        dplyr::mutate(
            primary_organ = as.character(primary_organ),
            record_id = as.character(.data[[id_col]])
        )

    ici_data <- data_step0 %>%
        dplyr::filter(has_ici == TRUE)

    mapped_data <- ici_data %>%
        dplyr::left_join(assignment_df, by = "primary_organ") %>%
        dplyr::filter(!is.na(pattern_group)) %>%
        dplyr::mutate(
            pattern_group = factor(pattern_group, levels = pattern_order)
        ) %>%
        dplyr::arrange(pattern_group, primary_organ, record_id) %>%
        dplyr::mutate(pattern_group = as.character(pattern_group))

    mapped_irae_data <- mapped_data %>%
        dplyr::filter(has_irae == TRUE)

    unmatched_organs <- ici_data %>%
        dplyr::anti_join(assignment_df, by = "primary_organ") %>%
        dplyr::mutate(primary_organ = dplyr::coalesce(primary_organ, "Missing")) %>%
        dplyr::count(primary_organ, name = "n") %>%
        dplyr::arrange(dplyr::desc(n), primary_organ)

    total_ici_n <- nrow(ici_data)
    mapped_ici_n <- nrow(mapped_data)
    mapped_irae_ici_n <- nrow(mapped_irae_data)
    unmatched_ici_n <- sum(unmatched_organs$n)

    flow_summary <- tibble::tribble(
        ~section, ~item, ~primary_organ, ~n, ~pct_total_ici, ~notes,
        "summary", "total_ici_n", NA_character_, total_ici_n, 100, "All ICI cases before pattern-group mapping",
        "summary", "mapped_ici_n", NA_character_, mapped_ici_n, ifelse(total_ici_n > 0, mapped_ici_n / total_ici_n * 100, NA_real_), "ICI cases with mappable primary_organ",
        "summary", "mapped_irae_ici_n", NA_character_, mapped_irae_ici_n, ifelse(total_ici_n > 0, mapped_irae_ici_n / total_ici_n * 100, NA_real_), "Mapped ICI cases with has_irae == TRUE",
        "summary", "unmatched_ici_n", NA_character_, unmatched_ici_n, ifelse(total_ici_n > 0, unmatched_ici_n / total_ici_n * 100, NA_real_), "ICI cases excluded because primary_organ is not in 6.1.1 assignment"
    )

    unmatched_flow <- unmatched_organs %>%
        dplyr::mutate(
            section = "unmatched_primary_organ",
            item = "excluded_primary_organ",
            pct_total_ici = ifelse(total_ici_n > 0, n / total_ici_n * 100, NA_real_),
            notes = "Excluded from downstream statistics"
        ) %>%
        dplyr::select(section, item, primary_organ, n, pct_total_ici, notes)

    detected_pattern_groups <- sort(unique(mapped_data$pattern_group))

    if (!all(detected_pattern_groups %in% pattern_order)) {
        stop(dataset_name, "存在非法pattern_group名称")
    }

    list(
        ici_data = ici_data,
        mapped_data = mapped_data,
        mapped_irae_data = mapped_irae_data,
        flow_table = dplyr::bind_rows(flow_summary, unmatched_flow)
    )
}

## 函数：summarise_pattern_group_n ####
## 功能：统计每个pattern group的mapped ICI样本量与mapped irAE-ICI样本量 ####
## 输入：mapped_data，已映射ICI病例；mapped_irae_data，已映射且有irAE的ICI病例；pattern_order，固定pattern顺序 ####
## 输出：返回pattern group样本量汇总表 ####
summarise_pattern_group_n <- function(mapped_data, mapped_irae_data, pattern_order) {
    mapped_ici_n_df <- mapped_data %>%
        dplyr::count(pattern_group, name = "mapped_ici_n")

    mapped_irae_n_df <- mapped_irae_data %>%
        dplyr::count(pattern_group, name = "mapped_irae_ici_n")

    tibble::tibble(pattern_group = pattern_order) %>%
        dplyr::left_join(mapped_ici_n_df, by = "pattern_group") %>%
        dplyr::left_join(mapped_irae_n_df, by = "pattern_group") %>%
        dplyr::mutate(
            mapped_ici_n = tidyr::replace_na(mapped_ici_n, 0L),
            mapped_irae_ici_n = tidyr::replace_na(mapped_irae_ici_n, 0L),
            any_irae_rp = dplyr::if_else(
                mapped_ici_n > 0,
                mapped_irae_ici_n / mapped_ici_n * 100,
                NA_real_
            )
        )
}

## 函数：build_soc_spectrum ####
## 功能：构建pattern group × SOC的reporting spectrum结果表 ####
## 输入：mapped_data，已映射ICI病例；soc_metadata，SOC元数据表；pattern_group_n_df，pattern group样本量表；pattern_order，固定pattern顺序 ####
## 输出：返回SOC谱汇总表 ####
build_soc_spectrum <- function(mapped_data, soc_metadata, pattern_group_n_df, pattern_order) {
    if (nrow(mapped_data) == 0) {
        return(tibble::tibble())
    }

    soc_count_df <- mapped_data %>%
        dplyr::group_by(pattern_group) %>%
        dplyr::summarise(
            dplyr::across(
                dplyr::all_of(soc_metadata$soc_var),
                ~ sum(.x == 1, na.rm = TRUE)
            ),
            .groups = "drop"
        ) %>%
        tidyr::pivot_longer(
            cols = dplyr::all_of(soc_metadata$soc_var),
            names_to = "soc_var",
            values_to = "n"
        )

    tidyr::expand_grid(
        pattern_group = pattern_order,
        soc_var = soc_metadata$soc_var
    ) %>%
        dplyr::left_join(soc_count_df, by = c("pattern_group", "soc_var")) %>%
        dplyr::left_join(soc_metadata, by = "soc_var") %>%
        dplyr::left_join(pattern_group_n_df, by = "pattern_group") %>%
        dplyr::mutate(
            n = tidyr::replace_na(n, 0L),
            RP = dplyr::if_else(mapped_ici_n > 0, n / mapped_ici_n * 100, 0)
        ) %>%
        dplyr::group_by(pattern_group) %>%
        dplyr::arrange(dplyr::desc(RP), dplyr::desc(n), soc_order, .by_group = TRUE) %>%
        dplyr::mutate(rank_within_pattern = dplyr::row_number()) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(factor(pattern_group, levels = pattern_order), soc_order, soc_abbrev) %>%
        dplyr::select(
            pattern_group, mapped_ici_n,
            soc_var, soc_full, soc_abbrev, soc_order,
            n, RP, rank_within_pattern
        )
}

## 函数：prepare_tto_case_level ####
## 功能：整理病例级TTO数据并执行99%分位数截尾 ####
## 输入：mapped_irae_data，已映射且有irAE的ICI病例；dataset_name，数据库名称 ####
## 输出：返回包含病例级TTO表、信息表和截尾值的列表 ####
prepare_tto_case_level <- function(mapped_irae_data, dataset_name) {
    tto_raw_df <- mapped_irae_data %>%
        dplyr::transmute(
            dataset = dataset_name,
            record_id = as.character(record_id),
            primary_organ = as.character(primary_organ),
            pattern_group = as.character(pattern_group),
            ici_tto_raw = as.character(ici_tto),
            ici_tto = suppressWarnings(as.numeric(as.character(ici_tto)))
        )

    valid_tto_df <- tto_raw_df %>%
        dplyr::filter(!is.na(ici_tto), ici_tto >= 0)

    tto_cutoff_99 <- if (nrow(valid_tto_df) > 0) {
        as.numeric(stats::quantile(valid_tto_df$ici_tto, probs = 0.99, na.rm = TRUE, names = FALSE))
    } else {
        NA_real_
    }

    tto_case_level_df <- valid_tto_df %>%
        dplyr::filter(is.na(tto_cutoff_99) | ici_tto <= tto_cutoff_99) %>%
        dplyr::mutate(tto_cutoff_99 = tto_cutoff_99) %>%
        dplyr::arrange(factor(pattern_group, levels = pattern_display_order), ici_tto, record_id)

    excluded_invalid_n <- nrow(tto_raw_df) - nrow(valid_tto_df)
    excluded_above_cutoff_n <- nrow(valid_tto_df) - nrow(tto_case_level_df)

    info_df <- tibble::tribble(
        ~section, ~metric, ~value, ~notes,
        "analysis_info", "input_mapped_irae_ici_n", nrow(tto_raw_df), "Mapped irAE-ICI cases before TTO filtering",
        "analysis_info", "valid_non_negative_tto_n", nrow(valid_tto_df), "Cases with numeric ici_tto >= 0",
        "analysis_info", "excluded_invalid_tto_n", excluded_invalid_n, "Excluded because ici_tto is No ICI, Unknown, non-numeric or negative",
        "analysis_info", "tto_99_percentile_cutoff", tto_cutoff_99, "Database-level 99th percentile cutoff for ici_tto",
        "analysis_info", "excluded_above_99pct_n", excluded_above_cutoff_n, "Excluded because ici_tto exceeds 99th percentile cutoff",
        "analysis_info", "final_tto_analysis_n", nrow(tto_case_level_df), "Final case count used for TTO comparison"
    )

    list(
        case_level = tto_case_level_df,
        info_df = info_df,
        tto_cutoff_99 = tto_cutoff_99
    )
}

## 函数：build_tto_summary ####
## 功能：汇总各群TTO描述统计并执行Kruskal-Wallis与Dunn事后检验 ####
## 输入：tto_case_level_df，病例级TTO表；info_df，TTO分析信息表；dataset_name，数据库名称；pattern_order，固定pattern顺序 ####
## 输出：返回TTO汇总结果表 ####
build_tto_summary <- function(tto_case_level_df, info_df, dataset_name, pattern_order) {
    group_summary_df <- tibble::tibble(pattern_group = pattern_order) %>%
        dplyr::left_join(
            tto_case_level_df %>%
                dplyr::group_by(pattern_group) %>%
                dplyr::summarise(
                    n = dplyr::n(),
                    median_tto = stats::median(ici_tto, na.rm = TRUE),
                    q1_tto = stats::quantile(ici_tto, probs = 0.25, na.rm = TRUE, names = FALSE),
                    q3_tto = stats::quantile(ici_tto, probs = 0.75, na.rm = TRUE, names = FALSE),
                    mean_tto = mean(ici_tto, na.rm = TRUE),
                    sd_tto = stats::sd(ici_tto, na.rm = TRUE),
                    .groups = "drop"
                ),
            by = "pattern_group"
        ) %>%
        dplyr::mutate(
            section = "group_summary",
            dataset = dataset_name
        ) %>%
        dplyr::select(
            dataset, section, pattern_group, n,
            median_tto, q1_tto, q3_tto, mean_tto, sd_tto
        )

    nonzero_groups_n <- sum(group_summary_df$n > 0, na.rm = TRUE)

    if (nrow(tto_case_level_df) > 0 && nonzero_groups_n >= 2) {
        kw_test <- stats::kruskal.test(ici_tto ~ pattern_group, data = tto_case_level_df)

        global_test_df <- tibble::tibble(
            dataset = dataset_name,
            section = "global_test",
            metric = "Kruskal-Wallis",
            test_name = "Kruskal-Wallis test",
            statistic = unname(kw_test$statistic),
            df = unname(kw_test$parameter),
            p_value = kw_test$p.value,
            significance = build_significance_label(kw_test$p.value),
            notes = "Overall TTO comparison across pattern groups"
        )

        dunn_result <- suppressWarnings(
            dunn.test::dunn.test(
                x = tto_case_level_df$ici_tto,
                g = tto_case_level_df$pattern_group,
                method = "bonferroni",
                kw = FALSE,
                list = TRUE,
                table = FALSE
            )
        )

        pairwise_df <- tibble::tibble(
            comparison = as.character(dunn_result$comparisons),
            z_statistic = as.numeric(dunn_result$Z),
            p_value = as.numeric(dunn_result$P.adjusted)
        ) %>%
            tidyr::separate(
                col = comparison,
                into = c("group1", "group2"),
                sep = "\\s+-\\s+"
            ) %>%
            dplyr::mutate(
                dataset = dataset_name,
                section = "pairwise_test",
                comparison = paste(group1, group2, sep = " vs "),
                test_name = "Dunn test",
                p_adjust_method = "bonferroni",
                significance = build_significance_label(p_value)
            ) %>%
            dplyr::select(
                dataset, section, group1, group2, comparison,
                test_name, z_statistic, p_value, p_adjust_method, significance
            )
    } else {
        global_test_df <- tibble::tibble(
            dataset = dataset_name,
            section = "global_test",
            metric = "Kruskal-Wallis",
            test_name = "Kruskal-Wallis test",
            statistic = NA_real_,
            df = NA_real_,
            p_value = NA_real_,
            significance = NA_character_,
            notes = "Skipped because fewer than two pattern groups have valid TTO values"
        )

        pairwise_df <- tibble::tibble(
            dataset = character(0),
            section = character(0),
            group1 = character(0),
            group2 = character(0),
            comparison = character(0),
            test_name = character(0),
            z_statistic = numeric(0),
            p_value = numeric(0),
            p_adjust_method = character(0),
            significance = character(0)
        )
    }

    info_output_df <- info_df %>%
        dplyr::mutate(dataset = dataset_name) %>%
        dplyr::select(dataset, dplyr::everything())

    dplyr::bind_rows(
        info_output_df,
        group_summary_df,
        global_test_df,
        pairwise_df
    )
}

## 函数：classify_weibull_hazard_type ####
## 功能：根据Weibull形状参数判断风险变化类型 ####
## 输入：shape_value，Weibull形状参数 ####
## 输出：返回风险类型字符值 ####
classify_weibull_hazard_type <- function(shape_value) {
    dplyr::case_when(
        is.na(shape_value) ~ NA_character_,
        abs(shape_value - 1) <= 0.1 ~ "Approximately constant hazard",
        shape_value < 1 ~ "Decreasing hazard",
        shape_value > 1 ~ "Increasing hazard",
        TRUE ~ NA_character_
    )
}

## 函数：build_tto_weibull_results ####
## 功能：对各pattern group的TTO进行Weibull拟合并整理结果表 ####
## 输入：tto_case_level_df，病例级TTO表；dataset_name，数据库名称；pattern_order，固定pattern顺序 ####
## 输出：返回包含Weibull汇总表、密度曲线表和QQ表的列表 ####
build_tto_weibull_results <- function(tto_case_level_df, dataset_name, pattern_order) {
    summary_list <- list()
    density_curve_list <- list()
    qq_list <- list()

    for (current_pattern_group in pattern_order) {
        current_group_df <- tto_case_level_df %>%
            dplyr::filter(pattern_group == current_pattern_group) %>%
            dplyr::arrange(ici_tto, record_id)

        original_n <- nrow(current_group_df)
        raw_tto_values <- current_group_df$ici_tto
        raw_tto_values <- raw_tto_values[!is.na(raw_tto_values) & raw_tto_values >= 0]
        zero_tto_n <- sum(raw_tto_values == 0, na.rm = TRUE)

        weibull_input_values <- raw_tto_values
        weibull_input_values[weibull_input_values == 0] <- 0.5

        if (length(weibull_input_values) < 2 || stats::sd(weibull_input_values) <= 0) {
            summary_list[[length(summary_list) + 1]] <- tibble::tibble(
                dataset = dataset_name,
                pattern_group = current_pattern_group,
                original_n = original_n,
                weibull_input_n = length(weibull_input_values),
                zero_tto_n = zero_tto_n,
                zero_tto_pct = ifelse(original_n > 0, zero_tto_n / original_n * 100, NA_real_),
                zero_replaced_value = 0.5,
                fit_status = "skipped",
                fit_note = "Insufficient variation for Weibull fitting",
                shape = NA_real_,
                scale = NA_real_,
                hazard_type = NA_character_,
                aic = NA_real_,
                bic = NA_real_,
                ks_p_value = NA_real_,
                theoretical_mean = NA_real_,
                theoretical_median = NA_real_,
                observed_mean = ifelse(length(weibull_input_values) > 0, mean(weibull_input_values), NA_real_),
                observed_median = ifelse(length(weibull_input_values) > 0, stats::median(weibull_input_values), NA_real_)
            )
            next
        }

        current_mean <- mean(weibull_input_values)
        current_var <- stats::var(weibull_input_values)

        shape_init <- ifelse(is.na(current_var) || current_var <= 0, 1, (current_mean^2) / current_var)
        scale_init <- current_mean

        shape_init <- max(0.5, min(shape_init, 10))
        scale_init <- max(0.5, scale_init)

        weibull_fit <- tryCatch(
            MASS::fitdistr(
                weibull_input_values,
                densfun = "weibull",
                start = list(shape = shape_init, scale = scale_init)
            ),
            error = function(e) e
        )

        if (inherits(weibull_fit, "error")) {
            summary_list[[length(summary_list) + 1]] <- tibble::tibble(
                dataset = dataset_name,
                pattern_group = current_pattern_group,
                original_n = original_n,
                weibull_input_n = length(weibull_input_values),
                zero_tto_n = zero_tto_n,
                zero_tto_pct = ifelse(original_n > 0, zero_tto_n / original_n * 100, NA_real_),
                zero_replaced_value = 0.5,
                fit_status = "failed",
                fit_note = conditionMessage(weibull_fit),
                shape = NA_real_,
                scale = NA_real_,
                hazard_type = NA_character_,
                aic = NA_real_,
                bic = NA_real_,
                ks_p_value = NA_real_,
                theoretical_mean = NA_real_,
                theoretical_median = NA_real_,
                observed_mean = mean(weibull_input_values),
                observed_median = stats::median(weibull_input_values)
            )
            next
        }

        shape_value <- unname(weibull_fit$estimate["shape"])
        scale_value <- unname(weibull_fit$estimate["scale"])
        loglik_value <- unname(weibull_fit$loglik)
        parameter_n <- 2

        aic_value <- -2 * loglik_value + 2 * parameter_n
        bic_value <- -2 * loglik_value + log(length(weibull_input_values)) * parameter_n

        ks_p_value <- tryCatch(
            suppressWarnings(
                stats::ks.test(
                    weibull_input_values,
                    "pweibull",
                    shape = shape_value,
                    scale = scale_value
                )$p.value
            ),
            error = function(e) NA_real_
        )

        theoretical_mean <- scale_value * gamma(1 + 1 / shape_value)
        theoretical_median <- scale_value * (log(2))^(1 / shape_value)
        hazard_type <- classify_weibull_hazard_type(shape_value)

        summary_list[[length(summary_list) + 1]] <- tibble::tibble(
            dataset = dataset_name,
            pattern_group = current_pattern_group,
            original_n = original_n,
            weibull_input_n = length(weibull_input_values),
            zero_tto_n = zero_tto_n,
            zero_tto_pct = ifelse(original_n > 0, zero_tto_n / original_n * 100, NA_real_),
            zero_replaced_value = 0.5,
            fit_status = "success",
            fit_note = NA_character_,
            shape = shape_value,
            scale = scale_value,
            hazard_type = hazard_type,
            aic = aic_value,
            bic = bic_value,
            ks_p_value = ks_p_value,
            theoretical_mean = theoretical_mean,
            theoretical_median = theoretical_median,
            observed_mean = mean(weibull_input_values),
            observed_median = stats::median(weibull_input_values)
        )

        curve_x_max <- max(raw_tto_values, na.rm = TRUE)
        if (!is.finite(curve_x_max) || curve_x_max <= 0) {
            curve_x_max <- max(weibull_input_values, na.rm = TRUE)
        }
        curve_x_max <- max(curve_x_max, 1)

        density_curve_df <- tibble::tibble(
            dataset = dataset_name,
            pattern_group = current_pattern_group,
            x = seq(0, curve_x_max * 1.1, length.out = 400)
        ) %>%
            dplyr::mutate(
                density = stats::dweibull(
                    pmax(x, 1e-06),
                    shape = shape_value,
                    scale = scale_value
                ),
                shape = shape_value,
                scale = scale_value,
                hazard_type = hazard_type
            )

        qq_df <- tibble::tibble(
            dataset = dataset_name,
            pattern_group = current_pattern_group,
            quantile_rank = seq_along(weibull_input_values),
            probability = stats::ppoints(length(weibull_input_values)),
            sample_quantile = sort(weibull_input_values),
            theoretical_quantile = stats::qweibull(
                stats::ppoints(length(weibull_input_values)),
                shape = shape_value,
                scale = scale_value
            ),
            shape = shape_value,
            scale = scale_value,
            hazard_type = hazard_type
        )

        density_curve_list[[length(density_curve_list) + 1]] <- density_curve_df
        qq_list[[length(qq_list) + 1]] <- qq_df
    }

    list(
        summary_df = dplyr::bind_rows(summary_list),
        density_curve_df = if (length(density_curve_list) > 0) {
            dplyr::bind_rows(density_curve_list)
        } else {
            tibble::tibble(
                dataset = character(0),
                pattern_group = character(0),
                x = numeric(0),
                density = numeric(0),
                shape = numeric(0),
                scale = numeric(0),
                hazard_type = character(0)
            )
        },
        qq_df = if (length(qq_list) > 0) {
            dplyr::bind_rows(qq_list)
        } else {
            tibble::tibble(
                dataset = character(0),
                pattern_group = character(0),
                quantile_rank = integer(0),
                probability = numeric(0),
                sample_quantile = numeric(0),
                theoretical_quantile = numeric(0),
                shape = numeric(0),
                scale = numeric(0),
                hazard_type = character(0)
            )
        }
    )
}

## 函数：build_faers_fatality_lookup ####
## 功能：按FAERS的OUTC表汇总病例级致死状态 ####
## 输入：faers_outc，FAERS结局数据框；target_record_ids，需要参与fatality分析的病例ID向量 ####
## 输出：返回包含record_id与fatality_status的结果表 ####
build_faers_fatality_lookup <- function(faers_outc, target_record_ids = NULL) {
    faers_outc_subset <- faers_outc

    if (!is.null(target_record_ids)) {
        target_record_ids <- unique(as.character(target_record_ids))
        faers_outc_subset <- faers_outc_subset %>%
            dplyr::filter(as.character(primaryid) %in% target_record_ids)
    }

    faers_outc_subset %>%
        dplyr::transmute(
            record_id = as.character(primaryid),
            outc_cod = as.character(outc_cod)
        ) %>%
        dplyr::group_by(record_id) %>%
        dplyr::summarise(
            fatality_status = dplyr::case_when(
                any(outc_cod == "DE", na.rm = TRUE) ~ "Fatal",
                any(!is.na(outc_cod) & outc_cod != "", na.rm = TRUE) ~ "Non-Fatal",
                TRUE ~ "Unknown"
            ),
            .groups = "drop"
        )
}

## 函数：build_jader_fatality_lookup ####
## 功能：按JADER的irAE反应记录汇总病例级致死状态 ####
## 输入：jader_reac，JADER反应数据框；irae_pt_keys，标准化后的irAE PT集合；target_record_ids，需要参与fatality分析的病例ID向量 ####
## 输出：返回包含record_id与fatality_status的结果表 ####
build_jader_fatality_lookup <- function(jader_reac, irae_pt_keys, target_record_ids = NULL) {
    jader_reac_subset <- jader_reac

    if (!is.null(target_record_ids)) {
        target_record_ids <- unique(as.character(target_record_ids))
        jader_reac_subset <- jader_reac_subset %>%
            dplyr::filter(as.character(primaryid) %in% target_record_ids)
    }

    pt_key_candidates <- intersect(c("pt_std", "pt_term_en", "pt_term"), names(jader_reac_subset))

    if (length(pt_key_candidates) == 0) {
        stop("jader_reac缺少可用于irAE匹配的PT字段")
    }

    pt_key <- standardize_medra_term(jader_reac_subset[[pt_key_candidates[[1]]]])

    if (length(pt_key_candidates) > 1) {
        for (candidate_col in pt_key_candidates[-1]) {
            pt_key <- dplyr::coalesce(
                pt_key,
                standardize_medra_term(jader_reac_subset[[candidate_col]])
            )
        }
    }

    tibble::tibble(
        record_id = as.character(jader_reac_subset$primaryid),
        pt_key = pt_key,
        outcome = standardize_jader_outcome(jader_reac_subset$outcome)
    ) %>%
        dplyr::filter(pt_key %in% irae_pt_keys) %>%
        dplyr::group_by(record_id) %>%
        dplyr::summarise(
            fatality_status = dplyr::case_when(
                any(outcome == "死亡", na.rm = TRUE) ~ "Fatal",
                any(!is.na(outcome), na.rm = TRUE) ~ "Non-Fatal",
                TRUE ~ "Unknown"
            ),
            .groups = "drop"
        )
}

## 函数：prepare_fatality_case_level ####
## 功能：合并病例级致死状态并生成fatality分析输入表 ####
## 输入：mapped_irae_data，已映射且有irAE的ICI病例；fatality_lookup_df，病例级致死状态表；dataset_name，数据库名称；
##       unknown_as_non_fatal，是否将结局未知的病例计为非致死（TRUE时分母为全部irAE病例） ####
## 输出：返回包含病例级fatality表和信息表的列表 ####
## 返修修改（一审，审稿人2）：原稿中FAERS排除了无结局编码的病例，而JADER将“不明”计为非致死，两库分母不一致。
##       返修后两库统一以全部irAE病例为分母（unknown_as_non_fatal = TRUE），即“记录到死亡的报告比例”；
##       仅以结局已知病例为分母的结果作为敏感性分析，见6.2.12。JADER在原定义下无结局未知病例，结果不变。
prepare_fatality_case_level <- function(mapped_irae_data, fatality_lookup_df, dataset_name, unknown_as_non_fatal = FALSE) {
    fatality_raw_df <- mapped_irae_data %>%
        dplyr::transmute(
            dataset = dataset_name,
            record_id = as.character(record_id),
            primary_organ = as.character(primary_organ),
            pattern_group = as.character(pattern_group)
        ) %>%
        dplyr::left_join(fatality_lookup_df, by = "record_id") %>%
        dplyr::mutate(
            fatality_status = dplyr::coalesce(fatality_status, "Unknown")
        )

    unknown_counted_as_non_fatal_n <- 0L
    if (unknown_as_non_fatal) {
        unknown_counted_as_non_fatal_n <- sum(fatality_raw_df$fatality_status == "Unknown")
        fatality_raw_df <- fatality_raw_df %>%
            dplyr::mutate(fatality_status = ifelse(fatality_status == "Unknown", "Non-Fatal", fatality_status))
    }

    fatality_case_level_df <- fatality_raw_df %>%
        dplyr::filter(fatality_status %in% c("Fatal", "Non-Fatal")) %>%
        dplyr::arrange(factor(pattern_group, levels = pattern_display_order), fatality_status, record_id)

    fatal_n <- sum(fatality_case_level_df$fatality_status == "Fatal", na.rm = TRUE)
    non_fatal_n <- sum(fatality_case_level_df$fatality_status == "Non-Fatal", na.rm = TRUE)
    excluded_unknown_n <- nrow(fatality_raw_df) - nrow(fatality_case_level_df)

    info_df <- tibble::tribble(
        ~section, ~metric, ~value, ~notes,
        "analysis_info", "input_mapped_irae_ici_n", nrow(fatality_raw_df), "Mapped irAE-ICI cases before fatality filtering",
        "analysis_info", "known_fatality_status_n", nrow(fatality_case_level_df), "Cases retained for fatality comparison",
        "analysis_info", "excluded_unknown_fatality_n", excluded_unknown_n, "Cases excluded because fatality status is unavailable",
        "analysis_info", "unknown_counted_as_non_fatal_n", unknown_counted_as_non_fatal_n, "Cases without recorded outcome retained in the denominator as non-fatal",
        "analysis_info", "fatal_n", fatal_n, "Cases classified as Fatal",
        "analysis_info", "non_fatal_n", non_fatal_n, "Cases classified as Non-Fatal"
    )

    list(
        case_level = fatality_case_level_df,
        info_df = info_df
    )
}

## 函数：extract_pairwise_prop_results ####
## 功能：将pairwise.prop.test结果矩阵展开为长表 ####
## 输入：pairwise_test，pairwise.prop.test返回对象；dataset_name，数据库名称 ####
## 输出：返回两两比较结果长表 ####
extract_pairwise_prop_results <- function(pairwise_test, dataset_name) {
    p_value_matrix <- pairwise_test$p.value

    if (is.null(p_value_matrix) || length(p_value_matrix) == 0) {
        return(tibble::tibble())
    }

    result_list <- list()

    for (row_index in seq_len(nrow(p_value_matrix))) {
        for (col_index in seq_len(ncol(p_value_matrix))) {
            current_p <- p_value_matrix[row_index, col_index]

            if (is.na(current_p)) {
                next
            }

            current_group1 <- rownames(p_value_matrix)[row_index]
            current_group2 <- colnames(p_value_matrix)[col_index]

            result_list[[length(result_list) + 1]] <- tibble::tibble(
                dataset = dataset_name,
                section = "pairwise_test",
                group1 = current_group1,
                group2 = current_group2,
                comparison = paste(current_group1, current_group2, sep = " vs "),
                test_name = "pairwise.prop.test",
                p_value = as.numeric(current_p),
                p_adjust_method = "bonferroni",
                significance = build_significance_label(current_p)
            )
        }
    }

    dplyr::bind_rows(result_list)
}

## 函数：build_fatality_summary ####
## 功能：汇总各群fatality比例并执行总体与两两比较 ####
## 输入：fatality_case_level_df，病例级fatality表；info_df，fatality分析信息表；dataset_name，数据库名称；pattern_order，固定pattern顺序 ####
## 输出：返回fatality汇总结果表 ####
build_fatality_summary <- function(fatality_case_level_df, info_df, dataset_name, pattern_order) {
    group_summary_df <- tibble::tibble(pattern_group = pattern_order) %>%
        dplyr::left_join(
            fatality_case_level_df %>%
                dplyr::group_by(pattern_group) %>%
                dplyr::summarise(
                    fatal_n = sum(fatality_status == "Fatal", na.rm = TRUE),
                    non_fatal_n = sum(fatality_status == "Non-Fatal", na.rm = TRUE),
                    total_n = dplyr::n(),
                    .groups = "drop"
                ),
            by = "pattern_group"
        ) %>%
        dplyr::mutate(
            fatal_n = tidyr::replace_na(fatal_n, 0L),
            non_fatal_n = tidyr::replace_na(non_fatal_n, 0L),
            total_n = tidyr::replace_na(total_n, 0L),
            fatality_pct = dplyr::if_else(total_n > 0, fatal_n / total_n * 100, NA_real_),
            section = "group_summary",
            dataset = dataset_name
        ) %>%
        dplyr::select(dataset, section, pattern_group, fatal_n, non_fatal_n, total_n, fatality_pct)

    analysis_groups_df <- group_summary_df %>%
        dplyr::filter(total_n > 0)

    if (nrow(analysis_groups_df) >= 2) {
        contingency_table <- fatality_case_level_df %>%
            dplyr::count(pattern_group, fatality_status) %>%
            tidyr::pivot_wider(
                names_from = fatality_status,
                values_from = n,
                values_fill = 0
            ) %>%
            dplyr::arrange(factor(pattern_group, levels = pattern_order))

        contingency_matrix <- contingency_table %>%
            tibble::column_to_rownames("pattern_group") %>%
            as.matrix()

        if (ncol(contingency_matrix) < 2) {
            global_test_df <- tibble::tibble(
                dataset = dataset_name,
                section = "global_test",
                metric = "Fatality",
                test_name = "No valid global test",
                statistic = NA_real_,
                df = NA_real_,
                p_value = NA_real_,
                significance = NA_character_,
                notes = "Skipped because fatality endpoint has only one observed status"
            )

            pairwise_df <- tibble::tibble(
                dataset = character(0),
                section = character(0),
                group1 = character(0),
                group2 = character(0),
                comparison = character(0),
                test_name = character(0),
                p_value = numeric(0),
                p_adjust_method = character(0),
                significance = character(0)
            )
        } else {
            chi_reference <- suppressWarnings(stats::chisq.test(contingency_matrix))

            if (any(chi_reference$expected < 5)) {
                global_test <- stats::fisher.test(contingency_matrix)
                global_test_df <- tibble::tibble(
                    dataset = dataset_name,
                    section = "global_test",
                    metric = "Fatality",
                    test_name = "Fisher's exact test",
                    statistic = NA_real_,
                    df = NA_real_,
                    p_value = global_test$p.value,
                    significance = build_significance_label(global_test$p.value),
                    notes = "Overall fatality comparison across pattern groups"
                )
            } else {
                global_test_df <- tibble::tibble(
                    dataset = dataset_name,
                    section = "global_test",
                    metric = "Fatality",
                    test_name = "Chi-squared test",
                    statistic = unname(chi_reference$statistic),
                    df = unname(chi_reference$parameter),
                    p_value = chi_reference$p.value,
                    significance = build_significance_label(chi_reference$p.value),
                    notes = "Overall fatality comparison across pattern groups"
                )
            }

            pairwise_test <- stats::pairwise.prop.test(
                x = stats::setNames(analysis_groups_df$fatal_n, analysis_groups_df$pattern_group),
                n = stats::setNames(analysis_groups_df$total_n, analysis_groups_df$pattern_group),
                p.adjust.method = "bonferroni"
            )

            pairwise_df <- extract_pairwise_prop_results(pairwise_test, dataset_name)
        }
    } else {
        global_test_df <- tibble::tibble(
            dataset = dataset_name,
            section = "global_test",
            metric = "Fatality",
            test_name = "No valid global test",
            statistic = NA_real_,
            df = NA_real_,
            p_value = NA_real_,
            significance = NA_character_,
            notes = "Skipped because fewer than two pattern groups have known fatality status"
        )

        pairwise_df <- tibble::tibble(
            dataset = character(0),
            section = character(0),
            group1 = character(0),
            group2 = character(0),
            comparison = character(0),
            test_name = character(0),
            p_value = numeric(0),
            p_adjust_method = character(0),
            significance = character(0)
        )
    }

    info_output_df <- info_df %>%
        dplyr::mutate(dataset = dataset_name) %>%
        dplyr::select(dataset, dplyr::everything())

    dplyr::bind_rows(
        info_output_df,
        group_summary_df,
        global_test_df,
        pairwise_df
    )
}

## 函数：write_dataset_workbook ####
## 功能：将单个数据库的全部分析结果写入Excel工作簿 ####
## 输入：output_file，输出文件路径；cohort_flow_df，队列流程表；pattern_group_n_df，样本量表；soc_spectrum_df，SOC谱表；tto_case_level_df，TTO病例级表；tto_summary_df，TTO汇总表；tto_weibull_summary_df，Weibull汇总表；tto_weibull_density_df，Weibull密度曲线表；tto_weibull_qq_df，Weibull QQ表；fatality_case_level_df，fatality病例级表；fatality_summary_df，fatality汇总表 ####
## 输出：无显式返回值，副作用为写出Excel文件 ####
write_dataset_workbook <- function(
    output_file,
    cohort_flow_df,
    pattern_group_n_df,
    soc_spectrum_df,
    tto_case_level_df,
    tto_summary_df,
    tto_weibull_summary_df,
    tto_weibull_density_df,
    tto_weibull_qq_df,
    fatality_case_level_df,
    fatality_summary_df
) {
    writexl::write_xlsx(
        list(
            "队列流程" = cohort_flow_df,
            "pattern_group_n" = pattern_group_n_df,
            "SOC_spectrum" = soc_spectrum_df,
            "TTO_case_level" = tto_case_level_df,
            "TTO_summary" = tto_summary_df,
            "TTO_weibull_summary" = tto_weibull_summary_df,
            "TTO_weibull_density" = tto_weibull_density_df,
            "TTO_weibull_qq" = tto_weibull_qq_df,
            "fatality_case_level" = fatality_case_level_df,
            "fatality_summary" = fatality_summary_df
        ),
        path = output_file
    )
}


# 4. 读取分群表与基础数据 ####
pattern_assignment_df <- load_pattern_group_assignment(
    input_file = pattern_assignment_file,
    pattern_order = pattern_display_order
)

soc_order_vec <- load_soc_order_from_w_matrix(pattern_assignment_file)

irae_dictionary <- readxl::read_xlsx(irae_dict_file) %>%
    dplyr::transmute(
        pt_key = standardize_medra_term(PT)
    ) %>%
    dplyr::filter(!is.na(pt_key), pt_key != "") %>%
    dplyr::distinct()

faers_solid <- load_named_rdata_object(faers_solid_file, "faers_solid")
faers_outc <- load_named_rdata_object(faers_outc_file, "faers_outc")
jader_solid <- load_named_rdata_object(jader_solid_file, "jader_solid")
jader_reac <- load_named_rdata_object(jader_reac_file, "jader_reac")


# 5. FAERS 分析 ####

## 5.1 队列映射 ####
faers_prepared <- prepare_mapped_ici_dataset(
    data = faers_solid,
    id_col = "primaryid",
    dataset_name = "FAERS",
    assignment_df = pattern_assignment_df,
    pattern_order = pattern_display_order
)

## 5.2 样本量汇总 ####
faers_pattern_group_n <- summarise_pattern_group_n(
    mapped_data = faers_prepared$mapped_data,
    mapped_irae_data = faers_prepared$mapped_irae_data,
    pattern_order = pattern_display_order
)

## 5.3 SOC谱 ####
faers_soc_metadata <- extract_ordered_soc_metadata(
    data = faers_prepared$mapped_data,
    soc_order_vec = soc_order_vec
)

faers_soc_spectrum <- build_soc_spectrum(
    mapped_data = faers_prepared$mapped_data,
    soc_metadata = faers_soc_metadata,
    pattern_group_n_df = faers_pattern_group_n,
    pattern_order = pattern_display_order
) %>%
    dplyr::mutate(dataset = "FAERS") %>%
    dplyr::select(dataset, dplyr::everything())

## 5.4 TTO ####
faers_tto_results <- prepare_tto_case_level(
    mapped_irae_data = faers_prepared$mapped_irae_data,
    dataset_name = "FAERS"
)

faers_tto_summary <- build_tto_summary(
    tto_case_level_df = faers_tto_results$case_level,
    info_df = faers_tto_results$info_df,
    dataset_name = "FAERS",
    pattern_order = pattern_display_order
)

faers_tto_weibull <- build_tto_weibull_results(
    tto_case_level_df = faers_tto_results$case_level,
    dataset_name = "FAERS",
    pattern_order = pattern_display_order
)

## 5.5 fatality ####
faers_fatality_lookup <- build_faers_fatality_lookup(
    faers_outc = faers_outc,
    target_record_ids = faers_prepared$mapped_irae_data$record_id
)

faers_fatality_results <- prepare_fatality_case_level(
    mapped_irae_data = faers_prepared$mapped_irae_data,
    fatality_lookup_df = faers_fatality_lookup,
    dataset_name = "FAERS",
    unknown_as_non_fatal = TRUE
)

faers_fatality_summary <- build_fatality_summary(
    fatality_case_level_df = faers_fatality_results$case_level,
    info_df = faers_fatality_results$info_df,
    dataset_name = "FAERS",
    pattern_order = pattern_display_order
)

## 5.6 导出结果 ####
write_dataset_workbook(
    output_file = faers_output_file,
    cohort_flow_df = faers_prepared$flow_table %>%
        dplyr::mutate(dataset = "FAERS") %>%
        dplyr::select(dataset, dplyr::everything()),
    pattern_group_n_df = faers_pattern_group_n %>%
        dplyr::mutate(dataset = "FAERS") %>%
        dplyr::select(dataset, dplyr::everything()),
    soc_spectrum_df = faers_soc_spectrum,
    tto_case_level_df = faers_tto_results$case_level,
    tto_summary_df = faers_tto_summary,
    tto_weibull_summary_df = faers_tto_weibull$summary_df,
    tto_weibull_density_df = faers_tto_weibull$density_curve_df,
    tto_weibull_qq_df = faers_tto_weibull$qq_df,
    fatality_case_level_df = faers_fatality_results$case_level,
    fatality_summary_df = faers_fatality_summary
)

message("FAERS结果已导出: ", faers_output_file)


# 6. JADER 分析 ####

## 6.1 队列映射 ####
jader_prepared <- prepare_mapped_ici_dataset(
    data = jader_solid,
    id_col = "primaryid",
    dataset_name = "JADER",
    assignment_df = pattern_assignment_df,
    pattern_order = pattern_display_order
)

## 6.2 样本量汇总 ####
jader_pattern_group_n <- summarise_pattern_group_n(
    mapped_data = jader_prepared$mapped_data,
    mapped_irae_data = jader_prepared$mapped_irae_data,
    pattern_order = pattern_display_order
)

## 6.3 SOC谱 ####
jader_soc_metadata <- extract_ordered_soc_metadata(
    data = jader_prepared$mapped_data,
    soc_order_vec = soc_order_vec
)

jader_soc_spectrum <- build_soc_spectrum(
    mapped_data = jader_prepared$mapped_data,
    soc_metadata = jader_soc_metadata,
    pattern_group_n_df = jader_pattern_group_n,
    pattern_order = pattern_display_order
) %>%
    dplyr::mutate(dataset = "JADER") %>%
    dplyr::select(dataset, dplyr::everything())

## 6.4 TTO ####
jader_tto_results <- prepare_tto_case_level(
    mapped_irae_data = jader_prepared$mapped_irae_data,
    dataset_name = "JADER"
)

jader_tto_summary <- build_tto_summary(
    tto_case_level_df = jader_tto_results$case_level,
    info_df = jader_tto_results$info_df,
    dataset_name = "JADER",
    pattern_order = pattern_display_order
)

jader_tto_weibull <- build_tto_weibull_results(
    tto_case_level_df = jader_tto_results$case_level,
    dataset_name = "JADER",
    pattern_order = pattern_display_order
)

## 6.5 fatality ####
jader_fatality_lookup <- build_jader_fatality_lookup(
    jader_reac = jader_reac,
    irae_pt_keys = irae_dictionary$pt_key,
    target_record_ids = jader_prepared$mapped_irae_data$record_id
)

jader_fatality_results <- prepare_fatality_case_level(
    mapped_irae_data = jader_prepared$mapped_irae_data,
    fatality_lookup_df = jader_fatality_lookup,
    dataset_name = "JADER",
    unknown_as_non_fatal = TRUE
)

jader_fatality_summary <- build_fatality_summary(
    fatality_case_level_df = jader_fatality_results$case_level,
    info_df = jader_fatality_results$info_df,
    dataset_name = "JADER",
    pattern_order = pattern_display_order
)

## 6.6 导出结果 ####
write_dataset_workbook(
    output_file = jader_output_file,
    cohort_flow_df = jader_prepared$flow_table %>%
        dplyr::mutate(dataset = "JADER") %>%
        dplyr::select(dataset, dplyr::everything()),
    pattern_group_n_df = jader_pattern_group_n %>%
        dplyr::mutate(dataset = "JADER") %>%
        dplyr::select(dataset, dplyr::everything()),
    soc_spectrum_df = jader_soc_spectrum,
    tto_case_level_df = jader_tto_results$case_level,
    tto_summary_df = jader_tto_summary,
    tto_weibull_summary_df = jader_tto_weibull$summary_df,
    tto_weibull_density_df = jader_tto_weibull$density_curve_df,
    tto_weibull_qq_df = jader_tto_weibull$qq_df,
    fatality_case_level_df = jader_fatality_results$case_level,
    fatality_summary_df = jader_fatality_summary
)

message("JADER结果已导出: ", jader_output_file)
message("========== 癌种群 irAE 表型分析完成 ==========")
