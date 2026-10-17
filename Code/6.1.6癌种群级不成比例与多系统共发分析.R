# 癌种群级不成比例与多系统共发分析（FAERS + JADER）####

rm(list = ls())
gc()

# 1. 加载配置与依赖 ####
source("Code/config.R")
source("Code/mapping_rules.R")

if (!requireNamespace("dunn.test", quietly = TRUE)) {
    install.packages("dunn.test", quiet = TRUE)
}


# 2. 结果目录与参数 ####
script_name <- "6.1.6癌种群级不成比例与多系统共发分析"
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
jader_solid_file <- file.path(output_path, "jader_solid.RData")

faers_output_file <- file.path(script_results_path, "FAERS_癌种群级ROR与多系统共发汇总.xlsx")
jader_output_file <- file.path(script_results_path, "JADER_癌种群级ROR与多系统共发汇总.xlsx")

pattern_display_order <- c(
    "Ocular",
    "Hematologic",
    "Neurologic",
    "Systemic-like",
    "Visceral-like"
)

ror_min_event_threshold <- 3L

## pattern主导SOC复合指标定义（来源于6.1.1的W矩阵：每个SOC归入其basis weight最大的pattern，且该weight≥0.20；如RENU在Visceral-like为0.34、Systemic-like为0.23，归入Visceral-like）####
## 用于将NMF的多SOC共同偏向在病例级合并为复合指标，与单SOC ROR形成互补层级 ####
pattern_dominant_soc_definitions <- list(
    "Ocular" = c("EYE"),
    "Hematologic" = c("HEME"),
    "Neurologic" = c("NEURO"),
    "Systemic-like" = c("CARD", "SKIN", "ENDO", "MSK"),
    "Visceral-like" = c("HEPA", "GI", "RESP", "RENU")
)

dominant_soc_w_threshold <- 0.20

message("========== 癌种群级不成比例与多系统共发分析开始 ==========")
message("结果目录: ", script_results_path)


# 3. 通用辅助函数 ####

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

    if (!setequal(unique(assignment_df$pattern_group), pattern_order)) {
        stop("分群表中的pattern_group与预设pattern_display_order不一致")
    }

    assignment_df %>%
        dplyr::arrange(factor(pattern_group, levels = pattern_order), primary_organ)
}

## 函数：prepare_mapped_ici_dataset ####
## 功能：筛选ICI病例并合并pattern group标签 ####
## 输入：data，原始分析集；id_col，主键列名；assignment_df，分群对照表；pattern_order，固定pattern顺序 ####
## 输出：返回mapped ICI数据框 ####
prepare_mapped_ici_dataset <- function(data, id_col, assignment_df, pattern_order) {
    data %>%
        dplyr::mutate(
            primary_organ = as.character(primary_organ),
            record_id = as.character(.data[[id_col]])
        ) %>%
        dplyr::filter(has_ici == TRUE) %>%
        dplyr::left_join(
            assignment_df %>% dplyr::select(primary_organ, pattern_group),
            by = "primary_organ"
        ) %>%
        dplyr::filter(!is.na(pattern_group)) %>%
        dplyr::mutate(
            pattern_group = factor(pattern_group, levels = pattern_order)
        ) %>%
        dplyr::arrange(pattern_group, primary_organ, record_id) %>%
        dplyr::mutate(pattern_group = as.character(pattern_group))
}

## 函数：summarise_pattern_group_n ####
## 功能：统计每个pattern group的mapped ICI样本量与mapped irAE-ICI样本量 ####
## 输入：mapped_data，已映射ICI病例；pattern_order，固定pattern顺序 ####
## 输出：返回pattern group样本量汇总表 ####
summarise_pattern_group_n <- function(mapped_data, pattern_order) {
    mapped_ici_n_df <- mapped_data %>%
        dplyr::count(pattern_group, name = "mapped_ici_n")

    mapped_irae_n_df <- mapped_data %>%
        dplyr::filter(has_irae == TRUE) %>%
        dplyr::count(pattern_group, name = "mapped_irae_ici_n")

    tibble::tibble(pattern_group = pattern_order) %>%
        dplyr::left_join(mapped_ici_n_df, by = "pattern_group") %>%
        dplyr::left_join(mapped_irae_n_df, by = "pattern_group") %>%
        dplyr::mutate(
            mapped_ici_n = tidyr::replace_na(mapped_ici_n, 0L),
            mapped_irae_ici_n = tidyr::replace_na(mapped_irae_ici_n, 0L)
        )
}

## 函数：build_target_pt_table ####
## 功能：构建固定纳入的12个重点PT清单 ####
## 输入：无 ####
## 输出：包含PT展示名、变量名、所属SOC的数据框 ####
build_target_pt_table <- function() {
    tibble::tribble(
        ~pt_order, ~irae_display, ~irae_var, ~parent_soc, ~parent_soc_abbrev,
        1, "Pneumonitis", "PT_PNEUMONITIS", "Respiratory, thoracic and mediastinal disorders", "RESP",
        2, "Colitis", "PT_COLITIS", "Gastrointestinal disorders", "GI",
        3, "Myocarditis", "PT_MYOCARDITIS", "Cardiac disorders", "CARD",
        4, "Hypophysitis", "PT_HYPOPHYSITIS", "Endocrine disorders", "ENDO",
        5, "Type 1 diabetes mellitus", "PT_TYPE_1_DIABETES_MELLITUS", "Endocrine disorders", "ENDO",
        6, "Hepatitis", "PT_HEPATITIS", "Hepatobiliary disorders", "HEPA",
        7, "Tubulointerstitial nephritis", "PT_TUBULOINTERSTITIAL_NEPHRITIS", "Renal and urinary disorders", "RENU",
        8, "Myositis", "PT_MYOSITIS", "Musculoskeletal and connective tissue disorders", "MSK",
        9, "Myasthenia gravis", "PT_MYASTHENIA_GRAVIS", "Nervous system disorders", "NEURO",
        10, "Uveitis", "PT_UVEITIS", "Eye disorders", "EYE",
        11, "Pemphigoid", "PT_PEMPHIGOID", "Skin and subcutaneous tissue disorders", "SKIN",
        12, "Stevens-johnson syndrome", "PT_STEVENS-JOHNSON_SYNDROME", "Skin and subcutaneous tissue disorders", "SKIN"
    )
}


# 4. ROR/IC 计算函数 ####

## 函数：compute_ror_per_indicator ####
## 功能：基于2x2列联表计算单个指标在指定pattern group与其余pattern group比较的ROR、95%CI和IC025 ####
## 输入：a，目标group内事件病例数；b，目标group内非事件病例数；c_count，对照组内事件病例数；d_count，对照组内非事件病例数 ####
## 输出：返回包含ROR与95%CI、IC025的列表 ####
compute_ror_per_indicator <- function(a, b, c_count, d_count) {
    a_adj <- a + 0.5
    b_adj <- b + 0.5
    c_adj <- c_count + 0.5
    d_adj <- d_count + 0.5

    log_ror <- log((a_adj * d_adj) / (b_adj * c_adj))
    se_log_ror <- sqrt(1 / a_adj + 1 / b_adj + 1 / c_adj + 1 / d_adj)
    ror <- exp(log_ror)
    ror_lower <- exp(log_ror - 1.96 * se_log_ror)
    ror_upper <- exp(log_ror + 1.96 * se_log_ror)

    if (is.finite(log_ror) && is.finite(se_log_ror) && se_log_ror > 0) {
        z_value <- log_ror / se_log_ror
        p_value <- 2 * stats::pnorm(-abs(z_value))
    } else {
        z_value <- NA_real_
        p_value <- NA_real_
    }

    n_total <- a + b + c_count + d_count
    expected <- (a + b) * (a + c_count) / max(n_total, 1)
    ic <- log2((a + 0.5) / max(expected + 0.5, 0.5))
    ic_var <- 1 / max(a + 0.5, 0.5)
    ic025 <- ic - 3.3 * sqrt(ic_var)

    list(
        a = a,
        b = b,
        c_count = c_count,
        d_count = d_count,
        ror = ror,
        ror_lower = ror_lower,
        ror_upper = ror_upper,
        log_ror = log_ror,
        se_log_ror = se_log_ror,
        z_value = z_value,
        p_value = p_value,
        ic = ic,
        ic025 = ic025
    )
}

## 函数：compute_pattern_group_ror ####
## 功能：对单一二元指标按每个pattern group以"非该group的mapped ICI"为对照计算ROR ####
## 输入：mapped_data，已映射ICI病例；indicator_var，二元指标列名；indicator_label，指标展示名；pattern_order，固定pattern顺序；dataset_name，数据库名称 ####
## 输出：返回ROR结果长表 ####
compute_pattern_group_ror <- function(mapped_data, indicator_var, indicator_label, pattern_order, dataset_name) {
    if (!indicator_var %in% names(mapped_data)) {
        stop("指标列不存在: ", indicator_var)
    }

    indicator_long <- mapped_data %>%
        dplyr::transmute(
            pattern_group = as.character(pattern_group),
            indicator = as.integer(.data[[indicator_var]] == 1)
        )

    result_list <- list()

    for (current_group in pattern_order) {
        current_target <- indicator_long %>%
            dplyr::filter(pattern_group == current_group)

        current_reference <- indicator_long %>%
            dplyr::filter(pattern_group != current_group)

        a_value <- sum(current_target$indicator == 1, na.rm = TRUE)
        b_value <- sum(current_target$indicator == 0, na.rm = TRUE)
        c_value <- sum(current_reference$indicator == 1, na.rm = TRUE)
        d_value <- sum(current_reference$indicator == 0, na.rm = TRUE)

        ror_metrics <- compute_ror_per_indicator(a_value, b_value, c_value, d_value)

        signal_logic <- a_value >= ror_min_event_threshold &
            !is.na(ror_metrics$ror_lower) &
            ror_metrics$ror_lower > 1

        result_list[[length(result_list) + 1]] <- tibble::tibble(
            dataset = dataset_name,
            pattern_group = current_group,
            indicator_var = indicator_var,
            indicator_label = indicator_label,
            target_event_n = a_value,
            target_non_event_n = b_value,
            reference_event_n = c_value,
            reference_non_event_n = d_value,
            ror = ror_metrics$ror,
            ror_lower = ror_metrics$ror_lower,
            ror_upper = ror_metrics$ror_upper,
            p_value = ror_metrics$p_value,
            ic025 = ror_metrics$ic025,
            signal_ror = signal_logic,
            signal_ic025 = !is.na(ror_metrics$ic025) & ror_metrics$ic025 > 0 & a_value >= ror_min_event_threshold
        )
    }

    dplyr::bind_rows(result_list) %>%
        dplyr::arrange(factor(pattern_group, levels = pattern_order), indicator_label)
}

## 函数：build_soc_metadata_for_data ####
## 功能：从分析集提取SOC列信息并附上缩写与展示名 ####
## 输入：data，分析数据框 ####
## 输出：返回SOC列表数据框 ####
build_soc_metadata_for_data <- function(data) {
    soc_vars <- grep("^SOC_", names(data), value = TRUE)

    if (length(soc_vars) == 0) {
        stop("数据中未检测到SOC列")
    }

    tibble::tibble(
        soc_var = soc_vars,
        soc_full = gsub("_", " ", sub("^SOC_", "", soc_vars))
    ) %>%
        dplyr::mutate(
            soc_abbrev = vapply(soc_full, get_soc_abbrev_safe, character(1))
        ) %>%
        dplyr::arrange(soc_full)
}

## 函数：build_soc_pattern_ror_table ####
## 功能：对全部SOC循环计算pattern group ROR ####
## 输入：mapped_data，已映射ICI病例；soc_metadata，SOC元数据；pattern_order，固定pattern顺序；dataset_name，数据库名称 ####
## 输出：返回SOC × pattern group的ROR长表 ####
build_soc_pattern_ror_table <- function(mapped_data, soc_metadata, pattern_order, dataset_name) {
    result_list <- list()

    for (current_index in seq_len(nrow(soc_metadata))) {
        current_var <- soc_metadata$soc_var[[current_index]]
        current_full <- soc_metadata$soc_full[[current_index]]
        current_abbrev <- soc_metadata$soc_abbrev[[current_index]]

        current_result <- compute_pattern_group_ror(
            mapped_data = mapped_data,
            indicator_var = current_var,
            indicator_label = current_abbrev,
            pattern_order = pattern_order,
            dataset_name = dataset_name
        )

        current_result$soc_var <- current_var
        current_result$soc_full <- current_full
        current_result$soc_abbrev <- current_abbrev

        result_list[[length(result_list) + 1]] <- current_result
    }

    dplyr::bind_rows(result_list) %>%
        dplyr::select(dataset, pattern_group, soc_var, soc_full, soc_abbrev, target_event_n, target_non_event_n, reference_event_n, reference_non_event_n, ror, ror_lower, ror_upper, p_value, ic025, signal_ror, signal_ic025)
}

## 函数：build_target_pt_pattern_ror_table ####
## 功能：对12个重点PT循环计算pattern group ROR ####
## 输入：mapped_data，已映射ICI病例；target_pt_table，目标PT清单；pattern_order，固定pattern顺序；dataset_name，数据库名称 ####
## 输出：返回PT × pattern group的ROR长表 ####
build_target_pt_pattern_ror_table <- function(mapped_data, target_pt_table, pattern_order, dataset_name) {
    missing_vars <- setdiff(target_pt_table$irae_var, names(mapped_data))

    if (length(missing_vars) > 0) {
        message(dataset_name, " 缺少PT列，按0补齐: ", paste(missing_vars, collapse = ", "))
        for (current_var in missing_vars) {
            mapped_data[[current_var]] <- 0L
        }
    }

    result_list <- list()

    for (current_index in seq_len(nrow(target_pt_table))) {
        current_var <- target_pt_table$irae_var[[current_index]]
        current_label <- target_pt_table$irae_display[[current_index]]
        current_order <- target_pt_table$pt_order[[current_index]]
        current_parent <- target_pt_table$parent_soc_abbrev[[current_index]]

        current_result <- compute_pattern_group_ror(
            mapped_data = mapped_data,
            indicator_var = current_var,
            indicator_label = current_label,
            pattern_order = pattern_order,
            dataset_name = dataset_name
        )

        current_result$pt_order <- current_order
        current_result$parent_soc_abbrev <- current_parent

        result_list[[length(result_list) + 1]] <- current_result
    }

    dplyr::bind_rows(result_list) %>%
        dplyr::select(dataset, pattern_group, pt_order, indicator_var, indicator_label, parent_soc_abbrev, target_event_n, target_non_event_n, reference_event_n, reference_non_event_n, ror, ror_lower, ror_upper, p_value, ic025, signal_ror, signal_ic025) %>%
        dplyr::rename(irae_var = indicator_var, irae_display = indicator_label) %>%
        dplyr::arrange(factor(pattern_group, levels = pattern_order), pt_order)
}

## 函数：build_pattern_dominant_soc_indicator ####
## 功能：根据pattern_dominant_soc_definitions为mapped ICI数据集追加每个pattern的主导SOC复合二元指标 ####
## 输入：mapped_data，已映射ICI病例；soc_metadata，SOC元数据；dominant_soc_def，pattern→SOC缩写向量列表 ####
## 输出：返回带有composite_<pattern>列的mapped_data及实际可用SOC映射表 ####
build_pattern_dominant_soc_indicator <- function(mapped_data, soc_metadata, dominant_soc_def) {
    soc_lookup <- soc_metadata %>%
        dplyr::select(soc_var, soc_abbrev) %>%
        dplyr::filter(!is.na(soc_abbrev), soc_abbrev != "")

    available_def_list <- list()
    for (current_pattern in names(dominant_soc_def)) {
        current_socs <- dominant_soc_def[[current_pattern]]
        current_match <- soc_lookup %>%
            dplyr::filter(soc_abbrev %in% current_socs)

        composite_col <- paste0("composite_dom_", make.names(current_pattern))

        if (nrow(current_match) == 0) {
            mapped_data[[composite_col]] <- 0L
            available_def_list[[length(available_def_list) + 1]] <- tibble::tibble(
                pattern_group = current_pattern,
                composite_var = composite_col,
                soc_abbrev_used = NA_character_,
                soc_var_used = NA_character_
            )
            next
        }

        current_var_set <- current_match$soc_var
        composite_values <- mapped_data %>%
            dplyr::select(dplyr::all_of(current_var_set)) %>%
            dplyr::mutate(dplyr::across(dplyr::everything(), ~ tidyr::replace_na(as.integer(.x), 0L))) %>%
            as.matrix() %>%
            (function(x) as.integer(rowSums(x) > 0))()

        mapped_data[[composite_col]] <- composite_values

        available_def_list[[length(available_def_list) + 1]] <- tibble::tibble(
            pattern_group = current_pattern,
            composite_var = composite_col,
            soc_abbrev_used = paste(current_match$soc_abbrev, collapse = " | "),
            soc_var_used = paste(current_match$soc_var, collapse = " | ")
        )
    }

    list(
        mapped_data = mapped_data,
        composite_def_df = dplyr::bind_rows(available_def_list)
    )
}

## 函数：compute_pattern_dominant_soc_ror ####
## 功能：对每个pattern以其自身的主导SOC复合指标为focal二元变量，计算该pattern相对其它pattern的ROR ####
## 输入：mapped_data，已映射ICI病例（含composite_dom_列）；composite_def_df，复合指标定义表；pattern_order，固定pattern顺序；dataset_name，数据库名称 ####
## 输出：返回pattern × 主导SOC复合的ROR长表 ####
compute_pattern_dominant_soc_ror <- function(mapped_data, composite_def_df, pattern_order, dataset_name) {
    result_list <- list()

    for (current_pattern in pattern_order) {
        current_def <- composite_def_df %>%
            dplyr::filter(pattern_group == current_pattern)

        if (nrow(current_def) == 0) {
            next
        }

        composite_col <- current_def$composite_var[[1]]
        soc_used <- current_def$soc_abbrev_used[[1]]

        if (!composite_col %in% names(mapped_data)) {
            next
        }

        target_df <- mapped_data %>%
            dplyr::filter(pattern_group == current_pattern)
        reference_df <- mapped_data %>%
            dplyr::filter(pattern_group != current_pattern)

        a_value <- sum(target_df[[composite_col]] == 1L, na.rm = TRUE)
        b_value <- sum(target_df[[composite_col]] == 0L, na.rm = TRUE)
        c_value <- sum(reference_df[[composite_col]] == 1L, na.rm = TRUE)
        d_value <- sum(reference_df[[composite_col]] == 0L, na.rm = TRUE)

        ror_metrics <- compute_ror_per_indicator(a_value, b_value, c_value, d_value)

        signal_logic <- a_value >= ror_min_event_threshold &
            !is.na(ror_metrics$ror_lower) &
            ror_metrics$ror_lower > 1

        target_total <- a_value + b_value
        reference_total <- c_value + d_value

        result_list[[length(result_list) + 1]] <- tibble::tibble(
            dataset = dataset_name,
            pattern_group = current_pattern,
            dominant_soc = soc_used,
            target_event_n = a_value,
            target_total_n = target_total,
            target_rate_pct = dplyr::if_else(target_total > 0, a_value / target_total * 100, NA_real_),
            reference_event_n = c_value,
            reference_total_n = reference_total,
            reference_rate_pct = dplyr::if_else(reference_total > 0, c_value / reference_total * 100, NA_real_),
            ror = ror_metrics$ror,
            ror_lower = ror_metrics$ror_lower,
            ror_upper = ror_metrics$ror_upper,
            p_value = ror_metrics$p_value,
            ic025 = ror_metrics$ic025,
            signal_ror = signal_logic,
            signal_ic025 = !is.na(ror_metrics$ic025) & ror_metrics$ic025 > 0 & a_value >= ror_min_event_threshold
        )
    }

    dplyr::bind_rows(result_list) %>%
        dplyr::arrange(factor(pattern_group, levels = pattern_order))
}


# 5. 多系统共发函数 ####

## 函数：build_case_level_organ_count ####
## 功能：对每个mapped ICI病例计算distinct irAE SOC数与distinct重点PT数 ####
## 输入：mapped_data，已映射ICI病例；soc_metadata，SOC元数据；target_pt_table，重点PT清单；dataset_name，数据库名称 ####
## 输出：返回病例级SOC/PT共发数数据框 ####
build_case_level_organ_count <- function(mapped_data, soc_metadata, target_pt_table, dataset_name) {
    soc_vars <- soc_metadata$soc_var
    pt_vars <- target_pt_table$irae_var

    missing_pt_vars <- setdiff(pt_vars, names(mapped_data))
    if (length(missing_pt_vars) > 0) {
        for (current_var in missing_pt_vars) {
            mapped_data[[current_var]] <- 0L
        }
    }

    soc_matrix <- as.matrix(mapped_data[, soc_vars, drop = FALSE])
    pt_matrix <- as.matrix(mapped_data[, pt_vars, drop = FALSE])
    soc_matrix[is.na(soc_matrix)] <- 0
    pt_matrix[is.na(pt_matrix)] <- 0

    tibble::tibble(
        dataset = dataset_name,
        record_id = as.character(mapped_data$record_id),
        primary_organ = as.character(mapped_data$primary_organ),
        pattern_group = as.character(mapped_data$pattern_group),
        has_irae = as.logical(mapped_data$has_irae),
        n_distinct_soc = as.integer(rowSums(soc_matrix > 0)),
        n_distinct_target_pt = as.integer(rowSums(pt_matrix > 0))
    )
}

## 函数：summarise_organ_count_by_pattern ####
## 功能：按pattern group汇总distinct SOC/PT数的描述统计 ####
## 输入：case_count_df，病例级共发数数据框；pattern_order，固定pattern顺序；dataset_name，数据库名称 ####
## 输出：返回共发数描述统计长表 ####
summarise_organ_count_by_pattern <- function(case_count_df, pattern_order, dataset_name) {
    irae_only_df <- case_count_df %>%
        dplyr::filter(has_irae == TRUE)

    summary_long_list <- list()

    for (current_metric in c("n_distinct_soc", "n_distinct_target_pt")) {
        current_metric_label <- ifelse(
            current_metric == "n_distinct_soc",
            "Distinct irAE SOC count (mapped irAE-ICI cases)",
            "Distinct target PT count (mapped irAE-ICI cases)"
        )

        current_summary <- tibble::tibble(pattern_group = pattern_order) %>%
            dplyr::left_join(
                irae_only_df %>%
                    dplyr::group_by(pattern_group) %>%
                    dplyr::summarise(
                        n_cases = dplyr::n(),
                        median_value = stats::median(.data[[current_metric]], na.rm = TRUE),
                        q1_value = stats::quantile(.data[[current_metric]], probs = 0.25, na.rm = TRUE, names = FALSE),
                        q3_value = stats::quantile(.data[[current_metric]], probs = 0.75, na.rm = TRUE, names = FALSE),
                        mean_value = mean(.data[[current_metric]], na.rm = TRUE),
                        sd_value = stats::sd(.data[[current_metric]], na.rm = TRUE),
                        max_value = suppressWarnings(max(.data[[current_metric]], na.rm = TRUE)),
                        .groups = "drop"
                    ),
                by = "pattern_group"
            ) %>%
            dplyr::mutate(
                dataset = dataset_name,
                metric = current_metric,
                metric_label = current_metric_label,
                n_cases = tidyr::replace_na(n_cases, 0L)
            ) %>%
            dplyr::select(dataset, metric, metric_label, pattern_group, n_cases, median_value, q1_value, q3_value, mean_value, sd_value, max_value)

        summary_long_list[[length(summary_long_list) + 1]] <- current_summary
    }

    dplyr::bind_rows(summary_long_list)
}

## 函数：test_organ_count_diff ####
## 功能：对每个共发数指标在pattern group之间执行Kruskal-Wallis与Dunn事后检验 ####
## 输入：case_count_df，病例级共发数数据框；pattern_order，固定pattern顺序；dataset_name，数据库名称 ####
## 输出：返回包含全局检验和两两比较的列表 ####
test_organ_count_diff <- function(case_count_df, pattern_order, dataset_name) {
    irae_only_df <- case_count_df %>%
        dplyr::filter(has_irae == TRUE) %>%
        dplyr::filter(pattern_group %in% pattern_order)

    global_list <- list()
    pairwise_list <- list()

    for (current_metric in c("n_distinct_soc", "n_distinct_target_pt")) {
        current_label <- ifelse(
            current_metric == "n_distinct_soc",
            "Distinct irAE SOC count",
            "Distinct target PT count"
        )

        current_df <- irae_only_df %>%
            dplyr::transmute(
                pattern_group = pattern_group,
                value = as.numeric(.data[[current_metric]])
            ) %>%
            dplyr::filter(!is.na(value))

        nonzero_groups <- current_df %>%
            dplyr::group_by(pattern_group) %>%
            dplyr::summarise(group_n = dplyr::n(), .groups = "drop") %>%
            dplyr::filter(group_n > 0)

        if (nrow(nonzero_groups) >= 2 && nrow(current_df) > 0) {
            kw_result <- stats::kruskal.test(value ~ pattern_group, data = current_df)

            global_list[[length(global_list) + 1]] <- tibble::tibble(
                dataset = dataset_name,
                metric = current_metric,
                metric_label = current_label,
                test_name = "Kruskal-Wallis",
                statistic = unname(kw_result$statistic),
                df = unname(kw_result$parameter),
                p_value = kw_result$p.value,
                significance = build_significance_label(kw_result$p.value),
                notes = "Overall comparison of distinct count across pattern groups"
            )

            dunn_result <- suppressWarnings(
                dunn.test::dunn.test(
                    x = current_df$value,
                    g = current_df$pattern_group,
                    method = "bonferroni",
                    kw = FALSE,
                    list = TRUE,
                    table = FALSE
                )
            )

            current_pairwise <- tibble::tibble(
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
                    metric = current_metric,
                    metric_label = current_label,
                    test_name = "Dunn test",
                    p_adjust_method = "bonferroni",
                    significance = build_significance_label(p_value),
                    comparison = paste(group1, group2, sep = " vs ")
                ) %>%
                dplyr::select(dataset, metric, metric_label, group1, group2, comparison, test_name, z_statistic, p_value, p_adjust_method, significance)

            pairwise_list[[length(pairwise_list) + 1]] <- current_pairwise
        } else {
            global_list[[length(global_list) + 1]] <- tibble::tibble(
                dataset = dataset_name,
                metric = current_metric,
                metric_label = current_label,
                test_name = "Skipped",
                statistic = NA_real_,
                df = NA_real_,
                p_value = NA_real_,
                significance = NA_character_,
                notes = "Insufficient pattern groups with data"
            )
        }
    }

    list(
        global_df = dplyr::bind_rows(global_list),
        pairwise_df = if (length(pairwise_list) > 0) dplyr::bind_rows(pairwise_list) else tibble::tibble()
    )
}


# 6. 工作簿写出 ####

## 函数：write_stage_b_workbook ####
## 功能：将单数据库的全部分析结果按中文sheet写出为Excel工作簿 ####
## 输入：output_file，输出路径；workbook_sheets，命名列表 ####
## 输出：无显式返回值 ####
write_stage_b_workbook <- function(output_file, workbook_sheets) {
    writexl::write_xlsx(workbook_sheets, path = output_file)
}


# 7. 读取分群表与基础数据 ####
pattern_assignment_df <- load_pattern_group_assignment(
    input_file = pattern_assignment_file,
    pattern_order = pattern_display_order
)

faers_solid <- load_named_rdata_object(faers_solid_file, "faers_solid")
jader_solid <- load_named_rdata_object(jader_solid_file, "jader_solid")

target_pt_table <- build_target_pt_table()


# 8. FAERS 分析 ####

## 8.1 队列映射 ####
faers_mapped <- prepare_mapped_ici_dataset(
    data = faers_solid,
    id_col = "primaryid",
    assignment_df = pattern_assignment_df,
    pattern_order = pattern_display_order
)

faers_pattern_group_n <- summarise_pattern_group_n(faers_mapped, pattern_display_order) %>%
    dplyr::mutate(dataset = "FAERS") %>%
    dplyr::select(dataset, dplyr::everything())

## 8.2 SOC 级 ROR ####
faers_soc_metadata <- build_soc_metadata_for_data(faers_mapped)

faers_soc_ror <- build_soc_pattern_ror_table(
    mapped_data = faers_mapped,
    soc_metadata = faers_soc_metadata,
    pattern_order = pattern_display_order,
    dataset_name = "FAERS"
)

## 8.3 重点PT 级 ROR ####
faers_target_pt_ror <- build_target_pt_pattern_ror_table(
    mapped_data = faers_mapped,
    target_pt_table = target_pt_table,
    pattern_order = pattern_display_order,
    dataset_name = "FAERS"
)

## 8.4 主导SOC复合指标ROR ####
faers_dom_indicator_result <- build_pattern_dominant_soc_indicator(
    mapped_data = faers_mapped,
    soc_metadata = faers_soc_metadata,
    dominant_soc_def = pattern_dominant_soc_definitions
)

faers_mapped <- faers_dom_indicator_result$mapped_data
faers_composite_def_df <- faers_dom_indicator_result$composite_def_df %>%
    dplyr::mutate(dataset = "FAERS") %>%
    dplyr::select(dataset, dplyr::everything())

faers_dominant_soc_ror <- compute_pattern_dominant_soc_ror(
    mapped_data = faers_mapped,
    composite_def_df = faers_dom_indicator_result$composite_def_df,
    pattern_order = pattern_display_order,
    dataset_name = "FAERS"
)

## 8.5 多系统共发 ####
faers_case_count <- build_case_level_organ_count(
    mapped_data = faers_mapped,
    soc_metadata = faers_soc_metadata,
    target_pt_table = target_pt_table,
    dataset_name = "FAERS"
)

faers_organ_count_summary <- summarise_organ_count_by_pattern(
    case_count_df = faers_case_count,
    pattern_order = pattern_display_order,
    dataset_name = "FAERS"
)

faers_organ_count_tests <- test_organ_count_diff(
    case_count_df = faers_case_count,
    pattern_order = pattern_display_order,
    dataset_name = "FAERS"
)

## 8.6 导出 ####
write_stage_b_workbook(
    output_file = faers_output_file,
    workbook_sheets = list(
        "pattern样本量" = faers_pattern_group_n,
        "SOC_ROR" = faers_soc_ror,
        "重点PT_ROR" = faers_target_pt_ror,
        "主导SOC复合定义" = faers_composite_def_df,
        "主导SOC复合ROR" = faers_dominant_soc_ror,
        "多系统共发_病例级" = faers_case_count,
        "多系统共发_汇总" = faers_organ_count_summary,
        "多系统共发_全局检验" = faers_organ_count_tests$global_df,
        "多系统共发_两两比较" = faers_organ_count_tests$pairwise_df
    )
)

message("FAERS结果已导出: ", faers_output_file)


# 9. JADER 分析 ####

## 9.1 队列映射 ####
jader_mapped <- prepare_mapped_ici_dataset(
    data = jader_solid,
    id_col = "primaryid",
    assignment_df = pattern_assignment_df,
    pattern_order = pattern_display_order
)

jader_pattern_group_n <- summarise_pattern_group_n(jader_mapped, pattern_display_order) %>%
    dplyr::mutate(dataset = "JADER") %>%
    dplyr::select(dataset, dplyr::everything())

## 9.2 SOC 级 ROR ####
jader_soc_metadata <- build_soc_metadata_for_data(jader_mapped)

jader_soc_ror <- build_soc_pattern_ror_table(
    mapped_data = jader_mapped,
    soc_metadata = jader_soc_metadata,
    pattern_order = pattern_display_order,
    dataset_name = "JADER"
)

## 9.3 重点PT 级 ROR ####
jader_target_pt_ror <- build_target_pt_pattern_ror_table(
    mapped_data = jader_mapped,
    target_pt_table = target_pt_table,
    pattern_order = pattern_display_order,
    dataset_name = "JADER"
)

## 9.4 主导SOC复合指标ROR ####
jader_dom_indicator_result <- build_pattern_dominant_soc_indicator(
    mapped_data = jader_mapped,
    soc_metadata = jader_soc_metadata,
    dominant_soc_def = pattern_dominant_soc_definitions
)

jader_mapped <- jader_dom_indicator_result$mapped_data
jader_composite_def_df <- jader_dom_indicator_result$composite_def_df %>%
    dplyr::mutate(dataset = "JADER") %>%
    dplyr::select(dataset, dplyr::everything())

jader_dominant_soc_ror <- compute_pattern_dominant_soc_ror(
    mapped_data = jader_mapped,
    composite_def_df = jader_dom_indicator_result$composite_def_df,
    pattern_order = pattern_display_order,
    dataset_name = "JADER"
)

## 9.5 多系统共发 ####
jader_case_count <- build_case_level_organ_count(
    mapped_data = jader_mapped,
    soc_metadata = jader_soc_metadata,
    target_pt_table = target_pt_table,
    dataset_name = "JADER"
)

jader_organ_count_summary <- summarise_organ_count_by_pattern(
    case_count_df = jader_case_count,
    pattern_order = pattern_display_order,
    dataset_name = "JADER"
)

jader_organ_count_tests <- test_organ_count_diff(
    case_count_df = jader_case_count,
    pattern_order = pattern_display_order,
    dataset_name = "JADER"
)

## 9.6 导出 ####
write_stage_b_workbook(
    output_file = jader_output_file,
    workbook_sheets = list(
        "pattern样本量" = jader_pattern_group_n,
        "SOC_ROR" = jader_soc_ror,
        "重点PT_ROR" = jader_target_pt_ror,
        "主导SOC复合定义" = jader_composite_def_df,
        "主导SOC复合ROR" = jader_dominant_soc_ror,
        "多系统共发_病例级" = jader_case_count,
        "多系统共发_汇总" = jader_organ_count_summary,
        "多系统共发_全局检验" = jader_organ_count_tests$global_df,
        "多系统共发_两两比较" = jader_organ_count_tests$pairwise_df
    )
)

message("JADER结果已导出: ", jader_output_file)
message("========== 癌种群级不成比例与多系统共发分析完成 ==========")
