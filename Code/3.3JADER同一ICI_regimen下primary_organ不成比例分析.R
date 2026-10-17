# 同一 ICI_regimen 下 primary_organ 不成比例分析（JADER）####

rm(list = ls())
gc()

source("Code/config.R")
source("Code/mapping_rules.R")
source("Code/disproportionality_analysis.R")

# 0. 并行线程设置 ####
Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1",
    NUMEXPR_NUM_THREADS = "1"
)

# 1. 结果目录与模型参数 ####
script_name <- "3.3JADER同一ICI_regimen下primary_organ不成比例分析"
script_results_path <- file.path(results_path, script_name)

if (!dir.exists(script_results_path)) {
    dir.create(script_results_path, recursive = TRUE, showWarnings = FALSE)
}

message("========== 同一 ICI_regimen 下 primary_organ 不成比例分析开始 ==========")
message(paste("结果目录:", script_results_path))

dataset_configs <- list(
    "JADER" = list(
        file_name = "jader_solid.RData",
        object_name = "jader_solid"
    )
)

model_specs <- list(
    "M0" = list(adjust_vars = NULL),
    "M1" = list(adjust_vars = c("treatment_strategy")),
    "M2" = list(adjust_vars = c("treatment_strategy", "Age", "Sex")),
    "M3" = list(adjust_vars = c("treatment_strategy", "Age", "Sex", "Chemotherapy", "Targeted_therapy"))
)

target_regimen_raw <- c(
    "anti-PD-L1",
    "anti-PD-1",
    "anti-CTLA-4+anti-PD-1"
)
excluded_primary_organs <- c("CUP", "Unmapped", "Unknown_Primary")
irae_level_to_analyze <- "SOC"
min_cases_target <- 3L
min_cases_reference <- 3L
expected_count_threshold <- 5L
within_regimen_fdr_method <- "BH"


# 2. 辅助函数 ####

## 确保目录存在，输入dir_path目录路径，输出为已存在目录的不可见返回值 ####
ensure_dir_exists <- function(dir_path) {
    if (dir.exists(dir_path)) {
        return(invisible(dir_path))
    }

    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)

    if (!dir.exists(dir_path)) {
        stop(paste("无法创建输出目录:", dir_path))
    }

    invisible(dir_path)
}

## 确保输出文件父目录存在，输入file_path文件路径，输出为父目录的不可见返回值 ####
ensure_parent_dir <- function(file_path) {
    ensure_dir_exists(dirname(file_path))
}

## 读取RData指定对象，输入input_path与object_name，输出为目标数据对象 ####
load_named_rdata <- function(input_path, object_name) {
    data_env <- new.env(parent = emptyenv())
    load(input_path, envir = data_env)

    if (!exists(object_name, envir = data_env, inherits = FALSE)) {
        stop(paste("RData中未找到对象:", object_name))
    }

    get(object_name, envir = data_env, inherits = FALSE)
}

## 标准化ICI regimen展示文本，输入x原始regimen字符向量，输出展示用字符向量 ####
format_ici_regimen_label <- function(x) {
    x %>%
        stringr::str_replace_all("anti-PD-1", "anti–PD-1") %>%
        stringr::str_replace_all("anti-PD-L1", "anti–PD-L1") %>%
        stringr::str_replace_all("anti-CTLA-4", "anti–CTLA-4") %>%
        stringr::str_replace_all("\\+", " + ")
}

## 标准化ICI regimen原始值，输入x原始regimen字符向量，输出按成分排序后的稳定raw字符向量 ####
normalize_ici_regimen_raw <- function(x) {
    vapply(x, function(current_regimen) {
        if (is.na(current_regimen) || current_regimen == "") {
            return(NA_character_)
        }

        current_parts <- strsplit(as.character(current_regimen), "\\+", fixed = FALSE)[[1]]
        current_parts <- trimws(current_parts)
        current_parts <- current_parts[current_parts != ""]

        if (length(current_parts) == 0) {
            return(NA_character_)
        }

        paste(sort(unique(current_parts)), collapse = "+")
    }, character(1))
}

## 构建目标regimen映射表，输入为空，输出为raw与display对照表 ####
get_target_regimen_table <- function() {
    tibble::tibble(
        ici_regimen_raw = normalize_ici_regimen_raw(target_regimen_raw),
        regimen_rank = seq_along(target_regimen_raw)
    ) %>%
        mutate(
            ici_regimen_display = format_ici_regimen_label(ici_regimen_raw)
        )
}

## 构建参数数据框，输入params_list参数列表，输出可写入Excel的数据框 ####
build_params_dataframe <- function(params_list) {
    data.frame(
        参数 = names(params_list),
        值 = vapply(params_list, function(param_value) {
            if (is.null(param_value) || length(param_value) == 0) {
                return(NA_character_)
            }

            paste(as.character(param_value), collapse = ", ")
        }, character(1)),
        stringsAsFactors = FALSE
    )
}

## 获取单个SOC缩写，输入irae_name名称，输出为SOC缩写或原始名称 ####
get_soc_abbrev_safe <- function(irae_name) {
    soc_full <- gsub("_", " ", irae_name)
    soc_abbrev <- get_soc_abbrev(soc_full)

    if (is.na(soc_abbrev) || soc_abbrev == "") {
        return(irae_name)
    }

    soc_abbrev
}

## 批量转换SOC缩写，输入irae_vec名称向量，输出为SOC缩写字符向量 ####
abbreviate_soc_vector <- function(irae_vec) {
    vapply(irae_vec, get_soc_abbrev_safe, character(1))
}

## 解析协变量别名，输入adjust_vars和data，输出数据中实际存在的列名 ####
resolve_adjust_var_aliases <- function(adjust_vars, data) {
    if (is.null(adjust_vars) || length(adjust_vars) == 0) {
        return(character(0))
    }

    data_names <- names(data)
    normalized_data_names <- tolower(gsub("[^[:alnum:]]+", "", data_names))

    resolved_vars <- vapply(adjust_vars, function(var_name) {
        if (var_name %in% data_names) {
            return(var_name)
        }

        normalized_var_name <- tolower(gsub("[^[:alnum:]]+", "", var_name))
        matched_index <- match(normalized_var_name, normalized_data_names)

        if (is.na(matched_index)) {
            return(NA_character_)
        }

        data_names[[matched_index]]
    }, character(1))

    resolved_vars
}

## 保留可用协变量，输入adjust_vars和data，输出当前子集中可用且有变异的协变量 ####
retain_varying_adjust_vars <- function(adjust_vars, data) {
    resolved_vars <- resolve_adjust_var_aliases(adjust_vars, data)
    available_vars <- unique(stats::na.omit(resolved_vars))

    retained_vars <- available_vars[purrr::map_lgl(available_vars, function(var_name) {
        current_values <- data[[var_name]]
        length(unique(current_values[!is.na(current_values)])) >= 2
    })]

    retained_vars
}

## 安全执行P值校正，输入p_values与method，输出与原长度一致的校正后P值 ####
safe_p_adjust <- function(p_values, method) {
    adjusted_p <- rep(NA_real_, length(p_values))
    valid_index <- which(!is.na(p_values))

    if (length(valid_index) == 0) {
        return(adjusted_p)
    }

    adjusted_p[valid_index] <- p.adjust(p_values[valid_index], method = method)
    adjusted_p
}

## 补齐调整模型字段，输入details结果表，输出字段完整的结果表 ####
ensure_adjusted_columns <- function(details) {
    numeric_cols <- c(
        "adjusted_OR", "adjusted_OR_lower", "adjusted_OR_upper",
        "adjusted_p_value", "adjusted_p_adjusted",
        "adjusted_p_adjusted_within_regimen"
    )
    logical_cols <- c(
        "adjusted_model_fitted",
        "adjusted_signal_detected",
        "adjusted_signal_detected_within_regimen"
    )

    for (current_col in numeric_cols) {
        if (!current_col %in% names(details)) {
            details[[current_col]] <- NA_real_
        }
    }

    for (current_col in logical_cols) {
        if (!current_col %in% names(details)) {
            details[[current_col]] <- NA
        }
    }

    details
}

## 构建比较标识，输入ici_regimen_display与primary_organ，输出comparison_id ####
build_comparison_id <- function(ici_regimen_display, primary_organ) {
    paste(ici_regimen_display, primary_organ, sep = " | ")
}

## 构建regimen元数据，输入data和dataset_name，输出分析数据与纳入统计信息 ####
build_regimen_metadata <- function(data, dataset_name) {
    target_regimen_table <- get_target_regimen_table()

    filtered_data <- data %>%
        mutate(
            primary_organ = as.character(primary_organ),
            ici_regimen_raw = normalize_ici_regimen_raw(as.character(ici_types))
        ) %>%
        filter(
            has_ici == TRUE,
            !is.na(primary_organ),
            primary_organ != "",
            !primary_organ %in% excluded_primary_organs,
            !is.na(ici_regimen_raw),
            ici_regimen_raw %in% target_regimen_table$ici_regimen_raw
        ) %>%
        left_join(target_regimen_table, by = "ici_regimen_raw") %>%
        mutate(
            comparison_id = build_comparison_id(ici_regimen_display, primary_organ)
        )

    regimen_status <- target_regimen_table %>%
        left_join(
            filtered_data %>%
                count(ici_regimen_raw, ici_regimen_display, regimen_rank, name = "regimen_total_n"),
            by = c("ici_regimen_raw", "ici_regimen_display", "regimen_rank")
        ) %>%
        left_join(
            filtered_data %>%
                distinct(ici_regimen_raw, ici_regimen_display, regimen_rank, primary_organ) %>%
                count(ici_regimen_raw, ici_regimen_display, regimen_rank, name = "n_primary_organs"),
            by = c("ici_regimen_raw", "ici_regimen_display", "regimen_rank")
        ) %>%
        mutate(
            regimen_total_n = tidyr::replace_na(regimen_total_n, 0L),
            n_primary_organs = tidyr::replace_na(n_primary_organs, 0L),
            included_in_main = n_primary_organs >= 2,
            regimen_status = dplyr::case_when(
                regimen_total_n == 0 ~ "排除：无病例",
                !included_in_main ~ "排除：器官数<2",
                TRUE ~ "主分析保留"
            )
        ) %>%
        arrange(regimen_rank)

    missing_regimens <- regimen_status %>%
        filter(regimen_total_n == 0) %>%
        pull(ici_regimen_display)

    if (length(missing_regimens) > 0) {
        stop(
            paste0(
                dataset_name,
                " 中以下ICI regimen无可分析病例: ",
                paste(missing_regimens, collapse = ", ")
            )
        )
    }

    insufficient_regimens <- regimen_status %>%
        filter(!included_in_main) %>%
        pull(ici_regimen_display)

    if (length(insufficient_regimens) > 0) {
        stop(
            paste0(
                dataset_name,
                " 中以下ICI regimen的primary_organ类别数不足2，无法进行不成比例分析: ",
                paste(insufficient_regimens, collapse = ", ")
            )
        )
    }

    organ_status <- filtered_data %>%
        count(ici_regimen_raw, ici_regimen_display, regimen_rank, primary_organ, name = "organ_n") %>%
        group_by(ici_regimen_raw, ici_regimen_display, regimen_rank) %>%
        arrange(desc(organ_n), primary_organ, .by_group = TRUE) %>%
        mutate(
            regimen_total_n = sum(organ_n),
            organ_pct = round(organ_n / regimen_total_n * 100, 2),
            organ_rank = row_number(),
            included_in_main = TRUE,
            comparison_id = build_comparison_id(ici_regimen_display, primary_organ)
        ) %>%
        ungroup() %>%
        arrange(regimen_rank, desc(organ_n), primary_organ)

    analysis_data <- filtered_data %>%
        left_join(
            organ_status %>%
                select(ici_regimen_raw, primary_organ, organ_n, organ_pct),
            by = c("ici_regimen_raw", "primary_organ")
        )

    list(
        analysis_data = analysis_data,
        regimen_status = regimen_status,
        organ_status = organ_status
    )
}

## 追加regimen内FDR结果，输入details结果表，输出追加regimen内FDR后的结果表 ####
append_within_regimen_fdr <- function(details) {
    if (nrow(details) == 0) {
        return(details)
    }

    details %>%
        group_by(ici_regimen_display) %>%
        mutate(
            p_adjusted_within_regimen = safe_p_adjust(p_value, method = within_regimen_fdr_method),
            signal_detected_within_regimen = !is.na(p_adjusted_within_regimen) &
                p_adjusted_within_regimen < 0.05,
            adjusted_p_adjusted_within_regimen = safe_p_adjust(
                adjusted_p_value,
                method = within_regimen_fdr_method
            ),
            adjusted_signal_detected_within_regimen = !is.na(adjusted_p_adjusted_within_regimen) &
                adjusted_p_adjusted_within_regimen < 0.05
        ) %>%
        ungroup() %>%
        arrange(regimen_rank, desc(organ_n), primary_organ, irae)
}

## 获取效应值字段，输入model_name模型名，输出该模型使用的效应值与P值列名 ####
get_effect_columns_for_model <- function(model_name) {
    if (model_name == "M0") {
        return(list(
            effect = "ROR",
            lower = "ROR_lower",
            upper = "ROR_upper",
            p_value = "p_value",
            p_report = "p_adjusted_within_regimen"
        ))
    }

    list(
        effect = "adjusted_OR",
        lower = "adjusted_OR_lower",
        upper = "adjusted_OR_upper",
        p_value = "adjusted_p_value",
        p_report = "adjusted_p_adjusted_within_regimen"
    )
}

## 获取显著性字段，输入model_name模型名，输出该模型使用的显著性列名 ####
get_signal_column_for_model <- function(model_name) {
    if (model_name == "M0") {
        return("signal_detected_within_regimen")
    }

    "adjusted_signal_detected_within_regimen"
}

## 构建宽表矩阵，输入details结果表和value_column列名，输出宽格式数据框 ####
build_wide_matrix_dataframe <- function(details, value_column) {
    row_meta <- details %>%
        distinct(
            ici_regimen_raw,
            ici_regimen_display,
            regimen_rank,
            primary_organ,
            organ_n,
            comparison_id
        )

    irae_levels <- sort(unique(details$irae_abbrev))
    matrix_template <- matrix(
        NA,
        nrow = nrow(row_meta),
        ncol = length(irae_levels),
        dimnames = list(row_meta$comparison_id, irae_levels)
    )

    for (row_index in seq_len(nrow(details))) {
        current_row <- details[row_index, ]
        matrix_template[current_row$comparison_id, current_row$irae_abbrev] <- current_row[[value_column]]
    }

    cbind(
        row_meta,
        as.data.frame(matrix_template, check.names = FALSE)
    )
}

## 汇总模型概览，输入model_results模型结果列表，输出模型级汇总数据框 ####
build_model_overview <- function(model_results) {
    purrr::map_dfr(names(model_results), function(model_name) {
        current_result <- model_results[[model_name]]
        signal_column <- get_signal_column_for_model(model_name)

        tibble::tibble(
            model = model_name,
            adjust_vars_requested = if (is.null(current_result$params$adjust_vars_requested)) {
                "未调整"
            } else {
                paste(current_result$params$adjust_vars_requested, collapse = ", ")
            },
            adjust_vars_used = if (is.null(current_result$params$adjust_vars_used)) {
                "未调整"
            } else {
                paste(current_result$params$adjust_vars_used, collapse = ", ")
            },
            n_regimens = dplyr::n_distinct(current_result$details$ici_regimen_display),
            n_regimen_primary_organs = dplyr::n_distinct(current_result$details$comparison_id),
            n_comparisons = nrow(current_result$details),
            n_fisher_test = sum(current_result$details$fisher_test_used, na.rm = TRUE),
            signal_count = sum(current_result$details[[signal_column]], na.rm = TRUE),
            within_regimen_fdr_method = current_result$params$within_regimen_fdr_method
        )
    })
}

## 汇总regimen信号数，输入model_results模型结果列表，输出按regimen与模型汇总的数据框 ####
build_regimen_signal_summary <- function(model_results) {
    regimen_order <- get_target_regimen_table()$ici_regimen_display

    signal_summary_long <- purrr::map_dfr(names(model_results), function(model_name) {
        signal_column <- get_signal_column_for_model(model_name)

        model_results[[model_name]]$details %>%
            group_by(ici_regimen_display) %>%
            summarise(signal_count = sum(.data[[signal_column]], na.rm = TRUE), .groups = "drop") %>%
            mutate(model = model_name)
    })

    signal_summary_long %>%
        mutate(ici_regimen_display = factor(ici_regimen_display, levels = regimen_order)) %>%
        arrange(ici_regimen_display) %>%
        tidyr::pivot_wider(names_from = model, values_from = signal_count) %>%
        mutate(across(where(is.numeric), ~ tidyr::replace_na(.x, 0)))
}

## 提取Top信号，输入result、model_name和top_n，输出显著SOC前N位数据框 ####
build_top_signals <- function(result, model_name, top_n = 20) {
    signal_column <- get_signal_column_for_model(model_name)
    effect_columns <- get_effect_columns_for_model(model_name)

    result$details %>%
        filter(.data[[signal_column]] %in% TRUE) %>%
        arrange(desc(.data[[effect_columns$effect]]), ici_regimen_display, primary_organ, irae) %>%
        transmute(
            ici_regimen_display = ici_regimen_display,
            primary_organ = primary_organ,
            comparison_id = comparison_id,
            irae_full = irae_full,
            irae = irae_abbrev,
            effect = .data[[effect_columns$effect]],
            lower = .data[[effect_columns$lower]],
            upper = .data[[effect_columns$upper]],
            raw_p_value = .data[[effect_columns$p_value]],
            within_regimen_fdr = .data[[effect_columns$p_report]]
        ) %>%
        slice_head(n = top_n)
}

## 导出模型工作簿，输入result、model_name和output_file，输出写入磁盘的Excel工作簿 ####
export_model_workbook <- function(result, model_name, output_file) {
    details_to_export <- result$details %>%
        arrange(regimen_rank, desc(organ_n), primary_organ, irae)

    fisher_notes <- details_to_export %>%
        filter(fisher_test_used %in% TRUE) %>%
        select(any_of(c(
            "database", "ici_regimen_raw", "ici_regimen_display", "regimen_rank",
            "primary_organ", "comparison_id", "irae", "irae_full", "irae_abbrev",
            "a", "b", "c", "d", "min_expected_count",
            "ROR", "ROR_lower", "ROR_upper",
            "adjusted_OR", "adjusted_OR_lower", "adjusted_OR_upper",
            "p_value", "p_adjusted_within_regimen",
            "adjusted_p_value", "adjusted_p_adjusted_within_regimen",
            "analysis_method", "p_value_method", "result_note"
        )))

    workbook_sheets <- list(
        "详细结果" = details_to_export,
        "ROR矩阵" = build_wide_matrix_dataframe(details_to_export, "ROR"),
        "regimen内FDR信号矩阵" = build_wide_matrix_dataframe(
            details_to_export,
            get_signal_column_for_model(model_name)
        ),
        "Fisher备注" = fisher_notes,
        "分析参数" = build_params_dataframe(result$params)
    )

    if (model_name != "M0") {
        workbook_sheets[["调整后OR矩阵"]] <- build_wide_matrix_dataframe(details_to_export, "adjusted_OR")
        workbook_sheets[["调整后regimen内FDR信号矩阵"]] <- build_wide_matrix_dataframe(
            details_to_export,
            "adjusted_signal_detected_within_regimen"
        )
    }

    ensure_parent_dir(output_file)
    writexl::write_xlsx(workbook_sheets, path = output_file)
}

## 导出汇总工作簿，输入dataset_name、metadata、model_results和output_file，输出summary工作簿 ####
export_summary_workbook <- function(dataset_name, metadata, model_results, output_file) {
    summary_sheets <- list(
        "regimen纳入情况" = metadata$regimen_status,
        "regimen内器官样本量" = metadata$organ_status,
        "模型概览" = build_model_overview(model_results),
        "regimen信号数" = build_regimen_signal_summary(model_results),
        "M0_top_signals" = build_top_signals(model_results$M0, "M0", top_n = 20),
        "M3_top_signals" = build_top_signals(model_results$M3, "M3", top_n = 20),
        "分析参数" = build_params_dataframe(list(
            dataset = dataset_name,
            stratification_var = "ici_regimen",
            target_regimens_raw = target_regimen_raw,
            target_regimens_display = get_target_regimen_table()$ici_regimen_display,
            irae_level = irae_level_to_analyze,
            min_cases_target = min_cases_target,
            min_cases_reference = min_cases_reference,
            expected_count_threshold = expected_count_threshold,
            within_regimen_fdr_method = within_regimen_fdr_method,
            excluded_primary_organs = excluded_primary_organs
        ))
    )

    ensure_parent_dir(output_file)
    writexl::write_xlsx(summary_sheets, path = output_file)
    message(paste(dataset_name, "汇总工作簿已保存:", output_file))
}

## 构建结果对象，输入dataset_name、model_name、details、协变量信息和metadata，输出统一结果列表 ####
build_result_object <- function(dataset_name, model_name, details, requested_adjust_vars, used_adjust_vars, metadata) {
    list(
        details = details,
        params = list(
            dataset = dataset_name,
            model = model_name,
            cancer_var = "primary_organ",
            stratification_var = "ici_regimen",
            irae_level = irae_level_to_analyze,
            filter_ici = TRUE,
            target_regimens_raw = target_regimen_raw,
            target_regimens_display = get_target_regimen_table()$ici_regimen_display,
            excluded_primary_organs = excluded_primary_organs,
            adjust_vars_requested = requested_adjust_vars,
            adjust_vars_used = used_adjust_vars,
            min_cases_target = min_cases_target,
            min_cases_reference = min_cases_reference,
            expected_count_threshold = expected_count_threshold,
            within_regimen_fdr_method = within_regimen_fdr_method,
            n_regimens = nrow(metadata$regimen_status),
            n_regimen_primary_organs = nrow(metadata$organ_status),
            n_comparisons = nrow(details),
            n_fisher_test = sum(details$fisher_test_used, na.rm = TRUE)
        )
    )
}

## 运行单个regimen模型，输入dataset_name、regimen_data、organ_status和requested_adjust_vars，输出当前regimen的结果表 ####
run_single_regimen_model <- function(dataset_name, regimen_data, organ_status, requested_adjust_vars) {
    regimen_meta <- regimen_data %>%
        distinct(ici_regimen_raw, ici_regimen_display, regimen_rank)

    if (nrow(regimen_meta) != 1) {
        stop("run_single_regimen_model 仅支持单个ICI regimen输入")
    }

    current_adjust_vars <- retain_varying_adjust_vars(requested_adjust_vars, regimen_data)

    current_result <- run_disproportionality_analysis(
        data = regimen_data,
        cancer_var = "primary_organ",
        irae_level = irae_level_to_analyze,
        min_cases = min_cases_target,
        min_cases_target = min_cases_target,
        min_cases_reference = min_cases_reference,
        expected_count_threshold = expected_count_threshold,
        filter_ici = FALSE,
        exclude_cancer = NULL,
        adjust_vars = if (length(current_adjust_vars) == 0) NULL else current_adjust_vars,
        p_adjust_method = "none"
    )

    if (is.null(current_result)) {
        stop(paste(dataset_name, regimen_meta$ici_regimen_display, "未返回有效分析结果"))
    }

    current_result$details %>%
        mutate(
            database = dataset_name,
            ici_regimen_raw = regimen_meta$ici_regimen_raw,
            ici_regimen_display = regimen_meta$ici_regimen_display,
            regimen_rank = regimen_meta$regimen_rank,
            primary_organ = cancer_type,
            comparison_id = build_comparison_id(regimen_meta$ici_regimen_display, cancer_type),
            irae_full = gsub("_", " ", irae),
            irae_abbrev = abbreviate_soc_vector(irae),
            adjust_vars_used = if (length(current_adjust_vars) == 0) {
                NA_character_
            } else {
                paste(current_adjust_vars, collapse = ", ")
            }
        ) %>%
        left_join(
            organ_status %>%
                select(primary_organ, organ_n, organ_pct, regimen_total_n),
            by = "primary_organ"
        ) %>%
        ensure_adjusted_columns()
}

## 运行单个数据库全流程，输入dataset_name和dataset_config，输出写入磁盘的分析结果 ####
run_dataset_analysis <- function(dataset_name, dataset_config) {
    message(sprintf("\n========== %s 同一 ICI_regimen 下 primary_organ 分析开始 ==========", dataset_name))

    analysis_data <- load_named_rdata(
        input_path = file.path(output_path, dataset_config$file_name),
        object_name = dataset_config$object_name
    )

    metadata <- build_regimen_metadata(analysis_data, dataset_name)

    message(paste(dataset_name, "目标regimen纳入情况:"))
    print(metadata$regimen_status)

    message(paste(dataset_name, "regimen内器官样本量（前20行）:"))
    print(head(metadata$organ_status, 20))

    model_results <- list()

    for (model_name in names(model_specs)) {
        message(sprintf("\n%s 运行 %s ...", dataset_name, model_name))
        requested_adjust_vars <- model_specs[[model_name]]$adjust_vars
        used_adjust_vars_by_regimen <- list()

        regimen_details <- lapply(target_regimen_raw, function(current_regimen_raw) {
            current_regimen_data <- metadata$analysis_data %>%
                filter(ici_regimen_raw == current_regimen_raw)
            current_organ_status <- metadata$organ_status %>%
                filter(ici_regimen_raw == current_regimen_raw)
            current_used_adjust_vars <- retain_varying_adjust_vars(requested_adjust_vars, current_regimen_data)
            used_adjust_vars_by_regimen[[current_regimen_raw]] <<- current_used_adjust_vars

            run_single_regimen_model(
                dataset_name = dataset_name,
                regimen_data = current_regimen_data,
                organ_status = current_organ_status,
                requested_adjust_vars = requested_adjust_vars
            )
        })

        combined_details <- bind_rows(regimen_details) %>%
            append_within_regimen_fdr()

        current_result <- build_result_object(
            dataset_name = dataset_name,
            model_name = model_name,
            details = combined_details,
            requested_adjust_vars = requested_adjust_vars,
            used_adjust_vars = unique(unlist(used_adjust_vars_by_regimen, use.names = FALSE)),
            metadata = metadata
        )

        output_file <- file.path(
            script_results_path,
            paste0(dataset_name, "_within_ICI_regimen_primary_organ_SOC_", model_name, ".xlsx")
        )

        export_model_workbook(
            result = current_result,
            model_name = model_name,
            output_file = output_file
        )

        message(paste(dataset_name, model_name, "结果已保存:", output_file))
        model_results[[model_name]] <- current_result
    }

    summary_output_file <- file.path(
        script_results_path,
        paste0(dataset_name, "_within_ICI_regimen_primary_organ_SOC_summary.xlsx")
    )

    export_summary_workbook(
        dataset_name = dataset_name,
        metadata = metadata,
        model_results = model_results,
        output_file = summary_output_file
    )

    message(sprintf("========== %s 同一 ICI_regimen 下 primary_organ 分析完成 ==========", dataset_name))
}


# 3. 统计分析与结果导出 ####
purrr::iwalk(dataset_configs, function(dataset_config, dataset_name) {
    run_dataset_analysis(
        dataset_name = dataset_name,
        dataset_config = dataset_config
    )
})

message("\n========== 同一 ICI_regimen 下 primary_organ 不成比例分析全部完成 ==========")
