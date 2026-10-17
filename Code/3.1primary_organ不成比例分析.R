# primary_organ不成比例分析（FAERS + JADER）####

## 加载配置与函数####
source("Code/config.R")
source("Code/mapping_rules.R")

## 并行线程设置####
Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1",
    NUMEXPR_NUM_THREADS = "1"
)

## 结果目录设置####
script_name <- "3.1primary_organ不成比例分析"
script_results_path <- file.path(results_path, script_name)

if (!dir.exists(script_results_path)) {
    dir.create(script_results_path, recursive = TRUE, showWarnings = FALSE)
}

message("========== primary_organ不成比例分析开始 ==========")
message(paste("结果目录:", script_results_path))

## 数据集配置####
dataset_configs <- list(
    "FAERS" = list(
        file_name = "faers_malignant_final.RData",
        object_name = "faers_malignant_final"
    ),
    "JADER" = list(
        file_name = "jader_malignant_final.RData",
        object_name = "jader_malignant_final"
    )
)

## 模型配置####
model_specs <- list(
    "M0" = list(adjust_vars = NULL),
    "M1" = list(adjust_vars = c("treatment_strategy")),
    "M2" = list(adjust_vars = c("treatment_strategy", "Age", "Sex")),
    "M3" = list(adjust_vars = c("treatment_strategy", "Age", "Sex", "Chemotherapy", "Targeted_therapy"))
)


## 辅助函数####

### 确保目录存在，必要时递归创建并校验创建结果####
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

### 确保输出文件的父目录存在####
ensure_parent_dir <- function(file_path) {
    ensure_dir_exists(dirname(file_path))
}


### 获取SOC缩写，无法映射时返回原始名称####
get_soc_abbrev_safe <- function(irae_name) {
    soc_full <- gsub("_", " ", irae_name)
    soc_abbrev <- get_soc_abbrev(soc_full)

    if (is.na(soc_abbrev)) {
        return(irae_name)
    }

    soc_abbrev
}

### 将SOC名称向量批量转换为缩写####
abbreviate_soc_vector <- function(irae_vec) {
    vapply(irae_vec, get_soc_abbrev_safe, character(1))
}

### 获取详情表中用于展示的SOC缩写列####
get_detail_irae_display <- function(details) {
    if ("irae_abbrev" %in% names(details)) {
        return(details$irae_abbrev)
    }

    abbreviate_soc_vector(details$irae)
}

### 获取详情表中的SOC全称列，缺失时由irae字段回推####
get_detail_irae_full <- function(details) {
    if ("irae_full" %in% names(details)) {
        return(details$irae_full)
    }

    gsub("_", " ", details$irae)
}

### 构建各模型结果工作簿路径####
get_model_result_paths <- function(model_specs, results_path, prefix = "FAERS_primary_organ_SOC_") {
    stats::setNames(
        file.path(
            results_path,
            paste0(prefix, names(model_specs), ".xlsx")
        ),
        names(model_specs)
    )
}

### 从参数工作表读取指定参数值，缺失时返回默认值####
read_param_value <- function(params_df, param_name, default = NA_character_) {
    matched_value <- params_df %>%
        filter(参数 == param_name) %>%
        pull(值)

    if (length(matched_value) == 0) {
        return(default)
    }

    matched_value[[1]]
}

### 将参数值转换为数值，无法转换时返回默认值####
parse_numeric_param <- function(value, default = NA_real_) {
    if (length(value) == 0 || is.na(value) || value == "") {
        return(default)
    }

    parsed_value <- suppressWarnings(as.numeric(value))

    if (is.na(parsed_value)) {
        return(default)
    }

    parsed_value
}

### 将参数值转换为逻辑值，无法转换时返回默认值####
parse_logical_param <- function(value, default = NA) {
    if (length(value) == 0 || is.na(value) || value == "") {
        return(default)
    }

    value_lower <- tolower(as.character(value))

    if (value_lower %in% c("true", "t", "1")) {
        return(TRUE)
    }

    if (value_lower %in% c("false", "f", "0")) {
        return(FALSE)
    }

    default
}

### 将逗号分隔的参数文本还原为字符向量####
parse_character_vector_param <- function(value) {
    if (length(value) == 0 || is.na(value) || value == "") {
        return(NULL)
    }

    trimws(strsplit(as.character(value), ",", fixed = TRUE)[[1]])
}

### 从模型结果Excel读取details和params，供后续汇总与绘图复用####
load_result_from_excel <- function(input_file) {
    details <- readxl::read_xlsx(input_file, sheet = "详细结果")
    params_df <- readxl::read_xlsx(input_file, sheet = "分析参数")

    character_cols <- intersect(
        c(
            "cancer_type", "irae", "irae_var", "irae_full", "irae_abbrev",
            "analysis_method", "p_value_method", "result_note"
        ),
        names(details)
    )
    logical_cols <- intersect(
        c(
            "signal_detected", "adjusted_signal_detected", "converged",
            "fisher_test_used", "adjusted_model_fitted"
        ),
        names(details)
    )
    numeric_cols <- setdiff(names(details), c(character_cols, logical_cols))

    if (length(character_cols) > 0) {
        details[character_cols] <- lapply(details[character_cols], as.character)
    }

    if (length(logical_cols) > 0) {
        details[logical_cols] <- lapply(details[logical_cols], as.logical)
    }

    if (length(numeric_cols) > 0) {
        details[numeric_cols] <- lapply(details[numeric_cols], function(column) {
            suppressWarnings(as.numeric(column))
        })
    }

    min_cases <- parse_numeric_param(read_param_value(params_df, "min_cases"))

    params <- list(
        cancer_var = read_param_value(params_df, "cancer_var"),
        irae_level = read_param_value(params_df, "irae_level"),
        n_cancers = parse_numeric_param(read_param_value(params_df, "n_cancers")),
        n_iraes = parse_numeric_param(read_param_value(params_df, "n_iraes")),
        filter_ici = parse_logical_param(read_param_value(params_df, "filter_ici"), default = FALSE),
        exclude_cancer = parse_character_vector_param(read_param_value(params_df, "exclude_cancer", default = "")),
        adjust_vars = parse_character_vector_param(read_param_value(params_df, "adjust_vars", default = "")),
        min_cases = min_cases,
        expected_count_threshold = parse_numeric_param(
            read_param_value(params_df, "expected_count_threshold", default = ""),
            default = NA_real_
        ),
        standard_test_method = read_param_value(params_df, "standard_test_method", default = NA_character_),
        fisher_test_method = read_param_value(params_df, "fisher_test_method", default = NA_character_),
        p_adjust_method = read_param_value(params_df, "p_adjust_method"),
        n_processed = parse_numeric_param(read_param_value(params_df, "n_processed")),
        n_skipped = parse_numeric_param(read_param_value(params_df, "n_skipped")),
        n_fisher_test = parse_numeric_param(read_param_value(params_df, "n_fisher_test"), default = 0)
    )

    list(
        details = tibble::as_tibble(details),
        params = params
    )
}

### 检查模型结果Excel是否已包含当前期望频数分流规则所需字段####
model_result_excel_is_compatible <- function(input_file) {
    if (!file.exists(input_file)) {
        return(FALSE)
    }

    tryCatch(
        {
            detail_columns <- names(readxl::read_xlsx(input_file, sheet = "详细结果", n_max = 0))
            params_df <- readxl::read_xlsx(input_file, sheet = "分析参数")
            param_names <- params_df$参数

            required_detail_columns <- c(
                "analysis_method", "p_value_method", "result_note",
                "fisher_test_used", "min_expected_count"
            )
            required_param_names <- c(
                "expected_count_threshold", "n_fisher_test",
                "standard_test_method", "fisher_test_method"
            )

            all(required_detail_columns %in% detail_columns) &&
                all(required_param_names %in% param_names)
        },
        error = function(e) {
            FALSE
        }
    )
}

### 从Excel工作簿载入各模型结果，供后续汇总与绘图使用####
load_model_results <- function(model_specs, results_path, prefix = "FAERS_primary_organ_SOC_") {
    model_result_paths <- get_model_result_paths(model_specs, results_path, prefix)
    all_model_files_exist <- all(file.exists(unname(model_result_paths)))
    all_model_files_compatible <- all_model_files_exist &&
        all(purrr::map_lgl(model_result_paths, model_result_excel_is_compatible))

    if (!all_model_files_compatible) {
        stop("模型结果Excel工作簿缺失或与当前字段要求不兼容，请先准备好可用的Excel结果文件")
    }

    model_results_loaded <- purrr::map(model_result_paths, load_result_from_excel)
    message("模型结果已从Excel载入，后续汇总和绘图仅使用Excel结果")

    model_results_loaded
}

### 将协变量名称解析为数据集中实际存在的列名，兼容大小写与下划线差异####
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

        if (!is.na(matched_index)) {
            return(data_names[[matched_index]])
        }

        NA_character_
    }, character(1))

    resolved_vars
}

### 根据数据集中可用且有变异的协变量解析模型配置####
resolve_model_specs_for_data <- function(model_specs, data, dataset_name) {
    purrr::imap(model_specs, function(model_spec, model_name) {
        adjust_vars <- model_spec$adjust_vars

        if (is.null(adjust_vars) || length(adjust_vars) == 0) {
            return(model_spec)
        }

        resolved_vars <- resolve_adjust_var_aliases(adjust_vars, data)
        available_vars <- unique(stats::na.omit(resolved_vars))
        missing_vars <- adjust_vars[is.na(resolved_vars)]

        retained_vars <- available_vars[purrr::map_lgl(available_vars, function(var_name) {
            values <- data[[var_name]]
            length(unique(values[!is.na(values)])) >= 2
        })]
        constant_vars <- setdiff(available_vars, retained_vars)

        removed_vars <- c(missing_vars, constant_vars)
        if (length(removed_vars) > 0) {
            message(sprintf(
                "%s %s 移除协变量: %s",
                dataset_name,
                model_name,
                paste(removed_vars, collapse = ", ")
            ))
        }

        model_spec$adjust_vars <- if (length(retained_vars) == 0) NULL else retained_vars
        model_spec
    })
}


## 汇总与绘图函数####

### 根据模型名称返回显著性标记字段####
get_signal_column <- function(model_name) {
    if (model_name == "M0") {
        return("signal_detected")
    }
    "adjusted_signal_detected"
}

### 根据模型名称返回效应值和P值字段####
get_effect_columns <- function(model_name) {
    if (model_name == "M0") {
        return(list(
            effect = "ROR",
            lower = "ROR_lower",
            upper = "ROR_upper",
            p_value = "p_value",
            p_report = "p_adjusted"
        ))
    }

    list(
        effect = "adjusted_OR",
        lower = "adjusted_OR_lower",
        upper = "adjusted_OR_upper",
        p_value = "adjusted_p_value",
        p_report = "adjusted_p_adjusted"
    )
}

### 汇总每个模型的运行规模和显著信号数量####
build_model_overview <- function(model_results, model_specs) {
    purrr::map_dfr(names(model_results), function(model_name) {
        result <- model_results[[model_name]]
        signal_col <- get_signal_column(model_name)
        adjust_vars <- model_specs[[model_name]]$adjust_vars

        tibble(
            model = model_name,
            adjust_vars = if (is.null(adjust_vars)) "未调整" else paste(adjust_vars, collapse = ", "),
            n_processed = result$params$n_processed,
            n_skipped = result$params$n_skipped,
            n_fisher_test = result$params$n_fisher_test,
            n_cancers = result$params$n_cancers,
            n_iraes = result$params$n_iraes,
            signal_count = sum(result$details[[signal_col]], na.rm = TRUE),
            p_adjust_method = result$params$p_adjust_method
        )
    })
}

### 汇总不同分类在各模型中的显著SOC信号数量####
build_signal_summary <- function(model_results, category_levels) {
    signal_summary_long <- purrr::map_dfr(names(model_results), function(model_name) {
        result <- model_results[[model_name]]
        signal_col <- get_signal_column(model_name)

        result$details %>%
            group_by(cancer_type) %>%
            summarise(signal_count = sum(.data[[signal_col]], na.rm = TRUE), .groups = "drop") %>%
            mutate(model = model_name)
    })

    signal_summary_long %>%
        mutate(cancer_type = factor(cancer_type, levels = category_levels)) %>%
        arrange(cancer_type) %>%
        tidyr::pivot_wider(names_from = model, values_from = signal_count)
}

### 汇总不同SOC在各模型中的显著信号数量####
build_soc_summary <- function(model_results, soc_levels) {
    soc_summary_long <- purrr::map_dfr(names(model_results), function(model_name) {
        result <- model_results[[model_name]]
        signal_col <- get_signal_column(model_name)

        result$details %>%
            mutate(irae = get_detail_irae_display(result$details)) %>%
            group_by(irae) %>%
            summarise(signal_count = sum(.data[[signal_col]], na.rm = TRUE), .groups = "drop") %>%
            mutate(model = model_name)
    })

    soc_summary_long %>%
        mutate(irae = factor(irae, levels = abbreviate_soc_vector(soc_levels))) %>%
        arrange(irae) %>%
        tidyr::pivot_wider(names_from = model, values_from = signal_count)
}

### 提取指定模型中效应值最高的显著信号####
build_top_signals <- function(result, model_name, top_n = 20) {
    signal_col <- get_signal_column(model_name)
    effect_cols <- get_effect_columns(model_name)
    irae_display <- get_detail_irae_display(result$details)
    irae_full <- get_detail_irae_full(result$details)

    result$details %>%
        mutate(
            irae_display = irae_display,
            irae_full_label = irae_full
        ) %>%
        filter(.data[[signal_col]]) %>%
        arrange(desc(.data[[effect_cols$effect]])) %>%
        transmute(
            cancer_type = cancer_type,
            irae_full = irae_full_label,
            irae = irae_display,
            effect = .data[[effect_cols$effect]],
            lower = .data[[effect_cols$lower]],
            upper = .data[[effect_cols$upper]],
            raw_p_value = .data[[effect_cols$p_value]],
            reported_p_value = .data[[effect_cols$p_report]]
        ) %>%
        slice_head(n = top_n)
}

# 分析执行####
source("Code/disproportionality_analysis.R")
### 运行单个数据库的primary_organ不成比例分析并输出结果####
run_primary_organ_analysis_for_dataset <- function(dataset_name, dataset_config) {
    message(sprintf("\n========== %s primary_organ不成比例分析开始 ==========", dataset_name))

    load(file.path(output_path, dataset_config$file_name))

    if (!exists(dataset_config$object_name, inherits = FALSE)) {
        stop(sprintf("%s 数据对象不存在: %s", dataset_name, dataset_config$object_name))
    }

    malignant_final <- get(dataset_config$object_name, inherits = FALSE)

    excluded_organs <- intersect(
        c("CUP", "Hematologic", "Unmapped", "Unknown_Primary"),
        unique(malignant_final$primary_organ)
    )

    primary_organ_analysis <- malignant_final %>%
        filter(has_ici == TRUE) %>%
        filter(!primary_organ %in% excluded_organs)

    organ_counts <- primary_organ_analysis %>%
        count(primary_organ, name = "n") %>%
        arrange(desc(n))


    soc_vars <- names(primary_organ_analysis)[grep("^SOC_", names(primary_organ_analysis))]

    message(paste(dataset_name, "排除的primary_organ:", paste(excluded_organs, collapse = ", ")))
    message(paste(dataset_name, "primary_organ分析样本量（仅ICI）:", nrow(primary_organ_analysis)))
    message(paste(dataset_name, "保留器官数:", nrow(organ_counts)))
    message(paste(dataset_name, "SOC变量数:", length(soc_vars)))
    print(organ_counts)

    model_specs_current <- resolve_model_specs_for_data(
        model_specs = model_specs,
        data = primary_organ_analysis,
        dataset_name = dataset_name
    )

    # primary_organ主分析：器官间irAE差异（器官 vs 其他器官，仅ICI患者）
    for (model_name in names(model_specs_current)) {
        adjust_vars <- model_specs_current[[model_name]]$adjust_vars

        output_file <- file.path(
            script_results_path,
            paste0(dataset_name, "_primary_organ_SOC_", model_name, ".xlsx")
        )

        if (file.exists(output_file) && model_result_excel_is_compatible(output_file)) {
            message(paste(dataset_name, model_name, "结果已存在，跳过"))
            next
        }

        message(paste(dataset_name, "运行", model_name, "..."))

        result <- tryCatch(
            {
                run_disproportionality_analysis(
                    data = primary_organ_analysis,
                    cancer_var = "primary_organ",
                    irae_level = "SOC",
                    min_cases = 3,
                    expected_count_threshold = 5,
                    filter_ici = FALSE,
                    exclude_cancer = NULL,
                    adjust_vars = adjust_vars,
                    p_adjust_method = "none"
                )
            },
            error = function(e) {
                message(paste(dataset_name, model_name, "分析失败:", conditionMessage(e)))
                NULL
            }
        )

        if (!is.null(result)) {
            ensure_parent_dir(output_file)
            export_to_excel(result, output_file)
            message(paste(dataset_name, model_name, "结果已保存:", output_file))
        }
    }

    message(paste(dataset_name, "primary_organ不成比例分析完成"))

    model_results <- load_model_results(
        model_specs = model_specs_current,
        results_path = script_results_path,
        prefix = paste0(dataset_name, "_primary_organ_SOC_")
    )

    soc_levels <- model_results$M0$details %>%
        group_by(irae) %>%
        summarise(signal_count = sum(signal_detected, na.rm = TRUE), .groups = "drop") %>%
        arrange(desc(signal_count), irae) %>%
        pull(irae)

    last_model_analysis1 <- tail(names(model_specs_current), 1)
    model_overview <- build_model_overview(model_results, model_specs_current)
    organ_signal_summary <- build_signal_summary(model_results, organ_counts$primary_organ)
    soc_signal_summary <- build_soc_summary(model_results, soc_levels)
    top_signals_m0 <- build_top_signals(model_results$M0, "M0", top_n = 20)
    top_signals_last <- build_top_signals(model_results[[last_model_analysis1]], last_model_analysis1, top_n = 20)

    summary_output_file <- file.path(
        script_results_path,
        paste0(dataset_name, "_primary_organ_SOC_summary.xlsx")
    )

    summary_sheets <- list(
        "模型概览" = model_overview,
        "器官样本量" = organ_counts,
        "器官信号数" = organ_signal_summary,
        "SOC信号数" = soc_signal_summary,
        "M0_top_signals" = top_signals_m0
    )
    summary_sheets[[paste0(last_model_analysis1, "_top_signals")]] <- top_signals_last

    ensure_parent_dir(summary_output_file)
    writexl::write_xlsx(summary_sheets, path = summary_output_file)
    message(paste(dataset_name, "总览工作簿已保存:", summary_output_file))


    gc()
    message(sprintf("\n========== %s primary_organ不成比例分析完成 ==========", dataset_name))
}

purrr::iwalk(dataset_configs, function(dataset_config, dataset_name) {
    run_primary_organ_analysis_for_dataset(dataset_name, dataset_config)
})

message("\n========== primary_organ不成比例分析全部完成 ==========")
message("热图绘制请运行 Code/3.1.1热图绘制.R")
