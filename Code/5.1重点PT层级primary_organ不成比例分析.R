# 重点PT层级primary_organ不成比例分析（FAERS + JADER）####

## 1. 加载配置与函数####
source("Code/config.R")
source("Code/mapping_rules.R")
source("Code/disproportionality_analysis.R")

## 2. 并行线程设置####
Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1",
    NUMEXPR_NUM_THREADS = "1"
)

## 3. 结果目录与模型参数####
script_name <- "5.1重点PT层级primary_organ不成比例分析"
script_results_path <- file.path(results_path, script_name)

if (!dir.exists(script_results_path)) {
    dir.create(script_results_path, recursive = TRUE, showWarnings = FALSE)
}

message("========== 重点PT层级primary_organ不成比例分析开始 ==========")
message(paste("结果目录:", script_results_path))

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

model_specs <- list(
    "M0" = list(adjust_vars = NULL),
    "M1" = list(adjust_vars = c("treatment_strategy")),
    "M2" = list(adjust_vars = c("treatment_strategy", "Age", "Sex")),
    "M3" = list(adjust_vars = c("treatment_strategy", "Age", "Sex", "Chemotherapy", "Targeted_therapy"))
)
available_parallel_workers <- max(1L, future::availableCores())
dataset_parallel_workers <- max(1L, min(length(dataset_configs), available_parallel_workers))
model_parallel_workers <- max(
    1L,
    min(
        length(model_specs),
        floor(available_parallel_workers / dataset_parallel_workers)
    )
)


## 4. 辅助函数####

## 函数：ensure_dir_exists ####
# 功能：确保目录存在，必要时递归创建并校验创建结果####
# 输入：dir_path，字符型目录路径####
# 输出：不可见返回已存在的目录路径####
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

## 函数：ensure_parent_dir ####
# 功能：确保输出文件的父目录存在####
# 输入：file_path，字符型文件路径####
# 输出：不可见返回父目录路径####
ensure_parent_dir <- function(file_path) {
    ensure_dir_exists(dirname(file_path))
}

## 函数：load_named_rdata ####
# 功能：从RData文件中读取指定对象且避免污染全局环境####
# 输入：input_path，RData路径；object_name，对象名称####
# 输出：目标R对象####
load_named_rdata <- function(input_path, object_name) {
    data_env <- new.env(parent = emptyenv())
    load(input_path, envir = data_env)

    if (!exists(object_name, envir = data_env, inherits = FALSE)) {
        stop(paste("RData中未找到对象:", object_name))
    }

    get(object_name, envir = data_env, inherits = FALSE)
}

## 函数：build_target_pt_table ####
# 功能：构建固定纳入的12个重点PT清单及其所属SOC信息####
# 输入：无####
# 输出：包含PT展示名、变量名、所属SOC及固定顺序的数据框####
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
    ) %>%
        dplyr::mutate(
            irae_raw = sub("^PT_", "", irae_var),
            irae_full = irae_display
        )
}

## 函数：subset_to_target_pt_vars ####
# 功能：保留全部非PT列以及固定12个目标PT列，缺失目标PT时自动补零列####
# 输入：data，原始分析数据框；target_pt_table，目标PT清单；dataset_name，数据库名称####
# 输出：仅保留目标PT列的数据框####
subset_to_target_pt_vars <- function(data, target_pt_table, dataset_name) {
    all_pt_vars <- names(data)[grepl("^PT_", names(data))]
    target_pt_vars <- target_pt_table$irae_var
    missing_pt_vars <- setdiff(target_pt_vars, names(data))

    if (length(missing_pt_vars) > 0) {
        message(
            dataset_name,
            " 缺少以下目标PT变量，按0事件补零列后继续分析: ",
            paste(missing_pt_vars, collapse = ", ")
        )

        for (current_pt_var in missing_pt_vars) {
            data[[current_pt_var]] <- 0L
        }
    }

    keep_columns <- c(
        setdiff(names(data), all_pt_vars),
        target_pt_vars
    )

    data[, unique(keep_columns), drop = FALSE]
}

## 函数：enhance_pt_result_for_export ####
# 功能：为PT层分析结果补充展示字段并按固定顺序重命名矩阵列####
# 输入：result，run_disproportionality_analysis返回对象；target_pt_table，目标PT清单####
# 输出：补充展示字段且矩阵列名已重命名的结果对象####
enhance_pt_result_for_export <- function(result, target_pt_table) {
    pt_lookup <- target_pt_table %>%
        dplyr::select(
            irae_var,
            irae_raw,
            irae_display,
            irae_full,
            parent_soc,
            parent_soc_abbrev,
            pt_order
        )

    result$details <- result$details %>%
        dplyr::left_join(pt_lookup, by = "irae_var") %>%
        dplyr::mutate(
            irae_display = dplyr::coalesce(as.character(irae_display), as.character(irae)),
            irae_full = dplyr::coalesce(as.character(irae_full), as.character(irae_display)),
            parent_soc = as.character(parent_soc),
            parent_soc_abbrev = as.character(parent_soc_abbrev),
            pt_order = suppressWarnings(as.numeric(pt_order))
        )

    if (any(is.na(result$details$parent_soc_abbrev))) {
        missing_pt_labels <- result$details %>%
            dplyr::filter(is.na(parent_soc_abbrev)) %>%
            dplyr::distinct(irae_var) %>%
            dplyr::pull(irae_var)

        stop(
            paste0(
                "以下PT结果未能映射回固定清单: ",
                paste(missing_pt_labels, collapse = ", ")
            )
        )
    }

    ordered_front_columns <- c(
        "cancer_type", "irae_display", "irae_full", "parent_soc", "parent_soc_abbrev",
        "pt_order", "irae", "irae_var"
    )

    result$details <- result$details %>%
        dplyr::arrange(cancer_type, pt_order, irae_display) %>%
        dplyr::select(dplyr::any_of(ordered_front_columns), dplyr::everything())

    matrix_column_lookup <- stats::setNames(
        target_pt_table$irae_display,
        target_pt_table$irae_raw
    )
    ordered_matrix_columns <- target_pt_table$irae_raw

    relabel_matrix_columns <- function(current_matrix) {
        if (is.null(current_matrix)) {
            return(NULL)
        }

        missing_columns <- setdiff(ordered_matrix_columns, colnames(current_matrix))
        if (length(missing_columns) > 0) {
            stop(
                paste0(
                    "结果矩阵缺少以下目标PT列: ",
                    paste(missing_columns, collapse = ", ")
                )
            )
        }

        current_matrix <- current_matrix[, ordered_matrix_columns, drop = FALSE]
        colnames(current_matrix) <- unname(matrix_column_lookup[colnames(current_matrix)])
        current_matrix
    }

    result$ror_matrix <- relabel_matrix_columns(result$ror_matrix)
    result$signal_matrix <- relabel_matrix_columns(result$signal_matrix)
    result$adjusted_or_matrix <- relabel_matrix_columns(result$adjusted_or_matrix)
    result$adjusted_signal_matrix <- relabel_matrix_columns(result$adjusted_signal_matrix)

    result
}

## 函数：read_param_value ####
# 功能：从分析参数工作表中读取指定参数值，缺失时返回默认值####
# 输入：params_df，参数数据框；param_name，参数名；default，默认值####
# 输出：字符型参数值或默认值####
read_param_value <- function(params_df, param_name, default = NA_character_) {
    matched_value <- params_df %>%
        dplyr::filter(参数 == param_name) %>%
        dplyr::pull(值)

    if (length(matched_value) == 0) {
        return(default)
    }

    matched_value[[1]]
}

## 函数：parse_numeric_param ####
# 功能：将参数文本转换为数值，无法转换时返回默认值####
# 输入：value，字符型参数值；default，默认数值####
# 输出：数值型参数####
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

## 函数：parse_logical_param ####
# 功能：将参数文本转换为逻辑值，无法转换时返回默认值####
# 输入：value，字符型参数值；default，默认逻辑值####
# 输出：逻辑型参数####
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

## 函数：parse_character_vector_param ####
# 功能：将逗号分隔的参数文本还原为字符向量####
# 输入：value，字符型参数值####
# 输出：字符向量或NULL####
parse_character_vector_param <- function(value) {
    if (length(value) == 0 || is.na(value) || value == "") {
        return(NULL)
    }

    trimws(strsplit(as.character(value), ",", fixed = TRUE)[[1]])
}

## 函数：load_result_from_excel ####
# 功能：从单个模型结果Excel读取明细表和参数表，供后续汇总使用####
# 输入：input_file，模型结果Excel路径####
# 输出：包含details与params的列表####
load_result_from_excel <- function(input_file) {
    details <- readxl::read_xlsx(input_file, sheet = "详细结果")
    params_df <- readxl::read_xlsx(input_file, sheet = "分析参数")

    character_cols <- intersect(
        c(
            "cancer_type", "irae", "irae_var", "irae_display", "irae_full",
            "parent_soc", "parent_soc_abbrev", "analysis_method",
            "p_value_method", "result_note"
        ),
        names(details)
    )
    logical_cols <- intersect(
        c(
            "signal_detected", "adjusted_signal_detected",
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

    params <- list(
        cancer_var = read_param_value(params_df, "cancer_var"),
        irae_level = read_param_value(params_df, "irae_level"),
        n_cancers = parse_numeric_param(read_param_value(params_df, "n_cancers")),
        n_iraes = parse_numeric_param(read_param_value(params_df, "n_iraes")),
        filter_ici = parse_logical_param(read_param_value(params_df, "filter_ici"), default = FALSE),
        exclude_cancer = parse_character_vector_param(read_param_value(params_df, "exclude_cancer", default = "")),
        adjust_vars = parse_character_vector_param(read_param_value(params_df, "adjust_vars", default = "")),
        min_cases = parse_numeric_param(read_param_value(params_df, "min_cases")),
        expected_count_threshold = parse_numeric_param(
            read_param_value(params_df, "expected_count_threshold", default = ""),
            default = NA_real_
        ),
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

## 函数：model_result_excel_is_compatible ####
# 功能：检查模型结果Excel是否已包含当前5.1分析所需字段####
# 输入：input_file，模型结果Excel路径####
# 输出：逻辑值，TRUE表示结果文件可直接复用####
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
                "fisher_test_used", "min_expected_count", "irae_display",
                "irae_full", "parent_soc", "parent_soc_abbrev", "pt_order"
            )
            required_param_names <- c(
                "expected_count_threshold", "n_fisher_test",
                "n_processed", "n_skipped"
            )

            all(required_detail_columns %in% detail_columns) &&
                all(required_param_names %in% param_names)
        },
        error = function(e) {
            FALSE
        }
    )
}

## 函数：get_model_result_paths ####
# 功能：构建各模型结果工作簿路径####
# 输入：model_specs，模型配置列表；results_path，结果目录；prefix，文件名前缀####
# 输出：按模型名命名的文件路径列表####
get_model_result_paths <- function(model_specs, results_path, prefix) {
    stats::setNames(
        file.path(
            results_path,
            paste0(prefix, names(model_specs), ".xlsx")
        ),
        names(model_specs)
    )
}

## 函数：load_model_results ####
# 功能：从Excel工作簿批量载入各模型结果，供汇总使用####
# 输入：model_specs，模型配置列表；results_path，结果目录；prefix，文件名前缀####
# 输出：按模型名命名的结果列表####
load_model_results <- function(model_specs, results_path, prefix) {
    model_result_paths <- get_model_result_paths(model_specs, results_path, prefix)
    all_model_files_exist <- all(file.exists(unname(model_result_paths)))
    all_model_files_compatible <- all_model_files_exist &&
        all(purrr::map_lgl(model_result_paths, model_result_excel_is_compatible))

    if (!all_model_files_compatible) {
        stop("模型结果Excel工作簿缺失或与当前字段要求不兼容，请先准备好5.1分析结果")
    }

    purrr::map(model_result_paths, load_result_from_excel)
}

## 函数：resolve_adjust_var_aliases ####
# 功能：将协变量名称解析为数据集中实际存在的列名，兼容大小写与下划线差异####
# 输入：adjust_vars，协变量名称向量；data，分析数据框####
# 输出：与adjust_vars等长的解析结果向量####
resolve_adjust_var_aliases <- function(adjust_vars, data) {
    if (is.null(adjust_vars) || length(adjust_vars) == 0) {
        return(character(0))
    }

    data_names <- names(data)
    normalized_data_names <- tolower(gsub("[^[:alnum:]]+", "", data_names))

    vapply(adjust_vars, function(var_name) {
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
}

## 函数：resolve_model_specs_for_data ####
# 功能：根据数据集中可用且有变异的协变量解析模型配置####
# 输入：model_specs，模型配置列表；data，分析数据框；dataset_name，数据库名称####
# 输出：剔除缺失或无变异协变量后的模型配置列表####
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

## 函数：run_model_tasks_in_parallel ####
# 功能：并行运行同一数据库下待执行的全部模型分析，并在结束后恢复future配置####
# 输入：data，分析数据框；model_specs_current，当前模型配置；model_names_to_run，待运行模型名向量；dataset_name，数据库名称####
# 输出：按模型名命名的分析结果列表####
run_model_tasks_in_parallel <- function(data, model_specs_current, model_names_to_run, dataset_name) {
    if (length(model_names_to_run) == 0) {
        return(list())
    }

    worker_count <- min(
        length(model_names_to_run),
        max(1L, as.integer(model_parallel_workers))
    )
    backend <- if (future::supportsMulticore()) future::multicore else future::multisession

    old_plan <- future::plan()
    old_future_globals_max_size <- getOption("future.globals.maxSize")

    on.exit(future::plan(old_plan), add = TRUE)
    on.exit(options(future.globals.maxSize = old_future_globals_max_size), add = TRUE)

    options(future.globals.maxSize = Inf)

    if (worker_count > 1L) {
        future::plan(backend, workers = worker_count)
        message(dataset_name, " 模型并行运行，worker数: ", worker_count)
    } else {
        future::plan(future::sequential)
        message(dataset_name, " 模型顺序运行，worker数: 1")
    }

    model_results <- future.apply::future_lapply(
        model_names_to_run,
        function(model_name) {
            adjust_vars <- model_specs_current[[model_name]]$adjust_vars

            tryCatch(
                {
                    run_disproportionality_analysis(
                        data = data,
                        cancer_var = "primary_organ",
                        irae_level = "PT",
                        min_cases = 3,
                        expected_count_threshold = 5,
                        filter_ici = FALSE,
                        exclude_cancer = NULL,
                        adjust_vars = adjust_vars,
                        p_adjust_method = "none"
                    )
                },
                error = function(e) {
                    structure(
                        list(
                            model_name = model_name,
                            error_message = conditionMessage(e)
                        ),
                        class = "pt_model_error"
                    )
                }
            )
        },
        future.seed = TRUE
    )

    stats::setNames(model_results, model_names_to_run)
}

## 函数：run_dataset_tasks_in_parallel ####
# 功能：并行运行双数据库任务，并在结束后恢复future配置####
# 输入：dataset_configs，数据库配置列表；target_pt_table，目标PT清单####
# 输出：无返回值，副作用为生成各数据库结果文件####
run_dataset_tasks_in_parallel <- function(dataset_configs, target_pt_table) {
    dataset_names <- names(dataset_configs)

    if (length(dataset_names) == 0) {
        return(invisible(NULL))
    }

    worker_count <- min(
        length(dataset_names),
        max(1L, as.integer(dataset_parallel_workers))
    )

    old_plan <- future::plan()
    old_future_globals_max_size <- getOption("future.globals.maxSize")

    on.exit(future::plan(old_plan), add = TRUE)
    on.exit(options(future.globals.maxSize = old_future_globals_max_size), add = TRUE)

    options(future.globals.maxSize = Inf)

    if (worker_count > 1L) {
        future::plan(future::multisession, workers = worker_count)
        message("数据库层并行运行，worker数: ", worker_count)
    } else {
        future::plan(future::sequential)
        message("数据库层顺序运行，worker数: 1")
    }

    future.apply::future_lapply(
        dataset_names,
        function(dataset_name) {
            run_primary_organ_pt_analysis_for_dataset(
                dataset_name = dataset_name,
                dataset_config = dataset_configs[[dataset_name]],
                target_pt_table = target_pt_table
            )
        },
        future.seed = TRUE
    )

    invisible(NULL)
}

## 函数：get_signal_column ####
# 功能：根据模型名称返回显著性标记字段名####
# 输入：model_name，模型名称####
# 输出：字符型显著性字段名####
get_signal_column <- function(model_name) {
    if (model_name == "M0") {
        return("signal_detected")
    }

    "adjusted_signal_detected"
}

## 函数：get_effect_columns ####
# 功能：根据模型名称返回效应值和P值字段名####
# 输入：model_name，模型名称####
# 输出：包含效应值与P值字段名的列表####
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

## 函数：get_detail_irae_display ####
# 功能：从结果明细表中获取PT展示标签####
# 输入：details，结果明细表####
# 输出：字符型PT展示标签向量####
get_detail_irae_display <- function(details) {
    if ("irae_display" %in% names(details)) {
        return(as.character(details$irae_display))
    }

    as.character(details$irae)
}

## 函数：get_detail_parent_soc_abbrev ####
# 功能：从结果明细表中获取所属SOC缩写####
# 输入：details，结果明细表####
# 输出：字符型SOC缩写向量####
get_detail_parent_soc_abbrev <- function(details) {
    if ("parent_soc_abbrev" %in% names(details)) {
        return(as.character(details$parent_soc_abbrev))
    }

    rep(NA_character_, nrow(details))
}

## 函数：build_model_overview ####
# 功能：汇总各模型运行规模和显著信号数量####
# 输入：model_results，模型结果列表；model_specs，模型配置列表####
# 输出：模型概览数据框####
build_model_overview <- function(model_results, model_specs) {
    purrr::map_dfr(names(model_results), function(model_name) {
        result <- model_results[[model_name]]
        signal_col <- get_signal_column(model_name)
        adjust_vars <- model_specs[[model_name]]$adjust_vars

        tibble::tibble(
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

## 函数：build_signal_summary ####
# 功能：汇总不同器官在各模型中的显著PT信号数量####
# 输入：model_results，模型结果列表；category_levels，器官顺序向量####
# 输出：器官信号数宽表####
build_signal_summary <- function(model_results, category_levels) {
    signal_summary_long <- purrr::map_dfr(names(model_results), function(model_name) {
        result <- model_results[[model_name]]
        signal_col <- get_signal_column(model_name)

        result$details %>%
            dplyr::group_by(cancer_type) %>%
            dplyr::summarise(signal_count = sum(.data[[signal_col]], na.rm = TRUE), .groups = "drop") %>%
            dplyr::mutate(model = model_name)
    })

    signal_summary_long %>%
        dplyr::mutate(cancer_type = factor(cancer_type, levels = category_levels)) %>%
        dplyr::arrange(cancer_type) %>%
        tidyr::pivot_wider(names_from = model, values_from = signal_count)
}

## 函数：build_pt_summary ####
# 功能：汇总不同重点PT在各模型中的显著器官数量####
# 输入：model_results，模型结果列表；target_pt_table，目标PT清单####
# 输出：PT信号数宽表####
build_pt_summary <- function(model_results, target_pt_table) {
    pt_summary_long <- purrr::map_dfr(names(model_results), function(model_name) {
        result <- model_results[[model_name]]
        signal_col <- get_signal_column(model_name)

        result$details %>%
            dplyr::mutate(
                irae_display = get_detail_irae_display(result$details),
                parent_soc_abbrev = get_detail_parent_soc_abbrev(result$details),
                pt_order = suppressWarnings(as.numeric(pt_order))
            ) %>%
            dplyr::group_by(irae_display, parent_soc_abbrev, pt_order) %>%
            dplyr::summarise(signal_count = sum(.data[[signal_col]], na.rm = TRUE), .groups = "drop") %>%
            dplyr::mutate(model = model_name)
    })

    target_pt_table %>%
        dplyr::select(pt_order, irae_display, parent_soc, parent_soc_abbrev) %>%
        dplyr::left_join(
            pt_summary_long %>%
                tidyr::pivot_wider(names_from = model, values_from = signal_count),
            by = c("pt_order", "irae_display", "parent_soc_abbrev")
        ) %>%
        dplyr::arrange(pt_order)
}

## 函数：build_top_signals ####
# 功能：提取指定模型中效应值最高的显著PT信号####
# 输入：result，单个模型结果；model_name，模型名称；top_n，返回条目数####
# 输出：Top信号数据框####
build_top_signals <- function(result, model_name, top_n = 20) {
    signal_col <- get_signal_column(model_name)
    effect_cols <- get_effect_columns(model_name)

    result$details %>%
        dplyr::mutate(
            irae_display = get_detail_irae_display(result$details),
            parent_soc_abbrev = get_detail_parent_soc_abbrev(result$details)
        ) %>%
        dplyr::filter(.data[[signal_col]]) %>%
        dplyr::arrange(desc(.data[[effect_cols$effect]]), pt_order, cancer_type) %>%
        dplyr::transmute(
            cancer_type = cancer_type,
            parent_soc = parent_soc_abbrev,
            irae_full = irae_display,
            irae = irae_display,
            effect = .data[[effect_cols$effect]],
            lower = .data[[effect_cols$lower]],
            upper = .data[[effect_cols$upper]],
            raw_p_value = .data[[effect_cols$p_value]],
            reported_p_value = .data[[effect_cols$p_report]]
        ) %>%
        dplyr::slice_head(n = top_n)
}

## 函数：run_primary_organ_pt_analysis_for_dataset ####
# 功能：运行单个数据库的重点PT层级primary_organ不成比例分析并导出结果####
# 输入：dataset_name，数据库名称；dataset_config，数据对象配置；target_pt_table，目标PT清单####
# 输出：无返回值，副作用为生成模型结果工作簿与汇总工作簿####
run_primary_organ_pt_analysis_for_dataset <- function(dataset_name, dataset_config, target_pt_table) {
    message(sprintf("\n========== %s 重点PT层级primary_organ不成比例分析开始 ==========", dataset_name))

    malignant_final <- load_named_rdata(
        input_path = file.path(output_path, dataset_config$file_name),
        object_name = dataset_config$object_name
    )

    excluded_organs <- intersect(
        c("CUP", "Hematologic", "Unmapped", "Unknown_Primary"),
        unique(malignant_final$primary_organ)
    )

    primary_organ_analysis <- malignant_final %>%
        dplyr::filter(has_ici == TRUE) %>%
        dplyr::filter(!primary_organ %in% excluded_organs)

    primary_organ_analysis <- subset_to_target_pt_vars(
        data = primary_organ_analysis,
        target_pt_table = target_pt_table,
        dataset_name = dataset_name
    )

    organ_counts <- primary_organ_analysis %>%
        dplyr::count(primary_organ, name = "n") %>%
        dplyr::arrange(desc(n))

    message(paste(dataset_name, "排除的primary_organ:", paste(excluded_organs, collapse = ", ")))
    message(paste(dataset_name, "重点PT层分析样本量（仅ICI）:", nrow(primary_organ_analysis)))
    message(paste(dataset_name, "保留器官数:", nrow(organ_counts)))
    message(paste(dataset_name, "固定纳入PT数:", nrow(target_pt_table)))
    print(target_pt_table %>% dplyr::select(pt_order, irae_display, parent_soc_abbrev))

    model_specs_current <- resolve_model_specs_for_data(
        model_specs = model_specs,
        data = primary_organ_analysis,
        dataset_name = dataset_name
    )

    model_output_paths <- stats::setNames(
        file.path(
            script_results_path,
            paste0(dataset_name, "_primary_organ_PT_", names(model_specs_current), ".xlsx")
        ),
        names(model_specs_current)
    )
    model_names_to_run <- names(model_output_paths)[!vapply(model_output_paths, model_result_excel_is_compatible, logical(1))]

    for (model_name in setdiff(names(model_specs_current), model_names_to_run)) {
        message(paste(dataset_name, model_name, "结果已存在，跳过"))
    }

    if (length(model_names_to_run) > 0) {
        message(dataset_name, " 待运行模型: ", paste(model_names_to_run, collapse = ", "))
    }

    parallel_model_results <- run_model_tasks_in_parallel(
        data = primary_organ_analysis,
        model_specs_current = model_specs_current,
        model_names_to_run = model_names_to_run,
        dataset_name = dataset_name
    )

    for (model_name in names(parallel_model_results)) {
        current_result <- parallel_model_results[[model_name]]
        output_file <- model_output_paths[[model_name]]

        if (inherits(current_result, "pt_model_error")) {
            message(paste(dataset_name, model_name, "分析失败:", current_result$error_message))
            next
        }

        current_result <- enhance_pt_result_for_export(
            result = current_result,
            target_pt_table = target_pt_table
        )
        ensure_parent_dir(output_file)
        export_to_excel(current_result, output_file)
        message(paste(dataset_name, model_name, "结果已保存:", output_file))
    }

    model_results <- load_model_results(
        model_specs = model_specs_current,
        results_path = script_results_path,
        prefix = paste0(dataset_name, "_primary_organ_PT_")
    )

    model_overview <- build_model_overview(model_results, model_specs_current)
    organ_signal_summary <- build_signal_summary(model_results, organ_counts$primary_organ)
    pt_signal_summary <- build_pt_summary(model_results, target_pt_table)
    top_signals_m0 <- build_top_signals(model_results$M0, "M0", top_n = 20)
    top_signals_m3 <- build_top_signals(model_results$M3, "M3", top_n = 20)

    summary_output_file <- file.path(
        script_results_path,
        paste0(dataset_name, "_primary_organ_PT_summary.xlsx")
    )

    summary_sheets <- list(
        "模型概览" = model_overview,
        "器官样本量" = organ_counts,
        "器官信号数" = organ_signal_summary,
        "PT信号数" = pt_signal_summary,
        "纳入PT清单" = target_pt_table %>%
            dplyr::select(pt_order, irae_display, irae_var, parent_soc, parent_soc_abbrev),
        "M0_top_signals" = top_signals_m0,
        "M3_top_signals" = top_signals_m3
    )

    ensure_parent_dir(summary_output_file)
    writexl::write_xlsx(summary_sheets, path = summary_output_file)
    message(paste(dataset_name, "总览工作簿已保存:", summary_output_file))

    gc()
    message(sprintf("\n========== %s 重点PT层级primary_organ不成比例分析完成 ==========", dataset_name))
}


## 5. 构建固定PT清单####
target_pt_table <- build_target_pt_table()
message("可用并行worker总数: ", available_parallel_workers)
message("数据库层worker数: ", dataset_parallel_workers)
message("每个数据库内模型层worker数: ", model_parallel_workers)


## 6. 运行双库重点PT层分析####
run_dataset_tasks_in_parallel(
    dataset_configs = dataset_configs,
    target_pt_table = target_pt_table
)

message("\n========== 重点PT层级primary_organ不成比例分析全部完成 ==========")
message("热图绘制请运行 Code/5.1.1热图绘制.R")
