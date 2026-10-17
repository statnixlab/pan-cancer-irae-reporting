# primary_organ 热图绘制（FAERS + JADER）####

# 1. 加载配置 ####
source("Code/config.R")
source("Code/mapping_rules.R")

# 2. 结果目录与模型配置 ####
script_name <- "3.1primary_organ不成比例分析"
script_results_path <- file.path(results_path, script_name)
jader_heatmap_min_cases <- 100

model_specs <- list(
    "M0" = list(adjust_vars = NULL),
    "M1" = list(adjust_vars = c("treatment_strategy")),
    "M2" = list(adjust_vars = c("treatment_strategy", "Age", "Sex")),
    "M3" = list(adjust_vars = c("treatment_strategy", "Age", "Sex", "Chemotherapy", "Targeted_therapy"))
)

message("========== primary_organ 热图绘制开始 ==========")
message(paste("结果目录:", script_results_path))

# 3. 结果读取辅助函数 ####

# 函数：get_soc_abbrev_safe ####
## 功能：将SOC名称转换为缩写，无法映射时返回原始名称 ####
## 输入：irae_name，字符型SOC名称 ####
## 输出：字符型SOC缩写或原始名称 ####
get_soc_abbrev_safe <- function(irae_name) {
    soc_full <- gsub("_", " ", irae_name)
    soc_abbrev <- get_soc_abbrev(soc_full)

    if (is.na(soc_abbrev)) {
        return(irae_name)
    }

    soc_abbrev
}

# 函数：abbreviate_soc_vector ####
## 功能：批量将SOC名称向量转换为缩写 ####
## 输入：irae_vec，字符型SOC名称向量 ####
## 输出：字符型SOC缩写向量 ####
abbreviate_soc_vector <- function(irae_vec) {
    vapply(irae_vec, get_soc_abbrev_safe, character(1))
}

# 函数：get_detail_irae_display ####
## 功能：为结果明细表提供绘图使用的SOC展示标签 ####
## 输入：details，不成比例分析结果明细表 ####
## 输出：字符型SOC展示标签向量 ####
get_detail_irae_display <- function(details) {
    if ("irae_abbrev" %in% names(details)) {
        return(details$irae_abbrev)
    }

    abbreviate_soc_vector(details$irae)
}

# 函数：get_model_result_paths ####
## 功能：构建各模型结果Excel文件路径 ####
## 输入：model_specs，模型配置列表；results_path，结果目录；prefix，文件名前缀 ####
## 输出：按模型名称命名的文件路径列表 ####
get_model_result_paths <- function(model_specs, results_path, prefix = "FAERS_primary_organ_SOC_") {
    stats::setNames(
        file.path(
            results_path,
            paste0(prefix, names(model_specs), ".xlsx")
        ),
        names(model_specs)
    )
}

# 函数：load_result_from_excel ####
## 功能：从单个模型结果Excel读取热图所需明细表 ####
## 输入：input_file，模型结果Excel路径 ####
## 输出：包含details数据框的列表 ####
load_result_from_excel <- function(input_file) {
    details <- readxl::read_xlsx(input_file, sheet = "详细结果")

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

    list(
        details = tibble::as_tibble(details)
    )
}

# 函数：model_result_excel_is_compatible ####
## 功能：检查模型结果Excel是否包含热图绘制所需字段 ####
## 输入：input_file，模型结果Excel路径 ####
## 输出：逻辑值，TRUE表示可直接用于热图绘制 ####
model_result_excel_is_compatible <- function(input_file) {
    if (!file.exists(input_file)) {
        return(FALSE)
    }

    tryCatch(
        {
            detail_columns <- names(readxl::read_xlsx(input_file, sheet = "详细结果", n_max = 0))
            required_detail_columns <- c("cancer_type", "irae", "ROR", "signal_detected")
            all(required_detail_columns %in% detail_columns)
        },
        error = function(e) {
            FALSE
        }
    )
}

# 函数：load_model_results ####
## 功能：批量读取单个数据库全部模型结果 ####
## 输入：model_specs，模型配置列表；results_path，结果目录；prefix，文件名前缀 ####
## 输出：按模型名称命名的结果列表 ####
load_model_results <- function(model_specs, results_path, prefix = "FAERS_primary_organ_SOC_") {
    model_result_paths <- get_model_result_paths(model_specs, results_path, prefix)
    all_model_files_exist <- all(file.exists(unname(model_result_paths)))
    all_model_files_compatible <- all_model_files_exist &&
        all(purrr::map_lgl(model_result_paths, model_result_excel_is_compatible))

    if (!all_model_files_compatible) {
        stop("模型结果Excel工作簿缺失或与当前热图字段要求不兼容，请先运行3.1primary_organ不成比例分析.R")
    }

    model_results_loaded <- purrr::map(model_result_paths, load_result_from_excel)
    message("模型结果已从Excel载入，热图仅使用已保存结果文件")

    model_results_loaded
}

# 函数：load_heatmap_organ_counts ####
## 功能：读取热图汇总表中的器官样本量工作表并识别计数列 ####
## 输入：summary_file，summary工作簿路径 ####
## 输出：包含器官样本量数据和计数列名称的列表 ####
load_heatmap_organ_counts <- function(summary_file) {
    if (!file.exists(summary_file)) {
        stop(paste("热图汇总文件不存在:", summary_file))
    }

    organ_counts <- readxl::read_xlsx(summary_file, sheet = "器官样本量") %>%
        dplyr::mutate(primary_organ = as.character(primary_organ))

    count_column_candidates <- c("n", "total", "ici_n", "non_ici_n")
    matched_count_columns <- intersect(count_column_candidates, names(organ_counts))

    if (length(matched_count_columns) > 0) {
        count_column <- matched_count_columns[[1]]
    } else {
        fallback_columns <- setdiff(names(organ_counts), "primary_organ")

        if (length(fallback_columns) == 0) {
            stop("器官样本量工作表缺少可用的计数列")
        }

        count_column <- fallback_columns[[1]]
    }

    organ_counts[[count_column]] <- suppressWarnings(as.numeric(organ_counts[[count_column]]))

    list(
        organ_counts = organ_counts,
        count_column = count_column
    )
}

# 函数：get_heatmap_masked_organs ####
## 功能：识别热图中需要整行置灰的低样本量器官 ####
## 输入：summary_file，summary工作簿路径；min_cases_threshold，样本量阈值 ####
## 输出：包含需置灰器官名称和计数列名称的列表 ####
get_heatmap_masked_organs <- function(summary_file, min_cases_threshold) {
    organ_count_info <- load_heatmap_organ_counts(summary_file)
    organ_counts <- organ_count_info$organ_counts
    count_column <- organ_count_info$count_column

    masked_organs <- organ_counts %>%
        dplyr::filter(is.na(.data[[count_column]]) | .data[[count_column]] < min_cases_threshold) %>%
        dplyr::pull(primary_organ)

    list(
        organ_names = unique(masked_organs),
        count_column = count_column
    )
}

# 函数：model_results_are_ready ####
## 功能：检查指定数据库全部模型结果文件是否已准备完成 ####
## 输入：model_specs，模型配置列表；results_path，结果目录；prefix，文件名前缀 ####
## 输出：逻辑值，TRUE表示全部模型结果可直接读取 ####
model_results_are_ready <- function(model_specs, results_path, prefix) {
    model_result_paths <- get_model_result_paths(model_specs, results_path, prefix)

    all(file.exists(unname(model_result_paths))) &&
        all(purrr::map_lgl(model_result_paths, model_result_excel_is_compatible))
}

# 4. 热图数据整理函数 ####

# 函数：classify_consistency_status ####
## 功能：比较参照模型与当前模型的显著信号方向是否一致 ####
## 输入：reference_effect，参照模型效应值；sensitivity_effect，当前模型效应值；reference_signal，参照模型是否为显著信号；sensitivity_signal，当前模型是否为显著信号 ####
## 输出：字符值same、opposite或unavailable ####
classify_consistency_status <- function(reference_effect, sensitivity_effect, reference_signal, sensitivity_signal) {
    dplyr::case_when(
        is.na(reference_signal) | !reference_signal | is.na(sensitivity_signal) | !sensitivity_signal ~ "unavailable",
        !is.finite(reference_effect) | !is.finite(sensitivity_effect) ~ "unavailable",
        reference_effect <= 0 | sensitivity_effect <= 0 ~ "unavailable",
        reference_effect == 1 | sensitivity_effect == 1 ~ "unavailable",
        (reference_effect > 1 & sensitivity_effect > 1) |
            (reference_effect < 1 & sensitivity_effect < 1) ~ "same",
        TRUE ~ "opposite"
    )
}

# 函数：get_model_signal_column ####
## 功能：根据模型结果中的效应值列选择对应的显著性信号列 ####
## 输入：model_details，模型结果明细表；effect_column，当前使用的效应值列名 ####
## 输出：字符型显著性信号列名 ####
get_model_signal_column <- function(model_details, effect_column) {
    if (identical(effect_column, "adjusted_OR") && "adjusted_signal_detected" %in% names(model_details)) {
        return("adjusted_signal_detected")
    }

    "signal_detected"
}

# 函数：cluster_heatmap_axes ####
## 功能：基于M0的log2(ROR)结果生成热图行列排序 ####
## 输入：details，M0结果明细表；row_var，行变量名；col_var，列变量名；value_var，聚类数值列名 ####
## 输出：包含行列顺序、聚类矩阵和hclust对象的列表 ####
cluster_heatmap_axes <- function(details, row_var = "cancer_type", col_var = "irae", value_var = "ROR") {
    clustering_df <- details %>%
        dplyr::transmute(
            row_name = as.character(.data[[row_var]]),
            col_name = as.character(.data[[col_var]]),
            cluster_value = dplyr::case_when(
                is.na(.data[[value_var]]) | !is.finite(.data[[value_var]]) | .data[[value_var]] <= 0 ~ NA_real_,
                TRUE ~ log2(.data[[value_var]])
            )
        ) %>%
        dplyr::distinct(row_name, col_name, .keep_all = TRUE)

    if (nrow(clustering_df) == 0) {
        stop("热图聚类失败：缺少可用的M0结果数据")
    }

    heatmap_matrix <- clustering_df %>%
        tidyr::pivot_wider(names_from = col_name, values_from = cluster_value) %>%
        tibble::column_to_rownames("row_name") %>%
        as.matrix()

    matrix_for_cluster <- heatmap_matrix
    matrix_for_cluster[is.na(matrix_for_cluster)] <- 0

    row_levels <- rownames(matrix_for_cluster)
    row_hclust <- NULL
    if (nrow(matrix_for_cluster) > 1) {
        row_hclust <- stats::hclust(stats::dist(matrix_for_cluster))
        row_levels <- rownames(matrix_for_cluster)[row_hclust$order]
    }

    col_levels <- colnames(matrix_for_cluster)
    col_hclust <- NULL
    if (ncol(matrix_for_cluster) > 1) {
        col_hclust <- stats::hclust(stats::dist(t(matrix_for_cluster)))
        col_levels <- colnames(matrix_for_cluster)[col_hclust$order]
    }

    list(
        row_levels = row_levels,
        col_levels = col_levels,
        heatmap_matrix = heatmap_matrix,
        row_hclust = row_hclust,
        col_hclust = col_hclust
    )
}

# 函数：build_sensitivity_heatmap_matrices ####
## 功能：将M0与敏感性模型结果整理为热图矩阵和一致性矩阵 ####
## 输入：reference_details，M0结果；model_results，全部模型结果；sensitivity_models，敏感性模型名；row_levels，行顺序；col_levels，列顺序；blank_rows，需置灰行 ####
## 输出：包含ROR矩阵、标签矩阵、一致性数组和置灰标记的列表 ####
build_sensitivity_heatmap_matrices <- function(reference_details, model_results, sensitivity_models, row_levels, col_levels, blank_rows = NULL) {
    row_index <- stats::setNames(seq_along(row_levels), row_levels)
    col_index <- stats::setNames(seq_along(col_levels), col_levels)
    blank_rows <- unique(as.character(blank_rows))
    blank_rows <- blank_rows[!is.na(blank_rows)]
    blank_row_flags <- row_levels %in% blank_rows
    blank_row_lookup <- stats::setNames(blank_row_flags, row_levels)

    log2_ror_matrix <- matrix(
        NA_real_,
        nrow = length(row_levels),
        ncol = length(col_levels),
        dimnames = list(row_levels, col_levels)
    )
    label_matrix <- matrix(
        "",
        nrow = length(row_levels),
        ncol = length(col_levels),
        dimnames = list(row_levels, col_levels)
    )
    label_color_matrix <- matrix(
        "#1F1F1F",
        nrow = length(row_levels),
        ncol = length(col_levels),
        dimnames = list(row_levels, col_levels)
    )
    consistency_array <- array(
        "unavailable",
        dim = c(length(row_levels), length(col_levels), length(sensitivity_models)),
        dimnames = list(row_levels, col_levels, sensitivity_models)
    )

    reference_irae_display <- get_detail_irae_display(reference_details)
    reference_details_prepared <- reference_details %>%
        dplyr::mutate(
            cancer_type = as.character(cancer_type),
            irae = reference_irae_display
        ) %>%
        dplyr::distinct(cancer_type, irae, .keep_all = TRUE)

    for (i in seq_len(nrow(reference_details_prepared))) {
        current_cancer_type <- reference_details_prepared$cancer_type[i]
        if (isTRUE(unname(blank_row_lookup[current_cancer_type]))) {
            next
        }

        row_position <- unname(row_index[current_cancer_type])
        col_position <- unname(col_index[reference_details_prepared$irae[i]])

        if (is.na(row_position) || is.na(col_position)) {
            next
        }

        current_ror <- reference_details_prepared$ROR[i]
        log2_ror_matrix[row_position, col_position] <- dplyr::case_when(
            is.na(current_ror) | !is.finite(current_ror) | current_ror <= 0 ~ NA_real_,
            TRUE ~ pmax(pmin(log2(current_ror), 2), -2)
        )
        label_matrix[row_position, col_position] <- dplyr::case_when(
            is.na(current_ror) | !is.finite(current_ror) ~ "",
            isTRUE(reference_details_prepared$signal_detected[i]) ~ sprintf("%.2f*", current_ror),
            TRUE ~ sprintf("%.2f", current_ror)
        )
        label_color_matrix[row_position, col_position] <- dplyr::case_when(
            is.na(current_ror) | !is.finite(current_ror) ~ NA_character_,
            isTRUE(reference_details_prepared$signal_detected[i]) ~ "#1F1F1F",
            TRUE ~ "#1F1F1F"
        )
    }

    reference_lookup <- reference_details_prepared %>%
        dplyr::transmute(
            cancer_type = cancer_type,
            irae = irae,
            reference_effect = ROR,
            reference_signal = signal_detected
        )

    for (model_name in sensitivity_models) {
        model_details_raw <- model_results[[model_name]]$details
        model_irae_display <- get_detail_irae_display(model_details_raw)
        effect_column <- if ("adjusted_OR" %in% names(model_details_raw)) "adjusted_OR" else "ROR"
        signal_column <- get_model_signal_column(model_details_raw, effect_column)

        model_details <- model_details_raw %>%
            dplyr::mutate(
                cancer_type = as.character(cancer_type),
                irae = model_irae_display
            ) %>%
            dplyr::left_join(reference_lookup, by = c("cancer_type", "irae")) %>%
            dplyr::mutate(
                sensitivity_effect = .data[[effect_column]],
                sensitivity_signal = .data[[signal_column]],
                consistency_status = classify_consistency_status(
                    reference_effect = reference_effect,
                    sensitivity_effect = sensitivity_effect,
                    reference_signal = reference_signal,
                    sensitivity_signal = sensitivity_signal
                )
            )

        if (nrow(model_details) == 0) {
            next
        }

        for (i in seq_len(nrow(model_details))) {
            current_cancer_type <- model_details$cancer_type[i]
            if (isTRUE(unname(blank_row_lookup[current_cancer_type]))) {
                next
            }

            row_position <- unname(row_index[current_cancer_type])
            col_position <- unname(col_index[model_details$irae[i]])

            if (is.na(row_position) || is.na(col_position)) {
                next
            }

            consistency_value <- model_details$consistency_status[i]
            if (is.na(consistency_value) || consistency_value == "") {
                consistency_value <- "unavailable"
            }

            consistency_array[row_position, col_position, model_name] <- consistency_value
        }
    }

    list(
        log2_ror_matrix = log2_ror_matrix,
        label_matrix = label_matrix,
        label_color_matrix = label_color_matrix,
        consistency_array = consistency_array,
        blank_row_flags = unname(blank_row_flags)
    )
}

# 5. 热图对象构建与保存函数 ####

# 函数：extract_drawn_heatmap_order ####
## 功能：提取ComplexHeatmap实际绘制后的行列显示顺序 ####
## 输入：drawn_heatmap，draw后的ComplexHeatmap对象；row_levels，绘制前矩阵行名顺序；col_levels，绘制前矩阵列名顺序 ####
## 输出：包含实际显示行列顺序的列表 ####
extract_drawn_heatmap_order <- function(drawn_heatmap, row_levels, col_levels) {
    row_order_index <- ComplexHeatmap::row_order(drawn_heatmap)
    col_order_index <- ComplexHeatmap::column_order(drawn_heatmap)

    if (is.list(row_order_index)) {
        row_order_index <- unlist(row_order_index, use.names = FALSE)
    }

    if (is.list(col_order_index)) {
        col_order_index <- unlist(col_order_index, use.names = FALSE)
    }

    row_levels_drawn <- if (is.character(row_order_index)) {
        row_order_index
    } else {
        row_levels[row_order_index]
    }

    col_levels_drawn <- if (is.character(col_order_index)) {
        col_order_index
    } else {
        col_levels[col_order_index]
    }

    list(
        row_levels = row_levels_drawn,
        col_levels = col_levels_drawn,
        heatmap_matrix = NULL,
        row_hclust = NULL,
        col_hclust = NULL
    )
}

# 函数：build_complex_sensitivity_heatmap ####
## 功能：构建带M0数值和敏感性一致性圆点的ComplexHeatmap对象 ####
## 输入：reference_details，M0结果；model_results，全部模型结果；sensitivity_models，敏感性模型名；row_title，行标题；fill_colors，填充色配置；status_colors，一致性色配置；predefined_heatmap_order，预设排序；blank_rows，需置灰行 ####
## 输出：包含热图对象、图例对象和热图排序信息的列表 ####
build_complex_sensitivity_heatmap <- function(reference_details, model_results, sensitivity_models, row_title, fill_colors, status_colors, predefined_heatmap_order = NULL, blank_rows = NULL) {
    reference_irae_display <- get_detail_irae_display(reference_details)
    reference_details_prepared <- reference_details %>%
        dplyr::mutate(
            cancer_type = as.character(cancer_type),
            irae = reference_irae_display
        )

    if (is.null(predefined_heatmap_order)) {
        heatmap_order <- cluster_heatmap_axes(
            details = reference_details_prepared,
            row_var = "cancer_type",
            col_var = "irae",
            value_var = "ROR"
        )
    } else {
        heatmap_order <- predefined_heatmap_order
    }

    if (!is.null(predefined_heatmap_order)) {
        blank_rows <- unique(c(
            blank_rows,
            setdiff(heatmap_order$row_levels, unique(reference_details_prepared$cancer_type))
        ))
    }

    heatmap_matrices <- build_sensitivity_heatmap_matrices(
        reference_details = reference_details,
        model_results = model_results,
        sensitivity_models = sensitivity_models,
        row_levels = heatmap_order$row_levels,
        col_levels = heatmap_order$col_levels,
        blank_rows = blank_rows
    )

    fill_function <- circlize::colorRamp2(
        c(-2, 0, 2),
        c(fill_colors$low, fill_colors$mid, fill_colors$high)
    )

    heatmap_object <- ComplexHeatmap::Heatmap(
        heatmap_matrices$log2_ror_matrix,
        name = "log2(ROR)",
        col = fill_function,
        na_col = fill_colors$na,
        show_heatmap_legend = FALSE,
        cluster_rows = if (is.null(heatmap_order$row_hclust)) FALSE else as.dendrogram(heatmap_order$row_hclust),
        cluster_columns = if (is.null(heatmap_order$col_hclust)) FALSE else as.dendrogram(heatmap_order$col_hclust),
        show_row_dend = !is.null(heatmap_order$row_hclust),
        show_column_dend = !is.null(heatmap_order$col_hclust),
        row_dend_side = "left",
        column_dend_side = "top",
        row_dend_width = grid::unit(22, "mm"),
        column_dend_height = grid::unit(22, "mm"),
        rect_gp = grid::gpar(col = "black", lwd = 0.7),
        row_names_gp = grid::gpar(fontsize = 16, col = "black"),
        column_names_gp = grid::gpar(fontsize = 16, col = "black"),
        column_names_rot = 45,
        row_title = row_title,
        row_title_gp = grid::gpar(fontsize = 18, col = "black"),
        column_title = "SOC",
        column_title_side = "bottom",
        column_title_gp = grid::gpar(fontsize = 18, col = "black"),
        heatmap_legend_param = list(
            at = c(-2, -1, 0, 1, 2),
            direction = "horizontal",
            title_position = "leftcenter",
            title_gp = grid::gpar(fontsize = 18, col = "black"),
            labels_gp = grid::gpar(fontsize = 16, col = "black"),
            legend_width = grid::unit(42, "mm")
        ),
        cell_fun = function(j, i, x, y, width, height, fill) {
            if (isTRUE(heatmap_matrices$blank_row_flags[i])) {
                return(invisible(NULL))
            }

            current_label <- heatmap_matrices$label_matrix[i, j]
            current_label_color <- heatmap_matrices$label_color_matrix[i, j]

            if (!is.na(current_label) && nzchar(current_label)) {
                grid::grid.text(
                    label = current_label,
                    x = x,
                    y = y,
                    gp = grid::gpar(fontsize = 11, fontface = "bold", col = current_label_color)
                )
            }
        }
    )

    heatmap_legend <- ComplexHeatmap::Legend(
        title = "log2(ROR)",
        col_fun = fill_function,
        at = c(-2, -1, 0, 1, 2),
        direction = "horizontal",
        title_position = "leftcenter",
        legend_width = grid::unit(42, "mm"),
        title_gp = grid::gpar(fontsize = 18, col = "black"),
        labels_gp = grid::gpar(fontsize = 16, col = "black")
    )

    list(
        heatmap = heatmap_object,
        heatmap_legend = heatmap_legend,
        heatmap_order = heatmap_order
    )
}

# 函数：build_model_robustness_heatmap ####
## 功能：构建单独展示敏感性模型与M0一致性的模型稳健性热图 ####
## 输入：reference_details，参照M0结果；model_results，待比较模型结果；sensitivity_models，待比较模型名；row_title，行标题；status_colors，一致性色配置；predefined_heatmap_order，预设排序；blank_rows，需置灰行；plot_title，可选图标题；legend_title，图例标题；column_names_rot，列名旋转角度 ####
## 输出：包含模型稳健性ComplexHeatmap对象和热图排序信息的列表 ####
build_model_robustness_heatmap <- function(reference_details, model_results, sensitivity_models, row_title, status_colors, predefined_heatmap_order = NULL, blank_rows = NULL, plot_title = NULL, legend_title = "Concordance with FAERS primary analysis", column_names_rot = 45) {
    if (length(sensitivity_models) == 0) {
        stop("模型稳健性图缺少可展示的敏感性模型")
    }

    reference_irae_display <- get_detail_irae_display(reference_details)
    reference_details_prepared <- reference_details %>%
        dplyr::mutate(
            cancer_type = as.character(cancer_type),
            irae = reference_irae_display
        )

    if (is.null(predefined_heatmap_order)) {
        heatmap_order <- cluster_heatmap_axes(
            details = reference_details_prepared,
            row_var = "cancer_type",
            col_var = "irae",
            value_var = "ROR"
        )
    } else {
        heatmap_order <- predefined_heatmap_order
    }

    if (!is.null(predefined_heatmap_order)) {
        blank_rows <- unique(c(
            blank_rows,
            setdiff(heatmap_order$row_levels, unique(reference_details_prepared$cancer_type))
        ))
    }

    heatmap_matrices <- build_sensitivity_heatmap_matrices(
        reference_details = reference_details,
        model_results = model_results,
        sensitivity_models = sensitivity_models,
        row_levels = heatmap_order$row_levels,
        col_levels = heatmap_order$col_levels,
        blank_rows = blank_rows
    )

    heatmap_list <- list()
    for (model_name in sensitivity_models) {
        model_matrix <- heatmap_matrices$consistency_array[, , model_name, drop = FALSE]
        model_matrix <- model_matrix[, , 1]

        heatmap_list[[model_name]] <- ComplexHeatmap::Heatmap(
            model_matrix,
            name = paste0(model_name, "_concordance"),
            col = status_colors,
            na_col = status_colors[["unavailable"]],
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            show_row_dend = FALSE,
            show_column_dend = FALSE,
            show_row_names = model_name == sensitivity_models[[1]],
            row_names_side = "left",
            row_names_gp = grid::gpar(fontsize = 7, col = "black"),
            column_names_gp = grid::gpar(fontsize = 7, col = "black"),
            column_names_rot = column_names_rot,
            row_title = if (model_name == sensitivity_models[[1]]) row_title else NULL,
            row_title_gp = grid::gpar(fontsize = 7, col = "black"),
            column_title = model_name,
            column_title_gp = grid::gpar(fontsize = 7, col = "black"),
            rect_gp = grid::gpar(col = "black", lwd = 0.7),
            show_heatmap_legend = FALSE
        )
    }

    robustness_heatmap <- heatmap_list[[1]]
    if (length(heatmap_list) > 1) {
        for (i in 2:length(heatmap_list)) {
            robustness_heatmap <- robustness_heatmap + heatmap_list[[i]]
        }
    }

    robustness_legend <- ComplexHeatmap::Legend(
        labels = c("Concordant", "Discordant", "Non-significant or Unavailable"),
        title = legend_title,
        type = "grid",
        direction = "horizontal",
        nrow = 1,
        title_position = "topcenter",
        legend_gp = grid::gpar(
            fill = unname(status_colors[c("same", "opposite", "unavailable")]),
            col = "black"
        ),
        title_gp = grid::gpar(fontsize = 7, col = "black"),
        labels_gp = grid::gpar(fontsize = 7, col = "black")
    )

    list(
        heatmap = robustness_heatmap,
        heatmap_order = heatmap_order,
        robustness_legend = robustness_legend,
        plot_title = plot_title
    )
}

# 函数：save_complex_heatmap_pdf ####
## 功能：将ComplexHeatmap对象和图例保存为PDF ####
## 输入：heatmap_object，热图对象列表；output_file，输出PDF路径；width，宽度；height，高度 ####
## 输出：实际绘制后的行列顺序，副作用为写出PDF文件 ####
save_complex_heatmap_pdf <- function(heatmap_object, output_file, width = 11, height = 12) {
    grDevices::cairo_pdf(output_file, width = width, height = height)
    grid::grid.newpage()

    legend_area_height <- grid::unit(22, "mm")
    page_layout <- grid::grid.layout(
        nrow = 2,
        ncol = 1,
        heights = grid::unit.c(
            grid::unit(1, "npc") - legend_area_height,
            legend_area_height
        )
    )

    packed_legends <- ComplexHeatmap::packLegend(
        heatmap_object$heatmap_legend,
        direction = "vertical",
        row_gap = grid::unit(2, "mm")
    )

    grid::pushViewport(grid::viewport(layout = page_layout))
    grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))

    drawn_heatmap <- ComplexHeatmap::draw(
        heatmap_object$heatmap,
        newpage = FALSE,
        padding = grid::unit(c(2, 4, 2, 4), "mm")
    )
    drawn_heatmap_order <- extract_drawn_heatmap_order(
        drawn_heatmap = drawn_heatmap,
        row_levels = heatmap_object$heatmap_order$row_levels,
        col_levels = heatmap_object$heatmap_order$col_levels
    )

    grid::upViewport()
    grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
    legend_box_width_mm <- max(80, width * 25.4 - 20)
    grid::pushViewport(grid::viewport(
        x = grid::unit(0.5, "npc"),
        y = grid::unit(0.5, "npc"),
        width = grid::unit(legend_box_width_mm, "mm"),
        height = grid::unit(1, "npc"),
        just = c("center", "center")
    ))

    ComplexHeatmap::draw(
        packed_legends,
        x = grid::unit(0.5, "npc"),
        y = grid::unit(0.5, "npc"),
        just = c("center", "center")
    )

    grid::upViewport(3)
    grDevices::dev.off()

    invisible(drawn_heatmap_order)
}

# 函数：save_model_robustness_heatmap_pdf ####
## 功能：将模型稳健性ComplexHeatmap对象保存为PDF ####
## 输入：heatmap_object，模型稳健性热图对象列表；output_file，输出PDF路径；width，宽度（cm）；height，高度（cm） ####
## 输出：无返回值，副作用为写出PDF文件 ####
save_model_robustness_heatmap_pdf <- function(heatmap_object, output_file, width = 12, height = 12) {
    grDevices::cairo_pdf(output_file, width = width / 2.54, height = height / 2.54)

    if (is.null(heatmap_object$plot_title) || !nzchar(heatmap_object$plot_title)) {
        ComplexHeatmap::draw(
            heatmap_object$heatmap,
            heatmap_legend_side = "bottom",
            merge_legends = TRUE,
            heatmap_legend_list = list(heatmap_object$robustness_legend),
            padding = grid::unit(c(2, 4, 8, 4), "mm")
        )
    } else {
        ComplexHeatmap::draw(
            heatmap_object$heatmap,
            heatmap_legend_side = "bottom",
            merge_legends = TRUE,
            heatmap_legend_list = list(heatmap_object$robustness_legend),
            column_title = heatmap_object$plot_title,
            column_title_gp = grid::gpar(fontsize = 7, col = "black"),
            padding = grid::unit(c(2, 4, 8, 4), "mm")
        )
    }

    grDevices::dev.off()
}

# 6. 1图：FAERS primary_organ 热图 ####
message("\n========== FAERS primary_organ 热图开始 ==========")

## 1. 读取并处理绘图数据 ####
faers_analysis_prefix <- "FAERS_primary_organ_SOC_"
faers_summary_file <- file.path(
    script_results_path,
    "FAERS_primary_organ_SOC_summary.xlsx"
)
faers_heatmap_object <- NULL
faers_m0_sensitivity_display_order <- NULL

if (!model_results_are_ready(model_specs, script_results_path, faers_analysis_prefix)) {
    message("FAERS primary_organ 模型结果未就绪，跳过绘图")
} else {
    faers_model_results <- load_model_results(
        model_specs = model_specs,
        results_path = script_results_path,
        prefix = faers_analysis_prefix
    )

    ## 筛选ICI使用例数 > 100 的癌种用于聚类 ####
    faers_organ_count_info <- load_heatmap_organ_counts(faers_summary_file)
    faers_organ_counts <- faers_organ_count_info$organ_counts
    faers_count_column <- faers_organ_count_info$count_column
    faers_qualified_organs <- faers_organ_counts %>%
        dplyr::filter(
            !is.na(.data[[faers_count_column]]),
            .data[[faers_count_column]] >= 100
        ) %>%
        dplyr::pull(primary_organ) %>%
        unique()
    message(
        paste0(
            "FAERS 热图聚类：按器官样本量工作表中的 ",
            faers_count_column,
            " 列筛选 >= 100 例癌种，共 ",
            length(faers_qualified_organs),
            " 个参与聚类"
        )
    )

    faers_reference_details <- faers_model_results$M0$details %>%
        dplyr::filter(as.character(cancer_type) %in% faers_qualified_organs)
    faers_sensitivity_models <- intersect(c("M1", "M2", "M3"), names(faers_model_results))
    faers_fill_colors <- list(
        low = "#8FB6D9",
        mid = "#F6F1E8",
        high = "#E5A089",
        na = "#AFAFAF"
    )
    faers_status_colors <- c(
        "same" = "#2F7D5A",
        "opposite" = "#C2553D",
        "unavailable" = "#AFAFAF"
    )

    ## 同步过滤各敏感性模型的明细表 ####
    faers_model_results_filtered <- faers_model_results
    for (m in faers_sensitivity_models) {
        faers_model_results_filtered[[m]]$details <- faers_model_results[[m]]$details %>%
            dplyr::filter(as.character(cancer_type) %in% faers_qualified_organs)
    }

    ## 2. 绘图 ####
    faers_heatmap_object <- build_complex_sensitivity_heatmap(
        reference_details = faers_reference_details,
        model_results = faers_model_results_filtered,
        sensitivity_models = faers_sensitivity_models,
        row_title = "Primary organ",
        fill_colors = faers_fill_colors,
        status_colors = faers_status_colors
    )

    ## 3. 保存图形 ####
    faers_output_file <- file.path(
        script_results_path,
        "FAERS_primary_organ_SOC_heatmap_M0_sensitivity.pdf"
    )
    faers_m0_sensitivity_display_order <- save_complex_heatmap_pdf(
        heatmap_object = faers_heatmap_object,
        output_file = faers_output_file,
        width = 8,
        height = 12
    )
    message(paste("FAERS primary_organ 热图已保存:", faers_output_file))

    ## 4. 模型稳健性图 ####
    faers_robustness_status_colors <- c(
        "same" = "#8CB4C5",
        "opposite" = "#D28A78",
        "unavailable" = "#BDBDBD"
    )
    faers_robustness_heatmap_object <- build_model_robustness_heatmap(
        reference_details = faers_reference_details,
        model_results = faers_model_results_filtered,
        sensitivity_models = faers_sensitivity_models,
        row_title = NULL,
        status_colors = faers_robustness_status_colors,
        predefined_heatmap_order = faers_m0_sensitivity_display_order,
        column_names_rot = 90
    )

    faers_robustness_output_file <- file.path(
        script_results_path,
        "FAERS_原发器官_SOC模型稳健性图.pdf"
    )
    save_model_robustness_heatmap_pdf(
        heatmap_object = faers_robustness_heatmap_object,
        output_file = faers_robustness_output_file,
        width = 14.65,
        height = 12
    )
    message(paste("FAERS primary_organ 模型稳健性图已保存:", faers_robustness_output_file))
}

# 7. 2图：JADER primary_organ 热图 ####
message("\n========== JADER primary_organ 热图开始 ==========")

## 1. 读取并处理绘图数据 ####
jader_analysis_prefix <- "JADER_primary_organ_SOC_"
jader_summary_file <- file.path(
    script_results_path,
    "JADER_primary_organ_SOC_summary.xlsx"
)

if (is.null(faers_heatmap_object)) {
    message("FAERS primary_organ 热图未生成，无法为 JADER 复用聚类顺序，跳过绘图")
} else if (!model_results_are_ready(model_specs, script_results_path, jader_analysis_prefix)) {
    message("JADER primary_organ 模型结果未就绪，跳过绘图")
} else {
    jader_model_results <- load_model_results(
        model_specs = model_specs,
        results_path = script_results_path,
        prefix = jader_analysis_prefix
    )
    jader_masked_organ_info <- get_heatmap_masked_organs(
        summary_file = jader_summary_file,
        min_cases_threshold = jader_heatmap_min_cases
    )
    jader_blank_rows <- jader_masked_organ_info$organ_names
    message(
        paste0(
            "JADER 热图展示：按器官样本量工作表中的 ",
            jader_masked_organ_info$count_column,
            " 列将 <",
            jader_heatmap_min_cases,
            " 例器官整行置灰，涉及 ",
            length(jader_blank_rows),
            " 个"
        )
    )

    ## JADER 热图复用 FAERS 聚类后的癌种与 SOC 排序 ####
    ## faers_heatmap_object$heatmap_order 包含筛选后参与聚类的 row_levels 与 col_levels，
    ## JADER 不在该列表中的癌种行会因 blank_rows 置灰逻辑自动处理（原有设计）
    jader_reference_details <- jader_model_results$M0$details
    jader_sensitivity_models <- intersect(c("M1", "M2", "M3"), names(jader_model_results))
    jader_robustness_models <- intersect(c("M0", "M1", "M2", "M3"), names(jader_model_results))
    jader_fill_colors <- list(
        low = "#8FB6D9",
        mid = "#F6F1E8",
        high = "#E5A089",
        na = "#AFAFAF"
    )
    jader_status_colors <- c(
        "same" = "#2F7D5A",
        "opposite" = "#C2553D",
        "unavailable" = "#AFAFAF"
    )

    ## 2. 绘图 ####
    jader_heatmap_object <- build_complex_sensitivity_heatmap(
        reference_details = jader_reference_details,
        model_results = jader_model_results,
        sensitivity_models = jader_sensitivity_models,
        row_title = "Primary organ",
        fill_colors = jader_fill_colors,
        status_colors = jader_status_colors,
        predefined_heatmap_order = faers_heatmap_object$heatmap_order,
        blank_rows = jader_blank_rows
    )

    ## 3. 保存图形 ####
    jader_output_file <- file.path(
        script_results_path,
        "JADER_primary_organ_SOC_heatmap_M0_sensitivity.pdf"
    )
    save_complex_heatmap_pdf(
        heatmap_object = jader_heatmap_object,
        output_file = jader_output_file,
        width = 8,
        height = 12
    )
    message(paste("JADER primary_organ 热图已保存:", jader_output_file))

    ## 4. 模型稳健性图 ####
    jader_robustness_status_colors <- c(
        "same" = "#8CB4C5",
        "opposite" = "#D28A78",
        "unavailable" = "#BDBDBD"
    )
    jader_robustness_heatmap_object <- build_model_robustness_heatmap(
        reference_details = faers_reference_details,
        model_results = jader_model_results,
        sensitivity_models = jader_robustness_models,
        row_title = NULL,
        status_colors = jader_robustness_status_colors,
        predefined_heatmap_order = faers_m0_sensitivity_display_order,
        blank_rows = jader_blank_rows,
        legend_title = "Concordance with FAERS primary analysis",
        column_names_rot = 90
    )

    jader_robustness_output_file <- file.path(
        script_results_path,
        "JADER_原发器官_SOC模型稳健性图.pdf"
    )
    save_model_robustness_heatmap_pdf(
        heatmap_object = jader_robustness_heatmap_object,
        output_file = jader_robustness_output_file,
        width = 14.65,
        height = 12
    )
    message(paste("JADER primary_organ 模型稳健性图已保存:", jader_robustness_output_file))
}

message("\n========== primary_organ 热图绘制全部完成 ==========")
