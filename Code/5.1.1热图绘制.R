# 重点PT层级primary_organ热图绘制（FAERS + JADER）####

## 1. 加载配置####
source("Code/config.R")
source("Code/mapping_rules.R")

## 2. 结果目录与模型配置####
analysis_script_name <- "5.1重点PT层级primary_organ不成比例分析"
analysis_script_results_path <- file.path(results_path, analysis_script_name)
reference_script_name <- "3.1primary_organ不成比例分析"
reference_script_results_path <- file.path(results_path, reference_script_name)
script_name <- "5.1.1热图绘制"
script_results_path <- file.path(results_path, script_name)
heatmap_min_cases <- 100

model_specs <- list(
    "M0" = list(adjust_vars = NULL),
    "M1" = list(adjust_vars = c("treatment_strategy")),
    "M2" = list(adjust_vars = c("treatment_strategy", "Age", "Sex")),
    "M3" = list(adjust_vars = c("treatment_strategy", "Age", "Sex", "Chemotherapy", "Targeted_therapy"))
)

if (!dir.exists(script_results_path)) {
    dir.create(script_results_path, recursive = TRUE, showWarnings = FALSE)
}

message("========== 重点PT层级primary_organ热图绘制开始 ==========")
message(paste("结果目录:", script_results_path))


## 3. 辅助函数####

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

## 函数：format_case_n ####
# 功能：格式化样本量，便于图中标注####
# 输入：x，数值型样本量####
# 输出：带千分位的字符型样本量####
format_case_n <- function(x) {
    format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
}

## 函数：get_soc_abbrev_safe ####
# 功能：将SOC名称转换为缩写，无法映射时返回原始名称####
# 输入：irae_name，字符型SOC名称####
# 输出：字符型SOC缩写或原始名称####
get_soc_abbrev_safe <- function(irae_name) {
    soc_full <- gsub("_", " ", irae_name)
    soc_abbrev <- get_soc_abbrev(soc_full)

    if (is.na(soc_abbrev) || soc_abbrev == "") {
        return(irae_name)
    }

    soc_abbrev
}

## 函数：build_target_pt_table ####
# 功能：构建固定纳入的12个重点PT清单，供列顺序与标签映射使用####
# 输入：无####
# 输出：包含PT展示名、变量名、所属SOC及顺序的数据框####
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

## 函数：wrap_pt_labels ####
# 功能：为PT列标签添加自动换行，避免热图顶端文字过长####
# 输入：pt_labels，字符型PT名称向量；width，换行宽度####
# 输出：命名字符向量，名称为原始PT标签，值为换行后的标签####
wrap_pt_labels <- function(pt_labels, width = 16) {
    wrapped_labels <- stringr::str_wrap(as.character(pt_labels), width = width)
    stats::setNames(wrapped_labels, pt_labels)
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
# 功能：从PT结果明细表中获取所属SOC缩写####
# 输入：details，结果明细表####
# 输出：字符型SOC缩写向量####
get_detail_parent_soc_abbrev <- function(details) {
    if ("parent_soc_abbrev" %in% names(details)) {
        return(as.character(details$parent_soc_abbrev))
    }

    rep(NA_character_, nrow(details))
}

## 函数：get_reference_soc_abbrev ####
# 功能：从3.1主分析SOC结果明细表中获取SOC缩写####
# 输入：details，SOC结果明细表####
# 输出：字符型SOC缩写向量####
get_reference_soc_abbrev <- function(details) {
    if ("irae_abbrev" %in% names(details)) {
        return(as.character(details$irae_abbrev))
    }

    vapply(as.character(details$irae), get_soc_abbrev_safe, character(1))
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

## 函数：load_result_from_excel ####
# 功能：从单个模型结果Excel读取热图所需明细表####
# 输入：input_file，模型结果Excel路径####
# 输出：包含details数据框的列表####
load_result_from_excel <- function(input_file) {
    details <- readxl::read_xlsx(input_file, sheet = "详细结果")

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

    list(
        details = tibble::as_tibble(details)
    )
}

## 函数：model_result_excel_is_compatible ####
# 功能：检查结果Excel是否包含PT热图绘制所需字段####
# 输入：input_file，模型结果Excel路径####
# 输出：逻辑值，TRUE表示结果文件可直接用于5.1.1####
model_result_excel_is_compatible <- function(input_file) {
    if (!file.exists(input_file)) {
        return(FALSE)
    }

    tryCatch(
        {
            detail_columns <- names(readxl::read_xlsx(input_file, sheet = "详细结果", n_max = 0))
            required_detail_columns <- c(
                "cancer_type", "irae", "irae_var", "irae_display",
                "parent_soc", "parent_soc_abbrev", "pt_order",
                "ROR", "signal_detected"
            )
            all(required_detail_columns %in% detail_columns)
        },
        error = function(e) {
            FALSE
        }
    )
}

## 函数：load_model_results ####
# 功能：批量读取单个数据库的全部模型结果####
# 输入：model_specs，模型配置列表；results_path，结果目录；prefix，文件名前缀####
# 输出：按模型名称命名的结果列表####
load_model_results <- function(model_specs, results_path, prefix) {
    model_result_paths <- get_model_result_paths(model_specs, results_path, prefix)
    all_model_files_exist <- all(file.exists(unname(model_result_paths)))
    all_model_files_compatible <- all_model_files_exist &&
        all(purrr::map_lgl(model_result_paths, model_result_excel_is_compatible))

    if (!all_model_files_compatible) {
        stop("5.1 PT模型结果Excel工作簿缺失或与当前热图字段要求不兼容，请先运行5.1分析脚本")
    }

    purrr::map(model_result_paths, load_result_from_excel)
}

## 函数：load_reference_soc_details ####
# 功能：读取3.1主分析M0 SOC结果明细表，供5类一致性热图做参照####
# 输入：input_file，3.1的M0结果Excel路径####
# 输出：SOC结果明细表####
load_reference_soc_details <- function(input_file) {
    if (!file.exists(input_file)) {
        stop(paste("3.1参考M0结果不存在:", input_file))
    }

    load_result_from_excel(input_file)$details
}

## 函数：build_reference_organ_row_order ####
# 功能：基于3.1 FAERS主分析热图顺序构建PT热图固定行顺序####
# 输入：reference_summary_file，3.1 summary工作簿路径；reference_m0_model_file，3.1 M0结果路径；min_cases_threshold，参与排序的最小器官样本量####
# 输出：包含固定行顺序且不显示聚类树的热图顺序列表####
build_reference_organ_row_order <- function(reference_summary_file, reference_m0_model_file, min_cases_threshold) {
    if (!file.exists(reference_summary_file)) {
        stop(paste("3.1参考summary文件不存在，请先准备3.1结果:", reference_summary_file))
    }

    if (!file.exists(reference_m0_model_file)) {
        stop(paste("3.1参考M0文件不存在，请先准备3.1结果:", reference_m0_model_file))
    }

    organ_count_info <- load_heatmap_organ_counts(reference_summary_file)
    reference_organs <- organ_count_info$organ_counts %>%
        dplyr::filter(
            !is.na(.data[[organ_count_info$count_column]]),
            .data[[organ_count_info$count_column]] >= min_cases_threshold
        ) %>%
        dplyr::pull(primary_organ) %>%
        unique()

    if (length(reference_organs) == 0) {
        stop("3.1参考热图无可用于排序的癌种")
    }

    reference_details <- load_reference_soc_details(reference_m0_model_file) %>%
        dplyr::mutate(
            cancer_type = as.character(cancer_type),
            irae_display = get_reference_soc_abbrev(.)
        ) %>%
        dplyr::filter(cancer_type %in% reference_organs)

    if (nrow(reference_details) == 0) {
        stop("3.1参考M0结果中未找到可用于构建顺序的数据")
    }

    soc_heatmap_order <- cluster_heatmap_axes(
        details = reference_details,
        row_var = "cancer_type",
        col_var = "irae_display",
        value_var = "ROR"
    )

    reference_heatmap_matrix <- reference_details %>%
        dplyr::transmute(
            cancer_type = as.character(cancer_type),
            irae_display = as.character(irae_display),
            log2_ror = dplyr::case_when(
                is.na(ROR) | !is.finite(ROR) | ROR <= 0 ~ NA_real_,
                TRUE ~ pmax(pmin(log2(ROR), 2), -2)
            )
        ) %>%
        tidyr::pivot_wider(names_from = irae_display, values_from = log2_ror) %>%
        tibble::column_to_rownames("cancer_type") %>%
        as.matrix()

    reference_heatmap_matrix <- reference_heatmap_matrix[
        soc_heatmap_order$row_levels,
        soc_heatmap_order$col_levels,
        drop = FALSE
    ]

    reference_heatmap_object <- ComplexHeatmap::Heatmap(
        reference_heatmap_matrix,
        name = "M0 log2(ROR)",
        col = circlize::colorRamp2(
            c(-2, 0, 2),
            c("#8FB6D9", "#F6F1E8", "#E5A089")
        ),
        na_col = "#AFAFAF",
        show_heatmap_legend = FALSE,
        cluster_rows = if (is.null(soc_heatmap_order$row_hclust)) FALSE else as.dendrogram(soc_heatmap_order$row_hclust),
        cluster_columns = if (is.null(soc_heatmap_order$col_hclust)) FALSE else as.dendrogram(soc_heatmap_order$col_hclust),
        show_row_dend = !is.null(soc_heatmap_order$row_hclust),
        show_column_dend = !is.null(soc_heatmap_order$col_hclust),
        row_dend_side = "left",
        column_dend_side = "top",
        row_dend_width = grid::unit(22, "mm"),
        column_dend_height = grid::unit(22, "mm"),
        rect_gp = grid::gpar(col = "black", lwd = 0.7),
        row_names_gp = grid::gpar(fontsize = 7, col = "black"),
        column_names_gp = grid::gpar(fontsize = 7, col = "black"),
        column_names_rot = 90
    )

    reference_display_order <- extract_complex_heatmap_display_order(
        heatmap_object = reference_heatmap_object,
        row_names = rownames(reference_heatmap_matrix),
        col_names = colnames(reference_heatmap_matrix)
    )

    list(
        row_levels = reference_display_order$row_order,
        row_hclust = NULL,
        col_hclust = NULL
    )
}

## 函数：build_faers_pt_column_order ####
# 功能：基于FAERS重点PT M0结果为PT热图构建固定列顺序####
# 输入：details，FAERS重点PT M0结果明细表；target_pt_table，固定PT清单####
# 输出：按FAERS结果生成的PT列顺序向量####
build_faers_pt_column_order <- function(details, target_pt_table) {
    details_prepared <- details %>%
        dplyr::mutate(
            cancer_type = as.character(cancer_type),
            irae_display = get_detail_irae_display(details)
        ) %>%
        dplyr::distinct(cancer_type, irae_display, .keep_all = TRUE)

    if (nrow(details_prepared) == 0) {
        stop("FAERS重点PT结果为空，无法构建PT列顺序")
    }

    pt_heatmap_order <- cluster_heatmap_axes(
        details = details_prepared,
        row_var = "cancer_type",
        col_var = "irae_display",
        value_var = "ROR"
    )

    ordered_pt <- pt_heatmap_order$col_levels
    missing_pt <- setdiff(target_pt_table$irae_display, ordered_pt)

    c(
        ordered_pt,
        target_pt_table$irae_display[target_pt_table$irae_display %in% missing_pt]
    )
}

## 函数：load_heatmap_organ_counts ####
# 功能：读取汇总工作簿中的器官样本量工作表并识别计数列####
# 输入：summary_file，summary工作簿路径####
# 输出：包含器官样本量数据和计数列名称的列表####
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

## 函数：get_heatmap_masked_organs ####
# 功能：识别热图中需要整行置灰的低样本量器官####
# 输入：summary_file，summary工作簿路径；min_cases_threshold，样本量阈值####
# 输出：包含低样本器官名称和计数列名称的列表####
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

## 函数：model_results_are_ready ####
# 功能：检查指定数据库全部PT模型结果文件是否已准备完成####
# 输入：model_specs，模型配置列表；results_path，结果目录；prefix，文件名前缀####
# 输出：逻辑值，TRUE表示全部模型结果可直接读取####
model_results_are_ready <- function(model_specs, results_path, prefix) {
    model_result_paths <- get_model_result_paths(model_specs, results_path, prefix)

    all(file.exists(unname(model_result_paths))) &&
        all(purrr::map_lgl(model_result_paths, model_result_excel_is_compatible))
}

## 函数：cluster_heatmap_axes ####
# 功能：基于FAERS M0 PT结果的log2(ROR)生成热图行列聚类顺序####
# 输入：details，M0结果明细表；row_var，行变量名；col_var，列变量名；value_var，聚类数值列名####
# 输出：包含行列顺序与聚类对象的列表####
cluster_heatmap_axes <- function(details, row_var = "cancer_type", col_var = "irae_display", value_var = "ROR") {
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
        row_hclust = row_hclust,
        col_hclust = col_hclust
    )
}

## 函数：extract_complex_heatmap_display_order ####
# 功能：提取ComplexHeatmap对象实际绘制后的行列顺序####
# 输入：heatmap_object，ComplexHeatmap对象；row_names，矩阵行名；col_names，矩阵列名####
# 输出：包含实际显示行列顺序的列表####
extract_complex_heatmap_display_order <- function(heatmap_object, row_names, col_names) {
    draw_env <- new.env(parent = emptyenv())

    grid::grid.grabExpr(
        {
            draw_env$drawn_heatmap <- ComplexHeatmap::draw(
                heatmap_object,
                heatmap_legend_side = "bottom",
                annotation_legend_side = "bottom",
                padding = grid::unit(c(2, 4, 2, 4), "mm"),
                newpage = FALSE
            )
        },
        wrap = FALSE
    )

    row_order_index <- ComplexHeatmap::row_order(draw_env$drawn_heatmap)
    col_order_index <- ComplexHeatmap::column_order(draw_env$drawn_heatmap)

    if (is.list(row_order_index)) {
        row_order_index <- unlist(row_order_index, use.names = FALSE)
    }

    if (is.list(col_order_index)) {
        col_order_index <- unlist(col_order_index, use.names = FALSE)
    }

    list(
        row_order = row_names[row_order_index],
        col_order = col_names[col_order_index]
    )
}

## 函数：classify_consistency_status ####
# 功能：比较M0与敏感性模型效应值方向是否一致####
# 输入：reference_effect，M0效应值；sensitivity_effect，敏感性模型效应值####
# 输出：字符值same、opposite或unavailable####
classify_consistency_status <- function(reference_effect, sensitivity_effect) {
    dplyr::case_when(
        !is.finite(reference_effect) | !is.finite(sensitivity_effect) ~ "unavailable",
        reference_effect <= 0 | sensitivity_effect <= 0 ~ "unavailable",
        reference_effect == 1 | sensitivity_effect == 1 ~ "unavailable",
        (reference_effect > 1 & sensitivity_effect > 1) |
            (reference_effect < 1 & sensitivity_effect < 1) ~ "same",
        TRUE ~ "opposite"
    )
}

## 函数：build_sensitivity_heatmap_matrices ####
# 功能：将M0与M1/M2/M3结果整理为PT主热图所需矩阵####
# 输入：reference_details，M0结果；model_results，全部模型结果；sensitivity_models，敏感性模型名；row_levels，行顺序；col_levels，列顺序；blank_rows，需置灰行####
# 输出：包含底图矩阵、标签矩阵、一致性数组与置灰标记的列表####
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
        "#9A9A9A",
        nrow = length(row_levels),
        ncol = length(col_levels),
        dimnames = list(row_levels, col_levels)
    )
    consistency_array <- array(
        "unavailable",
        dim = c(length(row_levels), length(col_levels), length(sensitivity_models)),
        dimnames = list(row_levels, col_levels, sensitivity_models)
    )

    reference_details_prepared <- reference_details %>%
        dplyr::mutate(
            cancer_type = as.character(cancer_type),
            irae_display = get_detail_irae_display(reference_details)
        ) %>%
        dplyr::distinct(cancer_type, irae_display, .keep_all = TRUE)

    for (i in seq_len(nrow(reference_details_prepared))) {
        current_cancer_type <- reference_details_prepared$cancer_type[i]
        if (isTRUE(unname(blank_row_lookup[current_cancer_type]))) {
            next
        }

        row_position <- unname(row_index[current_cancer_type])
        col_position <- unname(col_index[reference_details_prepared$irae_display[i]])

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
            TRUE ~ "#9A9A9A"
        )
    }

    reference_lookup <- reference_details_prepared %>%
        dplyr::transmute(
            cancer_type = cancer_type,
            irae_display = irae_display,
            reference_effect = ROR
        )

    for (model_name in sensitivity_models) {
        model_details_raw <- model_results[[model_name]]$details
        effect_column <- if ("adjusted_OR" %in% names(model_details_raw)) "adjusted_OR" else "ROR"

        model_details <- model_details_raw %>%
            dplyr::mutate(
                cancer_type = as.character(cancer_type),
                irae_display = get_detail_irae_display(model_details_raw)
            ) %>%
            dplyr::left_join(reference_lookup, by = c("cancer_type", "irae_display")) %>%
            dplyr::mutate(
                sensitivity_effect = .data[[effect_column]],
                consistency_status = classify_consistency_status(reference_effect, sensitivity_effect)
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
            col_position <- unname(col_index[model_details$irae_display[i]])

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

## 函数：build_complex_pt_heatmap ####
# 功能：构建带M0数值和M1/M2/M3方向一致性圆点的PT主热图####
# 输入：reference_details，M0结果；model_results，全部模型结果；sensitivity_models，敏感性模型名；fill_colors，底图颜色配置；status_colors，一致性点颜色配置；predefined_heatmap_order，预设排序；blank_rows，需置灰行；wrapped_pt_labels，换行后的PT标签####
# 输出：包含热图对象、图例对象和排序信息的列表####
build_complex_pt_heatmap <- function(reference_details, model_results, sensitivity_models, fill_colors, status_colors, predefined_heatmap_order = NULL, blank_rows = NULL, wrapped_pt_labels = NULL) {
    reference_details_prepared <- reference_details %>%
        dplyr::mutate(
            cancer_type = as.character(cancer_type),
            irae_display = get_detail_irae_display(reference_details)
        )

    if (is.null(predefined_heatmap_order)) {
        heatmap_order <- cluster_heatmap_axes(
            details = reference_details_prepared,
            row_var = "cancer_type",
            col_var = "irae_display",
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
    point_offsets <- if (length(sensitivity_models) <= 1) {
        0
    } else {
        seq(-0.30, 0.30, length.out = length(sensitivity_models))
    }

    column_labels <- if (is.null(wrapped_pt_labels)) {
        heatmap_order$col_levels
    } else {
        unname(wrapped_pt_labels[heatmap_order$col_levels])
    }

    heatmap_object <- ComplexHeatmap::Heatmap(
        heatmap_matrices$log2_ror_matrix,
        name = "M0 log2(ROR)",
        col = fill_function,
        na_col = fill_colors$na,
        show_heatmap_legend = FALSE,
        cluster_rows = if (is.null(heatmap_order$row_hclust)) FALSE else as.dendrogram(heatmap_order$row_hclust),
        cluster_columns = if (is.null(heatmap_order$col_hclust)) FALSE else as.dendrogram(heatmap_order$col_hclust),
        row_order = if (is.null(heatmap_order$row_hclust)) match(heatmap_order$row_levels, rownames(heatmap_matrices$log2_ror_matrix)) else NULL,
        column_order = if (is.null(heatmap_order$col_hclust)) match(heatmap_order$col_levels, colnames(heatmap_matrices$log2_ror_matrix)) else NULL,
        show_row_dend = !is.null(heatmap_order$row_hclust),
        show_column_dend = !is.null(heatmap_order$col_hclust),
        row_dend_side = "left",
        column_dend_side = "top",
        row_dend_width = grid::unit(22, "mm"),
        column_dend_height = grid::unit(22, "mm"),
        rect_gp = grid::gpar(col = "black", lwd = 0.7),
        row_names_gp = grid::gpar(fontsize = 7, col = "black"),
        column_names_gp = grid::gpar(fontsize = 7, col = "black"),
        column_names_rot = 90,
        column_labels = column_labels,
        row_title = "Primary organ",
        row_title_gp = grid::gpar(fontsize = 7, col = "black"),
        column_title = "PT",
        column_title_side = "bottom",
        column_title_gp = grid::gpar(fontsize = 7, col = "black"),
        heatmap_legend_param = list(
            at = c(-2, -1, 0, 1, 2),
            direction = "horizontal",
            title_position = "leftcenter",
            title_gp = grid::gpar(fontsize = 7, col = "black"),
            labels_gp = grid::gpar(fontsize = 7, col = "black"),
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
                    y = y + height * 0.17,
                    gp = grid::gpar(fontsize = 7, fontface = "bold", col = current_label_color)
                )
            }

            if (length(sensitivity_models) > 0) {
                for (k in seq_along(sensitivity_models)) {
                    current_status <- heatmap_matrices$consistency_array[i, j, sensitivity_models[k]]
                    point_color <- status_colors[[current_status]]

                    if (is.null(point_color) || is.na(point_color)) {
                        point_color <- status_colors[["unavailable"]]
                    }

                    grid::grid.points(
                        x = x + width * point_offsets[k],
                        y = y - height * 0.22,
                        pch = 16,
                        size = grid::unit(1.9, "mm"),
                        gp = grid::gpar(col = point_color, fill = point_color)
                    )
                }
            }
        }
    )

    heatmap_legend <- ComplexHeatmap::Legend(
        title = "M0 log2(ROR)",
        col_fun = fill_function,
        at = c(-2, -1, 0, 1, 2),
        direction = "horizontal",
        title_position = "leftcenter",
        legend_width = grid::unit(42, "mm"),
        title_gp = grid::gpar(fontsize = 7, col = "black"),
        labels_gp = grid::gpar(fontsize = 7, col = "black")
    )

    consistency_legend <- ComplexHeatmap::Legend(
        labels = c("Concordant", "Discordant", "Unavailable"),
        title = "Concordance with M0",
        type = "points",
        direction = "horizontal",
        nrow = 1,
        title_position = "leftcenter",
        pch = 16,
        legend_gp = grid::gpar(
            col = unname(status_colors[c("same", "opposite", "unavailable")]),
            fill = unname(status_colors[c("same", "opposite", "unavailable")])
        ),
        title_gp = grid::gpar(fontsize = 7, col = "black"),
        labels_gp = grid::gpar(fontsize = 7, col = "black")
    )

    list(
        heatmap = heatmap_object,
        legend_list = list(heatmap_legend, consistency_legend),
        heatmap_order = heatmap_order
    )
}

## 函数：classify_pt_soc_concordance ####
# 功能：根据3.1对应SOC主分析与5.1对应PT结果判定5类一致性标签####
# 输入：soc_effect，SOC的ROR；soc_p_value，SOC原始P值；pt_effect，PT的ROR；pt_p_value，PT原始P值####
# 输出：字符型5类标签####
classify_pt_soc_concordance <- function(soc_effect, soc_p_value, pt_effect, pt_p_value) {
    soc_positive <- is.finite(soc_effect) && soc_effect > 1 && is.finite(soc_p_value) && soc_p_value < 0.05
    soc_negative <- is.finite(soc_effect) && soc_effect < 1 && is.finite(soc_p_value) && soc_p_value < 0.05
    pt_positive <- is.finite(pt_effect) && pt_effect > 1 && is.finite(pt_p_value) && pt_p_value < 0.05
    pt_negative <- is.finite(pt_effect) && pt_effect < 1 && is.finite(pt_p_value) && pt_p_value < 0.05

    dplyr::case_when(
        soc_positive && pt_positive ~ "SOC+&PT+",
        soc_positive && pt_negative ~ "SOC+&PT-",
        soc_negative && pt_positive ~ "SOC-&PT+",
        soc_negative && pt_negative ~ "SOC-&PT-",
        TRUE ~ "Non-significant or Unavailable"
    )
}

## 函数：build_pt_soc_concordance_plot_data ####
# 功能：构建5类一致性热图所需状态矩阵####
# 输入：pt_details，5.1的M0 PT结果；reference_soc_details，3.1的M0 SOC结果；heatmap_order，主热图顺序；blank_rows，需置灰行####
# 输出：包含状态矩阵、行列顺序和换行标签的列表####
build_pt_soc_concordance_plot_data <- function(pt_details, reference_soc_details, heatmap_order, blank_rows = NULL, wrapped_pt_labels = NULL) {
    pt_details_prepared <- pt_details %>%
        dplyr::mutate(
            cancer_type = as.character(cancer_type),
            irae_display = get_detail_irae_display(pt_details),
            parent_soc_abbrev = get_detail_parent_soc_abbrev(pt_details),
            ROR = suppressWarnings(as.numeric(ROR)),
            p_value = suppressWarnings(as.numeric(p_value))
        ) %>%
        dplyr::distinct(cancer_type, irae_display, parent_soc_abbrev, .keep_all = TRUE)

    reference_lookup <- reference_soc_details %>%
        dplyr::mutate(
            cancer_type = as.character(cancer_type),
            parent_soc_abbrev = get_reference_soc_abbrev(reference_soc_details),
            soc_ROR = suppressWarnings(as.numeric(ROR)),
            soc_p_value = suppressWarnings(as.numeric(p_value))
        ) %>%
        dplyr::distinct(cancer_type, parent_soc_abbrev, .keep_all = TRUE) %>%
        dplyr::select(cancer_type, parent_soc_abbrev, soc_ROR, soc_p_value)

    pt_details_joined <- pt_details_prepared %>%
        dplyr::left_join(reference_lookup, by = c("cancer_type", "parent_soc_abbrev"))

    unmatched_rows <- pt_details_joined %>%
        dplyr::filter(is.na(soc_ROR) | is.na(soc_p_value))

    if (nrow(unmatched_rows) > 0) {
        stop(
            paste0(
                "3.1与5.1结果按 primary_organ + parent_soc_abbrev 匹配失败，未映射记录数: ",
                nrow(unmatched_rows)
            )
        )
    }

    row_levels <- heatmap_order$row_levels
    col_levels <- heatmap_order$col_levels
    blank_rows <- unique(as.character(blank_rows))
    blank_rows <- blank_rows[!is.na(blank_rows)]
    blank_row_flags <- row_levels %in% blank_rows
    blank_row_lookup <- stats::setNames(blank_row_flags, row_levels)

    status_matrix <- matrix(
        "Non-significant or Unavailable",
        nrow = length(row_levels),
        ncol = length(col_levels),
        dimnames = list(row_levels, col_levels)
    )

    for (i in seq_len(nrow(pt_details_joined))) {
        current_row <- pt_details_joined[i, ]

        if (!current_row$cancer_type %in% row_levels) {
            next
        }

        if (!current_row$irae_display %in% col_levels) {
            next
        }

        if (isTRUE(unname(blank_row_lookup[current_row$cancer_type]))) {
            next
        }

        status_matrix[current_row$cancer_type, current_row$irae_display] <- classify_pt_soc_concordance(
            soc_effect = current_row$soc_ROR,
            soc_p_value = current_row$soc_p_value,
            pt_effect = current_row$ROR,
            pt_p_value = current_row$p_value
        )
    }

    column_labels <- if (is.null(wrapped_pt_labels)) {
        col_levels
    } else {
        unname(wrapped_pt_labels[col_levels])
    }

    list(
        status_matrix = status_matrix,
        row_order = row_levels,
        col_order = col_levels,
        column_labels = column_labels,
        blank_row_flags = blank_row_flags
    )
}

## 函数：build_pt_soc_concordance_heatmap ####
# 功能：构建5类标签一致性热图对象####
# 输入：plot_data，状态矩阵与顺序信息；status_colors，标签颜色映射；column_title，列标题；column_title_side，列标题位置；show_row_names，是否显示行名；row_names_side，行名位置；row_title，行标题；heatmap_name，热图名称####
# 输出：包含热图对象和图例对象的列表####
build_pt_soc_concordance_heatmap <- function(plot_data, status_colors, column_title = "PT", column_title_side = "bottom", show_row_names = TRUE, row_names_side = "right", row_title = "Primary organ", heatmap_name = "Concordance") {
    status_levels <- c("SOC+&PT+", "SOC+&PT-", "SOC-&PT+", "SOC-&PT-", "Non-significant or Unavailable")
    status_matrix_factor <- matrix(
        factor(plot_data$status_matrix, levels = status_levels),
        nrow = nrow(plot_data$status_matrix),
        ncol = ncol(plot_data$status_matrix),
        dimnames = dimnames(plot_data$status_matrix)
    )

    heatmap_object <- ComplexHeatmap::Heatmap(
        status_matrix_factor,
        name = heatmap_name,
        col = status_colors,
        show_heatmap_legend = FALSE,
        rect_gp = grid::gpar(col = "black", lwd = 0.7),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_order = match(plot_data$row_order, rownames(status_matrix_factor)),
        column_order = match(plot_data$col_order, colnames(status_matrix_factor)),
        show_row_names = show_row_names,
        row_names_side = row_names_side,
        row_names_gp = grid::gpar(fontsize = 7, col = "black"),
        column_names_gp = grid::gpar(fontsize = 7, col = "black"),
        column_names_rot = 90,
        column_labels = plot_data$column_labels,
        row_title = row_title,
        row_title_gp = grid::gpar(fontsize = 7, col = "black"),
        column_title = column_title,
        column_title_side = column_title_side,
        column_title_gp = grid::gpar(fontsize = 7, col = "black")
    )

    state_legend <- ComplexHeatmap::Legend(
        labels = names(status_colors),
        title = "Concordance with FAERS primary analysis",
        type = "grid",
        nrow = 2,
        legend_gp = grid::gpar(fill = unname(status_colors), col = "black"),
        title_gp = grid::gpar(fontsize = 7, col = "black"),
        labels_gp = grid::gpar(fontsize = 7, col = "black")
    )

    list(
        heatmap = heatmap_object,
        legend_list = list(state_legend)
    )
}

## 函数：save_complex_heatmap_pdf ####
# 功能：将热图对象及图例保存为PDF####
# 输入：heatmap_object，包含热图对象和图例列表的列表；output_file，输出PDF路径；width，宽度（英寸）；height，高度（英寸）####
# 输出：无返回值，副作用为写出PDF文件####
save_complex_heatmap_pdf <- function(heatmap_object, output_file, width = 7 / 2.54, height = 12) {
    ensure_parent_dir(output_file)

    grDevices::cairo_pdf(output_file, width = width, height = height)
    grid::grid.newpage()

    legend_area_height <- grid::unit(34, "mm")
    page_layout <- grid::grid.layout(
        nrow = 2,
        ncol = 1,
        heights = grid::unit.c(
            grid::unit(1, "npc") - legend_area_height,
            legend_area_height
        )
    )

    packed_legends <- do.call(
        ComplexHeatmap::packLegend,
        c(
            heatmap_object$legend_list,
            list(
                direction = "vertical",
                row_gap = grid::unit(2, "mm")
            )
        )
    )

    grid::pushViewport(grid::viewport(layout = page_layout))
    grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))

    ComplexHeatmap::draw(
        heatmap_object$heatmap,
        newpage = FALSE,
        padding = grid::unit(c(2, 4, 2, 4), "mm")
    )

    grid::upViewport()
    grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
    legend_box_width_mm <- max(40, width * 25.4 - 8)
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
}

## 函数：save_combined_pt_soc_concordance_pdf ####
# 功能：将FAERS和JADER 5类标签一致性热图横向合并保存为一个PDF，并共用一套底部图例####
# 输入：faers_heatmap_object，FAERS热图与图例对象；jader_heatmap_object，JADER热图与图例对象；output_file，输出PDF路径；width，宽度（英寸）；height，高度（英寸）####
# 输出：无返回值，副作用为写出PDF文件####
save_combined_pt_soc_concordance_pdf <- function(faers_heatmap_object, jader_heatmap_object, output_file, width = 14.65 / 2.54, height = 14 / 2.54) {
    ensure_parent_dir(output_file)

    combined_heatmap <- faers_heatmap_object$heatmap + jader_heatmap_object$heatmap

    grDevices::cairo_pdf(output_file, width = width, height = height)
    grid::grid.newpage()

    legend_area_height <- grid::unit(28, "mm")
    page_layout <- grid::grid.layout(
        nrow = 2,
        ncol = 1,
        heights = grid::unit.c(
            grid::unit(1, "npc") - legend_area_height,
            legend_area_height
        )
    )

    packed_legends <- do.call(
        ComplexHeatmap::packLegend,
        c(
            faers_heatmap_object$legend_list,
            list(
                direction = "vertical",
                row_gap = grid::unit(2, "mm")
            )
        )
    )

    grid::pushViewport(grid::viewport(layout = page_layout))
    grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))

    ComplexHeatmap::draw(
        combined_heatmap,
        newpage = FALSE,
        padding = grid::unit(c(2, 4, 2, 4), "mm")
    )

    grid::upViewport()
    grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
    legend_box_width_mm <- max(80, width * 25.4 - 8)
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
}


## 4. 固定PT清单与标签配置####
target_pt_table <- build_target_pt_table()
wrapped_pt_labels <- stats::setNames(
    as.character(target_pt_table$irae_display),
    as.character(target_pt_table$irae_display)
)
reference_organ_row_order <- build_reference_organ_row_order(
    reference_summary_file = file.path(
        reference_script_results_path,
        "FAERS_primary_organ_SOC_summary.xlsx"
    ),
    reference_m0_model_file = file.path(
        reference_script_results_path,
        "FAERS_primary_organ_SOC_M0.xlsx"
    ),
    min_cases_threshold = heatmap_min_cases
)
pt_fill_colors <- list(
    low = "#8FB6D9",
    mid = "#F6F1E8",
    high = "#E5A089",
    na = "#AFAFAF"
)
pt_status_colors <- c(
    "same" = "#2F7D5A",
    "opposite" = "#C2553D",
    "unavailable" = "#AFAFAF"
)
five_class_colors <- c(
    "SOC+&PT+" = "#2F7D5A",
    "SOC+&PT-" = "#C2553D",
    "SOC-&PT+" = "#D98E04",
    "SOC-&PT-" = "#4C78A8",
    "Non-significant or Unavailable" = "#AFAFAF"
)


## 5. 图1：FAERS重点PT primary_organ主热图####
message("\n========== 图1：FAERS重点PT primary_organ主热图开始 ==========")

## 1. 读取并处理绘图数据####
faers_analysis_prefix <- "FAERS_primary_organ_PT_"
faers_summary_file <- file.path(
    analysis_script_results_path,
    "FAERS_primary_organ_PT_summary.xlsx"
)
faers_main_heatmap_object <- NULL

if (!model_results_are_ready(model_specs, analysis_script_results_path, faers_analysis_prefix)) {
    stop("FAERS 5.1 PT模型结果未就绪，请先运行5.1重点PT层分析脚本")
}

faers_model_results <- load_model_results(
    model_specs = model_specs,
    results_path = analysis_script_results_path,
    prefix = faers_analysis_prefix
)

faers_organ_count_info <- load_heatmap_organ_counts(faers_summary_file)
faers_qualified_organs <- faers_organ_count_info$organ_counts %>%
    dplyr::filter(
        !is.na(.data[[faers_organ_count_info$count_column]]),
        .data[[faers_organ_count_info$count_column]] >= heatmap_min_cases
    ) %>%
    dplyr::pull(primary_organ) %>%
    unique()

faers_reference_details <- faers_model_results$M0$details %>%
    dplyr::filter(as.character(cancer_type) %in% faers_qualified_organs)

faers_model_results_filtered <- faers_model_results
for (model_name in intersect(c("M1", "M2", "M3"), names(faers_model_results_filtered))) {
    faers_model_results_filtered[[model_name]]$details <- faers_model_results_filtered[[model_name]]$details %>%
        dplyr::filter(as.character(cancer_type) %in% faers_qualified_organs)
}

faers_pt_column_order <- build_faers_pt_column_order(
    details = faers_reference_details,
    target_pt_table = target_pt_table
)
faers_reference_pt_heatmap_order <- list(
    row_levels = reference_organ_row_order$row_levels,
    col_levels = faers_pt_column_order,
    row_hclust = NULL,
    col_hclust = NULL
)

## 2. 绘图####
faers_main_heatmap_object <- build_complex_pt_heatmap(
    reference_details = faers_reference_details,
    model_results = faers_model_results_filtered,
    sensitivity_models = intersect(c("M1", "M2", "M3"), names(faers_model_results_filtered)),
    fill_colors = pt_fill_colors,
    status_colors = pt_status_colors,
    predefined_heatmap_order = faers_reference_pt_heatmap_order,
    blank_rows = NULL,
    wrapped_pt_labels = wrapped_pt_labels
)

## 3. 保存图形####
save_complex_heatmap_pdf(
    heatmap_object = faers_main_heatmap_object,
    output_file = file.path(
        script_results_path,
        "FAERS_primary_organ_PT_heatmap_M0_sensitivity.pdf"
    ),
    width = 7 / 2.54,
    height = 12
)


## 6. 图2：JADER重点PT primary_organ主热图####
message("\n========== 图2：JADER重点PT primary_organ主热图开始 ==========")

## 1. 读取并处理绘图数据####
jader_analysis_prefix <- "JADER_primary_organ_PT_"
jader_summary_file <- file.path(
    analysis_script_results_path,
    "JADER_primary_organ_PT_summary.xlsx"
)

if (!model_results_are_ready(model_specs, analysis_script_results_path, jader_analysis_prefix)) {
    stop("JADER 5.1 PT模型结果未就绪，请先运行5.1重点PT层分析脚本")
}

jader_model_results <- load_model_results(
    model_specs = model_specs,
    results_path = analysis_script_results_path,
    prefix = jader_analysis_prefix
)

jader_masked_organ_info <- get_heatmap_masked_organs(
    summary_file = jader_summary_file,
    min_cases_threshold = heatmap_min_cases
)
jader_blank_rows <- jader_masked_organ_info$organ_names

## 2. 绘图####
jader_main_heatmap_object <- build_complex_pt_heatmap(
    reference_details = jader_model_results$M0$details,
    model_results = jader_model_results,
    sensitivity_models = intersect(c("M1", "M2", "M3"), names(jader_model_results)),
    fill_colors = pt_fill_colors,
    status_colors = pt_status_colors,
    predefined_heatmap_order = faers_reference_pt_heatmap_order,
    blank_rows = jader_blank_rows,
    wrapped_pt_labels = wrapped_pt_labels
)

## 3. 保存图形####
save_complex_heatmap_pdf(
    heatmap_object = jader_main_heatmap_object,
    output_file = file.path(
        script_results_path,
        "JADER_primary_organ_PT_heatmap_M0_sensitivity.pdf"
    ),
    width = 7 / 2.54,
    height = 12
)


## 7. 图3：FAERS重点PT与3.1主分析5类一致性热图####
message("\n========== 图3：FAERS重点PT与3.1主分析5类一致性热图开始 ==========")

## 1. 读取并处理绘图数据####
faers_reference_soc_details <- load_reference_soc_details(
    file.path(
        reference_script_results_path,
        "FAERS_primary_organ_SOC_M0.xlsx"
    )
)

faers_concordance_plot_data <- build_pt_soc_concordance_plot_data(
    pt_details = faers_model_results$M0$details,
    reference_soc_details = faers_reference_soc_details,
    heatmap_order = faers_main_heatmap_object$heatmap_order,
    blank_rows = NULL,
    wrapped_pt_labels = wrapped_pt_labels
)

## 2. 绘图####
faers_concordance_heatmap_object <- build_pt_soc_concordance_heatmap(
    plot_data = faers_concordance_plot_data,
    status_colors = five_class_colors
)

## 3. 保存图形####
save_complex_heatmap_pdf(
    heatmap_object = faers_concordance_heatmap_object,
    output_file = file.path(
        script_results_path,
        "FAERS_primary_organ_PT_SOC_concordance_heatmap_5class.pdf"
    ),
    width = 7 / 2.54,
    height = 12
)


## 8. 图4：JADER重点PT与3.1主分析5类一致性热图####
message("\n========== 图4：JADER重点PT与3.1主分析5类一致性热图开始 ==========")

## 1. 读取并处理绘图数据####
jader_reference_soc_details <- load_reference_soc_details(
    file.path(
        reference_script_results_path,
        "JADER_primary_organ_SOC_M0.xlsx"
    )
)

jader_concordance_plot_data <- build_pt_soc_concordance_plot_data(
    pt_details = jader_model_results$M0$details,
    reference_soc_details = jader_reference_soc_details,
    heatmap_order = jader_main_heatmap_object$heatmap_order,
    blank_rows = jader_blank_rows,
    wrapped_pt_labels = wrapped_pt_labels
)

## 2. 绘图####
jader_concordance_heatmap_object <- build_pt_soc_concordance_heatmap(
    plot_data = jader_concordance_plot_data,
    status_colors = five_class_colors
)

## 3. 保存图形####
save_complex_heatmap_pdf(
    heatmap_object = jader_concordance_heatmap_object,
    output_file = file.path(
        script_results_path,
        "JADER_primary_organ_PT_SOC_concordance_heatmap_5class.pdf"
    ),
    width = 7 / 2.54,
    height = 12
)


## 9. 图5：FAERS与JADER重点PT和SOC一致性5分类合并热图####
message("\n========== 图5：FAERS与JADER重点PT和SOC一致性5分类合并热图开始 ==========")

## 1. 读取并处理绘图数据####
faers_combined_concordance_plot_data <- faers_concordance_plot_data
jader_combined_concordance_plot_data <- jader_concordance_plot_data

## 2. 绘图####
faers_combined_concordance_heatmap_object <- build_pt_soc_concordance_heatmap(
    plot_data = faers_combined_concordance_plot_data,
    status_colors = five_class_colors,
    column_title = "FAERS",
    column_title_side = "top",
    show_row_names = TRUE,
    row_names_side = "left",
    row_title = NULL,
    heatmap_name = "FAERS_concordance"
)

jader_combined_concordance_heatmap_object <- build_pt_soc_concordance_heatmap(
    plot_data = jader_combined_concordance_plot_data,
    status_colors = five_class_colors,
    column_title = "JADER",
    column_title_side = "top",
    show_row_names = FALSE,
    row_title = NULL,
    heatmap_name = "JADER_concordance"
)

## 3. 保存图形####
save_combined_pt_soc_concordance_pdf(
    faers_heatmap_object = faers_combined_concordance_heatmap_object,
    jader_heatmap_object = jader_combined_concordance_heatmap_object,
    output_file = file.path(
        script_results_path,
        "FAERS和JADER重点PT与SOC一致性5分类合并热图.pdf"
    ),
    width = 14.65 / 2.54,
    height = 14 / 2.54
)

message("\n========== 重点PT层级primary_organ热图绘制完成 ==========")
