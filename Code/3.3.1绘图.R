# 3.3.1 绘图（FAERS）####

rm(list = ls())
gc()

source("Code/config.R")
source("Code/mapping_rules.R")

# 0. 结果目录与输入文件 ####
script_name <- "3.3同一ICI_regimen下primary_organ不成比例分析"
script_results_path <- file.path(results_path, script_name)
dataset_name <- "FAERS"
summary_file <- file.path(
    script_results_path,
    paste0(dataset_name, "_within_ICI_regimen_primary_organ_SOC_summary.xlsx")
)
m0_model_file <- file.path(
    script_results_path,
    paste0(dataset_name, "_within_ICI_regimen_primary_organ_SOC_M0.xlsx")
)

reference_script_name <- "3.1primary_organ不成比例分析"
reference_script_results_path <- file.path(results_path, reference_script_name)
reference_summary_file <- file.path(
    reference_script_results_path,
    "FAERS_primary_organ_SOC_summary.xlsx"
)
reference_m0_model_file <- file.path(
    reference_script_results_path,
    "FAERS_primary_organ_SOC_M0.xlsx"
)
reference_heatmap_min_cases <- 100L

if (!dir.exists(script_results_path)) {
    dir.create(script_results_path, recursive = TRUE, showWarnings = FALSE)
}

message("========== 3.3.1 绘图开始 ==========")
message(paste("结果目录:", script_results_path))


# 1. 公共辅助函数 ####

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

## 格式化样本量，输入x样本量数值，输出带千分位字符值 ####
format_case_n <- function(x) {
    format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
}

## 构建带样本量的热图标题，输入display_label与total_n，输出标题字符串 ####
build_heatmap_title_with_n <- function(display_label, total_n) {
    paste0(display_label, "\n（N=", format_case_n(total_n), "）")
}

## 清洗文件名片段，输入x原始文本，输出适合文件名的标签 ####
sanitize_file_label <- function(x) {
    cleaned_label <- gsub("[^[:alnum:]]+", "_", x)
    cleaned_label <- gsub("^_+|_+$", "", cleaned_label)
    gsub("_+", "_", cleaned_label)
}

## 标准化ICI regimen展示文本，输入x原始regimen字符向量，输出展示用字符向量 ####
format_ici_regimen_label <- function(x) {
    x %>%
        stringr::str_replace_all("anti-PD-1", "anti–PD-1") %>%
        stringr::str_replace_all("anti-PD-L1", "anti–PD-L1") %>%
        stringr::str_replace_all("anti-CTLA-4", "anti–CTLA-4") %>%
        stringr::str_replace_all("\\+", " + ")
}

## 构建比较标识，输入ici_regimen_display与primary_organ，输出comparison_id ####
build_comparison_id <- function(ici_regimen_display, primary_organ) {
    paste(ici_regimen_display, primary_organ, sep = " | ")
}

## 获取单个SOC缩写，输入irae_name原始SOC名称，输出SOC缩写或原始名称 ####
get_soc_abbrev_safe <- function(irae_name) {
    soc_full <- gsub("_", " ", irae_name)
    soc_abbrev <- get_soc_abbrev(soc_full)

    if (is.na(soc_abbrev) || soc_abbrev == "") {
        return(irae_name)
    }

    soc_abbrev
}

## 获取绘图使用的SOC展示标签，输入details结果明细表，输出字符型SOC展示标签向量 ####
get_detail_irae_display <- function(details) {
    if ("irae_abbrev" %in% names(details)) {
        return(as.character(details$irae_abbrev))
    }

    vapply(as.character(details$irae), get_soc_abbrev_safe, character(1))
}

## 读取summary工作簿，输入summary_file路径，输出regimen与器官汇总表 ####
load_summary_tables <- function(summary_file) {
    if (!file.exists(summary_file)) {
        stop(paste("summary文件不存在:", summary_file))
    }

    regimen_status <- readxl::read_xlsx(summary_file, sheet = "regimen纳入情况") %>%
        mutate(
            ici_regimen_raw = as.character(ici_regimen_raw),
            ici_regimen_display = as.character(ici_regimen_display),
            regimen_rank = suppressWarnings(as.numeric(regimen_rank)),
            regimen_total_n = suppressWarnings(as.numeric(regimen_total_n)),
            n_primary_organs = suppressWarnings(as.numeric(n_primary_organs)),
            included_in_main = as.logical(included_in_main),
            regimen_status = as.character(regimen_status)
        )

    organ_status <- readxl::read_xlsx(summary_file, sheet = "regimen内器官样本量") %>%
        mutate(
            ici_regimen_raw = as.character(ici_regimen_raw),
            ici_regimen_display = as.character(ici_regimen_display),
            regimen_rank = suppressWarnings(as.numeric(regimen_rank)),
            primary_organ = as.character(primary_organ),
            organ_n = suppressWarnings(as.numeric(organ_n)),
            regimen_total_n = suppressWarnings(as.numeric(regimen_total_n)),
            organ_pct = suppressWarnings(as.numeric(organ_pct)),
            organ_rank = suppressWarnings(as.numeric(organ_rank)),
            included_in_main = as.logical(included_in_main),
            comparison_id = as.character(comparison_id)
        )

    list(
        regimen_status = regimen_status,
        organ_status = organ_status
    )
}

## 读取热图器官样本量工作表，输入summary_file路径，输出器官样本量数据与计数列名 ####
load_heatmap_organ_counts <- function(summary_file) {
    if (!file.exists(summary_file)) {
        stop(paste("热图汇总文件不存在:", summary_file))
    }

    organ_counts <- readxl::read_xlsx(summary_file, sheet = "器官样本量") %>%
        mutate(primary_organ = as.character(primary_organ))

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

## 读取M0详细结果，输入m0_model_file路径，输出标准化后的详细结果表 ####
load_m0_model_details <- function(m0_model_file) {
    if (!file.exists(m0_model_file)) {
        stop(paste("M0结果文件不存在:", m0_model_file))
    }

    details_df <- readxl::read_xlsx(m0_model_file, sheet = "详细结果")

    if (nrow(details_df) == 0) {
        stop(paste("M0详细结果为空:", m0_model_file))
    }

    character_cols <- intersect(
        c(
            "database", "cancer_type", "irae", "irae_var", "irae_full", "irae_abbrev",
            "analysis_method", "p_value_method", "result_note",
            "ici_regimen_raw", "ici_regimen_display", "primary_organ",
            "comparison_id", "adjust_vars_used"
        ),
        names(details_df)
    )
    logical_cols <- intersect(
        c(
            "signal_detected", "signal_detected_within_regimen", "fisher_test_used",
            "adjusted_model_fitted", "adjusted_signal_detected",
            "adjusted_signal_detected_within_regimen"
        ),
        names(details_df)
    )
    numeric_cols <- setdiff(names(details_df), c(character_cols, logical_cols))

    if (length(character_cols) > 0) {
        details_df[character_cols] <- lapply(details_df[character_cols], as.character)
    }

    if (length(logical_cols) > 0) {
        details_df[logical_cols] <- lapply(details_df[logical_cols], as.logical)
    }

    if (length(numeric_cols) > 0) {
        details_df[numeric_cols] <- lapply(details_df[numeric_cols], function(current_col) {
            suppressWarnings(as.numeric(current_col))
        })
    }

    tibble::as_tibble(details_df)
}

## 构建primary_organ配色，输入primary_organs器官向量，输出命名颜色向量 ####
build_primary_organ_color_map <- function(primary_organs) {
    primary_organs <- unique(as.character(primary_organs))
    primary_organs <- primary_organs[!is.na(primary_organs) & primary_organs != ""]

    color_values <- scales::hue_pal(l = 68, c = 100)(length(primary_organs))
    primary_organ_color_map <- stats::setNames(color_values, primary_organs)

    fixed_color_map <- c(
        "LUNG" = "#2F6DB2",
        "HEAD/NECK" = "#F4A259",
        "KIDNEY" = "#39A74A",
        "SKIN" = "#6EC1E4",
        "LIVER" = "#D44C8C",
        "ESOPHAGUS/STOMACH" = "#D6A23D",
        "BOWEL" = "#2E9C8E",
        "BREAST" = "#63C65E"
    )

    for (current_organ in names(fixed_color_map)) {
        current_index <- which(toupper(primary_organs) == current_organ)

        if (length(current_index) > 0) {
            primary_organ_color_map[current_index] <- fixed_color_map[[current_organ]]
        }
    }

    primary_organ_color_map
}

## 对热图行列执行层次聚类，输入details结果表和行列字段名，输出热图行列顺序与聚类对象 ####
cluster_heatmap_axes <- function(details, row_var = "cancer_type", col_var = "irae_display", value_var = "ROR") {
    clustering_df <- details %>%
        transmute(
            row_name = as.character(.data[[row_var]]),
            col_name = as.character(.data[[col_var]]),
            cluster_value = dplyr::case_when(
                is.na(.data[[value_var]]) | !is.finite(.data[[value_var]]) | .data[[value_var]] <= 0 ~ NA_real_,
                TRUE ~ log2(.data[[value_var]])
            )
        ) %>%
        distinct(row_name, col_name, .keep_all = TRUE)

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

## 构建3.1 FAERS参考热图顺序，输入参考summary和M0文件路径及阈值，输出参考癌种与SOC顺序 ####
build_reference_heatmap_order <- function(reference_summary_file, reference_m0_model_file, min_cases_threshold) {
    if (!file.exists(reference_summary_file)) {
        stop(paste("3.1参考summary文件不存在，请先准备3.1结果:", reference_summary_file))
    }

    if (!file.exists(reference_m0_model_file)) {
        stop(paste("3.1参考M0文件不存在，请先准备3.1结果:", reference_m0_model_file))
    }

    organ_count_info <- load_heatmap_organ_counts(reference_summary_file)
    reference_organs <- organ_count_info$organ_counts %>%
        filter(
            !is.na(.data[[organ_count_info$count_column]]),
            .data[[organ_count_info$count_column]] >= min_cases_threshold
        ) %>%
        pull(primary_organ) %>%
        unique()

    if (length(reference_organs) == 0) {
        stop("3.1参考热图无可用于聚类的癌种")
    }

    reference_details <- load_m0_model_details(reference_m0_model_file)
    reference_irae_display <- get_detail_irae_display(reference_details)
    reference_details_prepared <- reference_details %>%
        mutate(
            cancer_type = as.character(cancer_type),
            irae_display = as.character(reference_irae_display)
        ) %>%
        filter(cancer_type %in% reference_organs)

    if (nrow(reference_details_prepared) == 0) {
        stop("3.1参考M0结果中未找到可用于构建热图顺序的数据")
    }

    heatmap_order <- cluster_heatmap_axes(
        details = reference_details_prepared,
        row_var = "cancer_type",
        col_var = "irae_display",
        value_var = "ROR"
    )

    list(
        row_order = heatmap_order$row_levels,
        col_order = heatmap_order$col_levels,
        row_hclust = heatmap_order$row_hclust,
        col_hclust = heatmap_order$col_hclust,
        count_column = organ_count_info$count_column,
        min_cases_threshold = min_cases_threshold
    )
}

## 构建单个regimen热图数据，输入details、organ_status、regimen信息和参考顺序，输出热图矩阵与标签 ####
build_regimen_heatmap_plot_data <- function(
    details,
    organ_status,
    target_regimen_display,
    reference_heatmap_order,
    blank_threshold = 100L
) {
    regimen_details <- details %>%
        filter(ici_regimen_display == target_regimen_display) %>%
        mutate(
            primary_organ = as.character(primary_organ),
            irae_abbrev = as.character(irae_abbrev)
        ) %>%
        distinct(primary_organ, irae_abbrev, .keep_all = TRUE)

    regimen_organ_status <- organ_status %>%
        filter(ici_regimen_display == target_regimen_display) %>%
        mutate(
            primary_organ = as.character(primary_organ),
            organ_n = suppressWarnings(as.numeric(organ_n))
        ) %>%
        distinct(primary_organ, .keep_all = TRUE)

    if (nrow(regimen_details) == 0 && nrow(regimen_organ_status) == 0) {
        stop(paste("当前regimen缺少可绘图结果:", target_regimen_display))
    }

    row_meta <- tibble::tibble(primary_organ = reference_heatmap_order$row_order) %>%
        left_join(
            regimen_organ_status %>% select(primary_organ, organ_n),
            by = "primary_organ"
        ) %>%
        mutate(
            organ_n = tidyr::replace_na(organ_n, 0),
            blank_row = organ_n <= blank_threshold,
            row_label = paste0(primary_organ, " (N=", format_case_n(organ_n), ")"),
            annotation_key = if_else(blank_row, "masked_row", primary_organ),
            comparison_id = build_comparison_id(target_regimen_display, primary_organ)
        )

    log2_ror_matrix <- matrix(
        NA_real_,
        nrow = nrow(row_meta),
        ncol = length(reference_heatmap_order$col_order),
        dimnames = list(row_meta$primary_organ, reference_heatmap_order$col_order)
    )
    label_matrix <- matrix(
        "",
        nrow = nrow(row_meta),
        ncol = length(reference_heatmap_order$col_order),
        dimnames = list(row_meta$primary_organ, reference_heatmap_order$col_order)
    )
    label_color_matrix <- matrix(
        NA_character_,
        nrow = nrow(row_meta),
        ncol = length(reference_heatmap_order$col_order),
        dimnames = list(row_meta$primary_organ, reference_heatmap_order$col_order)
    )

    blank_row_lookup <- stats::setNames(row_meta$blank_row, row_meta$primary_organ)

    for (row_index in seq_len(nrow(regimen_details))) {
        current_row <- regimen_details[row_index, ]

        if (!current_row$primary_organ %in% row_meta$primary_organ) {
            next
        }

        if (!current_row$irae_abbrev %in% reference_heatmap_order$col_order) {
            next
        }

        if (isTRUE(unname(blank_row_lookup[current_row$primary_organ]))) {
            next
        }

        current_effect <- current_row$ROR
        log2_ror_matrix[current_row$primary_organ, current_row$irae_abbrev] <- dplyr::case_when(
            is.na(current_effect) | !is.finite(current_effect) | current_effect <= 0 ~ NA_real_,
            TRUE ~ pmax(pmin(log2(current_effect), 2), -2)
        )
        label_matrix[current_row$primary_organ, current_row$irae_abbrev] <- dplyr::case_when(
            is.na(current_effect) | !is.finite(current_effect) ~ "",
            isTRUE(current_row$signal_detected_within_regimen) ~ sprintf("%.2f*", current_effect),
            TRUE ~ sprintf("%.2f", current_effect)
        )
        label_color_matrix[current_row$primary_organ, current_row$irae_abbrev] <- dplyr::case_when(
            is.na(current_effect) | !is.finite(current_effect) ~ NA_character_,
            isTRUE(current_row$signal_detected_within_regimen) ~ "#1F1F1F",
            TRUE ~ "#9A9A9A"
        )
    }

    list(
        regimen_display = target_regimen_display,
        row_meta = row_meta,
        row_labels = stats::setNames(row_meta$row_label, row_meta$primary_organ),
        blank_row_flags = row_meta$blank_row,
        row_hclust = reference_heatmap_order$row_hclust,
        col_hclust = reference_heatmap_order$col_hclust,
        log2_ror_matrix = log2_ror_matrix[row_meta$primary_organ, reference_heatmap_order$col_order, drop = FALSE],
        label_matrix = label_matrix[row_meta$primary_organ, reference_heatmap_order$col_order, drop = FALSE],
        label_color_matrix = label_color_matrix[row_meta$primary_organ, reference_heatmap_order$col_order, drop = FALSE]
    )
}

## 构建单个regimen热图对象，输入plot_data、颜色映射和图标题，输出ComplexHeatmap对象 ####
build_regimen_heatmap_object <- function(plot_data, primary_organ_color_map, figure_title) {
    annotation_color_map <- c(primary_organ_color_map, "masked_row" = "#AFAFAF")

    ComplexHeatmap::Heatmap(
        plot_data$log2_ror_matrix,
        name = "M0 log2(ROR)",
        col = circlize::colorRamp2(
            c(-2, 0, 2),
            c("#8FB6D9", "#F6F1E8", "#E5A089")
        ),
        na_col = "#AFAFAF",
        rect_gp = grid::gpar(col = "black", lwd = 0.7),
        cluster_rows = if (is.null(plot_data$row_hclust)) FALSE else as.dendrogram(plot_data$row_hclust),
        cluster_columns = if (is.null(plot_data$col_hclust)) FALSE else as.dendrogram(plot_data$col_hclust),
        show_row_dend = !is.null(plot_data$row_hclust),
        show_column_dend = !is.null(plot_data$col_hclust),
        row_dend_side = "left",
        column_dend_side = "top",
        row_dend_width = grid::unit(22, "mm"),
        column_dend_height = grid::unit(22, "mm"),
        row_labels = plot_data$row_labels,
        row_names_gp = grid::gpar(fontsize = 12, col = "black"),
        column_names_gp = grid::gpar(fontsize = 12, col = "black"),
        column_names_rot = 45,
        left_annotation = ComplexHeatmap::rowAnnotation(
            primary_organ = ComplexHeatmap::anno_simple(
                x = plot_data$row_meta$annotation_key,
                col = annotation_color_map,
                border = TRUE,
                gp = grid::gpar(col = "black", lwd = 0.5)
            ),
            show_legend = FALSE,
            annotation_name_side = "top",
            annotation_name_rot = 0,
            annotation_name_gp = grid::gpar(fontsize = 14, col = "black")
        ),
        column_title = figure_title,
        column_title_gp = grid::gpar(fontsize = 14, col = "black", fontface = "bold"),
        heatmap_legend_param = list(
            at = c(-2, -1, 0, 1, 2),
            direction = "horizontal",
            title_position = "leftcenter",
            title_gp = grid::gpar(fontsize = 14, col = "black"),
            labels_gp = grid::gpar(fontsize = 12, col = "black"),
            legend_width = grid::unit(42, "mm")
        ),
        cell_fun = function(j, i, x, y, width, height, fill) {
            if (isTRUE(plot_data$blank_row_flags[i])) {
                return(invisible(NULL))
            }

            current_label <- plot_data$label_matrix[i, j]
            current_label_color <- plot_data$label_color_matrix[i, j]

            if (!is.na(current_label) && nzchar(current_label)) {
                grid::grid.text(
                    label = current_label,
                    x = x,
                    y = y,
                    gp = grid::gpar(fontsize = 9, fontface = "bold", col = current_label_color)
                )
            }
        }
    )
}

## 返回一致性热图状态标签，输入为空，输出命名字符向量 ####
get_concordance_status_labels <- function() {
    c(
        concordant = "Concordant",
        discordant = "Discordant",
        other = "Non-significant or Unavailable"
    )
}

## 返回一致性热图状态配色，输入为空，输出命名颜色向量 ####
get_concordance_status_color_map <- function() {
    c(
        "Concordant" = "#8CB4C5",
        "Discordant" = "#911F23",
        "Non-significant or Unavailable" = "#253056"
    )
}

## 构建一致性热图共享图例，输入font_size字号，输出ComplexHeatmap图例对象 ####
build_concordance_status_legend <- function(font_size = 7) {
    status_color_map <- get_concordance_status_color_map()

    ComplexHeatmap::Legend(
        labels = names(status_color_map),
        title = "Concordance with FAERS primary analysis",
        type = "grid",
        direction = "horizontal",
        nrow = 1,
        title_position = "topcenter",
        grid_width = grid::unit(3, "mm"),
        grid_height = grid::unit(3, "mm"),
        legend_gp = grid::gpar(fill = unname(status_color_map), col = "black"),
        title_gp = grid::gpar(fontsize = font_size, col = "black"),
        labels_gp = grid::gpar(fontsize = font_size, col = "black")
    )
}

## 判断当前regimen结果与FAERS主分析方向是否一致，输入reference_effect与current_effect，输出concordant、discordant或other ####
classify_concordance_status <- function(reference_effect, current_effect) {
    if (
        is.na(reference_effect) || !is.finite(reference_effect) ||
        is.na(current_effect) || !is.finite(current_effect) ||
        reference_effect <= 0 || current_effect <= 0 ||
        reference_effect == 1 || current_effect == 1
    ) {
        return("other")
    }

    if (
        (reference_effect > 1 && current_effect > 1) ||
        (reference_effect < 1 && current_effect < 1)
    ) {
        return("concordant")
    }

    "discordant"
}

## 判断FAERS主分析结果是否达到严格一致性判定要求，输入reference_p_value，输出是否保留为方向比较 ####
is_reference_signal_strict <- function(reference_p_value) {
    !is.na(reference_p_value) && is.finite(reference_p_value) && reference_p_value < 0.05
}

## 构建单个regimen三分类一致性热图数据，输入3.3结果、3.1主分析结果、organ_status、regimen信息和参考顺序，输出分类矩阵与标签 ####
build_regimen_concordance_heatmap_plot_data <- function(
    details,
    reference_main_details,
    organ_status,
    target_regimen_display,
    reference_heatmap_order,
    blank_threshold = 100L
) {
    status_labels <- get_concordance_status_labels()

    regimen_details <- details %>%
        filter(ici_regimen_display == target_regimen_display) %>%
        mutate(
            primary_organ = as.character(primary_organ),
            irae = as.character(irae),
            irae_abbrev = as.character(irae_abbrev),
            ROR = suppressWarnings(as.numeric(ROR)),
            p_value = suppressWarnings(as.numeric(p_value))
        ) %>%
        distinct(primary_organ, irae, irae_abbrev, .keep_all = TRUE)

    regimen_organ_status <- organ_status %>%
        filter(ici_regimen_display == target_regimen_display) %>%
        mutate(
            primary_organ = as.character(primary_organ),
            organ_n = suppressWarnings(as.numeric(organ_n))
        ) %>%
        distinct(primary_organ, .keep_all = TRUE)

    if (nrow(regimen_details) == 0 && nrow(regimen_organ_status) == 0) {
        stop(paste("当前regimen缺少可绘图结果:", target_regimen_display))
    }

    reference_lookup <- reference_main_details %>%
        mutate(
            primary_organ = as.character(cancer_type),
            irae = as.character(irae),
            faers_main_ROR = suppressWarnings(as.numeric(ROR)),
            faers_main_p_value = suppressWarnings(as.numeric(p_value))
        ) %>%
        distinct(primary_organ, irae, .keep_all = TRUE) %>%
        select(primary_organ, irae, faers_main_ROR, faers_main_p_value)

    regimen_details_joined <- regimen_details %>%
        left_join(reference_lookup, by = c("primary_organ", "irae"))

    regimen_total_n <- regimen_organ_status %>%
        dplyr::pull(regimen_total_n) %>%
        unique() %>%
        .[!is.na(.)]

    regimen_total_n <- if (length(regimen_total_n) == 0) NA_real_ else regimen_total_n[[1]]

    row_meta <- tibble::tibble(primary_organ = reference_heatmap_order$row_order) %>%
        left_join(
            regimen_organ_status %>% select(primary_organ, organ_n),
            by = "primary_organ"
        ) %>%
        mutate(
            organ_n = tidyr::replace_na(organ_n, 0),
            blank_row = organ_n <= blank_threshold,
            row_label = primary_organ,
            annotation_key = if_else(blank_row, "masked_row", primary_organ),
            comparison_id = build_comparison_id(target_regimen_display, primary_organ)
        )

    status_matrix <- matrix(
        status_labels[["other"]],
        nrow = nrow(row_meta),
        ncol = length(reference_heatmap_order$col_order),
        dimnames = list(row_meta$primary_organ, reference_heatmap_order$col_order)
    )

    blank_row_lookup <- stats::setNames(row_meta$blank_row, row_meta$primary_organ)

    for (row_index in seq_len(nrow(regimen_details_joined))) {
        current_row <- regimen_details_joined[row_index, ]

        if (!current_row$primary_organ %in% row_meta$primary_organ) {
            next
        }

        if (!current_row$irae_abbrev %in% reference_heatmap_order$col_order) {
            next
        }

        if (isTRUE(unname(blank_row_lookup[current_row$primary_organ]))) {
            next
        }

        if (is.na(current_row$p_value) || current_row$p_value >= 0.05) {
            next
        }

        if (!is_reference_signal_strict(current_row$faers_main_p_value)) {
            next
        }

        current_status <- classify_concordance_status(
            reference_effect = current_row$faers_main_ROR,
            current_effect = current_row$ROR
        )

        status_matrix[current_row$primary_organ, current_row$irae_abbrev] <- status_labels[[current_status]]
    }

    list(
        regimen_display = target_regimen_display,
        group_total_n = regimen_total_n,
        row_meta = row_meta,
        row_labels = stats::setNames(row_meta$row_label, row_meta$primary_organ),
        blank_row_flags = row_meta$blank_row,
        row_order = reference_heatmap_order$row_order,
        col_order = reference_heatmap_order$col_order,
        row_hclust = reference_heatmap_order$row_hclust,
        col_hclust = reference_heatmap_order$col_hclust,
        status_matrix = status_matrix[row_meta$primary_organ, reference_heatmap_order$col_order, drop = FALSE]
    )
}

## 构建单个regimen三分类一致性热图对象，输入plot_data、颜色映射、图标题、行名显示、热图名、字号和列名角度，输出ComplexHeatmap对象 ####
build_regimen_concordance_heatmap_object <- function(
    plot_data,
    primary_organ_color_map,
    figure_title,
    show_row_names = TRUE,
    heatmap_name = "signal_status",
    font_size = 7,
    column_names_rot = 90
) {
    status_levels <- unname(get_concordance_status_labels())
    status_color_map <- get_concordance_status_color_map()
    status_matrix_factor <- matrix(
        factor(plot_data$status_matrix, levels = status_levels),
        nrow = nrow(plot_data$status_matrix),
        ncol = ncol(plot_data$status_matrix),
        dimnames = dimnames(plot_data$status_matrix)
    )

    ComplexHeatmap::Heatmap(
        status_matrix_factor,
        name = heatmap_name,
        col = status_color_map,
        show_heatmap_legend = FALSE,
        rect_gp = grid::gpar(col = "black", lwd = 0.7),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_order = match(plot_data$row_order, rownames(plot_data$status_matrix)),
        column_order = match(plot_data$col_order, colnames(plot_data$status_matrix)),
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        row_dend_side = "left",
        column_dend_side = "top",
        show_row_names = show_row_names,
        row_labels = plot_data$row_labels,
        row_names_side = "left",
        row_names_gp = grid::gpar(fontsize = font_size, col = "black"),
        column_names_gp = grid::gpar(fontsize = font_size, col = "black"),
        column_names_rot = column_names_rot,
        column_title = figure_title,
        column_title_gp = grid::gpar(fontsize = font_size, col = "black")
    )
}

## 横向合并多个ComplexHeatmap对象，输入heatmap_objects列表，输出合并后的HeatmapList对象 ####
combine_heatmaps_horizontally <- function(heatmap_objects) {
    if (length(heatmap_objects) == 0) {
        stop("缺少可合并的热图对象")
    }

    combined_heatmap <- heatmap_objects[[1]]

    if (length(heatmap_objects) > 1) {
        for (heatmap_index in 2:length(heatmap_objects)) {
            combined_heatmap <- combined_heatmap + heatmap_objects[[heatmap_index]]
        }
    }

    combined_heatmap
}

## 提取ComplexHeatmap对象实际绘制后的行列顺序，输入heatmap_object、row_names和col_names，输出显示顺序列表 ####
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

## 将ComplexHeatmap对象记录为可复用grob，输入heatmap_object对象，输出可绘制grob ####
capture_complex_heatmap_grob <- function(heatmap_object) {
    grid::recordGrob(
        {
            ComplexHeatmap::draw(
                heatmap_object,
                heatmap_legend_side = "bottom",
                annotation_legend_side = "bottom",
                padding = grid::unit(c(2, 4, 2, 4), "mm"),
                newpage = FALSE
            )
        },
        list(heatmap_object = heatmap_object)
    )
}

## 绘制合并一致性热图和共享图例，输入combined_heatmap、legend_list、宽度和图例高度，输出绘制到当前设备 ####
draw_combined_concordance_heatmap <- function(
    combined_heatmap,
    legend_list,
    width_cm,
    legend_area_height_mm = 14
) {
    legend_area_height <- grid::unit(legend_area_height_mm, "mm")
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
            legend_list,
            list(
                direction = "horizontal",
                column_gap = grid::unit(3, "mm")
            )
        )
    )

    grid::pushViewport(grid::viewport(layout = page_layout))
    grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))

    ComplexHeatmap::draw(
        combined_heatmap,
        newpage = FALSE,
        padding = grid::unit(c(1, 2, 1, 2), "mm")
    )

    grid::upViewport()
    grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
    legend_box_width_mm <- max(40, width_cm * 10 - 4)
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
}

## 将合并一致性热图记录为可复用grob，输入combined_heatmap、legend_list和宽度，输出可绘制grob ####
capture_combined_concordance_heatmap_grob <- function(
    combined_heatmap,
    legend_list,
    width_cm,
    legend_area_height_mm = 14
) {
    grid::recordGrob(
        {
            draw_combined_concordance_heatmap(
                combined_heatmap = combined_heatmap,
                legend_list = legend_list,
                width_cm = width_cm,
                legend_area_height_mm = legend_area_height_mm
            )
        },
        list(
            combined_heatmap = combined_heatmap,
            legend_list = legend_list,
            width_cm = width_cm,
            legend_area_height_mm = legend_area_height_mm
        )
    )
}

## 将grob对象保存为PDF，输入plot_grob、output_file、width和height，输出写入磁盘的PDF ####
save_grid_grob_pdf <- function(plot_grob, output_file, width, height) {
    ensure_parent_dir(output_file)
    grDevices::cairo_pdf(output_file, width = width, height = height)
    grid::grid.newpage()
    grid::grid.draw(plot_grob)
    invisible(grDevices::dev.off())
}


# 1图：FAERS anti-PD-L1 regimen下primary_organ SOC热图 ####

## 1. 读取并处理绘图数据 ####
reference_heatmap_order_fig1 <- build_reference_heatmap_order(
    reference_summary_file = reference_summary_file,
    reference_m0_model_file = reference_m0_model_file,
    min_cases_threshold = reference_heatmap_min_cases
)
summary_tables_fig1 <- load_summary_tables(summary_file)
m0_details_fig1 <- load_m0_model_details(m0_model_file)
primary_organ_color_map_fig1 <- build_primary_organ_color_map(reference_heatmap_order_fig1$row_order)

plot_data_fig1 <- build_regimen_heatmap_plot_data(
    details = m0_details_fig1,
    organ_status = summary_tables_fig1$organ_status,
    target_regimen_display = format_ici_regimen_label("anti-PD-L1"),
    reference_heatmap_order = reference_heatmap_order_fig1,
    blank_threshold = reference_heatmap_min_cases
)

## 2. 绘图 ####
heatmap_object_fig1 <- build_regimen_heatmap_object(
    plot_data = plot_data_fig1,
    primary_organ_color_map = primary_organ_color_map_fig1,
    figure_title = "FAERS primary organ SOC heatmap within anti–PD-L1 regimen"
)
display_order_fig1 <- extract_complex_heatmap_display_order(
    heatmap_object = heatmap_object_fig1,
    row_names = rownames(plot_data_fig1$log2_ror_matrix),
    col_names = colnames(plot_data_fig1$log2_ror_matrix)
)
plot_grob_fig1 <- capture_complex_heatmap_grob(heatmap_object_fig1)
grid::grid.newpage()
grid::grid.draw(plot_grob_fig1)

## 3. 保存图形 ####
save_grid_grob_pdf(
    plot_grob = plot_grob_fig1,
    output_file = file.path(
        script_results_path,
        paste0(
            "图1_",
            dataset_name,
            "_",
            sanitize_file_label(plot_data_fig1$regimen_display),
            "_primary_organ_SOC热图.pdf"
        )
    ),
    width = 10,
    height = 12
)


# 2图：FAERS anti-PD-1 regimen下primary_organ SOC热图 ####

## 1. 读取并处理绘图数据 ####
reference_heatmap_order_fig2 <- build_reference_heatmap_order(
    reference_summary_file = reference_summary_file,
    reference_m0_model_file = reference_m0_model_file,
    min_cases_threshold = reference_heatmap_min_cases
)
summary_tables_fig2 <- load_summary_tables(summary_file)
m0_details_fig2 <- load_m0_model_details(m0_model_file)
primary_organ_color_map_fig2 <- build_primary_organ_color_map(reference_heatmap_order_fig2$row_order)

plot_data_fig2 <- build_regimen_heatmap_plot_data(
    details = m0_details_fig2,
    organ_status = summary_tables_fig2$organ_status,
    target_regimen_display = format_ici_regimen_label("anti-PD-1"),
    reference_heatmap_order = reference_heatmap_order_fig2,
    blank_threshold = reference_heatmap_min_cases
)

## 2. 绘图 ####
heatmap_object_fig2 <- build_regimen_heatmap_object(
    plot_data = plot_data_fig2,
    primary_organ_color_map = primary_organ_color_map_fig2,
    figure_title = "FAERS primary organ SOC heatmap within anti–PD-1 regimen"
)
display_order_fig2 <- extract_complex_heatmap_display_order(
    heatmap_object = heatmap_object_fig2,
    row_names = rownames(plot_data_fig2$log2_ror_matrix),
    col_names = colnames(plot_data_fig2$log2_ror_matrix)
)
plot_grob_fig2 <- capture_complex_heatmap_grob(heatmap_object_fig2)
grid::grid.newpage()
grid::grid.draw(plot_grob_fig2)

## 3. 保存图形 ####
save_grid_grob_pdf(
    plot_grob = plot_grob_fig2,
    output_file = file.path(
        script_results_path,
        paste0(
            "图2_",
            dataset_name,
            "_",
            sanitize_file_label(plot_data_fig2$regimen_display),
            "_primary_organ_SOC热图.pdf"
        )
    ),
    width = 10,
    height = 12
)


# 3图：FAERS anti-CTLA-4 + anti-PD-1 regimen下primary_organ SOC热图 ####

## 1. 读取并处理绘图数据 ####
reference_heatmap_order_fig3 <- build_reference_heatmap_order(
    reference_summary_file = reference_summary_file,
    reference_m0_model_file = reference_m0_model_file,
    min_cases_threshold = reference_heatmap_min_cases
)
summary_tables_fig3 <- load_summary_tables(summary_file)
m0_details_fig3 <- load_m0_model_details(m0_model_file)
primary_organ_color_map_fig3 <- build_primary_organ_color_map(reference_heatmap_order_fig3$row_order)

plot_data_fig3 <- build_regimen_heatmap_plot_data(
    details = m0_details_fig3,
    organ_status = summary_tables_fig3$organ_status,
    target_regimen_display = format_ici_regimen_label("anti-CTLA-4+anti-PD-1"),
    reference_heatmap_order = reference_heatmap_order_fig3,
    blank_threshold = reference_heatmap_min_cases
)

## 2. 绘图 ####
heatmap_object_fig3 <- build_regimen_heatmap_object(
    plot_data = plot_data_fig3,
    primary_organ_color_map = primary_organ_color_map_fig3,
    figure_title = "FAERS primary organ SOC heatmap within anti–CTLA-4 + anti–PD-1 regimen"
)
display_order_fig3 <- extract_complex_heatmap_display_order(
    heatmap_object = heatmap_object_fig3,
    row_names = rownames(plot_data_fig3$log2_ror_matrix),
    col_names = colnames(plot_data_fig3$log2_ror_matrix)
)
plot_grob_fig3 <- capture_complex_heatmap_grob(heatmap_object_fig3)
grid::grid.newpage()
grid::grid.draw(plot_grob_fig3)

## 3. 保存图形 ####
save_grid_grob_pdf(
    plot_grob = plot_grob_fig3,
    output_file = file.path(
        script_results_path,
        paste0(
            "图3_",
            dataset_name,
            "_",
            sanitize_file_label(plot_data_fig3$regimen_display),
            "_primary_organ_SOC热图.pdf"
        )
    ),
    width = 10,
    height = 12
)


# 4-6图：FAERS不同ICI regimen下primary_organ SOC一致性合并热图 ####

## 1. 读取并处理绘图数据 ####
combined_concordance_row_order <- display_order_fig1$row_order
combined_concordance_col_order <- display_order_fig1$col_order

reference_heatmap_order_fig4 <- build_reference_heatmap_order(
    reference_summary_file = reference_summary_file,
    reference_m0_model_file = reference_m0_model_file,
    min_cases_threshold = reference_heatmap_min_cases
)
summary_tables_fig4 <- load_summary_tables(summary_file)
m0_details_fig4 <- load_m0_model_details(m0_model_file)
reference_main_details_fig4 <- load_m0_model_details(reference_m0_model_file)
primary_organ_color_map_fig4 <- build_primary_organ_color_map(reference_heatmap_order_fig4$row_order)

plot_data_fig4 <- build_regimen_concordance_heatmap_plot_data(
    details = m0_details_fig4,
    reference_main_details = reference_main_details_fig4,
    organ_status = summary_tables_fig4$organ_status,
    target_regimen_display = format_ici_regimen_label("anti-PD-L1"),
    reference_heatmap_order = reference_heatmap_order_fig4,
    blank_threshold = reference_heatmap_min_cases
)
plot_data_fig4$row_order <- combined_concordance_row_order
plot_data_fig4$col_order <- combined_concordance_col_order
plot_data_fig4$row_hclust <- NULL
plot_data_fig4$col_hclust <- NULL

reference_heatmap_order_fig5 <- build_reference_heatmap_order(
    reference_summary_file = reference_summary_file,
    reference_m0_model_file = reference_m0_model_file,
    min_cases_threshold = reference_heatmap_min_cases
)
summary_tables_fig5 <- load_summary_tables(summary_file)
m0_details_fig5 <- load_m0_model_details(m0_model_file)
reference_main_details_fig5 <- load_m0_model_details(reference_m0_model_file)
primary_organ_color_map_fig5 <- build_primary_organ_color_map(reference_heatmap_order_fig5$row_order)

plot_data_fig5 <- build_regimen_concordance_heatmap_plot_data(
    details = m0_details_fig5,
    reference_main_details = reference_main_details_fig5,
    organ_status = summary_tables_fig5$organ_status,
    target_regimen_display = format_ici_regimen_label("anti-PD-1"),
    reference_heatmap_order = reference_heatmap_order_fig5,
    blank_threshold = reference_heatmap_min_cases
)
plot_data_fig5$row_order <- combined_concordance_row_order
plot_data_fig5$col_order <- combined_concordance_col_order
plot_data_fig5$row_hclust <- NULL
plot_data_fig5$col_hclust <- NULL

reference_heatmap_order_fig6 <- build_reference_heatmap_order(
    reference_summary_file = reference_summary_file,
    reference_m0_model_file = reference_m0_model_file,
    min_cases_threshold = reference_heatmap_min_cases
)
summary_tables_fig6 <- load_summary_tables(summary_file)
m0_details_fig6 <- load_m0_model_details(m0_model_file)
reference_main_details_fig6 <- load_m0_model_details(reference_m0_model_file)
primary_organ_color_map_fig6 <- build_primary_organ_color_map(reference_heatmap_order_fig6$row_order)

plot_data_fig6 <- build_regimen_concordance_heatmap_plot_data(
    details = m0_details_fig6,
    reference_main_details = reference_main_details_fig6,
    organ_status = summary_tables_fig6$organ_status,
    target_regimen_display = format_ici_regimen_label("anti-CTLA-4+anti-PD-1"),
    reference_heatmap_order = reference_heatmap_order_fig6,
    blank_threshold = reference_heatmap_min_cases
)
plot_data_fig6$row_order <- combined_concordance_row_order
plot_data_fig6$col_order <- combined_concordance_col_order
plot_data_fig6$row_hclust <- NULL
plot_data_fig6$col_hclust <- NULL

## 2. 绘图 ####
heatmap_object_fig4 <- build_regimen_concordance_heatmap_object(
    plot_data = plot_data_fig4,
    primary_organ_color_map = primary_organ_color_map_fig4,
    figure_title = build_heatmap_title_with_n(
        display_label = plot_data_fig4$regimen_display,
        total_n = plot_data_fig4$group_total_n
    ),
    show_row_names = TRUE,
    heatmap_name = "signal_status_fig4",
    font_size = 7,
    column_names_rot = 90
)
heatmap_object_fig5 <- build_regimen_concordance_heatmap_object(
    plot_data = plot_data_fig5,
    primary_organ_color_map = primary_organ_color_map_fig5,
    figure_title = build_heatmap_title_with_n(
        display_label = plot_data_fig5$regimen_display,
        total_n = plot_data_fig5$group_total_n
    ),
    show_row_names = FALSE,
    heatmap_name = "signal_status_fig5",
    font_size = 7,
    column_names_rot = 90
)
heatmap_object_fig6 <- build_regimen_concordance_heatmap_object(
    plot_data = plot_data_fig6,
    primary_organ_color_map = primary_organ_color_map_fig6,
    figure_title = build_heatmap_title_with_n(
        display_label = plot_data_fig6$regimen_display,
        total_n = plot_data_fig6$group_total_n
    ),
    show_row_names = FALSE,
    heatmap_name = "signal_status_fig6",
    font_size = 7,
    column_names_rot = 90
)
combined_heatmap_fig4_6 <- combine_heatmaps_horizontally(
    list(
        heatmap_object_fig4,
        heatmap_object_fig5,
        heatmap_object_fig6
    )
)
concordance_legend_fig4_6 <- build_concordance_status_legend(font_size = 7)
plot_grob_fig4_6 <- capture_combined_concordance_heatmap_grob(
    combined_heatmap = combined_heatmap_fig4_6,
    legend_list = list(concordance_legend_fig4_6),
    width_cm = 14.65,
    legend_area_height_mm = 14
)
grid::grid.newpage()
grid::grid.draw(plot_grob_fig4_6)

## 3. 保存图形 ####
save_grid_grob_pdf(
    plot_grob = plot_grob_fig4_6,
    output_file = file.path(
        script_results_path,
        paste0("图4-6_", dataset_name, "_ICI_regimen_primary_organ_SOC一致性合并热图.pdf")
    ),
    width = 14.65 / 2.54,
    height = 10 / 2.54
)

message("\n========== 3.3.1 绘图完成 ==========")
