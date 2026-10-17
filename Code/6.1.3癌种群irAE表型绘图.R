# 癌种群 irAE 表型绘图（FAERS + JADER）####

rm(list = ls())
gc()

# 1. 加载配置 ####
source("Code/config.R")


# 2. 结果目录与输入文件 ####
analysis_script_name <- "6.1.2癌种群irAE表型分析"
analysis_results_path <- file.path(results_path, analysis_script_name)

script_name <- "6.1.3癌种群irAE表型绘图"
script_results_path <- file.path(results_path, script_name)

if (!dir.exists(script_results_path)) {
    dir.create(script_results_path, recursive = TRUE, showWarnings = FALSE)
}

faers_summary_file <- file.path(analysis_results_path, "FAERS_癌种群irAE表型汇总.xlsx")
jader_summary_file <- file.path(analysis_results_path, "JADER_癌种群irAE表型汇总.xlsx")

faers_soc_plot_file <- file.path(script_results_path, "图1_FAERS_癌种群SOC谱热图.pdf")
faers_tto_plot_file <- file.path(script_results_path, "图2_FAERS_癌种群TTO比较图.pdf")
faers_fatality_plot_file <- file.path(script_results_path, "图3_FAERS_癌种群fatality柱状图.pdf")
jader_soc_plot_file <- file.path(script_results_path, "图4_JADER_癌种群SOC谱热图.pdf")
jader_tto_plot_file <- file.path(script_results_path, "图5_JADER_癌种群TTO比较图.pdf")
jader_fatality_plot_file <- file.path(script_results_path, "图6_JADER_癌种群fatality柱状图.pdf")
faers_weibull_plot_file <- file.path(script_results_path, "图7_FAERS_癌种群TTO_Weibull分析.pdf")
jader_weibull_plot_file <- file.path(script_results_path, "图8_JADER_癌种群TTO_Weibull分析.pdf")

pattern_display_order <- c(
    "Ocular",
    "Hematologic",
    "Neurologic",
    "Systemic-like",
    "Visceral-like"
)

pattern_color_map <- c(
    "Ocular" = "#DF7A5E",
    "Hematologic" = "#82B29A",
    "Neurologic" = "#A15F88",
    "Systemic-like" = "#F2CC8E",
    "Visceral-like" = "#BC4749"
)

message("========== 癌种群 irAE 表型绘图开始 ==========")
message("分析结果目录: ", analysis_results_path)
message("绘图结果目录: ", script_results_path)


# 3. 辅助函数 ####

## 函数：format_case_n ####
## 功能：格式化样本量为带千分位分隔符的字符值 ####
## 输入：x，数值型样本量 ####
## 输出：返回格式化后的字符值 ####
format_case_n <- function(x) {
    format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
}

## 函数：format_p_value_text ####
## 功能：格式化P值文本，便于图中展示 ####
## 输入：p_value，数值型P值 ####
## 输出：返回格式化后的P值字符值 ####
format_p_value_text <- function(p_value) {
    if (length(p_value) == 0 || is.na(p_value)) {
        return("p = NA")
    }

    if (p_value < 0.001) {
        return("p < 0.001")
    }

    paste0("p = ", sprintf("%.3f", p_value))
}

## 函数：format_overall_p_value_plotmath ####
## 功能：格式化整体检验P值为plotmath表达式，使p显示为斜体 ####
## 输入：p_value，数值型P值 ####
## 输出：返回可被ggplot2 parse的plotmath字符值 ####
format_overall_p_value_plotmath <- function(p_value) {
    if (length(p_value) == 0 || is.na(p_value)) {
        return("'Overall'~italic(p)~'='~'NA'")
    }

    if (p_value < 0.001) {
        return("'Overall'~italic(p)~'<'~0.001")
    }

    paste0("'Overall'~italic(p)~'='~", sprintf("%.3f", p_value))
}

## 函数：format_weibull_failure_display ####
## 功能：将Weibull风险类型转换为图中展示的失效类型 ####
## 输入：hazard_type，风险类型字符值 ####
## 输出：返回图中展示的失效类型字符值 ####
format_weibull_failure_display <- function(hazard_type) {
    dplyr::case_when(
        is.na(hazard_type) ~ "Unknown failure",
        hazard_type == "Decreasing hazard" ~ "Early failure",
        hazard_type == "Approximately constant hazard" ~ "Random failure",
        hazard_type == "Increasing hazard" ~ "Wear-out failure",
        TRUE ~ as.character(hazard_type)
    )
}

## 函数：build_weibull_display_start_df ####
## 功能：为Weibull拟合曲线计算可视化起点，避免shape < 1时0附近尖峰压缩直方图 ####
## 输入：case_level_df，TTO病例级数据框；summary_df，Weibull汇总结果数据框 ####
## 输出：返回包含pattern_group和display_x_min的数据框 ####
build_weibull_display_start_df <- function(case_level_df, summary_df) {
    positive_min_df <- case_level_df %>%
        dplyr::filter(!is.na(ici_tto), ici_tto > 0) %>%
        dplyr::group_by(pattern_group) %>%
        dplyr::summarise(
            min_positive_tto = min(ici_tto, na.rm = TRUE),
            .groups = "drop"
        )

    summary_df %>%
        dplyr::transmute(
            pattern_group = as.character(pattern_group),
            zero_tto_n = suppressWarnings(as.numeric(zero_tto_n)),
            zero_replaced_value = suppressWarnings(as.numeric(zero_replaced_value))
        ) %>%
        dplyr::left_join(positive_min_df, by = "pattern_group") %>%
        dplyr::mutate(
            display_x_min = dplyr::case_when(
                !is.na(zero_tto_n) & zero_tto_n > 0 & !is.na(zero_replaced_value) ~ zero_replaced_value,
                !is.na(min_positive_tto) ~ pmax(min_positive_tto, 0.5),
                TRUE ~ 0.5
            )
        ) %>%
        dplyr::select(pattern_group, display_x_min)
}

## 函数：build_tto_stats_panel_text ####
## 功能：整理TTO描述统计、总体检验和显著两两比较为图内文本 ####
## 输入：group_summary_df，TTO分组汇总表；global_test_df，总体检验结果；pairwise_df，两两比较结果；pattern_order，固定分组顺序 ####
## 输出：返回用于统计信息面板的字符值 ####
build_tto_stats_panel_text <- function(group_summary_df, global_test_df, pairwise_df, pattern_order) {
    ordered_group_summary <- tibble::tibble(pattern_group = pattern_order) %>%
        dplyr::left_join(group_summary_df, by = "pattern_group") %>%
        dplyr::mutate(
            n = suppressWarnings(as.numeric(n)),
            median_tto = suppressWarnings(as.numeric(median_tto)),
            q1_tto = suppressWarnings(as.numeric(q1_tto)),
            q3_tto = suppressWarnings(as.numeric(q3_tto))
        ) %>%
        dplyr::filter(!is.na(n), n > 0)

    group_lines <- if (nrow(ordered_group_summary) > 0) {
        paste0(
            ordered_group_summary$pattern_group,
            ": ",
            sprintf("%.1f", ordered_group_summary$median_tto),
            " [",
            sprintf("%.1f", ordered_group_summary$q1_tto),
            "-",
            sprintf("%.1f", ordered_group_summary$q3_tto),
            "], n=",
            format_case_n(ordered_group_summary$n)
        )
    } else {
        "No valid TTO summary available."
    }

    global_p_value <- suppressWarnings(as.numeric(global_test_df$p_value[[1]]))
    global_line <- paste0("Kruskal-Wallis: ", format_p_value_text(global_p_value))

    significant_pairwise_df <- pairwise_df %>%
        dplyr::filter(
            !is.na(group1),
            !is.na(group2),
            !is.na(p_value),
            !is.na(significance),
            significance != "ns"
        ) %>%
        dplyr::mutate(
            group1 = as.character(group1),
            group2 = as.character(group2),
            p_value = suppressWarnings(as.numeric(p_value)),
            group1_order = match(group1, pattern_order),
            group2_order = match(group2, pattern_order)
        ) %>%
        dplyr::filter(!is.na(group1_order), !is.na(group2_order)) %>%
        dplyr::arrange(group1_order, group2_order, p_value)

    pairwise_lines <- if (nrow(significant_pairwise_df) > 0) {
        paste0(
            significant_pairwise_df$group1,
            " vs ",
            significant_pairwise_df$group2,
            ": ",
            vapply(significant_pairwise_df$p_value, format_p_value_text, character(1)),
            " ",
            significant_pairwise_df$significance
        )
    } else {
        "No significant pairwise differences."
    }

    paste(
        c(
            "Median [IQR] in days",
            group_lines,
            "",
            global_line,
            "Dunn-Bonferroni significant pairs",
            pairwise_lines
        ),
        collapse = "\n"
    )
}

## 函数：build_pattern_group_label_map ####
## 功能：为pattern group构建带样本量的展示标签 ####
## 输入：summary_df，包含pattern_group和样本量列的数据框；n_col，样本量列名；pattern_order，固定pattern顺序；separator，组名与样本量之间的连接符 ####
## 输出：返回以pattern group命名的展示标签向量 ####
build_pattern_group_label_map <- function(summary_df, n_col, pattern_order, separator = "\n") {
    summary_df <- summary_df %>%
        dplyr::mutate(
            pattern_group = as.character(pattern_group),
            n_value = suppressWarnings(as.numeric(.data[[n_col]]))
        )

    label_df <- tibble::tibble(pattern_group = pattern_order) %>%
        dplyr::left_join(
            summary_df %>% dplyr::select(pattern_group, n_value),
            by = "pattern_group"
        ) %>%
        dplyr::mutate(
            n_value = tidyr::replace_na(n_value, 0),
            display_label = paste0(pattern_group, separator, "(n=", format_case_n(n_value), ")")
        )

    stats::setNames(label_df$display_label, label_df$pattern_group)
}

## 函数：build_descending_group_order ####
## 功能：按指定数值列从高到低生成分组排序，缺失值排在最后 ####
## 输入：summary_df，包含pattern_group和目标数值列的数据框；value_col，排序数值列名；fallback_order，缺失或并列时的后备顺序 ####
## 输出：返回排序后的pattern group字符向量 ####
build_descending_group_order <- function(summary_df, value_col, fallback_order) {
    summary_df %>%
        dplyr::mutate(
            pattern_group = as.character(pattern_group),
            sort_value = suppressWarnings(as.numeric(.data[[value_col]])),
            sort_value = dplyr::coalesce(sort_value, -Inf)
        ) %>%
        dplyr::arrange(
            dplyr::desc(sort_value),
            factor(pattern_group, levels = fallback_order)
        ) %>%
        dplyr::pull(pattern_group) %>%
        unique()
}

## 函数：build_fatality_annotation_df ####
## 功能：将fatality两两比较结果转换为显著性线段所需数据 ####
## 输入：pairwise_df，fatality两两比较结果表；plot_data，fatality柱状图数据；pattern_order，固定pattern顺序；line_gap_scale，显著性线段垂直间距倍率；line_base_scale，显著性线段整体起始高度倍率 ####
## 输出：返回显著性标注数据框 ####
build_fatality_annotation_df <- function(pairwise_df, plot_data, pattern_order, line_gap_scale = 1, line_base_scale = 1) {
    if (nrow(pairwise_df) == 0) {
        return(tibble::tibble())
    }

    valid_pairwise_df <- pairwise_df %>%
        dplyr::filter(
            !is.na(group1),
            !is.na(group2),
            !is.na(significance),
            significance != "ns"
        ) %>%
        dplyr::mutate(
            group1 = as.character(group1),
            group2 = as.character(group2),
            group1_order = match(group1, pattern_order),
            group2_order = match(group2, pattern_order)
        ) %>%
        dplyr::filter(!is.na(group1_order), !is.na(group2_order)) %>%
        dplyr::arrange(group1_order, group2_order, p_value)

    if (nrow(valid_pairwise_df) == 0) {
        return(tibble::tibble())
    }

    x_order <- plot_data %>%
        dplyr::distinct(pattern_group) %>%
        dplyr::pull(pattern_group) %>%
        as.character()

    x_positions <- stats::setNames(seq_along(x_order), x_order)
    max_y <- max(plot_data$fatality_pct, na.rm = TRUE)
    if (!is.finite(max_y)) {
        max_y <- 0
    }

    base_offset <- max(5, max_y * 0.24) * line_base_scale
    y_step <- max(2, max_y * 0.08) * line_gap_scale
    tip_height <- max(0.5, y_step * 0.3)

    valid_pairwise_df %>%
        dplyr::mutate(
            xmin = pmin(x_positions[group1], x_positions[group2]),
            xmax = pmax(x_positions[group1], x_positions[group2]),
            annotation_rank = dplyr::row_number(),
            y_position = max_y + base_offset + (annotation_rank - 1) * y_step,
            label_x = (xmin + xmax) / 2,
            label_y = y_position + tip_height * 0.4,
            tip_height = tip_height
        )
}

## 函数：build_tto_boxplot_label_df ####
## 功能：整理TTO箱线图中位数和四分位数的数值标注数据 ####
## 输入：group_summary_df，TTO分组汇总表；label_map，pattern group到展示标签的映射；pattern_order，固定pattern顺序 ####
## 输出：返回包含箱线图数值标注位置和文本的数据框 ####
build_tto_boxplot_label_df <- function(group_summary_df, label_map, pattern_order) {
    detail_label_levels <- rev(unname(label_map[pattern_order]))

    tibble::tibble(pattern_group = pattern_order) %>%
        dplyr::left_join(group_summary_df, by = "pattern_group") %>%
        dplyr::mutate(
            pattern_group = as.character(pattern_group),
            q1_tto = suppressWarnings(as.numeric(q1_tto)),
            median_tto = suppressWarnings(as.numeric(median_tto)),
            q3_tto = suppressWarnings(as.numeric(q3_tto)),
            pattern_group_label_detail = factor(
                label_map[pattern_group],
                levels = detail_label_levels
            ),
            pattern_group_y = match(
                as.character(pattern_group_label_detail),
                detail_label_levels
            )
        ) %>%
        dplyr::select(pattern_group_label_detail, pattern_group_y, q1_tto, median_tto, q3_tto) %>%
        tidyr::pivot_longer(
            cols = c(q1_tto, median_tto, q3_tto),
            names_to = "stat_type",
            values_to = "stat_value"
        ) %>%
        dplyr::filter(!is.na(stat_value)) %>%
        dplyr::mutate(label = sprintf("%.0f", stat_value))
}

## 函数：build_tto_pairwise_annotation_df ####
## 功能：将TTO两两比较结果转换为箱线图右侧显著性括号所需数据 ####
## 输入：pairwise_df，TTO两两比较结果表；label_map，pattern group到展示标签的映射；pattern_order，固定pattern顺序；detail_label_levels，箱线图y轴标签顺序；plot_max_x，箱线图主数据最大x值；x_gap_scale，右侧显著性括号横向间距倍率 ####
## 输出：返回TTO显著性标注数据框 ####
build_tto_pairwise_annotation_df <- function(pairwise_df, label_map, pattern_order, detail_label_levels, plot_max_x, x_gap_scale = 1) {
    if (nrow(pairwise_df) == 0) {
        return(tibble::tibble())
    }

    if (!is.finite(plot_max_x) || plot_max_x <= 0) {
        plot_max_x <- 1
    }

    valid_pairwise_df <- pairwise_df %>%
        dplyr::filter(
            !is.na(group1),
            !is.na(group2),
            !is.na(significance),
            significance != "ns"
        ) %>%
        dplyr::mutate(
            group1 = as.character(group1),
            group2 = as.character(group2),
            p_value = suppressWarnings(as.numeric(p_value)),
            group1_order = match(group1, pattern_order),
            group2_order = match(group2, pattern_order),
            group1_label = unname(label_map[group1]),
            group2_label = unname(label_map[group2])
        ) %>%
        dplyr::filter(!is.na(group1_order), !is.na(group2_order)) %>%
        dplyr::arrange(group1_order, group2_order, p_value)

    if (nrow(valid_pairwise_df) == 0) {
        return(tibble::tibble())
    }

    y_positions <- stats::setNames(seq_along(detail_label_levels), detail_label_levels)
    base_x <- plot_max_x + max(15, plot_max_x * 0.03)
    x_step <- max(12, plot_max_x * 0.045) * x_gap_scale
    tip_width <- max(6, x_step * 0.35)
    label_offset <- max(4, x_step * 0.25)

    valid_pairwise_df %>%
        dplyr::mutate(
            y1 = y_positions[group1_label],
            y2 = y_positions[group2_label]
        ) %>%
        dplyr::filter(!is.na(y1), !is.na(y2)) %>%
        dplyr::mutate(
            ymin = pmin(y1, y2),
            ymax = pmax(y1, y2),
            annotation_rank = dplyr::row_number(),
            x_position = base_x + annotation_rank * x_step,
            label_x = x_position + label_offset,
            label_y = (ymin + ymax) / 2,
            tip_width = tip_width
        )
}


# 4. 1图：FAERS SOC reporting spectrum heatmap ####

## 1. 读取并处理绘图数据 ####
faers_soc_plot_data <- readxl::read_xlsx(faers_summary_file, sheet = "SOC_spectrum") %>%
    dplyr::mutate(
        pattern_group = as.character(pattern_group),
        soc_abbrev = as.character(soc_abbrev),
        soc_order = suppressWarnings(as.numeric(soc_order)),
        RP = suppressWarnings(as.numeric(RP))
    )

faers_pattern_group_n <- readxl::read_xlsx(faers_summary_file, sheet = "pattern_group_n") %>%
    dplyr::mutate(
        pattern_group = as.character(pattern_group),
        mapped_ici_n = suppressWarnings(as.numeric(mapped_ici_n))
    )

faers_soc_label_map <- build_pattern_group_label_map(
    summary_df = faers_pattern_group_n,
    n_col = "mapped_ici_n",
    pattern_order = pattern_display_order
)

faers_soc_levels <- faers_soc_plot_data %>%
    dplyr::arrange(soc_order, soc_abbrev) %>%
    dplyr::pull(soc_abbrev) %>%
    unique()

faers_soc_plot_data <- faers_soc_plot_data %>%
    dplyr::mutate(
        label_text = sprintf("%.1f", RP),
        label_color = ifelse(
            RP >= max(RP, na.rm = TRUE) * 0.55,
            "white",
            "black"
        ),
        pattern_group_label = factor(
            faers_soc_label_map[pattern_group],
            levels = unname(faers_soc_label_map[pattern_display_order])
        ),
        soc_abbrev = factor(
            soc_abbrev,
            levels = rev(faers_soc_levels)
        )
    )

## 2. 绘图 ####
p <- ggplot(faers_soc_plot_data, aes(x = pattern_group_label, y = soc_abbrev, fill = RP)) +
    geom_tile(color = "white", linewidth = 0.25) +
    geom_text(
        aes(label = label_text, color = label_color),
        size = 10 / ggplot2::.pt
    ) +
    scale_fill_gradient(
        low = "#F6F3EE",
        high = "#B23A48",
        name = "RP (%)"
    ) +
    guides(
        fill = guide_colorbar(
            direction = "horizontal",
            title.position = "left",
            title.theme = element_text(
                size = 14,
                color = "black",
                angle = 0,
                margin = ggplot2::margin(r = 8)
            ),
            barwidth = grid::unit(5.2, "cm"),
            barheight = grid::unit(0.35, "cm")
        )
    ) +
    scale_color_identity() +
    labs(
        x = NULL,
        y = NULL
    ) +
    theme_bw() +
    theme(
        axis.title.x = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14, color = "black"),
        axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 14, color = "black"),
        legend.text = element_text(size = 12, color = "black"),
        legend.position = "top",
        panel.grid = element_blank()
    )

print(p)

## 3. 保存图形 ####
ggsave(
    filename = faers_soc_plot_file,
    plot = p,
    width = 5,
    height = 4,
    device = cairo_pdf
)


# 5. 2图：FAERS Time to Onset cumulative incidence ####

## 1. 读取并处理绘图数据 ####
faers_tto_plot_data <- readxl::read_xlsx(faers_summary_file, sheet = "TTO_case_level") %>%
    dplyr::mutate(
        pattern_group = as.character(pattern_group),
        ici_tto = suppressWarnings(as.numeric(ici_tto))
    ) %>%
    dplyr::filter(!is.na(ici_tto))

faers_tto_summary <- readxl::read_xlsx(faers_summary_file, sheet = "TTO_summary") %>%
    dplyr::mutate(section = as.character(section))

faers_tto_group_summary <- faers_tto_summary %>%
    dplyr::filter(section == "group_summary") %>%
    dplyr::mutate(
        pattern_group = as.character(pattern_group),
        n = suppressWarnings(as.numeric(n)),
        median_tto = suppressWarnings(as.numeric(median_tto)),
        q1_tto = suppressWarnings(as.numeric(q1_tto)),
        q3_tto = suppressWarnings(as.numeric(q3_tto))
    )

faers_tto_pairwise <- faers_tto_summary %>%
    dplyr::filter(section == "pairwise_test") %>%
    dplyr::mutate(
        group1 = as.character(group1),
        group2 = as.character(group2),
        significance = as.character(significance),
        p_value = suppressWarnings(as.numeric(p_value))
    )

faers_tto_global_test <- faers_tto_summary %>%
    dplyr::filter(section == "global_test") %>%
    dplyr::slice(1)

faers_tto_label_map <- build_pattern_group_label_map(
    summary_df = faers_tto_group_summary,
    n_col = "n",
    pattern_order = pattern_display_order,
    separator = " "
)

faers_tto_plot_data <- faers_tto_plot_data %>%
    dplyr::mutate(
        pattern_group_label = factor(
            faers_tto_label_map[pattern_group],
            levels = unname(faers_tto_label_map[pattern_display_order])
        ),
        pattern_group_label_detail = factor(
            faers_tto_label_map[pattern_group],
            levels = rev(unname(faers_tto_label_map[pattern_display_order]))
        )
    )

faers_tto_color_map <- stats::setNames(
    unname(pattern_color_map[pattern_display_order]),
    unname(faers_tto_label_map[pattern_display_order])
)

faers_tto_max_x <- if (nrow(faers_tto_plot_data) > 0) {
    max(faers_tto_plot_data$ici_tto, na.rm = TRUE)
} else {
    1
}

faers_tto_detail_levels <- rev(unname(faers_tto_label_map[pattern_display_order]))
faers_tto_detail_levels_plain <- rev(unname(pattern_display_order))

faers_tto_plot_data <- faers_tto_plot_data %>%
    dplyr::mutate(
        pattern_group_y = match(
            as.character(pattern_group_label_detail),
            faers_tto_detail_levels
        )
    )

faers_tto_boxplot_labels <- build_tto_boxplot_label_df(
    group_summary_df = faers_tto_group_summary,
    label_map = faers_tto_label_map,
    pattern_order = pattern_display_order
) %>%
    dplyr::mutate(
        label_x = dplyr::case_when(
            stat_type == "q1_tto" ~ pmax(stat_value - 3, 0),
            stat_type == "q3_tto" ~ stat_value + 3,
            TRUE ~ stat_value
        ),
        label_y = pattern_group_y + dplyr::case_when(
            stat_type == "q1_tto" ~ -0.23,
            stat_type == "median_tto" ~ 0,
            stat_type == "q3_tto" ~ 0.23,
            TRUE ~ 0
        ),
        label_hjust = dplyr::case_when(
            stat_type == "q1_tto" ~ 1,
            stat_type == "q3_tto" ~ 0,
            TRUE ~ 0.5
        ),
        fontface = dplyr::if_else(stat_type == "median_tto", "bold", "plain")
    )

faers_tto_pairwise_annotations <- build_tto_pairwise_annotation_df(
    pairwise_df = faers_tto_pairwise,
    label_map = faers_tto_label_map,
    pattern_order = pattern_display_order,
    detail_label_levels = faers_tto_detail_levels,
    plot_max_x = faers_tto_max_x,
    x_gap_scale = 1.75
)

faers_tto_pairwise_annotations <- faers_tto_pairwise_annotations %>%
    dplyr::mutate(
        significance_vertical = vapply(
            strsplit(significance, ""),
            paste,
            collapse = "\n",
            FUN.VALUE = character(1)
        )
    )

faers_tto_global_p <- suppressWarnings(as.numeric(faers_tto_global_test$p_value[[1]]))
faers_tto_global_label <- format_overall_p_value_plotmath(faers_tto_global_p)

faers_tto_detail_xmax <- faers_tto_max_x
if (nrow(faers_tto_pairwise_annotations) > 0) {
    faers_tto_detail_xmax <- max(
        faers_tto_detail_xmax,
        faers_tto_pairwise_annotations$label_x
    )
}
if (!is.finite(faers_tto_detail_xmax) || faers_tto_detail_xmax <= 0) {
    faers_tto_detail_xmax <- 1
}
faers_tto_detail_xmax <- faers_tto_detail_xmax + max(2, faers_tto_max_x * 0.005)

## 2. 绘图 ####
p_detail <- ggplot(
    faers_tto_plot_data,
    aes(
        x = ici_tto,
        y = pattern_group_y,
        group = pattern_group_label_detail,
        color = pattern_group_label,
        fill = pattern_group_label
    )
) +
    geom_boxplot(
        orientation = "y",
        outlier.shape = NA,
        width = 0.6,
        alpha = 0.18,
        linewidth = 0.4,
        show.legend = FALSE
    ) +
    geom_jitter(
        height = 0.16,
        width = 0,
        alpha = 0.18,
        size = 0.7,
        show.legend = FALSE
    ) +
    geom_label(
        data = faers_tto_boxplot_labels,
        aes(
            x = label_x,
            y = label_y,
            label = label,
            hjust = label_hjust,
            fontface = fontface
        ),
        inherit.aes = FALSE,
        size = 7 / ggplot2::.pt,
        linewidth = 0.15,
        label.padding = grid::unit(0.12, "lines"),
        fill = "white",
        color = "black"
    ) +
    scale_color_manual(values = faers_tto_color_map, guide = "none") +
    scale_fill_manual(values = faers_tto_color_map, guide = "none") +
    scale_y_continuous(
        breaks = seq_along(faers_tto_detail_levels),
        labels = faers_tto_detail_levels_plain,
        expand = ggplot2::expansion(add = c(0.25, 0.75))
    ) +
    coord_cartesian(xlim = c(0, faers_tto_detail_xmax), clip = "off") +
    labs(
        x = "Time to Onset (Days)",
        y = NULL
    ) +
    annotate(
        geom = "text",
        x = faers_tto_detail_xmax * 0.82,
        y = length(faers_tto_detail_levels) + 0.45,
        label = faers_tto_global_label,
        parse = TRUE,
        size = 8 / ggplot2::.pt,
        color = "black"
    ) +
    theme_bw() +
    theme(
        axis.title.x = element_text(size = 8.5, color = "black"),
        axis.title.y = element_text(size = 8.5, color = "black"),
        axis.text.x = element_text(size = 8, color = "black"),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.line = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
        plot.margin = ggplot2::margin(5, 10, 2, 0)
    )

if (nrow(faers_tto_pairwise_annotations) > 0) {
    p_detail <- p_detail +
        geom_segment(
            data = faers_tto_pairwise_annotations,
            aes(x = x_position, xend = x_position, y = ymin, yend = ymax),
            inherit.aes = FALSE,
            color = "black",
            linewidth = 0.4
        ) +
        geom_segment(
            data = faers_tto_pairwise_annotations,
            aes(x = x_position - tip_width, xend = x_position, y = ymin, yend = ymin),
            inherit.aes = FALSE,
            color = "black",
            linewidth = 0.4
        ) +
        geom_segment(
            data = faers_tto_pairwise_annotations,
            aes(x = x_position - tip_width, xend = x_position, y = ymax, yend = ymax),
            inherit.aes = FALSE,
            color = "black",
            linewidth = 0.4
        ) +
        geom_text(
            data = faers_tto_pairwise_annotations,
            aes(x = label_x, y = label_y, label = significance_vertical),
            inherit.aes = FALSE,
            hjust = 0,
            lineheight = 0.7,
            size = 8 / ggplot2::.pt,
            color = "black"
        )
}

p <- p_detail

print(p)

## 3. 保存图形 ####
ggsave(
    filename = faers_tto_plot_file,
    plot = p,
    width = 8.6,
    height = 5,
    units = "cm",
    device = cairo_pdf
)


# 6. 3图：FAERS fatality comparison ####

## 1. 读取并处理绘图数据 ####
faers_fatality_summary <- readxl::read_xlsx(faers_summary_file, sheet = "fatality_summary") %>%
    dplyr::mutate(section = as.character(section))

faers_fatality_plot_data <- faers_fatality_summary %>%
    dplyr::filter(section == "group_summary") %>%
    dplyr::mutate(
        pattern_group = as.character(pattern_group),
        fatality_pct = suppressWarnings(as.numeric(fatality_pct)),
        total_n = suppressWarnings(as.numeric(total_n))
    )

faers_fatality_pairwise <- faers_fatality_summary %>%
    dplyr::filter(section == "pairwise_test") %>%
    dplyr::mutate(
        group1 = as.character(group1),
        group2 = as.character(group2),
        significance = as.character(significance),
        p_value = suppressWarnings(as.numeric(p_value))
    )

faers_fatality_global_test <- faers_fatality_summary %>%
    dplyr::filter(section == "global_test") %>%
    dplyr::slice(1)

faers_fatality_order <- build_descending_group_order(
    summary_df = faers_fatality_plot_data,
    value_col = "fatality_pct",
    fallback_order = pattern_display_order
)

faers_fatality_plot_data <- faers_fatality_plot_data %>%
    dplyr::mutate(
        pattern_group = factor(pattern_group, levels = faers_fatality_order),
        pattern_group_label = factor(
            as.character(pattern_group),
            levels = faers_fatality_order
        )
    ) %>%
    dplyr::arrange(pattern_group)

faers_fatality_annotations <- build_fatality_annotation_df(
    pairwise_df = faers_fatality_pairwise,
    plot_data = faers_fatality_plot_data,
    pattern_order = faers_fatality_order,
    line_gap_scale = 1.4,
    line_base_scale = 1.3
)

faers_fatality_global_p <- suppressWarnings(as.numeric(faers_fatality_global_test$p_value[[1]]))
faers_fatality_global_label <- format_overall_p_value_plotmath(faers_fatality_global_p)
faers_fatality_global_x <- max(1.5, (length(faers_fatality_order) + 1) / 2 - 0.45)

faers_fatality_ymax <- max(faers_fatality_plot_data$fatality_pct, na.rm = TRUE)
if (!is.finite(faers_fatality_ymax)) {
    faers_fatality_ymax <- 0
}

if (nrow(faers_fatality_annotations) > 0) {
    faers_fatality_ymax <- max(faers_fatality_ymax, faers_fatality_annotations$label_y) + 8
} else {
    faers_fatality_ymax <- faers_fatality_ymax + 10
}

## 2. 绘图 ####
p <- ggplot(
    faers_fatality_plot_data,
    aes(x = pattern_group_label, y = fatality_pct, fill = pattern_group)
) +
    geom_col(width = 0.7) +
    geom_text(
        aes(label = sprintf("%.1f", fatality_pct)),
        vjust = -0.4,
        size = 8 / ggplot2::.pt,
        color = "black"
    ) +
    scale_fill_manual(values = pattern_color_map, guide = "none") +
    scale_y_continuous(limits = c(0, faers_fatality_ymax)) +
    labs(
        #title = "FAERS Fatality by Cancer Pattern Group",
        x = NULL,
        y = "Fatality (%)"
    ) +
    annotate(
        geom = "text",
        x = faers_fatality_global_x,
        y = faers_fatality_ymax - 2,
        label = faers_fatality_global_label,
        parse = TRUE,
        size = 8 / ggplot2::.pt,
        color = "black"
    ) +
    theme_bw(base_family = "Helvetica") +
    theme(
        plot.title = element_text(size = 8, color = "black", face = "bold", hjust = 0.5, family = "Helvetica"),
        axis.title.x = element_text(size = 8, color = "black", family = "Helvetica"),
        axis.title.y = element_text(size = 8, color = "black", family = "Helvetica"),
        axis.text.x = element_text(size = 8, color = "black", angle = 45, hjust = 1, family = "Helvetica"),
        axis.text.y = element_text(size = 8, color = "black", family = "Helvetica"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(t = 2, r = 2, b = 0, l = 2)
    )

if (nrow(faers_fatality_annotations) > 0) {
    p <- p +
        geom_segment(
            data = faers_fatality_annotations,
            aes(x = xmin, xend = xmax, y = y_position, yend = y_position),
            inherit.aes = FALSE,
            color = "black",
            linewidth = 0.4
        ) +
        geom_segment(
            data = faers_fatality_annotations,
            aes(x = xmin, xend = xmin, y = y_position - tip_height, yend = y_position),
            inherit.aes = FALSE,
            color = "black",
            linewidth = 0.4
        ) +
        geom_segment(
            data = faers_fatality_annotations,
            aes(x = xmax, xend = xmax, y = y_position - tip_height, yend = y_position),
            inherit.aes = FALSE,
            color = "black",
            linewidth = 0.4
        ) +
        geom_text(
            data = faers_fatality_annotations,
            aes(x = label_x, y = label_y, label = significance),
            inherit.aes = FALSE,
            size = 8 / ggplot2::.pt,
            color = "black"
        )
}

print(p)

## 3. 保存图形 ####
ggsave(
    filename = faers_fatality_plot_file,
    plot = p,
    width = 6,
    height = 6,
    units = "cm",
    device = cairo_pdf
)


# 7. 4图：JADER SOC reporting spectrum heatmap ####

## 1. 读取并处理绘图数据 ####
jader_soc_plot_data <- readxl::read_xlsx(jader_summary_file, sheet = "SOC_spectrum") %>%
    dplyr::mutate(
        pattern_group = as.character(pattern_group),
        soc_abbrev = as.character(soc_abbrev),
        soc_order = suppressWarnings(as.numeric(soc_order)),
        RP = suppressWarnings(as.numeric(RP))
    )

jader_pattern_group_n <- readxl::read_xlsx(jader_summary_file, sheet = "pattern_group_n") %>%
    dplyr::mutate(
        pattern_group = as.character(pattern_group),
        mapped_ici_n = suppressWarnings(as.numeric(mapped_ici_n))
    )

jader_soc_label_map <- build_pattern_group_label_map(
    summary_df = jader_pattern_group_n,
    n_col = "mapped_ici_n",
    pattern_order = pattern_display_order
)

jader_soc_levels <- jader_soc_plot_data %>%
    dplyr::arrange(soc_order, soc_abbrev) %>%
    dplyr::pull(soc_abbrev) %>%
    unique()

jader_soc_plot_data <- jader_soc_plot_data %>%
    dplyr::mutate(
        label_text = sprintf("%.1f", RP),
        label_color = ifelse(
            RP >= max(RP, na.rm = TRUE) * 0.55,
            "white",
            "black"
        ),
        pattern_group_label = factor(
            jader_soc_label_map[pattern_group],
            levels = unname(jader_soc_label_map[pattern_display_order])
        ),
        soc_abbrev = factor(
            soc_abbrev,
            levels = rev(jader_soc_levels)
        )
    )

## 2. 绘图 ####
p <- ggplot(jader_soc_plot_data, aes(x = pattern_group_label, y = soc_abbrev, fill = RP)) +
    geom_tile(color = "white", linewidth = 0.25) +
    geom_text(
        aes(label = label_text, color = label_color),
        size = 10 / ggplot2::.pt
    ) +
    scale_fill_gradient(
        low = "#F6F3EE",
        high = "#B23A48",
        name = "RP (%)"
    ) +
    guides(
        fill = guide_colorbar(
            direction = "horizontal",
            title.position = "left",
            title.theme = element_text(
                size = 14,
                color = "black",
                angle = 0,
                margin = ggplot2::margin(r = 8)
            ),
            barwidth = grid::unit(5.2, "cm"),
            barheight = grid::unit(0.35, "cm")
        )
    ) +
    scale_color_identity() +
    labs(
        x = NULL,
        y = NULL
    ) +
    theme_bw() +
    theme(
        axis.title.x = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14, color = "black"),
        axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 14, color = "black"),
        legend.text = element_text(size = 12, color = "black"),
        legend.position = "top",
        panel.grid = element_blank()
    )

print(p)

## 3. 保存图形 ####
ggsave(
    filename = jader_soc_plot_file,
    plot = p,
    width = 5,
    height = 4,
    device = cairo_pdf
)


# 8. 5图：JADER Time to Onset cumulative incidence ####

## 1. 读取并处理绘图数据 ####
jader_tto_plot_data <- readxl::read_xlsx(jader_summary_file, sheet = "TTO_case_level") %>%
    dplyr::mutate(
        pattern_group = as.character(pattern_group),
        ici_tto = suppressWarnings(as.numeric(ici_tto))
    ) %>%
    dplyr::filter(!is.na(ici_tto))

jader_tto_summary <- readxl::read_xlsx(jader_summary_file, sheet = "TTO_summary") %>%
    dplyr::mutate(section = as.character(section))

jader_tto_group_summary <- jader_tto_summary %>%
    dplyr::filter(section == "group_summary") %>%
    dplyr::mutate(
        pattern_group = as.character(pattern_group),
        n = suppressWarnings(as.numeric(n)),
        median_tto = suppressWarnings(as.numeric(median_tto)),
        q1_tto = suppressWarnings(as.numeric(q1_tto)),
        q3_tto = suppressWarnings(as.numeric(q3_tto))
    )

jader_tto_pairwise <- jader_tto_summary %>%
    dplyr::filter(section == "pairwise_test") %>%
    dplyr::mutate(
        group1 = as.character(group1),
        group2 = as.character(group2),
        significance = as.character(significance),
        p_value = suppressWarnings(as.numeric(p_value))
    )

jader_tto_global_test <- jader_tto_summary %>%
    dplyr::filter(section == "global_test") %>%
    dplyr::slice(1)

jader_tto_label_map <- build_pattern_group_label_map(
    summary_df = jader_tto_group_summary,
    n_col = "n",
    pattern_order = pattern_display_order,
    separator = " "
)

jader_tto_plot_data <- jader_tto_plot_data %>%
    dplyr::mutate(
        pattern_group_label = factor(
            jader_tto_label_map[pattern_group],
            levels = unname(jader_tto_label_map[pattern_display_order])
        ),
        pattern_group_label_detail = factor(
            jader_tto_label_map[pattern_group],
            levels = rev(unname(jader_tto_label_map[pattern_display_order]))
        )
    )

jader_tto_color_map <- stats::setNames(
    unname(pattern_color_map[pattern_display_order]),
    unname(jader_tto_label_map[pattern_display_order])
)

jader_tto_max_x <- if (nrow(jader_tto_plot_data) > 0) {
    max(jader_tto_plot_data$ici_tto, na.rm = TRUE)
} else {
    1
}

jader_tto_detail_levels <- rev(unname(jader_tto_label_map[pattern_display_order]))
jader_tto_detail_levels_plain <- rev(unname(pattern_display_order))

jader_tto_plot_data <- jader_tto_plot_data %>%
    dplyr::mutate(
        pattern_group_y = match(
            as.character(pattern_group_label_detail),
            jader_tto_detail_levels
        )
    )

jader_tto_boxplot_labels <- build_tto_boxplot_label_df(
    group_summary_df = jader_tto_group_summary,
    label_map = jader_tto_label_map,
    pattern_order = pattern_display_order
) %>%
    dplyr::mutate(
        label_y = pattern_group_y + dplyr::case_when(
            stat_type == "q1_tto" ~ -0.23,
            stat_type == "median_tto" ~ 0,
            stat_type == "q3_tto" ~ 0.23,
            TRUE ~ 0
        ),
        fontface = dplyr::if_else(stat_type == "median_tto", "bold", "plain")
    )

jader_tto_pairwise_annotations <- build_tto_pairwise_annotation_df(
    pairwise_df = jader_tto_pairwise,
    label_map = jader_tto_label_map,
    pattern_order = pattern_display_order,
    detail_label_levels = jader_tto_detail_levels,
    plot_max_x = jader_tto_max_x
)

jader_tto_global_p <- suppressWarnings(as.numeric(jader_tto_global_test$p_value[[1]]))
jader_tto_global_label <- paste("Overall", format_p_value_text(jader_tto_global_p))

jader_tto_detail_xmax <- jader_tto_max_x
if (nrow(jader_tto_pairwise_annotations) > 0) {
    jader_tto_detail_xmax <- max(
        jader_tto_detail_xmax,
        jader_tto_pairwise_annotations$label_x + jader_tto_pairwise_annotations$tip_width
    )
}
if (!is.finite(jader_tto_detail_xmax) || jader_tto_detail_xmax <= 0) {
    jader_tto_detail_xmax <- 1
}
jader_tto_detail_xmax <- jader_tto_detail_xmax + max(8, jader_tto_max_x * 0.02)

## 2. 绘图 ####
p_curve <- ggplot(jader_tto_plot_data, aes(x = ici_tto, color = pattern_group_label)) +
    stat_ecdf(geom = "step", pad = FALSE, linewidth = 0.9) +
    scale_color_manual(values = jader_tto_color_map, name = "Cancer group") +
    scale_y_continuous(
        labels = scales::percent_format(accuracy = 1),
        limits = c(0, 1)
    ) +
    coord_cartesian(xlim = c(0, jader_tto_max_x)) +
    labs(
        x = "Time to Onset (Days)",
        y = "Cumulative Incidence"
    ) +
    annotate(
        geom = "text",
        x = jader_tto_max_x * 0.5,
        y = 0.98,
        label = jader_tto_global_label,
        size = 12 / ggplot2::.pt,
        color = "black"
    ) +
    theme_bw() +
    theme(
        axis.title.x = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14, color = "black", margin = ggplot2::margin(r = 4)),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.line = element_blank(),
        legend.title = element_text(size = 14, color = "black"),
        legend.text = element_text(size = 12, color = "black"),
        legend.position = c(0.98, 0.98),
        legend.justification = c(1, 1),
        legend.background = element_rect(fill = scales::alpha("white", 0.9), color = "black", linewidth = 0.3),
        legend.key = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
        plot.margin = ggplot2::margin(10, 10, 5, 10)
    )

p_detail <- ggplot(
    jader_tto_plot_data,
    aes(
        x = ici_tto,
        y = pattern_group_y,
        group = pattern_group_label_detail,
        color = pattern_group_label,
        fill = pattern_group_label
    )
) +
    geom_boxplot(
        orientation = "y",
        outlier.shape = NA,
        width = 0.6,
        alpha = 0.18,
        linewidth = 0.4,
        show.legend = FALSE
    ) +
    geom_jitter(
        height = 0.16,
        width = 0,
        alpha = 0.18,
        size = 0.7,
        show.legend = FALSE
    ) +
    geom_label(
        data = jader_tto_boxplot_labels,
        aes(
            x = stat_value,
            y = label_y,
            label = label,
            fontface = fontface
        ),
        inherit.aes = FALSE,
        size = 9 / ggplot2::.pt,
        linewidth = 0.15,
        label.padding = grid::unit(0.12, "lines"),
        fill = "white",
        color = "black"
    ) +
    scale_color_manual(values = jader_tto_color_map, guide = "none") +
    scale_fill_manual(values = jader_tto_color_map, guide = "none") +
    scale_y_continuous(
        breaks = seq_along(jader_tto_detail_levels),
        labels = jader_tto_detail_levels_plain,
        expand = ggplot2::expansion(mult = c(0.05, 0.08))
    ) +
    coord_cartesian(xlim = c(0, jader_tto_detail_xmax), clip = "off") +
    labs(
        x = "Time to Onset (Days)",
        y = NULL
    ) +
    theme_bw() +
    theme(
        axis.title.x = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.line = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
        plot.margin = ggplot2::margin(5, 10, 10, 10)
    )

if (nrow(jader_tto_pairwise_annotations) > 0) {
    p_detail <- p_detail +
        geom_segment(
            data = jader_tto_pairwise_annotations,
            aes(x = x_position, xend = x_position, y = ymin, yend = ymax),
            inherit.aes = FALSE,
            color = "black",
            linewidth = 0.4
        ) +
        geom_segment(
            data = jader_tto_pairwise_annotations,
            aes(x = x_position - tip_width, xend = x_position, y = ymin, yend = ymin),
            inherit.aes = FALSE,
            color = "black",
            linewidth = 0.4
        ) +
        geom_segment(
            data = jader_tto_pairwise_annotations,
            aes(x = x_position - tip_width, xend = x_position, y = ymax, yend = ymax),
            inherit.aes = FALSE,
            color = "black",
            linewidth = 0.4
        ) +
        geom_text(
            data = jader_tto_pairwise_annotations,
            aes(x = label_x, y = label_y, label = significance),
            inherit.aes = FALSE,
            hjust = 0,
            size = 11 / ggplot2::.pt,
            color = "black"
        )
}

p <- (p_curve / p_detail) +
    patchwork::plot_layout(heights = c(1.1, 1.4))

print(p)

## 3. 保存图形 ####
ggsave(
    filename = jader_tto_plot_file,
    plot = p,
    width = 10,
    height = 8.5,
    device = cairo_pdf
)


# 9. 6图：JADER fatality comparison ####

## 1. 读取并处理绘图数据 ####
jader_fatality_summary <- readxl::read_xlsx(jader_summary_file, sheet = "fatality_summary") %>%
    dplyr::mutate(section = as.character(section))

jader_fatality_plot_data <- jader_fatality_summary %>%
    dplyr::filter(section == "group_summary") %>%
    dplyr::mutate(
        pattern_group = as.character(pattern_group),
        fatality_pct = suppressWarnings(as.numeric(fatality_pct)),
        total_n = suppressWarnings(as.numeric(total_n))
    )

jader_fatality_pairwise <- jader_fatality_summary %>%
    dplyr::filter(section == "pairwise_test") %>%
    dplyr::mutate(
        group1 = as.character(group1),
        group2 = as.character(group2),
        significance = as.character(significance),
        p_value = suppressWarnings(as.numeric(p_value))
    )

jader_fatality_global_test <- jader_fatality_summary %>%
    dplyr::filter(section == "global_test") %>%
    dplyr::slice(1)

jader_fatality_order <- build_descending_group_order(
    summary_df = jader_fatality_plot_data,
    value_col = "fatality_pct",
    fallback_order = pattern_display_order
)

jader_fatality_plot_data <- jader_fatality_plot_data %>%
    dplyr::mutate(
        pattern_group = factor(pattern_group, levels = jader_fatality_order),
        pattern_group_label = factor(
            as.character(pattern_group),
            levels = jader_fatality_order
        )
    ) %>%
    dplyr::arrange(pattern_group)

jader_fatality_annotations <- build_fatality_annotation_df(
    pairwise_df = jader_fatality_pairwise,
    plot_data = jader_fatality_plot_data,
    pattern_order = jader_fatality_order
)

jader_fatality_global_p <- suppressWarnings(as.numeric(jader_fatality_global_test$p_value[[1]]))
jader_fatality_global_label <- format_overall_p_value_plotmath(jader_fatality_global_p)
jader_fatality_global_x <- max(1.5, (length(jader_fatality_order) + 1) / 2 - 0.45)

jader_fatality_ymax <- max(jader_fatality_plot_data$fatality_pct, na.rm = TRUE)
if (!is.finite(jader_fatality_ymax)) {
    jader_fatality_ymax <- 0
}

if (nrow(jader_fatality_annotations) > 0) {
    jader_fatality_ymax <- max(jader_fatality_ymax, jader_fatality_annotations$label_y) + 8
} else {
    jader_fatality_ymax <- jader_fatality_ymax + 10
}

## 2. 绘图 ####
p <- ggplot(
    jader_fatality_plot_data,
    aes(x = pattern_group_label, y = fatality_pct, fill = pattern_group)
) +
    geom_col(width = 0.7) +
    geom_text(
        aes(label = sprintf("%.1f", fatality_pct)),
        vjust = -0.4,
        size = 8 / ggplot2::.pt,
        color = "black"
    ) +
    scale_fill_manual(values = pattern_color_map, guide = "none") +
    scale_y_continuous(limits = c(0, jader_fatality_ymax)) +
    labs(
        #title = "JADER Fatality by Cancer Pattern Group",
        x = NULL,
        y = "Fatality (%)"
    ) +
    annotate(
        geom = "text",
        x = jader_fatality_global_x,
        y = jader_fatality_ymax - 2,
        label = jader_fatality_global_label,
        parse = TRUE,
        size = 8 / ggplot2::.pt,
        color = "black"
    ) +
    theme_bw(base_family = "Helvetica") +
    theme(
        plot.title = element_text(size = 8, color = "black", face = "bold", hjust = 0.5, family = "Helvetica"),
        axis.title.x = element_text(size = 8, color = "black", family = "Helvetica"),
        axis.title.y = element_text(size = 8, color = "black", family = "Helvetica"),
        axis.text.x = element_text(size = 8, color = "black", angle = 45, hjust = 1, family = "Helvetica"),
        axis.text.y = element_text(size = 8, color = "black", family = "Helvetica"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(t = 2, r = 2, b = 0, l = 2)
    )

if (nrow(jader_fatality_annotations) > 0) {
    p <- p +
        geom_segment(
            data = jader_fatality_annotations,
            aes(x = xmin, xend = xmax, y = y_position, yend = y_position),
            inherit.aes = FALSE,
            color = "black",
            linewidth = 0.4
        ) +
        geom_segment(
            data = jader_fatality_annotations,
            aes(x = xmin, xend = xmin, y = y_position - tip_height, yend = y_position),
            inherit.aes = FALSE,
            color = "black",
            linewidth = 0.4
        ) +
        geom_segment(
            data = jader_fatality_annotations,
            aes(x = xmax, xend = xmax, y = y_position - tip_height, yend = y_position),
            inherit.aes = FALSE,
            color = "black",
            linewidth = 0.4
        ) +
        geom_text(
            data = jader_fatality_annotations,
            aes(x = label_x, y = label_y, label = significance),
            inherit.aes = FALSE,
            size = 8 / ggplot2::.pt,
            color = "black"
        )
}

print(p)

## 3. 保存图形 ####
ggsave(
    filename = jader_fatality_plot_file,
    plot = p,
    width = 6,
    height = 6,
    units = "cm",
    device = cairo_pdf
)


# 10. 7图：FAERS TTO Weibull analysis ####

## 1. 读取并处理绘图数据 ####
faers_weibull_case_level <- readxl::read_xlsx(faers_summary_file, sheet = "TTO_case_level") %>%
    dplyr::mutate(
        pattern_group = as.character(pattern_group),
        ici_tto = suppressWarnings(as.numeric(ici_tto))
    ) %>%
    dplyr::filter(!is.na(ici_tto))

faers_weibull_summary <- readxl::read_xlsx(faers_summary_file, sheet = "TTO_weibull_summary") %>%
    dplyr::mutate(
        pattern_group = as.character(pattern_group),
        shape = suppressWarnings(as.numeric(shape)),
        scale = suppressWarnings(as.numeric(scale)),
        zero_tto_n = suppressWarnings(as.numeric(zero_tto_n)),
        zero_replaced_value = suppressWarnings(as.numeric(zero_replaced_value)),
        fit_status = as.character(fit_status),
        hazard_type = as.character(hazard_type)
    )

faers_weibull_density <- readxl::read_xlsx(faers_summary_file, sheet = "TTO_weibull_density") %>%
    dplyr::mutate(
        pattern_group = as.character(pattern_group),
        x = suppressWarnings(as.numeric(x)),
        density = suppressWarnings(as.numeric(density))
    )

faers_weibull_display_start <- build_weibull_display_start_df(
    case_level_df = faers_weibull_case_level,
    summary_df = faers_weibull_summary
)

faers_weibull_success_summary <- faers_weibull_summary %>%
    dplyr::filter(fit_status == "success") %>%
    dplyr::left_join(faers_weibull_display_start, by = "pattern_group") %>%
    dplyr::mutate(
        pattern_order = match(pattern_group, pattern_display_order)
    ) %>%
    dplyr::filter(!is.na(pattern_order)) %>%
    dplyr::arrange(pattern_order) %>%
    dplyr::mutate(
        failure_display = format_weibull_failure_display(hazard_type),
        panel_letter = LETTERS[dplyr::row_number()],
        facet_label = paste0(
            pattern_group,
            " - ",
            failure_display,
            "\nShape = ",
            sprintf("%.3f", shape),
            ", Scale = ",
            sprintf("%.3f", scale)
        )
    )

faers_weibull_facet_levels <- faers_weibull_success_summary$facet_label

faers_weibull_case_level <- faers_weibull_case_level %>%
    dplyr::filter(pattern_group %in% faers_weibull_success_summary$pattern_group) %>%
    dplyr::left_join(
        faers_weibull_success_summary %>%
            dplyr::select(pattern_group, facet_label, panel_letter),
        by = "pattern_group"
    ) %>%
    dplyr::mutate(
        facet_label = factor(facet_label, levels = faers_weibull_facet_levels)
    )

faers_weibull_density <- faers_weibull_density %>%
    dplyr::filter(pattern_group %in% faers_weibull_success_summary$pattern_group) %>%
    dplyr::left_join(
        faers_weibull_success_summary %>%
            dplyr::select(pattern_group, facet_label, panel_letter, display_x_min),
        by = "pattern_group"
    ) %>%
    dplyr::filter(is.na(display_x_min) | x >= display_x_min) %>%
    dplyr::mutate(
        facet_label = factor(facet_label, levels = faers_weibull_facet_levels)
    )

## 2. 绘图 ####
if (nrow(faers_weibull_success_summary) > 0) {
    p <- ggplot(faers_weibull_case_level, aes(x = ici_tto)) +
        geom_histogram(
            aes(y = after_stat(density)),
            bins = 25,
            fill = "#D9D9D9",
            alpha = 0.9,
            color = "black",
            linewidth = 0.5,
            show.legend = FALSE
        ) +
        geom_line(
            data = faers_weibull_density,
            aes(x = x, y = density),
            inherit.aes = FALSE,
            color = "red",
            linewidth = 0.9,
            show.legend = FALSE
        ) +
        facet_wrap(~ facet_label, ncol = 3, scales = "free") +
        coord_cartesian(clip = "off") +
        labs(
            x = "Time to Onset (Days)",
            y = "Density"
        ) +
        theme_bw() +
        theme(
            axis.title.x = element_text(size = 14, color = "black"),
            axis.title.y = element_text(size = 14, color = "black"),
            axis.text.x = element_text(size = 12, color = "black"),
            axis.text.y = element_text(size = 12, color = "black"),
            strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
            strip.text = element_text(size = 12, color = "black", face = "bold", lineheight = 1.15),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
            panel.spacing = grid::unit(1.2, "lines"),
            plot.margin = ggplot2::margin(12, 12, 10, 18)
        )
} else {
    p <- ggplot() +
        annotate(
            geom = "text",
            x = 0.5,
            y = 0.5,
            label = "No successful Weibull fit available.",
            size = 12 / ggplot2::.pt,
            color = "black"
        ) +
        theme_void()
}

print(p)

## 3. 保存图形 ####
ggsave(
    filename = faers_weibull_plot_file,
    plot = p,
    width = 15,
    height = 11,
    device = cairo_pdf
)


# 11. 8图：JADER TTO Weibull analysis ####

## 1. 读取并处理绘图数据 ####
jader_weibull_case_level <- readxl::read_xlsx(jader_summary_file, sheet = "TTO_case_level") %>%
    dplyr::mutate(
        pattern_group = as.character(pattern_group),
        ici_tto = suppressWarnings(as.numeric(ici_tto))
    ) %>%
    dplyr::filter(!is.na(ici_tto))

jader_weibull_summary <- readxl::read_xlsx(jader_summary_file, sheet = "TTO_weibull_summary") %>%
    dplyr::mutate(
        pattern_group = as.character(pattern_group),
        shape = suppressWarnings(as.numeric(shape)),
        scale = suppressWarnings(as.numeric(scale)),
        zero_tto_n = suppressWarnings(as.numeric(zero_tto_n)),
        zero_replaced_value = suppressWarnings(as.numeric(zero_replaced_value)),
        fit_status = as.character(fit_status),
        hazard_type = as.character(hazard_type)
    )

jader_weibull_density <- readxl::read_xlsx(jader_summary_file, sheet = "TTO_weibull_density") %>%
    dplyr::mutate(
        pattern_group = as.character(pattern_group),
        x = suppressWarnings(as.numeric(x)),
        density = suppressWarnings(as.numeric(density))
    )

jader_weibull_display_start <- build_weibull_display_start_df(
    case_level_df = jader_weibull_case_level,
    summary_df = jader_weibull_summary
)

jader_weibull_success_summary <- jader_weibull_summary %>%
    dplyr::filter(fit_status == "success") %>%
    dplyr::left_join(jader_weibull_display_start, by = "pattern_group") %>%
    dplyr::mutate(
        pattern_order = match(pattern_group, pattern_display_order)
    ) %>%
    dplyr::filter(!is.na(pattern_order)) %>%
    dplyr::arrange(pattern_order) %>%
    dplyr::mutate(
        failure_display = format_weibull_failure_display(hazard_type),
        panel_letter = LETTERS[dplyr::row_number()],
        facet_label = paste0(
            pattern_group,
            " - ",
            failure_display,
            "\nShape = ",
            sprintf("%.3f", shape),
            ", Scale = ",
            sprintf("%.3f", scale)
        )
    )

jader_weibull_facet_levels <- jader_weibull_success_summary$facet_label

jader_weibull_case_level <- jader_weibull_case_level %>%
    dplyr::filter(pattern_group %in% jader_weibull_success_summary$pattern_group) %>%
    dplyr::left_join(
        jader_weibull_success_summary %>%
            dplyr::select(pattern_group, facet_label, panel_letter),
        by = "pattern_group"
    ) %>%
    dplyr::mutate(
        facet_label = factor(facet_label, levels = jader_weibull_facet_levels)
    )

jader_weibull_density <- jader_weibull_density %>%
    dplyr::filter(pattern_group %in% jader_weibull_success_summary$pattern_group) %>%
    dplyr::left_join(
        jader_weibull_success_summary %>%
            dplyr::select(pattern_group, facet_label, panel_letter, display_x_min),
        by = "pattern_group"
    ) %>%
    dplyr::filter(is.na(display_x_min) | x >= display_x_min) %>%
    dplyr::mutate(
        facet_label = factor(facet_label, levels = jader_weibull_facet_levels)
    )

## 2. 绘图 ####
if (nrow(jader_weibull_success_summary) > 0) {
    p <- ggplot(jader_weibull_case_level, aes(x = ici_tto)) +
        geom_histogram(
            aes(y = after_stat(density)),
            bins = 25,
            fill = "#D9D9D9",
            alpha = 0.9,
            color = "black",
            linewidth = 0.5,
            show.legend = FALSE
        ) +
        geom_line(
            data = jader_weibull_density,
            aes(x = x, y = density),
            inherit.aes = FALSE,
            color = "red",
            linewidth = 0.9,
            show.legend = FALSE
        ) +
        facet_wrap(~ facet_label, ncol = 3, scales = "free") +
        coord_cartesian(clip = "off") +
        labs(
            x = "Time to Onset (Days)",
            y = "Density"
        ) +
        theme_bw() +
        theme(
            axis.title.x = element_text(size = 14, color = "black"),
            axis.title.y = element_text(size = 14, color = "black"),
            axis.text.x = element_text(size = 12, color = "black"),
            axis.text.y = element_text(size = 12, color = "black"),
            strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
            strip.text = element_text(size = 12, color = "black", face = "bold", lineheight = 1.15),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
            panel.spacing = grid::unit(1.2, "lines"),
            plot.margin = ggplot2::margin(12, 12, 10, 18)
        )
} else {
    p <- ggplot() +
        annotate(
            geom = "text",
            x = 0.5,
            y = 0.5,
            label = "No successful Weibull fit available.",
            size = 12 / ggplot2::.pt,
            color = "black"
        ) +
        theme_void()
}

print(p)

## 3. 保存图形 ####
ggsave(
    filename = jader_weibull_plot_file,
    plot = p,
    width = 15,
    height = 11,
    device = cairo_pdf
)

message("========== 癌种群 irAE 表型绘图完成 ==========")
