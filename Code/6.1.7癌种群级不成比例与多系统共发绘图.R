# 癌种群级不成比例与多系统共发绘图（FAERS + JADER）####

rm(list = ls())
gc()

# 1. 加载配置 ####
source("Code/config.R")


# 2. 结果目录与输入文件 ####
analysis_script_name <- "6.1.6癌种群级不成比例与多系统共发分析"
analysis_results_path <- file.path(results_path, analysis_script_name)

script_name <- "6.1.7癌种群级不成比例与多系统共发绘图"
script_results_path <- file.path(results_path, script_name)

if (!dir.exists(script_results_path)) {
    dir.create(script_results_path, recursive = TRUE, showWarnings = FALSE)
}

faers_summary_file <- file.path(analysis_results_path, "FAERS_癌种群级ROR与多系统共发汇总.xlsx")
jader_summary_file <- file.path(analysis_results_path, "JADER_癌种群级ROR与多系统共发汇总.xlsx")

faers_soc_ror_plot_file <- file.path(script_results_path, "图1_FAERS_癌种群级SOC_ROR热图.pdf")
faers_pt_ror_plot_file <- file.path(script_results_path, "图2_FAERS_癌种群级重点PT_ROR热图.pdf")
faers_organ_count_plot_file <- file.path(script_results_path, "图3_FAERS_多系统共发箱线图.pdf")
jader_soc_ror_plot_file <- file.path(script_results_path, "图4_JADER_癌种群级SOC_ROR热图.pdf")
jader_pt_ror_plot_file <- file.path(script_results_path, "图5_JADER_癌种群级重点PT_ROR热图.pdf")
jader_organ_count_plot_file <- file.path(script_results_path, "图6_JADER_多系统共发箱线图.pdf")
faers_dominant_soc_ror_plot_file <- file.path(script_results_path, "图7_FAERS_主导SOC复合ROR森林图.pdf")
jader_dominant_soc_ror_plot_file <- file.path(script_results_path, "图8_JADER_主导SOC复合ROR森林图.pdf")
faers_dominant_soc_ror_simple_plot_file <- file.path(script_results_path, "图7简化版_FAERS_主导SOC复合ROR森林图.pdf")
jader_dominant_soc_ror_simple_plot_file <- file.path(script_results_path, "图8简化版_JADER_主导SOC复合ROR森林图.pdf")
combined_dominant_soc_ror_plot_file <- file.path(script_results_path, "图9_双库主导SOC复合ROR森林图.pdf")

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

message("========== 癌种群级不成比例与多系统共发绘图开始 ==========")
message("分析结果目录: ", analysis_results_path)
message("绘图结果目录: ", script_results_path)


# 3. 辅助函数 ####

## 函数：format_case_n ####
## 功能：格式化样本量为带千分位分隔符 ####
## 输入：x，数值型样本量 ####
## 输出：返回格式化字符值 ####
format_case_n <- function(x) {
    format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
}

## 函数：format_log2_ror_for_fill ####
## 功能：将ROR转为log2刻度并截尾，方便颜色映射 ####
## 输入：ror，ROR数值向量；clip，截尾绝对值 ####
## 输出：返回截尾后的log2 ROR数值向量 ####
format_log2_ror_for_fill <- function(ror, clip = 3) {
    log2_ror <- log2(pmax(ror, 1e-06))
    log2_ror[is.na(ror)] <- NA_real_
    pmax(pmin(log2_ror, clip), -clip)
}

## 函数：build_signal_label ####
## 功能：根据ROR信号布尔生成显著性标记 ####
## 输入：signal_ror，逻辑型ROR信号向量 ####
## 输出：返回字符型标记向量 ####
build_signal_label <- function(signal_ror) {
    dplyr::case_when(
        is.na(signal_ror) ~ "",
        signal_ror == TRUE ~ "*",
        TRUE ~ ""
    )
}

## 函数：format_p_value_compact ####
## 功能：将P值格式化为森林图右侧文本使用的紧凑展示 ####
## 输入：p，数值型P值向量 ####
## 输出：返回字符型格式化P值向量 ####
format_p_value_compact <- function(p) {
    vapply(p, function(value) {
        if (is.na(value)) {
            return("P = NA")
        }
        if (value < 0.001) {
            return("P < 0.001")
        }
        if (value < 0.01) {
            return(sprintf("P = %.3f", value))
        }
        sprintf("P = %.3f", value)
    }, character(1))
}

## 函数：build_pattern_group_label_map ####
## 功能：为pattern group构建带样本量的展示标签 ####
## 输入：summary_df，包含pattern_group与样本量的数据框；n_col，样本量列名；pattern_order，固定pattern顺序；separator，连接符 ####
## 输出：返回pattern group标签命名向量 ####
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

## 函数：build_ror_heatmap ####
## 功能：构建pattern group × 指标的ROR热图（log2 fill + 文本数值 + 信号星号） ####
## 输入：plot_data，含pattern_group_label / row_label / ror / ror_lower / signal_ror的绘图数据；fill_low/mid/high，配色 ####
## 输出：返回ggplot对象 ####
build_ror_heatmap <- function(plot_data, fill_low = "#3B6E8F", fill_mid = "#F6F3EE", fill_high = "#7E1D1D") {
    ggplot(plot_data, aes(x = pattern_group_label, y = row_label, fill = log2_ror_clip)) +
        geom_tile(color = "white", linewidth = 0.25) +
        geom_text(
            aes(label = ror_text),
            size = 9 / ggplot2::.pt,
            color = "black"
        ) +
        geom_text(
            aes(label = signal_label),
            nudge_y = 0.32,
            size = 11 / ggplot2::.pt,
            color = "black",
            fontface = "bold"
        ) +
        scale_fill_gradient2(
            low = fill_low,
            mid = fill_mid,
            high = fill_high,
            midpoint = 0,
            limits = c(-3, 3),
            name = "log2(ROR)",
            oob = scales::squish
        ) +
        scale_y_discrete(limits = rev) +
        guides(
            fill = guide_colorbar(
                direction = "horizontal",
                title.position = "left",
                title.theme = element_text(size = 14, color = "black", angle = 0, margin = ggplot2::margin(r = 8)),
                barwidth = grid::unit(5.2, "cm"),
                barheight = grid::unit(0.35, "cm")
            )
        ) +
        labs(x = NULL, y = NULL) +
        theme_bw() +
        theme(
            axis.title.x = element_text(size = 14, color = "black"),
            axis.title.y = element_text(size = 14, color = "black"),
            axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1),
            axis.text.y = element_text(size = 12, color = "black", lineheight = 1.05),
            legend.title = element_text(size = 14, color = "black"),
            legend.text = element_text(size = 12, color = "black"),
            legend.position = "top",
            panel.grid = element_blank()
        )
}


# 4. 1图：FAERS 癌种群级 SOC ROR 热图 ####

## 1. 读取并处理绘图数据 ####
faers_soc_ror <- readxl::read_xlsx(faers_summary_file, sheet = "SOC_ROR") %>%
    dplyr::mutate(
        pattern_group = as.character(pattern_group),
        soc_abbrev = as.character(soc_abbrev),
        ror = suppressWarnings(as.numeric(ror)),
        ror_lower = suppressWarnings(as.numeric(ror_lower)),
        ror_upper = suppressWarnings(as.numeric(ror_upper)),
        target_event_n = suppressWarnings(as.numeric(target_event_n)),
        signal_ror = as.logical(signal_ror)
    )

faers_pattern_n <- readxl::read_xlsx(faers_summary_file, sheet = "pattern样本量") %>%
    dplyr::mutate(
        pattern_group = as.character(pattern_group),
        mapped_ici_n = suppressWarnings(as.numeric(mapped_ici_n))
    )

faers_label_map <- build_pattern_group_label_map(
    summary_df = faers_pattern_n,
    n_col = "mapped_ici_n",
    pattern_order = pattern_display_order
)

faers_soc_levels <- faers_soc_ror %>%
    dplyr::distinct(soc_abbrev) %>%
    dplyr::arrange(soc_abbrev) %>%
    dplyr::pull(soc_abbrev)

faers_soc_ror_plot_data <- faers_soc_ror %>%
    dplyr::mutate(
        pattern_group_label = factor(
            faers_label_map[pattern_group],
            levels = unname(faers_label_map[pattern_display_order])
        ),
        row_label = factor(soc_abbrev, levels = faers_soc_levels),
        log2_ror_clip = format_log2_ror_for_fill(ror),
        ror_text = ifelse(is.na(ror), "NA", sprintf("%.2f", ror)),
        signal_label = build_signal_label(signal_ror)
    )

## 2. 绘图 ####
p <- build_ror_heatmap(faers_soc_ror_plot_data) +
    labs(y = "irAE SOC")

print(p)

## 3. 保存图形 ####
ggsave(
    filename = faers_soc_ror_plot_file,
    plot = p,
    width = 7,
    height = 11,
    device = cairo_pdf
)


# 5. 2图：FAERS 癌种群级 重点PT ROR 热图 ####

## 1. 读取并处理绘图数据 ####
faers_pt_ror <- readxl::read_xlsx(faers_summary_file, sheet = "重点PT_ROR") %>%
    dplyr::mutate(
        pattern_group = as.character(pattern_group),
        irae_display = as.character(irae_display),
        pt_order = suppressWarnings(as.numeric(pt_order)),
        ror = suppressWarnings(as.numeric(ror)),
        ror_lower = suppressWarnings(as.numeric(ror_lower)),
        ror_upper = suppressWarnings(as.numeric(ror_upper)),
        target_event_n = suppressWarnings(as.numeric(target_event_n)),
        signal_ror = as.logical(signal_ror)
    )

faers_pt_levels <- faers_pt_ror %>%
    dplyr::distinct(pt_order, irae_display) %>%
    dplyr::arrange(pt_order) %>%
    dplyr::pull(irae_display)

faers_pt_ror_plot_data <- faers_pt_ror %>%
    dplyr::mutate(
        pattern_group_label = factor(
            faers_label_map[pattern_group],
            levels = unname(faers_label_map[pattern_display_order])
        ),
        row_label = factor(irae_display, levels = faers_pt_levels),
        log2_ror_clip = format_log2_ror_for_fill(ror),
        ror_text = ifelse(is.na(ror), "NA", sprintf("%.2f", ror)),
        signal_label = build_signal_label(signal_ror)
    )

## 2. 绘图 ####
p <- build_ror_heatmap(faers_pt_ror_plot_data) +
    labs(y = "Target irAE PT")

print(p)

## 3. 保存图形 ####
ggsave(
    filename = faers_pt_ror_plot_file,
    plot = p,
    width = 6.5,
    height = 7.5,
    device = cairo_pdf
)


# 6. 3图：FAERS 多系统共发箱线图 ####

## 1. 读取并处理绘图数据 ####
faers_case_count <- readxl::read_xlsx(faers_summary_file, sheet = "多系统共发_病例级") %>%
    dplyr::mutate(
        pattern_group = as.character(pattern_group),
        has_irae = as.logical(has_irae),
        n_distinct_soc = suppressWarnings(as.numeric(n_distinct_soc)),
        n_distinct_target_pt = suppressWarnings(as.numeric(n_distinct_target_pt))
    ) %>%
    dplyr::filter(has_irae == TRUE)

faers_organ_count_global <- readxl::read_xlsx(faers_summary_file, sheet = "多系统共发_全局检验") %>%
    dplyr::mutate(
        metric = as.character(metric),
        p_value = suppressWarnings(as.numeric(p_value))
    )

faers_organ_count_long <- faers_case_count %>%
    tidyr::pivot_longer(
        cols = c(n_distinct_soc, n_distinct_target_pt),
        names_to = "metric",
        values_to = "value"
    ) %>%
    dplyr::filter(!is.na(value)) %>%
    dplyr::mutate(
        pattern_group = factor(pattern_group, levels = pattern_display_order),
        metric_label = dplyr::case_when(
            metric == "n_distinct_soc" ~ "Distinct irAE SOC count",
            metric == "n_distinct_target_pt" ~ "Distinct target PT count",
            TRUE ~ metric
        ),
        metric_label = factor(metric_label, levels = c("Distinct irAE SOC count", "Distinct target PT count"))
    )

faers_organ_count_global_label <- faers_organ_count_global %>%
    dplyr::transmute(
        metric_label = dplyr::case_when(
            metric == "n_distinct_soc" ~ "Distinct irAE SOC count",
            metric == "n_distinct_target_pt" ~ "Distinct target PT count",
            TRUE ~ metric
        ),
        global_p_text = ifelse(
            is.na(p_value),
            "Overall p = NA",
            ifelse(p_value < 0.001, "Overall p < 0.001", sprintf("Overall p = %.3f", p_value))
        )
    )

faers_organ_count_panel_max <- faers_organ_count_long %>%
    dplyr::group_by(metric_label) %>%
    dplyr::summarise(panel_max = max(value, na.rm = TRUE) + 1, .groups = "drop")

faers_organ_count_global_label <- faers_organ_count_global_label %>%
    dplyr::left_join(faers_organ_count_panel_max, by = "metric_label") %>%
    dplyr::mutate(
        x_position = (length(pattern_display_order) + 1) / 2,
        y_position = panel_max + 0.5
    )

## 2. 绘图 ####
p <- ggplot(
    faers_organ_count_long,
    aes(x = pattern_group, y = value, fill = pattern_group)
) +
    geom_boxplot(
        width = 0.6,
        outlier.size = 0.4,
        outlier.alpha = 0.5,
        linewidth = 0.4,
        color = "black",
        show.legend = FALSE
    ) +
    geom_text(
        data = faers_organ_count_global_label,
        aes(x = x_position, y = y_position, label = global_p_text),
        inherit.aes = FALSE,
        size = 12 / ggplot2::.pt,
        color = "black"
    ) +
    facet_wrap(~ metric_label, ncol = 2, scales = "free_y") +
    scale_fill_manual(values = pattern_color_map, guide = "none") +
    labs(x = NULL, y = "Distinct count per case (mapped irAE-ICI)") +
    theme_bw() +
    theme(
        axis.title.x = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14, color = "black"),
        axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12, color = "black"),
        strip.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
        strip.text = element_text(size = 12, color = "black", face = "bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4)
    )

print(p)

## 3. 保存图形 ####
ggsave(
    filename = faers_organ_count_plot_file,
    plot = p,
    width = 9,
    height = 5,
    device = cairo_pdf
)


# 7. 4图：JADER 癌种群级 SOC ROR 热图 ####

## 1. 读取并处理绘图数据 ####
jader_soc_ror <- readxl::read_xlsx(jader_summary_file, sheet = "SOC_ROR") %>%
    dplyr::mutate(
        pattern_group = as.character(pattern_group),
        soc_abbrev = as.character(soc_abbrev),
        ror = suppressWarnings(as.numeric(ror)),
        ror_lower = suppressWarnings(as.numeric(ror_lower)),
        ror_upper = suppressWarnings(as.numeric(ror_upper)),
        target_event_n = suppressWarnings(as.numeric(target_event_n)),
        signal_ror = as.logical(signal_ror)
    )

jader_pattern_n <- readxl::read_xlsx(jader_summary_file, sheet = "pattern样本量") %>%
    dplyr::mutate(
        pattern_group = as.character(pattern_group),
        mapped_ici_n = suppressWarnings(as.numeric(mapped_ici_n))
    )

jader_label_map <- build_pattern_group_label_map(
    summary_df = jader_pattern_n,
    n_col = "mapped_ici_n",
    pattern_order = pattern_display_order
)

jader_soc_levels <- jader_soc_ror %>%
    dplyr::distinct(soc_abbrev) %>%
    dplyr::arrange(soc_abbrev) %>%
    dplyr::pull(soc_abbrev)

jader_soc_ror_plot_data <- jader_soc_ror %>%
    dplyr::mutate(
        pattern_group_label = factor(
            jader_label_map[pattern_group],
            levels = unname(jader_label_map[pattern_display_order])
        ),
        row_label = factor(soc_abbrev, levels = jader_soc_levels),
        log2_ror_clip = format_log2_ror_for_fill(ror),
        ror_text = ifelse(is.na(ror), "NA", sprintf("%.2f", ror)),
        signal_label = build_signal_label(signal_ror)
    )

## 2. 绘图 ####
p <- build_ror_heatmap(jader_soc_ror_plot_data) +
    labs(y = "irAE SOC")

print(p)

## 3. 保存图形 ####
ggsave(
    filename = jader_soc_ror_plot_file,
    plot = p,
    width = 7,
    height = 11,
    device = cairo_pdf
)


# 8. 5图：JADER 癌种群级 重点PT ROR 热图 ####

## 1. 读取并处理绘图数据 ####
jader_pt_ror <- readxl::read_xlsx(jader_summary_file, sheet = "重点PT_ROR") %>%
    dplyr::mutate(
        pattern_group = as.character(pattern_group),
        irae_display = as.character(irae_display),
        pt_order = suppressWarnings(as.numeric(pt_order)),
        ror = suppressWarnings(as.numeric(ror)),
        ror_lower = suppressWarnings(as.numeric(ror_lower)),
        ror_upper = suppressWarnings(as.numeric(ror_upper)),
        target_event_n = suppressWarnings(as.numeric(target_event_n)),
        signal_ror = as.logical(signal_ror)
    )

jader_pt_levels <- jader_pt_ror %>%
    dplyr::distinct(pt_order, irae_display) %>%
    dplyr::arrange(pt_order) %>%
    dplyr::pull(irae_display)

jader_pt_ror_plot_data <- jader_pt_ror %>%
    dplyr::mutate(
        pattern_group_label = factor(
            jader_label_map[pattern_group],
            levels = unname(jader_label_map[pattern_display_order])
        ),
        row_label = factor(irae_display, levels = jader_pt_levels),
        log2_ror_clip = format_log2_ror_for_fill(ror),
        ror_text = ifelse(is.na(ror), "NA", sprintf("%.2f", ror)),
        signal_label = build_signal_label(signal_ror)
    )

## 2. 绘图 ####
p <- build_ror_heatmap(jader_pt_ror_plot_data) +
    labs(y = "Target irAE PT")

print(p)

## 3. 保存图形 ####
ggsave(
    filename = jader_pt_ror_plot_file,
    plot = p,
    width = 6.5,
    height = 7.5,
    device = cairo_pdf
)


# 9. 6图：JADER 多系统共发箱线图 ####

## 1. 读取并处理绘图数据 ####
jader_case_count <- readxl::read_xlsx(jader_summary_file, sheet = "多系统共发_病例级") %>%
    dplyr::mutate(
        pattern_group = as.character(pattern_group),
        has_irae = as.logical(has_irae),
        n_distinct_soc = suppressWarnings(as.numeric(n_distinct_soc)),
        n_distinct_target_pt = suppressWarnings(as.numeric(n_distinct_target_pt))
    ) %>%
    dplyr::filter(has_irae == TRUE)

jader_organ_count_global <- readxl::read_xlsx(jader_summary_file, sheet = "多系统共发_全局检验") %>%
    dplyr::mutate(
        metric = as.character(metric),
        p_value = suppressWarnings(as.numeric(p_value))
    )

jader_organ_count_long <- jader_case_count %>%
    tidyr::pivot_longer(
        cols = c(n_distinct_soc, n_distinct_target_pt),
        names_to = "metric",
        values_to = "value"
    ) %>%
    dplyr::filter(!is.na(value)) %>%
    dplyr::mutate(
        pattern_group = factor(pattern_group, levels = pattern_display_order),
        metric_label = dplyr::case_when(
            metric == "n_distinct_soc" ~ "Distinct irAE SOC count",
            metric == "n_distinct_target_pt" ~ "Distinct target PT count",
            TRUE ~ metric
        ),
        metric_label = factor(metric_label, levels = c("Distinct irAE SOC count", "Distinct target PT count"))
    )

jader_organ_count_global_label <- jader_organ_count_global %>%
    dplyr::transmute(
        metric_label = dplyr::case_when(
            metric == "n_distinct_soc" ~ "Distinct irAE SOC count",
            metric == "n_distinct_target_pt" ~ "Distinct target PT count",
            TRUE ~ metric
        ),
        global_p_text = ifelse(
            is.na(p_value),
            "Overall p = NA",
            ifelse(p_value < 0.001, "Overall p < 0.001", sprintf("Overall p = %.3f", p_value))
        )
    )

jader_organ_count_panel_max <- jader_organ_count_long %>%
    dplyr::group_by(metric_label) %>%
    dplyr::summarise(panel_max = max(value, na.rm = TRUE) + 1, .groups = "drop")

jader_organ_count_global_label <- jader_organ_count_global_label %>%
    dplyr::left_join(jader_organ_count_panel_max, by = "metric_label") %>%
    dplyr::mutate(
        x_position = (length(pattern_display_order) + 1) / 2,
        y_position = panel_max + 0.5
    )

## 2. 绘图 ####
p <- ggplot(
    jader_organ_count_long,
    aes(x = pattern_group, y = value, fill = pattern_group)
) +
    geom_boxplot(
        width = 0.6,
        outlier.size = 0.4,
        outlier.alpha = 0.5,
        linewidth = 0.4,
        color = "black",
        show.legend = FALSE
    ) +
    geom_text(
        data = jader_organ_count_global_label,
        aes(x = x_position, y = y_position, label = global_p_text),
        inherit.aes = FALSE,
        size = 12 / ggplot2::.pt,
        color = "black"
    ) +
    facet_wrap(~ metric_label, ncol = 2, scales = "free_y") +
    scale_fill_manual(values = pattern_color_map, guide = "none") +
    labs(x = NULL, y = "Distinct count per case (mapped irAE-ICI)") +
    theme_bw() +
    theme(
        axis.title.x = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14, color = "black"),
        axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12, color = "black"),
        strip.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
        strip.text = element_text(size = 12, color = "black", face = "bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4)
    )

print(p)

## 3. 保存图形 ####
ggsave(
    filename = jader_organ_count_plot_file,
    plot = p,
    width = 9,
    height = 5,
    device = cairo_pdf
)


# 9. 7图：FAERS 癌种群级主导SOC复合ROR森林图 ####

## 1. 读取并处理绘图数据 ####
faers_dominant_soc_ror <- readxl::read_xlsx(faers_summary_file, sheet = "主导SOC复合ROR") %>%
    dplyr::mutate(
        pattern_group = as.character(pattern_group),
        dominant_soc = as.character(dominant_soc),
        target_event_n = suppressWarnings(as.numeric(target_event_n)),
        target_total_n = suppressWarnings(as.numeric(target_total_n)),
        ror = suppressWarnings(as.numeric(ror)),
        ror_lower = suppressWarnings(as.numeric(ror_lower)),
        ror_upper = suppressWarnings(as.numeric(ror_upper)),
        p_value = suppressWarnings(as.numeric(p_value)),
        signal_ror = as.logical(signal_ror)
    ) %>%
    dplyr::filter(pattern_group %in% pattern_display_order) %>%
    dplyr::mutate(
        pattern_group = factor(pattern_group, levels = pattern_display_order),
        log2_ror = log2(pmax(ror, 1e-06)),
        log2_lower = log2(pmax(ror_lower, 1e-06)),
        log2_upper = log2(pmax(ror_upper, 1e-06)),
        ror_ci_text = sprintf("%.2f (%.2f-%.2f)", ror, ror_lower, ror_upper),
        p_text = format_p_value_compact(p_value),
        soc_label = paste0(
            pattern_group, " (n=", format(target_total_n, big.mark = ",", scientific = FALSE, trim = TRUE), ")\n",
            "[", dominant_soc, "]"
        )
    )

faers_forest_label_levels <- faers_dominant_soc_ror %>%
    dplyr::arrange(pattern_group) %>%
    dplyr::pull(soc_label)

faers_dominant_soc_ror <- faers_dominant_soc_ror %>%
    dplyr::mutate(soc_label = factor(soc_label, levels = rev(faers_forest_label_levels)))

faers_x_max <- max(faers_dominant_soc_ror$log2_upper, na.rm = TRUE)
faers_x_min <- min(faers_dominant_soc_ror$log2_lower, na.rm = TRUE)
faers_x_range <- faers_x_max - faers_x_min
faers_text_x_ror <- faers_x_max + faers_x_range * 0.10
faers_text_x_p <- faers_x_max + faers_x_range * 0.55

## 2. 绘图：左侧森林图 panel ####
p_faers_forest <- ggplot(
    faers_dominant_soc_ror,
    aes(x = log2_ror, y = soc_label, color = pattern_group)
) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
    geom_errorbarh(
        aes(xmin = log2_lower, xmax = log2_upper),
        height = 0.25,
        linewidth = 0.5
    ) +
    geom_point(aes(fill = pattern_group), shape = 22, size = 2.8, stroke = 0.5, color = "black") +
    scale_color_manual(values = pattern_color_map, guide = "none") +
    scale_fill_manual(values = pattern_color_map, guide = "none") +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) +
    scale_y_discrete(expand = expansion(add = c(0.5, 1.2))) +
    labs(
        x = "log2(ROR)",
        y = NULL
    ) +
    theme_bw() +
    theme(
        axis.title.x = element_text(size = 8, color = "black"),
        axis.text.x = element_text(size = 8, color = "black"),
        axis.text.y = element_text(size = 8, color = "black"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
        plot.margin = margin(t = 5, r = 2, b = 5, l = 5)
    )

## 3. 绘图：右侧文本表 panel ####
p_faers_table <- ggplot(faers_dominant_soc_ror, aes(y = soc_label)) +
    geom_text(
        aes(x = 0.5, label = ror_ci_text),
        hjust = 0.5, size = 8 / ggplot2::.pt, color = "black"
    ) +
    annotate(
        "text", x = 0.5, y = length(faers_forest_label_levels) + 0.7,
        label = "ROR (95% CI)",
        hjust = 0.5, size = 8 / ggplot2::.pt, color = "black", fontface = "bold"
    ) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_discrete(expand = expansion(add = c(0.5, 1.2))) +
    labs(x = " ", y = NULL) +
    theme_bw() +
    theme(
        axis.title.x = element_text(size = 8, color = "transparent"),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.margin = margin(t = 5, r = 0, b = 5, l = 0)
    )

## 4. 横向拼接 ####
p <- patchwork::wrap_plots(p_faers_forest, p_faers_table, nrow = 1, widths = c(1.35, 0.55))

print(p)

## 5. 保存图形 ####
ggsave(
    filename = faers_dominant_soc_ror_plot_file,
    plot = p,
    width = 6,
    height = 4,
    units = "cm",
    device = cairo_pdf
)


# 10. 8图：JADER 癌种群级主导SOC复合ROR森林图 ####

## 1. 读取并处理绘图数据 ####
jader_dominant_soc_ror <- readxl::read_xlsx(jader_summary_file, sheet = "主导SOC复合ROR") %>%
    dplyr::mutate(
        pattern_group = as.character(pattern_group),
        dominant_soc = as.character(dominant_soc),
        target_event_n = suppressWarnings(as.numeric(target_event_n)),
        target_total_n = suppressWarnings(as.numeric(target_total_n)),
        ror = suppressWarnings(as.numeric(ror)),
        ror_lower = suppressWarnings(as.numeric(ror_lower)),
        ror_upper = suppressWarnings(as.numeric(ror_upper)),
        p_value = suppressWarnings(as.numeric(p_value)),
        signal_ror = as.logical(signal_ror)
    ) %>%
    dplyr::filter(pattern_group %in% pattern_display_order) %>%
    dplyr::mutate(
        pattern_group = factor(pattern_group, levels = pattern_display_order),
        log2_ror = log2(pmax(ror, 1e-06)),
        log2_lower = log2(pmax(ror_lower, 1e-06)),
        log2_upper = log2(pmax(ror_upper, 1e-06)),
        ror_ci_text = sprintf("%.2f (%.2f-%.2f)", ror, ror_lower, ror_upper),
        p_text = format_p_value_compact(p_value),
        soc_label = paste0(
            pattern_group, " (n=", format(target_total_n, big.mark = ",", scientific = FALSE, trim = TRUE), ")\n",
            "[", dominant_soc, "]"
        )
    )

jader_forest_label_levels <- jader_dominant_soc_ror %>%
    dplyr::arrange(pattern_group) %>%
    dplyr::pull(soc_label)

jader_dominant_soc_ror <- jader_dominant_soc_ror %>%
    dplyr::mutate(soc_label = factor(soc_label, levels = rev(jader_forest_label_levels)))

## 2. 绘图：左侧森林图 panel ####
p_jader_forest <- ggplot(
    jader_dominant_soc_ror,
    aes(x = log2_ror, y = soc_label, color = pattern_group)
) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
    geom_errorbarh(
        aes(xmin = log2_lower, xmax = log2_upper),
        height = 0.25,
        linewidth = 0.5
    ) +
    geom_point(aes(fill = pattern_group), shape = 22, size = 2.8, stroke = 0.5, color = "black") +
    scale_color_manual(values = pattern_color_map, guide = "none") +
    scale_fill_manual(values = pattern_color_map, guide = "none") +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) +
    scale_y_discrete(expand = expansion(add = c(0.5, 1.2))) +
    labs(
        x = "log2(ROR)",
        y = NULL
    ) +
    theme_bw() +
    theme(
        axis.title.x = element_text(size = 8, color = "black"),
        axis.text.x = element_text(size = 8, color = "black"),
        axis.text.y = element_text(size = 8, color = "black"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
        plot.margin = margin(t = 5, r = 2, b = 5, l = 5)
    )

## 3. 绘图：右侧文本表 panel ####
p_jader_table <- ggplot(jader_dominant_soc_ror, aes(y = soc_label)) +
    geom_text(
        aes(x = 0.45, label = ror_ci_text),
        hjust = 0.5, size = 8 / ggplot2::.pt, color = "black"
    ) +
    annotate(
        "text", x = 0.45, y = length(jader_forest_label_levels) + 0.7,
        label = "ROR (95% CI)",
        hjust = 0.5, size = 8 / ggplot2::.pt, color = "black", fontface = "bold"
    ) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_discrete(expand = expansion(add = c(0.5, 1.2))) +
    labs(x = " ", y = NULL) +
    theme_bw() +
    theme(
        axis.title.x = element_text(size = 8, color = "transparent"),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.margin = margin(t = 5, r = 0, b = 5, l = 0)
    )

## 4. 横向拼接 ####
p <- patchwork::wrap_plots(p_jader_forest, p_jader_table, nrow = 1, widths = c(1.25, 0.65))

print(p)

## 5. 保存图形 ####
ggsave(
    filename = jader_dominant_soc_ror_plot_file,
    plot = p,
    width = 6,
    height = 4,
    units = "cm",
    device = cairo_pdf
)


# 11. 7图简化版：FAERS 癌种群级主导SOC复合ROR森林图 ####

## 1. 读取并处理绘图数据 ####
faers_dominant_soc_ror_simple <- readxl::read_xlsx(faers_summary_file, sheet = "主导SOC复合ROR") %>%
    dplyr::mutate(
        pattern_group = as.character(pattern_group),
        ror = suppressWarnings(as.numeric(ror)),
        ror_lower = suppressWarnings(as.numeric(ror_lower)),
        ror_upper = suppressWarnings(as.numeric(ror_upper))
    ) %>%
    dplyr::filter(pattern_group %in% pattern_display_order) %>%
    dplyr::mutate(
        pattern_group = factor(pattern_group, levels = pattern_display_order),
        pattern_label = factor(pattern_group, levels = rev(pattern_display_order)),
        log2_ror = log2(pmax(ror, 1e-06)),
        log2_lower = log2(pmax(ror_lower, 1e-06)),
        log2_upper = log2(pmax(ror_upper, 1e-06))
    )

faers_simple_x_min <- min(0, faers_dominant_soc_ror_simple$log2_lower, na.rm = TRUE)
faers_simple_x_max <- max(faers_dominant_soc_ror_simple$log2_upper, na.rm = TRUE)
faers_simple_x_range <- faers_simple_x_max - faers_simple_x_min

## 2. 绘图 ####
p <- ggplot(
    faers_dominant_soc_ror_simple,
    aes(x = log2_ror, y = pattern_label, color = pattern_group)
) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.25) +
    geom_errorbarh(
        aes(xmin = log2_lower, xmax = log2_upper),
        height = 0.22,
        linewidth = 0.35
    ) +
    geom_point(aes(fill = pattern_group), shape = 22, size = 1.8, stroke = 0.35, color = "black") +
    scale_color_manual(values = pattern_color_map, guide = "none") +
    scale_fill_manual(values = pattern_color_map, guide = "none") +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) +
    scale_y_discrete(expand = expansion(add = c(0.4, 0.4))) +
    labs(
        title = "FAERS",
        x = "log2(ROR)",
        y = NULL
    ) +
    coord_cartesian(
        xlim = c(
            faers_simple_x_min - faers_simple_x_range * 0.03,
            faers_simple_x_max + faers_simple_x_range * 0.05
        )
    ) +
    theme_bw(base_family = "Helvetica") +
    theme(
        plot.title = element_text(size = 8, color = "black", face = "bold", hjust = 0.5, family = "Helvetica"),
        axis.title.x = element_text(size = 8, color = "black", family = "Helvetica"),
        axis.text.x = element_text(size = 8, color = "black", family = "Helvetica"),
        axis.text.y = element_text(size = 8, color = "black", family = "Helvetica"),
        axis.ticks.length = grid::unit(1.5, "pt"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
        plot.margin = margin(t = 2, r = 2, b = 2, l = 2)
    )

print(p)

## 3. 保存图形 ####
ggsave(
    filename = faers_dominant_soc_ror_simple_plot_file,
    plot = p,
    width = 6,
    height = 4,
    units = "cm",
    device = cairo_pdf
)


# 12. 8图简化版：JADER 癌种群级主导SOC复合ROR森林图 ####

## 1. 读取并处理绘图数据 ####
jader_dominant_soc_ror_simple <- readxl::read_xlsx(jader_summary_file, sheet = "主导SOC复合ROR") %>%
    dplyr::mutate(
        pattern_group = as.character(pattern_group),
        ror = suppressWarnings(as.numeric(ror)),
        ror_lower = suppressWarnings(as.numeric(ror_lower)),
        ror_upper = suppressWarnings(as.numeric(ror_upper))
    ) %>%
    dplyr::filter(pattern_group %in% pattern_display_order) %>%
    dplyr::mutate(
        pattern_group = factor(pattern_group, levels = pattern_display_order),
        pattern_label = factor(pattern_group, levels = rev(pattern_display_order)),
        log2_ror = log2(pmax(ror, 1e-06)),
        log2_lower = log2(pmax(ror_lower, 1e-06)),
        log2_upper = log2(pmax(ror_upper, 1e-06))
    )

jader_simple_x_min <- min(0, jader_dominant_soc_ror_simple$log2_lower, na.rm = TRUE)
jader_simple_x_max <- max(jader_dominant_soc_ror_simple$log2_upper, na.rm = TRUE)
jader_simple_x_range <- jader_simple_x_max - jader_simple_x_min

## 2. 绘图 ####
p <- ggplot(
    jader_dominant_soc_ror_simple,
    aes(x = log2_ror, y = pattern_label, color = pattern_group)
) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.25) +
    geom_errorbarh(
        aes(xmin = log2_lower, xmax = log2_upper),
        height = 0.22,
        linewidth = 0.35
    ) +
    geom_point(aes(fill = pattern_group), shape = 22, size = 1.8, stroke = 0.35, color = "black") +
    scale_color_manual(values = pattern_color_map, guide = "none") +
    scale_fill_manual(values = pattern_color_map, guide = "none") +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) +
    scale_y_discrete(expand = expansion(add = c(0.4, 0.4))) +
    labs(
        title = "JADER",
        x = "log2(ROR)",
        y = NULL
    ) +
    coord_cartesian(
        xlim = c(
            jader_simple_x_min - jader_simple_x_range * 0.03,
            jader_simple_x_max + jader_simple_x_range * 0.05
        )
    ) +
    theme_bw(base_family = "Helvetica") +
    theme(
        plot.title = element_text(size = 8, color = "black", face = "bold", hjust = 0.5, family = "Helvetica"),
        axis.title.x = element_text(size = 8, color = "black", family = "Helvetica"),
        axis.text.x = element_text(size = 8, color = "black", family = "Helvetica"),
        axis.text.y = element_text(size = 8, color = "black", family = "Helvetica"),
        axis.ticks.length = grid::unit(1.5, "pt"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
        plot.margin = margin(t = 2, r = 2, b = 2, l = 2)
    )

print(p)

## 3. 保存图形 ####
ggsave(
    filename = jader_dominant_soc_ror_simple_plot_file,
    plot = p,
    width = 6,
    height = 4,
    units = "cm",
    device = cairo_pdf
)


# 13. 9图：FAERS+JADER 双库主导SOC复合ROR合并森林图（横向拼接） ####

## 1. 读取并处理两库绘图数据 ####
combined_faers <- readxl::read_xlsx(faers_summary_file, sheet = "主导SOC复合ROR") %>%
    dplyr::mutate(
        pattern_group = as.character(pattern_group),
        dominant_soc = as.character(dominant_soc),
        target_total_n = suppressWarnings(as.numeric(target_total_n)),
        ror = suppressWarnings(as.numeric(ror)),
        ror_lower = suppressWarnings(as.numeric(ror_lower)),
        ror_upper = suppressWarnings(as.numeric(ror_upper)),
        p_value = suppressWarnings(as.numeric(p_value))
    ) %>%
    dplyr::filter(pattern_group %in% pattern_display_order) %>%
    dplyr::mutate(
        pattern_group = factor(pattern_group, levels = pattern_display_order),
        log2_ror = log2(pmax(ror, 1e-06)),
        log2_lower = log2(pmax(ror_lower, 1e-06)),
        log2_upper = log2(pmax(ror_upper, 1e-06)),
        ror_ci_text = sprintf("%.2f (%.2f-%.2f)", ror, ror_lower, ror_upper),
        p_text = format_p_value_compact(p_value),
        soc_label = paste0(
            pattern_group, " (n=", format(target_total_n, big.mark = ",", scientific = FALSE, trim = TRUE), ")\n",
            "[", dominant_soc, "]"
        )
    )

combined_jader <- readxl::read_xlsx(jader_summary_file, sheet = "主导SOC复合ROR") %>%
    dplyr::mutate(
        pattern_group = as.character(pattern_group),
        dominant_soc = as.character(dominant_soc),
        target_total_n = suppressWarnings(as.numeric(target_total_n)),
        ror = suppressWarnings(as.numeric(ror)),
        ror_lower = suppressWarnings(as.numeric(ror_lower)),
        ror_upper = suppressWarnings(as.numeric(ror_upper)),
        p_value = suppressWarnings(as.numeric(p_value))
    ) %>%
    dplyr::filter(pattern_group %in% pattern_display_order) %>%
    dplyr::mutate(
        pattern_group = factor(pattern_group, levels = pattern_display_order),
        log2_ror = log2(pmax(ror, 1e-06)),
        log2_lower = log2(pmax(ror_lower, 1e-06)),
        log2_upper = log2(pmax(ror_upper, 1e-06)),
        ror_ci_text = sprintf("%.2f (%.2f-%.2f)", ror, ror_lower, ror_upper),
        p_text = format_p_value_compact(p_value),
        soc_label = paste0(
            pattern_group, " (n=", format(target_total_n, big.mark = ",", scientific = FALSE, trim = TRUE), ")\n",
            "[", dominant_soc, "]"
        )
    )

combined_faers_levels <- combined_faers %>% dplyr::arrange(pattern_group) %>% dplyr::pull(soc_label)
combined_jader_levels <- combined_jader %>% dplyr::arrange(pattern_group) %>% dplyr::pull(soc_label)

combined_faers <- combined_faers %>%
    dplyr::mutate(soc_label = factor(soc_label, levels = rev(combined_faers_levels)))
combined_jader <- combined_jader %>%
    dplyr::mutate(soc_label = factor(soc_label, levels = rev(combined_jader_levels)))

combined_faers_xmax <- max(combined_faers$log2_upper, na.rm = TRUE)
combined_faers_xmin <- min(combined_faers$log2_lower, na.rm = TRUE)
combined_faers_xrng <- combined_faers_xmax - combined_faers_xmin
combined_faers_text_ror <- combined_faers_xmax + combined_faers_xrng * 0.10
combined_faers_text_p <- combined_faers_xmax + combined_faers_xrng * 0.55

combined_jader_xmax <- max(combined_jader$log2_upper, na.rm = TRUE)
combined_jader_xmin <- min(combined_jader$log2_lower, na.rm = TRUE)
combined_jader_xrng <- combined_jader_xmax - combined_jader_xmin
combined_jader_text_ror <- combined_jader_xmax + combined_jader_xrng * 0.10
combined_jader_text_p <- combined_jader_xmax + combined_jader_xrng * 0.55

## 2. 绘图：FAERS 子图 ####
p_faers <- ggplot(
    combined_faers,
    aes(x = log2_ror, y = soc_label, color = pattern_group)
) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
    geom_errorbarh(
        aes(xmin = log2_lower, xmax = log2_upper),
        height = 0.25,
        linewidth = 0.6
    ) +
    geom_point(aes(fill = pattern_group), shape = 22, size = 3.8, stroke = 0.6, color = "black") +
    annotate(
        "text",
        x = combined_faers_text_ror, y = nrow(combined_faers) + 0.7,
        label = "ROR (95% CI)",
        hjust = 0, size = 11 / ggplot2::.pt, color = "black", fontface = "bold"
    ) +
    geom_text(
        aes(x = combined_faers_text_ror, label = ror_ci_text),
        hjust = 0, size = 11 / ggplot2::.pt, color = "black"
    ) +
    scale_color_manual(values = pattern_color_map, guide = "none") +
    scale_fill_manual(values = pattern_color_map, guide = "none") +
    scale_x_continuous(
        breaks = scales::pretty_breaks(n = 5),
        expand = c(0, 0)
    ) +
    scale_y_discrete(expand = expansion(add = c(0.5, 1.2))) +
    labs(
        x = "log2(ROR) — FAERS",
        y = NULL
    ) +
    theme_bw() +
    theme(
        axis.title.x = element_text(size = 13, color = "black"),
        axis.text.x = element_text(size = 11, color = "black"),
        axis.text.y = element_text(size = 11, color = "black"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
        plot.margin = margin(t = 5, r = 90, b = 5, l = 5)
    ) +
    coord_cartesian(
        xlim = c(combined_faers_xmin - combined_faers_xrng * 0.05, combined_faers_xmax + combined_faers_xrng * 0.05),
        clip = "off"
    )

## 3. 绘图：JADER 子图 ####
p_jader <- ggplot(
    combined_jader,
    aes(x = log2_ror, y = soc_label, color = pattern_group)
) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
    geom_errorbarh(
        aes(xmin = log2_lower, xmax = log2_upper),
        height = 0.25,
        linewidth = 0.6
    ) +
    geom_point(aes(fill = pattern_group), shape = 22, size = 3.8, stroke = 0.6, color = "black") +
    annotate(
        "text",
        x = combined_jader_text_ror, y = nrow(combined_jader) + 0.7,
        label = "ROR (95% CI)",
        hjust = 0, size = 11 / ggplot2::.pt, color = "black", fontface = "bold"
    ) +
    geom_text(
        aes(x = combined_jader_text_ror, label = ror_ci_text),
        hjust = 0, size = 11 / ggplot2::.pt, color = "black"
    ) +
    scale_color_manual(values = pattern_color_map, guide = "none") +
    scale_fill_manual(values = pattern_color_map, guide = "none") +
    scale_x_continuous(
        breaks = scales::pretty_breaks(n = 5),
        expand = c(0, 0)
    ) +
    scale_y_discrete(expand = expansion(add = c(0.5, 1.2))) +
    labs(
        x = "log2(ROR) — JADER",
        y = NULL
    ) +
    theme_bw() +
    theme(
        axis.title.x = element_text(size = 13, color = "black"),
        axis.text.x = element_text(size = 11, color = "black"),
        axis.text.y = element_text(size = 11, color = "black"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
        plot.margin = margin(t = 5, r = 90, b = 5, l = 15)
    ) +
    coord_cartesian(
        xlim = c(combined_jader_xmin - combined_jader_xrng * 0.05, combined_jader_xmax + combined_jader_xrng * 0.05),
        clip = "off"
    )

## 4. 横向拼接 ####
p <- patchwork::wrap_plots(p_faers, p_jader, nrow = 1, widths = c(1, 1))

print(p)

## 5. 保存图形 ####
ggsave(
    filename = combined_dominant_soc_ror_plot_file,
    plot = p,
    width = 18,
    height = 4.8,
    device = cairo_pdf
)

message("========== 癌种群级不成比例与多系统共发绘图完成 ==========")
