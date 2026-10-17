# ICI类别分层数据稀疏性绘图（返修：审稿人1第12条） ####
# 读取6.2.9保存的Excel结果绘图，不依赖分析阶段的内存对象
# 四个面板分别保存，供拼图为2×2组合图

rm(list = ls())
gc()

# 1. 加载配置 ####
source("Code/config.R")

# 2. 结果目录与输入文件 ####
results_dir_name <- "6.2.10ICI类别分层数据稀疏性绘图"
script_results_path <- file.path(results_path, results_dir_name)

if (!dir.exists(script_results_path)) {
    dir.create(script_results_path, recursive = TRUE, showWarnings = FALSE)
}

input_xlsx_file <- file.path(results_path, "6.2.9ICI类别分层数据稀疏性分析", "ICI类别分层数据稀疏性分析结果.xlsx")

if (!file.exists(input_xlsx_file)) {
    stop("6.2.9结果文件不存在，请先运行 Code/6.2.9ICI类别分层数据稀疏性分析.R: ", input_xlsx_file)
}

message("========== ICI类别分层数据稀疏性绘图开始 ==========")


# a图：Events per cancer-SOC cell ####

## 1. 读取并处理绘图数据 ####
events_plot_df <- readxl::read_xlsx(input_xlsx_file, sheet = "单元明细") %>%
    dplyr::mutate(
        dataset = factor(
            dataset,
            levels = c("All ICI (main analysis)", "anti-PD-1", "anti-PD-L1", "anti-CTLA-4+anti-PD-1"),
            labels = c("All ICI", "Anti-PD-1", "Anti-PD-L1", "Anti-CTLA-4\n+ anti-PD-1")
        )
    )

events_median_df <- events_plot_df %>%
    dplyr::group_by(dataset) %>%
    dplyr::summarise(median_events = stats::median(n_events), .groups = "drop") %>%
    dplyr::mutate(median_label = paste0("Median ", median_events))

## 2. 绘图 ####
p <- ggplot(events_plot_df, aes(x = dataset, y = n_events, fill = dataset)) +
    geom_hline(yintercept = 5, linetype = "dashed", colour = "grey40", linewidth = 0.3) +
    geom_boxplot(width = 0.6, outlier.size = 0.3, outlier.colour = "grey40", linewidth = 0.25) +
    geom_text(
        data = events_median_df,
        aes(x = dataset, y = Inf, label = median_label),
        vjust = 1.5, size = 7 / .pt, colour = "black", inherit.aes = FALSE
    ) +
    scale_y_continuous(
        trans = scales::pseudo_log_trans(base = 10),
        breaks = c(0, 1, 5, 20, 100, 1000, 10000),
        labels = scales::comma,
        expand = expansion(mult = c(0.03, 0.12))
    ) +
    scale_fill_manual(values = c("All ICI" = "grey55", "Anti-PD-1" = "#2166AC", "Anti-PD-L1" = "#6BAED6", "Anti-CTLA-4\n+ anti-PD-1" = "#C6DBEF"), guide = "none") +
    labs(
        x = NULL,
        y = "Events per cancer-SOC cell"
    ) +
    theme_bw(base_size = 7) +
    theme(
        text = element_text(size = 7, colour = "black"),
        axis.text = element_text(size = 6, colour = "black"),
        axis.text.x = element_text(lineheight = 0.9),
        panel.grid = element_blank()
    )

print(p)

## 3. 保存图形 ####
ggsave(
    filename = file.path(script_results_path, "图a_各数据集单元事件数.pdf"),
    plot = p,
    width = 7.3,
    height = 5.5,
    units = "cm",
    device = cairo_pdf
)


# b图：Width of the 95% confidence interval of the ROR ####

## 1. 读取并处理绘图数据 ####
ci_plot_df <- readxl::read_xlsx(input_xlsx_file, sheet = "单元明细") %>%
    dplyr::mutate(
        dataset = factor(
            dataset,
            levels = c("All ICI (main analysis)", "anti-PD-1", "anti-PD-L1", "anti-CTLA-4+anti-PD-1"),
            labels = c("All ICI", "Anti-PD-1", "Anti-PD-L1", "Anti-CTLA-4\n+ anti-PD-1")
        )
    )

ci_median_df <- ci_plot_df %>%
    dplyr::group_by(dataset) %>%
    dplyr::summarise(median_ci_ratio = stats::median(ci_ratio), .groups = "drop") %>%
    dplyr::mutate(median_label = paste0("Median ", sprintf("%.1f", median_ci_ratio)))

## 2. 绘图 ####
p <- ggplot(ci_plot_df, aes(x = dataset, y = ci_ratio, fill = dataset)) +
    geom_boxplot(width = 0.6, outlier.size = 0.3, outlier.colour = "grey40", linewidth = 0.25) +
    geom_text(
        data = ci_median_df,
        aes(x = dataset, y = Inf, label = median_label),
        vjust = 1.5, size = 7 / .pt, colour = "black", inherit.aes = FALSE
    ) +
    scale_y_log10(
        breaks = c(1, 2, 5, 10, 100, 1000),
        labels = scales::comma,
        expand = expansion(mult = c(0.03, 0.12))
    ) +
    scale_fill_manual(values = c("All ICI" = "grey55", "Anti-PD-1" = "#2166AC", "Anti-PD-L1" = "#6BAED6", "Anti-CTLA-4\n+ anti-PD-1" = "#C6DBEF"), guide = "none") +
    labs(
        x = NULL,
        y = "Upper / lower 95% CI of ROR"
    ) +
    theme_bw(base_size = 7) +
    theme(
        text = element_text(size = 7, colour = "black"),
        axis.text = element_text(size = 6, colour = "black"),
        axis.text.x = element_text(lineheight = 0.9),
        panel.grid = element_blank()
    )

print(p)

## 3. 保存图形 ####
ggsave(
    filename = file.path(script_results_path, "图b_各数据集ROR置信区间宽度.pdf"),
    plot = p,
    width = 7.3,
    height = 5.5,
    units = "cm",
    device = cairo_pdf
)


# c图：Deviation of stratum-specific RORs from the main analysis by number of events ####

## 1. 读取并处理绘图数据 ####
deviation_plot_df <- readxl::read_xlsx(input_xlsx_file, sheet = "与主分析偏离_按事件数") %>%
    dplyr::mutate(
        dataset = factor(
            dataset,
            levels = c("anti-PD-1", "anti-PD-L1", "anti-CTLA-4+anti-PD-1"),
            labels = c("Anti-PD-1", "Anti-PD-L1", "Anti-CTLA-4 + anti-PD-1")
        ),
        event_bin = factor(event_bin, levels = c("0", "1-4", "5-19", ">=20"), labels = c("0", "1–4", "5–19", "≥20"))
    )

## 2. 绘图 ####
p <- ggplot(deviation_plot_df, aes(x = event_bin, y = median_abs_log2_diff_vs_main, colour = dataset, group = dataset)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40", linewidth = 0.3) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 1.4) +
    annotate("text", x = 4.4, y = 1, label = "2-fold", hjust = 1, vjust = -0.5, size = 5 / .pt, colour = "grey30") +
    scale_colour_manual(values = c("Anti-PD-1" = "#2166AC", "Anti-PD-L1" = "#6BAED6", "Anti-CTLA-4 + anti-PD-1" = "#9ECAE1"), name = NULL) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
    labs(
        x = "Events in the stratum-specific cell",
        y = "Median |log2 ROR difference|\nvs main analysis"
    ) +
    theme_bw(base_size = 7) +
    theme(
        text = element_text(size = 7, colour = "black"),
        axis.text = element_text(size = 6, colour = "black"),
        axis.title.y = element_text(lineheight = 0.9),
        legend.text = element_text(size = 6, colour = "black"),
        legend.key.size = unit(3, "mm"),
        legend.position = c(0.98, 0.98),
        legend.justification = c(1, 1),
        legend.background = element_blank(),
        panel.grid = element_blank()
    )

print(p)

## 3. 保存图形 ####
ggsave(
    filename = file.path(script_results_path, "图c_分层ROR与主分析偏离_按事件数.pdf"),
    plot = p,
    width = 7.3,
    height = 5.5,
    units = "cm",
    device = cairo_pdf
)


# d图：Zero-event cells with ROR above 1 ####

## 1. 读取并处理绘图数据 ####
zero_event_plot_df <- readxl::read_xlsx(input_xlsx_file, sheet = "单元明细") %>%
    dplyr::filter(n_events == 0) %>%
    dplyr::group_by(dataset) %>%
    dplyr::summarise(
        n_zero_cells = dplyr::n(),
        n_ror_above_1 = sum(ror > 1),
        .groups = "drop"
    ) %>%
    dplyr::mutate(
        prop_ror_above_1 = n_ror_above_1 / n_zero_cells,
        bar_label = paste0(n_ror_above_1, "/", n_zero_cells),
        dataset = factor(
            dataset,
            levels = c("All ICI (main analysis)", "anti-PD-1", "anti-PD-L1", "anti-CTLA-4+anti-PD-1"),
            labels = c("All ICI", "Anti-PD-1", "Anti-PD-L1", "Anti-CTLA-4\n+ anti-PD-1")
        )
    )

## 2. 绘图 ####
p <- ggplot(zero_event_plot_df, aes(x = dataset, y = prop_ror_above_1, fill = dataset)) +
    geom_col(width = 0.6, colour = "black", linewidth = 0.25) +
    geom_text(aes(label = bar_label), vjust = -0.5, size = 5 / .pt, colour = "black") +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent, expand = expansion(mult = c(0, 0.08))) +
    scale_fill_manual(values = c("All ICI" = "grey55", "Anti-PD-1" = "#2166AC", "Anti-PD-L1" = "#6BAED6", "Anti-CTLA-4\n+ anti-PD-1" = "#C6DBEF"), guide = "none") +
    labs(
        x = NULL,
        y = "Zero-event cells with ROR > 1"
    ) +
    theme_bw(base_size = 7) +
    theme(
        text = element_text(size = 7, colour = "black"),
        axis.text = element_text(size = 6, colour = "black"),
        axis.text.x = element_text(lineheight = 0.9),
        panel.grid = element_blank()
    )

print(p)

## 3. 保存图形 ####
ggsave(
    filename = file.path(script_results_path, "图d_零事件单元ROR大于1比例.pdf"),
    plot = p,
    width = 7.3,
    height = 5.5,
    units = "cm",
    device = cairo_pdf
)

message("已生成图形文件，目录: ", script_results_path)
message("========== ICI类别分层数据稀疏性绘图结束 ==========")
