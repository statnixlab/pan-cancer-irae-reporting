# FAERS实体瘤irAE发生谱分析####

rm(list = ls())
gc()

source("Code/config.R")
source("Code/mapping_rules.R")
source("Code/irae_reporting_summary.R")

script_name <- "4.2FAERS实体瘤irAE发生谱分析"
script_results_path <- file.path(results_path, script_name)
if (!dir.exists(script_results_path)) {
    dir.create(script_results_path, recursive = TRUE, showWarnings = FALSE)
}

message("========== FAERS实体瘤irAE发生谱分析开始 ==========")
message("结果目录: ", script_results_path)


## 辅助函数####

### 导出ComplexHeatmap对象为PDF，宽高在调用处单独指定####
save_complex_heatmap_pdf <- function(heatmap_object, output_file, width = 11, height = 8) {
    grDevices::cairo_pdf(output_file, width = width, height = height)
    grid::grid.newpage()
    ComplexHeatmap::draw(heatmap_object, newpage = FALSE)
    invisible(grDevices::dev.off())
}


## 参数设置####

# 实体瘤识别已在faers_solid上游完成，本脚本仅构建ICI使用病例主分析队列
top_pt_for_main_plot <- 5
top_pt_for_bubble_plot <- 2
top_soc_for_annual_trend <- 5
annual_year_min <- 2015
annual_year_max <- NULL
compare_with_existing_primary_organ_summary <- TRUE


## 1. 加载数据与字典####

faers_solid <- load_named_rdata(
    file.path(output_path, "faers_solid.RData"),
    object_name = "faers_solid"
)

irae_dictionary <- load_irae_dictionary(
    file.path(data_path, "irAE列表.xlsx")
)

message("已加载 faers_solid")
message("总样本数: ", nrow(faers_solid))
message("变量数: ", ncol(faers_solid))
message("irAE字典条目数: ", nrow(irae_dictionary))


## 2. 生成主分析队列####

faers_solid_ici_reporting <- faers_solid %>%
    dplyr::filter(has_ici == TRUE)

filter_summary <- tibble::tibble(
    step = c("faers_solid", "仅保留ICI使用病例"),
    n = c(nrow(faers_solid), nrow(faers_solid_ici_reporting))
) %>%
    dplyr::mutate(retained_pct = round(n / dplyr::first(n) * 100, 2))

message("\n主分析队列构建流程:")
print(filter_summary)

message("\n主分析队列样本数: ", nrow(faers_solid_ici_reporting))
message("原发器官数: ", dplyr::n_distinct(faers_solid_ici_reporting$primary_organ))


## 3. 构建发生谱汇总表####

cohort_summary <- build_cohort_summary(
    data = faers_solid_ici_reporting,
    excluded_organs = character(0)
) %>%
    dplyr::filter(item != "排除原发器官") %>%
    dplyr::mutate(
        value = dplyr::case_when(
            section == "队列定义" & item == "数据源" ~ "faers_solid",
            section == "队列定义" & item == "纳入病例" ~ "仅保留 has_ici == TRUE 的 faers_solid 病例",
            TRUE ~ value
        )
    )

soc_count_distribution <- summarise_soc_count_distribution(
    faers_solid_ici_reporting
)

any_irae_by_organ <- summarise_any_irae_by_organ(
    faers_solid_ici_reporting
)

soc_overall <- summarise_soc_overall(
    faers_solid_ici_reporting
)

soc_by_organ <- summarise_soc_by_organ(
    faers_solid_ici_reporting
)

soc_top3_by_organ <- build_top_soc_by_organ(
    soc_by_organ_df = soc_by_organ,
    top_n = 3
)

pt_metadata <- extract_pt_metadata(
    data = faers_solid_ici_reporting,
    irae_dictionary = irae_dictionary
)

pt_overall <- summarise_pt_overall(
    data = faers_solid_ici_reporting,
    irae_dictionary = irae_dictionary
)

pt_top5_by_soc <- select_top_pt_by_soc(
    pt_summary = pt_overall,
    top_n = top_pt_for_main_plot
)

pt_top2_by_soc_for_plot <- select_top_pt_by_soc(
    pt_summary = pt_overall,
    top_n = top_pt_for_bubble_plot
)

pt_by_organ_selected <- summarise_pt_by_organ(
    data = faers_solid_ici_reporting,
    selected_pt_df = pt_top2_by_soc_for_plot
)

soc_cooccurrence <- build_soc_cooccurrence_table(
    faers_solid_ici_reporting
)

soc_upset_matrix <- build_soc_upset_matrix(
    faers_solid_ici_reporting
)

annual_trend <- build_annual_reporting_trend(
    data = faers_solid_ici_reporting,
    soc_overall = soc_overall,
    year_min = annual_year_min,
    year_max = annual_year_max,
    top_n_soc = top_soc_for_annual_trend
)

annual_trend_table <- annual_trend$combined

message("\n各器官Any-irAE汇总（前10）:")
print(head(any_irae_by_organ, 10))

message("\n整体SOC汇总:")
print(soc_overall)

message("\n整体PT汇总（前15）:")
print(head(pt_overall, 15))


## 4. 数据校验####

compare_with_existing_primary_organ_summary_flag <- if (exists(
    "compare_with_existing_primary_organ_summary",
    inherits = FALSE
)) {
    isTRUE(compare_with_existing_primary_organ_summary)
} else {
    message(
        "未检测到 compare_with_existing_primary_organ_summary，",
        "按默认值 TRUE 进行原发器官汇总结果比对"
    )
    TRUE
}

existing_primary_organ_summary_path <- if (compare_with_existing_primary_organ_summary_flag) {
    file.path(
        results_path,
        "3.1primary_organ不成比例分析",
        "FAERS_primary_organ_SOC_summary.xlsx"
    )
} else {
    NULL
}

validation_results <- validate_reporting_summary(
    filtered_data = faers_solid_ici_reporting,
    any_irae_by_organ = any_irae_by_organ,
    soc_overall = soc_overall,
    soc_by_organ = soc_by_organ,
    pt_overall = pt_overall,
    pt_by_organ_selected = pt_by_organ_selected,
    pt_metadata = pt_metadata,
    existing_summary_path = existing_primary_organ_summary_path
)

message("\n数据校验结果:")
print(validation_results)

if (any(validation_results$status == "FAIL")) {
    stop("发生谱分析校验失败，请根据数据校验结果排查问题")
}


## 5. 保存结果表与对象####

summary_output_path <- file.path(
    script_results_path,
    "FAERS实体瘤irAE发生谱汇总.xlsx"
)

writexl::write_xlsx(
    list(
        "队列汇总" = cohort_summary,
        "每例SOC数分布" = soc_count_distribution,
        "器官_any_irAE" = any_irae_by_organ,
        "SOC_总体" = soc_overall,
        "SOC_器官分层" = soc_by_organ,
        "SOC_top3_by_organ" = soc_top3_by_organ,
        "PT_总体" = pt_overall,
        "PT_top5_by_SOC" = pt_top5_by_soc,
        "PT_top2_by_SOC_for_plot" = pt_top2_by_soc_for_plot,
        "PT_器官分层_selected" = pt_by_organ_selected,
        "SOC_cooccurrence" = soc_cooccurrence,
        "年度趋势" = annual_trend_table,
        "数据校验" = validation_results
    ),
    path = summary_output_path
)

analysis_metadata <- list(
    script_name = script_name,
    script_results_path = script_results_path,
    top_pt_for_main_plot = top_pt_for_main_plot,
    top_pt_for_bubble_plot = top_pt_for_bubble_plot,
    top_soc_for_annual_trend = top_soc_for_annual_trend,
    annual_year_min = annual_year_min,
    annual_year_max = annual_year_max
)

save(
    analysis_metadata,
    filter_summary,
    cohort_summary,
    soc_count_distribution,
    any_irae_by_organ,
    soc_overall,
    soc_by_organ,
    soc_top3_by_organ,
    pt_metadata,
    pt_overall,
    pt_top5_by_soc,
    pt_top2_by_soc_for_plot,
    pt_by_organ_selected,
    soc_cooccurrence,
    soc_upset_matrix,
    annual_trend,
    validation_results,
    file = file.path(script_results_path, "FAERS实体瘤irAE发生谱对象.RData")
)

message("\n分析结果文件已保存:")
message("  汇总工作簿: ", summary_output_path)
message("  对象文件: ", file.path(script_results_path, "FAERS实体瘤irAE发生谱对象.RData"))


## 6. 图1：各器官Any-irAE lollipop图####

### 图1：数据准备####

primary_organ_color_map <- get_primary_organ_color_map(any_irae_by_organ$primary_organ)

any_irae_plot_data <- any_irae_by_organ %>%
    dplyr::arrange(desc(N), primary_organ) %>%
    dplyr::mutate(
        organ_label = primary_organ,
        organ_label = factor(organ_label, levels = rev(unique(organ_label))),
        label_text = paste0(format_case_n(n_any), "/", format_case_n(N))
    )

### 图1：绘图代码####

p_any_irae_by_organ <- ggplot(any_irae_plot_data, aes(x = RP, y = organ_label)) +
    geom_segment(
        aes(x = 0, xend = RP, yend = organ_label),
        linewidth = 0.9,
        colour = "#B8B8B8"
    ) +
    geom_point(
        aes(fill = primary_organ),
        shape = 21,
        size = 3.8,
        stroke = 0.35,
        colour = "#2B2B2B"
    ) +
    geom_text(
        aes(label = label_text),
        hjust = -0.1,
        size = 12 / ggplot2::.pt,
        colour = "black"
    ) +
    scale_fill_manual(values = primary_organ_color_map, guide = "none") +
    scale_x_continuous(
        limits = c(0, max(any_irae_plot_data$RP, na.rm = TRUE) * 1.28),
        labels = function(x) paste0(round(x, 1), "%"),
        expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
        x = "Any-irAE reporting proportion",
        y = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
        text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
        axis.line.x = element_line(colour = "black"),
        axis.ticks.y = element_blank()
    )

print(p_any_irae_by_organ)

### 图1：导出PDF####

ggsave(
    filename = file.path(script_results_path, "FAERS实体瘤_Any_irAE_lollipop.pdf"),
    plot = p_any_irae_by_organ,
    width = 9,
    height = 8,
    device = cairo_pdf
)


## 7. 图2：SOC气泡热图####

### 图2：数据准备####

soc_bubble_organ_order <- soc_by_organ %>%
    dplyr::distinct(primary_organ, N_organ) %>%
    dplyr::arrange(desc(N_organ), primary_organ) %>%
    dplyr::mutate(organ_label = primary_organ)

soc_bubble_levels <- soc_by_organ %>%
    dplyr::group_by(SOC) %>%
    dplyr::summarise(
        overall_n = sum(n),
        overall_RP = sum(n) / sum(N_organ) * 100,
        .groups = "drop"
    ) %>%
    dplyr::arrange(desc(overall_RP), desc(overall_n), SOC) %>%
    dplyr::pull(SOC)

soc_bubble_plot_data <- soc_by_organ %>%
    dplyr::left_join(soc_bubble_organ_order, by = c("primary_organ", "N_organ")) %>%
    dplyr::mutate(
        organ_label = factor(organ_label, levels = rev(soc_bubble_organ_order$organ_label)),
        SOC = factor(SOC, levels = soc_bubble_levels)
    )

### 图2：绘图代码####

p_soc_bubble_heatmap <- ggplot(soc_bubble_plot_data, aes(x = SOC, y = organ_label)) +
    geom_point(
        aes(size = n, fill = RP),
        shape = 21,
        colour = "#2B2B2B",
        stroke = 0.3
    ) +
    scale_size_area(
        max_size = 9,
        breaks = scales::pretty_breaks(n = 4)
    ) +
    scale_fill_gradientn(
        colours = c("#FFF4E6", "#F8C27A", "#E67E4A", "#A64B2A"),
        name = "RP (%)"
    ) +
    labs(
        x = "SOC",
        y = NULL,
        size = "Report count"
    ) +
    theme_bw(base_size = 12) +
    theme(
        text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        legend.title = element_text(size = 14, colour = "black"),
        legend.text = element_text(size = 14, colour = "black"),
        panel.grid.major = element_line(colour = "#E7E7E7", linewidth = 0.25),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
        axis.text.x = element_text(size = 14, colour = "black", angle = 45, hjust = 1),
        legend.position = "right"
    )

print(p_soc_bubble_heatmap)

### 图2：导出PDF####

ggsave(
    filename = file.path(script_results_path, "FAERS实体瘤_SOC气泡热图.pdf"),
    plot = p_soc_bubble_heatmap,
    width = 9,
    height = 8,
    device = cairo_pdf
)


## 8. 图3：PT分面lollipop图####

### 图3：数据准备####

pt_top5_plot_data <- pt_top5_by_soc %>%
    dplyr::arrange(SOC, desc(RP), desc(n), PT) %>%
    dplyr::group_by(SOC) %>%
    dplyr::mutate(
        pt_key = paste0(PT, "___", SOC),
        label_offset = dplyr::if_else(
            max(RP, na.rm = TRUE) > 0,
            max(RP, na.rm = TRUE) * 0.06,
            0.05
        ),
        label_x = RP + label_offset
    ) %>%
    dplyr::ungroup()

pt_top5_levels <- pt_top5_plot_data %>%
    dplyr::arrange(SOC, RP, n, PT) %>%
    dplyr::pull(pt_key) %>%
    unique()

pt_top5_plot_data <- pt_top5_plot_data %>%
    dplyr::mutate(
        pt_key = factor(pt_key, levels = pt_top5_levels),
        label_text = paste0("n=", format_case_n(n))
    )

### 图3：绘图代码####

p_pt_top5_lollipop <- ggplot(pt_top5_plot_data, aes(x = RP, y = pt_key)) +
    geom_segment(
        aes(x = 0, xend = RP, yend = pt_key),
        linewidth = 0.75,
        colour = "#BCBCBC"
    ) +
    geom_point(
        colour = "#8C2D04",
        fill = "#F16913",
        shape = 21,
        size = 3.2,
        stroke = 0.25
    ) +
    geom_text(
        aes(x = label_x, label = label_text),
        hjust = 0,
        size = 12 / ggplot2::.pt,
        colour = "black"
    ) +
    facet_wrap(~ SOC, scales = "free", ncol = 3) +
    scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
    scale_x_continuous(
        labels = function(x) paste0(round(x, 1), "%"),
        expand = expansion(mult = c(0, 0.45))
    ) +
    coord_cartesian(clip = "off") +
    labs(
        x = "Overall PT reporting proportion",
        y = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
        text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "#F4E6D7", colour = "#D9B38C"),
        strip.text = element_text(size = 12, colour = "black", face = "bold"),
        plot.margin = margin(t = 5.5, r = 14, b = 5.5, l = 5.5)
    )

print(p_pt_top5_lollipop)

### 图3：导出PDF####

ggsave(
    filename = file.path(script_results_path, "FAERS实体瘤_PT分面lollipop.pdf"),
    plot = p_pt_top5_lollipop,
    width = 18,
    height = 8,
    device = cairo_pdf
)


## 9. 图4：selected PT器官分层气泡图####

### 图4：数据准备####

pt_bubble_organ_order <- pt_by_organ_selected %>%
    dplyr::distinct(primary_organ, N_organ) %>%
    dplyr::arrange(desc(N_organ), primary_organ) %>%
    dplyr::mutate(organ_label = primary_organ)

pt_bubble_levels <- pt_by_organ_selected %>%
    dplyr::arrange(SOC, rank_within_soc, PT) %>%
    dplyr::mutate(pt_key = paste0(PT, "___", SOC)) %>%
    dplyr::pull(pt_key) %>%
    unique()

pt_bubble_plot_data <- pt_by_organ_selected %>%
    dplyr::left_join(pt_bubble_organ_order, by = c("primary_organ", "N_organ")) %>%
    dplyr::mutate(
        organ_label = factor(organ_label, levels = rev(pt_bubble_organ_order$organ_label)),
        pt_key = factor(paste0(PT, "___", SOC), levels = pt_bubble_levels)
    )

### 图4：绘图代码####

p_pt_bubble_by_organ <- ggplot(pt_bubble_plot_data, aes(x = pt_key, y = organ_label)) +
    geom_point(
        aes(size = n, fill = RP),
        shape = 21,
        colour = "#2B2B2B",
        stroke = 0.25
    ) +
    facet_wrap(~ SOC, scales = "free_x", ncol = 3) +
    scale_x_discrete(labels = function(x) sub("___.*$", "", x)) +
    scale_size_area(max_size = 10) +
    scale_fill_gradientn(
        colours = c("#F4FBF6", "#B8E0C2", "#59A96A", "#1F6B45"),
        name = "RP (%)"
    ) +
    labs(
        x = "PT",
        y = NULL,
        size = "Report count"
    ) +
    theme_bw(base_size = 12) +
    theme(
        text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        legend.title = element_text(size = 14, colour = "black"),
        legend.text = element_text(size = 14, colour = "black"),
        panel.grid.major = element_line(colour = "#E7E7E7", linewidth = 0.25),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 14, colour = "black", angle = 45, hjust = 1),
        strip.background = element_rect(fill = "#E8F2E4", colour = "#B9D3B0"),
        strip.text = element_text(size = 12, colour = "black", face = "bold")
    )

print(p_pt_bubble_by_organ)

### 图4：导出PDF####

ggsave(
    filename = file.path(script_results_path, "FAERS实体瘤_selected_PT器官分层气泡图.pdf"),
    plot = p_pt_bubble_by_organ,
    width = 14,
    height = 10,
    device = cairo_pdf
)


## 10. 图5：年度趋势组合图####

### 图5：数据准备####

annual_any_irae_plot_data <- annual_trend$any_irae
annual_top_soc_plot_data <- annual_trend$top_soc
annual_top_soc_color_map <- c(
    "SKIN" = "#D7301F",
    "ENDO" = "#238B45",
    "GI" = "#DD8D29",
    "NEURO" = "#2171B5",
    "MSK" = "#6A51A3",
    "HEPA" = "#7F2704",
    "HEME" = "#54278F",
    "RESP" = "#08519C",
    "RENU" = "#006D2C",
    "CARD" = "#B30000",
    "EYE" = "#756BB1"
)

### 图5：绘图代码####

if (nrow(annual_any_irae_plot_data) > 0 && nrow(annual_top_soc_plot_data) > 0) {
    p_annual_reporting_trend_any <- ggplot(annual_any_irae_plot_data, aes(x = year, y = RP)) +
        geom_line(linewidth = 0.9, colour = "#A63603") +
        geom_point(size = 2.4, colour = "#A63603") +
        geom_text(
            aes(label = paste0("n=", format_case_n(n))),
            vjust = -0.8,
            size = 12 / ggplot2::.pt,
            colour = "black"
        ) +
        scale_x_continuous(breaks = annual_any_irae_plot_data$year) +
        scale_y_continuous(labels = function(x) paste0(round(x, 1), "%")) +
        labs(
            x = NULL,
            y = "Any-irAE RP",
            title = "Annual Any-irAE reporting proportion"
        ) +
        theme_bw(base_size = 12) +
        theme(
            text = element_text(size = 12, colour = "black"),
            axis.title = element_text(size = 14, colour = "black"),
            axis.text = element_text(size = 14, colour = "black"),
            plot.title = element_text(size = 12, colour = "black", face = "bold"),
            panel.grid.minor = element_blank()
        )

    p_annual_reporting_trend_soc <- ggplot(
        annual_top_soc_plot_data,
        aes(x = year, y = RP, colour = series)
    ) +
        geom_line(linewidth = 0.85) +
        geom_point(size = 2.1) +
        scale_x_continuous(breaks = sort(unique(annual_top_soc_plot_data$year))) +
        scale_y_continuous(labels = function(x) paste0(round(x, 1), "%")) +
        scale_colour_manual(
            values = annual_top_soc_color_map,
            name = "SOC"
        ) +
        labs(
            x = "Year",
            y = "SOC RP",
            title = "Annual reporting proportion of overall top SOCs"
        ) +
        theme_bw(base_size = 12) +
        theme(
            text = element_text(size = 12, colour = "black"),
            axis.title = element_text(size = 14, colour = "black"),
            axis.text = element_text(size = 14, colour = "black"),
            legend.title = element_text(size = 14, colour = "black"),
            legend.text = element_text(size = 14, colour = "black"),
            plot.title = element_text(size = 12, colour = "black", face = "bold"),
            panel.grid.minor = element_blank(),
            legend.position = "bottom"
        )

    p_annual_reporting_trend <- p_annual_reporting_trend_any /
        p_annual_reporting_trend_soc +
        patchwork::plot_layout(heights = c(1, 1.2))
} else {
    p_annual_reporting_trend_any <- NULL
    p_annual_reporting_trend_soc <- NULL
    p_annual_reporting_trend <- ggplot() +
        annotate(
            "text",
            x = 0,
            y = 0,
            label = "No annual trend data",
            size = 12 / ggplot2::.pt,
            colour = "black"
        ) +
        theme_void() +
        theme(text = element_text(size = 12, colour = "black"))
}

print(p_annual_reporting_trend)

### 图5：导出PDF####

ggsave(
    filename = file.path(script_results_path, "FAERS实体瘤_年度趋势组合图.pdf"),
    plot = p_annual_reporting_trend,
    width = 12,
    height = 10,
    device = cairo_pdf
)


## 11. 图6：SOC共报告UpSet图####

### 图6：数据准备####

soc_upset_top_n_comb <- 20
if (is.null(soc_upset_matrix) || nrow(soc_upset_matrix) == 0 || ncol(soc_upset_matrix) == 0) {
    soc_upset_comb_matrix <- NULL
} else {
    soc_upset_comb_matrix <- ComplexHeatmap::make_comb_mat(soc_upset_matrix)
    soc_upset_comb_order <- order(ComplexHeatmap::comb_size(soc_upset_comb_matrix), decreasing = TRUE)

    if (length(soc_upset_comb_order) > soc_upset_top_n_comb) {
        soc_upset_comb_order <- soc_upset_comb_order[seq_len(soc_upset_top_n_comb)]
    }

    soc_upset_comb_matrix <- soc_upset_comb_matrix[, soc_upset_comb_order]
}

### 图6：绘图代码####

if (is.null(soc_upset_comb_matrix)) {
    p_soc_upset <- NULL
} else {
    p_soc_upset <- ComplexHeatmap::UpSet(
        soc_upset_comb_matrix,
        top_annotation = ComplexHeatmap::upset_top_annotation(
            soc_upset_comb_matrix,
            add_numbers = TRUE
        ),
        right_annotation = ComplexHeatmap::upset_right_annotation(soc_upset_comb_matrix),
        pt_size = grid::unit(2, "mm"),
        lwd = 1
    )
}

if (!is.null(p_soc_upset)) {
    ComplexHeatmap::draw(p_soc_upset)
}

### 图6：导出PDF####

if (!is.null(p_soc_upset)) {
    save_complex_heatmap_pdf(
        heatmap_object = p_soc_upset,
        output_file = file.path(script_results_path, "FAERS实体瘤_SOC共报告UpSet图.pdf"),
        width = 12,
        height = 8
    )
}


## 12. 保存图形对象与绘图清单####

### 保存绘图清单####

plot_manifest <- tibble::tribble(
    ~plot_name, ~object_name, ~output_file, ~recommended_width, ~recommended_height, ~note,
    "各器官Any-irAE lollipop", "p_any_irae_by_organ", file.path(script_results_path, "FAERS实体瘤_Any_irAE_lollipop.pdf"), 11, 8, "主图1",
    "SOC气泡热图", "p_soc_bubble_heatmap", file.path(script_results_path, "FAERS实体瘤_SOC气泡热图.pdf"), 13, 9, "主图2",
    "PT分面lollipop", "p_pt_top5_lollipop", file.path(script_results_path, "FAERS实体瘤_PT分面lollipop.pdf"), 13, 10, "主图3",
    "selected PT器官分层气泡图", "p_pt_bubble_by_organ", file.path(script_results_path, "FAERS实体瘤_selected_PT器官分层气泡图.pdf"), 14, 10, "主图4",
    "年度趋势组合图", "p_annual_reporting_trend", file.path(script_results_path, "FAERS实体瘤_年度趋势组合图.pdf"), 12, 10, "补充图",
    "SOC共报告UpSet图", "p_soc_upset", file.path(script_results_path, "FAERS实体瘤_SOC共报告UpSet图.pdf"), 12, 8, "补充图"
)

### 保存绘图对象与清单文件####

save(
    any_irae_plot_data,
    soc_bubble_plot_data,
    pt_top5_plot_data,
    pt_bubble_plot_data,
    annual_any_irae_plot_data,
    annual_top_soc_plot_data,
    annual_top_soc_color_map,
    soc_upset_comb_matrix,
    p_any_irae_by_organ,
    p_soc_bubble_heatmap,
    p_pt_top5_lollipop,
    p_pt_bubble_by_organ,
    p_annual_reporting_trend_any,
    p_annual_reporting_trend_soc,
    p_annual_reporting_trend,
    p_soc_upset,
    plot_manifest,
    file = file.path(script_results_path, "FAERS实体瘤irAE发生谱绘图对象.RData")
)

writexl::write_xlsx(
    list(
        "绘图清单" = plot_manifest
    ),
    path = file.path(script_results_path, "FAERS实体瘤irAE发生谱绘图清单.xlsx")
)

message("\n绘图对象已保存:")
message("  绘图对象文件: ", file.path(script_results_path, "FAERS实体瘤irAE发生谱绘图对象.RData"))
message("  绘图清单: ", file.path(script_results_path, "FAERS实体瘤irAE发生谱绘图清单.xlsx"))
message("  PDF已在各图对应的导出区分别保存")
message("  如需调整尺寸，直接修改各图下面对应的 ggsave()/save_complex_heatmap_pdf() 调用")
message("\n========== FAERS实体瘤irAE发生谱分析与绘图完成 ==========")
