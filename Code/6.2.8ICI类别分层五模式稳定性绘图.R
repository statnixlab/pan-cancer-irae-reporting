# ICI类别分层五模式稳定性绘图（返修：审稿人1第12条） ####
# 读取6.2.7保存的Excel结果与3.3保存的分层报告数绘图，不依赖分析阶段的内存对象

rm(list = ls())
gc()

# 1. 加载配置 ####
source("Code/config.R")

# 2. 结果目录与输入文件 ####
results_dir_name <- "6.2.8ICI类别分层五模式稳定性绘图"
script_results_path <- file.path(results_path, results_dir_name)

if (!dir.exists(script_results_path)) {
    dir.create(script_results_path, recursive = TRUE, showWarnings = FALSE)
}

input_xlsx_file <- file.path(results_path, "6.2.7ICI类别分层五模式稳定性分析", "ICI类别分层五模式稳定性分析结果.xlsx")
regimen_count_file <- file.path(results_path, "3.3同一ICI_regimen下primary_organ不成比例分析", "FAERS_within_ICI_regimen_primary_organ_SOC_summary.xlsx")

if (!file.exists(input_xlsx_file)) {
    stop("6.2.7结果文件不存在，请先运行 Code/6.2.7ICI类别分层五模式稳定性分析.R: ", input_xlsx_file)
}

message("========== ICI类别分层五模式稳定性绘图开始 ==========")


# 1图：Projected group of each cancer type within ICI classes ####

## 1. 读取并处理绘图数据 ####
reference_group_levels <- c("Ocular", "Hematologic", "Neurologic", "Systemic-like", "Visceral-like")
regimen_levels <- c("anti-PD-1", "anti-PD-L1", "anti-CTLA-4+anti-PD-1")
regimen_column_labels <- c(
    "anti-PD-1" = "Anti-PD-1\nmonotherapy",
    "anti-PD-L1" = "Anti-PD-L1\nmonotherapy",
    "anti-CTLA-4+anti-PD-1" = "Anti-CTLA-4 +\nanti-PD-1"
)

projection_df <- readxl::read_xlsx(input_xlsx_file, sheet = "分层投影癌种分组") %>%
    dplyr::filter(sparse_handling == "sparse_to_1", min_cancer_n == 100)

final_group_df <- readxl::read_xlsx(input_xlsx_file, sheet = "投影方法核查") %>%
    dplyr::select(cancer_type, final_group)

regimen_count_df <- readxl::read_xlsx(regimen_count_file, sheet = "regimen内器官样本量") %>%
    dplyr::filter(primary_organ %in% final_group_df$cancer_type, ici_regimen_raw %in% regimen_levels) %>%
    dplyr::select(ici_regimen = ici_regimen_raw, cancer_type = primary_organ, organ_n)

# 用连续坐标压缩参考色条的占位，并在五个最终分组之间增加留白 ####
row_layout_df <- final_group_df %>%
    dplyr::mutate(final_group = factor(final_group, levels = reference_group_levels)) %>%
    dplyr::arrange(final_group, cancer_type) %>%
    dplyr::mutate(
        group_index = as.integer(final_group),
        row_y = dplyr::n() - dplyr::row_number() + 1 +
            0.45 * (max(group_index) - group_index)
    )
cancer_levels <- row_layout_df$cancer_type

column_layout_df <- tibble::tibble(
    ici_regimen = regimen_levels,
    column_x = c(1.05, 2.40, 3.60)
) %>%
    dplyr::left_join(projection_df %>% dplyr::count(ici_regimen), by = "ici_regimen") %>%
    dplyr::mutate(
        n = dplyr::coalesce(n, 0L),
        column_label = paste0(unname(regimen_column_labels[ici_regimen]), "\n(", n, " cancer types)")
    )

header_top <- max(row_layout_df$row_y)

stratum_tile_df <- tidyr::expand_grid(ici_regimen = regimen_levels, cancer_type = cancer_levels) %>%
    dplyr::left_join(regimen_count_df, by = c("ici_regimen", "cancer_type")) %>%
    dplyr::left_join(
        projection_df %>% dplyr::select(ici_regimen, cancer_type, projected_group, same_as_final),
        by = c("ici_regimen", "cancer_type")
    ) %>%
    dplyr::left_join(row_layout_df %>% dplyr::select(cancer_type, row_y), by = "cancer_type") %>%
    dplyr::left_join(column_layout_df, by = "ici_regimen") %>%
    dplyr::mutate(
        tile_fill = ifelse(is.na(projected_group), "Not included", projected_group),
        tile_label = dplyr::case_when(
            is.na(organ_n) ~ "\u2014",
            is.na(projected_group) ~ paste0("(", format(organ_n, big.mark = ",", trim = TRUE), ")"),
            TRUE ~ format(organ_n, big.mark = ",", trim = TRUE)
        ),
        label_colour = dplyr::case_when(
            is.na(projected_group) ~ "grey45",
            projected_group %in% c("Hematologic", "Visceral-like") ~ "white",
            TRUE ~ "black"
        ),
        changed_group = !is.na(projected_group) & !is.na(same_as_final) & !same_as_final,
        tile_width = 1.08
    )

changed_group_df <- stratum_tile_df %>%
    dplyr::filter(changed_group)

final_strip_df <- row_layout_df %>%
    dplyr::mutate(
        column_x = 0.12,
        tile_fill = as.character(final_group),
        tile_label = "",
        tile_width = 0.22
    )

stratum_plot_df <- dplyr::bind_rows(
    final_strip_df %>% dplyr::select(column_x, row_y, tile_fill, tile_width),
    stratum_tile_df %>% dplyr::select(column_x, row_y, tile_fill, tile_width)
) %>%
    dplyr::mutate(
        tile_fill = factor(tile_fill, levels = c(reference_group_levels, "Not included"))
    )

## 2. 绘图 ####
p <- ggplot(stratum_plot_df, aes(x = column_x, y = row_y)) +
    annotate(
        "rect", xmin = 0.51, xmax = 1.59,
        ymin = header_top + 3.25, ymax = header_top + 4.45,
        fill = "grey95", colour = NA
    ) +
    annotate(
        "rect", xmin = 1.86, xmax = 4.14,
        ymin = header_top + 3.25, ymax = header_top + 4.45,
        fill = "grey80", colour = NA
    ) +
    annotate(
        "text", x = 1.05, y = header_top + 3.85,
        label = "Main stratum", size = 8 / .pt, fontface = "bold", colour = "black"
    ) +
    annotate(
        "text", x = 3, y = header_top + 3.85,
        label = "Sparse strata (not interpreted)", size = 8 / .pt,
        fontface = "bold", colour = "black"
    ) +
    annotate(
        "text", x = 0.12, y = header_top + 1.85,
        label = "Final\ngroup", size = 8 / .pt, lineheight = 1.05, colour = "black"
    ) +
    geom_text(
        data = column_layout_df,
        aes(x = column_x, y = header_top + 1.85, label = column_label),
        inherit.aes = FALSE, size = 8 / .pt, lineheight = 1.05, colour = "black"
    ) +
    annotate(
        "segment", x = 1.725, xend = 1.725,
        y = 0.5, yend = header_top + 4.45,
        colour = "#D7DCE1", linewidth = 0.25
    ) +
    geom_tile(
        aes(fill = tile_fill, width = tile_width), height = 0.94,
        colour = "black", linewidth = 0.2
    ) +
    geom_text(
        data = stratum_tile_df, aes(label = tile_label, colour = label_colour),
        size = 8 / .pt
    ) +
    geom_text(
        data = changed_group_df,
        aes(x = column_x + 0.43, colour = label_colour), label = "×",
        size = 8 / .pt, fontface = "bold"
    ) +
    scale_colour_identity() +
    scale_fill_manual(
        values = c(
            "Ocular" = "#DF7A5E", "Hematologic" = "#A15F88", "Neurologic" = "#82B29A",
            "Systemic-like" = "#F2CC8E", "Visceral-like" = "#BC4749", "Not included" = "white"
        ),
        name = NULL,
        drop = FALSE
    ) +
    scale_x_continuous(limits = c(-0.15, 4.17), breaks = NULL, expand = expansion(mult = 0)) +
    scale_y_continuous(
        breaks = row_layout_df$row_y, labels = row_layout_df$cancer_type,
        limits = c(0.45, header_top + 4.55), expand = expansion(mult = 0)
    ) +
    labs(
        title = NULL,
        x = NULL,
        y = NULL
    ) +
    theme_minimal(base_size = 8) +
    theme(
        text = element_text(size = 8, colour = "black"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 8, colour = "black", margin = margin(r = 4)),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 8, colour = "black", margin = margin(l = 2, r = 0)),
        legend.key.size = unit(2.8, "mm"),
        legend.key.spacing.x = unit(1.5, "mm"),
        legend.key = element_rect(fill = "white", colour = NA),
        legend.position = "bottom",
        legend.location = "plot",
        legend.justification = "center",
        legend.box.just = "center",
        legend.box.spacing = unit(2, "mm"),
        legend.margin = margin(0, 0, 0, 0),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(6, 6, 6, 6)
    ) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(colour = "black", linewidth = 0.2)))

print(p)

## 3. 保存图形 ####
ggsave(
    filename = file.path(script_results_path, "图1_ICI类别分层下各癌种投影分组.pdf"),
    plot = p,
    width = 14.5,
    height = 14.5,
    units = "cm",
    device = cairo_pdf
)

message("已生成图形文件，目录: ", script_results_path)
message("========== ICI类别分层五模式稳定性绘图结束 ==========")
