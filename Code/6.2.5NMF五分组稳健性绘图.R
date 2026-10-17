# NMF五分组稳健性绘图（返修：审稿人2第一段） ####
# 读取6.2.4保存的Excel结果绘图，不依赖分析阶段的内存对象

rm(list = ls())
gc()

# 1. 加载配置 ####
source("Code/config.R")

# 2. 结果目录与输入文件 ####
results_dir_name <- "6.2.5NMF五分组稳健性绘图"
script_results_path <- file.path(results_path, results_dir_name)

if (!dir.exists(script_results_path)) {
    dir.create(script_results_path, recursive = TRUE, showWarnings = FALSE)
}

input_xlsx_file <- file.path(results_path, "6.2.4NMF五分组稳健性分析", "NMF五分组稳健性分析结果.xlsx")

if (!file.exists(input_xlsx_file)) {
    stop("6.2.4结果文件不存在，请先运行 Code/6.2.4NMF五分组稳健性分析.R: ", input_xlsx_file)
}

message("========== NMF五分组稳健性绘图开始 ==========")


# 1图：Agreement with the final five groups across robustness conditions ####

## 1. 读取并处理绘图数据 ####
ari_condition_levels <- c("Random seed", "Profile-preserving transformation", "Inclusion threshold", "Leave one cancer out", "Excluding uncertain cancers", "Profile-altering transformation")

ari_plot_df <- dplyr::bind_rows(
    readxl::read_xlsx(input_xlsx_file, sheet = "稳健性运行汇总"),
    readxl::read_xlsx(input_xlsx_file, sheet = "剔除不确定癌种汇总") %>%
        dplyr::mutate(condition_type = "Excluding uncertain cancers")
) %>%
    dplyr::filter(status == "completed") %>%
    dplyr::mutate(
        condition_type = factor(condition_type, levels = ari_condition_levels),
        x_index = as.integer(condition_type),
        agreement = ifelse(ari_vs_final >= 0.8, "ARI \u2265 0.8", "ARI < 0.8"),
        point_label = dplyr::case_when(
            condition_id == "transform_log2_minmax" ~ "log2 ROR, min-max",
            condition_id == "transform_max_scaled" ~ "ROR / maximum",
            condition_id == "transform_zscore_shifted" ~ "z-score",
            condition_id == "transform_log2_channel" ~ "log2 ROR, two-channel",
            condition_id == "transform_log2_channel_sigonly" ~ "log2 ROR, sig. only",
            condition_id == "transform_rank" ~ "Within-cancer rank",
            condition_id == "exclude_low_silhouette" ~ "5 lowest-silhouette\ntypes removed",
            condition_type == "Leave one cancer out" & ari_vs_final < 1 ~ sub("^Without ", "Without ", condition_label),
            TRUE ~ NA_character_
        ),
        # 相对谱形变换列和剔除不确定癌种列的右侧较空，标签放右侧；其余放左侧 ####
        label_side = ifelse(condition_type %in% c("Profile-preserving transformation", "Excluding uncertain cancers"), "right", "left")
    ) %>%
    # 同一列中ARI相同的点均匀横向排开，替代随机抖动，使标签引线准确指向各点 ####
    dplyr::group_by(condition_type, ari_rounded = round(ari_vs_final, 3)) %>%
    dplyr::mutate(
        n_same = dplyr::n(),
        x_position = x_index + (dplyr::row_number() - (n_same + 1) / 2) * pmin(0.05, 0.5 / n_same)
    ) %>%
    dplyr::ungroup()

ari_axis_labels <- ari_plot_df %>%
    dplyr::count(condition_type, x_index) %>%
    dplyr::arrange(x_index) %>%
    dplyr::mutate(axis_label = paste0(sub(" ", "\n", as.character(condition_type)), "\n(n = ", n, ")"))

## 2. 绘图 ####
p <- ggplot(ari_plot_df, aes(x = x_position, y = ari_vs_final)) +
    annotate("rect", xmin = 5.5, xmax = 6.5, ymin = -Inf, ymax = Inf, fill = "grey93") +
    geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey40", linewidth = 0.3) +
    geom_point(aes(fill = agreement), shape = 21, size = 1.8, colour = "black", stroke = 0.25) +
    ggrepel::geom_text_repel(
        aes(label = point_label, hjust = ifelse(label_side == "left", 1, 0)),
        size = 7 / .pt, colour = "black", bg.color = "white", bg.r = 0.2, na.rm = TRUE,
        nudge_x = ifelse(ari_plot_df$label_side == "left", -0.33, 0.2),
        direction = "y", lineheight = 0.9,
        min.segment.length = 0, segment.size = 0.2, segment.colour = "grey40",
        box.padding = 0.3, max.overlaps = Inf, seed = 123456,
        max.iter = 20000, max.time = 2
    ) +
    scale_fill_manual(values = c("ARI \u2265 0.8" = "#2166AC", "ARI < 0.8" = "#BC4749"), name = NULL) +
    scale_x_continuous(breaks = ari_axis_labels$x_index, labels = ari_axis_labels$axis_label, limits = c(0.45, 6.55), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1.04), breaks = seq(0, 1, 0.2)) +
    labs(
        x = NULL,
        y = "ARI vs final five groups"
    ) +
    theme_bw(base_size = 7) +
    theme(
        text = element_text(size = 7, colour = "black"),
        axis.text = element_text(size = 7, colour = "black"),
        legend.text = element_text(size = 7, colour = "black"),
        legend.key.size = unit(3, "mm"),
        legend.position = "top",
        legend.margin = margin(0, 0, 0, 0),
        panel.grid = element_blank()
    )

print(p)

## 3. 保存图形 ####
ggsave(
    filename = file.path(script_results_path, "图1_各稳健性条件与最终五分组的一致性.pdf"),
    plot = p,
    width = 14.65,
    height = 7,
    units = "cm",
    device = cairo_pdf
)


# 2图：Group assignment of each cancer type across robustness conditions ####

## 1. 读取并处理绘图数据 ####
reference_group_levels <- c("Ocular", "Hematologic", "Neurologic", "Systemic-like", "Visceral-like")
assignment_condition_levels <- c("Final group", "Random seed", "Profile-preserving transformation", "Inclusion threshold", "Leave one cancer out", "Profile-altering transformation")

assignment_raw_df <- readxl::read_xlsx(input_xlsx_file, sheet = "各运行癌种分组") %>%
    dplyr::filter(!is.na(reference_group)) %>%
    dplyr::mutate(
        column_label = dplyr::case_when(
            condition_type == "Random seed" ~ sub("^seed_0?", "", condition_id),
            condition_id == "transform_log2_minmax" ~ "log2 ROR, min-max",
            condition_id == "transform_max_scaled" ~ "ROR / maximum",
            condition_id == "transform_zscore_shifted" ~ "z-score",
            condition_id == "transform_log2_channel" ~ "log2 ROR, two-channel",
            condition_id == "transform_log2_channel_sigonly" ~ "log2 ROR, sig. only",
            condition_id == "transform_rank" ~ "Within-cancer rank",
            condition_type == "Inclusion threshold" ~ sub("^>= ", "≥", sub(" reports$", "", condition_label)),
            condition_type == "Leave one cancer out" ~ sub("^Without ", "", condition_label),
            TRUE ~ condition_label
        )
    )

# 每个检验类别的条件数，用于分面标题 ####
condition_counts <- assignment_raw_df %>%
    dplyr::distinct(condition_type, condition_id) %>%
    dplyr::count(condition_type)

facet_labels <- c(
    "Final group" = "",
    "Random seed" = paste0("Random seed\n(n = ", condition_counts$n[condition_counts$condition_type == "Random seed"], ")"),
    "Profile-preserving transformation" = paste0("Profile-\npre-\nserving\ntransfor-\nmation\n(n = ", condition_counts$n[condition_counts$condition_type == "Profile-preserving transformation"], ")"),
    "Inclusion threshold" = paste0("Inclusion\nthreshold\n(n = ", condition_counts$n[condition_counts$condition_type == "Inclusion threshold"], ")"),
    "Leave one cancer out" = paste0("Leave one cancer out\n(n = ", condition_counts$n[condition_counts$condition_type == "Leave one cancer out"], ")"),
    "Profile-altering transformation" = paste0("Profile-\naltering\ntransfor-\nmation\n(n = ", condition_counts$n[condition_counts$condition_type == "Profile-altering transformation"], ")")
)

cancer_levels <- assignment_raw_df %>%
    dplyr::distinct(cancer_type, reference_group) %>%
    dplyr::mutate(reference_group = factor(reference_group, levels = reference_group_levels)) %>%
    dplyr::arrange(reference_group, cancer_type) %>%
    dplyr::pull(cancer_type)

column_index_df <- assignment_raw_df %>%
    dplyr::distinct(condition_type, condition_id, column_label) %>%
    dplyr::mutate(
        condition_type = factor(condition_type, levels = assignment_condition_levels),
        sort_key = suppressWarnings(as.numeric(gsub("\\D", "", condition_id)))
    ) %>%
    dplyr::arrange(condition_type, sort_key, condition_id)

column_levels <- c("Final group", unique(column_index_df$column_label))

# 补全癌种×条件网格：某癌种未纳入该条件（留一法、纳入阈值）时标记为Not included ####
assignment_heatmap_df <- tidyr::expand_grid(
    column_index_df %>% dplyr::select(condition_type, condition_id, column_label),
    cancer_type = cancer_levels
) %>%
    dplyr::left_join(
        assignment_raw_df %>% dplyr::select(condition_id, cancer_type, mapped_group, same_as_final),
        by = c("condition_id", "cancer_type")
    ) %>%
    dplyr::mutate(
        tile_fill = ifelse(is.na(mapped_group), "Not included", mapped_group),
        tile_mark = dplyr::case_when(
            is.na(mapped_group) ~ "",
            !same_as_final ~ "×",
            TRUE ~ ""
        ),
        mark_colour = dplyr::case_when(
            is.na(mapped_group) ~ "grey55",
            mapped_group %in% c("Hematologic", "Visceral-like") ~ "white",
            TRUE ~ "black"
        )
    )

# 左侧颜色条：每个癌种的最终分组 ####
final_group_strip_df <- assignment_raw_df %>%
    dplyr::distinct(cancer_type, reference_group) %>%
    dplyr::mutate(
        condition_type = factor("Final group", levels = assignment_condition_levels),
        column_label = "Final group",
        tile_fill = reference_group,
        tile_mark = "",
        mark_colour = "black"
    )

assignment_plot_df <- dplyr::bind_rows(
    final_group_strip_df,
    assignment_heatmap_df %>% dplyr::select(condition_type, column_label, cancer_type, tile_fill, tile_mark, mark_colour)
) %>%
    dplyr::mutate(
        condition_type = factor(as.character(condition_type), levels = assignment_condition_levels),
        column_label = factor(column_label, levels = column_levels),
        cancer_type = factor(cancer_type, levels = rev(cancer_levels)),
        tile_fill = factor(tile_fill, levels = c(reference_group_levels, "Mixed", "Not included"))
    )

## 2. 绘图 ####
p <- ggplot(assignment_plot_df, aes(x = column_label, y = cancer_type, fill = tile_fill)) +
    geom_tile(colour = "black", linewidth = 0.15) +
    geom_text(aes(label = tile_mark, colour = mark_colour), size = 6 / .pt) +
    scale_colour_identity() +
    ggh4x::facet_grid2(
        ~condition_type,
        scales = "free_x", space = "free_x",
        labeller = ggplot2::as_labeller(facet_labels),
        strip = ggh4x::strip_themed(
            background_x = list(
                element_blank(),
                element_rect(fill = "grey95", colour = "grey60", linewidth = 0.3),
                element_rect(fill = "grey95", colour = "grey60", linewidth = 0.3),
                element_rect(fill = "grey95", colour = "grey60", linewidth = 0.3),
                element_rect(fill = "grey95", colour = "grey60", linewidth = 0.3),
                element_rect(fill = "grey80", colour = "grey60", linewidth = 0.3)
            ),
            text_x = list(
                element_text(size = 5.5, colour = "black"),
                element_text(size = 5.5, colour = "black", lineheight = 0.9),
                element_text(size = 5, colour = "black", lineheight = 0.9),
                element_text(size = 5.5, colour = "black", lineheight = 0.9),
                element_text(size = 5.5, colour = "black", lineheight = 0.9),
                element_text(size = 5, colour = "black", lineheight = 0.9)
            )
        )
    ) +
    scale_fill_manual(
        values = c(
            "Ocular" = "#DF7A5E", "Hematologic" = "#A15F88", "Neurologic" = "#82B29A",
            "Systemic-like" = "#F2CC8E", "Visceral-like" = "#BC4749",
            "Mixed" = "grey80", "Not included" = "white"
        ),
        name = NULL,
        drop = FALSE
    ) +
    labs(
        x = NULL,
        y = NULL
    ) +
    theme_minimal(base_size = 7) +
    theme(
        text = element_text(size = 7, colour = "black"),
        axis.text.x = element_text(size = 5, colour = "black", angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 6, colour = "black"),
        strip.text = element_text(size = 5.5, colour = "black", lineheight = 0.9),
        legend.text = element_text(size = 6, colour = "black"),
        legend.key.size = unit(3, "mm"),
        legend.key = element_rect(colour = "black", linewidth = 0.2),
        legend.position = "bottom",
        legend.margin = margin(0, 0, 0, 0),
        panel.grid = element_blank(),
        panel.spacing.x = unit(1.2, "mm")
    ) +
    guides(fill = guide_legend(nrow = 1))

print(p)

## 3. 保存图形 ####
ggsave(
    filename = file.path(script_results_path, "图2_各稳健性条件下癌种分组.pdf"),
    plot = p,
    width = 14.65,
    height = 11.5,
    units = "cm",
    device = cairo_pdf
)


message("已生成图形文件，目录: ", script_results_path)
message("========== NMF五分组稳健性绘图结束 ==========")
