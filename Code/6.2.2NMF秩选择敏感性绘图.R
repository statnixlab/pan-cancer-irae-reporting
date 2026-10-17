# NMF秩选择敏感性绘图（返修：审稿人1第1条、审稿人2秩选择） ####
# 读取6.2.1保存的Excel结果绘图，不依赖分析阶段的内存对象

rm(list = ls())
gc()

# 1. 加载配置与绘图R包 ####
source("Code/config.R")

suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))

# 2. 结果目录与输入文件 ####
results_dir_name <- "6.2.2NMF秩选择敏感性绘图"
script_results_path <- file.path(results_path, results_dir_name)

if (!dir.exists(script_results_path)) {
    dir.create(script_results_path, recursive = TRUE, showWarnings = FALSE)
}

input_xlsx_file <- file.path(results_path, "6.2.1NMF秩选择敏感性分析", "NMF秩选择敏感性分析结果.xlsx")

if (!file.exists(input_xlsx_file)) {
    stop("6.2.1结果文件不存在，请先运行 Code/6.2.1NMF秩选择敏感性分析.R: ", input_xlsx_file)
}

message("========== NMF秩选择敏感性绘图开始 ==========")


# 1图：Rank selection metrics, observed vs permutation null ####

## 1. 读取并处理绘图数据 ####
metrics_plot_df <- readxl::read_xlsx(input_xlsx_file, sheet = "秩选择指标_实际与零分布") %>%
    dplyr::filter(metric %in% c("cophenetic", "dispersion", "silhouette.consensus", "rss")) %>%
    dplyr::mutate(
        metric = factor(
            metric,
            levels = c("cophenetic", "dispersion", "silhouette.consensus", "rss"),
            labels = c("Cophenetic correlation", "Dispersion", "Consensus silhouette", "Residual sum of squares")
        )
    )

## 2. 绘图 ####
p <- ggplot(metrics_plot_df, aes(x = rank)) +
    annotate("rect", xmin = 4.5, xmax = 5.5, ymin = -Inf, ymax = Inf, fill = "grey92") +
    geom_ribbon(aes(ymin = null_min, ymax = null_max, fill = "Permutation null (range)"), alpha = 0.6) +
    geom_line(aes(y = null_mean, colour = "Permutation null (mean)"), linewidth = 0.4, linetype = "dashed") +
    geom_line(aes(y = observed, colour = "Observed"), linewidth = 0.4) +
    geom_point(aes(y = observed, colour = "Observed"), size = 1.2) +
    facet_wrap(~metric, nrow = 1, scales = "free_y") +
    scale_x_continuous(breaks = 2:8) +
    scale_colour_manual(values = c("Observed" = "#BC4749", "Permutation null (mean)" = "grey35")) +
    scale_fill_manual(values = c("Permutation null (range)" = "grey80")) +
    labs(
        x = "Factorization rank",
        y = "Value",
        colour = NULL,
        fill = NULL
    ) +
    theme_bw(base_size = 7) +
    theme(
        text = element_text(size = 7, colour = "black"),
        axis.text = element_text(size = 6, colour = "black"),
        plot.title = element_text(size = 7, colour = "black", face = "bold"),
        strip.text = element_text(size = 6, colour = "black"),
        strip.background = element_rect(fill = "grey95", colour = "grey60"),
        legend.position = "top",
        legend.text = element_text(size = 6, colour = "black"),
        legend.key.size = unit(3, "mm"),
        legend.margin = margin(0, 0, 0, 0),
        panel.grid.minor = element_blank()
    )

print(p)

## 3. 保存图形 ####
ggsave(
    filename = file.path(script_results_path, "图1_秩选择指标_实际与置换零分布对比.pdf"),
    plot = p,
    width = 14.65,
    height = 5.5,
    units = "cm",
    device = cairo_pdf
)


# 2图：Consensus matrix at rank 4 ####

## 1. 读取并处理绘图数据 ####
consensus_rank4_df <- readxl::read_xlsx(input_xlsx_file, sheet = "共识矩阵_秩4")
consensus_rank4_matrix <- consensus_rank4_df %>%
    tibble::column_to_rownames("cancer_type") %>%
    as.matrix()

assignment_rank4_df <- readxl::read_xlsx(input_xlsx_file, sheet = "各秩癌种分组") %>%
    dplyr::filter(rank == 4)
assignment_rank4_df <- assignment_rank4_df[match(rownames(consensus_rank4_matrix), assignment_rank4_df$cancer_type), ]

consensus_rank4_hc <- stats::hclust(stats::as.dist(1 - consensus_rank4_matrix), method = "average")

rank4_cluster_levels <- sort(unique(assignment_rank4_df$assigned_pattern))
rank4_cluster_colors <- stats::setNames(
    c("#4E79A7", "#F28E2B", "#59A14F", "#E15759", "#76B7B2", "#EDC948", "#B07AA1")[seq_along(rank4_cluster_levels)],
    rank4_cluster_levels
)
rank4_cluster_labels <- assignment_rank4_df %>%
    dplyr::distinct(assigned_pattern, assigned_pattern_label) %>%
    dplyr::arrange(assigned_pattern)

rank4_top_annotation <- ComplexHeatmap::HeatmapAnnotation(
    "Final group" = factor(assignment_rank4_df$reference_group, levels = c("Ocular", "Hematologic", "Neurologic", "Systemic-like", "Visceral-like")),
    "Rank-4 cluster" = assignment_rank4_df$assigned_pattern,
    "Silhouette" = assignment_rank4_df$silhouette,
    col = list(
        "Final group" = c("Ocular" = "#DF7A5E", "Hematologic" = "#A15F88", "Neurologic" = "#82B29A", "Systemic-like" = "#F2CC8E", "Visceral-like" = "#BC4749"),
        "Rank-4 cluster" = rank4_cluster_colors,
        "Silhouette" = circlize::colorRamp2(c(-0.2, 1), c("#F1F1F1", "#009E73"))
    ),
    annotation_name_side = "right",
    annotation_name_gp = grid::gpar(fontsize = 6, fontface = "bold", col = "black"),
    simple_anno_size = grid::unit(1.6, "mm"),
    show_legend = FALSE
)

rank4_cluster_legend <- ComplexHeatmap::Legend(
    title = "Rank-4 cluster",
    at = rank4_cluster_labels$assigned_pattern,
    labels = rank4_cluster_labels$assigned_pattern_label,
    legend_gp = grid::gpar(fill = unname(rank4_cluster_colors[rank4_cluster_labels$assigned_pattern])),
    ncol = 2,
    by_row = TRUE,
    gap = grid::unit(2, "mm"),
    labels_gp = grid::gpar(fontsize = 5, col = "black"),
    title_gp = grid::gpar(fontsize = 6, fontface = "bold", col = "black"),
    grid_height = grid::unit(2, "mm"),
    grid_width = grid::unit(2, "mm")
)

consensus_rank4_legend <- ComplexHeatmap::Legend(
    title = "Consensus",
    col_fun = circlize::colorRamp2(c(0, 0.2, 0.4, 0.6, 0.8, 1), c("#313695", "#74ADD1", "#E0F3F8", "#FEE090", "#F46D43", "#A50026")),
    at = c(0, 0.5, 1),
    direction = "horizontal",
    title_position = "topleft",
    legend_width = grid::unit(1.5, "cm"),
    grid_height = grid::unit(2, "mm"),
    labels_gp = grid::gpar(fontsize = 5, col = "black"),
    title_gp = grid::gpar(fontsize = 6, fontface = "bold", col = "black")
)

## 2. 绘图 ####
grDevices::cairo_pdf(
    filename = file.path(script_results_path, "图2_秩4共识矩阵.pdf"),
    width = 7.3 / 2.54,
    height = 8.6 / 2.54,
    pointsize = 7
)

p <- ComplexHeatmap::Heatmap(
    consensus_rank4_matrix,
    name = "Consensus",
    col = circlize::colorRamp2(c(0, 0.2, 0.4, 0.6, 0.8, 1), c("#313695", "#74ADD1", "#E0F3F8", "#FEE090", "#F46D43", "#A50026")),
    top_annotation = rank4_top_annotation,
    cluster_rows = stats::as.dendrogram(consensus_rank4_hc),
    cluster_columns = stats::as.dendrogram(consensus_rank4_hc),
    row_dend_width = grid::unit(0.3, "cm"),
    column_dend_height = grid::unit(0.3, "cm"),
    row_names_gp = grid::gpar(fontsize = 5, col = "black"),
    column_names_gp = grid::gpar(fontsize = 5, col = "black"),
    width = grid::unit(3.9, "cm"),
    height = grid::unit(3.9, "cm"),
    rect_gp = grid::gpar(col = "#4B4B7A", lwd = 0.2),
    show_heatmap_legend = FALSE
)

ComplexHeatmap::draw(
    p,
    annotation_legend_list = list(
        ComplexHeatmap::packLegend(
            consensus_rank4_legend,
            rank4_cluster_legend,
            direction = "horizontal",
            gap = grid::unit(3, "mm")
        )
    ),
    annotation_legend_side = "bottom"
)

dev.off()


# 3图：Consensus matrix at rank 5 ####

## 1. 读取并处理绘图数据 ####
consensus_rank5_df <- readxl::read_xlsx(input_xlsx_file, sheet = "共识矩阵_秩5")
consensus_rank5_matrix <- consensus_rank5_df %>%
    tibble::column_to_rownames("cancer_type") %>%
    as.matrix()

assignment_rank5_df <- readxl::read_xlsx(input_xlsx_file, sheet = "各秩癌种分组") %>%
    dplyr::filter(rank == 5)
assignment_rank5_df <- assignment_rank5_df[match(rownames(consensus_rank5_matrix), assignment_rank5_df$cancer_type), ]

consensus_rank5_hc <- stats::hclust(stats::as.dist(1 - consensus_rank5_matrix), method = "average")

rank5_cluster_levels <- sort(unique(assignment_rank5_df$assigned_pattern))
rank5_cluster_colors <- stats::setNames(
    c("#4E79A7", "#F28E2B", "#59A14F", "#E15759", "#76B7B2", "#EDC948", "#B07AA1")[seq_along(rank5_cluster_levels)],
    rank5_cluster_levels
)
rank5_cluster_labels <- assignment_rank5_df %>%
    dplyr::distinct(assigned_pattern, assigned_pattern_label) %>%
    dplyr::arrange(assigned_pattern)

rank5_top_annotation <- ComplexHeatmap::HeatmapAnnotation(
    "Final group" = factor(assignment_rank5_df$reference_group, levels = c("Ocular", "Hematologic", "Neurologic", "Systemic-like", "Visceral-like")),
    "Rank-5 cluster" = assignment_rank5_df$assigned_pattern,
    "Silhouette" = assignment_rank5_df$silhouette,
    col = list(
        "Final group" = c("Ocular" = "#DF7A5E", "Hematologic" = "#A15F88", "Neurologic" = "#82B29A", "Systemic-like" = "#F2CC8E", "Visceral-like" = "#BC4749"),
        "Rank-5 cluster" = rank5_cluster_colors,
        "Silhouette" = circlize::colorRamp2(c(-0.2, 1), c("#F1F1F1", "#009E73"))
    ),
    annotation_name_side = "right",
    annotation_name_gp = grid::gpar(fontsize = 6, fontface = "bold", col = "black"),
    simple_anno_size = grid::unit(1.6, "mm"),
    show_legend = FALSE
)

rank5_cluster_legend <- ComplexHeatmap::Legend(
    title = "Rank-5 cluster",
    at = rank5_cluster_labels$assigned_pattern,
    labels = rank5_cluster_labels$assigned_pattern_label,
    legend_gp = grid::gpar(fill = unname(rank5_cluster_colors[rank5_cluster_labels$assigned_pattern])),
    ncol = 2,
    by_row = TRUE,
    gap = grid::unit(2, "mm"),
    labels_gp = grid::gpar(fontsize = 5, col = "black"),
    title_gp = grid::gpar(fontsize = 6, fontface = "bold", col = "black"),
    grid_height = grid::unit(2, "mm"),
    grid_width = grid::unit(2, "mm")
)

consensus_rank5_legend <- ComplexHeatmap::Legend(
    title = "Consensus",
    col_fun = circlize::colorRamp2(c(0, 0.2, 0.4, 0.6, 0.8, 1), c("#313695", "#74ADD1", "#E0F3F8", "#FEE090", "#F46D43", "#A50026")),
    at = c(0, 0.5, 1),
    direction = "horizontal",
    title_position = "topleft",
    legend_width = grid::unit(1.5, "cm"),
    grid_height = grid::unit(2, "mm"),
    labels_gp = grid::gpar(fontsize = 5, col = "black"),
    title_gp = grid::gpar(fontsize = 6, fontface = "bold", col = "black")
)

## 2. 绘图 ####
grDevices::cairo_pdf(
    filename = file.path(script_results_path, "图3_秩5共识矩阵.pdf"),
    width = 7.3 / 2.54,
    height = 8.6 / 2.54,
    pointsize = 7
)

p <- ComplexHeatmap::Heatmap(
    consensus_rank5_matrix,
    name = "Consensus",
    col = circlize::colorRamp2(c(0, 0.2, 0.4, 0.6, 0.8, 1), c("#313695", "#74ADD1", "#E0F3F8", "#FEE090", "#F46D43", "#A50026")),
    top_annotation = rank5_top_annotation,
    cluster_rows = stats::as.dendrogram(consensus_rank5_hc),
    cluster_columns = stats::as.dendrogram(consensus_rank5_hc),
    row_dend_width = grid::unit(0.3, "cm"),
    column_dend_height = grid::unit(0.3, "cm"),
    row_names_gp = grid::gpar(fontsize = 5, col = "black"),
    column_names_gp = grid::gpar(fontsize = 5, col = "black"),
    width = grid::unit(3.9, "cm"),
    height = grid::unit(3.9, "cm"),
    rect_gp = grid::gpar(col = "#4B4B7A", lwd = 0.2),
    show_heatmap_legend = FALSE
)

ComplexHeatmap::draw(
    p,
    annotation_legend_list = list(
        ComplexHeatmap::packLegend(
            consensus_rank5_legend,
            rank5_cluster_legend,
            direction = "horizontal",
            gap = grid::unit(3, "mm")
        )
    ),
    annotation_legend_side = "bottom"
)

dev.off()


# 4图：Consensus matrix at rank 6 ####

## 1. 读取并处理绘图数据 ####
consensus_rank6_df <- readxl::read_xlsx(input_xlsx_file, sheet = "共识矩阵_秩6")
consensus_rank6_matrix <- consensus_rank6_df %>%
    tibble::column_to_rownames("cancer_type") %>%
    as.matrix()

assignment_rank6_df <- readxl::read_xlsx(input_xlsx_file, sheet = "各秩癌种分组") %>%
    dplyr::filter(rank == 6)
assignment_rank6_df <- assignment_rank6_df[match(rownames(consensus_rank6_matrix), assignment_rank6_df$cancer_type), ]

consensus_rank6_hc <- stats::hclust(stats::as.dist(1 - consensus_rank6_matrix), method = "average")

rank6_cluster_levels <- sort(unique(assignment_rank6_df$assigned_pattern))
rank6_cluster_colors <- stats::setNames(
    c("#4E79A7", "#F28E2B", "#59A14F", "#E15759", "#76B7B2", "#EDC948", "#B07AA1")[seq_along(rank6_cluster_levels)],
    rank6_cluster_levels
)
rank6_cluster_labels <- assignment_rank6_df %>%
    dplyr::distinct(assigned_pattern, assigned_pattern_label) %>%
    dplyr::arrange(assigned_pattern)

rank6_top_annotation <- ComplexHeatmap::HeatmapAnnotation(
    "Final group" = factor(assignment_rank6_df$reference_group, levels = c("Ocular", "Hematologic", "Neurologic", "Systemic-like", "Visceral-like")),
    "Rank-6 cluster" = assignment_rank6_df$assigned_pattern,
    "Silhouette" = assignment_rank6_df$silhouette,
    col = list(
        "Final group" = c("Ocular" = "#DF7A5E", "Hematologic" = "#A15F88", "Neurologic" = "#82B29A", "Systemic-like" = "#F2CC8E", "Visceral-like" = "#BC4749"),
        "Rank-6 cluster" = rank6_cluster_colors,
        "Silhouette" = circlize::colorRamp2(c(-0.2, 1), c("#F1F1F1", "#009E73"))
    ),
    annotation_name_side = "right",
    annotation_name_gp = grid::gpar(fontsize = 6, fontface = "bold", col = "black"),
    simple_anno_size = grid::unit(1.6, "mm"),
    show_legend = FALSE
)

rank6_cluster_legend <- ComplexHeatmap::Legend(
    title = "Rank-6 cluster",
    at = rank6_cluster_labels$assigned_pattern,
    labels = rank6_cluster_labels$assigned_pattern_label,
    legend_gp = grid::gpar(fill = unname(rank6_cluster_colors[rank6_cluster_labels$assigned_pattern])),
    ncol = 2,
    by_row = TRUE,
    gap = grid::unit(2, "mm"),
    labels_gp = grid::gpar(fontsize = 5, col = "black"),
    title_gp = grid::gpar(fontsize = 6, fontface = "bold", col = "black"),
    grid_height = grid::unit(2, "mm"),
    grid_width = grid::unit(2, "mm")
)

consensus_rank6_legend <- ComplexHeatmap::Legend(
    title = "Consensus",
    col_fun = circlize::colorRamp2(c(0, 0.2, 0.4, 0.6, 0.8, 1), c("#313695", "#74ADD1", "#E0F3F8", "#FEE090", "#F46D43", "#A50026")),
    at = c(0, 0.5, 1),
    direction = "horizontal",
    title_position = "topleft",
    legend_width = grid::unit(1.5, "cm"),
    grid_height = grid::unit(2, "mm"),
    labels_gp = grid::gpar(fontsize = 5, col = "black"),
    title_gp = grid::gpar(fontsize = 6, fontface = "bold", col = "black")
)

## 2. 绘图 ####
grDevices::cairo_pdf(
    filename = file.path(script_results_path, "图4_秩6共识矩阵.pdf"),
    width = 7.3 / 2.54,
    height = 8.6 / 2.54,
    pointsize = 7
)

p <- ComplexHeatmap::Heatmap(
    consensus_rank6_matrix,
    name = "Consensus",
    col = circlize::colorRamp2(c(0, 0.2, 0.4, 0.6, 0.8, 1), c("#313695", "#74ADD1", "#E0F3F8", "#FEE090", "#F46D43", "#A50026")),
    top_annotation = rank6_top_annotation,
    cluster_rows = stats::as.dendrogram(consensus_rank6_hc),
    cluster_columns = stats::as.dendrogram(consensus_rank6_hc),
    row_dend_width = grid::unit(0.3, "cm"),
    column_dend_height = grid::unit(0.3, "cm"),
    row_names_gp = grid::gpar(fontsize = 5, col = "black"),
    column_names_gp = grid::gpar(fontsize = 5, col = "black"),
    width = grid::unit(3.9, "cm"),
    height = grid::unit(3.9, "cm"),
    rect_gp = grid::gpar(col = "#4B4B7A", lwd = 0.2),
    show_heatmap_legend = FALSE
)

ComplexHeatmap::draw(
    p,
    annotation_legend_list = list(
        ComplexHeatmap::packLegend(
            consensus_rank6_legend,
            rank6_cluster_legend,
            direction = "horizontal",
            gap = grid::unit(3, "mm")
        )
    ),
    annotation_legend_side = "bottom"
)

dev.off()


# 5图：Consensus matrix at rank 7 ####

## 1. 读取并处理绘图数据 ####
consensus_rank7_df <- readxl::read_xlsx(input_xlsx_file, sheet = "共识矩阵_秩7")
consensus_rank7_matrix <- consensus_rank7_df %>%
    tibble::column_to_rownames("cancer_type") %>%
    as.matrix()

assignment_rank7_df <- readxl::read_xlsx(input_xlsx_file, sheet = "各秩癌种分组") %>%
    dplyr::filter(rank == 7)
assignment_rank7_df <- assignment_rank7_df[match(rownames(consensus_rank7_matrix), assignment_rank7_df$cancer_type), ]

consensus_rank7_hc <- stats::hclust(stats::as.dist(1 - consensus_rank7_matrix), method = "average")

rank7_cluster_levels <- sort(unique(assignment_rank7_df$assigned_pattern))
rank7_cluster_colors <- stats::setNames(
    c("#4E79A7", "#F28E2B", "#59A14F", "#E15759", "#76B7B2", "#EDC948", "#B07AA1")[seq_along(rank7_cluster_levels)],
    rank7_cluster_levels
)
rank7_cluster_labels <- assignment_rank7_df %>%
    dplyr::distinct(assigned_pattern, assigned_pattern_label) %>%
    dplyr::arrange(assigned_pattern)

rank7_top_annotation <- ComplexHeatmap::HeatmapAnnotation(
    "Final group" = factor(assignment_rank7_df$reference_group, levels = c("Ocular", "Hematologic", "Neurologic", "Systemic-like", "Visceral-like")),
    "Rank-7 cluster" = assignment_rank7_df$assigned_pattern,
    "Silhouette" = assignment_rank7_df$silhouette,
    col = list(
        "Final group" = c("Ocular" = "#DF7A5E", "Hematologic" = "#A15F88", "Neurologic" = "#82B29A", "Systemic-like" = "#F2CC8E", "Visceral-like" = "#BC4749"),
        "Rank-7 cluster" = rank7_cluster_colors,
        "Silhouette" = circlize::colorRamp2(c(-0.2, 1), c("#F1F1F1", "#009E73"))
    ),
    annotation_name_side = "right",
    annotation_name_gp = grid::gpar(fontsize = 6, fontface = "bold", col = "black"),
    simple_anno_size = grid::unit(1.6, "mm"),
    show_legend = FALSE
)

rank7_cluster_legend <- ComplexHeatmap::Legend(
    title = "Rank-7 cluster",
    at = rank7_cluster_labels$assigned_pattern,
    labels = rank7_cluster_labels$assigned_pattern_label,
    legend_gp = grid::gpar(fill = unname(rank7_cluster_colors[rank7_cluster_labels$assigned_pattern])),
    ncol = 3,
    by_row = TRUE,
    gap = grid::unit(2, "mm"),
    labels_gp = grid::gpar(fontsize = 5, col = "black"),
    title_gp = grid::gpar(fontsize = 6, fontface = "bold", col = "black"),
    grid_height = grid::unit(2, "mm"),
    grid_width = grid::unit(2, "mm")
)

consensus_rank7_legend <- ComplexHeatmap::Legend(
    title = "Consensus",
    col_fun = circlize::colorRamp2(c(0, 0.2, 0.4, 0.6, 0.8, 1), c("#313695", "#74ADD1", "#E0F3F8", "#FEE090", "#F46D43", "#A50026")),
    at = c(0, 0.5, 1),
    direction = "horizontal",
    title_position = "topleft",
    legend_width = grid::unit(1.5, "cm"),
    grid_height = grid::unit(2, "mm"),
    labels_gp = grid::gpar(fontsize = 5, col = "black"),
    title_gp = grid::gpar(fontsize = 6, fontface = "bold", col = "black")
)

## 2. 绘图 ####
grDevices::cairo_pdf(
    filename = file.path(script_results_path, "图5_秩7共识矩阵.pdf"),
    width = 7.3 / 2.54,
    height = 8.6 / 2.54,
    pointsize = 7
)

p <- ComplexHeatmap::Heatmap(
    consensus_rank7_matrix,
    name = "Consensus",
    col = circlize::colorRamp2(c(0, 0.2, 0.4, 0.6, 0.8, 1), c("#313695", "#74ADD1", "#E0F3F8", "#FEE090", "#F46D43", "#A50026")),
    top_annotation = rank7_top_annotation,
    cluster_rows = stats::as.dendrogram(consensus_rank7_hc),
    cluster_columns = stats::as.dendrogram(consensus_rank7_hc),
    row_dend_width = grid::unit(0.3, "cm"),
    column_dend_height = grid::unit(0.3, "cm"),
    row_names_gp = grid::gpar(fontsize = 5, col = "black"),
    column_names_gp = grid::gpar(fontsize = 5, col = "black"),
    width = grid::unit(3.9, "cm"),
    height = grid::unit(3.9, "cm"),
    rect_gp = grid::gpar(col = "#4B4B7A", lwd = 0.2),
    show_heatmap_legend = FALSE
)

ComplexHeatmap::draw(
    p,
    annotation_legend_list = list(
        ComplexHeatmap::packLegend(
            consensus_rank7_legend,
            rank7_cluster_legend,
            direction = "horizontal",
            gap = grid::unit(3, "mm")
        )
    ),
    annotation_legend_side = "bottom"
)

dev.off()


# 11图：Consensus matrix at rank 2 ####

## 1. 读取并处理绘图数据 ####
consensus_rank2_df <- readxl::read_xlsx(input_xlsx_file, sheet = "共识矩阵_秩2")
consensus_rank2_matrix <- consensus_rank2_df %>%
    tibble::column_to_rownames("cancer_type") %>%
    as.matrix()

assignment_rank2_df <- readxl::read_xlsx(input_xlsx_file, sheet = "各秩癌种分组") %>%
    dplyr::filter(rank == 2)
assignment_rank2_df <- assignment_rank2_df[match(rownames(consensus_rank2_matrix), assignment_rank2_df$cancer_type), ]

consensus_rank2_hc <- stats::hclust(stats::as.dist(1 - consensus_rank2_matrix), method = "average")

rank2_cluster_levels <- sort(unique(assignment_rank2_df$assigned_pattern))
rank2_cluster_colors <- stats::setNames(
    c("#4E79A7", "#F28E2B", "#59A14F", "#E15759", "#76B7B2", "#EDC948", "#B07AA1")[seq_along(rank2_cluster_levels)],
    rank2_cluster_levels
)
rank2_cluster_labels <- assignment_rank2_df %>%
    dplyr::distinct(assigned_pattern, assigned_pattern_label) %>%
    dplyr::arrange(assigned_pattern)

rank2_top_annotation <- ComplexHeatmap::HeatmapAnnotation(
    "Final group" = factor(assignment_rank2_df$reference_group, levels = c("Ocular", "Hematologic", "Neurologic", "Systemic-like", "Visceral-like")),
    "Rank-2 cluster" = assignment_rank2_df$assigned_pattern,
    "Silhouette" = assignment_rank2_df$silhouette,
    col = list(
        "Final group" = c("Ocular" = "#DF7A5E", "Hematologic" = "#A15F88", "Neurologic" = "#82B29A", "Systemic-like" = "#F2CC8E", "Visceral-like" = "#BC4749"),
        "Rank-2 cluster" = rank2_cluster_colors,
        "Silhouette" = circlize::colorRamp2(c(-0.2, 1), c("#F1F1F1", "#009E73"))
    ),
    annotation_name_side = "right",
    annotation_name_gp = grid::gpar(fontsize = 6, fontface = "bold", col = "black"),
    simple_anno_size = grid::unit(1.6, "mm"),
    show_legend = FALSE
)

rank2_cluster_legend <- ComplexHeatmap::Legend(
    title = "Rank-2 cluster",
    at = rank2_cluster_labels$assigned_pattern,
    labels = rank2_cluster_labels$assigned_pattern_label,
    legend_gp = grid::gpar(fill = unname(rank2_cluster_colors[rank2_cluster_labels$assigned_pattern])),
    ncol = 1,
    by_row = TRUE,
    gap = grid::unit(2, "mm"),
    labels_gp = grid::gpar(fontsize = 5, col = "black"),
    title_gp = grid::gpar(fontsize = 6, fontface = "bold", col = "black"),
    grid_height = grid::unit(2, "mm"),
    grid_width = grid::unit(2, "mm")
)

consensus_rank2_legend <- ComplexHeatmap::Legend(
    title = "Consensus",
    col_fun = circlize::colorRamp2(c(0, 0.2, 0.4, 0.6, 0.8, 1), c("#313695", "#74ADD1", "#E0F3F8", "#FEE090", "#F46D43", "#A50026")),
    at = c(0, 0.5, 1),
    direction = "horizontal",
    title_position = "topleft",
    legend_width = grid::unit(1.5, "cm"),
    grid_height = grid::unit(2, "mm"),
    labels_gp = grid::gpar(fontsize = 5, col = "black"),
    title_gp = grid::gpar(fontsize = 6, fontface = "bold", col = "black")
)

## 2. 绘图 ####
grDevices::cairo_pdf(
    filename = file.path(script_results_path, "图11_秩2共识矩阵.pdf"),
    width = 7.3 / 2.54,
    height = 8.6 / 2.54,
    pointsize = 7
)

p <- ComplexHeatmap::Heatmap(
    consensus_rank2_matrix,
    name = "Consensus",
    col = circlize::colorRamp2(c(0, 0.2, 0.4, 0.6, 0.8, 1), c("#313695", "#74ADD1", "#E0F3F8", "#FEE090", "#F46D43", "#A50026")),
    top_annotation = rank2_top_annotation,
    cluster_rows = stats::as.dendrogram(consensus_rank2_hc),
    cluster_columns = stats::as.dendrogram(consensus_rank2_hc),
    row_dend_width = grid::unit(0.3, "cm"),
    column_dend_height = grid::unit(0.3, "cm"),
    row_names_gp = grid::gpar(fontsize = 5, col = "black"),
    column_names_gp = grid::gpar(fontsize = 5, col = "black"),
    width = grid::unit(3.9, "cm"),
    height = grid::unit(3.9, "cm"),
    rect_gp = grid::gpar(col = "#4B4B7A", lwd = 0.2),
    show_heatmap_legend = FALSE
)

ComplexHeatmap::draw(
    p,
    annotation_legend_list = list(
        ComplexHeatmap::packLegend(
            consensus_rank2_legend,
            rank2_cluster_legend,
            direction = "horizontal",
            gap = grid::unit(3, "mm")
        )
    ),
    annotation_legend_side = "bottom"
)

dev.off()


# 12图：Consensus matrix at rank 3 ####

## 1. 读取并处理绘图数据 ####
consensus_rank3_df <- readxl::read_xlsx(input_xlsx_file, sheet = "共识矩阵_秩3")
consensus_rank3_matrix <- consensus_rank3_df %>%
    tibble::column_to_rownames("cancer_type") %>%
    as.matrix()

assignment_rank3_df <- readxl::read_xlsx(input_xlsx_file, sheet = "各秩癌种分组") %>%
    dplyr::filter(rank == 3)
assignment_rank3_df <- assignment_rank3_df[match(rownames(consensus_rank3_matrix), assignment_rank3_df$cancer_type), ]

consensus_rank3_hc <- stats::hclust(stats::as.dist(1 - consensus_rank3_matrix), method = "average")

rank3_cluster_levels <- sort(unique(assignment_rank3_df$assigned_pattern))
rank3_cluster_colors <- stats::setNames(
    c("#4E79A7", "#F28E2B", "#59A14F", "#E15759", "#76B7B2", "#EDC948", "#B07AA1")[seq_along(rank3_cluster_levels)],
    rank3_cluster_levels
)
rank3_cluster_labels <- assignment_rank3_df %>%
    dplyr::distinct(assigned_pattern, assigned_pattern_label) %>%
    dplyr::arrange(assigned_pattern)

rank3_top_annotation <- ComplexHeatmap::HeatmapAnnotation(
    "Final group" = factor(assignment_rank3_df$reference_group, levels = c("Ocular", "Hematologic", "Neurologic", "Systemic-like", "Visceral-like")),
    "Rank-3 cluster" = assignment_rank3_df$assigned_pattern,
    "Silhouette" = assignment_rank3_df$silhouette,
    col = list(
        "Final group" = c("Ocular" = "#DF7A5E", "Hematologic" = "#A15F88", "Neurologic" = "#82B29A", "Systemic-like" = "#F2CC8E", "Visceral-like" = "#BC4749"),
        "Rank-3 cluster" = rank3_cluster_colors,
        "Silhouette" = circlize::colorRamp2(c(-0.2, 1), c("#F1F1F1", "#009E73"))
    ),
    annotation_name_side = "right",
    annotation_name_gp = grid::gpar(fontsize = 6, fontface = "bold", col = "black"),
    simple_anno_size = grid::unit(1.6, "mm"),
    show_legend = FALSE
)

rank3_cluster_legend <- ComplexHeatmap::Legend(
    title = "Rank-3 cluster",
    at = rank3_cluster_labels$assigned_pattern,
    labels = rank3_cluster_labels$assigned_pattern_label,
    legend_gp = grid::gpar(fill = unname(rank3_cluster_colors[rank3_cluster_labels$assigned_pattern])),
    ncol = 1,
    by_row = TRUE,
    gap = grid::unit(2, "mm"),
    labels_gp = grid::gpar(fontsize = 5, col = "black"),
    title_gp = grid::gpar(fontsize = 6, fontface = "bold", col = "black"),
    grid_height = grid::unit(2, "mm"),
    grid_width = grid::unit(2, "mm")
)

consensus_rank3_legend <- ComplexHeatmap::Legend(
    title = "Consensus",
    col_fun = circlize::colorRamp2(c(0, 0.2, 0.4, 0.6, 0.8, 1), c("#313695", "#74ADD1", "#E0F3F8", "#FEE090", "#F46D43", "#A50026")),
    at = c(0, 0.5, 1),
    direction = "horizontal",
    title_position = "topleft",
    legend_width = grid::unit(1.5, "cm"),
    grid_height = grid::unit(2, "mm"),
    labels_gp = grid::gpar(fontsize = 5, col = "black"),
    title_gp = grid::gpar(fontsize = 6, fontface = "bold", col = "black")
)

## 2. 绘图 ####
grDevices::cairo_pdf(
    filename = file.path(script_results_path, "图12_秩3共识矩阵.pdf"),
    width = 7.3 / 2.54,
    height = 8.6 / 2.54,
    pointsize = 7
)

p <- ComplexHeatmap::Heatmap(
    consensus_rank3_matrix,
    name = "Consensus",
    col = circlize::colorRamp2(c(0, 0.2, 0.4, 0.6, 0.8, 1), c("#313695", "#74ADD1", "#E0F3F8", "#FEE090", "#F46D43", "#A50026")),
    top_annotation = rank3_top_annotation,
    cluster_rows = stats::as.dendrogram(consensus_rank3_hc),
    cluster_columns = stats::as.dendrogram(consensus_rank3_hc),
    row_dend_width = grid::unit(0.3, "cm"),
    column_dend_height = grid::unit(0.3, "cm"),
    row_names_gp = grid::gpar(fontsize = 5, col = "black"),
    column_names_gp = grid::gpar(fontsize = 5, col = "black"),
    width = grid::unit(3.9, "cm"),
    height = grid::unit(3.9, "cm"),
    rect_gp = grid::gpar(col = "#4B4B7A", lwd = 0.2),
    show_heatmap_legend = FALSE
)

ComplexHeatmap::draw(
    p,
    annotation_legend_list = list(
        ComplexHeatmap::packLegend(
            consensus_rank3_legend,
            rank3_cluster_legend,
            direction = "horizontal",
            gap = grid::unit(3, "mm")
        )
    ),
    annotation_legend_side = "bottom"
)

dev.off()


# 13图：Consensus matrix at rank 8 ####

## 1. 读取并处理绘图数据 ####
consensus_rank8_df <- readxl::read_xlsx(input_xlsx_file, sheet = "共识矩阵_秩8")
consensus_rank8_matrix <- consensus_rank8_df %>%
    tibble::column_to_rownames("cancer_type") %>%
    as.matrix()

assignment_rank8_df <- readxl::read_xlsx(input_xlsx_file, sheet = "各秩癌种分组") %>%
    dplyr::filter(rank == 8)
assignment_rank8_df <- assignment_rank8_df[match(rownames(consensus_rank8_matrix), assignment_rank8_df$cancer_type), ]

consensus_rank8_hc <- stats::hclust(stats::as.dist(1 - consensus_rank8_matrix), method = "average")

rank8_cluster_levels <- sort(unique(assignment_rank8_df$assigned_pattern))
rank8_cluster_colors <- stats::setNames(
    c("#4E79A7", "#F28E2B", "#59A14F", "#E15759", "#76B7B2", "#EDC948", "#B07AA1", "#9C755F")[seq_along(rank8_cluster_levels)],
    rank8_cluster_levels
)
rank8_cluster_labels <- assignment_rank8_df %>%
    dplyr::distinct(assigned_pattern, assigned_pattern_label) %>%
    dplyr::arrange(assigned_pattern)

rank8_top_annotation <- ComplexHeatmap::HeatmapAnnotation(
    "Final group" = factor(assignment_rank8_df$reference_group, levels = c("Ocular", "Hematologic", "Neurologic", "Systemic-like", "Visceral-like")),
    "Rank-8 cluster" = assignment_rank8_df$assigned_pattern,
    "Silhouette" = assignment_rank8_df$silhouette,
    col = list(
        "Final group" = c("Ocular" = "#DF7A5E", "Hematologic" = "#A15F88", "Neurologic" = "#82B29A", "Systemic-like" = "#F2CC8E", "Visceral-like" = "#BC4749"),
        "Rank-8 cluster" = rank8_cluster_colors,
        "Silhouette" = circlize::colorRamp2(c(-0.2, 1), c("#F1F1F1", "#009E73"))
    ),
    annotation_name_side = "right",
    annotation_name_gp = grid::gpar(fontsize = 6, fontface = "bold", col = "black"),
    simple_anno_size = grid::unit(1.6, "mm"),
    show_legend = FALSE
)

rank8_cluster_legend <- ComplexHeatmap::Legend(
    title = "Rank-8 cluster",
    at = rank8_cluster_labels$assigned_pattern,
    labels = rank8_cluster_labels$assigned_pattern_label,
    legend_gp = grid::gpar(fill = unname(rank8_cluster_colors[rank8_cluster_labels$assigned_pattern])),
    ncol = 3,
    by_row = TRUE,
    gap = grid::unit(2, "mm"),
    labels_gp = grid::gpar(fontsize = 5, col = "black"),
    title_gp = grid::gpar(fontsize = 6, fontface = "bold", col = "black"),
    grid_height = grid::unit(2, "mm"),
    grid_width = grid::unit(2, "mm")
)

consensus_rank8_legend <- ComplexHeatmap::Legend(
    title = "Consensus",
    col_fun = circlize::colorRamp2(c(0, 0.2, 0.4, 0.6, 0.8, 1), c("#313695", "#74ADD1", "#E0F3F8", "#FEE090", "#F46D43", "#A50026")),
    at = c(0, 0.5, 1),
    direction = "horizontal",
    title_position = "topleft",
    legend_width = grid::unit(1.5, "cm"),
    grid_height = grid::unit(2, "mm"),
    labels_gp = grid::gpar(fontsize = 5, col = "black"),
    title_gp = grid::gpar(fontsize = 6, fontface = "bold", col = "black")
)

## 2. 绘图 ####
grDevices::cairo_pdf(
    filename = file.path(script_results_path, "图13_秩8共识矩阵.pdf"),
    width = 7.3 / 2.54,
    height = 8.6 / 2.54,
    pointsize = 7
)

p <- ComplexHeatmap::Heatmap(
    consensus_rank8_matrix,
    name = "Consensus",
    col = circlize::colorRamp2(c(0, 0.2, 0.4, 0.6, 0.8, 1), c("#313695", "#74ADD1", "#E0F3F8", "#FEE090", "#F46D43", "#A50026")),
    top_annotation = rank8_top_annotation,
    cluster_rows = stats::as.dendrogram(consensus_rank8_hc),
    cluster_columns = stats::as.dendrogram(consensus_rank8_hc),
    row_dend_width = grid::unit(0.3, "cm"),
    column_dend_height = grid::unit(0.3, "cm"),
    row_names_gp = grid::gpar(fontsize = 5, col = "black"),
    column_names_gp = grid::gpar(fontsize = 5, col = "black"),
    width = grid::unit(3.9, "cm"),
    height = grid::unit(3.9, "cm"),
    rect_gp = grid::gpar(col = "#4B4B7A", lwd = 0.2),
    show_heatmap_legend = FALSE
)

ComplexHeatmap::draw(
    p,
    annotation_legend_list = list(
        ComplexHeatmap::packLegend(
            consensus_rank8_legend,
            rank8_cluster_legend,
            direction = "horizontal",
            gap = grid::unit(3, "mm")
        )
    ),
    annotation_legend_side = "bottom"
)

dev.off()


# 7图：Cancer assignments across ranks 2-8 relative to the final five groups ####

## 1. 读取并处理绘图数据 ####
assignment_plot_df <- readxl::read_xlsx(input_xlsx_file, sheet = "各秩癌种分组")
ari_plot_df <- readxl::read_xlsx(input_xlsx_file, sheet = "分组一致性ARI") %>%
    dplyr::filter(comparison == "vs final rank-5 groups (500 runs)")

reference_group_levels <- c("Ocular", "Hematologic", "Neurologic", "Systemic-like", "Visceral-like")

# 群内最终分组并列最多时（无唯一多数），记为Mixed ####
cluster_majority_df <- assignment_plot_df %>%
    dplyr::count(rank, assigned_pattern, reference_group, name = "n") %>%
    dplyr::group_by(rank, assigned_pattern) %>%
    dplyr::filter(n == max(n)) %>%
    dplyr::summarise(
        majority_group = ifelse(dplyr::n() == 1, dplyr::first(reference_group), "Mixed (no majority)"),
        .groups = "drop"
    )

assignment_plot_df <- assignment_plot_df %>%
    dplyr::left_join(cluster_majority_df, by = c("rank", "assigned_pattern")) %>%
    dplyr::left_join(ari_plot_df %>% dplyr::select(rank, ari), by = "rank") %>%
    dplyr::mutate(
        rank_label = paste0("Rank ", rank, "\nARI = ", sprintf("%.2f", ari)),
        reference_group = factor(reference_group, levels = reference_group_levels),
        majority_group = factor(majority_group, levels = c(reference_group_levels, "Mixed (no majority)")),
        differs_from_majority = majority_group != "Mixed (no majority)" & as.character(reference_group) != as.character(majority_group),
        # 每行最多放2个SOC，超出部分在"/"后换行 ####
        tile_label = purrr::map_chr(
            strsplit(assigned_pattern_label, "/", fixed = TRUE),
            ~ paste(
                vapply(split(.x, ceiling(seq_along(.x) / 2)), paste, character(1), collapse = "/"),
                collapse = "/\n"
            )
        ),
        tile_label = paste0(tile_label, ifelse(differs_from_majority, " *", "")),
        tile_text_colour = ifelse(majority_group %in% c("Hematologic", "Visceral-like"), "white", "black")
    ) %>%
    dplyr::arrange(reference_group, cancer_type)

cancer_levels <- assignment_plot_df %>%
    dplyr::filter(rank == 5) %>%
    dplyr::arrange(reference_group, cancer_type) %>%
    dplyr::pull(cancer_type)

rank_label_levels <- c("Final\ngroup", unique(assignment_plot_df$rank_label[order(assignment_plot_df$rank)]))

assignment_plot_df <- assignment_plot_df %>%
    dplyr::mutate(
        cancer_type = factor(cancer_type, levels = rev(cancer_levels)),
        rank_label = factor(rank_label, levels = rank_label_levels)
    )

# 左侧颜色条：每个癌种自身的最终rank-5分组 ####
final_group_strip_df <- assignment_plot_df %>%
    dplyr::distinct(cancer_type, reference_group) %>%
    dplyr::mutate(
        rank_label = factor("Final\ngroup", levels = rank_label_levels),
        strip_fill = factor(as.character(reference_group), levels = c(reference_group_levels, "Mixed (no majority)"))
    )

## 2. 绘图 ####
# 使用数值横坐标，将最终分组色条靠近Rank 2，并压缩癌种标签旁的留白 ####
p <- ggplot(assignment_plot_df, aes(x = rank, y = cancer_type, fill = majority_group)) +
    geom_tile(width = 1, colour = "black", linewidth = 0.3) +
    geom_tile(
        data = final_group_strip_df,
        aes(x = 1.2, y = cancer_type, fill = strip_fill),
        width = 0.5, colour = "black", linewidth = 0.3, inherit.aes = FALSE
    ) +
    geom_text(aes(label = tile_label, colour = tile_text_colour), size = 5 / .pt, lineheight = 0.75) +
    scale_colour_identity() +
    scale_x_continuous(
        breaks = c(1.2, 2:8),
        labels = rank_label_levels,
        limits = c(0.9, 8.5),
        expand = expansion(mult = 0)
    ) +
    scale_fill_manual(
        values = c("Ocular" = "#DF7A5E", "Hematologic" = "#A15F88", "Neurologic" = "#82B29A", "Systemic-like" = "#F2CC8E", "Visceral-like" = "#BC4749", "Mixed (no majority)" = "grey80"),
        name = "Final rank-5 group",
        drop = FALSE
    ) +
    labs(
        x = NULL,
        y = NULL
    ) +
    theme_minimal(base_size = 7) +
    theme(
        text = element_text(size = 7, colour = "black"),
        axis.text.x = element_text(size = 6, colour = "black"),
        axis.text.y = element_text(size = 6, colour = "black", margin = margin(r = 2)),
        plot.title = element_text(size = 7, colour = "black", face = "bold"),
        plot.caption = element_text(size = 5, colour = "black", hjust = 0),
        legend.title = element_text(size = 6, colour = "black"),
        legend.text = element_text(size = 6, colour = "black"),
        legend.key.size = unit(3, "mm"),
        legend.position = "bottom",
        legend.margin = margin(0, 0, 0, 0),
        panel.grid = element_blank()
    )

print(p)

## 3. 保存图形 ####
ggsave(
    filename = file.path(script_results_path, "图7_秩2至8癌种分组对照.pdf"),
    plot = p,
    width = 14.65,
    height = 20,
    units = "cm",
    device = cairo_pdf
)


# 8图：Excess decrease in RSS when adding each factor, relative to permutation null ####

## 1. 读取并处理绘图数据 ####
rss_excess_plot_df <- readxl::read_xlsx(input_xlsx_file, sheet = "RSS下降_实际与零分布") %>%
    dplyr::mutate(
        factor_added = factor(paste0(rank - 1, "→", rank), levels = paste0(rank - 1, "→", rank)),
        excess_decrease = rss_drop_observed - rss_drop_null_mean,
        null_range_low = rss_drop_null_min - rss_drop_null_mean,
        null_range_high = rss_drop_null_max - rss_drop_null_mean,
        excess_direction = ifelse(excess_decrease > 0, "Above permutation mean", "Below permutation mean"),
        excess_label = sprintf("%+.2f", excess_decrease),
        label_vjust = ifelse(excess_decrease > 0, -0.6, 1.6)
    )

selected_rank_position <- which(levels(rss_excess_plot_df$factor_added) == "4→5") + 0.5

## 2. 绘图 ####
p <- ggplot(rss_excess_plot_df, aes(x = factor_added)) +
    geom_hline(yintercept = 0, colour = "black", linewidth = 0.3) +
    geom_col(aes(y = excess_decrease, fill = excess_direction), width = 0.45) +
    geom_errorbar(
        aes(ymin = null_range_low, ymax = null_range_high, colour = "Permutation range (min to max)"),
        width = 0.12, linewidth = 0.3, position = position_nudge(x = 0.33)
    ) +
    geom_text(aes(y = excess_decrease, label = excess_label, vjust = label_vjust), size = 7 / .pt, colour = "black") +
    geom_vline(xintercept = selected_rank_position, linetype = "dashed", colour = "#BC4749", linewidth = 0.4) +
    annotate(
        "text", x = selected_rank_position, y = Inf, label = "Selected rank 5",
        hjust = -0.05, vjust = 1.5, size = 7 / .pt, colour = "#BC4749"
    ) +
    scale_fill_manual(
        values = c("Above permutation mean" = "#BC4749", "Below permutation mean" = "grey55"),
        name = NULL
    ) +
    scale_colour_manual(values = c("Permutation range (min to max)" = "grey40"), name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0.12, 0.15))) +
    labs(
        x = "Factorization rank (from → to)",
        y = "Excess decrease in RSS\n(observed - permutation mean)"
    ) +
    theme_bw(base_size = 7) +
    theme(
        text = element_text(size = 7, colour = "black"),
        axis.text = element_text(size = 7, colour = "black"),
        legend.position = "bottom",
        legend.text = element_text(size = 7, colour = "black"),
        legend.key.size = unit(3, "mm"),
        legend.margin = margin(0, 0, 0, 0),
        panel.grid = element_blank()
    )

print(p)

## 3. 保存图形 ####
ggsave(
    filename = file.path(script_results_path, "图8_增加各因子的RSS额外下降量_相对置换零分布.pdf"),
    plot = p,
    width = 14.65,
    height = 6,
    units = "cm",
    device = cairo_pdf
)


# 9图：Cancer-level consensus silhouette widths at ranks 2-8 ####

## 1. 读取并处理绘图数据 ####
silhouette_plot_df <- readxl::read_xlsx(input_xlsx_file, sheet = "各秩癌种分组") %>%
    dplyr::mutate(
        rank_label = factor(paste0("Rank ", rank), levels = paste0("Rank ", 2:8)),
        reference_group = factor(reference_group, levels = c("Ocular", "Hematologic", "Neurologic", "Systemic-like", "Visceral-like")),
        low_silhouette_label = ifelse(silhouette < 0.5, cancer_type, NA_character_)
    )

silhouette_mean_df <- silhouette_plot_df %>%
    dplyr::group_by(rank_label) %>%
    dplyr::summarise(mean_silhouette = mean(silhouette), .groups = "drop")

## 2. 绘图 ####
p <- ggplot(silhouette_plot_df, aes(x = rank_label, y = silhouette)) +
    geom_vline(xintercept = seq(1.5, 6.5, by = 1), colour = "grey75", linewidth = 0.3) +
    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey40", linewidth = 0.3) +
    geom_crossbar(
        data = silhouette_mean_df,
        aes(x = rank_label, y = mean_silhouette, ymin = mean_silhouette, ymax = mean_silhouette),
        width = 0.5, colour = "black", linewidth = 0.3, inherit.aes = FALSE
    ) +
    geom_point(
        aes(fill = reference_group),
        shape = 21, size = 1.6, colour = "grey20", stroke = 0.2,
        position = position_jitter(width = 0.15, height = 0, seed = 123456)
    ) +
    ggrepel::geom_text_repel(
        aes(label = low_silhouette_label),
        size = 7 / .pt, colour = "black", na.rm = TRUE,
        position = position_jitter(width = 0.15, height = 0, seed = 123456),
        min.segment.length = 0, segment.size = 0.2, box.padding = 0.15, max.overlaps = Inf
    ) +
    scale_fill_manual(
        values = c("Ocular" = "#DF7A5E", "Hematologic" = "#A15F88", "Neurologic" = "#82B29A", "Systemic-like" = "#F2CC8E", "Visceral-like" = "#BC4749"),
        name = "Final group"
    ) +
    scale_y_continuous(limits = c(-0.05, 1.02), breaks = seq(0, 1, 0.25)) +
    labs(
        x = NULL,
        y = "Silhouette width"
    ) +
    theme_bw(base_size = 7) +
    theme(
        text = element_text(size = 7, colour = "black"),
        axis.text = element_text(size = 7, colour = "black"),
        plot.title = element_text(size = 7, colour = "black", face = "bold"),
        plot.subtitle = element_text(size = 7, colour = "black"),
        legend.title = element_text(size = 7, colour = "black"),
        legend.text = element_text(size = 7, colour = "black"),
        legend.key.size = unit(3, "mm"),
        legend.position = "bottom",
        legend.margin = margin(0, 0, 0, 0),
        panel.grid = element_blank()
    )

print(p)

## 3. 保存图形 ####
ggsave(
    filename = file.path(script_results_path, "图9_秩2至8癌种轮廓系数.pdf"),
    plot = p,
    width = 14.65,
    height = 7.5,
    units = "cm",
    device = cairo_pdf
)


# 10图：Summary of clustering quality metrics at ranks 2-8 ####

## 1. 读取并处理绘图数据 ####
quality_plot_df <- readxl::read_xlsx(input_xlsx_file, sheet = "各秩聚类质量汇总") %>%
    dplyr::select(
        rank, cophenetic, mean_silhouette, min_silhouette,
        n_cancers_silhouette_below_threshold, n_singleton_clusters
    ) %>%
    tidyr::pivot_longer(-rank, names_to = "metric", values_to = "value") %>%
    dplyr::mutate(
        metric = factor(
            metric,
            levels = c("cophenetic", "mean_silhouette", "min_silhouette", "n_cancers_silhouette_below_threshold", "n_singleton_clusters"),
            labels = c(
                "Cophenetic correlation\n(higher is better)",
                "Mean silhouette\n(higher is better)",
                "Minimum silhouette\n(higher is better)",
                "Cancer types with silhouette < 0.5\n(lower is better)",
                "Singleton clusters\n(lower is better)"
            )
        ),
        is_selected_rank = rank == 5,
        value_label = ifelse(value %% 1 == 0, sprintf("%.0f", value), sprintf("%.3f", value))
    ) %>%
    dplyr::group_by(metric) %>%
    dplyr::arrange(rank, .by_group = TRUE) %>%
    dplyr::mutate(
        neighbour_mean = (dplyr::coalesce(dplyr::lag(value), dplyr::lead(value)) + dplyr::coalesce(dplyr::lead(value), dplyr::lag(value))) / 2,
        previous_value = dplyr::lag(value),
        next_value = dplyr::lead(value),
        label_position = dplyr::case_when(
            !is.na(previous_value) & !is.na(next_value) & previous_value < value & value < next_value ~ "rising",
            !is.na(previous_value) & !is.na(next_value) & previous_value > value & value > next_value ~ "falling",
            value >= neighbour_mean ~ "above",
            TRUE ~ "below"
        ),
        label_vjust = dplyr::case_when(
            label_position %in% c("rising", "falling") ~ -0.4,
            label_position == "above" ~ -0.9,
            TRUE ~ 1.9
        ),
        label_hjust = dplyr::case_when(
            label_position == "rising" ~ 1.1,
            label_position == "falling" ~ -0.1,
            TRUE ~ 0.5
        )
    ) %>%
    dplyr::ungroup()

## 2. 绘图 ####
p <- ggplot(quality_plot_df, aes(x = rank, y = value)) +
    geom_line(colour = "grey50", linewidth = 0.4) +
    geom_point(aes(colour = is_selected_rank), size = 1.6) +
    geom_text(aes(label = value_label, colour = is_selected_rank, vjust = label_vjust, hjust = label_hjust), size = 7 / .pt) +
    facet_wrap(~metric, nrow = 2, scales = "free_y", axes = "all_x") +
    scale_x_continuous(breaks = 2:8, expand = expansion(mult = 0.15)) +
    scale_y_continuous(
        expand = expansion(mult = c(0.25, 0.3)),
        breaks = function(limits) {
            if (max(limits) > 1.5) seq(max(0, ceiling(min(limits))), floor(max(limits)), by = 1) else scales::breaks_pretty()(limits)
        }
    ) +
    scale_colour_manual(values = c("TRUE" = "#BC4749", "FALSE" = "grey30"), guide = "none") +
    labs(
        x = "Factorization rank",
        y = NULL
    ) +
    coord_cartesian(clip = "off") +
    theme_bw(base_size = 7) +
    theme(
        text = element_text(size = 7, colour = "black"),
        axis.text = element_text(size = 7, colour = "black"),
        plot.title = element_text(size = 7, colour = "black", face = "bold"),
        strip.text = element_text(size = 7, colour = "black"),
        strip.background = element_rect(fill = "grey95", colour = "grey60"),
        panel.grid = element_blank()
    )

print(p)

## 3. 保存图形 ####
ggsave(
    filename = file.path(script_results_path, "图10_秩2至8聚类质量指标汇总.pdf"),
    plot = p,
    width = 14.65,
    height = 9,
    units = "cm",
    device = cairo_pdf
)

message("已生成图形文件，目录: ", script_results_path)
message("========== NMF秩选择敏感性绘图结束 ==========")
