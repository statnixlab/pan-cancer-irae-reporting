# 7.4.1 肿瘤免疫微环境绘图 ####
# 用途：读取7.4导出的免疫特征分析结果，按癌种群展示feature-axis heatmap、核心特征、方向稳定性、effect size和混合模型敏感性结果。
#       所有绘图数据均从已保存的Excel结果表读取，不依赖7.4分析阶段的内存对象。 ####

rm(list = ls())
gc()


# 1. 加载配置 ####
source("Code/config.R")
library(ggh4x)


# 2. 参数与路径 ####
script_name <- "7.4.1肿瘤免疫微环境绘图"
script_results_path <- file.path(results_path, script_name)
figure_output_path <- file.path(script_results_path, "figures")

dir.create(script_results_path, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_output_path, recursive = TRUE, showWarnings = FALSE)

input_results_path <- file.path(results_path, "7.4肿瘤免疫微环境分析")
table_output_path <- file.path(input_results_path, "tables")

immune_result_xlsx <- file.path(table_output_path, "TCGA肿瘤免疫微环境分析结果.xlsx")
stopifnot(file.exists(immune_result_xlsx))

pattern_group_order <- c(
  "Ocular",
  "Hematologic",
  "Neurologic",
  "Systemic-like",
  "Visceral-like"
)

pattern_group_color_map <- c(
  "Ocular"        = "#DF7A5E",
  "Hematologic"   = "#82B29A",
  "Neurologic"    = "#A15F88",
  "Systemic-like" = "#F2CC8E",
  "Visceral-like" = "#BC4749"
)

heatmap_color_low <- "#2C7BB6"
heatmap_color_mid <- "#F7F7F7"
heatmap_color_high <- "#D7191C"

core_boxplot_features <- c(
  "TMB",
  "TCR_Shannon",
  "Hallmark_HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "Hallmark_HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "CYT_score",
  "T_cell_inflamed_GEP",
  "Hallmark_HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "Antigen_Presentation",
  "PDCD1_expression",
  "CD274_expression",
  "CTLA4_expression"
)

checkpoint_features <- c(
  "PDCD1_expression",
  "CD274_expression",
  "CTLA4_expression"
)


# 3. 绘图辅助函数 ####

## 函数：prepare_group_heatmap_data ####
## 功能：将pattern group特征汇总表转换为heatmap数据，并对每个feature做行z-score ####
## 输入：summary_df，模式群特征汇总表；feature_info_df，特征信息表；target_axis，目标feature axis，可为空 ####
## 输出：用于ggplot heatmap的数据框 ####
prepare_group_heatmap_data <- function(summary_df, feature_info_df, target_axis = NULL) {
  feature_order_df <- feature_info_df %>%
    dplyr::filter(main_display == TRUE) %>%
    dplyr::arrange(axis_order, feature_order)

  if (!is.null(target_axis)) {
    feature_order_df <- feature_order_df %>%
      dplyr::filter(feature_axis %in% target_axis)
  }

  available_features <- intersect(feature_order_df$feature_name, summary_df$feature_name)

  summary_df %>%
    dplyr::filter(feature_name %in% available_features) %>%
    dplyr::left_join(
      feature_order_df %>%
        dplyr::select(feature_name, feature_label, feature_axis, axis_order, feature_order),
      by = "feature_name",
      suffix = c("", ".feature_info")
    ) %>%
    dplyr::mutate(
      feature_label = dplyr::coalesce(feature_label.feature_info, feature_label),
      feature_axis = dplyr::coalesce(feature_axis.feature_info, feature_axis),
      axis_order = dplyr::coalesce(axis_order.feature_info, axis_order),
      feature_order = dplyr::coalesce(feature_order.feature_info, feature_order)
    ) %>%
    dplyr::select(-dplyr::any_of(c(
      "feature_label.feature_info",
      "feature_axis.feature_info",
      "axis_order.feature_info",
      "feature_order.feature_info"
    ))) %>%
    dplyr::group_by(feature_name) %>%
    dplyr::mutate(
      feature_mean = mean(median_score, na.rm = TRUE),
      feature_sd = stats::sd(median_score, na.rm = TRUE),
      z_score = dplyr::case_when(
        is.na(feature_sd) | feature_sd == 0 ~ NA_real_,
        TRUE ~ (median_score - feature_mean) / feature_sd
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-feature_mean, -feature_sd) %>%
    dplyr::mutate(
      pattern_group = factor(pattern_group, levels = pattern_group_order),
      feature_axis = factor(feature_axis, levels = unique(feature_order_df$feature_axis)),
      feature_label = factor(
        feature_label,
        levels = rev(feature_order_df$feature_label[match(available_features, feature_order_df$feature_name)])
      )
    )
}

## 函数：prepare_primary_organ_long_data ####
## 功能：将primary organ级宽表整理为按癌种群展示的长表 ####
## 输入：score_df，原发部位级特征表；feature_info_df，特征信息表；feature_order，目标特征顺序 ####
## 输出：用于箱线图或点图的数据框 ####
prepare_primary_organ_long_data <- function(score_df, feature_info_df, feature_order) {
  available_features <- intersect(feature_order, colnames(score_df))

  score_df %>%
    dplyr::select(pattern_group, primary_organ, n_tumor_samples, dplyr::any_of(available_features)) %>%
    tidyr::pivot_longer(
      cols = dplyr::any_of(available_features),
      names_to = "feature_name",
      values_to = "score"
    ) %>%
    dplyr::left_join(feature_info_df, by = "feature_name") %>%
    dplyr::mutate(
      pattern_group = factor(pattern_group, levels = pattern_group_order),
      feature_label = factor(
        feature_label,
        levels = feature_info_df$feature_label[match(available_features, feature_info_df$feature_name)]
      )
    )
}


# 1图：癌种群Feature-axis Heatmap ####

## 1. 读取并处理绘图数据 ####
pattern_group_summary_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "模式群特征汇总"
)

feature_info_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "特征信息表"
)

heatmap_data <- prepare_group_heatmap_data(
  summary_df = pattern_group_summary_df,
  feature_info_df = feature_info_df
)

## 单 feature 一组的轴 strip 透明（含边框），多 feature 的轴显示灰底黑字竖排 strip ####
axis_levels_in_data <- levels(heatmap_data$feature_axis)

axis_feature_count <- heatmap_data %>%
  dplyr::distinct(feature_axis, feature_label) %>%
  dplyr::count(feature_axis, name = "n_features") %>%
  dplyr::mutate(feature_axis = as.character(feature_axis))

multi_feature_axes <- axis_feature_count$feature_axis[axis_feature_count$n_features > 1]

strip_fill_vec <- ifelse(
  axis_levels_in_data %in% multi_feature_axes,
  "grey85",
  "transparent"
)

strip_color_vec <- ifelse(
  axis_levels_in_data %in% multi_feature_axes,
  "grey50",
  "transparent"
)

strip_text_color_vec <- ifelse(
  axis_levels_in_data %in% multi_feature_axes,
  "grey20",
  "transparent"
)

heatmap_strip <- ggh4x::strip_themed(
  background_y = ggh4x::elem_list_rect(fill = strip_fill_vec, color = strip_color_vec),
  text_y = ggh4x::elem_list_text(color = strip_text_color_vec, face = "bold", size = 4, angle = 90)
)

## 手动指定每个生物轴 label：feature 多的轴保留单行，feature 少的主动换行避免截断 ####
axis_label_lookup <- c(
  "Tumor Antigenicity"                  = "Antigenicity",
  "TCR Repertoire"                      = "TCR\nRepertoire",
  "Cytokine and Inflammatory Signaling" = "Cytokine\nSignaling",
  "T Cell Activation and Cytotoxicity"  = "T cell\nCytotoxicity",
  "Antigen Presentation Machinery"      = "Antigen\nPresentation",
  "Key Immune Cell Infiltration"        = "Key Immune Cell Infiltration",
  "Checkpoint Gene Expression"          = "Checkpoint\nExpression",
  "TME Context"                         = "TME Context"
)

## 2. 绘图 ####
p <- ggplot(heatmap_data, aes(x = pattern_group, y = feature_label, fill = z_score)) +
  geom_tile(color = "white", linewidth = 0.25) +
  ggh4x::facet_grid2(
    feature_axis ~ .,
    scales = "free_y",
    space = "free_y",
    strip = heatmap_strip,
    labeller = labeller(feature_axis = function(x) unname(axis_label_lookup[as.character(x)]))
  ) +
  scale_fill_gradient2(
    low = heatmap_color_low,
    mid = heatmap_color_mid,
    high = heatmap_color_high,
    midpoint = 0,
    na.value = "grey90",
    name = "Row z-score"
  ) +
  guides(
    fill = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      barwidth = unit(2.6, "cm"),
      barheight = unit(0.22, "cm")
    )
  ) +
  labs(
    #title = "Immune Feature Patterns across Cancer Groups",
    x = "Cancer Group",
    y = NULL
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 7, face = "bold", color = "black"),
    axis.title = element_text(size = 7, color = "black"),
    axis.text.x = element_text(size = 7, color = "black", angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 7, color = "black"),
    strip.placement = "outside",
    strip.text.y = element_text(size = 4, face = "bold"),
    legend.position = "top",
    legend.justification = "center",
    legend.title = element_text(size = 7, color = "black"),
    legend.text = element_text(size = 7, color = "black"),
    legend.key.width = unit(0.55, "cm"),
    legend.key.height = unit(0.18, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    panel.grid = element_blank()
  )

print(p)

## 3. 保存图形 ####
ggsave(
  filename = file.path(figure_output_path, "图1_癌种群免疫特征生物轴热图.pdf"),
  plot = p,
  width = 7.5,
  height = 12,
  units = "cm",
  device = cairo_pdf
)


# 2图：核心免疫特征癌种群箱线图 ####

## 1. 读取并处理绘图数据 ####
primary_organ_score_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "原发部位级免疫特征"
)

feature_info_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "特征信息表"
)

all_main_display_features <- feature_info_df %>%
  dplyr::filter(main_display == TRUE) %>%
  dplyr::arrange(axis_order, feature_order) %>%
  dplyr::pull(feature_name)

plot_data <- prepare_primary_organ_long_data(
  score_df = primary_organ_score_df,
  feature_info_df = feature_info_df,
  feature_order = all_main_display_features
)

## 2. 绘图 ####
p <- ggplot(plot_data, aes(x = pattern_group, y = score, fill = pattern_group)) +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.75) +
  geom_jitter(aes(size = n_tumor_samples), width = 0.16, alpha = 0.8, color = "black") +
  facet_wrap(~ feature_label, scales = "free_y", ncol = 4) +
  scale_fill_manual(name = "Pattern group", values = pattern_group_color_map, drop = FALSE) +
  scale_size_continuous(range = c(1.2, 3.5), name = "Tumor samples") +
  labs(
    #title = "Immune Features across Cancer Groups",
    x = NULL,
    y = "Primary Organ-level Score"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold", color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 14, color = "black"),
    legend.text = element_text(size = 12, color = "black"),
    strip.text = element_text(size = 12, face = "bold", color = "black")
  )

print(p)

## 3. 保存图形 ####
ggsave(
  filename = file.path(figure_output_path, "图2_核心免疫特征癌种群箱线图.pdf"),
  plot = p,
  width = 18,
  height = 18,
  device = cairo_pdf
)


# 3图：Checkpoint基因表达癌种群箱线图 ####

## 1. 读取并处理绘图数据 ####
primary_organ_score_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "原发部位级免疫特征"
)

feature_info_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "特征信息表"
)

plot_data <- prepare_primary_organ_long_data(
  score_df = primary_organ_score_df,
  feature_info_df = feature_info_df,
  feature_order = checkpoint_features
)

## 2. 绘图 ####
p <- ggplot(plot_data, aes(x = pattern_group, y = score, fill = pattern_group)) +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.75) +
  geom_jitter(width = 0.16, size = 1.8, alpha = 0.85, color = "black") +
  facet_wrap(~ feature_label, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = pattern_group_color_map, drop = FALSE) +
  labs(
    title = "Checkpoint Gene Expression across Cancer Groups",
    x = "Cancer Group",
    y = "Primary Organ-level log2(TPM + 1)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold", color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 12, color = "black"),
    legend.position = "none",
    strip.text = element_text(size = 12, face = "bold", color = "black")
  )

print(p)

## 3. 保存图形 ####
ggsave(
  filename = file.path(figure_output_path, "图3_免疫检查点基因表达癌种群箱线图.pdf"),
  plot = p,
  width = 11,
  height = 5.5,
  device = cairo_pdf
)


# 4图：关键免疫细胞浸润癌种群热图 ####

## 1. 读取并处理绘图数据 ####
pattern_group_summary_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "模式群特征汇总"
)

feature_info_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "特征信息表"
)

immune_cell_heatmap_data <- prepare_group_heatmap_data(
  summary_df = pattern_group_summary_df,
  feature_info_df = feature_info_df,
  target_axis = "Key Immune Cell Infiltration"
)

## 2. 绘图 ####
p <- ggplot(immune_cell_heatmap_data, aes(x = pattern_group, y = feature_label, fill = z_score)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient2(
    low = heatmap_color_low,
    mid = heatmap_color_mid,
    high = heatmap_color_high,
    midpoint = 0,
    na.value = "grey90",
    name = "Row z-score"
  ) +
  labs(
    title = "Key Immune Cell Infiltration across Cancer Groups",
    x = "Cancer Group",
    y = "Immune Cell"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold", color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 14, color = "black"),
    legend.text = element_text(size = 12, color = "black"),
    panel.grid = element_blank()
  )

print(p)

## 3. 保存图形 ####
ggsave(
  filename = file.path(figure_output_path, "图4_关键免疫细胞浸润癌种群热图.pdf"),
  plot = p,
  width = 8,
  height = 5.5,
  device = cairo_pdf
)


# 5图：Bootstrap最高癌种群稳定性Dotplot ####

## 1. 读取并处理绘图数据 ####
bootstrap_stability_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "抽样特征稳定性"
)

feature_info_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "特征信息表"
)

plot_data <- bootstrap_stability_df %>%
  dplyr::left_join(
    feature_info_df %>%
      dplyr::select(feature_name, feature_label, feature_axis, axis_order, feature_order, main_display),
    by = "feature_name",
    suffix = c("", ".feature_info")
  ) %>%
  dplyr::mutate(
    feature_label = dplyr::coalesce(feature_label.feature_info, feature_label),
    feature_axis = dplyr::coalesce(feature_axis.feature_info, feature_axis),
    axis_order = dplyr::coalesce(axis_order.feature_info, axis_order),
    feature_order = dplyr::coalesce(feature_order.feature_info, feature_order),
    main_display = dplyr::coalesce(main_display.feature_info, main_display)
  ) %>%
  dplyr::filter(main_display == TRUE, !is.na(highest_group_probability)) %>%
  dplyr::arrange(axis_order, feature_order) %>%
  dplyr::mutate(
    feature_label = factor(feature_label, levels = rev(unique(feature_label))),
    feature_axis = factor(feature_axis, levels = unique(feature_axis)),
    most_often_highest_group = factor(most_often_highest_group, levels = pattern_group_order)
  )

## 2. 绘图 ####
p <- ggplot(plot_data, aes(x = most_often_highest_group, y = feature_label)) +
  geom_point(aes(size = highest_group_probability, color = bootstrap_median_range), alpha = 0.85) +
  scale_color_gradient(low = "#4A90D9", high = "#D94A4A", name = "Bootstrap Median Range") +
  scale_size_continuous(range = c(2.5, 8), limits = c(0, 1), name = "Probability Highest") +
  labs(
    title = "Bootstrap Stability of Highest Cancer Group",
    x = "Most Often Highest Cancer Group",
    y = "Feature"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold", color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 14, color = "black"),
    legend.text = element_text(size = 12, color = "black")
  )

print(p)

## 3. 保存图形 ####
ggsave(
  filename = file.path(figure_output_path, "图5_Bootstrap最高癌种群稳定性点图.pdf"),
  plot = p,
  width = 10,
  height = 10,
  device = cairo_pdf
)


# 6图：TCGA项目级补充癌种群点图 ####

## 1. 读取并处理绘图数据 ####
project_score_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "TCGA项目级补充特征"
)

feature_info_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "特征信息表"
)

available_features <- intersect(core_boxplot_features, colnames(project_score_df))

plot_data <- project_score_df %>%
  dplyr::select(pattern_group, primary_organ, project_id, tcga_code, n_tumor_samples, dplyr::any_of(available_features)) %>%
  tidyr::pivot_longer(
    cols = dplyr::any_of(available_features),
    names_to = "feature_name",
    values_to = "score"
  ) %>%
  dplyr::left_join(feature_info_df, by = "feature_name") %>%
  dplyr::mutate(
    pattern_group = factor(pattern_group, levels = pattern_group_order),
    feature_label = factor(
      feature_label,
      levels = feature_info_df$feature_label[match(available_features, feature_info_df$feature_name)]
    )
  )

## 2. 绘图 ####
p <- ggplot(plot_data, aes(x = pattern_group, y = score, color = pattern_group)) +
  geom_point(aes(size = n_tumor_samples), position = position_jitter(width = 0.16, height = 0), alpha = 0.85) +
  facet_wrap(~ feature_label, scales = "free_y", ncol = 3) +
  scale_color_manual(values = pattern_group_color_map, drop = FALSE) +
  scale_size_continuous(range = c(1.3, 3.8), name = "Tumor samples") +
  labs(
    title = "TCGA Project-level Supplementary Immune Features",
    x = "Cancer Group",
    y = "TCGA Project-level Score"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold", color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 14, color = "black"),
    legend.text = element_text(size = 12, color = "black"),
    strip.text = element_text(size = 12, face = "bold", color = "black")
  )

print(p)

## 3. 保存图形 ####
ggsave(
  filename = file.path(figure_output_path, "图6_TCGA项目级补充癌种群点图.pdf"),
  plot = p,
  width = 14,
  height = 12,
  device = cairo_pdf
)


# 7图：Bootstrap最高癌种群与效应量描述点图 ####

## 1. 读取并处理绘图数据 ####
direction_effect_size_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "方向效应量汇总"
)

plot_data <- direction_effect_size_df %>%
  dplyr::filter(
    main_display == TRUE,
    !is.na(most_often_highest_group),
    !is.na(bootstrap_median_range),
    !is.na(primary_organ_score_sd),
    primary_organ_score_sd > 0
  ) %>%
  dplyr::arrange(axis_order, feature_order) %>%
  dplyr::mutate(
    feature_label = factor(feature_label, levels = rev(unique(feature_label))),
    most_often_highest_group = factor(most_often_highest_group, levels = pattern_group_order),
    bootstrap_std_median_range = bootstrap_median_range / primary_organ_score_sd,
    interpretation_priority = factor(
      interpretation_priority,
      levels = c("Primary interpretation", "Secondary interpretation", "Direction only", "Exploratory only")
    )
  )

## 2. 绘图 ####
p <- ggplot(plot_data, aes(x = most_often_highest_group, y = feature_label)) +
  geom_point(aes(size = bootstrap_std_median_range, color = most_often_highest_group), alpha = 0.85) +
  scale_color_manual(
    values = pattern_group_color_map,
    drop = FALSE,
    name = "Highest\nCancer Group"
  ) +
  scale_size_continuous(
    range = c(1.2, 3.2),
    name = "Bootstrap\nStandardized Range"
  ) +
  guides(
    color = "none",
    size = guide_legend(order = 1, nrow = 1)
  ) +
  labs(
    #title = "Direction Stability and Effect Size across Cancer Groups",
    x = "Most Often Highest Cancer Group\n(Bootstrap-estimated)",
    y = NULL
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 7, face = "bold", color = "black"),
    axis.title = element_text(size = 7, color = "black"),
    axis.title.x = element_text(size = 7, color = "black", hjust = 0.35),
    axis.text.x = element_text(size = 7, color = "black", angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 7, color = "black"),
    legend.position = "top",
    legend.justification = "left",
    legend.direction = "horizontal",
    legend.box = "vertical",
    legend.box.just = "left",
    legend.background = element_blank(),
    legend.title = element_text(size = 7, color = "black"),
    legend.text = element_text(size = 7, color = "black"),
    legend.key.size = unit(0.22, "cm"),
    legend.spacing.x = unit(0.06, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, -60),
    plot.margin = margin(2, 8, 2, 2)
  )

print(p)

## 3. 保存图形 ####
ggsave(
  filename = file.path(figure_output_path, "图7_Bootstrap最高癌种群与效应量描述点图.pdf"),
  plot = p,
  width = 7.5,
  height = 12,
  units = "cm",
  device = cairo_pdf
)


# 8图：LMM最高癌种群与Bootstrap一致性点图 ####

## 1. 读取并处理绘图数据 ####
direction_effect_size_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "方向效应量汇总"
)

plot_data <- direction_effect_size_df %>%
  dplyr::filter(
    main_display == TRUE,
    mixed_model_status == "success",
    !is.na(mixed_highest_group)
  ) %>%
  dplyr::arrange(axis_order, feature_order) %>%
  dplyr::mutate(
    feature_label = factor(feature_label, levels = rev(unique(feature_label))),
    mixed_highest_group = factor(mixed_highest_group, levels = pattern_group_order),
    mixed_model_direction_agreement = factor(
      mixed_model_direction_agreement,
      levels = c("Mixed model supported", "Mixed model discordant"),
      labels = c("Same", "Different")
    )
  )

## 2. 绘图 ####
p <- ggplot(plot_data, aes(x = mixed_highest_group, y = feature_label)) +
  geom_point(aes(size = mixed_standardized_range, color = mixed_model_direction_agreement), alpha = 0.85) +
  scale_color_manual(
    values = c(
      "Same" = "#911F23",
      "Different" = "#8CB4C5"
    ),
    drop = FALSE,
    name = "LMM vs Bootstrap Highest Group"
  ) +
  scale_size_continuous(range = c(1.2, 3.2), name = "LMM Std. EMM Range") +
  guides(
    color = guide_legend(order = 1, ncol = 1, override.aes = list(size = 1.6)),
    size = guide_legend(order = 2)
  ) +
  labs(
    #title = "Mixed Model Direction Agreement",
    x = "Highest Cancer Group\n(LMM-estimated)",
    y = NULL
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 7, face = "bold", color = "black"),
    axis.title = element_text(size = 7, color = "black"),
    axis.title.x = element_text(size = 7, color = "black", hjust = 0.35),
    axis.text.x = element_text(size = 7, color = "black", angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 7, color = "black"),
    legend.position = "right",
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.background = element_blank(),
    legend.title = element_text(size = 7, color = "black"),
    legend.text = element_text(size = 7, color = "black"),
    legend.key.size = unit(0.22, "cm"),
    legend.spacing.x = unit(0.06, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    plot.margin = margin(2, 2, 2, 2)
  )

print(p)

## 3. 保存图形 ####
ggsave(
  filename = file.path(figure_output_path, "图8_LMM最高癌种群与Bootstrap一致性点图.pdf"),
  plot = p,
  width = 14.65,
  height = 12,
  units = "cm",
  device = cairo_pdf
)


# 8.1图：癌种群免疫特征生物轴热图、Bootstrap效应量和LMM一致性横向拼接图 ####

## 1. 读取并处理绘图数据 ####
combined_main_font_size <- 7
combined_axis_font_size <- 7
combined_legend_font_size <- 7

pattern_group_summary_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "模式群特征汇总"
)

feature_info_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "特征信息表"
)

combined_heatmap_data <- prepare_group_heatmap_data(
  summary_df = pattern_group_summary_df,
  feature_info_df = feature_info_df
)

combined_axis_levels_in_data <- levels(combined_heatmap_data$feature_axis)
combined_feature_label_levels <- levels(combined_heatmap_data$feature_label)

combined_axis_feature_count <- combined_heatmap_data %>%
  dplyr::distinct(feature_axis, feature_label) %>%
  dplyr::count(feature_axis, name = "n_features") %>%
  dplyr::mutate(feature_axis = as.character(feature_axis))

combined_multi_feature_axes <- combined_axis_feature_count$feature_axis[
  combined_axis_feature_count$n_features > 1
]

combined_strip_fill_vec <- ifelse(
  combined_axis_levels_in_data %in% combined_multi_feature_axes,
  "grey85",
  "transparent"
)

combined_strip_color_vec <- ifelse(
  combined_axis_levels_in_data %in% combined_multi_feature_axes,
  "grey50",
  "transparent"
)

combined_strip_text_color_vec <- ifelse(
  combined_axis_levels_in_data %in% combined_multi_feature_axes,
  "grey20",
  "transparent"
)

combined_heatmap_strip <- ggh4x::strip_themed(
  background_y = ggh4x::elem_list_rect(fill = combined_strip_fill_vec, color = combined_strip_color_vec),
  text_y = ggh4x::elem_list_text(
    color = combined_strip_text_color_vec,
    face = "bold",
    size = combined_axis_font_size,
    angle = 90
  )
)

combined_axis_label_lookup <- c(
  "Tumor Antigenicity"                  = "Antigenicity",
  "TCR Repertoire"                      = "TCR\nRepertoire",
  "Cytokine and Inflammatory Signaling" = "Inflammatory\nSignaling",
  "T Cell Activation and Cytotoxicity"  = "T Cell\nCytotoxicity",
  "Antigen Presentation Machinery"      = "Antigen\nPresentation",
  "Key Immune Cell Infiltration"        = "Immune Cell\nInfiltration",
  "Checkpoint Gene Expression"          = "Checkpoint\nExpression",
  "TME Context"                         = "TME Context"
)

direction_effect_size_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "方向效应量汇总"
)

combined_bootstrap_plot_data <- direction_effect_size_df %>%
  dplyr::filter(
    main_display == TRUE,
    !is.na(most_often_highest_group),
    !is.na(bootstrap_median_range),
    !is.na(primary_organ_score_sd),
    primary_organ_score_sd > 0
  ) %>%
  dplyr::arrange(axis_order, feature_order) %>%
  dplyr::mutate(
    feature_axis = factor(feature_axis, levels = combined_axis_levels_in_data),
    feature_label = factor(feature_label, levels = combined_feature_label_levels),
    most_often_highest_group = factor(most_often_highest_group, levels = pattern_group_order),
    bootstrap_std_median_range = bootstrap_median_range / primary_organ_score_sd,
    interpretation_priority = factor(
      interpretation_priority,
      levels = c("Primary interpretation", "Secondary interpretation", "Direction only", "Exploratory only")
    )
  )

combined_lmm_plot_data <- direction_effect_size_df %>%
  dplyr::filter(
    main_display == TRUE,
    mixed_model_status == "success",
    !is.na(mixed_highest_group)
  ) %>%
  dplyr::arrange(axis_order, feature_order) %>%
  dplyr::mutate(
    feature_axis = factor(feature_axis, levels = combined_axis_levels_in_data),
    feature_label = factor(feature_label, levels = combined_feature_label_levels),
    mixed_highest_group = factor(mixed_highest_group, levels = pattern_group_order),
    mixed_model_direction_agreement = factor(
      mixed_model_direction_agreement,
      levels = c("Mixed model supported", "Mixed model discordant"),
      labels = c("Same", "Different")
    )
  )

## 2. 绘图 ####
p_heatmap <- ggplot(combined_heatmap_data, aes(x = pattern_group, y = feature_label, fill = z_score)) +
  geom_tile(color = "white", linewidth = 0.18) +
  ggh4x::facet_grid2(
    feature_axis ~ .,
    scales = "free_y",
    space = "free_y",
    strip = combined_heatmap_strip,
    labeller = labeller(feature_axis = function(x) unname(combined_axis_label_lookup[as.character(x)]))
  ) +
  scale_fill_gradient2(
    low = heatmap_color_low,
    mid = heatmap_color_mid,
    high = heatmap_color_high,
    midpoint = 0,
    na.value = "grey90",
    name = "Row z-score"
  ) +
  guides(
    fill = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      barwidth = unit(2.2, "cm"),
      barheight = unit(0.22, "cm")
    )
  ) +
  labs(
    x = "Cancer Group",
    y = "Feature"
  ) +
  theme_bw(base_size = combined_main_font_size) +
  theme(
    plot.title = element_text(size = combined_main_font_size, face = "bold", color = "black"),
    axis.title = element_text(size = combined_main_font_size, color = "black"),
    axis.text.x = element_text(size = combined_axis_font_size, color = "black", angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = combined_axis_font_size, color = "black"),
    strip.placement = "outside",
    strip.text.y = element_text(size = combined_axis_font_size, face = "bold", color = "black"),
    legend.position = "bottom",
    legend.title = element_text(size = combined_legend_font_size, color = "black"),
    legend.text = element_text(size = combined_legend_font_size, color = "black"),
    legend.key.size = unit(0.28, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    plot.margin = margin(2, 2, 2, 2),
    panel.spacing.y = unit(0.45, "mm"),
    panel.grid = element_blank()
  )

p_bootstrap <- ggplot(combined_bootstrap_plot_data, aes(x = most_often_highest_group, y = feature_label)) +
  geom_point(aes(size = bootstrap_std_median_range, color = most_often_highest_group), alpha = 0.85) +
  ggh4x::facet_grid2(
    feature_axis ~ .,
    scales = "free_y",
    space = "free_y",
    drop = FALSE
  ) +
  scale_color_manual(
    values = pattern_group_color_map,
    drop = FALSE,
    name = "Highest\nCancer Group"
  ) +
  scale_size_continuous(
    range = c(0.9, 2.6),
    name = "Bootstrap\nStd. Range"
  ) +
  guides(
    color = "none",
    size = guide_legend(order = 1, override.aes = list(color = "grey20"))
  ) +
  labs(
    x = "Bootstrap highest\ngroup",
    y = NULL
  ) +
  theme_bw(base_size = combined_main_font_size) +
  theme(
    plot.title = element_text(size = combined_main_font_size, face = "bold", color = "black"),
    axis.title = element_text(size = combined_main_font_size, color = "black"),
    axis.text.x = element_text(size = combined_axis_font_size, color = "black", angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.background = element_blank(),
    strip.background.y = element_blank(),
    strip.text = element_blank(),
    strip.text.y = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = combined_legend_font_size, color = "black"),
    legend.text = element_text(size = combined_legend_font_size, color = "black"),
    legend.key.size = unit(0.28, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    plot.margin = margin(2, 2, 2, 2),
    panel.spacing.y = unit(0.45, "mm")
  )

p_lmm <- ggplot(combined_lmm_plot_data, aes(x = mixed_highest_group, y = feature_label)) +
  geom_point(aes(size = mixed_standardized_range, color = mixed_model_direction_agreement), alpha = 0.85) +
  ggh4x::facet_grid2(
    feature_axis ~ .,
    scales = "free_y",
    space = "free_y",
    drop = FALSE
  ) +
  scale_color_manual(
    values = c(
      "Same" = "#911F23",
      "Different" = "#8CB4C5"
    ),
    drop = FALSE,
    name = "LMM vs Bootstrap\nHighest Group"
  ) +
  scale_size_continuous(
    range = c(0.9, 2.6),
    name = "LMM Std.\nEMM Range"
  ) +
  guides(
    color = guide_legend(order = 1, ncol = 1, override.aes = list(size = 1.8)),
    size = guide_legend(order = 2)
  ) +
  labs(
    x = "LMM highest\ngroup",
    y = NULL
  ) +
  theme_bw(base_size = combined_main_font_size) +
  theme(
    plot.title = element_text(size = combined_main_font_size, face = "bold", color = "black"),
    axis.title = element_text(size = combined_main_font_size, color = "black"),
    axis.text.x = element_text(size = combined_axis_font_size, color = "black", angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.background = element_blank(),
    strip.background.y = element_blank(),
    strip.text = element_blank(),
    strip.text.y = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = combined_legend_font_size, color = "black"),
    legend.text = element_text(size = combined_legend_font_size, color = "black"),
    legend.key.size = unit(0.28, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    plot.margin = margin(2, 2, 2, 2),
    panel.spacing.y = unit(0.45, "mm")
  )

p <- (
  patchwork::wrap_plots(
    p_heatmap,
    p_bootstrap,
    p_lmm,
    nrow = 1,
    widths = c(1.45, 0.9, 0.9)
  ) +
    patchwork::plot_annotation(tag_levels = "A")
) &
  theme(
    plot.tag = element_text(size = combined_main_font_size, face = "bold", color = "black"),
    plot.margin = margin(2, 2, 2, 2)
  )

print(p)

## 3. 保存图形 ####
ggsave(
  filename = file.path(figure_output_path, "图8.1_图1图7图8横向拼接图.pdf"),
  plot = p,
  width = 14.65,
  height = 10,
  units = "cm",
  device = cairo_pdf
)


# 9图：生物轴方向汇总热图 ####

## 1. 读取并处理绘图数据 ####
direction_effect_size_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "方向效应量汇总"
)

feature_axis_levels <- direction_effect_size_df %>%
  dplyr::filter(main_display == TRUE) %>%
  dplyr::arrange(axis_order) %>%
  dplyr::pull(feature_axis) %>%
  unique()

plot_data <- direction_effect_size_df %>%
  dplyr::filter(
    main_display == TRUE,
    direction_stability == "Stable",
    effect_size_class %in% c("Moderate", "Large"),
    !is.na(most_often_highest_group)
  ) %>%
  dplyr::count(feature_axis, most_often_highest_group, name = "stable_feature_count") %>%
  tidyr::complete(
    feature_axis = feature_axis_levels,
    most_often_highest_group = pattern_group_order,
    fill = list(stable_feature_count = 0)
  ) %>%
  dplyr::mutate(
    feature_axis = factor(feature_axis, levels = rev(feature_axis_levels)),
    most_often_highest_group = factor(most_often_highest_group, levels = pattern_group_order)
  )

## 2. 绘图 ####
p <- ggplot(plot_data, aes(x = most_often_highest_group, y = feature_axis, fill = stable_feature_count)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(aes(label = stable_feature_count), size = 3.5, color = "black") +
  scale_fill_gradient(low = "#F7F7F7", high = "#D94A4A", name = "Stable Features") +
  labs(
    title = "Stable Direction Patterns by Immune Feature Axis",
    x = "Most Often Highest Cancer Group",
    y = "Feature Axis"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold", color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 14, color = "black"),
    legend.text = element_text(size = 12, color = "black"),
    panel.grid = element_blank()
  )

print(p)

## 3. 保存图形 ####
ggsave(
  filename = file.path(figure_output_path, "图9_生物轴方向汇总热图.pdf"),
  plot = p,
  width = 10,
  height = 5.8,
  device = cairo_pdf
)


# 10图：最高最低组标准化效应量森林图 ####

## 1. 读取并处理绘图数据 ####
direction_effect_size_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "方向效应量汇总"
)

plot_data <- direction_effect_size_df %>%
  dplyr::filter(
    main_display == TRUE,
    !is.na(standardized_median_range),
    !is.na(most_often_highest_group)
  ) %>%
  dplyr::arrange(axis_order, feature_order) %>%
  dplyr::mutate(
    feature_label = factor(feature_label, levels = rev(unique(feature_label))),
    most_often_highest_group = factor(most_often_highest_group, levels = pattern_group_order),
    direction_stability = factor(direction_stability, levels = c("Stable", "Trend", "Unstable", "Not evaluated")),
    comparison_label = paste0(observed_highest_group, " vs ", observed_lowest_group)
  )

## 2. 绘图 ####
p <- ggplot(plot_data, aes(x = standardized_median_range, y = feature_label)) +
  geom_vline(xintercept = c(0.50, 0.80), linetype = "dashed", color = "grey65", linewidth = 0.35) +
  geom_segment(aes(x = 0, xend = standardized_median_range, yend = feature_label), color = "grey70", linewidth = 0.45) +
  geom_point(aes(color = most_often_highest_group, shape = direction_stability), size = 3.2, alpha = 0.9) +
  geom_text(aes(label = comparison_label), hjust = -0.05, size = 3.1, color = "black") +
  scale_color_manual(values = pattern_group_color_map, drop = FALSE, name = "Highest Group") +
  scale_shape_manual(
    values = c("Stable" = 16, "Trend" = 17, "Unstable" = 1, "Not evaluated" = 4),
    drop = FALSE,
    name = "Direction Stability"
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.35))) +
  labs(
    title = "Highest-Lowest Standardized Effect Size",
    x = "Standardized Median Range",
    y = "Feature"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold", color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 14, color = "black"),
    legend.text = element_text(size = 12, color = "black")
  )

print(p)

## 3. 保存图形 ####
ggsave(
  filename = file.path(figure_output_path, "图10_最高最低组标准化效应量森林图.pdf"),
  plot = p,
  width = 13,
  height = 10,
  device = cairo_pdf
)


# 11图：混合模型与Bootstrap方向一致性Tile图 ####

## 1. 读取并处理绘图数据 ####
direction_effect_size_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "方向效应量汇总"
)

plot_data <- direction_effect_size_df %>%
  dplyr::filter(main_display == TRUE) %>%
  dplyr::arrange(axis_order, feature_order) %>%
  dplyr::mutate(
    feature_label = factor(feature_label, levels = rev(unique(feature_label))),
    agreement_panel = "Mixed model vs bootstrap",
    mixed_model_direction_agreement = factor(
      mixed_model_direction_agreement,
      levels = c("Mixed model supported", "Mixed model discordant", "Mixed model unavailable")
    ),
    direction_label = dplyr::case_when(
      is.na(most_often_highest_group) & is.na(mixed_highest_group) ~ "NA / NA",
      is.na(mixed_highest_group) ~ paste0(most_often_highest_group, " / NA"),
      is.na(most_often_highest_group) ~ paste0("NA / ", mixed_highest_group),
      TRUE ~ paste0(most_often_highest_group, " / ", mixed_highest_group)
    )
  )

## 2. 绘图 ####
p <- ggplot(plot_data, aes(x = agreement_panel, y = feature_label, fill = mixed_model_direction_agreement)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(aes(label = direction_label), size = 2.8, color = "black") +
  scale_fill_manual(
    values = c(
      "Mixed model supported" = "#8FD19E",
      "Mixed model discordant" = "#F1A3A3",
      "Mixed model unavailable" = "grey85"
    ),
    drop = FALSE,
    name = "Direction Agreement"
  ) +
  labs(
    title = "Bootstrap and Mixed Model Direction Agreement",
    x = NULL,
    y = "Feature"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold", color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 14, color = "black"),
    legend.text = element_text(size = 12, color = "black"),
    panel.grid = element_blank()
  )

print(p)

## 3. 保存图形 ####
ggsave(
  filename = file.path(figure_output_path, "图11_混合模型与Bootstrap方向一致性Tile图.pdf"),
  plot = p,
  width = 9,
  height = 10,
  device = cairo_pdf
)


# 12图：癌种群原始中位数热图 ####

## 1. 读取并处理绘图数据 ####
pattern_group_summary_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "模式群特征汇总"
)

feature_info_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "特征信息表"
)

feature_order_df <- feature_info_df %>%
  dplyr::filter(main_display == TRUE) %>%
  dplyr::arrange(axis_order, feature_order)

plot_data <- pattern_group_summary_df %>%
  dplyr::filter(feature_name %in% feature_order_df$feature_name) %>%
  dplyr::left_join(
    feature_order_df %>%
      dplyr::select(feature_name, feature_label, feature_axis, axis_order, feature_order),
    by = "feature_name",
    suffix = c("", ".feature_info")
  ) %>%
  dplyr::mutate(
    feature_label = dplyr::coalesce(feature_label.feature_info, feature_label),
    feature_axis = dplyr::coalesce(feature_axis.feature_info, feature_axis),
    pattern_group = factor(pattern_group, levels = pattern_group_order),
    feature_axis = factor(feature_axis, levels = unique(feature_order_df$feature_axis)),
    feature_label = factor(
      feature_label,
      levels = rev(feature_order_df$feature_label[match(feature_order_df$feature_name, feature_order_df$feature_name)])
    ),
    median_score_label = ifelse(is.na(median_score), "NA", sprintf("%.2f", median_score))
  )

## 2. 绘图 ####
p <- ggplot(plot_data, aes(x = pattern_group, y = feature_label, fill = median_score)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(aes(label = median_score_label), size = 2.4, color = "black") +
  facet_grid(feature_axis ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient(low = "#F7F7F7", high = "#D94A4A", na.value = "grey90", name = "Median Score") +
  labs(
    title = "Raw Median Immune Feature Scores across Cancer Groups",
    x = "Cancer Group",
    y = "Feature"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold", color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 12, color = "black"),
    strip.text.y = element_text(size = 12, face = "bold", color = "black", angle = 0),
    legend.title = element_text(size = 14, color = "black"),
    legend.text = element_text(size = 12, color = "black"),
    panel.grid = element_blank()
  )

print(p)

## 3. 保存图形 ####
ggsave(
  filename = file.path(figure_output_path, "图12_癌种群原始中位数热图.pdf"),
  plot = p,
  width = 11,
  height = 12,
  device = cairo_pdf
)


# 13图：两两方向一致性热图 ####

## 1. 读取并处理绘图数据 ####
pairwise_direction_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "抽样两两方向"
)

feature_info_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "特征信息表"
)

feature_order_df <- feature_info_df %>%
  dplyr::filter(main_display == TRUE) %>%
  dplyr::arrange(axis_order, feature_order)

plot_data <- pairwise_direction_df %>%
  dplyr::filter(feature_name %in% feature_order_df$feature_name) %>%
  dplyr::left_join(
    feature_order_df %>%
      dplyr::select(feature_name, feature_label, feature_axis, axis_order, feature_order),
    by = "feature_name",
    suffix = c("", ".feature_info")
  ) %>%
  dplyr::mutate(
    feature_label = dplyr::coalesce(feature_label.feature_info, feature_label),
    feature_axis = dplyr::coalesce(feature_axis.feature_info, feature_axis),
    pairwise_comparison = paste0(group2, " vs ", group1),
    feature_label = factor(
      feature_label,
      levels = rev(feature_order_df$feature_label[match(feature_order_df$feature_name, feature_order_df$feature_name)])
    ),
    direction_label = ifelse(
      is.na(bootstrap_median_diff_group2_minus_group1),
      "NA",
      sprintf("%+.2f", bootstrap_median_diff_group2_minus_group1)
    )
  )

comparison_order <- pairwise_direction_df %>%
  dplyr::distinct(group1, group2) %>%
  dplyr::mutate(pairwise_comparison = paste0(group2, " vs ", group1)) %>%
  dplyr::pull(pairwise_comparison)

plot_data <- plot_data %>%
  dplyr::mutate(pairwise_comparison = factor(pairwise_comparison, levels = comparison_order))

## 2. 绘图 ####
p <- ggplot(plot_data, aes(x = pairwise_comparison, y = feature_label, fill = direction_consistency)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(aes(label = direction_label), size = 2.3, color = "black") +
  scale_fill_gradient(low = "#F7F7F7", high = "#D94A4A", limits = c(0.5, 1), name = "Direction Consistency") +
  labs(
    title = "Pairwise Bootstrap Direction Consistency",
    x = "Pairwise Comparison",
    y = "Feature"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold", color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 10, color = "black", angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 10, color = "black"),
    legend.title = element_text(size = 14, color = "black"),
    legend.text = element_text(size = 12, color = "black"),
    panel.grid = element_blank()
  )

print(p)

## 3. 保存图形 ####
ggsave(
  filename = file.path(figure_output_path, "图13_两两方向一致性热图.pdf"),
  plot = p,
  width = 14,
  height = 12,
  device = cairo_pdf
)


# 14图：癌种级免疫Feature Heatmap ####

## 1. 读取并处理绘图数据 ####
primary_organ_score_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "原发部位级免疫特征"
)

feature_info_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "特征信息表"
)

feature_order_df <- feature_info_df %>%
  dplyr::filter(main_display == TRUE) %>%
  dplyr::arrange(axis_order, feature_order)

available_features <- intersect(feature_order_df$feature_name, colnames(primary_organ_score_df))

available_feature_order_df <- feature_order_df %>%
  dplyr::filter(feature_name %in% available_features)

primary_organ_order_df <- primary_organ_score_df %>%
  dplyr::mutate(pattern_group = factor(pattern_group, levels = pattern_group_order)) %>%
  dplyr::arrange(pattern_group, primary_organ) %>%
  dplyr::distinct(pattern_group, primary_organ, .keep_all = TRUE) %>%
  dplyr::mutate(
    primary_organ_label = paste0(
      gsub("_", " ", primary_organ),
      "\n(N=",
      n_tumor_samples,
      ")"
    )
  )

plot_data <- primary_organ_score_df %>%
  dplyr::select(pattern_group, primary_organ, n_tumor_samples, dplyr::any_of(available_features)) %>%
  tidyr::pivot_longer(
    cols = dplyr::any_of(available_features),
    names_to = "feature_name",
    values_to = "score"
  ) %>%
  dplyr::left_join(
    available_feature_order_df %>%
      dplyr::select(feature_name, feature_label, feature_axis, axis_order, feature_order),
    by = "feature_name"
  ) %>%
  dplyr::group_by(feature_name) %>%
  dplyr::mutate(
    feature_mean = mean(score, na.rm = TRUE),
    feature_sd = stats::sd(score, na.rm = TRUE),
    z_score = dplyr::case_when(
      is.na(feature_sd) | feature_sd == 0 ~ NA_real_,
      TRUE ~ (score - feature_mean) / feature_sd
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(-feature_mean, -feature_sd) %>%
  dplyr::left_join(
    primary_organ_order_df %>%
      dplyr::select(primary_organ, primary_organ_label),
    by = "primary_organ"
  ) %>%
  dplyr::mutate(
    pattern_group = factor(pattern_group, levels = pattern_group_order),
    primary_organ = factor(primary_organ, levels = primary_organ_order_df$primary_organ),
    primary_organ_label = factor(primary_organ_label, levels = primary_organ_order_df$primary_organ_label),
    feature_axis = factor(feature_axis, levels = unique(available_feature_order_df$feature_axis)),
    feature_label = factor(feature_label, levels = rev(available_feature_order_df$feature_label))
  )

## 2. 绘图 ####
p <- ggplot(plot_data, aes(x = primary_organ_label, y = feature_label, fill = z_score)) +
  geom_tile(color = "white", linewidth = 0.2) +
  facet_grid(
    feature_axis ~ pattern_group,
    scales = "free",
    space = "free"
  ) +
  scale_x_discrete(expand = expansion(add = 0)) +
  scale_y_discrete(expand = expansion(add = 0)) +
  scale_fill_gradient2(
    low = heatmap_color_low,
    mid = heatmap_color_mid,
    high = heatmap_color_high,
    midpoint = 0,
    na.value = "grey90",
    name = "Row z-score"
  ) +
  guides(
    fill = guide_colorbar(
      direction = "horizontal",
      title.position = "left",
      title.vjust = 1,
      barwidth = unit(5, "cm"),
      barheight = unit(0.4, "cm")
    )
  ) +
  labs(
    title = "Immune Feature Patterns across Cancer Types",
    x = "Cancer Type",
    y = "Feature"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold", color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 12, color = "black", angle = 60, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 12, color = "black"),
    strip.text.x = element_text(size = 12, face = "bold", color = "black"),
    strip.text.y = element_text(size = 12, face = "bold", color = "black", angle = 0),
    legend.position = "top",
    legend.justification = "center",
    legend.title = element_text(size = 14, color = "black"),
    legend.text = element_text(size = 12, color = "black"),
    panel.grid = element_blank(),
    panel.spacing.x = unit(0.25, "lines"),
    panel.spacing.y = unit(0.15, "lines")
  )

print(p)

## 3. 保存图形 ####
ggsave(
  filename = file.path(figure_output_path, "图14_癌种级免疫Feature热图.pdf"),
  plot = p,
  width = 18,
  height = 12,
  device = cairo_pdf
)


# 15图：LMM癌种随机截距Forest图 ####

## 1. 读取并处理绘图数据 ####
lmm_blup_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "LMM癌种随机截距"
)

feature_info_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "特征信息表"
)

feature_order_df <- feature_info_df %>%
  dplyr::filter(main_display == TRUE) %>%
  dplyr::arrange(axis_order, feature_order)

primary_organ_order_df <- lmm_blup_df %>%
  dplyr::group_by(primary_organ, pattern_group) %>%
  dplyr::summarise(n_samples_in_organ = max(n_samples_in_organ, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(pattern_group = factor(pattern_group, levels = pattern_group_order)) %>%
  dplyr::arrange(pattern_group, primary_organ) %>%
  dplyr::mutate(
    primary_organ_label = paste0(
      gsub("_", " ", primary_organ),
      " (N=",
      round(n_samples_in_organ),
      ")"
    )
  )

plot_data <- lmm_blup_df %>%
  dplyr::filter(feature_name %in% feature_order_df$feature_name) %>%
  dplyr::left_join(
    primary_organ_order_df %>% dplyr::select(primary_organ, primary_organ_label),
    by = "primary_organ"
  ) %>%
  dplyr::mutate(
    pattern_group = factor(pattern_group, levels = pattern_group_order),
    primary_organ_label = factor(primary_organ_label, levels = rev(primary_organ_order_df$primary_organ_label)),
    feature_label = factor(feature_label, levels = feature_order_df$feature_label)
  )

## 2. 绘图 ####
p <- ggplot(plot_data, aes(x = blup_estimate, y = primary_organ_label, color = pattern_group)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.5) +
  geom_errorbarh(
    aes(xmin = blup_ci_lower, xmax = blup_ci_upper),
    height = 0,
    linewidth = 0.6
  ) +
  geom_point(size = 1.8) +
  facet_wrap(~ feature_label, ncol = 4, scales = "free_x") +
  scale_color_manual(values = pattern_group_color_map, name = "Pattern group") +
  labs(
    title = "Cancer-type random intercepts (BLUP) from main LMM",
    subtitle = "feature_z ~ pattern_group + (1 | primary_organ); each point = primary_organ random intercept (z-score units), error bars = 95% CI",
    x = "Random intercept (BLUP, z-score units)",
    y = "Cancer Type"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold", color = "black"),
    plot.subtitle = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    strip.text = element_text(size = 12, face = "bold", color = "black"),
    legend.position = "top",
    legend.title = element_text(size = 14, color = "black"),
    legend.text = element_text(size = 12, color = "black"),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.5, "lines")
  )

print(p)

## 3. 保存图形 ####
ggsave(
  filename = file.path(figure_output_path, "图15_LMM癌种随机截距Forest图.pdf"),
  plot = p,
  width = 22,
  height = 30,
  device = cairo_pdf
)


# 16图：LMM癌种群固定效应Forest图_相对Ocular ####

## 1. 读取并处理绘图数据 ####
mixed_pairwise_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "混合模型两两比较"
)

feature_info_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "特征信息表"
)

feature_order_df <- feature_info_df %>%
  dplyr::filter(main_display == TRUE) %>%
  dplyr::arrange(axis_order, feature_order)

reference_group <- "Ocular"
other_groups <- pattern_group_order[pattern_group_order != reference_group]

plot_data <- mixed_pairwise_df %>%
  dplyr::filter(group1 == reference_group, mixed_pairwise_status == "success") %>%
  dplyr::filter(feature_name %in% feature_order_df$feature_name) %>%
  dplyr::mutate(
    estimate = mixed_emmean_diff_group2_minus_group1,
    ci_lower = estimate - 1.96 * mixed_pairwise_se,
    ci_upper = estimate + 1.96 * mixed_pairwise_se,
    group_label = factor(group2, levels = rev(other_groups)),
    feature_label = factor(feature_label, levels = feature_order_df$feature_label),
    significance_label = dplyr::case_when(
      is.na(mixed_pairwise_fdr) ~ "",
      mixed_pairwise_fdr < 0.001 ~ "***",
      mixed_pairwise_fdr < 0.01 ~ "**",
      mixed_pairwise_fdr < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

## 2. 绘图 ####
p <- ggplot(plot_data, aes(x = estimate, y = group_label, color = group2)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.5) +
  geom_errorbarh(
    aes(xmin = ci_lower, xmax = ci_upper),
    height = 0,
    linewidth = 0.6
  ) +
  geom_point(size = 2.5) +
  geom_text(
    aes(label = significance_label),
    color = "black",
    vjust = -0.6,
    hjust = 0.5,
    size = 4
  ) +
  facet_wrap(~ feature_label, ncol = 4, scales = "free_x") +
  scale_color_manual(values = pattern_group_color_map[other_groups], guide = "none") +
  labs(
    title = paste0("Pattern group fixed effects vs reference (", reference_group, ") from main LMM"),
    subtitle = paste0(
      "Estimate = mean(group) - mean(", reference_group,
      ") in feature z-score units; error bar = 95% CI; ",
      "FDR-adjusted significance: * < 0.05, ** < 0.01, *** < 0.001"
    ),
    x = paste0("Effect vs ", reference_group, " (z-score units)"),
    y = "Pattern group"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold", color = "black"),
    plot.subtitle = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    strip.text = element_text(size = 12, face = "bold", color = "black"),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.5, "lines")
  )

print(p)

## 3. 保存图形 ####
ggsave(
  filename = file.path(figure_output_path, "图16_LMM癌种群固定效应Forest图_相对Ocular.pdf"),
  plot = p,
  width = 18,
  height = 14,
  device = cairo_pdf
)


# 17图：加权LMM最高癌种群与Bootstrap一致性点图 ####

## 1. 读取并处理绘图数据 ####
direction_effect_size_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "方向效应量汇总_加权"
)

plot_data <- direction_effect_size_df %>%
  dplyr::filter(
    main_display == TRUE,
    mixed_model_status == "success",
    !is.na(mixed_highest_group)
  ) %>%
  dplyr::arrange(axis_order, feature_order) %>%
  dplyr::mutate(
    feature_label = factor(feature_label, levels = rev(unique(feature_label))),
    mixed_highest_group = factor(mixed_highest_group, levels = pattern_group_order),
    mixed_model_direction_agreement = factor(
      mixed_model_direction_agreement,
      levels = c("Mixed model supported", "Mixed model discordant"),
      labels = c("Same", "Different")
    )
  )

## 2. 绘图 ####
p <- ggplot(plot_data, aes(x = mixed_highest_group, y = feature_label)) +
  geom_point(aes(size = mixed_standardized_range, color = mixed_model_direction_agreement), alpha = 0.85) +
  scale_color_manual(
    values = c(
      "Same" = "#911F23",
      "Different" = "#8CB4C5"
    ),
    drop = FALSE,
    name = "Weighted LMM vs Bootstrap\nHighest Group"
  ) +
  scale_size_continuous(range = c(2.5, 8), name = "Weighted LMM Standardized\nEMM Range") +
  guides(
    color = guide_legend(order = 1, ncol = 1),
    size = guide_legend(order = 2)
  ) +
  labs(
    x = "Highest Cancer Group\n(Weighted LMM-estimated, weights = 1 / n_organ)",
    y = "Feature"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold", color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 12, color = "black"),
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 12, color = "black"),
    legend.text = element_text(size = 11, color = "black")
  )

print(p)

## 3. 保存图形 ####
ggsave(
  filename = file.path(figure_output_path, "图17_加权LMM最高癌种群与Bootstrap一致性点图.pdf"),
  plot = p,
  width = 6,
  height = 10,
  device = cairo_pdf
)


# 18图：加权混合模型与Bootstrap方向一致性Tile图 ####

## 1. 读取并处理绘图数据 ####
direction_effect_size_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "方向效应量汇总_加权"
)

plot_data <- direction_effect_size_df %>%
  dplyr::filter(main_display == TRUE) %>%
  dplyr::arrange(axis_order, feature_order) %>%
  dplyr::mutate(
    feature_label = factor(feature_label, levels = rev(unique(feature_label))),
    agreement_panel = "Weighted LMM vs bootstrap",
    mixed_model_direction_agreement = factor(
      mixed_model_direction_agreement,
      levels = c("Mixed model supported", "Mixed model discordant", "Mixed model unavailable")
    ),
    direction_label = dplyr::case_when(
      is.na(most_often_highest_group) & is.na(mixed_highest_group) ~ "NA / NA",
      is.na(mixed_highest_group) ~ paste0(most_often_highest_group, " / NA"),
      is.na(most_often_highest_group) ~ paste0("NA / ", mixed_highest_group),
      TRUE ~ paste0(most_often_highest_group, " / ", mixed_highest_group)
    )
  )

## 2. 绘图 ####
p <- ggplot(plot_data, aes(x = agreement_panel, y = feature_label, fill = mixed_model_direction_agreement)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(aes(label = direction_label), size = 2.8, color = "black") +
  scale_fill_manual(
    values = c(
      "Mixed model supported" = "#8FD19E",
      "Mixed model discordant" = "#F1A3A3",
      "Mixed model unavailable" = "grey85"
    ),
    drop = FALSE,
    name = "Direction Agreement"
  ) +
  labs(
    title = "Bootstrap and Weighted Mixed Model Direction Agreement",
    x = NULL,
    y = "Feature"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold", color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 14, color = "black"),
    legend.text = element_text(size = 12, color = "black"),
    panel.grid = element_blank()
  )

print(p)

## 3. 保存图形 ####
ggsave(
  filename = file.path(figure_output_path, "图18_加权混合模型与Bootstrap方向一致性Tile图.pdf"),
  plot = p,
  width = 9,
  height = 10,
  device = cairo_pdf
)


# 19图：加权LMM癌种随机截距Forest图 ####

## 1. 读取并处理绘图数据 ####
lmm_blup_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "LMM癌种随机截距_加权"
)

feature_info_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "特征信息表"
)

feature_order_df <- feature_info_df %>%
  dplyr::filter(main_display == TRUE) %>%
  dplyr::arrange(axis_order, feature_order)

primary_organ_order_df <- lmm_blup_df %>%
  dplyr::group_by(primary_organ, pattern_group) %>%
  dplyr::summarise(n_samples_in_organ = max(n_samples_in_organ, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(pattern_group = factor(pattern_group, levels = pattern_group_order)) %>%
  dplyr::arrange(pattern_group, primary_organ) %>%
  dplyr::mutate(
    primary_organ_label = paste0(
      gsub("_", " ", primary_organ),
      " (N=",
      round(n_samples_in_organ),
      ")"
    )
  )

plot_data <- lmm_blup_df %>%
  dplyr::filter(feature_name %in% feature_order_df$feature_name) %>%
  dplyr::left_join(
    primary_organ_order_df %>% dplyr::select(primary_organ, primary_organ_label),
    by = "primary_organ"
  ) %>%
  dplyr::mutate(
    pattern_group = factor(pattern_group, levels = pattern_group_order),
    primary_organ_label = factor(primary_organ_label, levels = rev(primary_organ_order_df$primary_organ_label)),
    feature_label = factor(feature_label, levels = feature_order_df$feature_label)
  )

## 2. 绘图 ####
p <- ggplot(plot_data, aes(x = blup_estimate, y = primary_organ_label, color = pattern_group)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.5) +
  geom_errorbarh(
    aes(xmin = blup_ci_lower, xmax = blup_ci_upper),
    height = 0,
    linewidth = 0.6
  ) +
  geom_point(size = 1.8) +
  facet_wrap(~ feature_label, ncol = 4, scales = "free_x") +
  scale_color_manual(values = pattern_group_color_map, name = "Pattern group") +
  labs(
    title = "Cancer-type random intercepts (BLUP) from weighted LMM",
    subtitle = "feature_z ~ pattern_group + (1 | primary_organ) with weights = 1 / n_samples_in_organ; error bars = 95% CI",
    x = "Random intercept (BLUP, z-score units)",
    y = "Cancer Type"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold", color = "black"),
    plot.subtitle = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    strip.text = element_text(size = 12, face = "bold", color = "black"),
    legend.position = "top",
    legend.title = element_text(size = 14, color = "black"),
    legend.text = element_text(size = 12, color = "black"),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.5, "lines")
  )

print(p)

## 3. 保存图形 ####
ggsave(
  filename = file.path(figure_output_path, "图19_加权LMM癌种随机截距Forest图.pdf"),
  plot = p,
  width = 22,
  height = 30,
  device = cairo_pdf
)


# 20图：加权LMM癌种群固定效应Forest图_相对Ocular ####

## 1. 读取并处理绘图数据 ####
weighted_pairwise_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "加权混合模型两两"
)

feature_info_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "特征信息表"
)

feature_order_df <- feature_info_df %>%
  dplyr::filter(main_display == TRUE) %>%
  dplyr::arrange(axis_order, feature_order)

reference_group <- "Ocular"
other_groups <- pattern_group_order[pattern_group_order != reference_group]

plot_data <- weighted_pairwise_df %>%
  dplyr::filter(group1 == reference_group, weighted_pairwise_status == "success") %>%
  dplyr::filter(feature_name %in% feature_order_df$feature_name) %>%
  dplyr::mutate(
    estimate = weighted_emmean_diff_group2_minus_group1,
    ci_lower = estimate - 1.96 * weighted_pairwise_se,
    ci_upper = estimate + 1.96 * weighted_pairwise_se,
    group_label = factor(group2, levels = rev(other_groups)),
    feature_label = factor(feature_label, levels = feature_order_df$feature_label),
    significance_label = dplyr::case_when(
      is.na(weighted_pairwise_fdr) ~ "",
      weighted_pairwise_fdr < 0.001 ~ "***",
      weighted_pairwise_fdr < 0.01 ~ "**",
      weighted_pairwise_fdr < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

## 2. 绘图 ####
p <- ggplot(plot_data, aes(x = estimate, y = group_label, color = group2)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.5) +
  geom_errorbarh(
    aes(xmin = ci_lower, xmax = ci_upper),
    height = 0,
    linewidth = 0.6
  ) +
  geom_point(size = 2.5) +
  geom_text(
    aes(label = significance_label),
    color = "black",
    vjust = -0.6,
    hjust = 0.5,
    size = 4
  ) +
  facet_wrap(~ feature_label, ncol = 4, scales = "free_x") +
  scale_color_manual(values = pattern_group_color_map[other_groups], guide = "none") +
  labs(
    title = paste0("Pattern group fixed effects vs reference (", reference_group, ") from weighted LMM"),
    subtitle = paste0(
      "Estimate = mean(group) - mean(", reference_group,
      ") in feature z-score units (weights = 1 / n_organ); error bar = 95% CI; ",
      "FDR-adjusted significance: * < 0.05, ** < 0.01, *** < 0.001"
    ),
    x = paste0("Effect vs ", reference_group, " (z-score units)"),
    y = "Pattern group"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold", color = "black"),
    plot.subtitle = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    strip.text = element_text(size = 12, face = "bold", color = "black"),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.5, "lines")
  )

print(p)

## 3. 保存图形 ####
ggsave(
  filename = file.path(figure_output_path, "图20_加权LMM癌种群固定效应Forest图_相对Ocular.pdf"),
  plot = p,
  width = 18,
  height = 14,
  device = cairo_pdf
)


# 21图：分层抽样LMM最高癌种群与Bootstrap一致性点图 ####

## 1. 读取并处理绘图数据 ####
direction_effect_size_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "方向效应量汇总_分层抽样"
)

plot_data <- direction_effect_size_df %>%
  dplyr::filter(
    main_display == TRUE,
    mixed_model_status == "success",
    !is.na(mixed_highest_group)
  ) %>%
  dplyr::arrange(axis_order, feature_order) %>%
  dplyr::mutate(
    feature_label = factor(feature_label, levels = rev(unique(feature_label))),
    mixed_highest_group = factor(mixed_highest_group, levels = pattern_group_order),
    mixed_model_direction_agreement = factor(
      mixed_model_direction_agreement,
      levels = c("Mixed model supported", "Mixed model discordant"),
      labels = c("Same", "Different")
    )
  )

## 2. 绘图 ####
p <- ggplot(plot_data, aes(x = mixed_highest_group, y = feature_label)) +
  geom_point(aes(size = mixed_standardized_range, color = mixed_model_direction_agreement), alpha = 0.85) +
  scale_color_manual(
    values = c(
      "Same" = "#911F23",
      "Different" = "#8CB4C5"
    ),
    drop = FALSE,
    name = "Subsample LMM vs Bootstrap\nHighest Group"
  ) +
  scale_size_continuous(range = c(2.5, 8), name = "Subsample LMM Median\nStandardized EMM Range") +
  guides(
    color = guide_legend(order = 1, ncol = 1),
    size = guide_legend(order = 2)
  ) +
  labs(
    x = "Most Often Highest Cancer Group\n(Stratified subsampling LMM-estimated)",
    y = "Feature"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold", color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 12, color = "black"),
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 12, color = "black"),
    legend.text = element_text(size = 11, color = "black")
  )

print(p)

## 3. 保存图形 ####
ggsave(
  filename = file.path(figure_output_path, "图21_分层抽样LMM最高癌种群与Bootstrap一致性点图.pdf"),
  plot = p,
  width = 6,
  height = 10,
  device = cairo_pdf
)


# 22图：分层抽样混合模型与Bootstrap方向一致性Tile图 ####

## 1. 读取并处理绘图数据 ####
direction_effect_size_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "方向效应量汇总_分层抽样"
)

plot_data <- direction_effect_size_df %>%
  dplyr::filter(main_display == TRUE) %>%
  dplyr::arrange(axis_order, feature_order) %>%
  dplyr::mutate(
    feature_label = factor(feature_label, levels = rev(unique(feature_label))),
    agreement_panel = "Subsample LMM vs bootstrap",
    mixed_model_direction_agreement = factor(
      mixed_model_direction_agreement,
      levels = c("Mixed model supported", "Mixed model discordant", "Mixed model unavailable")
    ),
    direction_label = dplyr::case_when(
      is.na(most_often_highest_group) & is.na(mixed_highest_group) ~ "NA / NA",
      is.na(mixed_highest_group) ~ paste0(most_often_highest_group, " / NA"),
      is.na(most_often_highest_group) ~ paste0("NA / ", mixed_highest_group),
      TRUE ~ paste0(most_often_highest_group, " / ", mixed_highest_group)
    )
  )

## 2. 绘图 ####
p <- ggplot(plot_data, aes(x = agreement_panel, y = feature_label, fill = mixed_model_direction_agreement)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(aes(label = direction_label), size = 2.8, color = "black") +
  scale_fill_manual(
    values = c(
      "Mixed model supported" = "#8FD19E",
      "Mixed model discordant" = "#F1A3A3",
      "Mixed model unavailable" = "grey85"
    ),
    drop = FALSE,
    name = "Direction Agreement"
  ) +
  labs(
    title = "Bootstrap and Stratified-subsampling LMM Direction Agreement",
    x = NULL,
    y = "Feature"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold", color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 14, color = "black"),
    legend.text = element_text(size = 12, color = "black"),
    panel.grid = element_blank()
  )

print(p)

## 3. 保存图形 ####
ggsave(
  filename = file.path(figure_output_path, "图22_分层抽样混合模型与Bootstrap方向一致性Tile图.pdf"),
  plot = p,
  width = 9,
  height = 10,
  device = cairo_pdf
)


# 23图：分层抽样LMM癌种随机截距Forest图 ####

## 1. 读取并处理绘图数据 ####
lmm_blup_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "LMM癌种随机截距_分层抽样"
)

feature_info_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "特征信息表"
)

feature_order_df <- feature_info_df %>%
  dplyr::filter(main_display == TRUE) %>%
  dplyr::arrange(axis_order, feature_order)

primary_organ_order_df <- lmm_blup_df %>%
  dplyr::group_by(primary_organ, pattern_group) %>%
  dplyr::summarise(n_samples_in_organ = max(n_samples_in_organ, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(pattern_group = factor(pattern_group, levels = pattern_group_order)) %>%
  dplyr::arrange(pattern_group, primary_organ) %>%
  dplyr::mutate(
    primary_organ_label = paste0(
      gsub("_", " ", primary_organ),
      " (N=",
      round(n_samples_in_organ),
      ")"
    )
  )

plot_data <- lmm_blup_df %>%
  dplyr::filter(feature_name %in% feature_order_df$feature_name) %>%
  dplyr::left_join(
    primary_organ_order_df %>% dplyr::select(primary_organ, primary_organ_label),
    by = "primary_organ"
  ) %>%
  dplyr::mutate(
    pattern_group = factor(pattern_group, levels = pattern_group_order),
    primary_organ_label = factor(primary_organ_label, levels = rev(primary_organ_order_df$primary_organ_label)),
    feature_label = factor(feature_label, levels = feature_order_df$feature_label)
  )

## 2. 绘图 ####
p <- ggplot(plot_data, aes(x = blup_estimate, y = primary_organ_label, color = pattern_group)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.5) +
  geom_errorbarh(
    aes(xmin = blup_ci_lower, xmax = blup_ci_upper),
    height = 0,
    linewidth = 0.6
  ) +
  geom_point(size = 1.8) +
  facet_wrap(~ feature_label, ncol = 4, scales = "free_x") +
  scale_color_manual(values = pattern_group_color_map, name = "Pattern group") +
  labs(
    title = "Cancer-type random intercepts (BLUP) from stratified-subsampling LMM",
    subtitle = "Median BLUP across 200 iterations (each iteration subsamples per primary_organ); error bars = 2.5%-97.5% quantile of BLUP distribution",
    x = "Random intercept (median BLUP, z-score units)",
    y = "Cancer Type"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold", color = "black"),
    plot.subtitle = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    strip.text = element_text(size = 12, face = "bold", color = "black"),
    legend.position = "top",
    legend.title = element_text(size = 14, color = "black"),
    legend.text = element_text(size = 12, color = "black"),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.5, "lines")
  )

print(p)

## 3. 保存图形 ####
ggsave(
  filename = file.path(figure_output_path, "图23_分层抽样LMM癌种随机截距Forest图.pdf"),
  plot = p,
  width = 22,
  height = 30,
  device = cairo_pdf
)


# 24图：分层抽样LMM癌种群固定效应Forest图_相对Ocular ####

## 1. 读取并处理绘图数据 ####
subsample_pairwise_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "分层抽样LMM两两"
)

feature_info_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "特征信息表"
)

feature_order_df <- feature_info_df %>%
  dplyr::filter(main_display == TRUE) %>%
  dplyr::arrange(axis_order, feature_order)

reference_group <- "Ocular"
other_groups <- pattern_group_order[pattern_group_order != reference_group]

plot_data <- subsample_pairwise_df %>%
  dplyr::filter(group1 == reference_group) %>%
  dplyr::filter(feature_name %in% feature_order_df$feature_name) %>%
  dplyr::mutate(
    estimate = subsample_pairwise_median_estimate_g2_minus_g1,
    ci_lower = subsample_pairwise_estimate_ci_lower,
    ci_upper = subsample_pairwise_estimate_ci_upper,
    group_label = factor(group2, levels = rev(other_groups)),
    feature_label = factor(feature_label, levels = feature_order_df$feature_label),
    significance_label = dplyr::case_when(
      is.na(subsample_pairwise_median_p_value) ~ "",
      subsample_pairwise_median_p_value < 0.001 ~ "***",
      subsample_pairwise_median_p_value < 0.01 ~ "**",
      subsample_pairwise_median_p_value < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

## 2. 绘图 ####
p <- ggplot(plot_data, aes(x = estimate, y = group_label, color = group2)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.5) +
  geom_errorbarh(
    aes(xmin = ci_lower, xmax = ci_upper),
    height = 0,
    linewidth = 0.6
  ) +
  geom_point(size = 2.5) +
  geom_text(
    aes(label = significance_label),
    color = "black",
    vjust = -0.6,
    hjust = 0.5,
    size = 4
  ) +
  facet_wrap(~ feature_label, ncol = 4, scales = "free_x") +
  scale_color_manual(values = pattern_group_color_map[other_groups], guide = "none") +
  labs(
    title = paste0("Pattern group fixed effects vs reference (", reference_group, ") from stratified-subsampling LMM"),
    subtitle = paste0(
      "Estimate = median(group - ", reference_group,
      ") across 200 subsample iterations; error bar = 2.5%-97.5% quantile; ",
      "significance based on median p-value: * < 0.05, ** < 0.01, *** < 0.001"
    ),
    x = paste0("Effect vs ", reference_group, " (z-score units)"),
    y = "Pattern group"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold", color = "black"),
    plot.subtitle = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    strip.text = element_text(size = 12, face = "bold", color = "black"),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.5, "lines")
  )

print(p)

## 3. 保存图形 ####
ggsave(
  filename = file.path(figure_output_path, "图24_分层抽样LMM癌种群固定效应Forest图_相对Ocular.pdf"),
  plot = p,
  width = 18,
  height = 14,
  device = cairo_pdf
)


# 25图：主LMM癌种群emmean Forest图 ####

## 1. 读取并处理绘图数据 ####
emmean_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "LMM癌种群emmean"
)

feature_info_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "特征信息表"
)

feature_order_df <- feature_info_df %>%
  dplyr::filter(main_display == TRUE) %>%
  dplyr::arrange(axis_order, feature_order)

plot_data <- emmean_df %>%
  dplyr::filter(mixed_emmean_status == "success") %>%
  dplyr::filter(feature_name %in% feature_order_df$feature_name) %>%
  dplyr::mutate(
    estimate = mixed_emmean_estimate,
    ci_lower = mixed_emmean_ci_lower,
    ci_upper = mixed_emmean_ci_upper,
    group_label = factor(pattern_group, levels = rev(pattern_group_order)),
    feature_label = factor(feature_label, levels = feature_order_df$feature_label)
  )

## 2. 绘图 ####
p <- ggplot(plot_data, aes(x = estimate, y = group_label, color = pattern_group)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.5) +
  geom_errorbarh(
    aes(xmin = ci_lower, xmax = ci_upper),
    height = 0,
    linewidth = 0.6
  ) +
  geom_point(size = 2.5) +
  facet_wrap(~ feature_label, ncol = 4, scales = "free_x") +
  scale_color_manual(values = pattern_group_color_map, guide = "none") +
  labs(
    title = "Pattern group emmean (model-adjusted mean) from main LMM",
    subtitle = "Estimate = emmean in feature z-score units; error bar = 95% CI (asymptotic z); 0 = global mean",
    x = "Emmean (z-score units)",
    y = "Pattern group"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold", color = "black"),
    plot.subtitle = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    strip.text = element_text(size = 12, face = "bold", color = "black"),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.5, "lines")
  )

print(p)

## 3. 保存图形 ####
ggsave(
  filename = file.path(figure_output_path, "图25_主LMM癌种群emmean_Forest图.pdf"),
  plot = p,
  width = 18,
  height = 16,
  device = cairo_pdf
)


# 26图：加权LMM癌种群emmean Forest图 ####

## 1. 读取并处理绘图数据 ####
weighted_emmean_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "LMM癌种群emmean_加权"
)

feature_info_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "特征信息表"
)

feature_order_df <- feature_info_df %>%
  dplyr::filter(main_display == TRUE) %>%
  dplyr::arrange(axis_order, feature_order)

plot_data <- weighted_emmean_df %>%
  dplyr::filter(weighted_emmean_status == "success") %>%
  dplyr::filter(feature_name %in% feature_order_df$feature_name) %>%
  dplyr::mutate(
    estimate = weighted_emmean_estimate,
    ci_lower = weighted_emmean_ci_lower,
    ci_upper = weighted_emmean_ci_upper,
    group_label = factor(pattern_group, levels = rev(pattern_group_order)),
    feature_label = factor(feature_label, levels = feature_order_df$feature_label)
  )

## 2. 绘图 ####
p <- ggplot(plot_data, aes(x = estimate, y = group_label, color = pattern_group)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.5) +
  geom_errorbarh(
    aes(xmin = ci_lower, xmax = ci_upper),
    height = 0,
    linewidth = 0.6
  ) +
  geom_point(size = 2.5) +
  facet_wrap(~ feature_label, ncol = 4, scales = "free_x") +
  scale_color_manual(values = pattern_group_color_map, guide = "none") +
  labs(
    title = "Pattern group emmean (model-adjusted mean) from equal-cancer-type-weight LMM",
    subtitle = "Estimate = emmean in feature z-score units; error bar = 95% CI (asymptotic z); weights = 1 / n_organ",
    x = "Emmean (z-score units)",
    y = "Pattern group"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold", color = "black"),
    plot.subtitle = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    strip.text = element_text(size = 12, face = "bold", color = "black"),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.5, "lines")
  )

print(p)

## 3. 保存图形 ####
ggsave(
  filename = file.path(figure_output_path, "图26_加权LMM癌种群emmean_Forest图.pdf"),
  plot = p,
  width = 18,
  height = 16,
  device = cairo_pdf
)


# 27图：分层抽样LMM癌种群emmean Forest图 ####

## 1. 读取并处理绘图数据 ####
subsample_emmean_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "LMM癌种群emmean_分层抽样"
)

feature_info_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "特征信息表"
)

feature_order_df <- feature_info_df %>%
  dplyr::filter(main_display == TRUE) %>%
  dplyr::arrange(axis_order, feature_order)

plot_data <- subsample_emmean_df %>%
  dplyr::filter(feature_name %in% feature_order_df$feature_name) %>%
  dplyr::mutate(
    estimate = subsample_emmean_median_estimate,
    ci_lower = subsample_emmean_ci_lower,
    ci_upper = subsample_emmean_ci_upper,
    group_label = factor(pattern_group, levels = rev(pattern_group_order)),
    feature_label = factor(feature_label, levels = feature_order_df$feature_label)
  )

## 2. 绘图 ####
p <- ggplot(plot_data, aes(x = estimate, y = group_label, color = pattern_group)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.5) +
  geom_errorbarh(
    aes(xmin = ci_lower, xmax = ci_upper),
    height = 0,
    linewidth = 0.6
  ) +
  geom_point(size = 2.5) +
  facet_wrap(~ feature_label, ncol = 4, scales = "free_x") +
  scale_color_manual(values = pattern_group_color_map, guide = "none") +
  labs(
    title = "Pattern group emmean from stratified-subsampling LMM",
    subtitle = "Estimate = median emmean across 200 subsample iterations; error bar = 2.5%-97.5% quantile",
    x = "Emmean (z-score units)",
    y = "Pattern group"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold", color = "black"),
    plot.subtitle = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    strip.text = element_text(size = 12, face = "bold", color = "black"),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.5, "lines")
  )

print(p)

## 3. 保存图形 ####
ggsave(
  filename = file.path(figure_output_path, "图27_分层抽样LMM癌种群emmean_Forest图.pdf"),
  plot = p,
  width = 18,
  height = 16,
  device = cairo_pdf
)


# 28图：主LMM癌种群one-vs-rest Forest图 ####

## 1. 读取并处理绘图数据 ####
one_vs_rest_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "LMM癌种群one_vs_rest"
)

feature_info_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "特征信息表"
)

feature_order_df <- feature_info_df %>%
  dplyr::filter(main_display == TRUE) %>%
  dplyr::arrange(axis_order, feature_order)

plot_data <- one_vs_rest_df %>%
  dplyr::filter(mixed_one_vs_rest_status == "success") %>%
  dplyr::filter(feature_name %in% feature_order_df$feature_name) %>%
  dplyr::mutate(
    estimate = mixed_one_vs_rest_estimate,
    ci_lower = mixed_one_vs_rest_ci_lower,
    ci_upper = mixed_one_vs_rest_ci_upper,
    group_label = factor(pattern_group, levels = rev(pattern_group_order)),
    feature_label = factor(feature_label, levels = feature_order_df$feature_label),
    significance_label = dplyr::case_when(
      is.na(mixed_one_vs_rest_fdr) ~ "",
      mixed_one_vs_rest_fdr < 0.001 ~ "***",
      mixed_one_vs_rest_fdr < 0.01 ~ "**",
      mixed_one_vs_rest_fdr < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

## 2. 绘图 ####
p <- ggplot(plot_data, aes(x = estimate, y = group_label, color = pattern_group)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.5) +
  geom_errorbarh(
    aes(xmin = ci_lower, xmax = ci_upper),
    height = 0,
    linewidth = 0.6
  ) +
  geom_point(size = 2.5) +
  geom_text(
    aes(label = significance_label),
    color = "black",
    vjust = -0.6,
    hjust = 0.5,
    size = 4
  ) +
  facet_wrap(~ feature_label, ncol = 4, scales = "free_x") +
  scale_color_manual(values = pattern_group_color_map, guide = "none") +
  labs(
    title = "Pattern group one-vs-rest contrast from main LMM",
    subtitle = "Estimate = emmean(group) - mean(rest 4 groups); error bar = 95% CI; FDR-adjusted: * < 0.05, ** < 0.01, *** < 0.001",
    x = "Effect vs rest (z-score units)",
    y = "Pattern group"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold", color = "black"),
    plot.subtitle = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    strip.text = element_text(size = 12, face = "bold", color = "black"),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.5, "lines")
  )

print(p)

## 3. 保存图形 ####
ggsave(
  filename = file.path(figure_output_path, "图28_主LMM癌种群one_vs_rest_Forest图.pdf"),
  plot = p,
  width = 18,
  height = 16,
  device = cairo_pdf
)


# 29图：加权LMM癌种群one-vs-rest Forest图 ####

## 1. 读取并处理绘图数据 ####
weighted_one_vs_rest_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "LMM癌种群one_vs_rest_加权"
)

feature_info_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "特征信息表"
)

feature_order_df <- feature_info_df %>%
  dplyr::filter(main_display == TRUE) %>%
  dplyr::arrange(axis_order, feature_order)

plot_data <- weighted_one_vs_rest_df %>%
  dplyr::filter(weighted_one_vs_rest_status == "success") %>%
  dplyr::filter(feature_name %in% feature_order_df$feature_name) %>%
  dplyr::mutate(
    estimate = weighted_one_vs_rest_estimate,
    ci_lower = weighted_one_vs_rest_ci_lower,
    ci_upper = weighted_one_vs_rest_ci_upper,
    group_label = factor(pattern_group, levels = rev(pattern_group_order)),
    feature_label = factor(feature_label, levels = feature_order_df$feature_label),
    significance_label = dplyr::case_when(
      is.na(weighted_one_vs_rest_fdr) ~ "",
      weighted_one_vs_rest_fdr < 0.001 ~ "***",
      weighted_one_vs_rest_fdr < 0.01 ~ "**",
      weighted_one_vs_rest_fdr < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

## 2. 绘图 ####
p <- ggplot(plot_data, aes(x = estimate, y = group_label, color = pattern_group)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.5) +
  geom_errorbarh(
    aes(xmin = ci_lower, xmax = ci_upper),
    height = 0,
    linewidth = 0.6
  ) +
  geom_point(size = 2.5) +
  geom_text(
    aes(label = significance_label),
    color = "black",
    vjust = -0.6,
    hjust = 0.5,
    size = 4
  ) +
  facet_wrap(~ feature_label, ncol = 4, scales = "free_x") +
  scale_color_manual(values = pattern_group_color_map, guide = "none") +
  labs(
    title = "Pattern group one-vs-rest contrast from equal-cancer-type-weight LMM",
    subtitle = "Estimate = emmean(group) - mean(rest 4 groups) in weighted LMM; error bar = 95% CI; FDR-adjusted: * < 0.05, ** < 0.01, *** < 0.001",
    x = "Effect vs rest (z-score units)",
    y = "Pattern group"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold", color = "black"),
    plot.subtitle = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    strip.text = element_text(size = 12, face = "bold", color = "black"),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.5, "lines")
  )

print(p)

## 3. 保存图形 ####
ggsave(
  filename = file.path(figure_output_path, "图29_加权LMM癌种群one_vs_rest_Forest图.pdf"),
  plot = p,
  width = 18,
  height = 16,
  device = cairo_pdf
)


# 30图：分层抽样LMM癌种群one-vs-rest Forest图 ####

## 1. 读取并处理绘图数据 ####
subsample_one_vs_rest_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "LMM癌种群one_vs_rest_分层抽样"
)

feature_info_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "特征信息表"
)

feature_order_df <- feature_info_df %>%
  dplyr::filter(main_display == TRUE) %>%
  dplyr::arrange(axis_order, feature_order)

plot_data <- subsample_one_vs_rest_df %>%
  dplyr::filter(feature_name %in% feature_order_df$feature_name) %>%
  dplyr::mutate(
    estimate = subsample_one_vs_rest_median_estimate,
    ci_lower = subsample_one_vs_rest_ci_lower,
    ci_upper = subsample_one_vs_rest_ci_upper,
    group_label = factor(pattern_group, levels = rev(pattern_group_order)),
    feature_label = factor(feature_label, levels = feature_order_df$feature_label),
    significance_label = dplyr::case_when(
      is.na(subsample_one_vs_rest_median_p_value) ~ "",
      subsample_one_vs_rest_median_p_value < 0.001 ~ "***",
      subsample_one_vs_rest_median_p_value < 0.01 ~ "**",
      subsample_one_vs_rest_median_p_value < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

## 2. 绘图 ####
p <- ggplot(plot_data, aes(x = estimate, y = group_label, color = pattern_group)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.5) +
  geom_errorbarh(
    aes(xmin = ci_lower, xmax = ci_upper),
    height = 0,
    linewidth = 0.6
  ) +
  geom_point(size = 2.5) +
  geom_text(
    aes(label = significance_label),
    color = "black",
    vjust = -0.6,
    hjust = 0.5,
    size = 4
  ) +
  facet_wrap(~ feature_label, ncol = 4, scales = "free_x") +
  scale_color_manual(values = pattern_group_color_map, guide = "none") +
  labs(
    title = "Pattern group one-vs-rest contrast from stratified-subsampling LMM",
    subtitle = "Estimate = median(group - rest) across 200 iterations; error bar = 2.5%-97.5% quantile; significance based on median p: * < 0.05, ** < 0.01, *** < 0.001",
    x = "Effect vs rest (z-score units)",
    y = "Pattern group"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold", color = "black"),
    plot.subtitle = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    strip.text = element_text(size = 12, face = "bold", color = "black"),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.5, "lines")
  )

print(p)

## 3. 保存图形 ####
ggsave(
  filename = file.path(figure_output_path, "图30_分层抽样LMM癌种群one_vs_rest_Forest图.pdf"),
  plot = p,
  width = 18,
  height = 16,
  device = cairo_pdf
)


# 31图：CEM匹配LMM癌种群固定效应Forest图_相对Ocular ####

## 1. 读取并处理绘图数据 ####
cem_pairwise_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "CEM匹配LMM两两"
)

feature_info_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "特征信息表"
)

feature_order_df <- feature_info_df %>%
  dplyr::filter(main_display == TRUE) %>%
  dplyr::arrange(axis_order, feature_order)

reference_group <- "Ocular"
other_groups <- pattern_group_order[pattern_group_order != reference_group]

plot_data <- cem_pairwise_df %>%
  dplyr::filter(group1 == reference_group, cem_mixed_pairwise_status == "success") %>%
  dplyr::filter(feature_name %in% feature_order_df$feature_name) %>%
  dplyr::mutate(
    estimate = cem_mixed_emmean_diff_group2_minus_group1,
    ci_lower = estimate - 1.96 * cem_mixed_pairwise_se,
    ci_upper = estimate + 1.96 * cem_mixed_pairwise_se,
    group_label = factor(group2, levels = rev(other_groups)),
    feature_label = factor(feature_label, levels = feature_order_df$feature_label),
    significance_label = dplyr::case_when(
      is.na(cem_mixed_pairwise_fdr) ~ "",
      cem_mixed_pairwise_fdr < 0.001 ~ "***",
      cem_mixed_pairwise_fdr < 0.01 ~ "**",
      cem_mixed_pairwise_fdr < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

## 2. 绘图 ####
p <- ggplot(plot_data, aes(x = estimate, y = group_label, color = group2)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.5) +
  geom_errorbarh(
    aes(xmin = ci_lower, xmax = ci_upper),
    height = 0,
    linewidth = 0.6
  ) +
  geom_point(size = 2.5) +
  geom_text(
    aes(label = significance_label),
    color = "black",
    vjust = -0.6,
    hjust = 0.5,
    size = 4
  ) +
  facet_wrap(~ feature_label, ncol = 4, scales = "free_x") +
  scale_color_manual(values = pattern_group_color_map[other_groups], guide = "none") +
  labs(
    title = paste0("Pattern group fixed effects vs reference (", reference_group, ") from CEM-matched LMM"),
    subtitle = paste0(
      "Per-organ coarsened exact matching on age + tumor purity, N = 35 per primary_organ; ",
      "error bar = 95% CI; FDR-adjusted significance: * < 0.05, ** < 0.01, *** < 0.001"
    ),
    x = paste0("Effect vs ", reference_group, " (z-score units)"),
    y = "Pattern group"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold", color = "black"),
    plot.subtitle = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    strip.text = element_text(size = 12, face = "bold", color = "black"),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.5, "lines")
  )

print(p)

## 3. 保存图形 ####
ggsave(
  filename = file.path(figure_output_path, "图31_CEM匹配LMM癌种群固定效应Forest图_相对Ocular.pdf"),
  plot = p,
  width = 18,
  height = 14,
  device = cairo_pdf
)


# 32图：CEM匹配LMM癌种群emmean_Forest图 ####

## 1. 读取并处理绘图数据 ####
cem_emmean_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "CEM匹配LMM癌种群emmean"
)

feature_info_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "特征信息表"
)

feature_order_df <- feature_info_df %>%
  dplyr::filter(main_display == TRUE) %>%
  dplyr::arrange(axis_order, feature_order)

plot_data <- cem_emmean_df %>%
  dplyr::filter(cem_mixed_emmean_status == "success") %>%
  dplyr::filter(feature_name %in% feature_order_df$feature_name) %>%
  dplyr::mutate(
    estimate = cem_mixed_emmean_estimate,
    ci_lower = cem_mixed_emmean_ci_lower,
    ci_upper = cem_mixed_emmean_ci_upper,
    group_label = factor(pattern_group, levels = rev(pattern_group_order)),
    feature_label = factor(feature_label, levels = feature_order_df$feature_label)
  )

## 2. 绘图 ####
p <- ggplot(plot_data, aes(x = estimate, y = group_label, color = pattern_group)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.5) +
  geom_errorbarh(
    aes(xmin = ci_lower, xmax = ci_upper),
    height = 0,
    linewidth = 0.6
  ) +
  geom_point(size = 2.5) +
  facet_wrap(~ feature_label, ncol = 4, scales = "free_x") +
  scale_color_manual(values = pattern_group_color_map, guide = "none") +
  labs(
    title = "Pattern group estimated marginal means from CEM-matched LMM",
    subtitle = "Per-organ coarsened exact matching on age + tumor purity, N = 35 per primary_organ; error bar = 95% CI (z-score units)",
    x = "Adjusted mean (z-score units)",
    y = "Pattern group"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold", color = "black"),
    plot.subtitle = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    strip.text = element_text(size = 12, face = "bold", color = "black"),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.5, "lines")
  )

print(p)

## 3. 保存图形 ####
ggsave(
  filename = file.path(figure_output_path, "图32_CEM匹配LMM癌种群emmean_Forest图.pdf"),
  plot = p,
  width = 18,
  height = 16,
  device = cairo_pdf
)


# 33图：CEM匹配LMM癌种群one_vs_rest_Forest图 ####

## 1. 读取并处理绘图数据 ####
cem_one_vs_rest_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "CEM匹配LMM癌种群one_vs_rest"
)

feature_info_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "特征信息表"
)

feature_order_df <- feature_info_df %>%
  dplyr::filter(main_display == TRUE) %>%
  dplyr::arrange(axis_order, feature_order)

plot_data <- cem_one_vs_rest_df %>%
  dplyr::filter(cem_mixed_one_vs_rest_status == "success") %>%
  dplyr::filter(feature_name %in% feature_order_df$feature_name) %>%
  dplyr::mutate(
    estimate = cem_mixed_one_vs_rest_estimate,
    ci_lower = cem_mixed_one_vs_rest_ci_lower,
    ci_upper = cem_mixed_one_vs_rest_ci_upper,
    group_label = factor(pattern_group, levels = rev(pattern_group_order)),
    feature_label = factor(feature_label, levels = feature_order_df$feature_label),
    significance_label = dplyr::case_when(
      is.na(cem_mixed_one_vs_rest_fdr) ~ "",
      cem_mixed_one_vs_rest_fdr < 0.001 ~ "***",
      cem_mixed_one_vs_rest_fdr < 0.01 ~ "**",
      cem_mixed_one_vs_rest_fdr < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

## 2. 绘图 ####
p <- ggplot(plot_data, aes(x = estimate, y = group_label, color = pattern_group)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.5) +
  geom_errorbarh(
    aes(xmin = ci_lower, xmax = ci_upper),
    height = 0,
    linewidth = 0.6
  ) +
  geom_point(size = 2.5) +
  geom_text(
    aes(label = significance_label),
    color = "black",
    vjust = -0.6,
    hjust = 0.5,
    size = 4
  ) +
  facet_wrap(~ feature_label, ncol = 4, scales = "free_x") +
  scale_color_manual(values = pattern_group_color_map, guide = "none") +
  labs(
    title = "Pattern group one-vs-rest contrast from CEM-matched LMM",
    subtitle = "Per-organ coarsened exact matching on age + tumor purity, N = 35 per primary_organ; error bar = 95% CI; FDR-adjusted: * < 0.05, ** < 0.01, *** < 0.001",
    x = "Effect vs rest (z-score units)",
    y = "Pattern group"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold", color = "black"),
    plot.subtitle = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    strip.text = element_text(size = 12, face = "bold", color = "black"),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.5, "lines")
  )

print(p)

## 3. 保存图形 ####
ggsave(
  filename = file.path(figure_output_path, "图33_CEM匹配LMM癌种群one_vs_rest_Forest图.pdf"),
  plot = p,
  width = 18,
  height = 16,
  device = cairo_pdf
)


# 34图：CEM匹配前后年龄均值对比图 ####

## 1. 读取并处理绘图数据 ####
cem_balance_df <- readxl::read_xlsx(
  path = immune_result_xlsx,
  sheet = "CEM匹配协变量平衡"
)

plot_data <- cem_balance_df %>%
  dplyr::filter(primary_organ != "ALL", covariate == "age_at_index") %>%
  dplyr::select(primary_organ, mean_before, mean_after) %>%
  tidyr::pivot_longer(
    cols = c(mean_before, mean_after),
    names_to = "stage",
    values_to = "mean_value"
  ) %>%
  dplyr::mutate(
    stage = dplyr::recode(stage, mean_before = "Before matching", mean_after = "After matching"),
    stage = factor(stage, levels = c("Before matching", "After matching"))
  )

## 2. 绘图 ####
p <- ggplot(plot_data, aes(x = stage, y = mean_value, group = primary_organ)) +
  geom_line(color = "grey70", linewidth = 0.4, alpha = 0.8) +
  geom_point(aes(color = stage), size = 2) +
  scale_color_manual(values = c("Before matching" = "#4C72B0", "After matching" = "#C44E52"), guide = "none") +
  labs(
    title = "Per-organ age-at-index means before vs after CEM matching",
    subtitle = "Each line = one primary_organ; lines converging toward a common level indicate improved cross-organ age balance",
    x = NULL,
    y = "Per-organ mean age at index (years)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold", color = "black"),
    plot.subtitle = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    strip.text = element_text(size = 12, face = "bold", color = "black"),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.8, "lines")
  )

print(p)

## 3. 保存图形 ####
ggsave(
  filename = file.path(figure_output_path, "图34_CEM匹配前后年龄均值对比图.pdf"),
  plot = p,
  width = 8,
  height = 6,
  device = cairo_pdf
)


message("========== 7.4.1 肿瘤免疫微环境绘图全部完成 ==========")
