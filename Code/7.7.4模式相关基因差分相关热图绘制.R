# 7.7.4 模式相关基因差分相关热图绘制 ####
# 用途: 从 7.7 模式相关基因差分相关识别结果中筛选五个模式的关键基因,
#       绘制五模式关键基因 Δρ dot-heatmap。 ####
# 输入: Results/7.7模式相关基因差分相关识别/tables/模式相关基因差分相关分析结果.xlsx ####
# 输出:
#   1. Results/7.7.4模式相关基因差分相关热图绘制/tables/五模式关键基因绘图数据.xlsx
#   2. Results/7.7.4模式相关基因差分相关热图绘制/figures/五模式关键基因Δρ点热图.pdf
#   3. Results/7.7.4模式相关基因差分相关热图绘制/figures/Δρ热图_主筛选hit.pdf ####

rm(list = ls())
gc()

## 设置默认 graphics device 为 cairo_pdf ####
options(device = function() cairo_pdf(file = tempfile(fileext = ".pdf")))


# 1. 加载配置与包 ####
source("Code/config.R")
library(dplyr)
library(tidyr)
library(readxl)
library(openxlsx)
library(ggplot2)
library(grid)
library(patchwork)


# 2. 路径与参数 ####
script_name <- "7.7.4模式相关基因差分相关热图绘制"
script_results_path <- file.path(results_path, script_name)
table_output_path <- file.path(script_results_path, "tables")
figure_output_path <- file.path(script_results_path, "figures")
log_output_path <- file.path(script_results_path, "logs")

dir.create(script_results_path, recursive = TRUE, showWarnings = FALSE)
dir.create(table_output_path, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_output_path, recursive = TRUE, showWarnings = FALSE)
dir.create(log_output_path, recursive = TRUE, showWarnings = FALSE)

input_xlsx_file <- file.path(
  results_path, "7.7模式相关基因差分相关识别",
  "tables", "模式相关基因差分相关分析结果.xlsx"
)
if (!file.exists(input_xlsx_file)) {
  stop("依赖文件不存在: ", input_xlsx_file, " (请先运行 7.7)")
}

plot_data_xlsx_file <- file.path(table_output_path, "五模式关键基因绘图数据.xlsx")
plot_pdf_file <- file.path(figure_output_path, "五模式关键基因Δρ点热图.pdf")
legacy_plot_pdf_file <- file.path(figure_output_path, "Δρ热图_主筛选hit.pdf")

run_timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_file <- file.path(log_output_path, paste0(script_name, "_", run_timestamp, ".log"))

pattern_order <- c(
  "Ocular",
  "Hematologic",
  "Neurologic",
  "Systemic-like",
  "Visceral-like"
)

pattern_color_map <- c(
  "Systemic-like" = "#F2CC8E",
  "Visceral-like" = "#BC4749",
  "Neurologic"    = "#A15F88",
  "Ocular"        = "#DF7A5E",
  "Hematologic"   = "#82B29A"
)

top_n_per_pattern_direction <- 5
max_neg_log10_perm_p <- 4


# 3. 日志与检查函数 ####

## 函数: log_message ####
## 功能: 同时向控制台和日志文件写入消息 ####
## 输入: message_text,待写入日志的文本 ####
## 输出: 无返回值 ####
log_message <- function(message_text) {
  timestamp_text <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  final_text <- paste0("[", timestamp_text, "] ", message_text)
  message(final_text)
  cat(final_text, "\n", file = log_file, append = TRUE)
}

## 函数: check_required_columns ####
## 功能: 检查数据框是否包含必要列,缺失则终止脚本 ####
## 输入: data_df,待检查数据框; required_columns,必要列名; data_label,数据说明 ####
## 输出: 无返回值 ####
check_required_columns <- function(data_df, required_columns, data_label) {
  missing_columns <- setdiff(required_columns, colnames(data_df))
  if (length(missing_columns) > 0) {
    stop(data_label, " 缺少必要列: ", paste(missing_columns, collapse = ", "))
  }
}

log_message(paste0("========== ", script_name, " 开始 =========="))


# 4. 统计分析: 读取 7.7 结果并筛选关键基因 ####

## 1. 读取数据 ####
log_message("读取 7.7 模式相关基因差分相关识别结果...")
hit_long_df <- read_excel(input_xlsx_file, sheet = "主筛选命中列表_主阈值")
merged_metric_df <- read_excel(input_xlsx_file, sheet = "合并指标全表")

check_required_columns(
  hit_long_df,
  c("pattern", "gene", "hit_direction", "hit_high_confidence",
    "delta_rho", "delta_perm_p_value"),
  "主筛选命中列表_主阈值"
)
check_required_columns(
  merged_metric_df,
  c("pattern", "gene", "rho_tumor", "rho_normal", "delta_rho",
    "fisher_p_value", "delta_perm_p_value"),
  "合并指标全表"
)

hit_long_df <- hit_long_df %>%
  mutate(
    pattern = as.character(pattern),
    gene = as.character(gene),
    hit_direction = as.character(hit_direction),
    hit_high_confidence = as.character(hit_high_confidence) %in% c("TRUE", "True", "true", "1")
  ) %>%
  filter(pattern %in% pattern_order)

merged_metric_df <- merged_metric_df %>%
  mutate(
    pattern = as.character(pattern),
    gene = as.character(gene)
  ) %>%
  filter(pattern %in% pattern_order)

## 2. 每模式按正负方向筛选 Top 基因 ####
top_gene_df <- hit_long_df %>%
  filter(
    hit_direction %in% c("positive", "negative"),
    !is.na(delta_rho)
  ) %>%
  mutate(
    pattern = factor(pattern, levels = pattern_order),
    hit_direction = factor(hit_direction, levels = c("positive", "negative"))
  ) %>%
  group_by(pattern, hit_direction) %>%
  arrange(desc(abs(delta_rho)), delta_perm_p_value, .by_group = TRUE) %>%
  slice_head(n = top_n_per_pattern_direction) %>%
  ungroup() %>%
  arrange(pattern, hit_direction, desc(abs(delta_rho))) %>%
  group_by(pattern) %>%
  mutate(pattern_rank = row_number()) %>%
  ungroup() %>%
  mutate(pattern = as.character(pattern), hit_direction = as.character(hit_direction))

if (nrow(top_gene_df) == 0) {
  stop("未筛选到任何关键基因,请检查 7.7 主筛选命中列表。")
}

## 3. 生成基因排序与绘图长表 ####
gene_order_df <- top_gene_df %>%
  group_by(gene) %>%
  arrange(desc(abs(delta_rho)), delta_perm_p_value, .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  transmute(
    gene,
    primary_pattern = pattern,
    primary_direction = hit_direction,
    primary_delta_rho = delta_rho,
    primary_abs_delta_rho = abs(delta_rho),
    primary_delta_perm_p_value = delta_perm_p_value
  ) %>%
  mutate(
    primary_pattern = factor(primary_pattern, levels = pattern_order),
    primary_direction = factor(primary_direction, levels = c("positive", "negative"))
  ) %>%
  arrange(primary_pattern, primary_direction, desc(primary_abs_delta_rho), gene) %>%
  mutate(gene_order = row_number()) %>%
  mutate(
    primary_pattern = as.character(primary_pattern),
    primary_direction = as.character(primary_direction)
  )

hit_annotation_df <- hit_long_df %>%
  select(gene, pattern, hit_direction, hit_high_confidence) %>%
  distinct()

plot_long_df <- tidyr::expand_grid(
  gene = gene_order_df$gene,
  pattern = pattern_order
) %>%
  left_join(merged_metric_df, by = c("gene", "pattern")) %>%
  left_join(hit_annotation_df, by = c("gene", "pattern")) %>%
  left_join(gene_order_df, by = "gene") %>%
  mutate(
    hit_flag = !is.na(hit_direction),
    hit_direction = if_else(is.na(hit_direction), "not_hit", hit_direction),
    hit_high_confidence = tidyr::replace_na(hit_high_confidence, FALSE),
    neg_log10_perm_p = -log10(pmax(delta_perm_p_value, 1e-300)),
    neg_log10_perm_p_capped = pmin(neg_log10_perm_p, max_neg_log10_perm_p),
    point_alpha = if_else(hit_flag, 0.95, 0.25)
  ) %>%
  arrange(gene_order, factor(pattern, levels = pattern_order))

delta_rho_matrix_df <- plot_long_df %>%
  select(gene, gene_order, pattern, delta_rho) %>%
  arrange(gene_order, factor(pattern, levels = pattern_order)) %>%
  select(-gene_order) %>%
  tidyr::pivot_wider(names_from = pattern, values_from = delta_rho)

plot_parameter_df <- data.frame(
  参数 = c(
    "top_n_per_pattern_direction",
    "pattern_order",
    "max_neg_log10_perm_p"
  ),
  值 = c(
    top_n_per_pattern_direction,
    paste(pattern_order, collapse = " -> "),
    max_neg_log10_perm_p
  )
)


# 5. 结果保存: 保存绘图数据 Excel ####
log_message("保存五模式关键基因绘图数据 Excel...")
openxlsx::write.xlsx(
  x = list(
    "关键基因长表" = plot_long_df,
    "每模式关键基因" = top_gene_df,
    "基因排序" = gene_order_df,
    "Δρ矩阵" = delta_rho_matrix_df,
    "绘图参数" = plot_parameter_df
  ),
  file = plot_data_xlsx_file,
  overwrite = TRUE
)
log_message(paste0("已保存: ", plot_data_xlsx_file))


# 6. 1图: 五模式关键基因 Δρ 点热图 ####

## 1. 读取并处理绘图数据 ####
plot_data <- read_excel(plot_data_xlsx_file, sheet = "关键基因长表")
gene_order_plot_df <- read_excel(plot_data_xlsx_file, sheet = "基因排序")

delta_rho_limit <- max(abs(plot_data$delta_rho), na.rm = TRUE)
if (!is.finite(delta_rho_limit) || delta_rho_limit <= 0) {
  delta_rho_limit <- 0.5
}
delta_rho_limit <- ceiling(delta_rho_limit * 10) / 10

max_gene_order <- max(gene_order_plot_df$gene_order, na.rm = TRUE)
pattern_label_plot_width <- 0.10
pattern_x_spacing <- 0.72
pattern_x_positions <- seq_along(pattern_order) * pattern_x_spacing

gene_axis_df <- gene_order_plot_df %>%
  mutate(
    primary_pattern = factor(primary_pattern, levels = pattern_order),
    y_pos = max_gene_order - gene_order + 1
  ) %>%
  arrange(y_pos)

pattern_label_data <- gene_axis_df %>%
  group_by(primary_pattern) %>%
  summarise(
    ymin = min(y_pos) - 0.5,
    ymax = max(y_pos) + 0.5,
    y = (min(y_pos) + max(y_pos)) / 2,
    .groups = "drop"
  ) %>%
  mutate(
    primary_pattern = as.character(primary_pattern),
    xmin = 0.5,
    xmax = 1.5,
    x = 1
  )

pattern_separator_data <- pattern_label_data %>%
  filter(ymin > 0.5)

plot_data <- plot_data %>%
  mutate(
    pattern = factor(pattern, levels = pattern_order),
    pattern_x = as.numeric(pattern) * pattern_x_spacing,
    y_pos = max_gene_order - gene_order + 1,
    primary_pattern = factor(primary_pattern, levels = pattern_order),
    hit_high_confidence = as.character(hit_high_confidence) %in% c("TRUE", "True", "true", "1")
  )

## 2. 绘图 ####
p_pattern <- ggplot(pattern_label_data) +
  geom_rect(data = pattern_label_data,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                fill = primary_pattern),
            inherit.aes = FALSE, color = "black", linewidth = 0.45, alpha = 0.95,
            show.legend = FALSE) +
  geom_text(data = pattern_label_data,
            aes(x = x, y = y, label = primary_pattern),
            inherit.aes = FALSE, angle = 90, size = 5.2, fontface = "bold",
            color = "black", family = "Arial") +
  scale_fill_manual(values = pattern_color_map, guide = "none") +
  scale_x_continuous(
    position = "top",
    breaks = 1,
    labels = "Systemic-like",
    limits = c(0.5, 1.5),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    limits = c(0.5, max_gene_order + 0.5),
    expand = expansion(mult = c(0, 0))
  ) +
  theme_void(base_family = "Arial", base_size = 12) +
  theme(
    axis.text.x = element_text(size = 12, color = NA, family = "Arial",
                               angle = 45, hjust = 0, vjust = 0),
    plot.margin = margin(t = 10, r = 4, b = 10, l = 4)
  )

p_heatmap <- ggplot(plot_data, aes(x = pattern_x, y = y_pos)) +
  geom_hline(data = pattern_separator_data,
             aes(yintercept = ymin),
             inherit.aes = FALSE,
             color = "grey70", linewidth = 0.45) +
  geom_point(
    aes(size = neg_log10_perm_p_capped, fill = delta_rho, alpha = point_alpha),
    shape = 21, color = "black", stroke = 0.25
  ) +
  geom_point(
    data = filter(plot_data, hit_high_confidence),
    aes(size = neg_log10_perm_p_capped),
    shape = 21, fill = NA, color = "black", stroke = 0.9
  ) +
  scale_fill_gradient2(
    name = "Δρ",
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0,
    limits = c(-delta_rho_limit, delta_rho_limit)
  ) +
  scale_size_continuous(
    name = expression(-log[10](P[perm])),
    range = c(1.8, 6.2),
    breaks = c(1, 2, 3, 4),
    limits = c(0, max_neg_log10_perm_p)
  ) +
  scale_alpha_identity() +
  scale_x_continuous(
    position = "top",
    breaks = pattern_x_positions,
    labels = pattern_order,
    limits = range(pattern_x_positions) + c(-0.28, 0.28),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    breaks = gene_axis_df$y_pos,
    labels = gene_axis_df$gene,
    limits = c(0.5, max_gene_order + 0.5),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_family = "Arial", base_size = 12) +
  theme(
    plot.title = element_text(size = 16, face = "bold", color = "black", family = "Arial"),
    plot.subtitle = element_text(size = 12, color = "black", family = "Arial"),
    axis.text.x = element_text(size = 12, color = "black", family = "Arial",
                               angle = 45, hjust = 0, vjust = 0),
    axis.text.y = element_text(size = 10, color = "black", family = "Arial",
                               margin = margin(r = 4)),
    axis.ticks = element_line(color = "black"),
    panel.grid.major = element_line(color = "grey88", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = 14, color = "black", family = "Arial"),
    legend.text = element_text(size = 12, color = "black", family = "Arial"),
    legend.position = "right",
    plot.margin = margin(t = 10, r = 16, b = 10, l = 4)
  )

p <- p_pattern + p_heatmap +
  plot_layout(widths = c(pattern_label_plot_width, 1))

print(p)

## 3. 保存图形 ####
plot_height <- max(8, nrow(gene_axis_df) * 0.22 + 3)
ggsave(
  filename = plot_pdf_file,
  plot = p,
  width = 6,
  height = 11,
  device = cairo_pdf,
  limitsize = FALSE
)
log_message(paste0("1图 五模式关键基因 Δρ 点热图 已保存: ", plot_pdf_file))

ggsave(
  filename = legacy_plot_pdf_file,
  plot = p,
  width = 6,
  height = plot_height,
  device = cairo_pdf,
  limitsize = FALSE
)
log_message(paste0("1图 兼容旧文件名已保存: ", legacy_plot_pdf_file))


log_message(paste0("========== ", script_name, " 完成 =========="))
