# 7.7.5 模式相关基因识别概览图绘制 ####
# 用途: 总结 7.7 主筛选 hit 基因数与多模式分类,作为 7.7.3 通路富集与 7.7.4 Δρ 热图
#       之间的桥接面板 (识别基因总览,可填补拼图空隙) ####
# 输入: Results/7.7模式相关基因差分相关识别/tables/模式相关基因差分相关分析结果.xlsx
#       sheet: 命中计数_主与高置信_全阈值, 多模式基因分类_主阈值 ####
# 输出:
#   1. Results/7.7.5模式相关基因识别概览图绘制/figures/识别基因数_双向柱状图.pdf
#   2. Results/7.7.5模式相关基因识别概览图绘制/figures/识别基因多模式分类_堆叠条.pdf
#   3. Results/7.7.5模式相关基因识别概览图绘制/figures/模式相关基因识别概览图.pdf ####

rm(list = ls())
gc()

## 设置默认 graphics device 为 cairo_pdf ####
## 用途: nohup/Rscript headless 模式下 print(p) 默认 PDF 设备对 UTF-8 字符 (Δρ 等) 支持不全 ####
options(device = function() cairo_pdf(file = tempfile(fileext = ".pdf")))


# 1. 加载配置与包 ####
source("Code/config.R")
library(gground)
library(dplyr)
library(tidyr)
library(readxl)
library(ggplot2)
library(patchwork)
library(scales)
library(grid)


# 2. 路径与参数 ####
script_name <- "7.7.5模式相关基因识别概览图绘制"
script_results_path <- file.path(results_path, script_name)
figure_output_path <- file.path(script_results_path, "figures")
log_output_path <- file.path(script_results_path, "logs")

dir.create(script_results_path, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_output_path, recursive = TRUE, showWarnings = FALSE)
dir.create(log_output_path, recursive = TRUE, showWarnings = FALSE)

input_xlsx_file <- file.path(
  results_path, "7.7模式相关基因差分相关识别",
  "tables", "模式相关基因差分相关分析结果.xlsx"
)
if (!file.exists(input_xlsx_file)) {
  stop("依赖文件不存在: ", input_xlsx_file, " (请先运行 7.7)")
}

panel_a_pdf_file  <- file.path(figure_output_path, "识别基因数_双向柱状图.pdf")
panel_b_pdf_file  <- file.path(figure_output_path, "识别基因多模式分类_堆叠条.pdf")
overview_pdf_file <- file.path(figure_output_path, "模式相关基因识别概览图.pdf")

run_timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_file <- file.path(log_output_path, paste0(script_name, "_", run_timestamp, ".log"))

## pattern 顺序与配色 (与 7.7.3/7.7.4 一致) ####
pattern_order <- c(
  "Ocular",
  "Hematologic",
  "Neurologic",
  "Systemic-like",
  "Visceral-like"
)
pattern_color_map <- c(
  "Ocular"        = "#DF7A5E",
  "Hematologic"   = "#82B29A",
  "Neurologic"    = "#A15F88",
  "Systemic-like" = "#F2CC8E",
  "Visceral-like" = "#BC4749"
)

## Δρ 方向配色 (与 7.7.4 热图 Δρ 色系一致) ####
positive_color <- "#B2182B"
negative_color <- "#2166AC"

## 多模式分类 (gene_category) 顺序、显示标签与配色 ####
category_order <- c(
  "uni-pattern",
  "multi-pattern-same-direction",
  "anti-correlated"
)
category_label_map <- c(
  "uni-pattern"                  = "Uni-pattern",
  "multi-pattern-same-direction" = "Multi-pattern (same direction)",
  "anti-correlated"              = "Multi-pattern (opposite direction)"
)
category_color_map <- c(
  "uni-pattern"                  = "#8DA0CB",
  "multi-pattern-same-direction" = "#66C2A5",
  "anti-correlated"              = "#FC8D62"
)

## 主阈值 (与 7.7 一致) ####
delta_rho_threshold_main <- 0.3


# 3. 日志辅助函数 ####

## 函数: log_message ####
## 功能: 同时向控制台和日志文件写入消息 ####
## 输入: msg, 待写入日志的文本 ####
## 输出: 无返回值 ####
log_message <- function(msg) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- paste0("[", ts, "] ", msg)
  message(line)
  cat(line, "\n", file = log_file, append = TRUE)
}
log_message(paste0("========== ", script_name, " 开始 =========="))


# 1图: 识别基因数_双向柱状图 ####

## 1. 读取并处理绘图数据 ####
log_message("读取主筛选 hit 计数表, 构造双向柱状图绘图数据...")
hit_count_df <- read_excel(input_xlsx_file, sheet = "命中计数_主与高置信_全阈值") %>%
  mutate(pattern = as.character(pattern)) %>%
  filter(delta_rho_threshold == delta_rho_threshold_main,
         pattern %in% pattern_order)

panel_a_wide_df <- hit_count_df %>%
  transmute(
    pattern        = factor(pattern, levels = rev(pattern_order)),
    positive_count = as.integer(n_positive_hits),
    negative_count = as.integer(n_negative_hits),
    total_count    = positive_count + negative_count
  )

panel_a_long_df <- panel_a_wide_df %>%
  select(pattern, positive_count, negative_count) %>%
  pivot_longer(
    cols      = c(positive_count, negative_count),
    names_to  = "direction_raw",
    values_to = "count"
  ) %>%
  mutate(
    direction = recode(direction_raw,
      "positive_count" = "Positive Δρ",
      "negative_count" = "Negative Δρ"
    ),
    direction = factor(direction, levels = c("Negative Δρ", "Positive Δρ")),
    direction_fill = recode(direction,
      "Negative Δρ" = negative_color,
      "Positive Δρ" = positive_color
    ),
    signed_count = ifelse(direction == "Negative Δρ", -count, count),
    label_hjust  = ifelse(direction == "Negative Δρ", 1.18, -0.18)
  )

## 左右两侧分别取最大: negative 最大 1161, positive 最大 2112, 数据本身不对称 ####
panel_a_neg_max <- max(panel_a_wide_df$negative_count, na.rm = TRUE)
panel_a_pos_max <- max(panel_a_wide_df$positive_count, na.rm = TRUE)

## 不对称 xlim: 左侧 padding 需够容纳 "1161" 标签 (hjust=1.18 向左延伸),
##              右侧 padding 略大以保证 "2112" 标签不贴边 ####
panel_a_xlim <- c(-panel_a_neg_max * 1.60, panel_a_pos_max * 1.38)

## 不对称 breaks: 左 1 个 (-1000), 右 4 个 (500/1000/1500/2000), 0 居中 ####
panel_a_breaks <- c(-1000, 0, 500, 1000, 1500, 2000)

## 模式标签背景框位置: 参考 7.7.3, 在主图左侧使用模式色圆角背景 ####
panel_a_x_range <- diff(panel_a_xlim)
panel_a_pattern_label_gap <- panel_a_x_range * 0.040
panel_a_pattern_label_width <- panel_a_x_range * 0.33
panel_a_pattern_label_xmax <- panel_a_xlim[1] - panel_a_pattern_label_gap
panel_a_pattern_label_xmin <- panel_a_pattern_label_xmax - panel_a_pattern_label_width

panel_a_pattern_label_df <- panel_a_wide_df %>%
  distinct(pattern) %>%
  mutate(
    pattern_label = as.character(pattern),
    y_pos = as.numeric(pattern),
    xmin = panel_a_pattern_label_xmin,
    xmax = panel_a_pattern_label_xmax,
    ymin = y_pos - 0.36,
    ymax = y_pos + 0.36,
    x = (xmin + xmax) / 2,
    y = y_pos,
    pattern_fill = unname(pattern_color_map[pattern_label])
  )

panel_a_long_df <- panel_a_long_df %>%
  mutate(y_pos = as.numeric(pattern))

## 2. 绘图 ####
p_a <- ggplot(panel_a_long_df,
              aes(x = signed_count, y = y_pos)) +
  geom_round_rect(
    data = panel_a_pattern_label_df,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
        fill = pattern_fill),
    radius = unit(2, "mm"), inherit.aes = FALSE,
    color = NA, alpha = 0.95
  ) +
  geom_text(
    data = panel_a_pattern_label_df,
    aes(x = x, y = y, label = pattern_label),
    inherit.aes = FALSE,
    size = 4.8, fontface = "bold", color = "black"
  ) +
  geom_col(aes(fill = direction_fill),
           width = 0.68, color = "black", linewidth = 0.3,
           orientation = "y") +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.45) +
  geom_text(
    aes(label = count, hjust = label_hjust),
    size = 4.0, color = "black"
  ) +
  scale_fill_identity(
    name   = "Δρ direction",
    breaks = c(negative_color, positive_color),
    labels = c("Negative Δρ", "Positive Δρ"),
    guide  = guide_legend(override.aes = list(color = "black", alpha = 1))
  ) +
  scale_x_continuous(
    name   = "Number of genes",
    breaks = panel_a_breaks,
    labels = function(x) abs(x),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    name   = NULL,
    breaks = panel_a_pattern_label_df$y_pos,
    labels = panel_a_pattern_label_df$pattern_label,
    limits = c(0.5, length(pattern_order) + 0.5),
    expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(xlim = panel_a_xlim, clip = "off") +
  theme_bw(base_size = 12) +
  theme(
    plot.title         = element_text(size = 14, face = "bold", color = "black"),
    plot.subtitle      = element_text(size = 12, color = "black"),
    axis.title.x       = element_text(size = 14, color = "black"),
    axis.title.y       = element_text(size = 14, color = "black"),
    axis.text.x        = element_text(size = 12, color = "black"),
    axis.text.y        = element_blank(),
    axis.ticks.y       = element_blank(),
    legend.title       = element_text(size = 14, color = "black"),
    legend.text        = element_text(size = 12, color = "black"),
    legend.position    = "top",
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(color = "grey88", linewidth = 0.25),
    plot.margin        = margin(t = 6, r = 18, b = 6, l = 115)
  )

print(p_a)

## 3. 保存图形 ####
ggsave(
  filename = panel_a_pdf_file,
  plot     = p_a,
  width    = 6,
  height   = 4,
  device   = cairo_pdf,
  limitsize = FALSE
)
log_message(paste0("1图 识别基因数_双向柱状图 已保存: ", panel_a_pdf_file))


# 2图: 识别基因多模式分类_堆叠条 ####

## 1. 读取并处理绘图数据 ####
log_message("读取多模式基因分类表, 构造堆叠条绘图数据...")
multi_pattern_classification_df <- read_excel(
  input_xlsx_file,
  sheet = "多模式基因分类_主阈值"
)

## 注: position_stack(reverse = TRUE) 把第一级 (uni-pattern) 摆在堆叠左端,
##     所以 cumsum 按 gene_category 升序, label_x 才能与实际渲染段中心对齐 ####
panel_b_df <- multi_pattern_classification_df %>%
  filter(gene_category %in% category_order) %>%
  count(gene_category, name = "n_genes") %>%
  mutate(gene_category = factor(gene_category, levels = category_order)) %>%
  arrange(gene_category) %>%
  mutate(
    proportion = n_genes / sum(n_genes),
    cum_end    = cumsum(n_genes),
    cum_start  = cum_end - n_genes,
    label_x    = (cum_start + cum_end) / 2,
    is_big     = n_genes >= 1000,
    label_text = sprintf("%d (%.1f%%)", n_genes, proportion * 100)
  )

panel_b_total       <- sum(panel_b_df$n_genes)
panel_b_inside_df   <- panel_b_df %>% filter(is_big)
panel_b_outside_df  <- panel_b_df %>% filter(!is_big) %>%
  arrange(label_x) %>%
  mutate(
    row_idx    = row_number(),
    ## 统一拐点高度, 使两标签文字对齐 ####
    bend_y     = 1.40,
    ## 第 1 行 (127, multi-pattern): 折线左移; 第 2 行 (568): 直线, 不移 ####
    text_x     = ifelse(row_idx == 1,
                        label_x - panel_b_total * 0.08,
                        label_x),
    text_hjust = ifelse(row_idx == 1, 1, 0.5),
    label_y    = bend_y,
    has_elbow  = (row_idx == 1)
  )

## 手动图例坐标 (相对 panel_b_total) ####
## 用 annotate() 绕开 ggplot 网格图例的列填充, 精确控制 row 1 两项之间的间距 ####
manual_box_w     <- panel_b_total * 0.022   # 颜色方块 x 宽度
manual_box_yh    <- 0.13                    # 颜色方块 y 半高
manual_text_gap  <- panel_b_total * 0.006   # 颜色方块右侧到文字的间距

## 两行布局: Row 1 = Title + Uni + Opp; Row 2 = Same (与 Uni 列对齐) ####
y_legend_row1    <- 2.00
y_legend_row2    <- 1.00
y_legend_title   <- y_legend_row1

## x 坐标 (标题靠左, Uni 在标题右侧, Opp 再右侧; Same 与 Uni 列对齐) ####
x_title          <- 0
x_r1_uni_box1    <- panel_b_total * 0.32
x_r1_uni_box2    <- x_r1_uni_box1 + manual_box_w
x_r1_uni_text    <- x_r1_uni_box2 + manual_text_gap
x_r1_opp_box1    <- panel_b_total * 0.55
x_r1_opp_box2    <- x_r1_opp_box1 + manual_box_w
x_r1_opp_text    <- x_r1_opp_box2 + manual_text_gap

## Row 2: Multi-pattern (same direction) 与 Row 1 Uni 对齐 ####
x_r2_same_box1   <- x_r1_uni_box1
x_r2_same_box2   <- x_r1_uni_box2
x_r2_same_text   <- x_r1_uni_text

## 2. 绘图 ####
## 2.1 手动图例 (独立 ggplot, theme_void, 与柱状图共享 x 轴范围) ####
p_b_legend <- ggplot() +
  ### 标题 "Gene specificity" (14pt) ####
  annotate(
    "text", x = x_title, y = y_legend_title,
    label = "Gene specificity", hjust = 0, vjust = 0.5,
    size  = 14 / .pt, color = "black"
  ) +
  ### Row 1 Item 1: Uni-pattern ####
  annotate(
    "rect",
    xmin = x_r1_uni_box1, xmax = x_r1_uni_box2,
    ymin = y_legend_row1 - manual_box_yh, ymax = y_legend_row1 + manual_box_yh,
    fill = category_color_map[["uni-pattern"]],
    color = "black", linewidth = 0.3
  ) +
  annotate(
    "text", x = x_r1_uni_text, y = y_legend_row1,
    label = "Uni-pattern", hjust = 0, vjust = 0.5,
    size  = 12 / .pt, color = "black"
  ) +
  ### Row 1 Item 2: Multi-pattern (opposite dir.) ####
  annotate(
    "rect",
    xmin = x_r1_opp_box1, xmax = x_r1_opp_box2,
    ymin = y_legend_row1 - manual_box_yh, ymax = y_legend_row1 + manual_box_yh,
    fill = category_color_map[["anti-correlated"]],
    color = "black", linewidth = 0.3
  ) +
  annotate(
    "text", x = x_r1_opp_text, y = y_legend_row1,
    label = "Multi-pattern (opposite direction)", hjust = 0, vjust = 0.5,
    size  = 12 / .pt, color = "black"
  ) +
  ### Row 2: Multi-pattern (same dir.) ####
  annotate(
    "rect",
    xmin = x_r2_same_box1, xmax = x_r2_same_box2,
    ymin = y_legend_row2 - manual_box_yh, ymax = y_legend_row2 + manual_box_yh,
    fill = category_color_map[["multi-pattern-same-direction"]],
    color = "black", linewidth = 0.3
  ) +
  annotate(
    "text", x = x_r2_same_text, y = y_legend_row2,
    label = "Multi-pattern (same direction)", hjust = 0, vjust = 0.5,
    size  = 12 / .pt, color = "black"
  ) +
  scale_x_continuous(
    limits = c(0, panel_b_total),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  scale_y_continuous(
    limits = c(y_legend_row2 - manual_box_yh - 0.2,
               y_legend_row1 + manual_box_yh + 0.2),
    expand = c(0, 0)
  ) +
  theme_void() +
  theme(plot.margin = margin(t = 4, r = 18, b = 8, l = 6))

## 2.2 柱状图 (无图例, 共享 x 轴) ####
p_b_bar <- ggplot(
  panel_b_df,
  aes(x = n_genes, y = "Main-hit genes", fill = gene_category)
) +
  geom_col(
    width     = 0.5,
    color     = "black",
    linewidth = 0.3,
    position  = position_stack(reverse = TRUE)
  ) +
  ### 大段 inside label ####
  geom_text(
    data = panel_b_inside_df,
    aes(x = label_x, label = label_text),
    y = 1, color = "black", size = 3.8,
    inherit.aes = FALSE
  ) +
  ### 垂直引导线 (两个标签都有) ####
  geom_segment(
    data = panel_b_outside_df,
    aes(x = label_x, xend = label_x, yend = bend_y),
    y = 1.25,
    linewidth = 0.3, color = "black",
    inherit.aes = FALSE
  ) +
  ### 折线 - 水平段 (仅 127 multi-pattern) ####
  geom_segment(
    data = panel_b_outside_df %>% filter(has_elbow),
    aes(x = label_x, xend = text_x, y = bend_y, yend = bend_y),
    linewidth = 0.3, color = "black",
    inherit.aes = FALSE
  ) +
  ### 标签文字 (统一高度) ####
  geom_text(
    data = panel_b_outside_df,
    aes(x = text_x, y = label_y, label = label_text, hjust = text_hjust),
    vjust = -0.2,
    color = "black", size = 3.4,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(
    values = category_color_map,
    breaks = category_order,
    labels = category_label_map[category_order],
    guide  = "none"
  ) +
  scale_x_continuous(
    name   = sprintf("Gene-level specificity (n = %d genes)", panel_b_total),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  scale_y_discrete(name = NULL, expand = expansion(add = c(0.6, 1.0))) +
  coord_cartesian(clip = "off") +
  theme_bw(base_size = 12) +
  theme(
    axis.title.x       = element_text(size = 14, color = "black"),
    axis.text.x        = element_text(size = 12, color = "black"),
    axis.text.y        = element_blank(),
    axis.ticks.y       = element_blank(),
    legend.position    = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(color = "grey88", linewidth = 0.25),
    plot.margin        = margin(t = 4, r = 18, b = 10, l = 6)
  )

## 2.3 patchwork 拼接 (图例在上, 柱状图在下, 共享 x 轴宽度) ####
p_b <- p_b_legend / p_b_bar + plot_layout(heights = c(0.55, 1))

print(p_b)

## 3. 保存图形 ####
ggsave(
  filename = panel_b_pdf_file,
  plot     = p_b,
  width    = 6,
  height   = 2,
  device   = cairo_pdf,
  limitsize = FALSE
)
log_message(paste0("2图 识别基因多模式分类_堆叠条 已保存: ", panel_b_pdf_file))


# 3图: 模式相关基因识别概览拼接图 ####

## 1. 拼接绘图数据 (复用前两张图) ####
log_message("拼接 panel A 与 panel B 为模式相关基因识别概览图...")

## 2. 绘图 ####
p_overview <- p_a / p_b + plot_layout(heights = c(2.6, 1))

print(p_overview)

## 3. 保存图形 ####
ggsave(
  filename = overview_pdf_file,
  plot     = p_overview,
  width    = 6,
  height   = 6,
  device   = cairo_pdf,
  limitsize = FALSE
)
log_message(paste0("3图 模式相关基因识别概览拼接图 已保存: ", overview_pdf_file))


log_message(paste0("========== ", script_name, " 完成 =========="))
