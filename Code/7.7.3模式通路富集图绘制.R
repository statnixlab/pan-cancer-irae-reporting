# 7.7.3 模式通路富集图绘制 (蝴蝶图 / Diverging Bar Chart) ####
# 用途: 为 5 个 NMF irAE 模式各绘制一张蝴蝶图,Positive Δρ 柱向右、Negative Δρ 柱向左,
#       共享同一 y 轴。图中元素:
#         - 横向圆角柱 (geom_round_col, orientation = "y", 通路数据库来源色填充)
#         - 通路名 + Top基因名斜体 (Pos 行右展开,Neg 行左展开)
#         - 计数圆圈 (固定在左侧,模式色填充)
#         - 数据库类别作为柱子颜色,并在图例中展示来源名称 ####
# 输入: Results/7.7.2模式通路富集图数据准备/tables/模式通路富集图数据.xlsx ####
# 输出:
#   1. Results/7.7.3模式通路富集图绘制/figures/<pattern>_通路富集.pdf × 5
#   2. Results/7.7.3模式通路富集图绘制/figures/模式通路富集图_拼接版.pdf ####

rm(list = ls())
gc()

## 设置默认 graphics device 为 cairo_pdf ####
## 用途: nohup/Rscript headless 模式下 print(p) 默认 PDF 设备不支持 Arial,会报 "invalid font type" ####
## 该选项让 print(p) 也走 cairo_pdf,与 ggsave(device = cairo_pdf) 一致 ####
options(device = function() cairo_pdf(file = tempfile(fileext = ".pdf")))


# 1. 加载配置与包 ####
source("Code/config.R")
library(gground)
library(ggprism)
library(patchwork)
library(dplyr)
library(tidyr)
library(readxl)
library(ggplot2)
library(grid)


# 2. 路径与参数 ####
script_name <- "7.7.3模式通路富集图绘制"
script_results_path <- file.path(results_path, script_name)
figure_output_path <- file.path(script_results_path, "figures")
log_output_path <- file.path(script_results_path, "logs")

dir.create(script_results_path, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_output_path, recursive = TRUE, showWarnings = FALSE)
dir.create(log_output_path, recursive = TRUE, showWarnings = FALSE)

input_xlsx_file <- file.path(
  results_path, "7.7.2模式通路富集图数据准备",
  "tables", "模式通路富集图数据.xlsx"
)
if (!file.exists(input_xlsx_file)) {
  stop("依赖文件不存在: ", input_xlsx_file, " (请先运行 7.7.2)")
}

run_timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_file <- file.path(log_output_path, paste0(script_name, "_", run_timestamp, ".log"))


# 3. 日志辅助函数 ####

## 函数: log_message ####
## 功能: 同时打印到控制台和日志文件 ####
## 输入: msg,字符串 ####
## 输出: 无返回值 ####
log_message <- function(msg) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- paste0("[", ts, "] ", msg)
  cat(line, "\n")
  write(line, log_file, append = TRUE)
}
log_message(paste0("========== ", script_name, " 开始 =========="))


# 4. 全局视觉参数 ####
## 数据库 5 类柔和色调 (用于通路柱子颜色与图例) ####
pal <- c(
  "GO_BP"          = "#498F99",
  "Hallmark"       = "#95AC66",
  "KEGG"           = "#A9494A",
  "Reactome"       = "#DEA361",
  "C7_IMMUNESIGDB" = "#2f4858"
)
## 5 个 NMF irAE 模式色板 (用于柱填充与计数圆圈填充) ####
pattern_color_map <- c(
  "Ocular"        = "#DF7A5E",
  "Hematologic"   = "#82B29A",
  "Neurologic"    = "#A15F88",
  "Systemic-like" = "#F2CC8E",
  "Visceral-like" = "#BC4749"
)
## 数据库因子顺序: rev 后,rank 1 (GO_BP) 在图顶,rank 5 (C7-IS) 在图底 ####
database_levels_rev <- c("C7_IMMUNESIGDB", "Reactome", "KEGG", "Hallmark", "GO_BP")


# 5. 读取数据 ####
log_message("读取 7.7.2 Excel...")
curated_final_df <- read_excel(input_xlsx_file, sheet = "通路富集长表")
plot_order_df <- read_excel(input_xlsx_file, sheet = "图形排序")
database_order_tbl <- read_excel(input_xlsx_file, sheet = "数据库顺序颜色")

curated_final_df <- curated_final_df %>%
  mutate(
    display_name = case_when(
      display_name %in% c(
        "Oxidative Phosphorylation (KEGG)",
        "Oxidative Phosphorylation (GO)",
        "Oxidative Phosphorylation (Hallmark)"
      ) ~ "Oxidative Phosphorylation",
      TRUE ~ display_name
    )
  )

log_message(sprintf("curated_final_df 行数: %d; plot_order_df 行数: %d; database_order_tbl 行数: %d",
                    nrow(curated_final_df), nrow(plot_order_df), nrow(database_order_tbl)))

## 5.1 五张富集图统一绝对 x 位置参数 ####
## 用途: 固定主图、计数圆圈和左侧模式标签的 x 坐标,保证五张图视觉位置一致 ####
## positive_label_right_padding 控制 x=0 右侧额外留白,数值越小右侧越短 ####
## negative_axis_left_padding 控制 x=0 左侧额外留白,数值越大左侧越长 ####
positive_label_right_padding <- 1.0
negative_axis_left_padding <- 1.0
absolute_pattern_label_width <- 2.3
absolute_count_circle_left_offset <- 1.50
absolute_pattern_label_gap <- 1.65

absolute_xaxis_max_pos <- max(
  curated_final_df$log10_padj[curated_final_df$hit_direction == "positive"],
  na.rm = TRUE
) + positive_label_right_padding
absolute_xaxis_max_neg <- max(
  curated_final_df$log10_padj[curated_final_df$hit_direction == "negative"],
  na.rm = TRUE
) + negative_axis_left_padding

absolute_row_xmin <- -absolute_xaxis_max_neg
absolute_plot_xmax <- absolute_xaxis_max_pos
absolute_circle_x <- absolute_row_xmin - absolute_count_circle_left_offset
absolute_pattern_label_xmax <- absolute_circle_x - absolute_pattern_label_gap
absolute_pattern_label_xmin <- absolute_pattern_label_xmax - absolute_pattern_label_width


# 6. 文本辅助函数 (非绘图函数) ####

## 函数: build_gene_display_text ####
## 功能: 将 "/" 分隔的 gene_top_list 折行,避免单行过长 ####
## 输入: gene_str 字符串,max_line_len 每行最大字符数 ####
## 输出: 折行后字符串 ####
build_gene_display_text <- function(gene_str, max_line_len = 60) {
  if (is.na(gene_str) || gene_str == "") return("")
  genes <- strsplit(gene_str, "/", fixed = TRUE)[[1]]
  lines <- character()
  current_line <- ""
  for (g in genes) {
    candidate <- if (current_line == "") g else paste(current_line, g, sep = "/")
    if (nchar(candidate) > max_line_len && current_line != "") {
      lines <- c(lines, current_line)
      current_line <- g
    } else {
      current_line <- candidate
    }
  }
  if (nzchar(current_line)) lines <- c(lines, current_line)
  paste(lines, collapse = "\n")
}


# 7. 1图: Systemic-like — 蝴蝶图 (Positive 向右 / Negative 向左) ####

## 1. 读取并合并 Pos + Neg 绘图数据 ####
focus_pattern       <- "Systemic-like"
focus_pattern_label <- "Systemic-like"
pattern_fill        <- pattern_color_map[[focus_pattern]]

## 1.1 本图布局参数 ####
## 本图使用 5.1 中统一的绝对 x 位置参数 ####

plot_data <- bind_rows(
  curated_final_df %>%
    filter(pattern == focus_pattern, hit_direction == "positive") %>%
    mutate(direction = "positive"),
  curated_final_df %>%
    filter(pattern == focus_pattern, hit_direction == "negative") %>%
    mutate(direction = "negative")
) %>%
  left_join(database_order_tbl %>% select(database, database_rank, database_label),
            by = "database") %>%
  mutate(direction_f = factor(direction, levels = c("positive", "negative"))) %>%
  arrange(direction_f, database_rank, desc(log10_padj)) %>%
  mutate(
    log10_padj_signed = if_else(direction == "positive", log10_padj, -log10_padj),
    index = rev(seq_len(n())),
    database = factor(database, levels = database_levels_rev)
  )

n_pos <- sum(plot_data$direction == "positive")
n_neg <- sum(plot_data$direction == "negative")

xaxis_max_pos <- absolute_xaxis_max_pos
xaxis_max_neg <- absolute_xaxis_max_neg

## 2. 计算模式标签与计数圆圈几何位置 ####
row_xmin <- absolute_row_xmin
circle_x <- absolute_circle_x
pattern_label_xmin <- absolute_pattern_label_xmin
pattern_label_xmax <- absolute_pattern_label_xmax
plot_xmin <- row_xmin
plot_xmax <- absolute_plot_xmax

pattern_label_data <- tibble::tibble(
  xmin = pattern_label_xmin,
  xmax = pattern_label_xmax,
  ymin = 0.45,
  ymax = max(plot_data$index) + 0.55,
  label = focus_pattern_label
)

## 3. 绘图 ####
p_systemic <- ggplot(plot_data) +
  ## 模式标签 — 放在主图左侧并使用模式色填充 ####
  geom_round_rect(data = pattern_label_data,
                  aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                  radius = unit(2, "mm"), inherit.aes = FALSE,
                  fill = pattern_fill, color = NA, alpha = 0.95) +
  geom_text(data = pattern_label_data,
            aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2,
                label = label),
            inherit.aes = FALSE, angle = 90, size = 5.2,
            fontface = "bold", color = "black", family = "Arial") +
  ## 横向圆角柱 (有符号 x): pos 向右 / neg 向左,按通路数据库上色 ####
  geom_round_col(aes(x = log10_padj_signed, y = index, fill = database),
                 width = 0.78, alpha = 0.78,
                 orientation = "y") +
  ## 通路名 — Positive 行,从 x=0 向右展开 ####
  geom_text(data = filter(plot_data, direction == "positive"),
            aes(x = xaxis_max_pos * 0.005, y = index, label = display_name),
            hjust = 0, vjust = 0.5, size = 4.0, family = "Arial", color = "black") +
  ## 通路名 — Negative 行,从 x=0 向左展开 ####
  geom_text(data = filter(plot_data, direction == "negative"),
            aes(x = -xaxis_max_neg * 0.005, y = index, label = display_name),
            hjust = 1, vjust = 0.5, size = 4.0, family = "Arial", color = "black") +
  ## 计数圆圈 — 全部固定在最左侧 ####
  geom_point(aes(x = circle_x, y = index, size = n_hit_in_term),
             shape = 21, fill = pattern_fill, color = "black") +
  geom_text(aes(x = circle_x, y = index, label = n_hit_in_term),
            size = 3.0, fontface = "bold", color = "black", family = "Arial") +
  ## x = 0 中线 (虚线区分正负方向) ####
  annotate("segment", x = 0, xend = 0, y = 0.3,
           yend = max(plot_data$index) + 0.7,
           linewidth = 0.4, color = "grey50", linetype = "dashed") +
  scale_fill_manual(name = "Gene-set collection", values = pal,
                    breaks = c("GO_BP", "Hallmark", "KEGG", "Reactome", "C7_IMMUNESIGDB"),
                    labels = c("GO BP", "MSigDB Hallmark",
                               "KEGG pathways", "Reactome pathways",
                               "C7-IS"),
                    drop = FALSE) +
  scale_size_continuous(name = "Count", range = c(4, 10),
                        breaks = c(5, 10, 20, 50)) +
  guides(fill = guide_legend(order = 1, override.aes = list(alpha = 0.78)),
         size = guide_legend(order = 2)) +
  scale_x_continuous(
    breaks = pretty(c(-xaxis_max_neg, xaxis_max_pos), n = 8),
    labels = function(b) abs(round(b, 1)),
    expand = expansion(c(0.02, 0.02))
  ) +
  scale_y_continuous(limits = c(0, max(plot_data$index) + 1), expand = c(0, 0)) +
  coord_cartesian(xlim = c(plot_xmin, plot_xmax), clip = "off") +
  labs(x = "← Negative              -log10(p.adjust)              Positive →",
       y = NULL) +
  theme_prism(base_family = "Arial", base_size = 12) +
  theme(
    plot.title   = element_text(size = 16, face = "bold", color = "black", family = "Arial"),
    axis.title.x = element_text(size = 14, color = "black", family = "Arial",
                                margin = margin(t = 8)),
    axis.text.x  = element_text(size = 12, color = "black", family = "Arial",
                                margin = margin(t = 4)),
    axis.line.x  = element_line(linewidth = 1.2, color = "black"),
    axis.ticks.x = element_line(linewidth = 1.0, color = "black"),
    axis.ticks.length = unit(3, "mm"),
    axis.text.y  = element_blank(),
    axis.line.y  = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(size = 14, color = "black", family = "Arial"),
    legend.text  = element_text(size = 12, color = "black", family = "Arial"),
    legend.position = "right",
    plot.margin = margin(t = 8, r = 20, b = 8, l = 50)
  )

print(p_systemic)

## 4. 保存图形 ####
panel_height <- 10  ## 单图保存高度 (英寸),按需手动调节 ####
ggsave(
  filename = file.path(figure_output_path,
                       sprintf("%s_通路富集.pdf", focus_pattern)),
  plot = p_systemic, width = 9, height = 4, device = cairo_pdf
)
log_message(sprintf("1图 %s 已保存 (h=%.1f)", focus_pattern, panel_height))


# 8. 2图: Visceral-like — 蝴蝶图 (Positive 向右 / Negative 向左) ####

## 1. 读取并合并 Pos + Neg 绘图数据 ####
focus_pattern       <- "Visceral-like"
focus_pattern_label <- "Visceral-like"
pattern_fill        <- pattern_color_map[[focus_pattern]]

## 1.1 本图布局参数 ####
## 本图使用 5.1 中统一的绝对 x 位置参数 ####

plot_data <- bind_rows(
  curated_final_df %>%
    filter(pattern == focus_pattern, hit_direction == "positive") %>%
    mutate(direction = "positive"),
  curated_final_df %>%
    filter(pattern == focus_pattern, hit_direction == "negative") %>%
    mutate(direction = "negative")
) %>%
  left_join(database_order_tbl %>% select(database, database_rank, database_label),
            by = "database") %>%
  mutate(direction_f = factor(direction, levels = c("positive", "negative"))) %>%
  arrange(direction_f, database_rank, desc(log10_padj)) %>%
  mutate(
    log10_padj_signed = if_else(direction == "positive", log10_padj, -log10_padj),
    index = rev(seq_len(n())),
    database = factor(database, levels = database_levels_rev)
  )

n_pos <- sum(plot_data$direction == "positive")
n_neg <- sum(plot_data$direction == "negative")

xaxis_max_pos <- absolute_xaxis_max_pos
xaxis_max_neg <- absolute_xaxis_max_neg

## 2. 计算模式标签与计数圆圈几何位置 ####
row_xmin <- absolute_row_xmin
circle_x <- absolute_circle_x
pattern_label_xmin <- absolute_pattern_label_xmin
pattern_label_xmax <- absolute_pattern_label_xmax
plot_xmin <- row_xmin
plot_xmax <- absolute_plot_xmax

pattern_label_data <- tibble::tibble(
  xmin = pattern_label_xmin,
  xmax = pattern_label_xmax,
  ymin = 0.45,
  ymax = max(plot_data$index) + 0.55,
  label = focus_pattern_label
)

## 3. 绘图 ####
p_visceral <- ggplot(plot_data) +
  ## 模式标签 — 放在主图左侧并使用模式色填充 ####
  geom_round_rect(data = pattern_label_data,
                  aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                  radius = unit(2, "mm"), inherit.aes = FALSE,
                  fill = pattern_fill, color = NA, alpha = 0.95) +
  geom_text(data = pattern_label_data,
            aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2,
                label = label),
            inherit.aes = FALSE, angle = 90, size = 5.2,
            fontface = "bold", color = "black", family = "Arial") +
  ## 横向圆角柱 (有符号 x): pos 向右 / neg 向左,按通路数据库上色 ####
  geom_round_col(aes(x = log10_padj_signed, y = index, fill = database),
                 width = 0.78, alpha = 0.78,
                 orientation = "y") +
  ## 通路名 — Positive 行,从 x=0 向右展开 ####
  geom_text(data = filter(plot_data, direction == "positive"),
            aes(x = xaxis_max_pos * 0.005, y = index, label = display_name),
            hjust = 0, vjust = 0.5, size = 4.0, family = "Arial", color = "black") +
  ## 通路名 — Negative 行,从 x=0 向左展开 ####
  geom_text(data = filter(plot_data, direction == "negative"),
            aes(x = -xaxis_max_neg * 0.005, y = index, label = display_name),
            hjust = 1, vjust = 0.5, size = 4.0, family = "Arial", color = "black") +
  ## 计数圆圈 — 全部固定在最左侧 ####
  geom_point(aes(x = circle_x, y = index, size = n_hit_in_term),
             shape = 21, fill = pattern_fill, color = "black") +
  geom_text(aes(x = circle_x, y = index, label = n_hit_in_term),
            size = 3.0, fontface = "bold", color = "black", family = "Arial") +
  ## x = 0 中线 (虚线区分正负方向) ####
  annotate("segment", x = 0, xend = 0, y = 0.3,
           yend = max(plot_data$index) + 0.7,
           linewidth = 0.4, color = "grey50", linetype = "dashed") +
  scale_fill_manual(name = "Gene-set collection", values = pal,
                    breaks = c("GO_BP", "Hallmark", "KEGG", "Reactome", "C7_IMMUNESIGDB"),
                    labels = c("GO BP", "MSigDB Hallmark",
                               "KEGG pathways", "Reactome pathways",
                               "C7-IS"),
                    drop = FALSE) +
  scale_size_continuous(name = "Count", range = c(4, 10),
                        breaks = c(5, 10, 20, 50)) +
  guides(fill = guide_legend(order = 1, override.aes = list(alpha = 0.78)),
         size = guide_legend(order = 2)) +
  scale_x_continuous(
    breaks = pretty(c(-xaxis_max_neg, xaxis_max_pos), n = 8),
    labels = function(b) abs(round(b, 1)),
    expand = expansion(c(0.02, 0.02))
  ) +
  scale_y_continuous(limits = c(0, max(plot_data$index) + 1), expand = c(0, 0)) +
  coord_cartesian(xlim = c(plot_xmin, plot_xmax), clip = "off") +
  labs(x = "← Negative              -log10(p.adjust)              Positive →",
       y = NULL) +
  theme_prism(base_family = "Arial", base_size = 12) +
  theme(
    plot.title   = element_text(size = 16, face = "bold", color = "black", family = "Arial"),
    axis.title.x = element_text(size = 14, color = "black", family = "Arial",
                                margin = margin(t = 8)),
    axis.text.x  = element_text(size = 12, color = "black", family = "Arial",
                                margin = margin(t = 4)),
    axis.line.x  = element_line(linewidth = 1.2, color = "black"),
    axis.ticks.x = element_line(linewidth = 1.0, color = "black"),
    axis.ticks.length = unit(3, "mm"),
    axis.text.y  = element_blank(),
    axis.line.y  = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(size = 14, color = "black", family = "Arial"),
    legend.text  = element_text(size = 12, color = "black", family = "Arial"),
    legend.position = "right",
    plot.margin = margin(t = 8, r = 20, b = 8, l = 50)
  )

print(p_visceral)

## 4. 保存图形 ####
panel_height <- 10  ## 单图保存高度 (英寸),按需手动调节 ####
ggsave(
  filename = file.path(figure_output_path,
                       sprintf("%s_通路富集.pdf", focus_pattern)),
  plot = p_visceral, width = 9, height = 4, device = cairo_pdf
)
log_message(sprintf("2图 %s 已保存 (h=%.1f)", focus_pattern, 4))


# 9. 3图: Neurologic — 蝴蝶图 (Positive 向右 / Negative 向左) ####

## 1. 读取并合并 Pos + Neg 绘图数据 ####
focus_pattern       <- "Neurologic"
focus_pattern_label <- "Neurologic"
pattern_fill        <- pattern_color_map[[focus_pattern]]

## 1.1 本图布局参数 ####
## 本图使用 5.1 中统一的绝对 x 位置参数 ####

plot_data <- bind_rows(
  curated_final_df %>%
    filter(pattern == focus_pattern, hit_direction == "positive") %>%
    mutate(direction = "positive"),
  curated_final_df %>%
    filter(pattern == focus_pattern, hit_direction == "negative") %>%
    mutate(direction = "negative")
) %>%
  left_join(database_order_tbl %>% select(database, database_rank, database_label),
            by = "database") %>%
  mutate(direction_f = factor(direction, levels = c("positive", "negative"))) %>%
  arrange(direction_f, database_rank, desc(log10_padj)) %>%
  mutate(
    log10_padj_signed = if_else(direction == "positive", log10_padj, -log10_padj),
    index = rev(seq_len(n())),
    database = factor(database, levels = database_levels_rev)
  )

n_pos <- sum(plot_data$direction == "positive")
n_neg <- sum(plot_data$direction == "negative")

xaxis_max_pos <- absolute_xaxis_max_pos
xaxis_max_neg <- absolute_xaxis_max_neg

## 2. 计算模式标签与计数圆圈几何位置 ####
row_xmin <- absolute_row_xmin
circle_x <- absolute_circle_x
pattern_label_xmin <- absolute_pattern_label_xmin
pattern_label_xmax <- absolute_pattern_label_xmax
plot_xmin <- row_xmin
plot_xmax <- absolute_plot_xmax

pattern_label_data <- tibble::tibble(
  xmin = pattern_label_xmin,
  xmax = pattern_label_xmax,
  ymin = 0.45,
  ymax = max(plot_data$index) + 0.55,
  label = focus_pattern_label
)

## 3. 绘图 ####
p_neuro <- ggplot(plot_data) +
  ## 模式标签 — 放在主图左侧并使用模式色填充 ####
  geom_round_rect(data = pattern_label_data,
                  aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                  radius = unit(2, "mm"), inherit.aes = FALSE,
                  fill = pattern_fill, color = NA, alpha = 0.95) +
  geom_text(data = pattern_label_data,
            aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2,
                label = label),
            inherit.aes = FALSE, angle = 90, size = 5.2,
            fontface = "bold", color = "black", family = "Arial") +
  ## 横向圆角柱 (有符号 x): pos 向右 / neg 向左,按通路数据库上色 ####
  geom_round_col(aes(x = log10_padj_signed, y = index, fill = database),
                 width = 0.78, alpha = 0.78,
                 orientation = "y") +
  ## 通路名 — Positive 行,从 x=0 向右展开 ####
  geom_text(data = filter(plot_data, direction == "positive"),
            aes(x = xaxis_max_pos * 0.005, y = index, label = display_name),
            hjust = 0, vjust = 0.5, size = 4.0, family = "Arial", color = "black") +
  ## 通路名 — Negative 行,从 x=0 向左展开 ####
  geom_text(data = filter(plot_data, direction == "negative"),
            aes(x = -xaxis_max_neg * 0.005, y = index, label = display_name),
            hjust = 1, vjust = 0.5, size = 4.0, family = "Arial", color = "black") +
  ## 计数圆圈 — 全部固定在最左侧 ####
  geom_point(aes(x = circle_x, y = index, size = n_hit_in_term),
             shape = 21, fill = pattern_fill, color = "black") +
  geom_text(aes(x = circle_x, y = index, label = n_hit_in_term),
            size = 3.0, fontface = "bold", color = "black", family = "Arial") +
  ## x = 0 中线 (虚线区分正负方向) ####
  annotate("segment", x = 0, xend = 0, y = 0.3,
           yend = max(plot_data$index) + 0.7,
           linewidth = 0.4, color = "grey50", linetype = "dashed") +
  scale_fill_manual(name = "Gene-set collection", values = pal,
                    breaks = c("GO_BP", "Hallmark", "KEGG", "Reactome"),
                    labels = c("GO BP", "MSigDB Hallmark",
                               "KEGG pathways", "Reactome pathways"),
                    drop = FALSE) +
  scale_size_continuous(name = "Count", range = c(4, 10),
                        breaks = c(5, 10, 20, 50)) +
  guides(fill = guide_legend(order = 1, override.aes = list(alpha = 0.78)),
         size = guide_legend(order = 2)) +
  scale_x_continuous(
    breaks = pretty(c(-xaxis_max_neg, xaxis_max_pos), n = 8),
    labels = function(b) abs(round(b, 1)),
    expand = expansion(c(0.02, 0.02))
  ) +
  scale_y_continuous(limits = c(0, max(plot_data$index) + 1), expand = c(0, 0)) +
  coord_cartesian(xlim = c(plot_xmin, plot_xmax), clip = "off") +
  labs(x = "← Negative              -log10(p.adjust)              Positive →",
       y = NULL) +
  theme_prism(base_family = "Arial", base_size = 12) +
  theme(
    plot.title   = element_text(size = 16, face = "bold", color = "black", family = "Arial"),
    axis.title.x = element_text(size = 14, color = "black", family = "Arial",
                                margin = margin(t = 8)),
    axis.text.x  = element_text(size = 12, color = "black", family = "Arial",
                                margin = margin(t = 4)),
    axis.line.x  = element_line(linewidth = 1.2, color = "black"),
    axis.ticks.x = element_line(linewidth = 1.0, color = "black"),
    axis.ticks.length = unit(3, "mm"),
    axis.text.y  = element_blank(),
    axis.line.y  = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(size = 14, color = "black", family = "Arial"),
    legend.text  = element_text(size = 12, color = "black", family = "Arial"),
    legend.position = "right",
    legend.box = "vertical",
    legend.box.margin = margin(t = 34, r = 0, b = 0, l = 0),
    legend.spacing.y = unit(2, "pt"),
    plot.margin = margin(t = 8, r = 20, b = 8, l = 50)
  )

print(p_neuro)

## 4. 保存图形 ####
panel_height <- 10  ## 单图保存高度 (英寸),按需手动调节 ####
ggsave(
  filename = file.path(figure_output_path,
                       sprintf("%s_通路富集.pdf", focus_pattern)),
  plot = p_neuro, width = 9, height = 3, device = cairo_pdf
)
log_message(sprintf("3图 %s 已保存 (h=%.1f)", focus_pattern, panel_height))


# 10. 4图: Ocular — 蝴蝶图 (Positive 向右 / Negative 向左) ####

## 1. 读取并合并 Pos + Neg 绘图数据 ####
focus_pattern       <- "Ocular"
focus_pattern_label <- "Ocular"
pattern_fill        <- pattern_color_map[[focus_pattern]]

## 1.1 本图布局参数 ####
## 本图使用 5.1 中统一的绝对 x 位置参数 ####

plot_data <- bind_rows(
  curated_final_df %>%
    filter(pattern == focus_pattern, hit_direction == "positive") %>%
    mutate(direction = "positive"),
  curated_final_df %>%
    filter(pattern == focus_pattern, hit_direction == "negative") %>%
    mutate(direction = "negative")
) %>%
  left_join(database_order_tbl %>% select(database, database_rank, database_label),
            by = "database") %>%
  mutate(direction_f = factor(direction, levels = c("positive", "negative"))) %>%
  arrange(direction_f, database_rank, desc(log10_padj)) %>%
  mutate(
    log10_padj_signed = if_else(direction == "positive", log10_padj, -log10_padj),
    index = rev(seq_len(n())),
    database = factor(database, levels = database_levels_rev)
  )

n_pos <- sum(plot_data$direction == "positive")
n_neg <- sum(plot_data$direction == "negative")

xaxis_max_pos <- max(
  max(plot_data$log10_padj[plot_data$direction == "positive"], na.rm = TRUE) + 0.8,
  absolute_xaxis_max_pos - 4.0
)
xaxis_max_neg <- absolute_xaxis_max_neg + 4.0

## 2. 计算模式标签与计数圆圈几何位置 ####
row_xmin <- -xaxis_max_neg
circle_x <- row_xmin - absolute_count_circle_left_offset
pattern_label_xmax <- circle_x - absolute_pattern_label_gap
pattern_label_xmin <- pattern_label_xmax - absolute_pattern_label_width
plot_xmin <- row_xmin
plot_xmax <- xaxis_max_pos

pattern_label_data <- tibble::tibble(
  xmin = pattern_label_xmin,
  xmax = pattern_label_xmax,
  ymin = 0.45,
  ymax = max(plot_data$index) + 0.55,
  label = focus_pattern_label
)

## 3. 绘图 ####
p_ocular <- ggplot(plot_data) +
  ## 模式标签 — 放在主图左侧并使用模式色填充 ####
  geom_round_rect(data = pattern_label_data,
                  aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                  radius = unit(2, "mm"), inherit.aes = FALSE,
                  fill = pattern_fill, color = NA, alpha = 0.95) +
  geom_text(data = pattern_label_data,
            aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2,
                label = label),
            inherit.aes = FALSE, angle = 90, size = 5.2,
            fontface = "bold", color = "black", family = "Arial") +
  ## 横向圆角柱 (有符号 x): pos 向右 / neg 向左,按通路数据库上色 ####
  geom_round_col(aes(x = log10_padj_signed, y = index, fill = database),
                 width = 0.78, alpha = 0.78,
                 orientation = "y") +
  ## 通路名 — Positive 行,从 x=0 向右展开 ####
  geom_text(data = filter(plot_data, direction == "positive"),
            aes(x = xaxis_max_pos * 0.005, y = index, label = display_name),
            hjust = 0, vjust = 0.5, size = 4.0, family = "Arial", color = "black") +
  ## 通路名 — Negative 行,从 x=0 向左展开 ####
  geom_text(data = filter(plot_data, direction == "negative"),
            aes(x = -xaxis_max_neg * 0.005, y = index, label = display_name),
            hjust = 1, vjust = 0.5, size = 4.0, family = "Arial", color = "black") +
  ## 计数圆圈 — 全部固定在最左侧 ####
  geom_point(aes(x = circle_x, y = index, size = n_hit_in_term),
             shape = 21, fill = pattern_fill, color = "black") +
  geom_text(aes(x = circle_x, y = index, label = n_hit_in_term),
            size = 3.0, fontface = "bold", color = "black", family = "Arial") +
  ## x = 0 中线 (虚线区分正负方向) ####
  annotate("segment", x = 0, xend = 0, y = 0.3,
           yend = max(plot_data$index) + 0.7,
           linewidth = 0.4, color = "grey50", linetype = "dashed") +
  scale_fill_manual(name = "Gene-set collection", values = pal,
                    breaks = c("GO_BP", "KEGG", "Reactome"),
                    labels = c("GO BP", "KEGG pathways", "Reactome pathways"),
                    drop = FALSE) +
  scale_size_continuous(name = "Count", range = c(4, 10),
                        breaks = c(5, 10, 20, 50)) +
  guides(fill = guide_legend(order = 1, override.aes = list(alpha = 0.78)),
         size = guide_legend(order = 2)) +
  scale_x_continuous(
    breaks = pretty(c(-xaxis_max_neg, xaxis_max_pos), n = 8),
    labels = function(b) abs(round(b, 1)),
    expand = expansion(c(0.02, 0.02))
  ) +
  scale_y_continuous(limits = c(0, max(plot_data$index) + 1), expand = c(0, 0)) +
  coord_cartesian(xlim = c(plot_xmin, plot_xmax), clip = "off") +
  labs(x = "← Negative              -log10(p.adjust)              Positive →",
       y = NULL) +
  theme_prism(base_family = "Arial", base_size = 12) +
  theme(
    plot.title   = element_text(size = 16, face = "bold", color = "black", family = "Arial"),
    axis.title.x = element_text(size = 14, color = "black", family = "Arial",
                                margin = margin(t = 8)),
    axis.text.x  = element_text(size = 12, color = "black", family = "Arial",
                                margin = margin(t = 4)),
    axis.line.x  = element_line(linewidth = 1.2, color = "black"),
    axis.ticks.x = element_line(linewidth = 1.0, color = "black"),
    axis.ticks.length = unit(3, "mm"),
    axis.text.y  = element_blank(),
    axis.line.y  = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(size = 14, color = "black", family = "Arial"),
    legend.text  = element_text(size = 12, color = "black", family = "Arial"),
    legend.position = "right",
    legend.box = "vertical",
    legend.box.margin = margin(t = 28, r = 0, b = 0, l = 0),
    legend.spacing.y = unit(2, "pt"),
    plot.margin = margin(t = 8, r = 20, b = 8, l = 50)
  )

print(p_ocular)

## 4. 保存图形 ####
panel_height <- 10  ## 单图保存高度 (英寸),按需手动调节 ####
ggsave(
  filename = file.path(figure_output_path,
                       sprintf("%s_通路富集.pdf", focus_pattern)),
  plot = p_ocular, width = 9, height = 3, device = cairo_pdf
)
log_message(sprintf("4图 %s 已保存 (h=%.1f)", focus_pattern, panel_height))


# 11. 5图: Hematologic — 蝴蝶图 (Positive 向右 / Negative 向左) ####

## 1. 读取并合并 Pos + Neg 绘图数据 ####
focus_pattern       <- "Hematologic"
focus_pattern_label <- "Hematologic"
pattern_fill        <- pattern_color_map[[focus_pattern]]

## 1.1 本图布局参数 ####
## 本图使用 5.1 中统一的绝对 x 位置参数 ####

plot_data <- bind_rows(
  curated_final_df %>%
    filter(pattern == focus_pattern, hit_direction == "positive") %>%
    mutate(direction = "positive"),
  curated_final_df %>%
    filter(pattern == focus_pattern, hit_direction == "negative") %>%
    mutate(direction = "negative")
) %>%
  left_join(database_order_tbl %>% select(database, database_rank, database_label),
            by = "database") %>%
  mutate(direction_f = factor(direction, levels = c("positive", "negative"))) %>%
  arrange(direction_f, database_rank, desc(log10_padj)) %>%
  mutate(
    log10_padj_signed = if_else(direction == "positive", log10_padj, -log10_padj),
    index = rev(seq_len(n())),
    database = factor(database, levels = database_levels_rev)
  )

n_pos <- sum(plot_data$direction == "positive")
n_neg <- sum(plot_data$direction == "negative")

xaxis_max_pos <- absolute_xaxis_max_pos
xaxis_max_neg <- absolute_xaxis_max_neg

## 2. 计算模式标签与计数圆圈几何位置 ####
row_xmin <- absolute_row_xmin
circle_x <- absolute_circle_x
pattern_label_xmin <- absolute_pattern_label_xmin
pattern_label_xmax <- absolute_pattern_label_xmax
plot_xmin <- row_xmin
plot_xmax <- absolute_plot_xmax

pattern_label_data <- tibble::tibble(
  xmin = pattern_label_xmin,
  xmax = pattern_label_xmax,
  ymin = 0.45,
  ymax = max(plot_data$index) + 0.55,
  label = focus_pattern_label
)

## 3. 绘图 ####
p_hema <- ggplot(plot_data) +
  ## 模式标签 — 放在主图左侧并使用模式色填充 ####
  geom_round_rect(data = pattern_label_data,
                  aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                  radius = unit(2, "mm"), inherit.aes = FALSE,
                  fill = pattern_fill, color = NA, alpha = 0.95) +
  geom_text(data = pattern_label_data,
            aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2,
                label = label),
            inherit.aes = FALSE, angle = 90, size = 5.2,
            fontface = "bold", color = "black", family = "Arial") +
  ## 横向圆角柱 (有符号 x): pos 向右 / neg 向左,按通路数据库上色 ####
  geom_round_col(aes(x = log10_padj_signed, y = index, fill = database),
                 width = 0.78, alpha = 0.78,
                 orientation = "y") +
  ## 通路名 — Positive 行,从 x=0 向右展开 ####
  geom_text(data = filter(plot_data, direction == "positive"),
            aes(x = xaxis_max_pos * 0.005, y = index, label = display_name),
            hjust = 0, vjust = 0.5, size = 4.0, family = "Arial", color = "black") +
  ## 通路名 — Negative 行,从 x=0 向左展开 ####
  geom_text(data = filter(plot_data, direction == "negative"),
            aes(x = -xaxis_max_neg * 0.005, y = index, label = display_name),
            hjust = 1, vjust = 0.5, size = 4.0, family = "Arial", color = "black") +
  ## 计数圆圈 — 全部固定在最左侧 ####
  geom_point(aes(x = circle_x, y = index, size = n_hit_in_term),
             shape = 21, fill = pattern_fill, color = "black") +
  geom_text(aes(x = circle_x, y = index, label = n_hit_in_term),
            size = 3.0, fontface = "bold", color = "black", family = "Arial") +
  ## x = 0 中线 (虚线区分正负方向) ####
  annotate("segment", x = 0, xend = 0, y = 0.3,
           yend = max(plot_data$index) + 0.7,
           linewidth = 0.4, color = "grey50", linetype = "dashed") +
  scale_fill_manual(name = "Gene-set collection", values = pal,
                    breaks = c("GO_BP", "Hallmark", "KEGG", "Reactome", "C7_IMMUNESIGDB"),
                    labels = c("GO BP", "MSigDB Hallmark",
                               "KEGG pathways", "Reactome pathways",
                               "C7-IS"),
                    drop = FALSE) +
  scale_size_continuous(name = "Count", range = c(4, 10),
                        breaks = c(5, 10, 20, 50)) +
  guides(fill = guide_legend(order = 1, override.aes = list(alpha = 0.78)),
         size = guide_legend(order = 2)) +
  scale_x_continuous(
    breaks = pretty(c(-xaxis_max_neg, xaxis_max_pos), n = 8),
    labels = function(b) abs(round(b, 1)),
    expand = expansion(c(0.02, 0.02))
  ) +
  scale_y_continuous(limits = c(0, max(plot_data$index) + 1), expand = c(0, 0)) +
  coord_cartesian(xlim = c(plot_xmin, plot_xmax), clip = "off") +
  labs(x = "← Negative              -log10(p.adjust)              Positive →",
       y = NULL) +
  theme_prism(base_family = "Arial", base_size = 12) +
  theme(
    plot.title   = element_text(size = 16, face = "bold", color = "black", family = "Arial"),
    axis.title.x = element_text(size = 14, color = "black", family = "Arial",
                                margin = margin(t = 8)),
    axis.text.x  = element_text(size = 12, color = "black", family = "Arial",
                                margin = margin(t = 4)),
    axis.line.x  = element_line(linewidth = 1.2, color = "black"),
    axis.ticks.x = element_line(linewidth = 1.0, color = "black"),
    axis.ticks.length = unit(3, "mm"),
    axis.text.y  = element_blank(),
    axis.line.y  = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(size = 14, color = "black", family = "Arial"),
    legend.text  = element_text(size = 12, color = "black", family = "Arial"),
    legend.position = "right",
    legend.box = "vertical",
    legend.box.margin = margin(t = 28, r = 0, b = 0, l = 0),
    legend.spacing.y = unit(2, "pt"),
    plot.margin = margin(t = 8, r = 20, b = 8, l = 50)
  )

print(p_hema)

## 4. 保存图形 ####
panel_height <- 10  ## 单图保存高度 (英寸),按需手动调节 ####
ggsave(
  filename = file.path(figure_output_path,
                       sprintf("%s_通路富集.pdf", focus_pattern)),
  plot = p_hema, width = 9, height = 3, device = cairo_pdf
)
log_message(sprintf("5图 %s 已保存 (h=%.1f)", focus_pattern, panel_height))


# 12. 6图: 拼接版 (5 patterns 蝴蝶图按 TTO 升序垂直拼接) ####

## 1. 计算每个模式的行数 (用于 patchwork heights 加权) ####
n_rows_per_pattern <- curated_final_df %>%
  group_by(pattern) %>%
  summarise(n_rows = n(), .groups = "drop")

get_n_rows <- function(focus) {
  n_rows_per_pattern$n_rows[n_rows_per_pattern$pattern == focus]
}

heights_vec <- c(
  get_n_rows("Systemic-like"),
  get_n_rows("Visceral-like"),
  get_n_rows("Neurologic"),
  get_n_rows("Ocular"),
  get_n_rows("Hematologic")
)

## 2. 拼接 5 张蝴蝶图 (按 TTO 升序: Systemic-like → Hematologic) ####
panel_grid <- (
  (p_systemic + theme(plot.title = element_text(size = 12))) /
  (p_visceral + theme(plot.title = element_text(size = 12))) /
  (p_neuro    + theme(plot.title = element_text(size = 12))) /
  (p_ocular   + theme(plot.title = element_text(size = 12))) /
  (p_hema     + theme(plot.title = element_text(size = 12)))
) +
  plot_layout(heights = heights_vec, guides = "collect")

## 3. 绘图 ####
p_combined <- panel_grid +
  plot_annotation(
    title = "Pattern-specific pathway enrichment landscape (5 patterns, butterfly view)",
    subtitle = "Bars filled with pathway database color; pattern labels and count circles use NMF irAE pattern color. Bar direction: positive Δρ → right; negative Δρ → left. Patterns ordered top-to-bottom by FAERS median TTO (Systemic-like 32d → Hematologic 43d).",
    theme = theme(
      plot.title    = element_text(size = 18, face = "bold", color = "black",
                                   family = "Arial"),
      plot.subtitle = element_text(size = 13, color = "black", family = "Arial",
                                   margin = margin(b = 10))
    )
  ) &
  theme(legend.position = "right")

print(p_combined)

## 4. 保存图形 ####
combined_pdf_path <- file.path(figure_output_path, "模式通路富集图_拼接版.pdf")
combined_height <- 30  ## 拼接版高度 (英寸),按需手动调节 ####
ggsave(
  filename = combined_pdf_path,
  plot = p_combined,
  width = 18, height = combined_height,
  device = cairo_pdf, limitsize = FALSE
)
log_message(sprintf("6图 拼接版 已保存: %s (h=%.1f)", combined_pdf_path, combined_height))


log_message(paste0("========== ", script_name, " 完成 =========="))
