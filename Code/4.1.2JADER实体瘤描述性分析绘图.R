# JADER实体瘤描述性分析绘图####

rm(list = ls())
gc()

source("Code/config.R")

script_name <- "4.1.2JADER实体瘤描述性分析绘图"
script_results_path <- file.path(results_path, script_name)
if (!dir.exists(script_results_path)) {
  dir.create(script_results_path, recursive = TRUE, showWarnings = FALSE)
}

message("========== JADER实体瘤描述性分析绘图开始 ==========")
message("结果目录: ", script_results_path)


## 1. 参数设置####

jader_input_path <- file.path(output_path, "jader_solid.RData")

summary_output_path <- file.path(
  script_results_path,
  "JADER实体瘤描述性分析汇总.xlsx"
)

object_output_path <- file.path(
  script_results_path,
  "JADER实体瘤描述性分析对象.RData"
)

annual_year_min <- 2015
annual_year_max <- NULL

gtex_mapping_path <- file.path(
  results_path,
  "7.1TCGA基础数据准备",
  "tables",
  "模式群-原发部位-TCGA项目-GTEx器官映射表（含样本量）.xlsx"
)


## 2. 辅助函数####

### 读取4.1脚本中的描述性分析通用辅助函数####
# 功能：复用4.1描述性分析脚本中的通用数据整理、配色和导出辅助函数，不执行FAERS主体分析####
# 输入：helper_script_path，4.1描述性分析脚本路径####
# 输出：在当前环境中注册4.1脚本“辅助函数”代码块内定义的函数####
load_descriptive_helper_functions <- function(helper_script_path) {
  if (!file.exists(helper_script_path)) {
    stop("未找到4.1描述性分析脚本: ", helper_script_path)
  }

  helper_lines <- readLines(helper_script_path, encoding = "UTF-8", warn = FALSE)
  helper_start <- grep("^## 辅助函数####$", helper_lines)
  helper_end <- grep("^## 1\\. 加载实体瘤分析集####$", helper_lines)

  if (length(helper_start) != 1 || length(helper_end) != 1 || helper_start >= helper_end) {
    stop("无法定位4.1脚本中的辅助函数代码块")
  }

  helper_code <- paste(helper_lines[helper_start:(helper_end - 1)], collapse = "\n")
  eval(parse(text = helper_code), envir = parent.frame())

  invisible(TRUE)
}

### 从JADER数据中提取报告年份####
# 功能：优先读取report_year字段，缺失时兼容quarter字段提取年份####
# 输入：data，JADER实体瘤分析集数据框####
# 输出：整数型报告年份向量####
extract_jader_report_year <- function(data) {
  if ("report_year" %in% names(data)) {
    return(suppressWarnings(as.integer(as.character(data$report_year))))
  }

  if ("quarter" %in% names(data)) {
    return(extract_report_year(data$quarter))
  }

  stop("JADER数据中未找到 report_year 或 quarter 字段，无法生成年度趋势")
}

### 从JADER数据中提取报告季度####
# 功能：优先使用report_year、report_q和year_quarter字段生成季度标签，缺失时从report_quarter文本解析####
# 输入：data，JADER实体瘤分析集数据框####
# 输出：包含report_period_plot和report_period_index的数据框####
extract_jader_report_quarter_data <- function(data) {
  if (all(c("report_year", "report_q") %in% names(data))) {
    report_year <- suppressWarnings(as.integer(as.character(data$report_year)))
    report_q <- suppressWarnings(as.integer(as.character(data$report_q)))

    if ("year_quarter" %in% names(data)) {
      report_period_index <- suppressWarnings(as.integer(as.character(data$year_quarter)))
    } else {
      report_period_index <- report_year * 10L + report_q
    }

    return(
      tibble(
        report_period_plot = if_else(
          !is.na(report_year) & !is.na(report_q),
          paste0(report_year, "Q", report_q),
          NA_character_
        ),
        report_period_index = report_period_index
      )
    )
  }

  if ("report_quarter" %in% names(data)) {
    report_quarter <- as.character(data$report_quarter)
    report_year <- suppressWarnings(as.integer(str_extract(report_quarter, "\\d{4}")))
    report_q <- case_when(
      str_detect(report_quarter, "第一|Q1|1Q") ~ 1L,
      str_detect(report_quarter, "第二|Q2|2Q") ~ 2L,
      str_detect(report_quarter, "第三|Q3|3Q") ~ 3L,
      str_detect(report_quarter, "第四|Q4|4Q") ~ 4L,
      TRUE ~ NA_integer_
    )

    return(
      tibble(
        report_period_plot = if_else(
          !is.na(report_year) & !is.na(report_q),
          paste0(report_year, "Q", report_q),
          NA_character_
        ),
        report_period_index = report_year * 10L + report_q
      )
    )
  }

  if ("quarter" %in% names(data)) {
    quarter_text <- normalize_report_quarter(data$quarter)

    return(
      tibble(
        report_period_plot = quarter_text,
        report_period_index = extract_report_quarter_index(quarter_text)
      )
    )
  }

  stop("JADER数据中未找到季度字段，无法生成季度粒度SOC景观图")
}

### 构建JADER年度实体瘤报告与ICI使用趋势数据####
# 功能：按报告年份汇总JADER实体瘤总报告数和ICI治疗策略分层报告数####
# 输入：data，JADER实体瘤分析集；year_min，最小展示年份；year_max，最大展示年份####
# 输出：包含展示年份、实体瘤总数、ICI分层数和ICI总数的数据列表####
build_jader_annual_trend_data <- function(data, year_min = 2015, year_max = NULL) {
  data_with_year <- data %>%
    mutate(year = extract_jader_report_year(data))

  observed_years <- data_with_year %>%
    filter(!is.na(year)) %>%
    pull(year) %>%
    unique() %>%
    sort()

  if (length(observed_years) == 0) {
    stop("JADER实体瘤数据中不存在可用report_year，无法生成年度趋势")
  }

  max_available_year <- max(observed_years, na.rm = TRUE)
  max_display_year <- if (is.null(year_max)) max_available_year else min(year_max, max_available_year)

  if (max_display_year < year_min) {
    available_years <- observed_years[observed_years <= max_display_year]
  } else {
    available_years <- seq(year_min, max_display_year)
  }

  if (!is.null(year_max)) {
    available_years <- available_years[available_years <= year_max]
  }

  if (length(available_years) == 0) {
    stop("年度趋势展示年份为空，请检查annual_year_min和annual_year_max")
  }

  strategy_levels <- format_therapy_strategy_bar_label(
    c("anti-PD-L1", "anti-PD-1", "anti-CTLA-4", "Combination Therapy")
  )

  solid_yearly_total <- data_with_year %>%
    filter(year %in% available_years) %>%
    count(year, name = "solid_total") %>%
    tidyr::complete(
      year = available_years,
      fill = list(solid_total = 0)
    ) %>%
    mutate(x_pos = year - 0.18)

  ici_yearly_by_strategy <- data_with_year %>%
    filter(year %in% available_years, has_ici == TRUE, treatment_strategy != "No ICI") %>%
    mutate(therapy_strategy = format_therapy_strategy_bar_label(treatment_strategy)) %>%
    count(year, therapy_strategy, name = "n") %>%
    tidyr::complete(
      year = available_years,
      therapy_strategy = strategy_levels,
      fill = list(n = 0)
    ) %>%
    mutate(
      therapy_strategy = factor(therapy_strategy, levels = strategy_levels),
      x_pos = year + 0.18
    )

  ici_yearly_total <- ici_yearly_by_strategy %>%
    group_by(year) %>%
    summarise(ici_total = sum(n), .groups = "drop") %>%
    mutate(x_pos = year + 0.18)

  list(
    available_years = available_years,
    solid_yearly_total = solid_yearly_total,
    ici_yearly_by_strategy = ici_yearly_by_strategy,
    ici_yearly_total = ici_yearly_total
  )
}

### 构建JADER实体瘤ICI报告SOC景观图数据####
# 功能：按癌种层级、报告季度和SOC汇总JADER ICI报告比例，并生成景观图注释配色####
# 输入：data，JADER实体瘤分析集；gtex_mapping_path，GTEx映射表路径####
# 输出：包含列元数据、SOC长表、SOC汇总表和注释配色的数据列表####
build_jader_ici_soc_landscape_data <- function(data, gtex_mapping_path = NULL) {
  soc_metadata <- extract_soc_metadata_simple(data)

  jader_quarter_data <- extract_jader_report_quarter_data(data)
  data_with_quarter <- bind_cols(data, jader_quarter_data)

  ici_data <- data_with_quarter %>%
    filter(has_ici == TRUE) %>%
    mutate(
      primary_organ_plot = standardize_cancer_display_label(primary_organ),
      histologic_category_plot = standardize_cancer_display_label(
        histologic_category,
        unknown_patterns = c("^other$", "^others$", "^unknown/others$")
      ),
      assigned_type_plot = standardize_cancer_display_label(
        assigned_type,
        unknown_patterns = c("^other$", "^others$", "^na$", "^n/a$", "^unknown/others$", "^unknown$")
      )
    ) %>%
    filter(!is.na(report_period_index))

  if (nrow(ici_data) == 0) {
    stop("JADER ICI病例中不存在可用报告季度，无法绘制SOC景观图")
  }

  primary_order <- ici_data %>%
    count(primary_organ_plot, name = "primary_total")

  histologic_order <- ici_data %>%
    count(primary_organ_plot, histologic_category_plot, name = "histologic_total") %>%
    left_join(primary_order, by = "primary_organ_plot")

  tcga_order <- ici_data %>%
    count(primary_organ_plot, histologic_category_plot, assigned_type_plot, name = "tcga_total") %>%
    left_join(primary_order, by = "primary_organ_plot") %>%
    left_join(
      histologic_order %>% select(primary_organ_plot, histologic_category_plot, histologic_total),
      by = c("primary_organ_plot", "histologic_category_plot")
    )

  column_meta <- ici_data %>%
    count(
      primary_organ_plot,
      histologic_category_plot,
      assigned_type_plot,
      report_period_plot,
      report_period_index,
      name = "report_cases"
    ) %>%
    left_join(primary_order, by = "primary_organ_plot") %>%
    left_join(
      histologic_order %>% select(primary_organ_plot, histologic_category_plot, histologic_total),
      by = c("primary_organ_plot", "histologic_category_plot")
    ) %>%
    left_join(
      tcga_order %>% select(primary_organ_plot, histologic_category_plot, assigned_type_plot, tcga_total),
      by = c("primary_organ_plot", "histologic_category_plot", "assigned_type_plot")
    ) %>%
    arrange(
      desc(primary_total),
      primary_organ_plot,
      desc(histologic_total),
      histologic_category_plot,
      desc(tcga_total),
      assigned_type_plot,
      report_period_index,
      report_period_plot
    ) %>%
    mutate(column_id = row_number())

  gtex_mapping_df <- NULL
  if (!is.null(gtex_mapping_path)) {
    gtex_mapping_df <- read_gtex_landscape_mapping(gtex_mapping_path)
    column_meta <- add_gtex_mapping_to_landscape_column_meta(
      column_meta = column_meta,
      gtex_mapping_df = gtex_mapping_df
    )
  }

  soc_long <- ici_data %>%
    group_by(
      primary_organ_plot,
      histologic_category_plot,
      assigned_type_plot,
      report_period_plot,
      report_period_index
    ) %>%
    summarise(
      report_cases = n(),
      across(all_of(soc_metadata$soc_var), ~ sum(.x == 1, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    pivot_longer(
      cols = all_of(soc_metadata$soc_var),
      names_to = "soc_var",
      values_to = "soc_reports"
    ) %>%
    left_join(soc_metadata, by = "soc_var") %>%
    left_join(
      column_meta %>%
        select(
          primary_organ_plot,
          histologic_category_plot,
          assigned_type_plot,
          report_period_plot,
          report_period_index,
          column_id
        ),
      by = c(
        "primary_organ_plot",
        "histologic_category_plot",
        "assigned_type_plot",
        "report_period_plot",
        "report_period_index"
      )
    ) %>%
    mutate(report_pct = if_else(report_cases > 0, soc_reports / report_cases * 100, 0))

  soc_summary <- soc_long %>%
    group_by(SOC_full) %>%
    summarise(total_reports = sum(soc_reports), .groups = "drop") %>%
    arrange(desc(total_reports), SOC_full)

  soc_abbrev_map <- build_soc_abbrev_map()
  soc_levels <- soc_summary$SOC_full
  soc_display_levels <- format_soc_display_label(soc_levels, soc_abbrev_map)

  soc_long <- soc_long %>%
    mutate(
      SOC_full = factor(SOC_full, levels = rev(soc_levels)),
      SOC_display = factor(format_soc_display_label(SOC_full, soc_abbrev_map), levels = rev(soc_display_levels))
    )

  soc_summary <- soc_summary %>%
    mutate(
      SOC_full = factor(SOC_full, levels = rev(soc_levels)),
      SOC_display = factor(format_soc_display_label(SOC_full, soc_abbrev_map), levels = rev(soc_display_levels))
    )

  cancer_color_results <- build_hierarchical_cancer_color_maps(column_meta)
  column_meta <- cancer_color_results$column_meta
  primary_color_map <- cancer_color_results$primary_color_map
  histologic_color_map <- cancer_color_results$histologic_color_map
  tcga_color_map <- cancer_color_results$tcga_color_map
  color_mapping <- cancer_color_results$color_mapping
  gtex_color_map <- build_gtex_tissue_color_map(
    column_meta = column_meta,
    primary_color_map = primary_color_map
  )
  gtex_color_mapping <- tibble(
    level_name = "GTEx",
    parent_category = NA_character_,
    color_key = names(gtex_color_map),
    category = names(gtex_color_map),
    color = unname(gtex_color_map)
  )
  color_mapping <- bind_rows(color_mapping, gtex_color_mapping)

  heatmap_max_pct <- max(soc_long$report_pct, na.rm = TRUE)
  heatmap_max_pct <- if (is.finite(heatmap_max_pct) && heatmap_max_pct > 0) heatmap_max_pct else 1

  list(
    ici_data = ici_data,
    column_meta = column_meta,
    soc_long = soc_long,
    soc_summary = soc_summary,
    soc_metadata = soc_metadata,
    soc_abbrev_map = soc_abbrev_map,
    gtex_mapping = gtex_mapping_df,
    gtex_color_map = gtex_color_map,
    gtex_color_mapping = gtex_color_mapping,
    primary_color_map = primary_color_map,
    histologic_color_map = histologic_color_map,
    tcga_color_map = tcga_color_map,
    color_mapping = color_mapping,
    heatmap_max_pct = heatmap_max_pct
  )
}

### 读取描述性分析结果Excel指定sheet####
# 功能：从已保存的JADER描述性分析汇总Excel中读取指定sheet####
# 输入：workbook_path，Excel路径；sheet_name，sheet名称####
# 输出：tibble格式数据框####
read_result_sheet <- function(workbook_path, sheet_name) {
  if (!file.exists(workbook_path)) {
    stop("未找到结果Excel: ", workbook_path)
  }

  readxl::read_xlsx(
    path = workbook_path,
    sheet = sheet_name
  ) %>%
    tibble::as_tibble()
}

## 3. 统计分析与结果导出####

load_descriptive_helper_functions("Code/4.1描述性分析.R")

jader_solid <- load_named_rdata_object(jader_input_path, "jader_solid")
jader_solid_tumor_final <- jader_solid

filter_summary <- tibble(
  step = "jader_solid",
  n = nrow(jader_solid_tumor_final),
  retained_pct = 100
)

message("已加载 jader_solid")
message("总样本数: ", nrow(jader_solid_tumor_final))
message("变量数: ", ncol(jader_solid_tumor_final))

primary_organ_summary <- build_distribution_summary(jader_solid_tumor_final, primary_organ)
histologic_summary <- build_distribution_summary(jader_solid_tumor_final, histologic_category)
assigned_type_summary <- build_distribution_summary(jader_solid_tumor_final, assigned_type)

ici_sankey_results <- build_ici_sankey_source_data(jader_solid_tumor_final)
ici_type_summary <- ici_sankey_results$ici_type_counts
ici_agent_summary <- ici_sankey_results$agent_counts
ici_regimen_summary <- ici_sankey_results$regimen_counts
ici_strategy_summary <- ici_sankey_results$strategy_counts
ici_sankey_flow_summary <- ici_sankey_results$sankey_data %>%
  count(ici_class, drug_name, ici_regimen, therapy_strategy, name = "n", sort = TRUE)
ici_sankey_plot_data <- ici_sankey_results$sankey_long
ici_reference_color_map <- build_ici_reference_color_map()
ici_sankey_node_color_map <- build_ici_sankey_node_color_map(
  sankey_data = ici_sankey_results$sankey_data,
  reference_color_map = ici_reference_color_map
)
ici_sankey_node_color_table <- tibble(
  node = names(ici_sankey_node_color_map),
  color = unname(ici_sankey_node_color_map)
)

annual_trend_results <- build_jader_annual_trend_data(
  jader_solid_tumor_final,
  year_min = annual_year_min,
  year_max = annual_year_max
)
annual_trend_years <- annual_trend_results$available_years
solid_yearly_total <- annual_trend_results$solid_yearly_total
ici_yearly_by_strategy <- annual_trend_results$ici_yearly_by_strategy
ici_yearly_total <- annual_trend_results$ici_yearly_total

ici_soc_landscape_results <- build_jader_ici_soc_landscape_data(
  jader_solid_tumor_final,
  gtex_mapping_path = gtex_mapping_path
)
jader_solid_ici_landscape <- ici_soc_landscape_results$ici_data
ici_landscape_column_meta <- ici_soc_landscape_results$column_meta
ici_landscape_soc_long <- ici_soc_landscape_results$soc_long
ici_landscape_soc_summary <- ici_soc_landscape_results$soc_summary
ici_landscape_soc_abbrev_map <- ici_soc_landscape_results$soc_abbrev_map
ici_landscape_gtex_mapping <- ici_soc_landscape_results$gtex_mapping
ici_landscape_gtex_color_mapping <- ici_soc_landscape_results$gtex_color_mapping
ici_landscape_color_mapping <- ici_soc_landscape_results$color_mapping
ici_landscape_primary_color_map <- ici_soc_landscape_results$primary_color_map
ici_landscape_histologic_color_map <- ici_soc_landscape_results$histologic_color_map
ici_landscape_tcga_color_map <- ici_soc_landscape_results$tcga_color_map
ici_landscape_heatmap_max_pct <- ici_soc_landscape_results$heatmap_max_pct

cancer_mapping_path <- file.path(results_path, "癌种映射.xlsx")
valid_organ_histology <- build_valid_organ_histology_mapping(cancer_mapping_path)
jader_hierarchy_export <- build_three_layer_hierarchy_export(
  all_data = jader_solid_tumor_final,
  ici_data = jader_solid_tumor_final %>% filter(has_ici == TRUE),
  valid_organ_histology = valid_organ_histology
)

jader_hierarchy_export_path <- file.path(script_results_path, "JADER图3_三层结构.xlsx")
writexl::write_xlsx(
  jader_hierarchy_export,
  path = jader_hierarchy_export_path
)

if (is.null(ici_landscape_gtex_mapping)) {
  ici_landscape_gtex_mapping <- tibble()
}

writexl::write_xlsx(
  list(
    "筛选流程" = filter_summary,
    "原发器官分布" = primary_organ_summary,
    "组织学分类分布" = histologic_summary,
    "癌种实体分布" = assigned_type_summary,
    "ICI类别分布" = ici_type_summary,
    "ICI药物分布" = ici_agent_summary,
    "ICI方案分布" = ici_regimen_summary,
    "治疗策略分布" = ici_strategy_summary,
    "ICI桑基图长表" = ici_sankey_plot_data,
    "ICI桑基图节点配色" = ici_sankey_node_color_table,
    "ICI桑基图路径统计" = ici_sankey_flow_summary,
    "年度年份" = tibble(year = annual_trend_years),
    "年度实体瘤总数" = solid_yearly_total,
    "年度ICI_按策略堆叠" = ici_yearly_by_strategy,
    "年度ICI总数" = ici_yearly_total,
    "第三图_列元数据" = ici_landscape_column_meta,
    "第三图_SOC百分比长表" = ici_landscape_soc_long,
    "第三图_SOC总报告数" = ici_landscape_soc_summary,
    "第三图_GTEx映射" = ici_landscape_gtex_mapping,
    "第三图_颜色映射" = ici_landscape_color_mapping,
    "第三图_GTEx颜色映射" = ici_landscape_gtex_color_mapping,
    "图3_三层结构" = jader_hierarchy_export
  ),
  path = summary_output_path
)

message("JADER描述性分析汇总Excel已保存: ", summary_output_path)
message("JADER图3三层结构已保存: ", jader_hierarchy_export_path)


# 1图：JADER solid tumor ICI therapy Sankey diagram ####

### 1. 读取并处理绘图数据####

jader_sankey_plot_data <- read_result_sheet(
  summary_output_path,
  "ICI桑基图长表"
) %>%
  mutate(
    x_new = as.numeric(x_new),
    next_x_new = as.numeric(next_x_new)
  )

jader_sankey_node_color_table <- read_result_sheet(
  summary_output_path,
  "ICI桑基图节点配色"
)
jader_sankey_node_color_map <- stats::setNames(
  jader_sankey_node_color_table$color,
  jader_sankey_node_color_table$node
)
jader_sankey_x_breaks <- c(0.9, 1.35, 1.95, 2.9)

### 2. 绘图####

p_jader_solid_sankey <- ggplot(
  jader_sankey_plot_data,
  aes(
    x = x_new,
    next_x = next_x_new,
    node = node,
    next_node = next_node,
    fill = factor(node),
    label = node
  )
) +
  geom_sankey(
    flow.alpha = 0.5,
    node.color = "gray30",
    width = 0.16,
    smooth = 8
  ) +
  geom_sankey_label(
    size = 14 / ggplot2::.pt,
    color = "black",
    fill = "white",
    family = "sans",
    hjust = 0.5
  ) +
  scale_fill_manual(values = jader_sankey_node_color_map, guide = "none") +
  scale_x_continuous(
    breaks = jader_sankey_x_breaks,
    labels = c("ICI class", "Agent", "ICI regimen", "Therapy strategy"),
    expand = expansion(mult = c(0.05, 0.17))
  ) +
  coord_cartesian(clip = "off") +
  theme_sankey(base_size = 14) +
  theme(
    legend.position = "none",
    text = element_text(size = 14, color = "black", family = "sans"),
    axis.title = element_text(size = 15, color = "black", family = "sans"),
    axis.text = element_text(size = 15, color = "black", family = "sans"),
    axis.text.x = element_text(size = 15, color = "black", family = "sans", margin = margin(t = -6)),
    plot.margin = margin(2, 18, 2, 8)
  ) +
  labs(x = NULL, y = NULL)

print(p_jader_solid_sankey)

### 3. 保存图形####

ggsave(
  filename = file.path(script_results_path, "JADER实体瘤_ICI情况桑基图.pdf"),
  plot = p_jader_solid_sankey,
  width = 11,
  height = 5,
  device = cairo_pdf,
  family = "sans"
)


# 2图：JADER annual solid tumor reports and ICI trend ####

### 1. 读取并处理绘图数据####

jader_annual_years <- read_result_sheet(
  summary_output_path,
  "年度年份"
) %>%
  pull(year) %>%
  as.integer()

jader_solid_yearly_total <- read_result_sheet(
  summary_output_path,
  "年度实体瘤总数"
) %>%
  mutate(
    year = as.integer(year),
    solid_total = as.numeric(solid_total),
    x_pos = as.numeric(x_pos)
  )

jader_ici_yearly_by_strategy <- read_result_sheet(
  summary_output_path,
  "年度ICI_按策略堆叠"
) %>%
  mutate(
    year = as.integer(year),
    n = as.numeric(n),
    x_pos = as.numeric(x_pos)
  )

jader_ici_yearly_total <- read_result_sheet(
  summary_output_path,
  "年度ICI总数"
) %>%
  mutate(
    year = as.integer(year),
    ici_total = as.numeric(ici_total),
    x_pos = as.numeric(x_pos)
  )

jader_annual_strategy_levels <- format_therapy_strategy_bar_label(
  c("anti-PD-L1", "anti-PD-1", "anti-CTLA-4", "Combination Therapy")
)
jader_ici_yearly_by_strategy <- jader_ici_yearly_by_strategy %>%
  mutate(
    therapy_strategy = factor(as.character(therapy_strategy), levels = jader_annual_strategy_levels)
  )
jader_annual_strategy_color_map <- build_ici_reference_color_map()[jader_annual_strategy_levels]
jader_annual_label_offset <- max(
  c(jader_solid_yearly_total$solid_total, jader_ici_yearly_total$ici_total),
  na.rm = TRUE
) * 0.03
jader_annual_label_offset <- ifelse(
  is.finite(jader_annual_label_offset) && jader_annual_label_offset > 0,
  jader_annual_label_offset,
  1
)

### 2. 绘图####

p_jader_solid_annual_trend_bar <- ggplot() +
  geom_col(
    data = jader_solid_yearly_total,
    aes(x = x_pos, y = solid_total),
    width = 0.32,
    fill = "#1F3B73",
    show.legend = FALSE
  ) +
  geom_text(
    data = jader_solid_yearly_total,
    aes(x = x_pos, y = solid_total, label = format_n(solid_total)),
    vjust = -0.3,
    size = 12 / ggplot2::.pt,
    color = "black"
  ) +
  geom_col(
    data = jader_ici_yearly_by_strategy,
    aes(x = x_pos, y = n, fill = therapy_strategy),
    width = 0.32
  ) +
  geom_text(
    data = jader_ici_yearly_total,
    aes(
      x = x_pos - 0.05,
      y = ici_total + jader_annual_label_offset,
      label = format_n(ici_total)
    ),
    vjust = 0,
    hjust = 0,
    angle = 45,
    size = 12 / ggplot2::.pt,
    color = "black"
  ) +
  scale_fill_manual(
    values = jader_annual_strategy_color_map,
    name = "Therapy strategy",
    drop = FALSE,
    guide = guide_legend(order = 2, nrow = 2, byrow = TRUE)
  ) +
  scale_x_continuous(
    breaks = jader_annual_years,
    labels = jader_annual_years,
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  scale_y_continuous(
    labels = scales::comma,
    expand = expansion(mult = c(0, 0.15))
  ) +
  labs(x = NULL, y = "Number of cases") +
  theme_classic(base_size = 14) +
  theme(
    text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 15, color = "black"),
    axis.text = element_text(size = 14, color = "black"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.margin = margin(t = -6, r = 0, b = 0, l = 0),
    legend.box.margin = margin(t = -8, r = 0, b = 0, l = 0),
    legend.title = element_text(size = 15, color = "black"),
    legend.text = element_text(size = 14, color = "black"),
    plot.margin = margin(6, 12, 2, 6)
  )

print(p_jader_solid_annual_trend_bar)

### 3. 保存图形####

ggsave(
  filename = file.path(script_results_path, "JADER实体瘤_年度报告与ICI趋势图.pdf"),
  plot = p_jader_solid_annual_trend_bar,
  width = 11,
  height = 3,
  device = cairo_pdf
)


# 3图：JADER ICI report SOC landscape ####

### 1. 读取并处理绘图数据####

jader_landscape_column_meta <- read_result_sheet(
  summary_output_path,
  "第三图_列元数据"
) %>%
  mutate(
    column_id = as.integer(column_id),
    report_cases = as.numeric(report_cases),
    primary_total = as.numeric(primary_total),
    report_period_plot = as.character(report_period_plot),
    histologic_color_key = as.character(histologic_color_key),
    tcga_color_key = as.character(tcga_color_key),
    gtex_tissue_plot = as.character(gtex_tissue_plot)
  )

jader_landscape_soc_long <- read_result_sheet(
  summary_output_path,
  "第三图_SOC百分比长表"
) %>%
  mutate(
    column_id = as.integer(column_id),
    report_cases = as.numeric(report_cases),
    soc_reports = as.numeric(soc_reports),
    report_pct = as.numeric(report_pct)
  )

jader_landscape_soc_summary <- read_result_sheet(
  summary_output_path,
  "第三图_SOC总报告数"
) %>%
  mutate(total_reports = as.numeric(total_reports))

jader_landscape_color_mapping <- read_result_sheet(
  summary_output_path,
  "第三图_颜色映射"
) %>%
  mutate(
    level_name = as.character(level_name),
    color_key = as.character(color_key),
    color = as.character(color)
  ) %>%
  distinct(level_name, color_key, .keep_all = TRUE)

jader_landscape_soc_abbrev_map <- build_soc_abbrev_map()
jader_landscape_soc_summary <- ensure_soc_display_column(
  jader_landscape_soc_summary,
  soc_abbrev_map = jader_landscape_soc_abbrev_map
)
jader_landscape_soc_long <- ensure_soc_display_column(
  jader_landscape_soc_long,
  soc_abbrev_map = jader_landscape_soc_abbrev_map
)

jader_landscape_primary_color_map <- jader_landscape_color_mapping %>%
  filter(level_name == "Primary organ") %>%
  {stats::setNames(.$color, .$color_key)}
jader_landscape_histologic_color_map <- jader_landscape_color_mapping %>%
  filter(level_name == "Histologic category") %>%
  {stats::setNames(.$color, .$color_key)}
jader_landscape_tcga_color_map <- jader_landscape_color_mapping %>%
  filter(level_name == "TCGA") %>%
  {stats::setNames(.$color, .$color_key)}
jader_landscape_gtex_color_map <- jader_landscape_color_mapping %>%
  filter(level_name == "GTEx") %>%
  {stats::setNames(.$color, .$color_key)}

jader_landscape_primary_annotation <- factor(
  jader_landscape_column_meta$primary_organ_plot,
  levels = unique(jader_landscape_column_meta$primary_organ_plot)
)
jader_landscape_histologic_annotation <- factor(
  jader_landscape_column_meta$histologic_color_key,
  levels = unique(jader_landscape_column_meta$histologic_color_key)
)
jader_landscape_tcga_annotation <- factor(
  jader_landscape_column_meta$tcga_color_key,
  levels = unique(jader_landscape_column_meta$tcga_color_key)
)
jader_landscape_gtex_annotation <- factor(
  jader_landscape_column_meta$gtex_tissue_plot,
  levels = unique(jader_landscape_column_meta$gtex_tissue_plot)
)
jader_landscape_column_order <- as.character(jader_landscape_column_meta$column_id)

jader_landscape_primary_legend_data <- jader_landscape_column_meta %>%
  distinct(primary_organ_plot, primary_total) %>%
  arrange(desc(primary_total), primary_organ_plot) %>%
  mutate(
    primary_label = if_else(
      primary_organ_plot == "Unknown/Others",
      "Unknown/Others",
      paste0(primary_organ_plot, "(", format_n(primary_total), ")")
    ),
    primary_color = unname(jader_landscape_primary_color_map[primary_organ_plot])
  )

if (!("Unknown/Others" %in% jader_landscape_primary_legend_data$primary_organ_plot)) {
  jader_landscape_primary_legend_data <- bind_rows(
    jader_landscape_primary_legend_data,
    tibble(
      primary_organ_plot = "Unknown/Others",
      primary_total = NA_real_,
      primary_label = "Unknown/Others",
      primary_color = "#BDBDBD"
    )
  )
}

jader_landscape_primary_legend_ncol <- 5L

jader_landscape_soc_summary_plot <- jader_landscape_soc_summary %>%
  mutate(
    SOC_full = as.character(SOC_full),
    SOC_display = as.character(SOC_display),
    SOC_plot_label = format_soc_abbrev_label(SOC_full, jader_landscape_soc_abbrev_map)
  ) %>%
  arrange(desc(total_reports), SOC_full)
jader_landscape_soc_order <- jader_landscape_soc_summary_plot$SOC_display
jader_landscape_soc_plot_order <- jader_landscape_soc_summary_plot$SOC_plot_label

jader_landscape_heatmap_matrix <- jader_landscape_soc_long %>%
  mutate(
    SOC_display = as.character(SOC_display),
    column_id = as.character(column_id)
  ) %>%
  select(SOC_display, column_id, report_pct) %>%
  tidyr::pivot_wider(
    names_from = column_id,
    values_from = report_pct,
    values_fill = 0
  ) %>%
  as.data.frame()

rownames(jader_landscape_heatmap_matrix) <- jader_landscape_heatmap_matrix$SOC_display
jader_landscape_heatmap_matrix$SOC_display <- NULL
jader_landscape_heatmap_matrix <- as.matrix(
  jader_landscape_heatmap_matrix[
    jader_landscape_soc_order,
    jader_landscape_column_order,
    drop = FALSE
  ]
)
rownames(jader_landscape_heatmap_matrix) <- jader_landscape_soc_plot_order
storage.mode(jader_landscape_heatmap_matrix) <- "double"

jader_landscape_report_cases_vector <- jader_landscape_column_meta$report_cases
names(jader_landscape_report_cases_vector) <- jader_landscape_column_order

jader_landscape_soc_total_vector <- jader_landscape_soc_summary_plot$total_reports
names(jader_landscape_soc_total_vector) <- jader_landscape_soc_plot_order

jader_landscape_heatmap_max_pct <- compute_heatmap_upper_bound(jader_landscape_heatmap_matrix)
jader_landscape_heatmap_breaks <- build_heatmap_color_breaks(
  mat = jader_landscape_heatmap_matrix,
  upper_bound = jader_landscape_heatmap_max_pct
)
jader_landscape_heatmap_legend_breaks <- c(0, 25, 50, 75, 100)
jader_landscape_heatmap_legend_labels <- as.character(jader_landscape_heatmap_legend_breaks)

jader_landscape_soc_bar_breaks <- pretty(
  c(0, max(jader_landscape_soc_total_vector, na.rm = TRUE)),
  n = 4
)
jader_landscape_soc_bar_breaks <- unique(
  jader_landscape_soc_bar_breaks[
    jader_landscape_soc_bar_breaks >= 0 &
      jader_landscape_soc_bar_breaks <= max(jader_landscape_soc_total_vector, na.rm = TRUE)
  ]
)
if (length(jader_landscape_soc_bar_breaks) < 2) {
  jader_landscape_soc_bar_breaks <- unique(c(0, max(jader_landscape_soc_total_vector, na.rm = TRUE)))
}

jader_landscape_heatmap_color_fun <- circlize::colorRamp2(
  jader_landscape_heatmap_breaks,
  grDevices::colorRampPalette(
    c("#FFFFFF", "#EAEAEA", "#D0D0D0", "#A8A8A8", "#737373", "#303030", "#000000")
  )(length(jader_landscape_heatmap_breaks))
)
jader_landscape_row_names_max_width <- ComplexHeatmap::max_text_width(
  jader_landscape_soc_plot_order,
  gp = grid::gpar(fontsize = 13, col = "black", fontfamily = "sans")
) + grid::unit(4, "mm")

### 2. 绘图####

p_jader_solid_ici_landscape_top_annotation <- ComplexHeatmap::HeatmapAnnotation(
  `Report cases per quarter` = ComplexHeatmap::anno_barplot(
    jader_landscape_report_cases_vector,
    gp = grid::gpar(fill = "#394B59", col = NA),
    border = TRUE,
    baseline = 0,
    bar_width = 1,
    height = grid::unit(24, "mm"),
    axis_param = list(
      side = "left",
      gp = grid::gpar(fontsize = 13, col = "black", fontfamily = "sans"),
      labels_rot = 0
    )
  ),
  `Primary organ` = jader_landscape_primary_annotation,
  `Histologic category` = jader_landscape_histologic_annotation,
  `TCGA` = jader_landscape_tcga_annotation,
  `GTEx` = jader_landscape_gtex_annotation,
  col = list(
    `Primary organ` = jader_landscape_primary_color_map,
    `Histologic category` = jader_landscape_histologic_color_map,
    `TCGA` = jader_landscape_tcga_color_map,
    `GTEx` = jader_landscape_gtex_color_map
  ),
  show_legend = FALSE,
  show_annotation_name = TRUE,
  annotation_name_side = "right",
  annotation_name_rot = 0,
  annotation_name_gp = grid::gpar(fontsize = 15, col = "black", fontfamily = "sans"),
  annotation_height = grid::unit.c(
    grid::unit(24, "mm"),
    grid::unit(4.2, "mm"),
    grid::unit(4.2, "mm"),
    grid::unit(4.2, "mm"),
    grid::unit(4.2, "mm")
  ),
  gap = grid::unit(c(2, 0.6, 0.6, 0.6), "mm"),
  gp = grid::gpar(col = NA)
)

p_jader_solid_ici_landscape_soc_bar <- ComplexHeatmap::rowAnnotation(
  `SOC reports` = ComplexHeatmap::anno_barplot(
    jader_landscape_soc_total_vector,
    gp = grid::gpar(fill = "#5C6969", col = NA),
    border = TRUE,
    baseline = 0,
    bar_width = 0.56,
    width = grid::unit(32, "mm"),
    axis_param = list(
      side = "bottom",
      at = jader_landscape_soc_bar_breaks,
      labels = scales::comma(jader_landscape_soc_bar_breaks),
      gp = grid::gpar(fontsize = 13, col = "black", fontfamily = "sans"),
      labels_rot = 90
    )
  ),
  show_annotation_name = FALSE,
  annotation_name_side = "bottom",
  annotation_name_rot = 0,
  annotation_name_gp = grid::gpar(fontsize = 15, col = "black", fontfamily = "sans")
)

p_jader_solid_ici_landscape_heatmap <- ComplexHeatmap::Heatmap(
  jader_landscape_heatmap_matrix,
  name = "% of report cases",
  col = jader_landscape_heatmap_color_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  row_names_side = "right",
  row_names_gp = grid::gpar(fontsize = 13, col = "black", fontfamily = "sans"),
  row_names_max_width = jader_landscape_row_names_max_width,
  row_title = NULL,
  column_title = NULL,
  height = grid::unit(58, "mm"),
  top_annotation = p_jader_solid_ici_landscape_top_annotation,
  right_annotation = p_jader_solid_ici_landscape_soc_bar,
  show_heatmap_legend = FALSE,
  rect_gp = grid::gpar(col = "#F2F2F2", lwd = 0.25),
  border = "#4D4D4D",
  use_raster = TRUE,
  raster_quality = 3
)

p_jader_solid_ici_landscape_heatmap_legend <- ComplexHeatmap::Legend(
  title = "% of report cases",
  col_fun = jader_landscape_heatmap_color_fun,
  at = jader_landscape_heatmap_legend_breaks,
  labels = jader_landscape_heatmap_legend_labels,
  direction = "horizontal",
  legend_width = grid::unit(42, "mm"),
  grid_height = grid::unit(4, "mm"),
  title_position = "topcenter",
  title_gp = grid::gpar(fontsize = 15, col = "black", fontfamily = "sans"),
  labels_gp = grid::gpar(fontsize = 13, col = "black", fontfamily = "sans")
)

p_jader_solid_ici_landscape_primary_legend <- ComplexHeatmap::Legend(
  title = "Primary organ",
  labels = jader_landscape_primary_legend_data$primary_label,
  legend_gp = grid::gpar(fill = jader_landscape_primary_legend_data$primary_color, col = NA),
  labels_gp = grid::gpar(fontsize = 12, col = "black", fontfamily = "sans"),
  title_gp = grid::gpar(fontsize = 15, col = "black", fontface = "bold", fontfamily = "sans"),
  ncol = jader_landscape_primary_legend_ncol,
  by_row = TRUE,
  gap = grid::unit(2.5, "mm"),
  row_gap = grid::unit(1.5, "mm"),
  grid_width = grid::unit(4.2, "mm"),
  grid_height = grid::unit(4.2, "mm")
)

jader_landscape_main_viewport_y <- grid::unit(0.62, "npc")
jader_landscape_main_viewport_height <- grid::unit(0.60, "npc")
jader_landscape_heatmap_padding <- grid::unit(c(3, 14, 6, 3), "mm")
jader_landscape_heatmap_legend_y <- grid::unit(0.45, "npc")
jader_landscape_heatmap_legend_height <- grid::unit(46, "mm")
jader_landscape_soc_bar_title_x <- grid::unit(0.88, "npc")
jader_landscape_soc_bar_title_y <- grid::unit(0.29, "npc")
jader_landscape_primary_legend_y <- grid::unit(0.04, "npc")

p_jader_solid_ici_landscape_heatmap_grob <- grid::grid.grabExpr(
  {
    ComplexHeatmap::draw(
      p_jader_solid_ici_landscape_heatmap,
      heatmap_legend_side = "left",
      annotation_legend_side = "bottom",
      merge_legends = FALSE,
      padding = jader_landscape_heatmap_padding
    )
  },
  width = 10,
  height = 5.4
)

p_jader_solid_ici_landscape_heatmap_legend_grob <- grid::grid.grabExpr(
  {
    ComplexHeatmap::draw(
      p_jader_solid_ici_landscape_heatmap_legend,
      x = grid::unit(0.5, "npc"),
      y = grid::unit(0.5, "npc"),
      just = c("center", "center")
    )
  },
  width = 1.4,
  height = 3
)

p_jader_solid_ici_landscape_primary_legend_grob <- grid::grid.grabExpr(
  {
    ComplexHeatmap::draw(
      p_jader_solid_ici_landscape_primary_legend,
      x = grid::unit(0, "npc"),
      y = grid::unit(0, "npc"),
      just = c("left", "bottom")
    )
  },
  width = 15,
  height = 2
)

p_jader_solid_ici_reporting_landscape <- grid::grid.grabExpr(
  {
    grid::pushViewport(
      grid::viewport(
        x = grid::unit(0.5, "npc"),
        y = jader_landscape_main_viewport_y,
        width = grid::unit(1, "npc"),
        height = jader_landscape_main_viewport_height,
        just = c("center", "center")
      )
    )
    grid::grid.draw(p_jader_solid_ici_landscape_heatmap_grob)
    grid::popViewport()

    grid::pushViewport(
      grid::viewport(
        x = grid::unit(0.055, "npc"),
        y = jader_landscape_heatmap_legend_y,
        width = grid::unit(18, "mm"),
        height = jader_landscape_heatmap_legend_height,
        just = c("left", "center"),
        angle = 90
      )
    )
    grid::grid.draw(p_jader_solid_ici_landscape_heatmap_legend_grob)
    grid::popViewport()

    grid::grid.text(
      "SOC reports",
      x = jader_landscape_soc_bar_title_x,
      y = jader_landscape_soc_bar_title_y,
      just = c("left", "center"),
      gp = grid::gpar(fontsize = 15, col = "black", fontfamily = "sans")
    )

    grid::pushViewport(
      grid::viewport(
        x = grid::unit(0.085, "npc"),
        y = jader_landscape_primary_legend_y,
        width = grid::unit(0.88, "npc"),
        height = grid::unit(0.18, "npc"),
        just = c("left", "bottom")
      )
    )
    grid::grid.draw(p_jader_solid_ici_landscape_primary_legend_grob)
    grid::popViewport()
  },
  width = 11,
  height = 10
)

print(p_jader_solid_ici_reporting_landscape)
grid::grid.newpage()
grid::grid.draw(p_jader_solid_ici_reporting_landscape)

### 3. 保存图形####

grDevices::cairo_pdf(
  file = file.path(script_results_path, "JADER实体瘤_ICI报告SOC景观图.pdf"),
  width = 11,
  height = 7,
  family = "sans"
)
grid::grid.newpage()
grid::grid.draw(p_jader_solid_ici_reporting_landscape)
grDevices::dev.off()


## 7. 保存R对象####

save(
  jader_solid_tumor_final,
  filter_summary,
  primary_organ_summary,
  histologic_summary,
  assigned_type_summary,
  ici_type_summary,
  ici_agent_summary,
  ici_regimen_summary,
  ici_strategy_summary,
  ici_sankey_flow_summary,
  ici_sankey_plot_data,
  ici_sankey_node_color_table,
  annual_trend_years,
  solid_yearly_total,
  ici_yearly_by_strategy,
  ici_yearly_total,
  jader_solid_ici_landscape,
  ici_landscape_column_meta,
  ici_landscape_soc_long,
  ici_landscape_soc_summary,
  ici_landscape_color_mapping,
  valid_organ_histology,
  jader_hierarchy_export,
  p_jader_solid_sankey,
  p_jader_solid_annual_trend_bar,
  p_jader_solid_ici_reporting_landscape,
  file = object_output_path
)

message("\n结果文件已保存:")
message("  汇总表: ", summary_output_path)
message("  对象: ", object_output_path)
message("  JADER三层结构: ", jader_hierarchy_export_path)
message("  桑基图: ", file.path(script_results_path, "JADER实体瘤_ICI情况桑基图.pdf"))
message("  年度趋势图: ", file.path(script_results_path, "JADER实体瘤_年度报告与ICI趋势图.pdf"))
message("  SOC景观图: ", file.path(script_results_path, "JADER实体瘤_ICI报告SOC景观图.pdf"))

message("\n========== JADER实体瘤描述性分析绘图完成 ==========")
