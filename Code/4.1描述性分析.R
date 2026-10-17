# FAERS实体瘤描述性分析####

rm(list = ls())
gc()

source("Code/config.R")

script_name <- "4.1描述性分析"
script_results_path <- file.path(results_path, script_name)
if (!dir.exists(script_results_path)) {
  dir.create(script_results_path, recursive = TRUE)
}

message("========== FAERS实体瘤描述性分析开始 ==========")
message("结果目录: ", script_results_path)


## 参数设置####

# 实体瘤识别已在faers_solid上游完成，本脚本不再二次筛选实体瘤

# 年度趋势图展示年份范围，max_year为NULL时自动使用数据中的最大年份
annual_year_min <- 2015
annual_year_max <- NULL

# GTEx映射表来自7.1脚本，景观图顶部注释按primary_organ接入该映射
gtex_mapping_path <- file.path(
  results_path,
  "7.1TCGA基础数据准备",
  "tables",
  "模式群-原发部位-TCGA项目-GTEx器官映射表（含样本量）.xlsx"
)


## 辅助函数####

### 格式化样本量，便于图表标注####
format_n <- function(x) {
  format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
}

### 从RData文件中读取指定对象####
# 功能：从指定RData文件中读取目标对象，并避免污染全局环境####
# 输入：input_path，RData文件路径；object_name，对象名称####
# 输出：返回读取后的R对象####
load_named_rdata_object <- function(input_path, object_name) {
  input_env <- new.env(parent = emptyenv())
  load(input_path, envir = input_env)

  if (!exists(object_name, envir = input_env, inherits = FALSE)) {
    stop(paste0("RData中未找到对象 ", object_name))
  }

  get(object_name, envir = input_env, inherits = FALSE)
}

### 加载FAERS实体瘤分析集####
# 功能：从RData文件中读取faers_solid对象####
# 输入：input_path，faers_solid.RData文件路径####
# 输出：返回faers_solid数据框####
load_faers_solid <- function(input_path) {
  load_named_rdata_object(input_path, "faers_solid")
}

### 生成指定层级的频数和构成比汇总表####
build_distribution_summary <- function(data, group_var) {
  group_var <- rlang::ensym(group_var)

  data %>%
    count(!!group_var, sort = TRUE, name = "n") %>%
    mutate(
      pct = n / sum(n) * 100,
      pct_label = paste0(round(pct, 1), "%")
    )
}

### 标准化ICI regimen显示文本####
format_ici_regimen_label <- function(x) {
  x %>%
    str_replace_all("anti-PD-1", "anti–PD-1") %>%
    str_replace_all("anti-PD-L1", "anti–PD-L1") %>%
    str_replace_all("anti-CTLA-4", "anti–CTLA-4") %>%
    str_replace_all("\\+", " + ")
}

### 标准化therapy strategy显示文本####
format_therapy_strategy_label <- function(x) {
  case_when(
    x == "anti-PD-1" ~ "anti–PD-1 monotherapy",
    x == "anti-PD-L1" ~ "anti–PD-L1 monotherapy",
    x == "anti-CTLA-4" ~ "anti–CTLA-4 monotherapy",
    x == "Combination Therapy" ~ "Combination therapy",
    TRUE ~ as.character(x)
  )
}

### 标准化年度趋势图中的therapy strategy文本####
format_therapy_strategy_bar_label <- function(x) {
  format_therapy_strategy_label(x)
}

### 从quarter字段中提取报告年份####
extract_report_year <- function(x) {
  suppressWarnings(as.integer(substr(as.character(x), 1, 4)))
}

### 标准化quarter文本，便于排序和展示####
normalize_report_quarter <- function(x) {
  x %>%
    as.character() %>%
    str_trim() %>%
    toupper() %>%
    str_replace_all("\\s+", "")
}

### 从quarter字段中提取可排序的季度索引####
extract_report_quarter_index <- function(x) {
  quarter_text <- normalize_report_quarter(x)
  report_year <- suppressWarnings(as.integer(str_sub(quarter_text, 1, 4)))
  report_quarter <- suppressWarnings(as.integer(str_extract(quarter_text, "(?<=Q)[1-4]")))

  if_else(
    is.na(report_year) | is.na(report_quarter),
    NA_integer_,
    report_year * 10L + report_quarter
  )
}

### 标准化癌种层级标签，统一缺失值和Other类表述####
standardize_cancer_display_label <- function(x, unknown_patterns = character(0), default_label = "Unknown/Others") {
  x_chr <- as.character(x)
  x_chr <- str_trim(x_chr)
  x_chr[is.na(x_chr) | x_chr == ""] <- default_label

  if (length(unknown_patterns) > 0) {
    pattern_text <- paste(unknown_patterns, collapse = "|")
    x_chr[str_detect(x_chr, regex(pattern_text, ignore_case = TRUE))] <- default_label
  }

  x_chr
}

### 为分类变量生成稳定配色，并将Unknown/Others固定为灰色####
build_categorical_color_map <- function(labels, unknown_label = "Unknown/Others", unknown_color = "#BDBDBD") {
  labels <- unique(as.character(labels))
  labels <- labels[!is.na(labels) & labels != ""]

  non_unknown_labels <- setdiff(labels, unknown_label)
  color_values <- character(0)

  if (length(non_unknown_labels) > 0) {
    color_values <- scales::hue_pal(l = 68, c = 110)(length(non_unknown_labels))
    names(color_values) <- non_unknown_labels
  }

  if (unknown_label %in% labels) {
    color_values <- c(color_values, stats::setNames(unknown_color, unknown_label))
  }

  color_values[labels]
}

### 按比例混合两种颜色，生成层级相近色####
blend_hex_color <- function(color, target_color = "#FFFFFF", target_ratio = 0.25) {
  target_ratio <- max(min(target_ratio, 1), 0)

  color_rgb <- grDevices::col2rgb(color)
  target_rgb <- grDevices::col2rgb(target_color)
  blended_rgb <- round(color_rgb * (1 - target_ratio) + target_rgb * target_ratio)

  grDevices::rgb(
    red = blended_rgb[1, ] / 255,
    green = blended_rgb[2, ] / 255,
    blue = blended_rgb[3, ] / 255
  )
}

### 标准化癌种配色匹配键，兼容空格斜线与下划线####
normalize_cancer_color_key <- function(x) {
  x %>%
    as.character() %>%
    str_to_lower() %>%
    str_replace_all("&", "and") %>%
    str_replace_all("[^a-z0-9]+", "")
}

### 构建高辨识度主层癌种固定色板####
build_reference_primary_palette <- function() {
  c(
    "lung" = "#2F6DB2",
    "skin" = "#6EC1E4",
    "melanoma" = "#6EC1E4",
    "eye" = "#6EC1E4",
    "uvealmelanoma" = "#6EC1E4",
    "headandneck" = "#F4A259",
    "headneck" = "#F4A259",
    "headandneckscc" = "#F4A259",
    "kidney" = "#39A74A",
    "breast" = "#63C65E",
    "bladder" = "#B89220",
    "esophagus" = "#F3C54D",
    "esophageal" = "#F3C54D",
    "stomach" = "#595959",
    "gastric" = "#595959",
    "esophagusstomach" = "#D6A23D",
    "colorectal" = "#2E9C8E",
    "bowel" = "#2E9C8E",
    "colon" = "#2E9C8E",
    "rectum" = "#5BBFB5",
    "prostate" = "#5BC0BE",
    "ovary" = "#FF4B4B",
    "ovarian" = "#FF4B4B",
    "lymphoma" = "#FF7682",
    "hematologic" = "#C77DD4",
    "leukemia" = "#C88AD8",
    "glioblastoma" = "#9B63A8",
    "cns" = "#9B63A8",
    "brain" = "#9B63A8",
    "nasopharyngeal" = "#8B5A3C",
    "liver" = "#D44C8C",
    "hepatocellular" = "#D44C8C",
    "pancreas" = "#7A9E1F",
    "pancreatic" = "#7A9E1F",
    "biliary" = "#C8A638",
    "biliarytract" = "#1F78B4",
    "adrenal" = "#9C6B4E",
    "adrenalgland" = "#E31A1C",
    "cervix" = "#F08C67",
    "pleura" = "#7F7F7F",
    "mesothelioma" = "#7F7F7F",
    "thyroid" = "#F39AAE",
    "testis" = "#94A83D",
    "penis" = "#6A3D9A",
    "thymus" = "#A58BC7",
    "uterus" = "#EC8FB4",
    "vulvavagina" = "#33A02C",
    "softtissue" = "#C06C84",
    "sarcoma" = "#C06C84",
    "bone" = "#FF7F00",
    "ampullaofvater" = "#A6CEE3",
    "peritoneum" = "#8E7C72",
    "rare" = "#7F6AAD",
    "pns" = "#4F7F8A",
    "peripheralnervoussystem" = "#4F7F8A",
    "others" = "#F0A5C5",
    "other" = "#F0A5C5",
    "unknownothers" = "#BDBDBD"
  )
}

### 为未预设的癌种补充高对比备用颜色####
build_fallback_distinct_palette <- function(n) {
  fallback_palette <- c(
    "#1F78B4", "#E31A1C", "#33A02C", "#FF7F00", "#6A3D9A", "#A6CEE3",
    "#B15928", "#E6AB02", "#FB9A99", "#1B9E77", "#D95F02", "#7570B3",
    "#66A61E", "#E7298A", "#A6761D", "#666666", "#8DD3C7", "#B3DE69"
  )

  rep(fallback_palette, length.out = n)
}

### 为第一层癌种生成高区分度主色####
build_primary_layer_color_map <- function(labels, unknown_label = "Unknown/Others", unknown_color = "#BDBDBD") {
  labels <- unique(as.character(labels))
  labels <- labels[!is.na(labels) & labels != ""]

  non_unknown_labels <- setdiff(labels, unknown_label)
  color_values <- character(0)

  if (length(non_unknown_labels) > 0) {
    normalized_keys <- normalize_cancer_color_key(non_unknown_labels)
    reference_palette <- build_reference_primary_palette()
    matched_colors <- unname(reference_palette[normalized_keys])

    if (any(is.na(matched_colors))) {
      fallback_labels <- non_unknown_labels[is.na(matched_colors)]
      fallback_colors <- build_fallback_distinct_palette(length(fallback_labels))
      matched_colors[is.na(matched_colors)] <- fallback_colors
    }

    color_values <- matched_colors
    names(color_values) <- non_unknown_labels
  }

  if (unknown_label %in% labels) {
    color_values <- c(color_values, stats::setNames(unknown_color, unknown_label))
  }

  color_values[labels]
}

### 基于父层颜色为子层生成同色系且更易区分的颜色####
build_child_layer_color_map <- function(
  relation_data,
  parent_col,
  child_key_col,
  child_label_col,
  parent_color_map,
  unknown_label = "Unknown/Others",
  unknown_color = "#BDBDBD",
  min_target_ratio = 0.12,
  max_target_ratio = 0.4
) {
  relation_data <- tibble::as_tibble(relation_data) %>%
    distinct(
      !!rlang::sym(parent_col),
      !!rlang::sym(child_key_col),
      !!rlang::sym(child_label_col)
    )

  parent_values <- unique(as.character(relation_data[[parent_col]]))
  child_color_map <- character(0)

  for (parent_value in parent_values) {
    child_subset <- relation_data %>%
      filter(as.character(.data[[parent_col]]) == parent_value)

    child_keys <- as.character(child_subset[[child_key_col]])
    child_labels <- as.character(child_subset[[child_label_col]])
    parent_color <- unname(parent_color_map[parent_value])

    if (is.na(parent_color) || length(parent_color) == 0) {
      parent_color <- "#808080"
    }

    is_unknown <- is.na(child_labels) | child_labels == unknown_label
    non_unknown_n <- sum(!is_unknown)

    if (non_unknown_n > 0) {
      if (non_unknown_n == 1) {
        child_colors <- blend_hex_color(parent_color, "#FFFFFF", min_target_ratio + 0.02)
      } else {
        child_colors <- character(non_unknown_n)
        child_index <- seq_len(non_unknown_n)
        light_index <- child_index %% 2 == 1
        dark_index <- !light_index

        light_n <- sum(light_index)
        dark_n <- sum(dark_index)

        if (light_n > 0) {
          light_ratios <- seq(min_target_ratio, max_target_ratio, length.out = light_n)
          child_colors[light_index] <- vapply(
            light_ratios,
            function(current_ratio) blend_hex_color(parent_color, "#FFFFFF", current_ratio),
            character(1)
          )
        }

        if (dark_n > 0) {
          dark_ratios <- seq(min_target_ratio * 0.55, max_target_ratio * 0.55, length.out = dark_n)
          child_colors[dark_index] <- vapply(
            dark_ratios,
            function(current_ratio) blend_hex_color(parent_color, "#000000", current_ratio),
            character(1)
          )
        }
      }

      names(child_colors) <- child_keys[!is_unknown]
      child_color_map <- c(child_color_map, child_colors)
    }

    if (any(is_unknown)) {
      unknown_colors <- rep(unknown_color, sum(is_unknown))
      names(unknown_colors) <- child_keys[is_unknown]
      child_color_map <- c(child_color_map, unknown_colors)
    }
  }

  child_color_map[!duplicated(names(child_color_map))]
}

### 构建癌种三层层级继承配色####
build_hierarchical_cancer_color_maps <- function(
  column_meta,
  unknown_label = "Unknown/Others",
  unknown_color = "#BDBDBD"
) {
  column_meta <- column_meta %>%
    mutate(
      histologic_color_key = paste(primary_organ_plot, histologic_category_plot, sep = "|||"),
      tcga_color_key = paste(primary_organ_plot, histologic_category_plot, assigned_type_plot, sep = "|||")
    )

  primary_levels <- column_meta %>%
    distinct(primary_organ_plot) %>%
    pull(primary_organ_plot)

  primary_color_map <- build_primary_layer_color_map(
    labels = primary_levels,
    unknown_label = unknown_label,
    unknown_color = unknown_color
  )

  histologic_relation <- column_meta %>%
    distinct(primary_organ_plot, histologic_category_plot, histologic_color_key)

  histologic_color_map <- build_child_layer_color_map(
    relation_data = histologic_relation,
    parent_col = "primary_organ_plot",
    child_key_col = "histologic_color_key",
    child_label_col = "histologic_category_plot",
    parent_color_map = primary_color_map,
    unknown_label = unknown_label,
    unknown_color = unknown_color,
    min_target_ratio = 0.08,
    max_target_ratio = 0.2
  )

  tcga_relation <- column_meta %>%
    distinct(primary_organ_plot, histologic_category_plot, histologic_color_key, assigned_type_plot, tcga_color_key)

  tcga_color_map <- build_child_layer_color_map(
    relation_data = tcga_relation,
    parent_col = "histologic_color_key",
    child_key_col = "tcga_color_key",
    child_label_col = "assigned_type_plot",
    parent_color_map = histologic_color_map,
    unknown_label = unknown_label,
    unknown_color = unknown_color,
    min_target_ratio = 0.06,
    max_target_ratio = 0.16
  )

  color_mapping <- bind_rows(
    tibble(
      level_name = "Primary organ",
      parent_category = NA_character_,
      color_key = primary_levels,
      category = primary_levels,
      color = unname(primary_color_map[primary_levels])
    ),
    histologic_relation %>%
      transmute(
        level_name = "Histologic category",
        parent_category = primary_organ_plot,
        color_key = histologic_color_key,
        category = histologic_category_plot,
        color = unname(histologic_color_map[histologic_color_key])
      ),
    tcga_relation %>%
      transmute(
        level_name = "TCGA",
        parent_category = paste(primary_organ_plot, histologic_category_plot, sep = " | "),
        color_key = tcga_color_key,
        category = assigned_type_plot,
        color = unname(tcga_color_map[tcga_color_key])
      )
  )

  list(
    column_meta = column_meta,
    primary_color_map = primary_color_map,
    histologic_color_map = histologic_color_map,
    tcga_color_map = tcga_color_map,
    color_mapping = color_mapping
  )
}

### 构建与2.1分析一致的ICI固定配色####
build_ici_reference_color_map <- function() {
  c(
    "anti–PD-L1" = "#6ECCC8",
    "anti–PD-L1 monotherapy" = "#6ECCC8",
    "anti–PD-1" = "#31859C",
    "anti–PD-1 monotherapy" = "#31859C",
    "anti–CTLA-4" = "#1F4E79",
    "anti–CTLA-4 monotherapy" = "#1F4E79",
    "Combination Therapy" = "#C5E033",
    "Combination therapy" = "#C5E033",
    "Other" = "#BDBDBD",
    "Unknown/Others" = "#BDBDBD"
  )
}

### 根据2.1固定配色生成桑基图节点颜色####
build_ici_sankey_node_color_map <- function(sankey_data, reference_color_map) {
  node_color_data <- bind_rows(
    sankey_data %>%
      distinct(ici_class_label, ici_class) %>%
      transmute(
        node = ici_class_label,
        node_color = unname(reference_color_map[ici_class])
      ),
    sankey_data %>%
      distinct(agent_label, ici_class) %>%
      transmute(
        node = agent_label,
        node_color = unname(reference_color_map[ici_class])
      ),
    sankey_data %>%
      distinct(regimen_label, therapy_strategy) %>%
      transmute(
        node = regimen_label,
        node_color = unname(reference_color_map[therapy_strategy])
      ),
    sankey_data %>%
      distinct(strategy_label, therapy_strategy) %>%
      transmute(
        node = strategy_label,
        node_color = unname(reference_color_map[therapy_strategy])
      )
  ) %>%
    mutate(node_color = replace_na(node_color, unname(reference_color_map["Other"]))) %>%
    distinct(node, .keep_all = TRUE)

  stats::setNames(node_color_data$node_color, node_color_data$node)
}

### 提取数据中的SOC元数据，便于构建发生谱矩阵####
extract_soc_metadata_simple <- function(data) {
  soc_vars <- grep("^SOC_", names(data), value = TRUE)

  if (length(soc_vars) == 0) {
    stop("数据中未检测到SOC列，无法绘制第三张图")
  }

  tibble(
    soc_var = soc_vars,
    SOC_full = str_replace_all(sub("^SOC_", "", soc_vars), "_", " ")
  )
}

### 构建SOC全称到缩写的映射表####
build_soc_abbrev_map <- function() {
  c(
    "Blood and lymphatic system disorders" = "HEME",
    "Cardiac disorders" = "CARD",
    "Endocrine disorders" = "ENDO",
    "Eye disorders" = "EYE",
    "Gastrointestinal disorders" = "GI",
    "Hepatobiliary disorders" = "HEPA",
    "Musculoskeletal and connective tissue disorders" = "MSK",
    "Nervous system disorders" = "NEURO",
    "Renal and urinary disorders" = "RENU",
    "Respiratory, thoracic and mediastinal disorders" = "RESP",
    "Skin and subcutaneous tissue disorders" = "SKIN"
  )
}

### 将SOC全称格式化为“全称（缩写）”展示文本####
format_soc_display_label <- function(soc_full, soc_abbrev_map = build_soc_abbrev_map()) {
  soc_full_chr <- as.character(soc_full)
  soc_abbrev <- unname(soc_abbrev_map[soc_full_chr])

  ifelse(
    !is.na(soc_abbrev) & soc_abbrev != "",
    paste0(soc_full_chr, " (", soc_abbrev, ")"),
    soc_full_chr
  )
}

### 将SOC全称格式化为缩写展示文本####
# 功能：优先返回SOC缩写，缺少缩写映射时保留SOC全称####
# 输入：soc_full，SOC全称向量；soc_abbrev_map，SOC全称到缩写的命名向量####
# 输出：用于绘图行名的SOC缩写向量####
format_soc_abbrev_label <- function(soc_full, soc_abbrev_map = build_soc_abbrev_map()) {
  soc_full_chr <- as.character(soc_full)
  soc_abbrev <- unname(soc_abbrev_map[soc_full_chr])

  ifelse(
    !is.na(soc_abbrev) & soc_abbrev != "",
    soc_abbrev,
    soc_full_chr
  )
}

### 格式化百分比标签，避免最大值被四舍五入为100%####
format_pct_label <- function(x) {
  label_raw <- ifelse(
    abs(x - round(x)) < 1e-8,
    as.character(round(x)),
    scales::number(x, accuracy = 0.1, trim = TRUE)
  )

  paste0(label_raw, "%")
}

### 计算热图颜色映射上限，保留数据中的真实最大比例####
compute_heatmap_upper_bound <- function(mat) {
  mat_values <- as.numeric(mat)
  mat_values <- mat_values[is.finite(mat_values) & !is.na(mat_values)]

  if (length(mat_values) == 0) {
    return(1)
  }

  upper_bound <- max(mat_values)

  if (!is.finite(upper_bound) || upper_bound <= 0) {
    return(1)
  }

  upper_bound
}

### 构建热图分段颜色断点，增强低比例区间的层次####
build_heatmap_color_breaks <- function(mat, upper_bound) {
  mat_values <- as.numeric(mat)
  mat_values <- mat_values[is.finite(mat_values) & !is.na(mat_values)]
  positive_values <- mat_values[mat_values > 0]

  if (length(positive_values) == 0) {
    return(c(0, upper_bound))
  }

  quantile_breaks <- stats::quantile(
    positive_values,
    probs = c(0, 0.5, 0.75, 0.9, 0.95, 0.99, 1),
    na.rm = TRUE,
    names = FALSE
  )

  heatmap_breaks <- c(0, quantile_breaks, upper_bound)
  heatmap_breaks <- unique(round(heatmap_breaks, 4))
  heatmap_breaks <- heatmap_breaks[
    heatmap_breaks >= 0 &
      heatmap_breaks <= upper_bound
  ]
  heatmap_breaks <- sort(unique(heatmap_breaks))

  if (length(heatmap_breaks) < 2) {
    heatmap_breaks <- c(0, upper_bound)
  }

  heatmap_breaks
}

### 确保SOC展示列存在，兼容旧对象####
ensure_soc_display_column <- function(data, soc_abbrev_map = build_soc_abbrev_map()) {
  data <- tibble::as_tibble(data)

  if (!"SOC_full" %in% names(data)) {
    stop("数据中缺少 SOC_full 列，无法构建SOC展示标签")
  }

  data <- data %>%
    mutate(SOC_full = as.character(SOC_full))

  if ("SOC_display" %in% names(data)) {
    data <- data %>%
      mutate(SOC_display = as.character(SOC_display))
  } else {
    data <- data %>%
      mutate(SOC_display = format_soc_display_label(SOC_full, soc_abbrev_map))
  }

  data
}

### 解析Excel读取后的逻辑标记####
# 功能：将Excel中的0/1、TRUE/FALSE或文本标记统一转换为逻辑值####
# 输入：input_vector，待转换的标记向量####
# 输出：逻辑向量####
parse_excel_logical_flag <- function(input_vector) {
  normalized_vector <- tolower(trimws(as.character(input_vector)))

  case_when(
    normalized_vector %in% c("true", "t", "1", "yes", "y") ~ TRUE,
    normalized_vector %in% c("false", "f", "0", "no", "n", "") ~ FALSE,
    is.na(normalized_vector) ~ FALSE,
    TRUE ~ FALSE
  )
}

### 读取景观图GTEx注释映射####
# 功能：读取7.1生成的GTEx映射表，并整理为可与景观图列元数据合并的格式####
# 输入：mapping_path，7.1输出的模式群-原发部位-TCGA项目-GTEx器官映射表路径####
# 输出：包含primary_organ_plot、GTEx组织标签和样本量信息的数据框####
read_gtex_landscape_mapping <- function(mapping_path) {
  if (!file.exists(mapping_path)) {
    stop("未找到GTEx映射表: ", mapping_path)
  }

  available_sheets <- readxl::excel_sheets(mapping_path)
  preferred_sheet <- "原发部位映射总表（含样本量）"
  fallback_sheet <- "原发部位映射总表"

  if (preferred_sheet %in% available_sheets) {
    mapping_sheet <- preferred_sheet
  } else if (fallback_sheet %in% available_sheets) {
    mapping_sheet <- fallback_sheet
  } else {
    stop("GTEx映射表中未找到原发部位映射总表sheet: ", mapping_path)
  }

  gtex_mapping_raw <- readxl::read_xlsx(
    path = mapping_path,
    sheet = mapping_sheet
  )

  required_columns <- c("primary_organ", "selected_gtex_general_tissue", "has_gtex_mapping")
  missing_columns <- setdiff(required_columns, names(gtex_mapping_raw))

  if (length(missing_columns) > 0) {
    stop("GTEx映射表缺少必要列: ", paste(missing_columns, collapse = ", "))
  }

  if (!"any_proxy_tissue" %in% names(gtex_mapping_raw)) {
    gtex_mapping_raw$any_proxy_tissue <- FALSE
  }

  if (!"n_gtex_normal" %in% names(gtex_mapping_raw)) {
    gtex_mapping_raw$n_gtex_normal <- NA_integer_
  }

  gtex_mapping_raw %>%
    transmute(
      primary_organ = as.character(primary_organ),
      primary_organ_plot = standardize_cancer_display_label(primary_organ),
      selected_gtex_general_tissue = str_squish(as.character(selected_gtex_general_tissue)),
      gtex_has_mapping = parse_excel_logical_flag(has_gtex_mapping),
      gtex_proxy_tissue = parse_excel_logical_flag(any_proxy_tissue),
      n_gtex_normal = suppressWarnings(as.integer(n_gtex_normal))
    ) %>%
    mutate(
      selected_gtex_general_tissue = na_if(selected_gtex_general_tissue, ""),
      gtex_has_mapping = gtex_has_mapping & !is.na(selected_gtex_general_tissue),
      gtex_tissue_plot = case_when(
        !gtex_has_mapping ~ "No GTEx mapping",
        gtex_proxy_tissue ~ paste0(selected_gtex_general_tissue, " (proxy)"),
        TRUE ~ selected_gtex_general_tissue
      )
    ) %>%
    distinct(primary_organ_plot, .keep_all = TRUE)
}

### 合并GTEx映射到景观图列元数据####
# 功能：按primary_organ为景观图每个列单元补充GTEx注释信息####
# 输入：column_meta，景观图列元数据；gtex_mapping_df，GTEx映射表####
# 输出：补充GTEx组织标签和样本量后的列元数据####
add_gtex_mapping_to_landscape_column_meta <- function(column_meta, gtex_mapping_df) {
  no_gtex_mapping_label <- "No GTEx mapping"

  if (is.null(gtex_mapping_df) || nrow(gtex_mapping_df) == 0) {
    return(
      column_meta %>%
        mutate(
          selected_gtex_general_tissue = NA_character_,
          gtex_has_mapping = FALSE,
          gtex_proxy_tissue = FALSE,
          n_gtex_normal = 0L,
          gtex_tissue_plot = no_gtex_mapping_label
        )
    )
  }

  column_meta %>%
    left_join(
      gtex_mapping_df %>%
        select(
          primary_organ_plot,
          selected_gtex_general_tissue,
          gtex_has_mapping,
          gtex_proxy_tissue,
          n_gtex_normal,
          gtex_tissue_plot
        ),
      by = "primary_organ_plot"
    ) %>%
    mutate(
      selected_gtex_general_tissue = if_else(
        is.na(selected_gtex_general_tissue),
        NA_character_,
        selected_gtex_general_tissue
      ),
      gtex_has_mapping = coalesce(gtex_has_mapping, FALSE),
      gtex_proxy_tissue = coalesce(gtex_proxy_tissue, FALSE),
      n_gtex_normal = coalesce(n_gtex_normal, 0L),
      gtex_tissue_plot = coalesce(gtex_tissue_plot, no_gtex_mapping_label)
    )
}

### 构建GTEx注释配色####
# 功能：让GTEx组织注释继承对应Primary organ颜色，并将无GTEx映射固定为灰色####
# 输入：column_meta，景观图列元数据；primary_color_map，Primary organ配色；no_mapping_label，无映射标签####
# 输出：命名颜色向量####
build_gtex_tissue_color_map <- function(column_meta, primary_color_map, no_mapping_label = "No GTEx mapping") {
  required_columns <- c("primary_organ_plot", "gtex_tissue_plot")
  missing_columns <- setdiff(required_columns, names(column_meta))

  if (length(missing_columns) > 0) {
    stop("GTEx注释配色缺少必要列: ", paste(missing_columns, collapse = ", "))
  }

  gtex_color_relation <- column_meta %>%
    transmute(
      primary_organ_plot = as.character(primary_organ_plot),
      gtex_tissue_plot = as.character(gtex_tissue_plot)
    ) %>%
    filter(!is.na(gtex_tissue_plot), gtex_tissue_plot != "") %>%
    distinct(gtex_tissue_plot, primary_organ_plot) %>%
    group_by(gtex_tissue_plot) %>%
    summarise(
      primary_organ_plot = first(primary_organ_plot),
      .groups = "drop"
    )

  color_values <- unname(primary_color_map[gtex_color_relation$primary_organ_plot])
  color_values[gtex_color_relation$gtex_tissue_plot == no_mapping_label] <- "#D0D0D0"
  color_values[is.na(color_values) | color_values == ""] <- "#808080"
  names(color_values) <- gtex_color_relation$gtex_tissue_plot

  color_values[gtex_color_relation$gtex_tissue_plot]
}

### 构建第三张图所需的ICI报告景观图数据####
# 功能：按癌种层级、季度和SOC汇总ICI报告比例，并补充GTEx注释映射####
# 输入：data，实体瘤分析集；gtex_mapping_path，GTEx映射表路径####
# 输出：包含景观图列元数据、SOC长表、SOC汇总表和注释配色的列表####
build_ici_soc_landscape_data <- function(data, gtex_mapping_path = NULL) {
  soc_metadata <- extract_soc_metadata_simple(data)

  ici_data <- data %>%
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
      ),
      quarter_plot = normalize_report_quarter(quarter),
      quarter_index = extract_report_quarter_index(quarter)
    ) %>%
    filter(!is.na(quarter_index))

  if (nrow(ici_data) == 0) {
    stop("ICI病例中不存在可用的quarter信息，无法绘制第三张图")
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
      quarter_plot,
      quarter_index,
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
      desc(report_cases),
      quarter_index,
      quarter_plot
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
      quarter_plot,
      quarter_index
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
          quarter_plot,
          quarter_index,
          column_id
        ),
      by = c(
        "primary_organ_plot",
        "histologic_category_plot",
        "assigned_type_plot",
        "quarter_plot",
        "quarter_index"
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

### 生成实体瘤ICI用药桑基图源数据####
build_ici_sankey_source_data <- function(data) {
  ici_class_map <- tibble(
    drug_name_upper = c(
      "NIVOLUMAB", "PEMBROLIZUMAB", "CEMIPLIMAB", "TISLELIZUMAB",
      "DOSTARLIMAB", "TORIPALIMAB", "RETIFANLIMAB",
      "ATEZOLIZUMAB", "DURVALUMAB", "AVELUMAB",
      "IPILIMUMAB", "TREMELIMUMAB"
    ),
    ici_class = c(
      rep("anti–PD-1", 7),
      rep("anti–PD-L1", 3),
      rep("anti–CTLA-4", 2)
    )
  )

  ici_detail <- data %>%
    filter(has_ici == TRUE, !is.na(ici_drugs), ici_drugs != "None") %>%
    select(primaryid, ici_drugs, ici_types, treatment_strategy) %>%
    distinct()

  drug_records <- ici_detail %>%
    separate_rows(ici_drugs, sep = "\\+") %>%
    transmute(
      primaryid = primaryid,
      drug_name_upper = toupper(trimws(ici_drugs))
    ) %>%
    left_join(ici_class_map, by = "drug_name_upper") %>%
    mutate(
      ici_class = replace_na(ici_class, "Other"),
      drug_name = str_to_title(str_to_lower(drug_name_upper))
    )

  ici_type_counts <- drug_records %>%
    count(ici_class, name = "type_n", sort = TRUE)

  agent_counts <- drug_records %>%
    count(drug_name, name = "agent_n", sort = TRUE)

  patient_info <- ici_detail %>%
    transmute(
      primaryid = primaryid,
      ici_regimen = format_ici_regimen_label(ici_types),
      therapy_strategy = format_therapy_strategy_label(treatment_strategy)
    )

  regimen_counts <- patient_info %>%
    count(ici_regimen, name = "regimen_n", sort = TRUE)

  strategy_counts <- patient_info %>%
    count(therapy_strategy, name = "strategy_n", sort = TRUE)

  sankey_data <- drug_records %>%
    left_join(ici_type_counts, by = "ici_class") %>%
    left_join(agent_counts, by = "drug_name") %>%
    left_join(patient_info, by = "primaryid") %>%
    left_join(regimen_counts, by = "ici_regimen") %>%
    left_join(strategy_counts, by = "therapy_strategy") %>%
    mutate(
      ici_class_label = ici_class,
      agent_label = drug_name,
      regimen_label = ici_regimen,
      strategy_label = paste0(therapy_strategy, " (n = ", format_n(strategy_n), ")")
    )

  sankey_long <- ggsankey::make_long(
    sankey_data,
    ici_class_label,
    agent_label,
    regimen_label,
    strategy_label
  ) %>%
    mutate(
      x_new = c(
        "ici_class_label" = 0.9,
        "agent_label" = 1.35,
        "regimen_label" = 1.95,
        "strategy_label" = 2.9
      )[as.character(x)],
      next_x_new = c(
        "ici_class_label" = 0.9,
        "agent_label" = 1.35,
        "regimen_label" = 1.95,
        "strategy_label" = 2.9
      )[as.character(next_x)]
    )

  list(
    ici_detail = ici_detail,
    drug_records = drug_records,
    ici_type_counts = ici_type_counts,
    agent_counts = agent_counts,
    regimen_counts = regimen_counts,
    strategy_counts = strategy_counts,
    sankey_data = sankey_data,
    sankey_long = sankey_long
  )
}

### 构建年度实体瘤报告与ICI使用趋势数据####
build_annual_trend_data <- function(data, year_min = 2015, year_max = NULL) {
  data_with_year <- data %>%
    mutate(year = extract_report_year(quarter))

  available_years <- data_with_year %>%
    filter(!is.na(year)) %>%
    pull(year) %>%
    unique() %>%
    sort()

  max_available_year <- max(available_years, na.rm = TRUE)
  max_display_year <- if (is.null(year_max)) max_available_year else min(year_max, max_available_year)
  available_years <- seq(year_min, max_display_year)

  if (!is.null(year_max)) {
    available_years <- available_years[available_years <= year_max]
  }

  strategy_levels <- format_therapy_strategy_bar_label(
    c("anti-PD-1", "anti-PD-L1", "anti-CTLA-4", "Combination Therapy")
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

### 读取合法的第一层与第二层癌种映射关系####
# 功能：读取癌种映射表，并生成三层结构导出时使用的合法器官-组织学对照####
# 输入：mapping_path，癌种映射Excel文件路径####
# 输出：返回包含primary_organ_plot、histologic_category_plot和.valid_histology列的tibble####
build_valid_organ_histology_mapping <- function(mapping_path) {
  if (!file.exists(mapping_path)) {
    stop("未找到癌种映射文件: ", mapping_path)
  }

  cancer_mapping_ref <- openxlsx::read.xlsx(mapping_path)

  cancer_mapping_ref %>%
    transmute(
      primary_organ_plot = standardize_cancer_display_label(Primary_organs_tissues),
      histologic_category_plot = standardize_cancer_display_label(
        Histologic_category,
        unknown_patterns = c("^other$", "^others$", "^unknown$", "^unknown/others$")
      )
    ) %>%
    filter(histologic_category_plot != "Unknown/Others") %>%
    distinct() %>%
    mutate(.valid_histology = TRUE)
}

### 标准化三层结构统计底表####
# 功能：统一癌种标签，并将映射表未覆盖的第二层和第三层归入Unknown/Others####
# 输入：data，待统计的实体瘤数据框；valid_organ_histology，合法器官-组织学映射表####
# 输出：返回可直接用于三层结构计数的数据框####
prepare_three_layer_hierarchy_base <- function(data, valid_organ_histology) {
  data %>%
    mutate(
      primary_organ_plot = standardize_cancer_display_label(primary_organ),
      histologic_category_plot = standardize_cancer_display_label(
        histologic_category,
        unknown_patterns = c("^other$", "^others$", "^unknown$", "^unknown/others$")
      ),
      assigned_type_plot = standardize_cancer_display_label(
        assigned_type,
        unknown_patterns = c("^other$", "^others$", "^na$", "^n/a$", "^unknown/others$", "^unknown$")
      )
    ) %>%
    left_join(valid_organ_histology, by = c("primary_organ_plot", "histologic_category_plot")) %>%
    mutate(
      histologic_category_plot = if_else(
        is.na(.valid_histology) & histologic_category_plot != "Unknown/Others",
        "Unknown/Others",
        histologic_category_plot
      ),
      assigned_type_plot = if_else(
        is.na(.valid_histology),
        "Unknown/Others",
        assigned_type_plot
      )
    ) %>%
    select(-.valid_histology)
}

### 构建三层结构导出表####
# 功能：基于全量实体瘤与ICI实体瘤数据，生成第一层至第三层的汇总导出表####
# 输入：all_data，全量实体瘤数据框；ici_data，ICI实体瘤数据框；valid_organ_histology，合法器官-组织学映射表####
# 输出：返回格式化后的三层结构导出数据框####
build_three_layer_hierarchy_export <- function(all_data, ici_data, valid_organ_histology) {
  all_hierarchy_base <- prepare_three_layer_hierarchy_base(all_data, valid_organ_histology)
  ici_hierarchy_base <- prepare_three_layer_hierarchy_base(ici_data, valid_organ_histology)

  all_primary_n <- all_hierarchy_base %>%
    count(primary_organ_plot, name = "all_primary_n")

  all_histologic_n_merged <- all_hierarchy_base %>%
    count(primary_organ_plot, histologic_category_plot, name = "all_histologic_n")

  all_tcga_n_merged <- all_hierarchy_base %>%
    count(primary_organ_plot, histologic_category_plot, assigned_type_plot, name = "all_tcga_n")

  ici_hierarchy_counts <- ici_hierarchy_base %>%
    count(primary_organ_plot, histologic_category_plot, assigned_type_plot, name = "tcga_total")

  ici_histologic_counts <- ici_hierarchy_counts %>%
    group_by(primary_organ_plot, histologic_category_plot) %>%
    summarise(
      histologic_total = sum(tcga_total),
      .groups = "drop"
    )

  ici_primary_counts <- ici_histologic_counts %>%
    group_by(primary_organ_plot) %>%
    summarise(
      primary_total = sum(histologic_total),
      .groups = "drop"
    )

  all_tcga_n_merged %>%
    left_join(all_histologic_n_merged, by = c("primary_organ_plot", "histologic_category_plot")) %>%
    left_join(all_primary_n, by = "primary_organ_plot") %>%
    left_join(ici_hierarchy_counts, by = c("primary_organ_plot", "histologic_category_plot", "assigned_type_plot")) %>%
    left_join(ici_histologic_counts, by = c("primary_organ_plot", "histologic_category_plot")) %>%
    left_join(ici_primary_counts, by = "primary_organ_plot") %>%
    mutate(
      primary_total = coalesce(primary_total, 0L),
      histologic_total = coalesce(histologic_total, 0L),
      tcga_total = coalesce(tcga_total, 0L),
      primary_ici_pct = round(if_else(all_primary_n > 0, primary_total / all_primary_n * 100, 0), 2),
      histologic_ici_pct = round(if_else(all_histologic_n > 0, histologic_total / all_histologic_n * 100, 0), 2),
      tcga_ici_pct = round(if_else(all_tcga_n > 0, tcga_total / all_tcga_n * 100, 0), 2),
      primary_is_uo = primary_organ_plot == "Unknown/Others",
      histologic_is_uo = histologic_category_plot == "Unknown/Others",
      tcga_is_uo = assigned_type_plot == "Unknown/Others"
    ) %>%
    arrange(
      primary_is_uo,
      desc(all_primary_n),
      primary_organ_plot,
      histologic_is_uo,
      desc(all_histologic_n),
      histologic_category_plot,
      tcga_is_uo,
      desc(all_tcga_n),
      assigned_type_plot
    ) %>%
    mutate(
      primary_organ = paste0(primary_organ_plot, "(", format_n(all_primary_n), ")"),
      primary_organ_ici = paste0(format_n(primary_total), "(", primary_ici_pct, "%)"),
      histologic_category = paste0(histologic_category_plot, "(", format_n(all_histologic_n), ")"),
      histologic_category_ici = paste0(format_n(histologic_total), "(", histologic_ici_pct, "%)"),
      assigned_type = paste0(assigned_type_plot, "(", format_n(all_tcga_n), ")"),
      assigned_type_ici = paste0(format_n(tcga_total), "(", tcga_ici_pct, "%)")
    ) %>%
    select(
      primary_organ,
      primary_organ_ici,
      histologic_category,
      histologic_category_ici,
      assigned_type,
      assigned_type_ici
    )
}


## 1. 加载实体瘤分析集####

faers_solid_path <- file.path(output_path, "faers_solid.RData")
faers_solid <- load_faers_solid(faers_solid_path)

message("已加载 faers_solid")
message("总样本数: ", nrow(faers_solid))
message("变量数: ", ncol(faers_solid))


## 2. 直接使用实体瘤分析集####

faers_solid_tumor_final <- faers_solid
filter_summary <- tibble(
  step = "faers_solid",
  n = nrow(faers_solid_tumor_final),
  retained_pct = 100
)

message("\n当前脚本直接使用 faers_solid，不再进行实体瘤识别筛选:")
print(filter_summary)

message("\n实体瘤患者数: ", nrow(faers_solid_tumor_final))
message("原发器官数: ", dplyr::n_distinct(faers_solid_tumor_final$primary_organ))
message("组织学分类数: ", dplyr::n_distinct(faers_solid_tumor_final$histologic_category))
message("癌种实体数: ", dplyr::n_distinct(faers_solid_tumor_final$assigned_type))


## 3. 汇总表准备####

primary_organ_summary <- build_distribution_summary(faers_solid_tumor_final, primary_organ)
histologic_summary <- build_distribution_summary(faers_solid_tumor_final, histologic_category)
assigned_type_summary <- build_distribution_summary(faers_solid_tumor_final, assigned_type)

ici_sankey_results <- build_ici_sankey_source_data(faers_solid_tumor_final)
ici_type_summary <- ici_sankey_results$ici_type_counts
ici_agent_summary <- ici_sankey_results$agent_counts
ici_regimen_summary <- ici_sankey_results$regimen_counts
ici_strategy_summary <- ici_sankey_results$strategy_counts
ici_sankey_flow_summary <- ici_sankey_results$sankey_data %>%
  count(ici_class, drug_name, ici_regimen, therapy_strategy, name = "n", sort = TRUE)
annual_trend_results <- build_annual_trend_data(
  faers_solid_tumor_final,
  year_min = annual_year_min,
  year_max = annual_year_max
)
ici_soc_landscape_results <- build_ici_soc_landscape_data(
  faers_solid_tumor_final,
  gtex_mapping_path = gtex_mapping_path
)
solid_yearly_total <- annual_trend_results$solid_yearly_total
ici_yearly_by_strategy <- annual_trend_results$ici_yearly_by_strategy
ici_yearly_total <- annual_trend_results$ici_yearly_total
faers_solid_ici_landscape <- ici_soc_landscape_results$ici_data
ici_landscape_column_meta <- ici_soc_landscape_results$column_meta
ici_landscape_soc_long <- ici_soc_landscape_results$soc_long
ici_landscape_soc_summary <- ici_soc_landscape_results$soc_summary
ici_landscape_soc_abbrev_map <- ici_soc_landscape_results$soc_abbrev_map
ici_landscape_gtex_mapping <- ici_soc_landscape_results$gtex_mapping
ici_landscape_gtex_color_map <- ici_soc_landscape_results$gtex_color_map
ici_landscape_gtex_color_mapping <- ici_soc_landscape_results$gtex_color_mapping
ici_landscape_color_mapping <- ici_soc_landscape_results$color_mapping
ici_landscape_primary_color_map <- ici_soc_landscape_results$primary_color_map
ici_landscape_histologic_color_map <- ici_soc_landscape_results$histologic_color_map
ici_landscape_tcga_color_map <- ici_soc_landscape_results$tcga_color_map
ici_landscape_heatmap_max_pct <- ici_soc_landscape_results$heatmap_max_pct

message("\n原发器官分布（前10）:")
print(head(primary_organ_summary, 10))

message("\n组织学分类分布（前10）:")
print(head(histologic_summary, 10))

message("\n癌种实体分布（前10）:")
print(head(assigned_type_summary, 10))

message("\nICI class分布:")
print(ici_type_summary)

message("\nICI agent分布（前10）:")
print(head(ici_agent_summary, 10))

message("\n年度实体瘤报告数:")
print(solid_yearly_total)

message("\n年度ICI报告数:")
print(ici_yearly_total)

message("\n第三张图列元数据（前10）:")
print(head(ici_landscape_column_meta, 10))

message("\n第三张图SOC汇总（前10）:")
print(head(ici_landscape_soc_summary, 10))

message("\n第三张图GTEx映射（前10）:")
print(head(ici_landscape_gtex_mapping, 10))


## 4. 图1：ICI情况桑基图####

### 图1：数据准备####

ici_sankey_plot_data <- ici_sankey_results$sankey_long
ici_sankey_x_breaks <- c(0.9, 1.35, 1.95, 2.9)
ici_reference_color_map <- build_ici_reference_color_map()
ici_sankey_node_color_map <- build_ici_sankey_node_color_map(
  sankey_data = ici_sankey_results$sankey_data,
  reference_color_map = ici_reference_color_map
)

### 图1：绘图代码####

p_faers_solid_sankey <- ggplot(
  ici_sankey_plot_data,
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
  scale_fill_manual(values = ici_sankey_node_color_map, guide = "none") +
  scale_x_continuous(
    breaks = ici_sankey_x_breaks,
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

print(p_faers_solid_sankey)

### 图1：导出PDF####

ggsave(
  filename = file.path(script_results_path, "FAERS实体瘤_ICI情况桑基图.pdf"),
  plot = p_faers_solid_sankey,
  width = 11,
  height = 5,
  device = cairo_pdf,
  family = "sans"
)


## 5. 图2：年度报告与ICI趋势图####

### 图2：数据准备####

annual_trend_years <- annual_trend_results$available_years
annual_strategy_levels <- format_therapy_strategy_bar_label(
  c("anti-PD-L1", "anti-PD-1", "anti-CTLA-4", "Combination Therapy")
)
ici_yearly_by_strategy <- ici_yearly_by_strategy %>%
  dplyr::mutate(
    therapy_strategy = factor(as.character(therapy_strategy), levels = annual_strategy_levels)
  )
annual_strategy_color_map <- ici_reference_color_map[annual_strategy_levels]

### 图2：绘图代码####

p_faers_solid_annual_trend_bar <- ggplot() +
  geom_col(
    data = solid_yearly_total,
    aes(x = x_pos, y = solid_total),
    width = 0.32,
    fill = "#1F3B73",
    show.legend = FALSE
  ) +
  geom_text(
    data = solid_yearly_total,
    aes(x = x_pos, y = solid_total, label = format_n(solid_total)),
    vjust = -0.3,
    size = 12 / ggplot2::.pt,
    color = "black"
  ) +
  geom_col(
    data = ici_yearly_by_strategy,
    aes(x = x_pos, y = n, fill = therapy_strategy),
    width = 0.32
  ) +
  geom_text(
    data = ici_yearly_total,
    aes(
      x = x_pos - 0.05,
      y = ici_total + pmax(800, ici_total * 0.04),
      label = format_n(ici_total)
    ),
    vjust = 0,
    hjust = 0,
    angle = 45,
    size = 12 / ggplot2::.pt,
    color = "black"
  ) +
  scale_fill_manual(
    values = annual_strategy_color_map,
    name = "Therapy strategy",
    drop = FALSE,
    guide = guide_legend(order = 2, nrow = 2, byrow = TRUE)
  ) +
  scale_x_continuous(
    breaks = annual_trend_years,
    labels = annual_trend_years,
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

print(p_faers_solid_annual_trend_bar)

### 图2：导出PDF####

ggsave(
  filename = file.path(script_results_path, "FAERS实体瘤_年度报告与ICI趋势图.pdf"),
  plot = p_faers_solid_annual_trend_bar,
  width = 11,
  height = 3,
  device = cairo_pdf
)


## 6. 图3：ICI报告SOC景观图####

### 图3：数据准备####

if (!exists("ici_landscape_soc_abbrev_map")) {
  ici_landscape_soc_abbrev_map <- build_soc_abbrev_map()
}

if (!all(c("histologic_color_key", "tcga_color_key") %in% names(ici_landscape_column_meta))) {
  cancer_color_results <- build_hierarchical_cancer_color_maps(ici_landscape_column_meta)
  ici_landscape_column_meta <- cancer_color_results$column_meta
  ici_landscape_primary_color_map <- cancer_color_results$primary_color_map
  ici_landscape_histologic_color_map <- cancer_color_results$histologic_color_map
  ici_landscape_tcga_color_map <- cancer_color_results$tcga_color_map
  ici_landscape_color_mapping <- cancer_color_results$color_mapping
}

if (!exists("ici_landscape_gtex_mapping") || is.null(ici_landscape_gtex_mapping)) {
  ici_landscape_gtex_mapping <- read_gtex_landscape_mapping(gtex_mapping_path)
}

if (!"gtex_tissue_plot" %in% names(ici_landscape_column_meta)) {
  ici_landscape_column_meta <- add_gtex_mapping_to_landscape_column_meta(
    column_meta = ici_landscape_column_meta,
    gtex_mapping_df = ici_landscape_gtex_mapping
  )
}

if (!exists("ici_landscape_gtex_color_map") || is.null(ici_landscape_gtex_color_map)) {
  ici_landscape_gtex_color_map <- build_gtex_tissue_color_map(
    column_meta = ici_landscape_column_meta,
    primary_color_map = ici_landscape_primary_color_map
  )
}

ici_landscape_soc_summary <- ensure_soc_display_column(
  ici_landscape_soc_summary,
  soc_abbrev_map = ici_landscape_soc_abbrev_map
)
ici_landscape_soc_long <- ensure_soc_display_column(
  ici_landscape_soc_long,
  soc_abbrev_map = ici_landscape_soc_abbrev_map
)

ici_landscape_primary_annotation <- factor(
  ici_landscape_column_meta$primary_organ_plot,
  levels = unique(ici_landscape_column_meta$primary_organ_plot)
)
ici_landscape_histologic_annotation <- factor(
  ici_landscape_column_meta$histologic_color_key,
  levels = unique(ici_landscape_column_meta$histologic_color_key)
)
ici_landscape_tcga_annotation <- factor(
  ici_landscape_column_meta$tcga_color_key,
  levels = unique(ici_landscape_column_meta$tcga_color_key)
)
ici_landscape_gtex_annotation <- factor(
  ici_landscape_column_meta$gtex_tissue_plot,
  levels = unique(ici_landscape_column_meta$gtex_tissue_plot)
)
ici_landscape_column_split <- factor(
  ici_landscape_column_meta$primary_organ_plot,
  levels = unique(ici_landscape_column_meta$primary_organ_plot)
)
ici_landscape_column_order <- as.character(ici_landscape_column_meta$column_id)
ici_landscape_primary_legend_data <- ici_landscape_column_meta %>%
  distinct(primary_organ_plot, primary_total) %>%
  arrange(desc(primary_total), primary_organ_plot) %>%
  mutate(
    primary_label = if_else(
      primary_organ_plot == "Unknown/Others",
      "Unknown/Others",
      paste0(primary_organ_plot, "(", format_n(primary_total), ")")
    ),
    primary_color = unname(ici_landscape_primary_color_map[primary_organ_plot])
  )

if (!("Unknown/Others" %in% ici_landscape_primary_legend_data$primary_organ_plot)) {
  ici_landscape_primary_legend_data <- bind_rows(
    ici_landscape_primary_legend_data,
    tibble(
      primary_organ_plot = "Unknown/Others",
      primary_total = NA_real_,
      primary_label = "Unknown/Others",
      primary_color = "#BDBDBD"
    )
  )
}

ici_landscape_primary_legend_ncol <- 5L

ici_landscape_soc_summary_plot <- ici_landscape_soc_summary %>%
  mutate(
    SOC_full = as.character(SOC_full),
    SOC_display = as.character(SOC_display),
    SOC_plot_label = format_soc_abbrev_label(SOC_full, ici_landscape_soc_abbrev_map)
  ) %>%
  arrange(desc(total_reports), SOC_full)
ici_landscape_soc_order <- ici_landscape_soc_summary_plot$SOC_display
ici_landscape_soc_plot_order <- ici_landscape_soc_summary_plot$SOC_plot_label

ici_landscape_heatmap_matrix <- ici_landscape_soc_long %>%
  mutate(
    SOC_full = as.character(SOC_full),
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

rownames(ici_landscape_heatmap_matrix) <- ici_landscape_heatmap_matrix$SOC_display
ici_landscape_heatmap_matrix$SOC_display <- NULL
ici_landscape_heatmap_matrix <- as.matrix(
  ici_landscape_heatmap_matrix[ici_landscape_soc_order, ici_landscape_column_order, drop = FALSE]
)
rownames(ici_landscape_heatmap_matrix) <- ici_landscape_soc_plot_order
storage.mode(ici_landscape_heatmap_matrix) <- "double"
ici_landscape_heatmap_max_pct <- compute_heatmap_upper_bound(ici_landscape_heatmap_matrix)

ici_landscape_report_cases_vector <- ici_landscape_column_meta$report_cases
names(ici_landscape_report_cases_vector) <- ici_landscape_column_order

ici_landscape_soc_total_vector <- ici_landscape_soc_summary_plot$total_reports
names(ici_landscape_soc_total_vector) <- ici_landscape_soc_plot_order

ici_landscape_heatmap_breaks <- build_heatmap_color_breaks(
  mat = ici_landscape_heatmap_matrix,
  upper_bound = ici_landscape_heatmap_max_pct
)
ici_landscape_heatmap_break_labels <- ifelse(
  abs(ici_landscape_heatmap_breaks - round(ici_landscape_heatmap_breaks)) < 1e-8,
  as.character(round(ici_landscape_heatmap_breaks)),
  scales::number(ici_landscape_heatmap_breaks, accuracy = 0.1, trim = TRUE)
)
ici_landscape_heatmap_legend_breaks <- c(0, 25, 50, 75, 100)
ici_landscape_heatmap_legend_labels <- as.character(ici_landscape_heatmap_legend_breaks)

ici_landscape_soc_bar_breaks <- pretty(c(0, max(ici_landscape_soc_total_vector, na.rm = TRUE)), n = 4)
ici_landscape_soc_bar_breaks <- unique(
  ici_landscape_soc_bar_breaks[
    ici_landscape_soc_bar_breaks >= 0 &
      ici_landscape_soc_bar_breaks <= max(ici_landscape_soc_total_vector, na.rm = TRUE)
  ]
)
if (length(ici_landscape_soc_bar_breaks) < 2) {
  ici_landscape_soc_bar_breaks <- unique(c(0, max(ici_landscape_soc_total_vector, na.rm = TRUE)))
}

ici_landscape_heatmap_color_fun <- circlize::colorRamp2(
  ici_landscape_heatmap_breaks,
  grDevices::colorRampPalette(
    c("#FFFFFF", "#EAEAEA", "#D0D0D0", "#A8A8A8", "#737373", "#303030", "#000000")
  )(length(ici_landscape_heatmap_breaks))
)
ici_landscape_row_names_max_width <- ComplexHeatmap::max_text_width(
  ici_landscape_soc_plot_order,
  gp = grid::gpar(fontsize = 13, col = "black", fontfamily = "sans")
) + grid::unit(4, "mm")

### 图3：绘图代码####

p_faers_solid_ici_landscape_top_annotation <- ComplexHeatmap::HeatmapAnnotation(
  `Report cases per quarter` = ComplexHeatmap::anno_barplot(
    ici_landscape_report_cases_vector,
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
  `Primary organ` = ici_landscape_primary_annotation,
  `Histologic category` = ici_landscape_histologic_annotation,
  `TCGA` = ici_landscape_tcga_annotation,
  `GTEx` = ici_landscape_gtex_annotation,
  col = list(
    `Primary organ` = ici_landscape_primary_color_map,
    `Histologic category` = ici_landscape_histologic_color_map,
    `TCGA` = ici_landscape_tcga_color_map,
    `GTEx` = ici_landscape_gtex_color_map
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

p_faers_solid_ici_landscape_soc_bar <- ComplexHeatmap::rowAnnotation(
  `SOC reports` = ComplexHeatmap::anno_barplot(
    ici_landscape_soc_total_vector,
    gp = grid::gpar(fill = "#5C6969", col = NA),
    border = TRUE,
    baseline = 0,
    bar_width = 0.56,
    width = grid::unit(32, "mm"),
    axis_param = list(
      side = "bottom",
      at = ici_landscape_soc_bar_breaks,
      labels = scales::comma(ici_landscape_soc_bar_breaks),
      gp = grid::gpar(fontsize = 13, col = "black", fontfamily = "sans"),
      labels_rot = 90
    )
  ),
  show_annotation_name = FALSE,
  annotation_name_side = "bottom",
  annotation_name_rot = 0,
  annotation_name_gp = grid::gpar(fontsize = 15, col = "black", fontfamily = "sans")
)

p_faers_solid_ici_landscape_heatmap <- ComplexHeatmap::Heatmap(
  ici_landscape_heatmap_matrix,
  name = "% of report cases",
  col = ici_landscape_heatmap_color_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  row_names_side = "right",
  row_names_gp = grid::gpar(fontsize = 13, col = "black", fontfamily = "sans"),
  row_names_max_width = ici_landscape_row_names_max_width,
  row_title = NULL,
  column_title = NULL,
  height = grid::unit(58, "mm"),
  top_annotation = p_faers_solid_ici_landscape_top_annotation,
  right_annotation = p_faers_solid_ici_landscape_soc_bar,
  show_heatmap_legend = FALSE,
  rect_gp = grid::gpar(col = "#F2F2F2", lwd = 0.25),
  border = "#4D4D4D",
  use_raster = TRUE,
  raster_quality = 3
)

p_faers_solid_ici_landscape_heatmap_legend <- ComplexHeatmap::Legend(
  title = "% of report cases",
  col_fun = ici_landscape_heatmap_color_fun,
  at = ici_landscape_heatmap_legend_breaks,
  labels = ici_landscape_heatmap_legend_labels,
  direction = "horizontal",
  legend_width = grid::unit(42, "mm"),
  grid_height = grid::unit(4, "mm"),
  title_position = "topcenter",
  title_gp = grid::gpar(fontsize = 15, col = "black", fontfamily = "sans"),
  labels_gp = grid::gpar(fontsize = 13, col = "black", fontfamily = "sans")
)

p_faers_solid_ici_landscape_primary_legend <- ComplexHeatmap::Legend(
  title = "Primary organ",
  labels = ici_landscape_primary_legend_data$primary_label,
  legend_gp = grid::gpar(fill = ici_landscape_primary_legend_data$primary_color, col = NA),
  labels_gp = grid::gpar(fontsize = 12, col = "black", fontfamily = "sans"),
  title_gp = grid::gpar(fontsize = 15, col = "black", fontface = "bold", fontfamily = "sans"),
  ncol = ici_landscape_primary_legend_ncol,
  by_row = TRUE,
  gap = grid::unit(2.5, "mm"),
  row_gap = grid::unit(1.5, "mm"),
  grid_width = grid::unit(4.2, "mm"),
  grid_height = grid::unit(4.2, "mm")
)

ici_landscape_main_viewport_y <- grid::unit(0.62, "npc")
ici_landscape_main_viewport_height <- grid::unit(0.60, "npc")
ici_landscape_heatmap_padding <- grid::unit(c(3, 14, 6, 3), "mm")
ici_landscape_heatmap_legend_y <- grid::unit(0.45, "npc")
ici_landscape_heatmap_legend_height <- grid::unit(46, "mm")
ici_landscape_soc_bar_title_x <- grid::unit(0.88, "npc")
ici_landscape_soc_bar_title_y <- grid::unit(0.29, "npc")
ici_landscape_primary_legend_y <- grid::unit(0.04, "npc")

p_faers_solid_ici_landscape_heatmap_grob <- grid::grid.grabExpr(
  {
    ComplexHeatmap::draw(
      p_faers_solid_ici_landscape_heatmap,
      heatmap_legend_side = "left",
      annotation_legend_side = "bottom",
      merge_legends = FALSE,
      padding = ici_landscape_heatmap_padding
    )
  },
  width = 10,
  height = 5.4
)

p_faers_solid_ici_landscape_heatmap_legend_grob <- grid::grid.grabExpr(
  {
    ComplexHeatmap::draw(
      p_faers_solid_ici_landscape_heatmap_legend,
      x = grid::unit(0.5, "npc"),
      y = grid::unit(0.5, "npc"),
      just = c("center", "center")
    )
  },
  width = 1.4,
  height = 3
)

p_faers_solid_ici_landscape_primary_legend_grob <- grid::grid.grabExpr(
  {
    ComplexHeatmap::draw(
      p_faers_solid_ici_landscape_primary_legend,
      x = grid::unit(0, "npc"),
      y = grid::unit(0, "npc"),
      just = c("left", "bottom")
    )
  },
  width = 15,
  height = 2
)

p_faers_solid_ici_reporting_landscape <- grid::grid.grabExpr(
  {
    grid::pushViewport(
      grid::viewport(
        x = grid::unit(0.5, "npc"),
        y = ici_landscape_main_viewport_y,
        width = grid::unit(1, "npc"),
        height = ici_landscape_main_viewport_height,
        just = c("center", "center")
      )
    )
    grid::grid.draw(p_faers_solid_ici_landscape_heatmap_grob)
    grid::popViewport()

    grid::pushViewport(
      grid::viewport(
        x = grid::unit(0.055, "npc"),
        y = ici_landscape_heatmap_legend_y,
        width = grid::unit(18, "mm"),
        height = ici_landscape_heatmap_legend_height,
        just = c("left", "center"),
        angle = 90
      )
    )
    grid::grid.draw(p_faers_solid_ici_landscape_heatmap_legend_grob)
    grid::popViewport()

    grid::grid.text(
      "SOC reports",
      x = ici_landscape_soc_bar_title_x,
      y = ici_landscape_soc_bar_title_y,
      just = c("left", "center"),
      gp = grid::gpar(fontsize = 15, col = "black", fontfamily = "sans")
    )

    grid::pushViewport(
      grid::viewport(
        x = grid::unit(0.085, "npc"),
        y = ici_landscape_primary_legend_y,
        width = grid::unit(0.88, "npc"),
        height = grid::unit(0.18, "npc"),
        just = c("left", "bottom")
      )
    )
    grid::grid.draw(p_faers_solid_ici_landscape_primary_legend_grob)
    grid::popViewport()
  },
  width = 11,
  height = 10
)

grid::grid.newpage()
grid::grid.draw(p_faers_solid_ici_reporting_landscape)

### 图3：导出PDF####

grDevices::cairo_pdf(
  file = file.path(script_results_path, "FAERS实体瘤_ICI报告SOC景观图.pdf"),
  width = 11,
  height = 7,
  family = "sans"
)
grid::grid.newpage()
grid::grid.draw(p_faers_solid_ici_reporting_landscape)
grDevices::dev.off()


## 图3：FAERS三层结构导出Excel####

cancer_mapping_path <- file.path(results_path, "癌种映射.xlsx")
valid_organ_histology <- build_valid_organ_histology_mapping(cancer_mapping_path)
faers_hierarchy_export_path <- file.path(script_results_path, "图3_三层结构.xlsx")

faers_hierarchy_export <- build_three_layer_hierarchy_export(
  all_data = faers_solid_tumor_final,
  ici_data = faers_solid_tumor_final %>% filter(has_ici == TRUE),
  valid_organ_histology = valid_organ_histology
)

writexl::write_xlsx(
  faers_hierarchy_export,
  path = faers_hierarchy_export_path
)

message("\nFAERS图3三层结构已导出: ", faers_hierarchy_export_path)


## JADER：三层结构导出Excel####

jader_solid_path <- file.path(output_path, "jader_solid.RData")
jader_solid_tumor_final <- load_named_rdata_object(jader_solid_path, "jader_solid")
jader_hierarchy_export_path <- file.path(script_results_path, "JADER图3_三层结构.xlsx")

jader_hierarchy_export <- build_three_layer_hierarchy_export(
  all_data = jader_solid_tumor_final,
  ici_data = jader_solid_tumor_final %>% filter(has_ici == TRUE),
  valid_organ_histology = valid_organ_histology
)

writexl::write_xlsx(
  jader_hierarchy_export,
  path = jader_hierarchy_export_path
)

message("已加载 jader_solid 用于三层结构导出，样本数: ", nrow(jader_solid_tumor_final))
message("JADER图3三层结构已导出: ", jader_hierarchy_export_path)


## 7. 保存结果####

summary_output_path <- file.path(script_results_path, "FAERS实体瘤描述性分析汇总.xlsx")
writexl::write_xlsx(
  list(
    "筛选流程" = filter_summary,
    "原发器官分布" = primary_organ_summary,
    "组织学分类分布" = histologic_summary,
    "癌种实体分布" = assigned_type_summary,
    "ICI class分布" = ici_type_summary,
    "ICI agent分布" = ici_agent_summary,
    "ICI regimen分布" = ici_regimen_summary,
    "Therapy strategy分布" = ici_strategy_summary,
    "ICI桑基图路径统计" = ici_sankey_flow_summary,
    "年度实体瘤总数" = solid_yearly_total,
    "年度ICI_按strategy堆叠" = ici_yearly_by_strategy,
    "年度ICI总数" = ici_yearly_total,
    "第三图_列元数据" = ici_landscape_column_meta,
    "第三图_SOC百分比长表" = ici_landscape_soc_long,
    "第三图_SOC总报告数" = ici_landscape_soc_summary,
    "第三图_GTEx映射" = ici_landscape_gtex_mapping,
    "第三图_颜色映射" = ici_landscape_color_mapping
  ),
  path = summary_output_path
)

save(
  faers_solid_tumor_final,
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
  ici_sankey_node_color_map,
  annual_trend_years,
  annual_strategy_color_map,
  solid_yearly_total,
  ici_yearly_by_strategy,
  ici_yearly_total,
  faers_solid_ici_landscape,
  ici_landscape_column_meta,
  ici_landscape_soc_long,
  ici_landscape_soc_summary,
  ici_landscape_soc_abbrev_map,
  ici_landscape_gtex_mapping,
  ici_landscape_gtex_color_map,
  ici_landscape_gtex_color_mapping,
  ici_landscape_color_mapping,
  ici_landscape_primary_color_map,
  ici_landscape_histologic_color_map,
  ici_landscape_tcga_color_map,
  ici_landscape_primary_annotation,
  ici_landscape_histologic_annotation,
  ici_landscape_tcga_annotation,
  ici_landscape_gtex_annotation,
  ici_landscape_column_split,
  ici_landscape_primary_legend_data,
  ici_landscape_primary_legend_ncol,
  ici_landscape_soc_summary_plot,
  ici_landscape_soc_plot_order,
  ici_landscape_heatmap_matrix,
  ici_landscape_report_cases_vector,
  ici_landscape_soc_total_vector,
  ici_landscape_heatmap_breaks,
  ici_landscape_heatmap_break_labels,
  ici_landscape_heatmap_legend_breaks,
  ici_landscape_heatmap_legend_labels,
  ici_landscape_soc_bar_breaks,
  ici_landscape_heatmap_color_fun,
  ici_landscape_row_names_max_width,
  valid_organ_histology,
  faers_hierarchy_export,
  jader_hierarchy_export,
  p_faers_solid_sankey,
  p_faers_solid_annual_trend_bar,
  p_faers_solid_ici_landscape_top_annotation,
  p_faers_solid_ici_landscape_heatmap,
  p_faers_solid_ici_landscape_soc_bar,
  p_faers_solid_ici_landscape_primary_legend,
  p_faers_solid_ici_reporting_landscape,
  file = file.path(script_results_path, "FAERS实体瘤描述性分析对象.RData")
)

old_sankey_path <- file.path(script_results_path, "FAERS实体瘤_癌种层级桑基图.pdf")
old_bar_path <- file.path(script_results_path, "FAERS实体瘤_原发器官柱状图.pdf")
if (file.exists(old_sankey_path)) {
  invisible(file.remove(old_sankey_path))
}
if (file.exists(old_bar_path)) {
  invisible(file.remove(old_bar_path))
}

message("\n结果文件已保存:")
message("  汇总表: ", summary_output_path)
message("  对象: ", file.path(script_results_path, "FAERS实体瘤描述性分析对象.RData"))
message("  桑基图: ", file.path(script_results_path, "FAERS实体瘤_ICI情况桑基图.pdf"))
message("  柱状图: ", file.path(script_results_path, "FAERS实体瘤_年度报告与ICI趋势图.pdf"))
message("  第三张图: ", file.path(script_results_path, "FAERS实体瘤_ICI报告SOC景观图.pdf"))
message("  FAERS三层结构: ", faers_hierarchy_export_path)
message("  JADER三层结构: ", jader_hierarchy_export_path)

message("\n========== FAERS实体瘤描述性分析完成 ==========")
