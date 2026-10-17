# 7.1 TCGA基础数据准备 ####
# 用途：基于6.1最终癌种分群结果，建立模式群-原发部位-TCGA项目-GTEx器官映射，
#       并整理GTEx正常组织、TCGA正常组织、TCGA肿瘤组织三类表达矩阵及相应样本注释，
#       最终跨平台统一基因维度（以Ensembl ID去版本号为桥梁），输出gene symbol统一的矩阵 ####

rm(list = ls())
gc()

# 1. 加载配置 ####
source("Code/config.R")


# 2. 参数与路径 ####
script_name <- "7.1TCGA基础数据准备"
script_results_path <- file.path(results_path, script_name)
table_output_path <- file.path(script_results_path, "tables")
rdata_output_path <- file.path(script_results_path, "rdata")
log_output_path <- file.path(script_results_path, "logs")

dir.create(script_results_path, recursive = TRUE, showWarnings = FALSE)
dir.create(table_output_path, recursive = TRUE, showWarnings = FALSE)
dir.create(rdata_output_path, recursive = TRUE, showWarnings = FALSE)
dir.create(log_output_path, recursive = TRUE, showWarnings = FALSE)

pattern_group_assignment_file <- file.path(
  results_path,
  "6.1癌种NMF分群分析_rank5",
  "FAERS_NMF_WH_matrix.xlsx"
)

gtex_sample_annotation_file <- file.path(
  data_path,
  "TCGA_GTEx",
  "GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt"
)

gtex_expression_file <- file.path(
  data_path,
  "TCGA_GTEx",
  "GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_tpm.gct.gz"
)

mapping_output_xlsx <- file.path(
  table_output_path,
  "模式群-原发部位-TCGA项目-GTEx器官映射表.xlsx"
)

gtex_tissue_summary_output_xlsx <- file.path(
  table_output_path,
  "GTEx器官可用性汇总表.xlsx"
)

gtex_sample_annotation_output_xlsx <- file.path(
  table_output_path,
  "gtex_normal_sample_annotation.xlsx"
)

gtex_gene_annotation_output_xlsx <- file.path(
  table_output_path,
  "gtex_expression_gene_annotation.xlsx"
)

gtex_sample_annotation_output_rdata <- file.path(
  rdata_output_path,
  "gtex_normal_sample_annotation.Rdata"
)

gtex_expression_matrix_output_rdata <- file.path(
  rdata_output_path,
  "gtex_expression_normal_matrix.Rdata"
)

gtex_gene_annotation_output_rdata <- file.path(
  rdata_output_path,
  "gtex_expression_gene_annotation.Rdata"
)

run_timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_file <- file.path(log_output_path, paste0("7.1TCGA基础数据准备_", run_timestamp, ".log"))
script_file <- file.path("Code", "7.1TCGA基础数据准备.R")

pattern_group_order <- c(
  "Ocular",
  "Hematologic",
  "Neurologic",
  "Systemic-like",
  "Visceral-like"
)


# 3. 日志与通用辅助函数 ####

## 函数：log_message ####
## 功能：同时向控制台和日志文件写入消息 ####
## 输入：message_text，待写入日志的文本 ####
## 输出：无返回值 ####
log_message <- function(message_text) {
  timestamp_text <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  final_text <- paste0("[", timestamp_text, "] ", message_text)

  message(final_text)
  cat(final_text, "\n", file = log_file, append = TRUE)
}

## 函数：assert_required_files ####
## 功能：检查必要输入文件是否存在，不存在则终止脚本 ####
## 输入：file_paths，文件路径向量；file_labels，文件说明向量 ####
## 输出：无返回值 ####
assert_required_files <- function(file_paths, file_labels) {
  missing_index <- which(!file.exists(file_paths))

  if (length(missing_index) > 0) {
    missing_messages <- paste0(file_labels[missing_index], ": ", file_paths[missing_index])
    stop("缺少必要输入文件:\n", paste(missing_messages, collapse = "\n"))
  }
}

## 函数：are_result_files_up_to_date ####
## 功能：检查一组结果文件是否都已存在，且生成时间不早于依赖文件 ####
## 输入：result_files，结果文件路径向量；dependency_files，依赖文件路径向量 ####
## 输出：逻辑值，TRUE表示结果文件可直接复用 ####
are_result_files_up_to_date <- function(result_files, dependency_files) {
  result_files <- unique(result_files)
  dependency_files <- unique(dependency_files)

  if (length(result_files) == 0 || any(!file.exists(result_files))) {
    return(FALSE)
  }

  dependency_files <- dependency_files[file.exists(dependency_files)]

  if (length(dependency_files) == 0) {
    return(TRUE)
  }

  min_result_mtime <- min(file.info(result_files)$mtime)
  max_dependency_mtime <- max(file.info(dependency_files)$mtime)

  min_result_mtime >= max_dependency_mtime
}

## 函数：load_single_object_from_rdata ####
## 功能：从单对象Rdata文件中读取指定对象 ####
## 输入：rdata_file，Rdata文件路径；object_name，对象名 ####
## 输出：读取出的对象 ####
load_single_object_from_rdata <- function(rdata_file, object_name) {
  temp_env <- new.env(parent = emptyenv())
  load(rdata_file, envir = temp_env)

  if (!exists(object_name, envir = temp_env, inherits = FALSE)) {
    stop("Rdata文件中不存在目标对象: ", object_name, " | ", rdata_file)
  }

  get(object_name, envir = temp_env, inherits = FALSE)
}

## 函数：try_reuse_cached_gtex_expression_result ####
## 功能：在当前样本集合未变化时，复用已保存的GTEx表达矩阵和基因注释 ####
## 输入：current_sample_ids，当前目标样本ID向量；sample_annotation_rdata_file，样本注释Rdata路径；expression_matrix_rdata_file，表达矩阵Rdata路径；gene_annotation_rdata_file，基因注释Rdata路径 ####
## 输出：列表，包含can_reuse、reuse_message、expression_matrix和gene_annotation_df ####
try_reuse_cached_gtex_expression_result <- function(
    current_sample_ids,
    sample_annotation_rdata_file,
    expression_matrix_rdata_file,
    gene_annotation_rdata_file
) {
  required_cache_files <- c(
    sample_annotation_rdata_file,
    expression_matrix_rdata_file,
    gene_annotation_rdata_file
  )

  if (any(!file.exists(required_cache_files))) {
    return(list(
      can_reuse = FALSE,
      reuse_message = "缓存Rdata文件不完整",
      expression_matrix = NULL,
      gene_annotation_df = NULL
    ))
  }

  cached_sample_annotation_df <- tryCatch(
    load_single_object_from_rdata(
      rdata_file = sample_annotation_rdata_file,
      object_name = "gtex_final_sample_annotation_df"
    ),
    error = function(e) e
  )

  if (inherits(cached_sample_annotation_df, "error")) {
    return(list(
      can_reuse = FALSE,
      reuse_message = paste0("读取缓存样本注释失败: ", conditionMessage(cached_sample_annotation_df)),
      expression_matrix = NULL,
      gene_annotation_df = NULL
    ))
  }

  if (!"sample_id" %in% colnames(cached_sample_annotation_df)) {
    return(list(
      can_reuse = FALSE,
      reuse_message = "缓存样本注释缺少sample_id列",
      expression_matrix = NULL,
      gene_annotation_df = NULL
    ))
  }

  current_sample_ids <- as.character(current_sample_ids)
  cached_sample_ids <- as.character(cached_sample_annotation_df$sample_id)

  if (!setequal(current_sample_ids, cached_sample_ids)) {
    return(list(
      can_reuse = FALSE,
      reuse_message = "缓存样本集合与当前映射结果不一致",
      expression_matrix = NULL,
      gene_annotation_df = NULL
    ))
  }

  cached_expression_matrix <- tryCatch(
    load_single_object_from_rdata(
      rdata_file = expression_matrix_rdata_file,
      object_name = "gtex_expression_normal_matrix"
    ),
    error = function(e) e
  )

  if (inherits(cached_expression_matrix, "error")) {
    return(list(
      can_reuse = FALSE,
      reuse_message = paste0("读取缓存表达矩阵失败: ", conditionMessage(cached_expression_matrix)),
      expression_matrix = NULL,
      gene_annotation_df = NULL
    ))
  }

  cached_gene_annotation_df <- tryCatch(
    load_single_object_from_rdata(
      rdata_file = gene_annotation_rdata_file,
      object_name = "gtex_gene_annotation_df"
    ),
    error = function(e) e
  )

  if (inherits(cached_gene_annotation_df, "error")) {
    return(list(
      can_reuse = FALSE,
      reuse_message = paste0("读取缓存基因注释失败: ", conditionMessage(cached_gene_annotation_df)),
      expression_matrix = NULL,
      gene_annotation_df = NULL
    ))
  }

  cached_expression_matrix <- as.matrix(cached_expression_matrix)

  if (is.null(colnames(cached_expression_matrix))) {
    return(list(
      can_reuse = FALSE,
      reuse_message = "缓存表达矩阵缺少样本列名",
      expression_matrix = NULL,
      gene_annotation_df = NULL
    ))
  }

  if (!all(current_sample_ids %in% colnames(cached_expression_matrix))) {
    return(list(
      can_reuse = FALSE,
      reuse_message = "缓存表达矩阵缺少当前目标样本列",
      expression_matrix = NULL,
      gene_annotation_df = NULL
    ))
  }

  cached_expression_matrix <- cached_expression_matrix[, current_sample_ids, drop = FALSE]

  if (nrow(cached_gene_annotation_df) != nrow(cached_expression_matrix)) {
    return(list(
      can_reuse = FALSE,
      reuse_message = "缓存基因注释与表达矩阵行数不一致",
      expression_matrix = NULL,
      gene_annotation_df = NULL
    ))
  }

  reuse_message <- if (identical(current_sample_ids, cached_sample_ids)) {
    "缓存样本顺序与当前结果一致"
  } else {
    "缓存样本集合一致，已按当前样本顺序重排列"
  }

  list(
    can_reuse = TRUE,
    reuse_message = reuse_message,
    expression_matrix = cached_expression_matrix,
    gene_annotation_df = cached_gene_annotation_df
  )
}

## 函数：collapse_unique_values ####
## 功能：将向量中的唯一非缺失值合并为单个分号分隔字符串 ####
## 输入：input_vector，待合并向量 ####
## 输出：单个字符串；若无有效值则返回NA ####
collapse_unique_values <- function(input_vector) {
  output_vector <- unique(as.character(input_vector))
  output_vector <- output_vector[!is.na(output_vector) & output_vector != ""]

  if (length(output_vector) == 0) {
    return(NA_character_)
  }

  paste(output_vector, collapse = "; ")
}


# 4. 映射规则辅助函数 ####

## 函数：read_pattern_group_assignment ####
## 功能：读取6.1 NMF结果中的原发部位-模式群对应表 ####
## 输入：xlsx_file，6.1输出的WH矩阵Excel路径 ####
## 输出：包含primary_organ和pattern_group的数据框 ####
read_pattern_group_assignment <- function(xlsx_file) {
  readxl::read_xlsx(
    path = xlsx_file,
    sheet = "Pattern group assignment"
  ) %>%
    dplyr::rename(
      primary_organ = cancer_type,
      pattern_group = `Pattern group`
    ) %>%
    dplyr::mutate(
      primary_organ = as.character(primary_organ),
      pattern_group = factor(as.character(pattern_group), levels = pattern_group_order)
    ) %>%
    dplyr::arrange(pattern_group, primary_organ) %>%
    dplyr::mutate(pattern_group = as.character(pattern_group))
}

## 函数：build_tcga_project_rule_table ####
## 功能：为7.1当前纳入的primary_organ建立TCGA项目映射规则 ####
## 输入：无 ####
## 输出：包含primary_organ、tcga_project和tcga_code的数据框 ####
build_tcga_project_rule_table <- function() {
  tibble::tribble(
    ~primary_organ,       ~tcga_project, ~tcga_code,
    "Adrenal gland",      "TCGA-ACC",    "ACC",
    "Biliary tract",      "TCGA-CHOL",   "CHOL",
    "Bladder",            "TCGA-BLCA",   "BLCA",
    "Bowel",              "TCGA-COAD",   "COAD",
    "Bowel",              "TCGA-READ",   "READ",
    "Brain",              "TCGA-GBM",    "GBM",
    "Brain",              "TCGA-LGG",    "LGG",
    "Breast",             "TCGA-BRCA",   "BRCA",
    "Cervix",             "TCGA-CESC",   "CESC",
    "Esophagus/stomach",  "TCGA-ESCA",   "ESCA",
    "Esophagus/stomach",  "TCGA-STAD",   "STAD",
    "Eye",                "TCGA-UVM",    "UVM",
    "Head/neck",          "TCGA-HNSC",   "HNSC",
    "Kidney",             "TCGA-KICH",   "KICH",
    "Kidney",             "TCGA-KIRC",   "KIRC",
    "Kidney",             "TCGA-KIRP",   "KIRP",
    "Liver",              "TCGA-LIHC",   "LIHC",
    "Lung",               "TCGA-LUAD",   "LUAD",
    "Lung",               "TCGA-LUSC",   "LUSC",
    "Ovary",              "TCGA-OV",     "OV",
    "Pancreas",           "TCGA-PAAD",   "PAAD",
    "PNS",                "TCGA-PCPG",   "PCPG",
    "Pleura",             "TCGA-MESO",   "MESO",
    "Prostate",           "TCGA-PRAD",   "PRAD",
    "Skin",               "TCGA-SKCM",   "SKCM",
    "Soft tissue",        "TCGA-SARC",   "SARC",
    "Thymus",             "TCGA-THYM",   "THYM",
    "Thyroid",            "TCGA-THCA",   "THCA",
    "Uterus",             "TCGA-UCEC",   "UCEC",
    "Uterus",             "TCGA-UCS",    "UCS"
  )
}

## 函数：build_gtex_bulk_tissue_candidate_table ####
## 功能：将GTEx注释中所有可用于本项目的bulk tissue整理为候选表 ####
## 输入：tissue_summary_df，GTEx组织汇总表 ####
## 输出：GTEx bulk tissue候选表 ####
build_gtex_bulk_tissue_candidate_table <- function(tissue_summary_df) {
  tissue_summary_df %>%
    dplyr::filter(n_bulk_samples > 0) %>%
    dplyr::mutate(
      candidate_status = "available_bulk_tissue",
      tissue_granularity = dplyr::if_else(
        gtex_general_tissue == gtex_specific_tissue,
        "general_equals_specific",
        "specific_subtissue"
      )
    ) %>%
    dplyr::arrange(dplyr::desc(n_bulk_samples), gtex_general_tissue, gtex_specific_tissue)
}

## 函数：build_gtex_primary_organ_review_table ####
## 功能：基于GTEx候选组织表，逐个primary_organ定义最终纳入或排除的GTEx映射决策 ####
## 输入：无 ####
## 输出：包含逐器官审阅决策的数据框 ####
build_gtex_primary_organ_review_table <- function() {
  tibble::tribble(
    ~primary_organ,       ~decision_status, ~match_column,           ~match_value,                              ~selected_gtex_general_tissue, ~is_proxy_tissue, ~review_note,
    "Adrenal gland",      "include_direct", "gtex_specific_tissue",  "Adrenal Gland",                          "Adrenal Gland",               FALSE,            "direct_match",
    "Biliary tract",      "no_match",       NA_character_,           NA_character_,                            NA_character_,                 FALSE,            "no_direct_bulk_tissue_in_gtex_v11",
    "Bladder",            "include_direct", "gtex_specific_tissue",  "Bladder",                                "Bladder",                     FALSE,            "direct_match",
    "Bone",               "no_match",       NA_character_,           NA_character_,                            NA_character_,                 FALSE,            "no_direct_bulk_tissue_in_gtex_v11",
    "Bowel",              "include_direct", "gtex_specific_tissue",  "Colon - Sigmoid",                        "Colon",                       FALSE,            "final_use_colorectal_normal_tissue",
    "Bowel",              "include_direct", "gtex_specific_tissue",  "Colon - Transverse",                     "Colon",                       FALSE,            "final_use_colorectal_normal_tissue",
    "Brain",              "include_direct", "gtex_general_tissue",   "Brain",                                  "Brain",                       FALSE,            "final_use_all_brain_bulk_tissues",
    "Breast",             "include_direct", "gtex_specific_tissue",  "Breast - Mammary Tissue",                "Breast",                      FALSE,            "direct_match",
    "Cervix",             "include_direct", "gtex_specific_tissue",  "Cervix - Endocervix",                    "Cervix Uteri",                FALSE,            "direct_match",
    "Cervix",             "include_direct", "gtex_specific_tissue",  "Cervix - Ectocervix",                    "Cervix Uteri",                FALSE,            "direct_match",
    "Esophagus/stomach",  "include_direct", "gtex_specific_tissue",  "Esophagus - Gastroesophageal Junction",  "Esophagus",                   FALSE,            "combined_organ_direct_match",
    "Esophagus/stomach",  "include_direct", "gtex_specific_tissue",  "Esophagus - Mucosa",                     "Esophagus",                   FALSE,            "combined_organ_direct_match",
    "Esophagus/stomach",  "include_direct", "gtex_specific_tissue",  "Esophagus - Muscularis",                 "Esophagus",                   FALSE,            "combined_organ_direct_match",
    "Esophagus/stomach",  "include_direct", "gtex_specific_tissue",  "Stomach",                                "Stomach",                     FALSE,            "combined_organ_direct_match",
    "Esophagus/stomach",  "include_direct", "gtex_specific_tissue",  "Stomach - Mucosa",                       "Stomach",                     FALSE,            "combined_organ_direct_match",
    "Esophagus/stomach",  "include_direct", "gtex_specific_tissue",  "Stomach - Muscularis",                   "Stomach",                     FALSE,            "combined_organ_direct_match",
    "Eye",                "no_match",       NA_character_,           NA_character_,                            NA_character_,                 FALSE,            "no_direct_bulk_tissue_in_gtex_v11",
    "Head/neck",          "include_proxy",  "gtex_specific_tissue",  "Minor Salivary Gland",                   "Salivary Gland",              TRUE,             "proxy_for_head_neck",
    "Kidney",             "include_direct", "gtex_specific_tissue",  "Kidney - Cortex",                        "Kidney",                      FALSE,            "final_use_kidney_cortex_bulk_tissue",
    "Liver",              "include_direct", "gtex_specific_tissue",  "Liver",                                  "Liver",                       FALSE,            "direct_match",
    "Lung",               "include_direct", "gtex_specific_tissue",  "Lung",                                   "Lung",                        FALSE,            "direct_match",
    "Ovary",              "include_direct", "gtex_specific_tissue",  "Ovary",                                  "Ovary",                       FALSE,            "direct_match",
    "Pancreas",           "include_direct", "gtex_specific_tissue",  "Pancreas",                               "Pancreas",                    FALSE,            "final_use_pancreas_bulk_tissue",
    "PNS",                "include_proxy",  "gtex_specific_tissue",  "Nerve - Tibial",                         "Nerve",                       TRUE,             "proxy_for_peripheral_nervous_system",
    "Pleura",             "no_match",       NA_character_,           NA_character_,                            NA_character_,                 FALSE,            "no_direct_bulk_tissue_in_gtex_v11",
    "Prostate",           "include_direct", "gtex_specific_tissue",  "Prostate",                               "Prostate",                    FALSE,            "direct_match",
    "Skin",               "include_direct", "gtex_specific_tissue",  "Skin - Not Sun Exposed (Suprapubic)",    "Skin",                        FALSE,            "direct_match",
    "Skin",               "include_direct", "gtex_specific_tissue",  "Skin - Sun Exposed (Lower leg)",         "Skin",                        FALSE,            "direct_match",
    "Soft tissue",        "include_proxy",  "gtex_specific_tissue",  "Muscle - Skeletal",                      "Muscle",                      TRUE,             "proxy_for_soft_tissue",
    "Thymus",             "no_match",       NA_character_,           NA_character_,                            NA_character_,                 FALSE,            "no_direct_bulk_tissue_in_gtex_v11",
    "Thyroid",            "include_direct", "gtex_specific_tissue",  "Thyroid",                                "Thyroid",                     FALSE,            "direct_match",
    "Uterus",             "include_direct", "gtex_specific_tissue",  "Uterus",                                 "Uterus",                      FALSE,            "direct_match",
    "Vulva/Vagina",       "include_proxy",  "gtex_specific_tissue",  "Vagina",                                 "Vagina",                      TRUE,             "proxy_for_vulva_vagina"
  )
}

## 函数：parse_excel_logical ####
## 功能：将Excel读取出的0/1、TRUE/FALSE或文本逻辑值稳健转换为逻辑型 ####
## 输入：input_vector，待转换向量 ####
## 输出：逻辑向量 ####
parse_excel_logical <- function(input_vector) {
  normalized_vector <- tolower(trimws(as.character(input_vector)))

  dplyr::case_when(
    normalized_vector %in% c("true", "t", "1", "yes", "y") ~ TRUE,
    normalized_vector %in% c("false", "f", "0", "no", "n", "") ~ FALSE,
    is.na(normalized_vector) ~ FALSE,
    TRUE ~ FALSE
  )
}

## 函数：read_or_build_gtex_primary_organ_review_table ####
## 功能：优先读取用户已修改的GTEx逐器官审阅表；若不存在则使用脚本默认规则 ####
## 输入：review_file，映射工作簿路径 ####
## 输出：GTEx逐器官审阅决策表 ####
read_or_build_gtex_primary_organ_review_table <- function(review_file) {
  default_review_df <- build_gtex_primary_organ_review_table()

  if (!file.exists(review_file)) {
    attr(default_review_df, "rule_source") <- "default_builtin_rules"
    return(default_review_df)
  }

  available_sheets <- readxl::excel_sheets(review_file)
  candidate_sheet_names <- c(
    "GTEx器官审阅表",
    "gtex_organ_review",
    "gtex_primary_organ_review"
  )
  matched_sheet_name <- candidate_sheet_names[candidate_sheet_names %in% available_sheets][[1]]

  if (is.null(matched_sheet_name) || is.na(matched_sheet_name)) {
    attr(default_review_df, "rule_source") <- "default_builtin_rules"
    return(default_review_df)
  }

  review_df <- readxl::read_xlsx(
    path = review_file,
    sheet = matched_sheet_name
  )

  required_columns <- c(
    "primary_organ",
    "decision_status",
    "match_column",
    "match_value",
    "selected_gtex_general_tissue",
    "is_proxy_tissue",
    "review_note"
  )

  missing_columns <- setdiff(required_columns, colnames(review_df))

  if (length(missing_columns) > 0) {
    stop(
      "用户修改的GTEx审阅表缺少必要列: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  output_df <- review_df %>%
    dplyr::select(dplyr::all_of(required_columns)) %>%
    dplyr::mutate(
      primary_organ = as.character(primary_organ),
      decision_status = as.character(decision_status),
      match_column = dplyr::if_else(is.na(match_column), NA_character_, as.character(match_column)),
      match_value = dplyr::if_else(is.na(match_value), NA_character_, as.character(match_value)),
      selected_gtex_general_tissue = dplyr::if_else(
        is.na(selected_gtex_general_tissue),
        NA_character_,
        as.character(selected_gtex_general_tissue)
      ),
      is_proxy_tissue = parse_excel_logical(is_proxy_tissue),
      review_note = dplyr::if_else(is.na(review_note), NA_character_, as.character(review_note))
    )

  attr(output_df, "rule_source") <- paste0("existing_mapping_workbook:", matched_sheet_name)

  output_df
}

## 函数：extract_final_gtex_mapping_rule_table ####
## 功能：从逐器官审阅表中提取最终用于筛样本的GTEx映射规则 ####
## 输入：review_df，逐器官审阅决策表 ####
## 输出：最终GTEx映射规则表 ####
extract_final_gtex_mapping_rule_table <- function(review_df) {
  review_df %>%
    dplyr::filter(decision_status %in% c("include_direct", "include_proxy")) %>%
    dplyr::transmute(
      primary_organ = primary_organ,
      match_column = match_column,
      match_value = match_value,
      selected_gtex_general_tissue = selected_gtex_general_tissue,
      is_proxy_tissue = is_proxy_tissue,
      mapping_note = review_note
    ) %>%
    dplyr::arrange(primary_organ, match_column, match_value)
}

## 函数：validate_gtex_review_coverage ####
## 功能：检查GTEx逐器官审阅表是否完整覆盖当前纳入研究的primary_organ ####
## 输入：pattern_group_df，模式群原发部位表；review_df，GTEx逐器官审阅表 ####
## 输出：无返回值；若覆盖不完整则报错 ####
validate_gtex_review_coverage <- function(pattern_group_df, review_df) {
  expected_primary_organs <- unique(as.character(pattern_group_df$primary_organ))
  reviewed_primary_organs <- unique(as.character(review_df$primary_organ))

  missing_primary_organs <- setdiff(expected_primary_organs, reviewed_primary_organs)
  extra_primary_organs <- setdiff(reviewed_primary_organs, expected_primary_organs)

  if (length(missing_primary_organs) > 0) {
    stop(
      "GTEx逐器官审阅表缺少以下primary_organ: ",
      paste(missing_primary_organs, collapse = ", ")
    )
  }

  if (length(extra_primary_organs) > 0) {
    stop(
      "GTEx逐器官审阅表包含未纳入本轮研究的primary_organ: ",
      paste(extra_primary_organs, collapse = ", ")
    )
  }
}


# 5. GTEx读取与筛选辅助函数 ####

## 函数：read_gtex_sample_annotation ####
## 功能：读取GTEx样本注释，并保留当前步骤所需列 ####
## 输入：annotation_file，GTEx样本注释文件路径 ####
## 输出：包含sample_id、SMTS和SMTSD等字段的数据框 ####
read_gtex_sample_annotation <- function(annotation_file) {
  data.table::fread(
    input = annotation_file,
    sep = "\t",
    data.table = FALSE,
    showProgress = FALSE
  ) %>%
    dplyr::select(SAMPID, SMTS, SMTSD) %>%
    dplyr::rename(
      sample_id = SAMPID,
      gtex_general_tissue = SMTS,
      gtex_specific_tissue = SMTSD
    ) %>%
    dplyr::mutate(
      sample_id = as.character(sample_id),
      gtex_general_tissue = as.character(gtex_general_tissue),
      gtex_specific_tissue = as.character(gtex_specific_tissue)
    )
}

## 函数：read_gtex_expression_header_sample_ids ####
## 功能：仅读取GTEx表达文件表头中的样本ID，用于约束后续注释筛选范围 ####
## 输入：expression_file，GTEx表达文件路径 ####
## 输出：表达矩阵表头中的样本ID字符向量 ####
read_gtex_expression_header_sample_ids <- function(expression_file) {
  header_df <- data.table::fread(
    cmd = paste("gzip -dc", shQuote(expression_file)),
    skip = 2,
    nrows = 0,
    data.table = FALSE,
    showProgress = FALSE
  )

  setdiff(colnames(header_df), c("Name", "Description"))
}

## 函数：split_gtex_annotation_by_expression_availability ####
## 功能：按GTEx表达矩阵表头样本ID拆分注释数据，区分有表达和无表达的注释样本 ####
## 输入：annotation_df，GTEx样本注释数据框；expression_sample_ids，表达矩阵表头样本ID向量 ####
## 输出：包含matched_annotation_df、unmatched_annotation_df和matching_summary_df的列表 ####
split_gtex_annotation_by_expression_availability <- function(annotation_df, expression_sample_ids) {
  annotation_with_flag_df <- annotation_df %>%
    dplyr::mutate(has_expression_in_matrix = sample_id %in% expression_sample_ids)

  matched_annotation_df <- annotation_with_flag_df %>%
    dplyr::filter(has_expression_in_matrix) %>%
    dplyr::select(-has_expression_in_matrix)

  unmatched_annotation_df <- annotation_with_flag_df %>%
    dplyr::filter(!has_expression_in_matrix) %>%
    dplyr::select(-has_expression_in_matrix)

  matching_summary_df <- tibble::tibble(
    metric = c(
      "annotation_sample_n",
      "expression_header_sample_n",
      "annotation_expression_overlap_n",
      "annotation_without_expression_n"
    ),
    value = c(
      nrow(annotation_df),
      length(expression_sample_ids),
      nrow(matched_annotation_df),
      nrow(unmatched_annotation_df)
    )
  )

  list(
    matched_annotation_df = matched_annotation_df,
    unmatched_annotation_df = unmatched_annotation_df,
    matching_summary_df = matching_summary_df
  )
}

## 函数：flag_bulk_gtex_tissues ####
## 功能：为GTEx样本标记是否属于可作为项目正常组织的bulk tissue ####
## 输入：annotation_df，GTEx样本注释数据框 ####
## 输出：新增bulk tissue标记列后的数据框 ####
flag_bulk_gtex_tissues <- function(annotation_df) {
  annotation_df %>%
    dplyr::mutate(
      is_cell_or_culture = grepl("^Cells - ", gtex_specific_tissue),
      is_microdissected_component = grepl(
        " - (Acini|Islets|Portal Tract|Hepatocyte|Lymphoid Aggregate|Mixed Cell)$",
        gtex_specific_tissue
      ),
      is_bulk_tissue_candidate = !is_cell_or_culture & !is_microdissected_component
    )
}

## 函数：summarize_gtex_tissues ####
## 功能：按GTEx组织层级汇总样本数，并统计bulk候选数量 ####
## 输入：annotation_df，带bulk标记的GTEx样本注释数据框 ####
## 输出：GTEx组织汇总表 ####
summarize_gtex_tissues <- function(annotation_df) {
  annotation_df %>%
    dplyr::group_by(gtex_general_tissue, gtex_specific_tissue) %>%
    dplyr::summarise(
      n_samples_total = dplyr::n(),
      n_bulk_samples = sum(is_bulk_tissue_candidate),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(n_bulk_samples), gtex_general_tissue, gtex_specific_tissue)
}

## 函数：validate_gtex_mapping_rules ####
## 功能：检查GTEx映射规则中的目标组织是否存在于当前bulk组织注释中 ####
## 输入：annotation_df，带bulk标记的GTEx样本注释数据框；mapping_rule_df，GTEx映射规则表 ####
## 输出：无返回值；若存在缺失规则则报错 ####
validate_gtex_mapping_rules <- function(annotation_df, mapping_rule_df) {
  bulk_annotation_df <- annotation_df %>%
    dplyr::filter(is_bulk_tissue_candidate)

  missing_rule_df <- mapping_rule_df %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      rule_exists = any(bulk_annotation_df[[match_column]] == match_value, na.rm = TRUE)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!rule_exists)

  if (nrow(missing_rule_df) > 0) {
    stop(
      "以下GTEx映射规则在当前bulk tissue注释中不存在，请先调整规则:\n",
      paste(
        paste0(
          missing_rule_df$primary_organ,
          " | ",
          missing_rule_df$match_column,
          " = ",
          missing_rule_df$match_value
        ),
        collapse = "\n"
      )
    )
  }
}

## 函数：apply_gtex_mapping_rules ####
## 功能：按预定义规则筛选进入研究的GTEx样本，并补充primary_organ和pattern_group信息 ####
## 输入：annotation_df，带bulk标记的GTEx样本注释数据框；mapping_rule_df，GTEx映射规则表；pattern_group_df，模式群映射表 ####
## 输出：最终纳入研究的GTEx样本注释数据框 ####
apply_gtex_mapping_rules <- function(annotation_df, mapping_rule_df, pattern_group_df) {
  bulk_annotation_df <- annotation_df %>%
    dplyr::filter(is_bulk_tissue_candidate)

  mapped_annotation_list <- lapply(
    seq_len(nrow(mapping_rule_df)),
    function(rule_index) {
      current_rule <- mapping_rule_df[rule_index, , drop = FALSE]

      bulk_annotation_df %>%
        dplyr::filter(.data[[current_rule$match_column[[1]]]] == current_rule$match_value[[1]]) %>%
        dplyr::mutate(
          primary_organ = current_rule$primary_organ[[1]],
          gtex_match_column = current_rule$match_column[[1]],
          gtex_match_value = current_rule$match_value[[1]],
          selected_gtex_general_tissue = current_rule$selected_gtex_general_tissue[[1]],
          is_proxy_tissue = current_rule$is_proxy_tissue[[1]],
          mapping_note = current_rule$mapping_note[[1]]
        )
    }
  )

  mapped_annotation_df <- dplyr::bind_rows(mapped_annotation_list) %>%
    dplyr::left_join(pattern_group_df, by = "primary_organ") %>%
    dplyr::mutate(
      pattern_group = factor(pattern_group, levels = pattern_group_order)
    ) %>%
    dplyr::arrange(pattern_group, primary_organ, gtex_general_tissue, gtex_specific_tissue, sample_id) %>%
    dplyr::mutate(pattern_group = as.character(pattern_group))

  duplicated_sample_df <- mapped_annotation_df %>%
    dplyr::count(sample_id, name = "mapping_n") %>%
    dplyr::filter(mapping_n > 1)

  if (nrow(duplicated_sample_df) > 0) {
    stop(
      "检测到GTEx样本被映射到多个primary_organ，请先调整规则: ",
      paste(duplicated_sample_df$sample_id, collapse = ", ")
    )
  }

  mapped_annotation_df
}

## 函数：read_gtex_expression_subset ####
## 功能：读取GTEx表达矩阵中指定样本列，并拆分为基因注释表和表达矩阵 ####
## 输入：expression_file，GTEx表达文件路径；sample_ids，目标样本ID向量 ####
## 输出：包含gene_annotation_df和expression_matrix的列表 ####
read_gtex_expression_subset <- function(expression_file, sample_ids) {
  header_df <- data.table::fread(
    cmd = paste("gzip -dc", shQuote(expression_file)),
    skip = 2,
    nrows = 0,
    data.table = FALSE,
    showProgress = FALSE
  )

  available_sample_ids <- setdiff(colnames(header_df), c("Name", "Description"))
  missing_sample_ids <- setdiff(sample_ids, available_sample_ids)

  if (length(missing_sample_ids) > 0) {
    stop(
      "以下GTEx样本在表达矩阵表头中不存在: ",
      paste(missing_sample_ids, collapse = ", ")
    )
  }

  selected_columns <- c("Name", "Description", sample_ids)

  expression_df <- data.table::fread(
    cmd = paste("gzip -dc", shQuote(expression_file)),
    skip = 2,
    select = selected_columns,
    data.table = FALSE,
    showProgress = TRUE
  )

  gene_annotation_df <- expression_df %>%
    dplyr::transmute(
      gene_id = as.character(Name),
      gene_symbol = as.character(Description)
    )

  expression_matrix <- as.matrix(expression_df[, sample_ids, drop = FALSE])
  storage.mode(expression_matrix) <- "numeric"
  rownames(expression_matrix) <- gene_annotation_df$gene_id

  list(
    gene_annotation_df = gene_annotation_df,
    expression_matrix = expression_matrix
  )
}


# 6. 读取输入数据 ####
assert_required_files(
  file_paths = c(
    pattern_group_assignment_file,
    gtex_sample_annotation_file,
    gtex_expression_file
  ),
  file_labels = c(
    "6.1模式群分配文件",
    "GTEx样本注释文件",
    "GTEx表达文件"
  )
)

log_message("========== 7.1 TCGA基础数据准备开始 ==========")
log_message(paste0("模式群文件: ", pattern_group_assignment_file))
log_message(paste0("GTEx样本注释文件: ", gtex_sample_annotation_file))
log_message(paste0("GTEx表达文件: ", gtex_expression_file))

full_result_files <- c(
  mapping_output_xlsx,
  gtex_tissue_summary_output_xlsx,
  gtex_sample_annotation_output_xlsx,
  gtex_gene_annotation_output_xlsx,
  gtex_sample_annotation_output_rdata,
  gtex_expression_matrix_output_rdata,
  gtex_gene_annotation_output_rdata
)

full_result_dependency_files <- c(
  pattern_group_assignment_file,
  gtex_sample_annotation_file,
  gtex_expression_file,
  script_file,
  mapping_output_xlsx
)

reuse_full_results <- are_result_files_up_to_date(
  result_files = full_result_files,
  dependency_files = full_result_dependency_files
)

full_result_cache_state <- list(
  is_loaded = FALSE,
  error_message = NULL
)

if (reuse_full_results) {
  tryCatch(
    {
      gtex_final_sample_annotation_df <- load_single_object_from_rdata(
        rdata_file = gtex_sample_annotation_output_rdata,
        object_name = "gtex_final_sample_annotation_df"
      )
      gtex_expression_normal_matrix <- load_single_object_from_rdata(
        rdata_file = gtex_expression_matrix_output_rdata,
        object_name = "gtex_expression_normal_matrix"
      )
      gtex_gene_annotation_df <- load_single_object_from_rdata(
        rdata_file = gtex_gene_annotation_output_rdata,
        object_name = "gtex_gene_annotation_df"
      )

      full_result_cache_state$is_loaded <- TRUE
    },
    error = function(e) {
      full_result_cache_state$error_message <- conditionMessage(e)
    }
  )
}

if (full_result_cache_state$is_loaded) {
  log_message("检测到已保存且为最新的7.1结果，直接读取缓存并跳过本次重算。")
  log_message(paste0("复用GTEx样本注释Rdata: ", gtex_sample_annotation_output_rdata))
  log_message(paste0("复用GTEx表达矩阵Rdata: ", gtex_expression_matrix_output_rdata))
  log_message(paste0("复用GTEx基因注释Rdata: ", gtex_gene_annotation_output_rdata))
  log_message(paste0("缓存GTEx样本数: ", nrow(gtex_final_sample_annotation_df)))
  log_message(
    paste0(
      "缓存GTEx表达矩阵维度: ",
      nrow(gtex_expression_normal_matrix),
      " × ",
      ncol(gtex_expression_normal_matrix)
    )
  )
} else {
  if (reuse_full_results && !is.null(full_result_cache_state$error_message)) {
    log_message(
      paste0(
        "检测到完整结果缓存，但读取失败，将重新计算: ",
        full_result_cache_state$error_message
      )
    )
  } else {
    log_message("未检测到可直接复用的完整结果，将重新整理GTEx映射和表达数据。")
  }

  pattern_group_assignment_df <- read_pattern_group_assignment(pattern_group_assignment_file)
  tcga_project_rule_df <- build_tcga_project_rule_table()

  gtex_annotation_raw_df <- read_gtex_sample_annotation(gtex_sample_annotation_file)
  gtex_expression_header_sample_ids <- read_gtex_expression_header_sample_ids(gtex_expression_file)
  gtex_annotation_matching_result <- split_gtex_annotation_by_expression_availability(
    annotation_df = gtex_annotation_raw_df,
    expression_sample_ids = gtex_expression_header_sample_ids
  )
  gtex_annotation_with_expression_df <- gtex_annotation_matching_result$matched_annotation_df
  gtex_annotation_without_expression_df <- gtex_annotation_matching_result$unmatched_annotation_df
  gtex_annotation_matching_summary_df <- gtex_annotation_matching_result$matching_summary_df
  gtex_annotation_flagged_df <- flag_bulk_gtex_tissues(gtex_annotation_with_expression_df)
  gtex_tissue_summary_df <- summarize_gtex_tissues(gtex_annotation_flagged_df)
  gtex_bulk_tissue_candidate_df <- build_gtex_bulk_tissue_candidate_table(gtex_tissue_summary_df)
  gtex_primary_organ_review_df <- read_or_build_gtex_primary_organ_review_table(mapping_output_xlsx)
  gtex_rule_source <- attr(gtex_primary_organ_review_df, "rule_source")
  validate_gtex_review_coverage(pattern_group_assignment_df, gtex_primary_organ_review_df)
  gtex_final_mapping_rule_df <- extract_final_gtex_mapping_rule_table(gtex_primary_organ_review_df)

  log_message(paste0("当前模式群纳入primary_organ数: ", nrow(pattern_group_assignment_df)))
  log_message(paste0("GTEx注释样本总数: ", nrow(gtex_annotation_raw_df)))
  log_message(paste0("GTEx表达表头样本数: ", length(gtex_expression_header_sample_ids)))
  log_message(paste0("GTEx注释-表达交集样本数: ", nrow(gtex_annotation_with_expression_df)))
  log_message(paste0("GTEx注释中无表达样本数: ", nrow(gtex_annotation_without_expression_df)))
  log_message(paste0("GTEx审阅规则来源: ", gtex_rule_source))
  log_message(
    paste0(
      "GTEx bulk tissue候选样本数: ",
      sum(gtex_annotation_flagged_df$is_bulk_tissue_candidate)
    )
  )
  log_message(paste0("GTEx bulk tissue候选组织数: ", nrow(gtex_bulk_tissue_candidate_df)))


  # 7. 构建项目映射表 ####

  ## 7.1 整理TCGA项目映射 ####
  tcga_project_summary_df <- tcga_project_rule_df %>%
    dplyr::group_by(primary_organ) %>%
    dplyr::summarise(
      tcga_project = collapse_unique_values(tcga_project),
      tcga_code = collapse_unique_values(tcga_code),
      n_tcga_projects = dplyr::n(),
      .groups = "drop"
    )

  ## 7.2 验证GTEx映射规则 ####
  validate_gtex_mapping_rules(
    annotation_df = gtex_annotation_flagged_df,
    mapping_rule_df = gtex_final_mapping_rule_df
  )

  ## 7.3 整理GTEx逐器官审阅结果 ####
  gtex_review_summary_df <- gtex_primary_organ_review_df %>%
    dplyr::group_by(primary_organ) %>%
    dplyr::summarise(
      review_status = collapse_unique_values(decision_status),
      gtex_match_column = collapse_unique_values(match_column[decision_status %in% c("include_direct", "include_proxy")]),
      gtex_match_value = collapse_unique_values(match_value[decision_status %in% c("include_direct", "include_proxy")]),
      selected_gtex_general_tissue = collapse_unique_values(selected_gtex_general_tissue[decision_status %in% c("include_direct", "include_proxy")]),
      any_proxy_tissue = any(is_proxy_tissue[decision_status %in% c("include_direct", "include_proxy")]),
      has_gtex_mapping = any(decision_status %in% c("include_direct", "include_proxy")),
      review_note = collapse_unique_values(review_note),
      n_gtex_rules = sum(decision_status %in% c("include_direct", "include_proxy")),
      .groups = "drop"
    )

  ## 7.4 生成总映射总表 ####
  primary_organ_mapping_df <- pattern_group_assignment_df %>%
    dplyr::left_join(tcga_project_summary_df, by = "primary_organ") %>%
    dplyr::left_join(gtex_review_summary_df, by = "primary_organ") %>%
    dplyr::mutate(
      has_tcga_project = !is.na(tcga_project)
    ) %>%
    dplyr::arrange(factor(pattern_group, levels = pattern_group_order), primary_organ)


  # 8. 根据映射整理GTEx样本注释 ####

  ## 8.1 筛选最终纳入研究的GTEx样本 ####
  gtex_final_sample_annotation_df <- apply_gtex_mapping_rules(
    annotation_df = gtex_annotation_flagged_df,
    mapping_rule_df = gtex_final_mapping_rule_df,
    pattern_group_df = pattern_group_assignment_df
  ) %>%
    dplyr::left_join(tcga_project_summary_df, by = "primary_organ") %>%
    dplyr::select(
      sample_id,
      pattern_group,
      primary_organ,
      tcga_project,
      tcga_code,
      selected_gtex_general_tissue,
      gtex_general_tissue,
      gtex_specific_tissue,
      gtex_match_column,
      gtex_match_value,
      is_proxy_tissue,
      mapping_note
    )

  ## 8.2 汇总纳入样本数 ####
  gtex_included_sample_summary_df <- gtex_final_sample_annotation_df %>%
    dplyr::group_by(pattern_group, primary_organ, tcga_project, tcga_code, gtex_general_tissue, gtex_specific_tissue) %>%
    dplyr::summarise(
      n_samples = dplyr::n(),
      any_proxy_tissue = any(is_proxy_tissue),
      .groups = "drop"
    ) %>%
    dplyr::arrange(factor(pattern_group, levels = pattern_group_order), primary_organ, gtex_general_tissue, gtex_specific_tissue)

  log_message(paste0("建立GTEx映射的primary_organ数: ", sum(primary_organ_mapping_df$has_gtex_mapping)))
  log_message(paste0("未建立GTEx映射的primary_organ数: ", sum(!primary_organ_mapping_df$has_gtex_mapping)))
  log_message(paste0("最终纳入GTEx样本数: ", nrow(gtex_final_sample_annotation_df)))


  # 9. 根据样本注释筛选GTEx表达矩阵 ####

  ## 9.1 优先复用已保存的GTEx表达结果 ####
  expression_result_files <- c(
    gtex_expression_matrix_output_rdata,
    gtex_gene_annotation_output_rdata
  )

  expression_dependency_files <- c(
    gtex_expression_file,
    script_file
  )

  reuse_expression_results <- are_result_files_up_to_date(
    result_files = expression_result_files,
    dependency_files = expression_dependency_files
  )

  if (reuse_expression_results) {
    gtex_cached_expression_result <- try_reuse_cached_gtex_expression_result(
      current_sample_ids = gtex_final_sample_annotation_df$sample_id,
      sample_annotation_rdata_file = gtex_sample_annotation_output_rdata,
      expression_matrix_rdata_file = gtex_expression_matrix_output_rdata,
      gene_annotation_rdata_file = gtex_gene_annotation_output_rdata
    )
  } else {
    gtex_cached_expression_result <- list(
      can_reuse = FALSE,
      reuse_message = "表达结果缓存不存在或已过期",
      expression_matrix = NULL,
      gene_annotation_df = NULL
    )
  }

  if (isTRUE(gtex_cached_expression_result$can_reuse)) {
    gtex_expression_normal_matrix <- gtex_cached_expression_result$expression_matrix
    gtex_gene_annotation_df <- gtex_cached_expression_result$gene_annotation_df

    log_message(
      paste0(
        "复用已保存的GTEx表达结果: ",
        gtex_cached_expression_result$reuse_message
      )
    )
  } else {
    log_message(
      paste0(
        "未复用GTEx表达结果，将重新读取表达矩阵。原因: ",
        gtex_cached_expression_result$reuse_message
      )
    )

    ## 9.2 读取目标样本表达矩阵 ####
    gtex_expression_result <- read_gtex_expression_subset(
      expression_file = gtex_expression_file,
      sample_ids = gtex_final_sample_annotation_df$sample_id
    )

    gtex_gene_annotation_df <- gtex_expression_result$gene_annotation_df
    gtex_expression_normal_matrix <- gtex_expression_result$expression_matrix
  }

  log_message(
    paste0(
      "GTEx表达矩阵维度: ",
      nrow(gtex_expression_normal_matrix),
      " × ",
      ncol(gtex_expression_normal_matrix)
    )
  )


  # 10. 结果保存 ####

  ## 10.1 保存表格 ####
  openxlsx::write.xlsx(
    x = list(
      "原发部位映射总表" = primary_organ_mapping_df,
      "TCGA项目规则表" = tcga_project_rule_df,
      "GTEx候选组织表" = gtex_bulk_tissue_candidate_df,
      "GTEx器官审阅表" = gtex_primary_organ_review_df,
      "GTEx最终规则表" = gtex_final_mapping_rule_df
    ),
    file = mapping_output_xlsx,
    overwrite = TRUE
  )

  openxlsx::write.xlsx(
    x = list(
      "注释表达匹配汇总" = gtex_annotation_matching_summary_df,
      "GTEx组织汇总表" = gtex_tissue_summary_df,
      "注释无表达样本表" = gtex_annotation_without_expression_df,
      "最终纳入样本汇总" = gtex_included_sample_summary_df
    ),
    file = gtex_tissue_summary_output_xlsx,
    overwrite = TRUE
  )

  openxlsx::write.xlsx(
    x = list(
      "GTEx样本注释表" = gtex_final_sample_annotation_df
    ),
    file = gtex_sample_annotation_output_xlsx,
    overwrite = TRUE
  )

  openxlsx::write.xlsx(
    x = list(
      "GTEx基因注释表" = gtex_gene_annotation_df
    ),
    file = gtex_gene_annotation_output_xlsx,
    overwrite = TRUE
  )

  ## 10.2 保存Rdata对象 ####
  save(gtex_final_sample_annotation_df, file = gtex_sample_annotation_output_rdata)
  save(gtex_expression_normal_matrix, file = gtex_expression_matrix_output_rdata)
  save(gtex_gene_annotation_df, file = gtex_gene_annotation_output_rdata)

  log_message(paste0("已导出映射表: ", mapping_output_xlsx))
  log_message(paste0("已导出GTEx器官汇总表: ", gtex_tissue_summary_output_xlsx))
  log_message(paste0("已导出GTEx样本注释表: ", gtex_sample_annotation_output_xlsx))
  log_message(paste0("已导出GTEx基因注释表: ", gtex_gene_annotation_output_xlsx))
  log_message(paste0("已保存GTEx样本注释Rdata: ", gtex_sample_annotation_output_rdata))
  log_message(paste0("已保存GTEx表达矩阵Rdata: ", gtex_expression_matrix_output_rdata))
  log_message(paste0("已保存GTEx基因注释Rdata: ", gtex_gene_annotation_output_rdata))
}

log_message("GTEx数据准备完成。")

# 11. TCGA参数与路径 ####

tcga_by_project_path <- file.path(
  data_path,
  "TCGA_GTEx",
  "by_project"
)

tcga_tumor_sample_annotation_output_xlsx <- file.path(
  table_output_path,
  "tcga_tumor_sample_annotation.xlsx"
)

tcga_normal_sample_annotation_output_xlsx <- file.path(
  table_output_path,
  "tcga_normal_sample_annotation.xlsx"
)

tcga_gene_annotation_output_xlsx <- file.path(
  table_output_path,
  "tcga_gene_annotation.xlsx"
)

tcga_tumor_sample_annotation_output_rdata <- file.path(
  rdata_output_path,
  "tcga_tumor_sample_annotation.Rdata"
)

tcga_normal_sample_annotation_output_rdata <- file.path(
  rdata_output_path,
  "tcga_normal_sample_annotation.Rdata"
)

tcga_tumor_expression_matrix_output_rdata <- file.path(
  rdata_output_path,
  "tcga_tumor_expression_matrix.Rdata"
)

tcga_normal_expression_matrix_output_rdata <- file.path(
  rdata_output_path,
  "tcga_normal_expression_matrix.Rdata"
)

tcga_gene_annotation_output_rdata <- file.path(
  rdata_output_path,
  "tcga_gene_annotation.Rdata"
)

tcga_processing_summary_output_xlsx <- file.path(
  table_output_path,
  "tcga_by_project处理汇总表.xlsx"
)


# 12. TCGA数据检查与缓存判断 ####

log_message("========== TCGA数据整理开始 ==========")
log_message(paste0("TCGA by_project数据路径: ", tcga_by_project_path))

tcga_result_files <- c(
  tcga_tumor_sample_annotation_output_rdata,
  tcga_normal_sample_annotation_output_rdata,
  tcga_tumor_expression_matrix_output_rdata,
  tcga_normal_expression_matrix_output_rdata,
  tcga_gene_annotation_output_rdata,
  tcga_tumor_sample_annotation_output_xlsx,
  tcga_normal_sample_annotation_output_xlsx,
  tcga_gene_annotation_output_xlsx,
  tcga_processing_summary_output_xlsx
)

tcga_dependency_files <- c(
  pattern_group_assignment_file,
  script_file
)

reuse_tcga_results <- are_result_files_up_to_date(
  result_files = tcga_result_files,
  dependency_files = tcga_dependency_files
)

if (reuse_tcga_results) {
  log_message("检测到已保存且为最新的TCGA结果，跳过TCGA数据重算。")
  log_message(paste0("复用TCGA肿瘤样本注释: ", tcga_tumor_sample_annotation_output_rdata))
  log_message(paste0("复用TCGA正常样本注释: ", tcga_normal_sample_annotation_output_rdata))
  log_message(paste0("复用TCGA肿瘤表达矩阵: ", tcga_tumor_expression_matrix_output_rdata))
  log_message(paste0("复用TCGA正常表达矩阵: ", tcga_normal_expression_matrix_output_rdata))
} else {
  log_message("未检测到可直接复用的TCGA结果，将重新整理TCGA数据。")


  # 13. 读取并整合TCGA各项目数据 ####

  ## 13.1 获取映射规则 ####
  if (!exists("pattern_group_assignment_df")) {
    pattern_group_assignment_df <- read_pattern_group_assignment(pattern_group_assignment_file)
  }
  if (!exists("tcga_project_rule_df")) {
    tcga_project_rule_df <- build_tcga_project_rule_table()
  }

  ## 13.2 检查各项目文件可用性 ####
  unique_tcga_projects <- unique(tcga_project_rule_df$tcga_project)

  tcga_project_file_check_list <- lapply(unique_tcga_projects, function(proj) {
    tpm_file <- file.path(tcga_by_project_path, paste0(proj, "_tpm_matrix.rds"))
    annotation_file <- file.path(tcga_by_project_path, paste0(proj, "_sample_annotation.csv.gz"))
    gene_annotation_file <- file.path(tcga_by_project_path, paste0(proj, "_gene_annotation.csv.gz"))

    tibble::tibble(
      tcga_project = proj,
      tpm_file = tpm_file,
      annotation_file = annotation_file,
      gene_annotation_file = gene_annotation_file,
      tpm_file_exists = file.exists(tpm_file),
      annotation_file_exists = file.exists(annotation_file),
      gene_annotation_file_exists = file.exists(gene_annotation_file),
      all_files_exist = file.exists(tpm_file) & file.exists(annotation_file) & file.exists(gene_annotation_file)
    )
  })
  tcga_project_file_check_df <- dplyr::bind_rows(tcga_project_file_check_list)

  available_projects <- tcga_project_file_check_df %>%
    dplyr::filter(all_files_exist) %>%
    dplyr::pull(tcga_project)

  missing_projects <- tcga_project_file_check_df %>%
    dplyr::filter(!all_files_exist) %>%
    dplyr::pull(tcga_project)

  log_message(paste0("TCGA映射规则涉及项目数: ", length(unique_tcga_projects)))
  log_message(paste0("文件齐全可处理的项目数: ", length(available_projects)))
  if (length(missing_projects) > 0) {
    log_message(paste0(
      "WARNING: 以下项目文件缺失，将跳过: ",
      paste(missing_projects, collapse = ", ")
    ))
  }

  ## 13.3 逐项目读取并拆分肿瘤/正常样本 ####
  tcga_tumor_sample_list <- list()
  tcga_normal_sample_list <- list()
  tcga_tumor_matrix_list <- list()
  tcga_normal_matrix_list <- list()
  tcga_project_summary_list <- list()
  tcga_reference_gene_annotation_df <- NULL

  for (proj in available_projects) {
    log_message(paste0("读取TCGA项目: ", proj))

    file_info <- tcga_project_file_check_df %>%
      dplyr::filter(tcga_project == proj)

    ## 读取样本注释 ####
    sample_annotation_df <- data.table::fread(
      input = file_info$annotation_file,
      data.table = FALSE,
      showProgress = FALSE
    ) %>%
      dplyr::mutate(
        project_id = as.character(project_id),
        sample_barcode = as.character(sample_barcode),
        sample_id = as.character(sample_id),
        patient_id = as.character(patient_id),
        sample_type = as.character(sample_type),
        tissue_type = as.character(tissue_type)
      )

    ## 读取TPM矩阵 ####
    tpm_matrix <- readRDS(file_info$tpm_file)

    ## 读取基因注释（首次读取时保存作为参考） ####
    if (is.null(tcga_reference_gene_annotation_df)) {
      tcga_reference_gene_annotation_df <- data.table::fread(
        input = file_info$gene_annotation_file,
        data.table = FALSE,
        showProgress = FALSE
      ) %>%
        dplyr::mutate(
          gene_name = as.character(gene_name),
          gene_id = as.character(gene_id)
        )
      log_message(paste0("参考基因注释来源: ", proj, " (", nrow(tcga_reference_gene_annotation_df), " genes)"))
    }

    ## 验证样本注释与矩阵的一致性 ####
    annotation_barcodes <- sample_annotation_df$sample_barcode
    matrix_barcodes <- colnames(tpm_matrix)
    common_barcodes <- intersect(annotation_barcodes, matrix_barcodes)

    if (length(common_barcodes) != length(annotation_barcodes)) {
      log_message(paste0(
        "WARNING: ", proj, " 样本注释与矩阵不完全匹配。",
        " 注释: ", length(annotation_barcodes),
        " 矩阵: ", length(matrix_barcodes),
        " 交集: ", length(common_barcodes)
      ))
    }

    ## 拆分肿瘤与正常样本 ####
    tumor_annotation_df <- sample_annotation_df %>%
      dplyr::filter(
        sample_type == "Primary Tumor",
        sample_barcode %in% matrix_barcodes
      )

    normal_annotation_df <- sample_annotation_df %>%
      dplyr::filter(
        sample_type == "Solid Tissue Normal",
        sample_barcode %in% matrix_barcodes
      )

    n_tumor <- nrow(tumor_annotation_df)
    n_normal <- nrow(normal_annotation_df)

    log_message(paste0(
      "  ", proj, ": 肿瘤样本 ", n_tumor, ", 正常样本 ", n_normal
    ))

    ## 提取对应子矩阵 ####
    if (n_tumor > 0) {
      tumor_barcodes <- tumor_annotation_df$sample_barcode
      tcga_tumor_matrix_list[[proj]] <- tpm_matrix[, tumor_barcodes, drop = FALSE]
      tcga_tumor_sample_list[[proj]] <- tumor_annotation_df
    }

    if (n_normal > 0) {
      normal_barcodes <- normal_annotation_df$sample_barcode
      tcga_normal_matrix_list[[proj]] <- tpm_matrix[, normal_barcodes, drop = FALSE]
      tcga_normal_sample_list[[proj]] <- normal_annotation_df
    }

    ## 记录汇总信息 ####
    tcga_project_summary_list[[proj]] <- tibble::tibble(
      tcga_project = proj,
      n_annotation_samples = length(annotation_barcodes),
      n_matrix_samples = length(matrix_barcodes),
      n_matched = length(common_barcodes),
      n_tumor = n_tumor,
      n_normal = n_normal,
      n_genes = nrow(tpm_matrix)
    )

    ## 及时释放大矩阵 ####
    rm(tpm_matrix, sample_annotation_df)
    gc()
  }

  ## 13.4 验证各项目基因维度一致性 ####
  gene_counts <- sapply(c(tcga_tumor_matrix_list, tcga_normal_matrix_list), nrow)
  if (length(unique(gene_counts)) > 1) {
    log_message("WARNING: 各项目基因数不一致，请检查数据源。")
    log_message(paste0("  基因数分布: ", paste(unique(gene_counts), collapse = ", ")))
  }

  gene_names_ref <- rownames(tcga_tumor_matrix_list[[1]])
  gene_name_consistent <- all(sapply(
    c(tcga_tumor_matrix_list, tcga_normal_matrix_list),
    function(m) identical(rownames(m), gene_names_ref)
  ))

  if (!gene_name_consistent) {
    log_message("WARNING: 各项目基因名称顺序不完全一致，将按参考基因名称重排。")
    tcga_tumor_matrix_list <- lapply(tcga_tumor_matrix_list, function(m) m[gene_names_ref, , drop = FALSE])
    tcga_normal_matrix_list <- lapply(tcga_normal_matrix_list, function(m) m[gene_names_ref, , drop = FALSE])
  } else {
    log_message("各项目基因名称顺序一致，无需重排。")
  }

  ## 13.5 合并为统一表达矩阵 ####

  ### 合并肿瘤表达矩阵 ####
  if (length(tcga_tumor_matrix_list) > 0) {
    tcga_tumor_expression_matrix <- do.call(cbind, tcga_tumor_matrix_list)
    storage.mode(tcga_tumor_expression_matrix) <- "numeric"
    log_message(paste0(
      "TCGA肿瘤表达矩阵维度: ",
      nrow(tcga_tumor_expression_matrix), " × ", ncol(tcga_tumor_expression_matrix)
    ))
  } else {
    tcga_tumor_expression_matrix <- matrix(nrow = 0, ncol = 0)
    log_message("WARNING: 未找到任何TCGA肿瘤样本。")
  }

  ### 合并正常组织表达矩阵 ####
  if (length(tcga_normal_matrix_list) > 0) {
    tcga_normal_expression_matrix <- do.call(cbind, tcga_normal_matrix_list)
    storage.mode(tcga_normal_expression_matrix) <- "numeric"
    log_message(paste0(
      "TCGA正常组织表达矩阵维度: ",
      nrow(tcga_normal_expression_matrix), " × ", ncol(tcga_normal_expression_matrix)
    ))
  } else {
    tcga_normal_expression_matrix <- matrix(nrow = 0, ncol = 0)
    log_message("WARNING: 未找到任何TCGA正常组织样本。")
  }

  ## 释放子矩阵列表 ####
  rm(tcga_tumor_matrix_list, tcga_normal_matrix_list)
  gc()

  ## 13.6 整合样本注释并添加primary_organ和pattern_group ####

  ### 合并肿瘤样本注释 ####
  if (length(tcga_tumor_sample_list) > 0) {
    tcga_tumor_sample_annotation_df <- dplyr::bind_rows(tcga_tumor_sample_list) %>%
      dplyr::left_join(
        tcga_project_rule_df %>% dplyr::select(tcga_project, primary_organ),
        by = c("project_id" = "tcga_project")
      ) %>%
      dplyr::left_join(
        pattern_group_assignment_df %>% dplyr::select(primary_organ, pattern_group),
        by = "primary_organ"
      ) %>%
      dplyr::select(
        sample_barcode,
        sample_id,
        patient_id,
        project_id,
        sample_type,
        tissue_type,
        primary_organ,
        pattern_group,
        sample_rank_within_patient
      ) %>%
      dplyr::arrange(
        factor(pattern_group, levels = pattern_group_order),
        primary_organ,
        project_id,
        sample_barcode
      ) %>%
      dplyr::mutate(pattern_group = as.character(pattern_group))
  } else {
    tcga_tumor_sample_annotation_df <- tibble::tibble(
      sample_barcode = character(),
      sample_id = character(),
      patient_id = character(),
      project_id = character(),
      sample_type = character(),
      tissue_type = character(),
      primary_organ = character(),
      pattern_group = character(),
      sample_rank_within_patient = integer()
    )
  }

  ### 合并正常组织样本注释 ####
  if (length(tcga_normal_sample_list) > 0) {
    tcga_normal_sample_annotation_df <- dplyr::bind_rows(tcga_normal_sample_list) %>%
      dplyr::left_join(
        tcga_project_rule_df %>% dplyr::select(tcga_project, primary_organ),
        by = c("project_id" = "tcga_project")
      ) %>%
      dplyr::left_join(
        pattern_group_assignment_df %>% dplyr::select(primary_organ, pattern_group),
        by = "primary_organ"
      ) %>%
      dplyr::select(
        sample_barcode,
        sample_id,
        patient_id,
        project_id,
        sample_type,
        tissue_type,
        primary_organ,
        pattern_group,
        sample_rank_within_patient
      ) %>%
      dplyr::arrange(
        factor(pattern_group, levels = pattern_group_order),
        primary_organ,
        project_id,
        sample_barcode
      ) %>%
      dplyr::mutate(pattern_group = as.character(pattern_group))
  } else {
    tcga_normal_sample_annotation_df <- tibble::tibble(
      sample_barcode = character(),
      sample_id = character(),
      patient_id = character(),
      project_id = character(),
      sample_type = character(),
      tissue_type = character(),
      primary_organ = character(),
      pattern_group = character(),
      sample_rank_within_patient = integer()
    )
  }

  ## 释放子列表 ####
  rm(tcga_tumor_sample_list, tcga_normal_sample_list)
  gc()

  ## 13.7 整理基因注释 ####
  tcga_gene_annotation_df <- tcga_reference_gene_annotation_df

  ## 13.8 汇总处理结果 ####
  tcga_project_processing_summary_df <- dplyr::bind_rows(tcga_project_summary_list)

  tcga_overall_summary_df <- tibble::tibble(
    metric = c(
      "映射规则涉及项目数",
      "文件齐全可处理项目数",
      "文件缺失跳过项目数",
      "TCGA肿瘤样本总数",
      "TCGA正常组织样本总数",
      "TCGA基因数",
      "TCGA肿瘤覆盖primary_organ数",
      "TCGA正常组织覆盖primary_organ数"
    ),
    value = c(
      length(unique_tcga_projects),
      length(available_projects),
      length(missing_projects),
      nrow(tcga_tumor_sample_annotation_df),
      nrow(tcga_normal_sample_annotation_df),
      nrow(tcga_gene_annotation_df),
      length(unique(tcga_tumor_sample_annotation_df$primary_organ)),
      length(unique(tcga_normal_sample_annotation_df$primary_organ))
    )
  )

  log_message(paste0("TCGA肿瘤样本总数: ", nrow(tcga_tumor_sample_annotation_df)))
  log_message(paste0("TCGA正常组织样本总数: ", nrow(tcga_normal_sample_annotation_df)))
  log_message(paste0("TCGA基因数: ", nrow(tcga_gene_annotation_df)))
  log_message(paste0(
    "TCGA肿瘤覆盖primary_organ: ",
    paste(sort(unique(tcga_tumor_sample_annotation_df$primary_organ)), collapse = ", ")
  ))
  log_message(paste0(
    "TCGA正常组织覆盖primary_organ: ",
    paste(sort(unique(tcga_normal_sample_annotation_df$primary_organ)), collapse = ", ")
  ))


  # 14. 保存TCGA结果 ####

  ## 14.1 保存表格 ####
  openxlsx::write.xlsx(
    x = list(
      "TCGA肿瘤样本注释" = tcga_tumor_sample_annotation_df
    ),
    file = tcga_tumor_sample_annotation_output_xlsx,
    overwrite = TRUE
  )

  openxlsx::write.xlsx(
    x = list(
      "TCGA正常样本注释" = tcga_normal_sample_annotation_df
    ),
    file = tcga_normal_sample_annotation_output_xlsx,
    overwrite = TRUE
  )

  openxlsx::write.xlsx(
    x = list(
      "TCGA基因注释" = tcga_gene_annotation_df
    ),
    file = tcga_gene_annotation_output_xlsx,
    overwrite = TRUE
  )

  openxlsx::write.xlsx(
    x = list(
      "逐项目处理汇总" = tcga_project_processing_summary_df,
      "文件可用性检查" = tcga_project_file_check_df %>%
        dplyr::select(tcga_project, tpm_file_exists, annotation_file_exists, gene_annotation_file_exists, all_files_exist),
      "总体汇总" = tcga_overall_summary_df
    ),
    file = tcga_processing_summary_output_xlsx,
    overwrite = TRUE
  )

  ## 14.2 保存Rdata对象 ####
  save(tcga_tumor_sample_annotation_df, file = tcga_tumor_sample_annotation_output_rdata)
  save(tcga_normal_sample_annotation_df, file = tcga_normal_sample_annotation_output_rdata)
  save(tcga_tumor_expression_matrix, file = tcga_tumor_expression_matrix_output_rdata)
  save(tcga_normal_expression_matrix, file = tcga_normal_expression_matrix_output_rdata)
  save(tcga_gene_annotation_df, file = tcga_gene_annotation_output_rdata)

  log_message(paste0("已导出TCGA肿瘤样本注释: ", tcga_tumor_sample_annotation_output_xlsx))
  log_message(paste0("已导出TCGA正常样本注释: ", tcga_normal_sample_annotation_output_xlsx))
  log_message(paste0("已导出TCGA基因注释: ", tcga_gene_annotation_output_xlsx))
  log_message(paste0("已导出TCGA处理汇总: ", tcga_processing_summary_output_xlsx))
  log_message(paste0("已保存TCGA肿瘤样本注释Rdata: ", tcga_tumor_sample_annotation_output_rdata))
  log_message(paste0("已保存TCGA正常样本注释Rdata: ", tcga_normal_sample_annotation_output_rdata))
  log_message(paste0("已保存TCGA肿瘤表达矩阵Rdata: ", tcga_tumor_expression_matrix_output_rdata))
  log_message(paste0("已保存TCGA正常表达矩阵Rdata: ", tcga_normal_expression_matrix_output_rdata))
  log_message(paste0("已保存TCGA基因注释Rdata: ", tcga_gene_annotation_output_rdata))
}

log_message("TCGA数据整理完成。")


# 15. 跨平台基因统一 ####
# 说明：GTEx矩阵行名为Ensembl ID(带版本号)，TCGA矩阵行名为gene symbol。
#       本节以Ensembl ID(去版本号)为桥梁取交集，统一三套矩阵的基因维度，
#       并将GTEx矩阵行名转换为gene symbol。

log_message("========== 跨平台基因统一开始 ==========")

## 15.1 输出路径 ####
unified_output_rdata_path <- file.path(rdata_output_path, "unified")
dir.create(unified_output_rdata_path, recursive = TRUE, showWarnings = FALSE)

unified_gtex_expression_rdata <- file.path(unified_output_rdata_path, "unified_gtex_normal_expression_matrix.Rdata")
unified_tcga_tumor_expression_rdata <- file.path(unified_output_rdata_path, "unified_tcga_tumor_expression_matrix.Rdata")
unified_tcga_normal_expression_rdata <- file.path(unified_output_rdata_path, "unified_tcga_normal_expression_matrix.Rdata")
unified_gene_annotation_rdata <- file.path(unified_output_rdata_path, "unified_gene_annotation.Rdata")
unified_gene_annotation_xlsx <- file.path(table_output_path, "unified_gene_annotation.xlsx")
unified_summary_xlsx <- file.path(table_output_path, "跨平台基因统一汇总表.xlsx")

unified_result_files <- c(
  unified_gtex_expression_rdata,
  unified_tcga_tumor_expression_rdata,
  unified_tcga_normal_expression_rdata,
  unified_gene_annotation_rdata,
  unified_gene_annotation_xlsx,
  unified_summary_xlsx
)

unified_dependency_files <- c(
  gtex_expression_matrix_output_rdata,
  gtex_gene_annotation_output_rdata,
  tcga_tumor_expression_matrix_output_rdata,
  tcga_normal_expression_matrix_output_rdata,
  tcga_gene_annotation_output_rdata,
  script_file
)

reuse_unified_results <- are_result_files_up_to_date(
  result_files = unified_result_files,
  dependency_files = unified_dependency_files
)

if (reuse_unified_results) {
  log_message("检测到已保存且为最新的跨平台统一结果，跳过重算。")
  log_message(paste0("复用统一GTEx矩阵: ", unified_gtex_expression_rdata))
  log_message(paste0("复用统一TCGA肿瘤矩阵: ", unified_tcga_tumor_expression_rdata))
  log_message(paste0("复用统一TCGA正常矩阵: ", unified_tcga_normal_expression_rdata))
} else {
  log_message("未检测到可直接复用的统一结果，将重新进行跨平台基因统一。")

  ## 15.2 加载所需数据 ####
  if (!exists("gtex_gene_annotation_df")) {
    gtex_gene_annotation_df <- load_single_object_from_rdata(
      rdata_file = gtex_gene_annotation_output_rdata,
      object_name = "gtex_gene_annotation_df"
    )
  }
  if (!exists("tcga_gene_annotation_df")) {
    tcga_gene_annotation_df <- load_single_object_from_rdata(
      rdata_file = tcga_gene_annotation_output_rdata,
      object_name = "tcga_gene_annotation_df"
    )
  }

  ## 15.3 构建Ensembl ID去版本号映射 ####
  gtex_gene_map <- data.frame(
    gtex_gene_id = gtex_gene_annotation_df$gene_id,
    gtex_gene_symbol = gtex_gene_annotation_df$gene_symbol,
    ensembl_base = gsub("[.].*", "", gtex_gene_annotation_df$gene_id),
    stringsAsFactors = FALSE
  )

  tcga_gene_map <- data.frame(
    tcga_gene_name = tcga_gene_annotation_df$gene_name,
    tcga_gene_id = tcga_gene_annotation_df$gene_id,
    ensembl_base = gsub("[.].*", "", tcga_gene_annotation_df$gene_id),
    stringsAsFactors = FALSE
  )

  ## 15.4 取交集（通过Ensembl base ID内连接） ####
  unified_gene_map <- merge(
    gtex_gene_map,
    tcga_gene_map,
    by = "ensembl_base",
    all = FALSE
  )

  ## 检查是否有重复的ensembl_base（理论上不应有） ####
  dup_ensembl <- unified_gene_map$ensembl_base[duplicated(unified_gene_map$ensembl_base)]
  if (length(dup_ensembl) > 0) {
    log_message(paste0(
      "WARNING: 合并后存在重复ensembl_base ", length(unique(dup_ensembl)),
      " 个，将保留首次出现项。"
    ))
    unified_gene_map <- unified_gene_map[!duplicated(unified_gene_map$ensembl_base), ]
  }

  ## 检查统一后gene symbol是否有重复 ####
  dup_symbols <- unified_gene_map$tcga_gene_name[duplicated(unified_gene_map$tcga_gene_name)]
  if (length(dup_symbols) > 0) {
    log_message(paste0(
      "WARNING: 统一后存在重复gene symbol ", length(unique(dup_symbols)),
      " 个，将保留首次出现项。"
    ))
    unified_gene_map <- unified_gene_map[!duplicated(unified_gene_map$tcga_gene_name), ]
  }

  ## 统一后的基因数 ####
  n_unified_genes <- nrow(unified_gene_map)
  log_message(paste0("GTEx基因数: ", nrow(gtex_gene_map)))
  log_message(paste0("TCGA基因数: ", nrow(tcga_gene_map)))
  log_message(paste0("统一后基因数: ", n_unified_genes))

  ## 15.5 制作统一基因注释表 ####
  unified_gene_annotation_df <- data.frame(
    gene_symbol = unified_gene_map$tcga_gene_name,
    ensembl_id = unified_gene_map$ensembl_base,
    gtex_gene_id = unified_gene_map$gtex_gene_id,
    tcga_gene_id = unified_gene_map$tcga_gene_id,
    stringsAsFactors = FALSE
  )

  ## 15.6 提取统一基因子集的GTEx表达矩阵 ####
  log_message("提取GTEx统一基因子集...")
  if (!exists("gtex_expression_normal_matrix")) {
    gtex_expression_normal_matrix <- load_single_object_from_rdata(
      rdata_file = gtex_expression_matrix_output_rdata,
      object_name = "gtex_expression_normal_matrix"
    )
  }

  unified_gtex_normal_expression_matrix <- gtex_expression_normal_matrix[
    unified_gene_map$gtex_gene_id, , drop = FALSE
  ]
  rownames(unified_gtex_normal_expression_matrix) <- unified_gene_map$tcga_gene_name
  storage.mode(unified_gtex_normal_expression_matrix) <- "numeric"

  log_message(paste0(
    "统一GTEx矩阵维度: ",
    nrow(unified_gtex_normal_expression_matrix), " × ",
    ncol(unified_gtex_normal_expression_matrix)
  ))

  ## 释放原始GTEx矩阵 ####
  rm(gtex_expression_normal_matrix)
  gc()

  ## 15.7 提取统一基因子集的TCGA肿瘤表达矩阵 ####
  log_message("提取TCGA肿瘤统一基因子集...")
  if (!exists("tcga_tumor_expression_matrix")) {
    tcga_tumor_expression_matrix <- load_single_object_from_rdata(
      rdata_file = tcga_tumor_expression_matrix_output_rdata,
      object_name = "tcga_tumor_expression_matrix"
    )
  }

  unified_tcga_tumor_expression_matrix <- tcga_tumor_expression_matrix[
    unified_gene_map$tcga_gene_name, , drop = FALSE
  ]
  storage.mode(unified_tcga_tumor_expression_matrix) <- "numeric"

  log_message(paste0(
    "统一TCGA肿瘤矩阵维度: ",
    nrow(unified_tcga_tumor_expression_matrix), " × ",
    ncol(unified_tcga_tumor_expression_matrix)
  ))

  rm(tcga_tumor_expression_matrix)
  gc()

  ## 15.8 提取统一基因子集的TCGA正常组织表达矩阵 ####
  log_message("提取TCGA正常组织统一基因子集...")
  if (!exists("tcga_normal_expression_matrix")) {
    tcga_normal_expression_matrix <- load_single_object_from_rdata(
      rdata_file = tcga_normal_expression_matrix_output_rdata,
      object_name = "tcga_normal_expression_matrix"
    )
  }

  unified_tcga_normal_expression_matrix <- tcga_normal_expression_matrix[
    unified_gene_map$tcga_gene_name, , drop = FALSE
  ]
  storage.mode(unified_tcga_normal_expression_matrix) <- "numeric"

  log_message(paste0(
    "统一TCGA正常矩阵维度: ",
    nrow(unified_tcga_normal_expression_matrix), " × ",
    ncol(unified_tcga_normal_expression_matrix)
  ))

  rm(tcga_normal_expression_matrix)
  gc()

  ## 15.9 汇总统一结果 ####
  unified_summary_df <- tibble::tibble(
    metric = c(
      "GTEx原始基因数",
      "TCGA原始基因数",
      "Ensembl ID(去版本号)交集数",
      "去重后统一基因数",
      "统一GTEx矩阵行数",
      "统一GTEx矩阵列数(样本)",
      "统一TCGA肿瘤矩阵行数",
      "统一TCGA肿瘤矩阵列数(样本)",
      "统一TCGA正常矩阵行数",
      "统一TCGA正常矩阵列数(样本)"
    ),
    value = c(
      nrow(gtex_gene_map),
      nrow(tcga_gene_map),
      length(intersect(gtex_gene_map$ensembl_base, tcga_gene_map$ensembl_base)),
      n_unified_genes,
      nrow(unified_gtex_normal_expression_matrix),
      ncol(unified_gtex_normal_expression_matrix),
      nrow(unified_tcga_tumor_expression_matrix),
      ncol(unified_tcga_tumor_expression_matrix),
      nrow(unified_tcga_normal_expression_matrix),
      ncol(unified_tcga_normal_expression_matrix)
    )
  )

  ## 15.10 保存统一结果 ####
  save(unified_gtex_normal_expression_matrix, file = unified_gtex_expression_rdata)
  save(unified_tcga_tumor_expression_matrix, file = unified_tcga_tumor_expression_rdata)
  save(unified_tcga_normal_expression_matrix, file = unified_tcga_normal_expression_rdata)
  save(unified_gene_annotation_df, file = unified_gene_annotation_rdata)

  openxlsx::write.xlsx(
    x = list(
      "统一基因注释" = unified_gene_annotation_df
    ),
    file = unified_gene_annotation_xlsx,
    overwrite = TRUE
  )

  openxlsx::write.xlsx(
    x = list(
      "统一汇总" = unified_summary_df
    ),
    file = unified_summary_xlsx,
    overwrite = TRUE
  )

  log_message(paste0("已保存统一GTEx矩阵: ", unified_gtex_expression_rdata))
  log_message(paste0("已保存统一TCGA肿瘤矩阵: ", unified_tcga_tumor_expression_rdata))
  log_message(paste0("已保存统一TCGA正常矩阵: ", unified_tcga_normal_expression_rdata))
  log_message(paste0("已保存统一基因注释: ", unified_gene_annotation_rdata))
  log_message(paste0("已导出统一基因注释表: ", unified_gene_annotation_xlsx))
  log_message(paste0("已导出统一汇总表: ", unified_summary_xlsx))
}

log_message("跨平台基因统一完成。")


# 16. 汇总映射表（含样本量） ####
# 说明：从已保存的样本注释文件中统计各primary_organ的三类样本量，更新映射总表 ####

log_message("========== 汇总映射表（含样本量）开始 ==========")

final_mapping_output_xlsx <- file.path(
  table_output_path,
  "模式群-原发部位-TCGA项目-GTEx器官映射表（含样本量）.xlsx"
)

## 16.1 读取映射总表 ####
if (!exists("primary_organ_mapping_df")) {
  primary_organ_mapping_df <- readxl::read_xlsx(
    path = mapping_output_xlsx,
    sheet = "原发部位映射总表"
  )
}

## 16.2 读取各来源样本注释 ####
if (!exists("gtex_final_sample_annotation_df")) {
  gtex_final_sample_annotation_df <- load_single_object_from_rdata(
    rdata_file = gtex_sample_annotation_output_rdata,
    object_name = "gtex_final_sample_annotation_df"
  )
}
if (!exists("tcga_tumor_sample_annotation_df")) {
  tcga_tumor_sample_annotation_df <- load_single_object_from_rdata(
    rdata_file = tcga_tumor_sample_annotation_output_rdata,
    object_name = "tcga_tumor_sample_annotation_df"
  )
}
if (!exists("tcga_normal_sample_annotation_df")) {
  tcga_normal_sample_annotation_df <- load_single_object_from_rdata(
    rdata_file = tcga_normal_sample_annotation_output_rdata,
    object_name = "tcga_normal_sample_annotation_df"
  )
}

## 16.3 统计各primary_organ的样本量 ####
gtex_sample_count_df <- gtex_final_sample_annotation_df %>%
  dplyr::count(primary_organ, name = "n_gtex_normal")

tcga_tumor_count_df <- tcga_tumor_sample_annotation_df %>%
  dplyr::count(primary_organ, name = "n_tcga_tumor")

tcga_normal_count_df <- tcga_normal_sample_annotation_df %>%
  dplyr::count(primary_organ, name = "n_tcga_normal")

## 16.4 合并到映射总表 ####
primary_organ_mapping_with_counts_df <- primary_organ_mapping_df %>%
  dplyr::left_join(gtex_sample_count_df, by = "primary_organ") %>%
  dplyr::left_join(tcga_tumor_count_df, by = "primary_organ") %>%
  dplyr::left_join(tcga_normal_count_df, by = "primary_organ") %>%
  dplyr::mutate(
    n_gtex_normal = dplyr::if_else(is.na(n_gtex_normal), 0L, as.integer(n_gtex_normal)),
    n_tcga_tumor = dplyr::if_else(is.na(n_tcga_tumor), 0L, as.integer(n_tcga_tumor)),
    n_tcga_normal = dplyr::if_else(is.na(n_tcga_normal), 0L, as.integer(n_tcga_normal)),
    n_total = n_gtex_normal + n_tcga_tumor + n_tcga_normal
  )

## 16.5 保存含样本量的映射表 ####
openxlsx::write.xlsx(
  x = list(
    "原发部位映射总表（含样本量）" = primary_organ_mapping_with_counts_df
  ),
  file = final_mapping_output_xlsx,
  overwrite = TRUE
)

log_message(paste0(
  "映射表样本量汇总: ",
  "GTEx正常 ", sum(primary_organ_mapping_with_counts_df$n_gtex_normal),
  " | TCGA肿瘤 ", sum(primary_organ_mapping_with_counts_df$n_tcga_tumor),
  " | TCGA正常 ", sum(primary_organ_mapping_with_counts_df$n_tcga_normal)
))
log_message(paste0("已导出含样本量映射表: ", final_mapping_output_xlsx))
log_message("汇总映射表完成。")

log_message("========== 7.1 TCGA基础数据准备全部完成 ==========")
