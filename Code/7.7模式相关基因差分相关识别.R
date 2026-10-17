# 7.7 模式相关基因差分相关识别 ####
# 用途：基于差分 Spearman 相关 (Δρ = ρ_tumor − ρ_normal) 框架，
#       在 primary_organ 级识别与 6.1 NMF 5 个 irAE 模式相关的基因表达特征。
#       本脚本采用纯 Δρ 框架: hit 由 |Δρ| 阈值 + matched-pair permutation p 决定 (主筛选),
#       高置信子集追加 Fisher r-to-z p (parametric robustness)。
#       正常端使用 7.2 ComBat 校正后的 GTEx + TCGA-normal 联合 log2 矩阵。
#       organ 按正常对照样本量分 Tier 1 / 1b / 2 / 3, Δρ 仅在 Tier1+1b+2 上计算。
#       Δρ 显著性以 matched organ 子集 organ 标签 same-seed permutation 零分布的双尾经验 p 为主检验,
#       Fisher r-to-z 变换作为参数化 robustness check 平行计算 (高置信子集追加 Fisher p < 0.05)。
#       BH FDR 在 N=20 organ × m=15388 gene 下数学上不可达,作为透明性指标输出, 不参与 hit 判定。 ####

rm(list = ls())
gc()


# 1. 加载配置 ####
source("Code/config.R")


# 2. 参数与路径 ####
script_name <- "7.7模式相关基因差分相关识别"
script_results_path <- file.path(results_path, script_name)
table_output_path <- file.path(script_results_path, "tables")
rdata_output_path <- file.path(script_results_path, "rdata")
log_output_path <- file.path(script_results_path, "logs")

dir.create(script_results_path, recursive = TRUE, showWarnings = FALSE)
dir.create(table_output_path, recursive = TRUE, showWarnings = FALSE)
dir.create(rdata_output_path, recursive = TRUE, showWarnings = FALSE)
dir.create(log_output_path, recursive = TRUE, showWarnings = FALSE)

## 输入文件路径 ####
input_71_path <- file.path(results_path, "7.1TCGA基础数据准备")
input_71_rdata_path <- file.path(input_71_path, "rdata")

tcga_tumor_expression_rdata <- file.path(
  input_71_rdata_path,
  "tcga_tumor_expression_matrix.Rdata"
)
tcga_tumor_sample_annotation_rdata <- file.path(
  input_71_rdata_path,
  "tcga_tumor_sample_annotation.Rdata"
)
tcga_normal_sample_annotation_rdata <- file.path(
  input_71_rdata_path,
  "tcga_normal_sample_annotation.Rdata"
)
gtex_normal_sample_annotation_rdata <- file.path(
  input_71_rdata_path,
  "gtex_normal_sample_annotation.Rdata"
)

input_72_path <- file.path(results_path, "7.2正常组织PCA分析")
input_72_rdata_path <- file.path(input_72_path, "rdata")
corrected_normal_log2_rdata <- file.path(
  input_72_rdata_path,
  "corrected_log2_matrix.Rdata"
)

input_61_path <- file.path(results_path, "6.1癌种NMF分群分析_rank5")
nmf_best_fit_rdata <- file.path(input_61_path, "FAERS_NMF_best_fit.rdata")
nmf_h_cluster_rdata <- file.path(input_61_path, "FAERS_NMF_H_cluster.rdata")

## 输出文件路径 ####
output_main_xlsx <- file.path(table_output_path, "模式相关基因差分相关分析结果.xlsx")
organ_median_log2_tumor_rdata <- file.path(rdata_output_path, "organ_median_log2_tumor_matrix.Rdata")
organ_median_log2_normal_rdata <- file.path(rdata_output_path, "organ_median_log2_normal_matrix.Rdata")
organ_tier_rdata <- file.path(rdata_output_path, "organ_tier_classification_df.Rdata")
gene_filter_summary_rdata <- file.path(rdata_output_path, "gene_filter_summary_df.Rdata")
rho_tumor_observed_rdata <- file.path(rdata_output_path, "rho_tumor_matched_observed_df.Rdata")
rho_normal_observed_rdata <- file.path(rdata_output_path, "rho_normal_matched_observed_df.Rdata")
delta_rho_df_rdata <- file.path(rdata_output_path, "delta_rho_df.Rdata")
permutation_null_array_rdata <- file.path(rdata_output_path, "delta_permutation_null_array.Rdata")
final_pattern_gene_list_rdata <- file.path(rdata_output_path, "final_pattern_gene_list.Rdata")

## 主要参数 ####
log2_offset <- 1
keep_first_tumor_per_patient <- TRUE

## organ 纳入与 Tier 分类参数 ####
min_samples_per_organ_tumor <- 10               # primary_organ 进入肿瘤端分析的最低肿瘤样本数
tcga_normal_min_for_tier1 <- 5                  # TCGA-normal ≥ 5 → Tier 1（联合 GTEx）
gtex_min_for_normal_baseline <- 30              # GTEx ≥ 30 → 可作为正常 organ-median 单源参照（Tier 1b/2）

## 基因预过滤参数（仅保留计算可行性过滤；tissue-of-origin 混杂由 Δρ 框架本身消除，无须额外预过滤） ####
min_organ_median_log2 <- 1
min_organ_median_mad <- 0.3

## 主分析参数 ####
n_permutation_iterations <- 1000
analysis_seed <- 20260506

## Δρ 筛选阈值 ####
## 主检验: |Δρ| ≥ 阈值 + matched-pair permutation 双尾经验 p < perm_p_threshold ####
## 高置信子集 (high-confidence subset): 主检验通过 + Fisher r-to-z p < fisher_p_threshold_high_confidence ####
## BH FDR (Fisher z 与 permutation) 仅作透明性输出, 不参与 hit 判定 (N=20 organ 下数学不可达) ####
delta_rho_threshold_main <- 0.3
delta_rho_threshold_supplementary <- c(0.2, 0.4)
perm_p_threshold <- 0.05
fisher_p_threshold_high_confidence <- 0.05

## pattern 顺序（与 6.1 / 7.6 一致） ####
pattern_group_order <- c(
  "Ocular",
  "Hematologic",
  "Neurologic",
  "Systemic-like",
  "Visceral-like"
)

## 并行设置 ####
detected_cores <- parallel::detectCores(logical = FALSE)
if (is.na(detected_cores)) {
  detected_cores <- 1
}
parallel_workers <- as.integer(max(1, min(8, detected_cores - 1)))
enable_parallel <- parallel_workers > 1 && .Platform$OS.type != "windows"

required_packages <- c(
  "matrixStats",
  "openxlsx",
  "data.table",
  "dplyr",
  "tidyr",
  "tibble",
  "purrr"
)

run_timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_file <- file.path(log_output_path, paste0("7.7模式相关基因差分相关识别_", run_timestamp, ".log"))


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

## 函数：load_single_object_from_rdata ####
## 功能：从单对象 Rdata 文件中读取指定对象 ####
## 输入：rdata_file，Rdata 文件路径；object_name，对象名 ####
## 输出：读取出的对象 ####
load_single_object_from_rdata <- function(rdata_file, object_name) {
  temp_env <- new.env(parent = emptyenv())
  load(rdata_file, envir = temp_env)

  if (!exists(object_name, envir = temp_env, inherits = FALSE)) {
    stop("Rdata 文件中不存在目标对象: ", object_name, " | ", rdata_file)
  }

  get(object_name, envir = temp_env, inherits = FALSE)
}

## 函数：load_named_objects_from_rdata ####
## 功能：从含多对象的 Rdata 文件中读取一个或多个对象，返回命名列表 ####
## 输入：rdata_file，Rdata 文件路径；object_names，对象名向量 ####
## 输出：命名列表，每个元素为对应对象 ####
load_named_objects_from_rdata <- function(rdata_file, object_names) {
  temp_env <- new.env(parent = emptyenv())
  load(rdata_file, envir = temp_env)

  missing_objects <- setdiff(object_names, ls(envir = temp_env))

  if (length(missing_objects) > 0) {
    stop(
      "Rdata 文件中不存在目标对象: ",
      paste(missing_objects, collapse = ", "),
      " | ",
      rdata_file
    )
  }

  result_list <- lapply(object_names, function(current_name) {
    get(current_name, envir = temp_env, inherits = FALSE)
  })
  names(result_list) <- object_names
  result_list
}

## 函数：check_required_columns ####
## 功能：检查数据框是否包含必要列，缺失则终止脚本 ####
## 输入：data_df，待检查数据框；required_columns，必要列名向量；data_label，数据说明 ####
## 输出：无返回值 ####
check_required_columns <- function(data_df, required_columns, data_label) {
  missing_columns <- setdiff(required_columns, colnames(data_df))

  if (length(missing_columns) > 0) {
    stop(data_label, " 缺少必要列: ", paste(missing_columns, collapse = ", "))
  }
}

## 函数：check_required_packages ####
## 功能：检查本脚本所需 R 包是否安装，缺失则终止脚本 ####
## 输入：packages，包名向量；output_path，依赖检查 Excel 输出路径 ####
## 输出：依赖检查数据框 ####
check_required_packages <- function(packages, output_path) {
  dependency_check_df <- data.frame(
    R包 = packages,
    已安装 = vapply(packages, requireNamespace, logical(1), quietly = TRUE),
    安装命令 = paste0("install.packages(\"", packages, "\")"),
    stringsAsFactors = FALSE
  )

  openxlsx::write.xlsx(
    x = list("依赖检查" = dependency_check_df),
    file = output_path,
    overwrite = TRUE
  )

  if (any(!dependency_check_df$已安装)) {
    missing_df <- dependency_check_df[!dependency_check_df$已安装, , drop = FALSE]
    missing_message <- paste0(
      missing_df$R包,
      ": ",
      missing_df$安装命令,
      collapse = "\n"
    )
    stop("缺少 7.7 分析必要 R 包，请手动安装后重新运行:\n", missing_message)
  }

  dependency_check_df
}

## 函数：run_parallel_lapply ####
## 功能：在 Linux/macOS 中用 mclapply 并行执行列表任务，无法并行时退回 lapply ####
## 输入：x，任务输入向量或列表；fun，任务函数；use_parallel；workers；...，传递给任务函数的其他参数 ####
## 输出：与 lapply 一致的列表结果 ####
run_parallel_lapply <- function(x, fun, use_parallel, workers, ...) {
  if (!use_parallel || workers < 2 || length(x) < 2 || .Platform$OS.type == "windows") {
    return(lapply(x, fun, ...))
  }

  parallel::mclapply(
    X = x,
    FUN = fun,
    ...,
    mc.cores = min(workers, length(x)),
    mc.preschedule = TRUE
  )
}

## 函数：safe_quantile ####
## 功能：在向量全为 NA 或长度为 0 时返回 NA，避免 quantile 报错 ####
## 输入：values，数值向量；probs，分位点 ####
## 输出：分位数（数值或 NA） ####
safe_quantile <- function(values, probs) {
  values <- as.numeric(values)
  values <- values[!is.na(values)]

  if (length(values) == 0) {
    return(NA_real_)
  }

  unname(stats::quantile(values, probs = probs, na.rm = TRUE, type = 7))
}


# 4. TCGA 肿瘤端输入加载 ####

## 函数：load_tcga_tumor_inputs ####
## 功能：读取 TCGA 肿瘤 TPM 表达矩阵和样本注释，并应用与 7.4/7.6 一致的样本筛选规则 ####
## 输入：expression_rdata；annotation_rdata ####
## 输出：含 tumor_tpm_matrix（gene × sample）与 tumor_meta_df 的列表 ####
load_tcga_tumor_inputs <- function(expression_rdata, annotation_rdata) {
  tumor_tpm_matrix <- load_single_object_from_rdata(
    rdata_file = expression_rdata,
    object_name = "tcga_tumor_expression_matrix"
  )
  tumor_annotation_df <- load_single_object_from_rdata(
    rdata_file = annotation_rdata,
    object_name = "tcga_tumor_sample_annotation_df"
  )

  check_required_columns(
    data_df = tumor_annotation_df,
    required_columns = c(
      "sample_barcode", "patient_id", "project_id",
      "primary_organ", "pattern_group", "sample_type"
    ),
    data_label = "TCGA 肿瘤样本注释"
  )

  tumor_meta_df <- tumor_annotation_df %>%
    dplyr::mutate(
      sample_col_id = as.character(sample_barcode),
      patient_id = as.character(patient_id),
      project_id = as.character(project_id),
      sample_type = as.character(sample_type),
      primary_organ = as.character(primary_organ),
      pattern_group = as.character(pattern_group)
    ) %>%
    dplyr::filter(
      sample_type == "Primary Tumor",
      sample_col_id %in% colnames(tumor_tpm_matrix),
      !is.na(primary_organ),
      !is.na(pattern_group)
    ) %>%
    dplyr::arrange(
      factor(pattern_group, levels = pattern_group_order),
      primary_organ, project_id, patient_id, sample_col_id
    )

  if (keep_first_tumor_per_patient && "sample_rank_within_patient" %in% colnames(tumor_meta_df)) {
    tumor_meta_df <- tumor_meta_df %>%
      dplyr::filter(is.na(sample_rank_within_patient) | sample_rank_within_patient == 1)
  }

  tumor_tpm_matrix <- tumor_tpm_matrix[, tumor_meta_df$sample_col_id, drop = FALSE]
  storage.mode(tumor_tpm_matrix) <- "numeric"

  list(
    tumor_tpm_matrix = tumor_tpm_matrix,
    tumor_meta_df = tumor_meta_df
  )
}


# 5. NMF H 矩阵加载 ####

## 函数：load_h_matrix_raw ####
## 功能：从 6.1 NMF best_fit 读取原始 H 矩阵，并整理为 primary_organ × pattern 数据矩阵 ####
## 输入：best_fit_rdata；h_cluster_rdata ####
## 输出：含 h_organ_pattern_matrix 与 cancer_pattern_groups 的列表 ####
## 说明：6.1.1 NMF 脚本只在内存中给 H 矩阵的行赋生物学名（"Hematologic"/"Ocular"/...），
##       并未将改名结果写回 best_fit；因此本函数需复刻 6.1.1 行 253-268 的硬编码命名约定，
##       避免直接读到 NMF 默认的空行名而 fallback 成 Pattern_1..5。 ####
load_h_matrix_raw <- function(best_fit_rdata, h_cluster_rdata) {
  best_fit <- load_single_object_from_rdata(
    rdata_file = best_fit_rdata,
    object_name = "best_fit"
  )

  if (!requireNamespace("NMF", quietly = TRUE)) {
    stop("缺少 NMF 包，无法从 best_fit 提取 coef(H)。")
  }

  h_raw_matrix <- NMF::coef(best_fit)
  h_raw_matrix <- as.matrix(h_raw_matrix)
  storage.mode(h_raw_matrix) <- "numeric"

  cluster_objects <- load_named_objects_from_rdata(
    rdata_file = h_cluster_rdata,
    object_names = c("cluster_to_pattern", "clusters", "cancer_pattern_groups")
  )

  cancer_pattern_groups <- cluster_objects$cancer_pattern_groups

  ## 复刻 6.1.1 的硬编码 pattern 命名约定（来源：Code/6.1.1癌种NMF分群分析.R 行 253-268） ####
  pattern_names_61_order <- c(
    "Hematologic",
    "Ocular",
    "Neurologic",
    "Systemic-like",
    "Visceral-like"
  )

  if (nrow(h_raw_matrix) != length(pattern_names_61_order)) {
    stop(
      "原始 H 矩阵行数 ", nrow(h_raw_matrix),
      " 与 6.1.1 预定义 pattern 数 ", length(pattern_names_61_order), " 不匹配 (rank=5 假设)"
    )
  }

  current_rownames <- rownames(h_raw_matrix)
  needs_relabel <- is.null(current_rownames) ||
    all(grepl("^Pattern_[0-9]+$", current_rownames)) ||
    all(grepl("^[0-9]+$", current_rownames)) ||
    !all(pattern_names_61_order %in% current_rownames)

  if (needs_relabel) {
    log_message(paste0(
      "H 矩阵行名为 ",
      if (is.null(current_rownames)) "NULL" else paste(current_rownames, collapse = ","),
      "; 按 6.1.1 硬编码顺序重命名为: ",
      paste(pattern_names_61_order, collapse = ", ")
    ))
    rownames(h_raw_matrix) <- pattern_names_61_order
  }

  pattern_row_names <- rownames(h_raw_matrix)
  pattern_match_index <- match(pattern_group_order, pattern_row_names)

  if (any(is.na(pattern_match_index))) {
    missing_patterns <- pattern_group_order[is.na(pattern_match_index)]
    stop(
      "原始 H 矩阵缺少预期模式名: ", paste(missing_patterns, collapse = ", "),
      "; 实际行名: ", paste(pattern_row_names, collapse = ", ")
    )
  }

  h_raw_matrix <- h_raw_matrix[pattern_match_index, , drop = FALSE]
  rownames(h_raw_matrix) <- pattern_group_order

  h_organ_pattern_matrix <- t(h_raw_matrix)
  colnames(h_organ_pattern_matrix) <- pattern_group_order

  list(
    h_organ_pattern_matrix = h_organ_pattern_matrix,
    cancer_pattern_groups = cancer_pattern_groups
  )
}


# 6. 正常端输入加载与 Tier 分类 ####

## 函数：load_normal_inputs ####
## 功能：读取 7.2 ComBat 校正后的正常组织 log2 矩阵 + 样本注释（GTEx + TCGA-normal） ####
## 输入：corrected_log2_rdata；gtex_annotation_rdata；tcga_normal_annotation_rdata ####
## 输出：含 normal_log2_matrix（gene × sample，已校正）与 normal_meta_df（含 data_source / primary_organ） 的列表 ####
load_normal_inputs <- function(corrected_log2_rdata,
                               gtex_annotation_rdata,
                               tcga_normal_annotation_rdata) {
  normal_log2_matrix <- load_single_object_from_rdata(
    rdata_file = corrected_log2_rdata,
    object_name = "corrected_log2_matrix"
  )
  storage.mode(normal_log2_matrix) <- "numeric"

  gtex_annotation_df <- load_single_object_from_rdata(
    rdata_file = gtex_annotation_rdata,
    object_name = "gtex_final_sample_annotation_df"
  )
  tcga_normal_annotation_df <- load_single_object_from_rdata(
    rdata_file = tcga_normal_annotation_rdata,
    object_name = "tcga_normal_sample_annotation_df"
  )

  gtex_meta_df <- data.frame(
    sample_col_id = as.character(gtex_annotation_df$sample_id),
    data_source = "GTEx",
    primary_organ = as.character(gtex_annotation_df$primary_organ),
    stringsAsFactors = FALSE
  )

  tcga_normal_meta_df <- data.frame(
    sample_col_id = as.character(tcga_normal_annotation_df$sample_barcode),
    data_source = "TCGA",
    primary_organ = as.character(tcga_normal_annotation_df$primary_organ),
    stringsAsFactors = FALSE
  )

  normal_meta_df <- rbind(gtex_meta_df, tcga_normal_meta_df) %>%
    dplyr::filter(
      !is.na(primary_organ),
      sample_col_id %in% colnames(normal_log2_matrix)
    ) %>%
    dplyr::arrange(primary_organ, data_source, sample_col_id)

  normal_log2_matrix <- normal_log2_matrix[, normal_meta_df$sample_col_id, drop = FALSE]

  list(
    normal_log2_matrix = normal_log2_matrix,
    normal_meta_df = normal_meta_df
  )
}

## 函数：classify_organ_tier ####
## 功能：按肿瘤侧样本量与正常侧 (TCGA-normal, GTEx) 样本量为每个 primary_organ 分配 Tier ####
## 输入：tumor_meta_df；normal_meta_df；min_tumor_n；tcga_normal_min；gtex_min ####
## 输出：data.frame，每行一个 organ，含 n_tumor / n_tcga_normal / n_gtex_normal / tier / normal_baseline_source ####
classify_organ_tier <- function(tumor_meta_df,
                                normal_meta_df,
                                min_tumor_n,
                                tcga_normal_min,
                                gtex_min) {
  tumor_count_df <- tumor_meta_df %>%
    dplyr::count(primary_organ, name = "n_tumor_samples")

  normal_count_df <- normal_meta_df %>%
    dplyr::count(primary_organ, data_source) %>%
    tidyr::pivot_wider(
      names_from = data_source,
      values_from = n,
      values_fill = 0L
    )

  for (col_name in c("TCGA", "GTEx")) {
    if (!col_name %in% colnames(normal_count_df)) {
      normal_count_df[[col_name]] <- 0L
    }
  }

  organ_tier_df <- tumor_count_df %>%
    dplyr::full_join(normal_count_df, by = "primary_organ") %>%
    dplyr::mutate(
      n_tumor_samples = ifelse(is.na(n_tumor_samples), 0L, as.integer(n_tumor_samples)),
      n_tcga_normal = as.integer(ifelse(is.na(TCGA), 0L, TCGA)),
      n_gtex_normal = as.integer(ifelse(is.na(GTEx), 0L, GTEx))
    ) %>%
    dplyr::select(primary_organ, n_tumor_samples, n_tcga_normal, n_gtex_normal) %>%
    dplyr::mutate(
      tumor_qualified = n_tumor_samples >= min_tumor_n,
      tier = dplyr::case_when(
        n_tcga_normal >= tcga_normal_min ~ "Tier1",
        n_tcga_normal > 0 & n_gtex_normal >= gtex_min ~ "Tier1b",
        n_tcga_normal == 0 & n_gtex_normal >= gtex_min ~ "Tier2",
        TRUE ~ "Tier3"
      ),
      normal_baseline_source = dplyr::case_when(
        tier == "Tier1"  ~ "TCGA-normal+GTEx",
        tier == "Tier1b" ~ "GTEx-only",
        tier == "Tier2"  ~ "GTEx-only",
        tier == "Tier3"  ~ "无（不计算ρ_normal）"
      ),
      enters_rho_tumor  = tumor_qualified,
      enters_rho_normal = tumor_qualified & tier %in% c("Tier1", "Tier1b", "Tier2")
    ) %>%
    dplyr::arrange(
      dplyr::desc(enters_rho_tumor),
      factor(tier, levels = c("Tier1", "Tier1b", "Tier2", "Tier3")),
      dplyr::desc(n_tumor_samples)
    )

  organ_tier_df
}


# 7. organ-median 表达矩阵计算 ####

## 函数：compute_tumor_organ_median_matrices ####
## 功能：将 TCGA 肿瘤 TPM 矩阵按 primary_organ 取中位，得到 organ-median TPM 与 organ-median log2(TPM+1) ####
## 输入：tumor_tpm_matrix；tumor_meta_df；included_primary_organs；log2_offset ####
## 输出：含 organ_median_tpm_matrix 与 organ_median_log2_matrix（gene × organ）的列表 ####
compute_tumor_organ_median_matrices <- function(tumor_tpm_matrix,
                                                tumor_meta_df,
                                                included_primary_organs,
                                                log2_offset) {
  gene_names <- rownames(tumor_tpm_matrix)

  organ_median_tpm_matrix <- matrix(
    NA_real_,
    nrow = length(gene_names),
    ncol = length(included_primary_organs),
    dimnames = list(gene_names, included_primary_organs)
  )
  organ_median_log2_matrix <- matrix(
    NA_real_,
    nrow = length(gene_names),
    ncol = length(included_primary_organs),
    dimnames = list(gene_names, included_primary_organs)
  )

  for (current_primary_organ in included_primary_organs) {
    current_samples <- tumor_meta_df$sample_col_id[tumor_meta_df$primary_organ == current_primary_organ]
    current_tpm_matrix <- tumor_tpm_matrix[, current_samples, drop = FALSE]
    storage.mode(current_tpm_matrix) <- "numeric"

    organ_median_tpm_matrix[, current_primary_organ] <- matrixStats::rowMedians(
      current_tpm_matrix, na.rm = TRUE
    )
    organ_median_log2_matrix[, current_primary_organ] <- matrixStats::rowMedians(
      log2(current_tpm_matrix + log2_offset), na.rm = TRUE
    )
  }

  list(
    organ_median_tpm_matrix = organ_median_tpm_matrix,
    organ_median_log2_matrix = organ_median_log2_matrix
  )
}

## 函数：compute_normal_organ_median_matrix ####
## 功能：基于 ComBat 校正后正常 log2 矩阵 + Tier 分类，按 organ 选取相应数据源（TCGA-normal+GTEx 或 GTEx-only）取 median ####
## 输入：normal_log2_matrix（gene × sample，已 log2，已 ComBat）；normal_meta_df；organ_tier_df；included_primary_organs ####
## 输出：含 organ_median_log2_normal_matrix（gene × organ）与 normal_organ_sample_summary_df（每 organ 实际使用的样本数与来源）的列表 ####
compute_normal_organ_median_matrix <- function(normal_log2_matrix,
                                               normal_meta_df,
                                               organ_tier_df,
                                               included_primary_organs) {
  gene_names <- rownames(normal_log2_matrix)
  organs_with_normal <- organ_tier_df$primary_organ[organ_tier_df$enters_rho_normal]
  organs_with_normal <- intersect(organs_with_normal, included_primary_organs)

  organ_median_log2_normal_matrix <- matrix(
    NA_real_,
    nrow = length(gene_names),
    ncol = length(organs_with_normal),
    dimnames = list(gene_names, organs_with_normal)
  )

  sample_summary_list <- vector("list", length(organs_with_normal))

  for (idx in seq_along(organs_with_normal)) {
    current_organ <- organs_with_normal[[idx]]
    current_tier_row <- organ_tier_df[organ_tier_df$primary_organ == current_organ, , drop = FALSE]
    current_tier <- as.character(current_tier_row$tier)

    if (current_tier == "Tier1") {
      keep_sources <- c("TCGA", "GTEx")
    } else {
      keep_sources <- "GTEx"
    }

    current_samples <- normal_meta_df$sample_col_id[
      normal_meta_df$primary_organ == current_organ &
        normal_meta_df$data_source %in% keep_sources
    ]
    current_samples <- intersect(current_samples, colnames(normal_log2_matrix))

    if (length(current_samples) == 0) {
      sample_summary_list[[idx]] <- data.frame(
        primary_organ = current_organ,
        tier = current_tier,
        normal_baseline_source = current_tier_row$normal_baseline_source,
        n_normal_samples_used = 0L,
        n_tcga_used = 0L,
        n_gtex_used = 0L,
        stringsAsFactors = FALSE
      )
      next
    }

    current_log2_matrix <- normal_log2_matrix[, current_samples, drop = FALSE]
    storage.mode(current_log2_matrix) <- "numeric"

    organ_median_log2_normal_matrix[, current_organ] <- matrixStats::rowMedians(
      current_log2_matrix, na.rm = TRUE
    )

    sample_summary_list[[idx]] <- data.frame(
      primary_organ = current_organ,
      tier = current_tier,
      normal_baseline_source = current_tier_row$normal_baseline_source,
      n_normal_samples_used = length(current_samples),
      n_tcga_used = sum(normal_meta_df$sample_col_id %in% current_samples &
                          normal_meta_df$data_source == "TCGA"),
      n_gtex_used = sum(normal_meta_df$sample_col_id %in% current_samples &
                          normal_meta_df$data_source == "GTEx"),
      stringsAsFactors = FALSE
    )
  }

  normal_organ_sample_summary_df <- dplyr::bind_rows(sample_summary_list)

  list(
    organ_median_log2_normal_matrix = organ_median_log2_normal_matrix,
    normal_organ_sample_summary_df = normal_organ_sample_summary_df
  )
}


# 8. 基因预过滤（仅保留计算可行性：低表达 + 低方差 + 正常端覆盖） ####

## 函数：apply_gene_prefilter ####
## 功能：对全基因依次应用低表达、低方差两道计算可行性过滤，并要求基因同时存在于肿瘤端与正常端矩阵以保证 Δρ 可计算；不再做组织特异性过滤（Δρ 框架本身已通过正常端阴性对照消除 tissue-of-origin 混杂） ####
## 输入：organ_median_log2_tumor_matrix（gene × tumor_organ）；organ_median_log2_normal_matrix（gene × normal_organ）；min_max_log2，最大 organ-median log2 阈值；min_mad_log2，肿瘤端 organ-median log2 的 MAD 阈值 ####
## 输出：含 retained_genes（向量）、filter_summary_df、filter_step_count_df 的列表 ####
apply_gene_prefilter <- function(organ_median_log2_tumor_matrix,
                                 organ_median_log2_normal_matrix,
                                 min_max_log2,
                                 min_mad_log2) {
  gene_names <- rownames(organ_median_log2_tumor_matrix)

  max_log2_per_gene <- matrixStats::rowMaxs(organ_median_log2_tumor_matrix, na.rm = TRUE)
  max_log2_per_gene[!is.finite(max_log2_per_gene)] <- NA_real_

  mad_log2_per_gene <- apply(organ_median_log2_tumor_matrix, 1, function(current_values) {
    valid_values <- current_values[!is.na(current_values)]
    if (length(valid_values) < 3) {
      return(NA_real_)
    }
    stats::mad(valid_values, constant = 1.4826, na.rm = TRUE)
  })

  in_normal_matrix <- gene_names %in% rownames(organ_median_log2_normal_matrix)

  filter_summary_df <- data.frame(
    gene = gene_names,
    max_organ_median_log2 = max_log2_per_gene,
    mad_organ_median_log2 = mad_log2_per_gene,
    in_normal_matrix = in_normal_matrix,
    stringsAsFactors = FALSE
  )

  filter_summary_df <- filter_summary_df %>%
    dplyr::mutate(
      fail_low_expression = is.na(max_organ_median_log2) | max_organ_median_log2 < min_max_log2,
      fail_low_variance = is.na(mad_organ_median_log2) | mad_organ_median_log2 < min_mad_log2,
      fail_not_in_normal = !in_normal_matrix,
      retained = !(fail_low_expression | fail_low_variance | fail_not_in_normal)
    )

  filter_step_count_df <- data.frame(
    步骤 = c(
      "全部基因",
      "失败_低表达(肿瘤端)",
      "失败_低方差(肿瘤端)",
      "失败_正常端无该基因",
      "保留基因"
    ),
    基因数 = c(
      nrow(filter_summary_df),
      sum(filter_summary_df$fail_low_expression, na.rm = TRUE),
      sum(filter_summary_df$fail_low_variance &
            !filter_summary_df$fail_low_expression, na.rm = TRUE),
      sum(filter_summary_df$fail_not_in_normal &
            !filter_summary_df$fail_low_expression &
            !filter_summary_df$fail_low_variance, na.rm = TRUE),
      sum(filter_summary_df$retained, na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
  )

  retained_genes <- filter_summary_df$gene[filter_summary_df$retained]

  list(
    retained_genes = retained_genes,
    filter_summary_df = filter_summary_df,
    filter_step_count_df = filter_step_count_df
  )
}


# 9. ρ_tumor / ρ_normal / Δρ 计算 ####

## 函数：compute_observed_correlation_long ####
## 功能：计算每个基因 × 每个模式的观察 Spearman 相关，返回长表（含 organ 共享数与 rho_label 列） ####
## 输入：organ_median_log2_matrix（gene × organ）；h_organ_pattern_matrix（organ × pattern）；side_label，"tumor" 或 "normal" ####
## 输出：长表数据框，含 gene、pattern、rho_observed、n_organ、side ####
compute_observed_correlation_long <- function(organ_median_log2_matrix,
                                              h_organ_pattern_matrix,
                                              side_label) {
  shared_organs <- intersect(colnames(organ_median_log2_matrix), rownames(h_organ_pattern_matrix))

  if (length(shared_organs) < 5) {
    stop("organ_median 与 H 矩阵共享的 primary_organ 少于 5 个，无法计算相关。 side: ", side_label)
  }

  organ_median_aligned_matrix <- organ_median_log2_matrix[, shared_organs, drop = FALSE]
  h_organ_pattern_aligned_matrix <- h_organ_pattern_matrix[shared_organs, , drop = FALSE]

  rho_matrix <- suppressWarnings(stats::cor(
    x = t(organ_median_aligned_matrix),
    y = h_organ_pattern_aligned_matrix,
    method = "spearman",
    use = "pairwise.complete.obs"
  ))

  observed_long_df <- as.data.frame(rho_matrix, check.names = FALSE) %>%
    tibble::rownames_to_column("gene") %>%
    tidyr::pivot_longer(
      cols = -gene,
      names_to = "pattern",
      values_to = "rho_observed"
    ) %>%
    dplyr::mutate(
      pattern = factor(pattern, levels = pattern_group_order),
      n_organ = length(shared_organs),
      side = side_label
    ) %>%
    dplyr::arrange(gene, pattern)

  observed_long_df
}

## 函数：compute_delta_rho_with_fisher_z ####
## 功能：合并 ρ_tumor 与 ρ_normal，计算 Δρ 并以 Fisher r-to-z 变换做独立样本相关比较 ####
## 输入：rho_tumor_df（含 gene/pattern/rho_observed/n_organ）；rho_normal_df（同）；fdr_method，p 校正方法 ####
## 输出：长表数据框，含 gene、pattern、rho_tumor、rho_normal、delta_rho、fisher_z、fisher_p_value、fisher_fdr、n_organ_tumor、n_organ_normal ####
compute_delta_rho_with_fisher_z <- function(rho_tumor_df, rho_normal_df, fdr_method = "BH") {
  rho_tumor_subset <- rho_tumor_df %>%
    dplyr::transmute(
      gene, pattern,
      rho_tumor = rho_observed,
      n_organ_tumor = n_organ
    )

  rho_normal_subset <- rho_normal_df %>%
    dplyr::transmute(
      gene, pattern,
      rho_normal = rho_observed,
      n_organ_normal = n_organ
    )

  delta_rho_long_df <- rho_tumor_subset %>%
    dplyr::left_join(rho_normal_subset, by = c("gene", "pattern")) %>%
    dplyr::mutate(
      delta_rho = rho_tumor - rho_normal,
      z_tumor = atanh(pmin(pmax(rho_tumor, -0.999999), 0.999999)),
      z_normal = atanh(pmin(pmax(rho_normal, -0.999999), 0.999999)),
      se_diff = sqrt(
        ifelse(n_organ_tumor > 3, 1 / (n_organ_tumor - 3), NA_real_) +
          ifelse(n_organ_normal > 3, 1 / (n_organ_normal - 3), NA_real_)
      ),
      fisher_z = ifelse(
        is.na(rho_tumor) | is.na(rho_normal) | is.na(se_diff) | se_diff == 0,
        NA_real_,
        (z_tumor - z_normal) / se_diff
      ),
      fisher_p_value = ifelse(
        is.na(fisher_z),
        NA_real_,
        2 * stats::pnorm(-abs(fisher_z))
      )
    )

  delta_rho_long_df <- delta_rho_long_df %>%
    dplyr::group_by(pattern) %>%
    dplyr::mutate(
      fisher_fdr = stats::p.adjust(fisher_p_value, method = fdr_method)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(
      gene, pattern,
      rho_tumor, rho_normal, delta_rho,
      n_organ_tumor, n_organ_normal,
      fisher_z, fisher_p_value, fisher_fdr
    ) %>%
    dplyr::arrange(gene, pattern)

  delta_rho_long_df
}


# 10. Permutation 灵敏性分析：matched 子集 Δρ 零分布 ####

## 函数：run_single_permutation_iteration ####
## 功能：单次 permutation 迭代：打乱 H 矩阵的 organ 行，重算 Spearman 相关 ####
## 输入：iteration_seed；organ_median_log2_matrix；h_organ_pattern_matrix；shared_organs，相关计算用到的 organ 顺序 ####
## 输出：n_gene × n_pattern 数值矩阵 ####
run_single_permutation_iteration <- function(iteration_seed,
                                             organ_median_log2_matrix,
                                             h_organ_pattern_matrix,
                                             shared_organs) {
  set.seed(iteration_seed)

  permuted_indices <- sample(seq_along(shared_organs))
  shared_organs_permuted <- shared_organs[permuted_indices]

  h_permuted_matrix <- h_organ_pattern_matrix[shared_organs_permuted, , drop = FALSE]

  organ_aligned_matrix <- organ_median_log2_matrix[, shared_organs, drop = FALSE]

  rho_matrix <- suppressWarnings(stats::cor(
    x = t(organ_aligned_matrix),
    y = h_permuted_matrix,
    method = "spearman",
    use = "pairwise.complete.obs"
  ))

  rho_matrix
}

## 函数：run_permutation_correlations ####
## 功能：执行 n_iterations 次 organ 标签 permutation，返回 iteration × gene × pattern 数组 ####
## 输入：organ_median_log2_matrix；h_organ_pattern_matrix；n_iterations；seed；use_parallel；workers ####
## 输出：3 维数组（H ⟂ expression 假设下的 rho 零分布） ####
run_permutation_correlations <- function(organ_median_log2_matrix,
                                         h_organ_pattern_matrix,
                                         n_iterations,
                                         seed,
                                         use_parallel,
                                         workers) {
  shared_organs <- intersect(colnames(organ_median_log2_matrix), rownames(h_organ_pattern_matrix))
  if (length(shared_organs) < 5) {
    stop("organ_median 与 H 矩阵共享 organ < 5，无法 permutation。")
  }

  set.seed(seed + 1L)
  iteration_seed_vector <- sample.int(.Machine$integer.max, n_iterations)
  iteration_indices <- seq_len(n_iterations)

  log_message(paste0(
    "执行 permutation 相关分析: ", n_iterations, " 次 organ 标签打乱, ",
    "并行: ", use_parallel, " (workers=", workers, ")"
  ))

  permutation_iteration_list <- run_parallel_lapply(
    x = iteration_indices,
    fun = function(iter_idx) {
      run_single_permutation_iteration(
        iteration_seed = iteration_seed_vector[iter_idx],
        organ_median_log2_matrix = organ_median_log2_matrix,
        h_organ_pattern_matrix = h_organ_pattern_matrix,
        shared_organs = shared_organs
      )
    },
    use_parallel = use_parallel,
    workers = workers
  )

  gene_names <- rownames(organ_median_log2_matrix)
  pattern_names_local <- colnames(h_organ_pattern_matrix)

  permutation_array <- array(
    NA_real_,
    dim = c(n_iterations, length(gene_names), length(pattern_names_local)),
    dimnames = list(
      iteration = paste0("perm_", iteration_indices),
      gene = gene_names,
      pattern = pattern_names_local
    )
  )

  for (i in iteration_indices) {
    permutation_array[i, , ] <- permutation_iteration_list[[i]]
  }

  permutation_array
}

## 函数：compute_delta_rho_permutation_pvalue ####
## 功能：基于双侧 permutation 零分布对 Δρ 做正交显著性验证（不强行依赖 Fisher z 独立性假设） ####
## 输入：observed_delta_df，含 gene/pattern/delta_rho；tumor_perm_array；normal_perm_array（同样 iteration 数与同样的种子打散 H）；fdr_method ####
## 输出：长表数据框，含 gene、pattern、delta_perm_p_value、delta_perm_fdr ####
compute_delta_rho_permutation_pvalue <- function(observed_delta_df,
                                                 tumor_perm_array,
                                                 normal_perm_array,
                                                 fdr_method = "BH") {
  gene_names <- dimnames(tumor_perm_array)[["gene"]]
  pattern_names_local <- dimnames(tumor_perm_array)[["pattern"]]

  if (!identical(dimnames(normal_perm_array)[["gene"]], gene_names)) {
    stop("肿瘤端与正常端 permutation 数组基因维度不一致。")
  }
  if (!identical(dim(tumor_perm_array)[1], dim(normal_perm_array)[1])) {
    stop("肿瘤端与正常端 permutation 迭代数不一致。")
  }

  n_iter <- dim(tumor_perm_array)[[1]]

  observed_wide_df <- observed_delta_df %>%
    dplyr::select(gene, pattern, delta_rho) %>%
    tidyr::pivot_wider(names_from = pattern, values_from = delta_rho)
  observed_wide_df <- observed_wide_df[match(gene_names, observed_wide_df$gene), , drop = FALSE]

  observed_matrix <- as.matrix(observed_wide_df[, pattern_names_local, drop = FALSE])
  rownames(observed_matrix) <- gene_names

  pvalue_list <- vector("list", length(pattern_names_local))

  for (pattern_index in seq_along(pattern_names_local)) {
    current_pattern <- pattern_names_local[[pattern_index]]
    current_observed <- observed_matrix[, current_pattern]

    delta_null_slice <- tumor_perm_array[, , pattern_index, drop = TRUE] -
      normal_perm_array[, , pattern_index, drop = TRUE]

    abs_observed <- abs(current_observed)
    abs_null <- abs(delta_null_slice)

    extreme_compare_matrix <- sweep(abs_null, 2, abs_observed, FUN = ">=")
    extreme_count_per_gene <- colSums(extreme_compare_matrix, na.rm = TRUE)

    delta_perm_p_value <- (extreme_count_per_gene + 1) / (n_iter + 1)
    delta_perm_fdr <- stats::p.adjust(delta_perm_p_value, method = fdr_method)

    pvalue_list[[pattern_index]] <- data.frame(
      gene = gene_names,
      pattern = current_pattern,
      delta_perm_p_value = delta_perm_p_value,
      delta_perm_fdr = delta_perm_fdr,
      stringsAsFactors = FALSE
    )
  }

  delta_pvalue_df <- dplyr::bind_rows(pvalue_list) %>%
    dplyr::mutate(pattern = factor(pattern, levels = pattern_group_order)) %>%
    dplyr::arrange(gene, pattern)

  delta_pvalue_df
}


# 11. 多阈值筛选与多模式分类 ####

## 函数：apply_delta_rho_filter ####
## 功能：在 Δρ 多阈值下做 hit 判定 ####
##       主筛选 hit_flag: |Δρ| ≥ 阈值 且 matched-pair permutation 双尾经验 p < perm_p_threshold ####
##       高置信 hit_flag_high_confidence: 主筛选通过 且 Fisher r-to-z p < fisher_p_threshold_high_confidence ####
## 输入：merged_long_df，已合并 Δρ / Fisher z / Δρ permutation 的长表；
##       delta_rho_thresholds；perm_p_threshold；fisher_p_threshold_high_confidence ####
## 输出：长表数据框，每行为 gene × pattern × delta_rho_threshold，含 hit_flag、hit_flag_high_confidence、hit_direction ####
apply_delta_rho_filter <- function(merged_long_df,
                                   delta_rho_thresholds,
                                   perm_p_threshold,
                                   fisher_p_threshold_high_confidence) {
  threshold_dfs <- lapply(delta_rho_thresholds, function(current_delta_t) {
    merged_long_df %>%
      dplyr::mutate(
        delta_rho_threshold = current_delta_t,
        passes_delta_rho_threshold = !is.na(delta_rho) & abs(delta_rho) >= current_delta_t,
        passes_perm_p = !is.na(delta_perm_p_value) & delta_perm_p_value < perm_p_threshold,
        passes_fisher_p = !is.na(fisher_p_value) & fisher_p_value < fisher_p_threshold_high_confidence,
        hit_flag = passes_delta_rho_threshold & passes_perm_p,
        hit_flag_high_confidence = hit_flag & passes_fisher_p,
        hit_direction = dplyr::case_when(
          !hit_flag ~ NA_character_,
          delta_rho > 0 ~ "positive",
          delta_rho < 0 ~ "negative",
          TRUE ~ NA_character_
        )
      )
  })

  dplyr::bind_rows(threshold_dfs) %>%
    dplyr::arrange(delta_rho_threshold, gene, pattern)
}

## 函数：build_per_pattern_hit_lists ####
## 功能：基于 apply_delta_rho_filter 输出，整理每模式 hit 基因列表与计数表 ####
##       同时输出主筛选 hit_long_df 与 high-confidence 子集 hit_long_high_confidence_df ####
##       计数表内列出 main / high-confidence 两套计数 ####
## 输入：threshold_filter_df ####
## 输出：含 hit_count_df、hit_long_df、hit_long_high_confidence_df 的列表 ####
build_per_pattern_hit_lists <- function(threshold_filter_df) {
  hit_long_df <- threshold_filter_df %>%
    dplyr::filter(hit_flag) %>%
    dplyr::transmute(
      delta_rho_threshold,
      pattern = as.character(pattern),
      gene, hit_direction,
      hit_high_confidence = hit_flag_high_confidence,
      rho_tumor, rho_normal, delta_rho,
      fisher_z, fisher_p_value, fisher_fdr,
      delta_perm_p_value, delta_perm_fdr,
      n_organ_tumor, n_organ_normal
    ) %>%
    dplyr::arrange(
      delta_rho_threshold,
      factor(pattern, levels = pattern_group_order),
      hit_direction,
      dplyr::desc(abs(delta_rho))
    )

  hit_long_high_confidence_df <- hit_long_df %>%
    dplyr::filter(hit_high_confidence)

  hit_count_df <- threshold_filter_df %>%
    dplyr::group_by(delta_rho_threshold, pattern) %>%
    dplyr::summarise(
      n_positive_hits = sum(hit_flag & hit_direction == "positive", na.rm = TRUE),
      n_negative_hits = sum(hit_flag & hit_direction == "negative", na.rm = TRUE),
      n_total_hits = n_positive_hits + n_negative_hits,
      n_positive_high_confidence = sum(hit_flag_high_confidence & hit_direction == "positive", na.rm = TRUE),
      n_negative_high_confidence = sum(hit_flag_high_confidence & hit_direction == "negative", na.rm = TRUE),
      n_total_high_confidence = n_positive_high_confidence + n_negative_high_confidence,
      .groups = "drop"
    ) %>%
    dplyr::arrange(
      delta_rho_threshold,
      factor(pattern, levels = pattern_group_order)
    )

  list(
    hit_long_df = hit_long_df,
    hit_long_high_confidence_df = hit_long_high_confidence_df,
    hit_count_df = hit_count_df
  )
}

## 函数：classify_multi_pattern_genes ####
## 功能：在主阈值下按 sign(Δρ) 对每个基因构建 5 维 hit 状态，分类为 uni-pattern / multi-pattern-same-direction / anti-correlated / non-specific ####
## 输入：threshold_filter_df；delta_rho_threshold_target ####
## 输出：每基因一行的数据框，含 5 个模式的 ρ_tumor / ρ_normal / Δρ / hit_status 与 gene_category ####
classify_multi_pattern_genes <- function(threshold_filter_df,
                                         delta_rho_threshold_target) {
  target_df <- threshold_filter_df %>%
    dplyr::filter(delta_rho_threshold == delta_rho_threshold_target) %>%
    dplyr::select(
      gene, pattern,
      rho_tumor, rho_normal, delta_rho,
      fisher_p_value, fisher_fdr,
      delta_perm_p_value, delta_perm_fdr,
      hit_flag, hit_flag_high_confidence, hit_direction
    )

  pivot_to_wide <- function(df, value_col, prefix) {
    wide_df <- df %>%
      dplyr::select(gene, pattern, dplyr::all_of(value_col)) %>%
      tidyr::pivot_wider(names_from = pattern, values_from = dplyr::all_of(value_col))
    colnames(wide_df)[-1] <- paste0(prefix, "_", colnames(wide_df)[-1])
    wide_df
  }

  rho_tumor_wide_df <- pivot_to_wide(target_df, "rho_tumor", "rho_tumor")
  rho_normal_wide_df <- pivot_to_wide(target_df, "rho_normal", "rho_normal")
  delta_rho_wide_df <- pivot_to_wide(target_df, "delta_rho", "delta_rho")
  perm_p_wide_df <- pivot_to_wide(target_df, "delta_perm_p_value", "perm_p")
  fisher_p_wide_df <- pivot_to_wide(target_df, "fisher_p_value", "fisher_p")
  fdr_delta_wide_df <- pivot_to_wide(target_df, "fisher_fdr", "fdr_delta_fisher")
  fdr_delta_perm_wide_df <- pivot_to_wide(target_df, "delta_perm_fdr", "fdr_delta_perm")

  hit_wide_df <- target_df %>%
    dplyr::mutate(
      hit_status = dplyr::case_when(
        !hit_flag ~ "none",
        hit_direction == "positive" ~ "positive",
        hit_direction == "negative" ~ "negative",
        TRUE ~ "none"
      )
    ) %>%
    dplyr::select(gene, pattern, hit_status) %>%
    tidyr::pivot_wider(names_from = pattern, values_from = hit_status)
  colnames(hit_wide_df)[-1] <- paste0("hit_", colnames(hit_wide_df)[-1])

  hit_high_confidence_wide_df <- target_df %>%
    dplyr::mutate(
      hit_high_confidence_status = dplyr::case_when(
        !hit_flag_high_confidence ~ "none",
        hit_direction == "positive" ~ "positive",
        hit_direction == "negative" ~ "negative",
        TRUE ~ "none"
      )
    ) %>%
    dplyr::select(gene, pattern, hit_high_confidence_status) %>%
    tidyr::pivot_wider(names_from = pattern, values_from = hit_high_confidence_status)
  colnames(hit_high_confidence_wide_df)[-1] <- paste0("hit_hc_", colnames(hit_high_confidence_wide_df)[-1])

  pattern_count_df <- target_df %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(
      n_positive_patterns = sum(hit_flag & hit_direction == "positive", na.rm = TRUE),
      n_negative_patterns = sum(hit_flag & hit_direction == "negative", na.rm = TRUE),
      n_high_confidence_patterns = sum(hit_flag_high_confidence, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      n_total_patterns = n_positive_patterns + n_negative_patterns,
      gene_category = dplyr::case_when(
        n_total_patterns == 0 ~ "non-specific",
        n_total_patterns == 1 ~ "uni-pattern",
        n_positive_patterns >= 1 & n_negative_patterns >= 1 ~ "anti-correlated",
        TRUE ~ "multi-pattern-same-direction"
      )
    )

  classification_df <- rho_tumor_wide_df %>%
    dplyr::left_join(rho_normal_wide_df, by = "gene") %>%
    dplyr::left_join(delta_rho_wide_df, by = "gene") %>%
    dplyr::left_join(perm_p_wide_df, by = "gene") %>%
    dplyr::left_join(fisher_p_wide_df, by = "gene") %>%
    dplyr::left_join(fdr_delta_wide_df, by = "gene") %>%
    dplyr::left_join(fdr_delta_perm_wide_df, by = "gene") %>%
    dplyr::left_join(hit_wide_df, by = "gene") %>%
    dplyr::left_join(hit_high_confidence_wide_df, by = "gene") %>%
    dplyr::left_join(pattern_count_df, by = "gene") %>%
    dplyr::arrange(
      factor(
        gene_category,
        levels = c("uni-pattern", "multi-pattern-same-direction", "anti-correlated", "non-specific")
      ),
      dplyr::desc(n_total_patterns),
      gene
    )

  classification_df
}


# 12. 主流程:数据加载 ####

log_message("========== 7.7 模式相关基因差分相关识别开始 ==========")
log_message(paste0(
  "并行设置: enable_parallel=", enable_parallel,
  "; parallel_workers=", parallel_workers
))

assert_required_files(
  file_paths = c(
    tcga_tumor_expression_rdata,
    tcga_tumor_sample_annotation_rdata,
    tcga_normal_sample_annotation_rdata,
    gtex_normal_sample_annotation_rdata,
    corrected_normal_log2_rdata,
    nmf_best_fit_rdata,
    nmf_h_cluster_rdata
  ),
  file_labels = c(
    "TCGA 肿瘤 TPM 矩阵(原始)",
    "TCGA 肿瘤样本注释",
    "TCGA 邻近正常样本注释",
    "GTEx 正常样本注释",
    "7.2 ComBat 校正后的正常 log2 矩阵",
    "6.1 NMF best_fit",
    "6.1 NMF H_cluster"
  )
)

dependency_check_df <- check_required_packages(
  packages = required_packages,
  output_path = file.path(table_output_path, "依赖检查.xlsx")
)
log_message("依赖检查通过。")

## 13.1 加载肿瘤端表达与注释 ####
log_message("读取 TCGA 肿瘤 TPM 矩阵和样本注释...")
tcga_inputs <- load_tcga_tumor_inputs(
  expression_rdata = tcga_tumor_expression_rdata,
  annotation_rdata = tcga_tumor_sample_annotation_rdata
)
tumor_tpm_matrix <- tcga_inputs$tumor_tpm_matrix
tumor_meta_df <- tcga_inputs$tumor_meta_df

log_message(paste0(
  "样本筛选后: ", nrow(tumor_meta_df), " 个 Primary Tumor 样本; ",
  length(unique(tumor_meta_df$primary_organ)), " 个 primary_organ"
))

## 13.2 加载 NMF H 矩阵 ####
log_message("读取 6.1 NMF 原始 H 矩阵...")
h_inputs <- load_h_matrix_raw(
  best_fit_rdata = nmf_best_fit_rdata,
  h_cluster_rdata = nmf_h_cluster_rdata
)
h_organ_pattern_matrix <- h_inputs$h_organ_pattern_matrix
cancer_pattern_groups <- h_inputs$cancer_pattern_groups

log_message(paste0(
  "原始 H 矩阵维度: ", nrow(h_organ_pattern_matrix), " primary_organ × ",
  ncol(h_organ_pattern_matrix), " pattern"
))

## 13.3 加载正常端表达与注释 ####
log_message("读取 7.2 ComBat 校正后的 GTEx + TCGA-normal 联合 log2 矩阵...")
normal_inputs <- load_normal_inputs(
  corrected_log2_rdata = corrected_normal_log2_rdata,
  gtex_annotation_rdata = gtex_normal_sample_annotation_rdata,
  tcga_normal_annotation_rdata = tcga_normal_sample_annotation_rdata
)
normal_log2_matrix <- normal_inputs$normal_log2_matrix
normal_meta_df <- normal_inputs$normal_meta_df

log_message(paste0(
  "正常端样本数: ", nrow(normal_meta_df),
  " (GTEx ", sum(normal_meta_df$data_source == "GTEx"),
  " + TCGA-normal ", sum(normal_meta_df$data_source == "TCGA"), ")"
))


# 13. organ Tier 分类与对齐 ####

log_message("按肿瘤侧样本量与正常侧 (TCGA + GTEx) 样本量对每个 primary_organ 分配 Tier...")
organ_tier_full_df <- classify_organ_tier(
  tumor_meta_df = tumor_meta_df,
  normal_meta_df = normal_meta_df,
  min_tumor_n = min_samples_per_organ_tumor,
  tcga_normal_min = tcga_normal_min_for_tier1,
  gtex_min = gtex_min_for_normal_baseline
)

## 仅保留肿瘤端达标且在 H 矩阵中存在的 organ ####
h_organ_levels <- rownames(h_organ_pattern_matrix)
organ_tier_full_df <- organ_tier_full_df %>%
  dplyr::mutate(
    in_h_matrix = primary_organ %in% h_organ_levels,
    enters_rho_tumor = enters_rho_tumor & in_h_matrix,
    enters_rho_normal = enters_rho_normal & in_h_matrix
  )

included_tumor_organs <- organ_tier_full_df %>%
  dplyr::filter(enters_rho_tumor) %>%
  dplyr::pull(primary_organ)

included_delta_organs <- organ_tier_full_df %>%
  dplyr::filter(enters_rho_normal) %>%
  dplyr::pull(primary_organ)

if (length(included_tumor_organs) < 5) {
  stop("肿瘤端纳入分析的 primary_organ 少于 5 个: ", length(included_tumor_organs))
}
if (length(included_delta_organs) < 5) {
  stop("Δρ 纳入分析(肿瘤+正常共有)的 primary_organ 少于 5 个: ", length(included_delta_organs))
}

log_message(paste0(
  "肿瘤端纳入 organ 数: ", length(included_tumor_organs),
  "; Δρ 纳入 organ 数 (Tier1+1b+2): ", length(included_delta_organs),
  "; Tier3 (无可用正常) organ 数: ", length(included_tumor_organs) - length(included_delta_organs)
))

save(organ_tier_full_df, file = organ_tier_rdata)


# 14. organ-median 表达矩阵 ####

## 14.1 肿瘤端 organ-median ####
log_message("计算肿瘤 organ-median TPM 与 organ-median log2(TPM+1) 矩阵 (全 organ 集)...")

tumor_meta_filtered_df <- tumor_meta_df %>%
  dplyr::filter(primary_organ %in% included_tumor_organs) %>%
  dplyr::arrange(
    factor(pattern_group, levels = pattern_group_order),
    primary_organ, project_id, patient_id, sample_col_id
  )

tumor_tpm_matrix <- tumor_tpm_matrix[, tumor_meta_filtered_df$sample_col_id, drop = FALSE]

tumor_organ_median_result <- compute_tumor_organ_median_matrices(
  tumor_tpm_matrix = tumor_tpm_matrix,
  tumor_meta_df = tumor_meta_filtered_df,
  included_primary_organs = included_tumor_organs,
  log2_offset = log2_offset
)

organ_median_tpm_tumor_matrix <- tumor_organ_median_result$organ_median_tpm_matrix
organ_median_log2_tumor_matrix <- tumor_organ_median_result$organ_median_log2_matrix

log_message(paste0(
  "肿瘤 organ-median 矩阵维度: ",
  nrow(organ_median_log2_tumor_matrix), " 基因 × ",
  ncol(organ_median_log2_tumor_matrix), " primary_organ"
))

## 14.2 正常端 organ-median (按 Tier 选择数据源) ####
log_message("按 Tier 规则计算正常 organ-median log2 矩阵 (Tier1=TCGA+GTEx, Tier1b/2=GTEx-only)...")

normal_organ_median_result <- compute_normal_organ_median_matrix(
  normal_log2_matrix = normal_log2_matrix,
  normal_meta_df = normal_meta_df,
  organ_tier_df = organ_tier_full_df,
  included_primary_organs = included_tumor_organs
)
organ_median_log2_normal_matrix <- normal_organ_median_result$organ_median_log2_normal_matrix
normal_organ_sample_summary_df <- normal_organ_median_result$normal_organ_sample_summary_df

log_message(paste0(
  "正常 organ-median 矩阵维度: ",
  nrow(organ_median_log2_normal_matrix), " 基因 × ",
  ncol(organ_median_log2_normal_matrix), " primary_organ"
))

## 14.3 落盘 organ-median 结果 ####
save(organ_median_log2_tumor_matrix, organ_median_tpm_tumor_matrix,
     file = organ_median_log2_tumor_rdata)
save(organ_median_log2_normal_matrix, normal_organ_sample_summary_df,
     file = organ_median_log2_normal_rdata)


# 15. 基因预过滤 ####

log_message("应用基因预过滤(低表达 + 低方差 + 正常端覆盖；不做 tissue-of-origin 预过滤，由 Δρ 框架本身消除组织来源混杂)...")

prefilter_result <- apply_gene_prefilter(
  organ_median_log2_tumor_matrix = organ_median_log2_tumor_matrix,
  organ_median_log2_normal_matrix = organ_median_log2_normal_matrix,
  min_max_log2 = min_organ_median_log2,
  min_mad_log2 = min_organ_median_mad
)

retained_genes <- prefilter_result$retained_genes
gene_filter_summary_df <- prefilter_result$filter_summary_df
gene_filter_step_count_df <- prefilter_result$filter_step_count_df

log_message("基因预过滤每步统计:")
for (i in seq_len(nrow(gene_filter_step_count_df))) {
  log_message(paste0(
    "  ", gene_filter_step_count_df$步骤[[i]], " : ",
    gene_filter_step_count_df$基因数[[i]]
  ))
}

if (length(retained_genes) < 100) {
  stop("预过滤后保留基因数 < 100,请检查输入或调整阈值。当前: ", length(retained_genes))
}

organ_median_log2_tumor_filtered_matrix <- organ_median_log2_tumor_matrix[retained_genes, , drop = FALSE]
organ_median_log2_normal_filtered_matrix <- organ_median_log2_normal_matrix[retained_genes, , drop = FALSE]

rm(tumor_tpm_matrix, organ_median_tpm_tumor_matrix,
   organ_median_log2_tumor_matrix, organ_median_log2_normal_matrix,
   normal_log2_matrix)
gc()

save(gene_filter_summary_df, gene_filter_step_count_df,
     file = gene_filter_summary_rdata)


# 16. 观察相关 ρ_tumor_matched / ρ_normal_matched 与 Δρ ####

## 16.1 ρ_tumor: matched 子集 (Tier1+1b+2 与 H 矩阵共有 organ) ####
log_message("计算肿瘤端观察 Spearman 相关 (Δρ matched 子集)...")
organ_median_log2_tumor_matched_matrix <- organ_median_log2_tumor_filtered_matrix[, included_delta_organs, drop = FALSE]
h_matched_matrix <- h_organ_pattern_matrix[included_delta_organs, , drop = FALSE]
rho_tumor_matched_df <- compute_observed_correlation_long(
  organ_median_log2_matrix = organ_median_log2_tumor_matched_matrix,
  h_organ_pattern_matrix = h_matched_matrix,
  side_label = "tumor_matched"
)

## 16.2 ρ_normal: 同 matched 子集 ####
log_message("计算正常端观察 Spearman 相关 (Δρ matched 子集)...")
organ_median_log2_normal_matched_matrix <- organ_median_log2_normal_filtered_matrix[, included_delta_organs, drop = FALSE]
rho_normal_df <- compute_observed_correlation_long(
  organ_median_log2_matrix = organ_median_log2_normal_matched_matrix,
  h_organ_pattern_matrix = h_matched_matrix,
  side_label = "normal_matched"
)

save(rho_tumor_matched_df, file = rho_tumor_observed_rdata)
save(rho_normal_df, file = rho_normal_observed_rdata)

## 16.3 Δρ + Fisher r-to-z ####
log_message("计算 Δρ = ρ_tumor_matched − ρ_normal 与 Fisher r-to-z 显著性...")
delta_rho_df <- compute_delta_rho_with_fisher_z(
  rho_tumor_df = rho_tumor_matched_df,
  rho_normal_df = rho_normal_df,
  fdr_method = "BH"
)

save(delta_rho_df, file = delta_rho_df_rdata)


# 17. Permutation 灵敏性分析 (matched 子集 Δρ 零分布,正交报告) ####

if (file.exists(permutation_null_array_rdata)) {
  log_message(paste0(
    "检测到已有 permutation 结果文件,直接读取: ", permutation_null_array_rdata
  ))
  permutation_objects <- load_named_objects_from_rdata(
    rdata_file = permutation_null_array_rdata,
    object_names = c("permutation_tumor_matched_array",
                     "permutation_normal_matched_array")
  )
  permutation_tumor_matched_array <- permutation_objects$permutation_tumor_matched_array
  permutation_normal_matched_array <- permutation_objects$permutation_normal_matched_array
} else {
  log_message("跑 matched 子集双侧 permutation (用于 Δρ orthogonal FDR,同种子打散保证肿瘤/正常匹配)...")
  ## 注意:matched 双侧使用相同的随机种子,使每次迭代的 H 打散映射在两侧完全一致
  matched_perm_seed <- analysis_seed + 100L
  permutation_tumor_matched_array <- run_permutation_correlations(
    organ_median_log2_matrix = organ_median_log2_tumor_matched_matrix,
    h_organ_pattern_matrix = h_matched_matrix,
    n_iterations = n_permutation_iterations,
    seed = matched_perm_seed,
    use_parallel = enable_parallel,
    workers = parallel_workers
  )
  permutation_normal_matched_array <- run_permutation_correlations(
    organ_median_log2_matrix = organ_median_log2_normal_filtered_matrix[, included_delta_organs, drop = FALSE],
    h_organ_pattern_matrix = h_matched_matrix,
    n_iterations = n_permutation_iterations,
    seed = matched_perm_seed,
    use_parallel = enable_parallel,
    workers = parallel_workers
  )

  save(permutation_tumor_matched_array,
       permutation_normal_matched_array,
       file = permutation_null_array_rdata)
  log_message(paste0("permutation 结果已保存: ", permutation_null_array_rdata))
}

permutation_pvalue_delta_df <- compute_delta_rho_permutation_pvalue(
  observed_delta_df = delta_rho_df,
  tumor_perm_array = permutation_tumor_matched_array,
  normal_perm_array = permutation_normal_matched_array,
  fdr_method = "BH"
)


# 18. 合并指标 + 纯 Δρ 多阈值筛选 + 多模式分类 ####

log_message("合并 ρ_tumor_matched / ρ_normal / Δρ / Fisher z / Δρ permutation 经验 p 与 FDR...")

merged_correlation_df <- delta_rho_df %>%
  dplyr::left_join(
    permutation_pvalue_delta_df,
    by = c("gene", "pattern")
  ) %>%
  dplyr::arrange(gene, pattern)

log_message(paste0("合并后行数: ", nrow(merged_correlation_df), " (基因 × 模式组合)"))

log_message(paste0(
  "应用 Δρ 多阈值筛选: |Δρ| ∈ {",
  paste(c(delta_rho_threshold_main, delta_rho_threshold_supplementary), collapse = ", "),
  "}; 主检验 perm_p < ", perm_p_threshold,
  "; 高置信子集追加 fisher_p < ", fisher_p_threshold_high_confidence,
  "; (BH FDR 仅作透明性输出, 不参与 hit 判定)"
))

threshold_filter_df <- apply_delta_rho_filter(
  merged_long_df = merged_correlation_df,
  delta_rho_thresholds = c(delta_rho_threshold_main, delta_rho_threshold_supplementary),
  perm_p_threshold = perm_p_threshold,
  fisher_p_threshold_high_confidence = fisher_p_threshold_high_confidence
)

hit_lists_result <- build_per_pattern_hit_lists(threshold_filter_df)
hit_long_df <- hit_lists_result$hit_long_df
hit_long_high_confidence_df <- hit_lists_result$hit_long_high_confidence_df
hit_count_df <- hit_lists_result$hit_count_df

log_message("主阈值下每模式 hit 计数 (主筛选 / high-confidence):")
hit_count_main_df <- hit_count_df %>%
  dplyr::filter(delta_rho_threshold == delta_rho_threshold_main)
for (i in seq_len(nrow(hit_count_main_df))) {
  log_message(paste0(
    "  ", as.character(hit_count_main_df$pattern[[i]]),
    " : positive=", hit_count_main_df$n_positive_hits[[i]],
    " (HC=", hit_count_main_df$n_positive_high_confidence[[i]], ")",
    "; negative=", hit_count_main_df$n_negative_hits[[i]],
    " (HC=", hit_count_main_df$n_negative_high_confidence[[i]], ")",
    "; total=", hit_count_main_df$n_total_hits[[i]],
    " (HC=", hit_count_main_df$n_total_high_confidence[[i]], ")"
  ))
}

multi_pattern_classification_df <- classify_multi_pattern_genes(
  threshold_filter_df = threshold_filter_df,
  delta_rho_threshold_target = delta_rho_threshold_main
)

save(merged_correlation_df, threshold_filter_df,
     hit_long_df, hit_long_high_confidence_df, hit_count_df,
     multi_pattern_classification_df,
     file = final_pattern_gene_list_rdata)


# 19. Excel 输出(中文 sheet 命名) ####

log_message("整理输出 Excel 工作簿...")

analysis_parameter_df <- data.frame(
  参数 = c(
    "log2_offset",
    "min_samples_per_organ_tumor",
    "tcga_normal_min_for_tier1",
    "gtex_min_for_normal_baseline",
    "min_organ_median_log2",
    "min_organ_median_mad",
    "n_permutation_iterations",
    "delta_rho_threshold_main",
    "delta_rho_threshold_supplementary",
    "perm_p_threshold (主检验, matched-pair permutation 双尾经验 p)",
    "fisher_p_threshold_high_confidence (高置信子集追加, Fisher r-to-z 平行检验)",
    "analysis_seed",
    "matched_permutation_seed (Δρ, 双侧同种子)"
  ),
  值 = c(
    log2_offset,
    min_samples_per_organ_tumor,
    tcga_normal_min_for_tier1,
    gtex_min_for_normal_baseline,
    min_organ_median_log2,
    min_organ_median_mad,
    n_permutation_iterations,
    delta_rho_threshold_main,
    paste(delta_rho_threshold_supplementary, collapse = ","),
    perm_p_threshold,
    fisher_p_threshold_high_confidence,
    analysis_seed,
    analysis_seed + 100L
  ),
  stringsAsFactors = FALSE
)

organ_tier_export_df <- organ_tier_full_df %>%
  dplyr::transmute(
    primary_organ,
    n_tumor_samples, n_tcga_normal, n_gtex_normal,
    tier, normal_baseline_source,
    in_h_matrix,
    enters_rho_tumor, enters_rho_normal
  )

merged_correlation_export_df <- merged_correlation_df %>%
  dplyr::mutate(pattern = as.character(pattern))

excel_sheets_list <- list(
  "依赖检查" = dependency_check_df,
  "分析参数" = analysis_parameter_df,
  "器官纳入与样本量分层" = organ_tier_export_df,
  "正常组织样本使用情况" = normal_organ_sample_summary_df,
  "基因预过滤每步统计" = gene_filter_step_count_df,
  "基因预过滤明细" = gene_filter_summary_df,
  "观察相关_肿瘤端配对子集" = rho_tumor_matched_df %>% dplyr::mutate(pattern = as.character(pattern)),
  "观察相关_正常端配对子集" = rho_normal_df %>% dplyr::mutate(pattern = as.character(pattern)),
  "差分相关全表_费舍尔变换" = delta_rho_df %>% dplyr::mutate(pattern = as.character(pattern)),
  "差分相关全表_置换检验" = permutation_pvalue_delta_df %>% dplyr::mutate(pattern = as.character(pattern)),
  "合并指标全表" = merged_correlation_export_df,
  "主筛选命中列表_主阈值" = hit_long_df %>%
    dplyr::filter(delta_rho_threshold == delta_rho_threshold_main),
  "高置信命中列表_主阈值" = hit_long_high_confidence_df %>%
    dplyr::filter(delta_rho_threshold == delta_rho_threshold_main),
  "命中计数_主与高置信_全阈值" = hit_count_df %>% dplyr::mutate(pattern = as.character(pattern)),
  "多模式基因分类_主阈值" = multi_pattern_classification_df
)

openxlsx::write.xlsx(
  x = excel_sheets_list,
  file = output_main_xlsx,
  overwrite = TRUE
)

log_message(paste0("已导出主结果 Excel: ", output_main_xlsx))


# 20. 完成 ####

log_message("========== 7.7 模式相关基因差分相关识别完成 ==========")
log_message(paste0(
  "主结果文件: ", output_main_xlsx,
  "; 中间 Rdata 路径: ", rdata_output_path
))
