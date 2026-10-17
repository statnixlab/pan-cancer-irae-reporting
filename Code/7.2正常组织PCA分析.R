# 7.2 正常组织PCA分析 ####
# 用途：对TCGA正常组织和GTEx正常组织的统一表达矩阵进行PCA分析，
#       评估两个平台正常组织样本在基因表达层面的异同。
#       第一阶段：校正前PCA（若已完成则跳过）
#       第二阶段：ComBat批次校正后PCA
#       纳入样本：TCGA Solid Tissue Normal + GTEx正常组织 ####

rm(list = ls())

# 1. 加载配置 ####
source("Code/config.R")


# 2. 参数与路径 ####
script_name <- "7.2正常组织PCA分析"
script_results_path <- file.path(results_path, script_name)
table_output_path <- file.path(script_results_path, "tables")
rdata_output_path <- file.path(script_results_path, "rdata")
figure_output_path <- file.path(script_results_path, "figures")
log_output_path <- file.path(script_results_path, "logs")

dir.create(script_results_path, recursive = TRUE, showWarnings = FALSE)
dir.create(table_output_path, recursive = TRUE, showWarnings = FALSE)
dir.create(rdata_output_path, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_output_path, recursive = TRUE, showWarnings = FALSE)
dir.create(log_output_path, recursive = TRUE, showWarnings = FALSE)

## 7.1输出路径（输入数据） ####
input_data_path <- file.path(results_path, "7.1TCGA基础数据准备")
input_rdata_path <- file.path(input_data_path, "rdata")
input_unified_path <- file.path(input_rdata_path, "unified")

unified_gtex_expression_rdata <- file.path(input_unified_path, "unified_gtex_normal_expression_matrix.Rdata")
unified_tcga_normal_expression_rdata <- file.path(input_unified_path, "unified_tcga_normal_expression_matrix.Rdata")
unified_gene_annotation_rdata <- file.path(input_unified_path, "unified_gene_annotation.Rdata")
gtex_sample_annotation_rdata <- file.path(input_rdata_path, "gtex_normal_sample_annotation.Rdata")
tcga_normal_sample_annotation_rdata <- file.path(input_rdata_path, "tcga_normal_sample_annotation.Rdata")

## PCA参数 ####
n_top_var_genes <- 5000  # 用于PCA的高变异基因数
log2_offset <- 1         # log2(TPM + offset)

## 日志 ####
run_timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_file <- file.path(log_output_path, paste0("7.2正常组织PCA分析_", run_timestamp, ".log"))


# 3. 日志与辅助函数 ####

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

log_message("========== 7.2 正常组织PCA分析开始 ==========")


# 4. 检查输入文件 ####
required_files <- c(
  unified_gtex_expression_rdata,
  unified_tcga_normal_expression_rdata,
  unified_gene_annotation_rdata,
  gtex_sample_annotation_rdata,
  tcga_normal_sample_annotation_rdata
)

required_labels <- c(
  "统一GTEx正常表达矩阵",
  "统一TCGA正常表达矩阵",
  "统一基因注释",
  "GTEx样本注释",
  "TCGA正常样本注释"
)

missing_idx <- which(!file.exists(required_files))
if (length(missing_idx) > 0) {
  stop("缺少必要输入文件:\n", paste(required_labels[missing_idx], ": ", required_files[missing_idx], collapse = "\n"))
}


# 5. 构建合并样本注释（两个阶段共用） ####
log_message("读取样本注释...")

gtex_sample_annotation_df <- load_single_object_from_rdata(
  rdata_file = gtex_sample_annotation_rdata,
  object_name = "gtex_final_sample_annotation_df"
)

tcga_normal_sample_annotation_df <- load_single_object_from_rdata(
  rdata_file = tcga_normal_sample_annotation_rdata,
  object_name = "tcga_normal_sample_annotation_df"
)

gtex_meta_df <- data.frame(
  sample_col_id = gtex_sample_annotation_df$sample_id,
  data_source = "GTEx",
  primary_organ = gtex_sample_annotation_df$primary_organ,
  pattern_group = gtex_sample_annotation_df$pattern_group,
  stringsAsFactors = FALSE
)

tcga_normal_meta_df <- data.frame(
  sample_col_id = tcga_normal_sample_annotation_df$sample_barcode,
  data_source = "TCGA",
  primary_organ = tcga_normal_sample_annotation_df$primary_organ,
  pattern_group = tcga_normal_sample_annotation_df$pattern_group,
  stringsAsFactors = FALSE
)

combined_meta_df <- rbind(gtex_meta_df, tcga_normal_meta_df)
log_message(paste0("合并样本注释总行数: ", nrow(combined_meta_df)))
log_message(paste0("数据来源分布: GTEx ", sum(combined_meta_df$data_source == "GTEx"),
                   " | TCGA ", sum(combined_meta_df$data_source == "TCGA")))


# ============================================================
# 第一阶段：校正前PCA（若结果已存在则跳过）
# ============================================================

# 6. 校正前PCA ####
pca_result_xlsx_path <- file.path(table_output_path, "PCA分析结果.xlsx")

if (file.exists(pca_result_xlsx_path)) {
  log_message("校正前PCA结果已存在，跳过第一阶段。")
} else {
  log_message("===== 开始校正前PCA =====")

  log_message(paste0("PCA高变异基因数: ", n_top_var_genes))
  log_message(paste0("log2变换偏移量: ", log2_offset))

  ## 读取表达矩阵 ####
  log_message("读取统一GTEx正常表达矩阵...")
  unified_gtex_normal_expression_matrix <- load_single_object_from_rdata(
    rdata_file = unified_gtex_expression_rdata,
    object_name = "unified_gtex_normal_expression_matrix"
  )
  log_message(paste0("GTEx矩阵维度: ", nrow(unified_gtex_normal_expression_matrix), " × ", ncol(unified_gtex_normal_expression_matrix)))

  log_message("读取统一TCGA正常表达矩阵...")
  unified_tcga_normal_expression_matrix <- load_single_object_from_rdata(
    rdata_file = unified_tcga_normal_expression_rdata,
    object_name = "unified_tcga_normal_expression_matrix"
  )
  log_message(paste0("TCGA正常矩阵维度: ", nrow(unified_tcga_normal_expression_matrix), " × ", ncol(unified_tcga_normal_expression_matrix)))

  ## 合并表达矩阵 ####
  log_message("合并GTEx与TCGA正常表达矩阵...")

  stopifnot(identical(rownames(unified_gtex_normal_expression_matrix),
                      rownames(unified_tcga_normal_expression_matrix)))

  combined_expression_matrix <- cbind(
    unified_gtex_normal_expression_matrix,
    unified_tcga_normal_expression_matrix
  )

  log_message(paste0("合并矩阵维度: ", nrow(combined_expression_matrix), " × ", ncol(combined_expression_matrix)))

  stopifnot(all(combined_meta_df$sample_col_id %in% colnames(combined_expression_matrix)))
  combined_expression_matrix <- combined_expression_matrix[, combined_meta_df$sample_col_id, drop = FALSE]

  ## log2变换 ####
  log_message("log2变换 (log2(TPM + 1))...")
  combined_log2_matrix <- log2(combined_expression_matrix + log2_offset)

  ## 计算基因方差，选取高变异基因 ####
  log_message("计算基因方差...")
  gene_variances <- apply(combined_log2_matrix, 1, var)

  valid_gene_mask <- !is.na(gene_variances) & gene_variances > 0
  n_valid_genes <- sum(valid_gene_mask)
  log_message(paste0("有效基因数（方差>0）: ", n_valid_genes))

  n_select <- min(n_top_var_genes, n_valid_genes)
  top_var_gene_names <- names(sort(gene_variances[valid_gene_mask], decreasing = TRUE))[1:n_select]
  log_message(paste0("选取高变异基因数: ", n_select))

  combined_log2_matrix_filtered <- combined_log2_matrix[top_var_gene_names, , drop = FALSE]

  ## PCA分析 ####
  log_message("执行PCA分析...")

  pca_input_matrix <- t(combined_log2_matrix_filtered)

  pca_result <- prcomp(pca_input_matrix, center = TRUE, scale. = TRUE)

  log_message("PCA分析完成。")

  ## 整理PCA结果 ####
  pca_variance_explained <- summary(pca_result)$importance
  variance_df <- data.frame(
    PC = paste0("PC", 1:ncol(pca_variance_explained)),
    standard_deviation = pca_variance_explained["Standard deviation", ],
    proportion_of_variance = pca_variance_explained["Proportion of Variance", ],
    cumulative_proportion = pca_variance_explained["Cumulative Proportion", ],
    stringsAsFactors = FALSE
  )

  n_pc_export <- min(20, nrow(variance_df))
  variance_export_df <- variance_df[1:n_pc_export, ]

  log_message(paste0("PC1方差解释比例: ", round(variance_df$proportion_of_variance[1] * 100, 2), "%"))
  log_message(paste0("PC2方差解释比例: ", round(variance_df$proportion_of_variance[2] * 100, 2), "%"))
  log_message(paste0("前10个PC累计方差解释: ", round(variance_df$cumulative_proportion[min(10, nrow(variance_df))] * 100, 2), "%"))

  n_pc_save <- min(20, ncol(pca_result$x))
  pca_coords_df <- as.data.frame(pca_result$x[, 1:n_pc_save, drop = FALSE])
  pca_coords_df$sample_col_id <- rownames(pca_coords_df)

  pca_plot_df <- merge(pca_coords_df, combined_meta_df, by = "sample_col_id", all.x = TRUE)
  log_message(paste0("PCA结果表行数: ", nrow(pca_plot_df)))

  ## 保存PCA分析结果 ####
  log_message("保存PCA分析结果...")

  openxlsx::write.xlsx(
    x = list(
      "PCA坐标与注释" = pca_plot_df,
      "方差解释比例" = variance_export_df
    ),
    file = file.path(table_output_path, "PCA分析结果.xlsx"),
    overwrite = TRUE
  )
  log_message(paste0("已导出PCA分析结果: ", file.path(table_output_path, "PCA分析结果.xlsx")))

  save(pca_result, file = file.path(rdata_output_path, "pca_result.Rdata"))
  save(pca_plot_df, file = file.path(rdata_output_path, "pca_plot_df.Rdata"))
  save(variance_df, file = file.path(rdata_output_path, "pca_variance_df.Rdata"))

  log_message(paste0("已保存PCA结果Rdata: ", file.path(rdata_output_path, "pca_result.Rdata")))

  organ_source_summary_df <- pca_plot_df %>%
    dplyr::group_by(primary_organ, data_source) %>%
    dplyr::summarise(
      n_samples = dplyr::n(),
      PC1_mean = mean(PC1, na.rm = TRUE),
      PC1_sd = sd(PC1, na.rm = TRUE),
      PC2_mean = mean(PC2, na.rm = TRUE),
      PC2_sd = sd(PC2, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(primary_organ, data_source)

  openxlsx::write.xlsx(
    x = list(
      "逐器官数据来源汇总" = organ_source_summary_df
    ),
    file = file.path(table_output_path, "PCA逐器官汇总.xlsx"),
    overwrite = TRUE
  )
  log_message(paste0("已导出PCA逐器官汇总: ", file.path(table_output_path, "PCA逐器官汇总.xlsx")))

  log_message("===== 校正前PCA分析完成 =====")
}  # end of skip logic


# ============================================================
# 第二阶段：批次校正后PCA分析
# 使用ComBat校正TCGA与GTEx之间的平台批次效应，
# 然后对校正后矩阵重新执行PCA。
# ============================================================


# 7. 加载sva包 ####
if (!requireNamespace("sva", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", quiet = TRUE)
  }
  BiocManager::install("sva", ask = FALSE, update = FALSE)
}
library(sva)

log_message("========== 批次校正后PCA分析开始 ==========")


# 8. 读取合并矩阵（用于批次校正） ####
log_message("读取表达矩阵用于批次校正...")

unified_gtex_normal_expression_matrix <- load_single_object_from_rdata(
  rdata_file = unified_gtex_expression_rdata,
  object_name = "unified_gtex_normal_expression_matrix"
)

unified_tcga_normal_expression_matrix <- load_single_object_from_rdata(
  rdata_file = unified_tcga_normal_expression_rdata,
  object_name = "unified_tcga_normal_expression_matrix"
)

stopifnot(identical(rownames(unified_gtex_normal_expression_matrix),
                    rownames(unified_tcga_normal_expression_matrix)))

combined_expression_matrix <- cbind(
  unified_gtex_normal_expression_matrix,
  unified_tcga_normal_expression_matrix
)

combined_expression_matrix <- combined_expression_matrix[, combined_meta_df$sample_col_id, drop = FALSE]

log_message(paste0("合并矩阵维度: ", nrow(combined_expression_matrix), " × ", ncol(combined_expression_matrix)))


# 9. log2变换 ####
log_message("log2变换 (log2(TPM + 1))...")
combined_log2_matrix <- log2(combined_expression_matrix + log2_offset)


# 10. ComBat批次校正 ####
log_message("执行ComBat批次校正（batch = data_source）...")

## 构建batch向量（与矩阵列顺序一致） ####
batch_vector <- combined_meta_df$data_source
names(batch_vector) <- combined_meta_df$sample_col_id
batch_vector <- batch_vector[colnames(combined_log2_matrix)]

## 构建协变量模型矩阵（保留primary_organ生物学差异） ####
mod <- model.matrix(~ primary_organ, data = combined_meta_df[match(colnames(combined_log2_matrix), combined_meta_df$sample_col_id), ])

## 过滤方差为0的基因（ComBat要求） ####
gene_vars_pre <- apply(combined_log2_matrix, 1, var)
nonzero_var_mask <- !is.na(gene_vars_pre) & gene_vars_pre > 0
log_message(paste0("ComBat前过滤零方差基因: 保留 ", sum(nonzero_var_mask), " / ", length(nonzero_var_mask)))

combined_log2_matrix_for_combat <- combined_log2_matrix[nonzero_var_mask, , drop = FALSE]

## 执行ComBat ####
corrected_log2_matrix <- ComBat(
  dat = combined_log2_matrix_for_combat,
  batch = batch_vector,
  mod = mod,
  par.prior = TRUE,
  prior.plots = FALSE
)

log_message(paste0("ComBat校正完成，校正后矩阵维度: ",
                   nrow(corrected_log2_matrix), " × ", ncol(corrected_log2_matrix)))

## 保存校正后矩阵 ####
save(corrected_log2_matrix, file = file.path(rdata_output_path, "corrected_log2_matrix.Rdata"))
log_message(paste0("已保存校正后矩阵: ", file.path(rdata_output_path, "corrected_log2_matrix.Rdata")))


# 11. 校正后选取高变异基因 + PCA ####
log_message("校正后计算基因方差...")

gene_variances_corrected <- apply(corrected_log2_matrix, 1, var)
valid_gene_mask_corrected <- !is.na(gene_variances_corrected) & gene_variances_corrected > 0
n_valid_genes_corrected <- sum(valid_gene_mask_corrected)
log_message(paste0("校正后有效基因数（方差>0）: ", n_valid_genes_corrected))

n_select_corrected <- min(n_top_var_genes, n_valid_genes_corrected)
top_var_gene_names_corrected <- names(sort(gene_variances_corrected[valid_gene_mask_corrected], decreasing = TRUE))[1:n_select_corrected]
log_message(paste0("校正后选取高变异基因数: ", n_select_corrected))

corrected_log2_matrix_filtered <- corrected_log2_matrix[top_var_gene_names_corrected, , drop = FALSE]

log_message("执行校正后PCA分析...")

pca_input_corrected <- t(corrected_log2_matrix_filtered)

pca_result_corrected <- prcomp(pca_input_corrected, center = TRUE, scale. = TRUE)

log_message("校正后PCA分析完成。")


# 12. 整理校正后PCA结果 ####

## 方差解释比例 ####
pca_var_corrected <- summary(pca_result_corrected)$importance
variance_corrected_df <- data.frame(
  PC = paste0("PC", 1:ncol(pca_var_corrected)),
  standard_deviation = pca_var_corrected["Standard deviation", ],
  proportion_of_variance = pca_var_corrected["Proportion of Variance", ],
  cumulative_proportion = pca_var_corrected["Cumulative Proportion", ],
  stringsAsFactors = FALSE
)

n_pc_export_c <- min(20, nrow(variance_corrected_df))
variance_corrected_export_df <- variance_corrected_df[1:n_pc_export_c, ]

log_message(paste0("校正后PC1方差解释比例: ", round(variance_corrected_df$proportion_of_variance[1] * 100, 2), "%"))
log_message(paste0("校正后PC2方差解释比例: ", round(variance_corrected_df$proportion_of_variance[2] * 100, 2), "%"))
log_message(paste0("校正后前10个PC累计方差解释: ", round(variance_corrected_df$cumulative_proportion[min(10, nrow(variance_corrected_df))] * 100, 2), "%"))

## PCA坐标 + 样本注释 ####
n_pc_save_c <- min(20, ncol(pca_result_corrected$x))
pca_coords_corrected_df <- as.data.frame(pca_result_corrected$x[, 1:n_pc_save_c, drop = FALSE])
pca_coords_corrected_df$sample_col_id <- rownames(pca_coords_corrected_df)

pca_plot_corrected_df <- merge(pca_coords_corrected_df, combined_meta_df, by = "sample_col_id", all.x = TRUE)

log_message(paste0("校正后PCA结果表行数: ", nrow(pca_plot_corrected_df)))


# 13. 保存校正后PCA结果 ####
log_message("保存校正后PCA分析结果...")

## 保存Excel ####
openxlsx::write.xlsx(
  x = list(
    "校正后PCA坐标与注释" = pca_plot_corrected_df,
    "校正后方差解释比例" = variance_corrected_export_df
  ),
  file = file.path(table_output_path, "PCA分析结果_批次校正后.xlsx"),
  overwrite = TRUE
)
log_message(paste0("已导出: ", file.path(table_output_path, "PCA分析结果_批次校正后.xlsx")))

## 保存Rdata ####
save(pca_result_corrected, file = file.path(rdata_output_path, "pca_result_corrected.Rdata"))
save(pca_plot_corrected_df, file = file.path(rdata_output_path, "pca_plot_corrected_df.Rdata"))
save(variance_corrected_df, file = file.path(rdata_output_path, "pca_variance_corrected_df.Rdata"))

log_message(paste0("已保存校正后PCA Rdata: ", file.path(rdata_output_path, "pca_result_corrected.Rdata")))

## 保存逐organ统计汇总 ####
organ_source_summary_corrected_df <- pca_plot_corrected_df %>%
  dplyr::group_by(primary_organ, data_source) %>%
  dplyr::summarise(
    n_samples = dplyr::n(),
    PC1_mean = mean(PC1, na.rm = TRUE),
    PC1_sd = sd(PC1, na.rm = TRUE),
    PC2_mean = mean(PC2, na.rm = TRUE),
    PC2_sd = sd(PC2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(primary_organ, data_source)

openxlsx::write.xlsx(
  x = list(
    "校正后逐器官数据来源汇总" = organ_source_summary_corrected_df
  ),
  file = file.path(table_output_path, "PCA逐器官汇总_批次校正后.xlsx"),
  overwrite = TRUE
)
log_message(paste0("已导出: ", file.path(table_output_path, "PCA逐器官汇总_批次校正后.xlsx")))

log_message("========== 7.2 正常组织PCA分析全部完成 ==========")
