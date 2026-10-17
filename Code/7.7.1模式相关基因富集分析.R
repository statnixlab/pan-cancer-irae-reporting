# 7.7.1 模式相关基因富集分析 ####
# 用途：读取 7.7 输出的差分 Spearman 相关分析结果，对每个 NMF irAE 模式 ×
#       (positive / negative) 方向 × (候选 hit 主轨 / 高置信子集 HC 辅轨) 组合，
#       在 GO BP / KEGG / Reactome / Hallmark / C7 ImmuneSigDB 五个 MSigDB 子集
#       上分别执行：
#         (1) ORA 主算法 (clusterProfiler::enricher, 超几何分布检验)
#         (2) GSEA 补充算法 (clusterProfiler::GSEA, Δρ 有符号排序)
#         (3) 生物学 positive control (每模式预定义 expected pathway sets) ####
# 备注：背景为 7.7 预过滤后保留的基因 (retained = TRUE)，与 hit 同一漏斗。 ####

rm(list = ls())
gc()


# 1. 加载配置 ####
source("Code/config.R")


# 2. 参数与路径 ####
script_name <- "7.7.1模式相关基因富集分析"
script_results_path <- file.path(results_path, script_name)
table_output_path <- file.path(script_results_path, "tables")
rdata_output_path <- file.path(script_results_path, "rdata")
log_output_path <- file.path(script_results_path, "logs")

dir.create(script_results_path, recursive = TRUE, showWarnings = FALSE)
dir.create(table_output_path, recursive = TRUE, showWarnings = FALSE)
dir.create(rdata_output_path, recursive = TRUE, showWarnings = FALSE)
dir.create(log_output_path, recursive = TRUE, showWarnings = FALSE)

## 输入文件路径 (来自 7.7 输出) ####
input_77_path <- file.path(results_path, "7.7模式相关基因差分相关识别")
input_77_rdata_path <- file.path(input_77_path, "rdata")

final_pattern_gene_input_rdata <- file.path(
  input_77_rdata_path,
  "final_pattern_gene_list.Rdata"
)
gene_filter_summary_input_rdata <- file.path(
  input_77_rdata_path,
  "gene_filter_summary_df.Rdata"
)

## 输出文件路径 ####
output_xlsx <- file.path(table_output_path, "模式相关基因富集分析结果.xlsx")
dependency_xlsx <- file.path(table_output_path, "依赖检查.xlsx")

database_term2gene_rdata <- file.path(rdata_output_path, "enrichment_database_term2gene_list.Rdata")
ora_long_df_rdata <- file.path(rdata_output_path, "ora_long_df_all.Rdata")
gsea_long_df_rdata <- file.path(rdata_output_path, "gsea_long_df_all.Rdata")
positive_control_rdata <- file.path(rdata_output_path, "positive_control_df.Rdata")
cross_pattern_shared_rdata <- file.path(rdata_output_path, "cross_pattern_shared_df.Rdata")
cross_pattern_unique_rdata <- file.path(rdata_output_path, "cross_pattern_unique_df.Rdata")
top_pathways_rdata <- file.path(rdata_output_path, "top_pathways_per_pattern_df.Rdata")
ora_gsea_agreement_rdata <- file.path(rdata_output_path, "ora_gsea_agreement_df.Rdata")

## 主参数 ####
msigdb_species <- "Homo sapiens"

## MSigDB 数据库配置 (v26+ 用 gs_collection / gs_subcollection) ####
msigdb_database_config <- list(
  GO_BP = list(
    label = "GO Biological Process",
    collection = "C5",
    subcollection = "GO:BP"
  ),
  KEGG = list(
    label = "KEGG (Legacy + Medicus)",
    collection = "C2",
    subcollection = c("CP:KEGG_LEGACY", "CP:KEGG_MEDICUS")
  ),
  Reactome = list(
    label = "Reactome",
    collection = "C2",
    subcollection = "CP:REACTOME"
  ),
  Hallmark = list(
    label = "MSigDB Hallmark",
    collection = "H",
    subcollection = NULL
  ),
  C7_IMMUNESIGDB = list(
    label = "MSigDB C7 ImmuneSigDB",
    collection = "C7",
    subcollection = "IMMUNESIGDB"
  )
)

## 富集分析参数 ####
enrichment_min_gs_size <- 10
enrichment_max_gs_size <- 500
enrichment_pvalue_cutoff <- 1            # enricher 内部不做 cutoff，下游统一 FDR<0.05 筛选
enrichment_qvalue_cutoff <- 1
enrichment_fdr_threshold <- 0.05         # 显著阈值
top_n_pathway_per_pattern <- 20          # 每模式每数据库 Top N 通路报告数
min_n_patterns_for_shared <- 2           # 跨模式共享通路最低模式数
min_hit_genes_for_enrichment <- 5        # 单 ORA 任务最低 hit 数 (< 5 跳过)

## GSEA 参数 ####
gsea_n_perm <- 10000                     # GSEA permutation 次数 (fgsea simple)
gsea_eps <- 0                            # fgsea 精度 (0 = 自动)
gsea_seed <- 20251208                    # GSEA 随机种子 (与 7.7 风格一致, 复现性)

## 阈值与方向 (与 7.7 一致) ####
delta_rho_threshold_main <- 0.3
delta_rho_thresholds_evaluated <- c(0.2, 0.3, 0.4)
hit_direction_levels <- c("positive", "negative")
hit_track_levels <- c("candidate", "high_confidence")

## pattern 顺序 (与 6.1 / 7.4 / 7.6 / 7.7 一致) ####
pattern_group_order <- c(
  "Ocular",
  "Hematologic",
  "Neurologic",
  "Systemic-like",
  "Visceral-like"
)

## 生物学 positive control: 每模式预定义 expected pathway gene set name (MSigDB term ID 关键词) ####
## 说明: 基于先前 irAE 文献先验的"expected biology"，用于验证方法学
##       是否能恢复模式特异性生物学信号。每模式给出多个 keyword (substring 匹配)，
##       从 5 库 MSigDB 中筛选 term name 包含任一 keyword 的 gene set 作为该模式的 expected set。
positive_control_keywords <- list(
  "Ocular" = c(
    "VISUAL_PERCEPTION",
    "EYE_DEVELOPMENT",
    "RETINA",
    "PIGMENTATION",
    "MELANIN",
    "MELANOCYTE",
    "PHOTORECEPTOR"
  ),
  "Hematologic" = c(
    "HEMATOPOIESIS",
    "HEMATOPOIETIC",
    "ERYTHROCYTE",
    "PLATELET",
    "LEUKOCYTE",
    "MYELOID",
    "LYMPHOID",
    "HEME_METABOLISM"
  ),
  "Neurologic" = c(
    "SYNAPTIC",
    "SYNAPSE",
    "NEURON",
    "NEUROTRANSMITTER",
    "AXON",
    "DENDRITE",
    "MYELIN",
    "NERVOUS_SYSTEM"
  ),
  "Systemic-like" = c(
    "INFLAMMATORY_RESPONSE",
    "TNFA_SIGNALING",
    "INTERFERON",
    "INTERLEUKIN",
    "CYTOKINE",
    "COMPLEMENT",
    "INNATE_IMMUNE",
    "ADAPTIVE_IMMUNE"
  ),
  "Visceral-like" = c(
    "BILE_ACID",
    "HEPATOCYTE",
    "LIVER",
    "THYROID",
    "ADRENAL",
    "GASTRIC",
    "INTESTINAL",
    "PANCREAS",
    "ENDOCRINE"
  )
)

## 并行设置 ####
detected_cores <- parallel::detectCores(logical = FALSE)
if (is.na(detected_cores)) {
  detected_cores <- 1
}
parallel_workers <- as.integer(max(1, min(8, detected_cores - 1)))
enable_parallel <- parallel_workers > 1 && .Platform$OS.type != "windows"

required_packages <- c(
  "msigdbr",
  "clusterProfiler",
  "fgsea",
  "openxlsx",
  "writexl",
  "dplyr",
  "tidyr",
  "tibble",
  "purrr"
)

run_timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_file <- file.path(log_output_path, paste0(script_name, "_", run_timestamp, ".log"))


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
## 功能：检查必要输入文件是否存在 ####
## 输入：file_paths，文件路径向量；file_labels，文件说明向量 ####
## 输出：无返回值 ####
assert_required_files <- function(file_paths, file_labels) {
  missing_index <- which(!file.exists(file_paths))
  if (length(missing_index) > 0) {
    missing_messages <- paste0(file_labels[missing_index], ": ", file_paths[missing_index])
    stop("缺少必要输入文件:\n", paste(missing_messages, collapse = "\n"))
  }
}

## 函数：load_named_objects_from_rdata ####
## 功能：从 Rdata 文件中读取一个或多个对象 ####
## 输入：rdata_file，Rdata 文件路径；object_names，对象名向量 ####
## 输出：命名列表 ####
load_named_objects_from_rdata <- function(rdata_file, object_names) {
  temp_env <- new.env(parent = emptyenv())
  load(rdata_file, envir = temp_env)

  missing_objects <- setdiff(object_names, ls(envir = temp_env))
  if (length(missing_objects) > 0) {
    stop("Rdata 文件中不存在目标对象: ", paste(missing_objects, collapse = ", "), " | ", rdata_file)
  }

  result_list <- lapply(object_names, function(current_name) {
    get(current_name, envir = temp_env, inherits = FALSE)
  })
  names(result_list) <- object_names
  result_list
}

## 函数：check_required_packages ####
## 功能：检查本脚本所需 R 包是否安装并写出依赖检查 Excel ####
## 输入：packages，包名向量；output_path，依赖检查结果保存路径 ####
## 输出：依赖检查数据框 ####
check_required_packages <- function(packages, output_path) {
  dependency_check_df <- data.frame(
    package = packages,
    installed = vapply(packages, requireNamespace, logical(1), quietly = TRUE),
    install_hint = dplyr::case_when(
      packages %in% c("msigdbr", "openxlsx", "writexl", "dplyr", "tidyr", "tibble", "purrr") ~
        paste0("install.packages(\"", packages, "\")"),
      packages %in% c("clusterProfiler", "fgsea") ~
        paste0("BiocManager::install(\"", packages, "\")"),
      TRUE ~ paste0("install.packages(\"", packages, "\")")
    ),
    stringsAsFactors = FALSE
  )

  openxlsx::write.xlsx(
    x = list("依赖检查" = dependency_check_df),
    file = output_path,
    overwrite = TRUE
  )

  if (any(!dependency_check_df$installed)) {
    missing_df <- dependency_check_df[!dependency_check_df$installed, , drop = FALSE]
    missing_message <- paste0(
      missing_df$package, ": ", missing_df$install_hint,
      collapse = "\n"
    )
    stop("缺少 7.7.1 富集分析必要 R 包:\n", missing_message)
  }

  dependency_check_df
}

## 函数：run_parallel_lapply ####
## 功能：在 Linux/macOS 中用 mclapply 并行执行任务，无法并行时退回 lapply ####
## 输入：x，任务输入；fun，任务函数；use_parallel；workers；... ####
## 输出：列表结果 ####
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

## 函数：truncate_excel_string_columns ####
## 功能：将数据框中所有字符列截断到 Excel 单元格上限以下 (32,767)；超长以 [TRUNCATED] 标记 ####
## 输入：df，待 sanitize 数据框；max_chars，单元格最大字符数 (默认 32000，留缓冲) ####
## 输出：sanitize 后的数据框 ####
truncate_excel_string_columns <- function(df, max_chars = 32000L) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0 || ncol(df) == 0) {
    return(df)
  }

  for (current_col in colnames(df)) {
    current_value <- df[[current_col]]
    if (is.character(current_value) || is.factor(current_value)) {
      current_value_chr <- as.character(current_value)
      n_char_vec <- nchar(current_value_chr, type = "chars", allowNA = TRUE)
      over_index <- which(!is.na(n_char_vec) & n_char_vec > max_chars)
      if (length(over_index) > 0) {
        truncated_chr <- current_value_chr
        for (idx in over_index) {
          truncated_chr[[idx]] <- paste0(
            substr(current_value_chr[[idx]], 1, max_chars - 80L),
            " ... [TRUNCATED, full content in Rdata]"
          )
        }
        df[[current_col]] <- truncated_chr
      } else if (is.factor(current_value)) {
        df[[current_col]] <- current_value_chr
      }
    }
  }

  df
}

## 函数：sanitize_sheet_list_for_excel ####
## 功能：对 Excel 写出前的 sheet 列表批量应用 truncate_excel_string_columns ####
## 输入：sheet_list，命名列表；max_chars，单元格最大字符数 ####
## 输出：sanitize 后的命名列表 ####
sanitize_sheet_list_for_excel <- function(sheet_list, max_chars = 32000L) {
  lapply(sheet_list, truncate_excel_string_columns, max_chars = max_chars)
}


# 4. MSigDB 数据库构建函数 ####

## 函数：fetch_msigdbr_with_fallback ####
## 功能：兼容 msigdbr v25+ (gs_collection / gs_subcollection) 和旧版 (category / subcategory) 两种 API ####
## 输入：species，物种名；collection，主类 (如 "H"、"C5")；subcollection，子类 (可为 NULL 或多值) ####
## 输出：tibble，列含 gs_name、gene_symbol、gs_description 等 ####
fetch_msigdbr_with_fallback <- function(species, collection, subcollection = NULL) {
  call_msigdbr <- function(use_new_api, current_subcollection) {
    args_new <- list(species = species, collection = collection)
    args_old <- list(species = species, category = collection)

    if (!is.null(current_subcollection)) {
      args_new$subcollection <- current_subcollection
      args_old$subcategory <- current_subcollection
    }

    args_to_use <- if (use_new_api) args_new else args_old

    tryCatch(
      do.call(msigdbr::msigdbr, args_to_use),
      error = function(e) NULL
    )
  }

  fetch_one <- function(current_subcollection) {
    result <- call_msigdbr(use_new_api = TRUE, current_subcollection = current_subcollection)
    if (is.null(result) || nrow(result) == 0) {
      result <- call_msigdbr(use_new_api = FALSE, current_subcollection = current_subcollection)
    }
    result
  }

  if (is.null(subcollection)) {
    return(fetch_one(NULL))
  }

  result_list <- lapply(subcollection, fetch_one)
  result_list <- result_list[!vapply(result_list, is.null, logical(1))]

  if (length(result_list) == 0) {
    return(NULL)
  }

  combined_df <- dplyr::bind_rows(result_list)
  combined_df %>% dplyr::distinct(gs_name, gene_symbol, .keep_all = TRUE)
}

## 函数：build_database_resources ####
## 功能：根据 msigdb_database_config 构建每个数据库的 TERM2GENE、TERM2NAME 与统计信息 ####
## 输入：species，物种名；database_config，数据库配置命名列表 ####
## 输出：命名列表，每元素含 label、term2gene_df、term2name_df、n_terms、n_unique_genes ####
build_database_resources <- function(species, database_config) {
  result_list <- lapply(names(database_config), function(current_database) {
    current_config <- database_config[[current_database]]

    log_message(paste0("读取 MSigDB ", current_database, " (", current_config$label, ")..."))

    raw_df <- fetch_msigdbr_with_fallback(
      species = species,
      collection = current_config$collection,
      subcollection = current_config$subcollection
    )

    if (is.null(raw_df) || nrow(raw_df) == 0) {
      log_message(paste0("警告: ", current_database, " 未取到任何 gene set。"))
      return(list(
        label = current_config$label,
        term2gene_df = NULL,
        term2name_df = NULL,
        n_terms = 0L,
        n_unique_genes = 0L
      ))
    }

    term2gene_df <- raw_df %>%
      dplyr::select(gs_name, gene_symbol) %>%
      dplyr::filter(!is.na(gs_name) & !is.na(gene_symbol) & nchar(gs_name) > 0 & nchar(gene_symbol) > 0) %>%
      dplyr::distinct() %>%
      dplyr::rename(term = gs_name, gene = gene_symbol)

    term2name_df <- raw_df %>%
      dplyr::select(gs_name, gs_description) %>%
      dplyr::distinct(gs_name, .keep_all = TRUE) %>%
      dplyr::mutate(
        gs_description = ifelse(
          is.na(gs_description) | nchar(gs_description) == 0,
          gs_name,
          gs_description
        )
      ) %>%
      dplyr::rename(term = gs_name, name = gs_description)

    list(
      label = current_config$label,
      term2gene_df = tibble::as_tibble(term2gene_df),
      term2name_df = tibble::as_tibble(term2name_df),
      n_terms = dplyr::n_distinct(term2gene_df$term),
      n_unique_genes = dplyr::n_distinct(term2gene_df$gene)
    )
  })

  names(result_list) <- names(database_config)
  result_list
}

## 函数：build_database_summary_df ####
## 功能：汇总每个数据库的 term 数、独立基因数与配置标签 ####
## 输入：database_resource_list，build_database_resources 的输出 ####
## 输出：数据框 ####
build_database_summary_df <- function(database_resource_list) {
  summary_rows <- lapply(names(database_resource_list), function(current_database) {
    current_resource <- database_resource_list[[current_database]]
    data.frame(
      database = current_database,
      label = current_resource$label,
      n_terms = current_resource$n_terms,
      n_unique_genes = current_resource$n_unique_genes,
      stringsAsFactors = FALSE
    )
  })

  dplyr::bind_rows(summary_rows)
}


# 5. ORA 富集分析函数 ####

## 函数：parse_gene_ratio ####
## 功能：解析 clusterProfiler::enricher 返回的 GeneRatio / BgRatio 字符串 (例 "5/100") ####
## 输入：ratio_string，比值字符串向量 ####
## 输出：list，含 numerator 与 denominator 两数值向量 ####
parse_gene_ratio <- function(ratio_string) {
  parts <- strsplit(as.character(ratio_string), "/", fixed = TRUE)
  numeric_parts <- vapply(parts, function(p) {
    if (length(p) >= 1) suppressWarnings(as.numeric(p[[1]])) else NA_real_
  }, numeric(1))
  denom_parts <- vapply(parts, function(p) {
    if (length(p) >= 2) suppressWarnings(as.numeric(p[[2]])) else NA_real_
  }, numeric(1))
  list(numerator = numeric_parts, denominator = denom_parts)
}

## 函数：run_single_ora_task ####
## 功能：对单个 hit 基因列表 + 单个数据库执行 ORA (clusterProfiler::enricher) ####
## 输入：hit_genes，待富集基因；background_genes，背景基因；term2gene_df、term2name_df，TERM2GENE/TERM2NAME 表；
##       min_gs_size、max_gs_size、pvalue_cutoff、qvalue_cutoff，富集参数 ####
## 输出：富集结果数据框 (可能为 NULL) ####
run_single_ora_task <- function(hit_genes,
                                background_genes,
                                term2gene_df,
                                term2name_df,
                                min_gs_size,
                                max_gs_size,
                                pvalue_cutoff,
                                qvalue_cutoff) {
  hit_genes <- unique(as.character(hit_genes))
  hit_genes <- hit_genes[!is.na(hit_genes) & nchar(hit_genes) > 0]

  if (length(hit_genes) < min_hit_genes_for_enrichment) {
    return(NULL)
  }

  background_genes <- unique(as.character(background_genes))
  background_genes <- background_genes[!is.na(background_genes) & nchar(background_genes) > 0]

  enrich_result_obj <- tryCatch(
    {
      clusterProfiler::enricher(
        gene = hit_genes,
        universe = background_genes,
        TERM2GENE = term2gene_df,
        TERM2NAME = term2name_df,
        pAdjustMethod = "BH",
        pvalueCutoff = pvalue_cutoff,
        qvalueCutoff = qvalue_cutoff,
        minGSSize = min_gs_size,
        maxGSSize = max_gs_size
      )
    },
    error = function(e) {
      log_message(paste0("ORA enricher 调用失败: ", conditionMessage(e)))
      NULL
    }
  )

  if (is.null(enrich_result_obj)) {
    return(NULL)
  }

  enrich_result_df <- tryCatch(
    as.data.frame(enrich_result_obj@result),
    error = function(e) NULL
  )

  if (is.null(enrich_result_df) || nrow(enrich_result_df) == 0) {
    return(NULL)
  }

  gene_ratio_parsed <- parse_gene_ratio(enrich_result_df$GeneRatio)
  bg_ratio_parsed <- parse_gene_ratio(enrich_result_df$BgRatio)

  fold_enrichment <- ifelse(
    !is.na(bg_ratio_parsed$numerator) & bg_ratio_parsed$numerator > 0 &
      !is.na(bg_ratio_parsed$denominator) & bg_ratio_parsed$denominator > 0 &
      !is.na(gene_ratio_parsed$denominator) & gene_ratio_parsed$denominator > 0,
    (gene_ratio_parsed$numerator / gene_ratio_parsed$denominator) /
      (bg_ratio_parsed$numerator / bg_ratio_parsed$denominator),
    NA_real_
  )

  result_clean_df <- tibble::tibble(
    term_id = as.character(enrich_result_df$ID),
    term_description = as.character(enrich_result_df$Description),
    n_hit_in_term = as.integer(gene_ratio_parsed$numerator),
    n_hit_total_in_universe = as.integer(gene_ratio_parsed$denominator),
    n_term_in_universe = as.integer(bg_ratio_parsed$numerator),
    n_universe = as.integer(bg_ratio_parsed$denominator),
    fold_enrichment = fold_enrichment,
    enrichment_pvalue = as.numeric(enrich_result_df$pvalue),
    enrichment_padj_BH = as.numeric(enrich_result_df$p.adjust),
    enrichment_qvalue = as.numeric(enrich_result_df$qvalue),
    term_gene_symbols = as.character(enrich_result_df$geneID)
  ) %>%
    dplyr::arrange(enrichment_padj_BH, enrichment_pvalue)

  result_clean_df
}

## 函数：build_ora_task_grid ####
## 功能：构建 ORA 任务网格: track × delta_rho_threshold × pattern × hit_direction × database ####
## 输入：tracks; thresholds; patterns; directions; databases ####
## 输出：tidyr::expand_grid 数据框 ####
build_ora_task_grid <- function(tracks, thresholds, patterns, directions, databases) {
  tidyr::expand_grid(
    track = tracks,
    delta_rho_threshold = thresholds,
    pattern = patterns,
    hit_direction = directions,
    database = databases
  )
}

## 函数：run_all_ora_tasks ####
## 功能：遍历所有 (track × threshold × pattern × direction × database) 组合执行 ORA，返回长表 ####
## 输入：hit_long_df_candidate，候选 hit 长表；hit_long_df_hc，HC hit 长表；background_genes；
##       database_resource_list；tracks；thresholds；patterns；directions ####
## 输出：长表数据框 ####
run_all_ora_tasks <- function(hit_long_df_candidate,
                              hit_long_df_hc,
                              background_genes,
                              database_resource_list,
                              tracks,
                              thresholds,
                              patterns,
                              directions) {
  task_grid <- build_ora_task_grid(
    tracks = tracks,
    thresholds = thresholds,
    patterns = patterns,
    directions = directions,
    databases = names(database_resource_list)
  )

  log_message(paste0(
    "执行 ORA 任务，共 ", nrow(task_grid),
    " 个 (track × delta_rho × pattern × direction × database) 组合; ",
    "并行: ", enable_parallel, " (workers=", parallel_workers, ")"
  ))

  task_indices <- seq_len(nrow(task_grid))

  result_list <- run_parallel_lapply(
    x = task_indices,
    fun = function(task_idx) {
      current_track <- task_grid$track[[task_idx]]
      current_threshold <- task_grid$delta_rho_threshold[[task_idx]]
      current_pattern <- task_grid$pattern[[task_idx]]
      current_direction <- task_grid$hit_direction[[task_idx]]
      current_database <- task_grid$database[[task_idx]]

      current_resource <- database_resource_list[[current_database]]

      if (is.null(current_resource$term2gene_df) || nrow(current_resource$term2gene_df) == 0) {
        return(NULL)
      }

      hit_source_df <- if (current_track == "candidate") {
        hit_long_df_candidate
      } else {
        hit_long_df_hc
      }

      current_hit_genes <- hit_source_df %>%
        dplyr::filter(
          delta_rho_threshold == current_threshold,
          as.character(pattern) == current_pattern,
          hit_direction == current_direction
        ) %>%
        dplyr::pull(gene) %>%
        unique()

      n_hit_input_value <- length(current_hit_genes)

      if (n_hit_input_value < min_hit_genes_for_enrichment) {
        return(NULL)
      }

      single_result_df <- run_single_ora_task(
        hit_genes = current_hit_genes,
        background_genes = background_genes,
        term2gene_df = current_resource$term2gene_df,
        term2name_df = current_resource$term2name_df,
        min_gs_size = enrichment_min_gs_size,
        max_gs_size = enrichment_max_gs_size,
        pvalue_cutoff = enrichment_pvalue_cutoff,
        qvalue_cutoff = enrichment_qvalue_cutoff
      )

      if (is.null(single_result_df) || nrow(single_result_df) == 0) {
        return(NULL)
      }

      single_result_df %>%
        dplyr::mutate(
          track = current_track,
          delta_rho_threshold = current_threshold,
          pattern = current_pattern,
          hit_direction = current_direction,
          database = current_database,
          n_hit_input = n_hit_input_value
        ) %>%
        dplyr::select(
          track,
          delta_rho_threshold,
          pattern,
          hit_direction,
          database,
          n_hit_input,
          dplyr::everything()
        )
    },
    use_parallel = enable_parallel,
    workers = parallel_workers
  )

  ora_long_df <- dplyr::bind_rows(result_list) %>%
    dplyr::arrange(
      track,
      delta_rho_threshold,
      factor(pattern, levels = pattern_group_order),
      hit_direction,
      database,
      enrichment_padj_BH,
      enrichment_pvalue
    )

  ora_long_df
}


# 6. GSEA 富集分析函数 ####

## 函数：build_gsea_ranking_vector ####
## 功能：基于 7.7 全表 (merged_correlation_df) 构建单 pattern 的 GSEA 排序向量 ####
##       排序指标：有符号 Δρ。命中正向通路 NES > 0；负向通路 NES < 0。 ####
## 输入：merged_correlation_df，7.7 全表；current_pattern，目标 pattern；background_genes，限定基因池 ####
## 输出：命名 numeric 向量 (gene → Δρ)，已按 Δρ 降序排序 ####
build_gsea_ranking_vector <- function(merged_correlation_df, current_pattern, background_genes) {
  background_set <- unique(as.character(background_genes))

  ranking_df <- merged_correlation_df %>%
    dplyr::filter(as.character(pattern) == current_pattern) %>%
    dplyr::filter(as.character(gene) %in% background_set) %>%
    dplyr::filter(!is.na(delta_rho)) %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(delta_rho = mean(delta_rho, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(delta_rho))

  ranking_vector <- ranking_df$delta_rho
  names(ranking_vector) <- ranking_df$gene
  ranking_vector
}

## 函数：run_single_gsea_task ####
## 功能：对单个 pattern × 单个数据库执行 GSEA (clusterProfiler::GSEA) ####
## 输入：ranking_vector，命名 numeric 向量；term2gene_df、term2name_df；min_gs_size、max_gs_size；
##       n_perm，permutation 次数；eps，fgsea 精度；seed，随机种子 ####
## 输出：GSEA 结果数据框 (可能为 NULL) ####
run_single_gsea_task <- function(ranking_vector,
                                 term2gene_df,
                                 term2name_df,
                                 min_gs_size,
                                 max_gs_size,
                                 n_perm,
                                 eps,
                                 seed) {
  if (length(ranking_vector) < min_gs_size) {
    return(NULL)
  }

  set.seed(seed)
  gsea_result_obj <- tryCatch(
    {
      clusterProfiler::GSEA(
        geneList = ranking_vector,
        TERM2GENE = term2gene_df,
        TERM2NAME = term2name_df,
        pAdjustMethod = "BH",
        pvalueCutoff = 1,
        minGSSize = min_gs_size,
        maxGSSize = max_gs_size,
        nPermSimple = n_perm,
        eps = eps,
        seed = TRUE,
        by = "fgsea",
        verbose = FALSE
      )
    },
    error = function(e) {
      log_message(paste0("GSEA 调用失败: ", conditionMessage(e)))
      NULL
    }
  )

  if (is.null(gsea_result_obj)) {
    return(NULL)
  }

  gsea_result_df <- tryCatch(
    as.data.frame(gsea_result_obj@result),
    error = function(e) NULL
  )

  if (is.null(gsea_result_df) || nrow(gsea_result_df) == 0) {
    return(NULL)
  }

  result_clean_df <- tibble::tibble(
    term_id = as.character(gsea_result_df$ID),
    term_description = as.character(gsea_result_df$Description),
    set_size = as.integer(gsea_result_df$setSize),
    enrichment_score = as.numeric(gsea_result_df$enrichmentScore),
    nes = as.numeric(gsea_result_df$NES),
    enrichment_pvalue = as.numeric(gsea_result_df$pvalue),
    enrichment_padj_BH = as.numeric(gsea_result_df$p.adjust),
    enrichment_qvalue = as.numeric(gsea_result_df$qvalue),
    rank_at_max = as.integer(gsea_result_df$rank),
    leading_edge_genes = as.character(gsea_result_df$core_enrichment),
    nes_direction = ifelse(as.numeric(gsea_result_df$NES) > 0, "positive", "negative")
  ) %>%
    dplyr::arrange(enrichment_padj_BH, enrichment_pvalue)

  result_clean_df
}

## 函数：run_all_gsea_tasks ####
## 功能：遍历 track × pattern × database 组合执行 GSEA ####
## 输入：merged_correlation_df，7.7 全表；hit_genes_per_track_pattern，list of list 用于过滤 ranking 池；
##       background_genes_candidate / background_genes_hc，每 track 的 ranking 池；
##       database_resource_list；tracks；patterns ####
## 输出：长表数据框 ####
run_all_gsea_tasks <- function(merged_correlation_df,
                               background_genes_candidate,
                               background_genes_hc,
                               database_resource_list,
                               tracks,
                               patterns) {
  task_grid <- tidyr::expand_grid(
    track = tracks,
    pattern = patterns,
    database = names(database_resource_list)
  )

  log_message(paste0(
    "执行 GSEA 任务，共 ", nrow(task_grid),
    " 个 (track × pattern × database) 组合; ",
    "并行: ", enable_parallel, " (workers=", parallel_workers, ")"
  ))

  task_indices <- seq_len(nrow(task_grid))

  result_list <- run_parallel_lapply(
    x = task_indices,
    fun = function(task_idx) {
      current_track <- task_grid$track[[task_idx]]
      current_pattern <- task_grid$pattern[[task_idx]]
      current_database <- task_grid$database[[task_idx]]

      current_resource <- database_resource_list[[current_database]]

      if (is.null(current_resource$term2gene_df) || nrow(current_resource$term2gene_df) == 0) {
        return(NULL)
      }

      ranking_pool <- if (current_track == "candidate") {
        background_genes_candidate
      } else {
        background_genes_hc
      }

      ranking_vector <- build_gsea_ranking_vector(
        merged_correlation_df = merged_correlation_df,
        current_pattern = current_pattern,
        background_genes = ranking_pool
      )

      n_ranking_input <- length(ranking_vector)

      if (n_ranking_input < enrichment_min_gs_size) {
        return(NULL)
      }

      single_result_df <- run_single_gsea_task(
        ranking_vector = ranking_vector,
        term2gene_df = current_resource$term2gene_df,
        term2name_df = current_resource$term2name_df,
        min_gs_size = enrichment_min_gs_size,
        max_gs_size = enrichment_max_gs_size,
        n_perm = gsea_n_perm,
        eps = gsea_eps,
        seed = gsea_seed
      )

      if (is.null(single_result_df) || nrow(single_result_df) == 0) {
        return(NULL)
      }

      single_result_df %>%
        dplyr::mutate(
          track = current_track,
          pattern = current_pattern,
          database = current_database,
          n_ranking_input = n_ranking_input
        ) %>%
        dplyr::select(
          track,
          pattern,
          database,
          n_ranking_input,
          dplyr::everything()
        )
    },
    use_parallel = enable_parallel,
    workers = parallel_workers
  )

  gsea_long_df <- dplyr::bind_rows(result_list) %>%
    dplyr::arrange(
      track,
      factor(pattern, levels = pattern_group_order),
      database,
      enrichment_padj_BH,
      enrichment_pvalue
    )

  gsea_long_df
}


# 7. 生物学 positive control 函数 ####

## 函数：collect_positive_control_term_ids ####
## 功能：根据每模式 keyword 列表，从 5 库 MSigDB term name 中匹配出 expected pathway gene set ####
## 输入：positive_control_keywords，每模式 keyword 列表；database_resource_list，5 库资源 ####
## 输出：长表数据框: pattern | database | term_id | matched_keyword | n_genes ####
collect_positive_control_term_ids <- function(positive_control_keywords, database_resource_list) {
  all_rows <- list()

  for (current_pattern in names(positive_control_keywords)) {
    current_keywords <- positive_control_keywords[[current_pattern]]

    for (current_database in names(database_resource_list)) {
      current_resource <- database_resource_list[[current_database]]
      term2gene_df <- current_resource$term2gene_df

      if (is.null(term2gene_df) || nrow(term2gene_df) == 0) {
        next
      }

      term_size_df <- term2gene_df %>%
        dplyr::group_by(term) %>%
        dplyr::summarise(n_genes = dplyr::n_distinct(gene), .groups = "drop")

      for (current_keyword in current_keywords) {
        matched_terms <- term_size_df$term[
          grepl(current_keyword, term_size_df$term, ignore.case = TRUE, fixed = FALSE)
        ]

        if (length(matched_terms) == 0) next

        for (one_term in matched_terms) {
          all_rows[[length(all_rows) + 1]] <- data.frame(
            pattern = current_pattern,
            database = current_database,
            term_id = one_term,
            matched_keyword = current_keyword,
            n_genes = term_size_df$n_genes[term_size_df$term == one_term][1],
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  if (length(all_rows) == 0) {
    return(data.frame(
      pattern = character(0),
      database = character(0),
      term_id = character(0),
      matched_keyword = character(0),
      n_genes = integer(0),
      stringsAsFactors = FALSE
    ))
  }

  dplyr::bind_rows(all_rows) %>%
    dplyr::distinct(pattern, database, term_id, .keep_all = TRUE) %>%
    dplyr::arrange(
      factor(pattern, levels = pattern_group_order),
      database,
      term_id
    )
}

## 函数：run_positive_control_test ####
## 功能：对每模式 expected gene set，用 fisher.test 检验 hit 是否优先富集于该集 ####
## 输入：hit_genes_by_pattern，命名 list (pattern → 该模式 hit 基因)；
##       background_genes，背景基因；positive_control_term_df，匹配的 expected term；
##       database_resource_list，5 库资源 ####
## 输出：长表数据框: pattern | database | term_id | matched_keyword | n_term | n_hit | n_overlap |
##                  fisher_odds_ratio | fisher_p | fisher_padj_BH ####
run_positive_control_test <- function(hit_genes_by_pattern,
                                      background_genes,
                                      positive_control_term_df,
                                      database_resource_list) {
  if (nrow(positive_control_term_df) == 0) {
    return(data.frame(
      pattern = character(0),
      database = character(0),
      term_id = character(0),
      matched_keyword = character(0),
      n_term_in_universe = integer(0),
      n_hit_input = integer(0),
      n_overlap = integer(0),
      n_universe = integer(0),
      fisher_odds_ratio = numeric(0),
      fisher_p_value = numeric(0),
      fisher_padj_BH = numeric(0),
      overlap_genes = character(0),
      stringsAsFactors = FALSE
    ))
  }

  background_set <- unique(as.character(background_genes))
  background_set <- background_set[!is.na(background_set) & nchar(background_set) > 0]
  n_universe <- length(background_set)

  result_rows <- list()

  for (i in seq_len(nrow(positive_control_term_df))) {
    current_pattern <- positive_control_term_df$pattern[[i]]
    current_database <- positive_control_term_df$database[[i]]
    current_term <- positive_control_term_df$term_id[[i]]
    current_keyword <- positive_control_term_df$matched_keyword[[i]]

    current_resource <- database_resource_list[[current_database]]
    term2gene_df <- current_resource$term2gene_df

    term_genes <- term2gene_df$gene[term2gene_df$term == current_term]
    term_genes <- intersect(unique(term_genes), background_set)
    n_term_in_universe <- length(term_genes)

    hit_genes <- unique(as.character(hit_genes_by_pattern[[current_pattern]]))
    hit_genes <- intersect(hit_genes, background_set)
    n_hit_input <- length(hit_genes)

    overlap_genes <- intersect(term_genes, hit_genes)
    n_overlap <- length(overlap_genes)

    if (n_term_in_universe < 1 || n_hit_input < 1) {
      next
    }

    a <- n_overlap
    b <- n_hit_input - n_overlap
    c <- n_term_in_universe - n_overlap
    d <- n_universe - a - b - c

    fisher_obj <- tryCatch(
      stats::fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater"),
      error = function(e) NULL
    )

    if (is.null(fisher_obj)) {
      odds_ratio <- NA_real_
      p_value <- NA_real_
    } else {
      odds_ratio <- as.numeric(fisher_obj$estimate)
      p_value <- as.numeric(fisher_obj$p.value)
    }

    result_rows[[length(result_rows) + 1]] <- data.frame(
      pattern = current_pattern,
      database = current_database,
      term_id = current_term,
      matched_keyword = current_keyword,
      n_term_in_universe = as.integer(n_term_in_universe),
      n_hit_input = as.integer(n_hit_input),
      n_overlap = as.integer(n_overlap),
      n_universe = as.integer(n_universe),
      fisher_odds_ratio = odds_ratio,
      fisher_p_value = p_value,
      overlap_genes = paste(sort(overlap_genes), collapse = "/"),
      stringsAsFactors = FALSE
    )
  }

  if (length(result_rows) == 0) {
    return(data.frame())
  }

  result_df <- dplyr::bind_rows(result_rows)

  result_df <- result_df %>%
    dplyr::group_by(pattern) %>%
    dplyr::mutate(fisher_padj_BH = stats::p.adjust(fisher_p_value, method = "BH")) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(
      factor(pattern, levels = pattern_group_order),
      fisher_padj_BH,
      database,
      term_id
    ) %>%
    dplyr::select(
      pattern,
      database,
      term_id,
      matched_keyword,
      n_term_in_universe,
      n_hit_input,
      n_overlap,
      n_universe,
      fisher_odds_ratio,
      fisher_p_value,
      fisher_padj_BH,
      overlap_genes
    )

  result_df
}


# 8. 跨模式综合分析函数 ####

## 函数：find_cross_pattern_shared_pathways ####
## 功能：在主阈值 candidate hit 主轨下，识别在 ≥ min_n_patterns_shared 个模式中同时显著的通路 ####
## 输入：ora_long_df，ORA 长表；track_target；threshold_target；fdr_threshold；min_n_patterns_shared ####
## 输出：跨模式共享通路数据框 ####
find_cross_pattern_shared_pathways <- function(ora_long_df,
                                               track_target,
                                               threshold_target,
                                               fdr_threshold,
                                               min_n_patterns_shared) {
  significant_df <- ora_long_df %>%
    dplyr::filter(
      track == track_target,
      delta_rho_threshold == threshold_target,
      !is.na(enrichment_padj_BH) & enrichment_padj_BH < fdr_threshold
    )

  if (nrow(significant_df) == 0) {
    return(data.frame(
      database = character(0),
      term_id = character(0),
      term_description = character(0),
      hit_direction = character(0),
      n_patterns_shared = integer(0),
      patterns_list = character(0),
      mean_padj = numeric(0),
      min_padj = numeric(0),
      max_fold_enrichment = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  shared_pathways_df <- significant_df %>%
    dplyr::group_by(database, term_id, term_description, hit_direction) %>%
    dplyr::summarise(
      n_patterns_shared = dplyr::n_distinct(pattern),
      patterns_list = paste(
        sort(unique(as.character(pattern))),
        collapse = ", "
      ),
      mean_padj = mean(enrichment_padj_BH, na.rm = TRUE),
      min_padj = min(enrichment_padj_BH, na.rm = TRUE),
      max_fold_enrichment = max(fold_enrichment, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::filter(n_patterns_shared >= min_n_patterns_shared) %>%
    dplyr::arrange(dplyr::desc(n_patterns_shared), mean_padj, database, term_id)

  shared_pathways_df
}

## 函数：find_cross_pattern_unique_pathways ####
## 功能：在主阈值 candidate hit 主轨下，识别仅在 1 个模式显著的"模式特异通路" ####
## 输入：ora_long_df；track_target；threshold_target；fdr_threshold ####
## 输出：模式特异通路长表 ####
find_cross_pattern_unique_pathways <- function(ora_long_df,
                                               track_target,
                                               threshold_target,
                                               fdr_threshold) {
  significant_df <- ora_long_df %>%
    dplyr::filter(
      track == track_target,
      delta_rho_threshold == threshold_target,
      !is.na(enrichment_padj_BH) & enrichment_padj_BH < fdr_threshold
    )

  if (nrow(significant_df) == 0) {
    return(significant_df)
  }

  unique_term_df <- significant_df %>%
    dplyr::group_by(database, term_id, term_description, hit_direction) %>%
    dplyr::summarise(
      n_patterns_significant = dplyr::n_distinct(pattern),
      .groups = "drop"
    ) %>%
    dplyr::filter(n_patterns_significant == 1)

  unique_pathways_df <- significant_df %>%
    dplyr::semi_join(unique_term_df, by = c("database", "term_id", "term_description", "hit_direction")) %>%
    dplyr::arrange(
      factor(pattern, levels = pattern_group_order),
      hit_direction,
      database,
      enrichment_padj_BH
    )

  unique_pathways_df
}

## 函数：build_top_pathways_per_pattern ####
## 功能：在主阈值 candidate 主轨下，为每个 (pattern × direction × database) 提取 Top N 显著通路 ####
## 输入：ora_long_df；track_target；threshold_target；fdr_threshold；top_n ####
## 输出：Top 通路长表 ####
build_top_pathways_per_pattern <- function(ora_long_df,
                                           track_target,
                                           threshold_target,
                                           fdr_threshold,
                                           top_n) {
  significant_df <- ora_long_df %>%
    dplyr::filter(
      track == track_target,
      delta_rho_threshold == threshold_target,
      !is.na(enrichment_padj_BH) & enrichment_padj_BH < fdr_threshold
    )

  if (nrow(significant_df) == 0) {
    return(significant_df)
  }

  top_df <- significant_df %>%
    dplyr::group_by(pattern, hit_direction, database) %>%
    dplyr::slice_min(
      order_by = enrichment_padj_BH,
      n = top_n,
      with_ties = FALSE
    ) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(
      factor(pattern, levels = pattern_group_order),
      hit_direction,
      database,
      enrichment_padj_BH
    )

  top_df
}

## 函数：build_per_pattern_database_summary_ora ####
## 功能：每 (track × threshold × pattern × direction × database) ORA hit 输入数与显著通路数汇总 ####
## 输入：ora_long_df；threshold_target；fdr_threshold ####
## 输出：宽表数据框 ####
build_per_pattern_database_summary_ora <- function(ora_long_df,
                                                   threshold_target,
                                                   fdr_threshold) {
  ora_main_df <- ora_long_df %>%
    dplyr::filter(delta_rho_threshold == threshold_target)

  if (nrow(ora_main_df) == 0) {
    return(data.frame())
  }

  summary_df <- ora_main_df %>%
    dplyr::group_by(track, pattern, hit_direction, database) %>%
    dplyr::summarise(
      n_hit_input = dplyr::first(n_hit_input),
      n_terms_tested = dplyr::n_distinct(term_id),
      n_terms_significant = sum(
        !is.na(enrichment_padj_BH) & enrichment_padj_BH < fdr_threshold,
        na.rm = TRUE
      ),
      min_padj = min(enrichment_padj_BH, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(
      track,
      factor(pattern, levels = pattern_group_order),
      hit_direction,
      database
    )

  summary_df
}


# 9. ORA-GSEA 一致性度量 ####

## 函数：compute_ora_gsea_agreement ####
## 功能：对每 (track × pattern × database)，计算 ORA 显著通路与 GSEA 显著通路的 Jaccard 系数 ####
## 输入：ora_long_df，主阈值过滤后的 ORA 长表；gsea_long_df，GSEA 长表；fdr_threshold ####
## 输出：长表数据框: track | pattern | database | n_ora_sig | n_gsea_sig | n_overlap | jaccard ####
compute_ora_gsea_agreement <- function(ora_long_df,
                                       gsea_long_df,
                                       threshold_target,
                                       fdr_threshold) {
  ora_sig_df <- ora_long_df %>%
    dplyr::filter(
      delta_rho_threshold == threshold_target,
      !is.na(enrichment_padj_BH) & enrichment_padj_BH < fdr_threshold
    ) %>%
    dplyr::group_by(track, pattern, database) %>%
    dplyr::summarise(
      ora_terms = list(unique(term_id)),
      .groups = "drop"
    )

  gsea_sig_df <- gsea_long_df %>%
    dplyr::filter(!is.na(enrichment_padj_BH) & enrichment_padj_BH < fdr_threshold) %>%
    dplyr::group_by(track, pattern, database) %>%
    dplyr::summarise(
      gsea_terms = list(unique(term_id)),
      .groups = "drop"
    )

  combined_df <- ora_sig_df %>%
    dplyr::full_join(gsea_sig_df, by = c("track", "pattern", "database"))

  combined_df$n_ora_sig <- vapply(combined_df$ora_terms, function(x) {
    if (is.null(x)) 0L else length(x)
  }, integer(1))

  combined_df$n_gsea_sig <- vapply(combined_df$gsea_terms, function(x) {
    if (is.null(x)) 0L else length(x)
  }, integer(1))

  combined_df$n_overlap <- mapply(function(a, b) {
    if (is.null(a) || is.null(b)) return(0L)
    length(intersect(a, b))
  }, combined_df$ora_terms, combined_df$gsea_terms)

  combined_df$n_union <- mapply(function(a, b) {
    a_v <- if (is.null(a)) character(0) else a
    b_v <- if (is.null(b)) character(0) else b
    length(union(a_v, b_v))
  }, combined_df$ora_terms, combined_df$gsea_terms)

  combined_df$jaccard <- ifelse(
    combined_df$n_union == 0,
    NA_real_,
    combined_df$n_overlap / combined_df$n_union
  )

  combined_df %>%
    dplyr::mutate(
      ora_terms = vapply(ora_terms, function(x) {
        if (is.null(x)) "" else paste(sort(unique(x)), collapse = "; ")
      }, character(1)),
      gsea_terms = vapply(gsea_terms, function(x) {
        if (is.null(x)) "" else paste(sort(unique(x)), collapse = "; ")
      }, character(1))
    ) %>%
    dplyr::select(
      track,
      pattern,
      database,
      n_ora_sig,
      n_gsea_sig,
      n_overlap,
      n_union,
      jaccard,
      ora_terms,
      gsea_terms
    ) %>%
    dplyr::arrange(
      track,
      factor(pattern, levels = pattern_group_order),
      database
    )
}


# 10. 主流程 ####

log_message("========== 7.7.1 模式相关基因富集分析开始 ==========")
log_message(paste0(
  "并行设置: enable_parallel=", enable_parallel,
  "; parallel_workers=", parallel_workers
))

assert_required_files(
  file_paths = c(
    final_pattern_gene_input_rdata,
    gene_filter_summary_input_rdata
  ),
  file_labels = c(
    "7.7 主筛选结果 (final_pattern_gene_list.Rdata)",
    "7.7 基因过滤汇总 (gene_filter_summary_df.Rdata)"
  )
)

dependency_check_df <- check_required_packages(
  packages = required_packages,
  output_path = dependency_xlsx
)
log_message("依赖检查通过。")


## 10.1 读取 7.7 输出 ####
log_message("读取 7.7 输出: hit 列表、全表 merged_correlation_df 与基因过滤汇总...")

final_pattern_gene_objects <- load_named_objects_from_rdata(
  rdata_file = final_pattern_gene_input_rdata,
  object_names = c(
    "hit_long_df",
    "hit_long_high_confidence_df",
    "hit_count_df",
    "merged_correlation_df",
    "multi_pattern_classification_df",
    "threshold_filter_df"
  )
)

hit_long_df_candidate <- final_pattern_gene_objects$hit_long_df
hit_long_df_hc <- final_pattern_gene_objects$hit_long_high_confidence_df
hit_count_df <- final_pattern_gene_objects$hit_count_df
merged_correlation_df <- final_pattern_gene_objects$merged_correlation_df

gene_filter_summary_objects <- load_named_objects_from_rdata(
  rdata_file = gene_filter_summary_input_rdata,
  object_names = c("gene_filter_summary_df")
)
gene_filter_summary_df <- gene_filter_summary_objects$gene_filter_summary_df


## 10.2 准备背景基因 (universe) ####
background_genes <- gene_filter_summary_df %>%
  dplyr::filter(retained) %>%
  dplyr::pull(gene) %>%
  unique() %>%
  as.character()

background_genes <- background_genes[!is.na(background_genes) & nchar(background_genes) > 0]

if (length(background_genes) < 100) {
  stop("背景基因数过少 (< 100)，无法做富集分析。当前: ", length(background_genes))
}

log_message(paste0("富集分析背景基因数: ", length(background_genes)))

bg_summary_df <- tibble::tibble(
  指标 = c(
    "7.7 输入基因总数",
    "通过预过滤基因数 (背景)",
    "主阈值 candidate hit 基因总数 (去重)",
    "主阈值 HC hit 基因总数 (去重)",
    "delta_rho_threshold_main"
  ),
  值 = c(
    as.character(nrow(gene_filter_summary_df)),
    as.character(length(background_genes)),
    as.character(dplyr::n_distinct(
      hit_long_df_candidate$gene[hit_long_df_candidate$delta_rho_threshold == delta_rho_threshold_main]
    )),
    as.character(dplyr::n_distinct(
      hit_long_df_hc$gene[hit_long_df_hc$delta_rho_threshold == delta_rho_threshold_main]
    )),
    as.character(delta_rho_threshold_main)
  )
)


## 10.3 hit 输入汇总 ####
hit_input_summary_candidate_df <- hit_long_df_candidate %>%
  dplyr::group_by(delta_rho_threshold, pattern, hit_direction) %>%
  dplyr::summarise(
    n_hit_genes = dplyr::n_distinct(gene),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    track = "candidate",
    pattern = as.character(pattern)
  ) %>%
  dplyr::select(track, dplyr::everything()) %>%
  dplyr::arrange(
    delta_rho_threshold,
    factor(pattern, levels = pattern_group_order),
    hit_direction
  )

hit_input_summary_hc_df <- hit_long_df_hc %>%
  dplyr::group_by(delta_rho_threshold, pattern, hit_direction) %>%
  dplyr::summarise(
    n_hit_genes = dplyr::n_distinct(gene),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    track = "high_confidence",
    pattern = as.character(pattern)
  ) %>%
  dplyr::select(track, dplyr::everything()) %>%
  dplyr::arrange(
    delta_rho_threshold,
    factor(pattern, levels = pattern_group_order),
    hit_direction
  )

hit_input_summary_df <- dplyr::bind_rows(
  hit_input_summary_candidate_df,
  hit_input_summary_hc_df
)

log_message("各 (track × 阈值 × pattern × 方向) hit 输入数量:")
for (i in seq_len(nrow(hit_input_summary_df))) {
  log_message(paste0(
    "  ", hit_input_summary_df$track[[i]],
    " | Δρ≥", hit_input_summary_df$delta_rho_threshold[[i]],
    " | ", hit_input_summary_df$pattern[[i]],
    " | ", hit_input_summary_df$hit_direction[[i]],
    " : ", hit_input_summary_df$n_hit_genes[[i]]
  ))
}


## 10.4 构建数据库资源 ####
log_message("构建 5 个 MSigDB 子集的 TERM2GENE / TERM2NAME...")

if (file.exists(database_term2gene_rdata)) {
  log_message(paste0(
    "检测到已有数据库资源文件，直接读取: ", database_term2gene_rdata
  ))
  loaded_db_objects <- load_named_objects_from_rdata(
    rdata_file = database_term2gene_rdata,
    object_names = c("database_resource_list", "database_summary_df")
  )
  database_resource_list <- loaded_db_objects$database_resource_list
  database_summary_df <- loaded_db_objects$database_summary_df
} else {
  database_resource_list <- build_database_resources(
    species = msigdb_species,
    database_config = msigdb_database_config
  )
  database_summary_df <- build_database_summary_df(database_resource_list)
  save(database_resource_list, database_summary_df, file = database_term2gene_rdata)
}

log_message("数据库统计:")
for (i in seq_len(nrow(database_summary_df))) {
  log_message(paste0(
    "  ", database_summary_df$database[[i]],
    " (", database_summary_df$label[[i]], ")",
    " : n_terms=", database_summary_df$n_terms[[i]],
    "; n_unique_genes=", database_summary_df$n_unique_genes[[i]]
  ))
}


## 10.5 执行 ORA 任务 ####
unique_thresholds_in_data <- sort(unique(hit_long_df_candidate$delta_rho_threshold))
thresholds_to_run <- intersect(delta_rho_thresholds_evaluated, unique_thresholds_in_data)
if (length(thresholds_to_run) == 0) {
  thresholds_to_run <- unique_thresholds_in_data
}
log_message(paste0(
  "ORA 待跑 Δρ 阈值: ", paste(thresholds_to_run, collapse = ", ")
))

if (file.exists(ora_long_df_rdata)) {
  log_message(paste0("检测到已有 ORA 结果文件，直接读取: ", ora_long_df_rdata))
  ora_long_df <- load_named_objects_from_rdata(
    rdata_file = ora_long_df_rdata,
    object_names = c("ora_long_df")
  )$ora_long_df
} else {
  ora_long_df <- run_all_ora_tasks(
    hit_long_df_candidate = hit_long_df_candidate,
    hit_long_df_hc = hit_long_df_hc,
    background_genes = background_genes,
    database_resource_list = database_resource_list,
    tracks = hit_track_levels,
    thresholds = thresholds_to_run,
    patterns = pattern_group_order,
    directions = hit_direction_levels
  )
  save(ora_long_df, file = ora_long_df_rdata)
  log_message(paste0("ORA 结果总表已保存: ", ora_long_df_rdata))
}

log_message(paste0("ORA 结果总表行数: ", nrow(ora_long_df)))


## 10.6 执行 GSEA 任务 ####
if (file.exists(gsea_long_df_rdata)) {
  log_message(paste0("检测到已有 GSEA 结果文件，直接读取: ", gsea_long_df_rdata))
  gsea_long_df <- load_named_objects_from_rdata(
    rdata_file = gsea_long_df_rdata,
    object_names = c("gsea_long_df")
  )$gsea_long_df
} else {
  ## GSEA 主轨 (candidate) 用全部 background_genes 作 ranking 池
  ## GSEA 辅轨 (high_confidence) 用 HC 主阈值任一方向至少有 hit 的基因 ∪ 当前模式背景作 ranking 池
  ##   (这里简化: 两轨道都用同一 background_genes 池; track 在 GSEA 主要起方向 sanity check 作用,
  ##    GSEA 不依赖 hit 集，依赖排序，故两轨结果实际相同；保留 track 维度便于未来差异化扩展)
  gsea_long_df <- run_all_gsea_tasks(
    merged_correlation_df = merged_correlation_df,
    background_genes_candidate = background_genes,
    background_genes_hc = background_genes,
    database_resource_list = database_resource_list,
    tracks = hit_track_levels,
    patterns = pattern_group_order
  )
  save(gsea_long_df, file = gsea_long_df_rdata)
  log_message(paste0("GSEA 结果总表已保存: ", gsea_long_df_rdata))
}

log_message(paste0("GSEA 结果总表行数: ", nrow(gsea_long_df)))


## 10.7 生物学 positive control ####
if (file.exists(positive_control_rdata)) {
  log_message(paste0("检测到已有 positive control 结果文件，直接读取: ", positive_control_rdata))
  loaded_pc_objects <- load_named_objects_from_rdata(
    rdata_file = positive_control_rdata,
    object_names = c("positive_control_term_df", "positive_control_test_df")
  )
  positive_control_term_df <- loaded_pc_objects$positive_control_term_df
  positive_control_test_df <- loaded_pc_objects$positive_control_test_df
} else {
  log_message("收集 positive control expected gene sets (按预定义 keyword 匹配 5 库 term)...")
  positive_control_term_df <- collect_positive_control_term_ids(
    positive_control_keywords = positive_control_keywords,
    database_resource_list = database_resource_list
  )
  log_message(paste0(
    "Positive control: 共匹配 ", nrow(positive_control_term_df),
    " 个 (pattern × term) 预定义集"
  ))

  ## 主阈值 candidate hit 按 pattern 聚合 (双向合并，positive control 关心命中富集而非方向)
  hit_genes_by_pattern_candidate <- list()
  for (current_pattern in pattern_group_order) {
    hit_genes_by_pattern_candidate[[current_pattern]] <- hit_long_df_candidate %>%
      dplyr::filter(
        delta_rho_threshold == delta_rho_threshold_main,
        as.character(pattern) == current_pattern
      ) %>%
      dplyr::pull(gene) %>%
      unique()
  }

  positive_control_test_df <- run_positive_control_test(
    hit_genes_by_pattern = hit_genes_by_pattern_candidate,
    background_genes = background_genes,
    positive_control_term_df = positive_control_term_df,
    database_resource_list = database_resource_list
  )

  save(positive_control_term_df, positive_control_test_df, file = positive_control_rdata)
  log_message(paste0("Positive control 检验结果已保存: ", positive_control_rdata))
}

log_message(paste0(
  "Positive control 检验通路数 (有效): ", nrow(positive_control_test_df)
))


## 10.8 跨模式综合 (基于 candidate 主轨主阈值) ####
log_message("跨模式综合分析: 共享 / 模式特异 / Top N (candidate 主轨主阈值)...")

shared_pathways_df <- find_cross_pattern_shared_pathways(
  ora_long_df = ora_long_df,
  track_target = "candidate",
  threshold_target = delta_rho_threshold_main,
  fdr_threshold = enrichment_fdr_threshold,
  min_n_patterns_shared = min_n_patterns_for_shared
)
save(shared_pathways_df, file = cross_pattern_shared_rdata)
log_message(paste0(
  "跨模式 (≥ ", min_n_patterns_for_shared, " 模式共享) 显著通路数: ",
  nrow(shared_pathways_df)
))

unique_pathways_df <- find_cross_pattern_unique_pathways(
  ora_long_df = ora_long_df,
  track_target = "candidate",
  threshold_target = delta_rho_threshold_main,
  fdr_threshold = enrichment_fdr_threshold
)
save(unique_pathways_df, file = cross_pattern_unique_rdata)
log_message(paste0("模式特异 (仅 1 模式显著) 通路条目数: ", nrow(unique_pathways_df)))

top_pathways_df <- build_top_pathways_per_pattern(
  ora_long_df = ora_long_df,
  track_target = "candidate",
  threshold_target = delta_rho_threshold_main,
  fdr_threshold = enrichment_fdr_threshold,
  top_n = top_n_pathway_per_pattern
)
save(top_pathways_df, file = top_pathways_rdata)
log_message(paste0(
  "每 (pattern × direction × database) Top ", top_n_pathway_per_pattern,
  " 通路汇总行数: ", nrow(top_pathways_df)
))


## 10.9 ORA-GSEA 一致性 ####
log_message("计算 ORA vs GSEA 显著通路 Jaccard 一致性...")
ora_gsea_agreement_df <- compute_ora_gsea_agreement(
  ora_long_df = ora_long_df,
  gsea_long_df = gsea_long_df,
  threshold_target = delta_rho_threshold_main,
  fdr_threshold = enrichment_fdr_threshold
)
save(ora_gsea_agreement_df, file = ora_gsea_agreement_rdata)
log_message(paste0("ORA-GSEA 一致性表行数: ", nrow(ora_gsea_agreement_df)))


## 10.10 主阈值显著通路提取 ####
significant_ora_main_candidate_df <- ora_long_df %>%
  dplyr::filter(
    track == "candidate",
    delta_rho_threshold == delta_rho_threshold_main,
    !is.na(enrichment_padj_BH) & enrichment_padj_BH < enrichment_fdr_threshold
  ) %>%
  dplyr::arrange(
    factor(pattern, levels = pattern_group_order),
    hit_direction,
    database,
    enrichment_padj_BH
  )

significant_ora_main_hc_df <- ora_long_df %>%
  dplyr::filter(
    track == "high_confidence",
    delta_rho_threshold == delta_rho_threshold_main,
    !is.na(enrichment_padj_BH) & enrichment_padj_BH < enrichment_fdr_threshold
  ) %>%
  dplyr::arrange(
    factor(pattern, levels = pattern_group_order),
    hit_direction,
    database,
    enrichment_padj_BH
  )

significant_gsea_candidate_df <- gsea_long_df %>%
  dplyr::filter(
    track == "candidate",
    !is.na(enrichment_padj_BH) & enrichment_padj_BH < enrichment_fdr_threshold
  ) %>%
  dplyr::arrange(
    factor(pattern, levels = pattern_group_order),
    database,
    enrichment_padj_BH
  )

significant_gsea_hc_df <- gsea_long_df %>%
  dplyr::filter(
    track == "high_confidence",
    !is.na(enrichment_padj_BH) & enrichment_padj_BH < enrichment_fdr_threshold
  ) %>%
  dplyr::arrange(
    factor(pattern, levels = pattern_group_order),
    database,
    enrichment_padj_BH
  )

log_message(paste0(
  "主阈值 candidate hit ORA 显著通路数 (FDR<", enrichment_fdr_threshold, "): ",
  nrow(significant_ora_main_candidate_df)
))
log_message(paste0(
  "主阈值 HC ORA 显著通路数 (FDR<", enrichment_fdr_threshold, "): ",
  nrow(significant_ora_main_hc_df)
))
log_message(paste0(
  "candidate GSEA 显著通路数: ", nrow(significant_gsea_candidate_df)
))
log_message(paste0(
  "HC GSEA 显著通路数: ", nrow(significant_gsea_hc_df)
))


## 10.11 阈值敏感性 (candidate ORA 在不同 Δρ 阈值下显著通路数) ####
threshold_sensitivity_df <- ora_long_df %>%
  dplyr::filter(track == "candidate") %>%
  dplyr::group_by(delta_rho_threshold, pattern, hit_direction, database) %>%
  dplyr::summarise(
    n_hit_input = dplyr::first(n_hit_input),
    n_terms_significant = sum(
      !is.na(enrichment_padj_BH) & enrichment_padj_BH < enrichment_fdr_threshold,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(pattern = as.character(pattern)) %>%
  dplyr::arrange(
    delta_rho_threshold,
    factor(pattern, levels = pattern_group_order),
    hit_direction,
    database
  )


## 10.12 每 track 每 db 统计 ####
per_pattern_db_summary_ora_df <- build_per_pattern_database_summary_ora(
  ora_long_df = ora_long_df,
  threshold_target = delta_rho_threshold_main,
  fdr_threshold = enrichment_fdr_threshold
) %>%
  dplyr::mutate(pattern = as.character(pattern))


# 11. 结果导出 (Excel, 全中文 sheet 名) ####

log_message(paste0("写出 7.7.1 富集分析 Excel: ", output_xlsx))

## 11.1 参数表 ####
parameter_summary_df <- tibble::tibble(
  参数 = c(
    "msigdb_species",
    "msigdb_databases",
    "enrichment_min_gs_size",
    "enrichment_max_gs_size",
    "enrichment_pvalue_cutoff",
    "enrichment_qvalue_cutoff",
    "enrichment_fdr_threshold",
    "top_n_pathway_per_pattern",
    "min_hit_genes_for_enrichment",
    "min_n_patterns_for_shared",
    "delta_rho_threshold_main",
    "delta_rho_thresholds_evaluated",
    "hit_direction_levels",
    "hit_track_levels",
    "gsea_n_perm",
    "gsea_eps",
    "gsea_seed",
    "gsea_ranking_metric",
    "pattern_group_order"
  ),
  值 = c(
    msigdb_species,
    paste(names(msigdb_database_config), collapse = ", "),
    as.character(enrichment_min_gs_size),
    as.character(enrichment_max_gs_size),
    as.character(enrichment_pvalue_cutoff),
    as.character(enrichment_qvalue_cutoff),
    as.character(enrichment_fdr_threshold),
    as.character(top_n_pathway_per_pattern),
    as.character(min_hit_genes_for_enrichment),
    as.character(min_n_patterns_for_shared),
    as.character(delta_rho_threshold_main),
    paste(thresholds_to_run, collapse = ", "),
    paste(hit_direction_levels, collapse = ", "),
    paste(hit_track_levels, collapse = ", "),
    as.character(gsea_n_perm),
    as.character(gsea_eps),
    as.character(gsea_seed),
    "signed delta_rho (per pattern)",
    paste(pattern_group_order, collapse = ", ")
  )
)

## 11.2 positive control keyword 表 ####
positive_control_keyword_df <- dplyr::bind_rows(lapply(
  names(positive_control_keywords),
  function(p) {
    data.frame(
      pattern = p,
      keyword = positive_control_keywords[[p]],
      stringsAsFactors = FALSE
    )
  }
)) %>%
  dplyr::arrange(
    factor(pattern, levels = pattern_group_order),
    keyword
  )

## 11.3 Excel sheet 列表 ####
output_sheet_list <- list(
  "依赖检查" = dependency_check_df,
  "分析参数" = parameter_summary_df,
  "数据库统计" = database_summary_df,
  "背景基因汇总" = bg_summary_df,
  "输入hit汇总" = hit_input_summary_df,
  "每模式数据库统计_主阈值" = per_pattern_db_summary_ora_df,

  "ORA结果总表" = ora_long_df %>% dplyr::mutate(pattern = as.character(pattern)),
  "ORA显著通路_主轨" = significant_ora_main_candidate_df %>% dplyr::mutate(pattern = as.character(pattern)),
  "ORA显著通路_HC子集" = significant_ora_main_hc_df %>% dplyr::mutate(pattern = as.character(pattern)),

  "GSEA结果总表" = gsea_long_df %>% dplyr::mutate(pattern = as.character(pattern)),
  "GSEA显著通路_主轨" = significant_gsea_candidate_df %>% dplyr::mutate(pattern = as.character(pattern)),
  "GSEA显著通路_HC子集" = significant_gsea_hc_df %>% dplyr::mutate(pattern = as.character(pattern)),

  "跨模式共享通路_主轨" = shared_pathways_df %>% dplyr::mutate(hit_direction = as.character(hit_direction)),
  "模式特异通路_主轨" = unique_pathways_df %>% dplyr::mutate(pattern = as.character(pattern)),
  "每模式Top通路_主轨" = top_pathways_df %>% dplyr::mutate(pattern = as.character(pattern)),

  "生物学positive_control_关键词" = positive_control_keyword_df,
  "生物学positive_control_匹配通路" = positive_control_term_df %>% dplyr::mutate(pattern = as.character(pattern)),
  "生物学positive_control_检验" = positive_control_test_df %>% dplyr::mutate(pattern = as.character(pattern)),

  "ORA_GSEA一致性" = ora_gsea_agreement_df %>% dplyr::mutate(pattern = as.character(pattern)),
  "阈值敏感性_汇总" = threshold_sensitivity_df
)

## 11.4 每模式独立 sheet (主阈值 candidate ORA 显著通路) ####
for (current_pattern in pattern_group_order) {
  pattern_safe_label <- gsub("-", "_", current_pattern, fixed = TRUE)
  pattern_significant_df <- significant_ora_main_candidate_df %>%
    dplyr::filter(as.character(pattern) == current_pattern) %>%
    dplyr::mutate(pattern = as.character(pattern))

  output_sheet_list[[paste0(pattern_safe_label, "_显著通路_主轨")]] <- pattern_significant_df
}

## sanitize: 截断任何单元字符串超过 Excel 32,767 上限的列 (完整数据已保留在 RData) ####
output_sheet_list_sanitized <- sanitize_sheet_list_for_excel(output_sheet_list, max_chars = 32000L)

writexl::write_xlsx(output_sheet_list_sanitized, path = output_xlsx)

log_message(paste0("Excel 富集结果已写出: ", output_xlsx))
log_message("========== 7.7.1 模式相关基因富集分析全部完成 ==========")
log_message(paste0(
  "主结果文件: ", output_xlsx,
  "; 中间 Rdata 路径: ", rdata_output_path
))
