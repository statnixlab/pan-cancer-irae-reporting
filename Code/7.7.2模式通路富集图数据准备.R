# 7.7.2 模式通路富集图数据准备 ####
# 用途: 为每模式 × 每方向 curate 一组临床表型驱动的代表性通路(覆盖 GO_BP/Hallmark/KEGG/Reactome/C7_IMMUNESIGDB 5 库),
#       从 7.7.1 ORA 结果中提取富集统计量,从 7.7 主筛选命中表中按 |Δρ| 排名选取 top 基因用于图形展示 ####
# 输入:
#   1. Results/7.7.1模式相关基因富集分析/rdata/ora_long_df_all.Rdata
#   2. Results/7.7模式相关基因差分相关识别/tables/模式相关基因差分相关分析结果.xlsx (sheet:主筛选命中列表_主阈值) ####
# 输出: Results/7.7.2模式通路富集图数据准备/tables/模式通路富集图数据.xlsx ####

rm(list = ls())
gc()


# 1. 加载配置 ####
source("Code/config.R")


# 2. 参数与路径 ####
script_name <- "7.7.2模式通路富集图数据准备"
script_results_path <- file.path(results_path, script_name)
table_output_path <- file.path(script_results_path, "tables")
log_output_path <- file.path(script_results_path, "logs")

dir.create(script_results_path, recursive = TRUE, showWarnings = FALSE)
dir.create(table_output_path, recursive = TRUE, showWarnings = FALSE)
dir.create(log_output_path, recursive = TRUE, showWarnings = FALSE)

## 输入文件路径 ####
ora_rdata_file <- file.path(results_path, "7.7.1模式相关基因富集分析", "rdata", "ora_long_df_all.Rdata")
hits_xlsx_file <- file.path(results_path, "7.7模式相关基因差分相关识别", "tables", "模式相关基因差分相关分析结果.xlsx")

## 主要参数 ####
main_track <- "candidate"
main_delta_rho_threshold <- 0.3
main_fdr_threshold <- 0.05
top_n_genes_to_display <- 10  # 每条通路展示前 N 个按 |Δρ| 排名的基因

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


# 4. 定义表型驱动 curate 通路清单 ####
# 每行对应一条图中展示的通路,字段: pattern / hit_direction / database / term_id /
#                                    display_name (图中文本) / clinical_note (表型支点说明) ####
log_message("定义 curate 通路清单...")

curated_pathways_df <- tibble::tribble(
  ~pattern,        ~hit_direction, ~database,        ~term_id,                                                             ~display_name,                          ~clinical_note,

  ## Ocular | Positive (视网膜代谢) - 3 条 ####
  "Ocular",        "positive",     "KEGG",           "KEGG_VALINE_LEUCINE_AND_ISOLEUCINE_DEGRADATION",                     "BCAA degradation",                     "retinal BCAA metabolism",
  "Ocular",        "positive",     "Reactome",       "REACTOME_MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION",                   "Mitochondrial β-oxidation",            "retinal mitochondrial energy",
  "Ocular",        "positive",     "GO_BP",          "GOBP_L_AMINO_ACID_METABOLIC_PROCESS",                                "L-amino acid metabolism",              "retinal amino acid metabolism",

  ## Ocular | Negative (神经/突触下调) - 4 条 ####
  "Ocular",        "negative",     "Reactome",       "REACTOME_NEURONAL_SYSTEM",                                           "Neuronal System",                      "eye as neural extension",
  "Ocular",        "negative",     "Reactome",       "REACTOME_TRANSMISSION_ACROSS_CHEMICAL_SYNAPSES",                     "Transmission across Chemical Synapses","retinal synaptic signaling",
  "Ocular",        "negative",     "Reactome",       "REACTOME_NEUROTRANSMITTER_RECEPTORS_AND_POSTSYNAPTIC_SIGNAL_TRANSMISSION", "Neurotransmitter receptors postsynaptic", "visual signal transduction",
  "Ocular",        "negative",     "GO_BP",          "GOBP_POSITIVE_REGULATION_OF_SYNAPSE_ASSEMBLY",                       "Positive reg. of synapse assembly",    "synapse formation",

  ## Hematologic | Positive (免疫激活/白细胞迁移) - 5 条 ####
  "Hematologic",   "positive",     "Hallmark",       "HALLMARK_TNFA_SIGNALING_VIA_NFKB",                                   "TNFA Signaling via NFKB",              "chronic inflammation",
  "Hematologic",   "positive",     "GO_BP",          "GOBP_LEUKOCYTE_MIGRATION",                                           "Leukocyte migration",                  "direct HEME phenotype",
  "Hematologic",   "positive",     "Hallmark",       "HALLMARK_ESTROGEN_RESPONSE_EARLY",                                   "Estrogen Response Early",              "chronic hormonal adaptation",
  "Hematologic",   "positive",     "C7_IMMUNESIGDB", "GSE19888_CTRL_VS_T_CELL_MEMBRANES_ACT_MAST_CELL_DN",                 "Mast cell: T-cell-stim DN",            "mast cell signature",
  "Hematologic",   "positive",     "C7_IMMUNESIGDB", "GSE46606_IRF4_KO_VS_WT_UNSTIM_BCELL_DN",                            "B cell: IRF4-KO DN",                   "B cell immunity",

  ## Hematologic | Negative (慢性脂代谢消耗) - 4 条 ####
  "Hematologic",   "negative",     "KEGG",           "KEGG_FATTY_ACID_METABOLISM",                                         "Fatty acid metabolism",                "chronic metabolic drain",
  "Hematologic",   "negative",     "Hallmark",       "HALLMARK_ADIPOGENESIS",                                              "Adipogenesis",                         "adipocyte differentiation",
  "Hematologic",   "negative",     "Reactome",       "REACTOME_MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION",                   "Mitochondrial β-oxidation",            "fatty acid oxidation loss",
  "Hematologic",   "negative",     "GO_BP",          "GOBP_LIPID_CATABOLIC_PROCESS",                                       "Lipid catabolic process",              "lipid breakdown",

  ## Neurologic | Positive (髓鞘合成上调) - 6 条 (5 + 1 补 GO_BP 以覆盖该库显著信号) ####
  "Neurologic",    "positive",     "Reactome",       "REACTOME_CHOLESTEROL_BIOSYNTHESIS",                                  "Cholesterol biosynthesis",             "myelin core",
  "Neurologic",    "positive",     "Reactome",       "REACTOME_LANOSTEROL_BIOSYNTHESIS",                                   "Lanosterol biosynthesis",              "myelin precursor",
  "Neurologic",    "positive",     "KEGG",           "KEGG_STEROID_BIOSYNTHESIS",                                          "Steroid biosynthesis",                 "myelin lipid synthesis",
  "Neurologic",    "positive",     "Reactome",       "REACTOME_ACTIVATION_OF_GENE_EXPRESSION_BY_SREBF_SREBP",              "SREBP gene expression",                "lipid biosynthesis regulation",
  "Neurologic",    "positive",     "Hallmark",       "HALLMARK_ANDROGEN_RESPONSE",                                         "Androgen Response",                    "steroid hormone signaling",
  "Neurologic",    "positive",     "GO_BP",          "GOBP_STEROL_BIOSYNTHETIC_PROCESS",                                   "Sterol biosynthetic process",          "myelin biosynthesis (GO)",

  ## Neurologic | Negative (神经退化) - 2 条 ####
  "Neurologic",    "negative",     "Hallmark",       "HALLMARK_MYOGENESIS",                                                "Myogenesis",                           "neuromuscular junction (MG)",
  "Neurologic",    "negative",     "GO_BP",          "GOBP_SENSORY_SYSTEM_DEVELOPMENT",                                    "Sensory system development",           "neural sensory loss",

  ## Systemic-like | Positive (急性炎症五重奏) - 6 条 ####
  "Systemic-like", "positive",     "Hallmark",       "HALLMARK_TNFA_SIGNALING_VIA_NFKB",                                   "TNFA Signaling via NFKB",              "acute systemic inflammation",
  "Systemic-like", "positive",     "Hallmark",       "HALLMARK_HYPOXIA",                                                   "Hypoxia",                              "acute stress response",
  "Systemic-like", "positive",     "Hallmark",       "HALLMARK_P53_PATHWAY",                                               "P53 Pathway",                          "DNA damage response",
  "Systemic-like", "positive",     "Hallmark",       "HALLMARK_APOPTOSIS",                                                 "Apoptosis",                            "cell death activation",
  "Systemic-like", "positive",     "Hallmark",       "HALLMARK_TGF_BETA_SIGNALING",                                        "TGF-β signaling",                      "immune modulation",
  "Systemic-like", "positive",     "C7_IMMUNESIGDB", "GSE37605_FOXP3_FUSION_GFP_VS_IRES_GFP_TREG_C57BL6_UP",              "Treg: FOXP3-fusion signature",         "Treg signature",

  ## Systemic-like | Negative (OXPHOS 抑制 - Warburg-like) - 5 条 ####
  "Systemic-like", "negative",     "KEGG",           "KEGG_OXIDATIVE_PHOSPHORYLATION",                                     "Oxidative Phosphorylation",            "Warburg-like shift",
  "Systemic-like", "negative",     "GO_BP",          "GOBP_OXIDATIVE_PHOSPHORYLATION",                                     "Oxidative Phosphorylation",            "mitochondrial suppression",
  "Systemic-like", "negative",     "Hallmark",       "HALLMARK_OXIDATIVE_PHOSPHORYLATION",                                 "Oxidative Phosphorylation",            "Warburg-like",
  "Systemic-like", "negative",     "Reactome",       "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT",                            "Respiratory electron transport",       "ETC suppression",
  "Systemic-like", "negative",     "KEGG",           "KEGG_PARKINSONS_DISEASE",                                            "Parkinson's disease",                  "mitochondrial disease map",

  ## Visceral-like | Positive (胝胆代谢 + 蛋白合成应激) - 7 条 (6 + 1 补 GO_BP 以覆盖该库显著信号) ####
  "Visceral-like", "positive",     "Hallmark",       "HALLMARK_BILE_ACID_METABOLISM",                                      "Bile acid metabolism",                 "HEPA direct",
  "Visceral-like", "positive",     "Hallmark",       "HALLMARK_UV_RESPONSE_DN",                                            "UV Response Down",                     "chronic stress",
  "Visceral-like", "positive",     "Hallmark",       "HALLMARK_IL2_STAT5_SIGNALING",                                       "IL2-STAT5 signaling",                  "immune activation",
  "Visceral-like", "positive",     "KEGG",           "KEGG_RIBOSOME",                                                      "Ribosome",                             "protein synthesis surge",
  "Visceral-like", "positive",     "Reactome",       "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION",                         "Eukaryotic Translation Elongation",    "global translation up",
  "Visceral-like", "positive",     "C7_IMMUNESIGDB", "GSE42088_UNINF_VS_LEISHMANIA_INF_DC_2H_DN",                          "DC: Leishmania-infected DN",           "dendritic cell signature",
  "Visceral-like", "positive",     "GO_BP",          "GOBP_CYTOPLASMIC_TRANSLATION",                                       "Cytoplasmic translation",              "global translation (GO)",

  ## Visceral-like | Negative (自身免疫 + 染色质调控) - 4 条 (3 + 1 补 GO_BP 以覆盖该库显著信号) ####
  "Visceral-like", "negative",     "KEGG",           "KEGG_SYSTEMIC_LUPUS_ERYTHEMATOSUS",                                  "Systemic lupus erythematosus",         "autoimmune mechanism",
  "Visceral-like", "negative",     "Reactome",       "REACTOME_HCMV_LATE_EVENTS",                                          "HCMV Late Events",                     "chronic viral response",
  "Visceral-like", "negative",     "Reactome",       "REACTOME_SIRT1_NEGATIVELY_REGULATES_RRNA_EXPRESSION",                "SIRT1 neg-regulates rRNA",             "chromatin regulation",
  "Visceral-like", "negative",     "GO_BP",          "GOBP_NUCLEOSOME_ORGANIZATION",                                       "Nucleosome organization",              "chromatin (GO)"
)

log_message(sprintf("curate 通路清单: %d 条", nrow(curated_pathways_df)))
log_message(sprintf("覆盖 %d pattern × dir 组合",
                    curated_pathways_df %>% dplyr::distinct(pattern, hit_direction) %>% nrow()))


# 5. 加载 ORA 富集数据 ####
log_message("加载 ORA 富集长表...")
ora_load_env <- new.env()
loaded_objects <- load(ora_rdata_file, envir = ora_load_env)
ora_long_df <- get(loaded_objects[1], envir = ora_load_env)
rm(ora_load_env)
log_message(sprintf("ORA 长表行数: %d (加载对象名: %s)", nrow(ora_long_df), loaded_objects[1]))


# 6. 加载主筛选命中表 (用于 |Δρ| 排序基因) ####
log_message("加载主筛选命中表...")
main_hits_df <- readxl::read_excel(hits_xlsx_file, sheet = "主筛选命中列表_主阈值") %>%
  dplyr::filter(delta_rho_threshold == main_delta_rho_threshold) %>%
  dplyr::select(pattern, hit_direction, gene, delta_rho)
log_message(sprintf("主筛选命中基因表行数: %d", nrow(main_hits_df)))


# 7. 提取 curate 通路的富集统计 ####
log_message("提取 curate 通路的富集统计...")

ora_candidate_df <- ora_long_df %>%
  dplyr::filter(track == main_track, delta_rho_threshold == main_delta_rho_threshold) %>%
  dplyr::select(
    pattern, hit_direction, database, term_id, term_description,
    n_hit_in_term, n_hit_total_in_universe, n_term_in_universe, n_universe,
    fold_enrichment, enrichment_pvalue, enrichment_padj_BH, enrichment_qvalue,
    term_gene_symbols
  )

curated_enrichment_df <- curated_pathways_df %>%
  dplyr::left_join(
    ora_candidate_df,
    by = c("pattern", "hit_direction", "database", "term_id")
  ) %>%
  dplyr::mutate(
    log10_padj = -log10(pmax(enrichment_padj_BH, 1e-30)),
    is_significant = !is.na(enrichment_padj_BH) & enrichment_padj_BH < main_fdr_threshold,
    term_found = !is.na(enrichment_padj_BH)
  )

## 检查是否有任何 curate 通路未找到对应的 ORA 结果 ####
missing_terms <- curated_enrichment_df %>% dplyr::filter(!term_found)
if (nrow(missing_terms) > 0) {
  log_message(sprintf("⚠ %d 条 curate 通路未在 ORA 结果中找到:", nrow(missing_terms)))
  for (i in seq_len(nrow(missing_terms))) {
    log_message(sprintf("  [%s | %s | %s] %s",
                        missing_terms$pattern[i], missing_terms$hit_direction[i],
                        missing_terms$database[i], missing_terms$term_id[i]))
  }
} else {
  log_message("✓ 所有 curate 通路均成功匹配 ORA 结果")
}


# 8. 为每条通路按 |Δρ| 排序提取 top 基因 ####
log_message(sprintf("为每条通路按 |Δρ| 排序提取 top %d 基因...", top_n_genes_to_display))

## 第一步: 展开 term_gene_symbols 逐行 ####
curated_with_row_id <- curated_enrichment_df %>%
  dplyr::mutate(curated_row_id = dplyr::row_number())

genes_long_df <- curated_with_row_id %>%
  dplyr::filter(!is.na(term_gene_symbols)) %>%
  dplyr::select(curated_row_id, pattern, hit_direction, term_gene_symbols) %>%
  tidyr::separate_rows(term_gene_symbols, sep = "/") %>%
  dplyr::rename(gene = term_gene_symbols) %>%
  dplyr::left_join(main_hits_df, by = c("pattern", "hit_direction", "gene"))

## 第二步: 按 |Δρ| 排名,取 top N ####
genes_top_df <- genes_long_df %>%
  dplyr::group_by(curated_row_id) %>%
  dplyr::arrange(dplyr::desc(abs(delta_rho)), .by_group = TRUE) %>%
  dplyr::mutate(gene_rank = dplyr::row_number()) %>%
  dplyr::ungroup()

gene_top_summary_df <- genes_top_df %>%
  dplyr::group_by(curated_row_id) %>%
  dplyr::summarise(
    n_genes_total = dplyr::n(),
    gene_top_list = paste(gene[gene_rank <= top_n_genes_to_display], collapse = "/"),
    gene_top_delta_rho = paste(sprintf("%.2f", delta_rho[gene_rank <= top_n_genes_to_display]), collapse = "/"),
    .groups = "drop"
  )


# 9. 合并 curate 通路 + top 基因 ####
curated_final_df <- curated_with_row_id %>%
  dplyr::left_join(gene_top_summary_df, by = "curated_row_id") %>%
  dplyr::select(-curated_row_id) %>%
  dplyr::mutate(
    n_genes_total = tidyr::replace_na(n_genes_total, 0L),
    gene_top_list = tidyr::replace_na(gene_top_list, ""),
    gene_top_delta_rho = tidyr::replace_na(gene_top_delta_rho, "")
  ) %>%
  dplyr::select(
    pattern, hit_direction, database, term_id, term_description, display_name, clinical_note,
    n_hit_in_term, n_hit_total_in_universe, n_term_in_universe, n_universe,
    fold_enrichment, enrichment_pvalue, enrichment_padj_BH, enrichment_qvalue,
    log10_padj, is_significant, term_found, n_genes_total, gene_top_list, gene_top_delta_rho
  )

log_message(sprintf("最终 curate 通路富集表行数: %d (期望 %d)",
                    nrow(curated_final_df), nrow(curated_pathways_df)))
log_message(sprintf("显著 (FDR<%.2f) 通路数: %d",
                    main_fdr_threshold, sum(curated_final_df$is_significant, na.rm = TRUE)))


# 10. Pattern × Direction 图形排序 (用于绘图顺序) ####
log_message("生成图形排序表...")

pattern_order_tbl <- tibble::tribble(
  ~pattern,        ~pattern_rank, ~pattern_label,                           ~median_TTO,
  "Systemic-like",  1,            "Systemic-like (TTO 32d)",                32,
  "Visceral-like",  2,            "Visceral-like (TTO 34d)",                34,
  "Neurologic",     3,            "Neurologic (TTO 37d)",                   37,
  "Ocular",         4,            "Ocular (TTO 39d)",                       39,
  "Hematologic",    5,            "Hematologic (TTO 43d)",                  43
)

dir_order_tbl <- tibble::tribble(
  ~hit_direction, ~dir_rank, ~dir_label,
  "positive",     1,         "Positive Δρ (tumor-up correlation)",
  "negative",     2,         "Negative Δρ (tumor-down correlation)"
)

plot_order_df <- tidyr::expand_grid(
  pattern = pattern_order_tbl$pattern,
  hit_direction = dir_order_tbl$hit_direction
) %>%
  dplyr::left_join(pattern_order_tbl, by = "pattern") %>%
  dplyr::left_join(dir_order_tbl, by = "hit_direction") %>%
  dplyr::arrange(pattern_rank, dir_rank) %>%
  dplyr::mutate(
    panel_id = seq_len(dplyr::n()),
    panel_label = sprintf("%s | %s", pattern, hit_direction)
  )


# 11. 数据库顺序 + 颜色映射 ####
log_message("生成数据库顺序与颜色映射...")

database_order_tbl <- tibble::tribble(
  ~database,        ~database_rank, ~database_color_hex, ~database_label,
  "GO_BP",          1,              "#3F88C5",           "BP",
  "Hallmark",       2,              "#D7263D",           "Hallmark",
  "KEGG",           3,              "#F6AE2D",           "KEGG",
  "Reactome",       4,              "#7E1E9C",           "Reactome",
  "C7_IMMUNESIGDB", 5,              "#1B998B",           "C7-IS"
)


# 12. 写出 Excel ####
log_message("写出 Excel...")

## 函数: truncate_excel_strings ####
## 功能: 将字符串列截断到 Excel 单元上限 (32767 字符) ####
## 输入: df,数据框 ####
## 输出: 截断后数据框 ####
truncate_excel_strings <- function(df) {
  df %>% dplyr::mutate(dplyr::across(
    dplyr::where(is.character),
    ~ ifelse(nchar(.x) > 32700, paste0(substr(.x, 1, 32700), "..."), .x)
  ))
}

output_sheet_list <- list(
  "通路富集长表" = truncate_excel_strings(curated_final_df),
  "通路清单" = truncate_excel_strings(curated_pathways_df),
  "图形排序" = plot_order_df,
  "数据库顺序颜色" = database_order_tbl
)

output_xlsx_path <- file.path(table_output_path, "模式通路富集图数据.xlsx")
writexl::write_xlsx(output_sheet_list, path = output_xlsx_path)
log_message(sprintf("Excel 已保存: %s", output_xlsx_path))

log_message(paste0("========== ", script_name, " 完成 =========="))
log_message(sprintf("输出: %s", output_xlsx_path))
