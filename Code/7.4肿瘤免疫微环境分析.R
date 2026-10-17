# 7.4 肿瘤免疫微环境分析 ####
# 用途：基于TCGA肿瘤组织表达矩阵计算整体微环境组成、免疫细胞浸润、免疫激活状态和checkpoint相关特征。
#       所有样本级评分先保存为结果表，再按TCGA project和primary_organ汇总。
#       主统计检验基于primary_organ级结果，避免将样本级差异误判为FAERS模式群差异。 ####

rm(list = ls())
gc()


# 1. 加载配置 ####
source("Code/config.R")


# 2. 参数与路径 ####
script_name <- "7.4肿瘤免疫微环境分析"
script_results_path <- file.path(results_path, script_name)
table_output_path <- file.path(script_results_path, "tables")
rdata_output_path <- file.path(script_results_path, "rdata")
log_output_path <- file.path(script_results_path, "logs")

dir.create(script_results_path, recursive = TRUE, showWarnings = FALSE)
dir.create(table_output_path, recursive = TRUE, showWarnings = FALSE)
dir.create(rdata_output_path, recursive = TRUE, showWarnings = FALSE)
dir.create(log_output_path, recursive = TRUE, showWarnings = FALSE)

input_data_path <- file.path(results_path, "7.1TCGA基础数据准备")
input_rdata_path <- file.path(input_data_path, "rdata")
input_unified_path <- file.path(input_rdata_path, "unified")

## 直接使用 7.1 生成的原始 TCGA 肿瘤表达矩阵，避免 unified inner-join 丢失 SOD2/SIGLEC5/VNN3 等 4796 个 TCGA 原生基因 ####
tcga_tumor_expression_rdata <- file.path(
  input_rdata_path,
  "tcga_tumor_expression_matrix.Rdata"
)

tcga_tumor_sample_annotation_rdata <- file.path(
  input_rdata_path,
  "tcga_tumor_sample_annotation.Rdata"
)

tcga_maf_dir <- file.path(data_path, "TCGA_GTEx", "by_project")
tcga_panimmune_tcr_file <- file.path(
  data_path,
  "TCGA_GTEx",
  "panimmune",
  "mitcr_sampleStatistics_20160714.tsv"
)

## 临床数据路径（用于CEM匹配LMM补充分析） ####
clinical_data_dir <- file.path(data_path, "TCGA_GTEx", "by_project")

output_xlsx <- file.path(table_output_path, "TCGA肿瘤免疫微环境分析结果.xlsx")
sample_score_rdata <- file.path(rdata_output_path, "sample_level_tme_score_df.Rdata")
project_score_rdata <- file.path(rdata_output_path, "tcga_project_level_tme_score_df.Rdata")
primary_organ_score_rdata <- file.path(rdata_output_path, "primary_organ_level_tme_score_df.Rdata")
bootstrap_result_rdata <- file.path(rdata_output_path, "bootstrap_tme_stability_result_list.Rdata")

log2_offset <- 1
keep_first_tumor_per_patient <- TRUE
min_genes_default <- 2
min_genes_bindea <- 1
exome_size_mb_for_tmb <- 38
bootstrap_n_iterations <- 1000
bootstrap_seed <- 20260427
bootstrap_samples_per_primary_organ <- 30
bootstrap_min_samples_per_primary_organ <- 10
bootstrap_primary_organs_per_pattern_group <- NA_integer_
subsample_lmm_n_iterations <- 200
subsample_lmm_n_per_organ <- 35
subsample_lmm_seed <- 20260512

## CEM匹配LMM参数：每organ抽样35、年龄按10年分箱（TumorPurity不再用作协变量，因在RNA-seq上不可靠且已从特征中移除） ####
cem_lmm_n_per_organ <- 35
cem_lmm_seed <- 20260513
cem_age_bin_breaks <- c(-Inf, 40, 50, 60, 70, 80, Inf)
detected_cores <- parallel::detectCores(logical = FALSE)
if (is.na(detected_cores)) {
  detected_cores <- 1
}
parallel_workers <- as.integer(max(1, min(8, detected_cores - 1)))
enable_parallel <- parallel_workers > 1

pattern_group_order <- c(
  "Ocular",
  "Hematologic",
  "Neurologic",
  "Systemic-like",
  "Visceral-like"
)

required_packages <- c(
  "tidyestimate",
  "msigdbr",
  "GSVA",
  "BiocParallel",
  "GSEABase",
  "matrixStats",
  "openxlsx",
  "data.table",
  "dplyr",
  "tidyr",
  "tibble",
  "purrr",
  "lme4",
  "lmerTest",
  "emmeans"
)

run_timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_file <- file.path(log_output_path, paste0("7.4肿瘤免疫微环境分析_", run_timestamp, ".log"))


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

## 函数：make_safe_feature_name ####
## 功能：将signature名称转换为适合作为列名的稳定英文名称 ####
## 输入：feature_name，原始特征名称 ####
## 输出：安全列名字符串 ####
make_safe_feature_name <- function(feature_name) {
  safe_name <- gsub("[^A-Za-z0-9]+", "_", feature_name)
  safe_name <- gsub("^_+|_+$", "", safe_name)
  safe_name <- gsub("_+", "_", safe_name)
  safe_name
}

## 函数：collapse_values ####
## 功能：将向量去重后折叠为分号分隔字符串 ####
## 输入：values，待折叠向量 ####
## 输出：折叠后的字符串 ####
collapse_values <- function(values) {
  values <- unique(as.character(values))
  values <- values[!is.na(values) & values != ""]

  if (length(values) == 0) {
    return(NA_character_)
  }

  paste(values, collapse = "; ")
}

## 函数：check_required_packages ####
## 功能：检查7.4分析所需R包是否安装，不自动安装缺失包 ####
## 输入：packages，包名向量；output_path，依赖检查结果保存路径 ####
## 输出：依赖检查数据框；若有缺失则终止脚本 ####
check_required_packages <- function(packages, output_path) {
  dependency_check_df <- data.frame(
    package = packages,
    installed = vapply(packages, requireNamespace, logical(1), quietly = TRUE),
    install_hint = dplyr::case_when(
      packages %in% c("tidyestimate", "msigdbr") ~ paste0("install.packages(\"", packages, "\")"),
      packages %in% c("GSVA", "GSEABase", "BiocParallel") ~ paste0("BiocManager::install(\"", packages, "\")"),
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
      missing_df$package,
      ": ",
      missing_df$install_hint,
      collapse = "\n"
    )
    stop("缺少7.4分析必要R包，请手动安装后重新运行:\n", missing_message)
  }

  dependency_check_df
}

## 函数：check_required_columns ####
## 功能：检查数据框是否包含必要列 ####
## 输入：data_df，待检查数据框；required_columns，必要列名向量；data_label，数据说明 ####
## 输出：无返回值 ####
check_required_columns <- function(data_df, required_columns, data_label) {
  missing_columns <- setdiff(required_columns, colnames(data_df))

  if (length(missing_columns) > 0) {
    stop(data_label, " 缺少必要列: ", paste(missing_columns, collapse = ", "))
  }
}

## 函数：make_biocparallel_param ####
## 功能：根据当前平台和并行参数生成BiocParallel后端 ####
## 输入：use_parallel，是否启用并行；workers，并行进程数 ####
## 输出：BiocParallel参数对象 ####
make_biocparallel_param <- function(use_parallel, workers) {
  if (!use_parallel || workers < 2) {
    return(BiocParallel::SerialParam(progressbar = FALSE))
  }

  if (.Platform$OS.type == "windows") {
    BiocParallel::SnowParam(workers = workers, type = "SOCK", progressbar = FALSE)
  } else {
    BiocParallel::MulticoreParam(workers = workers, progressbar = FALSE)
  }
}

## 函数：run_parallel_lapply ####
## 功能：在Linux/macOS中用mclapply并行执行列表任务，无法并行时自动退回lapply ####
## 输入：x，任务输入向量或列表；fun，任务函数；use_parallel，是否启用并行；workers，并行进程数；...，传递给任务函数的其他参数 ####
## 输出：与lapply一致的列表结果 ####
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

## 函数：run_gsva_with_parallel ####
## 功能：兼容不同GSVA版本调用gsva，并在可用时启用BiocParallel或parallel.sz并行 ####
## 输入：gsva_args，传给GSVA::gsva的参数列表；use_parallel，是否启用并行；workers，并行进程数 ####
## 输出：GSVA评分矩阵 ####
run_gsva_with_parallel <- function(gsva_args, use_parallel, workers) {
  if (!use_parallel || workers < 2) {
    return(do.call(GSVA::gsva, gsva_args))
  }

  bpparam <- make_biocparallel_param(use_parallel = use_parallel, workers = workers)
  gsva_args_with_bpparam <- c(gsva_args, list(BPPARAM = bpparam))

  tryCatch(
    do.call(GSVA::gsva, gsva_args_with_bpparam),
    error = function(e) {
      if (!grepl("unused argument|formal argument", conditionMessage(e))) {
        stop(e)
      }

      gsva_args_with_parallel_size <- c(gsva_args, list(parallel.sz = workers))

      tryCatch(
        do.call(GSVA::gsva, gsva_args_with_parallel_size),
        error = function(e2) {
          if (!grepl("unused argument|formal argument", conditionMessage(e2))) {
            stop(e2)
          }

          do.call(GSVA::gsva, gsva_args)
        }
      )
    }
  )
}


# 4. Signature构建函数 ####

## 函数：build_bindea_signature_list ####
## 功能：构建Bindea等报道的24类免疫细胞signature基因集 ####
## 输入：无 ####
## 输出：命名列表，每个元素为一个免疫细胞signature基因向量 ####
build_bindea_signature_list <- function() {
  list(
    "aDC" = c("CCL1", "EBI3", "INDO", "LAMP3", "OAS3"),
    "B cells" = c("ABCB4", "BACH2", "BCL11A", "BLK", "BLNK", "CCR9", "CD19", "CD72", "COCH", "CR2", "DTNB", "FCRL2", "GLDC", "GNG7", "HLA-DOB", "HLA-DQA1", "IGHA1", "IGHG1", "IGHM", "IGKC", "IGL", "KIAA0125", "MEF2C", "MICAL3", "MS4A1", "OSBPL10", "PNOC", "QRSL1", "SCN3A", "SLC15A2", "SPIB", "TCL1A", "TNFRSF17"),
    "CD8 T cells" = c("ABT1", "AES", "APBA2", "ARHGAP8", "C12orf47", "C19orf6", "C4orf15", "CAMLG", "CD8A", "CD8B", "CDKN2AIP", "DNAJB1", "FLT3LG", "GADD45A", "GZMM", "KLF9", "LEPROTL1", "LIME1", "MYST3", "PF4", "PPP1R2", "PRF1", "PRR5", "RBM3", "SF1", "SFRS7", "SLC16A7", "TBCC", "THUMPD1", "TMC6", "TSC22D3", "VAMP2", "ZEB1", "ZFP36L2", "ZNF22", "ZNF609", "ZNF91"),
    "Cytotoxic cells" = c("APBA2", "APOL3", "CTSW", "DUSP2", "GNLY", "GZMA", "GZMH", "KLRB1", "KLRD1", "KLRF1", "KLRK1", "NKG7", "RORA", "RUNX3", "SIGIRR", "WHAMMP3", "ZBTB16"),
    "DC" = c("CCL13", "CCL17", "CCL22", "CD209", "HSD11B1", "NPR1", "PPFIBP2"),
    "Eosinophils" = c("ABHD2", "ACACB", "C9orf156", "CAT", "CCR3", "CLC", "CYSLTR2", "EMR1", "EPN2", "GALC", "GPR44", "HES1", "HIST1H1C", "HRH4", "IGSF2", "IL5RA", "KBTBD11", "KCNH2", "LRP5L", "MYO15B", "RCOR3", "RNASE2", "RRP12", "SIAH1", "SMPD3", "SYNJ1", "TGIF1", "THBS1", "THBS4", "TIPARP", "TKTL1"),
    "iDC" = c("ABCG2", "BLVRB", "CARD9", "CD1A", "CD1B", "CD1C", "CD1E", "CH25H", "CLEC10A", "CSF1R", "CTNS", "F13A1", "FABP4", "FZD2", "GSTT1", "GUCA1A", "HS3ST2", "LMAN2L", "MMP12", "MS4A6A", "NUDT9", "PDXK", "PPARG", "PREP", "RAP1GAP", "SLC26A6", "SLC7A8", "SYT17", "TACSTD2", "TM7SF4", "VASH1"),
    "Macrophages" = c("APOE", "ATG7", "BCAT1", "CCL7", "CD163", "CD68", "CD84", "CHI3L1", "CHIT1", "CLEC5A", "COL8A2", "COLEC12", "CTSK", "CXCL5", "CYBB", "DNASE2B", "EMP1", "FDX1", "FN1", "GM2A", "GPC4", "KAL1", "MARCO", "ME1", "MS4A4A", "MSR1", "PCOLCE2", "PTGDS", "RAI14", "SCARB2", "SCG5", "SGMS1", "SULT1C2"),
    "Mast cells" = c("ABCC4", "ADCYAP1", "CALB2", "CEACAM8", "CMA1", "CPA3", "CTSG", "ELA2", "GATA2", "HDC", "HPGD", "HPGDS", "KIT", "LINC01140", "MAOB", "MLPH", "MPO", "MS4A2", "NR0B1", "PPM1H", "PRG2", "PTGS1", "SCG2", "SIGLEC6", "SLC18A2", "SLC24A3", "TAL1", "TPSAB1", "TPSB2", "VWA5A"),
    "Neutrophils" = c("ALPL", "BST1", "CD93", "CEACAM3", "CREB5", "CRISPLD2", "CSF3R", "CYP4F3", "DYSF", "FCAR", "FCGR3B", "FLJ11151", "FPR1", "FPRL1", "G0S2", "HIST1H2BC", "HPSE", "IL8RA", "IL8RB", "KCNJ15", "LILRB2", "MGAM", "MME", "PDE4B", "S100A12", "SIGLEC5", "SLC22A4", "SLC25A37", "TECPR2", "TNFRSF10C", "VNN3"),
    "NK CD56bright cells" = c("DUSP4", "FOXJ1", "LPCAT4", "MADD", "6-Mar", "MPPED1", "MUC3B", "NIBP", "PLA2G6", "RRAD", "XCL1"),
    "NK CD56dim cells" = c("EDG8", "FLJ20699", "GTF3C1", "GZMB", "IL21R", "KIR2DL3", "KIR2DS1", "KIR2DS2", "KIR2DS5", "KIR3DL1", "KIR3DL2", "KIR3DL3", "KIR3DS1", "SPON2", "TMEPAI"),
    "NK cells" = c("ADARB1", "AF107846", "ALDH1B1", "APBB2", "ATL2", "BCL2", "CDC5L", "FGF18", "FUT5", "FZR1", "GAGE2A", "IGFBP5", "KANK2", "LDB3", "MAPRE3", "MCM3AP", "MRC2", "NCR1", "PDLIM4", "PRX", "PSMD4", "RP5-886K2.1", "SGMS1", "SLC30A5", "SMEK1", "SPN", "TBXA2R", "TCTN2", "TINAGL1", "TRPV6", "XCL1", "XCL2", "ZNF205", "ZNF528", "ZNF747"),
    "pDC" = c("IL3RA"),
    "T cells" = c("BCL11B", "CD2", "CD28", "CD3D", "CD3E", "CD3G", "CD6", "CD96", "GIMAP5", "ITM2A", "LCK", "NCALD", "PRKCQ", "SH2D1A", "SKAP1", "TRA", "TRAC", "TRAT1", "TRBC1"),
    "T helper cells" = c("ANP32B", "ASF1A", "ATF2", "BATF", "C13orf34", "CD28", "DDX50", "FAM111A", "FRYL", "GOLGA8A", "ICOS", "ITM2A", "LRBA", "NAP1L4", "NUP107", "PHF10", "PPP2R5C", "RPA1", "SEC24C", "SLC25A12", "SRSF10", "TRA", "UBE2L3", "YME1L1"),
    "Tcm cells" = c("AQP3", "ATF7IP", "ATM", "CASP8", "CDC14A", "CEP68", "CLUAP1", "CREBZF", "CYLD", "DOCK9", "FAM153B", "FOXP1", "FYB", "HNRPH1", "INPP4B", "KLF12", "LOC441155", "MAP3K1", "MLL", "N4BP2L2-IT2", "NEFL", "NFATC3", "PCM1", "PCNX", "PDXDC2", "PHC3", "POLR2J2", "PSPC1", "REPS1", "RPP38", "SLC7A6", "SNRPN", "ST3GAL1", "STX16", "TIMM8A", "TRA", "TRAF3IP3", "TXK", "TXLNGY", "USP9Y"),
    "Tem cells" = c("AKT3", "C7orf54", "CCR2", "DDX17", "EWSR1", "FLI1", "GDPD5", "LTK", "MEFV", "NFATC4", "PRKY", "TBC1D5", "TBCD", "TRA", "VIL2"),
    "Tfh cells" = c("B3GAT1", "BLR1", "C18orf1", "CDK5R1", "CHGB", "CHI3L2", "CXCL13", "HEY1", "HIST1H4K", "ICA1", "KCNK5", "KIAA1324", "MAF", "MAGEH1", "MKL2", "MYO6", "MYO7A", "PASK", "PDCD1", "POMT1", "PTPN13", "PVALB", "SH3TC1", "SIRPG", "SLC7A10", "SMAD1", "ST8SIA1", "STK39", "THADA", "TOX", "TSHR", "ZNF764"),
    "Tgd cells" = c("C1orf61", "CD160", "FEZ1", "TARP", "TRD", "TRGV9"),
    "Th1 cells" = c("APBB2", "APOD", "ATP9A", "BST2", "BTG3", "CCL4", "CD38", "CD70", "CMAH", "CSF2", "CTLA4", "DGKI", "DOK5", "DPP4", "DUSP5", "EGFL6", "GGT1", "HBEGF", "IFNG", "IL12RB2", "IL22", "LRP8", "LRRN3", "LTA", "SGCB", "SYNGR3", "ZBTB32"),
    "Th17 cells" = c("IL17A", "IL17RA", "RORC"),
    "Th2 cells" = c("ADCY1", "AHI1", "ANK1", "BIRC5", "CDC25C", "CDC7", "CENPF", "CXCR6", "DHFR", "EVI5", "GATA3", "GSTA4", "HELLS", "IL26", "LAIR2", "LIMA1", "MB", "MICAL2", "NEIL3", "PHEX", "PMCH", "PTGIS", "SLC39A14", "SMAD2", "SNRPD1", "WDHD1"),
    "Treg cells" = c("FOXP3")
  )
}

## 函数：build_custom_signature_list ####
## 功能：构建小基因集signature，用于mean z-score计算 ####
## 输入：无 ####
## 输出：命名列表，每个元素为一个signature基因向量 ####
build_custom_signature_list <- function() {
  list(
    "T_cell_inflamed_GEP" = c(
      "CCL5", "CD27", "CD274", "CD276", "CD8A", "CMKLR1", "CXCL9", "CXCR6",
      "HLA-DQA1", "HLA-DRB1", "HLA-E", "IDO1", "LAG3", "NKG7",
      "PDCD1LG2", "PSMB10", "STAT1", "TIGIT"
    ),
    "Antigen_Presentation" = c(
      "HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "TAPBP",
      "PSMB8", "PSMB9", "PSMB10", "NLRC5"
    )
  )
}

## 函数：build_feature_info_table ####
## 功能：构建按生物解释轴分组的主分析特征信息表 ####
## 输入：无 ####
## 输出：特征信息数据框 ####
build_feature_info_table <- function() {
  tibble::tribble(
    ~feature_name, ~feature_label, ~feature_axis, ~axis_order, ~feature_order, ~feature_class, ~method, ~source, ~main_display,
    "TMB", "TMB", "Tumor Antigenicity", 1, 1, "mutation", "Nonsynonymous mutations per 38 Mb", "Ellrott et al., 2018; Chalmers et al., 2017", TRUE,
    "TCR_Shannon", "TCR Shannon diversity", "TCR Repertoire", 2, 1, "TCR", "Thorsson/TCGA PanImmune TCR Shannon diversity", "Thorsson et al., 2018", TRUE,
    "Hallmark_HALLMARK_INTERFERON_ALPHA_RESPONSE", "IFN-alpha response", "Cytokine and Inflammatory Signaling", 3, 1, "pathway", "MSigDB Hallmark + ssGSEA", "Liberzon et al., 2015; Hanzelmann et al., 2013", TRUE,
    "Hallmark_HALLMARK_TNFA_SIGNALING_VIA_NFKB", "TNF-alpha signaling via NF-kB", "Cytokine and Inflammatory Signaling", 3, 2, "pathway", "MSigDB Hallmark + ssGSEA", "Liberzon et al., 2015; Hanzelmann et al., 2013", TRUE,
    "CYT_score", "Cytolytic activity (GZMA/PRF1)", "T Cell Activation and Cytotoxicity", 4, 1, "gene_signature", "sqrt(GZMA_TPM * PRF1_TPM) (Rooney 2015 geometric mean on TPM scale)", "Rooney et al., 2015", TRUE,
    "T_cell_inflamed_GEP", "T cell-inflamed GEP", "T Cell Activation and Cytotoxicity", 4, 2, "gene_signature", "Mean z-score", "Ayers et al., 2017; Cristescu et al., 2018", TRUE,
    "Hallmark_HALLMARK_INTERFERON_GAMMA_RESPONSE", "IFN-gamma response", "T Cell Activation and Cytotoxicity", 4, 3, "pathway", "MSigDB Hallmark + ssGSEA", "Liberzon et al., 2015; Hanzelmann et al., 2013", TRUE,
    "Antigen_Presentation", "Antigen presentation", "Antigen Presentation Machinery", 5, 1, "gene_signature", "Mean z-score", "Thorsson et al., 2018; Charoentong et al., 2017", TRUE,
    "Bindea_CD8_T_cells", "CD8 T cells", "Key Immune Cell Infiltration", 6, 1, "cell_infiltration", "Bindea signatures + ssGSEA", "Bindea et al., 2013; Hanzelmann et al., 2013", TRUE,
    "Bindea_Cytotoxic_cells", "Cytotoxic cells", "Key Immune Cell Infiltration", 6, 2, "cell_infiltration", "Bindea signatures + ssGSEA", "Bindea et al., 2013; Hanzelmann et al., 2013", TRUE,
    "Bindea_Treg_cells", "Treg cells", "Key Immune Cell Infiltration", 6, 3, "cell_infiltration", "Bindea signatures (FOXP3-only, original Bindea 2013) + ssGSEA", "Bindea et al., 2013; Hanzelmann et al., 2013", TRUE,
    "Bindea_Macrophages", "Macrophages", "Key Immune Cell Infiltration", 6, 4, "cell_infiltration", "Bindea signatures + ssGSEA", "Bindea et al., 2013; Hanzelmann et al., 2013", TRUE,
    "Bindea_B_cells", "B cells", "Key Immune Cell Infiltration", 6, 5, "cell_infiltration", "Bindea signatures + ssGSEA", "Bindea et al., 2013; Hanzelmann et al., 2013", TRUE,
    "Bindea_NK_cells", "NK cells", "Key Immune Cell Infiltration", 6, 6, "cell_infiltration", "Bindea signatures + ssGSEA", "Bindea et al., 2013; Hanzelmann et al., 2013", TRUE,
    "Bindea_DC", "Dendritic cells", "Key Immune Cell Infiltration", 6, 7, "cell_infiltration", "Bindea signatures + ssGSEA", "Bindea et al., 2013; Hanzelmann et al., 2013", TRUE,
    "Bindea_Neutrophils", "Neutrophils", "Key Immune Cell Infiltration", 6, 8, "cell_infiltration", "Bindea signatures + ssGSEA", "Bindea et al., 2013; Hanzelmann et al., 2013", TRUE,
    "PDCD1_expression", "PD-1 (PDCD1)", "Checkpoint Gene Expression", 7, 1, "gene_expression", "log2(TPM + 1)", "Pardoll, 2012", TRUE,
    "CD274_expression", "PD-L1 (CD274)", "Checkpoint Gene Expression", 7, 2, "gene_expression", "log2(TPM + 1)", "Pardoll, 2012", TRUE,
    "CTLA4_expression", "CTLA4", "Checkpoint Gene Expression", 7, 3, "gene_expression", "log2(TPM + 1)", "Pardoll, 2012", TRUE,
    "ImmuneScore", "ImmuneScore", "TME Context", 8, 1, "TME_context", "tidyestimate::estimate_score (is_affymetrix=FALSE)", "Yoshihara et al., 2013; tidyestimate::estimate_score", TRUE,
    "StromalScore", "StromalScore", "TME Context", 8, 2, "TME_context", "tidyestimate::estimate_score (is_affymetrix=FALSE)", "Yoshihara et al., 2013; tidyestimate::estimate_score", TRUE,
    "ESTIMATEScore", "ESTIMATEScore", "TME Context", 8, 3, "TME_context", "ImmuneScore + StromalScore", "Yoshihara et al., 2013; tidyestimate::estimate_score", TRUE
  ) %>%
    dplyr::mutate(module = feature_axis)
}

## 函数：build_gene_alias_map ####
## 功能：提供部分历史基因符号到当前常用gene symbol的转换，用于提高signature基因匹配率 ####
## 输入：无 ####
## 输出：命名字符向量，名称为历史符号，值为替代符号 ####
build_gene_alias_map <- function() {
  c(
    ## 原有 11 条（已删除错误映射 GAGE2A->GAGE2C，因为 TPM 矩阵中本来就保留着 GAGE2A） ####
    "INDO" = "IDO1",
    "EMR1" = "ADGRE1",
    "FPRL1" = "FPR2",
    "IL8RA" = "CXCR1",
    "IL8RB" = "CXCR2",
    "BLR1" = "CXCR5",
    "EDG8" = "S1PR5",
    "TMEPAI" = "PMEPA1",
    "SFRS7" = "SRSF7",
    "VIL2" = "EZR",
    "MYST3" = "KAT6A",
    ## 新增 12 条：覆盖 Hallmark 与 Bindea signature 中已更名但未映射的基因 ####
    "RIGI" = "DDX58",
    "TMT1B" = "METTL7B",
    "AES" = "TLE5",
    "KAL1" = "ANOS1",
    "KIAA0125" = "FAM30A",
    "SMEK1" = "PPP4R3A",
    "HIST1H2BC" = "H2BC4",
    "C19orf6" = "MEX3D",
    "C12orf47" = "NEMP1",
    "FLJ11151" = "TMEM268",
    "6-Mar" = "MARCHF6",
    "MLL" = "KMT2A"
  )
}

## 函数：apply_gene_aliases ####
## 功能：根据表达矩阵可用基因符号对signature基因进行别名替换 ####
## 输入：genes，signature基因向量；available_genes，表达矩阵行名；alias_map，别名映射向量 ####
## 输出：去重后的基因向量 ####
apply_gene_aliases <- function(genes, available_genes, alias_map) {
  genes <- unique(as.character(genes))
  genes <- genes[!is.na(genes) & genes != ""]

  replaced_genes <- vapply(
    genes,
    function(gene_symbol) {
      if (gene_symbol %in% available_genes) {
        return(gene_symbol)
      }

      if (gene_symbol %in% names(alias_map) && alias_map[[gene_symbol]] %in% available_genes) {
        return(alias_map[[gene_symbol]])
      }

      gene_symbol
    },
    character(1)
  )

  unique(replaced_genes)
}

## 函数：get_hallmark_signature_list ####
## 功能：从msigdbr读取本分析需要的Hallmark通路基因集 ####
## 输入：target_gene_sets，目标Hallmark基因集名称向量 ####
## 输出：命名列表，每个元素为一个Hallmark基因集 ####
get_hallmark_signature_list <- function(target_gene_sets) {
  hallmark_df <- tryCatch(
    msigdbr::msigdbr(species = "Homo sapiens", collection = "H"),
    error = function(e) {
      msigdbr::msigdbr(species = "Homo sapiens", category = "H")
    }
  )

  hallmark_df %>%
    dplyr::filter(gs_name %in% target_gene_sets) %>%
    dplyr::select(gs_name, gene_symbol) %>%
    dplyr::distinct() %>%
    dplyr::group_by(gs_name) %>%
    dplyr::summarise(genes = list(unique(gene_symbol)), .groups = "drop") %>%
    tibble::deframe()
}

## 函数：format_signature_availability ####
## 功能：统计signature基因在表达矩阵中的可用性 ####
## 输入：signature_list，signature基因集列表；available_genes，表达矩阵行名；module_name，模块名；method_name，计算方法；min_genes，最低可用基因数 ####
## 输出：signature可用性数据框 ####
format_signature_availability <- function(signature_list, available_genes, module_name, method_name, min_genes) {
  availability_list <- lapply(names(signature_list), function(signature_name) {
    total_genes <- unique(as.character(signature_list[[signature_name]]))
    total_genes <- total_genes[!is.na(total_genes) & total_genes != ""]
    available_signature_genes <- intersect(total_genes, available_genes)
    missing_signature_genes <- setdiff(total_genes, available_genes)

    data.frame(
      signature_name = signature_name,
      feature_name = signature_name,
      module = module_name,
      method = method_name,
      n_total_genes = length(total_genes),
      n_available_genes = length(available_signature_genes),
      min_required_genes = min_genes,
      enough_genes = length(available_signature_genes) >= min_genes,
      available_genes = collapse_values(available_signature_genes),
      missing_genes = collapse_values(missing_signature_genes),
      stringsAsFactors = FALSE
    )
  })

  dplyr::bind_rows(availability_list)
}


# 5. 评分函数 ####

## 函数：run_ssgsea_scores ####
## 功能：使用GSVA包的ssGSEA方法计算每个样本的signature enrichment score ####
## 输入：expression_matrix，行为基因、列为样本的log2表达矩阵；gene_sets，signature列表；min_genes，最低可用基因数；use_parallel，是否启用并行；workers，并行进程数 ####
## 输出：signature为行、样本为列的ssGSEA分数矩阵 ####
run_ssgsea_scores <- function(expression_matrix, gene_sets, min_genes, use_parallel = FALSE, workers = 1) {
  available_genes <- rownames(expression_matrix)

  filtered_gene_sets <- lapply(gene_sets, function(genes) {
    intersect(unique(as.character(genes)), available_genes)
  })

  filtered_gene_sets <- filtered_gene_sets[lengths(filtered_gene_sets) >= min_genes]

  if (length(filtered_gene_sets) == 0) {
    stop("没有signature达到ssGSEA最低可用基因数要求。")
  }

  if (exists("ssgseaParam", envir = asNamespace("GSVA"), inherits = FALSE)) {
    ssgsea_param <- GSVA::ssgseaParam(
      exprData = expression_matrix,
      geneSets = filtered_gene_sets,
      minSize = min_genes,
      maxSize = Inf,
      normalize = TRUE
    )
    score_matrix <- run_gsva_with_parallel(
      gsva_args = list(ssgsea_param, verbose = FALSE),
      use_parallel = use_parallel,
      workers = workers
    )
  } else {
    score_matrix <- run_gsva_with_parallel(
      gsva_args = list(
        expr = expression_matrix,
        gset.idx.list = filtered_gene_sets,
        method = "ssgsea",
        kcdf = "Gaussian",
        min.sz = min_genes,
        max.sz = Inf,
        ssgsea.norm = TRUE,
        verbose = FALSE
      ),
      use_parallel = use_parallel,
      workers = workers
    )
  }

  as.matrix(score_matrix)
}

## 函数：transpose_score_matrix ####
## 功能：将signature×sample矩阵转换为sample级数据框并添加前缀 ####
## 输入：score_matrix，signature为行、样本为列的矩阵；feature_prefix，列名前缀；label_prefix，展示标签前缀 ####
## 输出：列表，包含样本级评分数据框和特征信息数据框 ####
transpose_score_matrix <- function(score_matrix, feature_prefix, label_prefix, module_name, method_name, source_name) {
  feature_names <- paste0(feature_prefix, make_safe_feature_name(rownames(score_matrix)))
  colnames_for_df <- feature_names

  score_df <- as.data.frame(t(score_matrix), check.names = FALSE)
  colnames(score_df) <- colnames_for_df
  score_df$sample_col_id <- rownames(score_df)
  score_df <- score_df[, c("sample_col_id", colnames_for_df), drop = FALSE]

  feature_info_df <- data.frame(
    feature_name = feature_names,
    feature_label = paste0(label_prefix, rownames(score_matrix)),
    module = module_name,
    method = method_name,
    source = source_name,
    stringsAsFactors = FALSE
  )

  list(score_df = score_df, feature_info_df = feature_info_df)
}

## 函数：calculate_mean_z_score ####
## 功能：对小基因集计算样本级mean z-score ####
## 输入：expression_matrix，行为基因、列为样本的log2表达矩阵；genes，signature基因向量；min_genes，最低可用基因数 ####
## 输出：样本级signature分数向量 ####
calculate_mean_z_score <- function(expression_matrix, genes, min_genes) {
  available_genes <- intersect(unique(as.character(genes)), rownames(expression_matrix))

  if (length(available_genes) < min_genes) {
    return(rep(NA_real_, ncol(expression_matrix)))
  }

  signature_matrix <- expression_matrix[available_genes, , drop = FALSE]
  gene_means <- rowMeans(signature_matrix, na.rm = TRUE)
  gene_sds <- matrixStats::rowSds(signature_matrix, na.rm = TRUE)
  valid_gene_mask <- !is.na(gene_sds) & gene_sds > 0

  if (sum(valid_gene_mask) < min_genes) {
    return(rep(NA_real_, ncol(expression_matrix)))
  }

  signature_matrix <- signature_matrix[valid_gene_mask, , drop = FALSE]
  gene_means <- gene_means[valid_gene_mask]
  gene_sds <- gene_sds[valid_gene_mask]

  z_matrix <- sweep(signature_matrix, 1, gene_means, "-")
  z_matrix <- sweep(z_matrix, 1, gene_sds, "/")

  colMeans(z_matrix, na.rm = TRUE)
}

## 函数：calculate_cyt_score ####
## 功能：按 Rooney et al. 2015 原公式 计算 cytolytic activity score = sqrt(GZMA_TPM × PRF1_TPM)（即 TPM 量纲下的几何平均），输入为 log2(TPM + log2_offset) 矩阵，函数内部先还原到 TPM 再算 ####
## 输入：expression_matrix，行为基因、列为样本的log2表达矩阵；log2_offset，log2 转换前加的偏移量 ####
## 输出：样本级 CYT score 向量（TPM 量纲，未取 log） ####
calculate_cyt_score <- function(expression_matrix, log2_offset = 1) {
  cyt_genes <- c("GZMA", "PRF1")
  available_genes <- intersect(cyt_genes, rownames(expression_matrix))

  if (length(available_genes) < 2) {
    return(rep(NA_real_, ncol(expression_matrix)))
  }

  ## 1. 反 log2 还原到 TPM 量纲：TPM = 2^x - log2_offset，截断为非负 ####
  log2_subset <- expression_matrix[available_genes, , drop = FALSE]
  tpm_subset <- pmax(2 ^ log2_subset - log2_offset, 0)

  ## 2. 几何平均：sqrt(GZMA_TPM × PRF1_TPM) ####
  sqrt(tpm_subset["GZMA", ] * tpm_subset["PRF1", ])
}

## 函数：calculate_single_gene_expression ####
## 功能：提取指定基因的样本级log2表达值 ####
## 输入：expression_matrix，行为基因、列为样本的log2表达矩阵；gene_symbol，目标基因符号 ####
## 输出：样本级基因表达向量 ####
calculate_single_gene_expression <- function(expression_matrix, gene_symbol) {
  if (!gene_symbol %in% rownames(expression_matrix)) {
    return(rep(NA_real_, ncol(expression_matrix)))
  }

  as.numeric(expression_matrix[gene_symbol, , drop = TRUE])
}

## 函数：calculate_tmb_scores ####
## 功能：基于本地TCGA MAF文件计算样本级TMB ####
## 输入：tumor_meta_df，肿瘤样本注释；maf_dir，MAF文件目录；exome_size_mb，外显子组大小Mb；use_parallel，是否启用并行；workers，并行进程数 ####
## 输出：列表，包含样本级TMB数据框和可用性汇总表 ####
calculate_tmb_scores <- function(tumor_meta_df, maf_dir, exome_size_mb, use_parallel = FALSE, workers = 1) {
  maf_files <- list.files(
    maf_dir,
    pattern = "_mutation\\.maf\\.gz$",
    full.names = TRUE
  )

  if (length(maf_files) == 0) {
    stop("未找到TCGA MAF文件目录或文件: ", maf_dir)
  }

  nonsynonymous_classes <- c(
    "Missense_Mutation",
    "Nonsense_Mutation",
    "Nonstop_Mutation",
    "Splice_Site",
    "Translation_Start_Site",
    "Frame_Shift_Del",
    "Frame_Shift_Ins",
    "In_Frame_Del",
    "In_Frame_Ins"
  )

  tmb_count_list <- run_parallel_lapply(maf_files, function(maf_file) {
    maf_df <- data.table::fread(
      maf_file,
      select = c(
        "Hugo_Symbol",
        "Chromosome",
        "Start_Position",
        "End_Position",
        "Reference_Allele",
        "Tumor_Seq_Allele2",
        "Variant_Classification",
        "Tumor_Sample_Barcode"
      ),
      data.table = FALSE,
      nThread = 1,
      showProgress = FALSE
    )

    maf_df %>%
      dplyr::filter(Variant_Classification %in% nonsynonymous_classes) %>%
      dplyr::mutate(sample_id = substr(Tumor_Sample_Barcode, 1, 16)) %>%
      dplyr::distinct(
        sample_id,
        Hugo_Symbol,
        Chromosome,
        Start_Position,
        End_Position,
        Reference_Allele,
        Tumor_Seq_Allele2
      ) %>%
      dplyr::count(sample_id, name = "nonsynonymous_mutation_count")
  }, use_parallel = use_parallel, workers = workers)

  tmb_count_df <- dplyr::bind_rows(tmb_count_list) %>%
    dplyr::group_by(sample_id) %>%
    dplyr::summarise(
      nonsynonymous_mutation_count = sum(nonsynonymous_mutation_count, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(TMB = nonsynonymous_mutation_count / exome_size_mb)

  tmb_score_df <- tumor_meta_df %>%
    dplyr::select(sample_col_id, sample_id, patient_id, project_id, primary_organ, pattern_group) %>%
    dplyr::left_join(tmb_count_df, by = "sample_id") %>%
    dplyr::mutate(has_tmb_data = !is.na(TMB)) %>%
    dplyr::select(sample_col_id, TMB, nonsynonymous_mutation_count, has_tmb_data)

  tmb_availability_df <- tmb_score_df %>%
    dplyr::summarise(
      n_total_samples = dplyr::n(),
      n_samples_with_tmb = sum(has_tmb_data, na.rm = TRUE),
      n_samples_without_tmb = sum(!has_tmb_data, na.rm = TRUE)
    )

  list(score_df = tmb_score_df, availability_df = tmb_availability_df)
}

## 函数：detect_tcr_columns ####
## 功能：自动识别TCR文件中的样本条码列和Shannon多样性列 ####
## 输入：tcr_df，TCR数据框 ####
## 输出：包含样本列名和Shannon列名的列表 ####
detect_tcr_columns <- function(tcr_df) {
  column_names <- colnames(tcr_df)
  sample_candidates <- column_names[grepl("sample|barcode|aliquot|patient", column_names, ignore.case = TRUE)]
  shannon_candidates <- column_names[grepl("shannon", column_names, ignore.case = TRUE)]

  if (length(sample_candidates) == 0 || length(shannon_candidates) == 0) {
    stop(
      "TCR文件无法自动识别样本列或Shannon列，当前列名: ",
      paste(column_names, collapse = ", ")
    )
  }

  list(
    sample_column = sample_candidates[[1]],
    shannon_column = shannon_candidates[[1]]
  )
}

## 函数：calculate_tcr_scores ####
## 功能：读取Thorsson/TCGA PanImmune TCR summary并提取样本级TCR Shannon多样性 ####
## 输入：tumor_meta_df，肿瘤样本注释；tcr_file，TCR summary文件路径 ####
## 输出：列表，包含样本级TCR数据框和可用性汇总表 ####
calculate_tcr_scores <- function(tumor_meta_df, tcr_file) {
  if (!file.exists(tcr_file)) {
    stop(
      "缺少TCR diversity输入文件: ",
      tcr_file,
      "\n请将Thorsson/TCGA PanImmune mitcr_sampleStatistics_20160714.tsv放到该路径后重新运行。"
    )
  }

  tcr_df <- data.table::fread(tcr_file, data.table = FALSE, showProgress = FALSE)
  tcr_columns <- detect_tcr_columns(tcr_df)

  tcr_score_source_df <- tcr_df %>%
    dplyr::transmute(
      sample_id = substr(as.character(.data[[tcr_columns$sample_column]]), 1, 16),
      TCR_Shannon = as.numeric(.data[[tcr_columns$shannon_column]])
    ) %>%
    dplyr::filter(!is.na(sample_id), sample_id != "") %>%
    dplyr::group_by(sample_id) %>%
    dplyr::summarise(TCR_Shannon = median(TCR_Shannon, na.rm = TRUE), .groups = "drop")

  tcr_score_df <- tumor_meta_df %>%
    dplyr::select(sample_col_id, sample_id, patient_id, project_id, primary_organ, pattern_group) %>%
    dplyr::left_join(tcr_score_source_df, by = "sample_id") %>%
    dplyr::mutate(has_tcr_data = !is.na(TCR_Shannon)) %>%
    dplyr::select(sample_col_id, TCR_Shannon, has_tcr_data)

  tcr_availability_df <- tcr_score_df %>%
    dplyr::summarise(
      n_total_samples = dplyr::n(),
      n_samples_with_tcr = sum(has_tcr_data, na.rm = TRUE),
      n_samples_without_tcr = sum(!has_tcr_data, na.rm = TRUE),
      detected_sample_column = tcr_columns$sample_column,
      detected_shannon_column = tcr_columns$shannon_column
    )

  list(score_df = tcr_score_df, availability_df = tcr_availability_df)
}

## 函数：calculate_estimate_scores ####
## 功能：调用tidyestimate原始ESTIMATE算法 (Yoshihara et al., 2013) 计算ImmuneScore/StromalScore/ESTIMATEScore；不输出TumorPurity，因原作者明确说cos(...)公式仅适用于Affymetrix数据，对RNA-seq不可作绝对purity解释 ####
## 输入：expression_matrix，行为基因(hgnc_symbol)、列为样本的log2表达矩阵；alias_map和并行参数已无须；保留参数签名以保持向后兼容 ####
## 输出：列表，包含样本级ESTIMATE分数、feature信息和signature可用性 ####
calculate_estimate_scores <- function(expression_matrix, alias_map = NULL, use_parallel = FALSE, workers = 1) {
  ## 1. 转换为 tidyestimate 期望的 data.frame：第一列 hgnc_symbol，其余为各样本表达 ####
  expr_df <- as.data.frame(expression_matrix, check.names = FALSE)
  expr_df <- tibble::rownames_to_column(expr_df, var = "hgnc_symbol")

  ## 2. tidyestimate::filter_common_genes：过滤为ESTIMATE训练时使用的常见基因，并启用别名解析 ####
  ## tidy=TRUE 表示df第1列已是基因ID，不要再次从rownames读取 ####
  filtered_df <- tidyestimate::filter_common_genes(
    expr_df,
    id = "hgnc_symbol",
    tidy = TRUE,
    tell_missing = FALSE,
    find_alias = TRUE
  )

  ## 3. tidyestimate::estimate_score：is_affymetrix=FALSE 适配 RNA-seq 数据 ####
  estimate_result_df <- tidyestimate::estimate_score(
    filtered_df,
    is_affymetrix = FALSE
  )

  ## 4. 整理为统一的样本级表 ####
  estimate_score_df <- data.frame(
    sample_col_id = as.character(estimate_result_df$sample),
    ImmuneScore = as.numeric(estimate_result_df$immune),
    StromalScore = as.numeric(estimate_result_df$stromal),
    ESTIMATEScore = as.numeric(estimate_result_df$estimate),
    stringsAsFactors = FALSE
  )

  ## 5. 可用性表：基于 ESTIMATE 训练时的 common gene panel，参考行展示有多少基因被覆盖 ####
  common_gene_panel <- as.character(filtered_df$hgnc_symbol)
  estimate_signature_list <- list(
    "ImmuneScore" = common_gene_panel,
    "StromalScore" = common_gene_panel
  )
  estimate_availability_df <- format_signature_availability(
    signature_list = estimate_signature_list,
    available_genes = rownames(expression_matrix),
    module_name = "Overall TME composition",
    method_name = "tidyestimate::estimate_score (RNA-seq mode)",
    min_genes = min_genes_default
  )

  estimate_feature_info_df <- data.frame(
    feature_name = c("ImmuneScore", "StromalScore", "ESTIMATEScore"),
    feature_label = c("ImmuneScore", "StromalScore", "ESTIMATEScore"),
    module = "Overall TME composition",
    method = c(
      "tidyestimate::estimate_score (Yoshihara et al., 2013), is_affymetrix=FALSE",
      "tidyestimate::estimate_score (Yoshihara et al., 2013), is_affymetrix=FALSE",
      "ImmuneScore + StromalScore"
    ),
    source = "Yoshihara et al., 2013; tidyestimate::estimate_score",
    stringsAsFactors = FALSE
  )

  list(
    score_df = estimate_score_df,
    feature_info_df = estimate_feature_info_df,
    availability_df = estimate_availability_df
  )
}

# 6. 汇总与统计函数 ####

## 函数：summarise_scores_by_group ####
## 功能：按指定分组变量汇总样本级免疫特征为中位数 ####
## 输入：sample_score_df，样本级评分表；group_columns，分组列名；feature_columns，特征列名 ####
## 输出：分组汇总后的评分表 ####
summarise_scores_by_group <- function(sample_score_df, group_columns, feature_columns) {
  sample_score_df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_columns))) %>%
    dplyr::summarise(
      n_tumor_samples = dplyr::n(),
      n_patients = dplyr::n_distinct(patient_id),
      dplyr::across(
        dplyr::all_of(feature_columns),
        ~ safe_quantile(.x, 0.5)
      ),
      .groups = "drop"
    )
}

## 函数：build_pattern_group_summary ####
## 功能：按pattern group和feature整理组内中位数、四分位数及样本数 ####
## 输入：primary_organ_score_df，原发部位级评分表；feature_columns，特征列名；feature_info_df，特征信息表 ####
## 输出：模式群特征汇总长表 ####
build_pattern_group_summary <- function(primary_organ_score_df, feature_columns, feature_info_df) {
  primary_organ_score_df %>%
    dplyr::select(pattern_group, primary_organ, dplyr::all_of(feature_columns)) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(feature_columns),
      names_to = "feature_name",
      values_to = "score"
    ) %>%
    dplyr::group_by(pattern_group, feature_name) %>%
    dplyr::summarise(
      n_primary_organs = sum(!is.na(score)),
      median_score = safe_quantile(score, 0.5),
      q1_score = safe_quantile(score, 0.25),
      q3_score = safe_quantile(score, 0.75),
      .groups = "drop"
    ) %>%
    dplyr::left_join(feature_info_df, by = "feature_name") %>%
    dplyr::arrange(axis_order, feature_order, factor(pattern_group, levels = pattern_group_order))
}

## 函数：run_kruskal_tests ####
## 功能：对每个TME特征进行pattern group间Kruskal-Wallis检验 ####
## 输入：primary_organ_score_df，原发部位级评分表；feature_columns，特征列名；feature_info_df，特征信息表 ####
## 输出：Kruskal-Wallis检验结果表 ####
run_kruskal_tests <- function(primary_organ_score_df, feature_columns, feature_info_df) {
  test_result_list <- lapply(feature_columns, function(current_feature) {
    test_df <- primary_organ_score_df %>%
      dplyr::transmute(pattern_group = pattern_group, score = .data[[current_feature]]) %>%
      dplyr::filter(!is.na(score), !is.na(pattern_group))

    group_median_df <- test_df %>%
      dplyr::group_by(pattern_group) %>%
      dplyr::summarise(
        n_primary_organs = dplyr::n(),
        median_score = median(score, na.rm = TRUE),
        .groups = "drop"
      )

    n_groups <- dplyr::n_distinct(test_df$pattern_group)
    n_values <- dplyr::n_distinct(test_df$score)

    if (n_groups < 2 || n_values < 2) {
      p_value <- NA_real_
      statistic <- NA_real_
    } else {
      test_result <- kruskal.test(score ~ pattern_group, data = test_df)
      p_value <- test_result$p.value
      statistic <- unname(test_result$statistic)
    }

    if (nrow(group_median_df) == 0) {
      high_group <- NA_character_
      low_group <- NA_character_
      median_range <- NA_real_
    } else {
      high_group <- group_median_df$pattern_group[which.max(group_median_df$median_score)]
      low_group <- group_median_df$pattern_group[which.min(group_median_df$median_score)]
      median_range <- max(group_median_df$median_score, na.rm = TRUE) - min(group_median_df$median_score, na.rm = TRUE)
    }

    data.frame(
      feature_name = current_feature,
      n_groups = n_groups,
      n_primary_organs = nrow(test_df),
      statistic = statistic,
      p_value = p_value,
      highest_pattern_group = high_group,
      lowest_pattern_group = low_group,
      median_range = median_range,
      stringsAsFactors = FALSE
    )
  })

  dplyr::bind_rows(test_result_list) %>%
    dplyr::mutate(fdr = p.adjust(p_value, method = "BH")) %>%
    dplyr::left_join(feature_info_df, by = "feature_name") %>%
    dplyr::arrange(fdr, p_value)
}

## 函数：run_pairwise_wilcox_tests ####
## 功能：对每个TME特征进行所有pattern group两两Wilcoxon秩和检验 ####
## 输入：primary_organ_score_df，原发部位级评分表；feature_columns，特征列名；feature_info_df，特征信息表 ####
## 输出：两两比较结果表 ####
run_pairwise_wilcox_tests <- function(primary_organ_score_df, feature_columns, feature_info_df) {
  group_levels <- pattern_group_order[pattern_group_order %in% unique(primary_organ_score_df$pattern_group)]
  group_pairs <- utils::combn(group_levels, 2, simplify = FALSE)

  pairwise_result_list <- lapply(feature_columns, function(current_feature) {
    feature_result_list <- lapply(group_pairs, function(current_pair) {
      group1 <- current_pair[[1]]
      group2 <- current_pair[[2]]

      group1_values <- primary_organ_score_df[
        primary_organ_score_df$pattern_group == group1,
        current_feature,
        drop = TRUE
      ]
      group2_values <- primary_organ_score_df[
        primary_organ_score_df$pattern_group == group2,
        current_feature,
        drop = TRUE
      ]

      group1_values <- group1_values[!is.na(group1_values)]
      group2_values <- group2_values[!is.na(group2_values)]

      if (length(group1_values) == 0 || length(group2_values) == 0 || length(unique(c(group1_values, group2_values))) < 2) {
        p_value <- NA_real_
        statistic <- NA_real_
      } else {
        test_result <- wilcox.test(group1_values, group2_values, exact = FALSE)
        p_value <- test_result$p.value
        statistic <- unname(test_result$statistic)
      }

      data.frame(
        feature_name = current_feature,
        group1 = group1,
        group2 = group2,
        n_group1 = length(group1_values),
        n_group2 = length(group2_values),
        median_group1 = median(group1_values, na.rm = TRUE),
        median_group2 = median(group2_values, na.rm = TRUE),
        median_diff_group2_minus_group1 = median(group2_values, na.rm = TRUE) - median(group1_values, na.rm = TRUE),
        statistic = statistic,
        p_value = p_value,
        stringsAsFactors = FALSE
      )
    })

    dplyr::bind_rows(feature_result_list)
  })

  dplyr::bind_rows(pairwise_result_list) %>%
    dplyr::mutate(fdr = p.adjust(p_value, method = "BH")) %>%
    dplyr::left_join(feature_info_df, by = "feature_name") %>%
    dplyr::arrange(feature_name, fdr, p_value)
}

## 函数：safe_quantile ####
## 功能：在允许全NA输入的情况下计算分位数 ####
## 输入：x，数值向量；prob，分位数概率 ####
## 输出：单个数值分位数或NA ####
safe_quantile <- function(x, prob) {
  valid_x <- x[!is.na(x)]

  if (length(valid_x) == 0) {
    return(NA_real_)
  }

  stats::quantile(valid_x, probs = prob, na.rm = TRUE, names = FALSE, type = 7)
}

## 函数：get_extreme_group ####
## 功能：从一个分组得分向量中提取最高或最低分组 ####
## 输入：group_scores，带分组名称的数值向量；direction，最高或最低 ####
## 输出：最高或最低分组名称 ####
get_extreme_group <- function(group_scores, direction = c("highest", "lowest")) {
  direction <- match.arg(direction)
  valid_index <- which(!is.na(group_scores))

  if (length(valid_index) == 0) {
    return(NA_character_)
  }

  valid_scores <- group_scores[valid_index]

  if (direction == "highest") {
    names(valid_scores)[which.max(valid_scores)]
  } else {
    names(valid_scores)[which.min(valid_scores)]
  }
}

## 函数：run_stratified_tme_bootstrap ####
## 功能：执行两层分层bootstrap，先在每个primary organ内等量抽样本，再在每个pattern group内等量抽primary organs ####
## 输入：sample_score_df，样本级评分表；feature_columns，特征列名；n_iterations，抽样次数；samples_per_primary_organ，每个癌种每轮抽取样本数；min_samples_per_primary_organ，纳入癌种最小样本数；primary_organs_per_pattern_group，每组每轮抽取癌种数；seed，随机种子；use_parallel，是否启用并行；workers，并行进程数 ####
## 输出：包含bootstrap三维数组、纳入癌种信息和参数表的列表 ####
run_stratified_tme_bootstrap <- function(sample_score_df,
                                         feature_columns,
                                         n_iterations,
                                         samples_per_primary_organ,
                                         min_samples_per_primary_organ,
                                         primary_organs_per_pattern_group,
                                         seed,
                                         use_parallel = FALSE,
                                         workers = 1) {
  set.seed(seed)

  bootstrap_input_df <- sample_score_df %>%
    dplyr::select(
      pattern_group,
      primary_organ,
      sample_barcode,
      patient_id,
      dplyr::all_of(feature_columns)
    ) %>%
    dplyr::filter(
      !is.na(pattern_group),
      !is.na(primary_organ),
      pattern_group %in% pattern_group_order
    )

  organ_info_df <- bootstrap_input_df %>%
    dplyr::group_by(pattern_group, primary_organ) %>%
    dplyr::summarise(
      n_tumor_samples = dplyr::n(),
      n_patients = dplyr::n_distinct(patient_id),
      .groups = "drop"
    ) %>%
    dplyr::filter(n_tumor_samples >= min_samples_per_primary_organ) %>%
    dplyr::arrange(factor(pattern_group, levels = pattern_group_order), primary_organ)

  available_group_levels <- pattern_group_order[pattern_group_order %in% unique(organ_info_df$pattern_group)]
  organs_by_group <- split(organ_info_df$primary_organ, organ_info_df$pattern_group)
  organs_by_group <- organs_by_group[available_group_levels]

  if (length(available_group_levels) < 2) {
    stop("分层bootstrap至少需要2个pattern group。")
  }

  group_organ_count <- vapply(organs_by_group, length, numeric(1))

  if (is.na(primary_organs_per_pattern_group)) {
    primary_organs_per_pattern_group <- min(group_organ_count)
  }

  if (primary_organs_per_pattern_group < 1) {
    stop("primary_organs_per_pattern_group必须大于等于1。")
  }

  included_primary_organs <- organ_info_df$primary_organ
  organ_sample_matrix_list <- lapply(included_primary_organs, function(current_primary_organ) {
    current_df <- bootstrap_input_df[bootstrap_input_df$primary_organ == current_primary_organ, feature_columns, drop = FALSE]
    current_matrix <- as.matrix(current_df)
    storage.mode(current_matrix) <- "numeric"
    current_matrix
  })
  names(organ_sample_matrix_list) <- included_primary_organs

  iteration_seed_vector <- sample.int(.Machine$integer.max, n_iterations)
  run_single_bootstrap_iteration <- function(iteration_index) {
    set.seed(iteration_seed_vector[[iteration_index]])

    organ_score_matrix <- matrix(
      NA_real_,
      nrow = length(included_primary_organs),
      ncol = length(feature_columns),
      dimnames = list(included_primary_organs, feature_columns)
    )

    group_score_matrix <- matrix(
      NA_real_,
      nrow = length(available_group_levels),
      ncol = length(feature_columns),
      dimnames = list(available_group_levels, feature_columns)
    )

    for (current_primary_organ in included_primary_organs) {
      current_matrix <- organ_sample_matrix_list[[current_primary_organ]]
      selected_sample_index <- sample(
        seq_len(nrow(current_matrix)),
        size = samples_per_primary_organ,
        replace = TRUE
      )
      current_organ_score <- matrixStats::colMedians(
        current_matrix[selected_sample_index, , drop = FALSE],
        na.rm = TRUE
      )
      current_organ_score[is.nan(current_organ_score)] <- NA_real_
      organ_score_matrix[current_primary_organ, ] <- current_organ_score
    }

    for (current_group in available_group_levels) {
      current_organs <- organs_by_group[[current_group]]
      selected_organs <- sample(
        current_organs,
        size = primary_organs_per_pattern_group,
        replace = TRUE
      )
      current_group_score <- matrixStats::colMedians(
        organ_score_matrix[selected_organs, , drop = FALSE],
        na.rm = TRUE
      )
      current_group_score[is.nan(current_group_score)] <- NA_real_
      group_score_matrix[current_group, ] <- current_group_score
    }

    group_score_matrix
  }

  bootstrap_iteration_list <- run_parallel_lapply(
    seq_len(n_iterations),
    run_single_bootstrap_iteration,
    use_parallel = use_parallel,
    workers = workers
  )

  bootstrap_group_median_array <- array(
    NA_real_,
    dim = c(n_iterations, length(available_group_levels), length(feature_columns)),
    dimnames = list(
      iteration = paste0("iteration_", seq_len(n_iterations)),
      pattern_group = available_group_levels,
      feature_name = feature_columns
    )
  )

  for (iteration_index in seq_len(n_iterations)) {
    bootstrap_group_median_array[iteration_index, , ] <- bootstrap_iteration_list[[iteration_index]]
  }

  parameter_df <- tibble::tibble(
    parameter = c(
      "algorithm",
      "n_iterations",
      "seed",
      "samples_per_primary_organ",
      "min_samples_per_primary_organ",
      "primary_organs_per_pattern_group",
      "included_primary_organs",
      "included_pattern_groups",
      "enable_parallel",
      "parallel_workers",
      "interpretation_note"
    ),
    value = c(
      "Two-stage balanced hierarchical bootstrap: samples within primary organs, then primary organs within pattern groups",
      as.character(n_iterations),
      as.character(seed),
      as.character(samples_per_primary_organ),
      as.character(min_samples_per_primary_organ),
      as.character(primary_organs_per_pattern_group),
      as.character(nrow(organ_info_df)),
      as.character(length(available_group_levels)),
      as.character(use_parallel),
      as.character(workers),
      "Bootstrap probabilities describe directional stability and should not be interpreted as independent-sample p values."
    )
  )

  list(
    bootstrap_group_median_array = bootstrap_group_median_array,
    organ_info_df = organ_info_df,
    parameter_df = parameter_df
  )
}

## 函数：summarise_bootstrap_group_scores ####
## 功能：汇总每个特征在每个pattern group中的bootstrap中位数、置信区间和最高/最低概率 ####
## 输入：bootstrap_group_median_array，bootstrap三维数组；feature_info_df，特征信息表 ####
## 输出：bootstrap模式群评分汇总表 ####
summarise_bootstrap_group_scores <- function(bootstrap_group_median_array, feature_info_df) {
  group_levels <- dimnames(bootstrap_group_median_array)$pattern_group
  feature_names <- dimnames(bootstrap_group_median_array)$feature_name

  summary_list <- lapply(feature_names, function(current_feature) {
    feature_matrix <- bootstrap_group_median_array[, , current_feature, drop = TRUE]

    feature_result_list <- lapply(group_levels, function(current_group) {
      current_values <- feature_matrix[, current_group]
      highest_group_vector <- apply(feature_matrix, 1, get_extreme_group, direction = "highest")
      lowest_group_vector <- apply(feature_matrix, 1, get_extreme_group, direction = "lowest")

      data.frame(
        feature_name = current_feature,
        pattern_group = current_group,
        n_iterations = sum(!is.na(current_values)),
        bootstrap_median_score = safe_quantile(current_values, 0.5),
        bootstrap_ci_lower = safe_quantile(current_values, 0.025),
        bootstrap_ci_upper = safe_quantile(current_values, 0.975),
        probability_highest = mean(highest_group_vector == current_group, na.rm = TRUE),
        probability_lowest = mean(lowest_group_vector == current_group, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    })

    dplyr::bind_rows(feature_result_list)
  })

  dplyr::bind_rows(summary_list) %>%
    dplyr::left_join(feature_info_df, by = "feature_name") %>%
    dplyr::arrange(feature_name, factor(pattern_group, levels = pattern_group_order))
}

## 函数：summarise_bootstrap_feature_stability ####
## 功能：汇总每个特征最高分组、最低分组和组间范围的bootstrap稳定性 ####
## 输入：bootstrap_group_median_array，bootstrap三维数组；feature_info_df，特征信息表 ####
## 输出：bootstrap特征稳定性汇总表 ####
summarise_bootstrap_feature_stability <- function(bootstrap_group_median_array, feature_info_df) {
  feature_names <- dimnames(bootstrap_group_median_array)$feature_name

  stability_list <- lapply(feature_names, function(current_feature) {
    feature_matrix <- bootstrap_group_median_array[, , current_feature, drop = TRUE]

    highest_group_vector <- apply(feature_matrix, 1, get_extreme_group, direction = "highest")
    lowest_group_vector <- apply(feature_matrix, 1, get_extreme_group, direction = "lowest")
    range_vector <- apply(feature_matrix, 1, function(current_scores) {
      valid_scores <- current_scores[!is.na(current_scores)]

      if (length(valid_scores) == 0) {
        return(NA_real_)
      }

      max(valid_scores) - min(valid_scores)
    })

    highest_group_table <- table(factor(highest_group_vector, levels = pattern_group_order))
    lowest_group_table <- table(factor(lowest_group_vector, levels = pattern_group_order))

    highest_probability <- if (sum(highest_group_table) == 0) {
      rep(NA_real_, length(pattern_group_order))
    } else {
      as.numeric(highest_group_table) / sum(highest_group_table)
    }

    lowest_probability <- if (sum(lowest_group_table) == 0) {
      rep(NA_real_, length(pattern_group_order))
    } else {
      as.numeric(lowest_group_table) / sum(lowest_group_table)
    }

    most_often_highest_group <- if (all(is.na(highest_probability))) {
      NA_character_
    } else {
      pattern_group_order[which.max(highest_probability)]
    }

    highest_group_probability <- if (all(is.na(highest_probability))) {
      NA_real_
    } else {
      max(highest_probability, na.rm = TRUE)
    }

    most_often_lowest_group <- if (all(is.na(lowest_probability))) {
      NA_character_
    } else {
      pattern_group_order[which.max(lowest_probability)]
    }

    lowest_group_probability <- if (all(is.na(lowest_probability))) {
      NA_real_
    } else {
      max(lowest_probability, na.rm = TRUE)
    }

    data.frame(
      feature_name = current_feature,
      most_often_highest_group = most_often_highest_group,
      highest_group_probability = highest_group_probability,
      most_often_lowest_group = most_often_lowest_group,
      lowest_group_probability = lowest_group_probability,
      bootstrap_median_range = safe_quantile(range_vector, 0.5),
      bootstrap_range_ci_lower = safe_quantile(range_vector, 0.025),
      bootstrap_range_ci_upper = safe_quantile(range_vector, 0.975),
      stringsAsFactors = FALSE
    )
  })

  dplyr::bind_rows(stability_list) %>%
    dplyr::left_join(feature_info_df, by = "feature_name") %>%
    dplyr::arrange(dplyr::desc(highest_group_probability), dplyr::desc(bootstrap_median_range))
}

## 函数：summarise_bootstrap_pairwise_direction ####
## 功能：汇总每个特征在两两pattern group之间差值方向的bootstrap稳定性 ####
## 输入：bootstrap_group_median_array，bootstrap三维数组；feature_info_df，特征信息表 ####
## 输出：bootstrap两两方向稳定性汇总表 ####
summarise_bootstrap_pairwise_direction <- function(bootstrap_group_median_array, feature_info_df) {
  group_levels <- dimnames(bootstrap_group_median_array)$pattern_group
  feature_names <- dimnames(bootstrap_group_median_array)$feature_name
  group_pairs <- utils::combn(group_levels, 2, simplify = FALSE)

  pairwise_list <- lapply(feature_names, function(current_feature) {
    feature_matrix <- bootstrap_group_median_array[, , current_feature, drop = TRUE]

    feature_pairwise_list <- lapply(group_pairs, function(current_pair) {
      group1 <- current_pair[[1]]
      group2 <- current_pair[[2]]
      diff_values <- feature_matrix[, group2] - feature_matrix[, group1]
      valid_diff_values <- diff_values[!is.na(diff_values)]

      if (length(valid_diff_values) == 0) {
        probability_group2_higher <- NA_real_
        probability_group1_higher <- NA_real_
        direction_consistency <- NA_real_
        median_diff <- NA_real_
      } else {
        probability_group2_higher <- mean(valid_diff_values > 0)
        probability_group1_higher <- mean(valid_diff_values < 0)
        direction_consistency <- max(probability_group2_higher, probability_group1_higher, na.rm = TRUE)
        median_diff <- median(valid_diff_values, na.rm = TRUE)
      }

      data.frame(
        feature_name = current_feature,
        group1 = group1,
        group2 = group2,
        n_iterations = length(valid_diff_values),
        bootstrap_median_diff_group2_minus_group1 = median_diff,
        bootstrap_diff_ci_lower = safe_quantile(valid_diff_values, 0.025),
        bootstrap_diff_ci_upper = safe_quantile(valid_diff_values, 0.975),
        probability_group2_higher = probability_group2_higher,
        probability_group1_higher = probability_group1_higher,
        direction_consistency = direction_consistency,
        stringsAsFactors = FALSE
      )
    })

    dplyr::bind_rows(feature_pairwise_list)
  })

  dplyr::bind_rows(pairwise_list) %>%
    dplyr::left_join(feature_info_df, by = "feature_name") %>%
    dplyr::arrange(dplyr::desc(direction_consistency), feature_name, group1, group2)
}

## 函数：classify_direction_stability ####
## 功能：根据bootstrap最高组概率划分方向稳定性等级 ####
## 输入：probability，最高组或最低组bootstrap概率 ####
## 输出：方向稳定性标签 ####
classify_direction_stability <- function(probability) {
  if (is.na(probability)) {
    return("Not evaluated")
  }

  if (probability >= 0.80) {
    return("Stable")
  }

  if (probability >= 0.60) {
    return("Trend")
  }

  "Unstable"
}

## 函数：classify_effect_size ####
## 功能：根据标准化组间中位数范围划分effect size等级 ####
## 输入：standardized_range，标准化后的组间中位数范围 ####
## 输出：effect size等级标签 ####
classify_effect_size <- function(standardized_range) {
  if (is.na(standardized_range)) {
    return("Not evaluated")
  }

  if (standardized_range >= 0.80) {
    return("Large")
  }

  if (standardized_range >= 0.50) {
    return("Moderate")
  }

  "Small"
}

## 函数：extract_model_table_value ####
## 功能：从模型结果表中按行名和候选列名安全提取数值 ####
## 输入：table_df，模型结果数据框；row_name，目标行名；candidate_columns，候选列名向量 ####
## 输出：单个数值或NA ####
extract_model_table_value <- function(table_df, row_name, candidate_columns) {
  if (is.null(table_df) || nrow(table_df) == 0) {
    return(NA_real_)
  }

  row_index <- which(rownames(table_df) == row_name)
  column_name <- intersect(candidate_columns, colnames(table_df))

  if (length(row_index) == 0 || length(column_name) == 0) {
    return(NA_real_)
  }

  as.numeric(table_df[row_index[[1]], column_name[[1]]])
}

## 函数：prepare_mixed_model_feature_data ####
## 功能：整理单个feature的样本级混合效应模型输入并执行标准化 ####
## 输入：sample_score_df，样本级免疫特征表；feature_name，目标特征名 ####
## 输出：包含标准化feature和分组变量的数据框 ####
prepare_mixed_model_feature_data <- function(sample_score_df, feature_name) {
  feature_values <- as.numeric(sample_score_df[[feature_name]])
  transformed_values <- if (feature_name == "TMB") {
    log1p(feature_values)
  } else {
    feature_values
  }

  feature_mean <- mean(transformed_values, na.rm = TRUE)
  feature_sd <- stats::sd(transformed_values, na.rm = TRUE)
  feature_z <- if (is.na(feature_sd) || feature_sd == 0) {
    rep(NA_real_, length(transformed_values))
  } else {
    (transformed_values - feature_mean) / feature_sd
  }

  data.frame(
    pattern_group = as.character(sample_score_df$pattern_group),
    primary_organ = as.character(sample_score_df$primary_organ),
    feature_value_raw = feature_values,
    feature_value_transformed = transformed_values,
    feature_value_z = feature_z,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::filter(
      !is.na(pattern_group),
      !is.na(primary_organ),
      !is.na(feature_value_z),
      pattern_group %in% pattern_group_order,
      primary_organ != ""
    ) %>%
    dplyr::mutate(
      pattern_group = factor(pattern_group, levels = pattern_group_order),
      primary_organ = factor(primary_organ)
    )
}

## 函数：build_mixed_model_failure_row ####
## 功能：为无法拟合的feature生成混合模型整体检验占位行 ####
## 输入：feature_name，目标特征名；model_data_df，模型输入数据；model_status，模型状态；model_message，状态说明 ####
## 输出：混合模型整体检验数据框单行 ####
build_mixed_model_failure_row <- function(feature_name, model_data_df, model_status, model_message) {
  data.frame(
    feature_name = feature_name,
    model_status = model_status,
    model_message = model_message,
    transformed_feature = ifelse(feature_name == "TMB", "log1p_then_z_score", "z_score"),
    n_samples = if (is.null(model_data_df)) 0 else nrow(model_data_df),
    n_primary_organs = if (is.null(model_data_df)) 0 else dplyr::n_distinct(model_data_df$primary_organ),
    n_pattern_groups = if (is.null(model_data_df)) 0 else dplyr::n_distinct(model_data_df$pattern_group),
    mixed_highest_group = NA_character_,
    mixed_lowest_group = NA_character_,
    mixed_highest_emmean = NA_real_,
    mixed_lowest_emmean = NA_real_,
    mixed_standardized_range = NA_real_,
    mixed_pattern_group_num_df = NA_real_,
    mixed_pattern_group_den_df = NA_real_,
    mixed_pattern_group_f_value = NA_real_,
    mixed_pattern_group_p_value = NA_real_,
    is_singular = NA,
    stringsAsFactors = FALSE
  )
}

## 函数：parse_emmeans_pairwise_contrast ####
## 功能：解析emmeans两两比较中的contrast标签，去除含特殊字符level(如"Systemic-like")被emmeans自动添加的外层括号 ####
## 输入：contrast_values，形如group1 - group2的contrast字符串向量 ####
## 输出：包含group1和group2的数据框 ####
parse_emmeans_pairwise_contrast <- function(contrast_values) {
  contrast_parts <- strsplit(as.character(contrast_values), " - ", fixed = TRUE)

  strip_outer_parens <- function(x) {
    if (is.na(x)) return(NA_character_)
    sub("^\\((.*)\\)$", "\\1", x)
  }

  data.frame(
    group1 = vapply(contrast_parts, function(x) strip_outer_parens(ifelse(length(x) >= 1, x[[1]], NA_character_)), character(1)),
    group2 = vapply(contrast_parts, function(x) strip_outer_parens(ifelse(length(x) >= 2, x[[2]], NA_character_)), character(1)),
    stringsAsFactors = FALSE
  )
}

## 函数：build_one_vs_rest_contrast_list ####
## 功能：根据emmeans可用pattern_group构造one-vs-rest对比列表(每群vs其它群平均)，用于emmeans::contrast ####
## 输入：available_group_levels，emmeans输出的pattern_group顺序(字符向量) ####
## 输出：named list，每个元素为一个one-vs-rest对比向量；groups<2时返回空list ####
build_one_vs_rest_contrast_list <- function(available_group_levels) {
  n_groups <- length(available_group_levels)
  if (n_groups < 2) return(list())

  contrast_list <- vector("list", n_groups)
  for (i in seq_len(n_groups)) {
    contrast_vec <- rep(-1 / (n_groups - 1), n_groups)
    contrast_vec[i] <- 1
    contrast_list[[i]] <- contrast_vec
  }
  names(contrast_list) <- paste0(available_group_levels, "_vs_rest")
  contrast_list
}

## 函数：fit_single_mixed_model_feature ####
## 功能：对单个feature拟合pattern_group固定效应和primary_organ随机截距模型 ####
## 输入：sample_score_df，样本级免疫特征表；feature_name，目标特征名 ####
## 输出：列表，包含混合模型整体检验和两两比较结果 ####
fit_single_mixed_model_feature <- function(sample_score_df, feature_name) {
  model_data_df <- prepare_mixed_model_feature_data(
    sample_score_df = sample_score_df,
    feature_name = feature_name
  )

  if (
    nrow(model_data_df) < 3 ||
      dplyr::n_distinct(model_data_df$pattern_group) < 2 ||
      dplyr::n_distinct(model_data_df$primary_organ) < 2 ||
      dplyr::n_distinct(model_data_df$feature_value_z) < 2
  ) {
    return(list(
      overall_df = build_mixed_model_failure_row(
        feature_name = feature_name,
        model_data_df = model_data_df,
        model_status = "insufficient_data",
        model_message = "Mixed model requires at least two pattern groups, two primary organs and non-constant feature values."
      ),
      pairwise_df = data.frame(),
      emmean_df = data.frame(),
      one_vs_rest_df = data.frame()
    ))
  }

  available_group_levels <- pattern_group_order[pattern_group_order %in% unique(as.character(model_data_df$pattern_group))]
  model_data_df$pattern_group <- factor(as.character(model_data_df$pattern_group), levels = available_group_levels)
  model_data_df$primary_organ <- factor(as.character(model_data_df$primary_organ))

  model_fit <- tryCatch(
    lmerTest::lmer(
      feature_value_z ~ pattern_group + (1 | primary_organ),
      data = model_data_df,
      REML = FALSE,
      control = lme4::lmerControl(optimizer = "bobyqa")
    ),
    error = function(e) e
  )

  if (inherits(model_fit, "error")) {
    return(list(
      overall_df = build_mixed_model_failure_row(
        feature_name = feature_name,
        model_data_df = model_data_df,
        model_status = "model_failed",
        model_message = conditionMessage(model_fit)
      ),
      pairwise_df = data.frame(),
      emmean_df = data.frame(),
      one_vs_rest_df = data.frame()
    ))
  }

  ## 大样本下 lmerTest 同样会在 anova 中打印 d.f. 停用 Note，统一抑制 ####
  model_anova <- tryCatch(
    suppressMessages(as.data.frame(stats::anova(model_fit, type = 3))),
    error = function(e) data.frame()
  )

  ## 大样本下显式使用渐近 z 检验，避免 emmeans 反复打印 d.f. 停用 Note ####
  emm_result <- tryCatch(
    emmeans::emmeans(model_fit, specs = "pattern_group", lmer.df = "asymptotic"),
    error = function(e) e
  )

  if (inherits(emm_result, "error")) {
    return(list(
      overall_df = build_mixed_model_failure_row(
        feature_name = feature_name,
        model_data_df = model_data_df,
        model_status = "emmeans_failed",
        model_message = conditionMessage(emm_result)
      ),
      pairwise_df = data.frame(),
      emmean_df = data.frame(),
      one_vs_rest_df = data.frame()
    ))
  }

  emm_df <- as.data.frame(emm_result)
  emm_df <- emm_df[!is.na(emm_df$emmean), , drop = FALSE]

  if (nrow(emm_df) == 0) {
    return(list(
      overall_df = build_mixed_model_failure_row(
        feature_name = feature_name,
        model_data_df = model_data_df,
        model_status = "emmeans_failed",
        model_message = "Estimated marginal means were all missing."
      ),
      pairwise_df = data.frame(),
      emmean_df = data.frame(),
      one_vs_rest_df = data.frame()
    ))
  }

  highest_index <- which.max(emm_df$emmean)
  lowest_index <- which.min(emm_df$emmean)
  mixed_highest_emmean <- emm_df$emmean[[highest_index]]
  mixed_lowest_emmean <- emm_df$emmean[[lowest_index]]

  overall_df <- data.frame(
    feature_name = feature_name,
    model_status = "success",
    model_message = NA_character_,
    transformed_feature = ifelse(feature_name == "TMB", "log1p_then_z_score", "z_score"),
    n_samples = nrow(model_data_df),
    n_primary_organs = dplyr::n_distinct(model_data_df$primary_organ),
    n_pattern_groups = dplyr::n_distinct(model_data_df$pattern_group),
    mixed_highest_group = as.character(emm_df$pattern_group[[highest_index]]),
    mixed_lowest_group = as.character(emm_df$pattern_group[[lowest_index]]),
    mixed_highest_emmean = mixed_highest_emmean,
    mixed_lowest_emmean = mixed_lowest_emmean,
    mixed_standardized_range = mixed_highest_emmean - mixed_lowest_emmean,
    mixed_pattern_group_num_df = extract_model_table_value(model_anova, "pattern_group", c("NumDF", "Num Df", "numDF")),
    mixed_pattern_group_den_df = extract_model_table_value(model_anova, "pattern_group", c("DenDF", "Den Df", "denDF")),
    mixed_pattern_group_f_value = extract_model_table_value(model_anova, "pattern_group", c("F value", "F.value", "F")),
    mixed_pattern_group_p_value = extract_model_table_value(model_anova, "pattern_group", c("Pr(>F)", "Pr..F.", "p.value")),
    is_singular = tryCatch(lme4::isSingular(model_fit, tol = 1e-04), error = function(e) NA),
    stringsAsFactors = FALSE
  )

  pairwise_result <- tryCatch(
    emmeans::contrast(emm_result, method = "pairwise", adjust = "none"),
    error = function(e) e
  )

  pairwise_df <- if (inherits(pairwise_result, "error")) {
    data.frame(
      feature_name = feature_name,
      group1 = NA_character_,
      group2 = NA_character_,
      mixed_pairwise_status = "pairwise_failed",
      mixed_pairwise_message = conditionMessage(pairwise_result),
      mixed_pairwise_test_statistic = NA_character_,
      mixed_emmean_diff_group1_minus_group2 = NA_real_,
      mixed_emmean_diff_group2_minus_group1 = NA_real_,
      mixed_pairwise_se = NA_real_,
      mixed_pairwise_df = NA_real_,
      mixed_pairwise_t_ratio = NA_real_,
      mixed_pairwise_p_value = NA_real_,
      stringsAsFactors = FALSE
    )
  } else {
    contrast_df <- as.data.frame(pairwise_result)
    parsed_contrast_df <- parse_emmeans_pairwise_contrast(contrast_df$contrast)
    n_contrasts <- nrow(parsed_contrast_df)

    ## 兼容大样本下 lmerTest 停用 d.f. 计算时 emmeans 返回 z.ratio 的情况 ####
    if ("t.ratio" %in% colnames(contrast_df)) {
      pairwise_test_statistic <- "t.ratio"
      pairwise_test_ratio_values <- as.numeric(contrast_df$t.ratio)
    } else if ("z.ratio" %in% colnames(contrast_df)) {
      pairwise_test_statistic <- "z.ratio"
      pairwise_test_ratio_values <- as.numeric(contrast_df$z.ratio)
    } else {
      pairwise_test_statistic <- NA_character_
      pairwise_test_ratio_values <- rep(NA_real_, n_contrasts)
    }

    pairwise_df_values <- if ("df" %in% colnames(contrast_df)) {
      as.numeric(contrast_df$df)
    } else {
      rep(NA_real_, n_contrasts)
    }

    pairwise_p_values <- if ("p.value" %in% colnames(contrast_df)) {
      as.numeric(contrast_df$p.value)
    } else {
      rep(NA_real_, n_contrasts)
    }

    data.frame(
      feature_name = feature_name,
      group1 = parsed_contrast_df$group1,
      group2 = parsed_contrast_df$group2,
      mixed_pairwise_status = "success",
      mixed_pairwise_message = NA_character_,
      mixed_pairwise_test_statistic = pairwise_test_statistic,
      mixed_emmean_diff_group1_minus_group2 = as.numeric(contrast_df$estimate),
      mixed_emmean_diff_group2_minus_group1 = -as.numeric(contrast_df$estimate),
      mixed_pairwise_se = as.numeric(contrast_df$SE),
      mixed_pairwise_df = pairwise_df_values,
      mixed_pairwise_t_ratio = pairwise_test_ratio_values,
      mixed_pairwise_p_value = pairwise_p_values,
      stringsAsFactors = FALSE
    )
  }

  ## 构造emmean_df：每pattern_group模型调节后均值 + 95% CI(asymptotic z) ####
  mixed_emmean_df <- data.frame(
    feature_name = feature_name,
    pattern_group = as.character(emm_df$pattern_group),
    mixed_emmean_status = "success",
    mixed_emmean_message = NA_character_,
    mixed_emmean_estimate = as.numeric(emm_df$emmean),
    mixed_emmean_se = as.numeric(emm_df$SE),
    mixed_emmean_ci_lower = as.numeric(emm_df$emmean) - 1.96 * as.numeric(emm_df$SE),
    mixed_emmean_ci_upper = as.numeric(emm_df$emmean) + 1.96 * as.numeric(emm_df$SE),
    stringsAsFactors = FALSE
  )

  ## 构造one_vs_rest对比：每pattern_group vs 其它群平均 ####
  one_vs_rest_contrast_list <- build_one_vs_rest_contrast_list(
    available_group_levels = as.character(emm_df$pattern_group)
  )

  one_vs_rest_result <- if (length(one_vs_rest_contrast_list) == 0) {
    NULL
  } else {
    tryCatch(
      emmeans::contrast(emm_result, method = one_vs_rest_contrast_list, adjust = "none"),
      error = function(e) e
    )
  }

  mixed_one_vs_rest_df <- if (is.null(one_vs_rest_result) || inherits(one_vs_rest_result, "error")) {
    data.frame(
      feature_name = feature_name,
      pattern_group = NA_character_,
      mixed_one_vs_rest_status = "one_vs_rest_failed",
      mixed_one_vs_rest_message = if (inherits(one_vs_rest_result, "error")) conditionMessage(one_vs_rest_result) else "insufficient_groups",
      mixed_one_vs_rest_estimate = NA_real_,
      mixed_one_vs_rest_se = NA_real_,
      mixed_one_vs_rest_ci_lower = NA_real_,
      mixed_one_vs_rest_ci_upper = NA_real_,
      mixed_one_vs_rest_t_ratio = NA_real_,
      mixed_one_vs_rest_p_value = NA_real_,
      stringsAsFactors = FALSE
    )
  } else {
    one_vs_rest_contrast_df <- as.data.frame(one_vs_rest_result)
    n_one_vs_rest <- nrow(one_vs_rest_contrast_df)
    pattern_group_names <- sub("_vs_rest$", "", as.character(one_vs_rest_contrast_df$contrast))

    if ("t.ratio" %in% colnames(one_vs_rest_contrast_df)) {
      one_vs_rest_test_ratio_values <- as.numeric(one_vs_rest_contrast_df$t.ratio)
    } else if ("z.ratio" %in% colnames(one_vs_rest_contrast_df)) {
      one_vs_rest_test_ratio_values <- as.numeric(one_vs_rest_contrast_df$z.ratio)
    } else {
      one_vs_rest_test_ratio_values <- rep(NA_real_, n_one_vs_rest)
    }

    one_vs_rest_p_values <- if ("p.value" %in% colnames(one_vs_rest_contrast_df)) {
      as.numeric(one_vs_rest_contrast_df$p.value)
    } else {
      rep(NA_real_, n_one_vs_rest)
    }

    data.frame(
      feature_name = feature_name,
      pattern_group = pattern_group_names,
      mixed_one_vs_rest_status = "success",
      mixed_one_vs_rest_message = NA_character_,
      mixed_one_vs_rest_estimate = as.numeric(one_vs_rest_contrast_df$estimate),
      mixed_one_vs_rest_se = as.numeric(one_vs_rest_contrast_df$SE),
      mixed_one_vs_rest_ci_lower = as.numeric(one_vs_rest_contrast_df$estimate) - 1.96 * as.numeric(one_vs_rest_contrast_df$SE),
      mixed_one_vs_rest_ci_upper = as.numeric(one_vs_rest_contrast_df$estimate) + 1.96 * as.numeric(one_vs_rest_contrast_df$SE),
      mixed_one_vs_rest_t_ratio = one_vs_rest_test_ratio_values,
      mixed_one_vs_rest_p_value = one_vs_rest_p_values,
      stringsAsFactors = FALSE
    )
  }

  list(
    overall_df = overall_df,
    pairwise_df = pairwise_df,
    emmean_df = mixed_emmean_df,
    one_vs_rest_df = mixed_one_vs_rest_df
  )
}

## 函数：run_mixed_effect_models ####
## 功能：对所有TME feature执行样本级混合效应模型敏感性分析 ####
## 输入：sample_score_df，样本级免疫特征表；feature_columns，特征列名；feature_info_df，特征信息表 ####
## 输出：列表，包含混合模型整体检验表和两两比较表 ####
run_mixed_effect_models <- function(sample_score_df, feature_columns, feature_info_df) {
  mixed_result_list <- lapply(feature_columns, function(current_feature) {
    fit_single_mixed_model_feature(
      sample_score_df = sample_score_df,
      feature_name = current_feature
    )
  })

  mixed_overall_df <- dplyr::bind_rows(lapply(mixed_result_list, `[[`, "overall_df")) %>%
    dplyr::mutate(mixed_pattern_group_fdr = p.adjust(mixed_pattern_group_p_value, method = "BH")) %>%
    dplyr::left_join(feature_info_df, by = "feature_name") %>%
    dplyr::arrange(axis_order, feature_order)

  mixed_pairwise_raw_df <- dplyr::bind_rows(lapply(mixed_result_list, `[[`, "pairwise_df"))
  mixed_pairwise_df <- if (nrow(mixed_pairwise_raw_df) == 0) {
    data.frame(
      feature_name = character(),
      group1 = character(),
      group2 = character(),
      mixed_pairwise_status = character(),
      mixed_pairwise_message = character(),
      mixed_pairwise_test_statistic = character(),
      mixed_emmean_diff_group1_minus_group2 = numeric(),
      mixed_emmean_diff_group2_minus_group1 = numeric(),
      mixed_pairwise_se = numeric(),
      mixed_pairwise_df = numeric(),
      mixed_pairwise_t_ratio = numeric(),
      mixed_pairwise_p_value = numeric(),
      mixed_pairwise_fdr = numeric(),
      stringsAsFactors = FALSE
    ) %>%
      dplyr::left_join(feature_info_df, by = "feature_name")
  } else {
    mixed_pairwise_raw_df %>%
      dplyr::mutate(mixed_pairwise_fdr = p.adjust(mixed_pairwise_p_value, method = "BH")) %>%
      dplyr::left_join(feature_info_df, by = "feature_name") %>%
      dplyr::arrange(axis_order, feature_order, group1, group2)
  }

  ## 汇总emmean_df：每pattern_group模型调节后均值 + 95% CI ####
  mixed_emmean_raw_df <- dplyr::bind_rows(lapply(mixed_result_list, `[[`, "emmean_df"))
  mixed_emmean_df <- if (nrow(mixed_emmean_raw_df) == 0) {
    data.frame()
  } else {
    mixed_emmean_raw_df %>%
      dplyr::left_join(feature_info_df, by = "feature_name") %>%
      dplyr::mutate(pattern_group = factor(pattern_group, levels = pattern_group_order)) %>%
      dplyr::arrange(axis_order, feature_order, pattern_group) %>%
      dplyr::mutate(pattern_group = as.character(pattern_group))
  }

  ## 汇总one_vs_rest_df：每pattern_group vs 其它群平均对比 + FDR ####
  mixed_one_vs_rest_raw_df <- dplyr::bind_rows(lapply(mixed_result_list, `[[`, "one_vs_rest_df"))
  mixed_one_vs_rest_df <- if (nrow(mixed_one_vs_rest_raw_df) == 0) {
    data.frame()
  } else {
    mixed_one_vs_rest_raw_df %>%
      dplyr::mutate(mixed_one_vs_rest_fdr = p.adjust(mixed_one_vs_rest_p_value, method = "BH")) %>%
      dplyr::left_join(feature_info_df, by = "feature_name") %>%
      dplyr::mutate(pattern_group = factor(pattern_group, levels = pattern_group_order)) %>%
      dplyr::arrange(axis_order, feature_order, pattern_group) %>%
      dplyr::mutate(pattern_group = as.character(pattern_group))
  }

  list(
    overall_df = mixed_overall_df,
    pairwise_df = mixed_pairwise_df,
    emmean_df = mixed_emmean_df,
    one_vs_rest_df = mixed_one_vs_rest_df
  )
}

## 函数：build_direction_effect_size_summary ####
## 功能：整合观察到的方向、bootstrap稳定性、effect size和混合模型方向一致性 ####
## 输入：pattern_group_summary_df，模式群汇总表；primary_organ_score_df，原发部位级特征表；feature_info_df，特征信息表；bootstrap_feature_stability_df，bootstrap稳定性表；mixed_model_overall_df，混合模型整体检验表 ####
## 输出：方向和effect size主解释汇总表 ####
build_direction_effect_size_summary <- function(pattern_group_summary_df,
                                                primary_organ_score_df,
                                                feature_info_df,
                                                bootstrap_feature_stability_df,
                                                mixed_model_overall_df) {
  feature_order_df <- feature_info_df %>%
    dplyr::filter(main_display == TRUE) %>%
    dplyr::arrange(axis_order, feature_order)

  observed_direction_df <- dplyr::bind_rows(lapply(feature_order_df$feature_name, function(current_feature) {
    current_group_df <- pattern_group_summary_df %>%
      dplyr::filter(feature_name == current_feature, !is.na(median_score))

    current_primary_values <- primary_organ_score_df[[current_feature]]
    current_primary_values <- current_primary_values[!is.na(current_primary_values)]
    primary_organ_score_sd <- stats::sd(current_primary_values, na.rm = TRUE)

    if (nrow(current_group_df) == 0) {
      return(data.frame(
        feature_name = current_feature,
        observed_highest_group = NA_character_,
        observed_lowest_group = NA_character_,
        highest_median_score = NA_real_,
        lowest_median_score = NA_real_,
        median_range_raw = NA_real_,
        primary_organ_score_sd = primary_organ_score_sd,
        standardized_median_range = NA_real_,
        n_primary_organs_total = length(current_primary_values),
        stringsAsFactors = FALSE
      ))
    }

    highest_index <- which.max(current_group_df$median_score)
    lowest_index <- which.min(current_group_df$median_score)
    median_range_raw <- current_group_df$median_score[[highest_index]] - current_group_df$median_score[[lowest_index]]
    standardized_median_range <- if (is.na(primary_organ_score_sd) || primary_organ_score_sd == 0) {
      NA_real_
    } else {
      median_range_raw / primary_organ_score_sd
    }

    data.frame(
      feature_name = current_feature,
      observed_highest_group = as.character(current_group_df$pattern_group[[highest_index]]),
      observed_lowest_group = as.character(current_group_df$pattern_group[[lowest_index]]),
      highest_median_score = current_group_df$median_score[[highest_index]],
      lowest_median_score = current_group_df$median_score[[lowest_index]],
      median_range_raw = median_range_raw,
      primary_organ_score_sd = primary_organ_score_sd,
      standardized_median_range = standardized_median_range,
      n_primary_organs_total = length(current_primary_values),
      stringsAsFactors = FALSE
    )
  }))

  observed_direction_df %>%
    dplyr::left_join(
      feature_order_df %>%
        dplyr::select(feature_name, feature_label, feature_axis, axis_order, feature_order, feature_class, method, source, main_display),
      by = "feature_name"
    ) %>%
    dplyr::left_join(
      bootstrap_feature_stability_df %>%
        dplyr::select(
          feature_name,
          most_often_highest_group,
          highest_group_probability,
          most_often_lowest_group,
          lowest_group_probability,
          bootstrap_median_range,
          bootstrap_range_ci_lower,
          bootstrap_range_ci_upper,
          observed_p_value,
          observed_fdr
        ),
      by = "feature_name"
    ) %>%
    dplyr::left_join(
      mixed_model_overall_df %>%
        dplyr::select(
          feature_name,
          mixed_model_status = model_status,
          mixed_model_message = model_message,
          mixed_highest_group,
          mixed_lowest_group,
          mixed_highest_emmean,
          mixed_lowest_emmean,
          mixed_standardized_range,
          mixed_pattern_group_f_value,
          mixed_pattern_group_p_value,
          mixed_pattern_group_fdr,
          mixed_is_singular = is_singular
        ),
      by = "feature_name"
    ) %>%
    dplyr::mutate(
      direction_stability = vapply(highest_group_probability, classify_direction_stability, character(1)),
      effect_size_class = vapply(standardized_median_range, classify_effect_size, character(1)),
      mixed_model_direction_agreement = dplyr::case_when(
        mixed_model_status != "success" ~ "Mixed model unavailable",
        is.na(mixed_highest_group) | is.na(most_often_highest_group) ~ "Mixed model unavailable",
        mixed_highest_group == most_often_highest_group ~ "Mixed model supported",
        TRUE ~ "Mixed model discordant"
      ),
      interpretation_priority = dplyr::case_when(
        direction_stability == "Stable" & effect_size_class %in% c("Moderate", "Large") ~ "Primary interpretation",
        direction_stability == "Trend" & effect_size_class %in% c("Moderate", "Large") ~ "Secondary interpretation",
        direction_stability == "Stable" & effect_size_class == "Small" ~ "Direction only",
        TRUE ~ "Exploratory only"
      )
    ) %>%
    dplyr::arrange(axis_order, feature_order)
}

## 函数：build_weighted_mixed_model_failure_row ####
## 功能：为加权混合模型拟合失败的feature生成整体检验占位行 ####
## 输入：feature_name，目标特征名；model_data_df，模型输入数据；model_status，模型状态；model_message，状态说明 ####
## 输出：加权混合模型整体检验数据框单行 ####
build_weighted_mixed_model_failure_row <- function(feature_name, model_data_df, model_status, model_message) {
  data.frame(
    feature_name = feature_name,
    weighted_mixed_status = model_status,
    weighted_mixed_message = model_message,
    transformed_feature = ifelse(feature_name == "TMB", "log1p_then_z_score", "z_score"),
    n_samples = if (is.null(model_data_df)) 0 else nrow(model_data_df),
    n_primary_organs = if (is.null(model_data_df)) 0 else dplyr::n_distinct(model_data_df$primary_organ),
    n_pattern_groups = if (is.null(model_data_df)) 0 else dplyr::n_distinct(model_data_df$pattern_group),
    weighted_mixed_highest_group = NA_character_,
    weighted_mixed_lowest_group = NA_character_,
    weighted_mixed_highest_emmean = NA_real_,
    weighted_mixed_lowest_emmean = NA_real_,
    weighted_mixed_standardized_range = NA_real_,
    weighted_mixed_pattern_group_num_df = NA_real_,
    weighted_mixed_pattern_group_den_df = NA_real_,
    weighted_mixed_pattern_group_f_value = NA_real_,
    weighted_mixed_pattern_group_p_value = NA_real_,
    weighted_mixed_is_singular = NA,
    stringsAsFactors = FALSE
  )
}

## 函数：fit_single_weighted_mixed_model_feature ####
## 功能：对单个feature拟合每个primary_organ等权重的混合效应模型，权重=1/该primary_organ样本数 ####
## 输入：sample_score_df，样本级免疫特征表；feature_name，目标特征名 ####
## 输出：列表，包含加权混合模型整体检验和两两比较结果 ####
fit_single_weighted_mixed_model_feature <- function(sample_score_df, feature_name) {
  model_data_df <- prepare_mixed_model_feature_data(
    sample_score_df = sample_score_df,
    feature_name = feature_name
  )

  if (
    nrow(model_data_df) < 3 ||
      dplyr::n_distinct(model_data_df$pattern_group) < 2 ||
      dplyr::n_distinct(model_data_df$primary_organ) < 2 ||
      dplyr::n_distinct(model_data_df$feature_value_z) < 2
  ) {
    return(list(
      overall_df = build_weighted_mixed_model_failure_row(
        feature_name = feature_name,
        model_data_df = model_data_df,
        model_status = "insufficient_data",
        model_message = "Weighted mixed model requires at least two pattern groups, two primary organs and non-constant feature values."
      ),
      pairwise_df = data.frame(),
      emmean_df = data.frame(),
      one_vs_rest_df = data.frame()
    ))
  }

  ## 计算每个样本的等癌种权重，每个primary_organ的总权重为1 ####
  model_data_df <- model_data_df %>%
    dplyr::group_by(primary_organ) %>%
    dplyr::mutate(weight = 1 / dplyr::n()) %>%
    dplyr::ungroup() %>%
    as.data.frame()

  available_group_levels <- pattern_group_order[pattern_group_order %in% unique(as.character(model_data_df$pattern_group))]
  model_data_df$pattern_group <- factor(as.character(model_data_df$pattern_group), levels = available_group_levels)
  model_data_df$primary_organ <- factor(as.character(model_data_df$primary_organ))

  model_fit <- tryCatch(
    lmerTest::lmer(
      feature_value_z ~ pattern_group + (1 | primary_organ),
      data = model_data_df,
      weights = model_data_df$weight,
      REML = FALSE,
      control = lme4::lmerControl(optimizer = "bobyqa")
    ),
    error = function(e) e
  )

  if (inherits(model_fit, "error")) {
    return(list(
      overall_df = build_weighted_mixed_model_failure_row(
        feature_name = feature_name,
        model_data_df = model_data_df,
        model_status = "model_failed",
        model_message = conditionMessage(model_fit)
      ),
      pairwise_df = data.frame(),
      emmean_df = data.frame(),
      one_vs_rest_df = data.frame()
    ))
  }

  model_anova <- tryCatch(
    suppressMessages(as.data.frame(stats::anova(model_fit, type = 3))),
    error = function(e) data.frame()
  )

  emm_result <- tryCatch(
    emmeans::emmeans(model_fit, specs = "pattern_group", lmer.df = "asymptotic"),
    error = function(e) e
  )

  if (inherits(emm_result, "error")) {
    return(list(
      overall_df = build_weighted_mixed_model_failure_row(
        feature_name = feature_name,
        model_data_df = model_data_df,
        model_status = "emmeans_failed",
        model_message = conditionMessage(emm_result)
      ),
      pairwise_df = data.frame(),
      emmean_df = data.frame(),
      one_vs_rest_df = data.frame()
    ))
  }

  emm_df <- as.data.frame(emm_result)
  emm_df <- emm_df[!is.na(emm_df$emmean), , drop = FALSE]

  if (nrow(emm_df) == 0) {
    return(list(
      overall_df = build_weighted_mixed_model_failure_row(
        feature_name = feature_name,
        model_data_df = model_data_df,
        model_status = "emmeans_failed",
        model_message = "Estimated marginal means were all missing."
      ),
      pairwise_df = data.frame(),
      emmean_df = data.frame(),
      one_vs_rest_df = data.frame()
    ))
  }

  highest_index <- which.max(emm_df$emmean)
  lowest_index <- which.min(emm_df$emmean)
  weighted_highest_emmean <- emm_df$emmean[[highest_index]]
  weighted_lowest_emmean <- emm_df$emmean[[lowest_index]]

  overall_df <- data.frame(
    feature_name = feature_name,
    weighted_mixed_status = "success",
    weighted_mixed_message = NA_character_,
    transformed_feature = ifelse(feature_name == "TMB", "log1p_then_z_score", "z_score"),
    n_samples = nrow(model_data_df),
    n_primary_organs = dplyr::n_distinct(model_data_df$primary_organ),
    n_pattern_groups = dplyr::n_distinct(model_data_df$pattern_group),
    weighted_mixed_highest_group = as.character(emm_df$pattern_group[[highest_index]]),
    weighted_mixed_lowest_group = as.character(emm_df$pattern_group[[lowest_index]]),
    weighted_mixed_highest_emmean = weighted_highest_emmean,
    weighted_mixed_lowest_emmean = weighted_lowest_emmean,
    weighted_mixed_standardized_range = weighted_highest_emmean - weighted_lowest_emmean,
    weighted_mixed_pattern_group_num_df = extract_model_table_value(model_anova, "pattern_group", c("NumDF", "Num Df", "numDF")),
    weighted_mixed_pattern_group_den_df = extract_model_table_value(model_anova, "pattern_group", c("DenDF", "Den Df", "denDF")),
    weighted_mixed_pattern_group_f_value = extract_model_table_value(model_anova, "pattern_group", c("F value", "F.value", "F")),
    weighted_mixed_pattern_group_p_value = extract_model_table_value(model_anova, "pattern_group", c("Pr(>F)", "Pr..F.", "p.value")),
    weighted_mixed_is_singular = tryCatch(lme4::isSingular(model_fit, tol = 1e-04), error = function(e) NA),
    stringsAsFactors = FALSE
  )

  pairwise_result <- tryCatch(
    emmeans::contrast(emm_result, method = "pairwise", adjust = "none"),
    error = function(e) e
  )

  pairwise_df <- if (inherits(pairwise_result, "error")) {
    data.frame(
      feature_name = feature_name,
      group1 = NA_character_,
      group2 = NA_character_,
      weighted_pairwise_status = "pairwise_failed",
      weighted_pairwise_message = conditionMessage(pairwise_result),
      weighted_pairwise_test_statistic = NA_character_,
      weighted_emmean_diff_group1_minus_group2 = NA_real_,
      weighted_emmean_diff_group2_minus_group1 = NA_real_,
      weighted_pairwise_se = NA_real_,
      weighted_pairwise_df = NA_real_,
      weighted_pairwise_t_ratio = NA_real_,
      weighted_pairwise_p_value = NA_real_,
      stringsAsFactors = FALSE
    )
  } else {
    contrast_df <- as.data.frame(pairwise_result)
    parsed_contrast_df <- parse_emmeans_pairwise_contrast(contrast_df$contrast)
    n_contrasts <- nrow(parsed_contrast_df)

    if ("t.ratio" %in% colnames(contrast_df)) {
      pairwise_test_statistic <- "t.ratio"
      pairwise_test_ratio_values <- as.numeric(contrast_df$t.ratio)
    } else if ("z.ratio" %in% colnames(contrast_df)) {
      pairwise_test_statistic <- "z.ratio"
      pairwise_test_ratio_values <- as.numeric(contrast_df$z.ratio)
    } else {
      pairwise_test_statistic <- NA_character_
      pairwise_test_ratio_values <- rep(NA_real_, n_contrasts)
    }

    pairwise_df_values <- if ("df" %in% colnames(contrast_df)) {
      as.numeric(contrast_df$df)
    } else {
      rep(NA_real_, n_contrasts)
    }

    pairwise_p_values <- if ("p.value" %in% colnames(contrast_df)) {
      as.numeric(contrast_df$p.value)
    } else {
      rep(NA_real_, n_contrasts)
    }

    data.frame(
      feature_name = feature_name,
      group1 = parsed_contrast_df$group1,
      group2 = parsed_contrast_df$group2,
      weighted_pairwise_status = "success",
      weighted_pairwise_message = NA_character_,
      weighted_pairwise_test_statistic = pairwise_test_statistic,
      weighted_emmean_diff_group1_minus_group2 = as.numeric(contrast_df$estimate),
      weighted_emmean_diff_group2_minus_group1 = -as.numeric(contrast_df$estimate),
      weighted_pairwise_se = as.numeric(contrast_df$SE),
      weighted_pairwise_df = pairwise_df_values,
      weighted_pairwise_t_ratio = pairwise_test_ratio_values,
      weighted_pairwise_p_value = pairwise_p_values,
      stringsAsFactors = FALSE
    )
  }

  ## 构造emmean_df：每pattern_group模型调节后均值 + 95% CI(asymptotic z)，加权LMM下每癌种等权 ####
  weighted_emmean_df <- data.frame(
    feature_name = feature_name,
    pattern_group = as.character(emm_df$pattern_group),
    weighted_emmean_status = "success",
    weighted_emmean_message = NA_character_,
    weighted_emmean_estimate = as.numeric(emm_df$emmean),
    weighted_emmean_se = as.numeric(emm_df$SE),
    weighted_emmean_ci_lower = as.numeric(emm_df$emmean) - 1.96 * as.numeric(emm_df$SE),
    weighted_emmean_ci_upper = as.numeric(emm_df$emmean) + 1.96 * as.numeric(emm_df$SE),
    stringsAsFactors = FALSE
  )

  ## 构造one_vs_rest对比：每pattern_group vs 其它群平均（加权） ####
  one_vs_rest_contrast_list <- build_one_vs_rest_contrast_list(
    available_group_levels = as.character(emm_df$pattern_group)
  )

  one_vs_rest_result <- if (length(one_vs_rest_contrast_list) == 0) {
    NULL
  } else {
    tryCatch(
      emmeans::contrast(emm_result, method = one_vs_rest_contrast_list, adjust = "none"),
      error = function(e) e
    )
  }

  weighted_one_vs_rest_df <- if (is.null(one_vs_rest_result) || inherits(one_vs_rest_result, "error")) {
    data.frame(
      feature_name = feature_name,
      pattern_group = NA_character_,
      weighted_one_vs_rest_status = "one_vs_rest_failed",
      weighted_one_vs_rest_message = if (inherits(one_vs_rest_result, "error")) conditionMessage(one_vs_rest_result) else "insufficient_groups",
      weighted_one_vs_rest_estimate = NA_real_,
      weighted_one_vs_rest_se = NA_real_,
      weighted_one_vs_rest_ci_lower = NA_real_,
      weighted_one_vs_rest_ci_upper = NA_real_,
      weighted_one_vs_rest_t_ratio = NA_real_,
      weighted_one_vs_rest_p_value = NA_real_,
      stringsAsFactors = FALSE
    )
  } else {
    one_vs_rest_contrast_df <- as.data.frame(one_vs_rest_result)
    n_one_vs_rest <- nrow(one_vs_rest_contrast_df)
    pattern_group_names <- sub("_vs_rest$", "", as.character(one_vs_rest_contrast_df$contrast))

    if ("t.ratio" %in% colnames(one_vs_rest_contrast_df)) {
      one_vs_rest_test_ratio_values <- as.numeric(one_vs_rest_contrast_df$t.ratio)
    } else if ("z.ratio" %in% colnames(one_vs_rest_contrast_df)) {
      one_vs_rest_test_ratio_values <- as.numeric(one_vs_rest_contrast_df$z.ratio)
    } else {
      one_vs_rest_test_ratio_values <- rep(NA_real_, n_one_vs_rest)
    }

    one_vs_rest_p_values <- if ("p.value" %in% colnames(one_vs_rest_contrast_df)) {
      as.numeric(one_vs_rest_contrast_df$p.value)
    } else {
      rep(NA_real_, n_one_vs_rest)
    }

    data.frame(
      feature_name = feature_name,
      pattern_group = pattern_group_names,
      weighted_one_vs_rest_status = "success",
      weighted_one_vs_rest_message = NA_character_,
      weighted_one_vs_rest_estimate = as.numeric(one_vs_rest_contrast_df$estimate),
      weighted_one_vs_rest_se = as.numeric(one_vs_rest_contrast_df$SE),
      weighted_one_vs_rest_ci_lower = as.numeric(one_vs_rest_contrast_df$estimate) - 1.96 * as.numeric(one_vs_rest_contrast_df$SE),
      weighted_one_vs_rest_ci_upper = as.numeric(one_vs_rest_contrast_df$estimate) + 1.96 * as.numeric(one_vs_rest_contrast_df$SE),
      weighted_one_vs_rest_t_ratio = one_vs_rest_test_ratio_values,
      weighted_one_vs_rest_p_value = one_vs_rest_p_values,
      stringsAsFactors = FALSE
    )
  }

  list(
    overall_df = overall_df,
    pairwise_df = pairwise_df,
    emmean_df = weighted_emmean_df,
    one_vs_rest_df = weighted_one_vs_rest_df
  )
}

## 函数：run_weighted_mixed_effect_models ####
## 功能：对所有TME feature执行每癌种等权重的混合效应模型 ####
## 输入：sample_score_df，样本级免疫特征表；feature_columns，特征列名；feature_info_df，特征信息表 ####
## 输出：列表，包含加权混合模型整体检验表和两两比较表 ####
run_weighted_mixed_effect_models <- function(sample_score_df, feature_columns, feature_info_df) {
  weighted_result_list <- lapply(feature_columns, function(current_feature) {
    fit_single_weighted_mixed_model_feature(
      sample_score_df = sample_score_df,
      feature_name = current_feature
    )
  })

  weighted_overall_df <- dplyr::bind_rows(lapply(weighted_result_list, `[[`, "overall_df")) %>%
    dplyr::mutate(weighted_mixed_pattern_group_fdr = p.adjust(weighted_mixed_pattern_group_p_value, method = "BH")) %>%
    dplyr::left_join(feature_info_df, by = "feature_name") %>%
    dplyr::arrange(axis_order, feature_order)

  weighted_pairwise_raw_df <- dplyr::bind_rows(lapply(weighted_result_list, `[[`, "pairwise_df"))
  weighted_pairwise_df <- if (nrow(weighted_pairwise_raw_df) == 0) {
    data.frame(
      feature_name = character(),
      group1 = character(),
      group2 = character(),
      weighted_pairwise_status = character(),
      weighted_pairwise_message = character(),
      weighted_pairwise_test_statistic = character(),
      weighted_emmean_diff_group1_minus_group2 = numeric(),
      weighted_emmean_diff_group2_minus_group1 = numeric(),
      weighted_pairwise_se = numeric(),
      weighted_pairwise_df = numeric(),
      weighted_pairwise_t_ratio = numeric(),
      weighted_pairwise_p_value = numeric(),
      weighted_pairwise_fdr = numeric(),
      stringsAsFactors = FALSE
    ) %>%
      dplyr::left_join(feature_info_df, by = "feature_name")
  } else {
    weighted_pairwise_raw_df %>%
      dplyr::mutate(weighted_pairwise_fdr = p.adjust(weighted_pairwise_p_value, method = "BH")) %>%
      dplyr::left_join(feature_info_df, by = "feature_name") %>%
      dplyr::arrange(axis_order, feature_order, group1, group2)
  }

  ## 汇总emmean_df：每pattern_group模型调节后均值 + 95% CI（加权LMM） ####
  weighted_emmean_raw_df <- dplyr::bind_rows(lapply(weighted_result_list, `[[`, "emmean_df"))
  weighted_emmean_df <- if (nrow(weighted_emmean_raw_df) == 0) {
    data.frame()
  } else {
    weighted_emmean_raw_df %>%
      dplyr::left_join(feature_info_df, by = "feature_name") %>%
      dplyr::mutate(pattern_group = factor(pattern_group, levels = pattern_group_order)) %>%
      dplyr::arrange(axis_order, feature_order, pattern_group) %>%
      dplyr::mutate(pattern_group = as.character(pattern_group))
  }

  ## 汇总one_vs_rest_df：每pattern_group vs 其它群平均 + FDR（加权LMM） ####
  weighted_one_vs_rest_raw_df <- dplyr::bind_rows(lapply(weighted_result_list, `[[`, "one_vs_rest_df"))
  weighted_one_vs_rest_df <- if (nrow(weighted_one_vs_rest_raw_df) == 0) {
    data.frame()
  } else {
    weighted_one_vs_rest_raw_df %>%
      dplyr::mutate(weighted_one_vs_rest_fdr = p.adjust(weighted_one_vs_rest_p_value, method = "BH")) %>%
      dplyr::left_join(feature_info_df, by = "feature_name") %>%
      dplyr::mutate(pattern_group = factor(pattern_group, levels = pattern_group_order)) %>%
      dplyr::arrange(axis_order, feature_order, pattern_group) %>%
      dplyr::mutate(pattern_group = as.character(pattern_group))
  }

  list(
    overall_df = weighted_overall_df,
    pairwise_df = weighted_pairwise_df,
    emmean_df = weighted_emmean_df,
    one_vs_rest_df = weighted_one_vs_rest_df
  )
}

## 函数：subsample_balanced_per_organ ####
## 功能：对每个primary_organ随机抽取相同样本量，得到等癌种平衡的样本子集 ####
## 输入：sample_score_df，样本级免疫特征表；n_per_organ，每个organ抽取的样本量上限 ####
## 输出：抽样后的样本级数据框，每个primary_organ最多保留n_per_organ条 ####
subsample_balanced_per_organ <- function(sample_score_df, n_per_organ) {
  sample_score_df %>%
    dplyr::filter(!is.na(primary_organ), primary_organ != "") %>%
    dplyr::group_by(primary_organ) %>%
    dplyr::group_modify(~ {
      n_avail <- nrow(.x)
      if (n_avail == 0) {
        return(.x[FALSE, , drop = FALSE])
      }
      n_take <- min(n_avail, n_per_organ)
      .x[sample.int(n_avail, n_take), , drop = FALSE]
    }) %>%
    dplyr::ungroup() %>%
    as.data.frame()
}

## 函数：run_stratified_subsample_mixed_models ####
## 功能：分层抽样平衡每个癌种样本量后多次拟合混合模型，输出每次迭代的整体、两两、BLUP、emmean和one-vs-rest结果 ####
## 输入：sample_score_df，样本级特征表；feature_columns，特征列；n_iterations，迭代次数；n_per_organ，每个organ抽取样本量；seed，随机种子；use_parallel，是否并行；workers，并行核数 ####
## 输出：列表，包含iteration_overall_df、iteration_pairwise_df、iteration_blup_df、iteration_emmean_df和iteration_one_vs_rest_df五张迭代级结果表 ####
run_stratified_subsample_mixed_models <- function(sample_score_df,
                                                  feature_columns,
                                                  n_iterations,
                                                  n_per_organ,
                                                  seed,
                                                  use_parallel,
                                                  workers) {
  iter_result_list <- run_parallel_lapply(
    x = seq_len(n_iterations),
    fun = function(iter_index) {
      set.seed(seed + iter_index)
      subsample_df <- subsample_balanced_per_organ(
        sample_score_df = sample_score_df,
        n_per_organ = n_per_organ
      )

      iter_results <- lapply(feature_columns, function(current_feature) {
        fit_single_mixed_model_feature(
          sample_score_df = subsample_df,
          feature_name = current_feature
        )
      })

      iter_overall_df <- dplyr::bind_rows(lapply(iter_results, `[[`, "overall_df")) %>%
        dplyr::mutate(iteration = iter_index)
      iter_pairwise_df <- dplyr::bind_rows(lapply(iter_results, `[[`, "pairwise_df"))
      if (nrow(iter_pairwise_df) > 0) {
        iter_pairwise_df$iteration <- iter_index
      }
      iter_emmean_df <- dplyr::bind_rows(lapply(iter_results, `[[`, "emmean_df"))
      if (nrow(iter_emmean_df) > 0) {
        iter_emmean_df$iteration <- iter_index
      }
      iter_one_vs_rest_df <- dplyr::bind_rows(lapply(iter_results, `[[`, "one_vs_rest_df"))
      if (nrow(iter_one_vs_rest_df) > 0) {
        iter_one_vs_rest_df$iteration <- iter_index
      }

      iter_blup_list <- lapply(feature_columns, function(current_feature) {
        extract_organ_blup_for_feature(
          sample_score_df = subsample_df,
          feature_name = current_feature,
          use_equal_organ_weights = FALSE
        )
      })
      iter_blup_df <- dplyr::bind_rows(iter_blup_list)
      if (nrow(iter_blup_df) > 0) {
        iter_blup_df$iteration <- iter_index
      }

      list(
        overall_df = iter_overall_df,
        pairwise_df = iter_pairwise_df,
        blup_df = iter_blup_df,
        emmean_df = iter_emmean_df,
        one_vs_rest_df = iter_one_vs_rest_df
      )
    },
    use_parallel = use_parallel,
    workers = workers
  )

  iteration_overall_df <- dplyr::bind_rows(lapply(iter_result_list, `[[`, "overall_df"))
  iteration_pairwise_df <- dplyr::bind_rows(lapply(iter_result_list, `[[`, "pairwise_df"))
  iteration_blup_df <- dplyr::bind_rows(lapply(iter_result_list, `[[`, "blup_df"))
  iteration_emmean_df <- dplyr::bind_rows(lapply(iter_result_list, `[[`, "emmean_df"))
  iteration_one_vs_rest_df <- dplyr::bind_rows(lapply(iter_result_list, `[[`, "one_vs_rest_df"))

  list(
    iteration_overall_df = iteration_overall_df,
    iteration_pairwise_df = iteration_pairwise_df,
    iteration_blup_df = iteration_blup_df,
    iteration_emmean_df = iteration_emmean_df,
    iteration_one_vs_rest_df = iteration_one_vs_rest_df
  )
}

## 函数：summarise_subsample_overall ####
## 功能：聚合多次分层抽样混合模型的整体检验结果，得到每个feature的稳健估计 ####
## 输入：iteration_overall_df，每次迭代整体检验结果；feature_info_df，特征信息表 ####
## 输出：每个feature的聚合整体检验数据框，包含%显著、中位p、中位effect size、95%CI、最常最高/最低组及概率 ####
summarise_subsample_overall <- function(iteration_overall_df, feature_info_df) {
  if (is.null(iteration_overall_df) || nrow(iteration_overall_df) == 0) {
    return(data.frame())
  }

  success_df <- iteration_overall_df %>%
    dplyr::filter(model_status == "success")

  iter_total_df <- iteration_overall_df %>%
    dplyr::group_by(feature_name) %>%
    dplyr::summarise(subsample_n_iterations_total = dplyr::n(), .groups = "drop")

  iter_success_df <- success_df %>%
    dplyr::group_by(feature_name) %>%
    dplyr::summarise(
      subsample_n_iterations_success = dplyr::n(),
      subsample_pct_p_below_005 = mean(mixed_pattern_group_p_value < 0.05, na.rm = TRUE),
      subsample_median_p_value = stats::median(mixed_pattern_group_p_value, na.rm = TRUE),
      subsample_median_standardized_range = stats::median(mixed_standardized_range, na.rm = TRUE),
      subsample_range_ci_lower = stats::quantile(mixed_standardized_range, 0.025, na.rm = TRUE),
      subsample_range_ci_upper = stats::quantile(mixed_standardized_range, 0.975, na.rm = TRUE),
      .groups = "drop"
    )

  highest_summary_df <- success_df %>%
    dplyr::filter(!is.na(mixed_highest_group)) %>%
    dplyr::count(feature_name, mixed_highest_group, name = "n_highest") %>%
    dplyr::group_by(feature_name) %>%
    dplyr::mutate(highest_total = sum(n_highest)) %>%
    dplyr::slice_max(n_highest, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      feature_name = feature_name,
      subsample_most_often_highest_group = mixed_highest_group,
      subsample_highest_group_probability = n_highest / highest_total
    )

  lowest_summary_df <- success_df %>%
    dplyr::filter(!is.na(mixed_lowest_group)) %>%
    dplyr::count(feature_name, mixed_lowest_group, name = "n_lowest") %>%
    dplyr::group_by(feature_name) %>%
    dplyr::mutate(lowest_total = sum(n_lowest)) %>%
    dplyr::slice_max(n_lowest, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      feature_name = feature_name,
      subsample_most_often_lowest_group = mixed_lowest_group,
      subsample_lowest_group_probability = n_lowest / lowest_total
    )

  iter_total_df %>%
    dplyr::left_join(iter_success_df, by = "feature_name") %>%
    dplyr::left_join(highest_summary_df, by = "feature_name") %>%
    dplyr::left_join(lowest_summary_df, by = "feature_name") %>%
    dplyr::left_join(feature_info_df, by = "feature_name") %>%
    dplyr::arrange(axis_order, feature_order)
}

## 函数：summarise_subsample_pairwise ####
## 功能：聚合多次分层抽样混合模型的两两比较结果，得到每个feature×pair的稳健估计 ####
## 输入：iteration_pairwise_df，每次迭代两两结果；feature_info_df，特征信息表 ####
## 输出：每个feature×pair的聚合两两比较数据框，包含%显著、中位p、中位estimate、95%CI、方向一致性 ####
summarise_subsample_pairwise <- function(iteration_pairwise_df, feature_info_df) {
  if (is.null(iteration_pairwise_df) || nrow(iteration_pairwise_df) == 0) {
    return(data.frame())
  }

  success_df <- iteration_pairwise_df %>%
    dplyr::filter(mixed_pairwise_status == "success")

  if (nrow(success_df) == 0) {
    return(data.frame())
  }

  success_df %>%
    dplyr::group_by(feature_name, group1, group2) %>%
    dplyr::summarise(
      subsample_n_iterations_success = dplyr::n(),
      subsample_pairwise_pct_p_below_005 = mean(mixed_pairwise_p_value < 0.05, na.rm = TRUE),
      subsample_pairwise_median_p_value = stats::median(mixed_pairwise_p_value, na.rm = TRUE),
      subsample_pairwise_median_estimate_g2_minus_g1 = stats::median(mixed_emmean_diff_group2_minus_group1, na.rm = TRUE),
      subsample_pairwise_estimate_ci_lower = stats::quantile(mixed_emmean_diff_group2_minus_group1, 0.025, na.rm = TRUE),
      subsample_pairwise_estimate_ci_upper = stats::quantile(mixed_emmean_diff_group2_minus_group1, 0.975, na.rm = TRUE),
      subsample_pairwise_direction_consistency = max(
        mean(mixed_emmean_diff_group2_minus_group1 > 0, na.rm = TRUE),
        mean(mixed_emmean_diff_group2_minus_group1 < 0, na.rm = TRUE)
      ),
      .groups = "drop"
    ) %>%
    dplyr::left_join(feature_info_df, by = "feature_name") %>%
    dplyr::arrange(axis_order, feature_order, group1, group2)
}

## 函数：extract_organ_blup_for_feature ####
## 功能：对单个feature拟合LMM(可选加权)并提取每个primary_organ的随机截距BLUP和条件SE ####
## 输入：sample_score_df，样本级特征表；feature_name，目标特征名；use_equal_organ_weights，是否使用加权LMM (weights = 1 / n_organ) ####
## 输出：每个primary_organ一行的BLUP数据框，含点估计、SE和95%CI；模型失败时返回空数据框 ####
extract_organ_blup_for_feature <- function(sample_score_df, feature_name, use_equal_organ_weights = FALSE) {
  model_data_df <- prepare_mixed_model_feature_data(
    sample_score_df = sample_score_df,
    feature_name = feature_name
  )

  if (
    nrow(model_data_df) < 3 ||
      dplyr::n_distinct(model_data_df$pattern_group) < 2 ||
      dplyr::n_distinct(model_data_df$primary_organ) < 2 ||
      dplyr::n_distinct(model_data_df$feature_value_z) < 2
  ) {
    return(data.frame())
  }

  if (isTRUE(use_equal_organ_weights)) {
    model_data_df <- model_data_df %>%
      dplyr::group_by(primary_organ) %>%
      dplyr::mutate(weight = 1 / dplyr::n()) %>%
      dplyr::ungroup() %>%
      as.data.frame()
  }

  available_group_levels <- pattern_group_order[pattern_group_order %in% unique(as.character(model_data_df$pattern_group))]
  model_data_df$pattern_group <- factor(as.character(model_data_df$pattern_group), levels = available_group_levels)
  model_data_df$primary_organ <- factor(as.character(model_data_df$primary_organ))

  model_fit <- tryCatch(
    if (isTRUE(use_equal_organ_weights)) {
      lmerTest::lmer(
        feature_value_z ~ pattern_group + (1 | primary_organ),
        data = model_data_df,
        weights = model_data_df$weight,
        REML = FALSE,
        control = lme4::lmerControl(optimizer = "bobyqa")
      )
    } else {
      lmerTest::lmer(
        feature_value_z ~ pattern_group + (1 | primary_organ),
        data = model_data_df,
        REML = FALSE,
        control = lme4::lmerControl(optimizer = "bobyqa")
      )
    },
    error = function(e) e
  )

  if (inherits(model_fit, "error")) {
    return(data.frame())
  }

  ranef_result <- tryCatch(
    lme4::ranef(model_fit, condVar = TRUE),
    error = function(e) e
  )

  if (inherits(ranef_result, "error") || is.null(ranef_result$primary_organ)) {
    return(data.frame())
  }

  organ_blup_df <- ranef_result$primary_organ
  post_var <- attr(organ_blup_df, "postVar")
  organ_se <- sqrt(as.numeric(c(post_var)))
  organ_blup <- organ_blup_df[, 1]
  organ_names <- rownames(organ_blup_df)

  organ_meta_df <- model_data_df %>%
    dplyr::group_by(primary_organ, pattern_group) %>%
    dplyr::summarise(n_samples_in_organ = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(
      primary_organ = as.character(primary_organ),
      pattern_group = as.character(pattern_group)
    )

  data.frame(
    feature_name = feature_name,
    primary_organ = organ_names,
    blup_estimate = organ_blup,
    blup_se = organ_se,
    blup_ci_lower = organ_blup - 1.96 * organ_se,
    blup_ci_upper = organ_blup + 1.96 * organ_se,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::left_join(organ_meta_df, by = "primary_organ")
}

## 函数：extract_lmm_organ_random_intercepts ####
## 功能：对所有TME feature提取LMM(可选加权)中primary_organ随机截距的BLUP及95%CI，便于绘制organ级森林图 ####
## 输入：sample_score_df，样本级特征表；feature_columns，特征列；feature_info_df，特征信息表；use_equal_organ_weights，是否使用加权LMM ####
## 输出：每个feature×organ一行的数据框，含BLUP估计、SE、95%CI、organ样本数和pattern_group ####
extract_lmm_organ_random_intercepts <- function(sample_score_df, feature_columns, feature_info_df, use_equal_organ_weights = FALSE) {
  blup_list <- lapply(feature_columns, function(current_feature) {
    extract_organ_blup_for_feature(
      sample_score_df = sample_score_df,
      feature_name = current_feature,
      use_equal_organ_weights = use_equal_organ_weights
    )
  })

  blup_df <- dplyr::bind_rows(blup_list)
  if (nrow(blup_df) == 0) {
    return(data.frame())
  }

  blup_df %>%
    dplyr::left_join(feature_info_df, by = "feature_name") %>%
    dplyr::mutate(
      pattern_group = factor(pattern_group, levels = pattern_group_order)
    ) %>%
    dplyr::arrange(axis_order, feature_order, pattern_group, primary_organ) %>%
    dplyr::mutate(pattern_group = as.character(pattern_group))
}

## 函数：summarise_subsample_blup ####
## 功能：聚合分层抽样LMM每次迭代提取的BLUP，按feature×primary_organ汇总median和95% quantile CI ####
## 输入：iteration_blup_df，每次迭代各feature各organ的BLUP估计；feature_info_df，特征信息表 ####
## 输出：每个feature×organ一行的聚合BLUP数据框，列结构与extract_lmm_organ_random_intercepts兼容 ####
summarise_subsample_blup <- function(iteration_blup_df, feature_info_df) {
  if (is.null(iteration_blup_df) || nrow(iteration_blup_df) == 0) {
    return(data.frame())
  }

  aggregated_df <- iteration_blup_df %>%
    dplyr::filter(!is.na(blup_estimate)) %>%
    dplyr::group_by(feature_name, primary_organ, pattern_group) %>%
    dplyr::summarise(
      n_iterations_covered = dplyr::n(),
      blup_estimate_median = stats::median(blup_estimate, na.rm = TRUE),
      blup_se_median = stats::median(blup_se, na.rm = TRUE),
      blup_ci_lower_q = stats::quantile(blup_estimate, 0.025, na.rm = TRUE),
      blup_ci_upper_q = stats::quantile(blup_estimate, 0.975, na.rm = TRUE),
      n_samples_in_organ_median = stats::median(n_samples_in_organ, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::transmute(
      feature_name = feature_name,
      primary_organ = primary_organ,
      pattern_group = pattern_group,
      n_iterations_covered = n_iterations_covered,
      blup_estimate = blup_estimate_median,
      blup_se = blup_se_median,
      blup_ci_lower = blup_ci_lower_q,
      blup_ci_upper = blup_ci_upper_q,
      n_samples_in_organ = n_samples_in_organ_median
    )

  aggregated_df %>%
    dplyr::left_join(feature_info_df, by = "feature_name") %>%
    dplyr::mutate(
      pattern_group = factor(pattern_group, levels = pattern_group_order)
    ) %>%
    dplyr::arrange(axis_order, feature_order, pattern_group, primary_organ) %>%
    dplyr::mutate(pattern_group = as.character(pattern_group))
}

## 函数：summarise_subsample_emmean ####
## 功能：聚合分层抽样LMM每次迭代各pattern_group的模型调节后均值(emmean)，按feature×pattern_group汇总median和95% quantile CI ####
## 输入：iteration_emmean_df，每次迭代各feature各pattern_group的emmean估计；feature_info_df，特征信息表 ####
## 输出：每个feature×pattern_group一行的聚合emmean数据框 ####
summarise_subsample_emmean <- function(iteration_emmean_df, feature_info_df) {
  if (is.null(iteration_emmean_df) || nrow(iteration_emmean_df) == 0) {
    return(data.frame())
  }

  success_df <- iteration_emmean_df %>%
    dplyr::filter(mixed_emmean_status == "success", !is.na(mixed_emmean_estimate))

  if (nrow(success_df) == 0) {
    return(data.frame())
  }

  success_df %>%
    dplyr::group_by(feature_name, pattern_group) %>%
    dplyr::summarise(
      subsample_emmean_n_iterations_success = dplyr::n(),
      subsample_emmean_median_estimate = stats::median(mixed_emmean_estimate, na.rm = TRUE),
      subsample_emmean_median_se = stats::median(mixed_emmean_se, na.rm = TRUE),
      subsample_emmean_ci_lower = stats::quantile(mixed_emmean_estimate, 0.025, na.rm = TRUE),
      subsample_emmean_ci_upper = stats::quantile(mixed_emmean_estimate, 0.975, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::left_join(feature_info_df, by = "feature_name") %>%
    dplyr::mutate(pattern_group = factor(pattern_group, levels = pattern_group_order)) %>%
    dplyr::arrange(axis_order, feature_order, pattern_group) %>%
    dplyr::mutate(pattern_group = as.character(pattern_group))
}

## 函数：summarise_subsample_one_vs_rest ####
## 功能：聚合分层抽样LMM每次迭代one-vs-rest对比结果，按feature×pattern_group汇总median、95% quantile CI、方向一致性和中位p ####
## 输入：iteration_one_vs_rest_df，每次迭代各feature各pattern_group的one-vs-rest对比；feature_info_df，特征信息表 ####
## 输出：每个feature×pattern_group一行的聚合one-vs-rest数据框 ####
summarise_subsample_one_vs_rest <- function(iteration_one_vs_rest_df, feature_info_df) {
  if (is.null(iteration_one_vs_rest_df) || nrow(iteration_one_vs_rest_df) == 0) {
    return(data.frame())
  }

  success_df <- iteration_one_vs_rest_df %>%
    dplyr::filter(mixed_one_vs_rest_status == "success", !is.na(mixed_one_vs_rest_estimate))

  if (nrow(success_df) == 0) {
    return(data.frame())
  }

  success_df %>%
    dplyr::group_by(feature_name, pattern_group) %>%
    dplyr::summarise(
      subsample_one_vs_rest_n_iterations_success = dplyr::n(),
      subsample_one_vs_rest_pct_p_below_005 = mean(mixed_one_vs_rest_p_value < 0.05, na.rm = TRUE),
      subsample_one_vs_rest_median_p_value = stats::median(mixed_one_vs_rest_p_value, na.rm = TRUE),
      subsample_one_vs_rest_median_estimate = stats::median(mixed_one_vs_rest_estimate, na.rm = TRUE),
      subsample_one_vs_rest_ci_lower = stats::quantile(mixed_one_vs_rest_estimate, 0.025, na.rm = TRUE),
      subsample_one_vs_rest_ci_upper = stats::quantile(mixed_one_vs_rest_estimate, 0.975, na.rm = TRUE),
      subsample_one_vs_rest_direction_consistency = max(
        mean(mixed_one_vs_rest_estimate > 0, na.rm = TRUE),
        mean(mixed_one_vs_rest_estimate < 0, na.rm = TRUE)
      ),
      .groups = "drop"
    ) %>%
    dplyr::left_join(feature_info_df, by = "feature_name") %>%
    dplyr::mutate(pattern_group = factor(pattern_group, levels = pattern_group_order)) %>%
    dplyr::arrange(axis_order, feature_order, pattern_group) %>%
    dplyr::mutate(pattern_group = as.character(pattern_group))
}


## 函数：load_tcga_patient_age ####
## 功能：按project_id批量读取 TCGA_GTEx/by_project/{project}_clinical_patient.csv.gz，提取每位患者的 age_at_index ####
## 输入：project_ids，TCGA项目ID向量（如"TCGA-BRCA"）；clinical_dir，临床数据根目录 ####
## 输出：两列数据框 patient_id (= bcr_patient_barcode) 和 age_at_index (数值) ####
load_tcga_patient_age <- function(project_ids, clinical_dir) {
  project_ids <- unique(as.character(project_ids))
  age_list <- lapply(project_ids, function(current_project) {
    clinical_file <- file.path(clinical_dir, paste0(current_project, "_clinical_patient.csv.gz"))
    if (!file.exists(clinical_file)) {
      warning(paste0("临床文件不存在，跳过: ", clinical_file))
      return(data.frame(patient_id = character(), age_at_index = numeric(), stringsAsFactors = FALSE))
    }

    current_df <- tryCatch(
      utils::read.csv(gzfile(clinical_file), stringsAsFactors = FALSE),
      error = function(e) {
        warning(paste0("读取临床文件失败: ", clinical_file, " - ", conditionMessage(e)))
        NULL
      }
    )
    if (is.null(current_df)) {
      return(data.frame(patient_id = character(), age_at_index = numeric(), stringsAsFactors = FALSE))
    }

    key_col <- if ("bcr_patient_barcode" %in% colnames(current_df)) {
      "bcr_patient_barcode"
    } else if ("submitter_id" %in% colnames(current_df)) {
      "submitter_id"
    } else {
      warning(paste0("临床文件缺少patient key列: ", clinical_file))
      return(data.frame(patient_id = character(), age_at_index = numeric(), stringsAsFactors = FALSE))
    }

    if (!"age_at_index" %in% colnames(current_df)) {
      warning(paste0("临床文件缺少age_at_index列: ", clinical_file))
      return(data.frame(patient_id = character(), age_at_index = numeric(), stringsAsFactors = FALSE))
    }

    data.frame(
      patient_id = as.character(current_df[[key_col]]),
      age_at_index = suppressWarnings(as.numeric(current_df$age_at_index)),
      stringsAsFactors = FALSE
    )
  })

  dplyr::bind_rows(age_list) %>%
    dplyr::filter(!is.na(patient_id), nzchar(patient_id)) %>%
    dplyr::distinct(patient_id, .keep_all = TRUE)
}


## 函数：attach_patient_age_to_sample_df ####
## 功能：把患者年龄合并到样本级特征表；若某sample缺失age_at_index，用其primary_organ内中位数填充 ####
## 输入：sample_score_df，样本级特征表（需含patient_id和primary_organ）；patient_age_df，患者年龄表 ####
## 输出：原表追加列 age_at_index 和 age_imputed_flag 后的数据框 ####
attach_patient_age_to_sample_df <- function(sample_score_df, patient_age_df) {
  merged_df <- sample_score_df %>%
    dplyr::left_join(patient_age_df, by = "patient_id")

  merged_df <- merged_df %>%
    dplyr::group_by(primary_organ) %>%
    dplyr::mutate(
      age_imputed_flag = is.na(age_at_index),
      age_at_index = dplyr::if_else(
        is.na(age_at_index),
        stats::median(age_at_index, na.rm = TRUE),
        age_at_index
      )
    ) %>%
    dplyr::ungroup() %>%
    as.data.frame()

  ## 极端情况：某organ内全部缺失，用整体中位数填充 ####
  global_median_age <- stats::median(merged_df$age_at_index, na.rm = TRUE)
  merged_df$age_at_index[is.na(merged_df$age_at_index)] <- global_median_age

  merged_df
}


## 函数：build_cem_matched_sample_set ####
## 功能：对每个primary_organ在 age_bin 分层结构上做CEM抽样，使各organ样本量相等且 age 分布接近共同目标 ####
## 输入：sample_score_df，含age_at_index的样本表；n_per_organ，每organ目标样本量；age_breaks，age 分箱切点；seed，随机种子 ####
## 输出：列表，包含matched_sample_df（匹配后的样本级表）、balance_df（匹配前后 age 分布对比） ####
build_cem_matched_sample_set <- function(sample_score_df,
                                         n_per_organ,
                                         age_breaks,
                                         seed) {
  set.seed(seed)

  work_df <- sample_score_df %>%
    dplyr::filter(!is.na(age_at_index)) %>%
    dplyr::mutate(
      age_bin = cut(age_at_index, breaks = age_breaks, include.lowest = TRUE, right = FALSE),
      cem_stratum = as.character(age_bin)
    )

  ## 全体分层目标比例：以所有organ汇总的 age_bin 分布作为共同目标 ####
  target_prop_df <- work_df %>%
    dplyr::count(cem_stratum, name = "n_overall") %>%
    dplyr::mutate(target_prop = n_overall / sum(n_overall))

  organ_levels <- sort(unique(as.character(work_df$primary_organ)))

  matched_list <- lapply(organ_levels, function(current_organ) {
    organ_df <- work_df %>% dplyr::filter(primary_organ == current_organ)
    if (nrow(organ_df) < n_per_organ) {
      ## 样本不足时保留所有样本 ####
      return(organ_df)
    }

    ## 按target_prop对每个stratum分配抽样量（向下取整后补随机） ####
    stratum_count_df <- organ_df %>%
      dplyr::count(cem_stratum, name = "n_available") %>%
      dplyr::left_join(target_prop_df, by = "cem_stratum") %>%
      dplyr::mutate(
        target_prop = dplyr::if_else(is.na(target_prop), 0, target_prop),
        n_target_raw = target_prop * n_per_organ,
        n_target_floor = pmin(floor(n_target_raw), n_available)
      )

    n_left <- n_per_organ - sum(stratum_count_df$n_target_floor)

    if (n_left > 0) {
      stratum_count_df <- stratum_count_df %>%
        dplyr::mutate(
          remainder = n_target_raw - floor(n_target_raw),
          capacity_remaining = n_available - n_target_floor
        )

      ordered_indices <- order(-stratum_count_df$remainder)
      for (stratum_idx in ordered_indices) {
        if (n_left <= 0) break
        if (stratum_count_df$capacity_remaining[stratum_idx] > 0) {
          take <- min(stratum_count_df$capacity_remaining[stratum_idx], n_left)
          stratum_count_df$n_target_floor[stratum_idx] <- stratum_count_df$n_target_floor[stratum_idx] + take
          stratum_count_df$capacity_remaining[stratum_idx] <- stratum_count_df$capacity_remaining[stratum_idx] - take
          n_left <- n_left - take
        }
      }
    }

    ## 若仍有剩余（target_prop在organ内全无），从organ内随机补齐 ####
    picked_df_list <- lapply(seq_len(nrow(stratum_count_df)), function(row_idx) {
      current_stratum <- stratum_count_df$cem_stratum[row_idx]
      n_take <- stratum_count_df$n_target_floor[row_idx]
      if (n_take <= 0) return(NULL)
      stratum_pool <- organ_df %>% dplyr::filter(cem_stratum == current_stratum)
      stratum_pool[sample(nrow(stratum_pool), n_take), , drop = FALSE]
    })
    picked_df <- dplyr::bind_rows(picked_df_list)

    still_left <- n_per_organ - nrow(picked_df)
    if (still_left > 0) {
      remaining_pool <- dplyr::anti_join(
        organ_df, picked_df,
        by = "sample_col_id"
      )
      if (nrow(remaining_pool) > 0) {
        extra_df <- remaining_pool[sample(nrow(remaining_pool), min(still_left, nrow(remaining_pool))), , drop = FALSE]
        picked_df <- dplyr::bind_rows(picked_df, extra_df)
      }
    }
    picked_df
  })

  matched_sample_df <- dplyr::bind_rows(matched_list) %>% as.data.frame()

  ## 匹配前后协变量分布对比：age的均值/SD/中位数（仅年龄） ####
  build_balance_row <- function(covariate_name, before_values, after_values) {
    data.frame(
      covariate = covariate_name,
      n_before = sum(!is.na(before_values)),
      mean_before = mean(before_values, na.rm = TRUE),
      sd_before = stats::sd(before_values, na.rm = TRUE),
      median_before = stats::median(before_values, na.rm = TRUE),
      n_after = sum(!is.na(after_values)),
      mean_after = mean(after_values, na.rm = TRUE),
      sd_after = stats::sd(after_values, na.rm = TRUE),
      median_after = stats::median(after_values, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }

  balance_overall_df <- build_balance_row("age_at_index", work_df$age_at_index, matched_sample_df$age_at_index) %>%
    dplyr::mutate(primary_organ = "ALL", .before = 1)

  balance_by_organ_df <- dplyr::bind_rows(lapply(organ_levels, function(current_organ) {
    before_organ <- work_df %>% dplyr::filter(primary_organ == current_organ)
    after_organ  <- matched_sample_df %>% dplyr::filter(primary_organ == current_organ)
    build_balance_row("age_at_index", before_organ$age_at_index, after_organ$age_at_index) %>%
      dplyr::mutate(primary_organ = current_organ, .before = 1)
  }))

  balance_df <- dplyr::bind_rows(balance_overall_df, balance_by_organ_df)

  list(
    matched_sample_df = matched_sample_df,
    balance_df = balance_df,
    target_prop_df = target_prop_df
  )
}


## 函数：run_cem_matched_mixed_models ####
## 功能：在CEM匹配后的样本子集上对所有feature拟合LMM，并把返回列前缀从mixed_重命名为cem_mixed_ ####
## 输入：matched_sample_df，CEM匹配后的样本表；feature_columns，特征列；feature_info_df，特征信息表 ####
## 输出：列表，包含cem_overall_df、cem_pairwise_df、cem_emmean_df、cem_one_vs_rest_df ####
run_cem_matched_mixed_models <- function(matched_sample_df, feature_columns, feature_info_df) {
  raw_result_list <- run_mixed_effect_models(
    sample_score_df = matched_sample_df,
    feature_columns = feature_columns,
    feature_info_df = feature_info_df
  )

  rename_mixed_cols <- function(df) {
    if (is.null(df) || nrow(df) == 0) return(df)
    old_names <- colnames(df)
    new_names <- gsub("^mixed_", "cem_mixed_", old_names)
    new_names <- ifelse(old_names == "is_singular", "cem_is_singular", new_names)
    colnames(df) <- new_names
    df
  }

  list(
    cem_overall_df = rename_mixed_cols(raw_result_list$overall_df),
    cem_pairwise_df = rename_mixed_cols(raw_result_list$pairwise_df),
    cem_emmean_df = rename_mixed_cols(raw_result_list$emmean_df),
    cem_one_vs_rest_df = rename_mixed_cols(raw_result_list$one_vs_rest_df)
  )
}


# 7. 检查输入与依赖 ####
log_message("========== 7.4 肿瘤免疫微环境分析开始 ==========")

assert_required_files(
  file_paths = c(tcga_tumor_expression_rdata, tcga_tumor_sample_annotation_rdata, tcga_panimmune_tcr_file),
  file_labels = c("TCGA肿瘤表达矩阵(原始)", "TCGA肿瘤样本注释", "Thorsson/TCGA PanImmune TCR summary")
)

if (!dir.exists(tcga_maf_dir)) {
  stop("缺少TCGA MAF目录: ", tcga_maf_dir)
}

dependency_check_df <- check_required_packages(
  packages = required_packages,
  output_path = file.path(table_output_path, "依赖检查.xlsx")
)

log_message("依赖检查通过。")
log_message(paste0("并行设置: enable_parallel=", enable_parallel, "; parallel_workers=", parallel_workers))


# 8. 读取并整理表达矩阵和样本注释 ####
log_message("读取TCGA肿瘤表达矩阵和样本注释...")

tcga_tumor_expression_matrix <- load_single_object_from_rdata(
  rdata_file = tcga_tumor_expression_rdata,
  object_name = "tcga_tumor_expression_matrix"
)

tcga_tumor_sample_annotation_df <- load_single_object_from_rdata(
  rdata_file = tcga_tumor_sample_annotation_rdata,
  object_name = "tcga_tumor_sample_annotation_df"
)

check_required_columns(
  data_df = tcga_tumor_sample_annotation_df,
  required_columns = c("sample_barcode", "patient_id", "project_id", "primary_organ", "pattern_group", "sample_type"),
  data_label = "TCGA肿瘤样本注释"
)

log_message(paste0("表达矩阵维度: ", nrow(tcga_tumor_expression_matrix), " × ", ncol(tcga_tumor_expression_matrix)))
log_message(paste0("样本注释行数: ", nrow(tcga_tumor_sample_annotation_df)))

tumor_meta_df <- tcga_tumor_sample_annotation_df %>%
  dplyr::mutate(
    sample_col_id = as.character(sample_barcode),
    sample_id = as.character(sample_id),
    patient_id = as.character(patient_id),
    project_id = as.character(project_id),
    tcga_code = gsub("^TCGA-", "", project_id),
    sample_type = as.character(sample_type),
    tissue_type = as.character(tissue_type),
    primary_organ = as.character(primary_organ),
    pattern_group = as.character(pattern_group)
  ) %>%
  dplyr::filter(
    sample_type == "Primary Tumor",
    sample_col_id %in% colnames(tcga_tumor_expression_matrix),
    !is.na(primary_organ),
    !is.na(pattern_group)
  ) %>%
  dplyr::arrange(
    factor(pattern_group, levels = pattern_group_order),
    primary_organ,
    project_id,
    patient_id,
    sample_col_id
  )

if (keep_first_tumor_per_patient && "sample_rank_within_patient" %in% colnames(tumor_meta_df)) {
  n_before_patient_filter <- nrow(tumor_meta_df)
  tumor_meta_df <- tumor_meta_df %>%
    dplyr::filter(is.na(sample_rank_within_patient) | sample_rank_within_patient == 1)
  log_message(paste0("按患者保留首个肿瘤样本: ", n_before_patient_filter, " -> ", nrow(tumor_meta_df)))
}

tcga_tumor_expression_matrix <- tcga_tumor_expression_matrix[, tumor_meta_df$sample_col_id, drop = FALSE]
stopifnot(identical(colnames(tcga_tumor_expression_matrix), tumor_meta_df$sample_col_id))

log_message(paste0("最终纳入样本数: ", nrow(tumor_meta_df)))
log_message(paste0("最终纳入TCGA项目数: ", length(unique(tumor_meta_df$project_id))))
log_message(paste0("最终纳入primary_organ数: ", length(unique(tumor_meta_df$primary_organ))))
log_message(paste0("最终纳入pattern_group数: ", length(unique(tumor_meta_df$pattern_group))))

log_message(paste0("log2转换: log2(TPM + ", log2_offset, ")"))
tcga_tumor_log2_matrix <- log2(tcga_tumor_expression_matrix + log2_offset)
storage.mode(tcga_tumor_log2_matrix) <- "numeric"

rm(tcga_tumor_expression_matrix)
gc()


# 9. 构建signature并检查基因可用性 ####
log_message("构建signature基因集...")

available_genes <- rownames(tcga_tumor_log2_matrix)
gene_alias_map <- build_gene_alias_map()
feature_info_df <- build_feature_info_table()

bindea_signature_list <- build_bindea_signature_list()
selected_bindea_signature_names <- c(
  "CD8 T cells",
  "Cytotoxic cells",
  "Treg cells",
  "Macrophages",
  "B cells",
  "NK cells",
  "DC",
  "Neutrophils"
)
bindea_signature_list <- bindea_signature_list[selected_bindea_signature_names]
bindea_signature_list <- lapply(
  bindea_signature_list,
  apply_gene_aliases,
  available_genes = available_genes,
  alias_map = gene_alias_map
)

hallmark_target_sets <- c(
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
)
hallmark_signature_list <- get_hallmark_signature_list(hallmark_target_sets)
hallmark_signature_list <- lapply(
  hallmark_signature_list,
  apply_gene_aliases,
  available_genes = available_genes,
  alias_map = gene_alias_map
)

custom_signature_list <- build_custom_signature_list()
custom_signature_list <- lapply(
  custom_signature_list,
  apply_gene_aliases,
  available_genes = available_genes,
  alias_map = gene_alias_map
)

bindea_availability_df <- format_signature_availability(
  signature_list = bindea_signature_list,
  available_genes = available_genes,
  module_name = "Key Immune Cell Infiltration",
  method_name = "Bindea signatures + ssGSEA",
  min_genes = min_genes_bindea
) %>%
  dplyr::mutate(feature_name = paste0("Bindea_", make_safe_feature_name(signature_name)))

hallmark_availability_df <- format_signature_availability(
  signature_list = hallmark_signature_list,
  available_genes = available_genes,
  module_name = "Hallmark pathway activity",
  method_name = "MSigDB Hallmark + ssGSEA",
  min_genes = min_genes_default
) %>%
  dplyr::mutate(feature_name = paste0("Hallmark_", make_safe_feature_name(signature_name)))

custom_availability_df <- format_signature_availability(
  signature_list = custom_signature_list,
  available_genes = available_genes,
  module_name = "Custom immune signatures",
  method_name = "Mean z-score",
  min_genes = min_genes_default
)

checkpoint_gene_availability_df <- format_signature_availability(
  signature_list = list(
    "PDCD1_expression" = "PDCD1",
    "CD274_expression" = "CD274",
    "CTLA4_expression" = "CTLA4"
  ),
  available_genes = available_genes,
  module_name = "Checkpoint Gene Expression",
  method_name = "log2(TPM + 1)",
  min_genes = 1
)


# 10. 计算样本级免疫特征 ####
log_message("计算TMB和TCR diversity...")
tmb_result <- calculate_tmb_scores(
  tumor_meta_df = tumor_meta_df,
  maf_dir = tcga_maf_dir,
  exome_size_mb = exome_size_mb_for_tmb,
  use_parallel = enable_parallel,
  workers = parallel_workers
)

tcr_result <- calculate_tcr_scores(
  tumor_meta_df = tumor_meta_df,
  tcr_file = tcga_panimmune_tcr_file
)

log_message("计算ESTIMATE相关分数...")
estimate_result <- calculate_estimate_scores(
  expression_matrix = tcga_tumor_log2_matrix,
  alias_map = gene_alias_map,
  use_parallel = enable_parallel,
  workers = parallel_workers
)

log_message("计算关键免疫细胞ssGSEA分数...")
bindea_score_matrix <- run_ssgsea_scores(
  expression_matrix = tcga_tumor_log2_matrix,
  gene_sets = bindea_signature_list,
  min_genes = min_genes_bindea,
  use_parallel = enable_parallel,
  workers = parallel_workers
)
bindea_result <- transpose_score_matrix(
  score_matrix = bindea_score_matrix,
  feature_prefix = "Bindea_",
  label_prefix = "",
  module_name = "Key Immune Cell Infiltration",
  method_name = "Bindea signatures + ssGSEA",
  source_name = "Bindea et al., 2013"
)
rm(bindea_score_matrix)
gc()

log_message("计算MSigDB Hallmark ssGSEA分数...")
hallmark_score_matrix <- run_ssgsea_scores(
  expression_matrix = tcga_tumor_log2_matrix,
  gene_sets = hallmark_signature_list,
  min_genes = min_genes_default,
  use_parallel = enable_parallel,
  workers = parallel_workers
)
hallmark_result <- transpose_score_matrix(
  score_matrix = hallmark_score_matrix,
  feature_prefix = "Hallmark_",
  label_prefix = "",
  module_name = "Hallmark pathway activity",
  method_name = "MSigDB Hallmark + ssGSEA",
  source_name = "Liberzon et al., 2015; msigdbr"
)
rm(hallmark_score_matrix)
gc()

log_message("计算小基因集mean z-score、CYT score与checkpoint基因表达...")
custom_score_df <- data.frame(sample_col_id = colnames(tcga_tumor_log2_matrix), stringsAsFactors = FALSE)
custom_score_df$CYT_score <- calculate_cyt_score(tcga_tumor_log2_matrix)

for (signature_name in names(custom_signature_list)) {
  custom_score_df[[signature_name]] <- calculate_mean_z_score(
    expression_matrix = tcga_tumor_log2_matrix,
    genes = custom_signature_list[[signature_name]],
    min_genes = min_genes_default
  )
}

custom_score_df$PDCD1_expression <- calculate_single_gene_expression(tcga_tumor_log2_matrix, "PDCD1")
custom_score_df$CD274_expression <- calculate_single_gene_expression(tcga_tumor_log2_matrix, "CD274")
custom_score_df$CTLA4_expression <- calculate_single_gene_expression(tcga_tumor_log2_matrix, "CTLA4")

cyt_availability_df <- format_signature_availability(
  signature_list = list("CYT_score" = c("GZMA", "PRF1")),
  available_genes = available_genes,
  module_name = "T Cell Activation and Cytotoxicity",
  method_name = "Mean log2 expression",
  min_genes = 2
)

signature_availability_df <- dplyr::bind_rows(
  estimate_result$availability_df,
  bindea_availability_df,
  hallmark_availability_df,
  cyt_availability_df,
  custom_availability_df,
  checkpoint_gene_availability_df
) %>%
  dplyr::left_join(
    feature_info_df %>% dplyr::select(feature_name, feature_label, feature_axis, source),
    by = "feature_name"
  ) %>%
  dplyr::arrange(module, method, feature_name)

sample_score_df <- tumor_meta_df %>%
  dplyr::select(
    sample_col_id,
    sample_barcode,
    sample_id,
    patient_id,
    project_id,
    tcga_code,
    sample_type,
    tissue_type,
    primary_organ,
    pattern_group,
    dplyr::any_of("sample_rank_within_patient")
  ) %>%
  dplyr::left_join(tmb_result$score_df %>% dplyr::select(sample_col_id, TMB), by = "sample_col_id") %>%
  dplyr::left_join(tcr_result$score_df %>% dplyr::select(sample_col_id, TCR_Shannon), by = "sample_col_id") %>%
  dplyr::left_join(estimate_result$score_df, by = "sample_col_id") %>%
  dplyr::left_join(bindea_result$score_df, by = "sample_col_id") %>%
  dplyr::left_join(hallmark_result$score_df, by = "sample_col_id") %>%
  dplyr::left_join(custom_score_df, by = "sample_col_id")

feature_columns <- feature_info_df$feature_name
missing_feature_columns <- setdiff(feature_columns, colnames(sample_score_df))
if (length(missing_feature_columns) > 0) {
  stop("样本级评分表缺少以下特征列: ", paste(missing_feature_columns, collapse = ", "))
}

log_message(paste0("样本级免疫特征表维度: ", nrow(sample_score_df), " × ", ncol(sample_score_df)))


# 11. 癌种层面汇总 ####
log_message("汇总TCGA project级和primary_organ级免疫特征...")

project_score_df <- summarise_scores_by_group(
  sample_score_df = sample_score_df,
  group_columns = c("pattern_group", "primary_organ", "project_id", "tcga_code"),
  feature_columns = feature_columns
) %>%
  dplyr::arrange(
    factor(pattern_group, levels = pattern_group_order),
    primary_organ,
    project_id
  )

primary_organ_score_df <- summarise_scores_by_group(
  sample_score_df = sample_score_df,
  group_columns = c("pattern_group", "primary_organ"),
  feature_columns = feature_columns
) %>%
  dplyr::arrange(
    factor(pattern_group, levels = pattern_group_order),
    primary_organ
  )

pattern_group_summary_df <- build_pattern_group_summary(
  primary_organ_score_df = primary_organ_score_df,
  feature_columns = feature_columns,
  feature_info_df = feature_info_df
)


# 12. 主统计检验 ####
log_message("执行primary_organ级Kruskal-Wallis和两两Wilcoxon检验...")

kruskal_result_df <- run_kruskal_tests(
  primary_organ_score_df = primary_organ_score_df,
  feature_columns = feature_columns,
  feature_info_df = feature_info_df
)

pairwise_result_df <- run_pairwise_wilcox_tests(
  primary_organ_score_df = primary_organ_score_df,
  feature_columns = feature_columns,
  feature_info_df = feature_info_df
)


# 13. 分层bootstrap稳定性分析 ####
log_message("执行样本-癌种两层分层bootstrap稳定性分析...")

bootstrap_result_list <- run_stratified_tme_bootstrap(
  sample_score_df = sample_score_df,
  feature_columns = feature_columns,
  n_iterations = bootstrap_n_iterations,
  samples_per_primary_organ = bootstrap_samples_per_primary_organ,
  min_samples_per_primary_organ = bootstrap_min_samples_per_primary_organ,
  primary_organs_per_pattern_group = bootstrap_primary_organs_per_pattern_group,
  seed = bootstrap_seed,
  use_parallel = enable_parallel,
  workers = parallel_workers
)

bootstrap_group_summary_df <- summarise_bootstrap_group_scores(
  bootstrap_group_median_array = bootstrap_result_list$bootstrap_group_median_array,
  feature_info_df = feature_info_df
)

bootstrap_feature_stability_df <- summarise_bootstrap_feature_stability(
  bootstrap_group_median_array = bootstrap_result_list$bootstrap_group_median_array,
  feature_info_df = feature_info_df
) %>%
  dplyr::left_join(
    kruskal_result_df %>%
      dplyr::select(
        feature_name,
        observed_highest_pattern_group = highest_pattern_group,
        observed_lowest_pattern_group = lowest_pattern_group,
        observed_median_range = median_range,
        observed_p_value = p_value,
        observed_fdr = fdr
      ),
    by = "feature_name"
  )

bootstrap_pairwise_direction_df <- summarise_bootstrap_pairwise_direction(
  bootstrap_group_median_array = bootstrap_result_list$bootstrap_group_median_array,
  feature_info_df = feature_info_df
)

bootstrap_parameter_df <- bootstrap_result_list$parameter_df


# 14. 混合效应模型敏感性分析 ####
log_message("执行样本级混合效应模型敏感性分析...")

mixed_model_result_list <- run_mixed_effect_models(
  sample_score_df = sample_score_df,
  feature_columns = feature_columns,
  feature_info_df = feature_info_df
)

mixed_model_overall_df <- mixed_model_result_list$overall_df
mixed_model_pairwise_df <- mixed_model_result_list$pairwise_df
mixed_model_emmean_df <- mixed_model_result_list$emmean_df
mixed_model_one_vs_rest_df <- mixed_model_result_list$one_vs_rest_df


# 14b. 等癌种权重混合效应模型补充分析 ####
log_message("执行每癌种等权重(weights = 1 / n_organ)的混合效应模型补充分析...")

weighted_mixed_model_result_list <- run_weighted_mixed_effect_models(
  sample_score_df = sample_score_df,
  feature_columns = feature_columns,
  feature_info_df = feature_info_df
)

weighted_mixed_model_overall_df <- weighted_mixed_model_result_list$overall_df
weighted_mixed_model_pairwise_df <- weighted_mixed_model_result_list$pairwise_df
weighted_mixed_model_emmean_df <- weighted_mixed_model_result_list$emmean_df
weighted_mixed_model_one_vs_rest_df <- weighted_mixed_model_result_list$one_vs_rest_df


# 14c. 分层抽样混合效应模型补充分析 ####
log_message(paste0(
  "执行分层抽样混合模型补充分析: 每个primary_organ抽 ",
  subsample_lmm_n_per_organ,
  " 样本，重复 ",
  subsample_lmm_n_iterations,
  " 次..."
))

subsample_mixed_model_iteration_list <- run_stratified_subsample_mixed_models(
  sample_score_df = sample_score_df,
  feature_columns = feature_columns,
  n_iterations = subsample_lmm_n_iterations,
  n_per_organ = subsample_lmm_n_per_organ,
  seed = subsample_lmm_seed,
  use_parallel = enable_parallel,
  workers = parallel_workers
)

subsample_mixed_model_iteration_overall_df <- subsample_mixed_model_iteration_list$iteration_overall_df
subsample_mixed_model_iteration_pairwise_df <- subsample_mixed_model_iteration_list$iteration_pairwise_df
subsample_mixed_model_iteration_blup_df <- subsample_mixed_model_iteration_list$iteration_blup_df
subsample_mixed_model_iteration_emmean_df <- subsample_mixed_model_iteration_list$iteration_emmean_df
subsample_mixed_model_iteration_one_vs_rest_df <- subsample_mixed_model_iteration_list$iteration_one_vs_rest_df

subsample_mixed_model_overall_df <- summarise_subsample_overall(
  iteration_overall_df = subsample_mixed_model_iteration_overall_df,
  feature_info_df = feature_info_df
)

subsample_mixed_model_pairwise_df <- summarise_subsample_pairwise(
  iteration_pairwise_df = subsample_mixed_model_iteration_pairwise_df,
  feature_info_df = feature_info_df
)

subsample_mixed_model_emmean_df <- summarise_subsample_emmean(
  iteration_emmean_df = subsample_mixed_model_iteration_emmean_df,
  feature_info_df = feature_info_df
)

subsample_mixed_model_one_vs_rest_df <- summarise_subsample_one_vs_rest(
  iteration_one_vs_rest_df = subsample_mixed_model_iteration_one_vs_rest_df,
  feature_info_df = feature_info_df
)


# 14d. LMM癌种随机截距(BLUP)提取 ####
log_message("提取主LMM中primary_organ的随机截距BLUP及条件SE...")

lmm_organ_random_intercept_df <- extract_lmm_organ_random_intercepts(
  sample_score_df = sample_score_df,
  feature_columns = feature_columns,
  feature_info_df = feature_info_df,
  use_equal_organ_weights = FALSE
)

log_message("提取加权LMM中primary_organ的随机截距BLUP及条件SE...")

weighted_lmm_organ_random_intercept_df <- extract_lmm_organ_random_intercepts(
  sample_score_df = sample_score_df,
  feature_columns = feature_columns,
  feature_info_df = feature_info_df,
  use_equal_organ_weights = TRUE
)

log_message("聚合分层抽样LMM每次迭代的primary_organ随机截距BLUP...")

subsample_lmm_organ_random_intercept_df <- summarise_subsample_blup(
  iteration_blup_df = subsample_mixed_model_iteration_blup_df,
  feature_info_df = feature_info_df
)


# 14e. CEM匹配LMM补充分析：基于年龄的粗糙精确匹配 ####
log_message(paste0(
  "执行CEM匹配LMM补充分析：每个primary_organ按 age_bin 匹配抽样 ",
  cem_lmm_n_per_organ, " 样本..."
))

log_message("读取TCGA各项目临床文件提取age_at_index...")
project_ids_used <- unique(sample_score_df$project_id)
patient_age_df <- load_tcga_patient_age(
  project_ids = project_ids_used,
  clinical_dir = clinical_data_dir
)
log_message(paste0("共读取到 ", nrow(patient_age_df), " 位患者的age_at_index"))

log_message("将age_at_index合并到样本级特征表...")
sample_score_with_age_df <- attach_patient_age_to_sample_df(
  sample_score_df = sample_score_df,
  patient_age_df = patient_age_df
)
log_message(paste0(
  "合并后样本总数: ", nrow(sample_score_with_age_df),
  "; age缺失填充样本数: ", sum(sample_score_with_age_df$age_imputed_flag)
))

log_message("构建CEM匹配后的样本子集...")
cem_match_result_list <- build_cem_matched_sample_set(
  sample_score_df = sample_score_with_age_df,
  n_per_organ = cem_lmm_n_per_organ,
  age_breaks = cem_age_bin_breaks,
  seed = cem_lmm_seed
)
cem_matched_sample_df <- cem_match_result_list$matched_sample_df
cem_covariate_balance_df <- cem_match_result_list$balance_df
cem_stratum_target_df <- cem_match_result_list$target_prop_df
log_message(paste0(
  "CEM匹配后样本总数: ", nrow(cem_matched_sample_df),
  "; 覆盖organ数: ", dplyr::n_distinct(cem_matched_sample_df$primary_organ)
))

log_message("在CEM匹配子集上拟合LMM...")
cem_mixed_model_result_list <- run_cem_matched_mixed_models(
  matched_sample_df = cem_matched_sample_df,
  feature_columns = feature_columns,
  feature_info_df = feature_info_df
)
cem_mixed_model_overall_df <- cem_mixed_model_result_list$cem_overall_df
cem_mixed_model_pairwise_df <- cem_mixed_model_result_list$cem_pairwise_df
cem_mixed_model_emmean_df <- cem_mixed_model_result_list$cem_emmean_df
cem_mixed_model_one_vs_rest_df <- cem_mixed_model_result_list$cem_one_vs_rest_df


# 15. 方向和effect size主解释汇总 ####
log_message("整理方向稳定性和effect size主解释汇总...")

direction_effect_size_summary_df <- build_direction_effect_size_summary(
  pattern_group_summary_df = pattern_group_summary_df,
  primary_organ_score_df = primary_organ_score_df,
  feature_info_df = feature_info_df,
  bootstrap_feature_stability_df = bootstrap_feature_stability_df,
  mixed_model_overall_df = mixed_model_overall_df
)

## 加权LMM版本方向效应量汇总：复用bootstrap方向，仅替换LMM来源 ####
weighted_overall_standardized_df <- weighted_mixed_model_overall_df %>%
  dplyr::transmute(
    feature_name = feature_name,
    model_status = weighted_mixed_status,
    model_message = weighted_mixed_message,
    mixed_highest_group = weighted_mixed_highest_group,
    mixed_lowest_group = weighted_mixed_lowest_group,
    mixed_highest_emmean = weighted_mixed_highest_emmean,
    mixed_lowest_emmean = weighted_mixed_lowest_emmean,
    mixed_standardized_range = weighted_mixed_standardized_range,
    mixed_pattern_group_f_value = weighted_mixed_pattern_group_f_value,
    mixed_pattern_group_p_value = weighted_mixed_pattern_group_p_value,
    mixed_pattern_group_fdr = weighted_mixed_pattern_group_fdr,
    is_singular = weighted_mixed_is_singular
  )

weighted_direction_effect_size_summary_df <- build_direction_effect_size_summary(
  pattern_group_summary_df = pattern_group_summary_df,
  primary_organ_score_df = primary_organ_score_df,
  feature_info_df = feature_info_df,
  bootstrap_feature_stability_df = bootstrap_feature_stability_df,
  mixed_model_overall_df = weighted_overall_standardized_df
)

## 分层抽样LMM版本方向效应量汇总：用聚合的median标准化指标 ####
subsample_overall_standardized_df <- subsample_mixed_model_overall_df %>%
  dplyr::transmute(
    feature_name = feature_name,
    model_status = dplyr::case_when(
      is.na(subsample_median_p_value) ~ "insufficient_data",
      TRUE ~ "success"
    ),
    model_message = NA_character_,
    mixed_highest_group = subsample_most_often_highest_group,
    mixed_lowest_group = subsample_most_often_lowest_group,
    mixed_highest_emmean = NA_real_,
    mixed_lowest_emmean = NA_real_,
    mixed_standardized_range = subsample_median_standardized_range,
    mixed_pattern_group_f_value = NA_real_,
    mixed_pattern_group_p_value = subsample_median_p_value,
    mixed_pattern_group_fdr = NA_real_,
    is_singular = NA
  )

subsample_direction_effect_size_summary_df <- build_direction_effect_size_summary(
  pattern_group_summary_df = pattern_group_summary_df,
  primary_organ_score_df = primary_organ_score_df,
  feature_info_df = feature_info_df,
  bootstrap_feature_stability_df = bootstrap_feature_stability_df,
  mixed_model_overall_df = subsample_overall_standardized_df
)


# 16. 参数汇总 ####
parameter_summary_df <- tibble::tibble(
  parameter = c(
    "log2_transform",
    "keep_first_tumor_per_patient",
    "main_statistical_unit",
    "supplemental_unit",
    "bootstrap_algorithm",
    "bootstrap_n_iterations",
    "bootstrap_samples_per_primary_organ",
    "bootstrap_primary_organs_per_pattern_group",
    "enable_parallel",
    "parallel_workers",
    "tmb_exome_size_mb",
    "tcr_input_file",
    "n_sample_level_samples",
    "n_tcga_projects",
    "n_primary_organs",
    "n_pattern_groups",
    "n_tme_features",
    "estimate_purity_note",
    "direction_stability_rule",
    "effect_size_rule",
    "mixed_model_role",
    "mixed_model_formula",
    "mixed_model_transform",
    "weighted_mixed_model_role",
    "weighted_mixed_model_weights",
    "subsample_mixed_model_role",
    "subsample_mixed_model_n_iterations",
    "subsample_mixed_model_n_per_organ",
    "subsample_mixed_model_seed",
    "cem_mixed_model_role",
    "cem_mixed_model_covariates",
    "cem_mixed_model_age_breaks",
    "cem_mixed_model_n_per_organ",
    "cem_mixed_model_seed",
    "cem_mixed_model_n_matched_samples",
    "cem_mixed_model_n_organs"
  ),
  value = c(
    paste0("log2(TPM + ", log2_offset, ")"),
    as.character(keep_first_tumor_per_patient),
    "primary_organ",
    "TCGA project",
    "Two-stage balanced hierarchical bootstrap",
    as.character(bootstrap_n_iterations),
    as.character(bootstrap_samples_per_primary_organ),
    bootstrap_parameter_df$value[bootstrap_parameter_df$parameter == "primary_organs_per_pattern_group"],
    as.character(enable_parallel),
    as.character(parallel_workers),
    as.character(exome_size_mb_for_tmb),
    tcga_panimmune_tcr_file,
    as.character(nrow(sample_score_df)),
    as.character(nrow(project_score_df)),
    as.character(nrow(primary_organ_score_df)),
    as.character(length(unique(primary_organ_score_df$pattern_group))),
    as.character(length(feature_columns)),
    "TumorPurity_ESTIMATE is formula-derived and used for relative comparison in TCGA TPM RNA-seq data.",
    "Stable: highest_group_probability >= 0.80; Trend: 0.60-0.80; Unstable: < 0.60.",
    "Large: standardized_median_range >= 0.80; Moderate: 0.50-0.80; Small: < 0.50.",
    "Sensitivity analysis only; p values and FDR do not replace direction stability and effect size interpretation.",
    "feature_z ~ pattern_group + (1 | primary_organ)",
    "TMB uses log1p(TMB) then z-score; all other features use sample-level z-score.",
    "Equal-cancer-type-weight LMM: each primary_organ contributes equal total weight; addresses dominance of large organs in pattern_group estimation.",
    "weights = 1 / n_samples_in_primary_organ (every primary_organ gets total weight = 1)",
    "Stratified-subsampling LMM: subsample to equal N per primary_organ then refit LMM; repeat to obtain robust estimates and 95% CI.",
    as.character(subsample_lmm_n_iterations),
    as.character(subsample_lmm_n_per_organ),
    as.character(subsample_lmm_seed),
    "CEM-matched LMM: per-organ coarsened exact matching on age strata to balance age distribution across primary_organs; single LMM fit on matched subset.",
    "age_at_index",
    paste(as.character(cem_age_bin_breaks), collapse = ","),
    as.character(cem_lmm_n_per_organ),
    as.character(cem_lmm_seed),
    as.character(nrow(cem_matched_sample_df)),
    as.character(dplyr::n_distinct(cem_matched_sample_df$primary_organ))
  )
)


# 17. 保存结果 ####
log_message("保存7.4分析结果...")

openxlsx::write.xlsx(
  x = list(
    "样本级免疫特征" = sample_score_df,
    "TCGA项目级补充特征" = project_score_df,
    "原发部位级免疫特征" = primary_organ_score_df,
    "模式群特征汇总" = pattern_group_summary_df,
    "方向效应量汇总" = direction_effect_size_summary_df,
    "方向效应量汇总_加权" = weighted_direction_effect_size_summary_df,
    "方向效应量汇总_分层抽样" = subsample_direction_effect_size_summary_df,
    "模式群差异检验" = kruskal_result_df,
    "两两比较结果" = pairwise_result_df,
    "混合模型整体检验" = mixed_model_overall_df,
    "混合模型两两比较" = mixed_model_pairwise_df,
    "加权混合模型整体" = weighted_mixed_model_overall_df,
    "加权混合模型两两" = weighted_mixed_model_pairwise_df,
    "分层抽样LMM整体" = subsample_mixed_model_overall_df,
    "分层抽样LMM两两" = subsample_mixed_model_pairwise_df,
    "LMM癌种随机截距" = lmm_organ_random_intercept_df,
    "LMM癌种随机截距_加权" = weighted_lmm_organ_random_intercept_df,
    "LMM癌种随机截距_分层抽样" = subsample_lmm_organ_random_intercept_df,
    "LMM癌种群emmean" = mixed_model_emmean_df,
    "LMM癌种群emmean_加权" = weighted_mixed_model_emmean_df,
    "LMM癌种群emmean_分层抽样" = subsample_mixed_model_emmean_df,
    "LMM癌种群one_vs_rest" = mixed_model_one_vs_rest_df,
    "LMM癌种群one_vs_rest_加权" = weighted_mixed_model_one_vs_rest_df,
    "LMM癌种群one_vs_rest_分层抽样" = subsample_mixed_model_one_vs_rest_df,
    "CEM匹配LMM整体" = cem_mixed_model_overall_df,
    "CEM匹配LMM两两" = cem_mixed_model_pairwise_df,
    "CEM匹配LMM癌种群emmean" = cem_mixed_model_emmean_df,
    "CEM匹配LMM癌种群one_vs_rest" = cem_mixed_model_one_vs_rest_df,
    "CEM匹配协变量平衡" = cem_covariate_balance_df,
    "CEM匹配分层目标比例" = cem_stratum_target_df,
    "抽样模式群评分" = bootstrap_group_summary_df,
    "抽样特征稳定性" = bootstrap_feature_stability_df,
    "抽样两两方向" = bootstrap_pairwise_direction_df,
    "抽样纳入癌种" = bootstrap_result_list$organ_info_df,
    "抽样参数" = bootstrap_parameter_df,
    "TMB可用性" = tmb_result$availability_df,
    "TCR可用性" = tcr_result$availability_df,
    "Signature基因可用性" = signature_availability_df,
    "特征信息表" = feature_info_df,
    "依赖检查" = dependency_check_df,
    "参数汇总" = parameter_summary_df
  ),
  file = output_xlsx,
  overwrite = TRUE
)

save(sample_score_df, file = sample_score_rdata)
save(project_score_df, file = project_score_rdata)
save(primary_organ_score_df, file = primary_organ_score_rdata)
save(
  bootstrap_group_summary_df,
  bootstrap_feature_stability_df,
  bootstrap_pairwise_direction_df,
  bootstrap_parameter_df,
  direction_effect_size_summary_df,
  weighted_direction_effect_size_summary_df,
  subsample_direction_effect_size_summary_df,
  mixed_model_overall_df,
  mixed_model_pairwise_df,
  weighted_mixed_model_overall_df,
  weighted_mixed_model_pairwise_df,
  subsample_mixed_model_overall_df,
  subsample_mixed_model_pairwise_df,
  subsample_mixed_model_iteration_overall_df,
  subsample_mixed_model_iteration_pairwise_df,
  subsample_mixed_model_iteration_blup_df,
  subsample_mixed_model_iteration_emmean_df,
  subsample_mixed_model_iteration_one_vs_rest_df,
  lmm_organ_random_intercept_df,
  weighted_lmm_organ_random_intercept_df,
  subsample_lmm_organ_random_intercept_df,
  mixed_model_emmean_df,
  weighted_mixed_model_emmean_df,
  subsample_mixed_model_emmean_df,
  mixed_model_one_vs_rest_df,
  weighted_mixed_model_one_vs_rest_df,
  subsample_mixed_model_one_vs_rest_df,
  cem_matched_sample_df,
  cem_covariate_balance_df,
  cem_stratum_target_df,
  cem_mixed_model_overall_df,
  cem_mixed_model_pairwise_df,
  cem_mixed_model_emmean_df,
  cem_mixed_model_one_vs_rest_df,
  file = bootstrap_result_rdata
)

log_message(paste0("已导出Excel结果: ", output_xlsx))
log_message(paste0("已保存样本级Rdata: ", sample_score_rdata))
log_message(paste0("已保存TCGA项目级Rdata: ", project_score_rdata))
log_message(paste0("已保存原发部位级Rdata: ", primary_organ_score_rdata))
log_message(paste0("已保存bootstrap稳定性Rdata: ", bootstrap_result_rdata))
log_message("========== 7.4 肿瘤免疫微环境分析全部完成 ==========")
