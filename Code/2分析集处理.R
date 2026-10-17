# 两个数据库分析集处理（FAERS + JADER）####
# 本脚本将2.1、2.2合并为一个脚本，只保留数据处理，不包含可视化
rm(list = ls())
gc()
# 加载通用配置
source("Code/config.R")
source("Code/mapping_rules.R")
## 1.1 加载FAERS数据####

message("加载FAERS数据...")
load(file.path(.output_path, "faers_demo.RData"))
load(file.path(.output_path, "faers_drug.RData"))
load(file.path(.output_path, "faers_reac.RData"))
load(file.path(.output_path, "faers_indi.RData"))

message(paste("FAERS Demo:", nrow(faers_demo), "条"))
message(paste("FAERS Drug:", nrow(faers_drug), "条"))
message(paste("FAERS Reac:", nrow(faers_reac), "条"))
message(paste("FAERS Indi:", nrow(faers_indi), "条"))

## 1.2 加载共用筛选规则（从mapping_rules.R）####
# 筛选规则现已移至 Code/mapping_rules.R 统一管理
## 1.3 癌种映射函数（已移至mapping_rules.R）####

## 1.4 从适应证表筛选癌症患者####

message("\n筛选癌症患者...")

# 使用mapping_rules.R中定义的排除PT列表（已在脚本开头加载）
# benign_tumor_pts: 良性/交界性肿瘤PT（误匹配BLASTOMA等关键词）
# non_malignant_pts: 非肿瘤疾病PT（误匹配METASTATIC关键词）
exclude_pts <- toupper(c(benign_tumor_pts, non_malignant_pts))
message(paste("  待排除的非恶性PT数量:", length(exclude_pts)))

# 查看indi表结构
head(faers_indi)
colnames(faers_indi)

# 按primaryid + indi_pt_std去重
message(paste("去重前适应证记录数:", nrow(faers_indi)))
faers_indi_dedup <- faers_indi %>%
  distinct(primaryid, indi_pt_std, .keep_all = TRUE)
message(paste("去重后适应证记录数:", nrow(faers_indi_dedup)))

# 设置并行计算（使用所有可用核心）
# Linux系统使用multicore（fork模式）更高效，避免socket传输大数据
n_cores <- parallel::detectCores() - 1
# options(future.globals.maxSize = 1024 * 1024 * 1024)  # 1GB
plan(multicore, workers = n_cores)
message(paste("使用", n_cores, "个核心进行并行计算"))

# 并行进行PT级标记 + 癌种映射
message("开始并行PT级标记和癌种映射...")
start_time <- Sys.time()

# 分别并行计算每个字段（避免内存爆炸）
message("  计算hit_malignant...")
hit_malignant_vec <- future_sapply(faers_indi_dedup$indi_pt_std, function(pt) {
  grepl(malignant_confirmed_pattern, pt, perl = TRUE)
}, USE.NAMES = FALSE)

message("  计算hit_benign...")
hit_benign_vec <- future_sapply(faers_indi_dedup$indi_pt_std, function(pt) {
  grepl(benign_pattern, pt, perl = TRUE)
}, USE.NAMES = FALSE)

message("  计算hit_uncertain...")
hit_uncertain_vec <- future_sapply(faers_indi_dedup$indi_pt_std, function(pt) {
  grepl(uncertain_pattern, pt, perl = TRUE)
}, USE.NAMES = FALSE)

message("  计算hit_metastasis...")
hit_metastasis_vec <- future_sapply(faers_indi_dedup$indi_pt_std, function(pt) {
  grepl(metastasis_pattern, pt, perl = TRUE)
}, USE.NAMES = FALSE)

message("  计算hit_generalized...")
hit_generalized_vec <- future_sapply(faers_indi_dedup$indi_pt_std, function(pt) {
  grepl(generalized_only_pattern, pt, perl = TRUE)
}, USE.NAMES = FALSE)

message("  计算hit_false_positive（假阳性排除）...")
hit_false_positive_vec <- future_sapply(faers_indi_dedup$indi_pt_std, function(pt) {
  grepl(false_positive_pattern, pt, perl = TRUE)
}, USE.NAMES = FALSE)

message("  计算hit_met_to_site（转移到某部位）...")
hit_met_to_site_vec <- future_sapply(faers_indi_dedup$indi_pt_std, function(pt) {
  grepl(met_to_site_pattern, pt, perl = TRUE)
}, USE.NAMES = FALSE)

message("  计算cancer_type...")
cancer_type_vec <- future_sapply(faers_indi_dedup$indi_pt_std, map_cancer_type, USE.NAMES = FALSE)

message("  计算is_excluded（排除非恶性PT）...")
is_excluded_vec <- future_sapply(faers_indi_dedup$indi_pt_std, function(pt) {
  toupper(pt) %in% exclude_pts
}, USE.NAMES = FALSE)

message("  计算indi_pt_normalized（PT标准化）...")
indi_pt_normalized_vec <- future_sapply(faers_indi_dedup$indi_pt_std, function(pt) {
  # 标准化PT名称：
  # 1. 保护HIGH/LOW GRADE
  # 2. 删除STAGE/GRADE及其后面内容
  # 3. 删除状态修饰词（METASTATIC、RECURRENT等），确保同一癌种的不同状态被识别为同一器官
  # 4. 恢复保护的内容
  pt_temp <- gsub("(HIGH|LOW)([ -])GRADE", "\\1\\2XPROTECTEDX", pt, ignore.case = TRUE)
  pt_temp <- trimws(sub("\\s*(STAGE|GRADE).*", "", pt_temp, ignore.case = TRUE))
  # 删除状态修饰词（确保"COLORECTAL CANCER METASTATIC"和"COLORECTAL CANCER"被识别为同一器官）
  pt_temp <- trimws(sub("\\s+(METASTATIC|RECURRENT|REFRACTORY|IN\\s+REMISSION|LOCALISED|LOCALIZED|REGIONAL|EXTENSIVE|LIMITED)$", "", pt_temp, ignore.case = TRUE))
  pt_temp <- gsub("XPROTECTEDX", "GRADE", pt_temp)
  return(pt_temp)
}, USE.NAMES = FALSE)

# 计算最终的hit_malignant（向量化操作，非常快）
# 排除非恶性PT（如Lipoblastoma、Metastatic Cutaneous Crohn's Disease等）
hit_malignant_final <- ifelse(is_excluded_vec, FALSE, hit_malignant_vec)
# 排除假阳性PT（如Metastatic Pulmonary Embolism等非癌症诊断）
hit_malignant_final <- ifelse(hit_false_positive_vec, FALSE, hit_malignant_final)
# 如果有癌种映射，则标记为恶性
hit_malignant_final <- ifelse(!is.na(cancer_type_vec), TRUE, hit_malignant_final)

# 预计算器官映射（用于未映射PT的合并判定）
# 关键修正：同一器官的多个PT应合并为单一原发癌，避免被归入other_malignancy
message("  计算primary_organ_for_grouping（用于合并同器官PT）...")
primary_organ_vec <- future_mapply(
  function(pt_norm, cancer) {
    map_to_primary_organ(pt_norm, cancer)
  },
  indi_pt_normalized_vec,
  cancer_type_vec,
  SIMPLIFY = TRUE,
  USE.NAMES = FALSE
)

message("  计算type_fine和type_coarse...")
# 并行计算type_fine和type_coarse
# 关键修正：type_fine使用Histologic_category（严格按照Results/癌种映射.xlsx）
# 同一Histologic_category下的PT会合并为同一type_fine
type_results <- future_mapply(
  function(cancer, hit_mal, hit_mts, hit_gen, pt_norm, organ) {
    # type_fine：用于癌种判定
    # 严格按照用户的三层映射体系：
    # 1. 通过TCGA代码或PT获取Histologic_category
    # 2. 无法映射但有转移关键词 → INDETERMINATE
    # 3. 无法映射但确认恶性 → 使用器官默认的Histologic_category
    type_fine <- if (hit_mal && hit_mts) {
      "INDETERMINATE"
    } else if (hit_mal && hit_gen) {
      "INDETERMINATE"
    } else if (hit_mal && !hit_gen && !hit_mts) {
      # 使用map_to_histologic_category获取Histologic_category
      histcat <- map_to_histologic_category(pt_norm, cancer, organ)
      if (!is.na(histcat)) {
        histcat
      } else if (!is.na(organ) && organ != "Unknown" && organ != "CUP") {
        paste0("UNMAPPED:", organ)
      } else {
        paste0("UNMAPPED_PRIMARY:", pt_norm)
      }
    } else {
      NA_character_
    }

    # type_coarse：最终癌种展示（使用Primary_organs_tissues）
    type_coarse <- if (!is.na(cancer)) {
      cancer
    } else if (hit_mal && hit_mts) {
      "INDETERMINATE"
    } else if (hit_mal && hit_gen) {
      "INDETERMINATE"
    } else if (hit_mal && !hit_gen && !hit_mts) {
      "Others"
    } else {
      NA_character_
    }

    c(type_fine, type_coarse)
  },
  cancer_type_vec,
  hit_malignant_final,
  hit_met_to_site_vec,
  hit_generalized_vec,
  indi_pt_normalized_vec,
  primary_organ_vec,
  SIMPLIFY = TRUE,
  USE.NAMES = FALSE
)

# 提取type_fine和type_coarse向量
type_fine_vec <- type_results[1, ]
type_coarse_vec <- type_results[2, ]

# 组装最终结果（简单赋值，非常快）
faers_indi_tagged <- faers_indi_dedup %>%
  mutate(
    hit_malignant = hit_malignant_final,
    hit_metastasis = hit_metastasis_vec,
    hit_met_to_site = hit_met_to_site_vec,
    hit_benign = hit_benign_vec,
    hit_uncertain = hit_uncertain_vec,
    hit_generalized = hit_generalized_vec,
    cancer_type = cancer_type_vec,
    indi_pt_normalized = indi_pt_normalized_vec,
    type_fine = type_fine_vec,
    type_coarse = type_coarse_vec,
    # 添加PT级器官映射（用于判断同器官多癌种情况）
    primary_organ_pt = primary_organ_vec
  )

# 清理临时向量
rm(
  hit_malignant_vec, hit_benign_vec, hit_uncertain_vec, hit_metastasis_vec,
  hit_generalized_vec, hit_false_positive_vec, hit_met_to_site_vec, cancer_type_vec,
  is_excluded_vec, indi_pt_normalized_vec, hit_malignant_final,
  type_results, type_fine_vec, type_coarse_vec, primary_organ_vec
)

end_time <- Sys.time()
message(paste("PT级标记和癌种映射完成，耗时:", round(difftime(end_time, start_time, units = "mins"), 2), "分钟"))

# 恢复顺序计划
plan(sequential)
gc()

# 统计PT级各分类
message(paste("确认恶性记录数:", sum(faers_indi_tagged$hit_malignant)))
message(paste("  其中转移相关:", sum(faers_indi_tagged$hit_metastasis)))
message(paste("明确良性记录数:", sum(faers_indi_tagged$hit_benign)))
message(paste("良恶性不明记录数:", sum(faers_indi_tagged$hit_uncertain)))
message(paste("可映射癌种记录数:", sum(!is.na(faers_indi_tagged$cancer_type))))

## 1.5 PRIMARYID级互斥归类（三大类：恶性/良性/良恶性不明）####

message("\n开始PRIMARYID级互斥归类...")

# 使用data.table进行高效并行分组聚合（data.table内置多线程优化）
library(data.table)
n_cores <- parallel::detectCores() - 1
setDTthreads(n_cores)
message(paste("使用", n_cores, "个data.table线程进行并行分组聚合"))

start_time <- Sys.time()

# 转换为data.table
indi_dt <- as.data.table(faers_indi_tagged)

# 分组聚合（data.table自动多线程）
primaryid_summary_dt <- indi_dt[, .(
  # 基础信号
  has_malignant = any(hit_malignant),
  has_benign = any(hit_benign),
  has_uncertain = any(hit_uncertain),
  # 纯不明信号（排除恶性重叠，用于统计）
  has_uncertain_only = any(hit_uncertain & !hit_malignant),
  # metastasis_only临时标记：所有适应证都包含转移关键词
  metastasis_only_flag = all(hit_metastasis),
  # 统计type_fine数量（排除INDETERMINATE，用于判断单癌种/多癌种）
  n_types_fine = uniqueN(type_fine[!is.na(type_fine) & type_fine != "INDETERMINATE"]),
  n_types_tcga = uniqueN(cancer_type[!is.na(cancer_type)]),
  # 是否只有INDETERMINATE（无明确原发癌种）
  has_only_indeterminate = {
    type_fine_nonNA <- type_fine[!is.na(type_fine)]
    length(type_fine_nonNA) > 0 && all(type_fine_nonNA == "INDETERMINATE")
  },
  # 获取唯一type_coarse（排除INDETERMINATE，临时存储）
  assigned_type_temp = {
    valid_types <- type_coarse[!is.na(type_coarse) & type_coarse != "INDETERMINATE"]
    if (length(valid_types) > 0) unique(valid_types)[1] else NA_character_
  },
  # 获取唯一type_fine（排除INDETERMINATE，临时存储）
  assigned_type_fine_temp = {
    valid_types <- type_fine[!is.na(type_fine) & type_fine != "INDETERMINATE"]
    if (length(valid_types) > 0) unique(valid_types)[1] else NA_character_
  },
  # 获取癌种列表（用于描述）
  types_list = paste(unique(type_fine[!is.na(type_fine)]), collapse = ";"),
  # 新增：计算恶性肿瘤PT的唯一器官数量（用于同器官多癌种判断）
  # 例如：GLIOBLASTOMA和BRAIN NEOPLASM MALIGNANT都属于Brain，则unique_organs=1
  n_unique_organs = {
    malignant_organs <- primary_organ_pt[hit_malignant & !is.na(primary_organ_pt)]
    uniqueN(malignant_organs)
  },
  # 新增：获取同器官情况下的器官名称
  same_organ_temp = {
    malignant_organs <- primary_organ_pt[hit_malignant & !is.na(primary_organ_pt)]
    unique_organs <- unique(malignant_organs)
    if (length(unique_organs) == 1) unique_organs[1] else NA_character_
  }
), by = primaryid]

# 计算派生字段（metastasis_only需要结合has_malignant）
primaryid_summary_dt[, `:=`(
  metastasis_only = has_malignant & metastasis_only_flag,
  assigned_type = fifelse(n_types_fine == 1, assigned_type_temp, NA_character_),
  assigned_type_fine = fifelse(n_types_fine == 1, assigned_type_fine_temp, NA_character_)
)]

# 丢弃判定（互斥，更清晰）
primaryid_summary_dt[, `:=`(
  benign_only = has_benign & !has_malignant & !has_uncertain,
  indeterminate_only = has_uncertain & !has_malignant
)]

# 最终分区（5类互斥，按优先级）使用data.table的fcase
primaryid_summary_dt[, cohort := fcase(
  # A. 非肿瘤组（排除）
  !has_malignant & !benign_only & !indeterminate_only, "non_tumor",
  # B. 排除：确定良性（纯良性，无恶性/不明信号）
  benign_only, "exclude_benign",
  # C. 描述性分析：良恶性不明（有borderline信号但无恶性）
  indeterminate_only, "descriptive_indeterminate",
  # D. 主分析：单一恶性肿瘤（n_types_fine=1，包括TCGA映射和Others）
  has_malignant & n_types_fine == 1, "single_malignancy",
  # E. 其他恶性肿瘤（多癌种共存或原发不可判定）
  has_malignant & n_types_fine != 1, "other_malignancy",
  # 其他（理论上不应出现）
  default = "other"
)]

# 删除临时列（保留same_organ_temp用于后续other_malignancy的器官分配）
primaryid_summary_dt[, c("metastasis_only_flag", "assigned_type_temp", "assigned_type_fine_temp") := NULL]

# 转换回tibble（保持与后续代码兼容）
primaryid_summary <- as_tibble(primaryid_summary_dt)

# 清理临时对象
rm(indi_dt, primaryid_summary_dt)

end_time <- Sys.time()
message(paste("PRIMARYID级归类完成，耗时:", round(difftime(end_time, start_time, units = "secs"), 2), "秒"))

# 保存中间数据集，方便后续调用（避免重复运行耗时的处理步骤）
message("\n保存中间数据集...")
save(primaryid_summary, file = file.path(.output_path, "各主键信号标记.RData"))
save(faers_indi_tagged, file = file.path(.output_path, "faers_indi_tagged.RData"))
message(paste("  primaryid_summary已保存:", file.path(.output_path, "各主键信号标记.RData")))
message(paste("  faers_indi_tagged已保存:", file.path(.output_path, "faers_indi_tagged.RData")))

# 统计各分区
message("\n========== PRIMARYID级分区统计 ==========")
message("\n分区分布:")
print(table(primaryid_summary$cohort))

## 1.55 第一层器官映射（进一步精细化处理）####
message("\n========== 第一层器官映射 ==========")

# 加载中间数据集（如果从此处开始运行）
if (!exists("faers_indi_tagged")) {
  message("加载中间数据集 faers_indi_tagged...")
  load(file.path(.output_path, "faers_indi_tagged.RData"))
}
if (!exists("primaryid_summary")) {
  message("加载中间数据集 primaryid_summary...")
  load(file.path(.output_path, "各主键信号标记.RData"))
}

# 确保映射规则已加载（包含select_primary_pt函数）
# 这对于分段运行或并行计算时函数可用性至关重要
if (!exists("select_primary_pt")) {
  message("加载映射规则（Code/mapping_rules.R）...")
  source("Code/mapping_rules.R")
}

# 从faers_indi_tagged获取肿瘤相关PT，构建first_pt
# 获取每个primaryid的肿瘤相关PT（只保留匹配malignant_confirmed_pattern的PT）
pt_tumor_by_primaryid_temp <- faers_indi_tagged %>%
  filter(hit_malignant) %>%
  group_by(primaryid) %>%
  summarise(indi_pt_list_tumor = paste(unique(indi_pt_std), collapse = " & "), .groups = "drop")

# 对single_malignancy的病例进行第一层器官映射
single_malignancy_for_mapping <- primaryid_summary %>%
  filter(cohort == "single_malignancy") %>%
  left_join(pt_tumor_by_primaryid_temp, by = "primaryid")

# 智能选择原发癌种PT（排除转移部位PT，如 "Metastases To Bone"）
# 使用 select_primary_pt 函数，已包含标准化处理（去除STAGE/GRADE，首字母大写）
single_malignancy_for_mapping <- single_malignancy_for_mapping %>%
  mutate(first_pt = sapply(indi_pt_list_tumor, select_primary_pt, USE.NAMES = FALSE))

a <- single_malignancy_for_mapping %>%
  filter(assigned_type == "Others") %>%
  count(first_pt, sort = TRUE)
write_delim(a, file.path(.results_path, "Others_first_pt.csv"), delim = ",")
# 应用第一层器官映射
message("应用第一层器官映射...")
message("  映射优先级: 1. assigned_type (TCGA代码) -> 2. first_pt (PT名称)")

single_malignancy_for_mapping <- single_malignancy_for_mapping %>%
  mutate(primary_organ = mapply(
    map_to_primary_organ,
    first_pt = first_pt,
    assigned_type = assigned_type,
    USE.NAMES = FALSE
  ))

# 统计映射结果
organ_dist <- single_malignancy_for_mapping %>%
  filter(!is.na(primary_organ)) %>%
  count(primary_organ, sort = TRUE) %>%
  mutate(pct = round(n / sum(n) * 100, 2))
message("\n已映射的器官分布:")
print(organ_dist)

# 未映射的first_pt统计
unmapped <- single_malignancy_for_mapping %>%
  filter(is.na(primary_organ)) %>%
  count(first_pt, sort = TRUE)
message(paste("\n未映射的first_pt数量:", nrow(unmapped)))
if (nrow(unmapped) > 0) {
  message("未映射的first_pt列表（按频次排序，前20个）:")
  print(head(unmapped, 20))
  # 将完整的未映射PT列表保存到文件，方便后续补充mapping_rules.R
  if (nrow(unmapped) > 20) {
    message(paste("完整未映射PT列表已保存至:", file.path(results_path, "FAERS_unmapped_pts.csv")))
    write.csv(unmapped, file.path(results_path, "FAERS_unmapped_pts.csv"), row.names = FALSE)
  }
}

# 将first_pt和primary_organ合并回primaryid_summary
primaryid_summary <- primaryid_summary %>%
  left_join(
    single_malignancy_for_mapping %>% select(primaryid, first_pt, primary_organ),
    by = "primaryid"
  ) %>%
  # 确保所有恶性肿瘤病例都有器官映射
  mutate(
    primary_organ = case_when(
      # other_malignancy：先检查是否所有PT属于同一器官
      # 例如：GLIOBLASTOMA和BRAIN NEOPLASM MALIGNANT都属于Brain
      # 此时primary_organ归类为Brain而不是Unknown_Primary，cohort仍为other_malignancy
      cohort == "other_malignancy" & !is.na(same_organ_temp) ~ same_organ_temp,
      # other_malignancy但多个器官或无法确认原发灶：标记为Unknown_Primary
      cohort == "other_malignancy" & is.na(same_organ_temp) ~ "Unknown_Primary",
      # single_malignancy但未能映射到具体器官：标记为Unmapped
      cohort == "single_malignancy" & is.na(primary_organ) ~ "Unmapped",
      # 已有映射的保持不变
      TRUE ~ primary_organ
    )
  )

# 清理临时对象
rm(pt_tumor_by_primaryid_temp, organ_dist, unmapped)
gc()

## 1.56 第二层组织学分类映射####
message("\n========== 第二层组织学分类映射 ==========")

# 应用第二层组织学分类映射
message("应用第二层组织学分类映射...")
message("  映射优先级: 1. first_pt (PT名称) -> 2. assigned_type (TCGA代码)")

single_malignancy_for_mapping <- single_malignancy_for_mapping %>%
  mutate(histologic_category = mapply(
    map_to_histologic_category,
    first_pt = first_pt,
    assigned_type = assigned_type,
    primary_organ = primary_organ,
    USE.NAMES = FALSE
  ))

# 统计第二层映射结果
histology_dist <- single_malignancy_for_mapping %>%
  filter(!is.na(histologic_category)) %>%
  count(histologic_category, sort = TRUE) %>%
  mutate(pct = round(n / sum(n) * 100, 2))
message("\n已映射的组织学分类分布:")
print(histology_dist)

# 未映射的first_pt统计（第二层）
unmapped_histology <- single_malignancy_for_mapping %>%
  filter(is.na(histologic_category)) %>%
  count(first_pt, sort = TRUE)
message(paste("\n未映射组织学分类的first_pt数量:", nrow(unmapped_histology)))
if (nrow(unmapped_histology) > 0) {
  message("未映射组织学分类的first_pt列表（按频次排序，前20个）:")
  print(head(unmapped_histology, 20))
}

# 将histologic_category合并回primaryid_summary
# 首先检查single_malignancy_for_mapping中是否有histologic_category列
if ("histologic_category" %in% names(single_malignancy_for_mapping)) {
  primaryid_summary <- primaryid_summary %>%
    left_join(
      single_malignancy_for_mapping %>% select(primaryid, histologic_category),
      by = "primaryid"
    )
} else {
  warning("single_malignancy_for_mapping中没有histologic_category列，可能需要先运行第1.56节的组织学分类映射代码")
  primaryid_summary <- primaryid_summary %>%
    mutate(histologic_category = NA_character_)
}

# 确保histologic_category列存在（即使left_join没有添加）
if (!"histologic_category" %in% names(primaryid_summary)) {
  primaryid_summary <- primaryid_summary %>%
    mutate(histologic_category = NA_character_)
}

# 对于非single_malignancy的病例，histologic_category设为Unknown
primaryid_summary <- primaryid_summary %>%
  mutate(
    histologic_category = case_when(
      cohort == "other_malignancy" ~ "Unknown",
      cohort == "single_malignancy" & is.na(histologic_category) ~ "Unknown",
      TRUE ~ histologic_category
    )
  )

# 清理临时对象
rm(single_malignancy_for_mapping, histology_dist, unmapped_histology)
gc()

## 1.6 保存最终映射表####

# 加载中间数据集（如果从此处开始运行）
if (!exists("faers_indi_tagged")) {
  message("加载中间数据集 faers_indi_tagged...")
  load(file.path(.output_path, "faers_indi_tagged.RData"))
}
if (!exists("primaryid_summary")) {
  message("加载中间数据集 primaryid_summary...")
  load(file.path(.output_path, "各主键信号标记.RData"))
}

# 获取每个primaryid的适应症PT（去重后用&合并，收集所有PT）
pt_by_primaryid <- faers_indi_tagged %>%
  group_by(primaryid) %>%
  summarise(indi_pt_list = paste(unique(indi_pt_std), collapse = " & "), .groups = "drop")

# 获取每个primaryid的肿瘤相关PT（只保留匹配malignant_confirmed_pattern的PT）
pt_tumor_by_primaryid <- faers_indi_tagged %>%
  filter(hit_malignant) %>%
  group_by(primaryid) %>%
  summarise(indi_pt_list_tumor = paste(unique(indi_pt_std), collapse = " & "), .groups = "drop")

# 保存主键+队列+癌种标注+适应症PT映射表（包含三层分类：primary_organ, histologic_category, assigned_type）
cohort_mapping <- primaryid_summary %>%
  select(primaryid, cohort, assigned_type, first_pt, primary_organ, histologic_category) %>%
  left_join(pt_by_primaryid, by = "primaryid") %>%
  left_join(pt_tumor_by_primaryid, by = "primaryid")
cohort_mapping <- cohort_mapping %>%
  mutate(primary_organ = if_else(cohort == "other_malignancy", "CUP", primary_organ))
save(cohort_mapping, file = file.path(.output_path, "cohort_mapping.RData"))
message(paste("\n映射表已保存:", file.path(.output_path, "cohort_mapping.RData")))
message(paste("共", nrow(cohort_mapping), "条记录"))

# 统计各分区
message("\n各分区统计:")
print(table(cohort_mapping$cohort))

# 释放内存
gc()


# 清理环境，准备第二部分
rm(list = ls())
gc()

# 第二部分：筛选ICI用药患者####
# 路径已在config.R中设置
# 加载必要数据（可独立运行）
if (!exists("faers_drug")) {
  load(file.path(.output_path, "faers_drug.RData"))
}
if (!exists("cohort_mapping")) {
  load(file.path(.output_path, "cohort_mapping.RData"))
}

# 只分析恶性肿瘤患者（排除non_tumor、exclude_benign、descriptive_indeterminate）
malignant_primaryids <- cohort_mapping %>%
  filter(!cohort %in% c("non_tumor", "exclude_benign", "descriptive_indeterminate")) %>%
  pull(primaryid)
message(paste("恶性肿瘤患者数:", length(malignant_primaryids)))

## 2.1 定义ICI药物列表####
# PD-1抑制剂
drug_PD_1 <- list(
  c("NIVOLUMAB", "OPDUALAG", "OPDIVO"),
  c("PEMBROLIZUMAB", "KEYTRUDA"),
  c("CEMIPLIMAB", "LIBTAYO"),
  c("DOSTARLIMAB", "JEMPERLI"),
  c("RETIFANLIMAB", "ZYNYZ"),
  c("TORIPALIMAB", "LOQTORZI"),
  c("TISLELIZUMAB", "TEVIMBRA")
)

# PD-L1抑制剂
drug_PD_L1 <- list(
  c("ATEZOLIZUMAB", "TECENTRIQ", "TECENTRIQ HYBREZA"),
  c("AVELUMAB", "BAVENCIO"),
  c("DURVALUMAB", "IMFINZI")
)

# CTLA-4抑制剂
drug_CTLA_4 <- list(
  c("IPILIMUMAB", "YERVOY"),
  c("TREMELIMUMAB", "IMJUDO")
)

# 合并所有ICI药物
drug_ICI_all <- c(drug_PD_1, drug_PD_L1, drug_CTLA_4)

## 2.2 ICI药物筛选函数####
# 匹配药物名称（drugname_std或prod_ai列）
match_ICI <- function(drugname_std, prod_ai, drug_list) {
  drugname_std <- toupper(drugname_std)
  prod_ai <- toupper(prod_ai)

  for (i in seq_along(drug_list)) {
    pattern <- paste(drug_list[[i]], collapse = "|")
    if (grepl(pattern, drugname_std, ignore.case = TRUE) ||
      grepl(pattern, prod_ai, ignore.case = TRUE)) {
      return(drug_list[[i]][1]) # 返回通用名
    }
  }
  return(NA_character_)
}

## 2.3 筛选ICI用药患者####
message("\n筛选ICI用药患者...")

# 只保留恶性肿瘤患者的用药记录
faers_drug_malignant <- faers_drug %>%
  filter(primaryid %in% malignant_primaryids)
message(paste("恶性肿瘤患者用药记录数:", nrow(faers_drug_malignant)))

# 设置并行计算（使用multicore避免socket传输问题）
library(future.apply)
n_cores <- parallel::detectCores() - 1
plan(multicore, workers = n_cores)
message(paste("使用", n_cores, "个核心进行并行计算"))

# 并行标记ICI药物
message("并行匹配ICI药物...")
start_time <- Sys.time()

ici_drug_vec <- future_mapply(
  function(d, p) match_ICI(d, p, drug_ICI_all),
  faers_drug_malignant$drugname_std,
  faers_drug_malignant$prod_ai,
  USE.NAMES = FALSE
)

end_time <- Sys.time()
message(paste("ICI匹配完成，耗时:", round(difftime(end_time, start_time, units = "mins"), 2), "分钟"))

# 恢复顺序计划
plan(sequential)

# 筛选ICI用药记录
faers_drug_ici <- faers_drug_malignant %>%
  mutate(ici_drug = ici_drug_vec) %>%
  filter(!is.na(ici_drug))

message(paste("ICI用药记录数:", nrow(faers_drug_ici)))

# 按primaryid去重（同一患者使用同一ICI药物只计一次）
faers_drug_ici_dedup <- faers_drug_ici %>%
  distinct(primaryid, ici_drug, .keep_all = TRUE)

message(paste("去重后ICI用药记录数:", nrow(faers_drug_ici_dedup)))

## 2.4 ICI类型分类####
# 定义ICI类型映射
ici_type_map <- data.frame(
  ici_drug = c(sapply(drug_PD_1, `[`, 1), sapply(drug_PD_L1, `[`, 1), sapply(drug_CTLA_4, `[`, 1)),
  ici_type = c(
    rep("anti-PD-1", length(drug_PD_1)),
    rep("anti-PD-L1", length(drug_PD_L1)),
    rep("anti-CTLA-4", length(drug_CTLA_4))
  )
)

faers_drug_ici_typed <- faers_drug_ici_dedup %>%
  left_join(ici_type_map, by = "ici_drug")

# 统计各ICI药物使用情况
message("\nICI药物使用统计:")
print(table(faers_drug_ici_typed$ici_drug))

message("\nICI类型统计:")
print(table(faers_drug_ici_typed$ici_type))

## 2.5 按primaryid汇总ICI用药####
# 首先计算每个患者的ICI TTO
# 规则：优先选择有有效TTO的记录，在有效TTO中再按怀疑优先级（PS > SS > 其他）排序
ici_tto_priority <- faers_drug_ici_typed %>%
  mutate(
    # 怀疑优先级：PS=1, SS=2, 其他=3
    role_priority = case_when(
      role_cod_std == "PS" ~ 1L,
      role_cod_std == "SS" ~ 2L,
      TRUE ~ 3L
    ),
    # TTO是否有效（非Unknown）
    tto_valid = (tto != "Unknown" & !is.na(tto)),
    # TTO数值（Unknown转为NA用于排序）
    tto_num = suppressWarnings(as.numeric(tto))
  ) %>%
  # 按primaryid分组，先按tto_valid降序（有效TTO优先），再按role_priority升序
  group_by(primaryid) %>%
  arrange(desc(tto_valid), role_priority, tto_num, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  select(primaryid, ici_tto = tto)

ici_by_primaryid <- faers_drug_ici_typed %>%
  group_by(primaryid) %>%
  summarise(
    ici_drugs = paste(sort(unique(ici_drug)), collapse = "+"),
    ici_types = paste(sort(unique(ici_type)), collapse = "+"),
    n_ici_drugs = n_distinct(ici_drug),
    n_ici_types = n_distinct(ici_type),
    .groups = "drop"
  ) %>%
  left_join(ici_tto_priority, by = "primaryid") %>%
  mutate(
    # 治疗方案分类
    treatment_strategy = case_when(
      n_ici_types > 1 ~ "Combination Therapy",
      TRUE ~ ici_types
    )
  )

# 统计ICI TTO情况
ici_tto_stats <- ici_by_primaryid %>%
  summarise(
    total = n(),
    known_tto = sum(ici_tto != "Unknown", na.rm = TRUE),
    unknown_tto = sum(ici_tto == "Unknown" | is.na(ici_tto))
  )
message(paste(
  "\nICI TTO统计 - 已知:", ici_tto_stats$known_tto,
  "未知:", ici_tto_stats$unknown_tto
))

message(paste("\nICI用药患者数:", nrow(ici_by_primaryid)))

message("\n治疗方案统计:")
print(table(ici_by_primaryid$treatment_strategy))

## 2.6 保存ICI用药映射表####
save(ici_by_primaryid, file = file.path(.output_path, "ici_mapping.RData"))
message(paste("\nICI映射表已保存:", file.path(.output_path, "ici_mapping.RData")))

# 第三部分：识别放疗、化疗和靶向治疗病例####

# 路径已在config.R中设置，如需独立运行请先执行: source("Code/config.R")

# 加载必要数据
if (!exists("faers_drug")) {
  load(file.path(.output_path, "faers_drug.RData"))
}
if (!exists("cohort_mapping")) {
  load(file.path(.output_path, "cohort_mapping.RData"))
}

# 只分析恶性肿瘤患者
malignant_primaryids <- cohort_mapping %>%
  filter(!cohort %in% c("non_tumor", "exclude_benign", "descriptive_indeterminate")) %>%
  pull(primaryid)
message(paste("恶性肿瘤患者数:", length(malignant_primaryids)))

# 筛选恶性肿瘤患者的用药记录
faers_drug_malignant <- faers_drug %>%
  filter(primaryid %in% malignant_primaryids)
message(paste("恶性肿瘤患者用药记录数:", nrow(faers_drug_malignant)))

## 3.1 定义放疗关键词####
# FAERS中放疗关键词分两组：强特异词 + 技术缩写词
# 强特异词（直接命中）
radiotherapy_strong <- c(
  "RADIATION( THERAPY| TREATMENT)?", "RADIOTHERAPY", "IRRADIATION",
  "BRACHYTHERAPY", "RADIOSURGERY",
  "EXTERNAL BEAM", "PROTON( THERAPY| BEAM)?",
  "GAMMA KNIFE", "CYBER ?KNIFE", "TOMOTHERAPY",
  "LINEAR ACCELERATOR", "\\bLINAC\\b",
  "STEREOTACTIC( RADIOTHERAPY| RADIOSURGERY)?",
  "INTRAOPERATIVE RADIOTHERAPY", "\\bIORT\\b",
  "WHOLE BRAIN", "CRANIAL IRRADIATION", "PELVIC IRRADIATION",
  "HEAVY ION", "CARBON ION", "PARTICLE THERAPY",
  "X-RAY THERAPY", "PHOTON THERAPY", "ELECTRON BEAM"
)
# 技术缩写词（加词边界减少误检）
radiotherapy_tech <- c("\\bIMRT\\b", "\\bIGRT\\b", "\\bVMAT\\b", "\\bSRS\\b", "\\bSBRT\\b")
# 合并为正则表达式
radiotherapy_pattern <- paste0(
  "(?i)(", paste(radiotherapy_strong, collapse = "|"), ")",
  "|(?i)(", paste(radiotherapy_tech, collapse = "|"), ")"
)

## 3.2 定义细胞毒性化疗药物####
chemo_drugs <- c(
  "DAUNORUBICIN", "CYTARABINE", "ARSENIC", "ASPARAGINASE", "AZACITIDINE",
  "BENDAMUSTINE", "BUSULFAN", "CAPECITABINE", "CARBOPLATIN", "CARMUSTINE",
  "CHLORAMBUCIL", "CISPLATIN", "CALDRIBINE", "CLOFARABINE", "CYCLOPHOSPHAMIDE",
  "DACARBAZINE", "DECITABINE", "DEGARELIX", "DOCETAXEL", "ERIBULIN",
  "ETOPOSIDE", "FLOXIRIDINE", "FLUDARABINE", "FLUOROURACIL", "GEMCITABINE",
  "GOSERELIN", "HISTRELIN", "HYDROXYUREA", "IFOSFAMIDE", "INTERFERON",
  "IRINOTECAN", "IXABEPILONE", "LANREOTIDE", "LEUCOVORIN", "LOMUSTINE",
  "MEDROXYPROGESTERONE", "MEGESTROL", "MELPHALAN", "MERCAPTOPURINE", "MESNA",
  "METHOTREXATE", "MITOXANTRONE", "NELARABINE", "OCTREOTIDE", "OMACETAXINE",
  "OXALIPLATIN", "PACLITAXEL", "PEGASPARGASE", "PEMETREXED", "PENTOSTATIN",
  "PROCARBAZINE", "SAMARIUM", "LEXIDRONAM", "TEMOZOLOMIDE", "TENIPOSIDE",
  "THIOTEPA", "TOPOTECAN", "TRABECTEDIN", "VINBLASTINE", "VINORELBINE",
  "DOXORUBICIN", "IDARUBICIN", "EPIRUBICIN", "BLEOMYCIN", "CABAZITAXEL",
  "DACTINOMYCIN", "MITOMYCIN", "STREPTOZOCIN", "5-FU", "FLUOROURACIL"
)
chemo_pattern <- paste(chemo_drugs, collapse = "|")

## 3.3 定义靶向治疗药物####
targeted_drugs <- c(
  "TRASTUZUMAB", "AFLIBERCEPT", "ALEMTUZUMAB", "ANASTROZOLE", "BELIMUMAB",
  "BELINOSTAT", "BEVACIZUMAB", "BLINATUMOMAB", "BORTEZOMIB", "BRENTUXIMAB",
  "CARFILZOMIB", "CETUXIMAB", "COPANLISIB", "DENOSUMAB", "ELOTUZUMAB",
  "EVEROLIMUS", "EXEMESTANE", "FULVESTRANT", "GEFITINIB", "GEMTUZUMAB",
  "IMATINIB", "INOTUZUMAB", "MOGAMULIZUMAB", "MOXETUMOMAB", "PASUDOTOX",
  "NATALIZUMAB", "NECITUMUMAB", "OBINUTUZUMAB", "OFATUMUMAB", "OLARATUMAB",
  "PANITUMUMAB", "PERTUZUMAB", "PRALATREXATE", "RADIUM", "RAMUCIRUMAB",
  "RITUXIMAB", "ROMIDEPSIN", "SILTUXIMAB", "TAGRAXOFUSP", "TALAZOPARIB",
  "TAMOXIFEN", "TEMSIROLIMUS", "VALRUBICIN", "ZIV-AFLIBERCEPT", "BICALUTAMIDE",
  "FLOXURIDINE", "FLUTAMIDE", "KETOCONAZOLE", "LANREOTIDE", "LENALIDOMIDE",
  "LEVOLEUCOVORIN", "POMALIDOMIDE", "PORFIMER", "PROPYLTHIOURACIL",
  "THALIDOMIDE", "TRIFLURIDINE", "VINCRISTINE", "ABEMACICLIB", "ABIRATERONE",
  "ACALABRUTINIB", "AFATINIB", "ALECTINIB", "ALPELISIB", "APALUTAMIDE",
  "AXITINIB", "BINIMETINIB", "BOSUTINIB", "BRIGATINIB", "CABOZANTINIB",
  "CERITINIB", "COBIMETINIB", "CRIZOTINIB", "DABRAFENIB", "DACOMITINIB",
  "DARATUMUMAB", "DAROLUTAMIDE", "DASATINIB", "DINUTUXIMAB", "DUVELISIB",
  "ELOTUZUMAB", "ENASIDENIB", "ENCORAFENIB", "ENFORTUMAB", "ENTRECTINIB",
  "ENZALUTAMIDE", "ERDAFITINIB", "ERLOTINIB", "FAM-TRASTUZUMAB", "DERUXTECA",
  "FEDRATINIB", "FULVESTRANT", "GILTERITINIB", "GLASDEGIB", "IBRUTINIB",
  "IDELALISIB", "ISATUXIMAB", "ISOTRETINOIN", "IVOSIDENIB", "IXAZOMIB",
  "LAPATINIB", "LAROTRECTINIB", "LENVATINIB", "LETROZOLE", "LORLATINIB",
  "NERATINIB", "NILOTINIB", "NIRAPARIB", "OBINUTUZUMAB", "OLAPARIB",
  "OLARATUMAB", "OSIMERTINIB", "PALBOCICLIB", "PANOBINOSTAT", "PAZOPANIB",
  "PEXIDARTINIB", "POLATUZUMAB", "PONATINIB", "REGORAFENIB", "RIBOCICLIB",
  "RUCAPARIB", "RUXOLITINIB", "MIDOSTAURIN", "SELINEXOR", "SONIDEGIB",
  "SORAFENIB", "SUNITINIB", "TALAZOPARIB", "TOSYKATE", "TAZEMETOSTAT",
  "TOREMIFENE", "TRAMETINIB", "TRETINOIN", "VANDETANIB", "VEMURAFENIB",
  "VENETOCLAX", "VISMODEGIB", "VORINOSTAT", "ZANUBRUTINIB"
)
targeted_pattern <- paste(targeted_drugs, collapse = "|")

## 3.4 标记治疗类型（并行计算）####
message("\n标记治疗类型...")

# 设置并行计算（使用multicore避免socket传输问题）
library(future.apply)
n_cores <- parallel::detectCores() - 1
plan(multicore, workers = n_cores)
message(paste("使用", n_cores, "个核心进行并行计算"))

start_time <- Sys.time()

# 并行标记放疗（检测drugname_std和prod_ai）
message("并行标记放疗...")
is_radiotherapy_vec <- future_mapply(
  function(d, p) {
    grepl(radiotherapy_pattern, d, perl = TRUE) |
      grepl(radiotherapy_pattern, p, perl = TRUE)
  },
  faers_drug_malignant$drugname_std,
  faers_drug_malignant$prod_ai,
  USE.NAMES = FALSE
)

# 并行标记化疗（检测drugname_std和prod_ai）
message("并行标记化疗...")
is_chemotherapy_vec <- future_mapply(
  function(d, p) {
    grepl(chemo_pattern, d, ignore.case = TRUE) |
      grepl(chemo_pattern, p, ignore.case = TRUE)
  },
  faers_drug_malignant$drugname_std,
  faers_drug_malignant$prod_ai,
  USE.NAMES = FALSE
)

# 并行标记靶向治疗（检测drugname_std和prod_ai）
message("并行标记靶向治疗...")
is_targeted_vec <- future_mapply(
  function(d, p) {
    grepl(targeted_pattern, d, ignore.case = TRUE) |
      grepl(targeted_pattern, p, ignore.case = TRUE)
  },
  faers_drug_malignant$drugname_std,
  faers_drug_malignant$prod_ai,
  USE.NAMES = FALSE
)

end_time <- Sys.time()
message(paste("标记完成，耗时:", round(difftime(end_time, start_time, units = "mins"), 2), "分钟"))

# 恢复顺序计划
plan(sequential)

# 合并标记结果
faers_drug_treatment <- faers_drug_malignant %>%
  mutate(
    is_radiotherapy = is_radiotherapy_vec,
    is_chemotherapy = is_chemotherapy_vec,
    is_targeted = is_targeted_vec
  )

# 按primaryid汇总治疗标注
treatment_by_primaryid <- faers_drug_treatment %>%
  group_by(primaryid) %>%
  summarise(
    has_radiotherapy = any(is_radiotherapy),
    has_chemotherapy = any(is_chemotherapy),
    has_targeted = any(is_targeted),
    .groups = "drop"
  ) %>%
  mutate(
    radiotherapy = ifelse(has_radiotherapy, "Yes", "No"),
    chemotherapy = ifelse(has_chemotherapy, "Yes", "No"),
    targeted_therapy = ifelse(has_targeted, "Yes", "No"),
    Targeted_therapy = targeted_therapy
  )

# 统计结果
message("\n放疗使用情况:")
print(table(treatment_by_primaryid$radiotherapy))

message("\n化疗使用情况:")
print(table(treatment_by_primaryid$chemotherapy))

message("\n靶向治疗使用情况:")
print(table(treatment_by_primaryid$targeted_therapy))

## 3.5 保存治疗标注结果####
save(treatment_by_primaryid, file = file.path(.output_path, "treatment_mapping.RData"))
message(paste("\n治疗标注表已保存:", file.path(.output_path, "treatment_mapping.RData")))

# 第四部分：标记irAE不良反应####

# 路径已在config.R中设置，如需独立运行请先执行: source("Code/config.R")

# 加载必要数据
load(file.path(.output_path, "faers_reac.RData"))

if (!exists("cohort_mapping")) {
  load(file.path(.output_path, "cohort_mapping.RData"))
}

## 4.1 加载irAE列表####
library(readxl)
irae_list <- read_excel(file.path(.data_path, "irAE列表.xlsx"))
message(paste("irAE PT数量:", nrow(irae_list)))
message("\nirAE SOC分布:")
print(table(irae_list$SOC))

# 提取irAE PT列表（标准化为大写）
irae_pts <- toupper(trimws(irae_list$PT))

## 4.2 筛选恶性肿瘤患者####
malignant_primaryids <- cohort_mapping %>%
  filter(!cohort %in% c("non_tumor", "exclude_benign", "descriptive_indeterminate")) %>%
  pull(primaryid)
message(paste("\n恶性肿瘤患者数:", length(malignant_primaryids)))

# 筛选恶性肿瘤患者的不良反应记录
faers_reac_malignant <- faers_reac %>%
  filter(primaryid %in% malignant_primaryids)
message(paste("恶性肿瘤患者不良反应记录数:", nrow(faers_reac_malignant)))

## 4.3 标记irAE（并行计算）####
message("\n标记irAE...")

library(future.apply)
n_cores <- parallel::detectCores() - 1
plan(multicore, workers = n_cores)
message(paste("使用", n_cores, "个核心进行并行计算"))

start_time <- Sys.time()

# 创建irAE匹配模式（精确匹配）
irae_pattern <- paste0("^(", paste(irae_pts, collapse = "|"), ")$")

# 并行标记irAE
is_irae_vec <- future_sapply(
  faers_reac_malignant$pt_std,
  function(x) x %in% irae_pts,
  USE.NAMES = FALSE
)

end_time <- Sys.time()
message(paste("标记完成，耗时:", round(difftime(end_time, start_time, units = "secs"), 2), "秒"))

plan(sequential)

# 合并标记结果
faers_reac_irae <- faers_reac_malignant %>%
  mutate(is_irae = is_irae_vec)

## 4.4 按primaryid汇总irAE标注（宽格式：每个PT和SOC为二分类）####
# 创建PT到SOC的映射
pt_to_soc <- setNames(irae_list$SOC, toupper(trimws(irae_list$PT)))

# 为每条irAE记录添加SOC
faers_reac_irae <- faers_reac_irae %>%
  mutate(irae_soc = ifelse(is_irae, pt_to_soc[pt_std], NA_character_))

# 获取所有唯一的PT和SOC
unique_pts <- sort(unique(irae_pts))
unique_socs <- sort(unique(irae_list$SOC))

message(paste("irAE PT数量:", length(unique_pts)))
message(paste("irAE SOC数量:", length(unique_socs)))

# 筛选命中irAE的记录
irae_hits <- faers_reac_irae %>%
  filter(is_irae) %>%
  select(primaryid, pt_std, irae_soc)

# 创建PT宽格式（每个PT一列，0/1）
# 列名中空格替换为下划线，确保代码兼容性
message("\n创建PT宽格式...")
irae_pt_wide <- irae_hits %>%
  distinct(primaryid, pt_std) %>%
  mutate(pt_clean = gsub(" ", "_", pt_std)) %>%
  select(primaryid, pt_clean) %>%
  distinct() %>%
  mutate(value = 1L) %>%
  pivot_wider(
    names_from = pt_clean,
    values_from = value,
    values_fill = 0L,
    names_prefix = "PT_"
  )

# 创建SOC宽格式（每个SOC一列，0/1）
message("创建SOC宽格式...")
irae_soc_wide <- irae_hits %>%
  distinct(primaryid, irae_soc) %>%
  mutate(soc_clean = gsub(" ", "_", irae_soc)) %>%
  select(primaryid, soc_clean) %>%
  distinct() %>%
  mutate(value = 1L) %>%
  pivot_wider(
    names_from = soc_clean,
    values_from = value,
    values_fill = 0L,
    names_prefix = "SOC_"
  )

# 获取所有恶性肿瘤患者primaryid（包括无irAE的）
all_malignant_ids <- data.frame(primaryid = malignant_primaryids)

# 合并PT和SOC宽格式，补全无irAE患者为0
irae_by_primaryid <- all_malignant_ids %>%
  left_join(irae_pt_wide, by = "primaryid") %>%
  left_join(irae_soc_wide, by = "primaryid")

# 将NA替换为0
irae_by_primaryid[is.na(irae_by_primaryid)] <- 0L

# 添加汇总列
irae_by_primaryid <- irae_by_primaryid %>%
  mutate(
    has_irae = rowSums(select(., starts_with("PT_"))) > 0,
    irae_pt_count = rowSums(select(., starts_with("PT_"))),
    irae_soc_count = rowSums(select(., starts_with("SOC_")))
  )

# 统计结果
message("\nirAE发生情况:")
print(table(irae_by_primaryid$has_irae))

message("\nirAE PT数量分布（有irAE的患者）:")
print(summary(irae_by_primaryid$irae_pt_count[irae_by_primaryid$has_irae]))

message("\n各SOC发生率（前10）:")
soc_cols <- names(irae_by_primaryid)[grepl("^SOC_", names(irae_by_primaryid))]
soc_rates <- colSums(irae_by_primaryid[soc_cols]) / nrow(irae_by_primaryid) * 100
print(head(sort(soc_rates, decreasing = TRUE), 10))

## 4.5 保存irAE标注结果####
save(irae_by_primaryid, file = file.path(.output_path, "irae_mapping.RData"))
message(paste("\nirAE标注表已保存:", file.path(.output_path, "irae_mapping.RData")))

# 第五部分：合并所有标记并保存恶性肿瘤样本####
library(tidyverse)

# 路径已在config.R中设置，如需独立运行请先执行: source("Code/config.R")

## 5.1 加载所有数据####
message("\n加载demo表格...")
load(file.path(.output_path, "faers_demo.RData"))
message(paste("demo记录数:", nrow(faers_demo)))

message("\n加载各标记表...")
# 队列映射
if (!exists("cohort_mapping")) {
  load(file.path(.output_path, "cohort_mapping.RData"))
}
# ICI映射
if (!exists("ici_by_primaryid")) {
  load(file.path(.output_path, "ici_mapping.RData"))
}
# 治疗类型映射
if (!exists("treatment_by_primaryid")) {
  load(file.path(.output_path, "treatment_mapping.RData"))
}
# irAE映射
if (!exists("irae_by_primaryid")) {
  load(file.path(.output_path, "irae_mapping.RData"))
}

## 5.2 筛选恶性肿瘤患者demo####
malignant_cohorts <- c("other_malignancy", "single_malignancy")
malignant_primaryids <- cohort_mapping %>%
  filter(cohort %in% malignant_cohorts) %>%
  pull(primaryid)

faers_demo_malignant <- faers_demo %>%
  filter(primaryid %in% malignant_primaryids)
message(paste("\n恶性肿瘤患者demo记录数:", nrow(faers_demo_malignant)))

## 5.3 合并队列信息####
message("\n合并队列信息...")
faers_malignant_final <- faers_demo_malignant %>%
  left_join(
    cohort_mapping,
    by = "primaryid"
  )

## 5.4 合并ICI信息####
message("合并ICI信息...")
# 简化ici_by_primaryid，只保留关键列（包含ici_tto）
ici_simple <- ici_by_primaryid %>%
  select(primaryid, ici_drugs, ici_types, treatment_strategy, ici_tto)

faers_malignant_final <- faers_malignant_final %>%
  left_join(ici_simple, by = "primaryid") %>%
  mutate(
    has_ici = !is.na(ici_drugs),
    ici_drugs = ifelse(is.na(ici_drugs), "None", ici_drugs),
    ici_types = ifelse(is.na(ici_types), "None", ici_types),
    treatment_strategy = ifelse(is.na(treatment_strategy), "No ICI", treatment_strategy),
    ici_tto = ifelse(is.na(ici_tto), "No ICI", ici_tto)
  )

## 5.5 合并治疗类型信息####
message("合并治疗类型信息...")
treatment_simple <- treatment_by_primaryid %>%
  select(primaryid, radiotherapy, chemotherapy, targeted_therapy)

faers_malignant_final <- faers_malignant_final %>%
  left_join(treatment_simple, by = "primaryid") %>%
  mutate(
    radiotherapy = ifelse(is.na(radiotherapy), "No", radiotherapy),
    chemotherapy = ifelse(is.na(chemotherapy), "No", chemotherapy),
    targeted_therapy = ifelse(is.na(targeted_therapy), "No", targeted_therapy)
  )

## 5.6 合并irAE信息####
message("合并irAE信息...")
faers_malignant_final <- faers_malignant_final %>%
  left_join(irae_by_primaryid, by = "primaryid")

# 将irAE相关NA替换为0
irae_cols <- names(faers_malignant_final)[grepl("^(PT_|SOC_)", names(faers_malignant_final))]
for (col in irae_cols) {
  faers_malignant_final[[col]][is.na(faers_malignant_final[[col]])] <- 0L
}
faers_malignant_final$has_irae[is.na(faers_malignant_final$has_irae)] <- FALSE
faers_malignant_final$irae_pt_count[is.na(faers_malignant_final$irae_pt_count)] <- 0L
faers_malignant_final$irae_soc_count[is.na(faers_malignant_final$irae_soc_count)] <- 0L

## 5.7 变量规范化####
message("\n规范化描述性分析变量...")

# 1. 年龄分组（Age）
faers_malignant_final <- faers_malignant_final %>%
  mutate(
    Age = case_when(
      is.na(age_years) ~ "Unknown",
      age_years < 40 ~ "<40",
      age_years >= 40 & age_years < 60 ~ "40–<60",
      age_years >= 60 & age_years < 80 ~ "60–<80",
      age_years >= 80 ~ "≥80",
      TRUE ~ "Unknown"
    ),
    Age = factor(Age, levels = c("<40", "40–<60", "60–<80", "≥80", "Unknown"))
  )

# 2. 性别（Sex）- 已经是 Male/Female/Unknown，转为factor
faers_malignant_final <- faers_malignant_final %>%
  mutate(
    Sex = factor(sex_std, levels = c("Male", "Female", "Unknown"))
  )

# 3. 体重分组（Weight）
faers_malignant_final <- faers_malignant_final %>%
  mutate(
    Weight = case_when(
      is.na(wt) ~ "Unknown",
      wt < 50 ~ "<50",
      wt >= 50 & wt < 70 ~ "50–<70",
      wt >= 70 & wt < 90 ~ "70–<90",
      wt >= 90 ~ "≥90",
      TRUE ~ "Unknown"
    ),
    Weight = factor(Weight, levels = c("<50", "50–<70", "70–<90", "≥90", "Unknown"))
  )

# 4. 报告者类型（Reporter）
faers_malignant_final <- faers_malignant_final %>%
  mutate(
    Reporter = case_when(
      is.na(occp_cod) ~ "Unknown",
      occp_cod %in% c("MD", "PH", "OT") ~ "Health-professional",
      TRUE ~ "Non-health-professional"
    ),
    Reporter = factor(Reporter, levels = c("Health-professional", "Non-health-professional", "Unknown"))
  )

# 5. 国家（Country）- 动态计算前5国家
# 统计国家频次（NA和"COUNTRY NOT SPECIFIED"归类为Unknown）
country_counts <- faers_malignant_final %>%
  mutate(
    country_clean = case_when(
      is.na(reporter_country) ~ "Unknown",
      reporter_country == "COUNTRY NOT SPECIFIED" ~ "Unknown",
      TRUE ~ reporter_country
    )
  ) %>%
  count(country_clean, sort = TRUE)

# 获取前5国家（排除Unknown）
top5_countries <- country_counts %>%
  filter(country_clean != "Unknown") %>%
  head(5) %>%
  pull(country_clean)
message("前5报告国家: ", paste(top5_countries, collapse = ", "))

# 国家代码到全称的映射
country_name_map <- c(
  "US" = "United States",
  "JP" = "Japan",
  "FR" = "France",
  "CA" = "Canada",
  "DE" = "Germany",
  "GB" = "United Kingdom",
  "IT" = "Italy",
  "CN" = "China",
  "AU" = "Australia",
  "ES" = "Spain"
)

# 将top5国家代码转为全称
top5_country_names <- country_name_map[top5_countries]
top5_country_names <- ifelse(is.na(top5_country_names), top5_countries, top5_country_names)

# 创建Country分组
faers_malignant_final <- faers_malignant_final %>%
  mutate(
    Country = case_when(
      is.na(reporter_country) | reporter_country == "COUNTRY NOT SPECIFIED" ~ "Unknown",
      reporter_country %in% top5_countries ~ reporter_country,
      TRUE ~ "Other Countries"
    ),
    # 转换国家代码为全称
    Country = case_when(
      Country == "US" ~ "United States",
      Country == "JP" ~ "Japan",
      Country == "FR" ~ "France",
      Country == "CA" ~ "Canada",
      Country == "DE" ~ "Germany",
      Country == "GB" ~ "United Kingdom",
      Country == "IT" ~ "Italy",
      Country == "CN" ~ "China",
      Country == "AU" ~ "Australia",
      Country == "ES" ~ "Spain",
      TRUE ~ Country # 保留 Unknown 和 Other Countries 不变
    ),
    Country = factor(Country, levels = c(top5_country_names, "Other Countries", "Unknown"))
  )

# 6. ICI治疗（ICI_therapy）
faers_malignant_final <- faers_malignant_final %>%
  mutate(
    ICI_therapy = ifelse(has_ici == TRUE, "Yes", "No"),
    ICI_therapy = factor(ICI_therapy, levels = c("Yes", "No"))
  )

# 7. 治疗变量转为factor（radiotherapy/chemotherapy/targeted_therapy已是Yes/No字符串）
faers_malignant_final <- faers_malignant_final %>%
  mutate(
    Radiotherapy = factor(radiotherapy, levels = c("Yes", "No")),
    Chemotherapy = factor(chemotherapy, levels = c("Yes", "No")),
    Targeted_therapy = factor(targeted_therapy, levels = c("Yes", "No"))
  )

message("变量规范化完成")

# 清理临时对象
rm(country_counts, top5_countries, top5_country_names, country_name_map)

## 5.8 数据概览####
message("\n=== 恶性肿瘤样本概览 ===")
message(paste("总样本数:", nrow(faers_malignant_final)))
message(paste("变量数:", ncol(faers_malignant_final)))

message("\n队列分布:")
print(table(faers_malignant_final$cohort))

message("\nICI使用情况:")
print(table(faers_malignant_final$has_ici))

message("\n治疗策略分布:")
print(table(faers_malignant_final$treatment_strategy))

message("\n放疗使用情况:")
print(table(faers_malignant_final$radiotherapy))

message("\n化疗使用情况:")
print(table(faers_malignant_final$chemotherapy))

message("\n靶向治疗使用情况:")
print(table(faers_malignant_final$targeted_therapy))

message("\nirAE发生情况:")
print(table(faers_malignant_final$has_irae))

## 5.86 按癌种映射表进一步标准化三层标签####
message("按癌种映射表进一步标准化三层标签...")

cancer_mapping_ref <- openxlsx::read.xlsx(file.path(results_path, "癌种映射.xlsx"))
valid_organ_histology <- cancer_mapping_ref %>%
  transmute(
    primary_organ_ref = case_when(
      is.na(Primary_organs_tissues) | str_trim(Primary_organs_tissues) == "" ~ "Unknown/Others",
      TRUE ~ str_trim(Primary_organs_tissues)
    ),
    histologic_category_ref = case_when(
      is.na(Histologic_category) | str_trim(Histologic_category) == "" ~ "Unknown/Others",
      str_detect(str_trim(Histologic_category), regex("^(other|others|unknown|unknown/others)$", ignore_case = TRUE)) ~ "Unknown/Others",
      TRUE ~ str_trim(Histologic_category)
    )
  ) %>%
  filter(histologic_category_ref != "Unknown/Others") %>%
  distinct() %>%
  mutate(.valid_histology = TRUE)

faers_malignant_final <- faers_malignant_final %>%
  mutate(
    histologic_category = case_when(
      is.na(histologic_category) | str_trim(histologic_category) == "" ~ "Unknown/Others",
      str_detect(str_trim(histologic_category), regex("^(other|others|unknown/others)$", ignore_case = TRUE)) ~ "Unknown/Others",
      TRUE ~ str_trim(histologic_category)
    ),
    assigned_type = case_when(
      is.na(assigned_type) | str_trim(assigned_type) == "" ~ "Unknown/Others",
      str_detect(str_trim(assigned_type), regex("^(other|others|na|n/a|unknown|unknown/others)$", ignore_case = TRUE)) ~ "Unknown/Others",
      TRUE ~ str_trim(assigned_type)
    )
  ) %>%
  left_join(
    valid_organ_histology,
    by = c("primary_organ" = "primary_organ_ref", "histologic_category" = "histologic_category_ref")
  ) %>%
  mutate(
    histologic_category = if_else(
      is.na(.valid_histology) & histologic_category != "Unknown/Others",
      "Unknown/Others", histologic_category
    ),
    assigned_type = if_else(
      is.na(.valid_histology),
      "Unknown/Others", assigned_type
    )
  ) %>%
  select(-.valid_histology)

rm(cancer_mapping_ref, valid_organ_histology)

## 5.87 按业务规则修正Eye第二层分类####
faers_malignant_final <- faers_malignant_final %>%
  mutate(
    histologic_category = case_when(
      primary_organ == "Eye" & assigned_type == "UVM" ~ "Uveal Melanoma",
      primary_organ == "Eye" ~ "Unknown/Others",
      TRUE ~ histologic_category
    )
  ) %>%
  mutate(
    histologic_category = case_when(
      primary_organ == "Eye" &
        histologic_category == "Unknown/Others" &
        str_detect(coalesce(indi_pt_list, ""), regex("RETINOBLASTOMA", ignore_case = TRUE)) ~ "Retinoblastoma",
      TRUE ~ histologic_category
    )
  )
## 5.88 按业务规则拆分Rare并修正Unmapped的第一层分类####
faers_malignant_final <- faers_malignant_final %>%
  mutate(
    primary_organ = case_when(
      primary_organ == "Rare" & assigned_type == "ACC" ~ "Adrenal gland",
      primary_organ == "Rare" & assigned_type == "TGCT" ~ "Testis",
      primary_organ == "Rare" & first_pt %in% c(
        "Adrenal Gland Cancer",
        "Adrenocortical Carcinoma"
      ) ~ "Adrenal gland",
      primary_organ == "Rare" & str_detect(
        coalesce(first_pt, ""),
        regex("testicular|testis", ignore_case = TRUE)
      ) ~ "Testis",
      primary_organ == "Rare" & first_pt %in% c(
        "Penile Cancer",
        "Penile Squamous Cell Carcinoma",
        "Penis Carcinoma",
        "Carcinoma In Situ Of Penis"
      ) ~ "Penis",
      primary_organ == "Rare" & first_pt %in% c(
        "Vulval Cancer",
        "Genital Neoplasm Malignant Female",
        "Vaginal Cancer",
        "Squamous Cell Carcinoma Of The Vulva",
        "Squamous Cell Carcinoma Of The Vagina",
        "Vulvar Adenocarcinoma",
        "Vaginal Adenocarcinoma",
        "Vulvar Basal Cell Carcinoma"
      ) ~ "Vulva/Vagina",
      primary_organ == "Rare" & first_pt %in% c(
        "Malignant Neoplasm Of Ampulla Of Vater",
        "Malignant Neoplasm Papilla Of Vater"
      ) ~ "Ampulla of Vater",
      primary_organ == "Rare" ~ "CUP",
      primary_organ == "Unmapped" ~ "Hematologic",
      TRUE ~ primary_organ
    )
  ) %>%
  mutate(histologic_category = case_when(
    primary_organ == "Testis" ~ "Testis Cancer",
    primary_organ == "Penis" ~ "Penile Cancer",
    primary_organ == "Ampulla of vater" ~ "Ampullary cancer",
    primary_organ == "Vulva/Vagina" & first_pt %in% c(
      "Vulval Cancer",
      "Squamous Cell Carcinoma Of The Vulva",
      "Vulvar Adenocarcinoma",
      "Vulvar Basal Cell Carcinoma"
    ) ~ "Vulval Cancer",
    primary_organ == "Vulva/Vagina" & first_pt %in% c(
      "Vaginal Cancer",
      "Squamous Cell Carcinoma Of The Vagina",
      "Vaginal Adenocarcinoma"
    ) ~ "Vaginal Cancer",
    primary_organ == "Vulva/Vagina" & first_pt %in% c(
      "Genital Neoplasm Malignant Female"
    ) ~ "Others",
    primary_organ == "Adrenal gland" & first_pt %in% c(
      "Adrenocortical Carcinoma"
    ) ~ "Adrenocortical Carcinoma",
    primary_organ == "Adrenal gland" & first_pt %in% c(
      "Adrenal Gland Cancer"
    ) ~ "Others",
    TRUE ~ histologic_category
  )) %>%
  mutate(assigned_type = case_when(
    histologic_category == "Testis Cancer" ~ "TGCT",
    TRUE ~ assigned_type
  ))
## 5.9 保存最终数据####
save(faers_malignant_final, file = file.path(.output_path, "faers_malignant_final.RData"))
faers_solid <- faers_malignant_final %>%
  filter(!primary_organ %in% c("CUP", "Hematologic", "Unmapped"))
save(faers_solid, file = file.path(.output_path, "faers_solid.RData"))
message(paste("\n恶性肿瘤样本已保存:", file.path(.output_path, "faers_malignant_final.RData")))

# 输出列名供参考
message("\n数据列名:")
print(names(faers_malignant_final))


# 清理环境，准备处理JADER数据####
rm(list = ls())
gc()


# 加载通用路径配置（自动设置 data_path, output_path, results_path）
source("Code/config.R")

# 加载共用筛选规则和映射规则
source("Code/mapping_rules.R")

# 第一部分：筛选JADER癌症患者####

## 1.1 加载JADER数据####
message("加载JADER数据...")
load(file.path(.output_path, "jader_demo.RData"))
load(file.path(.output_path, "jader_drug.RData"))
load(file.path(.output_path, "jader_reac.RData"))
load(file.path(.output_path, "jader_hist.RData"))

message(paste("JADER Demo:", nrow(jader_demo), "条"))
message(paste("JADER Drug:", nrow(jader_drug), "条"))
message(paste("JADER Reac:", nrow(jader_reac), "条"))
message(paste("JADER Hist:", nrow(jader_hist), "条"))

## 1.2 癌症筛选规则（已从mapping_rules.R加载）####
# 以下规则现已移至 Code/mapping_rules.R 统一管理：
# - primary_strong, generalized_histology, generalized_malignant
# - metastasis_pattern, met_to_site_pattern, malignant_confirmed_pattern
# - benign_pattern, uncertain_pattern, generalized_only_pattern

## 1.3 JADER使用hist表进行PT级标记（与FAERS一致）####
# JADER的hist表类似于FAERS的indi表，记录患者既往病史/原疾患
# 注意：数据处理.R中已将日文PT转换为英文，字段名为medical_history_std

message("\n从hist表进行PT级标记...")

# 检查是否已有英文PT字段
if (!"medical_history_std" %in% names(jader_hist)) {
  jader_hist <- jader_hist %>%
    mutate(medical_history_std = toupper(trimws(medical_history)))
  message("注意：使用日文medical_history，建议先运行数据处理.R进行日英转换")
}

# 按primaryid + medical_history_std去重
message(paste("去重前hist记录数:", nrow(jader_hist)))
jader_hist_dedup <- jader_hist %>%
  distinct(primaryid, medical_history_std, .keep_all = TRUE)
message(paste("去重后hist记录数:", nrow(jader_hist_dedup)))

# 使用mapping_rules.R中定义的排除PT列表（已在脚本开头加载）
# benign_tumor_pts: 良性/交界性肿瘤PT（误匹配BLASTOMA等关键词）
# non_malignant_pts: 非肿瘤疾病PT（误匹配METASTATIC关键词）
exclude_pts <- toupper(c(benign_tumor_pts, non_malignant_pts))
message(paste("  待排除的非恶性PT数量:", length(exclude_pts)))

# 设置并行计算（使用所有可用核心）
# Linux系统使用multicore（fork模式）更高效，避免socket传输大数据
n_cores <- parallel::detectCores() - 1
plan(multicore, workers = n_cores)
message(paste("使用", n_cores, "个核心进行并行计算"))

# 并行进行PT级标记 + 癌种映射
message("开始并行PT级标记和癌种映射...")
start_time <- Sys.time()

message("  计算hit_malignant...")
hit_malignant_vec <- future_sapply(jader_hist_dedup$medical_history_std, function(pt) {
  grepl(malignant_confirmed_pattern, pt, perl = TRUE)
}, USE.NAMES = FALSE)

message("  计算hit_benign...")
hit_benign_vec <- future_sapply(jader_hist_dedup$medical_history_std, function(pt) {
  grepl(benign_pattern, pt, perl = TRUE)
}, USE.NAMES = FALSE)

message("  计算hit_uncertain...")
hit_uncertain_vec <- future_sapply(jader_hist_dedup$medical_history_std, function(pt) {
  grepl(uncertain_pattern, pt, perl = TRUE)
}, USE.NAMES = FALSE)

message("  计算hit_metastasis...")
hit_metastasis_vec <- future_sapply(jader_hist_dedup$medical_history_std, function(pt) {
  grepl(metastasis_pattern, pt, perl = TRUE)
}, USE.NAMES = FALSE)

message("  计算hit_generalized...")
hit_generalized_vec <- future_sapply(jader_hist_dedup$medical_history_std, function(pt) {
  grepl(generalized_only_pattern, pt, perl = TRUE)
}, USE.NAMES = FALSE)

message("  计算hit_false_positive（假阳性排除）...")
hit_false_positive_vec <- future_sapply(jader_hist_dedup$medical_history_std, function(pt) {
  grepl(false_positive_pattern, pt, perl = TRUE)
}, USE.NAMES = FALSE)

message("  计算hit_met_to_site（转移到某部位）...")
hit_met_to_site_vec <- future_sapply(jader_hist_dedup$medical_history_std, function(pt) {
  grepl(met_to_site_pattern, pt, perl = TRUE)
}, USE.NAMES = FALSE)

message("  计算cancer_type...")
cancer_type_vec <- future_sapply(jader_hist_dedup$medical_history_std, map_cancer_type, USE.NAMES = FALSE)

message("  计算is_excluded（排除非恶性PT）...")
is_excluded_vec <- future_sapply(jader_hist_dedup$medical_history_std, function(pt) {
  toupper(pt) %in% exclude_pts
}, USE.NAMES = FALSE)

# 组装结果
# 重要修正：转移状态不等于原发不明
# - hit_metastasis：仅表示"有转移关键词"，作为分层变量
# - hit_met_to_site：表示"转移到某部位"模式，这类才更可能原发不明
# - cancer_type映射不受hit_metastasis影响，只受hit_met_to_site影响（已在map_cancer_type中处理）

# 并行计算 medical_history_normalized
message("  计算medical_history_normalized（PT标准化）...")
medical_history_normalized_vec <- future_sapply(jader_hist_dedup$medical_history_std, function(pt) {
  # 标准化PT名称：
  # 1. 保护HIGH/LOW GRADE
  # 2. 删除STAGE/GRADE及其后面内容
  # 3. 删除状态修饰词（METASTATIC、RECURRENT等），确保同一癌种的不同状态被识别为同一器官
  # 4. 恢复保护的内容
  pt_temp <- gsub("(HIGH|LOW)([ -])GRADE", "\\1\\2XPROTECTEDX", pt, ignore.case = TRUE)
  pt_temp <- trimws(sub("\\s*(STAGE|GRADE).*", "", pt_temp, ignore.case = TRUE))
  # 删除状态修饰词（确保"COLORECTAL CANCER METASTATIC"和"COLORECTAL CANCER"被识别为同一器官）
  pt_temp <- trimws(sub("\\s+(METASTATIC|RECURRENT|REFRACTORY|IN\\s+REMISSION|LOCALISED|LOCALIZED|REGIONAL|EXTENSIVE|LIMITED)$", "", pt_temp, ignore.case = TRUE))
  pt_temp <- gsub("XPROTECTEDX", "GRADE", pt_temp)
  return(pt_temp)
}, USE.NAMES = FALSE)

# 向量化计算 hit_malignant_final
# 排除非恶性PT（如Lipoblastoma、Metastatic Cutaneous Crohn's Disease等）
hit_malignant_final <- ifelse(is_excluded_vec, FALSE, hit_malignant_vec)
# 排除假阳性PT（如Metastatic Pulmonary Embolism等非癌症诊断）
hit_malignant_final <- ifelse(hit_false_positive_vec, FALSE, hit_malignant_final)
# 如果有癌种映射，则标记为恶性
hit_malignant_final <- ifelse(!is.na(cancer_type_vec), TRUE, hit_malignant_final)

# 预计算器官映射（用于未映射PT的合并判定）
# 关键修正：同一器官的多个PT应合并为单一原发癌，避免被归入other_malignancy
message("  计算primary_organ_for_grouping（用于合并同器官PT）...")
primary_organ_vec <- future_mapply(
  function(pt_norm, cancer) {
    map_to_primary_organ(pt_norm, cancer)
  },
  medical_history_normalized_vec,
  cancer_type_vec,
  SIMPLIFY = TRUE,
  USE.NAMES = FALSE
)

# 并行计算 type_fine 和 type_coarse
# 关键修正：type_fine使用Histologic_category（严格按照Results/癌种映射.xlsx）
# 同一Histologic_category下的PT会合并为同一type_fine
message("  计算type_fine和type_coarse...")
type_results <- future_mapply(
  function(cancer, hit_mal, hit_mts, hit_gen, pt_norm, organ) {
    # type_fine：用于癌种判定
    # 严格按照用户的三层映射体系：
    # 1. 通过TCGA代码或PT获取Histologic_category
    # 2. 无法映射但有转移关键词 → INDETERMINATE
    # 3. 无法映射但确认恶性 → 使用器官默认的Histologic_category
    type_fine <- if (hit_mal && hit_mts) {
      "INDETERMINATE"
    } else if (hit_mal && hit_gen) {
      "INDETERMINATE"
    } else if (hit_mal && !hit_gen && !hit_mts) {
      # 使用map_to_histologic_category获取Histologic_category
      histcat <- map_to_histologic_category(pt_norm, cancer, organ)
      if (!is.na(histcat)) {
        histcat
      } else if (!is.na(organ) && organ != "Unknown" && organ != "CUP") {
        paste0("UNMAPPED:", organ)
      } else {
        paste0("UNMAPPED_PRIMARY:", pt_norm)
      }
    } else {
      NA_character_
    }

    # type_coarse：最终癌种展示（使用Primary_organs_tissues）
    type_coarse <- if (!is.na(cancer)) {
      cancer
    } else if (hit_mal && hit_mts) {
      "INDETERMINATE"
    } else if (hit_mal && hit_gen) {
      "INDETERMINATE"
    } else if (hit_mal && !hit_gen && !hit_mts) {
      "Others"
    } else {
      NA_character_
    }

    return(c(type_fine, type_coarse))
  },
  cancer_type_vec,
  hit_malignant_final,
  hit_met_to_site_vec,
  hit_generalized_vec,
  medical_history_normalized_vec,
  primary_organ_vec,
  SIMPLIFY = TRUE,
  USE.NAMES = FALSE
)

type_fine_vec <- type_results[1, ]
type_coarse_vec <- type_results[2, ]

# 组装最终结果
jader_hist_tagged <- jader_hist_dedup %>%
  mutate(
    hit_malignant = hit_malignant_final,
    hit_metastasis = hit_metastasis_vec,
    hit_met_to_site = hit_met_to_site_vec,
    hit_benign = hit_benign_vec,
    hit_uncertain = hit_uncertain_vec,
    hit_generalized = hit_generalized_vec,
    cancer_type = cancer_type_vec,
    medical_history_normalized = medical_history_normalized_vec,
    type_fine = type_fine_vec,
    type_coarse = type_coarse_vec,
    # 添加PT级器官映射（用于判断同器官多癌种情况）
    primary_organ_pt = primary_organ_vec
  )

# 清理临时向量
rm(
  hit_malignant_vec, hit_benign_vec, hit_uncertain_vec, hit_metastasis_vec, hit_generalized_vec,
  hit_false_positive_vec, hit_met_to_site_vec, cancer_type_vec, is_excluded_vec,
  hit_malignant_final, medical_history_normalized_vec, type_results, type_fine_vec, type_coarse_vec,
  primary_organ_vec
)

end_time <- Sys.time()
message(paste("PT级标记和癌种映射完成，耗时:", round(difftime(end_time, start_time, units = "mins"), 2), "分钟"))

plan(sequential)
gc()

# 统计PT级各分类
message(paste("确认恶性记录数:", sum(jader_hist_tagged$hit_malignant)))
message(paste("  其中转移相关:", sum(jader_hist_tagged$hit_metastasis)))
message(paste("明确良性记录数:", sum(jader_hist_tagged$hit_benign)))
message(paste("良恶性不明记录数:", sum(jader_hist_tagged$hit_uncertain)))
message(paste("可映射癌种记录数:", sum(!is.na(jader_hist_tagged$cancer_type))))

## 1.5 PRIMARYID级互斥归类（与FAERS一致的6类分类）####

message("\n开始PRIMARYID级互斥归类...")

# 使用data.table进行高效并行分组聚合（data.table内置多线程优化）
n_cores <- parallel::detectCores() - 1
setDTthreads(n_cores)
message(paste("使用", n_cores, "个data.table线程进行并行分组聚合"))

start_time <- Sys.time()

# 转换为data.table
hist_dt <- as.data.table(jader_hist_tagged)

# 分组聚合（data.table自动多线程）
primaryid_summary_dt <- hist_dt[, .(
  # 基础信号
  has_malignant = any(hit_malignant),
  has_benign = any(hit_benign),
  has_uncertain = any(hit_uncertain),
  # 纯不明信号（排除恶性重叠）
  has_uncertain_only = any(hit_uncertain & !hit_malignant),
  # metastasis_only临时标记：所有适应证都包含转移关键词
  metastasis_only_flag = all(hit_metastasis),
  # 统计type_fine数量（排除INDETERMINATE）
  n_types_fine = uniqueN(type_fine[!is.na(type_fine) & type_fine != "INDETERMINATE"]),
  n_types_tcga = uniqueN(cancer_type[!is.na(cancer_type)]),
  # 是否只有INDETERMINATE
  has_only_indeterminate = {
    type_fine_nonNA <- type_fine[!is.na(type_fine)]
    length(type_fine_nonNA) > 0 && all(type_fine_nonNA == "INDETERMINATE")
  },
  # 获取唯一type_coarse（排除INDETERMINATE，临时存储）
  assigned_type_temp = {
    valid_types <- type_coarse[!is.na(type_coarse) & type_coarse != "INDETERMINATE"]
    if (length(valid_types) > 0) unique(valid_types)[1] else NA_character_
  },
  # 获取唯一type_fine（排除INDETERMINATE，临时存储）
  assigned_type_fine_temp = {
    valid_types <- type_fine[!is.na(type_fine) & type_fine != "INDETERMINATE"]
    if (length(valid_types) > 0) unique(valid_types)[1] else NA_character_
  },
  # 获取癌种列表（用于描述）
  types_list = paste(unique(type_fine[!is.na(type_fine)]), collapse = ";"),
  # 新增：计算恶性肿瘤PT的唯一器官数量（用于同器官多癌种判断）
  # 例如：GLIOBLASTOMA和BRAIN NEOPLASM MALIGNANT都属于Brain，则unique_organs=1
  n_unique_organs = {
    malignant_organs <- primary_organ_pt[hit_malignant & !is.na(primary_organ_pt)]
    uniqueN(malignant_organs)
  },
  # 新增：获取同器官情况下的器官名称
  same_organ_temp = {
    malignant_organs <- primary_organ_pt[hit_malignant & !is.na(primary_organ_pt)]
    unique_organs <- unique(malignant_organs)
    if (length(unique_organs) == 1) unique_organs[1] else NA_character_
  }
), by = primaryid]

# 计算派生字段（metastasis_only需要结合has_malignant）
primaryid_summary_dt[, `:=`(
  metastasis_only = has_malignant & metastasis_only_flag,
  assigned_type = fifelse(n_types_fine == 1, assigned_type_temp, NA_character_),
  assigned_type_fine = fifelse(n_types_fine == 1, assigned_type_fine_temp, NA_character_)
)]

# 丢弃判定（互斥）
primaryid_summary_dt[, `:=`(
  benign_only = has_benign & !has_malignant & !has_uncertain,
  indeterminate_only = has_uncertain & !has_malignant
)]

# 最终分区（5类互斥，按优先级）使用data.table的fcase
primaryid_summary_dt[, cohort := fcase(
  # A. 非肿瘤组（排除）
  !has_malignant & !benign_only & !indeterminate_only, "non_tumor",
  # B. 排除：确定良性（纯良性，无恶性/不明信号）
  benign_only, "exclude_benign",
  # C. 描述性分析：良恶性不明（有borderline信号但无恶性）
  indeterminate_only, "descriptive_indeterminate",
  # D. 主分析：单一恶性肿瘤（n_types_fine=1，包括TCGA映射和Others）
  has_malignant & n_types_fine == 1, "single_malignancy",
  # E. 其他恶性肿瘤（多癌种共存或原发不可判定）
  has_malignant & n_types_fine != 1, "other_malignancy",
  # 其他（理论上不应出现）
  default = "other"
)]

# 删除临时列（保留same_organ_temp用于后续other_malignancy的器官分配）
primaryid_summary_dt[, c("metastasis_only_flag", "assigned_type_temp", "assigned_type_fine_temp") := NULL]

# 转换回tibble（保持与后续代码兼容）
primaryid_summary <- as_tibble(primaryid_summary_dt)

# 清理临时对象
rm(hist_dt, primaryid_summary_dt)

end_time <- Sys.time()
message(paste("PRIMARYID级归类完成，耗时:", round(difftime(end_time, start_time, units = "secs"), 2), "秒"))

# 保存中间数据集，方便后续调用（避免重复运行耗时的处理步骤）
message("\n保存中间数据集...")
save(primaryid_summary, file = file.path(.output_path, "jader_primaryid_summary.RData"))
save(jader_hist_tagged, file = file.path(.output_path, "jader_hist_tagged.RData"))
message(paste("  primaryid_summary已保存:", file.path(.output_path, "jader_primaryid_summary.RData")))
message(paste("  jader_hist_tagged已保存:", file.path(.output_path, "jader_hist_tagged.RData")))

# 统计各分区（数据处理阶段已保证四表共有主键）
message("\n========== PRIMARYID级分区统计 ==========")
message("\n分区分布:")
print(table(primaryid_summary$cohort))

# 主分析癌种分布
message("\n主分析癌种分布:")
main_type_dist <- primaryid_summary %>%
  filter(cohort == "single_malignancy") %>%
  count(assigned_type, sort = TRUE)
print(main_type_dist)

## 1.55 第一层器官映射（进一步精细化处理）####
message("\n========== 第一层器官映射 ==========")

# 加载中间数据集（如果从此处开始运行）
if (!exists("jader_hist_tagged")) {
  message("加载中间数据集 jader_hist_tagged...")
  load(file.path(.output_path, "jader_hist_tagged.RData"))
}
if (!exists("primaryid_summary")) {
  message("加载中间数据集 primaryid_summary...")
  load(file.path(.output_path, "jader_primaryid_summary.RData"))
}

# 确保映射规则已加载（包含select_primary_pt函数）
# 这对于分段运行或并行计算时函数可用性至关重要
if (!exists("select_primary_pt")) {
  message("加载映射规则（Code/mapping_rules.R）...")
  source("Code/mapping_rules.R")
}

# 从jader_hist_tagged获取肿瘤相关PT，构建first_pt
# 获取每个primaryid的肿瘤相关PT（只保留匹配malignant_confirmed_pattern的PT）
# 注意：排除空字符串的PT，避免first_pt为空
pt_tumor_by_primaryid_temp <- jader_hist_tagged %>%
  filter(hit_malignant, !is.na(medical_history_std), medical_history_std != "") %>%
  group_by(primaryid) %>%
  summarise(hist_pt_list_tumor = paste(unique(medical_history_std), collapse = " & "), .groups = "drop")

# 对single_malignancy的病例进行第一层器官映射
single_malignancy_for_mapping <- primaryid_summary %>%
  filter(cohort == "single_malignancy") %>%
  left_join(pt_tumor_by_primaryid_temp, by = "primaryid")

# 检查是否有无法获取肿瘤PT的病例（通常是空medical_history_std导致）
n_missing_tumor_pt <- sum(is.na(single_malignancy_for_mapping$hist_pt_list_tumor))
if (n_missing_tumor_pt > 0) {
  message(paste("  警告：有", n_missing_tumor_pt, "个single_malignancy病例缺少肿瘤PT信息（medical_history_std为空）"))
}

# 使用智能PT选择函数（优先选择原发癌种PT，排除转移部位PT）
# 注意：由于select_primary_pt依赖多个全局函数，使用顺序sapply更稳定
message("计算first_pt...")
start_time_pt <- Sys.time()
single_malignancy_for_mapping <- single_malignancy_for_mapping %>%
  mutate(first_pt = sapply(hist_pt_list_tumor, select_primary_pt, USE.NAMES = FALSE))
end_time_pt <- Sys.time()
message(paste("first_pt计算完成，耗时:", round(difftime(end_time_pt, start_time_pt, units = "secs"), 2), "秒"))


# 应用第一层器官映射
message("应用第一层器官映射...")
message("  映射优先级: 1. assigned_type (TCGA代码) -> 2. first_pt (PT名称)")

single_malignancy_for_mapping <- single_malignancy_for_mapping %>%
  mutate(primary_organ = mapply(
    map_to_primary_organ,
    first_pt = first_pt,
    assigned_type = assigned_type,
    USE.NAMES = FALSE
  ))

# 统计映射结果
organ_dist <- single_malignancy_for_mapping %>%
  filter(!is.na(primary_organ)) %>%
  count(primary_organ, sort = TRUE) %>%
  mutate(pct = round(n / sum(n) * 100, 2))
message("\n已映射的器官分布:")
print(organ_dist)

# 未映射的first_pt统计
unmapped <- single_malignancy_for_mapping %>%
  filter(is.na(primary_organ)) %>%
  count(first_pt, sort = TRUE)
message(paste("\n未映射的first_pt数量:", nrow(unmapped)))
if (nrow(unmapped) > 0 & nrow(unmapped) <= 20) {
  message("未映射的first_pt列表（前20个）:")
  print(head(unmapped, 20))
}
if (nrow(unmapped) > 20) {
  message(paste("完整未映射PT列表已保存至:", file.path(.results_path, "JADER_unmapped_pts.csv")))
  write.csv(unmapped, file.path(.results_path, "JADER_unmapped_pts.csv"), row.names = FALSE)
}

# 将first_pt和primary_organ合并回primaryid_summary
primaryid_summary <- primaryid_summary %>%
  left_join(
    single_malignancy_for_mapping %>% select(primaryid, first_pt, primary_organ),
    by = "primaryid"
  ) %>%
  # 确保所有恶性肿瘤病例都有器官映射
  mutate(
    primary_organ = case_when(
      # other_malignancy：先检查是否所有PT属于同一器官
      # 例如：GLIOBLASTOMA和BRAIN NEOPLASM MALIGNANT都属于Brain
      # 此时primary_organ归类为Brain而不是Unknown_Primary，cohort仍为other_malignancy
      cohort == "other_malignancy" & !is.na(same_organ_temp) ~ same_organ_temp,
      # other_malignancy但多个器官或无法确认原发灶：标记为Unknown_Primary
      cohort == "other_malignancy" & is.na(same_organ_temp) ~ "Unknown_Primary",
      # single_malignancy但未能映射到具体器官：标记为Unmapped
      cohort == "single_malignancy" & is.na(primary_organ) ~ "Unmapped",
      # 已有映射的保持不变
      TRUE ~ primary_organ
    )
  )

# 将first_pt和primary_organ暂存用于后续映射
single_malignancy_for_mapping_temp <- single_malignancy_for_mapping %>%
  select(primaryid, first_pt, primary_organ, assigned_type)

# 清理临时对象
rm(pt_tumor_by_primaryid_temp, single_malignancy_for_mapping, organ_dist, unmapped)
gc()

## 1.56 第二层组织学分类映射####
message("\n========== 第二层组织学分类映射 ==========")

# 对single_malignancy的病例进行第二层组织学分类映射
single_malignancy_for_mapping <- single_malignancy_for_mapping_temp

single_malignancy_for_mapping <- single_malignancy_for_mapping %>%
  mutate(histologic_category = mapply(
    map_to_histologic_category,
    first_pt = first_pt,
    assigned_type = assigned_type,
    primary_organ = primary_organ,
    USE.NAMES = FALSE
  ))

# 统计映射结果
histology_dist <- single_malignancy_for_mapping %>%
  filter(!is.na(histologic_category)) %>%
  count(histologic_category, sort = TRUE) %>%
  mutate(pct = round(n / sum(n) * 100, 2))
message("\n已映射的组织学分类分布:")
print(histology_dist)

# 未映射的first_pt统计（第二层）
unmapped_histology <- single_malignancy_for_mapping %>%
  filter(is.na(histologic_category)) %>%
  count(first_pt, sort = TRUE)
message(paste("\n未映射组织学分类的first_pt数量:", nrow(unmapped_histology)))
if (nrow(unmapped_histology) > 0) {
  message("未映射组织学分类的first_pt列表（按频次排序，前20个）:")
  print(head(unmapped_histology, 20))
}

# 将histologic_category合并回primaryid_summary
# 首先检查single_malignancy_for_mapping中是否有histologic_category列
if ("histologic_category" %in% names(single_malignancy_for_mapping)) {
  primaryid_summary <- primaryid_summary %>%
    left_join(
      single_malignancy_for_mapping %>% select(primaryid, histologic_category),
      by = "primaryid"
    )
} else {
  warning("single_malignancy_for_mapping中没有histologic_category列，可能需要先运行第1.56节的组织学分类映射代码")
  primaryid_summary <- primaryid_summary %>%
    mutate(histologic_category = NA_character_)
}

# 确保histologic_category列存在（即使left_join没有添加）
if (!"histologic_category" %in% names(primaryid_summary)) {
  primaryid_summary <- primaryid_summary %>%
    mutate(histologic_category = NA_character_)
}

# 对于非single_malignancy的病例，histologic_category设为Unknown
primaryid_summary <- primaryid_summary %>%
  mutate(
    histologic_category = case_when(
      cohort == "other_malignancy" ~ "Unknown",
      cohort == "single_malignancy" & is.na(histologic_category) ~ "Unknown",
      TRUE ~ histologic_category
    )
  )

# 清理临时对象
rm(single_malignancy_for_mapping, single_malignancy_for_mapping_temp, histology_dist, unmapped_histology)
gc()

## 1.6 保存队列映射####

# 加载中间数据集（如果从此处开始运行）
if (!exists("jader_hist_tagged")) {
  message("加载中间数据集 jader_hist_tagged...")
  load(file.path(.output_path, "jader_hist_tagged.RData"))
}
if (!exists("primaryid_summary")) {
  message("加载中间数据集 primaryid_summary...")
  load(file.path(.output_path, "jader_primaryid_summary.RData"))
}

# 获取每个primaryid的原疾患PT（去重后用&合并，收集所有PT）
pt_by_primaryid <- jader_hist_tagged %>%
  group_by(primaryid) %>%
  summarise(hist_pt_list = paste(unique(medical_history_std), collapse = " & "), .groups = "drop")

# 获取每个primaryid的肿瘤相关PT（只保留匹配malignant_confirmed_pattern的PT）
pt_tumor_by_primaryid <- jader_hist_tagged %>%
  filter(hit_malignant) %>%
  group_by(primaryid) %>%
  summarise(hist_pt_list_tumor = paste(unique(medical_history_std), collapse = " & "), .groups = "drop")

# 保存主键+队列+癌种标注+原疾患PT映射表（包含first_pt、primary_organ和histologic_category）
jader_cohort_mapping <- primaryid_summary %>%
  select(primaryid, cohort, assigned_type, first_pt, primary_organ, histologic_category) %>%
  left_join(pt_by_primaryid, by = "primaryid") %>%
  left_join(pt_tumor_by_primaryid, by = "primaryid")
jader_cohort_mapping <- jader_cohort_mapping %>%
  mutate(primary_organ = if_else(cohort == "other_malignancy", "CUP", primary_organ))

save(jader_cohort_mapping, file = file.path(.output_path, "jader_cohort_mapping.RData"))
message(paste("\nJADER队列映射已保存:", file.path(.output_path, "jader_cohort_mapping.RData")))
message(paste("共", nrow(jader_cohort_mapping), "条记录"))

# 定义恶性肿瘤患者primaryid（用于后续ICI筛选）
tumor_primaryids <- primaryid_summary %>%
  filter(cohort %in% c("single_malignancy", "other_malignancy")) %>%
  pull(primaryid)
message(paste("恶性肿瘤患者数:", length(tumor_primaryids)))

# 第二部分：ICI患者筛选####

## 2.1 定义ICI药物列表（含日文药名）####
message("\n\n========== 第二部分：ICI患者筛选 ==========")

# PD-1抑制剂
drug_PD_1 <- list(
  NIVOLUMAB = c(
    "NIVOLUMAB", "ＮＩＶＯＬＵＭＡＢ", "ＢＬＩＮＤＥＤ　ＮＩＶＯＬＵＭＡＢ",
    "NIVOLUMAB;PLACEBO", "ニボルマブ", "ニボルマブ（遺伝子組換え）",
    "OPDIVO for intravenous injection", "ＯＰＤＩＶＯ", "OPDUALAG", "OPDIVO"
  ),
  PEMBROLIZUMAB = c(
    "PEMBROLIZUMAB", "ペンブロリズマブ", "ペムブロリズマブ",
    "PEMBROLIZUMAB (GENETICAL RECOMBINATION)", "ペムブロリズマブ(遺伝子組換え)注",
    "PEMBROLIZUMAB(GENETICAL RECOMBINATION)", "ペムブロリズマブ（遺伝子組換え）",
    "Keytruda", "KEYTRUDA"
  ),
  CEMIPLIMAB = c("CEMIPLIMAB", "セミプリマブ", "LIBTAYO"),
  DOSTARLIMAB = c("DOSTARLIMAB", "JEMPERLI"),
  RETIFANLIMAB = c("RETIFANLIMAB", "ZYNYZ"),
  TORIPALIMAB = c("TORIPALIMAB", "LOQTORZI"),
  TISLELIZUMAB = c("TISLELIZUMAB", "TEVIMBRA")
)

# PD-L1抑制剂
drug_PD_L1 <- list(
  ATEZOLIZUMAB = c(
    "ATEZOLIZUMAB", "ＡＴＥＺＯＬＩＺＵＭＡＢ", "ＥＸＴＥＲＮＡＬ−ＡＴＥＺＯＬＩＺＵＭＡＢ",
    "ATEZOLIZUMAB(GENETICAL RECOMBINATION)", "Atezolizumab (Genetical recombination)",
    "アテゾリズマブ（遺伝子組換え）", "アテゾリズマブ", "テセントリク（遺伝子組換え）",
    "テセントリク点滴静注", "TECENTRIQ", "TECENTRIQ HYBREZA", "抗PD-L1抗体"
  ),
  AVELUMAB = c(
    "AVELUMAB", "ＡＶＥＬＵＭＡＢ", "アベルマブ", "アベルマブ（遺伝子組換え）",
    "バベンシオ", "BAVENCIO"
  ),
  DURVALUMAB = c(
    "DURVALUMAB", "ＤＵＲＶＡＬＵＭＡＢ", "デュルバルマブ",
    "デュルバルマブ（遺伝子組換え）", "IMFINZI"
  )
)

# CTLA-4抑制剂
drug_CTLA_4 <- list(
  IPILIMUMAB = c(
    "IPILIMUMAB", "ＩＰＩＬＩＭＵＭＡＢ", "イピリムマブ",
    "イピリムマブ（遺伝子組換え）", "YERVOY"
  ),
  TREMELIMUMAB = c(
    "TREMELIMUMAB", "ＴＲＥＭＥＬＩＭＵＭＡＢ", "トレメリムマブ",
    "トレメリムマブ（遺伝子組換え）", "IMJUDO"
  )
)

# 合并所有ICI药物
drug_ICI_all <- c(drug_PD_1, drug_PD_L1, drug_CTLA_4)

## 2.2 ICI匹配函数（支持日文药名）####
match_ICI <- function(drugname_std, brand_name, drug_list) {
  drugname_std <- ifelse(is.na(drugname_std) || drugname_std == "", "", drugname_std)
  brand_name <- ifelse(is.na(brand_name) || brand_name == "", "", brand_name)

  drugname_std <- toupper(drugname_std)
  brand_name <- toupper(brand_name)

  for (ici_name in names(drug_list)) {
    patterns <- drug_list[[ici_name]]
    for (p in patterns) {
      # 日文药名使用fixed匹配，英文药名使用正则匹配
      if (grepl("[^\x01-\x7F]", p)) {
        # 包含非ASCII字符（日文），使用fixed匹配
        if (grepl(p, drugname_std, fixed = TRUE) ||
          grepl(p, brand_name, fixed = TRUE)) {
          return(ici_name)
        }
      } else {
        # 纯英文，使用正则边界匹配
        if (grepl(paste0("(?i)\\b", p, "\\b"), drugname_std, perl = TRUE) ||
          grepl(paste0("(?i)\\b", p, "\\b"), brand_name, perl = TRUE)) {
          return(ici_name)
        }
      }
    }
  }
  return(NA_character_)
}

## 2.3 筛选恶性肿瘤患者的ICI使用####
malignant_primaryids <- tumor_primaryids

# 筛选恶性肿瘤患者的用药记录
jader_drug_malignant <- jader_drug %>%
  filter(primaryid %in% malignant_primaryids)

message(paste("恶性肿瘤患者用药记录数:", nrow(jader_drug_malignant)))

# 设置并行计算（使用multicore避免socket传输问题）
library(future.apply)
n_cores <- parallel::detectCores() - 1
plan(multicore, workers = n_cores)
message(paste("使用", n_cores, "个核心进行并行计算"))

# 并行匹配ICI药物
message("并行匹配ICI药物...")
start_time <- Sys.time()

ici_drug_vec <- future_mapply(
  function(d, b) match_ICI(d, b, drug_ICI_all),
  jader_drug_malignant$generic_name_std,
  jader_drug_malignant$brand_name,
  USE.NAMES = FALSE
)

end_time <- Sys.time()
message(paste("ICI匹配耗时:", round(difftime(end_time, start_time, units = "mins"), 2), "分钟"))

# 恢复顺序计划
plan(sequential)

# 筛选ICI用药记录
jader_drug_ici <- jader_drug_malignant %>%
  mutate(ici_drug = ici_drug_vec) %>%
  filter(!is.na(ici_drug))

message(paste("ICI用药记录数:", nrow(jader_drug_ici)))
message(paste("使用ICI的患者数:", uniqueN(jader_drug_ici$primaryid)))

# 按primaryid去重（同一患者使用同一ICI药物只计一次）
jader_drug_ici_dedup <- jader_drug_ici %>%
  distinct(primaryid, ici_drug, .keep_all = TRUE)

message(paste("去重后ICI用药记录数:", nrow(jader_drug_ici_dedup)))

## 2.4 ICI类型分类####
# 定义ICI类型映射
ici_type_map <- data.frame(
  ici_drug = c(sapply(drug_PD_1, `[`, 1), sapply(drug_PD_L1, `[`, 1), sapply(drug_CTLA_4, `[`, 1)),
  ici_type = c(
    rep("anti-PD-1", length(drug_PD_1)),
    rep("anti-PD-L1", length(drug_PD_L1)),
    rep("anti-CTLA-4", length(drug_CTLA_4))
  )
)

jader_drug_ici_typed <- jader_drug_ici_dedup %>%
  left_join(ici_type_map, by = "ici_drug")

# 统计各ICI药物使用情况
message("\nICI药物使用统计:")
print(table(jader_drug_ici_typed$ici_drug))

message("\nICI类型统计:")
print(table(jader_drug_ici_typed$ici_type))

## 2.5 按primaryid汇总ICI用药####
message("\n汇总每个患者的ICI使用情况...")

jader_ici_by_primaryid <- jader_drug_ici_typed %>%
  group_by(primaryid) %>%
  summarise(
    ici_drugs = paste(sort(unique(ici_drug)), collapse = "+"),
    ici_types = paste(sort(unique(ici_type)), collapse = "+"),
    n_ici_drugs = n_distinct(ici_drug),
    n_ici_types = n_distinct(ici_type),
    .groups = "drop"
  ) %>%
  mutate(
    treatment_strategy = case_when(
      n_ici_types > 1 ~ "Combination Therapy",
      TRUE ~ ici_types
    ),
    ici_tto = "Unknown"
  )

message("\nICI治疗策略分布:")
print(table(jader_ici_by_primaryid$treatment_strategy))

## 2.6 保存ICI映射####
save(jader_ici_by_primaryid, file = file.path(.output_path, "jader_ici_mapping.RData"))
message(paste("\nJADER ICI映射已保存:", file.path(.output_path, "jader_ici_mapping.RData")))

# 第三部分：识别放疗、化疗和靶向治疗病例####

# 路径已在config.R中设置，如需独立运行请先执行: source("Code/config.R")

# 加载队列映射（如未加载）
if (!exists("jader_cohort_mapping") || is.null(jader_cohort_mapping)) {
  load(file.path(.output_path, "jader_cohort_mapping.RData"))
}

# 定义恶性肿瘤患者
malignant_primaryids <- jader_cohort_mapping %>%
  filter(cohort %in% c("single_malignancy", "other_malignancy")) %>%
  pull(primaryid)

# 加载drug数据（如未加载）
if (!exists("jader_drug") || is.null(jader_drug)) {
  load(file.path(.output_path, "jader_drug.RData"))
}

# 筛选恶性肿瘤患者的用药记录
jader_drug_malignant <- jader_drug %>%
  filter(primaryid %in% malignant_primaryids)
message(paste("恶性肿瘤患者用药记录数:", nrow(jader_drug_malignant)))

message("\n\n========== 第三部分：识别放疗、化疗和靶向治疗 ==========")

## 3.1 定义放疗关键词（英文+日文）####

radiotherapy_strong <- c(
  "RADIATION( THERAPY| TREATMENT)?", "RADIOTHERAPY", "IRRADIATION",
  "BRACHYTHERAPY", "RADIOSURGERY",
  "EXTERNAL BEAM", "PROTON( THERAPY| BEAM)?",
  "GAMMA KNIFE", "CYBER ?KNIFE", "TOMOTHERAPY",
  "LINEAR ACCELERATOR", "\\bLINAC\\b",
  "STEREOTACTIC( RADIOTHERAPY| RADIOSURGERY)?",
  "INTRAOPERATIVE RADIOTHERAPY", "\\bIORT\\b",
  "WHOLE BRAIN", "CRANIAL IRRADIATION", "PELVIC IRRADIATION",
  "HEAVY ION", "CARBON ION", "PARTICLE THERAPY",
  "X-RAY THERAPY", "PHOTON THERAPY", "ELECTRON BEAM"
)
radiotherapy_tech <- c("\\bIMRT\\b", "\\bIGRT\\b", "\\bVMAT\\b", "\\bSRS\\b", "\\bSBRT\\b")
# 日文放疗关键词
radiotherapy_jp <- c(
  "放射線療法", "放射線治療", "照射", "放射線",
  "小線源治療", "密封小線源", "ガンマナイフ", "サイバーナイフ",
  "陽子線", "重粒子線", "炭素イオン", "定位照射", "定位放射線"
)
radiotherapy_pattern <- paste0(
  "(?i)(", paste(radiotherapy_strong, collapse = "|"), ")",
  "|(?i)(", paste(radiotherapy_tech, collapse = "|"), ")",
  "|(", paste(radiotherapy_jp, collapse = "|"), ")"
)

# 化疗药物（英文+日文）
chemo_drugs <- c(
  "DAUNORUBICIN", "CYTARABINE", "ARSENIC", "ASPARAGINASE", "AZACITIDINE",
  "BENDAMUSTINE", "BUSULFAN", "CAPECITABINE", "CARBOPLATIN", "CARMUSTINE",
  "CHLORAMBUCIL", "CISPLATIN", "CALDRIBINE", "CLOFARABINE", "CYCLOPHOSPHAMIDE",
  "DACARBAZINE", "DECITABINE", "DEGARELIX", "DOCETAXEL", "ERIBULIN",
  "ETOPOSIDE", "FLOXIRIDINE", "FLUDARABINE", "FLUOROURACIL", "GEMCITABINE",
  "GOSERELIN", "HISTRELIN", "HYDROXYUREA", "IFOSFAMIDE", "INTERFERON",
  "IRINOTECAN", "IXABEPILONE", "LANREOTIDE", "LEUCOVORIN", "LOMUSTINE",
  "MEDROXYPROGESTERONE", "MEGESTROL", "MELPHALAN", "MERCAPTOPURINE", "MESNA",
  "METHOTREXATE", "MITOXANTRONE", "NELARABINE", "OCTREOTIDE", "OMACETAXINE",
  "OXALIPLATIN", "PACLITAXEL", "PEGASPARGASE", "PEMETREXED", "PENTOSTATIN",
  "PROCARBAZINE", "SAMARIUM", "LEXIDRONAM", "TEMOZOLOMIDE", "TENIPOSIDE",
  "THIOTEPA", "TOPOTECAN", "TRABECTEDIN", "VINBLASTINE", "VINORELBINE",
  "DOXORUBICIN", "IDARUBICIN", "EPIRUBICIN", "BLEOMYCIN", "CABAZITAXEL",
  "DACTINOMYCIN", "MITOMYCIN", "STREPTOZOCIN", "5-FU", "FLUOROURACIL"
)
# 日文化疗药名
chemo_drugs_jp <- c(
  "シスプラチン", "カルボプラチン", "オキサリプラチン",
  "パクリタキセル", "ドセタキセル", "ゲムシタビン",
  "フルオロウラシル", "カペシタビン", "イリノテカン",
  "エトポシド", "ドキソルビシン", "エピルビシン",
  "メトトレキサート", "シクロホスファミド", "イホスファミド",
  "ビンブラスチン", "ビノレルビン", "ビンクリスチン",
  "ブレオマイシン", "マイトマイシン", "ペメトレキセド",
  "テモゾロミド", "アザシチジン", "デシタビン"
)
chemo_pattern <- paste0(
  "(?i)(", paste(chemo_drugs, collapse = "|"), ")",
  "|(", paste(chemo_drugs_jp, collapse = "|"), ")"
)

# 靶向治疗药物（英文+日文）
targeted_drugs <- c(
  "TRASTUZUMAB", "AFLIBERCEPT", "ALEMTUZUMAB", "ANASTROZOLE", "BELIMUMAB",
  "BELINOSTAT", "BEVACIZUMAB", "BLINATUMOMAB", "BORTEZOMIB", "BRENTUXIMAB",
  "CARFILZOMIB", "CETUXIMAB", "COPANLISIB", "DENOSUMAB", "ELOTUZUMAB",
  "EVEROLIMUS", "EXEMESTANE", "FULVESTRANT", "GEFITINIB", "GEMTUZUMAB",
  "IMATINIB", "INOTUZUMAB", "MOGAMULIZUMAB", "MOXETUMOMAB", "PASUDOTOX",
  "NATALIZUMAB", "NECITUMUMAB", "OBINUTUZUMAB", "OFATUMUMAB", "OLARATUMAB",
  "PANITUMUMAB", "PERTUZUMAB", "PRALATREXATE", "RADIUM", "RAMUCIRUMAB",
  "RITUXIMAB", "ROMIDEPSIN", "SILTUXIMAB", "TAGRAXOFUSP", "TALAZOPARIB",
  "TAMOXIFEN", "TEMSIROLIMUS", "VALRUBICIN", "ZIV-AFLIBERCEPT", "BICALUTAMIDE",
  "FLOXURIDINE", "FLUTAMIDE", "KETOCONAZOLE", "LANREOTIDE", "LENALIDOMIDE",
  "LEVOLEUCOVORIN", "POMALIDOMIDE", "PORFIMER", "PROPYLTHIOURACIL",
  "THALIDOMIDE", "TRIFLURIDINE", "VINCRISTINE", "ABEMACICLIB", "ABIRATERONE",
  "ACALABRUTINIB", "AFATINIB", "ALECTINIB", "ALPELISIB", "APALUTAMIDE",
  "AXITINIB", "BINIMETINIB", "BOSUTINIB", "BRIGATINIB", "CABOZANTINIB",
  "CERITINIB", "COBIMETINIB", "CRIZOTINIB", "DABRAFENIB", "DACOMITINIB",
  "DARATUMUMAB", "DAROLUTAMIDE", "DASATINIB", "DINUTUXIMAB", "DUVELISIB",
  "ELOTUZUMAB", "ENASIDENIB", "ENCORAFENIB", "ENFORTUMAB", "ENTRECTINIB",
  "ENZALUTAMIDE", "ERDAFITINIB", "ERLOTINIB", "FAM-TRASTUZUMAB", "DERUXTECA",
  "FEDRATINIB", "FULVESTRANT", "GILTERITINIB", "GLASDEGIB", "IBRUTINIB",
  "IDELALISIB", "ISATUXIMAB", "ISOTRETINOIN", "IVOSIDENIB", "IXAZOMIB",
  "LAPATINIB", "LAROTRECTINIB", "LENVATINIB", "LETROZOLE", "LORLATINIB",
  "NERATINIB", "NILOTINIB", "NIRAPARIB", "OBINUTUZUMAB", "OLAPARIB",
  "OLARATUMAB", "OSIMERTINIB", "PALBOCICLIB", "PANOBINOSTAT", "PAZOPANIB",
  "PEXIDARTINIB", "POLATUZUMAB", "PONATINIB", "REGORAFENIB", "RIBOCICLIB",
  "RUCAPARIB", "RUXOLITINIB", "MIDOSTAURIN", "SELINEXOR", "SONIDEGIB",
  "SORAFENIB", "SUNITINIB", "TALAZOPARIB", "TOSYKATE", "TAZEMETOSTAT",
  "TOREMIFENE", "TRAMETINIB", "TRETINOIN", "VANDETANIB", "VEMURAFENIB",
  "VENETOCLAX", "VISMODEGIB", "VORINOSTAT", "ZANUBRUTINIB"
)
# 日文靶向药名
targeted_drugs_jp <- c(
  "トラスツズマブ", "ベバシズマブ", "リツキシマブ", "セツキシマブ",
  "パニツムマブ", "ペルツズマブ", "ラムシルマブ", "デノスマブ",
  "イマチニブ", "ゲフィチニブ", "エルロチニブ", "アファチニブ",
  "オシメルチニブ", "クリゾチニブ", "アレクチニブ", "ロルラチニブ",
  "ダブラフェニブ", "トラメチニブ", "ベムラフェニブ", "エンコラフェニブ",
  "ソラフェニブ", "スニチニブ", "レンバチニブ", "レゴラフェニブ",
  "パゾパニブ", "アキシチニブ", "カボザンチニブ",
  "ボルテゾミブ", "カルフィルゾミブ", "イキサゾミブ",
  "エベロリムス", "テムシロリムス",
  "オラパリブ", "ニラパリブ", "ルカパリブ", "タラゾパリブ",
  "パルボシクリブ", "リボシクリブ", "アベマシクリブ",
  "ベネトクラクス", "イブルチニブ", "アカラブルチニブ",
  "タモキシフェン", "レトロゾール", "アナストロゾール", "エキセメスタン",
  "フルベストラント", "エンザルタミド", "アビラテロン", "アパルタミド",
  "レナリドミド", "ポマリドミド", "サリドマイド"
)
targeted_pattern <- paste0(
  "(?i)(", paste(targeted_drugs, collapse = "|"), ")",
  "|(", paste(targeted_drugs_jp, collapse = "|"), ")"
)

# 并行标记治疗类型
message("并行标记治疗类型...")

# 设置并行计算（使用multicore避免socket传输问题）
library(future.apply)
n_cores <- parallel::detectCores() - 1
plan(multicore, workers = n_cores)
message(paste("使用", n_cores, "个核心进行并行计算"))

start_time <- Sys.time()

message("  标记放疗...")
is_radiotherapy_vec <- future_mapply(
  function(d, b) {
    grepl(radiotherapy_pattern, d, perl = TRUE) |
      grepl(radiotherapy_pattern, b, perl = TRUE)
  },
  jader_drug_malignant$generic_name_std,
  jader_drug_malignant$brand_name,
  USE.NAMES = FALSE
)

message("  标记化疗...")
is_chemotherapy_vec <- future_mapply(
  function(d, b) {
    grepl(chemo_pattern, d, perl = TRUE) |
      grepl(chemo_pattern, b, perl = TRUE)
  },
  jader_drug_malignant$generic_name_std,
  jader_drug_malignant$brand_name,
  USE.NAMES = FALSE
)

message("  标记靶向治疗...")
is_targeted_vec <- future_mapply(
  function(d, b) {
    grepl(targeted_pattern, d, perl = TRUE) |
      grepl(targeted_pattern, b, perl = TRUE)
  },
  jader_drug_malignant$generic_name_std,
  jader_drug_malignant$brand_name,
  USE.NAMES = FALSE
)

# 添加标记列
jader_drug_treatment <- jader_drug_malignant %>%
  mutate(
    is_radiotherapy = is_radiotherapy_vec,
    is_chemotherapy = is_chemotherapy_vec,
    is_targeted = is_targeted_vec
  )

end_time <- Sys.time()
message(paste("治疗类型标记完成，耗时:", round(difftime(end_time, start_time, units = "mins"), 2), "分钟"))

# 恢复顺序计划
plan(sequential)

# 汇总到primaryid级别
jader_treatment_by_primaryid <- jader_drug_treatment %>%
  group_by(primaryid) %>%
  summarise(
    has_radiotherapy = any(is_radiotherapy),
    has_chemotherapy = any(is_chemotherapy),
    has_targeted = any(is_targeted),
    .groups = "drop"
  ) %>%
  mutate(
    radiotherapy = ifelse(has_radiotherapy, "Yes", "No"),
    chemotherapy = ifelse(has_chemotherapy, "Yes", "No"),
    targeted_therapy = ifelse(has_targeted, "Yes", "No")
  )

message("\n治疗类型统计:")
message(paste("  放疗:", sum(jader_treatment_by_primaryid$has_radiotherapy)))
message(paste("  化疗:", sum(jader_treatment_by_primaryid$has_chemotherapy)))
message(paste("  靶向:", sum(jader_treatment_by_primaryid$has_targeted)))

## 3.2 保存治疗类型映射####
save(jader_treatment_by_primaryid, file = file.path(.output_path, "jader_treatment_mapping.RData"))
message(paste("\nJADER治疗类型映射已保存:", file.path(.output_path, "jader_treatment_mapping.RData")))

# 第四部分：标记irAE不良反应####

# 路径已在config.R中设置，如需独立运行请先执行: source("Code/config.R")

# 加载队列映射（如未加载）
if (!exists("jader_cohort_mapping") || is.null(jader_cohort_mapping)) {
  load(file.path(.output_path, "jader_cohort_mapping.RData"))
}

# 定义恶性肿瘤患者
malignant_primaryids <- jader_cohort_mapping %>%
  filter(cohort %in% c("single_malignancy", "other_malignancy")) %>%
  pull(primaryid)

message("\n\n========== 第四部分：标记irAE不良反应 ==========")

## 4.1 加载irAE列表####

library(readxl)
irae_list <- read_excel(file.path(.data_path, "irAE列表.xlsx"))
message(paste("irAE PT数量:", nrow(irae_list)))

# 提取irAE PT列表（标准化为大写）
irae_pts <- toupper(trimws(irae_list$PT))

# 筛选恶性肿瘤患者的不良反应记录
if (!exists("jader_reac") || is.null(jader_reac)) {
  load(file.path(.output_path, "jader_reac.RData"))
}

jader_reac_malignant <- jader_reac %>%
  filter(primaryid %in% malignant_primaryids)
message(paste("恶性肿瘤患者不良反应记录数:", nrow(jader_reac_malignant)))

# 并行标记irAE
message("并行标记irAE...")
n_cores <- parallel::detectCores() - 1
plan(multicore, workers = n_cores)
message(paste("使用", n_cores, "个核心进行并行计算"))
start_time <- Sys.time()

is_irae_vec <- future_sapply(
  jader_reac_malignant$pt_std,
  function(x) x %in% irae_pts,
  USE.NAMES = FALSE
)

end_time <- Sys.time()
message(paste("irAE标记完成，耗时:", round(difftime(end_time, start_time, units = "mins"), 2), "分钟"))

plan(sequential)

# 合并标记结果
jader_reac_irae <- jader_reac_malignant %>%
  mutate(is_irae = is_irae_vec)

# 创建PT到SOC的映射
pt_to_soc <- setNames(irae_list$SOC, toupper(trimws(irae_list$PT)))

# 为每条irAE记录添加SOC
jader_reac_irae <- jader_reac_irae %>%
  mutate(irae_soc = ifelse(is_irae, pt_to_soc[pt_std], NA_character_))

# 获取所有唯一的PT和SOC
unique_pts <- sort(unique(irae_pts))
unique_socs <- sort(unique(irae_list$SOC))

message(paste("irAE PT数量:", length(unique_pts)))
message(paste("irAE SOC数量:", length(unique_socs)))

# 筛选命中irAE的记录
irae_hits <- jader_reac_irae %>%
  filter(is_irae) %>%
  select(primaryid, pt_std, irae_soc)

# 创建PT宽格式（每个PT一列，0/1）
message("创建PT宽格式...")
irae_pt_wide <- irae_hits %>%
  distinct(primaryid, pt_std) %>%
  mutate(pt_clean = gsub(" ", "_", pt_std)) %>%
  select(primaryid, pt_clean) %>%
  distinct() %>%
  mutate(value = 1L) %>%
  pivot_wider(
    names_from = pt_clean,
    values_from = value,
    values_fill = 0L,
    names_prefix = "PT_"
  )

# 创建SOC宽格式（每个SOC一列，0/1）
message("创建SOC宽格式...")
irae_soc_wide <- irae_hits %>%
  distinct(primaryid, irae_soc) %>%
  mutate(soc_clean = gsub(" ", "_", irae_soc)) %>%
  select(primaryid, soc_clean) %>%
  distinct() %>%
  mutate(value = 1L) %>%
  pivot_wider(
    names_from = soc_clean,
    values_from = value,
    values_fill = 0L,
    names_prefix = "SOC_"
  )

# 获取所有恶性肿瘤患者primaryid（包括无irAE的）
all_malignant_ids <- data.frame(primaryid = malignant_primaryids)

# 合并PT和SOC宽格式，补全无irAE患者为0
jader_irae_by_primaryid <- all_malignant_ids %>%
  left_join(irae_pt_wide, by = "primaryid") %>%
  left_join(irae_soc_wide, by = "primaryid")

# 将NA替换为0
jader_irae_by_primaryid[is.na(jader_irae_by_primaryid)] <- 0L

# 添加汇总列
jader_irae_by_primaryid <- jader_irae_by_primaryid %>%
  mutate(
    has_irae = rowSums(select(., starts_with("PT_"))) > 0,
    irae_pt_count = rowSums(select(., starts_with("PT_"))),
    irae_soc_count = rowSums(select(., starts_with("SOC_")))
  )

# 统计结果
message("\nirAE发生情况:")
print(table(jader_irae_by_primaryid$has_irae))

message("\nirAE PT数量分布（有irAE的患者）:")
print(summary(jader_irae_by_primaryid$irae_pt_count[jader_irae_by_primaryid$has_irae]))

## 4.2 保存irAE标注结果####
save(jader_irae_by_primaryid, file = file.path(.output_path, "jader_irae_mapping.RData"))
message(paste("\nJADER irAE标注表已保存:", file.path(.output_path, "jader_irae_mapping.RData")))

# 第五部分：合并所有标记并保存恶性肿瘤样本####
message("\n\n========== 第五部分：合并所有标记并保存 ==========")

## 5.1 加载所有映射表####
# 队列映射
if (!exists("jader_cohort_mapping") || is.null(jader_cohort_mapping)) {
  load(file.path(.output_path, "jader_cohort_mapping.RData"))
}
# ICI映射
if (!exists("jader_ici_by_primaryid") || is.null(jader_ici_by_primaryid)) {
  load(file.path(.output_path, "jader_ici_mapping.RData"))
}
# 治疗类型映射
if (!exists("jader_treatment_by_primaryid") || is.null(jader_treatment_by_primaryid)) {
  load(file.path(.output_path, "jader_treatment_mapping.RData"))
}
# irAE映射
if (!exists("jader_irae_by_primaryid") || is.null(jader_irae_by_primaryid)) {
  load(file.path(.output_path, "jader_irae_mapping.RData"))
}

# 定义恶性肿瘤患者
malignant_primaryids <- jader_cohort_mapping %>%
  filter(cohort %in% c("other_malignancy", "single_malignancy")) %>%
  pull(primaryid)

## 5.2 筛选恶性肿瘤患者demo####
# 加载demo数据（如未加载）
if (!exists("jader_demo") || is.null(jader_demo)) {
  load(file.path(.output_path, "jader_demo.RData"))
}

jader_demo_malignant <- jader_demo %>%
  filter(primaryid %in% malignant_primaryids)

message(paste("恶性肿瘤患者demo记录数:", nrow(jader_demo_malignant)))

## 5.3 合并队列、ICI、治疗类型和irAE信息####
jader_malignant_final <- jader_demo_malignant %>%
  left_join(jader_cohort_mapping, by = "primaryid") %>%
  left_join(jader_ici_by_primaryid, by = "primaryid") %>%
  left_join(jader_treatment_by_primaryid, by = "primaryid") %>%
  left_join(jader_irae_by_primaryid, by = "primaryid") %>%
  mutate(
    has_ici = !is.na(ici_drugs),
    ici_drugs = ifelse(is.na(ici_drugs), "None", ici_drugs),
    ici_types = ifelse(is.na(ici_types), "None", ici_types),
    treatment_strategy = ifelse(is.na(treatment_strategy), "No ICI", treatment_strategy),
    ici_tto = ifelse(is.na(ici_tto), "No ICI", ici_tto),
    has_radiotherapy = ifelse(is.na(has_radiotherapy), FALSE, has_radiotherapy),
    has_chemotherapy = ifelse(is.na(has_chemotherapy), FALSE, has_chemotherapy),
    has_targeted = ifelse(is.na(has_targeted), FALSE, has_targeted),
    radiotherapy = ifelse(is.na(radiotherapy), "No", radiotherapy),
    chemotherapy = ifelse(is.na(chemotherapy), "No", chemotherapy),
    targeted_therapy = ifelse(is.na(targeted_therapy), "No", targeted_therapy),
    Targeted_therapy = ifelse(is.na(Targeted_therapy), targeted_therapy, Targeted_therapy),
    source = "JADER"
  )

# 将irAE相关NA替换为0
irae_cols <- names(jader_malignant_final)[grepl("^(PT_|SOC_)", names(jader_malignant_final))]
for (col in irae_cols) {
  jader_malignant_final[[col]][is.na(jader_malignant_final[[col]])] <- 0L
}
jader_malignant_final$has_irae[is.na(jader_malignant_final$has_irae)] <- FALSE
jader_malignant_final$irae_pt_count[is.na(jader_malignant_final$irae_pt_count)] <- 0L
jader_malignant_final$irae_soc_count[is.na(jader_malignant_final$irae_soc_count)] <- 0L

## 5.4 数据概览####
message("\n=== JADER恶性肿瘤样本概览 ===")
message(paste("总样本数:", nrow(jader_malignant_final)))

message("\n队列分布:")
print(table(jader_malignant_final$cohort, useNA = "ifany"))

message("\nICI使用情况:")
print(table(jader_malignant_final$has_ici))

message("\n治疗策略分布:")
print(table(jader_malignant_final$treatment_strategy))

message("\n治疗类型统计:")
message(paste("  放疗:", sum(jader_malignant_final$has_radiotherapy)))
message(paste("  化疗:", sum(jader_malignant_final$has_chemotherapy)))
message(paste("  靶向:", sum(jader_malignant_final$has_targeted)))

message("\nirAE发生情况:")
print(table(jader_malignant_final$has_irae))

## 5.5 变量规范化（与FAERS保持一致）####
message("\n规范化描述性分析变量...")

jader_malignant_final <- jader_malignant_final %>%
  mutate(
    # 年龄分组（JADER年龄为年代分类，如10歳代=10-19岁）
    Age = case_when(
      is.na(age) | age %in% c("不明", "") ~ "Unknown",
      # <40岁：10歳代, 20歳代, 30歳代, 10歳未満, 及特殊分类
      age %in% c(
        "10歳未満", "10歳代", "20歳代", "30歳代",
        "乳幼児", "新生児", "小児", "胎児", "青少年", "成人"
      ) ~ "<40",
      # 40-59岁：40歳代, 50歳代
      age %in% c("40歳代", "50歳代") ~ "40–<60",
      # 60-79岁：60歳代, 70歳代
      age %in% c("60歳代", "70歳代") ~ "60–<80",
      # ≥80岁：80歳代, 90歳代, 100歳代, 高齢者
      age %in% c("80歳代", "90歳代", "100歳代", "高齢者") ~ "≥80",
      TRUE ~ "Unknown"
    ),
    Age = factor(Age, levels = c("<40", "40–<60", "60–<80", "≥80", "Unknown")),
    # 性别
    Sex = case_when(
      is.na(sex) ~ "Unknown",
      sex %in% c("男性") ~ "Male",
      sex %in% c("女性") ~ "Female",
      TRUE ~ "Unknown"
    ),
    Sex = factor(Sex, levels = c("Male", "Female", "Unknown")),
    # 体重分组（JADER使用日文分类，如50kg台=50-59kg）
    Weight = case_when(
      is.na(weight) | weight == "" ~ "Unknown",
      # <50kg: 10kg未満, 10kg台, 20kg台, 30kg台, 40kg台
      weight %in% c("10kg未満", "10kg台", "20kg台", "30kg台", "40kg台") ~ "<50",
      # 50-69kg: 50kg台, 60kg台
      weight %in% c("50kg台", "60kg台") ~ "50–<70",
      # 70-89kg: 70kg台, 80kg台
      weight %in% c("70kg台", "80kg台") ~ "70–<90",
      # ≥90kg: 90kg台及以上
      weight %in% c(
        "90kg台", "100kg台", "110kg台", "120kg台", "130kg台",
        "140kg台", "150kg台", "160kg台", "170kg台", "180kg台",
        "200kg以上"
      ) ~ "≥90",
      TRUE ~ "Unknown"
    ),
    Weight = factor(Weight, levels = c("<50", "50–<70", "70–<90", "≥90", "Unknown")),
    # 报告者类型（JADER使用reporter_qual，包含医師、薬剤師等日文分类）
    # 包含医師/薬剤師/医療関係者的归为Health-professional
    Reporter = case_when(
      is.na(reporter_qual) | reporter_qual == "不明" ~ "Unknown",
      grepl("医師|薬剤師|医療関係者", reporter_qual) ~ "Health-professional",
      TRUE ~ "Non-health-professional"
    ),
    Reporter = factor(Reporter, levels = c("Health-professional", "Non-health-professional", "Unknown")),
    # 身高分组（JADER特有，FAERS无此变量）
    Height = case_when(
      is.na(height) | height == "" ~ "Unknown",
      # <150cm
      height %in% c(
        "10cm未満", "30cm台", "40cm台", "50cm台", "60cm台",
        "70cm台", "80cm台", "90cm台", "100cm台", "110cm台",
        "120cm台", "130cm台", "140cm台"
      ) ~ "<150",
      # 150-159cm
      height %in% c("150cm台") ~ "150–<160",
      # 160-169cm
      height %in% c("160cm台") ~ "160–<170",
      # ≥170cm
      height %in% c("170cm台", "180cm台", "190cm台") ~ "≥170",
      TRUE ~ "Unknown"
    ),
    Height = factor(Height, levels = c("<150", "150–<160", "160–<170", "≥170", "Unknown")),
    # 国家（JADER是日本国内数据库，所有记录都是Japan）
    Country = factor("Japan", levels = c("Japan")),
    # 治疗变量（直接使用已有的has_ici, has_radiotherapy等变量）
    # ICI治疗
    `ICI therapy` = ifelse(has_ici == TRUE, "Yes", "No"),
    `ICI therapy` = factor(`ICI therapy`, levels = c("Yes", "No")),
    # 放疗
    Radiotherapy = ifelse(has_radiotherapy == TRUE, "Yes", "No"),
    Radiotherapy = factor(Radiotherapy, levels = c("Yes", "No")),
    # 化疗
    Chemotherapy = ifelse(has_chemotherapy == TRUE, "Yes", "No"),
    Chemotherapy = factor(Chemotherapy, levels = c("Yes", "No")),
    # 靶向治疗
    `Targeted therapy` = ifelse(has_targeted == TRUE, "Yes", "No"),
    `Targeted therapy` = factor(`Targeted therapy`, levels = c("Yes", "No")),
    # 标准别名列，便于与FAERS共用分析脚本
    Targeted_therapy = factor(targeted_therapy, levels = c("Yes", "No"))
  )

message("变量规范化完成")

# 规范化后的变量统计
message("\n规范化变量统计:")
message("年龄分布:")
print(table(jader_malignant_final$Age))
message("\n性别分布:")
print(table(jader_malignant_final$Sex))
message("\n体重分布:")
print(table(jader_malignant_final$Weight))
message("\n身高分布:")
print(table(jader_malignant_final$Height))
message("\n报告者类型分布:")
print(table(jader_malignant_final$Reporter))
message("\nICI治疗分布:")
print(table(jader_malignant_final$`ICI therapy`))
message("\n靶向治疗分布:")
print(table(jader_malignant_final$`Targeted therapy`))

## 5.56 按癌种映射表进一步标准化三层标签####
message("按癌种映射表进一步标准化三层标签...")

cancer_mapping_ref <- openxlsx::read.xlsx(file.path(results_path, "癌种映射.xlsx"))
valid_organ_histology <- cancer_mapping_ref %>%
  transmute(
    primary_organ_ref = case_when(
      is.na(Primary_organs_tissues) | str_trim(Primary_organs_tissues) == "" ~ "Unknown/Others",
      TRUE ~ str_trim(Primary_organs_tissues)
    ),
    histologic_category_ref = case_when(
      is.na(Histologic_category) | str_trim(Histologic_category) == "" ~ "Unknown/Others",
      str_detect(str_trim(Histologic_category), regex("^(other|others|unknown|unknown/others)$", ignore_case = TRUE)) ~ "Unknown/Others",
      TRUE ~ str_trim(Histologic_category)
    )
  ) %>%
  filter(histologic_category_ref != "Unknown/Others") %>%
  distinct() %>%
  mutate(.valid_histology = TRUE)

jader_malignant_final <- jader_malignant_final %>%
  mutate(
    histologic_category = case_when(
      is.na(histologic_category) | str_trim(histologic_category) == "" ~ "Unknown/Others",
      str_detect(str_trim(histologic_category), regex("^(other|others|unknown/others)$", ignore_case = TRUE)) ~ "Unknown/Others",
      TRUE ~ str_trim(histologic_category)
    ),
    assigned_type = case_when(
      is.na(assigned_type) | str_trim(assigned_type) == "" ~ "Unknown/Others",
      str_detect(str_trim(assigned_type), regex("^(other|others|na|n/a|unknown|unknown/others)$", ignore_case = TRUE)) ~ "Unknown/Others",
      TRUE ~ str_trim(assigned_type)
    )
  ) %>%
  left_join(
    valid_organ_histology,
    by = c("primary_organ" = "primary_organ_ref", "histologic_category" = "histologic_category_ref")
  ) %>%
  mutate(
    histologic_category = if_else(
      is.na(.valid_histology) & histologic_category != "Unknown/Others",
      "Unknown/Others", histologic_category
    ),
    assigned_type = if_else(
      is.na(.valid_histology),
      "Unknown/Others", assigned_type
    )
  ) %>%
  select(-.valid_histology)

rm(cancer_mapping_ref, valid_organ_histology)

## 5.57 按业务规则修正Eye第二层分类####
jader_malignant_final <- jader_malignant_final %>%
  mutate(
    histologic_category = case_when(
      primary_organ == "Eye" & assigned_type == "UVM" ~ "Uveal Melanoma",
      primary_organ == "Eye" ~ "Unknown/Others",
      TRUE ~ histologic_category
    )
  ) %>%
  mutate(
    histologic_category = case_when(
      primary_organ == "Eye" &
        histologic_category == "Unknown/Others" &
        str_detect(coalesce(hist_pt_list, ""), regex("RETINOBLASTOMA", ignore_case = TRUE)) ~ "Retinoblastoma",
      TRUE ~ histologic_category
    )
  )

## 5.58 按业务规则拆分Rare并修正Unmapped的第一层分类####
jader_malignant_final <- jader_malignant_final %>%
  mutate(
    primary_organ = case_when(
      primary_organ == "Rare" & assigned_type == "ACC" ~ "Adrenal gland",
      primary_organ == "Rare" & assigned_type == "TGCT" ~ "Testis",
      primary_organ == "Rare" & first_pt %in% c(
        "Adrenal Gland Cancer",
        "Adrenocortical Carcinoma"
      ) ~ "Adrenal gland",
      primary_organ == "Rare" & str_detect(
        coalesce(first_pt, ""),
        regex("testicular|testis", ignore_case = TRUE)
      ) ~ "Testis",
      primary_organ == "Rare" & first_pt %in% c(
        "Penile Cancer",
        "Penile Squamous Cell Carcinoma",
        "Penis Carcinoma",
        "Carcinoma In Situ Of Penis"
      ) ~ "Penis",
      primary_organ == "Rare" & first_pt %in% c(
        "Vulval Cancer",
        "Genital Neoplasm Malignant Female",
        "Vaginal Cancer",
        "Squamous Cell Carcinoma Of The Vulva",
        "Squamous Cell Carcinoma Of The Vagina",
        "Vulvar Adenocarcinoma",
        "Vaginal Adenocarcinoma",
        "Vulvar Basal Cell Carcinoma"
      ) ~ "Vulva/Vagina",
      primary_organ == "Rare" & first_pt %in% c(
        "Malignant Neoplasm Of Ampulla Of Vater",
        "Malignant Neoplasm Papilla Of Vater"
      ) ~ "Ampulla of Vater",
      primary_organ == "Rare" ~ "CUP",
      primary_organ == "Unmapped" ~ "Hematologic",
      TRUE ~ primary_organ
    )
  ) %>%
  mutate(histologic_category = case_when(
    primary_organ == "Testis" ~ "Testis Cancer",
    primary_organ == "Penis" ~ "Penile Cancer",
    primary_organ == "Ampulla of Vater" ~ "Ampullary cancer",
    primary_organ == "Vulva/Vagina" & first_pt %in% c(
      "Vulval Cancer",
      "Squamous Cell Carcinoma Of The Vulva",
      "Vulvar Adenocarcinoma",
      "Vulvar Basal Cell Carcinoma"
    ) ~ "Vulval Cancer",
    primary_organ == "Vulva/Vagina" & first_pt %in% c(
      "Vaginal Cancer",
      "Squamous Cell Carcinoma Of The Vagina",
      "Vaginal Adenocarcinoma"
    ) ~ "Vaginal Cancer",
    primary_organ == "Vulva/Vagina" & first_pt %in% c(
      "Genital Neoplasm Malignant Female"
    ) ~ "Others",
    primary_organ == "Adrenal gland" & first_pt %in% c(
      "Adrenocortical Carcinoma"
    ) ~ "Adrenocortical Carcinoma",
    primary_organ == "Adrenal gland" & first_pt %in% c(
      "Adrenal Gland Cancer"
    ) ~ "Others",
    TRUE ~ histologic_category
  )) %>%
  mutate(assigned_type = case_when(
    histologic_category == "Testis Cancer" ~ "TGCT",
    TRUE ~ assigned_type
  ))

## 5.59 保存最终数据####
save(jader_malignant_final, file = file.path(.output_path, "jader_malignant_final.RData"))
jader_solid <- jader_malignant_final %>%
  filter(!primary_organ %in% c("CUP", "Hematologic", "Unmapped"))
save(jader_solid, file = file.path(.output_path, "jader_solid.RData"))
message(paste("\nJADER恶性肿瘤样本已保存:", file.path(.output_path, "jader_malignant_final.RData")))

# 输出列名供参考
message("\n数据列名:")
print(names(jader_malignant_final))

message("\n========== 两个数据库分析集处理全部完成 ==========")
