# 癌种映射流程与适应症映射表（返修：审稿人2，cancer-level mappings 与字典误分类） ####
# 目的：
#      1) 统计癌种识别流程（适应症去重 → 肿瘤信号 → 病例分类 → 原发器官映射 → 实体瘤 → ICI → 纳入癌种）各步的病例数，供补充图使用；
#      2) 列出FAERS中所有带肿瘤信号的适应症PT（标准MedDRA PT）及其归类（信号类型、TCGA代码、最终原发器官）和报告数，作为补充表；
#         JADER原疾患为日英复合的半自由文本，不逐条列出；
#      3) 按最终原发器官汇总映射到该器官的PT数和ICI报告数。
# 只读取2分析集处理.R、2.1分析集处理FAERS.R等已保存的中间数据，不改变映射规则与分析集。

rm(list = ls())
gc()

# 1. 加载配置 ####
source("Code/config.R")
suppressPackageStartupMessages(library(data.table))


# 2. 结果目录与输入文件 ####
results_dir_name <- "6.2.19癌种映射流程与适应症映射表"
script_results_path <- file.path(results_path, results_dir_name)

if (!dir.exists(script_results_path)) {
    dir.create(script_results_path, recursive = TRUE, showWarnings = FALSE)
}

input_files <- list(
    FAERS = list(
        summary = file.path(output_path, "各主键信号标记.RData"),
        tagged = file.path(output_path, "faers_indi_tagged.RData"),
        tagged_object = "faers_indi_tagged",
        tagged_pt_column = "indi_pt_std",
        malignant = file.path(output_path, "faers_malignant_final.RData"),
        malignant_object = "faers_malignant_final",
        solid = file.path(output_path, "faers_solid.RData"),
        solid_object = "faers_solid",
        solid_pt_column = "indi_pt_list_tumor"
    ),
    JADER = list(
        summary = file.path(output_path, "jader_primaryid_summary.RData"),
        tagged = file.path(output_path, "jader_hist_tagged.RData"),
        tagged_object = "jader_hist_tagged",
        tagged_pt_column = "medical_history_std",
        malignant = file.path(output_path, "jader_malignant_final.RData"),
        malignant_object = "jader_malignant_final",
        solid = file.path(output_path, "jader_solid.RData"),
        solid_object = "jader_solid",
        solid_pt_column = "hist_pt_list_tumor"
    )
)

pattern_assignment_file <- file.path(results_path, "6.1癌种NMF分群分析_rank5", "FAERS_NMF_WH_matrix.xlsx")
output_xlsx_file <- file.path(script_results_path, "癌种映射流程与适应症映射表.xlsx")

min_ici_reports_per_cancer <- 100L

message("========== 癌种映射流程与适应症映射表开始 ==========")


# 3. 辅助函数 ####

## 函数：load_named_rdata_object ####
# 功能：从RData文件中读取指定对象 ####
# 输入：input_path，RData文件路径；object_name，对象名称 ####
# 输出：返回读取出的对象 ####
load_named_rdata_object <- function(input_path, object_name) {
    input_env <- new.env(parent = emptyenv())
    load(input_path, envir = input_env)

    if (!exists(object_name, envir = input_env, inherits = FALSE)) {
        stop("RData中未找到对象: ", object_name)
    }

    input_env[[object_name]]
}

## 函数：split_tumor_pts ####
# 功能：将分析集中以" & "连接的肿瘤适应症PT拆成长表（报告 × PT），PT统一为大写 ####
# 输入：record_ids，报告ID；pt_strings，PT连接字符串 ####
# 输出：含record_id与pt两列的data.table ####
split_tumor_pts <- function(record_ids, pt_strings) {
    long_dt <- data.table(record_id = as.character(record_ids), pt = as.character(pt_strings))
    long_dt <- long_dt[!is.na(pt) & pt != ""]
    long_dt <- long_dt[, .(pt = trimws(unlist(strsplit(pt, " & ", fixed = TRUE)))), by = record_id]
    unique(long_dt[pt != "", .(record_id, pt = toupper(pt))])
}


# 4. 数据读取 ####
final_cancer_types <- readxl::read_xlsx(pattern_assignment_file, sheet = "Pattern group assignment")[[1]]

flow_list <- list()
pt_tag_list <- list()
pt_report_count_list <- list()
pt_ici_count_list <- list()
organ_summary_list <- list()

for (dataset_name in names(input_files)) {
    message("处理数据库: ", dataset_name)
    files <- input_files[[dataset_name]]

    summary_dt <- as.data.table(load_named_rdata_object(files$summary, "primaryid_summary"))
    tagged_dt <- as.data.table(load_named_rdata_object(files$tagged, files$tagged_object))
    setnames(tagged_dt, files$tagged_pt_column, "pt")
    tagged_dt[, `:=`(record_id = as.character(primaryid), pt = toupper(pt))]
    malignant_dt <- as.data.table(load_named_rdata_object(files$malignant, files$malignant_object))
    malignant_dt[, record_id := as.character(primaryid)]
    solid_dt <- as.data.table(load_named_rdata_object(files$solid, files$solid_object))
    solid_dt[, record_id := as.character(primaryid)]


    # 5. 流程各步病例数 ####
    single_dt <- malignant_dt[cohort == "single_malignancy"]
    ici_solid_dt <- solid_dt[has_ici == TRUE]
    ici_counts_by_organ <- ici_solid_dt[, .N, by = primary_organ]
    included_organs <- ici_counts_by_organ[N >= min_ici_reports_per_cancer, primary_organ]

    flow_list[[dataset_name]] <- data.table(
        dataset = dataset_name,
        step = c(
            "Deduplicated reports with indication terms",
            "Non-tumor (no tumor signal) - excluded",
            "Benign tumor only - excluded",
            "Uncertain/borderline tumor only - excluded",
            "Malignant tumor signal",
            "Multiple primary cancer types - excluded",
            "Primary cancer type indeterminate - excluded",
            "Single primary cancer type",
            "  of which mapped to a TCGA cancer type",
            "  of which not mapped to a TCGA cancer type (Others)",
            "Hematologic malignancy or no assignable organ - excluded",
            "Rare site without a defined cancer type - excluded",
            "Solid tumor reports",
            "ICI-exposed solid tumor reports",
            paste0("ICI-exposed reports in the 25 cancer types with >= ", min_ici_reports_per_cancer, " FAERS ICI-exposed reports (analysis set)"),
            paste0("ICI-exposed reports in cancer types with < ", min_ici_reports_per_cancer, " FAERS ICI-exposed reports - not included")
        ),
        n_reports = c(
            nrow(summary_dt),
            summary_dt[cohort == "non_tumor", .N],
            summary_dt[cohort == "exclude_benign", .N],
            summary_dt[cohort == "descriptive_indeterminate", .N],
            summary_dt[cohort %in% c("single_malignancy", "other_malignancy"), .N],
            summary_dt[cohort == "other_malignancy" & n_types_fine >= 2, .N],
            summary_dt[cohort == "other_malignancy" & n_types_fine < 2, .N],
            summary_dt[cohort == "single_malignancy", .N],
            single_dt[!is.na(assigned_type) & !assigned_type %in% c("Others", "Unknown/Others"), .N],
            single_dt[is.na(assigned_type) | assigned_type %in% c("Others", "Unknown/Others"), .N],
            single_dt[primary_organ %in% c("Hematologic", "Unmapped"), .N],
            single_dt[primary_organ == "CUP", .N],
            nrow(solid_dt),
            nrow(ici_solid_dt),
            ici_solid_dt[primary_organ %in% final_cancer_types, .N],
            ici_solid_dt[!primary_organ %in% final_cancer_types, .N]
        )
    )

    # 核查：FAERS中按 ≥100 例确定的癌种应正是最终模型的25个癌种 ####
    if (dataset_name == "FAERS" && !setequal(included_organs, final_cancer_types)) {
        stop("FAERS中ICI报告数 ≥ ", min_ici_reports_per_cancer, " 的癌种与最终模型的25个癌种不一致")
    }

    # 一致性核查：单一原发癌 = 血液/无器官 + 罕见部位 + 实体瘤 ####
    check_total <- single_dt[primary_organ %in% c("Hematologic", "Unmapped"), .N] + single_dt[primary_organ == "CUP", .N] + nrow(solid_dt)
    if (check_total != nrow(single_dt)) {
        stop(dataset_name, "：单一原发癌拆分后的病例数不相加（", check_total, " vs ", nrow(single_dt), "）")
    }


    # 6. PT层面的归类与报告数 ####

    ## 6.1 PT标记（每个PT一套标记） ####
    pt_tag_list[[dataset_name]] <- unique(tagged_dt[hit_malignant | hit_benign | hit_uncertain, .(
        pt, hit_malignant, hit_benign, hit_uncertain, hit_generalized, cancer_type, primary_organ_pt
    )])

    ## 6.2 全库中含该PT的报告数 ####
    pt_report_count_list[[dataset_name]] <- tagged_dt[hit_malignant | hit_benign | hit_uncertain, .(n_reports = uniqueN(record_id)), by = pt][, dataset := dataset_name]

    ## 6.3 分析集（ICI实体瘤，纳入的25个癌种）中含该PT的报告数及其最终原发器官 ####
    ici_pt_dt <- split_tumor_pts(ici_solid_dt$record_id, ici_solid_dt[[files$solid_pt_column]])
    ici_pt_dt <- merge(ici_pt_dt, ici_solid_dt[, .(record_id, primary_organ)], by = "record_id")
    pt_ici_count_list[[dataset_name]] <- ici_pt_dt[, .(n_ici_reports = uniqueN(record_id)), by = .(pt, primary_organ)][, dataset := dataset_name]

    # 罕见部位PT的最终器官（2分析集处理.R中按first_pt指定），用于PT映射表 ####
    rare_map <- unique(single_dt[, .(pt = toupper(first_pt), final_organ = primary_organ)])
    pt_tag_list[[paste0(dataset_name, "_rare")]] <- rare_map


    # 7. 原发器官层面的汇总 ####
    organ_summary_list[[dataset_name]] <- merge(
        ici_counts_by_organ[, .(primary_organ, n_ici_reports = N)],
        ici_pt_dt[, .(n_distinct_pts = uniqueN(pt)), by = primary_organ],
        by = "primary_organ", all.x = TRUE
    )[, `:=`(dataset = dataset_name, included_in_final_model = primary_organ %in% final_cancer_types)]

    rm(summary_dt, tagged_dt, malignant_dt, solid_dt, single_dt, ici_solid_dt, ici_pt_dt)
    gc()
}


# 8. 汇总PT映射表 ####

## 8.1 PT标记：只列FAERS的标准MedDRA PT（JADER原疾患为日英复合的半自由文本，不逐条列出） ####
pt_tags <- unique(pt_tag_list[["FAERS"]])
if (anyDuplicated(pt_tags$pt) > 0) {
    stop("FAERS中同一PT存在多套标记")
}

rare_map_all <- unique(pt_tag_list[["FAERS_rare"]])
rare_map_all <- rare_map_all[!final_organ %in% c("CUP", "Hematologic", "Unmapped")][, .SD[1], by = pt]

pt_table <- copy(pt_tags)
pt_table[, signal := fcase(
    hit_malignant & hit_generalized, "Malignant, primary site not specified",
    hit_malignant, "Malignant",
    hit_benign & !hit_uncertain, "Benign",
    default = "Uncertain/borderline"
)]
pt_table <- merge(pt_table, rare_map_all, by = "pt", all.x = TRUE)
pt_table[, final_primary_organ := fcase(
    !hit_malignant, "Not applicable (excluded)",
    hit_generalized | primary_organ_pt == "CUP", "Primary site not specified (excluded)",
    primary_organ_pt == "Hematologic", "Hematologic malignancy (excluded)",
    is.na(primary_organ_pt), "No assignable organ (excluded)",
    primary_organ_pt == "Rare" & !is.na(final_organ), final_organ,
    primary_organ_pt == "Rare", "Rare site (excluded)",
    default = primary_organ_pt
)]

## 8.2 报告数（FAERS） ####
report_counts <- pt_report_count_list[["FAERS"]][, .(pt, faers_reports = n_reports)]
ici_counts <- pt_ici_count_list[["FAERS"]][, .(faers_ici_reports = sum(n_ici_reports)), by = pt]
pt_table <- merge(pt_table, report_counts, by = "pt", all.x = TRUE)
pt_table <- merge(pt_table, ici_counts, by = "pt", all.x = TRUE)
for (col in c("faers_reports", "faers_ici_reports")) {
    set(pt_table, which(is.na(pt_table[[col]])), col, 0L)
}

# 核查：分析集中PT的最终器官应与映射表一致 ####
ici_organ_check <- pt_ici_count_list[["FAERS"]][, .(pt, primary_organ)]
ici_organ_check <- merge(unique(ici_organ_check), pt_table[, .(pt, final_primary_organ)], by = "pt", all.x = TRUE)
n_inconsistent <- ici_organ_check[primary_organ != final_primary_organ, uniqueN(pt)]
message("分析集中最终器官与PT映射表不一致的PT数（同一报告含多个同器官PT时可能出现）: ", n_inconsistent)

signal_order <- c("Malignant", "Malignant, primary site not specified", "Uncertain/borderline", "Benign")
pt_table[, signal := factor(signal, levels = signal_order)]
setorder(pt_table, signal, final_primary_organ, pt)

pt_mapping_table <- pt_table[, .(
    `Indication term (MedDRA PT)` = tools::toTitleCase(tolower(pt)),
    `Tumor signal` = as.character(signal),
    `TCGA cancer type` = fifelse(is.na(cancer_type), "", cancer_type),
    `Primary organ (cancer type in analysis)` = final_primary_organ,
    `FAERS reports` = faers_reports,
    `ICI-exposed reports in the analysis set` = faers_ici_reports
)]


# 9. 汇总表 ####
flow_df <- dcast(rbindlist(flow_list), step ~ dataset, value.var = "n_reports")
flow_df <- flow_df[match(flow_list$FAERS$step, step)]

organ_df <- dcast(rbindlist(organ_summary_list), primary_organ + included_in_final_model ~ dataset, value.var = c("n_ici_reports", "n_distinct_pts"), fill = 0)
setorder(organ_df, -included_in_final_model, -n_ici_reports_FAERS)

signal_summary <- pt_table[, .(
    n_pts = .N,
    n_pts_in_analysis_set = sum(faers_ici_reports > 0)
), by = .(signal, final_primary_organ_group = fifelse(grepl("excluded", final_primary_organ), final_primary_organ, "Mapped to a solid-tumor organ"))]

message("流程各步病例数:")
print(flow_df)
message("PT映射表行数: ", nrow(pt_mapping_table))
print(signal_summary)


# 10. 结果保存 ####
writexl::write_xlsx(
    list(
        "流程各步病例数" = flow_df,
        "适应症PT映射表" = pt_mapping_table,
        "PT归类汇总" = signal_summary,
        "原发器官汇总" = organ_df,
        "分析参数" = tibble::tibble(
            参数 = c("数据来源", "肿瘤信号", "病例分类", "原发器官", "纳入癌种", "报告数"),
            值 = c(
                "2分析集处理.R / 2.1分析集处理FAERS.R 保存的中间数据（FAERS适应症表；JADER原疾患表，经日英MedDRA对照转为英文PT）",
                "恶性 / 明确良性 / 良恶性不明（交界）；恶性中包含原发部位未指明的泛化表述",
                "非肿瘤、仅良性、仅良恶性不明、多原发癌、原发癌种不可判定均排除；仅单一原发癌进入分析",
                "单一原发癌映射到TCGA代码（不能确定时为Others），再映射到原发器官；血液肿瘤、无可映射器官、未定义罕见部位排除",
                paste0("ICI暴露实体瘤报告数 ≥ ", min_ici_reports_per_cancer, " 的癌种（最终模型25个）"),
                "FAERS全库中含该PT的去重报告数；分析集报告数为纳入25个癌种的ICI暴露实体瘤报告中含该PT的报告数（仅FAERS）"
            )
        )
    ),
    path = output_xlsx_file
)

message("已生成结果文件: ", output_xlsx_file)
message("========== 癌种映射流程与适应症映射表结束 ==========")
