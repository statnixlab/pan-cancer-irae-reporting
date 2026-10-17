# 通路富集的置换校准分析（返修：审稿人2第三段，基因与通路间的依赖） ####
# 目的：ORA（超几何检验）假设基因相互独立；共表达基因会成批进入候选名单，使富集P值偏小。
#      本脚本在7.7已保存的1,000次器官标签置换上，原样重复“筛选候选基因 → ORA”流程，
#      得到保留基因间相关结构的零分布，对Figure 7展示的全部46条通路计算校准后的经验P值。
# 预先规定：检验集合为Figure 7（7.7.2「通路富集长表」）的全部46条通路；校准P在46条内做BH，BH < 0.05视为校准后仍显著。
# 置换k下各基因的P值用其余999次置换作参照（留一法），候选基因规则与7.7一致（|Δρ| ≥ 0.3 且 P < 0.05，按Δρ方向分正负）。
# ORA口径与7.7.1一致：背景为7.7预过滤保留基因与该数据库基因的交集，通路基因取其与背景的交集。
# 只读取7.7、7.7.1、7.7.2已保存的结果，不修改原结果。

rm(list = ls())
gc()

# 1. 加载配置 ####
source("Code/config.R")


# 2. 结果目录与参数 ####
results_dir_name <- "6.2.16通路富集置换校准分析"
script_results_path <- file.path(results_path, results_dir_name)

if (!dir.exists(script_results_path)) {
    dir.create(script_results_path, recursive = TRUE, showWarnings = FALSE)
}

input_77_rdata_path <- file.path(results_path, "7.7模式相关基因差分相关识别", "rdata")
null_array_file <- file.path(input_77_rdata_path, "delta_permutation_null_array.Rdata")
final_gene_list_file <- file.path(input_77_rdata_path, "final_pattern_gene_list.Rdata")
gene_filter_file <- file.path(input_77_rdata_path, "gene_filter_summary_df.Rdata")
term2gene_file <- file.path(results_path, "7.7.1模式相关基因富集分析", "rdata", "enrichment_database_term2gene_list.Rdata")
figure_pathway_file <- file.path(results_path, "7.7.2模式通路富集图数据准备", "tables", "模式通路富集图数据.xlsx")

output_xlsx_file <- file.path(script_results_path, "通路富集置换校准结果.xlsx")

# 与7.7一致的候选基因规则 ####
delta_rho_threshold <- 0.3
perm_p_threshold <- 0.05
calibrated_fdr_threshold <- 0.05

pattern_levels <- c("Ocular", "Hematologic", "Neurologic", "Systemic-like", "Visceral-like")

message("========== 通路富集置换校准分析开始 ==========")
message("结果目录: ", script_results_path)


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

## 函数：calc_ora_p_value ####
# 功能：超几何检验的上尾P值（与clusterProfiler::enricher一致：P(X ≥ 重叠数)） ####
# 输入：n_overlap，候选基因与通路的重叠数；n_term，通路在背景中的基因数；n_universe，背景基因数；n_hit，候选基因在背景中的个数（均可为向量） ####
# 输出：数值向量P值 ####
calc_ora_p_value <- function(n_overlap, n_term, n_universe, n_hit) {
    stats::phyper(n_overlap - 1, n_term, n_universe - n_term, n_hit, lower.tail = FALSE)
}

## 函数：calc_leave_one_out_p ####
# 功能：每次置换下各基因的经验P值：以其余置换的|Δρ|作参照，P = (≥ 当前值的其余置换数 + 1) / 置换总数 ####
# 输入：abs_null_matrix，置换 × 基因的|Δρ|矩阵（可含NA） ####
# 输出：同维度的P值矩阵（NA位置保持NA） ####
calc_leave_one_out_p <- function(abs_null_matrix) {
    n_perm <- nrow(abs_null_matrix)
    p_matrix <- apply(abs_null_matrix, 2, function(x) {
        if (all(is.na(x))) {
            return(rep(NA_real_, length(x)))
        }
        n_greater_or_equal <- rank(-x, ties.method = "max", na.last = "keep")
        n_greater_or_equal / n_perm
    })
    p_matrix
}


# 4. 数据读取 ####
figure_pathway_df <- readxl::read_xlsx(figure_pathway_file, sheet = "通路富集长表") %>%
    dplyr::select(pattern, hit_direction, database, term_id, n_hit_in_term, n_term_in_universe, n_universe, n_hit_total_in_universe, fold_enrichment, enrichment_pvalue, enrichment_padj_BH)

term2gene_list <- load_named_rdata_object(term2gene_file, "database_resource_list")

background_genes <- load_named_rdata_object(gene_filter_file, "gene_filter_summary_df") %>%
    dplyr::filter(retained) %>%
    dplyr::pull(gene) %>%
    unique() %>%
    as.character()

hit_long_df <- load_named_rdata_object(final_gene_list_file, "hit_long_df") %>%
    dplyr::filter(delta_rho_threshold == !!delta_rho_threshold)

null_tumor_array <- load_named_rdata_object(null_array_file, "permutation_tumor_matched_array")
null_normal_array <- load_named_rdata_object(null_array_file, "permutation_normal_matched_array")
null_gene_names <- dimnames(null_tumor_array)[[2]]
n_permutations <- dim(null_tumor_array)[1]

if (!identical(dimnames(null_tumor_array), dimnames(null_normal_array))) {
    stop("肿瘤端与正常端置换数组的维度名不一致")
}

message("Figure 7通路数: ", nrow(figure_pathway_df), "；置换次数: ", n_permutations, "；置换基因数: ", length(null_gene_names))


# 5. 构建各数据库的背景与通路基因集（与7.7.1口径一致） ####
database_universe_list <- lapply(names(term2gene_list), function(db_name) {
    intersect(background_genes, unique(term2gene_list[[db_name]]$term2gene_df[[2]]))
})
names(database_universe_list) <- names(term2gene_list)

term_gene_list <- lapply(seq_len(nrow(figure_pathway_df)), function(i) {
    current_db <- figure_pathway_df$database[i]
    current_term2gene <- term2gene_list[[current_db]]$term2gene_df
    term_genes <- unique(current_term2gene[[2]][current_term2gene[[1]] == figure_pathway_df$term_id[i]])
    intersect(term_genes, database_universe_list[[current_db]])
})


# 6. 方法核查：用观察到的候选基因重算46条通路的ORA，须与7.7.1保存值一致 ####
reproduction_df <- figure_pathway_df %>%
    dplyr::mutate(
        recalc_n_universe = vapply(database, function(db) length(database_universe_list[[db]]), integer(1)),
        recalc_n_term = vapply(term_gene_list, length, integer(1)),
        recalc_n_hit = NA_integer_,
        recalc_n_overlap = NA_integer_
    )

for (i in seq_len(nrow(reproduction_df))) {
    observed_hits <- hit_long_df$gene[hit_long_df$pattern == reproduction_df$pattern[i] & hit_long_df$hit_direction == reproduction_df$hit_direction[i]]
    observed_hits_in_universe <- intersect(observed_hits, database_universe_list[[reproduction_df$database[i]]])
    reproduction_df$recalc_n_hit[i] <- length(observed_hits_in_universe)
    reproduction_df$recalc_n_overlap[i] <- length(intersect(observed_hits_in_universe, term_gene_list[[i]]))
}

reproduction_df <- reproduction_df %>%
    dplyr::mutate(
        recalc_p_value = calc_ora_p_value(recalc_n_overlap, recalc_n_term, recalc_n_universe, recalc_n_hit),
        counts_identical = recalc_n_universe == n_universe & recalc_n_term == n_term_in_universe &
            recalc_n_hit == n_hit_total_in_universe & recalc_n_overlap == n_hit_in_term,
        relative_p_difference = abs(recalc_p_value - enrichment_pvalue) / enrichment_pvalue
    )

message("计数完全一致的通路: ", sum(reproduction_df$counts_identical), "/", nrow(reproduction_df),
        "；P值最大相对差: ", signif(max(reproduction_df$relative_p_difference), 3))

if (!all(reproduction_df$counts_identical) || max(reproduction_df$relative_p_difference) > 1e-6) {
    stop("未能重现7.7.1保存的ORA结果，停止分析")
}


# 7. 置换零分布：每次置换下重复候选基因筛选与ORA ####
null_p_list <- vector("list", nrow(figure_pathway_df))
null_hit_count_list <- list()

for (current_pattern in pattern_levels) {
    message("处理模式: ", current_pattern)

    ## 7.1 置换下的Δρ与留一法P值 ####
    null_delta <- null_tumor_array[, , current_pattern] - null_normal_array[, , current_pattern]
    null_p <- calc_leave_one_out_p(abs(null_delta))

    null_hit_positive <- !is.na(null_delta) & null_delta >= delta_rho_threshold & null_p < perm_p_threshold
    null_hit_negative <- !is.na(null_delta) & null_delta <= -delta_rho_threshold & null_p < perm_p_threshold

    ## 7.2 置换下的候选基因数（与观察值比较） ####
    for (current_direction in c("positive", "negative")) {
        current_hits <- if (current_direction == "positive") null_hit_positive else null_hit_negative
        n_null_hits <- rowSums(current_hits)
        n_observed_hits <- sum(hit_long_df$pattern == current_pattern & hit_long_df$hit_direction == current_direction)

        null_hit_count_list[[paste(current_pattern, current_direction)]] <- tibble::tibble(
            pattern = current_pattern,
            hit_direction = current_direction,
            n_observed_candidates = n_observed_hits,
            null_median_candidates = stats::median(n_null_hits),
            null_q1_candidates = stats::quantile(n_null_hits, 0.25, names = FALSE),
            null_q3_candidates = stats::quantile(n_null_hits, 0.75, names = FALSE),
            null_p95_candidates = stats::quantile(n_null_hits, 0.95, names = FALSE),
            empirical_p_candidate_count = (sum(n_null_hits >= n_observed_hits) + 1) / (n_permutations + 1)
        )
    }

    ## 7.3 该模式下各通路的置换ORA P值 ####
    for (i in which(figure_pathway_df$pattern == current_pattern)) {
        current_hits <- if (figure_pathway_df$hit_direction[i] == "positive") null_hit_positive else null_hit_negative
        universe_index <- match(database_universe_list[[figure_pathway_df$database[i]]], null_gene_names)
        universe_index <- universe_index[!is.na(universe_index)]
        term_index <- match(term_gene_list[[i]], null_gene_names)
        term_index <- term_index[!is.na(term_index)]

        n_hit_in_universe <- rowSums(current_hits[, universe_index, drop = FALSE])
        n_overlap <- rowSums(current_hits[, term_index, drop = FALSE])

        null_p_list[[i]] <- ifelse(
            n_hit_in_universe > 0,
            calc_ora_p_value(n_overlap, length(term_index), length(universe_index), n_hit_in_universe),
            1
        )
    }

    rm(null_delta, null_p, null_hit_positive, null_hit_negative)
    gc()
}


# 8. 统计分析：校准后的经验P值与BH ####
calibration_df <- reproduction_df %>%
    dplyr::select(pattern, hit_direction, database, term_id, n_hit_in_term, n_term_in_universe, fold_enrichment, enrichment_pvalue, enrichment_padj_BH) %>%
    dplyr::mutate(
        null_median_ora_p = vapply(null_p_list, stats::median, numeric(1)),
        n_permutations_as_extreme = vapply(seq_along(null_p_list), function(i) sum(null_p_list[[i]] <= enrichment_pvalue[i]), integer(1)),
        calibrated_p = (n_permutations_as_extreme + 1) / (n_permutations + 1),
        calibrated_p_BH = stats::p.adjust(calibrated_p, method = "BH"),
        calibrated_significant = calibrated_p_BH < calibrated_fdr_threshold,
        pattern = factor(pattern, levels = pattern_levels)
    ) %>%
    dplyr::arrange(pattern, dplyr::desc(hit_direction), calibrated_p) %>%
    dplyr::mutate(pattern = as.character(pattern))

null_hit_count_df <- dplyr::bind_rows(null_hit_count_list)

message("校准后BH < 0.05的通路: ", sum(calibration_df$calibrated_significant), "/", nrow(calibration_df))
print(as.data.frame(calibration_df %>% dplyr::select(pattern, hit_direction, term_id, enrichment_padj_BH, calibrated_p, calibrated_p_BH)), digits = 3)
print(as.data.frame(null_hit_count_df), digits = 3)


# 9. 结果保存 ####
writexl::write_xlsx(
    list(
        "通路校准结果" = calibration_df,
        "置换下候选基因数" = null_hit_count_df,
        "ORA重现核查" = reproduction_df,
        "分析参数" = tibble::tibble(
            参数 = c("检验集合", "置换", "置换下基因P值", "候选基因规则", "ORA", "校准P值", "多重校正", "判读规则"),
            值 = c(
                "Figure 7（7.7.2「通路富集长表」）展示的全部通路，预先规定，不按结果增减",
                paste0("7.7已保存的", n_permutations, "次器官标签置换（肿瘤端与正常端共用置换，保留基因间相关结构）"),
                "留一法：以其余置换的|Δρ|为参照，P = (≥当前值的其余置换数 + 1) / 置换总数",
                paste0("|Δρ| ≥ ", delta_rho_threshold, " 且 P < ", perm_p_threshold, "，按Δρ方向分正负（与7.7一致）"),
                "超几何上尾检验；背景为7.7预过滤保留基因与该数据库基因的交集（与7.7.1 clusterProfiler::enricher一致）",
                "校准P = (置换ORA P ≤ 观察ORA P的置换数 + 1) / (置换次数 + 1)",
                "在全部Figure 7通路内做BH",
                "校准后BH < 0.05视为在考虑基因间相关后仍显著"
            )
        )
    ),
    path = output_xlsx_file
)

message("已生成结果文件: ", output_xlsx_file)
message("========== 通路富集置换校准分析结束 ==========")
