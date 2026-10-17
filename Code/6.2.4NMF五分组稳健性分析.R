# NMF五分组稳健性分析（返修：审稿人2第一段） ####
# 目的：检验rank 5的五分组对以下分析选择是否稳健，并识别归属不确定的癌种
#      1) 随机种子：10个不同种子；
#      2) 矩阵变换：
#         a. 保留癌种内相对谱形的变换：log2ROR癌种内min-max、ROR除以癌种内最大值、癌种内z分数平移为非负；
#         b. 其他表示方式：截尾log2ROR双通道、仅显著单元的log2ROR双通道、癌种内排名；
#      3) 癌种纳入阈值：报告数≥30、≥50、≥200、≥500；逐个剔除一个癌种（留一法）；
#      4) 归属不确定性：按1-3全部重跑中各癌种仍归入原组的比例定义不确定癌种（<80%），
#         剔除后重跑；另按最终模型silhouette剔除最低20%作为补充。
# 每个条件均以rank 5、200次运行进行NMF（PNMF），以一致性聚类得到分组，并与最终五分组比较。
# 每个条件的NMF结果单独缓存，中断后重新运行会跳过已完成的条件。

rm(list = ls())
gc()

# 1. 加载配置与NMF相关R包 ####
source("Code/config.R")
source("Code/mapping_rules.R")

suppressPackageStartupMessages(library(NMF))
suppressPackageStartupMessages(library(pNMF))
suppressPackageStartupMessages(library(cluster))

pnmf_registration_status <- tryCatch(
    {
        NMF::setNMFMethod("PNMF", pNMF::PNMF)
        "registered"
    },
    error = function(e) {
        if (grepl("key already exists", e$message, fixed = TRUE)) {
            "already_registered"
        } else {
            stop(e)
        }
    }
)


# 2. 结果目录与参数 ####
results_dir_name <- "6.2.4NMF五分组稳健性分析"
script_results_path <- file.path(results_path, results_dir_name)
run_cache_path <- file.path(script_results_path, "各条件NMF结果")

for (current_dir in c(script_results_path, run_cache_path)) {
    if (!dir.exists(current_dir)) {
        dir.create(current_dir, recursive = TRUE, showWarnings = FALSE)
    }
}

main_nmf_results_path <- file.path(results_path, "6.1癌种NMF分群分析_rank5")
ror_input_file <- file.path(results_path, "3.1primary_organ不成比例分析", "FAERS_primary_organ_SOC_M0.xlsx")
summary_input_file <- file.path(results_path, "3.1primary_organ不成比例分析", "FAERS_primary_organ_SOC_summary.xlsx")
main_wh_matrix_file <- file.path(main_nmf_results_path, "FAERS_NMF_WH_matrix.xlsx")
main_nmf_result_file <- file.path(main_nmf_results_path, "FAERS_NMF_result.rdata")
main_best_fit_file <- file.path(main_nmf_results_path, "FAERS_NMF_best_fit.rdata")
output_xlsx_file <- file.path(script_results_path, "NMF五分组稳健性分析结果.xlsx")

# 与主分析（6.1.1）一致的参数 ####
nmf_rank <- 5L
nmf_nrun <- 200L
main_seed <- 123456L
main_min_cancer_n <- 100L
log2_ror_cap <- 2

# 稳健性分析参数 ####
alternative_seeds <- main_seed + 1000L * (1:10)
alternative_min_cancer_n <- c(30L, 50L, 200L, 500L)
stability_threshold <- 0.80
silhouette_exclusion_prop <- 0.20

reference_group_levels <- c("Ocular", "Hematologic", "Neurologic", "Systemic-like", "Visceral-like")

message("========== NMF五分组稳健性分析开始 ==========")
message("结果目录: ", script_results_path)


# 3. 辅助函数 ####

## 函数：build_params_dataframe ####
# 功能：将参数列表转换为可写入Excel的数据框####
# 输入：params_list，命名列表####
# 输出：包含参数和值两列的数据框####
build_params_dataframe <- function(params_list) {
    data.frame(
        参数 = names(params_list),
        值 = vapply(params_list, function(param_value) {
            paste(as.character(param_value), collapse = ", ")
        }, character(1)),
        stringsAsFactors = FALSE
    )
}

## 函数：calc_adjusted_rand_index ####
# 功能：计算两个分组向量之间的调整兰德指数（ARI）####
# 输入：labels_a、labels_b，等长的分组向量####
# 输出：数值型ARI####
calc_adjusted_rand_index <- function(labels_a, labels_b) {
    contingency <- table(labels_a, labels_b)
    n_total <- sum(contingency)
    sum_cells <- sum(choose(contingency, 2))
    sum_rows <- sum(choose(rowSums(contingency), 2))
    sum_cols <- sum(choose(colSums(contingency), 2))
    expected_index <- sum_rows * sum_cols / choose(n_total, 2)
    max_index <- (sum_rows + sum_cols) / 2

    if (isTRUE(all.equal(max_index, expected_index))) {
        return(1)
    }

    (sum_cells - expected_index) / (max_index - expected_index)
}

## 函数：build_minmax_matrix ####
# 功能：主分析输入矩阵：ROR按癌种min-max标准化，转为SOC×癌种####
# 输入：ror_matrix，癌种×SOC的ROR矩阵####
# 输出：SOC×癌种的非负矩阵####
build_minmax_matrix <- function(ror_matrix) {
    apply(t(ror_matrix), 2, function(x) {
        if (isTRUE(all.equal(min(x), max(x)))) {
            rep(0, length(x))
        } else {
            (x - min(x)) / (max(x) - min(x))
        }
    })
}

## 函数：build_log2_channel_matrix ####
# 功能：截尾log2(ROR)拆为正向(up)与负向(down)两个非负通道（与6.1输入矩阵构建规则一致）####
# 输入：ror_matrix，癌种×SOC的ROR矩阵；cap_value，log2(ROR)截尾阈值；keep_matrix，可选，癌种×SOC逻辑矩阵，FALSE的单元置0####
# 输出：特征（SOC_up/SOC_down）×癌种的非负矩阵，已去除全零特征####
build_log2_channel_matrix <- function(ror_matrix, cap_value, keep_matrix = NULL) {
    clipped_log2_ror <- pmax(pmin(log2(ror_matrix), cap_value), -cap_value)

    if (!is.null(keep_matrix)) {
        clipped_log2_ror[!keep_matrix] <- 0
    }

    up_matrix <- pmax(clipped_log2_ror, 0)
    down_matrix <- pmax(-clipped_log2_ror, 0)
    colnames(up_matrix) <- paste0(colnames(ror_matrix), "_up")
    colnames(down_matrix) <- paste0(colnames(ror_matrix), "_down")

    channel_matrix <- t(cbind(up_matrix, down_matrix))
    channel_matrix[rowSums(channel_matrix) > 0, , drop = FALSE]
}

## 函数：build_rank_matrix ####
# 功能：癌种内按ROR大小对SOC排名，并缩放到0-1####
# 输入：ror_matrix，癌种×SOC的ROR矩阵####
# 输出：SOC×癌种的非负矩阵####
build_rank_matrix <- function(ror_matrix) {
    apply(t(ror_matrix), 2, function(x) {
        (rank(x, ties.method = "average") - 1) / (length(x) - 1)
    })
}

## 函数：build_log2_minmax_matrix ####
# 功能：log2(ROR)后按癌种min-max标准化，转为SOC×癌种####
# 输入：ror_matrix，癌种×SOC的ROR矩阵####
# 输出：SOC×癌种的非负矩阵####
build_log2_minmax_matrix <- function(ror_matrix) {
    build_minmax_matrix(log2(ror_matrix))
}

## 函数：build_max_scaled_matrix ####
# 功能：ROR除以该癌种的最大ROR（不减最小值），转为SOC×癌种####
# 输入：ror_matrix，癌种×SOC的ROR矩阵####
# 输出：SOC×癌种的非负矩阵####
build_max_scaled_matrix <- function(ror_matrix) {
    apply(t(ror_matrix), 2, function(x) x / max(x))
}

## 函数：build_zscore_shifted_matrix ####
# 功能：癌种内对ROR做z分数，再减去该癌种最小值使其非负，转为SOC×癌种####
# 输入：ror_matrix，癌种×SOC的ROR矩阵####
# 输出：SOC×癌种的非负矩阵####
build_zscore_shifted_matrix <- function(ror_matrix) {
    apply(t(ror_matrix), 2, function(x) {
        z_value <- (x - mean(x)) / stats::sd(x)
        z_value - min(z_value)
    })
}

## 函数：run_rank5_nmf_cached ####
# 功能：对输入矩阵运行rank 5 PNMF并做一致性聚类；结果缓存到rdata，已存在则直接读取####
# 输入：input_matrix，特征×癌种矩阵；seed_value，随机种子；condition_id，条件编号（用于缓存文件名）####
# 输出：命名整数向量（癌种→一致性聚类簇）####
run_rank5_nmf_cached <- function(input_matrix, seed_value, condition_id) {
    cache_file <- file.path(run_cache_path, paste0(condition_id, ".rdata"))

    if (file.exists(cache_file)) {
        load(cache_file)
        return(consensus_clusters)
    }

    nmf_fit <- nmf(
        input_matrix,
        rank = nmf_rank,
        nrun = nmf_nrun,
        seed = seed_value,
        method = "PNMF"
    )

    consensus_clusters <- predict(nmf_fit, what = "consensus")
    consensus_clusters <- stats::setNames(as.integer(as.character(consensus_clusters)), colnames(input_matrix))
    consensus_matrix <- NMF::consensus(nmf_fit)

    save(consensus_clusters, consensus_matrix, file = cache_file)
    consensus_clusters
}

## 函数：map_clusters_to_reference ####
# 功能：将新分组的每个簇映射为簇内多数癌种的最终分组；并列最多时记为Mixed####
# 输入：clusters，命名向量（癌种→簇）；reference_groups，命名向量（癌种→最终分组）####
# 输出：命名字符向量（癌种→映射后的分组）####
map_clusters_to_reference <- function(clusters, reference_groups) {
    mapped_group <- rep(NA_character_, length(clusters))
    names(mapped_group) <- names(clusters)

    for (current_cluster in unique(clusters)) {
        member_cancers <- names(clusters)[clusters == current_cluster]
        group_counts <- table(reference_groups[member_cancers])
        top_groups <- names(group_counts)[group_counts == max(group_counts)]
        mapped_group[member_cancers] <- if (length(top_groups) == 1) top_groups else "Mixed"
    }

    mapped_group
}


# 4. 数据读取 ####
ror_matrix_all <- readxl::read_xlsx(ror_input_file, sheet = "ROR矩阵") %>%
    dplyr::mutate(cancer_type = as.character(cancer_type)) %>%
    dplyr::mutate(dplyr::across(-cancer_type, ~ suppressWarnings(as.numeric(.x)))) %>%
    tibble::column_to_rownames("cancer_type") %>%
    as.matrix()
ror_matrix_all <- ror_matrix_all[, !colnames(ror_matrix_all) %in% c("Infections_and_infestations", "Infections and Infestations"), drop = FALSE]

signal_matrix_all <- readxl::read_xlsx(ror_input_file, sheet = "信号矩阵") %>%
    dplyr::mutate(cancer_type = as.character(cancer_type)) %>%
    tibble::column_to_rownames("cancer_type") %>%
    as.matrix()
signal_matrix_all <- signal_matrix_all[rownames(ror_matrix_all), colnames(ror_matrix_all), drop = FALSE]
storage.mode(signal_matrix_all) <- "logical"

soc_abbrev <- vapply(colnames(ror_matrix_all), function(x) {
    current_abbrev <- get_soc_abbrev(gsub("_", " ", x))
    if (is.na(current_abbrev) || current_abbrev == "") x else current_abbrev
}, character(1), USE.NAMES = FALSE)
colnames(ror_matrix_all) <- soc_abbrev
colnames(signal_matrix_all) <- soc_abbrev

organ_counts <- readxl::read_xlsx(summary_input_file, sheet = "器官样本量") %>%
    dplyr::transmute(cancer_type = as.character(primary_organ), n_reports = suppressWarnings(as.numeric(n)))

reference_assignment_df <- readxl::read_xlsx(main_wh_matrix_file, sheet = "Pattern group assignment") %>%
    dplyr::transmute(cancer_type = as.character(cancer_type), reference_group = as.character(`Pattern group`))
reference_groups <- stats::setNames(reference_assignment_df$reference_group, reference_assignment_df$cancer_type)

main_cancers <- organ_counts %>%
    dplyr::filter(n_reports >= main_min_cancer_n) %>%
    dplyr::pull(cancer_type)
main_cancers <- intersect(rownames(ror_matrix_all), main_cancers)

if (!setequal(main_cancers, names(reference_groups))) {
    stop("纳入癌种与主分析最终分组的癌种不一致，请检查上游结果")
}

message("ROR矩阵癌种数: ", nrow(ror_matrix_all), "；主分析癌种数: ", length(main_cancers))


# 5. 构建各稳健性条件 ####
condition_list <- list()

## 5.1 随机种子 ####
for (seed_index in seq_along(alternative_seeds)) {
    condition_list[[length(condition_list) + 1]] <- list(
        condition_id = sprintf("seed_%02d", seed_index),
        condition_type = "Random seed",
        condition_label = paste0("Seed ", alternative_seeds[seed_index]),
        cancers = main_cancers,
        matrix_type = "minmax",
        seed = alternative_seeds[seed_index]
    )
}

## 5.2 矩阵变换 ####
transformation_labels <- c(
    log2_minmax = "log2(ROR), within-cancer min-max",
    max_scaled = "ROR / within-cancer maximum",
    zscore_shifted = "Within-cancer z-score, shifted",
    log2_channel = "Capped log2(ROR), two channels",
    log2_channel_sigonly = "Capped log2(ROR), significant cells only",
    rank = "Within-cancer rank"
)

for (current_matrix_type in names(transformation_labels)) {
    condition_list[[length(condition_list) + 1]] <- list(
        condition_id = paste0("transform_", current_matrix_type),
        condition_type = ifelse(
            current_matrix_type %in% c("log2_minmax", "max_scaled", "zscore_shifted"),
            "Profile-preserving transformation",
            "Profile-altering transformation"
        ),
        condition_label = transformation_labels[[current_matrix_type]],
        cancers = main_cancers,
        matrix_type = current_matrix_type,
        seed = main_seed
    )
}

## 5.3 癌种纳入阈值 ####
for (current_min_n in alternative_min_cancer_n) {
    current_cancers <- organ_counts %>%
        dplyr::filter(n_reports >= current_min_n) %>%
        dplyr::pull(cancer_type)

    condition_list[[length(condition_list) + 1]] <- list(
        condition_id = paste0("threshold_", current_min_n),
        condition_type = "Inclusion threshold",
        condition_label = paste0(">= ", current_min_n, " reports"),
        cancers = intersect(rownames(ror_matrix_all), current_cancers),
        matrix_type = "minmax",
        seed = main_seed
    )
}

## 5.4 留一法 ####
for (left_out_cancer in main_cancers) {
    condition_list[[length(condition_list) + 1]] <- list(
        condition_id = paste0("leave_one_out_", gsub("[^A-Za-z0-9]", "_", left_out_cancer)),
        condition_type = "Leave one cancer out",
        condition_label = paste0("Without ", left_out_cancer),
        cancers = setdiff(main_cancers, left_out_cancer),
        matrix_type = "minmax",
        seed = main_seed
    )
}

message("稳健性条件数（1-3）: ", length(condition_list))


# 6. 运行各条件NMF并整理分组 ####
run_summary_list <- list()
run_assignment_list <- list()

for (current_condition in condition_list) {
    message("运行条件: ", current_condition$condition_id, "（癌种数 ", length(current_condition$cancers), "）")

    current_ror <- ror_matrix_all[current_condition$cancers, , drop = FALSE]

    # ROR存在缺失或非正值时无法构建输入矩阵，记录为失败 ####
    if (any(!is.finite(current_ror)) || any(current_ror <= 0)) {
        message("  ROR存在缺失或非正值，跳过")
        run_summary_list[[current_condition$condition_id]] <- tibble::tibble(
            condition_id = current_condition$condition_id,
            condition_type = current_condition$condition_type,
            condition_label = current_condition$condition_label,
            n_cancers = length(current_condition$cancers),
            status = "failed: invalid ROR"
        )
        next
    }

    current_input_matrix <- switch(
        current_condition$matrix_type,
        minmax = build_minmax_matrix(current_ror),
        log2_minmax = build_log2_minmax_matrix(current_ror),
        max_scaled = build_max_scaled_matrix(current_ror),
        zscore_shifted = build_zscore_shifted_matrix(current_ror),
        log2_channel = build_log2_channel_matrix(current_ror, log2_ror_cap),
        log2_channel_sigonly = build_log2_channel_matrix(current_ror, log2_ror_cap, signal_matrix_all[current_condition$cancers, , drop = FALSE]),
        rank = build_rank_matrix(current_ror)
    )

    # 全零癌种（如显著单元为空）无法参与NMF，记录后剔除 ####
    zero_cancers <- colnames(current_input_matrix)[colSums(current_input_matrix) == 0]
    current_input_matrix <- current_input_matrix[, colSums(current_input_matrix) > 0, drop = FALSE]

    current_clusters <- tryCatch(
        run_rank5_nmf_cached(current_input_matrix, current_condition$seed, current_condition$condition_id),
        error = function(e) {
            message("  条件运行失败: ", e$message)
            NULL
        }
    )

    if (is.null(current_clusters)) {
        run_summary_list[[current_condition$condition_id]] <- tibble::tibble(
            condition_id = current_condition$condition_id,
            condition_type = current_condition$condition_type,
            condition_label = current_condition$condition_label,
            n_cancers = ncol(current_input_matrix),
            status = "failed"
        )
        next
    }

    common_cancers <- intersect(names(current_clusters), names(reference_groups))
    mapped_groups <- map_clusters_to_reference(current_clusters[common_cancers], reference_groups[common_cancers])

    run_summary_list[[current_condition$condition_id]] <- tibble::tibble(
        condition_id = current_condition$condition_id,
        condition_type = current_condition$condition_type,
        condition_label = current_condition$condition_label,
        n_cancers = length(current_clusters),
        n_common_cancers = length(common_cancers),
        excluded_zero_cancers = paste(zero_cancers, collapse = ", "),
        added_cancers = paste(setdiff(names(current_clusters), names(reference_groups)), collapse = ", "),
        ari_vs_final = calc_adjusted_rand_index(current_clusters[common_cancers], reference_groups[common_cancers]),
        n_cancers_reassigned = sum(mapped_groups != reference_groups[common_cancers]),
        reassigned_cancers = paste(common_cancers[mapped_groups != reference_groups[common_cancers]], collapse = ", "),
        status = "completed"
    )

    run_assignment_list[[current_condition$condition_id]] <- tibble::tibble(
        condition_id = current_condition$condition_id,
        condition_type = current_condition$condition_type,
        condition_label = current_condition$condition_label,
        cancer_type = names(current_clusters),
        cluster = unname(current_clusters),
        reference_group = unname(reference_groups[names(current_clusters)]),
        mapped_group = unname(mapped_groups[names(current_clusters)]),
        same_as_final = mapped_group == reference_group
    )
}

run_summary_df <- dplyr::bind_rows(run_summary_list)
run_assignment_df <- dplyr::bind_rows(run_assignment_list)


# 7. 各癌种归属稳定性 ####

## 7.1 最终模型的silhouette与主导模式H得分占比 ####
load(main_nmf_result_file)
load(main_best_fit_file)

final_consensus <- NMF::consensus(nmf_result)
final_h <- coef(best_fit)
colnames(final_h) <- colnames(final_consensus)
final_consensus <- final_consensus[names(reference_groups), names(reference_groups)]
final_silhouette <- cluster::silhouette(
    as.integer(factor(reference_groups)),
    stats::as.dist(1 - final_consensus)
)
final_silhouette_values <- stats::setNames(final_silhouette[, "sil_width"], names(reference_groups))

dominant_h_share <- apply(final_h[, names(reference_groups), drop = FALSE], 2, function(x) max(x) / sum(x))

## 7.2 重跑中仍归入原组的比例 ####
cancer_stability_df <- run_assignment_df %>%
    dplyr::filter(!is.na(reference_group)) %>%
    dplyr::group_by(cancer_type, reference_group) %>%
    dplyr::summarise(
        n_runs = dplyr::n(),
        n_same_as_final = sum(same_as_final, na.rm = TRUE),
        prop_same_as_final = n_same_as_final / n_runs,
        conditions_reassigned = paste(condition_label[!same_as_final], collapse = "; "),
        .groups = "drop"
    ) %>%
    dplyr::mutate(
        final_silhouette = unname(final_silhouette_values[cancer_type]),
        dominant_h_share = unname(dominant_h_share[cancer_type]),
        uncertain_assignment = prop_same_as_final < stability_threshold,
        reference_group = factor(reference_group, levels = reference_group_levels)
    ) %>%
    dplyr::arrange(reference_group, prop_same_as_final) %>%
    dplyr::mutate(reference_group = as.character(reference_group))

## 7.3 按检验类别分别统计各癌种仍归入原组的比例 ####
cancer_stability_by_type_df <- run_assignment_df %>%
    dplyr::filter(!is.na(reference_group)) %>%
    dplyr::group_by(cancer_type, reference_group, condition_type) %>%
    dplyr::summarise(
        n_runs = dplyr::n(),
        n_same_as_final = sum(same_as_final, na.rm = TRUE),
        prop_same_as_final = n_same_as_final / n_runs,
        .groups = "drop"
    ) %>%
    dplyr::mutate(reference_group = factor(reference_group, levels = reference_group_levels)) %>%
    dplyr::arrange(reference_group, cancer_type, condition_type) %>%
    dplyr::mutate(reference_group = as.character(reference_group))

uncertain_cancers <- cancer_stability_df %>%
    dplyr::filter(uncertain_assignment) %>%
    dplyr::pull(cancer_type)

low_silhouette_cancers <- names(sort(final_silhouette_values))[seq_len(round(silhouette_exclusion_prop * length(final_silhouette_values)))]

message("归属不确定癌种（<", stability_threshold * 100, "%）: ", ifelse(length(uncertain_cancers) == 0, "无", paste(uncertain_cancers, collapse = ", ")))
message("silhouette最低", silhouette_exclusion_prop * 100, "%的癌种: ", paste(low_silhouette_cancers, collapse = ", "))


# 8. 剔除不确定癌种后重跑 ####
exclusion_conditions <- list(
    list(
        condition_id = "exclude_unstable",
        condition_label = paste0("Excluding cancers stable in <", stability_threshold * 100, "% of runs"),
        excluded = uncertain_cancers
    ),
    list(
        condition_id = "exclude_low_silhouette",
        condition_label = paste0("Excluding lowest ", silhouette_exclusion_prop * 100, "% by silhouette"),
        excluded = low_silhouette_cancers
    )
)

exclusion_summary_list <- list()
exclusion_assignment_list <- list()

for (current_exclusion in exclusion_conditions) {
    if (length(current_exclusion$excluded) == 0) {
        exclusion_summary_list[[current_exclusion$condition_id]] <- tibble::tibble(
            condition_id = current_exclusion$condition_id,
            condition_label = current_exclusion$condition_label,
            excluded_cancers = "",
            n_cancers = length(main_cancers),
            status = "not run (no cancer met the exclusion criterion)"
        )
        next
    }

    retained_cancers <- setdiff(main_cancers, current_exclusion$excluded)
    current_input_matrix <- build_minmax_matrix(ror_matrix_all[retained_cancers, , drop = FALSE])
    current_clusters <- run_rank5_nmf_cached(current_input_matrix, main_seed, current_exclusion$condition_id)
    mapped_groups <- map_clusters_to_reference(current_clusters, reference_groups[names(current_clusters)])

    retained_group_sizes <- table(factor(reference_groups[retained_cancers], levels = reference_group_levels))

    exclusion_summary_list[[current_exclusion$condition_id]] <- tibble::tibble(
        condition_id = current_exclusion$condition_id,
        condition_label = current_exclusion$condition_label,
        excluded_cancers = paste(current_exclusion$excluded, collapse = ", "),
        n_cancers = length(retained_cancers),
        retained_group_sizes = paste(paste0(names(retained_group_sizes), "=", retained_group_sizes), collapse = "; "),
        ari_vs_final = calc_adjusted_rand_index(current_clusters, reference_groups[names(current_clusters)]),
        n_cancers_reassigned = sum(mapped_groups != reference_groups[names(current_clusters)]),
        reassigned_cancers = paste(names(current_clusters)[mapped_groups != reference_groups[names(current_clusters)]], collapse = ", "),
        status = "completed"
    )

    exclusion_assignment_list[[current_exclusion$condition_id]] <- tibble::tibble(
        condition_id = current_exclusion$condition_id,
        condition_label = current_exclusion$condition_label,
        cancer_type = names(current_clusters),
        cluster = unname(current_clusters),
        reference_group = unname(reference_groups[names(current_clusters)]),
        mapped_group = unname(mapped_groups[names(current_clusters)]),
        same_as_final = mapped_group == reference_group
    )
}

exclusion_summary_df <- dplyr::bind_rows(exclusion_summary_list)
exclusion_assignment_df <- dplyr::bind_rows(exclusion_assignment_list)

message("各类条件ARI汇总:")
print(run_summary_df %>%
    dplyr::group_by(condition_type) %>%
    dplyr::summarise(n = dplyr::n(), min_ari = min(ari_vs_final, na.rm = TRUE), median_ari = stats::median(ari_vs_final, na.rm = TRUE), .groups = "drop"))
print(as.data.frame(exclusion_summary_df))


# 9. 结果保存 ####
writexl::write_xlsx(
    list(
        "稳健性运行汇总" = run_summary_df,
        "各运行癌种分组" = run_assignment_df,
        "癌种归属稳定性" = cancer_stability_df,
        "癌种归属稳定性_按检验类别" = cancer_stability_by_type_df,
        "剔除不确定癌种汇总" = exclusion_summary_df,
        "剔除不确定癌种分组" = exclusion_assignment_df,
        "分析参数" = build_params_dataframe(list(
            nmf_rank = nmf_rank,
            nmf_nrun = nmf_nrun,
            method = "PNMF",
            main_seed = main_seed,
            alternative_seeds = alternative_seeds,
            main_min_cancer_n = main_min_cancer_n,
            alternative_min_cancer_n = alternative_min_cancer_n,
            log2_ror_cap = log2_ror_cap,
            matrix_transformations = "minmax（主分析）；Profile-preserving（保留癌种内相对谱形）：log2_minmax、max_scaled、zscore_shifted；Profile-altering（改变相对谱形）：log2_channel、log2_channel_sigonly、rank",
            cluster_assignment = "一致性聚类（predict(what = 'consensus')）",
            mapping_to_final_groups = "新分组各簇映射为簇内多数癌种的最终分组，并列最多记为Mixed",
            stability_threshold = stability_threshold,
            silhouette_exclusion_prop = silhouette_exclusion_prop,
            ror_source = ror_input_file,
            reference_groups = main_wh_matrix_file
        ))
    ),
    path = output_xlsx_file
)

message("已生成结果文件: ", output_xlsx_file)
message("========== NMF五分组稳健性分析结束 ==========")
