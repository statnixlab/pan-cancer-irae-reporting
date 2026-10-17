# NMF秩选择敏感性分析（返修：审稿人1第1条、审稿人2秩选择） ####
# 目的：在主分析rank survey（rank 2-8，每个rank 200次运行，seed 123456）基础上，
#      1) 与20个随机置换矩阵的rank survey（零分布）对比，评估各指标是否超出随机结构；
#      2) 提取rank 2-8的共识矩阵、W矩阵、癌种分组和轮廓系数（审稿人点名rank 4-7，为完整起见扩展到主分析的全部候选秩）；
#      3) 用客观指标（最小群大小、单癌种群、轮廓系数、W模式冗余度）比较rank 2-8；
#      4) 量化rank 4/6/7与最终rank 5分组之间的合并/拆分关系（列联表与ARI）。
# 说明：主分析的rank survey结果直接读取，不重新运行；仅随机置换零分布需要新计算。

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
results_dir_name <- "6.2.1NMF秩选择敏感性分析"
script_results_path <- file.path(results_path, results_dir_name)

if (!dir.exists(script_results_path)) {
    dir.create(script_results_path, recursive = TRUE, showWarnings = FALSE)
}

main_nmf_results_path <- file.path(results_path, "6.1癌种NMF分群分析_rank5")

ror_input_file <- file.path(
    results_path,
    "3.1primary_organ不成比例分析",
    "FAERS_primary_organ_SOC_M0.xlsx"
)
summary_input_file <- file.path(
    results_path,
    "3.1primary_organ不成比例分析",
    "FAERS_primary_organ_SOC_summary.xlsx"
)
main_rank_estimate_file <- file.path(main_nmf_results_path, "FAERS_NMF_rank_estimate.rdata")
main_wh_matrix_file <- file.path(main_nmf_results_path, "FAERS_NMF_WH_matrix.xlsx")

permutation_results_path <- file.path(script_results_path, "随机置换基线")
output_xlsx_file <- file.path(script_results_path, "NMF秩选择敏感性分析结果.xlsx")

# 以下参数与主分析（6.1.1）保持一致 ####
rank_range <- 2:8
rank_nrun <- 200L
analysis_seed <- 123456L
min_cancer_n <- 100L
hclust_method <- "average"

# 随机置换基线参数：20个独立置换矩阵，每个rank 100次运行，第i个置换使用seed = analysis_seed + i ####
n_permutations <- 20L
permutation_nrun <- 100L

# 敏感性分析的比较对象 ####
compare_ranks <- 2:8
reference_rank <- 5L
dominant_w_threshold <- 0.20
low_silhouette_threshold <- 0.50

message("========== NMF秩选择敏感性分析开始 ==========")
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

## 函数：to_soc_abbrev ####
# 功能：将ROR矩阵中的SOC列名转换为SOC缩写####
# 输入：soc_names，字符型SOC列名向量####
# 输出：字符型SOC缩写向量####
to_soc_abbrev <- function(soc_names) {
    vapply(soc_names, function(x) {
        current_abbrev <- get_soc_abbrev(gsub("_", " ", x))

        if (is.na(current_abbrev) || current_abbrev == "") {
            x
        } else {
            current_abbrev
        }
    }, character(1), USE.NAMES = FALSE)
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

## 函数：normalize_h_by_cancer ####
# 功能：对H矩阵按癌种（列）做min-max标准化，并转置为癌种×模式矩阵####
# 输入：h_matrix，模式×癌种的H矩阵####
# 输出：癌种×模式的标准化H矩阵####
normalize_h_by_cancer <- function(h_matrix) {
    h_normalized <- apply(h_matrix, 2, function(x) {
        if (isTRUE(all.equal(min(x), max(x)))) {
            rep(0, length(x))
        } else {
            (x - min(x)) / (max(x) - min(x))
        }
    })

    h_normalized <- t(h_normalized)
    colnames(h_normalized) <- rownames(h_matrix)
    h_normalized
}

## 函数：map_clusters_to_patterns ####
# 功能：按主分析（6.1.1）相同规则，将一致性聚类簇一对一映射到W模式####
# 输入：clusters，命名整数向量（癌种→簇）；h_normalized，癌种×模式标准化H矩阵####
# 输出：命名字符向量（簇编号→模式编号）####
map_clusters_to_patterns <- function(clusters, h_normalized) {
    cluster_ids <- sort(unique(as.integer(clusters)))
    pattern_ids <- colnames(h_normalized)

    cluster_pattern_score <- t(vapply(cluster_ids, function(current_cluster_id) {
        colMeans(h_normalized[clusters == current_cluster_id, , drop = FALSE])
    }, numeric(length(pattern_ids))))
    colnames(cluster_pattern_score) <- pattern_ids

    cluster_to_pattern <- rep(NA_character_, length(cluster_ids))
    names(cluster_to_pattern) <- as.character(cluster_ids)
    available_patterns <- pattern_ids

    for (i in order(apply(cluster_pattern_score, 1, max), decreasing = TRUE)) {
        for (j in order(cluster_pattern_score[i, ], decreasing = TRUE)) {
            if (pattern_ids[j] %in% available_patterns) {
                cluster_to_pattern[as.character(cluster_ids[i])] <- pattern_ids[j]
                available_patterns <- setdiff(available_patterns, pattern_ids[j])
                break
            }
        }
    }

    remaining_clusters <- names(cluster_to_pattern)[is.na(cluster_to_pattern)]
    cluster_to_pattern[remaining_clusters] <- available_patterns[seq_along(remaining_clusters)]

    cluster_to_pattern
}

## 函数：label_patterns_by_w ####
# 功能：为每个模式生成主导SOC组合标签，规则与6.1.6主导SOC复合定义一致：####
#       每个SOC归入其W权重最大的模式，且该权重≥阈值；某模式未分到任何SOC时，以其权重最大的SOC作标签####
# 输入：w_matrix，SOC×模式的W矩阵；threshold，主导权重阈值####
# 输出：命名字符向量（模式编号→主导SOC组合，如"GI/HEPA/RENU/RESP"）####
label_patterns_by_w <- function(w_matrix, threshold) {
    soc_best_pattern <- colnames(w_matrix)[apply(w_matrix, 1, which.max)]
    soc_best_weight <- apply(w_matrix, 1, max)

    vapply(colnames(w_matrix), function(current_pattern) {
        dominant_socs <- rownames(w_matrix)[soc_best_pattern == current_pattern & soc_best_weight >= threshold]

        if (length(dominant_socs) == 0) {
            dominant_socs <- rownames(w_matrix)[which.max(w_matrix[, current_pattern])]
        }

        paste(sort(dominant_socs), collapse = "/")
    }, character(1))
}

## 函数：calc_max_w_cosine ####
# 功能：计算W矩阵任意两个模式之间余弦相似度的最大值，用于衡量模式冗余度####
# 输入：w_matrix，SOC×模式的W矩阵####
# 输出：数值型最大余弦相似度####
calc_max_w_cosine <- function(w_matrix) {
    w_unit <- sweep(w_matrix, 2, sqrt(colSums(w_matrix^2)), "/")
    cosine_matrix <- crossprod(w_unit)
    max(cosine_matrix[upper.tri(cosine_matrix)])
}


# 4. 读取数据并构建与主分析一致的NMF输入矩阵 ####
raw_ror_matrix_df <- readxl::read_xlsx(ror_input_file, sheet = "ROR矩阵")
organ_counts_df <- readxl::read_xlsx(summary_input_file, sheet = "器官样本量") %>%
    dplyr::mutate(
        primary_organ = as.character(primary_organ),
        n = suppressWarnings(as.numeric(n))
    )

qualified_cancers <- organ_counts_df %>%
    dplyr::filter(!is.na(n), n >= min_cancer_n) %>%
    dplyr::pull(primary_organ) %>%
    unique()

raw_ror_matrix <- raw_ror_matrix_df %>%
    dplyr::mutate(cancer_type = as.character(cancer_type)) %>%
    dplyr::filter(cancer_type %in% qualified_cancers) %>%
    dplyr::mutate(dplyr::across(-cancer_type, ~ suppressWarnings(as.numeric(.x)))) %>%
    tibble::column_to_rownames("cancer_type") %>%
    as.matrix()

mat <- t(raw_ror_matrix)
mat <- mat[!rownames(mat) %in% c("Infections_and_infestations", "Infections and Infestations"), , drop = FALSE]

mat_row_names <- rownames(mat)
mat_col_names <- colnames(mat)

mat <- apply(mat, 2, function(x) {
    if (isTRUE(all.equal(min(x), max(x)))) {
        rep(0, length(x))
    } else {
        (x - min(x)) / (max(x) - min(x))
    }
})
rownames(mat) <- to_soc_abbrev(mat_row_names)
colnames(mat) <- mat_col_names

message("NMF输入矩阵维度（SOC × 癌种）: ", nrow(mat), " × ", ncol(mat))


# 5. 读取主分析rank survey与最终分组 ####
load(main_rank_estimate_file)
observed_estimate <- estim.r
rm(estim.r)

reference_assignment_df <- readxl::read_xlsx(main_wh_matrix_file, sheet = "Pattern group assignment") %>%
    dplyr::transmute(cancer_type = as.character(cancer_type), reference_group = as.character(`Pattern group`))

reference_groups <- stats::setNames(reference_assignment_df$reference_group, reference_assignment_df$cancer_type)

if (!setequal(names(reference_groups), colnames(mat))) {
    stop("主分析分组中的癌种与当前NMF输入矩阵癌种不一致，请检查上游结果")
}

reference_group_levels <- c("Ocular", "Hematologic", "Neurologic", "Systemic-like", "Visceral-like")


# 6. 随机置换矩阵rank survey（20次置换的零分布） ####
# 每个癌种（列）内部独立置换SOC取值，保留各癌种取值分布、破坏SOC间协同结构 ####
# 每个置换结果单独保存，中断后重新运行会跳过已完成的置换 ####
if (!dir.exists(permutation_results_path)) {
    dir.create(permutation_results_path, recursive = TRUE, showWarnings = FALSE)
}

metric_columns <- c(
    "cophenetic", "dispersion", "silhouette.consensus",
    "silhouette.coef", "silhouette.basis", "rss", "evar",
    "sparseness.basis", "sparseness.coef"
)

permutation_metrics_list <- list()

for (permutation_id in seq_len(n_permutations)) {
    permutation_file <- file.path(
        permutation_results_path,
        sprintf("随机置换%02d_NMF秩估计.rdata", permutation_id)
    )

    if (file.exists(permutation_file)) {
        load(permutation_file)
        message("读取已有随机置换结果: ", permutation_file)
    } else {
        message("运行随机置换 ", permutation_id, "/", n_permutations)
        set.seed(analysis_seed + permutation_id)
        mat_randomized <- NMF::randomize(mat)

        random_estimate <- nmf(
            mat_randomized,
            rank_range,
            nrun = permutation_nrun,
            seed = analysis_seed + permutation_id,
            method = "PNMF"
        )

        save(mat_randomized, random_estimate, file = permutation_file)
    }

    permutation_metrics_list[[permutation_id]] <- as.data.frame(summary(random_estimate)) %>%
        dplyr::mutate(data_type = "Randomized", permutation = permutation_id)
}


# 7. 秩选择指标汇总 ####

## 7.1 实际数据与各置换的逐rank指标（长表） ####
rank_metrics_df <- dplyr::bind_rows(
    as.data.frame(summary(observed_estimate)) %>%
        dplyr::mutate(data_type = "Observed", permutation = 0L),
    dplyr::bind_rows(permutation_metrics_list)
) %>%
    dplyr::select(data_type, permutation, rank, dplyr::all_of(metric_columns)) %>%
    dplyr::arrange(data_type, permutation, rank)

## 7.2 实际值与零分布比较 ####
# 经验P值：cophenetic、dispersion、silhouette取零分布中≥实际值的比例；rss取≤实际值的比例；均按(b+1)/(n+1)计算 ####
observed_metrics_long <- rank_metrics_df %>%
    dplyr::filter(data_type == "Observed") %>%
    dplyr::select(rank, dplyr::all_of(metric_columns)) %>%
    tidyr::pivot_longer(-rank, names_to = "metric", values_to = "observed")

null_metrics_long <- rank_metrics_df %>%
    dplyr::filter(data_type == "Randomized") %>%
    dplyr::select(permutation, rank, dplyr::all_of(metric_columns)) %>%
    tidyr::pivot_longer(-c(permutation, rank), names_to = "metric", values_to = "null_value")

null_comparison_df <- null_metrics_long %>%
    dplyr::left_join(observed_metrics_long, by = c("rank", "metric")) %>%
    dplyr::group_by(metric, rank) %>%
    dplyr::summarise(
        observed = dplyr::first(observed),
        null_mean = mean(null_value),
        null_sd = stats::sd(null_value),
        null_min = min(null_value),
        null_max = max(null_value),
        n_permutations = dplyr::n(),
        n_null_as_extreme = if (dplyr::first(metric) == "rss") {
            sum(null_value <= dplyr::first(observed))
        } else {
            sum(null_value >= dplyr::first(observed))
        },
        .groups = "drop"
    ) %>%
    dplyr::mutate(
        z_score = (observed - null_mean) / null_sd,
        empirical_p = (n_null_as_extreme + 1) / (n_permutations + 1)
    ) %>%
    dplyr::filter(metric %in% c("cophenetic", "dispersion", "silhouette.consensus", "rss", "evar")) %>%
    dplyr::arrange(metric, rank)

## 7.3 每增加一个rank的RSS下降量：实际值与零分布比较 ####
rss_drop_long <- rank_metrics_df %>%
    dplyr::group_by(data_type, permutation) %>%
    dplyr::arrange(rank, .by_group = TRUE) %>%
    dplyr::mutate(rss_drop_from_previous_rank = dplyr::lag(rss) - rss) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!is.na(rss_drop_from_previous_rank))

rss_drop_compare_df <- rss_drop_long %>%
    dplyr::filter(data_type == "Randomized") %>%
    dplyr::group_by(rank) %>%
    dplyr::summarise(
        rss_drop_null_mean = mean(rss_drop_from_previous_rank),
        rss_drop_null_min = min(rss_drop_from_previous_rank),
        rss_drop_null_max = max(rss_drop_from_previous_rank),
        .groups = "drop"
    ) %>%
    dplyr::left_join(
        rss_drop_long %>%
            dplyr::filter(data_type == "Observed") %>%
            dplyr::select(rank, rss_drop_observed = rss_drop_from_previous_rank),
        by = "rank"
    ) %>%
    dplyr::mutate(observed_minus_null_mean = rss_drop_observed - rss_drop_null_mean) %>%
    dplyr::relocate(rank, rss_drop_observed)


# 8. 提取rank 2-8的共识矩阵、W矩阵和癌种分组 ####
rank_assignment_list <- list()
rank_w_list <- list()
rank_pattern_label_list <- list()
rank_cluster_summary_list <- list()
rank_consensus_sheets <- list()

for (current_rank in compare_ranks) {
    current_fit <- observed_estimate$fit[[as.character(current_rank)]]

    ## 8.1 W与H矩阵（最佳拟合） ####
    current_w <- basis(current_fit)
    current_h <- coef(current_fit)
    pattern_ids <- paste0("P", seq_len(current_rank))
    colnames(current_w) <- pattern_ids
    rownames(current_w) <- rownames(mat)
    rownames(current_h) <- pattern_ids
    colnames(current_h) <- colnames(mat)

    current_h_normalized <- normalize_h_by_cancer(current_h)
    pattern_labels <- label_patterns_by_w(current_w, dominant_w_threshold)

    ## 8.2 一致性聚类与轮廓系数 ####
    current_consensus <- NMF::consensus(current_fit)
    current_clusters <- predict(current_fit, what = "consensus")
    current_clusters <- stats::setNames(as.integer(as.character(current_clusters)), colnames(mat))

    current_hc <- NMF::consensushc(current_fit, what = "consensus", method = hclust_method, dendrogram = FALSE)
    leaf_order <- current_hc$labels[current_hc$order]

    current_silhouette <- cluster::silhouette(
        current_clusters,
        stats::as.dist(1 - current_consensus[colnames(mat), colnames(mat)])
    )
    silhouette_values <- stats::setNames(current_silhouette[, "sil_width"], colnames(mat))

    cluster_to_pattern <- map_clusters_to_patterns(current_clusters, current_h_normalized)
    assigned_pattern <- cluster_to_pattern[as.character(current_clusters)]

    ## 8.3 汇总表 ####
    rank_assignment_list[[as.character(current_rank)]] <- tibble::tibble(
        rank = current_rank,
        cancer_type = colnames(mat),
        consensus_cluster = current_clusters,
        assigned_pattern = unname(assigned_pattern),
        assigned_pattern_label = unname(pattern_labels[assigned_pattern]),
        silhouette = unname(silhouette_values),
        max_normalized_h = apply(current_h_normalized, 1, max),
        leaf_order = match(colnames(mat), leaf_order),
        reference_group = unname(reference_groups[colnames(mat)])
    )

    rank_w_list[[as.character(current_rank)]] <- as.data.frame(current_w) %>%
        tibble::rownames_to_column("soc_abbrev") %>%
        tidyr::pivot_longer(-soc_abbrev, names_to = "pattern", values_to = "w_value") %>%
        dplyr::mutate(rank = current_rank, pattern_label = unname(pattern_labels[pattern]), .before = 1)

    rank_pattern_label_list[[as.character(current_rank)]] <- tibble::tibble(
        rank = current_rank,
        pattern = names(pattern_labels),
        pattern_label = unname(pattern_labels),
        n_dominant_soc = lengths(strsplit(unname(pattern_labels), "/", fixed = TRUE)),
        n_assigned_cancers = vapply(names(pattern_labels), function(x) sum(assigned_pattern == x), integer(1))
    )

    cluster_sizes <- table(current_clusters)
    rank_cluster_summary_list[[as.character(current_rank)]] <- tibble::tibble(
        rank = current_rank,
        n_clusters = length(cluster_sizes),
        min_cluster_size = min(cluster_sizes),
        n_singleton_clusters = sum(cluster_sizes == 1),
        mean_silhouette = mean(silhouette_values),
        min_silhouette = min(silhouette_values),
        n_cancers_silhouette_below_threshold = sum(silhouette_values < low_silhouette_threshold),
        n_single_soc_patterns = sum(lengths(strsplit(unname(pattern_labels), "/", fixed = TRUE)) == 1),
        max_w_cosine_similarity = calc_max_w_cosine(current_w),
        n_socs_without_dominant_pattern = sum(apply(current_w, 1, max) < dominant_w_threshold)
    )

    rank_consensus_sheets[[paste0("共识矩阵_秩", current_rank)]] <- as.data.frame(current_consensus[colnames(mat), colnames(mat)]) %>%
        tibble::rownames_to_column("cancer_type")
}

rank_assignment_df <- dplyr::bind_rows(rank_assignment_list)
rank_w_df <- dplyr::bind_rows(rank_w_list)
rank_pattern_label_df <- dplyr::bind_rows(rank_pattern_label_list)

rank_cluster_summary_df <- dplyr::bind_rows(rank_cluster_summary_list) %>%
    dplyr::left_join(
        rank_metrics_df %>%
            dplyr::filter(data_type == "Observed") %>%
            dplyr::select(rank, cophenetic, dispersion, silhouette.consensus, rss, evar),
        by = "rank"
    )


# 9. 跨秩分组一致性：与最终rank 5分组的列联表与ARI ####
rank_vs_reference_df <- rank_assignment_df %>%
    dplyr::group_by(rank, assigned_pattern, assigned_pattern_label, reference_group) %>%
    dplyr::summarise(
        n_cancers = dplyr::n(),
        cancers = paste(sort(cancer_type), collapse = ", "),
        .groups = "drop"
    ) %>%
    dplyr::mutate(reference_group = factor(reference_group, levels = reference_group_levels)) %>%
    dplyr::arrange(rank, assigned_pattern, reference_group) %>%
    dplyr::mutate(reference_group = as.character(reference_group))

assignment_wide <- rank_assignment_df %>%
    dplyr::select(rank, cancer_type, consensus_cluster) %>%
    tidyr::pivot_wider(names_from = rank, values_from = consensus_cluster, names_prefix = "rank_")

ari_df <- dplyr::bind_rows(
    lapply(compare_ranks, function(current_rank) {
        current_clusters <- assignment_wide[[paste0("rank_", current_rank)]]

        tibble::tibble(
            rank = current_rank,
            comparison = "vs final rank-5 groups (500 runs)",
            ari = calc_adjusted_rand_index(current_clusters, reference_groups[assignment_wide$cancer_type])
        )
    }),
    lapply(compare_ranks[-1], function(current_rank) {
        tibble::tibble(
            rank = current_rank,
            comparison = paste0("vs rank ", current_rank - 1),
            ari = calc_adjusted_rand_index(
                assignment_wide[[paste0("rank_", current_rank)]],
                assignment_wide[[paste0("rank_", current_rank - 1)]]
            )
        )
    })
)

message("各秩与最终rank 5分组的ARI:")
print(ari_df)
message("各秩聚类质量汇总:")
print(as.data.frame(rank_cluster_summary_df))


# 10. 结果保存 ####
output_sheets <- c(
    list(
        "秩选择指标_实际与各置换" = rank_metrics_df,
        "秩选择指标_实际与零分布" = null_comparison_df,
        "RSS下降_实际与零分布" = rss_drop_compare_df,
        "各秩聚类质量汇总" = rank_cluster_summary_df,
        "各秩癌种分组" = rank_assignment_df,
        "各秩模式注释" = rank_pattern_label_df,
        "各秩W矩阵" = rank_w_df,
        "与最终五分组列联表" = rank_vs_reference_df,
        "分组一致性ARI" = ari_df
    ),
    rank_consensus_sheets,
    list(
        "分析参数" = build_params_dataframe(list(
            input_matrix = "6.1.1主分析相同输入（FAERS ROR矩阵，癌种内min-max标准化，SOC×癌种）",
            observed_rank_survey = main_rank_estimate_file,
            reference_groups = main_wh_matrix_file,
            rank_range = rank_range,
            observed_rank_nrun = rank_nrun,
            seed = analysis_seed,
            method = "PNMF",
            randomization = "NMF::randomize，各癌种列内独立置换",
            n_permutations = n_permutations,
            permutation_nrun = permutation_nrun,
            permutation_seeds = paste0(analysis_seed, " + 置换编号(1-", n_permutations, ")"),
            empirical_p = "(零分布中不弱于实际值的置换数 + 1) / (置换数 + 1)；rss取≤，其余取≥",
            compare_ranks = compare_ranks,
            consensus_hclust_method = hclust_method,
            dominant_w_threshold = dominant_w_threshold,
            low_silhouette_threshold = low_silhouette_threshold
        ))
    )
)

writexl::write_xlsx(output_sheets, path = output_xlsx_file)

message("已生成结果文件:")
message(" - ", permutation_results_path)
message(" - ", output_xlsx_file)
message("========== NMF秩选择敏感性分析结束 ==========")
