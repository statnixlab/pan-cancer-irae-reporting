# ICI类别分层五模式稳定性分析（返修：审稿人1第12条） ####
# 目的：检验五个潜在模式及癌种分组在不同ICI类别中是否稳定
#      1) 主分析：固定最终模型的W矩阵，将各ICI类别分层的ROR谱投影到W上，按得分最高的模式归组，与最终五分组比较；
#      2) 补充分析：在抗PD-1单药层内重新进行rank 5 NMF，比较新模式与最终模式的相似度及癌种分组。
# 分层ROR直接读取3.3（同一ICI_regimen下primary_organ不成比例分析）的结果，不重新计算。
# 分层后稀疏单元增多，事件数为0的单元经连续性校正会得到虚高的ROR；主分析中将分层内事件数<3的单元ROR设为1（视为无信息），
# 并同时报告使用原始ROR的结果。

rm(list = ls())
gc()

# 1. 加载配置与NMF相关R包 ####
source("Code/config.R")
source("Code/mapping_rules.R")

suppressPackageStartupMessages(library(NMF))
suppressPackageStartupMessages(library(pNMF))

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
results_dir_name <- "6.2.7ICI类别分层五模式稳定性分析"
script_results_path <- file.path(results_path, results_dir_name)

if (!dir.exists(script_results_path)) {
    dir.create(script_results_path, recursive = TRUE, showWarnings = FALSE)
}

regimen_ror_file <- file.path(results_path, "3.3同一ICI_regimen下primary_organ不成比例分析", "FAERS_within_ICI_regimen_primary_organ_SOC_M0.xlsx")
main_ror_file <- file.path(results_path, "3.1primary_organ不成比例分析", "FAERS_primary_organ_SOC_M0.xlsx")
main_wh_matrix_file <- file.path(results_path, "6.1癌种NMF分群分析_rank5", "FAERS_NMF_WH_matrix.xlsx")
denovo_cache_files <- c(
    raw = "抗PD1单药层重新NMF结果.rdata",
    sparse_to_1 = "抗PD1单药层重新NMF结果_稀疏单元置1.rdata"
)
output_xlsx_file <- file.path(script_results_path, "ICI类别分层五模式稳定性分析结果.xlsx")

# 与主分析一致的参数 ####
nmf_rank <- 5L
nmf_nrun <- 200L
main_seed <- 123456L

# 分层分析参数 ####
stratum_min_cancer_n <- c(100L, 50L)
sparse_event_threshold <- 3L
sparse_handling_levels <- c("sparse_to_1", "raw")
denovo_regimen <- "anti-PD-1"
denovo_min_cancer_n <- 100L

pattern_levels <- c("Ocular", "Hematologic", "Neurologic", "Systemic-like", "Visceral-like")
regimen_levels <- c("anti-PD-1", "anti-PD-L1", "anti-CTLA-4+anti-PD-1")

message("========== ICI类别分层五模式稳定性分析开始 ==========")
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

## 函数：to_soc_abbrev ####
# 功能：将ROR矩阵中的SOC列名转换为SOC缩写####
# 输入：soc_names，字符型SOC列名向量####
# 输出：字符型SOC缩写向量####
to_soc_abbrev <- function(soc_names) {
    vapply(soc_names, function(x) {
        current_abbrev <- get_soc_abbrev(gsub("_", " ", x))
        if (is.na(current_abbrev) || current_abbrev == "") x else current_abbrev
    }, character(1), USE.NAMES = FALSE)
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

## 函数：solve_nonnegative_h ####
# 功能：在固定W下，对输入矩阵每一列求非负最小二乘得分（optim L-BFGS-B，下界为0）####
# 输入：input_matrix，SOC×癌种矩阵（行顺序与W一致）；basis_matrix，SOC×模式的固定W####
# 输出：模式×癌种的非负得分矩阵H####
solve_nonnegative_h <- function(input_matrix, basis_matrix) {
    h_matrix <- vapply(seq_len(ncol(input_matrix)), function(j) {
        target_vector <- input_matrix[, j]
        fit <- stats::optim(
            par = rep(0.1, ncol(basis_matrix)),
            fn = function(b) sum((target_vector - basis_matrix %*% b)^2),
            gr = function(b) as.numeric(-2 * crossprod(basis_matrix, target_vector - basis_matrix %*% b)),
            method = "L-BFGS-B",
            lower = 0,
            control = list(maxit = 1000)
        )
        fit$par
    }, numeric(ncol(basis_matrix)))

    h_matrix <- matrix(h_matrix, nrow = ncol(basis_matrix))
    rownames(h_matrix) <- colnames(basis_matrix)
    colnames(h_matrix) <- colnames(input_matrix)
    h_matrix
}

## 函数：project_onto_fixed_basis ####
# 功能：在固定W下求非负最小二乘得分H，并按得分最高的模式归组####
# 输入：input_matrix，SOC×癌种矩阵；basis_matrix，SOC×模式的固定W####
# 输出：癌种层面的投影结果数据框（各模式得分、首位与次位模式、主导模式占比）####
project_onto_fixed_basis <- function(input_matrix, basis_matrix) {
    projected_h <- solve_nonnegative_h(input_matrix[rownames(basis_matrix), , drop = FALSE], basis_matrix)

    tibble::tibble(
        cancer_type = colnames(projected_h),
        projected_group = rownames(projected_h)[apply(projected_h, 2, which.max)],
        second_group = rownames(projected_h)[apply(projected_h, 2, function(x) order(x, decreasing = TRUE)[2])],
        dominant_share = apply(projected_h, 2, function(x) ifelse(sum(x) > 0, max(x) / sum(x), NA_real_))
    ) %>%
        dplyr::bind_cols(tibble::as_tibble(t(projected_h)))
}

## 函数：calc_explained_variance ####
# 功能：计算固定W投影对输入矩阵的解释方差（1 - RSS / 总平方和）####
# 输入：input_matrix，SOC×癌种矩阵；basis_matrix，固定W####
# 输出：数值型解释方差####
calc_explained_variance <- function(input_matrix, basis_matrix) {
    aligned_matrix <- input_matrix[rownames(basis_matrix), , drop = FALSE]
    projected_h <- solve_nonnegative_h(aligned_matrix, basis_matrix)
    1 - sum((aligned_matrix - basis_matrix %*% projected_h)^2) / sum(aligned_matrix^2)
}


# 4. 读取最终模型与分层ROR ####
final_w <- readxl::read_xlsx(main_wh_matrix_file, sheet = "W_matrix") %>%
    tibble::column_to_rownames("soc_abbrev") %>%
    as.matrix()
final_w <- final_w[, pattern_levels, drop = FALSE]

final_h <- readxl::read_xlsx(main_wh_matrix_file, sheet = "H_matrix") %>%
    tibble::column_to_rownames("pattern_name") %>%
    as.matrix()

reference_assignment_df <- readxl::read_xlsx(main_wh_matrix_file, sheet = "Pattern group assignment") %>%
    dplyr::transmute(cancer_type = as.character(cancer_type), final_group = as.character(`Pattern group`))
final_groups <- stats::setNames(reference_assignment_df$final_group, reference_assignment_df$cancer_type)

regimen_ror_df <- readxl::read_xlsx(regimen_ror_file, sheet = "ROR矩阵")
regimen_meta_columns <- c("ici_regimen_raw", "ici_regimen_display", "regimen_rank", "primary_organ", "organ_n", "regimen_total_n", "organ_pct", "organ_rank", "included_in_main", "comparison_id")
regimen_soc_columns <- setdiff(colnames(regimen_ror_df), regimen_meta_columns)

regimen_ror_df <- regimen_ror_df %>%
    dplyr::filter(primary_organ %in% names(final_groups)) %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(regimen_soc_columns), ~ suppressWarnings(as.numeric(.x))))

main_ror_matrix <- readxl::read_xlsx(main_ror_file, sheet = "ROR矩阵") %>%
    dplyr::mutate(cancer_type = as.character(cancer_type)) %>%
    dplyr::filter(cancer_type %in% names(final_groups)) %>%
    dplyr::mutate(dplyr::across(-cancer_type, ~ suppressWarnings(as.numeric(.x)))) %>%
    tibble::column_to_rownames("cancer_type") %>%
    as.matrix()
main_ror_matrix <- main_ror_matrix[, !colnames(main_ror_matrix) %in% c("Infections_and_infestations", "Infections and Infestations"), drop = FALSE]
colnames(main_ror_matrix) <- to_soc_abbrev(colnames(main_ror_matrix))


# 5. 投影方法核查：最终模型H按得分最高归组、全样本投影归组，是否与最终五分组一致 ####
final_h_top_group <- rownames(final_h)[apply(final_h[, names(final_groups), drop = FALSE], 2, which.max)]
names(final_h_top_group) <- names(final_groups)

full_projection_df <- project_onto_fixed_basis(build_minmax_matrix(main_ror_matrix), final_w)

method_check_df <- tibble::tibble(
    cancer_type = names(final_groups),
    final_group = unname(final_groups),
    final_h_top_group = unname(final_h_top_group[names(final_groups)]),
    full_sample_projected_group = full_projection_df$projected_group[match(names(final_groups), full_projection_df$cancer_type)]
) %>%
    dplyr::mutate(
        final_h_top_same = final_h_top_group == final_group,
        full_projection_same = full_sample_projected_group == final_group
    )

message("最终模型H得分最高归组与最终五分组一致的癌种数: ", sum(method_check_df$final_h_top_same), "/", nrow(method_check_df))
message("全样本投影归组与最终五分组一致的癌种数: ", sum(method_check_df$full_projection_same), "/", nrow(method_check_df))


# 6. 各ICI类别分层的固定W投影 ####

## 6.1 读取分层内各癌种-SOC的事件数（a），用于识别稀疏单元 ####
regimen_event_df <- readxl::read_xlsx(regimen_ror_file, sheet = "详细结果") %>%
    dplyr::filter(primary_organ %in% names(final_groups)) %>%
    dplyr::transmute(ici_regimen_raw, primary_organ, soc_abbrev = irae_abbrev, n_events = suppressWarnings(as.numeric(a)))

## 函数：build_stratum_ror_matrix ####
# 功能：构建某一ICI类别分层的癌种×SOC ROR矩阵；按需将事件数低于阈值的单元ROR设为1（视为无信息）####
# 输入：regimen_name，分层名称；min_n，纳入阈值；sparse_handling，"sparse_to_1"或"raw"####
# 输出：列表，含ROR矩阵、癌种报告数和被置1的单元数####
build_stratum_ror_matrix <- function(regimen_name, min_n, sparse_handling) {
    stratum_df <- regimen_ror_df %>%
        dplyr::filter(ici_regimen_raw == regimen_name, organ_n >= min_n)

    ror_matrix <- stratum_df %>%
        dplyr::select(primary_organ, dplyr::all_of(regimen_soc_columns)) %>%
        tibble::column_to_rownames("primary_organ") %>%
        as.matrix()
    colnames(ror_matrix) <- to_soc_abbrev(colnames(ror_matrix))
    ror_matrix <- ror_matrix[, rownames(final_w), drop = FALSE]

    if (any(!is.finite(ror_matrix)) || any(ror_matrix <= 0)) {
        stop("分层ROR存在缺失或非正值: ", regimen_name)
    }

    event_matrix <- regimen_event_df %>%
        dplyr::filter(ici_regimen_raw == regimen_name, primary_organ %in% rownames(ror_matrix)) %>%
        dplyr::select(primary_organ, soc_abbrev, n_events) %>%
        tidyr::pivot_wider(names_from = soc_abbrev, values_from = n_events) %>%
        tibble::column_to_rownames("primary_organ") %>%
        as.matrix()
    event_matrix <- event_matrix[rownames(ror_matrix), colnames(ror_matrix), drop = FALSE]

    sparse_cells <- event_matrix < sparse_event_threshold
    if (sparse_handling == "sparse_to_1") {
        ror_matrix[sparse_cells] <- 1
    }

    list(
        ror_matrix = ror_matrix,
        organ_n = stats::setNames(stratum_df$organ_n, stratum_df$primary_organ),
        n_sparse_cells = sum(sparse_cells)
    )
}

stratum_assignment_list <- list()
stratum_summary_list <- list()

for (current_handling in sparse_handling_levels) {
    for (current_regimen in regimen_levels) {
        for (current_min_n in stratum_min_cancer_n) {
            current_stratum <- build_stratum_ror_matrix(current_regimen, current_min_n, current_handling)
            current_input <- build_minmax_matrix(current_stratum$ror_matrix)

            current_projection <- project_onto_fixed_basis(current_input, final_w) %>%
                dplyr::mutate(
                    sparse_handling = current_handling,
                    ici_regimen = current_regimen,
                    min_cancer_n = current_min_n,
                    organ_n = unname(current_stratum$organ_n[cancer_type]),
                    final_group = unname(final_groups[cancer_type]),
                    same_as_final = projected_group == final_group,
                    .before = 1
                )

            current_key <- paste(current_handling, current_regimen, current_min_n)
            stratum_assignment_list[[current_key]] <- current_projection

            stratum_summary_list[[current_key]] <- tibble::tibble(
                sparse_handling = current_handling,
                ici_regimen = current_regimen,
                min_cancer_n = current_min_n,
                n_cancers = nrow(current_projection),
                n_sparse_cells = current_stratum$n_sparse_cells,
                prop_sparse_cells = current_stratum$n_sparse_cells / length(current_stratum$ror_matrix),
                n_final_groups_represented = dplyr::n_distinct(current_projection$final_group),
                excluded_cancers = paste(setdiff(names(final_groups), current_projection$cancer_type), collapse = ", "),
                ari_vs_final = calc_adjusted_rand_index(current_projection$projected_group, current_projection$final_group),
                n_cancers_retained = sum(current_projection$same_as_final),
                n_cancers_reassigned = sum(!current_projection$same_as_final),
                reassigned_cancers = paste(
                    paste0(current_projection$cancer_type, " (", current_projection$final_group, " -> ", current_projection$projected_group, ")")[!current_projection$same_as_final],
                    collapse = "; "
                ),
                explained_variance_fixed_w = calc_explained_variance(current_input, final_w)
            )
        }
    }
}

stratum_assignment_df <- dplyr::bind_rows(stratum_assignment_list)
stratum_summary_df <- dplyr::bind_rows(stratum_summary_list)

message("各分层固定W投影汇总:")
print(as.data.frame(stratum_summary_df %>% dplyr::select(sparse_handling, ici_regimen, min_cancer_n, n_cancers, n_sparse_cells, ari_vs_final, n_cancers_retained, explained_variance_fixed_w)))


# 7. 补充分析：抗PD-1单药层重新NMF（原始ROR与稀疏单元置1两个版本） ####
denovo_summary_list <- list()
denovo_assignment_list <- list()
denovo_match_list <- list()

for (current_handling in sparse_handling_levels) {
    denovo_stratum <- build_stratum_ror_matrix(denovo_regimen, denovo_min_cancer_n, current_handling)
    denovo_input <- build_minmax_matrix(denovo_stratum$ror_matrix)
    denovo_cache_file <- file.path(script_results_path, denovo_cache_files[[current_handling]])

    if (file.exists(denovo_cache_file)) {
        load(denovo_cache_file)
        message("读取已有抗PD-1层重新NMF结果: ", denovo_cache_file)
    } else {
        denovo_fit <- nmf(denovo_input, rank = nmf_rank, nrun = nmf_nrun, seed = main_seed, method = "PNMF")
        save(denovo_fit, file = denovo_cache_file)
    }

    denovo_clusters <- predict(denovo_fit, what = "consensus")
    denovo_clusters <- stats::setNames(as.integer(as.character(denovo_clusters)), colnames(denovo_input))

    ## 7.1 新模式与最终模式的余弦相似度，按相似度贪心一对一匹配 ####
    denovo_w <- basis(denovo_fit)
    rownames(denovo_w) <- rownames(denovo_input)
    colnames(denovo_w) <- paste0("P", seq_len(ncol(denovo_w)))

    cosine_matrix <- crossprod(
        sweep(final_w, 2, sqrt(colSums(final_w^2)), "/"),
        sweep(denovo_w, 2, sqrt(colSums(denovo_w^2)), "/")
    )

    current_match <- list()
    remaining_final <- rownames(cosine_matrix)
    remaining_denovo <- colnames(cosine_matrix)
    while (length(remaining_final) > 0) {
        sub_cosine <- cosine_matrix[remaining_final, remaining_denovo, drop = FALSE]
        best_index <- which(sub_cosine == max(sub_cosine), arr.ind = TRUE)[1, ]
        current_match[[length(current_match) + 1]] <- tibble::tibble(
            sparse_handling = current_handling,
            final_pattern = remaining_final[best_index[1]],
            denovo_pattern = remaining_denovo[best_index[2]],
            cosine_similarity = sub_cosine[best_index[1], best_index[2]]
        )
        remaining_final <- remaining_final[-best_index[1]]
        remaining_denovo <- remaining_denovo[-best_index[2]]
    }
    current_match_df <- dplyr::bind_rows(current_match) %>%
        dplyr::mutate(final_pattern = factor(final_pattern, levels = pattern_levels)) %>%
        dplyr::arrange(final_pattern) %>%
        dplyr::mutate(final_pattern = as.character(final_pattern))

    ## 7.2 新分组与最终五分组比较：各簇映射为簇内多数癌种的最终分组 ####
    denovo_mapped_group <- rep(NA_character_, length(denovo_clusters))
    names(denovo_mapped_group) <- names(denovo_clusters)
    for (current_cluster in unique(denovo_clusters)) {
        member_cancers <- names(denovo_clusters)[denovo_clusters == current_cluster]
        group_counts <- table(final_groups[member_cancers])
        top_groups <- names(group_counts)[group_counts == max(group_counts)]
        denovo_mapped_group[member_cancers] <- if (length(top_groups) == 1) top_groups else "Mixed"
    }

    current_assignment_df <- tibble::tibble(
        sparse_handling = current_handling,
        cancer_type = names(denovo_clusters),
        organ_n = unname(denovo_stratum$organ_n[names(denovo_clusters)]),
        final_group = unname(final_groups[names(denovo_clusters)]),
        denovo_cluster = unname(denovo_clusters),
        denovo_mapped_group = unname(denovo_mapped_group),
        same_as_final = denovo_mapped_group == final_group
    )

    denovo_summary_list[[current_handling]] <- tibble::tibble(
        sparse_handling = current_handling,
        ici_regimen = denovo_regimen,
        min_cancer_n = denovo_min_cancer_n,
        n_cancers = length(denovo_clusters),
        n_sparse_cells = denovo_stratum$n_sparse_cells,
        ari_vs_final = calc_adjusted_rand_index(denovo_clusters, final_groups[names(denovo_clusters)]),
        n_cancers_retained = sum(current_assignment_df$same_as_final),
        n_cancers_reassigned = sum(!current_assignment_df$same_as_final),
        reassigned_cancers = paste(current_assignment_df$cancer_type[!current_assignment_df$same_as_final], collapse = ", "),
        min_matched_cosine = min(current_match_df$cosine_similarity),
        cophenetic = cophcor(denovo_fit)
    )
    denovo_assignment_list[[current_handling]] <- current_assignment_df
    denovo_match_list[[current_handling]] <- current_match_df
}

denovo_summary_df <- dplyr::bind_rows(denovo_summary_list)
denovo_assignment_df <- dplyr::bind_rows(denovo_assignment_list)
denovo_pattern_match_df <- dplyr::bind_rows(denovo_match_list)

message("抗PD-1层重新NMF：")
print(as.data.frame(denovo_summary_df))
print(as.data.frame(denovo_pattern_match_df))


# 8. 结果保存 ####
writexl::write_xlsx(
    list(
        "投影方法核查" = method_check_df,
        "分层投影汇总" = stratum_summary_df,
        "分层投影癌种分组" = stratum_assignment_df,
        "抗PD1重新NMF汇总" = denovo_summary_df,
        "抗PD1重新NMF分组" = denovo_assignment_df,
        "抗PD1新模式与最终模式匹配" = denovo_pattern_match_df,
        "分析参数" = build_params_dataframe(list(
            regimen_ror_source = regimen_ror_file,
            final_model = main_wh_matrix_file,
            regimens = regimen_levels,
            stratum_min_cancer_n = stratum_min_cancer_n,
            sparse_handling = "sparse_to_1（主）：分层内事件数<阈值的单元ROR设为1，视为无信息；raw：使用原始ROR",
            sparse_event_threshold = sparse_event_threshold,
            preprocessing = "各分层内ROR按癌种min-max标准化（与主分析一致）",
            projection = "固定最终W，非负最小二乘（optim L-BFGS-B，下界0）求H，按H最大的模式归组",
            explained_variance = "1 - RSS / 输入矩阵平方和",
            denovo_regimen = denovo_regimen,
            denovo_min_cancer_n = denovo_min_cancer_n,
            denovo_nmf = paste0("rank ", nmf_rank, ", nrun ", nmf_nrun, ", seed ", main_seed, ", PNMF"),
            denovo_pattern_matching = "新W与最终W按余弦相似度贪心一对一匹配",
            denovo_cluster_mapping = "新分组各簇映射为簇内多数癌种的最终分组，并列最多记为Mixed"
        ))
    ),
    path = output_xlsx_file
)

message("已生成结果文件: ", output_xlsx_file)
message("========== ICI类别分层五模式稳定性分析结束 ==========")
