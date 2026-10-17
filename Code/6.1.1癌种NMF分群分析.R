# 癌种NMF分群分析（rank=5最终版） ####

rm(list = ls())
gc()

# 1. 加载配置与NMF相关R包 ####
source("Code/config.R")
source("Code/mapping_rules.R")

if (!requireNamespace("NMF", quietly = TRUE)) {
    install.packages("NMF", quiet = TRUE)
}

if (!requireNamespace("pNMF", quietly = TRUE)) {
    install.packages("pNMF", quiet = TRUE)
}

if (!requireNamespace("pheatmap", quietly = TRUE)) {
    install.packages("pheatmap", quiet = TRUE)
}

if (!requireNamespace("dendextend", quietly = TRUE)) {
    install.packages("dendextend", quiet = TRUE)
}

suppressPackageStartupMessages(library(NMF))
suppressPackageStartupMessages(library(pNMF))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(dendextend))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))

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

if (identical(pnmf_registration_status, "already_registered")) {
    message("PNMF算法已在当前R会话中注册，跳过重复注册")
}


# 2. 结果目录与参数 ####
results_dir_name <- "6.1癌种NMF分群分析_rank5"
legacy_results_dir_name <- "6.1癌种NMF分群分析_rank5尝试"
script_results_path <- file.path(results_path, results_dir_name)
legacy_results_path <- file.path(results_path, legacy_results_dir_name)

if (!dir.exists(script_results_path) && dir.exists(legacy_results_path)) {
    script_results_path <- legacy_results_path
}

if (!dir.exists(script_results_path)) {
    dir.create(script_results_path, recursive = TRUE, showWarnings = FALSE)
}

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

rank_estimate_file <- file.path(script_results_path, "FAERS_NMF_rank_estimate.rdata")
nmf_result_file <- file.path(script_results_path, "FAERS_NMF_result.rdata")
best_fit_file <- file.path(script_results_path, "FAERS_NMF_best_fit.rdata")
wh_matrix_file <- file.path(script_results_path, "FAERS_NMF_WH_matrix.xlsx")
h_cluster_file <- file.path(script_results_path, "FAERS_NMF_H_cluster.rdata")

rank_plot_file <- file.path(script_results_path, "图1_FAERS_NMF秩选择图.pdf")
consensus_map_file <- file.path(script_results_path, "图2_FAERS_NMF共识矩阵图.pdf")
w_heatmap_file <- file.path(script_results_path, "图3_FAERS_NMF_W热图.pdf")
h_tree_file <- file.path(script_results_path, "图4_FAERS_NMF_H层次聚类树.pdf")
h_heatmap_file <- file.path(script_results_path, "图5_FAERS_NMF_H热图.pdf")

rank_range <- 2:8
rank_nrun <- 200L
nmf_rank_manual <- 5L
nmf_final_nrun <- 500L
analysis_seed <- 123456L
hclust_method <- "average"
hclust_distance_method <- "euclidean"
cluster_k <- nmf_rank_manual
min_cancer_n <- 100L

message("========== 癌种NMF分群分析（rank=5最终版）开始 ==========")
message("结果目录: ", script_results_path)
message("输入文件: ", ror_input_file)
message("样本量汇总文件: ", summary_input_file)
message("纳入癌种最小样本量: ", min_cancer_n)
message("秩搜索范围: ", paste(rank_range, collapse = ", "))
message("当前手动rank: ", nmf_rank_manual)


# 3. 读取原始ROR矩阵 ####
if (!file.exists(ror_input_file)) {
    stop("原始ROR输入文件不存在: ", ror_input_file)
}

if (!file.exists(summary_input_file)) {
    stop("3.1汇总文件不存在: ", summary_input_file)
}

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

raw_ror_matrix_df <- raw_ror_matrix_df %>%
    dplyr::mutate(cancer_type = as.character(cancer_type)) %>%
    dplyr::filter(cancer_type %in% qualified_cancers)

raw_ror_matrix_df <- raw_ror_matrix_df %>%
    dplyr::mutate(
        cancer_type = as.character(cancer_type),
        dplyr::across(
            -cancer_type,
            ~ suppressWarnings(as.numeric(.x))
        )
    )

if (nrow(raw_ror_matrix_df) == 0) {
    stop("按样本量阈值筛选后没有可用于NMF的癌种")
}

message("纳入NMF的癌种数: ", nrow(raw_ror_matrix_df))
message("纳入NMF的癌种: ", paste(raw_ror_matrix_df$cancer_type, collapse = ", "))

raw_ror_matrix <- raw_ror_matrix_df %>%
    tibble::column_to_rownames("cancer_type") %>%
    as.matrix()
storage.mode(raw_ror_matrix) <- "double"

if (any(is.na(raw_ror_matrix))) {
    stop("原始ROR矩阵存在缺失值，请先检查上游结果表")
}

if (any(raw_ror_matrix <= 0)) {
    stop("原始ROR矩阵存在小于等于0的值，无法进行NMF分析")
}


# 4. 构建NMF输入矩阵 ####
mat <- t(raw_ror_matrix)

if ("Infections_and_infestations" %in% rownames(mat)) {
    mat <- mat[!rownames(mat) %in% c("Infections_and_infestations"), , drop = FALSE]
}

if ("Infections and Infestations" %in% rownames(mat)) {
    mat <- mat[!rownames(mat) %in% c("Infections and Infestations"), , drop = FALSE]
}

mat_row_names <- rownames(mat)
mat_col_names <- colnames(mat)

mat <- apply(mat, 2, function(x) {
    current_min <- min(x)
    current_max <- max(x)

    if (isTRUE(all.equal(current_min, current_max))) {
        rep(0, length(x))
    } else {
        (x - current_min) / (current_max - current_min)
    }
})

mat <- as.matrix(mat)
rownames(mat) <- mat_row_names
colnames(mat) <- mat_col_names


# 5. rank survey ####
if (file.exists(rank_estimate_file)) {
    load(rank_estimate_file)
    message("检测到已有rank survey结果，直接读取: ", rank_estimate_file)
} else {
    estim.r <- nmf(
        mat,
        rank_range,
        nrun = rank_nrun,
        seed = analysis_seed,
        method = "PNMF"
    )

    save(estim.r, file = rank_estimate_file)
}


# 6. 运行最终NMF模型 ####
# 看完图1后，如需修改rank，直接改上面的 nmf_rank_manual 再重新运行本脚本 ####
if (file.exists(nmf_result_file) && file.exists(best_fit_file)) {
    load(nmf_result_file)
    load(best_fit_file)
    message("检测到已有最终NMF结果，直接读取: ", nmf_result_file)
    message("检测到已有best_fit结果，直接读取: ", best_fit_file)
} else {
    nmf_result <- nmf(
        mat,
        rank = nmf_rank_manual,
        nrun = nmf_final_nrun,
        seed = analysis_seed,
        method = "PNMF"
    )

    best_fit <- minfit(nmf_result)

    save(nmf_result, file = nmf_result_file)
    save(best_fit, file = best_fit_file)
}


# 7. 提取W和H矩阵 ####
W <- basis(best_fit)
H <- coef(best_fit)

nmf_rank_used <- ncol(W)

if (!identical(as.integer(nmf_rank_manual), as.integer(nmf_rank_used))) {
    message(
        "当前脚本设置的 nmf_rank_manual = ",
        nmf_rank_manual,
        "，但已读取结果的实际rank = ",
        nmf_rank_used,
        "；后续输出将按实际rank解释"
    )
}

if (nmf_rank_used != 5) {
    stop("当前英文模式命名仅适用于 rank = 5，请先确认 best_fit 的实际rank")
}

pattern_names <- c(
    "Hematologic",
    "Ocular",
    "Neurologic",
    "Systemic-like",
    "Visceral-like"
)

# 如需手动调整5个主群从左到右的显示顺序，直接修改下面这一行 ####
pattern_display_order <- c(
    "Ocular",
    "Hematologic",
    "Neurologic",
    "Systemic-like",
    "Visceral-like"
)

if (!setequal(pattern_display_order, pattern_names) || length(pattern_display_order) != length(pattern_names)) {
    stop("pattern_display_order 必须与 pattern_names 完全一致，只允许调整顺序")
}

colnames(W) <- pattern_names
rownames(H) <- pattern_names

if (is.null(rownames(W))) {
    rownames(W) <- rownames(mat)
}

if (is.null(colnames(H))) {
    colnames(H) <- colnames(mat)
}

W <- W[, pattern_display_order, drop = FALSE]
H <- H[pattern_display_order, , drop = FALSE]
pattern_names <- pattern_display_order

soc_display_names <- rownames(W)
soc_display_names <- gsub("_", " ", soc_display_names)
soc_display_names <- vapply(
    soc_display_names,
    function(x) {
        current_abbrev <- get_soc_abbrev(x)

        if (is.na(current_abbrev) || current_abbrev == "") {
            x
        } else {
            current_abbrev
        }
    },
    character(1)
)
rownames(W) <- soc_display_names

W_df <- as.data.frame(W, check.names = FALSE) %>%
    tibble::rownames_to_column("soc_abbrev")

H_df <- as.data.frame(H, check.names = FALSE) %>%
    tibble::rownames_to_column("pattern_name")


# 8. 基于一致性聚类进行癌种分组 ####
H_normalized <- H

for (current_cancer in seq_len(ncol(H_normalized))) {
    current_values <- H_normalized[, current_cancer]
    current_min <- min(current_values)
    current_max <- max(current_values)

    if (isTRUE(all.equal(current_min, current_max))) {
        H_normalized[, current_cancer] <- 0
    } else {
        H_normalized[, current_cancer] <- (current_values - current_min) / (current_max - current_min)
    }
}

H_normalized <- t(H_normalized)

colnames(H_normalized) <- pattern_names
rownames(H_normalized) <- colnames(mat)

cluster_k <- nmf_rank_used

consensus_clusters <- predict(nmf_result, what = "consensus")
consensus_clusters <- factor(
    as.character(consensus_clusters),
    levels = levels(consensus_clusters)
)
names(consensus_clusters) <- rownames(H_normalized)

hc <- NMF::consensushc(
    nmf_result,
    what = "consensus",
    method = hclust_method,
    dendrogram = FALSE
)

clusters <- as.integer(consensus_clusters)
names(clusters) <- names(consensus_clusters)

cluster_ids <- sort(unique(as.integer(clusters)))
cluster_pattern_score <- matrix(
    0,
    nrow = length(cluster_ids),
    ncol = ncol(H_normalized),
    dimnames = list(
        paste0("Cluster_", cluster_ids),
        pattern_names
    )
)

for (i in seq_along(cluster_ids)) {
    current_cluster_id <- cluster_ids[i]
    current_cluster_matrix <- H_normalized[clusters == current_cluster_id, , drop = FALSE]

    if (nrow(current_cluster_matrix) == 1) {
        cluster_pattern_score[i, ] <- as.numeric(current_cluster_matrix[1, ])
    } else {
        cluster_pattern_score[i, ] <- colMeans(current_cluster_matrix)
    }
}

cluster_to_pattern <- rep(NA_character_, length(cluster_ids))
names(cluster_to_pattern) <- as.character(cluster_ids)
available_patterns <- pattern_names
cluster_order <- order(apply(cluster_pattern_score, 1, max), decreasing = TRUE)

for (i in cluster_order) {
    current_cluster_id <- cluster_ids[i]
    current_pattern_order <- order(cluster_pattern_score[i, ], decreasing = TRUE)

    for (j in current_pattern_order) {
        current_pattern_name <- pattern_names[j]

        if (current_pattern_name %in% available_patterns) {
            cluster_to_pattern[as.character(current_cluster_id)] <- current_pattern_name
            available_patterns <- setdiff(available_patterns, current_pattern_name)
            break
        }
    }
}

if (any(is.na(cluster_to_pattern))) {
    remaining_clusters <- names(cluster_to_pattern)[is.na(cluster_to_pattern)]

    for (i in seq_along(remaining_clusters)) {
        cluster_to_pattern[remaining_clusters[i]] <- available_patterns[i]
    }
}

cancer_pattern_groups <- cluster_to_pattern[as.character(clusters)]
names(cancer_pattern_groups) <- names(clusters)

current_leaf_order <- hc$labels[hc$order]
desired_leaf_order <- unlist(
    lapply(
        pattern_display_order,
        function(current_pattern_name) {
            current_leaf_order[cancer_pattern_groups[current_leaf_order] == current_pattern_name]
        }
    ),
    use.names = FALSE
)

if (length(desired_leaf_order) != length(current_leaf_order) || !setequal(desired_leaf_order, current_leaf_order)) {
    stop("无法根据 pattern_display_order 生成完整的癌种显示顺序")
}

hc <- stats::as.hclust(
    dendextend::rotate(
        as.dendrogram(hc),
        order = desired_leaf_order
    )
)

message("一致性聚类与模式名称对应关系:")
for (current_cluster_id in names(cluster_to_pattern)) {
    message(" - Cluster ", current_cluster_id, ": ", cluster_to_pattern[current_cluster_id])
}

save(
    H_normalized,
    hc,
    clusters,
    consensus_clusters,
    cluster_to_pattern,
    cancer_pattern_groups,
    file = h_cluster_file
)

H_normalized_df <- as.data.frame(H_normalized, check.names = FALSE) %>%
    tibble::rownames_to_column("cancer_type")

cluster_assignment_df <- data.frame(
    cancer_type = names(clusters),
    "Pattern group" = cancer_pattern_groups,
    stringsAsFactors = FALSE,
    check.names = FALSE
) %>%
    dplyr::mutate(`Pattern group` = factor(`Pattern group`, levels = pattern_names)) %>%
    dplyr::arrange(`Pattern group`, cancer_type) %>%
    dplyr::mutate(`Pattern group` = as.character(`Pattern group`))

writexl::write_xlsx(
    list(
        "W_matrix" = W_df,
        "H_matrix" = H_df,
        "H_normalized" = H_normalized_df,
        "Pattern group assignment" = cluster_assignment_df
    ),
    path = wh_matrix_file
)


# 9. 1图：NMF秩选择图 ####

## 1. 读取并处理绘图数据 ####
load(rank_estimate_file)

## 2. 绘图 ####
rank_plot_width_cm <- 4.85
rank_plot_height_cm <- rank_plot_width_cm

grDevices::cairo_pdf(
    filename = rank_plot_file,
    width = rank_plot_width_cm / 2.54,
    height = rank_plot_height_cm / 2.54,
    pointsize = 7
)

p_rank <- plot(
    estim.r,
    what = "cophenetic",
    main = NULL
) +
    ggplot2::theme(
        text = ggplot2::element_text(size = 7, colour = "black"),
        legend.position = "none",
        plot.title = ggplot2::element_blank()
    )

print(p_rank)

dev.off()


# 10. 2图：NMF共识矩阵图 ####

## 1. 读取并处理绘图数据 ####
load(nmf_result_file)
load(h_cluster_file)

consensus_plot_matrix <- NMF::consensus(nmf_result)
consensus_plot_matrix <- consensus_plot_matrix[hc$labels, hc$labels, drop = FALSE]

H_consensus_plot_df <- readxl::read_xlsx(wh_matrix_file, sheet = "H_matrix")
H_consensus_plot_matrix <- H_consensus_plot_df %>%
    as.data.frame(check.names = FALSE) %>%
    tibble::column_to_rownames("pattern_name") %>%
    as.matrix()

basis_id_lookup <- stats::setNames(as.character(seq_along(pattern_names)), pattern_names)
basis_pattern_assignment <- apply(
    H_consensus_plot_matrix[, colnames(consensus_plot_matrix), drop = FALSE],
    2,
    function(x) {
        pattern_names[which.max(x)]
    }
)
basis_annotation <- factor(
    basis_id_lookup[basis_pattern_assignment],
    levels = as.character(seq_along(pattern_names))
)

consensus_cluster_annotation <- factor(
    as.character(consensus_clusters[colnames(consensus_plot_matrix)]),
    levels = as.character(sort(unique(as.integer(consensus_clusters))))
)

consensus_annotation <- factor(
    cancer_pattern_groups[colnames(consensus_plot_matrix)],
    levels = pattern_names
)

consensus_distance <- stats::as.dist(1 - consensus_plot_matrix)
consensus_silhouette <- cluster::silhouette(
    as.integer(consensus_cluster_annotation),
    consensus_distance
)
silhouette_annotation <- consensus_silhouette[, "sil_width"]
names(silhouette_annotation) <- colnames(consensus_plot_matrix)

basis_colors <- c(
    "1" = "#9999FF",
    "2" = "#E76BF3",
    "3" = "#32B0FF",
    "4" = "#C77CFF",
    "5" = "#00BA38"
)
basis_colors <- basis_colors[levels(basis_annotation)]

consensus_colors <- c(
    "Ocular" = "#DF7A5E",
    "Hematologic" = "#A15F88",
    "Neurologic" = "#82B29A",
    "Systemic-like" = "#F2CC8E",
    "Visceral-like" = "#BC4749"
)
consensus_colors <- consensus_colors[levels(consensus_annotation)]

consensus_color_fun <- circlize::colorRamp2(
    c(0, 0.2, 0.4, 0.6, 0.8, 1),
    c("#313695", "#74ADD1", "#E0F3F8", "#FEE090", "#F46D43", "#A50026")
)

silhouette_min <- floor(min(silhouette_annotation, na.rm = TRUE) * 100) / 100
silhouette_color_fun <- circlize::colorRamp2(
    c(silhouette_min, 1),
    c("#F1F1F1", "#009E73")
)

consensus_top_annotation <- ComplexHeatmap::HeatmapAnnotation(
    basis = basis_annotation,
    consensus = consensus_annotation,
    silhouette = silhouette_annotation,
    col = list(
        basis = basis_colors,
        consensus = consensus_colors,
        silhouette = silhouette_color_fun
    ),
    show_annotation_name = TRUE,
    annotation_name_side = "right",
    annotation_name_rot = 0,
    annotation_name_gp = grid::gpar(fontsize = 7, fontface = "bold", col = "black"),
    show_legend = FALSE,
    simple_anno_size = grid::unit(1.8, "mm")
)

## 2. 绘图 ####
consensus_map_width_cm <- 9.7
consensus_map_height_cm <- consensus_map_width_cm

grDevices::cairo_pdf(
    filename = consensus_map_file,
    width = consensus_map_width_cm / 2.54,
    height = consensus_map_height_cm / 2.54,
    pointsize = 7
)

consensus_heatmap <- ComplexHeatmap::Heatmap(
    consensus_plot_matrix,
    name = "Consensus",
    col = consensus_color_fun,
    top_annotation = consensus_top_annotation,
    cluster_rows = stats::as.dendrogram(hc),
    cluster_columns = stats::as.dendrogram(hc),
    row_names_side = "right",
    show_row_names = TRUE,
    show_column_names = TRUE,
    column_names_rot = 90,
    width = grid::unit(6.0, "cm"),
    height = grid::unit(5.4, "cm"),
    row_dend_width = grid::unit(0.35, "cm"),
    column_dend_height = grid::unit(0.35, "cm"),
    row_names_gp = grid::gpar(fontsize = 7, col = "black"),
    column_names_gp = grid::gpar(fontsize = 7, col = "black"),
    column_names_max_height = grid::unit(2.4, "cm"),
    column_title = NULL,
    rect_gp = grid::gpar(col = "#4B4B7A", lwd = 0.2),
    show_heatmap_legend = TRUE,
    heatmap_legend_param = list(
        title = "Consensus",
        direction = "horizontal",
        title_position = "leftcenter",
        at = c(0, 0.5, 1),
        labels = c("0", "0.5", "1"),
        legend_width = grid::unit(2.6, "cm"),
        grid_height = grid::unit(1.4, "mm"),
        labels_gp = grid::gpar(fontsize = 7, col = "black"),
        title_gp = grid::gpar(fontsize = 7, fontface = "bold", col = "black")
    )
)

draw(
    consensus_heatmap,
    heatmap_legend_side = "top",
    annotation_legend_side = "top",
    padding = grid::unit(c(0, 0, 0, 0), "mm")
)

dev.off()


# 11. 3图：W矩阵热图 ####

## 1. 读取并处理绘图数据 ####
W_plot_df <- readxl::read_xlsx(wh_matrix_file, sheet = "W_matrix")
W_plot_matrix <- W_plot_df %>%
    as.data.frame(check.names = FALSE) %>%
    tibble::column_to_rownames("soc_abbrev") %>%
    as.matrix()

## 2. 绘图 ####
w_heatmap_width_cm <- 4.85
w_heatmap_height_cm <- w_heatmap_width_cm
w_heatmap_layout_scale <- w_heatmap_width_cm / (4 * 2.54)

grDevices::cairo_pdf(
    filename = w_heatmap_file,
    width = w_heatmap_width_cm / 2.54,
    height = w_heatmap_height_cm / 2.54,
    pointsize = 7
)

W_heatmap <- ComplexHeatmap::Heatmap(
    W_plot_matrix,
    name = "Basis weight",
    column_title = "irAE reporting pattern",
    column_title_side = "bottom",
    row_title_side = "left",
    row_names_side = "left",
    col = circlize::colorRamp2(
        seq(min(W_plot_matrix), max(W_plot_matrix), length.out = 5),
        c("#FFF6F3", "#FFCEBF", "#FF472F", "#DA2A19", "#741616")
    ),
    width = grid::unit(7.8 * w_heatmap_layout_scale, "cm"),
    height = grid::unit(5.8 * w_heatmap_layout_scale, "cm"),
    heatmap_legend_param = list(
        title = "Basis weight",
        direction = "horizontal",
        title_position = "leftcenter",
        legend_width = grid::unit(4 * w_heatmap_layout_scale, "cm"),
        grid_width = grid::unit(1 * w_heatmap_layout_scale, "cm"),
        grid_height = grid::unit(1.2, "mm"),
        labels_gp = grid::gpar(fontsize = 7, col = "black"),
        title_gp = grid::gpar(fontsize = 7, col = "black")
    ),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    cell_fun = function(j, i, x, y, width, height, fill) {
        grid::grid.text(
            sprintf("%.2f", W_plot_matrix[i, j]),
            x,
            y,
            gp = grid::gpar(fontsize = 6, col = "black")
        )
    },
    row_names_gp = grid::gpar(fontsize = 7, col = "black"),
    row_names_max_width = grid::unit(20 * w_heatmap_layout_scale, "mm"),
    column_names_gp = grid::gpar(fontsize = 7, col = "black"),
    column_names_max_height = grid::unit(46 * w_heatmap_layout_scale, "mm"),
    column_names_rot = 45,
    column_title_gp = grid::gpar(fontsize = 7, fontface = "bold", col = "black"),
    row_title_gp = grid::gpar(fontsize = 7, fontface = "bold", col = "black"),
    rect_gp = grid::gpar(col = "grey60", lwd = 1)
)

draw(
    W_heatmap,
    heatmap_legend_side = "top",
    annotation_legend_side = "top",
    padding = grid::unit(c(8, 12, 20, 12) * w_heatmap_layout_scale, "mm")
)

dev.off()


# 12. 4图：H矩阵层次聚类树 ####

## 1. 读取并处理绘图数据 ####
load(h_cluster_file)

## 2. 绘图 ####
grDevices::cairo_pdf(
    filename = h_tree_file,
    width = 9,
    height = 3
)
plot(
    hc,
    main = "Consensus Clustering Dendrogram of Cancer Types",
    xlab = "",
    ylab = "Height",
    hang = -1
)
dev.off()


# 13. 5图：H矩阵热图 ####

## 1. 读取并处理绘图数据 ####
H_normalized_plot_df <- readxl::read_xlsx(wh_matrix_file, sheet = "H_normalized")
cluster_assignment_plot_df <- readxl::read_xlsx(wh_matrix_file, sheet = "Pattern group assignment")

H_normalized_plot_matrix <- H_normalized_plot_df %>%
    as.data.frame(check.names = FALSE) %>%
    tibble::column_to_rownames("cancer_type") %>%
    as.matrix()

annotation_col_df <- cluster_assignment_plot_df %>%
    dplyr::mutate(`Cancer group` = factor(`Pattern group`, levels = pattern_names)) %>%
    dplyr::select(-`Pattern group`) %>%
    as.data.frame(stringsAsFactors = FALSE)

rownames(annotation_col_df) <- annotation_col_df$cancer_type
annotation_col_df$cancer_type <- NULL
annotation_col_df <- annotation_col_df[colnames(t(H_normalized_plot_matrix)), , drop = FALSE]

cluster_levels <- levels(annotation_col_df$`Cancer group`)
cluster_colors <- c("#DF7A5E", "#A15F88", "#82B29A", "#F2CC8E", "#BC4749", "#A15F88", "#4C956C", "#6C91BF", "#BC4749")
cluster_colors <- cluster_colors[seq_along(cluster_levels)]

annotation_colors <- list(
    "Cancer group" = stats::setNames(cluster_colors, cluster_levels)
)

load(h_cluster_file)

## 2. 绘图 ####
h_heatmap_width_cm <- 14.65
h_heatmap_height_cm <- 7

grDevices::cairo_pdf(
    filename = h_heatmap_file,
    width = h_heatmap_width_cm / 2.54,
    height = h_heatmap_height_cm / 2.54,
    pointsize = 7
)

h_top_annotation <- ComplexHeatmap::HeatmapAnnotation(
    df = annotation_col_df,
    col = annotation_colors,
    show_legend = FALSE,
    show_annotation_name = TRUE,
    annotation_name_side = "right",
    annotation_name_rot = 0,
    annotation_name_gp = grid::gpar(fontsize = 7, fontface = "bold", col = "black"),
    simple_anno_size = grid::unit(3, "mm")
)

H_heatmap <- ComplexHeatmap::Heatmap(
    t(H_normalized_plot_matrix),
    name = "Pattern score",
    row_title = "irAE reporting pattern",
    row_title_gp = grid::gpar(fontsize = 7, fontface = "bold", col = "black"),
    cluster_rows = FALSE,
    cluster_columns = stats::as.dendrogram(hc),
    top_annotation = h_top_annotation,
    show_row_names = TRUE,
    show_column_names = TRUE,
    col = rev(colorRampPalette(c("#14161E", "#4A1A4C", "#E22C56", "#FF5C3D", "#FFA079", "#FFE4D4"))(50)),
    column_dend_height = grid::unit(9, "mm"),
    row_names_side = "right",
    row_names_gp = grid::gpar(fontsize = 7, col = "black"),
    column_labels = colnames(t(H_normalized_plot_matrix)),
    column_names_rot = 45,
    column_names_gp = grid::gpar(fontsize = 7, col = "black"),
    column_names_max_height = grid::unit(2.3, "cm"),
    width = grid::unit(11.3, "cm"),
    height = grid::unit(2.0, "cm"),
    rect_gp = grid::gpar(col = "grey60", lwd = 0.5),
    heatmap_legend_param = list(
        title = "Pattern score",
        direction = "horizontal",
        title_position = "leftcenter",
        legend_width = grid::unit(4, "cm"),
        grid_height = grid::unit(1.4, "mm"),
        labels_gp = grid::gpar(fontsize = 7, col = "black"),
        title_gp = grid::gpar(fontsize = 7, col = "black")
    )
)

pattern_group_legend <- ComplexHeatmap::Legend(
    title = "Cancer group",
    at = pattern_names,
    legend_gp = grid::gpar(fill = unname(annotation_colors[["Cancer group"]][pattern_names])),
    direction = "horizontal",
    nrow = 1,
    title_position = "leftcenter",
    labels_gp = grid::gpar(fontsize = 7, col = "black"),
    title_gp = grid::gpar(fontsize = 7, fontface = "bold", col = "black")
)

draw(
    H_heatmap,
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom",
    annotation_legend_list = list(pattern_group_legend)
)

dev.off()

message("已生成结果文件:")
message(" - ", rank_estimate_file)
message(" - ", nmf_result_file)
message(" - ", best_fit_file)
message(" - ", wh_matrix_file)
message(" - ", h_cluster_file)
message(" - ", rank_plot_file)
message(" - ", consensus_map_file)
message(" - ", w_heatmap_file)
message(" - ", h_tree_file)
message(" - ", h_heatmap_file)
message("========== 癌种NMF分群分析（rank=5最终版）结束 ==========")
