# 6.2.15 差分相关高次数置换实验 ####
# 固定完成 1,000,000 次器官标签置换；各检查点使用同一序列的前缀。
# 复用 7.7 的输入、H 分数及原始 1,000 次种子，不运行原始分析脚本。
# 全部有效基因参与每模式 BH；无有效相关值的条目保持 NA。
# 结果写入独立运行目录，不覆盖 7.7 结果；不根据显著性提前停止。

# 1. 参数与依赖 ####
script_name <- "6.2.15差分相关高次数置换实验"
checkpoints <- c(1000L, 10000L, 100000L, 1000000L)
workers <- 8L
batch_size <- 100L
original_seed <- 20260607L  # 原 7.7: analysis_seed + 100 + 1
extension_seed <- 20261015L
pattern_order <- c("Ocular", "Hematologic", "Neurologic", "Systemic-like", "Visceral-like")
validate_only <- "--validate-only" %in% commandArgs(trailingOnly = TRUE)
if (validate_only) checkpoints <- 1000L
for (pkg in c("NMF", "readxl", "openxlsx", "Rcpp")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("缺少已安装依赖: ", pkg)
}
if (Sys.getenv("OPENBLAS_NUM_THREADS") != "1") {
  stop("请在启动 R 前设置 OPENBLAS_NUM_THREADS=1，避免多进程嵌套线程。")
}
input_dir <- file.path("Results", "7.7模式相关基因差分相关识别")
nmf_dir <- file.path("Results", "6.1癌种NMF分群分析_rank5")
run_label <- paste0(if (validate_only) "复现检查_" else "百万次置换_", format(Sys.time(), "%Y%m%d_%H%M%S"))
output_dir <- file.path("Results", script_name, run_label)
if (dir.exists(output_dir)) stop("运行目录已经存在: ", output_dir)
dir.create(output_dir, recursive = TRUE)
started_at <- Sys.time()

# 2. 输入读取辅助函数 ####
## 功能：从 Rdata 读取指定对象，不污染全局环境；输入：文件、对象名；输出：对象 ####
read_object <- function(path, object_name) {
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  if (!exists(object_name, env, inherits = FALSE)) stop("缺少对象: ", object_name)
  get(object_name, env, inherits = FALSE)
}

## 功能：记录时间与阶段；输入：文本；输出：控制台日志 ####
log_message <- function(text) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), text, "\n")
  flush.console()
}

## 功能：将每行秩中心化并归一；输入：有限值矩阵；输出：单位长度秩矩阵，常数行为 NA ####
rank_rows <- function(x) {
  ranks <- t(apply(x, 1, rank, ties.method = "average"))
  centered <- ranks - rowMeans(ranks)
  centered / sqrt(rowSums(centered^2))
}

# 3. 对齐原始输入与已保存的观察统计量 ####
log_message(paste("运行目录:", output_dir))
workbook <- file.path(input_dir, "tables", "模式相关基因差分相关分析结果.xlsx")
saved <- as.data.frame(readxl::read_excel(workbook, sheet = "合并指标全表"))
tiers <- as.data.frame(readxl::read_excel(workbook, sheet = "器官纳入与样本量分层"))
organs <- tiers$primary_organ[tiers$enters_rho_tumor & tiers$enters_rho_normal]
genes <- unique(saved$gene)
stopifnot(length(organs) == 20L, length(genes) == 15388L)
tumor_file <- file.path(input_dir, "rdata", "organ_median_log2_tumor_matrix.Rdata")
normal_file <- file.path(input_dir, "rdata", "organ_median_log2_normal_matrix.Rdata")
best_file <- file.path(nmf_dir, "FAERS_NMF_best_fit.rdata")
tumor <- read_object(tumor_file, "organ_median_log2_tumor_matrix")[genes, organs, drop = FALSE]
normal <- read_object(normal_file, "organ_median_log2_normal_matrix")[genes, organs, drop = FALSE]
h_raw <- as.matrix(NMF::coef(read_object(best_file, "best_fit")))
if (is.null(rownames(h_raw)) || !all(pattern_order %in% rownames(h_raw))) {
  stopifnot(nrow(h_raw) == 5L)
  rownames(h_raw) <- c("Hematologic", "Ocular", "Neurologic", "Systemic-like", "Visceral-like")
}
h <- t(h_raw[pattern_order, organs, drop = FALSE])
stopifnot(all(is.finite(tumor)), all(is.finite(normal)), all(is.finite(h)))
observed <- saved_p <- saved_q <- matrix(NA_real_, length(genes), 5L, dimnames = list(genes, pattern_order))
for (k in seq_along(pattern_order)) {
  d <- saved[saved$pattern == pattern_order[k], ]
  stopifnot(!anyDuplicated(d$gene), nrow(d) == length(genes))
  d <- d[match(genes, d$gene), ]
  observed[, k] <- d$delta_rho
  saved_p[, k] <- d$delta_perm_p_value
  saved_q[, k] <- d$delta_perm_fdr
}
tumor_rank <- rank_rows(tumor)
normal_rank <- rank_rows(normal)
h_rank <- t(rank_rows(t(h)))
delta_rank <- tumor_rank - normal_rank
observed_fast <- delta_rank %*% h_rank
stopifnot(identical(is.finite(observed_fast), is.finite(observed)))
observed_error <- max(abs(observed_fast - observed), na.rm = TRUE)
stopifnot(observed_error < 1e-12)
# Excel 导出会截断末位；>= 边界必须使用原始 cor 重新计算的双精度观察值。
observed_reference <- suppressWarnings(
  stats::cor(t(tumor), h, method = "spearman", use = "pairwise.complete.obs") -
    stats::cor(t(normal), h, method = "spearman", use = "pairwise.complete.obs"))
stopifnot(max(abs(observed_reference - observed), na.rm = TRUE) < 1e-12)
observed <- observed_reference
valid <- is.finite(observed)
valid_rows <- which(rowSums(valid) == ncol(valid))
stopifnot(all(rowSums(valid) %in% c(0L, 5L)))
excluded <- data.frame(基因 = genes[-valid_rows], 原因 = "原观察相关或差分相关缺失，置换P与BH均保持NA")
log_message(paste("有效基因:", length(valid_rows), "排除:", nrow(excluded), "观察差值最大误差:", observed_error))
tumor_valid <- tumor[valid_rows, , drop = FALSE]
normal_valid <- normal[valid_rows, , drop = FALSE]
tumor_raw_ranks <- t(apply(tumor_valid, 1, rank, ties.method = "average"))
normal_raw_ranks <- t(apply(normal_valid, 1, rank, ties.method = "average"))
h_raw_ranks <- apply(h, 2, rank, ties.method = "average")
delta_rank_valid <- delta_rank[valid_rows, , drop = FALSE]
observed_valid <- observed[valid_rows, , drop = FALSE]
abs_observed <- abs(observed_valid)

# 4. 固定随机序列与输入记录 ####
set.seed(original_seed)
original_iteration_seeds <- sample.int(.Machine$integer.max, 1000L)
set.seed(extension_seed)
extension_iteration_seeds <- sample.int(.Machine$integer.max, max(checkpoints) - 1000L)
iteration_seeds <- c(original_iteration_seeds, extension_iteration_seeds)
saveRDS(list(seeds = iteration_seeds, RNGkind = RNGkind(), checkpoints = checkpoints,
             organs = organs, genes = genes, patterns = pattern_order),
        file.path(output_dir, "固定置换序列与对齐信息.rds"))
input_files <- c(workbook, tumor_file, normal_file, best_file,
                 file.path("Code", paste0(script_name, ".R")))
input_manifest <- data.frame(文件 = input_files, MD5 = unname(tools::md5sum(input_files)))
openxlsx::write.xlsx(list("输入文件" = input_manifest, "排除基因" = excluded),
                    file.path(output_dir, "输入与排除记录.xlsx"), overwrite = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "运行环境.txt"))

# 5. 分块置换计数 ####
## 功能：单次扫描矩阵累计极端次数并列出浮点边界；输入：置换差值、观察绝对差；输出：计数与边界行 ####
## 编译循环减少大型临时逻辑矩阵的复制，不改变统计量或比较规则。
Rcpp::cppFunction(code = '
Rcpp::List scan_extremes(Rcpp::NumericMatrix delta, Rcpp::NumericMatrix observed_abs) {
  int ng = delta.nrow(), nc = delta.ncol(), nk = observed_abs.ncol();
  Rcpp::NumericMatrix counts(ng, nk);
  Rcpp::List boundaries(nc);
  for (int col = 0; col < nc; ++col) {
    int k = col % nk;
    std::vector<int> rows;
    for (int g = 0; g < ng; ++g) {
      double value = std::fabs(delta(g, col));
      double threshold = observed_abs(g, k);
      if (value >= threshold) counts(g, k) += 1;
      if (threshold > 0 && std::fabs(value - threshold) < 1e-12) rows.push_back(g + 1);
    }
    boundaries[col] = Rcpp::wrap(rows);
  }
  return Rcpp::List::create(Rcpp::Named("counts") = counts, Rcpp::Named("boundaries") = boundaries);
}')

## 功能：计算一块置换的极端次数；输入：迭代索引；输出：每个有效基因和模式的计数及边界复核次数 ####
## Spearman 等价于标准化秩向量内积；预计算秩避免每次重新排序。
## 接近比较边界时回到原始 stats::cor 计算，保持原 >= 比较的浮点行为。
count_batch <- function(indices) {
  permutations <- lapply(iteration_seeds[indices], function(seed) {
    set.seed(seed)
    sample(seq_along(organs))
  })
  permuted_h <- do.call(cbind, lapply(permutations, function(ix) h_rank[ix, , drop = FALSE]))
  delta <- delta_rank_valid %*% permuted_h
  scanned <- scan_extremes(delta, abs_observed)
  counts <- scanned$counts
  boundary_n <- 0L
  for (column in seq_len(ncol(delta))) {
    rows <- scanned$boundaries[[column]]
    if (length(rows)) {
        k <- (column - 1L) %% 5L + 1L
        j <- (column - 1L) %/% 5L + 1L
        hp <- h_raw_ranks[permutations[[j]], k, drop = FALSE]
        # 输入无缺失，预先求秩后的 Pearson 与原 Spearman 相同；仅回算当前模式。
        exact <- stats::cor(t(tumor_raw_ranks[rows, , drop = FALSE]), hp,
                            method = "pearson", use = "all.obs") -
          stats::cor(t(normal_raw_ranks[rows, , drop = FALSE]), hp,
                     method = "pearson", use = "all.obs")
        counts[rows, k] <- counts[rows, k] +
          as.integer(abs(exact[, 1]) >= abs_observed[rows, k]) -
          as.integer(abs(delta[rows, column]) >= abs_observed[rows, k])
        boundary_n <- boundary_n + length(rows)
    }
  }
  list(counts = counts, boundary_n = boundary_n)
}

# 6. 按固定检查点计算并导出 ####
counts <- matrix(0, length(valid_rows), 5L, dimnames = list(genes[valid_rows], pattern_order))
completed <- 0L
summaries <- list()
checkpoint_details <- list()
total_boundary <- 0
for (target in checkpoints) {
  stage_started <- Sys.time()
  starts <- seq(completed + 1L, target, by = batch_size)
  batches <- lapply(starts, function(first) seq.int(first, min(first + batch_size - 1L, target)))
  # 每波至多 8 个任务；累计计数，不保存百万次零分布数组。
  for (wave_start in seq(1L, length(batches), by = workers)) {
    selected <- seq.int(wave_start, min(wave_start + workers - 1L, length(batches)))
    result <- parallel::mclapply(batches[selected], count_batch,
                                mc.cores = min(workers, length(selected)), mc.set.seed = FALSE)
    for (item in result) {
      if (inherits(item, "try-error") || is.null(item$counts)) stop("置换任务失败")
      counts <- counts + item$counts
      total_boundary <- total_boundary + item$boundary_n
    }
    completed <- tail(batches[[tail(selected, 1L)]], 1L)
    if (wave_start == 1L || completed == target || (wave_start - 1L) %% (workers * 10L) == 0L) {
      log_message(paste("已完成", completed, "/", max(checkpoints),
                        "总用时(秒)", round(as.numeric(difftime(Sys.time(), started_at, units = "secs")), 1)))
    }
  }
  p <- (counts + 1) / (target + 1)
  q <- apply(p, 2, p.adjust, method = "BH")
  if (target == 1000L) {
    previous_counts <- round(saved_p[valid_rows, , drop = FALSE] * 1001 - 1)
    discrepancy <- counts - previous_counts
    validation <- data.frame(检查 = c("观察差值最大误差", "原1000次极端计数不一致条目数", "原1000次极端计数最大误差"),
                             数值 = c(observed_error, sum(discrepancy != 0), max(abs(discrepancy))))
    openxlsx::write.xlsx(list("复现检查" = validation, "排除基因" = excluded),
                        file.path(output_dir, "原1000次复现检查.xlsx"), overwrite = FALSE)
    if (any(discrepancy != 0)) {
      saveRDS(list(new_counts = counts, original_counts = previous_counts, difference = discrepancy),
              file.path(output_dir, "复现差异.rds"))
      stop("原1000次置换未完全复现；停止扩展并检查种子、对齐或数值边界。")
    }
    log_message("原1000次的全部有效基因极端计数完全复现。")
  }
  detail <- do.call(rbind, lapply(seq_len(5L), function(k) {
    count <- counts[, k]
    ci_low <- ifelse(count == 0, 0, stats::qbeta(0.025, count, target - count + 1))
    ci_high <- ifelse(count == target, 1, stats::qbeta(0.975, count + 1, target - count))
    data.frame(基因 = genes[valid_rows], 模式 = pattern_order[k], 置换次数 = target,
               差分相关 = observed_valid[, k], 极端次数 = count, 置换P = p[, k], BH校正值 = q[, k],
               尾部概率MC下限 = ci_low, 尾部概率MC上限 = ci_high,
               名义候选 = abs(observed_valid[, k]) >= 0.3 & p[, k] < 0.05,
               BH小于005 = q[, k] < 0.05, 原1000次P = saved_p[valid_rows, k],
               原1000次BH = saved_q[valid_rows, k], check.names = FALSE)
  }))
  checkpoint_details[[as.character(target)]] <- detail
  summary <- do.call(rbind, lapply(seq_len(5L), function(k) {
    significant <- q[, k] < 0.05
    data.frame(模式 = pattern_order[k], 置换次数 = target, 有效基因数 = length(valid_rows),
               最小可达P = 1 / (target + 1), 最小P = min(p[, k]), 最小BH = min(q[, k]),
               零极端次数基因数 = sum(counts[, k] == 0), BH显著基因数 = sum(significant),
               BH显著且效应达标数 = sum(significant & abs(observed_valid[, k]) >= 0.3),
               名义候选基因数 = sum(p[, k] < 0.05 & abs(observed_valid[, k]) >= 0.3))
  }))
  summaries[[as.character(target)]] <- summary
  saveRDS(list(completed = completed, counts = counts, detail = detail,
               boundary_recalculations = total_boundary),
          file.path(output_dir, paste0("检查点_", target, "次.rds")))
  ordered <- detail[order(detail$BH校正值, detail$置换P, -abs(detail$差分相关)), ]
  openxlsx::write.xlsx(list("阶段汇总" = do.call(rbind, summaries), "全部有效基因" = detail,
                           "BH显著基因" = ordered[ordered$BH小于005, ],
                           "每模式前二十" = do.call(rbind, lapply(pattern_order, function(name) head(ordered[ordered$模式 == name, ], 20))),
                           "排除基因" = excluded),
                      file.path(output_dir, paste0("置换实验_", target, "次.xlsx")), overwrite = FALSE)
  log_message(paste("检查点", target, "次完成；BH显著条目数:", sum(q < 0.05),
                    "阶段秒数:", round(as.numeric(difftime(Sys.time(), stage_started, units = "secs")), 1)))
  print(summary, row.names = FALSE)
}

# 7. 汇总比较与完成标记 ####
comparison <- data.frame(基因 = rep(genes[valid_rows], 5L), 模式 = rep(pattern_order, each = length(valid_rows)),
                         差分相关 = as.vector(observed_valid), check.names = FALSE)
for (target in checkpoints) {
  d <- checkpoint_details[[as.character(target)]]
  comparison[[paste0("P_", target, "次")]] <- d$置换P
  comparison[[paste0("BH_", target, "次")]] <- d$BH校正值
  comparison[[paste0("极端次数_", target, "次")]] <- d$极端次数
}
openxlsx::write.xlsx(list("阶段汇总" = do.call(rbind, summaries), "各次数逐基因比较" = comparison,
                         "排除基因" = excluded, "输入文件" = input_manifest),
                    file.path(output_dir, "高次数置换实验汇总.xlsx"), overwrite = FALSE)
stopifnot(identical(unname(tools::md5sum(input_files)), input_manifest$MD5))
writeLines(c(paste("完成置换次数:", completed), paste("有效基因数:", length(valid_rows)),
             paste("边界比较回算次数:", total_boundary), paste("结束时间:", Sys.time()),
             "Monte Carlo区间为逐检验二项分布区间，不是同时置信区间，也不是生物学效应置信区间。"),
           file.path(output_dir, "完成标记.txt"))
log_message("实验完成，输入文件MD5保持一致。")
