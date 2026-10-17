# 药物不良反应数据读取与清洗####
# 数据来源: FAERS (美国), JADER (日本)
# 可逐行运行查看结果

# 第一部分：环境设置####
# 加载必要的包
library(tidyverse)
library(data.table)
rm(list = ls())
gc()

# 加载通用路径配置（自动设置 data_path, output_path, results_path）
source("Code/config.R")

# 分析时间范围：统一用于各数据库的年份筛选
analysis_year_min <- 2015
analysis_year_max <- 2025

# 起始处理模块：FAERS / JADER / ALL
start_section <- toupper(Sys.getenv("PAN_CANCER_START", unset = "ALL"))

# 并行设置：在专用计算服务器上优先使用全部可用核心
total_parallel_workers <- max(1L, future::availableCores())
future_globals_max_gb <- suppressWarnings(as.numeric(Sys.getenv("PAN_CANCER_FUTURE_GLOBALS_GB", unset = "50")))
if (is.na(future_globals_max_gb) || future_globals_max_gb <= 0) {
  future_globals_max_gb <- 50
}
options(future.globals.maxSize = future_globals_max_gb * 1024^3)
data.table::setDTthreads(total_parallel_workers)
options(mc.cores = total_parallel_workers)

# 并行设置：优先并行独立任务，避免对大对象join强行并行
resolve_parallel_workers <- function(env_var, default_workers) {
  env_value <- Sys.getenv(env_var, unset = "")
  env_workers <- suppressWarnings(as.integer(env_value))

  max_workers <- total_parallel_workers

  if (!is.na(env_workers) && env_workers >= 1L) {
    return(min(env_workers, max_workers))
  }

  min(as.integer(default_workers), max_workers)
}

# 函数：执行受控并行lapply，结束后恢复原有future计划
run_parallel_lapply <- function(X, FUN, workers, worker_dt_threads = NULL, backend = NULL, ...) {
  if (length(X) == 0) {
    return(list())
  }

  workers <- min(as.integer(workers), length(X), total_parallel_workers)

  if (workers <= 1L) {
    data.table::setDTthreads(total_parallel_workers)
    return(lapply(X, FUN, ...))
  }

  if (is.null(worker_dt_threads)) {
    worker_dt_threads <- max(1L, floor(total_parallel_workers / workers))
  }

  if (is.null(backend)) {
    backend <- if (future::supportsMulticore()) future::multicore else future::multisession
  }
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

  future::plan(backend, workers = workers)
  future.apply::future_lapply(
    X,
    function(.x, ...) {
      data.table::setDTthreads(worker_dt_threads)
      FUN(.x, ...)
    },
    ...
  )
}

io_parallel_workers <- resolve_parallel_workers("PAN_CANCER_IO_WORKERS", total_parallel_workers)
faers_parallel_workers <- resolve_parallel_workers("PAN_CANCER_FAERS_WORKERS", total_parallel_workers)
jader_parallel_workers <- min(
  io_parallel_workers,
  max(1L, suppressWarnings(as.integer(Sys.getenv("PAN_CANCER_JADER_WORKERS", unset = "8"))))
)
message(sprintf(
  "并行设置: total cores=%d, data.table threads=%d, IO workers=%d, FAERS workers=%d",
  total_parallel_workers,
  data.table::getDTthreads(),
  io_parallel_workers,
  faers_parallel_workers
))
message(sprintf("future.globals.maxSize=%.0fGB", future_globals_max_gb))

# 查看数据文件夹结构
list.dirs(.data_path, recursive = FALSE)

# 第二部分：FAERS数据处理####
## 2.1 查看FAERS文件列表####
faers_path <- file.path(.data_path, "FAERS")

# 函数：从FAERS压缩包文件名中提取季度，并统一为YYYYQn格式
extract_faers_quarter <- function(zip_file) {
  zip_name <- basename(zip_file)
  quarter <- sub(
    "^faers_ascii_([0-9]{4}[Qq][1-4])\\.zip$",
    "\\1",
    zip_name,
    ignore.case = TRUE
  )

  if (identical(quarter, zip_name)) {
    return(NA_character_)
  }

  toupper(quarter)
}

# 函数：读取FAERS压缩包列表，并按分析年份范围与季度顺序排序
list_faers_zip_files <- function(faers_path) {
  faers_zips <- list.files(
    faers_path,
    pattern = "^faers_ascii_[0-9]{4}[Qq][1-4]\\.zip$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )

  if (length(faers_zips) == 0) {
    return(character(0))
  }

  quarter_labels <- vapply(faers_zips, extract_faers_quarter, character(1))
  valid_idx <- !is.na(quarter_labels)
  faers_zips <- faers_zips[valid_idx]
  quarter_labels <- quarter_labels[valid_idx]

  report_years <- as.integer(substr(quarter_labels, 1, 4))
  year_idx <- report_years >= analysis_year_min & report_years <= analysis_year_max
  faers_zips <- faers_zips[year_idx]
  quarter_labels <- quarter_labels[year_idx]

  order_idx <- order(
    as.integer(substr(quarter_labels, 1, 4)),
    as.integer(substr(quarter_labels, 6, 6))
  )

  faers_zips[order_idx]
}

# 函数：处理单个FAERS压缩包并返回指定表
process_faers_zip <- function(zip_file, table_pattern) {
  quarter <- extract_faers_quarter(zip_file)

  if (is.na(quarter)) {
    return(NULL)
  }

  tryCatch(
    {
      tmp_extract <- tempfile(pattern = paste0("faers_", table_pattern, "_"), tmpdir = tempdir())
      dir.create(tmp_extract, recursive = TRUE, showWarnings = FALSE)
      on.exit(unlink(tmp_extract, recursive = TRUE, force = TRUE), add = TRUE)

      unzip(zip_file, exdir = tmp_extract, overwrite = TRUE)
      txt_files <- list.files(tmp_extract, pattern = "\\.txt$", recursive = TRUE, full.names = TRUE)

      target_f <- txt_files[grepl(paste0("^", table_pattern), basename(txt_files), ignore.case = TRUE)]
      target_f <- target_f[!grepl("DELETE", basename(target_f), ignore.case = TRUE)]

      if (length(target_f) == 0) {
        return(NULL)
      }

      dt <- data.table::fread(target_f[1], sep = "$", quote = "", fill = TRUE, encoding = "Latin-1")
      dt$quarter <- quarter

      list(quarter = quarter, data = dt)
    },
    error = function(e) {
      message(sprintf("读取%s季度%s表失败: %s", quarter, table_pattern, e$message))
      NULL
    }
  )
}

## 2.2 内存优化：逐表处理FAERS数据####
# 策略：每次只处理一个表，完成后释放内存，避免同时加载7个表

# 函数：并行读取FAERS压缩包中的指定表
read_faers_table <- function(faers_zips, table_pattern, workers = faers_parallel_workers) {
  if (length(faers_zips) == 0) {
    return(list())
  }

  message(sprintf("开始并行读取FAERS %s表: %d个季度, workers=%d", table_pattern, length(faers_zips), workers))
  zip_results <- run_parallel_lapply(
    faers_zips,
    function(zip_file) process_faers_zip(zip_file, table_pattern),
    workers = workers
  )

  zip_results <- Filter(function(x) !is.null(x), zip_results)
  data_list <- lapply(zip_results, `[[`, "data")
  names(data_list) <- vapply(zip_results, `[[`, character(1), "quarter")

  message(sprintf("完成读取FAERS %s表: 成功%d个季度", table_pattern, length(data_list)))
  data_list
}


# 检查各季度变量是否一致
check_column_consistency <- function(data_list, table_name) {
  if (length(data_list) == 0) {
    return(NULL)
  }

  # 获取每个季度的列名
  col_list <- lapply(data_list, names)
  quarters <- names(data_list)

  # 找出所有唯一列名
  all_cols <- unique(unlist(col_list))

  # 创建列存在矩阵
  col_matrix <- sapply(col_list, function(x) all_cols %in% x)
  rownames(col_matrix) <- all_cols

  # 检查是否所有季度都有相同的列
  consistent <- all(rowSums(col_matrix) == length(data_list))

  if (!consistent) {
    message(paste0("\n[警告] ", table_name, " 表各季度变量不一致:"))
    # 显示不一致的列
    inconsistent_cols <- all_cols[rowSums(col_matrix) != length(data_list)]
    for (col in inconsistent_cols) {
      present_in <- quarters[col_matrix[col, ]]
      missing_in <- quarters[!col_matrix[col, ]]
      message(paste0("  - ", col, ": 存在于 ", length(present_in), " 个季度, 缺失于 ", length(missing_in), " 个季度"))
    }
  } else {
    message(paste0("[OK] ", table_name, " 表所有季度变量一致 (", length(all_cols), " 列)"))
  }

  return(list(consistent = consistent, all_cols = all_cols, col_matrix = col_matrix))
}

# 获取所有季度的完整列名集合
get_all_columns <- function(data_list) {
  unique(unlist(lapply(data_list, names)))
}

# 标准化单个数据表 (原地修改，减少内存复制)
standardize_table <- function(dt, standard_cols, rename_map = NULL) {
  # future worker 反序列化后，data.table 可能丢失列预分配信息
  # 在按引用新增列前显式恢复data.table状态
  data.table::setDT(dt)
  data.table::setalloccol(dt)

  # 1. 重命名列
  if (!is.null(rename_map)) {
    for (old_name in names(rename_map)) {
      new_name <- rename_map[[old_name]]
      if (old_name %in% names(dt)) {
        if (new_name %in% names(dt)) {
          dt[, (new_name) := fifelse(is.na(get(new_name)) | get(new_name) == "", get(old_name), get(new_name))]
          dt[, (old_name) := NULL]
        } else {
          setnames(dt, old_name, new_name)
        }
      }
    }
  }
  # 2. 转换为字符（注意：integer64和大数值numeric需要特殊处理避免精度丢失）
  for (col in names(dt)) {
    col_val <- dt[[col]]
    if (inherits(col_val, "integer64")) {
      # integer64 使用 bit64 专用转换，避免精度丢失
      set(dt, j = col, value = bit64::as.character.integer64(col_val))
    } else if (is.numeric(col_val)) {
      # numeric 使用 format 避免科学计数法
      set(dt, j = col, value = format(col_val, scientific = FALSE, trim = TRUE))
    } else {
      set(dt, j = col, value = as.character(col_val))
    }
  }
  # 3. 添加缺失列
  for (col in standard_cols) {
    if (!col %in% names(dt)) {
      set(dt, j = col, value = NA_character_)
    }
  }
  # 4. 删除多余列
  extra_cols <- setdiff(names(dt), standard_cols)
  if (length(extra_cols) > 0) dt[, (extra_cols) := NULL]
  # 5. 统一顺序
  setcolorder(dt, standard_cols)
  invisible(dt)
}

# 函数：并行标准化各季度数据表
standardize_table_list <- function(data_list, standard_cols, rename_map = NULL, workers = faers_parallel_workers) {
  if (length(data_list) == 0) {
    return(data_list)
  }

  quarter_names <- names(data_list)
  standardized_list <- run_parallel_lapply(
    data_list,
    function(dt) {
      standardize_table(dt, standard_cols, rename_map)
      dt
    },
    workers = workers
  )
  names(standardized_list) <- quarter_names
  standardized_list
}

# 函数：按列名稳健拼接各季度表，兼容季度间变量漂移
bind_standardized_tables <- function(data_list, table_name) {
  if (length(data_list) == 0) {
    return(data.table())
  }

  col_counts <- vapply(data_list, ncol, integer(1))
  if (length(unique(col_counts)) > 1) {
    message(sprintf(
      "[WARN] %s 标准化后列数仍不一致: %s；按列名补齐缺失列后拼接",
      table_name,
      paste(sort(unique(col_counts)), collapse = ", ")
    ))
  }

  rbindlist(data_list, use.names = TRUE, fill = TRUE)
}

# 函数：统一FAERS主键列类型，避免integer64/double/character混用导致筛选失败
normalize_faers_key_cols <- function(dt, key_cols) {
  if (nrow(dt) == 0) {
    return(invisible(dt))
  }

  data.table::setDT(dt)
  data.table::setalloccol(dt)

  existing_key_cols <- intersect(key_cols, names(dt))
  for (col in existing_key_cols) {
    col_val <- dt[[col]]
    if (inherits(col_val, "integer64")) {
      set(dt, j = col, value = bit64::as.character.integer64(col_val))
    } else if (is.numeric(col_val)) {
      set(dt, j = col, value = format(col_val, scientific = FALSE, trim = TRUE))
    } else {
      set(dt, j = col, value = as.character(col_val))
    }
  }

  invisible(dt)
}

## 2.3 逐表处理FAERS数据####
if (start_section %in% c("ALL", "FAERS")) {
  faers_zips <- list_faers_zip_files(faers_path)
  length(faers_zips) # 查看有多少个zip文件
  basename(faers_zips) # 查看文件名
  ### DEMO表####
  message("\n========== 处理 DEMO 表 ==========")
  all_faers_demo <- read_faers_table(faers_zips, "DEMO")
  check_column_consistency(all_faers_demo, "DEMO")

  demo_rename <- list("gndr_cod" = "sex")
  demo_all_cols <- get_all_columns(all_faers_demo)
  demo_all_cols <- setdiff(demo_all_cols, "gndr_cod")
  if (!"sex" %in% demo_all_cols) demo_all_cols <- c(demo_all_cols, "sex")

  all_faers_demo <- standardize_table_list(all_faers_demo, demo_all_cols, demo_rename)
  faers_demo_all <- bind_standardized_tables(all_faers_demo, "DEMO")
  normalize_faers_key_cols(faers_demo_all, c("primaryid", "caseid", "caseversion"))

  # ===== 时间筛选：只保留2015-2025的报告 =====
  message(paste("筛选前DEMO记录数:", nrow(faers_demo_all)))
  table(faers_demo_all$quarter)

  # 提取年份（quarter格式如 "2015Q1"）
  faers_demo_all[, report_year := as.numeric(substr(quarter, 1, 4))]
  faers_demo_all <- faers_demo_all[
    report_year >= analysis_year_min & report_year <= analysis_year_max
  ]

  message(paste("筛选后DEMO记录数:", nrow(faers_demo_all)))
  message(paste("时间范围:", min(faers_demo_all$quarter), "至", max(faers_demo_all$quarter)))

  # ===== FDA推荐去重策略（第一步去重） =====
  # 原始记录数
  message(paste("DEMO原始记录数:", nrow(faers_demo_all)))
  message(paste("DEMO原始病例数(caseid):", uniqueN(faers_demo_all$caseid)))

  # FDA去重：每个caseid只保留caseversion最大的记录
  faers_demo_all[, caseversion_num := as.numeric(caseversion)]
  faers_demo_dedup <- faers_demo_all[faers_demo_all[, .I[which.max(caseversion_num)], by = caseid]$V1]

  message(paste("caseid去重后记录数:", nrow(faers_demo_dedup)))

  # 检查primaryid是否有重复
  n_dup_primaryid <- sum(duplicated(faers_demo_dedup$primaryid))
  message(paste("primaryid重复数:", n_dup_primaryid))

  # 如果primaryid有重复，进一步去重（保留每个primaryid的最新版本）
  if (n_dup_primaryid > 0) {
    faers_demo_dedup <- faers_demo_dedup[faers_demo_dedup[, .I[which.max(caseversion_num)], by = primaryid]$V1]
    message(paste("primaryid去重后记录数:", nrow(faers_demo_dedup)))
  }

  message(paste("总去重比例:", round((1 - nrow(faers_demo_dedup) / nrow(faers_demo_all)) * 100, 2), "%"))

  # 保存去重后的primaryid列表，用于筛选其他表
  valid_primaryids <- as.character(unique(faers_demo_dedup$primaryid))
  message(paste("有效primaryid数:", length(valid_primaryids)))
  save(valid_primaryids, file = file.path(.output_path, "valid_primaryids.RData"))
  # 释放原始数据
  rm(faers_demo_all)
  gc()

  # 清洗DEMO
  faers_demo <- faers_demo_dedup %>%
    as_tibble() %>%
    mutate(
      across(where(is.character), ~ na_if(trimws(.), "")),
      sex_std = case_when(
        toupper(sex) %in% c("M", "MALE") ~ "Male",
        toupper(sex) %in% c("F", "FEMALE") ~ "Female",
        TRUE ~ "Unknown"
      ),
      age_numeric = suppressWarnings(as.numeric(age)),
      age_unit_std = toupper(age_cod),
      age_years = case_when(
        age_unit_std == "YR" ~ age_numeric,
        age_unit_std == "MON" ~ age_numeric / 12,
        age_unit_std == "WK" ~ age_numeric / 52,
        age_unit_std == "DY" ~ age_numeric / 365,
        age_unit_std == "DEC" ~ age_numeric * 10,
        TRUE ~ age_numeric
      ),
      event_date = as.Date(as.character(event_dt), format = "%Y%m%d"),
      mfr_date = as.Date(as.character(mfr_dt), format = "%Y%m%d"),
      init_fda_date = as.Date(as.character(init_fda_dt), format = "%Y%m%d"),
      fda_date = as.Date(as.character(fda_dt), format = "%Y%m%d"),
      report_delay = as.numeric(init_fda_date - event_date),
      source = "FAERS"
    ) %>%
    select(primaryid, caseid, caseversion, event_dt, age_years, sex_std, wt, reporter_country, occp_cod, quarter)
  dim(faers_demo)
  save(faers_demo, file = file.path(.output_path, "faers_demo.RData"))
  # 释放内存
  rm(all_faers_demo, faers_demo_dedup, faers_demo)
  gc()

  # --- DRUG表 ---
  message("\n========== 处理 DRUG 表 ==========")
  all_faers_drug <- read_faers_table(faers_zips, "DRUG")
  check_column_consistency(all_faers_drug, "DRUG")
  drug_cols <- get_all_columns(all_faers_drug)
  all_faers_drug <- standardize_table_list(all_faers_drug, drug_cols)
  faers_drug_all <- bind_standardized_tables(all_faers_drug, "DRUG")
  normalize_faers_key_cols(faers_drug_all, c("primaryid", "caseid", "drug_seq"))
  message(paste("DRUG原始记录数:", nrow(faers_drug_all)))
  # 用valid_primaryids筛选
  faers_drug_all <- faers_drug_all[primaryid %in% valid_primaryids]
  message(paste("DRUG筛选后记录数:", nrow(faers_drug_all)))

  faers_drug <- faers_drug_all %>%
    as_tibble() %>%
    mutate(
      across(where(is.character), ~ na_if(trimws(.), "")),
      drugname_std = toupper(trimws(drugname)),
      role_cod_std = toupper(trimws(role_cod)),
      source = "FAERS"
    ) %>%
    select(primaryid, caseid, drug_seq, drugname_std, prod_ai, role_cod_std, quarter)
  dim(faers_drug)
  save(faers_drug, file = file.path(.output_path, "faers_drug.RData"))
  # 释放内存
  rm(all_faers_drug, faers_drug_all, faers_drug)
  gc()

  # --- REAC表 ---
  message("\n========== 处理 REAC 表 ==========")
  all_faers_reac <- read_faers_table(faers_zips, "REAC")
  check_column_consistency(all_faers_reac, "REAC")
  reac_cols <- get_all_columns(all_faers_reac)
  all_faers_reac <- standardize_table_list(all_faers_reac, reac_cols)
  faers_reac_all <- bind_standardized_tables(all_faers_reac, "REAC")
  normalize_faers_key_cols(faers_reac_all, c("primaryid", "caseid"))
  message(paste("REAC原始记录数:", nrow(faers_reac_all)))
  faers_reac_all <- faers_reac_all[primaryid %in% valid_primaryids]
  message(paste("REAC筛选后记录数:", nrow(faers_reac_all)))

  faers_reac <- faers_reac_all %>%
    as_tibble() %>%
    mutate(
      across(where(is.character), ~ na_if(trimws(.), "")),
      pt_std = toupper(trimws(pt)),
      source = "FAERS"
    ) %>%
    select(primaryid, caseid, pt_std)
  dim(faers_reac)
  save(faers_reac, file = file.path(.output_path, "faers_reac.RData"))
  # 释放内存
  rm(all_faers_reac, faers_reac_all, faers_reac)
  gc()

  # --- INDI表 ---
  message("\n========== 处理 INDI 表 ==========")
  all_faers_indi <- read_faers_table(faers_zips, "INDI")
  check_column_consistency(all_faers_indi, "INDI")
  indi_cols <- get_all_columns(all_faers_indi)
  all_faers_indi <- standardize_table_list(all_faers_indi, indi_cols)
  faers_indi_all <- bind_standardized_tables(all_faers_indi, "INDI")
  normalize_faers_key_cols(faers_indi_all, c("primaryid", "caseid"))
  message(paste("INDI原始记录数:", nrow(faers_indi_all)))
  faers_indi_all <- faers_indi_all[primaryid %in% valid_primaryids]
  message(paste("INDI筛选后记录数:", nrow(faers_indi_all)))

  faers_indi <- faers_indi_all %>%
    as_tibble() %>%
    mutate(
      across(where(is.character), ~ na_if(trimws(.), "")),
      indi_pt_std = toupper(trimws(indi_pt)),
      source = "FAERS"
    ) %>%
    select(primaryid, caseid, indi_pt_std)
  dim(faers_indi)
  save(faers_indi, file = file.path(.output_path, "faers_indi.RData"))
  # 释放内存
  rm(all_faers_indi, faers_indi_all, faers_indi)
  gc()

  # --- OUTC表 ---
  message("\n========== 处理 OUTC 表 ==========")
  all_faers_outc <- read_faers_table(faers_zips, "OUTC")
  check_column_consistency(all_faers_outc, "OUTC")
  outc_cols <- get_all_columns(all_faers_outc)
  all_faers_outc <- standardize_table_list(all_faers_outc, outc_cols)
  faers_outc_all <- bind_standardized_tables(all_faers_outc, "OUTC")
  normalize_faers_key_cols(faers_outc_all, c("primaryid", "caseid"))
  message(paste("OUTC原始记录数:", nrow(faers_outc_all)))
  faers_outc_all <- faers_outc_all[primaryid %in% valid_primaryids]
  message(paste("OUTC筛选后记录数:", nrow(faers_outc_all)))

  faers_outc <- faers_outc_all %>%
    as_tibble() %>%
    mutate(
      across(where(is.character), ~ na_if(trimws(.), "")),
      source = "FAERS"
    ) %>%
    select(primaryid, caseid, outc_cod)
  dim(faers_outc)
  save(faers_outc, file = file.path(.output_path, "faers_outc.RData"))
  # 释放内存
  rm(all_faers_outc, faers_outc_all, faers_outc)
  gc()

  # --- RPSR表 ---
  message("\n========== 处理 RPSR 表 ==========")
  all_faers_rpsr <- read_faers_table(faers_zips, "RPSR")
  check_column_consistency(all_faers_rpsr, "RPSR")
  rpsr_cols <- get_all_columns(all_faers_rpsr)
  all_faers_rpsr <- standardize_table_list(all_faers_rpsr, rpsr_cols)
  faers_rpsr_all <- bind_standardized_tables(all_faers_rpsr, "RPSR")
  normalize_faers_key_cols(faers_rpsr_all, c("primaryid", "caseid"))
  message(paste("RPSR原始记录数:", nrow(faers_rpsr_all)))
  faers_rpsr_all <- faers_rpsr_all[primaryid %in% valid_primaryids]
  message(paste("RPSR筛选后记录数:", nrow(faers_rpsr_all)))

  faers_rpsr <- faers_rpsr_all %>%
    as_tibble() %>%
    mutate(
      across(where(is.character), ~ na_if(trimws(.), "")),
      source = "FAERS"
    ) %>%
    select(primaryid, caseid, rpsr_cod)
  dim(faers_rpsr)
  save(faers_rpsr, file = file.path(.output_path, "faers_rpsr.RData"))
  # 释放内存
  rm(all_faers_rpsr, faers_rpsr_all, faers_rpsr)
  gc()

  # --- THER表 ---
  message("\n========== 处理 THER 表 ==========")
  all_faers_ther <- read_faers_table(faers_zips, "THER")
  check_column_consistency(all_faers_ther, "THER")
  ther_cols <- get_all_columns(all_faers_ther)
  all_faers_ther <- standardize_table_list(all_faers_ther, ther_cols)
  faers_ther_all <- bind_standardized_tables(all_faers_ther, "THER")
  normalize_faers_key_cols(faers_ther_all, c("primaryid", "caseid", "dsg_drug_seq"))
  message(paste("THER原始记录数:", nrow(faers_ther_all)))
  faers_ther_all <- faers_ther_all[primaryid %in% valid_primaryids]
  message(paste("THER筛选后记录数:", nrow(faers_ther_all)))

  faers_ther <- faers_ther_all %>%
    as_tibble() %>%
    mutate(
      across(where(is.character), ~ na_if(trimws(.), "")),
      start_dt = as.Date(as.character(start_dt), format = "%Y%m%d"),
      end_dt = as.Date(as.character(end_dt), format = "%Y%m%d"),
      source = "FAERS"
    ) %>%
    select(primaryid, caseid, dsg_drug_seq, start_dt, end_dt, dur, dur_cod)
  dim(faers_ther)
  save(faers_ther, file = file.path(.output_path, "faers_ther.RData"))
  # 释放内存
  rm(all_faers_ther, faers_ther_all, faers_ther)
  gc()

  # 第二步去重：临床特征去重####
  # 年龄、性别、体重、报告国家、事件发生时间、用药、不良反应、适应证相同 → 同一病例，保留最新版本

  message("\n========== 第二步去重：临床特征去重 ==========")

  # 加载各表
  load(file.path(.output_path, "faers_demo.RData"))
  load(file.path(.output_path, "faers_drug.RData"))
  load(file.path(.output_path, "faers_reac.RData"))
  load(file.path(.output_path, "faers_indi.RData"))
  dim(faers_demo)
  dim(faers_drug)
  dim(faers_reac)
  dim(faers_indi)
  # 1. 从DEMO提取关键特征
  demo_key <- faers_demo %>%
    select(primaryid, caseid, caseversion, age_years, sex_std, wt, reporter_country, event_dt) %>%
    mutate(caseversion_num = as.numeric(caseversion))

  # 2. 从DRUG提取用药（合并同一报告的所有药物为一个字符串）
  drug_key <- faers_drug %>%
    group_by(primaryid) %>%
    summarise(drugs = paste(sort(unique(na.omit(drugname_std))), collapse = "|"), .groups = "drop")

  # 3. 从REAC提取不良反应（合并同一报告的所有PT为一个字符串）
  reac_key <- faers_reac %>%
    group_by(primaryid) %>%
    summarise(reacs = paste(sort(unique(na.omit(pt_std))), collapse = "|"), .groups = "drop")

  # 4. 从INDI提取适应证（合并同一报告的所有适应证为一个字符串）
  indi_key <- faers_indi %>%
    group_by(primaryid) %>%
    summarise(indis = paste(sort(unique(na.omit(indi_pt_std))), collapse = "|"), .groups = "drop")

  # 5. 合并所有关键字段
  dedup_df <- demo_key %>%
    left_join(drug_key, by = "primaryid") %>%
    left_join(reac_key, by = "primaryid") %>%
    left_join(indi_key, by = "primaryid") %>%
    mutate(
      # 创建去重键：年龄+性别+体重+报告国家+事件日期+用药+不良反应+适应证
      dedup_key = paste(age_years, sex_std, wt, reporter_country, event_dt, drugs, reacs, indis, sep = "||")
    )

  message(paste("去重前报告数:", nrow(dedup_df)))
  message(paste("去重键数:", uniqueN(dedup_df$dedup_key)))

  # 6. 每个去重键只保留caseversion最大的记录
  dedup_df <- as.data.table(dedup_df)
  dedup_final <- dedup_df[dedup_df[, .I[which.max(caseversion_num)], by = dedup_key]$V1]

  message(paste("去重后报告数:", nrow(dedup_final)))
  message(paste("第二步去重比例:", round((1 - nrow(dedup_final) / nrow(dedup_df)) * 100, 2), "%"))

  # 7. 更新valid_primaryids
  valid_primaryids_v2 <- unique(dedup_final$primaryid)
  message(paste("最终有效primaryid数:", length(valid_primaryids_v2)))
  save(valid_primaryids_v2, file = file.path(.output_path, "valid_primaryids_v2.RData"))

  # 8. 找出四表共有的primaryid
  demo_ids <- unique(faers_demo$primaryid)
  drug_ids <- unique(faers_drug$primaryid)
  reac_ids <- unique(faers_reac$primaryid)
  indi_ids <- unique(faers_indi$primaryid)

  # 取交集
  common_ids <- Reduce(intersect, list(valid_primaryids_v2, demo_ids, drug_ids, reac_ids, indi_ids))
  message(paste("DEMO主键数:", length(demo_ids)))
  message(paste("DRUG主键数:", length(drug_ids)))
  message(paste("REAC主键数:", length(reac_ids)))
  message(paste("INDI主键数:", length(indi_ids)))
  message(paste("四表共有主键数:", length(common_ids)))

  # 9. 用共有主键筛选并保存各表
  faers_demo <- faers_demo %>% filter(primaryid %in% common_ids)
  faers_drug <- faers_drug %>% filter(primaryid %in% common_ids)
  faers_reac <- faers_reac %>% filter(primaryid %in% common_ids)
  faers_indi <- faers_indi %>% filter(primaryid %in% common_ids)

  save(faers_demo, file = file.path(.output_path, "faers_demo.RData"))
  save(faers_drug, file = file.path(.output_path, "faers_drug.RData"))
  save(faers_reac, file = file.path(.output_path, "faers_reac.RData"))
  save(faers_indi, file = file.path(.output_path, "faers_indi.RData"))

  # 加载并筛选OUTC, RPSR, THER
  load(file.path(.output_path, "faers_outc.RData"))
  load(file.path(.output_path, "faers_rpsr.RData"))
  load(file.path(.output_path, "faers_ther.RData"))

  faers_outc <- faers_outc %>% filter(primaryid %in% common_ids)
  faers_rpsr <- faers_rpsr %>% filter(primaryid %in% common_ids)
  faers_ther <- faers_ther %>% filter(primaryid %in% common_ids)

  save(faers_outc, file = file.path(.output_path, "faers_outc.RData"))
  save(faers_rpsr, file = file.path(.output_path, "faers_rpsr.RData"))
  save(faers_ther, file = file.path(.output_path, "faers_ther.RData"))

  message("第二步去重完成，所有表已更新")

  # 第三步：计算TTO（Time to Onset）####
  # TTO = event_dt (demo) - start_dt (ther)，精确到天
  # 如果日期缺失或不精确到天，TTO记为"Unknown"

  message("\n========== 计算TTO (Time to Onset) ==========")

  # 检查event_dt是否精确到天（8位数字格式YYYYMMDD）
  check_date_precision <- function(date_str) {
    if (is.na(date_str) || date_str == "") {
      return(FALSE)
    }
    # 必须是8位数字且日部分不为00
    if (nchar(date_str) != 8) {
      return(FALSE)
    }
    day_part <- substr(date_str, 7, 8)
    if (day_part == "00" || day_part == "") {
      return(FALSE)
    }
    return(TRUE)
  }

  # 从demo表获取event_dt
  demo_event <- faers_demo %>%
    select(primaryid, event_dt) %>%
    mutate(
      event_dt_str = as.character(event_dt),
      event_date = as.Date(event_dt_str, format = "%Y%m%d"),
      event_precise = sapply(event_dt_str, check_date_precision)
    )

  # 从ther表获取start_dt（已转换为Date格式，但需要检查原始精度）
  # 重新加载ther原始数据检查精度
  load(file.path(.output_path, "faers_ther.RData"))

  # ther表可能存在同一primaryid + dsg_drug_seq的多条治疗记录
  # 先计算每条治疗记录对应的TTO候选值，再汇总为唯一键，避免join后放大drug表行数
  ther_tto <- faers_ther %>%
    select(primaryid, dsg_drug_seq, start_dt) %>%
    distinct() %>%
    left_join(demo_event %>% select(primaryid, event_date, event_precise), by = "primaryid") %>%
    mutate(
      start_precise = !is.na(start_dt),
      tto_days = as.numeric(event_date - start_dt),
      tto_valid = event_precise & start_precise & !is.na(tto_days) & tto_days >= 0
    ) %>%
    group_by(primaryid, dsg_drug_seq) %>%
    arrange(desc(tto_valid), tto_days, start_dt, .by_group = TRUE) %>%
    summarise(
      tto = {
        valid_tto_days <- tto_days[tto_valid]
        if (length(valid_tto_days) > 0) as.character(first(valid_tto_days)) else "Unknown"
      },
      .groups = "drop"
    )

  # 关联drug表，添加TTO
  # drug_seq与dsg_drug_seq对应
  faers_drug_with_tto <- faers_drug %>%
    left_join(ther_tto, by = c("primaryid" = "primaryid", "drug_seq" = "dsg_drug_seq")) %>%
    mutate(
      tto = coalesce(tto, "Unknown")
    ) %>%
    select(primaryid, caseid, drug_seq, drugname_std, prod_ai, role_cod_std, quarter, tto)

  # 统计TTO情况
  tto_stats <- faers_drug_with_tto %>%
    summarise(
      total_records = n(),
      known_tto = sum(tto != "Unknown"),
      unknown_tto = sum(tto == "Unknown"),
      known_pct = round(known_tto / total_records * 100, 2)
    )
  message(paste("TTO统计 - 总记录:", tto_stats$total_records))
  message(paste("TTO已知:", tto_stats$known_tto, "(", tto_stats$known_pct, "%)"))
  message(paste("TTO未知:", tto_stats$unknown_tto))

  # 对于已知TTO的记录，显示分布
  known_tto <- faers_drug_with_tto %>%
    filter(tto != "Unknown") %>%
    mutate(tto_num = as.numeric(tto))
  if (nrow(known_tto) > 0) {
    message(paste("TTO范围: ", min(known_tto$tto_num), "至", max(known_tto$tto_num), "天"))
    message(paste("TTO中位数:", median(known_tto$tto_num), "天"))
  }

  # 更新并保存drug表（包含TTO）
  faers_drug <- faers_drug_with_tto
  save(faers_drug, file = file.path(.output_path, "faers_drug.RData"))
  message("Drug表已更新，包含TTO字段")

  # 释放内存
  rm(
    demo_event, ther_tto, faers_drug_with_tto, known_tto, tto_stats,
    demo_key, drug_key, reac_key, indi_key, dedup_df, dedup_final,
    demo_ids, drug_ids, reac_ids, indi_ids, common_ids,
    faers_demo, faers_drug, faers_reac, faers_indi,
    faers_outc, faers_rpsr, faers_ther
  )
  gc()
}


# 第三部分：JADER数据处理####
if (start_section %in% c("ALL", "JADER")) {
  ## 3.1 查看JADER文件####

  jader_path <- file.path(.data_path, "JADER")
  list.files(jader_path)

  # 函数：从JADER目录或压缩包名称中提取批次号（YYYYMM）
  extract_jader_release <- function(path) {
    file_name <- basename(path)
    release_label <- sub("^pmdacasereport([0-9]{6})(\\.zip)?$", "\\1", file_name)

    if (identical(release_label, file_name)) {
      return(NA_character_)
    }

    release_label
  }

  # 函数：返回JADER中最新的已解压目录
  get_latest_jader_dir <- function(jader_path) {
    jader_dirs <- list.dirs(jader_path, recursive = FALSE)
    jader_dirs <- jader_dirs[grepl("^pmdacasereport[0-9]{6}$", basename(jader_dirs))]

    if (length(jader_dirs) == 0) {
      return(NA_character_)
    }

    release_labels <- vapply(jader_dirs, extract_jader_release, character(1))
    jader_dirs[order(release_labels, decreasing = TRUE)][1]
  }

  # 函数：返回JADER中最新的压缩包
  get_latest_jader_zip <- function(jader_path) {
    jader_zips <- list.files(
      jader_path,
      pattern = "^pmdacasereport[0-9]{6}\\.zip$",
      full.names = TRUE
    )

    if (length(jader_zips) == 0) {
      return(NA_character_)
    }

    release_labels <- vapply(jader_zips, extract_jader_release, character(1))
    jader_zips[order(release_labels, decreasing = TRUE)][1]
  }

  latest_jader_dir <- get_latest_jader_dir(jader_path)
  latest_jader_zip <- get_latest_jader_zip(jader_path)

  if (is.na(latest_jader_dir) && is.na(latest_jader_zip)) {
    stop("未找到JADER已解压目录或压缩包")
  }

  if (!is.na(latest_jader_zip)) {
    latest_zip_release <- extract_jader_release(latest_jader_zip)
    latest_dir_release <- if (!is.na(latest_jader_dir)) extract_jader_release(latest_jader_dir) else NA_character_

    if (is.na(latest_jader_dir) || latest_zip_release > latest_dir_release) {
      message(sprintf("解压最新JADER压缩包: %s", basename(latest_jader_zip)))
      unzip(latest_jader_zip, exdir = jader_path)
      latest_jader_dir <- get_latest_jader_dir(jader_path)
    }
  }

  jader_csv_dir <- latest_jader_dir
  message(sprintf("使用JADER数据目录: %s", basename(jader_csv_dir)))
  list.files(jader_csv_dir)

  ## 3.2 读取JADER各表####

  # 函数：读取并重命名单个JADER表
  read_jader_table <- function(jader_csv_dir, table_prefix, col_names) {
    table_file <- list.files(jader_csv_dir, pattern = paste0("^", table_prefix, ".*\\.csv$"), full.names = TRUE)

    if (length(table_file) == 0) {
      stop(sprintf("未找到JADER %s表文件", table_prefix))
    }

    # JADER优先使用CP932读取，避免大文件在Shift-JIS路径下触发非法多字节序列崩溃
    read_csv_with_encoding <- function(enc) {
      readr::read_csv(
        table_file[1],
        locale = readr::locale(encoding = enc),
        show_col_types = FALSE,
        progress = FALSE,
        col_types = readr::cols(.default = readr::col_character())
      )
    }

    dt <- tryCatch(
      read_csv_with_encoding("CP932"),
      error = function(e) {
        message(sprintf("JADER %s 用readr读取失败(CP932): %s，回退到base::read.csv", table_prefix, e$message))
        utils::read.csv(
          table_file[1],
          fileEncoding = "CP932",
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
      }
    ) %>%
      data.table::as.data.table()

    rename_n <- min(ncol(dt), length(col_names))
    if (rename_n > 0) {
      data.table::setnames(dt, seq_len(rename_n), col_names[seq_len(rename_n)])
    }

    dt
  }

  jader_demo_cols <- c(
    "primaryid", "report_num", "sex", "age", "weight", "height",
    "report_quarter", "status", "report_type", "reporter_qual", "e2b"
  )
  jader_drug_cols <- c(
    "primaryid", "report_num", "drug_seq", "drug_role", "generic_name",
    "brand_name", "route", "start_date", "end_date", "dose", "dose_unit",
    "frequency", "indication", "action", "rechallenge", "risk_category"
  )
  jader_reac_cols <- c("primaryid", "report_num", "reaction_seq", "pt_term", "outcome", "event_date")
  jader_hist_cols <- c("primaryid", "report_num", "hist_seq", "medical_history")

  jader_specs <- list(
    demo = list(prefix = "demo", cols = jader_demo_cols),
    drug = list(prefix = "drug", cols = jader_drug_cols),
    reac = list(prefix = "reac", cols = jader_reac_cols),
    hist = list(prefix = "hist", cols = jader_hist_cols)
  )

  message("开始串行读取JADER各表")
  jader_tables <- lapply(
    names(jader_specs),
    function(table_name) {
      spec <- jader_specs[[table_name]]
      read_jader_table(jader_csv_dir, spec$prefix, spec$cols)
    }
  )
  names(jader_tables) <- names(jader_specs)

  # JADER使用Shift-JIS编码（日语），需要用readr读取
  # 读取demo表 (症例一览)
  jader_demo <- jader_tables$demo
  dim(jader_demo)
  head(jader_demo)
  colnames(jader_demo)
  head(jader_demo)

  # 读取drug表 (医药品信息)
  jader_drug <- jader_tables$drug
  head(jader_drug)

  # 读取reac表 (副作用)
  jader_reac <- jader_tables$reac
  dim(jader_reac)
  head(jader_reac)
  head(jader_reac)

  # 读取hist表 (原疾患)
  jader_hist <- jader_tables$hist
  dim(jader_hist)
  head(jader_hist)

  ## 3.3 加载MedDRA词典并转换日文PT为英文####
  message("\n加载MedDRA日英对照词典...")
  load(file.path(.data_path, "JADER/DIC.RData"))
  message(paste("MedDRA词典条目数:", nrow(DIC)))

  # 创建日文PT到英文PT的映射表
  # DIC中pt列是日文PT名称，pt_name_en是英文PT名称
  meddra_jp_en <- DIC %>%
    as_tibble() %>%
    select(pt_jp = pt, pt_en = pt_name_en) %>%
    distinct() %>%
    filter(!is.na(pt_jp) & pt_jp != "" & !is.na(pt_en) & pt_en != "")

  message(paste("有效日英映射条目数:", nrow(meddra_jp_en)))

  # 函数：标准化JADER部分日期，去掉时间部分并保留年月日精度
  normalize_jader_partial_date <- function(x) {
    x <- as.character(x)
    x <- trimws(x)
    x[x == ""] <- NA_character_

    digits <- gsub("[^0-9]", "", x)
    digits[is.na(x)] <- NA_character_

    normalized <- ifelse(nchar(digits) >= 8L, substr(digits, 1L, 8L), digits)
    valid_len <- !is.na(normalized) & nchar(normalized) %in% c(4L, 6L, 8L)
    normalized[!valid_len] <- NA_character_

    normalized
  }

  # 函数：并行执行JADER日文到英文的MedDRA映射
  translate_jader_meddra <- function(dt, join_col, output_col, meddra_jp_en) {
    translated_dt <- dt %>%
      as_tibble() %>%
      left_join(meddra_jp_en, by = setNames("pt_jp", join_col)) %>%
      mutate(!!output_col := ifelse(is.na(pt_en), .data[[join_col]], pt_en)) %>%
      select(-pt_en)

    translated_dt
  }

  message("开始串行进行JADER MedDRA转换")
  jader_translation_results <- lapply(
    c("reac", "hist", "drug"),
    function(table_name) {
      switch(table_name,
        reac = translate_jader_meddra(jader_reac, "pt_term", "pt_term_en", meddra_jp_en),
        hist = translate_jader_meddra(jader_hist, "medical_history", "medical_history_en", meddra_jp_en),
        drug = translate_jader_meddra(jader_drug, "indication", "indication_en", meddra_jp_en)
      )
    }
  )
  names(jader_translation_results) <- c("reac", "hist", "drug")

  # 将reac表的日文pt_term转换为英文
  message("转换reac表PT为英文...")
  jader_reac <- jader_translation_results$reac

  # 统计转换情况
  n_matched_reac <- sum(!is.na(jader_reac$pt_term_en) & jader_reac$pt_term_en != jader_reac$pt_term)
  message(paste(
    "reac表成功转换:", n_matched_reac, "/", nrow(jader_reac),
    "(", round(n_matched_reac / nrow(jader_reac) * 100, 1), "%)"
  ))

  # 将hist表的日文medical_history转换为英文
  message("转换hist表PT为英文...")
  jader_hist <- jader_translation_results$hist

  n_matched_hist <- sum(!is.na(jader_hist$medical_history_en) & jader_hist$medical_history_en != jader_hist$medical_history)
  message(paste(
    "hist表成功转换:", n_matched_hist, "/", nrow(jader_hist),
    "(", round(n_matched_hist / nrow(jader_hist) * 100, 1), "%)"
  ))

  # 将drug表的日文indication转换为英文
  message("转换drug表indication为英文...")
  jader_drug <- jader_translation_results$drug

  n_matched_drug <- sum(!is.na(jader_drug$indication_en) & jader_drug$indication_en != jader_drug$indication, na.rm = TRUE)
  message(paste(
    "drug表indication成功转换:", n_matched_drug, "/", sum(!is.na(jader_drug$indication)),
    "(", round(n_matched_drug / sum(!is.na(jader_drug$indication)) * 100, 1), "%)"
  ))

  ## 3.4 清洗JADER数据####

  # 查看report_quarter格式
  table(jader_demo$report_quarter)

  # 清洗demo表
  # 格式为 "2014・第一"，需要把日文季度转为阿拉伯数字
  jader_demo_clean <- jader_demo %>%
    as_tibble() %>%
    mutate(
      across(where(is.character), ~ na_if(trimws(.), "")),
      # 提取年龄数值 (格式如 "60代")
      age_decade = suppressWarnings(as.numeric(gsub("[^0-9]", "", age))),
      # 提取报告年份（前4位数字）
      report_year = as.numeric(substr(report_quarter, 1, 4)),
      # 提取季度（中文数字转阿拉伯数字）
      report_q = case_when(
        grepl("第一", report_quarter) ~ 1L,
        grepl("第二", report_quarter) ~ 2L,
        grepl("第三", report_quarter) ~ 3L,
        grepl("第四", report_quarter) ~ 4L,
        TRUE ~ NA_integer_
      ),
      # 创建年度季度数值用于筛选
      year_quarter = report_year * 10 + report_q,
      source = "JADER"
    ) %>%
    # 只保留2015-2025年的报告（排除NA）
    filter(!is.na(report_year) & report_year >= analysis_year_min & report_year <= analysis_year_max)
  table(jader_demo_clean$report_quarter)
  message(paste("筛选后JADER报告数:", nrow(jader_demo_clean)))
  # 使用year_quarter数值找到正确的时间范围
  min_idx <- which.min(jader_demo_clean$year_quarter)
  max_idx <- which.max(jader_demo_clean$year_quarter)
  message(paste(
    "时间范围:", jader_demo_clean$report_quarter[min_idx],
    "至", jader_demo_clean$report_quarter[max_idx]
  ))

  dim(jader_demo_clean)
  head(jader_demo_clean)
  table(jader_demo_clean$sex)
  summary(jader_demo_clean$age_decade)

  # 加载日英药物名映射表
  message("加载日英药物名映射表...")
  library(readxl)
  drug_name_mapping <- read_excel(
    file.path(.data_path, "JADER/药物名查询表格.xlsx"),
    col_types = "text"
  )
  names(drug_name_mapping) <- c("jp_name", "en_name", "cn_name", "data")
  drug_name_mapping <- drug_name_mapping %>%
    select(jp_name, en_name) %>%
    filter(!is.na(jp_name) & !is.na(en_name)) %>%
    mutate(
      jp_name = trimws(jp_name),
      en_name = toupper(trimws(en_name))
    ) %>%
    distinct(jp_name, en_name) %>%
    group_by(jp_name) %>%
    summarise(
      en_name_n = n_distinct(en_name),
      en_name = first(en_name),
      .groups = "drop"
    )

  ambiguous_drug_name_mapping <- drug_name_mapping %>%
    filter(en_name_n > 1)

  if (nrow(ambiguous_drug_name_mapping) > 0) {
    message(sprintf(
      "[WARN] 药名映射表中有%d个日文名对应多个英文名；当前默认保留首个英文名",
      nrow(ambiguous_drug_name_mapping)
    ))
  }

  drug_name_mapping <- drug_name_mapping %>%
    select(jp_name, en_name)

  message(paste("日英药物映射数:", nrow(drug_name_mapping)))

  # 清洗drug表（使用英文indication，日文药名转英文）
  jader_clean_results <- lapply(
    c("drug", "reac", "hist"),
    function(table_name) {
      switch(table_name,
        drug = jader_drug %>%
          mutate(
            across(where(is.character), ~ na_if(trimws(.), "")),
            generic_name_clean = trimws(generic_name)
          ) %>%
          left_join(drug_name_mapping, by = c("generic_name_clean" = "jp_name")) %>%
          mutate(
            generic_name_mapped = !is.na(en_name),
            generic_name_std = toupper(ifelse(!is.na(en_name), en_name, generic_name_clean)),
            indication_std = toupper(trimws(indication_en)),
            source = "JADER"
          ) %>%
          select(-generic_name_clean, -en_name),
        reac = jader_reac %>%
          mutate(
            across(where(is.character), ~ na_if(trimws(.), "")),
            pt_std = toupper(trimws(pt_term_en)),
            source = "JADER"
          ),
        hist = jader_hist %>%
          mutate(
            across(where(is.character), ~ na_if(trimws(.), "")),
            medical_history_std = toupper(trimws(medical_history_en)),
            source = "JADER"
          )
      )
    }
  )
  names(jader_clean_results) <- c("drug", "reac", "hist")
  jader_drug_clean <- jader_clean_results$drug

  message(paste("药名映射成功数:", sum(jader_drug_clean$generic_name_mapped, na.rm = TRUE)))

  jader_drug_clean <- jader_drug_clean %>%
    select(-generic_name_mapped)

  dim(jader_drug_clean)
  head(jader_drug_clean)

  # 清洗reac表（使用英文PT）
  jader_reac_clean <- jader_clean_results$reac

  dim(jader_reac_clean)
  head(jader_reac_clean)
  jader_reac_clean %>%
    count(pt_std, sort = TRUE) %>%
    head(20)

  # 清洗hist表（使用英文medical_history）
  jader_hist_clean <- jader_clean_results$hist

  dim(jader_hist_clean)

  ## 3.5 JADER去重####

  message("\n========== JADER第一步去重：版本去重 ==========")
  message(paste("JADER原始记录数:", nrow(jader_demo_clean)))
  message(paste("JADER原始病例数(case_id):", uniqueN(jader_demo_clean$primaryid)))

  # 第一步：每个case_id只保留report_num最大的记录
  jader_demo_clean <- jader_demo_clean %>%
    mutate(report_num_num = suppressWarnings(as.integer(report_num))) %>%
    group_by(primaryid) %>%
    arrange(
      desc(!is.na(report_num_num)),
      desc(report_num_num),
      desc(report_num),
      .by_group = TRUE
    ) %>%
    slice(1) %>%
    ungroup() %>%
    select(-report_num_num)

  message(paste("版本去重后记录数:", nrow(jader_demo_clean)))

  # 创建有效的primaryid+report_num组合
  valid_jader_keys <- jader_demo_clean %>% select(primaryid, report_num)

  # 筛选其他表
  jader_drug_clean <- jader_drug_clean %>% semi_join(valid_jader_keys, by = c("primaryid", "report_num"))
  jader_reac_clean <- jader_reac_clean %>% semi_join(valid_jader_keys, by = c("primaryid", "report_num"))
  jader_hist_clean <- jader_hist_clean %>% semi_join(valid_jader_keys, by = c("primaryid", "report_num"))

  message("\n========== JADER补充：将drug的indication合并到hist表 ==========")
  # 从drug表提取indication作为hist表的补充记录
  # 这样hist表就能完整对应FAERS的indi表（包含原疾患+用药适应证）
  drug_indication_for_hist <- jader_drug_clean %>%
    filter(!is.na(indication_std) & indication_std != "") %>%
    select(primaryid, report_num, medical_history_std = indication_std) %>%
    distinct()

  message(paste("Drug表中有效indication记录数:", nrow(drug_indication_for_hist)))

  # 合并到hist表（去重避免重复）
  jader_hist_clean <- bind_rows(
    jader_hist_clean %>% select(primaryid, report_num, medical_history_std),
    drug_indication_for_hist
  ) %>%
    distinct(primaryid, report_num, medical_history_std, .keep_all = TRUE)

  message(paste("合并后hist表记录数:", nrow(jader_hist_clean)))
  message(paste("合并后hist表主键数:", n_distinct(paste(jader_hist_clean$primaryid, jader_hist_clean$report_num))))

  message("\n========== JADER第二步去重：临床特征去重 ==========")
  # 去重键：年龄、性别、体重、身高、发生时间、用药、不良反应、适应证

  # 从demo提取关键特征
  jader_demo_key <- jader_demo_clean %>%
    select(primaryid, report_num, age, sex, weight, height)

  # 从drug提取用药
  jader_drug_key <- jader_drug_clean %>%
    group_by(primaryid, report_num) %>%
    summarise(
      drugs = paste(sort(unique(na.omit(generic_name_std))), collapse = "|"),
      .groups = "drop"
    )

  # 从hist提取适应证/原疾患
  jader_hist_key <- jader_hist_clean %>%
    group_by(primaryid, report_num) %>%
    summarise(
      indis = paste(sort(unique(na.omit(medical_history_std))), collapse = "|"),
      .groups = "drop"
    )

  # 从reac提取不良反应和发生时间
  jader_reac_key <- jader_reac_clean %>%
    group_by(primaryid, report_num) %>%
    summarise(
      reacs = paste(sort(unique(na.omit(pt_std))), collapse = "|"),
      event_dates = paste(sort(unique(na.omit(normalize_jader_partial_date(event_date)))), collapse = "|"),
      .groups = "drop"
    )

  # 合并创建去重键
  jader_dedup_df <- jader_demo_key %>%
    left_join(jader_drug_key, by = c("primaryid", "report_num")) %>%
    left_join(jader_hist_key, by = c("primaryid", "report_num")) %>%
    left_join(jader_reac_key, by = c("primaryid", "report_num")) %>%
    mutate(
      dedup_key = paste(age, sex, weight, height, event_dates, drugs, reacs, indis, sep = "||")
    )

  message(paste("去重前报告数:", nrow(jader_dedup_df)))
  message(paste("去重键数:", uniqueN(jader_dedup_df$dedup_key)))

  # 每个去重键只保留一条记录
  jader_dedup_final <- jader_dedup_df %>%
    group_by(dedup_key) %>%
    slice(1) %>%
    ungroup()

  message(paste("去重后报告数:", nrow(jader_dedup_final)))
  message(paste("第二步去重比例:", round((1 - nrow(jader_dedup_final) / nrow(jader_dedup_df)) * 100, 2), "%"))

  # 更新有效keys
  valid_jader_keys <- jader_dedup_final %>% select(primaryid, report_num)

  message("\n========== JADER第三步：四表共有主键筛选 ==========")

  # 找出四表共有的case_id+report_num
  demo_keys <- jader_demo_clean %>%
    select(primaryid, report_num) %>%
    distinct()
  drug_keys <- jader_drug_clean %>%
    select(primaryid, report_num) %>%
    distinct()
  reac_keys <- jader_reac_clean %>%
    select(primaryid, report_num) %>%
    distinct()
  hist_keys <- jader_hist_clean %>%
    select(primaryid, report_num) %>%
    distinct()

  common_keys <- demo_keys %>%
    semi_join(drug_keys, by = c("primaryid", "report_num")) %>%
    semi_join(reac_keys, by = c("primaryid", "report_num")) %>%
    semi_join(hist_keys, by = c("primaryid", "report_num")) %>%
    semi_join(valid_jader_keys, by = c("primaryid", "report_num"))

  message(paste("DEMO主键数:", nrow(demo_keys)))
  message(paste("DRUG主键数:", nrow(drug_keys)))
  message(paste("REAC主键数:", nrow(reac_keys)))
  message(paste("HIST主键数:", nrow(hist_keys)))
  message(paste("四表共有主键数:", nrow(common_keys)))

  # 用共有主键筛选各表
  jader_demo <- jader_demo_clean %>% semi_join(common_keys, by = c("primaryid", "report_num"))
  jader_drug <- jader_drug_clean %>% semi_join(common_keys, by = c("primaryid", "report_num"))
  jader_reac <- jader_reac_clean %>% semi_join(common_keys, by = c("primaryid", "report_num"))
  jader_hist <- jader_hist_clean %>% semi_join(common_keys, by = c("primaryid", "report_num"))

  message(paste("JADER最终报告数:", nrow(jader_demo)))

  # -----------------------------------------------------------------------------
  # 3.5 JADER TTO说明
  # JADER的drug和reac是多对多关系，TTO计算需在后续ICI筛选时进行
  # 此处保留原始日期字段：drug表的start_date，reac表的event_date
  # -----------------------------------------------------------------------------
  message("\n========== JADER TTO说明 ==========")
  message("JADER TTO计算将在后续ICI药物筛选时进行")
  message("当前保留原始日期字段：drug.start_date, reac.event_date")

  # 保存JADER数据
  save(jader_demo, file = file.path(.output_path, "jader_demo.RData"))
  save(jader_drug, file = file.path(.output_path, "jader_drug.RData"))
  save(jader_reac, file = file.path(.output_path, "jader_reac.RData"))
  save(jader_hist, file = file.path(.output_path, "jader_hist.RData"))

  message("JADER去重完成，数据已保存")

  # 释放内存
  rm(
    jader_demo, jader_drug, jader_reac, jader_hist,
    jader_demo_clean, jader_drug_clean, jader_reac_clean, jader_hist_clean,
    jader_demo_key, jader_drug_key, jader_hist_key, jader_reac_key, jader_dedup_df, jader_dedup_final,
    valid_jader_keys, demo_keys, drug_keys, reac_keys, common_keys,
    jader_translation_results, jader_clean_results
  )
  gc()
}
