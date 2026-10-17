#不成比例分析函数 (Disproportionality Analysis)####

# 加载必要包
library(tidyverse)
library(data.table)

##1.核心计算函数####

#' 计算ROR（报告比值比）及其95%置信区间
#'
#' @param a 癌种A且irAE(+)的人数
#' @param b 癌种A且irAE(-)的人数
#' @param c 非癌种A且irAE(+)的人数
#' @param d 非癌种A且irAE(-)的人数
#' @param alpha 显著性水平，默认0.05
#' @param p_value_method P值计算方法，可选"fisher"或"chisq"
#' @param continuity_correction 连续性校正值，默认0.5
#' @param force_continuity_correction 是否强制对四格表全部格子应用连续性校正，默认FALSE
#' @return 包含ROR、CI、p值等的列表
calculate_ROR <- function(
  a,
  b,
  c,
  d,
  alpha = 0.05,
  p_value_method = c("fisher", "chisq"),
  continuity_correction = 0.5,
  force_continuity_correction = FALSE
) {
  p_value_method <- match.arg(p_value_method)

  # 根据组合类型计算P值
  contingency_table <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  p_value <- tryCatch({
    if (p_value_method == "fisher") {
      fisher.test(contingency_table)$p.value
    } else {
      suppressWarnings(chisq.test(contingency_table, correct = FALSE)$p.value)
    }
  }, error = function(e) {
    NA_real_
  })
  
  # 应用连续性校正（低计数兜底时可强制应用）
  correction_applied <- force_continuity_correction || any(c(a, b, c, d) == 0)
  if (correction_applied) {
    a <- a + continuity_correction
    b <- b + continuity_correction
    c <- c + continuity_correction
    d <- d + continuity_correction
  }
  
  # 计算ROR
  ROR <- (a * d) / (b * c)
  
  # 计算对数ROR的标准误
  log_ROR <- log(ROR)
  SE_log_ROR <- sqrt(1/a + 1/b + 1/c + 1/d)
  
  # 计算95%置信区间
  z_alpha <- qnorm(1 - alpha/2)
  ROR_lower <- exp(log_ROR - z_alpha * SE_log_ROR)
  ROR_upper <- exp(log_ROR + z_alpha * SE_log_ROR)
  
  # 信号检测（ROR_lower > 1 表示风险升高信号）
  signal <- !is.na(ROR_lower) && ROR_lower > 1
  
  return(list(
    ROR = ROR,
    ROR_lower = ROR_lower,
    ROR_upper = ROR_upper,
    log_ROR = log_ROR,
    SE_log_ROR = SE_log_ROR,
    p_value = p_value,
    p_value_method = p_value_method,
    signal = signal,
    continuity_correction_applied = correction_applied
  ))
}


#' 构建2×2列联表（内部函数）
#'
#' @param data 输入数据框
#' @param cancer_var 癌种分类变量名
#' @param cancer_value 目标癌种值
#' @param irae_var irAE变量名
#' @return 包含a, b, c, d的列表
.build_contingency_table <- function(data, cancer_var, cancer_value, irae_var) {
  # 创建癌种和irAE的二分类变量
  is_cancer <- data[[cancer_var]] == cancer_value
  is_irae <- data[[irae_var]] == 1
  
  # 计算四格表
  a <- sum(is_cancer & is_irae, na.rm = TRUE)    # 癌种A且irAE(+)
  b <- sum(is_cancer & !is_irae, na.rm = TRUE)   # 癌种A且irAE(-)
  c <- sum(!is_cancer & is_irae, na.rm = TRUE)   # 非癌种A且irAE(+)
  d <- sum(!is_cancer & !is_irae, na.rm = TRUE)  # 非癌种A且irAE(-)
  
  return(list(
    a = a, b = b, c = c, d = d,
    n_cancer = a + b,
    n_non_cancer = c + d,
    n_irae = a + c,
    n_total = a + b + c + d
  ))
}


#' 计算2×2列联表的期望频数（内部函数）
#'
#' @param a 癌种A且irAE(+)的人数
#' @param b 癌种A且irAE(-)的人数
#' @param c 非癌种A且irAE(+)的人数
#' @param d 非癌种A且irAE(-)的人数
#' @return 命名数值向量，包含a_expected、b_expected、c_expected、d_expected
.calculate_expected_counts <- function(a, b, c, d) {
  row_totals <- c(a + b, c + d)
  col_totals <- c(a + c, b + d)
  grand_total <- sum(row_totals)

  expected_matrix <- outer(row_totals, col_totals) / grand_total

  c(
    a_expected = expected_matrix[1, 1],
    b_expected = expected_matrix[1, 2],
    c_expected = expected_matrix[2, 1],
    d_expected = expected_matrix[2, 2]
  )
}


#' 计算调整后的OR（多变量逻辑回归）
#'
#' @param data 输入数据框
#' @param cancer_var 癌种分类变量名
#' @param cancer_value 目标癌种值
#' @param irae_var irAE变量名
#' @param adjust_vars 协变量向量
#' @return 包含调整后OR、CI、p值的列表
.calculate_adjusted_OR <- function(data, cancer_var, cancer_value, irae_var, adjust_vars) {
  # 创建模型数据
  model_data <- data %>%
    mutate(
      cancer_binary = as.integer(.data[[cancer_var]] == cancer_value),
      irae_binary = as.integer(.data[[irae_var]] == 1)
    )
  
  # 构建公式
  formula_str <- paste("irae_binary ~ cancer_binary +", paste(adjust_vars, collapse = " + "))
  formula_obj <- as.formula(formula_str)
  
  # 拟合逻辑回归模型（抑制完全分离警告）
  result <- tryCatch({
    model <- suppressWarnings(
      glm(formula_obj, data = model_data, family = binomial(link = "logit"))
    )
    summary_model <- summary(model)
    coef_table <- summary_model$coefficients
    
    # 提取cancer_binary的系数
    if ("cancer_binary" %in% rownames(coef_table)) {
      coef <- coef_table["cancer_binary", "Estimate"]
      se <- coef_table["cancer_binary", "Std. Error"]
      p_value <- coef_table["cancer_binary", "Pr(>|z|)"]
      
      # 计算OR和95%CI
      OR <- exp(coef)
      OR_lower <- exp(coef - 1.96 * se)
      OR_upper <- exp(coef + 1.96 * se)
      
      # 处理完全分离导致的非有限值
      make_finite <- function(x) ifelse(is.finite(x), x, NA_real_)
      OR <- make_finite(OR)
      OR_lower <- make_finite(OR_lower)
      OR_upper <- make_finite(OR_upper)
      
      list(
        adjusted_OR = OR,
        adjusted_OR_lower = OR_lower,
        adjusted_OR_upper = OR_upper,
        adjusted_p_value = p_value,
        converged = model$converged
      )
    } else {
      list(
        adjusted_OR = NA_real_,
        adjusted_OR_lower = NA_real_,
        adjusted_OR_upper = NA_real_,
        adjusted_p_value = NA_real_,
        converged = FALSE
      )
    }
  }, error = function(e) {
    list(
      adjusted_OR = NA_real_,
      adjusted_OR_lower = NA_real_,
      adjusted_OR_upper = NA_real_,
      adjusted_p_value = NA_real_,
      converged = FALSE,
      error = conditionMessage(e)
    )
  })
  
  return(result)
}


##2.主分析函数####

#' 运行不成比例分析（自动遍历所有癌种和irAE组合）
#'
#' @param data 输入数据框（如faers_malignant_final）
#' @param cancer_var 癌种分类变量名（如"primary_organ"、"assigned_type"）
#' @param irae_level irAE层级："SOC"或"PT"
#' @param min_cases 最小病例数通用阈值（向后兼容；若未单独指定target/reference阈值，则两组都使用该值），默认5
#' @param min_cases_target 目标癌种组最小病例数阈值（a < min_cases_target的组合跳过），默认NULL时继承min_cases
#' @param min_cases_reference 参照组最小病例数阈值（c < min_cases_reference的组合跳过），默认NULL时继承min_cases
#' @param low_count_target_threshold 旧版目标组低计数兜底阈值；保留向后兼容，默认NULL
#' @param expected_count_threshold 期望频数阈值；当2×2列联表任一期望频数 < 该阈值时，改用Fisher精确检验 + 连续性校正ROR，默认NULL表示不启用
#' @param filter_ici 是否只筛选ICI使用者，默认TRUE
#' @param exclude_cancer 排除的癌种（如c("CUP", "Unknown")），默认NULL
#' @param adjust_vars 协变量调整变量名向量（如c("age_group", "sex_std")），默认NULL不调整
#' @param p_adjust_method 多重比较校正方法，默认"none"
#' @param target_cancer 二分类比较时指定目标组类别；默认NULL表示按排序后第一类
#' @param reference_cancer 二分类比较时指定参照组类别；默认NULL表示按排序后第二类
#' @return 包含多个矩阵和详细结果的列表
run_disproportionality_analysis <- function(
    data,
    cancer_var,
    irae_level,
    min_cases = 5,
    min_cases_target = NULL,
    min_cases_reference = NULL,
    low_count_target_threshold = NULL,
    expected_count_threshold = NULL,
    filter_ici = TRUE,
    exclude_cancer = NULL,
    adjust_vars = NULL,
    p_adjust_method = "none",
    target_cancer = NULL,
    reference_cancer = NULL
) {
  

  # 参数验证

  message("========== 不成比例分析开始 ==========")
  message(paste("癌种分类变量:", cancer_var))
  message(paste("irAE层级:", irae_level))
  message(paste("只分析ICI使用者:", filter_ici))

  if (is.null(min_cases_target)) {
    min_cases_target <- min_cases
  }
  if (is.null(min_cases_reference)) {
    min_cases_reference <- min_cases
  }

  message(paste("目标组最小病例阈值(a):", min_cases_target))
  message(paste("参照组最小病例阈值(c):", min_cases_reference))
  if (is.null(expected_count_threshold)) {
    message("期望频数阈值: 未启用")
  } else {
    message(paste("期望频数阈值(任一期望频数<):", expected_count_threshold))
    message("已启用按期望频数分流；全部组合保留，不再按a/c阈值跳过")
  }
  message("常规组合P值检验: Pearson chi-squared test")
  message("低期望频数组合P值检验: Fisher exact test")
  
  # 验证irae_level
  if (!irae_level %in% c("SOC", "PT")) {
    stop("irae_level必须是'SOC'或'PT'")
  }
  
  # 验证cancer_var存在
  if (!cancer_var %in% names(data)) {
    stop(paste("癌种变量不存在:", cancer_var))
  }
  

  # 数据筛选

  analysis_data <- data
  
  # 筛选ICI使用者
  if (filter_ici) {
    if (!"has_ici" %in% names(analysis_data)) {
      stop("数据中不存在has_ici变量，无法筛选ICI使用者")
    }
    analysis_data <- analysis_data %>% filter(has_ici == TRUE)
    message(paste("ICI使用者数量:", nrow(analysis_data)))
  }
  
  # 排除指定癌种
  if (!is.null(exclude_cancer)) {
    analysis_data <- analysis_data %>%
      filter(!.data[[cancer_var]] %in% exclude_cancer)
    message(paste("排除癌种后剩余:", nrow(analysis_data)))
  }
  
  # 移除NA
  analysis_data <- analysis_data %>%
    filter(!is.na(.data[[cancer_var]]))
  

  # 获取癌种和irAE变量列表

  # 获取所有癌种类别
  cancer_values <- unique(analysis_data[[cancer_var]])
  cancer_values <- cancer_values[!is.na(cancer_values)]
  n_cancer_categories <- length(cancer_values)
  message(paste("癌种类别数:", n_cancer_categories))
  
  if ((!is.null(target_cancer) || !is.null(reference_cancer)) && n_cancer_categories != 2) {
    stop("target_cancer/reference_cancer 仅适用于二分类比较")
  }

  # 如果只有两类癌种，只分析目标组（另一类是互补的，ROR互为倒数）
  is_binary <- n_cancer_categories == 2
  if (is_binary) {
    cancer_values <- sort(cancer_values)

    if (!is.null(target_cancer) && !target_cancer %in% cancer_values) {
      stop(paste("target_cancer 不在当前二分类水平中:", target_cancer))
    }
    if (!is.null(reference_cancer) && !reference_cancer %in% cancer_values) {
      stop(paste("reference_cancer 不在当前二分类水平中:", reference_cancer))
    }

    if (!is.null(target_cancer) && !is.null(reference_cancer) && identical(target_cancer, reference_cancer)) {
      stop("target_cancer 与 reference_cancer 不能相同")
    }

    if (is.null(target_cancer) && is.null(reference_cancer)) {
      target_cancer <- cancer_values[1]
      reference_cancer <- cancer_values[2]
    } else if (is.null(target_cancer)) {
      target_candidates <- setdiff(cancer_values, reference_cancer)
      if (length(target_candidates) != 1) {
        stop("无法根据 reference_cancer 唯一确定 target_cancer")
      }
      target_cancer <- target_candidates[[1]]
    } else if (is.null(reference_cancer)) {
      reference_candidates <- setdiff(cancer_values, target_cancer)
      if (length(reference_candidates) != 1) {
        stop("无法根据 target_cancer 唯一确定 reference_cancer")
      }
      reference_cancer <- reference_candidates[[1]]
    }

    cancer_values <- target_cancer
    message(paste("检测到二分类癌种分组，只展示:", target_cancer, "（参照组:", reference_cancer, "）"))
  }
  
  # 获取指定层级的所有irAE变量
  prefix <- paste0("^", irae_level, "_")
  irae_vars <- names(analysis_data)[grep(prefix, names(analysis_data))]
  
  if (length(irae_vars) == 0) {
    stop(paste("未找到以", irae_level, "_为前缀的irAE变量"))
  }
  message(paste("irAE变量数:", length(irae_vars)))
  

  # 初始化结果收集器

  results_list <- list()
  total_combinations <- length(cancer_values) * length(irae_vars)
  processed <- 0
  skipped <- 0
  fisher_test_count <- 0
  skipped_target_count <- 0
  skipped_reference_count <- 0
  
  message(paste("总组合数:", total_combinations))
  message("开始计算ROR...")
  start_time <- Sys.time()
  

  # 遍历所有癌种×irAE组合

  for (cancer in cancer_values) {
    for (irae in irae_vars) {
      # 构建列联表
      table <- .build_contingency_table(analysis_data, cancer_var, cancer, irae)
      expected_counts <- .calculate_expected_counts(table$a, table$b, table$c, table$d)
      min_expected_count <- min(expected_counts, na.rm = TRUE)

      if (!is.null(expected_count_threshold)) {
        use_fisher_test <- any(expected_counts < expected_count_threshold)
      } else {
        use_fisher_test <- !is.null(low_count_target_threshold) &&
          table$a <= low_count_target_threshold

        # 旧版逻辑下，仍按目标组/参照组阈值决定是否跳过
        if (!use_fisher_test && table$a < min_cases_target) {
          skipped <- skipped + 1
          skipped_target_count <- skipped_target_count + 1
          next
        }

        if (table$c < min_cases_reference) {
          skipped <- skipped + 1
          skipped_reference_count <- skipped_reference_count + 1
          next
        }
      }
      
      # 计算ROR；低期望频数组合改用Fisher精确检验 + 连续性校正ROR
      ror_result <- calculate_ROR(
        table$a,
        table$b,
        table$c,
        table$d,
        p_value_method = if (use_fisher_test) "fisher" else "chisq",
        force_continuity_correction = use_fisher_test
      )

      analysis_method <- if (use_fisher_test) {
        "Fisher exact test + continuity-corrected ROR"
      } else {
        "Pearson chi-squared test + standard ROR"
      }

      result_note <- if (use_fisher_test) {
        if (!is.null(expected_count_threshold)) {
          paste0(
            "最小期望频数=",
            round(min_expected_count, 3),
            "，小于阈值",
            expected_count_threshold,
            "；采用Fisher精确检验 + 连续性校正ROR"
          )
        } else {
          paste0(
            "目标组病例数a=",
            table$a,
            "，不满足大于",
            low_count_target_threshold,
            "；采用Fisher精确检验 + 连续性校正ROR"
          )
        }
      } else {
        NA_character_
      }
      
      # 初始化结果行
      result_row <- data.frame(
        cancer_type = cancer,
        irae = sub(prefix, "", irae),  # 移除前缀
        irae_var = irae,
        a = table$a,
        b = table$b,
        c = table$c,
        d = table$d,
        n_cancer = table$n_cancer,
        n_irae = table$n_irae,
        n_total = table$n_total,
        ROR = ror_result$ROR,
        ROR_lower = ror_result$ROR_lower,
        ROR_upper = ror_result$ROR_upper,
        p_value = ror_result$p_value,
        min_expected_count = min_expected_count,
        p_value_method = if (ror_result$p_value_method == "fisher") "Fisher exact test" else "Pearson chi-squared test",
        fisher_test_used = use_fisher_test,
        analysis_method = analysis_method,
        result_note = result_note,
        stringsAsFactors = FALSE
      )
      
      # 如果需要调整
      if (!is.null(adjust_vars)) {
        if (use_fisher_test) {
          result_row$adjusted_OR <- NA_real_
          result_row$adjusted_OR_lower <- NA_real_
          result_row$adjusted_OR_upper <- NA_real_
          result_row$adjusted_p_value <- NA_real_
          result_row$adjusted_model_fitted <- FALSE
          result_row$result_note <- paste0(
            result_row$result_note,
            "；调整模型未拟合"
          )
        } else {
          adj_result <- .calculate_adjusted_OR(
            analysis_data, cancer_var, cancer, irae, adjust_vars
          )
          result_row$adjusted_OR <- adj_result$adjusted_OR
          result_row$adjusted_OR_lower <- adj_result$adjusted_OR_lower
          result_row$adjusted_OR_upper <- adj_result$adjusted_OR_upper
          result_row$adjusted_p_value <- adj_result$adjusted_p_value
          result_row$adjusted_model_fitted <- isTRUE(adj_result$converged)
        }
      }
      
      results_list[[length(results_list) + 1]] <- result_row
      processed <- processed + 1
      if (use_fisher_test) {
        fisher_test_count <- fisher_test_count + 1
      }
    }
  }
  
  end_time <- Sys.time()
  message(paste("计算完成，耗时:", round(difftime(end_time, start_time, units = "secs"), 2), "秒"))
  message(paste("有效组合:", processed, "/ 跳过:", skipped))
  message(paste("Fisher检验组合数:", fisher_test_count))
  if (skipped > 0) {
    message(paste("因目标组未达阈值而跳过:", skipped_target_count))
    message(paste("因参照组未达阈值而跳过:", skipped_reference_count))
  }
  

  # 合并结果并进行FDR校正

  if (length(results_list) == 0) {
    warning("没有符合条件的癌种×irAE组合，返回空结果")
    return(NULL)
  }
  
  details <- bind_rows(results_list)
  
  # FDR校正（未调整p值）
  details$p_adjusted <- p.adjust(details$p_value, method = p_adjust_method)
  # 信号检测：只要p_adjusted < 0.05即为显著差异（ROR>1风险升高，ROR<1风险降低都有意义）
  details$signal_detected <- !is.na(details$p_adjusted) &
    details$p_adjusted < 0.05
  
  # FDR校正（调整后p值，如果有）
  if (!is.null(adjust_vars) && "adjusted_p_value" %in% names(details)) {
    details$adjusted_p_adjusted <- p.adjust(details$adjusted_p_value, method = p_adjust_method)
    # 调整后信号检测：只要adjusted_p_adjusted < 0.05即为显著差异
    details$adjusted_signal_detected <- !is.na(details$adjusted_p_adjusted) &
      details$adjusted_p_adjusted < 0.05
  }
  
  message(paste("检出信号数（未调整）:", sum(details$signal_detected, na.rm = TRUE)))
  if (!is.null(adjust_vars)) {
    message(paste("检出信号数（调整后）:", sum(details$adjusted_signal_detected, na.rm = TRUE)))
  }
  

  # 转换为矩阵格式

  message("转换为矩阵格式...")
  
  # 获取唯一的癌种和irAE
  all_cancers <- sort(unique(details$cancer_type))
  all_iraes <- sort(unique(details$irae))
  
  # 创建矩阵模板
  create_matrix <- function(values, fill_value = NA) {
    mat <- matrix(fill_value, 
                  nrow = length(all_cancers), 
                  ncol = length(all_iraes),
                  dimnames = list(all_cancers, all_iraes))
    return(mat)
  }
  
  # 填充矩阵
  ror_matrix <- create_matrix(NA_real_)
  ci_lower_matrix <- create_matrix(NA_real_)
  ci_upper_matrix <- create_matrix(NA_real_)
  p_matrix <- create_matrix(NA_real_)
  p_adjusted_matrix <- create_matrix(NA_real_)
  signal_matrix <- create_matrix(FALSE)
  
  for (i in seq_len(nrow(details))) {
    row <- details[i, ]
    cancer <- row$cancer_type
    irae <- row$irae
    
    ror_matrix[cancer, irae] <- row$ROR
    ci_lower_matrix[cancer, irae] <- row$ROR_lower
    ci_upper_matrix[cancer, irae] <- row$ROR_upper
    p_matrix[cancer, irae] <- row$p_value
    p_adjusted_matrix[cancer, irae] <- row$p_adjusted
    signal_matrix[cancer, irae] <- row$signal_detected
  }
  
  # 调整后矩阵（如果有）
  if (!is.null(adjust_vars)) {
    adjusted_or_matrix <- create_matrix(NA_real_)
    adjusted_ci_lower_matrix <- create_matrix(NA_real_)
    adjusted_ci_upper_matrix <- create_matrix(NA_real_)
    adjusted_p_matrix <- create_matrix(NA_real_)
    adjusted_p_adjusted_matrix <- create_matrix(NA_real_)
    adjusted_signal_matrix <- create_matrix(FALSE)
    
    for (i in seq_len(nrow(details))) {
      row <- details[i, ]
      cancer <- row$cancer_type
      irae <- row$irae
      
      adjusted_or_matrix[cancer, irae] <- row$adjusted_OR
      adjusted_ci_lower_matrix[cancer, irae] <- row$adjusted_OR_lower
      adjusted_ci_upper_matrix[cancer, irae] <- row$adjusted_OR_upper
      adjusted_p_matrix[cancer, irae] <- row$adjusted_p_value
      adjusted_p_adjusted_matrix[cancer, irae] <- row$adjusted_p_adjusted
      adjusted_signal_matrix[cancer, irae] <- row$adjusted_signal_detected
    }
  }
  

  # 构建返回结果

  result <- list(
    # 未调整结果矩阵
    ror_matrix = ror_matrix,
    ci_lower_matrix = ci_lower_matrix,
    ci_upper_matrix = ci_upper_matrix,
    p_matrix = p_matrix,
    p_adjusted_matrix = p_adjusted_matrix,
    signal_matrix = signal_matrix,
    
    # 详细结果
    details = details,
    
    # 分析参数
    params = list(
      cancer_var = cancer_var,
      irae_level = irae_level,
      n_cancers = length(all_cancers),
      n_iraes = length(all_iraes),
      filter_ici = filter_ici,
      exclude_cancer = exclude_cancer,
      adjust_vars = adjust_vars,
      min_cases = min_cases,
      min_cases_target = min_cases_target,
      min_cases_reference = min_cases_reference,
      low_count_target_threshold = low_count_target_threshold,
      expected_count_threshold = expected_count_threshold,
      standard_test_method = "Pearson chi-squared test",
      fisher_test_method = "Fisher exact test",
      p_adjust_method = p_adjust_method,
      n_processed = processed,
      n_skipped = skipped,
      n_fisher_test = fisher_test_count,
      n_skipped_target = skipped_target_count,
      n_skipped_reference = skipped_reference_count,
      is_binary = is_binary,
      target_cancer = if (is_binary) target_cancer else NULL,
      reference_cancer = if (is_binary) reference_cancer else NULL
    )
  )
  
  # 添加调整后结果（如果有）
  if (!is.null(adjust_vars)) {
    result$adjusted_or_matrix <- adjusted_or_matrix
    result$adjusted_ci_lower_matrix <- adjusted_ci_lower_matrix
    result$adjusted_ci_upper_matrix <- adjusted_ci_upper_matrix
    result$adjusted_p_matrix <- adjusted_p_matrix
    result$adjusted_p_adjusted_matrix <- adjusted_p_adjusted_matrix
    result$adjusted_signal_matrix <- adjusted_signal_matrix
  }
  
  message("========== 不成比例分析完成 ==========")
  
  return(result)
}


##3.辅助函数####

#' 打印分析结果摘要
#'
#' @param result run_disproportionality_analysis()的返回结果
print_summary <- function(result) {
  if (is.null(result)) {
    message("结果为空")
    return(invisible(NULL))
  }
  
  params <- result$params
  details <- result$details
  
  cat("\n========== 不成比例分析结果摘要 ==========\n")
  cat(paste("癌种分类变量:", params$cancer_var, "\n"))
  cat(paste("irAE层级:", params$irae_level, "\n"))
  if (isTRUE(params$is_binary)) {
    target_label <- if (!is.null(params$target_cancer)) params$target_cancer else unique(details$cancer_type)[1]
    cat(paste("二分类分组: 是（", target_label, " vs ", params$reference_cancer, "）\n", sep = ""))
  }
  cat(paste("癌种类别数:", params$n_cancers, "\n"))
  cat(paste("irAE变量数:", params$n_iraes, "\n"))
  cat(paste("有效组合数:", params$n_processed, "\n"))
  cat(paste("跳过组合数:", params$n_skipped, "\n"))
  cat(paste("FDR校正方法:", params$p_adjust_method, "\n"))
  
  cat("\n--- 未调整分析 ---\n")
  cat(paste("检出信号数:", sum(details$signal_detected, na.rm = TRUE), "\n"))
  
  # 显示Top 10信号
  top_signals <- details %>%
    filter(signal_detected) %>%
    arrange(desc(ROR)) %>%
    head(10)
  
  if (nrow(top_signals) > 0) {
    cat("\nTop 10 信号（按ROR降序）:\n")
    print(top_signals %>% 
            select(cancer_type, irae, ROR, ROR_lower, ROR_upper, p_adjusted))
  }
  
  # 调整后结果（如果有）
  if (!is.null(params$adjust_vars)) {
    cat("\n--- 调整后分析 ---\n")
    cat(paste("调整变量:", paste(params$adjust_vars, collapse = ", "), "\n"))
    cat(paste("检出信号数:", sum(details$adjusted_signal_detected, na.rm = TRUE), "\n"))
  }
  
  cat("\n=========================================\n")
}


#' 导出分析结果到Excel
#'
#' @param result run_disproportionality_analysis()的返回结果
#' @param output_path 输出文件路径
export_to_excel <- function(result, output_path) {
  if (!require(writexl)) {
    install.packages("writexl")
    library(writexl)
  }
  
  low_count_notes <- result$details %>%
    dplyr::filter(.data$fisher_test_used %in% TRUE) %>%
    dplyr::select(any_of(c(
      "cancer_type", "irae", "irae_full", "a", "b", "c", "d",
      "min_expected_count", "ROR", "ROR_lower", "ROR_upper", "p_value",
      "p_value_method", "analysis_method", "result_note"
    )))

  sheets <- list(
    "详细结果" = result$details,
    "ROR矩阵" = as.data.frame(result$ror_matrix) %>% 
      tibble::rownames_to_column("cancer_type"),
    "信号矩阵" = as.data.frame(result$signal_matrix) %>% 
      tibble::rownames_to_column("cancer_type"),
    "Fisher备注" = low_count_notes,
    "分析参数" = data.frame(
      参数 = names(result$params),
      值 = sapply(result$params, function(x) paste(x, collapse = ", "))
    )
  )
  
  # 添加调整后矩阵（如果有）
  if (!is.null(result$adjusted_or_matrix)) {
    sheets$"调整后OR矩阵" <- as.data.frame(result$adjusted_or_matrix) %>%
      tibble::rownames_to_column("cancer_type")
    sheets$"调整后信号矩阵" <- as.data.frame(result$adjusted_signal_matrix) %>%
      tibble::rownames_to_column("cancer_type")
  }
  
  write_xlsx(sheets, output_path)
  message(paste("结果已导出至:", output_path))
}


##4.使用示例（注释状态）####

# # 加载数据
# load("Data/Cleaned/faers_malignant_final.RData")
# 
# # 示例1：SOC层级，按器官分类，不调整
# result_soc <- run_disproportionality_analysis(
#   data = faers_malignant_final,
#   cancer_var = "primary_organ",
#   irae_level = "SOC",
#   filter_ici = TRUE,
#   exclude_cancer = c("CUP", "Unknown")
# )
# print_summary(result_soc)
# 
# # 示例2：PT层级，按TCGA分类，协变量调整
# result_pt <- run_disproportionality_analysis(
#   data = faers_malignant_final,
#   cancer_var = "assigned_type",
#   irae_level = "PT",
#   filter_ici = TRUE,
#   exclude_cancer = c("Others"),
#   adjust_vars = c("age_group", "sex_std", "treatment_strategy")
# )
# print_summary(result_pt)
# 
# # 导出结果
# export_to_excel(result_soc, "Results/disproportionality_SOC.xlsx")
