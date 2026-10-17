# FAERS irAE reporting summary helper functions####
# 用途：为4.2脚本提供FAERS实体瘤irAE发生谱汇总与绘图函数

### 标准化MedDRA术语，便于PT和SOC名称匹配####
standardize_medra_term <- function(x) {
    x <- as.character(x)
    x <- stringr::str_replace_all(x, "_", " ")
    x <- stringr::str_replace_all(x, "-", " ")
    x <- stringr::str_replace_all(x, "\\s+", " ")
    x <- trimws(x)
    toupper(x)
}

### 将变量名转换为更友好的显示文本####
format_display_term <- function(x) {
    x <- as.character(x)
    x <- stringr::str_replace_all(x, "_", " ")
    x <- stringr::str_replace_all(x, "\\s+", " ")
    x <- trimws(x)
    stringr::str_to_title(tolower(x))
}

### 格式化样本量，便于表格和图中标注####
format_case_n <- function(x) {
    format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
}

### 从quarter字段提取年份####
extract_report_year <- function(x) {
    suppressWarnings(as.integer(substr(as.character(x), 1, 4)))
}

### 为primary organ提供固定配色，输入primary_organs可选器官向量，输出命名颜色向量####
get_primary_organ_color_map <- function(primary_organs = NULL) {
    primary_organ_color_map <- c(
        "Breast" = "#63C65E",
        "Lung" = "#2F6DB2",
        "Prostate" = "#5BC0BE",
        "Kidney" = "#39A74A",
        "Bowel" = "#2E9C8E",
        "Skin" = "#6EC1E4",
        "Ovary" = "#FF4B4B",
        "Liver" = "#D44C8C",
        "Esophagus/stomach" = "#D6A23D",
        "Pancreas" = "#7A9E1F",
        "Bladder" = "#B89220",
        "Head/neck" = "#F4A259",
        "Brain" = "#9B63A8",
        "Uterus" = "#EC8FB4",
        "Soft tissue" = "#C06C84",
        "Thyroid" = "#F39AAE",
        "Biliary tract" = "#1F78B4",
        "Cervix" = "#F08C67",
        "PNS" = "#4F7F8A",
        "Pleura" = "#7F7F7F",
        "Bone" = "#FF7F00",
        "Eye" = "#6EC1E4",
        "Rare" = "#7F6AAD",
        "Thymus" = "#A58BC7",
        "Adrenal gland" = "#E31A1C",
        "Vulva/Vagina" = "#33A02C",
        "Penis" = "#6A3D9A",
        "Testis" = "#94A83D",
        "Ampulla of Vater" = "#A6CEE3",
        "Peritoneum" = "#8E7C72"
    )

    if (is.null(primary_organs)) {
        return(primary_organ_color_map)
    }

    primary_organs <- unique(as.character(primary_organs))
    primary_organs <- primary_organs[!is.na(primary_organs) & primary_organs != ""]

    missing_organs <- setdiff(primary_organs, names(primary_organ_color_map))

    if (length(missing_organs) > 0) {
        fallback_colors <- scales::hue_pal(l = 68, c = 100)(length(missing_organs))
        names(fallback_colors) <- missing_organs
        primary_organ_color_map <- c(primary_organ_color_map, fallback_colors)
    }

    primary_organ_color_map[primary_organs]
}

### 从RData中加载指定对象####
load_named_rdata <- function(input_path, object_name) {
    input_env <- new.env(parent = emptyenv())
    load(input_path, envir = input_env)

    if (!exists(object_name, envir = input_env)) {
        stop(paste0("RData中未找到对象: ", object_name))
    }

    input_env[[object_name]]
}

### 准备ICI使用实体瘤发生谱分析集####
prepare_faers_solid_ici_reporting_data <- function(
    data,
    exclude_primary_organs = c("CUP", "Hematologic", "Unmapped", "Unknown_Primary", "Unknown")
) {
    data_step0 <- data %>%
        mutate(primary_organ = as.character(primary_organ))

    data_step1 <- data_step0 %>%
        filter(has_ici == TRUE)

    data_step2 <- data_step1 %>%
        filter(
            !is.na(primary_organ),
            !primary_organ %in% exclude_primary_organs
        )

    filter_log <- tibble::tibble(
        step = c(
            "faers_malignant_final",
            "仅保留ICI使用病例",
            "排除原发器官未知/不纳入类别"
        ),
        n = c(
            nrow(data_step0),
            nrow(data_step1),
            nrow(data_step2)
        )
    ) %>%
        mutate(retained_pct = round(n / first(n) * 100, 2))

    list(
        filtered_data = data_step2,
        filter_log = filter_log,
        excluded_organs = exclude_primary_organs
    )
}

### 读取PT与SOC映射字典####
load_irae_dictionary <- function(dict_path) {
    irae_dictionary <- readxl::read_xlsx(dict_path) %>%
        mutate(
            PT = as.character(PT),
            SOC_full = as.character(SOC),
            pt_key = standardize_medra_term(PT),
            SOC = vapply(SOC_full, get_soc_abbrev, character(1))
        ) %>%
        distinct(pt_key, .keep_all = TRUE)

    if (!all(c("PT", "SOC_full", "SOC", "pt_key") %in% names(irae_dictionary))) {
        stop("irAE字典缺少必要字段")
    }

    irae_dictionary
}

### 从数据集中提取SOC元数据####
extract_soc_metadata <- function(data) {
    soc_vars <- grep("^SOC_", names(data), value = TRUE)

    if (length(soc_vars) == 0) {
        stop("数据中未检测到SOC列")
    }

    tibble::tibble(
        soc_var = soc_vars,
        SOC_full = gsub("_", " ", sub("^SOC_", "", soc_vars)),
        SOC = vapply(gsub("_", " ", sub("^SOC_", "", soc_vars)), get_soc_abbrev, character(1))
    ) %>%
        mutate(
            SOC = dplyr::coalesce(SOC, SOC_full)
        )
}

### 从数据集中提取PT元数据并关联SOC####
extract_pt_metadata <- function(data, irae_dictionary) {
    pt_vars <- grep("^PT_", names(data), value = TRUE)

    if (length(pt_vars) == 0) {
        stop("数据中未检测到PT列")
    }

    tibble::tibble(
        pt_var = pt_vars,
        pt_raw = sub("^PT_", "", pt_vars),
        pt_key = standardize_medra_term(pt_raw),
        PT_fallback = format_display_term(pt_raw)
    ) %>%
        left_join(
            irae_dictionary %>% select(pt_key, PT, SOC, SOC_full),
            by = "pt_key"
        ) %>%
        mutate(
            PT = dplyr::coalesce(as.character(PT), PT_fallback),
            SOC = as.character(SOC),
            SOC_full = as.character(SOC_full)
        ) %>%
        select(pt_var, PT, SOC, SOC_full, pt_key)
}

### 汇总队列定义、核心指标与说明文字####
build_cohort_summary <- function(data, excluded_organs) {
    total_n <- nrow(data)
    any_irae_n <- sum(data$has_irae == TRUE, na.rm = TRUE)
    any_irae_rp <- any_irae_n / total_n * 100

    tibble::tribble(
        ~section, ~item, ~value,
        "队列定义", "数据源", "faers_malignant_final",
        "队列定义", "纳入病例", "仅保留 has_ici == TRUE 的实体瘤病例",
        "队列定义", "排除原发器官", paste(excluded_organs, collapse = ", "),
        "队列定义", "分析单位", "每行对应单个FAERS report/case",
        "队列定义", "主分层变量", "primary_organ",
        "样本概览", "总样本量", format_case_n(total_n),
        "样本概览", "原发器官数", as.character(dplyr::n_distinct(data$primary_organ)),
        "样本概览", "Any-irAE报告数", format_case_n(any_irae_n),
        "样本概览", "Any-irAE reporting proportion (%)", sprintf("%.2f", any_irae_rp),
        "说明", "多重报告说明", "同一病例可同时命中多个SOC/PT，因此各SOC/PT的RP之和可超过100%",
        "说明", "术语说明", "全文统一使用report count、reporting proportion和reporting spectrum，不使用incidence/发生率",
        "说明", "小样本器官提示", "主图保留全部24个器官；N<300的器官解释需谨慎"
    )
}

### 汇总每例命中的SOC个数分布####
summarise_soc_count_distribution <- function(data) {
    soc_vars <- grep("^SOC_", names(data), value = TRUE)
    soc_count <- rowSums(data[, soc_vars, drop = FALSE] == 1, na.rm = TRUE)

    tibble::tibble(soc_count = soc_count) %>%
        mutate(
            soc_count_group = dplyr::case_when(
                soc_count == 0 ~ "0",
                soc_count == 1 ~ "1",
                soc_count == 2 ~ "2",
                soc_count >= 3 ~ ">=3",
                TRUE ~ NA_character_
            ),
            soc_count_group = factor(soc_count_group, levels = c("0", "1", "2", ">=3"))
        ) %>%
        count(soc_count_group, name = "n") %>%
        tidyr::complete(
            soc_count_group = factor(c("0", "1", "2", ">=3"), levels = c("0", "1", "2", ">=3")),
            fill = list(n = 0L)
        ) %>%
        mutate(
            RP = n / sum(n) * 100
        )
}

### 汇总各器官Any-irAE报告比例####
summarise_any_irae_by_organ <- function(data) {
    data %>%
        group_by(primary_organ) %>%
        summarise(
            N = n(),
            n_any = sum(has_irae == TRUE, na.rm = TRUE),
            RP = n_any / N * 100,
            .groups = "drop"
        ) %>%
        arrange(desc(N), primary_organ)
}

### 汇总整体SOC报告情况####
summarise_soc_overall <- function(data) {
    soc_metadata <- extract_soc_metadata(data)
    soc_counts <- colSums(data[, soc_metadata$soc_var, drop = FALSE] == 1, na.rm = TRUE)
    total_n <- nrow(data)

    soc_metadata %>%
        mutate(
            n = as.integer(soc_counts[soc_var]),
            RP = n / total_n * 100
        ) %>%
        select(SOC, SOC_full, n, RP, soc_var) %>%
        arrange(desc(RP), desc(n), SOC)
}

### 汇总各器官SOC报告情况####
summarise_soc_by_organ <- function(data) {
    soc_metadata <- extract_soc_metadata(data)
    organ_summary <- data %>%
        count(primary_organ, name = "N_organ") %>%
        arrange(desc(N_organ), primary_organ)

    data %>%
        group_by(primary_organ) %>%
        summarise(
            across(
                all_of(soc_metadata$soc_var),
                ~ sum(.x == 1, na.rm = TRUE)
            ),
            .groups = "drop"
        ) %>%
        pivot_longer(
            cols = all_of(soc_metadata$soc_var),
            names_to = "soc_var",
            values_to = "n"
        ) %>%
        left_join(soc_metadata, by = "soc_var") %>%
        left_join(organ_summary, by = "primary_organ") %>%
        mutate(
            RP = n / N_organ * 100
        ) %>%
        select(primary_organ, N_organ, SOC, SOC_full, n, RP, soc_var) %>%
        arrange(desc(N_organ), primary_organ, desc(RP), desc(n), SOC)
}

### 提取每个器官中排名靠前的SOC####
build_top_soc_by_organ <- function(soc_by_organ_df, top_n = 3) {
    soc_by_organ_df %>%
        group_by(primary_organ) %>%
        arrange(desc(n), desc(RP), SOC, .by_group = TRUE) %>%
        mutate(rank_within_organ = row_number()) %>%
        slice_head(n = top_n) %>%
        ungroup()
}

### 汇总整体PT报告情况####
summarise_pt_overall <- function(data, irae_dictionary) {
    pt_metadata <- extract_pt_metadata(data, irae_dictionary)
    pt_counts <- colSums(data[, pt_metadata$pt_var, drop = FALSE] == 1, na.rm = TRUE)
    total_n <- nrow(data)

    pt_metadata %>%
        mutate(
            n = as.integer(pt_counts[pt_var]),
            RP = n / total_n * 100
        ) %>%
        arrange(desc(n), desc(RP), PT) %>%
        group_by(SOC) %>%
        mutate(rank_within_soc = row_number()) %>%
        ungroup() %>%
        mutate(rank_overall = row_number()) %>%
        select(PT, SOC, SOC_full, n, RP, rank_overall, rank_within_soc, pt_var)
}

### 按SOC提取排名靠前的PT####
select_top_pt_by_soc <- function(pt_summary, top_n = 5) {
    pt_summary %>%
        filter(!is.na(SOC), SOC != "") %>%
        group_by(SOC) %>%
        arrange(desc(n), desc(RP), PT, .by_group = TRUE) %>%
        mutate(rank_within_soc = row_number()) %>%
        slice_head(n = top_n) %>%
        ungroup()
}

### 汇总各器官在selected PT上的报告情况####
summarise_pt_by_organ <- function(data, selected_pt_df) {
    if (nrow(selected_pt_df) == 0) {
        return(tibble::tibble())
    }

    organ_summary <- data %>%
        count(primary_organ, name = "N_organ") %>%
        arrange(desc(N_organ), primary_organ)

    selected_pt_metadata <- selected_pt_df %>%
        select(PT, SOC, SOC_full, rank_within_soc, pt_var) %>%
        distinct()

    pt_counts_by_organ <- data %>%
        group_by(primary_organ) %>%
        summarise(
            across(
                all_of(selected_pt_metadata$pt_var),
                ~ sum(.x == 1, na.rm = TRUE)
            ),
            .groups = "drop"
        ) %>%
        pivot_longer(
            cols = all_of(selected_pt_metadata$pt_var),
            names_to = "pt_var",
            values_to = "n"
        )

    tidyr::expand_grid(
        primary_organ = organ_summary$primary_organ,
        pt_var = selected_pt_metadata$pt_var
    ) %>%
        left_join(pt_counts_by_organ, by = c("primary_organ", "pt_var")) %>%
        mutate(n = tidyr::replace_na(n, 0L)) %>%
        left_join(selected_pt_metadata, by = "pt_var") %>%
        left_join(organ_summary, by = "primary_organ") %>%
        mutate(
            RP = n / N_organ * 100
        ) %>%
        select(primary_organ, N_organ, SOC, SOC_full, PT, rank_within_soc, n, RP, pt_var) %>%
        arrange(desc(N_organ), primary_organ, SOC, rank_within_soc, PT)
}

### 汇总SOC共报告模式####
build_soc_cooccurrence_table <- function(data) {
    soc_binary_matrix <- build_soc_upset_matrix(data)

    combination_label <- apply(soc_binary_matrix, 1, function(x) {
        active_soc <- colnames(soc_binary_matrix)[x]
        if (length(active_soc) == 0) {
            return("None")
        }
        paste(active_soc, collapse = " + ")
    })

    tibble::tibble(
        combination = combination_label,
        soc_count = rowSums(soc_binary_matrix)
    ) %>%
        count(combination, soc_count, name = "n", sort = TRUE) %>%
        mutate(
            RP = n / sum(n) * 100
        )
}

### 构建供UpSet绘图使用的SOC二值矩阵####
build_soc_upset_matrix <- function(data) {
    soc_metadata <- extract_soc_metadata(data)
    upset_matrix <- as.matrix(data[data$has_irae == TRUE, soc_metadata$soc_var, drop = FALSE] == 1)
    colnames(upset_matrix) <- soc_metadata$SOC
    storage.mode(upset_matrix) <- "logical"
    upset_matrix
}

### 汇总年度Any-irAE与Top SOC报告趋势####
build_annual_reporting_trend <- function(
    data,
    soc_overall,
    year_min = 2014,
    year_max = NULL,
    top_n_soc = 5
) {
    if (!"quarter" %in% names(data)) {
        return(list(
            any_irae = tibble::tibble(),
            top_soc = tibble::tibble(),
            combined = tibble::tibble()
        ))
    }

    data_with_year <- data %>%
        mutate(year = extract_report_year(quarter)) %>%
        filter(!is.na(year))

    if (nrow(data_with_year) == 0) {
        return(list(
            any_irae = tibble::tibble(),
            top_soc = tibble::tibble(),
            combined = tibble::tibble()
        ))
    }

    min_available_year <- min(data_with_year$year, na.rm = TRUE)
    max_available_year <- max(data_with_year$year, na.rm = TRUE)
    year_lower <- max(year_min, min_available_year)
    year_upper <- if (is.null(year_max)) max_available_year else min(year_max, max_available_year)
    year_sequence <- seq(year_lower, year_upper)

    data_with_year <- data_with_year %>%
        filter(year %in% year_sequence)

    year_denominator <- tibble::tibble(year = year_sequence) %>%
        left_join(
            data_with_year %>%
                count(year, name = "N_year"),
            by = "year"
        ) %>%
        mutate(N_year = tidyr::replace_na(N_year, 0L))

    any_irae_trend <- tibble::tibble(year = year_sequence) %>%
        left_join(
            data_with_year %>%
                group_by(year) %>%
                summarise(
                    n = sum(has_irae == TRUE, na.rm = TRUE),
                    .groups = "drop"
                ),
            by = "year"
        ) %>%
        left_join(year_denominator, by = "year") %>%
        mutate(
            n = tidyr::replace_na(n, 0L),
            RP = dplyr::if_else(N_year > 0, n / N_year * 100, 0)
        ) %>%
        transmute(
            year = year,
            trend_type = "Any-irAE",
            series = "Any-irAE",
            n = n,
            N = N_year,
            RP = RP
        )

    selected_soc <- soc_overall %>%
        arrange(desc(RP), desc(n), SOC) %>%
        slice_head(n = top_n_soc) %>%
        pull(SOC)

    soc_metadata <- extract_soc_metadata(data_with_year) %>%
        filter(SOC %in% selected_soc)

    top_soc_trend_raw <- data_with_year %>%
        group_by(year) %>%
        summarise(
            across(
                all_of(soc_metadata$soc_var),
                ~ sum(.x == 1, na.rm = TRUE)
            ),
            .groups = "drop"
        ) %>%
        pivot_longer(
            cols = all_of(soc_metadata$soc_var),
            names_to = "soc_var",
            values_to = "n"
        )

    top_soc_trend <- tidyr::expand_grid(
        year = year_sequence,
        soc_var = soc_metadata$soc_var
    ) %>%
        left_join(top_soc_trend_raw, by = c("year", "soc_var")) %>%
        left_join(soc_metadata, by = "soc_var") %>%
        left_join(year_denominator, by = "year") %>%
        mutate(
            n = tidyr::replace_na(n, 0L),
            RP = dplyr::if_else(N_year > 0, n / N_year * 100, 0)
        ) %>%
        transmute(
            year = year,
            trend_type = "SOC",
            series = SOC,
            n = n,
            N = N_year,
            RP = RP,
            SOC_full = SOC_full
        ) %>%
        arrange(year, series)

    combined_trend <- bind_rows(
        any_irae_trend %>% mutate(SOC_full = NA_character_),
        top_soc_trend
    )

    list(
        any_irae = any_irae_trend,
        top_soc = top_soc_trend,
        combined = combined_trend
    )
}

### 构建一行校验结果####
build_validation_row <- function(check_name, status, details) {
    tibble::tibble(
        check = check_name,
        status = status,
        details = details
    )
}

### 校验发生谱汇总结果是否满足预设规则####
validate_reporting_summary <- function(
    filtered_data,
    any_irae_by_organ,
    soc_overall,
    soc_by_organ,
    pt_overall,
    pt_by_organ_selected,
    pt_metadata,
    existing_summary_path = NULL
) {
    validation_rows <- list()

    add_row <- function(check_name, status, details) {
        validation_rows[[length(validation_rows) + 1]] <<- build_validation_row(
            check_name = check_name,
            status = status,
            details = details
        )
    }

    raw_any_n <- sum(filtered_data$has_irae == TRUE, na.rm = TRUE)
    aggregated_any_n <- sum(any_irae_by_organ$n_any, na.rm = TRUE)
    add_row(
        "Any-irAE RP与has_irae原始统计一致",
        ifelse(raw_any_n == aggregated_any_n, "PASS", "FAIL"),
        paste0("raw_any_n=", raw_any_n, "; aggregated_any_n=", aggregated_any_n)
    )

    soc_rp_valid <- all(dplyr::between(soc_overall$RP, 0, 100), na.rm = TRUE) &&
        all(dplyr::between(soc_by_organ$RP, 0, 100), na.rm = TRUE)
    add_row(
        "SOC RP位于0到100之间",
        ifelse(soc_rp_valid, "PASS", "FAIL"),
        paste0("SOC overall max RP=", round(max(soc_overall$RP, na.rm = TRUE), 2))
    )

    pt_rp_valid <- all(dplyr::between(pt_overall$RP, 0, 100), na.rm = TRUE) &&
        all(dplyr::between(pt_by_organ_selected$RP, 0, 100), na.rm = TRUE)
    add_row(
        "PT RP位于0到100之间",
        ifelse(pt_rp_valid, "PASS", "FAIL"),
        paste0("PT overall max RP=", round(max(pt_overall$RP, na.rm = TRUE), 2))
    )

    soc_n_valid <- all(soc_by_organ$n <= soc_by_organ$N_organ, na.rm = TRUE)
    add_row(
        "SOC报告数不超过对应分母",
        ifelse(soc_n_valid, "PASS", "FAIL"),
        paste0("max n/N=", max(soc_by_organ$n / soc_by_organ$N_organ, na.rm = TRUE))
    )

    pt_n_valid <- all(pt_overall$n <= nrow(filtered_data), na.rm = TRUE) &&
        all(pt_by_organ_selected$n <= pt_by_organ_selected$N_organ, na.rm = TRUE)
    add_row(
        "PT报告数不超过对应分母",
        ifelse(pt_n_valid, "PASS", "FAIL"),
        paste0("max organ-level n/N=", max(pt_by_organ_selected$n / pt_by_organ_selected$N_organ, na.rm = TRUE))
    )

    mapped_pt_n <- sum(!is.na(pt_metadata$SOC) & pt_metadata$SOC != "")
    add_row(
        "PT到SOC映射完整",
        ifelse(mapped_pt_n == nrow(pt_metadata), "PASS", "FAIL"),
        paste0("mapped=", mapped_pt_n, "/", nrow(pt_metadata))
    )

    add_row(
        "多SOC/PT累计RP可超过100%的说明已纳入汇总",
        "PASS",
        paste0(
            "overall SOC RP sum=", round(sum(soc_overall$RP, na.rm = TRUE), 2),
            "; overall PT RP sum=", round(sum(pt_overall$RP, na.rm = TRUE), 2),
            "；说明文字已写入队列汇总"
        )
    )

    if (!is.null(existing_summary_path) && file.exists(existing_summary_path)) {
        existing_organ_counts <- readxl::read_xlsx(existing_summary_path, sheet = "器官样本量") %>%
            mutate(
                primary_organ = as.character(primary_organ),
                n = as.integer(n)
            ) %>%
            arrange(desc(n), primary_organ)

        current_organ_counts <- any_irae_by_organ %>%
            transmute(primary_organ = as.character(primary_organ), n = as.integer(N)) %>%
            arrange(desc(n), primary_organ)

        organ_counts_match <- identical(existing_organ_counts, current_organ_counts)
        add_row(
            "当前primary_organ样本量与3.1主分析一致",
            ifelse(organ_counts_match, "PASS", "FAIL"),
            paste0("existing_total=", sum(existing_organ_counts$n), "; current_total=", sum(current_organ_counts$n))
        )
    } else {
        add_row(
            "当前primary_organ样本量与3.1主分析一致",
            "SKIPPED",
            "未找到3.1主分析汇总工作簿，跳过对照"
        )
    }

    bind_rows(validation_rows)
}

### 绘制各器官Any-irAE报告比例lollipop图####
plot_any_irae_by_organ <- function(any_irae_by_organ_df) {
    organ_colors <- get_primary_organ_color_map(any_irae_by_organ_df$primary_organ)

    plot_df <- any_irae_by_organ_df %>%
        arrange(desc(N), primary_organ) %>%
        mutate(
            organ_label = primary_organ,
            organ_label = factor(organ_label, levels = rev(unique(organ_label))),
            label_text = paste0(format_case_n(n_any), "/", format_case_n(N))
        )

    ggplot(plot_df, aes(x = RP, y = organ_label)) +
        geom_segment(
            aes(x = 0, xend = RP, yend = organ_label),
            linewidth = 0.9,
            colour = "#B8B8B8"
        ) +
        geom_point(
            aes(fill = primary_organ),
            shape = 21,
            size = 3.8,
            stroke = 0.35,
            colour = "#2B2B2B"
        ) +
        geom_text(
            aes(label = label_text),
            hjust = -0.1,
            size = pt_to_geom_text_size(12),
            colour = "black"
        ) +
        scale_fill_manual(values = organ_colors, guide = "none") +
        scale_x_continuous(
            limits = c(0, max(plot_df$RP, na.rm = TRUE) * 1.18),
            labels = function(x) paste0(round(x, 1), "%"),
            expand = expansion(mult = c(0, 0.02))
        ) +
        labs(
            x = "Any-irAE reporting proportion",
            y = NULL
        ) +
        get_pan_cancer_theme("bw") +
        theme(
            panel.grid.major.y = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.line.x = element_line(colour = "black"),
            axis.ticks.y = element_blank()
        )
}

### 绘制SOC气泡热图####
plot_soc_bubble_heatmap <- function(soc_by_organ_df) {
    organ_order <- soc_by_organ_df %>%
        distinct(primary_organ, N_organ) %>%
        arrange(desc(N_organ), primary_organ) %>%
        mutate(organ_label = primary_organ)

    soc_levels <- soc_by_organ_df %>%
        group_by(SOC) %>%
        summarise(overall_n = sum(n), overall_RP = sum(n) / sum(N_organ) * 100, .groups = "drop") %>%
        arrange(desc(overall_RP), desc(overall_n), SOC) %>%
        pull(SOC)

    plot_df <- soc_by_organ_df %>%
        left_join(organ_order, by = c("primary_organ", "N_organ")) %>%
        mutate(
            organ_label = factor(organ_label, levels = rev(organ_order$organ_label)),
            SOC = factor(SOC, levels = soc_levels)
        )

    ggplot(plot_df, aes(x = SOC, y = organ_label)) +
        geom_point(
            aes(size = n, fill = RP),
            shape = 21,
            colour = "#2B2B2B",
            stroke = 0.3
        ) +
        scale_size_area(
            max_size = 12,
            breaks = scales::pretty_breaks(n = 4)
        ) +
        scale_fill_gradientn(
            colours = get_pan_cancer_gradient_palette("warm"),
            name = "RP (%)"
        ) +
        labs(
            x = "SOC",
            y = NULL,
            size = "Report count"
        ) +
        get_pan_cancer_theme("bw") +
        theme(
            panel.grid.major = element_line(colour = "#E7E7E7", linewidth = 0.25),
            panel.grid.minor = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "right"
        )
}

### 绘制按SOC分面的Top PT lollipop图####
plot_pt_lollipop <- function(pt_top_df) {
    plot_df <- pt_top_df %>%
        arrange(SOC, desc(RP), desc(n), PT) %>%
        group_by(SOC) %>%
        mutate(
            pt_key = paste0(PT, "___", SOC),
            label_offset = dplyr::if_else(
                max(RP, na.rm = TRUE) > 0,
                max(RP, na.rm = TRUE) * 0.06,
                0.05
            ),
            label_x = RP + label_offset
        ) %>%
        ungroup()

    plot_levels <- plot_df %>%
        arrange(SOC, RP, n, PT) %>%
        pull(pt_key) %>%
        unique()

    plot_df <- plot_df %>%
        mutate(
            pt_key = factor(pt_key, levels = plot_levels),
            label_text = paste0("n=", format_case_n(n))
        )

    ggplot(plot_df, aes(x = RP, y = pt_key)) +
        geom_segment(
            aes(x = 0, xend = RP, yend = pt_key),
            linewidth = 0.75,
            colour = "#BCBCBC"
        ) +
        geom_point(
            colour = "#8C2D04",
            fill = "#F16913",
            shape = 21,
            size = 3.2,
            stroke = 0.25
        ) +
        geom_text(
            aes(x = label_x, label = label_text),
            hjust = 0,
            size = pt_to_geom_text_size(12),
            colour = "black"
        ) +
        facet_wrap(~ SOC, scales = "free", ncol = 3) +
        scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
        scale_x_continuous(
            labels = function(x) paste0(round(x, 1), "%"),
            expand = expansion(mult = c(0, 0.08))
        ) +
        labs(
            x = "Overall PT reporting proportion",
            y = NULL
        ) +
        get_pan_cancer_theme("bw") +
        theme(
            panel.grid.major.y = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_rect(fill = "#F4E6D7", colour = "#D9B38C"),
            strip.text = element_text(face = "bold")
        )
}

### 绘制selected PT在各器官中的气泡图####
plot_pt_bubble_by_organ <- function(pt_by_organ_df) {
    organ_order <- pt_by_organ_df %>%
        distinct(primary_organ, N_organ) %>%
        arrange(desc(N_organ), primary_organ) %>%
        mutate(organ_label = primary_organ)

    pt_levels <- pt_by_organ_df %>%
        arrange(SOC, rank_within_soc, PT) %>%
        mutate(pt_key = paste0(PT, "___", SOC)) %>%
        pull(pt_key) %>%
        unique()

    plot_df <- pt_by_organ_df %>%
        left_join(organ_order, by = c("primary_organ", "N_organ")) %>%
        mutate(
            organ_label = factor(organ_label, levels = rev(organ_order$organ_label)),
            pt_key = factor(paste0(PT, "___", SOC), levels = pt_levels)
        )

    ggplot(plot_df, aes(x = pt_key, y = organ_label)) +
        geom_point(
            aes(size = n, fill = RP),
            shape = 21,
            colour = "#2B2B2B",
            stroke = 0.25
        ) +
        facet_wrap(~ SOC, scales = "free_x", ncol = 3) +
        scale_x_discrete(labels = function(x) sub("___.*$", "", x)) +
        scale_size_area(max_size = 10) +
        scale_fill_gradientn(
            colours = get_pan_cancer_gradient_palette("green"),
            name = "RP (%)"
        ) +
        labs(
            x = "PT",
            y = NULL,
            size = "Report count"
        ) +
        get_pan_cancer_theme("bw") +
        theme(
            panel.grid.major = element_line(colour = "#E7E7E7", linewidth = 0.25),
            panel.grid.minor = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            strip.background = element_rect(fill = "#E8F2E4", colour = "#B9D3B0"),
            strip.text = element_text(face = "bold")
        )
}

### 绘制SOC共报告UpSet图####
plot_soc_upset <- function(soc_binary_matrix, top_n_comb = 20) {
    if (is.null(soc_binary_matrix) || nrow(soc_binary_matrix) == 0 || ncol(soc_binary_matrix) == 0) {
        return(NULL)
    }

    comb_matrix <- ComplexHeatmap::make_comb_mat(soc_binary_matrix)
    comb_order <- order(ComplexHeatmap::comb_size(comb_matrix), decreasing = TRUE)

    if (length(comb_order) > top_n_comb) {
        comb_order <- comb_order[seq_len(top_n_comb)]
    }

    comb_matrix <- comb_matrix[, comb_order]

    ComplexHeatmap::UpSet(
        comb_matrix,
        top_annotation = ComplexHeatmap::upset_top_annotation(comb_matrix, add_numbers = TRUE),
        right_annotation = ComplexHeatmap::upset_right_annotation(comb_matrix),
        pt_size = grid::unit(2, "mm"),
        lwd = 1
    )
}

### 绘制年度Any-irAE与Top SOC趋势组合图####
plot_annual_reporting_trend <- function(annual_any_irae_df, annual_top_soc_df) {
    p_any <- ggplot(annual_any_irae_df, aes(x = year, y = RP)) +
        geom_line(linewidth = 0.9, colour = "#A63603") +
        geom_point(size = 2.4, colour = "#A63603") +
        geom_text(
            aes(label = paste0("n=", format_case_n(n))),
            vjust = -0.8,
            size = pt_to_geom_text_size(12),
            colour = "black"
        ) +
        scale_x_continuous(breaks = annual_any_irae_df$year) +
        scale_y_continuous(labels = function(x) paste0(round(x, 1), "%")) +
        labs(
            x = NULL,
            y = "Any-irAE RP",
            title = "Annual Any-irAE reporting proportion"
        ) +
        get_pan_cancer_theme("bw") +
        theme(
            panel.grid.minor = element_blank()
        )

    p_soc <- ggplot(annual_top_soc_df, aes(x = year, y = RP, colour = series)) +
        geom_line(linewidth = 0.85) +
        geom_point(size = 2.1) +
        scale_x_continuous(breaks = sort(unique(annual_top_soc_df$year))) +
        scale_y_continuous(labels = function(x) paste0(round(x, 1), "%")) +
        scale_colour_manual(
            values = c(
                "SKIN" = "#D7301F",
                "ENDO" = "#238B45",
                "GI" = "#DD8D29",
                "NEURO" = "#2171B5",
                "MSK" = "#6A51A3",
                "HEPA" = "#7F2704",
                "HEME" = "#54278F",
                "RESP" = "#08519C",
                "RENU" = "#006D2C",
                "CARD" = "#B30000",
                "EYE" = "#756BB1"
            ),
            name = "SOC"
        ) +
        labs(
            x = "Year",
            y = "SOC RP",
            title = "Annual reporting proportion of overall top SOCs"
        ) +
        get_pan_cancer_theme("bw") +
        theme(
            panel.grid.minor = element_blank(),
            legend.position = "bottom"
        )

    p_any / p_soc + patchwork::plot_layout(heights = c(1, 1.2))
}
