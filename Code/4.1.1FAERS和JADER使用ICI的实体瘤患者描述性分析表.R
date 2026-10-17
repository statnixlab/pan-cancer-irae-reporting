# FAERS和JADER使用ICI的实体瘤患者描述性分析表####

rm(list = ls())
gc()

source("Code/config.R")

script_name <- "4.1.1FAERS和JADER使用ICI的实体瘤患者描述性分析表"
script_results_path <- file.path(results_path, script_name)
if (!dir.exists(script_results_path)) {
    dir.create(script_results_path, recursive = TRUE, showWarnings = FALSE)
}

message("========== FAERS和JADER使用ICI的实体瘤患者描述性分析表开始 ==========")
message("结果目录: ", script_results_path)


## 1. 参数设置####

faers_input_path <- file.path(output_path, "faers_solid.RData")
jader_input_path <- file.path(output_path, "jader_solid.RData")

table_output_path <- file.path(
    script_results_path,
    "FAERS和JADER使用ICI实体瘤患者描述性分析表.xlsx"
)
object_output_path <- file.path(
    script_results_path,
    "FAERS和JADER使用ICI实体瘤患者描述性分析对象.RData"
)


## 2. 辅助函数####

### 从RData文件读取指定对象####
# 功能：从指定RData文件中读取目标对象，并避免污染全局环境####
# 输入：input_path，RData文件路径；object_name，对象名称####
# 输出：返回读取后的R对象####
load_named_rdata <- function(input_path, object_name) {
    data_env <- new.env(parent = emptyenv())
    load(input_path, envir = data_env)

    if (!exists(object_name, envir = data_env, inherits = FALSE)) {
        stop(paste("RData中未找到对象:", object_name))
    }

    get(object_name, envir = data_env, inherits = FALSE)
}

### 将TableOne结果转换为数据框####
# 功能：将CreateTableOne生成的结果整理为可直接导出的数据框####
# 输入：tableone_object，CreateTableOne返回的对象####
# 输出：返回包含Variable列的data.frame####
tableone_to_dataframe <- function(tableone_object) {
    tableone_df <- print(
        tableone_object,
        showAllLevels = TRUE,
        quote = FALSE,
        noSpaces = TRUE,
        printToggle = FALSE
    )

    tableone_df <- as.data.frame(tableone_df)
    tableone_df <- cbind(Variable = rownames(tableone_df), tableone_df)
    rownames(tableone_df) <- NULL

    tableone_df
}

### 为单个数字文本添加千分位逗号####
# 功能：仅对绝对值大于等于1000的数字文本添加千分位分隔符####
# 输入：number_text，单个数字对应的字符文本####
# 输出：返回添加千分位逗号后的字符文本####
format_single_number_text <- function(number_text) {
    number_value <- suppressWarnings(as.numeric(number_text))

    if (is.na(number_value) || abs(number_value) < 1000) {
        return(number_text)
    }

    prettyNum(
        number_text,
        big.mark = ",",
        scientific = FALSE,
        preserve.width = "none"
    )
}

### 为字符向量中的数字添加千分位逗号####
# 功能：识别字符中的数字片段，并为其中的大数添加千分位逗号####
# 输入：text_vector，字符向量####
# 输出：返回格式化后的字符向量####
format_numbers_in_text <- function(text_vector) {
    vapply(text_vector, function(single_text) {
        if (is.na(single_text)) {
            return(NA_character_)
        }

        match_positions <- gregexpr("\\d+(?:\\.\\d+)?", single_text, perl = TRUE)
        matched_tokens <- regmatches(single_text, match_positions)[[1]]

        if (length(matched_tokens) == 0) {
            return(single_text)
        }

        formatted_tokens <- vapply(
            matched_tokens,
            format_single_number_text,
            character(1)
        )

        regmatches(single_text, match_positions) <- list(formatted_tokens)
        single_text
    }, character(1))
}

### 为导出的TableOne数据框添加千分位逗号####
# 功能：对TableOne结果中的所有字符列执行千分位格式化####
# 输入：data，待导出的TableOne数据框####
# 输出：返回格式化后的数据框####
format_tableone_export <- function(data) {
    data %>%
        dplyr::mutate(
            dplyr::across(
                .cols = everything(),
                .fns = ~ format_numbers_in_text(as.character(.x))
            )
        )
}

### 为队列筛选汇总表添加千分位逗号####
# 功能：对筛选汇总表中的样本量列添加千分位逗号，保留其他列原始结构####
# 输入：data，队列筛选汇总数据框####
# 输出：返回格式化后的数据框####
format_filter_summary_export <- function(data) {
    data %>%
        dplyr::mutate(
            n = format_numbers_in_text(as.character(n))
        )
}

### 统计FAERS前5报告国家并生成显示名称####
# 功能：为FAERS描述性表构建稳定的国家分组水平####
# 输入：data，FAERS ICI实体瘤病例数据框####
# 输出：返回包含top5_countries和top5_country_names的列表####
build_faers_top5_country_levels <- function(data) {
    country_counts <- data %>%
        dplyr::mutate(
            country_clean = dplyr::case_when(
                is.na(reporter_country) ~ "Unknown",
                reporter_country == "COUNTRY NOT SPECIFIED" ~ "Unknown",
                TRUE ~ reporter_country
            )
        ) %>%
        dplyr::count(country_clean, sort = TRUE)

    top5_countries <- country_counts %>%
        dplyr::filter(country_clean != "Unknown") %>%
        dplyr::slice_head(n = 5) %>%
        dplyr::pull(country_clean)

    country_name_map <- c(
        "US" = "United States",
        "JP" = "Japan",
        "FR" = "France",
        "CA" = "Canada",
        "DE" = "Germany",
        "GB" = "United Kingdom",
        "IT" = "Italy",
        "CN" = "China",
        "AU" = "Australia",
        "ES" = "Spain"
    )

    top5_country_names <- country_name_map[top5_countries]
    top5_country_names <- ifelse(
        is.na(top5_country_names),
        top5_countries,
        top5_country_names
    )

    list(
        top5_countries = top5_countries,
        top5_country_names = top5_country_names
    )
}

### 构建FAERS描述性分析输入表####
# 功能：基于FAERS ICI实体瘤队列生成TableOne输入数据####
# 输入：data，已筛选的FAERS ICI实体瘤病例数据框####
# 输出：返回包含table_data与vars的列表####
prepare_faers_tableone_input <- function(data) {
    country_levels <- build_faers_top5_country_levels(data)

    table_data <- data %>%
        dplyr::mutate(
            Age = dplyr::case_when(
                is.na(age_years) ~ "Unknown",
                age_years < 40 ~ "<40",
                age_years >= 40 & age_years < 60 ~ "40–<60",
                age_years >= 60 & age_years < 80 ~ "60–<80",
                age_years >= 80 ~ "≥80",
                TRUE ~ "Unknown"
            ),
            Age = factor(Age, levels = c("<40", "40–<60", "60–<80", "≥80", "Unknown")),
            Sex = factor(sex_std, levels = c("Male", "Female", "Unknown")),
            Weight = dplyr::case_when(
                is.na(wt) ~ "Unknown",
                wt < 50 ~ "<50",
                wt >= 50 & wt < 70 ~ "50–<70",
                wt >= 70 & wt < 90 ~ "70–<90",
                wt >= 90 ~ "≥90",
                TRUE ~ "Unknown"
            ),
            Weight = factor(Weight, levels = c("<50", "50–<70", "70–<90", "≥90", "Unknown")),
            Reporter = dplyr::case_when(
                is.na(occp_cod) ~ "Unknown",
                occp_cod %in% c("MD", "PH", "OT") ~ "Health-professional",
                TRUE ~ "Non-health-professional"
            ),
            Reporter = factor(
                Reporter,
                levels = c("Health-professional", "Non-health-professional", "Unknown")
            ),
            Country = dplyr::case_when(
                is.na(reporter_country) | reporter_country == "COUNTRY NOT SPECIFIED" ~ "Unknown",
                reporter_country %in% country_levels$top5_countries ~ reporter_country,
                TRUE ~ "Other Countries"
            ),
            Country = dplyr::case_when(
                Country == "US" ~ "United States",
                Country == "JP" ~ "Japan",
                Country == "FR" ~ "France",
                Country == "CA" ~ "Canada",
                Country == "DE" ~ "Germany",
                TRUE ~ Country
            ),
            Country = factor(
                Country,
                levels = unique(c(
                    country_levels$top5_country_names,
                    "Other Countries",
                    "Unknown"
                ))
            ),
            `ICI therapy` = ifelse(has_ici == TRUE, "Yes", "No"),
            `ICI therapy` = factor(`ICI therapy`, levels = c("Yes", "No")),
            Chemotherapy = factor(chemotherapy, levels = c("Yes", "No")),
            `Targeted therapy` = factor(targeted_therapy, levels = c("Yes", "No"))
        ) %>%
        dplyr::select(
            Age,
            Sex,
            Weight,
            Reporter,
            Country,
            `ICI therapy`,
            Chemotherapy,
            `Targeted therapy`
        )

    vars <- c(
        "Age", "Sex", "Weight", "Reporter", "Country",
        "ICI therapy", "Chemotherapy", "Targeted therapy"
    )

    list(
        table_data = table_data,
        vars = vars
    )
}

### 构建JADER描述性分析输入表####
# 功能：基于JADER ICI实体瘤队列生成TableOne输入数据####
# 输入：data，已筛选的JADER ICI实体瘤病例数据框####
# 输出：返回包含table_data与vars的列表####
prepare_jader_tableone_input <- function(data) {
    table_data <- data %>%
        dplyr::mutate(
            Age = dplyr::case_when(
                is.na(age) | age %in% c("不明", "") ~ "Unknown",
                age %in% c(
                    "10歳未満", "10歳代", "20歳代", "30歳代",
                    "乳幼児", "新生児", "小児", "胎児", "青少年", "成人"
                ) ~ "<40",
                age %in% c("40歳代", "50歳代") ~ "40–<60",
                age %in% c("60歳代", "70歳代") ~ "60–<80",
                age %in% c("80歳代", "90歳代", "100歳代", "高齢者") ~ "≥80",
                TRUE ~ "Unknown"
            ),
            Age = factor(Age, levels = c("<40", "40–<60", "60–<80", "≥80", "Unknown")),
            Sex = dplyr::case_when(
                is.na(sex) ~ "Unknown",
                sex %in% c("男性") ~ "Male",
                sex %in% c("女性") ~ "Female",
                TRUE ~ "Unknown"
            ),
            Sex = factor(Sex, levels = c("Male", "Female", "Unknown")),
            Weight = dplyr::case_when(
                is.na(weight) | weight == "" ~ "Unknown",
                weight %in% c("10kg未満", "10kg台", "20kg台", "30kg台", "40kg台") ~ "<50",
                weight %in% c("50kg台", "60kg台") ~ "50–<70",
                weight %in% c("70kg台", "80kg台") ~ "70–<90",
                weight %in% c(
                    "90kg台", "100kg台", "110kg台", "120kg台", "130kg台",
                    "140kg台", "150kg台", "160kg台", "170kg台", "180kg台",
                    "200kg以上"
                ) ~ "≥90",
                TRUE ~ "Unknown"
            ),
            Weight = factor(Weight, levels = c("<50", "50–<70", "70–<90", "≥90", "Unknown")),
            Reporter = dplyr::case_when(
                is.na(reporter_qual) | reporter_qual == "不明" ~ "Unknown",
                grepl("医師|薬剤師|医療関係者", reporter_qual) ~ "Health-professional",
                TRUE ~ "Non-health-professional"
            ),
            Reporter = factor(
                Reporter,
                levels = c("Health-professional", "Non-health-professional", "Unknown")
            ),
            Height = dplyr::case_when(
                is.na(height) | height == "" ~ "Unknown",
                height %in% c(
                    "10cm未満", "30cm台", "40cm台", "50cm台", "60cm台",
                    "70cm台", "80cm台", "90cm台", "100cm台", "110cm台",
                    "120cm台", "130cm台", "140cm台"
                ) ~ "<150",
                height %in% c("150cm台") ~ "150–<160",
                height %in% c("160cm台") ~ "160–<170",
                height %in% c("170cm台", "180cm台", "190cm台") ~ "≥170",
                TRUE ~ "Unknown"
            ),
            Height = factor(Height, levels = c("<150", "150–<160", "160–<170", "≥170", "Unknown")),
            `ICI therapy` = ifelse(has_ici == TRUE, "Yes", "No"),
            `ICI therapy` = factor(`ICI therapy`, levels = c("Yes", "No")),
            Chemotherapy = ifelse(has_chemotherapy == TRUE, "Yes", "No"),
            Chemotherapy = factor(Chemotherapy, levels = c("Yes", "No")),
            `Targeted therapy` = ifelse(has_targeted == TRUE, "Yes", "No"),
            `Targeted therapy` = factor(`Targeted therapy`, levels = c("Yes", "No"))
        ) %>%
        dplyr::select(
            Age,
            Sex,
            Weight,
            Height,
            Reporter,
            `ICI therapy`,
            Chemotherapy,
            `Targeted therapy`
        )

    vars <- c(
        "Age", "Sex", "Weight", "Height", "Reporter",
        "ICI therapy", "Chemotherapy", "Targeted therapy"
    )

    list(
        table_data = table_data,
        vars = vars
    )
}


## 3. 加载实体瘤分析集####

faers_solid <- load_named_rdata(
    input_path = faers_input_path,
    object_name = "faers_solid"
)

jader_solid <- load_named_rdata(
    input_path = jader_input_path,
    object_name = "jader_solid"
)

message("已加载 faers_solid: ", nrow(faers_solid), " 例")
message("已加载 jader_solid: ", nrow(jader_solid), " 例")


## 4. 生成ICI实体瘤主分析队列####

faers_ici_solid <- faers_solid %>%
    dplyr::filter(has_ici == TRUE)

jader_ici_solid <- jader_solid %>%
    dplyr::filter(has_ici == TRUE)

if (nrow(faers_ici_solid) == 0) {
    stop("faers_solid 中不存在 has_ici == TRUE 的病例")
}

if (nrow(jader_ici_solid) == 0) {
    stop("jader_solid 中不存在 has_ici == TRUE 的病例")
}

faers_filter_summary <- tibble::tibble(
    step = c("faers_solid", "faers_solid_ici"),
    n = c(nrow(faers_solid), nrow(faers_ici_solid))
) %>%
    dplyr::mutate(retained_pct = round(n / dplyr::first(n) * 100, 2))

jader_filter_summary <- tibble::tibble(
    step = c("jader_solid", "jader_solid_ici"),
    n = c(nrow(jader_solid), nrow(jader_ici_solid))
) %>%
    dplyr::mutate(retained_pct = round(n / dplyr::first(n) * 100, 2))

message("FAERS ICI实体瘤病例数: ", nrow(faers_ici_solid))
message("JADER ICI实体瘤病例数: ", nrow(jader_ici_solid))


## 5. 构建描述性分析表####

faers_table_input <- prepare_faers_tableone_input(faers_ici_solid)
faers_tableone <- tableone::CreateTableOne(
    vars = faers_table_input$vars,
    data = faers_table_input$table_data
)
faers_tableone_df <- tableone_to_dataframe(faers_tableone) %>%
    format_tableone_export()

jader_table_input <- prepare_jader_tableone_input(jader_ici_solid)
jader_tableone <- tableone::CreateTableOne(
    vars = jader_table_input$vars,
    data = jader_table_input$table_data
)
jader_tableone_df <- tableone_to_dataframe(jader_tableone) %>%
    format_tableone_export()


## 6. 保存结果表####

workbook_sheets <- list(
    "FAERS_队列筛选" = format_filter_summary_export(faers_filter_summary),
    "FAERS_TableOne" = faers_tableone_df,
    "JADER_队列筛选" = format_filter_summary_export(jader_filter_summary),
    "JADER_TableOne" = jader_tableone_df
)

writexl::write_xlsx(workbook_sheets, path = table_output_path)

save(
    faers_filter_summary,
    jader_filter_summary,
    faers_tableone_df,
    jader_tableone_df,
    file = object_output_path
)

message("描述性分析表已保存: ", table_output_path)
message("描述性分析对象已保存: ", object_output_path)
message("========== FAERS和JADER使用ICI的实体瘤患者描述性分析表完成 ==========")
