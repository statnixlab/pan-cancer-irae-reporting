required_packages <- c(
  "dplyr",
  "future",
  "future.apply",
  "limma",
  "readr",
  "tidyr",
  "tibble",
  "TCGAbiolinks",
  "SummarizedExperiment",
  "data.table"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Missing required packages: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages(library(dplyr))

options(timeout = 1000, download.file.method = "libcurl")

# 在项目根目录下运行
base_path <- getwd()
output_root <- file.path(base_path, "Data", "TCGA_GTEx")
cache_root <- file.path(output_root, "cache")
project_dir <- file.path(output_root, "by_project")
n_workers <- min(4L, max(1L, parallel::detectCores() - 1L))
overwrite <- FALSE
download_method <- "api"
default_tcga_solid_projects <- c(
  "TCGA-ACC", "TCGA-BLCA", "TCGA-BRCA", "TCGA-CESC", "TCGA-CHOL",
  "TCGA-COAD", "TCGA-ESCA", "TCGA-GBM", "TCGA-HNSC", "TCGA-KICH",
  "TCGA-KIRC", "TCGA-KIRP", "TCGA-LGG", "TCGA-LIHC", "TCGA-LUAD",
  "TCGA-LUSC", "TCGA-MESO", "TCGA-OV", "TCGA-PAAD", "TCGA-PCPG",
  "TCGA-PRAD", "TCGA-READ", "TCGA-SARC", "TCGA-SKCM", "TCGA-STAD",
  "TCGA-TGCT", "TCGA-THCA", "TCGA-THYM", "TCGA-UCEC", "TCGA-UCS",
  "TCGA-UVM"
)

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(path)
}

read_csv_quiet <- function(path) {
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}

read_csv_character <- function(path) {
  readr::read_csv(
    path,
    show_col_types = FALSE,
    progress = FALSE,
    col_types = readr::cols(.default = readr::col_character())
  )
}

read_tsv_quiet <- function(path) {
  readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
}

read_tsv_character <- function(path) {
  readr::read_tsv(
    path,
    show_col_types = FALSE,
    progress = FALSE,
    col_types = readr::cols(.default = readr::col_character())
  )
}

write_csv_gz <- function(x, path) {
  ensure_dir(dirname(path))
  readr::write_csv(x, path)
}

write_tsv_gz <- function(x, path) {
  ensure_dir(dirname(path))
  readr::write_tsv(x, path)
}

get_default_projects <- function() {
  tryCatch(
    {
      all_projects <- TCGAbiolinks::getGDCprojects()$project_id
      tcga_projects <- sort(all_projects[grepl("^TCGA-", all_projects)])
      setdiff(tcga_projects, c("TCGA-LAML", "TCGA-DLBC"))
    },
    error = function(e) {
      message("Falling back to built-in TCGA solid tumor project list: ", conditionMessage(e))
      default_tcga_solid_projects
    }
  )
}

parse_project_env <- function(default_projects) {
  env_projects <- trimws(Sys.getenv("TCGA_PROJECTS", unset = ""))
  if (!nzchar(env_projects)) {
    return(default_projects)
  }

  projects <- trimws(strsplit(env_projects, ",", fixed = TRUE)[[1]])
  projects <- projects[nzchar(projects)]

  if (length(projects) == 0L) {
    stop("`TCGA_PROJECTS` was set but no valid project IDs were provided.", call. = FALSE)
  }

  invalid_projects <- setdiff(projects, default_projects)
  if (length(invalid_projects) > 0L) {
    stop(
      "Invalid or non-solid-tumor projects in `TCGA_PROJECTS`: ",
      paste(invalid_projects, collapse = ", "),
      call. = FALSE
    )
  }

  unique(projects)
}

parse_worker_env <- function(default_workers) {
  env_workers <- suppressWarnings(as.integer(Sys.getenv("TCGA_N_WORKERS", unset = "")))
  if (is.na(env_workers) || env_workers < 1L) {
    return(default_workers)
  }
  env_workers
}

parse_overwrite_env <- function(default_value) {
  env_overwrite <- tolower(trimws(Sys.getenv("TCGA_OVERWRITE", unset = "")))
  if (!nzchar(env_overwrite)) {
    return(default_value)
  }
  env_overwrite %in% c("1", "true", "yes", "y")
}

parse_cache_root_env <- function(default_cache_root) {
  env_cache_root <- trimws(Sys.getenv("TCGA_CACHE_ROOT", unset = ""))
  if (!nzchar(env_cache_root)) {
    return(default_cache_root)
  }
  env_cache_root
}

resolve_gdc_client_bin <- function() {
  candidate_bins <- c(
    trimws(Sys.getenv("TCGA_GDC_CLIENT_BIN", unset = "")),
    Sys.which("gdc-client.exe"),
    Sys.which("gdc-client")
  )
  candidate_bins <- unique(candidate_bins[nzchar(candidate_bins)])

  for (candidate in candidate_bins) {
    if (file.exists(candidate)) {
      return(normalizePath(candidate, winslash = "/", mustWork = FALSE))
    }
  }

  ""
}

configure_gdc_client_env <- function() {
  gdc_client_bin <- resolve_gdc_client_bin()
  if (!nzchar(gdc_client_bin)) {
    return("")
  }

  gdc_client_dir <- normalizePath(dirname(gdc_client_bin), winslash = "/", mustWork = FALSE)
  current_path <- Sys.getenv("PATH", unset = "")
  path_parts <- strsplit(current_path, .Platform$path.sep, fixed = TRUE)[[1]]

  if (!gdc_client_dir %in% path_parts) {
    Sys.setenv(PATH = paste(c(gdc_client_dir, current_path), collapse = .Platform$path.sep))
  }

  gdc_client_bin
}

parse_download_method_env <- function(default_method) {
  env_method <- tolower(trimws(Sys.getenv("TCGA_DOWNLOAD_METHOD", unset = "")))
  if (!nzchar(env_method)) {
    return(default_method)
  }

  if (!env_method %in% c("api", "client")) {
    stop("`TCGA_DOWNLOAD_METHOD` must be either 'api' or 'client'.", call. = FALSE)
  }

  gdc_client_bin <- configure_gdc_client_env()

  if (identical(env_method, "client") && !nzchar(gdc_client_bin)) {
    stop(
      "`TCGA_DOWNLOAD_METHOD=client` was requested, but `gdc-client` was not found.",
      call. = FALSE
    )
  }

  env_method
}

parse_proxy_env <- function() {
  first_nonempty_env <- function(keys) {
    for (key in keys) {
      value <- trimws(Sys.getenv(key, unset = ""))
      if (nzchar(value)) {
        return(value)
      }
    }
    ""
  }

  proxy_value <- first_nonempty_env(c("TCGA_PROXY"))
  http_proxy <- first_nonempty_env(c("TCGA_HTTP_PROXY", "http_proxy", "HTTP_PROXY"))
  https_proxy <- first_nonempty_env(c("TCGA_HTTPS_PROXY", "https_proxy", "HTTPS_PROXY"))
  no_proxy <- first_nonempty_env(c("TCGA_NO_PROXY", "no_proxy", "NO_PROXY"))

  if (!nzchar(http_proxy) && nzchar(proxy_value)) {
    http_proxy <- proxy_value
  }
  if (!nzchar(https_proxy) && nzchar(proxy_value)) {
    https_proxy <- proxy_value
  }

  list(
    http_proxy = http_proxy,
    https_proxy = https_proxy,
    no_proxy = no_proxy
  )
}

configure_proxy_env <- function(proxy_config) {
  if (nzchar(proxy_config$http_proxy)) {
    Sys.setenv(
      http_proxy = proxy_config$http_proxy,
      HTTP_PROXY = proxy_config$http_proxy
    )
  }

  if (nzchar(proxy_config$https_proxy)) {
    Sys.setenv(
      https_proxy = proxy_config$https_proxy,
      HTTPS_PROXY = proxy_config$https_proxy
    )
  }

  if (nzchar(proxy_config$no_proxy)) {
    Sys.setenv(
      no_proxy = proxy_config$no_proxy,
      NO_PROXY = proxy_config$no_proxy
    )
  }

  invisible(proxy_config)
}

project_files <- function(project_id) {
  prefix <- file.path(project_dir, project_id)
  list(
    cache_dir = file.path(cache_root, project_id),
    expr = paste0(prefix, "_tpm_matrix.rds"),
    anno = paste0(prefix, "_sample_annotation.csv.gz"),
    clinical = paste0(prefix, "_clinical_patient.csv.gz"),
    mutation = paste0(prefix, "_mutation.maf.gz"),
    gene = paste0(prefix, "_gene_annotation.csv.gz"),
    summary = paste0(prefix, "_summary.rds")
  )
}

has_project_outputs <- function(files) {
  required_paths <- unlist(files[c("expr", "anno", "clinical", "mutation", "gene", "summary")])
  all(file.exists(required_paths))
}

sample_code_from_barcode <- function(barcode) {
  substr(barcode, 14L, 15L)
}

standardize_gene_annotation <- function(se_object) {
  assay_names <- SummarizedExperiment::assayNames(se_object)
  if (!"tpm_unstrand" %in% assay_names) {
    stop("`tpm_unstrand` assay was not found in expression object.")
  }

  expr_matrix <- as.matrix(SummarizedExperiment::assay(se_object, "tpm_unstrand"))
  gene_info <- as.data.frame(SummarizedExperiment::rowData(se_object), stringsAsFactors = FALSE)

  gene_name <- if ("gene_name" %in% names(gene_info)) {
    as.character(gene_info$gene_name)
  } else {
    rep(NA_character_, nrow(expr_matrix))
  }

  gene_id <- if ("gene_id" %in% names(gene_info)) {
    as.character(gene_info$gene_id)
  } else {
    rownames(expr_matrix)
  }

  missing_gene_name <- is.na(gene_name) | !nzchar(gene_name) | gene_name == "?"
  gene_name[missing_gene_name] <- gene_id[missing_gene_name]

  rownames(expr_matrix) <- gene_name
  expr_matrix <- limma::avereps(expr_matrix, ID = rownames(expr_matrix))

  gene_map <- tibble::tibble(
    gene_name = gene_name,
    gene_id = gene_id
  ) %>%
    dplyr::group_by(gene_name) %>%
    dplyr::summarise(
      gene_id = paste(
        sort(unique(gene_id[!is.na(gene_id) & nzchar(gene_id)])),
        collapse = ";"
      ),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      gene_id = dplyr::if_else(!nzchar(gene_id), gene_name, gene_id)
    )

  final_gene_annotation <- tibble::tibble(gene_name = rownames(expr_matrix)) %>%
    dplyr::left_join(gene_map, by = "gene_name") %>%
    dplyr::mutate(
      gene_id = dplyr::if_else(is.na(gene_id) | !nzchar(gene_id), gene_name, gene_id)
    )

  list(
    expr_matrix = expr_matrix,
    gene_annotation = final_gene_annotation
  )
}

run_gdc_download <- function(query, cache_dir, download_method) {
  if (identical(download_method, "client")) {
    client_workdir <- file.path(
      cache_dir,
      sprintf(".gdc_client_work_%s_%s", Sys.getpid(), format(Sys.time(), "%Y%m%d%H%M%S"))
    )
    ensure_dir(client_workdir)

    old_wd <- getwd()
    on.exit({
      setwd(old_wd)
      if (dir.exists(client_workdir)) {
        unlink(client_workdir, recursive = TRUE, force = TRUE)
      }
    }, add = TRUE)

    setwd(client_workdir)
    TCGAbiolinks::GDCdownload(
      query = query,
      method = download_method,
      files.per.chunk = 1L,
      directory = normalizePath(cache_dir, winslash = "/", mustWork = FALSE)
    )
  } else {
    TCGAbiolinks::GDCdownload(
      query = query,
      method = download_method,
      files.per.chunk = 1L,
      directory = cache_dir
    )
  }
}

build_sample_annotation <- function(project_id, sample_barcodes) {
  sample_tbl <- tibble::tibble(
    project_id = project_id,
    sample_barcode = sample_barcodes,
    sample_id = substr(sample_barcodes, 1L, 16L),
    patient_id = substr(sample_barcodes, 1L, 12L),
    sample_code = sample_code_from_barcode(sample_barcodes),
    original_order = seq_along(sample_barcodes)
  ) %>%
    dplyr::mutate(
      sample_type = dplyr::case_when(
        sample_code == "01" ~ "Primary Tumor",
        sample_code == "11" ~ "Solid Tissue Normal",
        TRUE ~ NA_character_
      ),
      tissue_type = dplyr::case_when(
        sample_type == "Primary Tumor" ~ "Tumor",
        sample_type == "Solid Tissue Normal" ~ "Normal",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(sample_type))

  tumor_tbl_raw <- sample_tbl %>%
    dplyr::filter(sample_type == "Primary Tumor", sample_code == "01") %>%
    dplyr::arrange(patient_id, sample_barcode) %>%
    dplyr::group_by(patient_id) %>%
    dplyr::mutate(sample_rank_within_patient = dplyr::row_number()) %>%
    dplyr::ungroup()

  tumor_tbl_final <- tumor_tbl_raw %>%
    dplyr::filter(sample_rank_within_patient == 1L)

  normal_tbl <- sample_tbl %>%
    dplyr::filter(sample_type == "Solid Tissue Normal") %>%
    dplyr::arrange(patient_id, sample_barcode) %>%
    dplyr::group_by(patient_id) %>%
    dplyr::mutate(sample_rank_within_patient = dplyr::row_number()) %>%
    dplyr::ungroup()

  final_tbl <- dplyr::bind_rows(tumor_tbl_final, normal_tbl) %>%
    dplyr::arrange(original_order) %>%
    dplyr::select(
      project_id,
      sample_barcode,
      sample_id,
      patient_id,
      sample_type,
      tissue_type,
      sample_rank_within_patient
    )

  list(
    sample_annotation = final_tbl,
    raw_tumor_samples = nrow(tumor_tbl_raw),
    retained_tumor_samples = nrow(tumor_tbl_final),
    normal_samples = nrow(normal_tbl),
    selected_tumor_barcodes = tumor_tbl_final$sample_barcode
  )
}

download_expression_data <- function(project_id, cache_dir, download_method = download_method) {
  query <- TCGAbiolinks::GDCquery(
    project = project_id,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    sample.type = c("Primary Tumor", "Solid Tissue Normal")
  )

  query_results <- query[["results"]][[1]]
  if (nrow(query_results) == 0L) {
    stop("No expression files returned for ", project_id, ".")
  }

  run_gdc_download(
    query = query,
    cache_dir = cache_dir,
    download_method = download_method
  )

  se_object <- TCGAbiolinks::GDCprepare(
    query = query,
    save = FALSE,
    directory = cache_dir
  )

  expr_processed <- standardize_gene_annotation(se_object)
  annotation_processed <- build_sample_annotation(project_id, colnames(expr_processed$expr_matrix))

  expr_matrix <- expr_processed$expr_matrix[
    ,
    annotation_processed$sample_annotation$sample_barcode,
    drop = FALSE
  ]

  if (!identical(colnames(expr_matrix), annotation_processed$sample_annotation$sample_barcode)) {
    stop("Expression matrix columns and sample annotation order do not match for ", project_id, ".")
  }

  list(
    expr_matrix = expr_matrix,
    sample_annotation = annotation_processed$sample_annotation,
    gene_annotation = expr_processed$gene_annotation,
    raw_tumor_samples = annotation_processed$raw_tumor_samples,
    retained_tumor_samples = annotation_processed$retained_tumor_samples,
    normal_samples = annotation_processed$normal_samples,
    selected_tumor_barcodes = annotation_processed$selected_tumor_barcodes
  )
}

download_clinical_data <- function(project_id, patient_ids) {
  clinical_tbl <- TCGAbiolinks::GDCquery_clinic(
    project = project_id,
    type = "clinical",
    save.csv = FALSE
  )

  clinical_tbl <- as.data.frame(clinical_tbl, stringsAsFactors = FALSE)
  clinical_tbl$project_id <- project_id
  clinical_tbl$patient_id <- NA_character_

  if ("submitter_id" %in% names(clinical_tbl)) {
    clinical_tbl$patient_id <- clinical_tbl$submitter_id
  }
  if ("bcr_patient_barcode" %in% names(clinical_tbl)) {
    missing_patient_id <- is.na(clinical_tbl$patient_id) | !nzchar(clinical_tbl$patient_id)
    clinical_tbl$patient_id[missing_patient_id] <- clinical_tbl$bcr_patient_barcode[missing_patient_id]
  }

  clinical_tbl <- clinical_tbl %>%
    dplyr::filter(!is.na(patient_id), patient_id %in% unique(patient_ids)) %>%
    dplyr::distinct(patient_id, .keep_all = TRUE)

  clinical_tbl
}

download_mutation_data <- function(
  project_id,
  cache_dir,
  selected_tumor_barcodes,
  download_method = download_method
) {
  selected_sample_ids <- unique(substr(selected_tumor_barcodes, 1L, 16L))

  query <- TCGAbiolinks::GDCquery(
    project = project_id,
    data.category = "Simple Nucleotide Variation",
    data.type = "Masked Somatic Mutation"
  )

  query_results <- query[["results"]][[1]]
  if (nrow(query_results) == 0L) {
    return(
      tibble::tibble(
        project_id = character(),
        patient_id = character(),
        sample_id = character()
      )
    )
  }

  run_gdc_download(
    query = query,
    cache_dir = cache_dir,
    download_method = download_method
  )

  mutation_tbl <- TCGAbiolinks::GDCprepare(
    query = query,
    directory = cache_dir
  )

  mutation_tbl <- tibble::as_tibble(mutation_tbl)
  if (!"Tumor_Sample_Barcode" %in% names(mutation_tbl)) {
    stop("`Tumor_Sample_Barcode` was not found in mutation data for ", project_id, ".")
  }

  mutation_tbl %>%
    dplyr::mutate(
      project_id = project_id,
      patient_id = substr(Tumor_Sample_Barcode, 1L, 12L),
      sample_id = substr(Tumor_Sample_Barcode, 1L, 16L)
    ) %>%
    dplyr::filter(sample_id %in% selected_sample_ids)
}

build_project_summary_from_disk <- function(project_id, status = "skipped", error_message = NA_character_) {
  files <- project_files(project_id)
  if (file.exists(files$summary)) {
    summary_tbl <- readRDS(files$summary)
    summary_tbl$status <- status
    summary_tbl$error_message <- error_message
    return(summary_tbl)
  }

  expr_matrix <- readRDS(files$expr)
  sample_annotation <- read_csv_quiet(files$anno)
  clinical_tbl <- read_csv_quiet(files$clinical)
  mutation_tbl <- read_tsv_quiet(files$mutation)

  tibble::tibble(
    project_id = project_id,
    status = status,
    raw_tumor_samples = sum(sample_annotation$sample_type == "Primary Tumor"),
    retained_tumor_samples = sum(sample_annotation$tissue_type == "Tumor"),
    normal_samples = sum(sample_annotation$tissue_type == "Normal"),
    clinical_patients = nrow(clinical_tbl),
    mutation_rows = nrow(mutation_tbl),
    expression_genes = nrow(expr_matrix),
    sample_total = ncol(expr_matrix),
    error_message = error_message
  )
}

process_project <- function(project_id, overwrite = FALSE) {
  data.table::setDTthreads(1L)
  files <- project_files(project_id)
  current_download_method <- download_method

  ensure_dir(project_dir)
  ensure_dir(cache_root)

  if (!overwrite && has_project_outputs(files)) {
    message("[", project_id, "] existing project outputs detected, skipping download.")
    return(build_project_summary_from_disk(project_id, status = "skipped"))
  }

  if (dir.exists(files$cache_dir)) {
    unlink(files$cache_dir, recursive = TRUE, force = TRUE)
  }
  ensure_dir(files$cache_dir)

  message("[", project_id, "] downloading expression, clinical, and mutation data.")

  tryCatch(
    {
      expression_result <- download_expression_data(
        project_id = project_id,
        cache_dir = files$cache_dir,
        download_method = current_download_method
      )
      clinical_tbl <- download_clinical_data(project_id, expression_result$sample_annotation$patient_id)
      mutation_tbl <- download_mutation_data(
        project_id = project_id,
        cache_dir = files$cache_dir,
        selected_tumor_barcodes = expression_result$selected_tumor_barcodes,
        download_method = current_download_method
      )

      saveRDS(expression_result$expr_matrix, files$expr)
      write_csv_gz(expression_result$sample_annotation, files$anno)
      write_csv_gz(clinical_tbl, files$clinical)
      write_tsv_gz(mutation_tbl, files$mutation)
      write_csv_gz(expression_result$gene_annotation, files$gene)

      summary_tbl <- tibble::tibble(
        project_id = project_id,
        status = "completed",
        raw_tumor_samples = expression_result$raw_tumor_samples,
        retained_tumor_samples = expression_result$retained_tumor_samples,
        normal_samples = expression_result$normal_samples,
        clinical_patients = nrow(clinical_tbl),
        mutation_rows = nrow(mutation_tbl),
        expression_genes = nrow(expression_result$expr_matrix),
        sample_total = ncol(expression_result$expr_matrix),
        error_message = NA_character_
      )

      saveRDS(summary_tbl, files$summary)
      summary_tbl
    },
    error = function(e) {
      message("[", project_id, "] failed: ", conditionMessage(e))
      tibble::tibble(
        project_id = project_id,
        status = "failed",
        raw_tumor_samples = NA_integer_,
        retained_tumor_samples = NA_integer_,
        normal_samples = NA_integer_,
        clinical_patients = NA_integer_,
        mutation_rows = NA_integer_,
        expression_genes = NA_integer_,
        sample_total = NA_integer_,
        error_message = conditionMessage(e)
      )
    }
  )
}

combine_expression_matrices <- function(success_projects) {
  expr_list <- lapply(success_projects, function(project_id) {
    readRDS(project_files(project_id)$expr)
  })
  names(expr_list) <- success_projects

  common_genes <- Reduce(intersect, lapply(expr_list, rownames))
  if (length(common_genes) == 0L) {
    stop("No common genes were found across successful projects.")
  }

  expr_list <- lapply(expr_list, function(expr_matrix) {
    expr_matrix[common_genes, , drop = FALSE]
  })

  tcga_tpm_matrix <- do.call(cbind, expr_list)

  gene_annotation <- lapply(success_projects, function(project_id) {
    read_csv_quiet(project_files(project_id)$gene)
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::group_by(gene_name) %>%
    dplyr::summarise(
      gene_id = paste(sort(unique(gene_id[!is.na(gene_id) & nzchar(gene_id)])), collapse = ";"),
      .groups = "drop"
    )

  tcga_gene_annotation <- tibble::tibble(gene_name = common_genes) %>%
    dplyr::left_join(gene_annotation, by = "gene_name") %>%
    dplyr::mutate(
      gene_id = dplyr::if_else(is.na(gene_id) | !nzchar(gene_id), gene_name, gene_id)
    )

  list(
    matrix = tcga_tpm_matrix,
    gene_annotation = tcga_gene_annotation
  )
}

combine_sample_annotation <- function(success_projects, matrix_colnames) {
  sample_annotation <- lapply(success_projects, function(project_id) {
    read_csv_quiet(project_files(project_id)$anno)
  }) %>%
    dplyr::bind_rows()

  sample_annotation <- sample_annotation[match(matrix_colnames, sample_annotation$sample_barcode), , drop = FALSE]
  if (!identical(matrix_colnames, sample_annotation$sample_barcode)) {
    stop("Combined sample annotation order does not match the expression matrix.")
  }

  sample_annotation
}

combine_clinical_data <- function(success_projects, sample_annotation) {
  patient_ids <- unique(sample_annotation$patient_id)

  clinical_tbl <- lapply(success_projects, function(project_id) {
    read_csv_character(project_files(project_id)$clinical)
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::filter(patient_id %in% patient_ids) %>%
    dplyr::distinct(patient_id, .keep_all = TRUE)

  clinical_tbl
}

combine_mutation_data <- function(success_projects) {
  mutation_list <- lapply(success_projects, function(project_id) {
    read_tsv_character(project_files(project_id)$mutation)
  })

  dplyr::bind_rows(mutation_list)
}

validate_outputs <- function(
  tcga_tpm_matrix,
  tcga_sample_annotation,
  tcga_clinical_patient,
  tcga_mutation,
  tcga_download_manifest,
  expected_projects
) {
  if (!identical(colnames(tcga_tpm_matrix), tcga_sample_annotation$sample_barcode)) {
    stop("Validation failed: matrix columns and sample annotation are misaligned.")
  }

  sample_types <- unique(tcga_sample_annotation$sample_type)
  if (!all(sample_types %in% c("Primary Tumor", "Solid Tissue Normal"))) {
    stop("Validation failed: unexpected sample types detected.")
  }

  tumor_sample_codes <- sample_code_from_barcode(
    tcga_sample_annotation$sample_barcode[tcga_sample_annotation$tissue_type == "Tumor"]
  )
  if (!all(tumor_sample_codes == "01")) {
    stop("Validation failed: non-01 tumor sample codes were retained.")
  }

  tumor_patient_counts <- tcga_sample_annotation %>%
    dplyr::filter(tissue_type == "Tumor") %>%
    dplyr::count(patient_id)
  if (any(tumor_patient_counts$n > 1L)) {
    stop("Validation failed: one or more patients retained multiple tumor samples.")
  }

  if (anyDuplicated(tcga_clinical_patient$patient_id) > 0L) {
    stop("Validation failed: clinical patient_id values are not unique.")
  }

  annotation_patients <- unique(tcga_sample_annotation$patient_id)
  if (!all(annotation_patients %in% tcga_clinical_patient$patient_id)) {
    missing_patients <- setdiff(annotation_patients, tcga_clinical_patient$patient_id)
    stop(
      "Validation failed: clinical data do not cover all retained patients. Missing: ",
      paste(head(missing_patients, 10L), collapse = ", ")
    )
  }

  if (nrow(tcga_mutation) > 0L) {
    tumor_sample_ids <- tcga_sample_annotation$sample_id[tcga_sample_annotation$tissue_type == "Tumor"]
    if (!all(tcga_mutation$sample_id %in% tumor_sample_ids)) {
      stop("Validation failed: mutation table contains tumor sample_id values outside the retained tumor set.")
    }
  }

  if (length(expected_projects) != 31L || any(expected_projects %in% c("TCGA-LAML", "TCGA-DLBC"))) {
    stop("Validation failed: default project set is incorrect.")
  }

  manifest_counts <- tcga_download_manifest %>%
    dplyr::filter(status %in% c("completed", "skipped")) %>%
    dplyr::select(project_id, retained_tumor_samples, normal_samples)

  annotation_counts <- tcga_sample_annotation %>%
    dplyr::count(project_id, tissue_type, name = "n") %>%
    tidyr::pivot_wider(
      names_from = tissue_type,
      values_from = n,
      values_fill = 0
    )

  if (!"Tumor" %in% names(annotation_counts)) {
    annotation_counts$Tumor <- 0L
  }
  if (!"Normal" %in% names(annotation_counts)) {
    annotation_counts$Normal <- 0L
  }

  annotation_counts <- annotation_counts %>%
    dplyr::rename(
      retained_tumor_samples = Tumor,
      normal_samples = Normal
    )

  count_check <- manifest_counts %>%
    dplyr::left_join(annotation_counts, by = "project_id", suffix = c("_manifest", "_annotation"))

  tumor_match <- identical(
    count_check$retained_tumor_samples_manifest,
    count_check$retained_tumor_samples_annotation
  )
  normal_match <- identical(
    count_check$normal_samples_manifest,
    count_check$normal_samples_annotation
  )
  if (!tumor_match || !normal_match) {
    stop("Validation failed: manifest sample counts do not match sample annotation counts.")
  }

  TRUE
}

ensure_dir(output_root)
ensure_dir(cache_root)
ensure_dir(project_dir)

main <- function() {
  default_projects <- get_default_projects()
  projects <- parse_project_env(default_projects)
  n_workers_local <- parse_worker_env(n_workers)
  overwrite_local <- parse_overwrite_env(overwrite)
  cache_root_local <- parse_cache_root_env(cache_root)
  gdc_client_bin <- configure_gdc_client_env()
  download_method_local <- parse_download_method_env(download_method)
  proxy_config <- parse_proxy_env()

  cache_root <<- cache_root_local
  download_method <<- download_method_local
  ensure_dir(cache_root)
  configure_proxy_env(proxy_config)

  message("TCGA solid tumor projects selected: ", length(projects))
  message("Worker count: ", n_workers_local)
  message("Overwrite existing outputs: ", overwrite_local)
  message("Cache root: ", cache_root_local)
  message("Download method: ", download_method_local)
  if (nzchar(gdc_client_bin)) {
    message("gdc-client: ", gdc_client_bin)
  }
  message("HTTP proxy: ", if (nzchar(proxy_config$http_proxy)) proxy_config$http_proxy else "<none>")
  message("HTTPS proxy: ", if (nzchar(proxy_config$https_proxy)) proxy_config$https_proxy else "<none>")
  if (nzchar(proxy_config$no_proxy)) {
    message("NO_PROXY: ", proxy_config$no_proxy)
  }

  future::plan(future::multisession, workers = n_workers_local)
  on.exit(future::plan(future::sequential), add = TRUE)

  project_summaries <- future.apply::future_lapply(
    projects,
    process_project,
    overwrite = overwrite_local,
    future.seed = TRUE
  )

  tcga_download_manifest <- dplyr::bind_rows(project_summaries) %>%
    dplyr::arrange(project_id)

  successful_projects <- tcga_download_manifest %>%
    dplyr::filter(status %in% c("completed", "skipped")) %>%
    dplyr::pull(project_id)

  if (length(successful_projects) == 0L) {
    write_csv_gz(tcga_download_manifest, file.path(output_root, "tcga_download_manifest.csv"))
    stop("No projects completed successfully. Manifest was written to disk.", call. = FALSE)
  }

  combined_expression <- combine_expression_matrices(successful_projects)
  tcga_tpm_matrix <- combined_expression$matrix
  tcga_gene_annotation <- combined_expression$gene_annotation
  tcga_sample_annotation <- combine_sample_annotation(successful_projects, colnames(tcga_tpm_matrix))
  tcga_clinical_patient <- combine_clinical_data(successful_projects, tcga_sample_annotation)
  tcga_mutation <- combine_mutation_data(successful_projects)

  validate_outputs(
    tcga_tpm_matrix = tcga_tpm_matrix,
    tcga_sample_annotation = tcga_sample_annotation,
    tcga_clinical_patient = tcga_clinical_patient,
    tcga_mutation = tcga_mutation,
    tcga_download_manifest = tcga_download_manifest,
    expected_projects = default_projects
  )

  saveRDS(tcga_tpm_matrix, file.path(output_root, "tcga_tpm_matrix.rds"))
  write_csv_gz(tcga_sample_annotation, file.path(output_root, "tcga_sample_annotation.csv.gz"))
  write_csv_gz(tcga_clinical_patient, file.path(output_root, "tcga_clinical_patient.csv.gz"))
  write_tsv_gz(tcga_mutation, file.path(output_root, "tcga_mutation.maf.gz"))
  write_csv_gz(tcga_gene_annotation, file.path(output_root, "tcga_gene_annotation.csv.gz"))
  write_csv_gz(tcga_download_manifest, file.path(output_root, "tcga_download_manifest.csv"))

  message("TCGA download pipeline completed successfully.")
  message("Successful projects: ", length(successful_projects), "/", length(projects))
  message("Final matrix dimension: ", nrow(tcga_tpm_matrix), " genes x ", ncol(tcga_tpm_matrix), " samples")
}

if (sys.nframe() == 0L) {
  main()
}
