# 通用配置文件####
# 用途：统一管理项目路径、R包加载和绑图设置，支持本地和服务器通用
# 使用方法：在每个脚本开头添加 source("Code/config.R")
# 注意：请确保在"泛癌种"文件夹下运行代码

## R镜像设置####

# CRAN镜像设置（使用国内镜像加速包下载）
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

# Bioconductor镜像设置
options(BioC_mirror = "https://mirrors.westlake.edu.cn/bioconductor")

# 下载设置
options(download.file.method = "libcurl") # 解决部分SSL问题
options(timeout = 300) # 增加下载超时时间（秒）

## 加载R包####

# 定义需要加载的包列表
packages_to_load <- c(
  # 核心数据处理
  "tidyverse", # 数据处理核心包（包含ggplot2, dplyr, tidyr等）
  "data.table", # 高效数据处理
  "lubridate", # 日期时间处理

  # 并行计算
  "future.apply", # 并行计算

  # 数据读写
  "readxl", # 读取Excel文件
  "writexl", # 写入Excel文件
  "openxlsx", # Excel文件操作（含格式化）

  # 统计分析
  "tableone" # 描述性统计表
)

# 可视化包（单独列出，部分可能需要手动安装）
viz_packages <- c(
  "ggplot2", # tidyverse已包含，此处为明确
  "patchwork", # 图形组合
  "ggrepel", # 标签避让
  "ggsankey", # 桑基图
  "ComplexHeatmap", # 复杂热图
  "circlize", # 圆形可视化和颜色映射
  "scales" # 坐标轴刻度格式化
)

# 合并所有包
all_packages <- c(packages_to_load, viz_packages)

# 加载包的函数
load_packages <- function(packages, quietly = TRUE) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste0("包 '", pkg, "' 未安装，尝试安装..."))
      # 尝试从CRAN安装
      tryCatch(
        {
          install.packages(pkg, quiet = TRUE)
        },
        error = function(e) {
          # 如果CRAN安装失败，尝试Bioconductor（如ComplexHeatmap）
          if (pkg %in% c("ComplexHeatmap", "circlize")) {
            if (!requireNamespace("BiocManager", quietly = TRUE)) {
              install.packages("BiocManager", quiet = TRUE)
            }
            BiocManager::install(pkg, ask = FALSE, update = FALSE)
          } else {
            message(paste0("无法安装包 '", pkg, "': ", e$message))
          }
        }
      )
    }

    # 加载包
    tryCatch(
      {
        suppressPackageStartupMessages(library(pkg, character.only = TRUE))
      },
      error = function(e) {
        message(paste0("无法加载包 '", pkg, "': ", e$message))
      }
    )
  }
}

# 加载所有包
message("=== 加载R包 ===")
load_packages(all_packages)
message("R包加载完成")
message("===============")

## 工作目录检查####

# 检查并设置工作目录
# 如果当前目录不是项目根目录（泛癌种），尝试自动检测
check_and_set_wd <- function() {
  current_wd <- getwd()

  # 检查当前目录是否包含必要的子文件夹
  required_dirs <- c("Code", "Data", "Results")
  has_required_dirs <- all(dir.exists(required_dirs))

  if (!has_required_dirs) {
    # 尝试查找项目根目录
    # 检查是否在Code子目录中
    if (basename(current_wd) == "Code" && dir.exists("../Data")) {
      setwd("..")
      message("已自动切换到项目根目录: ", getwd())
    } else {
      warning("当前目录可能不是项目根目录，请确保在'泛癌种'文件夹下运行代码")
      message("当前工作目录: ", current_wd)
    }
  }
}

# 运行目录检查
check_and_set_wd()

## 路径设置（使用相对路径）####
# 注意：使用以.开头的隐藏变量，防止被rm(list = ls())清除
# 同时提供普通变量别名方便使用

# 项目根目录（当前工作目录）
.base_path <- getwd()

# 数据目录
.data_path <- file.path(.base_path, "Data")

# 清洗后数据输出目录
.output_path <- file.path(.data_path, "Cleaned")

# 结果输出目录
.results_path <- file.path(.base_path, "Results")

# 创建普通变量别名（方便使用，但会被rm清除）
base_path <- .base_path
data_path <- .data_path
output_path <- .output_path
results_path <- .results_path

## 创建必要的目录（如果不存在）####

# 创建输出目录
if (!dir.exists(.output_path)) {
  dir.create(.output_path, recursive = TRUE)
  message("已创建输出目录: ", .output_path)
}

# 创建结果目录
if (!dir.exists(.results_path)) {
  dir.create(.results_path, recursive = TRUE)
  message("已创建结果目录: ", .results_path)
}

## 输出路径信息（用于确认）####

message("=== 路径配置信息 ===")
message("项目根目录: ", .base_path)
message("数据目录: ", .data_path)
message("清洗数据目录: ", .output_path)
message("结果目录: ", .results_path)
message("====================")

## 路径恢复函数####
# 如果执行了rm(list = ls())后需要恢复路径变量，调用此函数
.restore_paths <- function() {
  base_path <<- .base_path
  data_path <<- .data_path
  output_path <<- .output_path
  results_path <<- .results_path
  message("路径变量已恢复")
}

## 绑图主题设置####

# 加载ggplot2（如果尚未加载）
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  message("ggplot2 未安装，跳过主题设置")
} else {
  library(ggplot2)

  # 统一绘图字体与字号参数
  pan_cancer_plot_font_family <- ""
  pan_cancer_plot_text_size <- 12
  pan_cancer_plot_axis_size <- 14
  pan_cancer_plot_title_size <- 16

  ### 将pt字号转换为geom_text使用的size####
  pt_to_geom_text_size <- function(pt_size) {
    pt_size / ggplot2::.pt
  }

  ### 返回项目统一的ggplot主题####
  get_pan_cancer_theme <- function(theme_name = c("bw", "classic", "minimal")) {
    theme_name <- match.arg(theme_name)

    base_theme <- switch(
      theme_name,
      "bw" = ggplot2::theme_bw(base_family = pan_cancer_plot_font_family),
      "classic" = ggplot2::theme_classic(base_family = pan_cancer_plot_font_family),
      "minimal" = ggplot2::theme_minimal(base_family = pan_cancer_plot_font_family)
    )

    base_theme +
      ggplot2::theme(
        text = element_text(
          family = pan_cancer_plot_font_family,
          size = pan_cancer_plot_text_size,
          color = "black"
        ),
        axis.title = element_text(
          family = pan_cancer_plot_font_family,
          size = pan_cancer_plot_axis_size,
          color = "black"
        ),
        axis.text = element_text(
          family = pan_cancer_plot_font_family,
          size = pan_cancer_plot_axis_size,
          color = "black"
        ),
        legend.title = element_text(
          family = pan_cancer_plot_font_family,
          size = pan_cancer_plot_axis_size,
          color = "black"
        ),
        legend.text = element_text(
          family = pan_cancer_plot_font_family,
          size = pan_cancer_plot_axis_size,
          color = "black"
        ),
        plot.title = element_text(
          family = pan_cancer_plot_font_family,
          size = pan_cancer_plot_title_size,
          color = "black",
          face = "bold"
        ),
        plot.subtitle = element_text(
          family = pan_cancer_plot_font_family,
          size = pan_cancer_plot_text_size,
          color = "black"
        ),
        strip.text = element_text(
          family = pan_cancer_plot_font_family,
          size = pan_cancer_plot_text_size,
          color = "black",
          face = "bold"
        ),
        plot.caption = element_text(
          family = pan_cancer_plot_font_family,
          size = pan_cancer_plot_text_size,
          color = "black"
        )
      )
  }

  ### 返回年度ICI治疗策略统一配色####
  get_ici_strategy_color_map <- function() {
    c(
      "anti-PD-1" = "#2F5D8A",
      "anti-PD-L1" = "#5B8FB9",
      "anti-CTLA-4" = "#C97B2A",
      "Combination Therapy" = "#5B8F72"
    )
  }

  ### 返回项目统一的连续型配色####
  get_pan_cancer_gradient_palette <- function(
    palette_name = c("warm", "green", "blue", "neutral")
  ) {
    palette_name <- match.arg(palette_name)

    switch(
      palette_name,
      "warm" = c("#FFF4E6", "#F8C27A", "#E67E4A", "#A64B2A"),
      "green" = c("#F4FBF6", "#B8E0C2", "#59A96A", "#1F6B45"),
      "blue" = c("#F3F7FB", "#A9C4E2", "#5D89B3", "#1F4E79"),
      "neutral" = c("#F5F5F5", "#D0D7DE", "#7B8794", "#334E68")
    )
  }

  # 设置全局主题，使用统一字体和字号规范
  theme_set(get_pan_cancer_theme("bw"))

  message("=== 绑图主题配置 ===")
  message("字体: 设备默认字体")
  message("普通文字: ", pan_cancer_plot_text_size, "号，黑色")
  message("坐标轴/图例文字: ", pan_cancer_plot_axis_size, "号，黑色")
  message("====================")
}

## 清理环境（可选，移除辅助函数）####
rm(check_and_set_wd)
