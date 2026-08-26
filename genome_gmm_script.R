#!/usr/bin/env Rscript
# ============================================================================
# genome_gmm_script.R —— 网页端包装脚本（适配老师新算法 genome_gmm.R）
#
# 用法:
#   Rscript genome_gmm_script.R <data_file> <species_name> <n> <output_dir>
#
# 成功（退出码 0，stdout 含算法自身打印的 "Algorithm execution completed!"）:
#   <output_dir>/<样本名>.histo_hap_genome_size_est.pdf   拟合图（矢量，论文用）
#   <output_dir>/<样本名>.histo_hap_genome_size_est.png   拟合图（网页显示用）
#   <output_dir>/<样本名>.histo_haploid_size.csv          结果表
#   <output_dir>/<样本名>.histo_search.csv                自动调参记录（备查）
#
# 失败（退出码 1，stdout 含 "GENOME_GMM_ERROR: <原因>"，不写出任何结果文件）
# ============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4L) {
  cat("Usage: Rscript genome_gmm_script.R <data_file> <species_name> <n> <output_dir>\n")
  quit(save = "no", status = 2)
}

data_file    <- args[1]
species_name <- args[2]
n            <- suppressWarnings(as.integer(args[3]))
output_dir   <- args[4]

fail <- function(msg) {
  cat("GENOME_GMM_ERROR:", msg, "\n")
  quit(save = "no", status = 1)
}

if (is.na(n) || n < 1L)      fail("n (ploidy) must be a positive integer.")
if (!file.exists(data_file)) fail(paste("data file not found:", data_file))
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

# 定位本脚本所在目录，找到老师的算法文件 genome_gmm.R（两者须放在同一目录）
cmd_args    <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", cmd_args[grepl("^--file=", cmd_args)][1])
script_dir  <- dirname(normalizePath(script_path))
gmm_source  <- file.path(script_dir, "genome_gmm.R")
if (!file.exists(gmm_source)) {
  fail(paste("genome_gmm.R not found next to this script:", gmm_source))
}

data_file_abs  <- normalizePath(data_file)
output_dir_abs <- normalizePath(output_dir)
sample_base    <- tools::file_path_sans_ext(basename(data_file_abs))

# 算法把 PDF 和调参表写到“当前工作目录”，所以先切到输出目录再运行
setwd(output_dir_abs)
source(gmm_source)

fit <- NULL
err_msg <- NULL
set.seed(1)
tryCatch(
  fit <- GenomeGMM(data_file = data_file_abs, n = n, species_name = species_name),
  error = function(e) err_msg <<- conditionMessage(e)
)
if (is.null(fit)) {
  fail(if (is.null(err_msg)) "unknown error." else err_msg)
}
if (!is.finite(fit$genome_size_bp) || fit$genome_size_bp <= 0) {
  fail("fit finished but did not produce a valid genome size estimate.")
}

# ---- 整理输出文件名（去掉时间戳，保持网页端文件契约稳定）----
pdf_out <- paste0(sample_base, ".histo_hap_genome_size_est.pdf")
if (!is.null(fit$pdf_file) && file.exists(fit$pdf_file)) {
  if (file.exists(pdf_out)) file.remove(pdf_out)
  file.rename(fit$pdf_file, pdf_out)
}
if (!is.null(fit$search_file) && file.exists(fit$search_file)) {
  search_out <- paste0(sample_base, ".histo_search.csv")
  if (file.exists(search_out)) file.remove(search_out)
  file.rename(fit$search_file, search_out)
}

# ---- PDF -> PNG（网页内嵌显示用；依次尝试 pdftoppm / ImageMagick / R 包 pdftools）----
png_out <- paste0(sample_base, ".histo_hap_genome_size_est.png")
if (file.exists(pdf_out)) {
  if (file.exists(png_out)) file.remove(png_out)
  if (nzchar(Sys.which("pdftoppm"))) {
    system2("pdftoppm", c("-png", "-r", "150", "-singlefile",
                          shQuote(pdf_out), shQuote(sub("\\.png$", "", png_out))))
  } else if (nzchar(Sys.which("magick"))) {
    system2("magick", c("-density", "150", shQuote(paste0(pdf_out, "[0]")),
                        shQuote(png_out)))
  } else if (nzchar(Sys.which("convert"))) {
    system2("convert", c("-density", "150", shQuote(paste0(pdf_out, "[0]")),
                         shQuote(png_out)))
  } else if (requireNamespace("pdftools", quietly = TRUE)) {
    try(pdftools::pdf_convert(pdf_out, format = "png", pages = 1, dpi = 150,
                              filenames = png_out), silent = TRUE)
  }
  if (!file.exists(png_out)) {
    cat("GENOME_GMM_WARN: PNG conversion failed; install poppler-utils (pdftoppm)",
        "to enable in-page plot display. The PDF result is still available.\n")
  }
}

# ---- 结果表（网页表格与下载用）----
result <- data.frame(
  ploidy          = n,
  genome_size_bp  = round(fit$genome_size_bp),
  haploid_size_bp = round(fit$genome_size_bp / n),
  kmer_depth      = round(min(fit$miu), 2)
)
write.csv(result, paste0(sample_base, ".histo_haploid_size.csv"),
          row.names = FALSE)

quit(save = "no", status = 0)
