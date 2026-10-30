#!/usr/bin/env Rscript
# ============================================================================
# genome_gmm_script.R —— 网页端包装脚本（适配老师新算法 genome_gmm.R）
#
# 用法:
#   Rscript genome_gmm_script.R <data_file> <species_name> <n|auto> <output_dir>
#                               [ploidy_type] [mono_low_mb mono_high_mb]
#
#   第 3 个参数传 "auto"   表示不指定倍性，由算法自动推断并生成 TXT 报告
#   第 5 个参数 ploidy_type 仅在新版 genome_gmm.R + 指定 n 时需要：
#                          "autopolyploid" / "allopolyploid"
#   第 5/6（或 6/7）个参数为可选 Monoploid genome size 范围（Mb），仅自动模式可用
#
#   参数顺序容错：是否传了 ploidy_type 由"第 5 个参数是否为合法字符串"判定，
#   纯数字则视为 monoploid 范围，从而兼容旧版调用方式（老 server.R）。
#
# 成功（退出码 0，stdout 由本脚本打印 "Algorithm execution completed!"，
# 因为新版 genome_gmm.R 的 GenomeGMM() 用 capture.output 吞掉了算法内部的打印）:
#   <output_dir>/<样本名>.histo_hap_genome_size_est.pdf   拟合图（矢量，论文用）
#   <output_dir>/<样本名>.histo_hap_genome_size_est.png   拟合图（网页显示用）
#   <output_dir>/<样本名>.histo_haploid_size.csv          结果表
#   <output_dir>/<样本名>.histo_search.csv                自动调参记录（备查，若有）
#   <output_dir>/<样本名>.histo_result.txt                文字报告（自动/异源模式）
#
# 失败（退出码 1，stdout 含 "GENOME_GMM_ERROR: <原因>"，不写出任何结果 CSV）
#   新版算法的两类“不抛错但无结果”路径（无可靠基因组峰、异源倍性无法分解为
#   平衡基础模式）也会被识别，并把算法给出的说明原样上报给网页端。
# ============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4L) {
  cat("Usage: Rscript genome_gmm_script.R <data_file> <species_name> <n|auto> <output_dir>",
      "[ploidy_type] [mono_low_mb mono_high_mb]\n")
  quit(save = "no", status = 2)
}

data_file    <- args[1]
species_name <- args[2]
# 第 3 个参数传 "auto" 时，以 n = NULL 调用 GenomeGMM：算法自动推断倍性，
# 并额外生成 <样本名>.histo_result.txt 推断报告
n <- if (identical(tolower(trimws(args[3])), "auto")) {
  NULL
} else {
  suppressWarnings(as.integer(args[3]))
}
output_dir   <- args[4]

# ---- 第 5 个参数：ploidy_type（新版）或 monoploid 下限（旧调用方式）----------
valid_ploidy_types <- c("autopolyploid", "allopolyploid")
ploidy_type <- NULL
extra_args  <- character(0)
if (length(args) >= 5L) {
  candidate <- tolower(trimws(args[5L]))
  if (candidate %in% valid_ploidy_types) {
    ploidy_type <- candidate
    extra_args  <- args[-seq_len(5L)]
  } else {
    extra_args  <- args[-seq_len(4L)]
  }
}

# ---- 可选 Monoploid genome size 范围（Mb，仅自动模式）------------------------
numeric_extra <- suppressWarnings(as.numeric(extra_args))
numeric_extra <- numeric_extra[is.finite(numeric_extra)]
mono_low  <- if (length(numeric_extra) >= 1L) numeric_extra[1L] else NA_real_
mono_high <- if (length(numeric_extra) >= 2L) numeric_extra[2L] else NA_real_

fail <- function(msg) {
  cat("GENOME_GMM_ERROR:", msg, "\n")
  quit(save = "no", status = 1)
}

if (!is.null(n) && (is.na(n) || n < 2L)) {
  fail("n (ploidy) must be an integer >= 2, or 'auto' for automatic inference.")
}
if (!is.null(ploidy_type) && is.null(n)) {
  fail("ploidy type can only be used together with a supplied n.")
}

monoploid_range <- NULL
if (!is.na(mono_low) || !is.na(mono_high)) {
  if (!is.null(n)) {
    fail("monoploid size range can only be used in auto mode (without n).")
  }
  if (is.na(mono_low) || is.na(mono_high) || mono_low <= 0 || mono_high <= 0 ||
      mono_low >= mono_high) {
    fail("monoploid size range must be two positive values with lower < upper (Mb).")
  }
  monoploid_range <- c(mono_low, mono_high)
}
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

# ---- 后端能力探测：GenomeGMM() 是否支持 ploidy_type 形参 ---------------------
# 旧版 genome_gmm.R 没有该形参，传了会报 unused argument；新版在“指定 n”时强制要求。
#
# 注意：必须 source 到全局环境（与老师原本的用法一致）。genome_gmm.R 内部的
# 辅助函数（如 .GenomeGMM_integrate_red_then_observed）是按“同一文件同处一个
# 全局命名空间”的前提互相调用的；若 source 进私有 environment，闭包会找不到
# 这些同文件函数，报“没有XXX这个函数”，进而让所有候选拟合失败。
source(gmm_source)
backend_supports_ploidy_type <- "ploidy_type" %in% names(formals(GenomeGMM))
if (!is.null(ploidy_type) && !backend_supports_ploidy_type) {
  fail(paste("the installed genome_gmm.R does not support ploidy_type;",
             "update the backend or choose the legacy behaviour."))
}
if (!is.null(n) && backend_supports_ploidy_type && is.null(ploidy_type)) {
  fail(paste("this genome_gmm.R requires a ploidy type (autopolyploid / allopolyploid)",
             "when n is supplied."))
}

# 算法把 PDF 和调参表写到“当前工作目录”，所以先切到输出目录再运行
setwd(output_dir_abs)

fit <- NULL
err_msg <- NULL
set.seed(1)
tryCatch(
  {
    if (backend_supports_ploidy_type) {
      fit <- GenomeGMM(
        data_file = data_file_abs, n = n, species_name = species_name,
        monoploid_range = monoploid_range, ploidy_type = ploidy_type
      )
    } else {
      fit <- GenomeGMM(
        data_file = data_file_abs, n = n, species_name = species_name,
        monoploid_range = monoploid_range
      )
    }
  },
  error = function(e) err_msg <<- conditionMessage(e)
)
if (is.null(fit)) {
  fail(if (is.null(err_msg)) "unknown error." else err_msg)
}

# ---- 标量取值辅助 ------------------------------------------------------------
scalar_num <- function(value) {
  value <- suppressWarnings(as.numeric(value))
  if (length(value) == 1L && is.finite(value)) value else NA_real_
}
scalar_chr <- function(value) {
  if (is.null(value) || length(value) != 1L || is.na(value)) return(NA_character_)
  value <- as.character(value)
  if (nzchar(value)) value else NA_character_
}

# ---- 提取算法给出的文字说明（用于失败上报）-----------------------------------
report_message <- function(fit) {
  message <- scalar_chr(fit$message)
  if (!is.na(message)) return(message)
  result_file <- scalar_chr(fit$result_file)
  if (!is.na(result_file) && file.exists(result_file)) {
    lines <- tryCatch(readLines(result_file, warn = FALSE, encoding = "UTF-8"),
                      error = function(e) character(0))
    lines <- trimws(lines)
    lines <- lines[nzchar(lines)]
    # 跳过通用标题行，保留最有信息量的说明（通常为第 2 段之后）
    if (length(lines) > 3L) lines <- lines[-(1:3)]
    if (length(lines)) return(paste(lines, collapse = " "))
  }
  status <- scalar_chr(fit$status)
  if (!is.na(status)) return(paste("backend status:", status))
  "fit finished but did not produce a valid genome size estimate."
}

# ---- 无有效基因组大小的三类退出路径统一处理 ----------------------------------
genome_size_bp <- scalar_num(fit$genome_size_bp)
if (is.na(genome_size_bp) || genome_size_bp <= 0) {
  # 算法可能已写好一份 TXT 报告；收编为网页约定的文件名后保留，便于人工核对，
  # 但按“失败”处理：不产出 CSV，网页端不会显示任何结果。
  result_file <- scalar_chr(fit$result_file)
  if (!is.na(result_file) && file.exists(result_file)) {
    txt_out <- paste0(sample_base, ".histo_result.txt")
    if (file.exists(txt_out)) file.remove(txt_out)
    file.rename(result_file, txt_out)
  }
  fail(report_message(fit))
}

# ---- 整理输出文件名（去掉时间戳，保持网页端文件契约稳定）--------------------
pdf_out <- paste0(sample_base, ".histo_hap_genome_size_est.pdf")
if (!is.null(fit$pdf_file) && length(fit$pdf_file) == 1L &&
    !is.na(fit$pdf_file) && file.exists(fit$pdf_file)) {
  if (file.exists(pdf_out)) file.remove(pdf_out)
  file.rename(fit$pdf_file, pdf_out)
}
search_out <- paste0(sample_base, ".histo_search.csv")
if (!is.null(fit$search_file) && length(fit$search_file) == 1L &&
    !is.na(fit$search_file) && file.exists(fit$search_file)) {
  if (file.exists(search_out)) file.remove(search_out)
  file.rename(fit$search_file, search_out)
}
# 自动倍性模式 / 异源模式的文字报告（新版已知 n 的自交模式不生成报告文件）
txt_final <- paste0(sample_base, ".histo_result.txt")
backend_report <- scalar_chr(fit$result_file)
backend_report_kept <- FALSE
if (!is.na(backend_report) && file.exists(backend_report)) {
  if (file.exists(txt_final)) file.remove(txt_final)
  if (file.rename(backend_report, txt_final)) backend_report_kept <- TRUE
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

# ---- 倍性与基因组大小语义（新版口径）-----------------------------------------
# 新版 GenomeGMM() 传入 n 时：n 是“生物倍性”，算法内部再解析出用于拟合的
# “基础模式数”（base_pattern，异源多倍体可能小于 n）。返回字段中：
#   actual_ploidy          用户给定的生物倍性（可能为 NA：自动模式）
#   base_pattern           实际用于拟合曲线的基础模式数
#   subgenome_unit_size_bp 单个基础单元的基因组大小（异源=亚基因组单元）
#   g_hap_bp               Monoploid genome size（单套染色体组 DNA 含量）
#   genome_size_bp         Holoploid genome size（该倍性下总 DNA 含量）
# 旧版只有 g_hap_bp / genome_size_bp，且 genome_size_bp = n * g_hap_bp。
actual_ploidy <- scalar_num(fit$actual_ploidy)
base_pattern  <- scalar_num(fit$base_pattern)
fit_harmonic  <- scalar_num(fit$fit_harmonic_count)
g_hap_bp      <- scalar_num(fit$g_hap_bp)
subgenome_bp  <- scalar_num(fit$subgenome_unit_size_bp)
minimum_size  <- scalar_num(fit$minimum_pattern_genome_size_bp)
homozygous_depth <- scalar_num(fit$homozygous_peak_depth)

# 展示用倍性：已知 n 用用户给的生物倍性；自动模式用算法推断的最小模式数
fitted_ploidy <- if (!is.na(actual_ploidy) && actual_ploidy >= 1) {
  as.integer(actual_ploidy)
} else if (!is.na(fit_harmonic) && fit_harmonic >= 1) {
  as.integer(fit_harmonic)
} else if (!is.na(base_pattern) && base_pattern >= 1) {
  as.integer(base_pattern)
} else {
  n
}

# Monoploid（单套）基因组大小：优先用算法明确给出的 g_hap_bp，
# 避免旧脚本用 genome_size_bp / ploidy 导致异源多倍体口径错误。
if (is.na(g_hap_bp) || g_hap_bp <= 0) {
  g_hap_bp <- if (!is.na(genome_size_bp) && !is.na(fitted_ploidy) &&
                  fitted_ploidy >= 1) {
    genome_size_bp / fitted_ploidy
  } else {
    NA_real_
  }
}
if (is.na(homozygous_depth)) {
  miu <- suppressWarnings(as.numeric(fit$miu))
  miu <- miu[is.finite(miu) & miu > 0]
  homozygous_depth <- if (length(miu)) min(miu) else NA_real_
}

pattern_confidence <- if (!is.na(scalar_chr(fit$pattern_confidence))) {
  scalar_chr(fit$pattern_confidence)
} else if (!is.null(fit$known_ploidy_resolution$base_pattern_confidence)) {
  scalar_chr(fit$known_ploidy_resolution$base_pattern_confidence)
} else {
  NA_character_
}

# ---- 结果表（网页表格与下载用）----------------------------------------------
# 列名保持向后兼容（ploidy / genome_size_bp / haploid_size_bp / kmer_depth），
# 其后追加倍性口径相关的新列；旧列语义不变，前端表格无需改动即可显示。
result <- data.frame(
  ploidy          = fitted_ploidy,
  genome_size_bp  = round(genome_size_bp),
  haploid_size_bp = if (is.na(g_hap_bp)) NA_real_ else round(g_hap_bp),
  kmer_depth      = if (is.na(homozygous_depth)) NA_real_ else round(homozygous_depth, 2),
  stringsAsFactors = FALSE
)
if (isTRUE(backend_supports_ploidy_type)) {
  result$base_pattern <- if (is.na(base_pattern)) NA_integer_ else as.integer(base_pattern)
  # 自动模式没有 ploidy type，该列会是全 NA，索性不加，免得表格多一列空值
  if (!is.null(ploidy_type)) result$ploidy_type <- ploidy_type
  result$subgenome_size_bp <- if (is.na(subgenome_bp)) NA_real_ else round(subgenome_bp)
  if (!is.na(pattern_confidence)) result$pattern_confidence <- pattern_confidence
}
write.csv(result, paste0(sample_base, ".histo_haploid_size.csv"), row.names = FALSE)

# ---- TXT 报告 ----------------------------------------------------------------
# 自动模式 / 异源模式：算法自己生成，前面已收编为 <样本名>.histo_result.txt，保留原文。
# 已知 n 的自交模式：算法按设计不生成报告，由本脚本补写一份结果摘要。
# 注意必须每次都重写（若未保留后端报告）：同名旧 TXT 残留会让网页端 is_fresh
# 时间戳校验判定它不是本次运行的结果，下载按钮不显示。
if (!backend_report_kept) {
  writeLines(enc2utf8(c(
    "GenomeGMM ploidy and genome-size inference",
    paste("Sample:", sample_base),
    if (nzchar(species_name)) paste("Species:", species_name) else NULL,
    "",
    if (is.null(n)) {
      "Mode: automatic ploidy inference (n not supplied)"
    } else {
      "Mode: known n supplied by user (no automatic ploidy inference)"
    },
    if (!is.null(n)) {
      sprintf("Supplied biological ploidy n: %d", as.integer(n))
    } else NULL,
    if (!is.null(ploidy_type)) {
      sprintf("Supplied ploidy type: %s", ploidy_type)
    } else NULL,
    if (!is.na(base_pattern)) {
      sprintf("Base pattern used for curve fitting: %dx", as.integer(base_pattern))
    } else NULL,
    if (!is.na(pattern_confidence)) {
      sprintf("Pattern confidence: %s", pattern_confidence)
    } else NULL,
    "",
    sprintf("Holoploid genome size estimate: %.2f Mb", genome_size_bp / 1e6),
    if (!is.na(g_hap_bp)) {
      sprintf("Monoploid genome size (g-hap, 1Cx): %.2f Mb", g_hap_bp / 1e6)
    } else NULL,
    if (!is.na(subgenome_bp)) {
      sprintf("Subgenome unit size: %.2f Mb", subgenome_bp / 1e6)
    } else NULL,
    if (!is.na(minimum_size)) {
      sprintf("Minimum-pattern genome size: %.2f Mb", minimum_size / 1e6)
    } else NULL,
    if (!is.na(homozygous_depth)) {
      sprintf("Homozygous peak depth: %.2f", homozygous_depth)
    } else NULL
  )), con = txt_final, useBytes = TRUE)
}

# 网页端 server.R 靠这句判断成功；新版算法内部打印被 capture.output 吞掉，故由本脚本打印
cat("Algorithm execution completed!\n")

quit(save = "no", status = 0)
