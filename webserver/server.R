# ============================================================================
# server.R —— findGSEP Webserver（GMM 单算法版，v1.6）
#
# v1.6 改动（按老师要求重构倍性输入，取代 v1.5 的“复选框 + 常驻单选组”）：
#   1. 倍性输入改为二选一（互斥单选），不再用复选框：
#        Option 1  Estimate ploidy n  —— 不输入 n，由算法推断
#        Option 2  Provide ploidy n   —— 输入 n + 倍性状态（Autoploidy / Alloploidy）
#      选择 Provide 时，n 与倍性状态通过 conditionalPanel 才出现。
#   2. n 一律指 ploidy，因此界面上不再用长句解释“基础模式数”：
#      留一行简短提示 “n is the ploidy number (2 = diploid, 4 = tetraploid, ...)”。
#   3. 异源（Alloploidy）时前端不再自行把 n 除以 2/3/… 后分多次调用：
#      新版 genome_gmm.R 收到 n + allopolyploid 后，内部会枚举 n 的全部整数分拆
#      （等价于按 n/2、n/3… 多次尝试），用峰证据 + 曲线质量打分选出基础模式。
#      若前端先除算再传，后端会把它当成另一个“生物倍性”，结果反而错。
#   4. 移除 v1.5 的动态 n 标签（n_label）与启动时禁用 ploidy_type 的逻辑。
#   5. Monoploid 范围先验仅在 Estimate 模式下可用，切到 Provide 自动置灰。
#
# 以下为 v1.5 的历史说明：
#
# v1.5 改动（适配新版 genome_gmm.R：ploidy_type / 异源基础模式解析）：
#   1. 新增 ploidy type 控件：新版 GenomeGMM() 在“指定 n”时必须显式传入
#      autopolyploid / allopolyploid，否则直接 stop()。仍按旧契约调用的前端
#      在新版后端下会 100% 失败，这是本次适配的核心破坏性变更。
#   2. 后端能力启动探测 detect_backend_caps()：只有确认 GenomeGMM() 具备
#      ploidy_type 形参时才传该参数；旧版后端自动回退老行为，控件置灰。
#   3. n 的语义随后端版本切换：新版 n = 生物倍性（曲线基础模式由后端解析），
#      旧版 n = 最小 k-mer 基础模式数。
#   4. 示例数据同时设定 n 与倍性状态。
#   5. 结果 CSV 由包装脚本追加 base_pattern / ploidy_type / subgenome_size_bp /
#      pattern_confidence 列；旧的 ploidy / genome_size_bp / haploid_size_bp /
#      kmer_depth 列名与语义保持不变，表格与下载逻辑无需改动。
#   6. 失败提示原样传递后端说明（含“无可靠基因组峰”“异源不可分解”两类
#      不抛错但无结果的提前退出）。
#
# 以下为 v1.4 的历史说明：
#
# 本版改动：
#   1. 去掉算法二选一：网页只保留老师新算法 GenomeGMM（genome_gmm.R，
#      通过包装脚本 genome_gmm_script.R 调用），只需 species_name + n。
#      原 findGSEP 算法及其参数（exp_hom / sizek / ploidy / xlimit / ylimit）
#      相关的界面、调用、下载逻辑全部移除。
#   2. Rscript 直接写死为 GSEGMM 环境：
#      /home/shiny/miniconda3/envs/GSEGMM/bin/Rscript
#   3. 失败时不显示任何结果（图表区保持为空），只弹通知告知原因。
#      run_state 标志 + 结果文件时间戳校验（is_fresh）双重保证，
#      杜绝“失败后显示上次结果”。旧结果文件全部保留留档，不做删除。
#   4. 成功判定：退出码 0 + 输出含 "Algorithm execution completed!"
#      + 结果 CSV 为本次运行新生成。
#   5. 新增“自动倍性估计”模式（老师 2026-09 要求，对应新版 genome_gmm.R）：
#      勾选 auto_n 后无需输入 n，包装脚本以 n = NULL 调用 GenomeGMM，
#      由算法在 2–10 中自动推断倍性，并额外生成 TXT 推断报告供下载。
#   6. 适配最新版 genome_gmm.R（19466 行）：
#      a. n 的含义改为“最低 k-mer 基础模式数”（AABB 填 2，AABBCC 填 3）；
#      b. 自动模式支持可选 monoploid_range（Mb 范围），界面加两个可选输入框；
#      c. 算法新增“无可靠基因组峰”提前退出路径（只写 TXT 说明、不抛错），
#         由包装脚本识别并按失败处理、把原因展示给用户。
# ============================================================================
if (file.exists("prepare_findGSE.R")) source("prepare_findGSE.R")
library(shiny)
library(shinyjs)
library(shinyBS)
library(DT)
library(ggplot2)
library(stringr)
library(seqinr)

# GMM 算法所在的 R 环境（genome_gmm.R 只用 base 包，此环境可裸跑）
GMM_RSCRIPT <- "Rscript"   # 本机版：走 PATH；若本机 R 不在 PATH，改成你电脑上 Rscript 的完整路径

# ---- 后端能力探测（向后兼容）----------------------------------------------
# 新版 genome_gmm.R 的 GenomeGMM() 增加了 ploidy_type 形参，并在“指定 n”时必须显式
# 给出 autopolyploid / allopolyploid；同时新增“无可靠峰”“异源不可分解”等不抛错的
# 提前退出路径。旧版没有这些概念。这里在启动时读一次函数签名：
#   支持新契约：只有确定支持时才让包装脚本传 ploidy_type，否则回退旧行为。
# 探测失败一律按“不支持”处理（安全侧），避免误把 ploidy_type 传给旧版而报 unused argument。
#
# 注意：探测放在**独立 Rscript 子进程**里执行，不把 873 KB 的后端 source 进 Shiny
# 的全局环境（后端内部函数依赖同一全局命名空间互相调用，若在本进程 source，
# 既会污染 app 变量，也可能与其内部辅助函数解析产生干扰）。
detect_backend_caps <- function() {
  caps <- list(ploidy_type = FALSE, version = "unknown")
  if (!nzchar(Sys.which(GMM_RSCRIPT)) ||
      !file.exists("genome_gmm.R")) return(caps)
  expr <- paste(
    "source('genome_gmm.R');",
    "cat(if ('ploidy_type' %in% names(formals(GenomeGMM))) 'NEW' else 'OLD')"
  )
  out <- tryCatch(
    suppressWarnings(system2(GMM_RSCRIPT, c("-e", shQuote(expr)),
      stdout = TRUE, stderr = FALSE, timeout = 300)),
    error = function(e) character(0)
  )
  out <- trimws(paste(out, collapse = ""))
  if (identical(out, "NEW")) {
    caps$ploidy_type <- TRUE
    caps$version <- "new (ploidy_type supported)"
  } else if (identical(out, "OLD")) {
    caps$version <- "legacy (no ploidy_type)"
  }
  caps
}

# 新版不传 n 时不允许同时给 ploidy_type；这里只为“指定 n”模式做校验
ploidy_type_problem <- function(ploidy_type, caps, ploidy_requested) {
  if (!isTRUE(ploidy_requested)) return(NULL)
  if (!isTRUE(caps$ploidy_type)) return(NULL)
  if (is.null(ploidy_type) || length(ploidy_type) != 1L || is.na(ploidy_type) ||
      !ploidy_type %in% c("autopolyploid", "allopolyploid")) {
    return("Please choose a ploidy type (autopolyploid or allopolyploid) when you enter n.")
  }
  NULL
}

# NULL 合并运算符（shiny 旧版不保证存在）
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a

# Define server logic required to draw a histogram
shinyServer(function(input, output,session) {
  shinyjs::disable("ResultDown")
  shinyjs::hide("seq_name")
  shinyjs::hide("dVisualize")
  shinyjs::hide("dTxt")   # TXT 报告下载：成功运行后出现
  # Monoploid size 范围仅 Estimate 模式可用；默认就是 Estimate，故启动时保持可用，
  # 切到 Provide 时由 observeEvent 关闭
  shinyjs::enable("mono_low")
  shinyjs::enable("mono_high")
  shinyjs::hide("seq_table")
  shinyjs::disable("draw_plot")
  shinyjs::hideElement(id= "panel_forna")

  # 本次会话是否已有“可用的新结果”（失败/未运行时一律 FALSE，结果区不显示任何内容）
  run_state <- reactiveVal(FALSE)

  # 后端能力：决定是否传 ploidy_type（旧版后端没有该形参，传了会报 unused argument）
  backend_caps <- detect_backend_caps()

  # 本次提交时间戳：判断结果文件是否为“本次运行新生成”。
  # 旧结果文件全部保留在服务器上留档，不做任何删除
  submit_time <- reactiveVal(NULL)
  is_fresh <- function(f) {
    !is.null(submit_time()) && file.exists(f) && file.mtime(f) > submit_time()
  }

  # 结果文件路径
  out_dir <- "input_histo/outfile"
  result_paths <- reactive({
    base <- file.path(out_dir, eff_outfile())
    list(
      png = paste0(base, ".histo_hap_genome_size_est.png"),
      pdf = paste0(base, ".histo_hap_genome_size_est.pdf"),
      csv = paste0(base, ".histo_haploid_size.csv"),
      txt = paste0(base, ".histo_result.txt")   # 自动倍性估计的文字报告（仅自动模式生成）
    )
  })

  # ---- 输出文件名自动生成：数据文件名（去扩展名）+ "_ploidy<n>" ----
  # 例如 this_is_filename.histo、n=4 → this_is_filename_ploidy4
  # 示例数据同样按示例原文件名（物种名.histo）生成；输出名留空时运行前自动补默认名
  auto_name <- reactiveVal("result")        # 当前自动生成的默认值
  outfile_edited <- reactiveVal(FALSE)      # 用户是否手动改过输出名（改过就不再覆盖）
  src_file <- reactiveVal(NULL)             # 当前数据来源文件名（上传文件或示例数据）
  eff_outfile <- reactiveVal("result")      # 本次运行实际生效的输出名（提交时确定）

  make_default_name <- function(fname, n_val) {
    base <- sub("\\.[^.]*$", "", fname)     # 去掉最后一个扩展名（.histo/.txt 等）
    paste0(base, "_ploidy", n_val)
  }

  # 按当前数据来源 + 倍性输入方式生成默认名；估计模式用 "_ploidyAuto"
  apply_default_name <- function() {
    sf <- src_file()
    if (is.null(sf)) return(invisible())
    nm <- if (identical(input$ploidy_mode, "estimate")) {
      paste0(sub("\\.[^.]*$", "", sf), "_ploidyAuto")
    } else {
      if (is.na(input$n)) return(invisible())
      make_default_name(sf, input$n)
    }
    auto_name(nm)
    updateTextInput(session, "outfile", value = nm)
  }

  # 倍性输入方式二选一：Estimate（估计） vs Provide（提供）
  #   Estimate → 不输入 n，Monoploid 范围先验可用
  #   Provide  → 输入 n（ploidy）+ 倍性状态；Monoploid 范围不适用
  observeEvent(input$ploidy_mode, {
    if (identical(input$ploidy_mode, "estimate")) {
      shinyjs::enable("mono_low")
      shinyjs::enable("mono_high")
    } else {
      shinyjs::disable("mono_low")
      shinyjs::disable("mono_high")
    }
    if (!is.null(src_file()) && !outfile_edited()) apply_default_name()
  })

  # 后端是旧版（无 ploidy_type 形参）：倍性状态控件不参与计算，仅提示
  observe({
    if (!isTRUE(backend_caps$ploidy_type)) {
      shinyjs::disable("ploidy_type")
    }
  })

  # ---- species name 同理：默认取数据文件名（去扩展名），留空时运行前自动补 ----
  auto_species <- reactiveVal("")          # 当前自动生成的默认物种名
  species_edited <- reactiveVal(FALSE)     # 用户是否手动改过物种名

  make_default_species <- function(fname) sub("\\.[^.]*$", "", fname)

  apply_default_species <- function() {
    sf <- src_file()
    if (is.null(sf)) return(invisible())
    sp <- make_default_species(sf)
    auto_species(sp)
    updateTextInput(session, "species_name", value = sp)
  }

  # 检测用户手动编辑（自动更新写入的值等于 auto_species，不算手动）
  observeEvent(input$species_name, {
    if (input$species_name != auto_species()) species_edited(TRUE)
  }, ignoreInit = TRUE)

  # 上传新文件：总是重置为自动生成（换数据 = 新的默认值）
  observeEvent(input$file_input, {
    req(input$file_input)
    src_file(input$file_input$name)
    outfile_edited(FALSE)
    species_edited(FALSE)
    apply_default_name()
    apply_default_species()
  })

  # n 变化：已有数据来源且用户没手动改过名字时，联动更新
  observeEvent(input$n, {
    if (!is.null(src_file()) && !outfile_edited()) apply_default_name()
  })

  # 检测用户手动编辑（自动更新写入的值等于 auto_name，不算手动）
  observeEvent(input$outfile, {
    if (input$outfile != auto_name()) outfile_edited(TRUE)
  }, ignoreInit = TRUE)

  observeEvent(input$example,{
    selected_species <- input$example_sel
    sel_file_name <- switch(selected_species,
      "javanica_tetraploid" = "demo_data/javanica.histo",
      "potato_tetraploid" = "demo_data/O_21mer_new20231012.histo",
      "wheat_hexaploid" = "demo_data/wheat.histo",
      "strawberry_octoploid" = "demo_data/strawberry.histo",
      "Chinese_sturgeon_octoploid" = "demo_data/zhonghuaxun_newgenerate.histo"
    )
    shinyjs::disable("example")
    shinyjs::disable("your_text")
    data_example <- read.table(sel_file_name,header =FALSE)
    text_example <- paste(data_example$V1, data_example$V2, sep = " ",collapse="\n")
    updateTextAreaInput(session, "your_text", value = text_example)
    shinyjs::enable("example")
    shinyjs::enable("your_text")
    # 示例数据：按示例原文件名（物种名.histo）+ 当前 n 生成默认输出名；
    # species name 同理，默认填示例物种名
    src_file(basename(sel_file_name))
    outfile_edited(FALSE)
    species_edited(FALSE)
    apply_default_name()
    apply_default_species()
  })

  observeEvent(input$submit, {
    # 校验提示只在点击 Submit 的瞬间出现一次：
    # 用 isolate 断开响应式依赖，避免跑完一组回到 Home 修改输入时灰字自动冒出
    output$form_table <- renderTable({
      sp <- isolate(input$species_name)
      n_val <- isolate(input$n)
      estimate_mode <- identical(isolate(input$ploidy_mode), "estimate")
      m_lo <- isolate(input$mono_low)
      m_hi <- isolate(input$mono_high)
      p_type <- isolate(input$ploidy_type)
      range_ok <- !estimate_mode || ((is.na(m_lo) && is.na(m_hi)) ||
        (!is.na(m_lo) && !is.na(m_hi) && m_lo > 0 && m_hi > m_lo))
      has_data <- isolate(!is.null(input$file_input) || str_length(input$your_text) > 0)
      type_problem <- ploidy_type_problem(p_type, backend_caps, !estimate_mode)
      validate(
        need(estimate_mode || (!is.na(n_val) && n_val >= 2), "Ploidy n should be an integer >= 2 (or choose Estimate ploidy n)!!!"),
        need(is.null(type_problem), type_problem),
        need(range_ok, "Monoploid size range: fill both bounds (Mb), lower < upper; or leave both empty!!!"),
        need(has_data, "You should specify either file or text")
      )
      data.frame()

    })
  })


  dt <- reactive({
    if(is.null(input$file_input)){
      return()
    }else{
      read.table(input$file_input$datapath,header =FALSE)
    }
  })

  observeEvent(input$submit, {
    estimate_mode <- identical(input$ploidy_mode, "estimate")
    range_ok <- !estimate_mode ||
      ((is.na(input$mono_low) && is.na(input$mono_high)) ||
       (!is.na(input$mono_low) && !is.na(input$mono_high) &&
        input$mono_low > 0 && input$mono_high > input$mono_low))
    # 新版后端要求：提供倍性时必须显式给出倍性状态，否则 GenomeGMM() 直接报错
    p_type <- input$ploidy_type
    type_problem <- ploidy_type_problem(p_type, backend_caps, !estimate_mode)
    validate(
      need(estimate_mode || (!is.na(input$n) && input$n >= 2), "Ploidy n should be an integer >= 2 (or choose Estimate ploidy n)!!!"),
      need(is.null(type_problem), type_problem),
      need(range_ok, "Monoploid size range: fill both bounds (Mb), lower < upper; or leave both empty!!!"),
      need(!is.null(input$file_input) || str_length(input$your_text) > 0, "You should specify either file or text")
    )

    # species name 留空 → 按规则补默认名（数据文件名去扩展名；无数据来源用 "unknown_species"）
    sp_val <- input$species_name
    if (str_length(sp_val) == 0) {
      sp_val <- if (!is.null(src_file())) make_default_species(src_file()) else "unknown_species"
      species_edited(FALSE)
      auto_species(sp_val)
      updateTextInput(session, "species_name", value = sp_val)
    }

    # 输出名留空 → 自动补默认名（文件名_ploidy<n>；无数据来源时用 "result"）
    out_val <- input$outfile
    if (str_length(out_val) == 0) {
      out_val <- auto_name()
      outfile_edited(FALSE)
      updateTextInput(session, "outfile", value = out_val)
    }
    eff_outfile(out_val)

    showNotification("Start processing, please wait...",duration = 3)
    shinyjs::disable("submit")
    shinyjs::disable("draw_plot")
    shinyjs::disable("ResultDown")

    shinyjs::hide("myImage")
    shinyjs::hide("seq_name")
    shinyjs::hide("dVisualize")
    shinyjs::hide("dTxt")
    shinyjs::hide("table")

    # ---- 第 1 重保证：提交即失效旧结果（界面层）----
    run_state(FALSE)
    # 记录本次提交时间：第 2 重保证用它判断结果文件是不是“这次跑出来的”
    submit_time(Sys.time())

    if(is.null(input$file_input)){
      write.table(input$your_text, paste0('input_histo/',out_val,'.histo'), sep = ' ', col.names = F, row.names = F, quote = F)
    }else{
      write.table(dt(), paste0('input_histo/',out_val,'.histo'), sep = ' ', col.names = F, row.names = F, quote = F)
    }

    full_file_path <- paste0("input_histo/", out_val, ".histo")
    withProgress(message = 'Running findGSEP', value = 0, {
      incProgress(1/2, detail = c("Please wait..."))
      Sys.sleep(0.1)

      # ---- 老师新算法 genome_gmm.R（包装脚本负责 PDF 转 PNG、写结果 CSV）----
      # Estimate 模式：传 "auto"，包装脚本以 n = NULL 调用 GenomeGMM，
      # 由算法自动推断倍性并额外生成 TXT 推断报告
      # Provide 模式：传 n（ploidy）+ 倍性状态。异源时由后端内部枚举 n 的整数分拆
      # （等价于老师要求的“按 n/2、n/3… 多次估计、按峰数挑最优”），前端不再自行除算，
      # 否则后段会把除后的值误当作另一个生物倍性。
      n_arg <- if (estimate_mode) "auto" else as.character(input$n)
      # 可选 Monoploid genome size 范围（Mb，仅 Estimate 模式；两个都填才传）
      range_args <- character(0)
      if (estimate_mode && !is.na(input$mono_low) && !is.na(input$mono_high)) {
        range_args <- c(as.character(input$mono_low), as.character(input$mono_high))
      }
      # 倍性状态：仅在新版后端且 Provide 模式时传递；Estimate/旧版后端一律不传，
      # 保持与 GenomeGMM() 的形参契约一致（旧版没有该形参，传了会报 unused argument）
      type_arg <- character(0)
      if (isTRUE(backend_caps$ploidy_type) && !estimate_mode) {
        type_arg <- as.character(p_type)
      }
      output_GSE <- system2(GMM_RSCRIPT,
                            c("genome_gmm_script.R",
                              shQuote(full_file_path),
                              shQuote(sp_val),
                              n_arg,
                              shQuote(out_dir),
                              type_arg,
                              range_args
                            ),
                            stdout = TRUE, stderr = TRUE)
      incProgress(1/3, detail = c("Please wait..."))
      cat(output_GSE)
      incProgress(1/2, detail = c("Task finished."))
      Sys.sleep(1)
    })

    # 退出码：system2 捕获输出时，非零退出码会记录在 "status" 属性里
    run_status <- attr(output_GSE, "status")
    if (is.null(run_status)) run_status <- 0L
    output_GSE <- paste(output_GSE, collapse = "\n")

    # ---- 成功判定：成功标记 + 结果文件必须是本次运行新生成的 ----
    # is_fresh()：文件存在且修改时间晚于本次提交时间 —— 同名重跑失败时，
    # 磁盘上的旧文件不会被误判为本次结果（旧文件保留，不删除）
    paths <- result_paths()
    success <- identical(as.integer(run_status), 0L) &&
      grepl("Algorithm execution completed", output_GSE) &&
      is_fresh(paths$csv)

    if (!success) {
      # ================= 失败分支：不显示任何结果，只提示 =================
      run_state(FALSE)
      # 包装脚本会把后端原因统一成一行 "GENOME_GMM_ERROR: <原因>" 输出。
      # 新版后端新增两类“不抛错但无结果”的提前退出（无可靠基因组峰、异源倍性
      # 无法分解为平衡基础模式），其说明文字也经此通道原样送达用户。
      err_line <- regmatches(output_GSE,
                             regexpr("GENOME_GMM_ERROR:[^\n]*", output_GSE))
      if (length(err_line) > 0 && nchar(err_line) > 0) {
        err_msg <- trimws(gsub("GENOME_GMM_ERROR:\\s*", "", err_line))
      } else {
        r_err <- grep("Error", strsplit(output_GSE, "\n", fixed = TRUE)[[1]],
                      value = TRUE)
        err_msg <- if (length(r_err) > 0) tail(r_err, 1) else "please check your data"
      }
      if (!nzchar(err_msg)) err_msg <- "please check your data"
      showNotification(
        paste0("Run failed: ", err_msg,
               " No result is shown. Please check your data and parameters, then try again."),
        type = "error", duration = 30
      )
      session$sendCustomMessage(type = 'testmessage',
        message = 'The task failed. Please check the error message and try again.')
      shinyjs::enable("submit")
      return()
    }

    # ================= 成功分支 =================
    run_state(TRUE)

    # 自动倍性估计的 TXT 报告：本次运行新生成才显示下载按钮
    if (is_fresh(paths$txt)) {
      shinyjs::show("dTxt")
    }

    # 若 PNG 转换失败（服务器缺 poppler），结果表照常，仅提示
    if (!is_fresh(paths$png) && is_fresh(paths$pdf)) {
      showNotification(
        "Result table is ready, but the in-page plot is unavailable on this server (PDF-to-PNG converter missing). You can still download the PDF figure.",
        type = "warning", duration = 15
      )
    }

    t_count <- try(system("echo $((`cat user_num`+1)) > user_num"))
    session$sendCustomMessage(type = 'testmessage',
                              message = paste0('The task is done!',' Thank you for waiting.'))
    shinyjs::enable("submit")
    shinyjs::enable("draw_plot")
    updateTabsetPanel(session,"tabs", selected = "Download")
    showNotification(paste("Congrates, you are the current user number: ",system("echo $((`cat user_num`))",intern=TRUE)),duration = 3)
  })

  observeEvent(input$example_sel, {
      selected_species <- input$example_sel
      # 示例数据同时给 n 和倍性状态：新版后端在“提供倍性”时强制要求 ploidy type。
      # 各示例的分类依据（均为公开文献，非推测）：
      #   wheat    Triticum aestivum  异源六倍体 AABBDD（教科书结论）
      #   potato   Solanum tuberosum  同源四倍体（减数分裂形成四价体）
      #   strawberry Fragaria × ananassa 异源八倍体
      #   javanica Meloidogyne javanica 4n 异源，基因组结构 A1A2B1B2
      #            Nat Commun 2023, doi:10.1038/s41467-023-42700-w
      #   Chinese_sturgeon Acipenser sinensis 同源八倍体（autooctoploid）
      #            Wang et al. 2024, Genom Proteom Bioinform, doi:10.1093/gpbjnl/qzad002
      example_setup <- switch(selected_species,
        "javanica_tetraploid"          = list(n = 4, type = "allopolyploid"),
        "potato_tetraploid"            = list(n = 4, type = "autopolyploid"),
        "wheat_hexaploid"              = list(n = 6, type = "allopolyploid"),
        "strawberry_octoploid"         = list(n = 8, type = "allopolyploid"),
        "Chinese_sturgeon_octoploid"   = list(n = 8, type = "autopolyploid")
      )
      if (!is.null(example_setup)) {
        updateNumericInput(session, "n", value = example_setup$n)
        if (isTRUE(backend_caps$ploidy_type)) {
          updateRadioButtons(session, "ploidy_type", selected = example_setup$type)
        }
        updateActionButton(session, "example", label = paste("Use", selected_species, "data"))
      }

  })

  observeEvent(input$draw_plot,{
    # 失败或未运行时：不给图，只提示
    if (!isTRUE(run_state())) {
      showNotification(
        "No result to display: the last run failed or has not finished. Please re-run successfully first.",
        type = "warning", duration = 8
      )
      return()
    }
    shinyjs::enable("ResultDown")
    shinyjs::show("dVisualize")
    shinyjs::show("table")

    png_file <- result_paths()$png
    if(is_fresh(png_file)){
      shinyjs::show("myImage")
      output$myImage <-renderImage({
        list(src = png_file, width = "100%", align = "center")
      })
    }
  })

  out_haploid_size <- reactive({
    # 只有本次运行成功才读结果文件；失败/未运行一律返回空，杜绝残留旧结果
    if (!input$draw_plot || !isTRUE(run_state())) {
      return(NULL)
    }
    paths <- result_paths()
    if (file.exists(paths$csv)) {
      read.csv(paths$csv, header = TRUE, sep = ',', row.names = NULL)
    } else {
      return(NULL)
    }
  })

  output$table <- DT::renderDataTable(DT::datatable(
	    out_haploid_size(),
	    options = list(
                          searching = FALSE,  # Disable searching
                          paging = FALSE       # Disable paging
                      ),
	    caption = tags$caption(
	        style = 'font-weight: bold; font-size: 18px;',
	        'Table 1: Size estimation results.'
	    ),
	    escape = FALSE
	))

  # 图下载：优先 PNG；若服务器缺 PDF 转换工具则退化为 PDF
  # 注：磁盘上的结果文件带 .histo 中缀（算法侧命名，不动）；下载名清理掉，更干净
  output$dVisualize <- downloadHandler(
    filename = function() {
      if (file.exists(result_paths()$png)) {
        paste0(eff_outfile(), ".png")
      } else {
        paste0(eff_outfile(), ".pdf")
      }
    },
    content = function(file) {
      png_file <- result_paths()$png
      pdf_file <- result_paths()$pdf
      if (file.exists(png_file)) {
        file.copy(png_file, file)
      } else if (file.exists(pdf_file)) {
        file.copy(pdf_file, file)
      }
    },
    contentType = "application/octet-stream"
  )

  # 自动倍性估计的 TXT 推断报告下载（仅自动模式产生该文件）
  output$dTxt <- downloadHandler(
    filename = function() {
      paste0(eff_outfile(), "_result.txt")
    },
    content = function(file) {
      txt_file <- result_paths()$txt
      if (file.exists(txt_file)) {
        file.copy(txt_file, file)
      }
    },
    contentType = "text/plain"
  )

  # 结果下载：单个 CSV
  output$ResultDown <- downloadHandler(
    filename = function() {
      paste0(eff_outfile(), "_genome_size_est.csv")
    },
    content = function(file) {
      csv_file <- result_paths()$csv
      if (file.exists(csv_file)) {
        file.copy(csv_file, file)
      }
    },
    contentType = "text/csv"
  )

})
