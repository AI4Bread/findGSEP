source("prepare_findGSE.R")
library(shiny)
library(shinyjs)
library(shinyBS)
library(DT)
library(ggplot2)
library(stringr)
library(seqinr)

GMM_RSCRIPT <- "/home/shiny/miniconda3/envs/GSEGMM/bin/Rscript"

shinyServer(function(input, output,session) {
  shinyjs::disable("ResultDown")
  shinyjs::hide("seq_name")
  shinyjs::hide("dVisualize")
  shinyjs::hide("seq_table")
  shinyjs::disable("draw_plot")
  shinyjs::hideElement(id= "panel_forna")

  run_state <- reactiveVal(FALSE)

  submit_time <- reactiveVal(NULL)
  is_fresh <- function(f) {
    !is.null(submit_time()) && file.exists(f) && file.mtime(f) > submit_time()
  }

  out_dir <- "input_histo/outfile"
  result_paths <- reactive({
    base <- file.path(out_dir, eff_outfile())
    list(
      png = paste0(base, ".histo_hap_genome_size_est.png"),
      pdf = paste0(base, ".histo_hap_genome_size_est.pdf"),
      csv = paste0(base, ".histo_haploid_size.csv")
    )
  })

  auto_name <- reactiveVal("result")        
  outfile_edited <- reactiveVal(FALSE)      
  src_file <- reactiveVal(NULL)         
  eff_outfile <- reactiveVal("result")    

  make_default_name <- function(fname, n_val) {
    base <- sub("\\.[^.]*$", "", fname)     
    paste0(base, "_ploidy", n_val)
  }

  apply_default_name <- function() {
    sf <- src_file()
    if (is.null(sf) || is.na(input$n)) return(invisible())
    nm <- make_default_name(sf, input$n)
    auto_name(nm)
    updateTextInput(session, "outfile", value = nm)
  }

  auto_species <- reactiveVal("")          
  species_edited <- reactiveVal(FALSE)     

  make_default_species <- function(fname) sub("\\.[^.]*$", "", fname)

  apply_default_species <- function() {
    sf <- src_file()
    if (is.null(sf)) return(invisible())
    sp <- make_default_species(sf)
    auto_species(sp)
    updateTextInput(session, "species_name", value = sp)
  }

  observeEvent(input$species_name, {
    if (input$species_name != auto_species()) species_edited(TRUE)
  }, ignoreInit = TRUE)

  observeEvent(input$file_input, {
    req(input$file_input)
    src_file(input$file_input$name)
    outfile_edited(FALSE)
    species_edited(FALSE)
    apply_default_name()
    apply_default_species()
  })

  observeEvent(input$n, {
    if (!is.null(src_file()) && !outfile_edited()) apply_default_name()
  })

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
    src_file(basename(sel_file_name))
    outfile_edited(FALSE)
    species_edited(FALSE)
    apply_default_name()
    apply_default_species()
  })

  observeEvent(input$submit, {
    output$form_table <- renderTable({
      sp <- isolate(input$species_name)
      n_val <- isolate(input$n)
      has_data <- isolate(!is.null(input$file_input) || str_length(input$your_text) > 0)
      validate(
        need(!is.na(n_val) && n_val >= 1, "n should be an integer >= 1!!!"),
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
    validate(
      need(!is.na(input$n) && input$n >= 1, "n should be an integer >= 1!!!"),
      need(!is.null(input$file_input) || str_length(input$your_text) > 0, "You should specify either file or text")
    )

    sp_val <- input$species_name
    if (str_length(sp_val) == 0) {
      sp_val <- if (!is.null(src_file())) make_default_species(src_file()) else "unknown_species"
      species_edited(FALSE)
      auto_species(sp_val)
      updateTextInput(session, "species_name", value = sp_val)
    }

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
    shinyjs::hide("table")

    run_state(FALSE)
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

      output_GSE <- system2(GMM_RSCRIPT,
                            c("genome_gmm_script.R",
                              shQuote(full_file_path),
                              shQuote(sp_val),
                              input$n,
                              shQuote(out_dir)
                            ),
                            stdout = TRUE, stderr = TRUE)
      incProgress(1/3, detail = c("Please wait..."))
      cat(output_GSE)
      incProgress(1/2, detail = c("Task finished."))
      Sys.sleep(1)
    })

    run_status <- attr(output_GSE, "status")
    if (is.null(run_status)) run_status <- 0L
    output_GSE <- paste(output_GSE, collapse = "\n")

    paths <- result_paths()
    success <- identical(as.integer(run_status), 0L) &&
      grepl("Algorithm execution completed", output_GSE) &&
      is_fresh(paths$csv)

    if (!success) {
      run_state(FALSE)
      err_line <- regmatches(output_GSE,
                             regexpr("GENOME_GMM_ERROR:[^\n]*", output_GSE))
      if (length(err_line) > 0 && nchar(err_line) > 0) {
        err_msg <- gsub("GENOME_GMM_ERROR:\\s*", "", err_line)
      } else {
        r_err <- grep("Error", strsplit(output_GSE, "\n", fixed = TRUE)[[1]],
                      value = TRUE)
        err_msg <- if (length(r_err) > 0) tail(r_err, 1) else "please check your data"
      }
      showNotification(
        paste0("Run failed: ", err_msg,
               " No result is shown. Please check your data and parameters, then try again."),
        type = "error", duration = 20
      )
      session$sendCustomMessage(type = 'testmessage',
        message = 'The task failed. Please check the error message and try again.')
      shinyjs::enable("submit")
      return()
    }


    run_state(TRUE)

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
      n_val <- switch(selected_species,
        "javanica_tetraploid" = 4,
        "potato_tetraploid" = 4,
        "wheat_hexaploid" = 6,
        "strawberry_octoploid" = 8,
        "Chinese_sturgeon_octoploid" = 8
      )
      if (!is.null(n_val)) {
        updateNumericInput(session, "n", value = n_val)
        updateActionButton(session, "example", label = paste("Use", selected_species, "data"))
      }

  })

  observeEvent(input$draw_plot,{
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
