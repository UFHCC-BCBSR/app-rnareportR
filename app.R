# app.R
library(shiny)
library(shinydashboard)
library(shinyFiles)
library(shinyjs)
library(DT)

# UI
ui <- fluidPage(
  useShinyjs(),
  tags$head(
    tags$style(HTML("
      body { background-color: #f4f4f4; }
      .title-section {
        text-align: center;
        padding: 15px 0;
        background-color: #2c3e50;
        color: white;
      }
      .step-section {
        background-color: white;
        margin: 10px;
        padding: 15px;
        border-radius: 8px;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
      }
      .btn-primary { background-color: #3498db; border-color: #3498db; }
      .btn-success { background-color: #27ae60; border-color: #27ae60; }
      .param-group {
        border: 1px solid #dee2e6;
        border-radius: 8px;
        padding: 10px;
        margin-bottom: 10px;
        background-color: #f8f9fa;
      }
      .param-group h4 {
        margin-top: 0;
        margin-bottom: 8px;
        color: #495057;
        border-bottom: 1px solid #dee2e6;
        padding-bottom: 5px;
      }
      .auth-status {
        padding: 8px;
        border-radius: 4px;
        margin-bottom: 10px;
      }
      .auth-success { background-color: #d4edda; border: 1px solid #c3e6cb; color: #155724; }
      .auth-needed { background-color: #fff3cd; border: 1px solid #ffeeba; color: #856404; }
      .selected-file {
        background-color: #f8f9fa;
        border: 1px solid #dee2e6;
        border-radius: 4px;
        padding: 6px 10px;
        margin-top: 6px;
        font-family: monospace;
        font-size: 11px;
        color: #495057;
        max-height: 60px;
        overflow-y: auto;
      }
      .validation-message {
        padding: 8px;
        border-radius: 4px;
        margin: 8px 0;
      }
      .validation-success { background-color: #d4edda; border: 1px solid #c3e6cb; color: #155724; }
      .validation-error { background-color: #f8d7da; border: 1px solid #f5c6cb; color: #721c24; }
      .info-box {
        background-color: #d1ecf1;
        border: 1px solid #bee5eb;
        border-radius: 4px;
        padding: 10px;
        margin: 10px 0;
        font-size: 13px;
        color: #0c5460;
      }
    "))
  ),
  
  div(class = "title-section",
      h1("RNA-seq Reporter", style = "margin: 0; font-size: 48px; font-weight: 300;"),
      p("Configure and run differential expression analysis reports",
        style = "margin: 10px 0 0 0; font-size: 14px; opacity: 0.8;")
  ),
  
  # Load Existing Parameters section
  div(class = "step-section",
      h3("Load Existing Parameters (Optional)", style = "text-align: center; margin-bottom: 20px;"),
      p("Start from an existing params.txt file to make modifications:",
        style = "text-align: center; color: #6c757d;"),
      
      # Option 1: Browse for existing params file
      div(style = "border: 1px solid #dee2e6; padding: 15px; margin: 10px 0; border-radius: 5px;",
          h5("Browse HiPerGator for Existing params.txt"),
          conditionalPanel(
            condition = "output.authenticated",
            textInput("custom_path_existing_params", "Directory Path:", value = "",
                      placeholder = "Enter path relative to volume..."),
            shinyFilesButton("browse_existing_params", "Browse Existing params.txt",
                             "Select params file", class = "btn-secondary", multiple = FALSE),
            uiOutput("selected_existing_params_file")
          ),
          conditionalPanel(
            condition = "!output.authenticated",
            div(style = "padding: 10px; text-align: center; color: #856404; font-size: 12px;",
                tags$i(class = "fa fa-lock"), " Login below to browse HiPerGator files"
            )
          )
      ),
      
      # Option 2: Upload params file
      div(style = "border: 1px solid #dee2e6; padding: 15px; margin: 10px 0; border-radius: 5px;",
          h5("Upload Existing params.txt"),
          fileInput("upload_existing_params", "Upload params.txt file", accept = ".txt"),
          uiOutput("uploaded_params_status")
      ),
      
      # Load button and status
      div(style = "text-align: center; margin-top: 20px;",
          actionButton("load_existing_params", "Load Parameters",
                       class = "btn-warning btn-lg", disabled = TRUE),
          br(), br(),
          uiOutput("params_load_status")
      )
  ),
  
  # Authentication section
  div(class = "step-section",
      h3("HiPerGator File Access", style = "text-align: center; margin-bottom: 20px;"),
      uiOutput("auth_status"),
      conditionalPanel(
        condition = "!output.authenticated",
        div(class = "login-box",
            h4("Login for HiPerGator File Access"),
            p("Enter your group credentials to browse HiPerGator files:"),
            fluidRow(
              column(6, textInput("group_name", "HiPerGator Group", placeholder = "e.g., cancercenter-dept")),
              column(6, passwordInput("group_password", "Password"))
            ),
            actionButton("login_btn", "Login", class = "btn-primary")
        )
      ),
      conditionalPanel(
        condition = "output.authenticated",
        div(style = "text-align: right;",
            actionButton("logout_btn", "Logout", class = "btn-secondary btn-sm")
        )
      )
  ),
  
  # Required Parameters
  div(class = "step-section",
      h2("Required Parameters", style = "text-align: center; margin-bottom: 30px;"),
      
      # Basic Configuration
      div(class = "param-group",
          h4("Basic Configuration"),
          fluidRow(
            column(6, textInput("sample_id", "Sample/Project ID", placeholder = "CX_liver_test")),
            column(6, textInput("report_title", "Report Title",
                                value = "RNA-seq Differential Expression Report"))
          ),
          fluidRow(
            column(4, selectInput("organism", "Organism",
                                  choices = list("Mouse (mmu)" = "mmu", "Human (hsa)" = "hsa"),
                                  selected = "mmu")),
            column(4, selectInput("annotation_db", "Annotation Database",
                                  choices = list("org.Mm.eg.db" = "org.Mm.eg.db",
                                                 "org.Hs.eg.db" = "org.Hs.eg.db"),
                                  selected = "org.Mm.eg.db")),
            column(4, textInput("hipergator_group", "HiPerGator Group",
                                placeholder = "e.g., cancercenter-dept"))
          ),
          textInput("output_path", "Output Path",
                    placeholder = "/blue/your-group/path/to/output"),
          textInput("user_email", "Email Address",
                    placeholder = "your.email@ufl.edu")
      ),
      
      # RNA-seq Specific Files
      div(class = "param-group",
          h4("RNA-seq Data Files"),
          
          # RSEM Directory
          tags$p(tags$strong("RSEM Directory"),
                 " - Directory containing .genes.results files from STAR RSEM"),
          conditionalPanel(
            condition = "output.authenticated",
            shinyDirButton("browse_rsem_dir", "Browse RSEM Directory",
                           "Select directory", class = "btn-info"),
            uiOutput("selected_rsem_dir")
          ),
          tags$br(),
          
          # Sample Sheet
          tags$p(tags$strong("Sample Sheet"), " - CSV file with sample metadata"),
          conditionalPanel(
            condition = "output.authenticated",
            fluidRow(
              column(6,
                     shinyFilesButton("browse_sample_sheet", "Browse Sample Sheet",
                                      "Select CSV file", class = "btn-info", multiple = FALSE),
                     uiOutput("selected_sample_sheet")
              ),
              column(6,
                     downloadButton("download_sample_sheet", "Download Selected File",
                                    class = "btn-secondary")
              )
            )
          ),
          fileInput("upload_sample_sheet", "OR Upload Sample Sheet", accept = ".csv"),
          tags$br(),
          
          # Contrasts File
          tags$p(tags$strong("Contrasts File"), " - Text file listing comparisons"),
          conditionalPanel(
            condition = "output.authenticated",
            shinyFilesButton("browse_contrasts", "Browse Contrasts File",
                             "Select file", class = "btn-info", multiple = FALSE),
            uiOutput("selected_contrasts")
          ),
          fileInput("upload_contrasts", "OR Upload Contrasts File", accept = ".txt")
      ),
      
      # Analysis Parameters - UPDATED SECTION
      div(class = "param-group",
          h4("Analysis Parameters"),
          
          # Differential Expression Method
          div(class = "info-box",
              tags$i(class = "fa fa-info-circle"),
              " Select the statistical method for differential expression analysis"
          ),
          selectInput("DE_tool", "Differential Expression Method",
                      choices = list(
                        "limma-voom (default)" = "limma_voom",
                        "DESeq2" = "deseq2",
                        "edgeR GLM" = "edger_GLM"
                      ),
                      selected = "limma_voom"),
          
          tags$br(),
          
          # Grouping and Batch Variables
          fluidRow(
            column(6, textInput("group_var", "Grouping Variable",
                                value = "Condition",
                                placeholder = "Column name in sample sheet")),
            column(6, textInput("batch_var", "Batch Variable (optional)",
                                value = "",
                                placeholder = "Leave empty if no batch correction"))
          ),
          
          tags$br(),
          
          # Gene Filtering Parameters
          div(class = "info-box",
              tags$i(class = "fa fa-filter"),
              " Gene filtering removes lowly expressed genes before analysis"
          ),
          h5("Gene Filtering Parameters", style = "margin-top: 15px;"),
          
          # Method selection
          fluidRow(
            column(12,
                   radioButtons("filter_method", "Filtering Method:",
                                choices = c("edgeR filterByExpr" = "edgeR",
                                            "NOISeq filtered.data" = "NOISeq"),
                                selected = "edgeR",
                                inline = TRUE)
            )
          ),
          
          # edgeR parameters (conditional)
          conditionalPanel(
            condition = "input.filter_method == 'edgeR'",
            fluidRow(
              column(6,
                     numericInput("filter_min_count", "Minimum count per sample",
                                  value = 10, min = 1, max = 100, step = 1)
              ),
              column(6,
                     numericInput("filter_min_prop", "Minimum sample proportion",
                                  value = 0.7, min = 0.1, max = 1.0, step = 0.1)
              )
            ),
            div(class = "help-text", style = "margin-top: 5px; font-size: 0.9em; color: #666;",
                "Genes must have ≥ min_count in ≥ min_prop of samples in at least one group")
          ),
          
          # NOISeq parameters (conditional)
          conditionalPanel(
            condition = "input.filter_method == 'NOISeq'",
            fluidRow(
              column(6,
                     numericInput("noiseq_method", "NOISeq Method",
                                  value = 1, min = 1, max = 3, step = 1)
              ),
              column(6,
                     numericInput("cv_cutoff", "CV Cutoff (%)",
                                  value = 100, min = 0, max = 500, step = 10)
              )
            ),
            fluidRow(
              column(6,
                     numericInput("cpm", "CPM Threshold",
                                  value = 1, min = 0, max = 10, step = 0.5)
              ),
              column(6,
                     selectInput("p_adj", "P-value Adjustment",
                                 choices = c("fdr", "bonferroni", "holm", "hochberg", "hommel", "BH", "BY"),
                                 selected = "fdr")
              )
            ),
            div(class = "help-text", style = "margin-top: 5px; font-size: 0.9em; color: #666;",
                "Method 1: CPM + CV filter, Method 2: Wilcoxon test, Method 3: Proportion test")
          ),
          
          tags$br(),
          
          # NEW: CPU cores for enrichment
          h5("Performance Settings", style = "margin-top: 15px;"),
          div(class = "info-box",
              tags$i(class = "fa fa-microchip"),
              " Number of CPU cores for GO/KEGG enrichment analysis (must match your SLURM allocation)"
          ),
          fluidRow(
            column(6,
                   numericInput("n_cores", "CPU cores",
                                value = 4, min = 1, max = 16, step = 1)
            ),
            column(6,
                   checkboxInput("filter_samples", "Filter low-quality samples", value = FALSE)
            )
          )
      )
  ),
  
  # Report Metadata (Optional)
  div(class = "step-section",
      h2("Report Metadata (Optional)", style = "text-align: center; margin-bottom: 30px;"),
      
      div(class = "param-group",
          h4("Project Information"),
          fluidRow(
            column(6, textInput("PI", "Principal Investigator", placeholder = "Dr. Jane Smith")),
            column(6, textInput("Institution", "Institution", placeholder = "University of Florida"))
          ),
          fluidRow(
            column(6, textInput("Department", "Department", placeholder = "Department of Medicine")),
            column(6, textInput("Study_Contact", "Study Contact", placeholder = "jane.smith@ufl.edu"))
          ),
          textInput("Project_Title", "Project Title", placeholder = "RNA-seq analysis of...")
      ),
      
      div(class = "param-group",
          h4("Study Details"),
          textAreaInput("Study_Summary", "Study Summary",
                        placeholder = "Brief description of the study objectives...",
                        rows = 3),
          fluidRow(
            column(6, textInput("Sample_Types", "Sample Type(s)",
                                placeholder = "Cell lines, tissue samples, etc.")),
            column(6, textInput("Analysis_Goals", "Analysis Goal(s)",
                                placeholder = "DE analysis, pathway analysis, etc."))
          )
      ),
      
      div(class = "param-group",
          h4("Report Credits"),
          fluidRow(
            column(6, textInput("Report_Prepared_By", "Report Prepared By",
                                placeholder = "Analyst Name")),
            column(6, textInput("Report_Reviewed_By", "Report Reviewed By",
                                placeholder = "Reviewer Name"))
          )
      )
  ),
  
  # Validation and Generation
  div(class = "step-section",
      h2("Generate Parameters & Submit", style = "text-align: center; margin-bottom: 30px;"),
      
      div(style = "text-align: center;",
          actionButton("validate_params", "Validate Parameters", class = "btn-secondary btn-lg"),
          br(), br(),
          actionButton("generate_params", "Generate params.txt", class = "btn-success btn-lg",
                       disabled = TRUE),
          br(), br(),
          downloadButton("download_params", "Download params.txt", class = "btn-info btn-lg",
                         disabled = TRUE),
          br(), br(),
          actionButton("submit_job", "Submit Analysis Job", class = "btn-primary btn-lg",
                       disabled = TRUE)
      ),
      br(),
      uiOutput("validation_status"),
      br(),
      verbatimTextOutput("params_preview"),
      br(),
      uiOutput("job_submission_status")
  )
)

# Server
server <- function(input, output, session) {
  
  values <- reactiveValues(
    authenticated = FALSE,
    current_group = NULL,
    selected_files = list(),
    params_valid = FALSE,
    params_generated = FALSE,
    validation_messages = c(),
    existing_params_file = NULL,
    params_loaded = FALSE         
  )
  
  # Group passwords
  group_passwords <- list(
    "cancercenter-dept" = "Abc123",
    "licht" = "licht123"
  )
  
  # Authentication
  output$authenticated <- reactive({ values$authenticated })
  outputOptions(output, "authenticated", suspendWhenHidden = FALSE)
  
  output$auth_status <- renderUI({
    if (values$authenticated) {
      div(class = "auth-status auth-success",
          tags$i(class = "fa fa-check-circle"),
          strong("Authenticated: "),
          paste("Logged in as", values$current_group),
          " - HiPerGator file browsing enabled"
      )
    } else {
      div(class = "auth-status auth-needed",
          tags$i(class = "fa fa-info-circle"),
          strong("HiPerGator Access: "),
          "Login required to browse HiPerGator files."
      )
    }
  })
  
  # Login handler
  observeEvent(input$login_btn, {
    req(input$group_name, input$group_password)
    
    if (input$group_name %in% names(group_passwords) &&
        input$group_password == group_passwords[[input$group_name]]) {
      
      values$authenticated <- TRUE
      values$current_group <- input$group_name
      
      # Set up file browser roots
      group_path <- paste0("/blue/", values$current_group)
      if (dir.exists(group_path)) {
        group_volumes <- setNames(group_path, values$current_group)
        
        # Initialize file browsers
        shinyDirChoose(input, "browse_rsem_dir", roots = group_volumes, session = session)
        shinyFileChoose(input, "browse_sample_sheet", roots = group_volumes,
                        session = session, filetypes = c("", "csv"))
        shinyFileChoose(input, "browse_contrasts", roots = group_volumes,
                        session = session, filetypes = c("", "txt"))
        shinyFileChoose(input, "browse_existing_params", roots = group_volumes,
                        session = session, filetypes = c("", "txt"))
      }
      
      showNotification("Login successful! HiPerGator file browsing enabled.", type = "message")
    } else {
      showNotification("Invalid credentials", type = "error")
    }
  })
  
  # Logout handler
  observeEvent(input$logout_btn, {
    values$authenticated <- FALSE
    values$current_group <- NULL
    values$selected_files <- list()
    showNotification("Logged out - file browsing disabled", type = "message")
  })
  
  # Parse params file function
  parse_params_file <- function(params_file) {
    if (!file.exists(params_file)) {
      stop("Parameters file not found")
    }
    
    lines <- readLines(params_file)
    params <- list()
    
    for (line in lines) {
      line <- trimws(line)
      if (line == "" || startsWith(line, "#")) {
        next
      }
      
      # Parse --param_name "value" or --param_name value
      if (startsWith(line, "--")) {
        space_pos <- regexpr("\\s", line)
        if (space_pos > 0) {
          param_name <- gsub("^--", "", substr(line, 1, space_pos - 1))
          param_value <- trimws(substr(line, space_pos + 1, nchar(line)))
          
          # Remove surrounding quotes if present
          if ((startsWith(param_value, '"') && endsWith(param_value, '"')) ||
              (startsWith(param_value, "'") && endsWith(param_value, "'"))) {
            param_value <- substr(param_value, 2, nchar(param_value) - 1)
          }
          
          params[[param_name]] <- param_value
        }
      }
    }
    
    return(params)
  }
  
  `%||%` <- function(x, y) if (is.null(x)) y else x
  
  # File browser for existing params
  observeEvent(input$browse_existing_params, {
    if (values$authenticated) {
      group_path <- paste0("/blue/", values$current_group)
      group_volumes <- setNames(group_path, values$current_group)
      
      tryCatch({
        if (!is.null(input$browse_existing_params) && length(input$browse_existing_params) > 0 &&
            !is.integer(input$browse_existing_params)) {
          selected_files <- parseFilePaths(group_volumes, input$browse_existing_params)
          if (!is.null(selected_files) && nrow(selected_files) > 0) {
            values$existing_params_file <- selected_files$datapath[1]
            shinyjs::enable("load_existing_params")
            showNotification("Params file selected - click Load Parameters to populate fields",
                             type = "message")
          }
        }
      }, error = function(e) {
        # Handle errors silently
      })
    }
  }, ignoreInit = TRUE)
  
  # Upload existing params
  observeEvent(input$upload_existing_params, {
    if (!is.null(input$upload_existing_params)) {
      permanent_dir <- file.path(dirname(getwd()), "uploaded_files")
      dir.create(permanent_dir, recursive = TRUE, showWarnings = FALSE)
      
      permanent_path <- file.path(permanent_dir, input$upload_existing_params$name)
      file.copy(input$upload_existing_params$datapath, permanent_path, overwrite = TRUE)
      
      values$existing_params_file <- permanent_path
      shinyjs::enable("load_existing_params")
      showNotification("Params file uploaded - click Load Parameters to populate fields",
                       type = "message")
    }
  })
  
  # Load existing parameters
  observeEvent(input$load_existing_params, {
    req(values$existing_params_file)
    
    tryCatch({
      if (!file.exists(values$existing_params_file)) {
        stop("File not found: ", values$existing_params_file)
      }
      
      params <- parse_params_file(values$existing_params_file)
      
      # ===== BASIC CONFIGURATION =====
      if ("sample_id" %in% names(params)) {
        updateTextInput(session, "sample_id", value = params$sample_id)
      }
      if ("report_title" %in% names(params)) {
        updateTextInput(session, "report_title", value = params$report_title)
      }
      if ("organism" %in% names(params)) {
        updateSelectInput(session, "organism", selected = params$organism)
      }
      if ("annotation_db" %in% names(params)) {
        updateSelectInput(session, "annotation_db", selected = params$annotation_db)
      }
      
      hipergator_group_key <- if ("hipergator-group" %in% names(params)) "hipergator-group" else "hipergator_group"
      if (hipergator_group_key %in% names(params)) {
        updateTextInput(session, "hipergator_group", value = params[[hipergator_group_key]])
      }
      
      output_path_key <- if ("output-path" %in% names(params)) "output-path" else "output_path"
      if (output_path_key %in% names(params)) {
        updateTextInput(session, "output_path", value = params[[output_path_key]])
      }
      
      if ("user_email" %in% names(params)) {
        updateTextInput(session, "user_email", value = params$user_email)
      }
      
      # ===== RNA-SEQ SPECIFIC =====
      if ("rsem_dir" %in% names(params)) {
        values$selected_files$rsem_dir <- params$rsem_dir
      }
      
      if ("group_var" %in% names(params)) {
        updateTextInput(session, "group_var", value = params$group_var)
      }
      
      if ("batch_var" %in% names(params)) {
        updateTextInput(session, "batch_var", value = params$batch_var)
      }
      
      # NEW: Load DE tool parameter
      if ("DE_tool" %in% names(params)) {
        updateSelectInput(session, "DE_tool", selected = params$DE_tool)
      }
      
      # NEW: Load filtering parameters
      # Filter method selection
      if ("filter_method" %in% names(params)) {
        updateRadioButtons(session, "filter_method", selected = params$filter_method)
      }
      
      # edgeR parameters
      if ("filter_min_count" %in% names(params)) {
        updateNumericInput(session, "filter_min_count", value = as.numeric(params$filter_min_count))
      }
      if ("filter_min_prop" %in% names(params)) {
        updateNumericInput(session, "filter_min_prop", value = as.numeric(params$filter_min_prop))
      }
      
      # NOISeq parameters
      if ("noiseq_method" %in% names(params)) {
        updateNumericInput(session, "noiseq_method", value = as.numeric(params$noiseq_method))
      }
      if ("cv_cutoff" %in% names(params)) {
        updateNumericInput(session, "cv_cutoff", value = as.numeric(params$cv_cutoff))
      }
      if ("cpm" %in% names(params)) {
        updateNumericInput(session, "cpm", value = as.numeric(params$cpm))
      }
      if ("p_adj" %in% names(params)) {
        updateSelectInput(session, "p_adj", selected = params$p_adj)
      }
      if ("n_cores" %in% names(params)) {
        updateNumericInput(session, "n_cores", value = as.numeric(params$n_cores))
      }
      
      if ("filter_samples" %in% names(params)) {
        updateCheckboxInput(session, "filter_samples", value = TRUE)
      }
      
      # ===== FILES =====
      if ("sample_data" %in% names(params)) {
        values$selected_files$sample_sheet <- params$sample_data
      }
      
      if ("contrasts" %in% names(params)) {
        if (file.exists(params$contrasts)) {
          values$selected_files$contrasts <- params$contrasts
        } else {
          updateTextInput(session, "contrasts_text", value = params$contrasts)
        }
      }
      
      # ===== OPTIONAL URLs =====
      if ("raw_seq_URL" %in% names(params)) {
        updateTextInput(session, "raw_seq_URL", value = params$raw_seq_URL)
      }
      if ("multiqc_url" %in% names(params)) {
        updateTextInput(session, "multiqc_url", value = params$multiqc_url)
      }
      
      # ===== REPORT METADATA =====
      if ("PI" %in% names(params)) {
        updateTextInput(session, "PI", value = params$PI)
      }
      if ("Institution" %in% names(params)) {
        updateTextInput(session, "Institution", value = params$Institution)
      }
      if ("Department" %in% names(params)) {
        updateTextInput(session, "Department", value = params$Department)
      }
      if ("Study_Contact" %in% names(params)) {
        updateTextInput(session, "Study_Contact", value = params$Study_Contact)
      }
      if ("Project_Title" %in% names(params)) {
        updateTextInput(session, "Project_Title", value = params$Project_Title)
      }
      if ("Study_Summary" %in% names(params)) {
        updateTextAreaInput(session, "Study_Summary", value = params$Study_Summary)
      }
      if ("Sample_Types" %in% names(params)) {
        updateTextInput(session, "Sample_Types", value = params$Sample_Types)
      }
      if ("Analysis_Goals" %in% names(params)) {
        updateTextInput(session, "Analysis_Goals", value = params$Analysis_Goals)
      }
      if ("Report_Prepared_By" %in% names(params)) {
        updateTextInput(session, "Report_Prepared_By", value = params$Report_Prepared_By)
      }
      if ("Report_Reviewed_By" %in% names(params)) {
        updateTextInput(session, "Report_Reviewed_By", value = params$Report_Reviewed_By)
      }
      
      values$params_loaded <- TRUE
      showNotification("Parameters loaded successfully from file!", type = "message")
      
    }, error = function(e) {
      showNotification(paste("Error loading parameters:", e$message), type = "error")
    })
  })
  
  # Display outputs for loaded params
  output$selected_existing_params_file <- renderUI({
    if (!is.null(values$selected_files$existing_params)) {
      div(class = "selected-file",
          strong("Selected: "), basename(values$selected_files$existing_params),
          tags$br(),
          tags$span(style = "font-size: 10px;", values$selected_files$existing_params)
      )
    }
  })
  
  output$uploaded_params_status <- renderUI({
    if (!is.null(values$existing_params_file)) {
      div(style = "background-color: #d4edda; padding: 10px; border-radius: 5px; margin-top: 10px;",
          tags$i(class = "fa fa-check-circle", style = "color: #155724;"),
          strong(" File ready: "), basename(values$existing_params_file)
      )
    }
  })
  
  output$params_load_status <- renderUI({
    if (!is.null(input$load_existing_params) && input$load_existing_params > 0) {
      div(style = "background-color: #d4edda; padding: 10px; border-radius: 5px; margin-top: 10px;",
          tags$i(class = "fa fa-info-circle", style = "color: #155724;"),
          strong("Parameters loaded! "), "You can now modify fields and regenerate."
      )
    }
  })
  
  # File selection handlers
  observeEvent(input$browse_rsem_dir, {
    if (values$authenticated) {
      group_path <- paste0("/blue/", values$current_group)
      group_volumes <- setNames(group_path, values$current_group)
      
      tryCatch({
        if (!is.null(input$browse_rsem_dir) && length(input$browse_rsem_dir) > 0) {
          selected_dir <- parseDirPath(group_volumes, input$browse_rsem_dir)
          if (length(selected_dir) > 0) {
            values$selected_files$rsem_dir <- selected_dir
            showNotification(paste("Selected RSEM directory:", basename(selected_dir)),
                             type = "message")
          }
        }
      }, error = function(e) {
        # Handle errors silently
      })
    }
  }, ignoreInit = TRUE)
  
  observeEvent(input$browse_sample_sheet, {
    if (values$authenticated) {
      group_path <- paste0("/blue/", values$current_group)
      group_volumes <- setNames(group_path, values$current_group)
      
      tryCatch({
        if (!is.null(input$browse_sample_sheet) && length(input$browse_sample_sheet) > 0 &&
            !is.integer(input$browse_sample_sheet)) {
          selected_files <- parseFilePaths(group_volumes, input$browse_sample_sheet)
          if (!is.null(selected_files) && nrow(selected_files) > 0) {
            values$selected_files$sample_sheet <- selected_files$datapath[1]
            showNotification(paste("Selected:", basename(selected_files$name[1])),
                             type = "message")
          }
        }
      }, error = function(e) {
        # Handle errors silently
      })
    }
  }, ignoreInit = TRUE)
  
  observeEvent(input$browse_contrasts, {
    if (values$authenticated) {
      group_path <- paste0("/blue/", values$current_group)
      group_volumes <- setNames(group_path, values$current_group)
      
      tryCatch({
        if (!is.null(input$browse_contrasts) && length(input$browse_contrasts) > 0 &&
            !is.integer(input$browse_contrasts)) {
          selected_files <- parseFilePaths(group_volumes, input$browse_contrasts)
          if (!is.null(selected_files) && nrow(selected_files) > 0) {
            values$selected_files$contrasts <- selected_files$datapath[1]
            showNotification(paste("Selected:", basename(selected_files$name[1])),
                             type = "message")
          }
        }
      }, error = function(e) {
        # Handle errors silently
      })
    }
  }, ignoreInit = TRUE)
  
  # Display selected files
  output$selected_rsem_dir <- renderUI({
    if (!is.null(values$selected_files$rsem_dir)) {
      div(class = "selected-file",
          strong("Selected: "), basename(values$selected_files$rsem_dir),
          tags$br(),
          tags$span(style = "font-size: 10px;", values$selected_files$rsem_dir)
      )
    }
  })
  
  output$selected_sample_sheet <- renderUI({
    if (!is.null(values$selected_files$sample_sheet)) {
      div(class = "selected-file",
          strong("Selected: "), basename(values$selected_files$sample_sheet),
          tags$br(),
          tags$span(style = "font-size: 10px;", values$selected_files$sample_sheet)
      )
    }
  })
  
  output$selected_contrasts <- renderUI({
    if (!is.null(values$selected_files$contrasts)) {
      div(class = "selected-file",
          strong("Selected: "), basename(values$selected_files$contrasts),
          tags$br(),
          tags$span(style = "font-size: 10px;", values$selected_files$contrasts)
      )
    }
  })
  
  # Download handler for sample sheet
  output$download_sample_sheet <- downloadHandler(
    filename = function() {
      if (!is.null(values$selected_files$sample_sheet)) {
        basename(values$selected_files$sample_sheet)
      } else {
        "sample_sheet.csv"
      }
    },
    content = function(file) {
      if (!is.null(values$selected_files$sample_sheet)) {
        file.copy(values$selected_files$sample_sheet, file)
      }
    }
  )
  
  # Handle uploaded sample sheet
  observeEvent(input$upload_sample_sheet, {
    if (!is.null(input$upload_sample_sheet)) {
      output_path <- input$output_path
      if (is.null(output_path) || output_path == "") {
        showNotification("Please set an output path before uploading files", type = "error")
        return()
      }
      
      if (!dir.exists(output_path)) {
        dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
      }
      
      output_file_path <- file.path(output_path, input$upload_sample_sheet$name)
      file.copy(input$upload_sample_sheet$datapath, output_file_path, overwrite = TRUE)
      values$selected_files$sample_sheet <- output_file_path
      showNotification(paste("File uploaded to:", output_file_path), type = "message")
    }
  })
  
  # Handle uploaded contrasts
  observeEvent(input$upload_contrasts, {
    if (!is.null(input$upload_contrasts)) {
      output_path <- input$output_path
      if (is.null(output_path) || output_path == "") {
        showNotification("Please set an output path before uploading files", type = "error")
        return()
      }
      
      if (!dir.exists(output_path)) {
        dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
      }
      
      output_file_path <- file.path(output_path, input$upload_contrasts$name)
      file.copy(input$upload_contrasts$datapath, output_file_path, overwrite = TRUE)
      values$selected_files$contrasts <- output_file_path
      showNotification(paste("File uploaded to:", output_file_path), type = "message")
    }
  })
  
  # Validation function
  validate_parameters <- function() {
    messages <- c()
    
    # Check required fields
    if (is.null(input$sample_id) || input$sample_id == "") {
      messages <- c(messages, "❌ Sample/Project ID is required")
    }
    if (is.null(input$hipergator_group) || input$hipergator_group == "") {
      messages <- c(messages, "❌ HiPerGator Group is required")
    }
    if (is.null(input$output_path) || input$output_path == "") {
      messages <- c(messages, "❌ Output path is required")
    }
    if (is.null(input$user_email) || input$user_email == "") {
      messages <- c(messages, "❌ User email is required")
    }
    if (is.null(input$group_var) || input$group_var == "") {
      messages <- c(messages, "❌ Grouping variable is required")
    }
    
    # Check file inputs
    if (is.null(values$selected_files$rsem_dir)) {
      messages <- c(messages, "❌ RSEM directory is required")
    }
    if (is.null(values$selected_files$sample_sheet)) {
      messages <- c(messages, "❌ Sample sheet is required")
    }
    if (is.null(values$selected_files$contrasts)) {
      messages <- c(messages, "❌ Contrasts file is required")
    }
    
    # If all required fields are present, add success message
    if (length(messages) == 0) {
      messages <- c(messages, "✅ All required parameters validated successfully")
    }
    
    values$validation_messages <- messages
    values$params_valid <- !any(grepl("^❌", messages))
    return(messages)
  }
  
  # Validation event
  observeEvent(input$validate_params, {
    messages <- validate_parameters()
    
    if (values$params_valid) {
      shinyjs::enable("generate_params")
    } else {
      shinyjs::disable("generate_params")
      values$params_generated <- FALSE
    }
  })
  
  # Display validation status
  output$validation_status <- renderUI({
    if (length(values$validation_messages) > 0) {
      message_class <- if (values$params_valid) "validation-success" else "validation-error"
      div(class = paste("validation-message", message_class),
          HTML(paste(values$validation_messages, collapse = "<br>"))
      )
    }
  })
  
  # Generate params.txt - UPDATED FUNCTION
  generate_params_content <- function() {
    lines <- c()
    # Basic parameters
    lines <- c(lines, paste("--report_title", shQuote(input$report_title)))
    lines <- c(lines, paste("--rsem_dir", shQuote(values$selected_files$rsem_dir)))
    lines <- c(lines, paste("--group_var", shQuote(input$group_var)))
    lines <- c(lines, paste("--output-path", shQuote(input$output_path)))
    # Optional batch variable
    if (!is.null(input$batch_var) && input$batch_var != "") {
      lines <- c(lines, paste("--batch_var", shQuote(input$batch_var)))
    }
    # NEW: DE tool parameter
    lines <- c(lines, paste("--DE_tool", shQuote(input$DE_tool)))
    
    # NEW: Filtering method and parameters
    lines <- c(lines, paste("--filter_method", shQuote(input$filter_method)))
    
    # edgeR filtering parameters
    lines <- c(lines, paste("--filter_min_count", input$filter_min_count))
    lines <- c(lines, paste("--filter_min_prop", input$filter_min_prop))
    
    # NOISeq filtering parameters
    lines <- c(lines, paste("--noiseq_method", input$noiseq_method))
    lines <- c(lines, paste("--cv_cutoff", input$cv_cutoff))
    lines <- c(lines, paste("--cpm", input$cpm))
    lines <- c(lines, paste("--p_adj", shQuote(input$p_adj)))
    
    lines <- c(lines, paste("--n_cores", input$n_cores))
    # Filter samples flag
    if (input$filter_samples) {
      lines <- c(lines, "--filter_samples")
    }
    # Sample sheet and contrasts
    lines <- c(lines, paste("--sample_data", shQuote(values$selected_files$sample_sheet)))
    lines <- c(lines, paste("--contrasts", shQuote(values$selected_files$contrasts)))
    # Organism and annotation
    lines <- c(lines, paste("--organism", shQuote(input$organism)))
    lines <- c(lines, paste("--annotation_db", shQuote(input$annotation_db)))
    # Report metadata (only add if not empty)
    if (!is.null(input$PI) && input$PI != "") {
      lines <- c(lines, paste("--PI", shQuote(input$PI)))
    }
    if (!is.null(input$Institution) && input$Institution != "") {
      lines <- c(lines, paste("--Institution", shQuote(input$Institution)))
    }
    if (!is.null(input$Department) && input$Department != "") {
      lines <- c(lines, paste("--Department", shQuote(input$Department)))
    }
    if (!is.null(input$Study_Contact) && input$Study_Contact != "") {
      lines <- c(lines, paste("--Study_Contact", shQuote(input$Study_Contact)))
    }
    if (!is.null(input$Project_Title) && input$Project_Title != "") {
      lines <- c(lines, paste("--Project_Title", shQuote(input$Project_Title)))
    }
    if (!is.null(input$Study_Summary) && input$Study_Summary != "") {
      lines <- c(lines, paste("--Study_Summary", shQuote(input$Study_Summary)))
    }
    if (!is.null(input$Sample_Types) && input$Sample_Types != "") {
      lines <- c(lines, paste("--Sample_Types", shQuote(input$Sample_Types)))
    }
    if (!is.null(input$Analysis_Goals) && input$Analysis_Goals != "") {
      lines <- c(lines, paste("--Analysis_Goals", shQuote(input$Analysis_Goals)))
    }
    if (!is.null(input$Report_Prepared_By) && input$Report_Prepared_By != "") {
      lines <- c(lines, paste("--Report_Prepared_By", shQuote(input$Report_Prepared_By)))
    }
    if (!is.null(input$Report_Reviewed_By) && input$Report_Reviewed_By != "") {
      lines <- c(lines, paste("--Report_Reviewed_By", shQuote(input$Report_Reviewed_By)))
    }
    # Optional URLs (empty strings by default)
    lines <- c(lines, '--raw_seq_URL ""')
    lines <- c(lines, '--multiqc_url ""')
    return(lines)
  }
  
  # Generate params file
  observeEvent(input$generate_params, {
    req(values$params_valid)
    
    tryCatch({
      # Ensure output directory exists
      if (!dir.exists(input$output_path)) {
        dir.create(input$output_path, recursive = TRUE, showWarnings = FALSE)
      }
      
      # Generate params file in output directory
      params_file <- file.path(input$output_path, paste0(input$sample_id, "_params.txt"))
      params_content <- generate_params_content()
      writeLines(params_content, params_file)
      
      # Store the path
      values$params_file_path <- params_file
      values$params_generated <- TRUE
      
      shinyjs::enable("submit_job")
      shinyjs::enable("download_params")
      
      showNotification(paste("Parameters file saved to HiPerGator:", params_file),
                       type = "message")
      
    }, error = function(e) {
      showNotification(paste("Error generating params file:", e$message), type = "error")
    })
  })
  
  # Display params preview
  output$params_preview <- renderText({
    if (values$params_generated && !is.null(values$params_file_path)) {
      paste("Generated params.txt:\n\n",
            paste(readLines(values$params_file_path), collapse = "\n"))
    } else {
      "Click 'Validate Parameters' then 'Generate params.txt' to see preview"
    }
  })
  
  # Download params.txt file
  output$download_params <- downloadHandler(
    filename = function() {
      paste0(input$sample_id, "_params.txt")
    },
    content = function(file) {
      if (values$params_generated && !is.null(values$params_file_path)) {
        file.copy(values$params_file_path, file)
      }
    }
  )
  
  # Enable/disable buttons based on validation state
  observe({
    if (values$params_valid && values$params_generated) {
      shinyjs::enable("download_params")
      shinyjs::enable("submit_job")
    } else {
      shinyjs::disable("download_params")
      shinyjs::disable("submit_job")
    }
  })
  
  # Submit job
  observeEvent(input$submit_job, {
    req(values$params_valid, values$params_generated)
    
    tryCatch({
      # Construct sbatch command
      sbatch_script <- "render-rnaseq-report.sbatch"
      sbatch_args <- c(
        sbatch_script,
        "--params-file", values$params_file_path,
        "--title", shQuote(input$report_title),
        "--output-dir", input$output_path,
        "--email", input$user_email
      )
      
      full_command <- paste("sbatch", paste(sbatch_args, collapse = " "))
      
      # Submit the job
      result <- system2("sbatch", args = sbatch_args, stdout = TRUE, stderr = TRUE)
      
      if (attr(result, "status") == 0 || is.null(attr(result, "status"))) {
        job_output <- paste(result, collapse = "\n")
        job_id_match <- regexpr("Submitted batch job ([0-9]+)", job_output)
        job_id <- if (job_id_match > 0) {
          gsub("Submitted batch job ", "", regmatches(job_output, job_id_match))
        } else {
          "Unknown"
        }
        
        # Show success notification
        showNotification(paste("RNA-seq analysis job submitted! Job ID:", job_id),
                         type = "message", duration = 10)
        
        # Update permanent UI
        output$job_submission_status <- renderUI({
          div(style = "background-color: #d4edda; padding: 15px; border-radius: 5px; margin-top: 15px;",
              tags$h4("✅ Job Submitted Successfully!", style = "color: #155724; margin-top: 0;"),
              tags$p(strong("Job ID: "), job_id),
              tags$p(strong("DE Method: "), input$DE_tool),
              tags$p(strong("Command: "), tags$code(full_command)),
              tags$p(strong("Parameters file: "), tags$code(values$params_file_path)),
              tags$p(strong("Expected output: "),
                     tags$code(file.path(input$output_path, paste0(input$sample_id, "_Report.html")))),
              tags$p(strong("Check job status: "), tags$code(paste0("squeue -j ", job_id))),
              tags$p(strong("View logs: "), tags$code(paste0("tail -f logs/rnaseq-report_", job_id, ".log")))
          )
        })
      } else {
        error_msg <- paste("SLURM submission failed:", paste(result, collapse = "\n"))
        
        output$job_submission_status <- renderUI({
          div(style = "background-color: #f8d7da; padding: 15px; border-radius: 5px; margin-top: 15px;",
              tags$h4("❌ Job Submission Failed", style = "color: #721c24; margin-top: 0;"),
              tags$p(strong("Error: "), error_msg),
              tags$p(strong("Command attempted: "), tags$code(full_command))
          )
        })
        
        showNotification(error_msg, type = "error", duration = 15)
      }
      
    }, error = function(e) {
      error_msg <- paste("Error submitting SLURM job:", e$message)
      
      output$job_submission_status <- renderUI({
        div(style = "background-color: #f8d7da; padding: 15px; border-radius: 5px; margin-top: 15px;",
            tags$h4("❌ Job Submission Error", style = "color: #721c24; margin-top: 0;"),
            tags$p(strong("Error: "), error_msg)
        )
      })
      
      showNotification(error_msg, type = "error")
    })
  })
}

shinyApp(ui = ui, server = server)