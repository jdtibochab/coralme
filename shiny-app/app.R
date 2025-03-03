library(shiny)
library(shinyFiles)
library(zip)
# Set maximum upload size (in bytes)
options(shiny.maxRequestSize = 50 * 1024^2)  # 50 MB
# Define UI for application
ui <- fluidPage(
  
  # Application title
  titlePanel("coralME: COmprehensive Reconstruction ALgorithm for ME-models"),
  
  # Sidebar with file upload inputs
  sidebarLayout(
    tabsetPanel(
      sidebarPanel("Mandatory Inputs",
                   style = "font-size: 18px;",
                   fileInput("organism_json", "Choose organism.json", multiple = FALSE, accept = c('.json'), width = NULL, buttonLabel = "Browse...", placeholder = "No file selected", capture = NULL),
                   fileInput("file1", "Choose M-model", multiple = FALSE, accept = c('.json', '.xml'), width = NULL, buttonLabel = "Browse...", placeholder = "No file selected", capture = NULL),
                   fileInput("file2", "Choose GenBank file", multiple = FALSE, accept = c('.gb', '.gbff'), width = NULL, buttonLabel = "Browse...", placeholder = "No file selected", capture = NULL),
                   fluidRow(
                     style = "display: flex; align-items: center;",
                     column(8, p("Run BLASTp:", style = "font-weight: bold;")),
                     column(4, uiOutput("toggleButtonUI1"))
                   ),
                   
                   # Added E-value cutoff input
                   fluidRow(
                     style = "display: flex; align-items: center;",
                     column(8, p("E-value cutoff:", style = "font-weight: bold;")),
                     column(4, numericInput("e_value_cutoff", label = NULL, value = 0.001, min = 0, max = 1, step = 0.0001, width = '750px'))
                   ),
                   
                   # Modified locus tag selector
                   fluidRow(
                     style = "display: flex; align-items: center;",
                     column(8, p("Locus tag format:", style = "font-weight: bold;")),
                     column(4, selectInput("locus_tag_format", label = NULL, choices = c("new", "old"), selected = "new", width = '750px'))
                   ),
                   
                   fluidRow(
                     style = "display: flex; align-items: center;",
                     column(8, p("Number of cores:", style = "font-weight: bold;")),
                     column(4, numericInput("num_cores_input", label = NULL, value = 1, min = 1, max = parallel::detectCores(), step = 1, width = '750px'))
                   ),
                   
                   # Changed "Reference:" from toggle to file upload
                   fileInput("reference_file", "Reference:", multiple = FALSE, accept = NULL, width = NULL, buttonLabel = "Browse...", placeholder = "No file selected", capture = NULL),
                   
                   fluidRow(
                     style = "display: flex; align-items: center;",
                     column(8, p("Include pseudogenes:", style = "font-weight: bold;")),
                     column(4, uiOutput("toggleButtonUI3"))
                   ),
                   fluidRow(
                     style = "display: flex; align-items: center;",
                     column(8, p("Estimate Keffs:", style = "font-weight: bold;")),
                     column(4, uiOutput("toggleButtonUI4"))
                   ),
                   fluidRow(
                     style = "display: flex; align-items: center;",
                     column(8, p("Add lipoproteins:", style = "font-weight: bold;")),
                     column(4, uiOutput("toggleButtonUI5"))
                   ),
                   
                   br(), # Insert a line break
                   tags$div(style = "font-size: 18px;", "Mandatory Outputs"),
                   
                   # Changed to directory selection
                   tags$div(style = "font-size: 18px; font-weight: bold;", "Logging directory"),
                   shinyDirButton("log_directory", 'Select Directory', 'Select'),
                   verbatimTextOutput("log_dir_path"),
                   
                   # Changed to directory selection
                   tags$div(style = "font-size: 18px; font-weight: bold;", "Output directory"),
                   shinyDirButton("out_directory", 'Select Directory', 'Select'),
                   verbatimTextOutput("out_dir_path"),
                   
                   # Changed to file upload
                   fileInput("organism_matrix", "Organism-Specific Matrix", multiple = FALSE, accept = c('.xlsx'), width = NULL, buttonLabel = "Browse...", placeholder = "No file selected", capture = NULL)
      ),
      sidebarPanel("Optional inputs: Manual Curation",
                   style = "font-size: 18px;",
                   fileInput("file3", "Choose Transcription Units file", multiple = FALSE, accept = c('.txt'), width = NULL, buttonLabel = "Browse...", placeholder = "No file selected", capture = NULL),
                   fileInput("file4", "Choose Reaction file", multiple = FALSE, accept = c('.txt'), width = NULL, buttonLabel = "Browse...", placeholder = "No file selected", capture = NULL),
                   fileInput("file5", "Choose Subreactions file", multiple = FALSE, accept = c('.txt'), width = NULL, buttonLabel = "Browse...", placeholder = "No file selected", capture = NULL),
                   fileInput("file6", "Choose Reactions metadata file", multiple = FALSE, accept = c('.txt'), width = NULL, buttonLabel = "Browse...", placeholder = "No file selected", capture = NULL),
                   fileInput("file7", "Choose Metabolites metadata file", multiple = FALSE, accept = c('.txt'), width = NULL, buttonLabel = "Browse...", placeholder = "No file selected", capture = NULL)
      ),
      sidebarPanel("Optional inputs: BioCyc",
                   style = "font-size: 18px;",
                   fileInput("file8", "Choose BioCyc genes file", multiple = FALSE, accept = c('.txt'), width = NULL, buttonLabel = "Browse...", placeholder = "No file selected", capture = NULL),
                   fileInput("file9", "Choose BioCyc proteins file", multiple = FALSE, accept = c('.txt'), width = NULL, buttonLabel = "Browse...", placeholder = "No file selected", capture = NULL),
                   fileInput("file10", "Choose BioCyc TU file", multiple = FALSE, accept = c('.txt'), width = NULL, buttonLabel = "Browse...", placeholder = "No file selected", capture = NULL),
                   fileInput("file11", "Choose BioCyc RNA file", multiple = FALSE, accept = c('.txt'), width = NULL, buttonLabel = "Browse...", placeholder = "No file selected", capture = NULL),
                   fileInput("file12", "Choose BioCyc sequences file", multiple = FALSE, accept = c('.txt'), width = NULL, buttonLabel = "Browse...", placeholder = "No file selected", capture = NULL)
      )
    ),
    
    # Main panel to display output and download button
    mainPanel(
      actionButton("run_script", "Run coralME"),
      
      # Added output display area with scrolling
      tags$div(
        style = "height: 400px; overflow-y: auto; margin-top: 20px; border: 1px solid #ccc; padding: 10px;",
        verbatimTextOutput("output_text")
      ),
      
      downloadButton('downloadData', 'Download Output')
    )
  )
)
# helper functions
ToggleTrueFalseLabelAndBackground <- function(toggleState, label) {
  if (toggleState) {
    actionButton(label, "True", style = "background-color: #337ab7; width: 75px; height: 30px;")
  } else {
    actionButton(label, "False", style = "background-color: #cccccc; width: 75px; height: 30px;")
  }
}
# Define server logic
server <- function(input, output, session) {
  # Define paths for shinyFiles
  volumes <- c(Home = fs::path_home(), WD = getwd())
  
  # Setup directory choosers
  shinyDirChoose(input, 'log_directory', roots = volumes, session = session)
  shinyDirChoose(input, 'out_directory', roots = volumes, session = session)
  
  # Display selected logging directory path
  output$log_dir_path <- renderText({
    if (is.null(input$log_directory)) {
      return("No directory selected")
    }
    parseDirPath(volumes, input$log_directory)
  })
  
  # Display selected output directory path
  output$out_dir_path <- renderText({
    if (is.null(input$out_directory)) {
      return("No directory selected")
    }
    parseDirPath(volumes, input$out_directory)
  })
  
  # Reactive value to store script output
  script_output_result <- reactiveVal("")
  
  # Options that can be true or false
  # Reactive value to store the toggle state
  toggleState1 <- reactiveVal(FALSE)
  toggleState3 <- reactiveVal(FALSE)
  toggleState4 <- reactiveVal(FALSE)
  toggleState5 <- reactiveVal(FALSE)
  
  # Observe event for "Toggle" button
  observeEvent(input$toggleButton1, { toggleState1(!toggleState1()) })
  observeEvent(input$toggleButton3, { toggleState3(!toggleState3()) })
  observeEvent(input$toggleButton4, { toggleState4(!toggleState4()) })
  observeEvent(input$toggleButton5, { toggleState5(!toggleState5()) })
  
  # Dynamically render the button with updated label
  output$toggleButtonUI1 <- renderUI(ToggleTrueFalseLabelAndBackground(toggleState1(), "toggleButton1"))
  output$toggleButtonUI3 <- renderUI(ToggleTrueFalseLabelAndBackground(toggleState3(), "toggleButton3"))
  output$toggleButtonUI4 <- renderUI(ToggleTrueFalseLabelAndBackground(toggleState4(), "toggleButton4"))
  output$toggleButtonUI5 <- renderUI(ToggleTrueFalseLabelAndBackground(toggleState5(), "toggleButton5"))
  
  # Helper function to get file path or NULL
  get_file_path <- function(file_input) {
    if (is.null(file_input)) {
      return(NULL)
    }
    return(file_input$datapath)
  }
  
  # Observe event for "Run Python Script" button
  observeEvent(input$run_script, {
    # Initialize output
    script_output_result("Running coralME...\n")
    # Get file paths (mandatory inputs)
    organism_json_path <- get_file_path(input$organism_json)
    file1_path <- get_file_path(input$file1)
    file2_path <- get_file_path(input$file2)
    
    # Check if mandatory files are provided
    if (is.null(organism_json_path) || is.null(file1_path) || is.null(file2_path)) {
      script_output_result("Error: organism.json, M-model, and GenBank files are required.")
      return()
    }
    
    # Get directory paths
    log_directory <- parseDirPath(volumes, input$log_directory)
    out_directory <- parseDirPath(volumes, input$out_directory)
    
    # Check if directories are selected
    if (length(log_directory) == 0 || length(out_directory) == 0) {
      script_output_result("Error: Log directory and output directory must be selected.")
      return()
    }
    
    # Construct command to run Python script
    if (.Platform$OS.type == "windows") { 
      command <- "C:\\Python312\\python.exe"
    } else {
      command <- "python3"
    }
    #command <- "C:\\Python312\\python.exe"
    args <- c(
      "cli.py"
    )
    
    # Add mandatory parameters
    args <- c(args, paste0("--organism-json=", organism_json_path))
    args <- c(args, paste0("--m_model_path=", file1_path))
    args <- c(args, paste0("--genbank_path=", file2_path))
    
    args <- c(args, paste0("--e-value=", input$e_value_cutoff))
    
    # Modified locus tag handling
    if (input$locus_tag_format == "new") {
      args <- c(args, "--locus-tag=locus_tag")
    } else {
      args <- c(args, "--locus-tag=old_locus_tag")
    }
    
    args <- c(args, paste0("--cores=", input$num_cores_input))
    
    # Add toggle states
    if(toggleState1()) args <- c(args, "--run-blastp")
    if(toggleState3()) args <- c(args, "--include-pseudogenes")
    if(toggleState4()) args <- c(args, "--estimate-keffs")
    if(toggleState5()) args <- c(args, "--add-lipoproteins")
    
    # Add directory paths
    args <- c(args, paste0("--log-directory=", log_directory))
    args <- c(args, paste0("--out-directory=", out_directory))
    
    # Add file paths if provided
    if (!is.null(input$reference_file)) {
      args <- c(args, paste0("--reference=", input$reference_file$datapath))
    }
    
    if (!is.null(input$organism_matrix)) {
      args <- c(args, paste0("--organism-matrix=", input$organism_matrix$datapath))
    }
    
    # Add optional manual curation files if provided
    if (!is.null(input$file3)) {
      args <- c(args, paste0("--tu-file=", input$file3$datapath))
    }
    
    if (!is.null(input$file4)) {
      args <- c(args, paste0("--reaction-file=", input$file4$datapath))
    }
    
    if (!is.null(input$file5)) {
      args <- c(args, paste0("--subreactions-file=", input$file5$datapath))
    }
    
    if (!is.null(input$file6)) {
      args <- c(args, paste0("--reactions-metadata=", input$file6$datapath))
    }
    
    if (!is.null(input$file7)) {
      args <- c(args, paste0("--metabolites-metadata=", input$file7$datapath))
    }
    
    # Add BioCyc files if provided
    if (!is.null(input$file8)) {
      args <- c(args, paste0("--biocyc-genes=", input$file8$datapath))
    }
    
    if (!is.null(input$file9)) {
      args <- c(args, paste0("--biocyc-proteins=", input$file9$datapath))
    }
    
    if (!is.null(input$file10)) {
      args <- c(args, paste0("--biocyc-tu=", input$file10$datapath))
    }
    
    if (!is.null(input$file11)) {
      args <- c(args, paste0("--biocyc-rna=", input$file11$datapath))
    }
    
    if (!is.null(input$file12)) {
      args <- c(args, paste0("--biocyc-sequences=", input$file12$datapath))
    }
    
    # Debug: Show command and arguments
    command_str <- paste(command, paste(args, collapse = " "))
    script_output_result(paste(script_output_result(), "Running command:\n", command_str, "\n\n"))
    
    tryCatch({
      # Execute the command with system2
      result <- system2(command, args = args, stdout = TRUE, stderr = TRUE)
      print(args)
      
      # Directly display the output without filtering
      script_output_result(paste(script_output_result(), paste(result, collapse = "\n")))
      
    }, error = function(e) {
      script_output_result(paste(script_output_result(), "\n\nError executing coralME:", e$message))
    })
  })
  
  # Render output text
  output$output_text <- renderText({
    script_output_result()
  })
  
  # Download handler for output
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("coralME-output-", format(Sys.time(), "%Y%m%d-%H%M%S"), ".zip", sep = "")
    },
    content = function(file) {
      # Get selected directories
      log_dir <- parseDirPath(volumes, input$log_directory)
      out_dir <- parseDirPath(volumes, input$out_directory)
      
      # Validate directories
      if (length(log_dir) == 0 || length(out_dir) == 0) {
        stop("Please select both log and output directories before downloading.")
      }
      if (!dir.exists(log_dir) || !dir.exists(out_dir)) {
        stop("Selected directories do not exist.")
      }
      
      # Create temporary staging area
      temp_dir <- tempfile()
      dir.create(temp_dir)
      on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)
      
      # Copy directories with preserved structure
      dir.create(file.path(temp_dir, "log"))
      dir.create(file.path(temp_dir, "output"))
      
      # Copy all contents including hidden files
      copy_results <- c(
        file.copy(list.files(log_dir, full.names = TRUE, all.files = TRUE, no.. = TRUE),
                  file.path(temp_dir, "log"), 
                  recursive = TRUE),
        file.copy(list.files(out_dir, full.names = TRUE, all.files = TRUE, no.. = TRUE),
                  file.path(temp_dir, "output"), 
                  recursive = TRUE)
      )
      
      if (!all(copy_results)) {
        stop("Error copying directory contents")
      }
      
      # Create zip file
      zip_file <- tempfile(fileext = ".zip")
      zip::zipr(zip_file, files = c("log", "output"), root = temp_dir)
      
      # Move to final download location
      file.rename(zip_file, file)
    },
    contentType = "application/zip"
  )
}
# Run the application
shinyApp(ui = ui, server = server)