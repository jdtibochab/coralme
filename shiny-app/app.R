library(shiny)
library(shinyFiles)

# Set maximum upload size (in bytes)
options(shiny.maxRequestSize = 50 * 1024^2)  # 50 MB

# Define UI for application
library(shiny)

ui <- fluidPage(

    # Application title
    titlePanel("coralME: COmprehensive Reconstruction ALgorithm for ME-models"),

    # Sidebar with file upload inputs
    sidebarLayout(
        tabsetPanel(
            sidebarPanel("Mandatory Inputs",
                     style = "font-size: 18px;",
                     fileInput("file1", "Choose M-model", multiple = FALSE, accept = c('.json', '.xml'), width = NULL, buttonLabel = "Browse...", placeholder = "No file selected", capture = NULL),
                     fileInput("file2", "Choose GenBank file", multiple = FALSE, accept = c('.gb', '.gbff'), width = NULL, buttonLabel = "Browse...", placeholder = "No file selected", capture = NULL),
                     fluidRow(
                         style = "display: flex; align-items: center;",
                         column(8, p("Run BLASTp:", style = "font-weight: bold;")),
                         column(4, uiOutput("toggleButtonUI1"))
                     ),
                     fluidRow(
                         style = "display: flex; align-items: center;",
                         column(8, p("Number of cores:", style = "font-weight: bold;")),
                         column(4, numericInput("num_cores_input", label = NULL, value = 1, min = 1, max = parallel::detectCores(), step = 1, width = '750px'))
                     ),
                     fluidRow(
                         style = "display: flex; align-items: center;",
                         column(8, p("Reference:", style = "font-weight: bold;")),
                         column(4, uiOutput("toggleButtonUI2"))
                     ),

                     shinyDirButton("UserReferenceDirectory", 'Select Directory', 'Upload'),

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
                     tags$div(style = "font-size: 18px; font-weight: bold;", "Logging directory"),
                     shinyDirButton("log_directory", 'Select Directory', 'Upload'),
                     tags$div(style = "font-size: 18px; font-weight: bold;", "Output directory"),
                     shinyDirButton("out_directory", 'Select Directory', 'Upload'),
                     tags$div(style = "font-size: 18px; font-weight: bold;", "Organism-Specific Matrix"),
                     shinySaveButton('saveData', 'Select Filename', 'Save', filetype = c('.xlsx'))

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

            textOutput("output_text"),
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
    shinyFileSave(input, 'saveData', roots=volumes, session=session)

    shinyFiles::shinyDirChoose(input, 'dir', roots=volumes, session=session)

    # Display the selected directory path
    output$log_directory <- renderText({
        if (is.null(input$log_directory)) {
            return("No directory selected")
        }

        parseDirPath(volumes, input$log_directory)
    })

    # Reactive value to store output file path
    output_file <- reactive({
        tempfile(fileext = ".txt")
    })

    # Options that can be true or false
    # Reactive value to store the toggle state
    toggleState1 <- reactiveVal(FALSE)
    toggleState2 <- reactiveVal(FALSE)
    toggleState3 <- reactiveVal(FALSE)
    toggleState4 <- reactiveVal(FALSE)
    toggleState5 <- reactiveVal(FALSE)

    # Observe event for "Toggle" button
    observeEvent(input$toggleButton1, { toggleState1(!toggleState1()) })
    observeEvent(input$toggleButton2, { toggleState2(!toggleState2()) })
    observeEvent(input$toggleButton3, { toggleState3(!toggleState3()) })
    observeEvent(input$toggleButton4, { toggleState4(!toggleState4()) })
    observeEvent(input$toggleButton5, { toggleState5(!toggleState5()) })

    # Dynamically render the button with updated label
    output$toggleButtonUI1 <- renderUI(ToggleTrueFalseLabelAndBackground(toggleState1(), "toggleButton1"))
    output$toggleButtonUI2 <- renderUI(ToggleTrueFalseLabelAndBackground(toggleState2(), "toggleButton2"))
    output$toggleButtonUI3 <- renderUI(ToggleTrueFalseLabelAndBackground(toggleState3(), "toggleButton3"))
    output$toggleButtonUI4 <- renderUI(ToggleTrueFalseLabelAndBackground(toggleState4(), "toggleButton4"))
    output$toggleButtonUI5 <- renderUI(ToggleTrueFalseLabelAndBackground(toggleState5(), "toggleButton5"))







    # Observe event for "Run Python Script" button
    observeEvent(input$run_script, {
        # Get file paths
        file1_path <- input$file1$datapath
        file2_path <- input$file2$datapath

        # Construct command to run Python script
        command <- "python3"
        args <- c("cli.py", file1_path, file2_path)

        # Execute the command and capture output
        script_output <- system2(command, args = args, stdout = TRUE)
        print(script_output)

        # Check for successful execution (assuming your script returns 0 on success)
        if (length(script_output) > 0 && script_output[length(script_output)] == "0") {
            output$output_text <- renderText("Python script executed successfully.")
        } else {
            output$output_text <- renderText("Python script execution failed.")
        }
    })

    # Download handler
    output$downloadData <- downloadHandler(
        filename = function() {
            basename(output_file())
        },
        content = function(file) {
            file.copy(output_file(), file)
        }
    )
}

# Run the application
shinyApp(ui = ui, server = server)
