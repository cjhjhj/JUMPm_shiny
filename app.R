library(shiny)
library(readxl)


## ui.R
ui = fluidPage(
  pageWithSidebar(
    headerPanel("Analaysis and visualization of JUMP quantification results"),
    sidebarPanel(
      h2("About"),
      p("This tool draws a volcano plot- typically used for showing relationships between fold change and p-values of genes (probes) in microarray studies.",
        "Points are sized according to their X-value, and colour intensity is by their Y-value"), 
      p("Which genes to labeSl are chosen based on Euclidean distance from the origin, but can be biased towards the X or Y dimension. All other settings should be self explanatory.",
        "Editing any settings will cause the plot to be redrawn, label positions are fixed above the point currently."),
      fileInput("inputFile", "Choose File"),
      # accept=c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
      numericInput("nGroups", "Number of groups in samples", value = 2, min = 2)
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("test", uiOutput("table")),
        # tabPanel("Select samples for each group",
        #         uiOutput("sampleGroups"),
        #         actionButton("submitGroups", "Update"),
        #         uiOutput("DE")),
        tabPanel("Volcano Plot", plotOutput(outputId = 'volcano', width = "800px", height = "800px")),
        tabPanel("Histogram", plotOutput(outputId = 'histogram', width = "100%")),
        tabPanel("Data Preview", tableOutput("contents"))
      )
    )
  )
)


server = function (input, output) {
  ## Load JUMP -q output file (either id_uni_pep_quan.xlsx or id_uni_prot_quan.xlsx)
  inputData = eventReactive(input$inputFile, {
    inFileName = input$inputFile$name
    if (length(grep("pep", inFileName))) {
      tbl = read_excel(inFileName, skip = 4) # Peptide publication table
    } else {
      tbl = read_excel(inFileName, skip = 1) # Protein publication table
    }
    as.data.frame(tbl)
  })
  
  # eventReactive(input$inputFile, {
  #   inputNumberGroups = reactive(as.integer(input$nGroups))
  #   print(inputNumberGroups)
  # })
  
  # ## Specificiation of sample groups
  # eventReactive(!is.null(input$inputFile), {
  #   inputNumberGroups = reactive(as.integer(input$nGroups))
  #   output$sampleGroups = renderUI({
  #     nGroups = inputNumberGroups()
  #     data = inputData()
  #     colSampleNames = grep('sig', colnames(data))
  #     sampleNames = colnames(data)[colSampleNames]
  #     lapply (1:nGroups, function(i) {
  #       checkboxGroupInput(inputId = paste0("Group", i), label = paste("Group", i),
  #                          choiceNames = as.list(sampleNames), choiceValues = as.list(sampleNames))
  #     })
  #   })
  # })
  
  # output$DE = renderText({
  #   nGroups = inputNumberGroups()
  #   comparison = as.character()
  #   for (g in 1:nGroups) {
  #     groupName = paste0("Group", g)
  #     comparison[g] = paste(input[[groupName]], collapse = ",")
  #   }
  #   paste0(comparison)
  # })
  
  # # Comparison between groups specified by the user input
  # DE = eventReactive(input$submitGroups, {
  #   ## Load the data
  #   data = inputData()
  #   nGroups = inputNumberGroups()
  # 
  #   ## Get the information of samples to be compared
  #   comparison = as.character()
  #   for (g in 1:nGroups) {
  #     groupName = paste0("Group", g)
  #     comparison[g] = paste(input[[groupName]], collapse = ",")
  #   }
  #   groups = list()
  #   nSamples = 0
  #   compSamples = NULL
  #   for (g in 1:nGroups) {
  #     groups[[g]] = unlist(strsplit(comparison[g], ","))
  #     nSamples = nSamples + length(groups[[g]])
  #     compSamples = c(compSamples, groups[[g]])
  #   }
  # 
  #   ## Generate a design matrix (which contains the information of comparison)
  #   subColNames = colnames(data)[which(colnames(data) %in% compSamples)]
  #   design = matrix(0, nrow = nSamples, ncol = nGroups)
  #   for (g in 1:nGroups) {
  #     design[which(subColNames %in% groups[[g]]), g] = 1
  #   }
  #   colnames(design) = paste("group", seq(1, nGroups), sep = "")
  # 
  #   ## Generate a contrast matrix and new column names for the LIMMA result table
  #   contVec = NULL
  #   newColNames = NULL
  #   combMatrix = combn(seq(1, nGroups), 2)
  #   for (j in 1:ncol(combMatrix)) {
  #     contVec = c(contVec, paste(paste("group", combMatrix[1, j], sep = ""), paste("group", combMatrix[2, j], sep = ""), sep = "-"))
  #     newColNames = c(newColNames, paste(comparison[combMatrix[1, j]], "/", comparison[combMatrix[2, j]], sep = ""))
  #   }
  #   contMatrix = makeContrasts(contrasts = contVec, levels = design)
  # })
  
  # output$selectedSamples = renderText({
  #   groupName = paste0("Group", 1)
  #   paste(input[[groupName]])
  # })
  
  # output$txt <- renderText({
  #   icons <- paste(input$icons, collapse = ", ")
  #   paste("You chose", icons)
  # })
  
  
  output$table = renderTable(head(inputData()))
  
}

shinyApp(ui = ui, server = server)