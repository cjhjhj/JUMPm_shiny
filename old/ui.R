## ui.R
library(shiny)
fluidPage(
    # h1("Analysis and visualization of JUMP quantification results", align = "center"),
    # headerPanel("Analysis and visualization of JUMP quantification results"),
    navbarPage(
        "JUMPm: metabolomics data analysis",
        tabPanel(
            "Exploratory data analysis",
            ## Sidebar panel controlling an exploratory data analysis
            sidebarPanel(
                width = 3,
                tags$head(tags$style(HTML('h5 {margin-bottom:0px; margin-top:0px;}'))),
                p("In this panel, you can perform explorative data analyses for your metabolomics dataset.", br(),
                  "1. PCA (principal component analysis) plot of samples", br(),
                  "2. Hierarchical clustering result of samples and features", br(),
                  "3. Data table of highly variant features"),
                br(),
                fileInput("inputFile1", label = HTML("Choose a file<h5>e.g. example_fully_aligned.feature</h5>")),
                numericInput("variant1",
                             label = HTML("Proportion of highly variant elements<h5>(e.g. if you choose 10, top 10% of highly variant features will be used)</h5>"),
                             value = 10),
                selectInput("metric1", label = "Select the measure of variation",
                            choice = list("Coefficient of variation (CV)" = 1, "Median absolute deviation (MAD)" = 2),
                            selected = 1),
                actionButton("submit1", "Submit")
            ),
            
            ## Main panel showing the results of the exploratory data analysis
            mainPanel(
                tabsetPanel(
                    tabPanel("Principal component analysis (PCA)", br(), plotOutput("pcaPlot", height = "600px")),
                    tabPanel("Heatmap of the subset of features", br(), plotOutput("hclustPlot", height = "700px", width = "500px")),
                    tabPanel("Data table", br(), DT::dataTableOutput("dataTable1"), br(), downloadButton("download1", "Download"), br(), plotOutput("plotDataTable1"))
                )
            )
        ),
        
        ## Differential expression analysis
        tabPanel(
            "Differential expression",
            ## Sidebar panel controlling an exploratory data analysis
            sidebarPanel(
                width = 3,
                p("This panel provides differential expression analysis results for your metabolomics dataset.", br(),
                  "1. For two group comparison, a volcano plot and heatmap will be provided", br(),
                  "2. For three or more group comparison, a heatmap will be provided", br(),
                  "3. Data table of differentially expressed features"),
                br(),
                fileInput("inputFile2", label = HTML("Choose a file<h5>e.g. example_fully_aligned.feature</h5>")),
                numericInput("numGroups2", label = "Number of groups", value = 2),
                uiOutput("groups2"),
                selectInput("metric2", label = "Select the measure of significance",
                            choice = list("p-value" = "p-value", "FDR" = "FDR"), selected = 1),
                numericInput("cutoff2", label = "Significance level", min = 0, max = 1, step = 0.01, value = 0.05),
                numericInput("logfc2", label = "Log2-fold cutoff", value = 1),
                actionButton("submit2", "Submit")
            ),
            
            ## Main panel showing the results of the differential expression analysis
            mainPanel(
                tabsetPanel(
                    tabPanel("Volcano plot", br(), plotOutput("volcanoPlot", height = "600px")),
                    tabPanel("Heatmap of differentially expressed features", br(), plotOutput("hclustDE", height = "700px", width = "500px")),
                    tabPanel("Data table", br(), DT::dataTableOutput("dataTable2"), br(), downloadButton("download2", "Download"), br(), plotOutput("plotDataTable2"))
                )
            )
        )
    )
)
