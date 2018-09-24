## ui.R
fluidPage(
  # h1("Analysis and visualization of JUMP quantification results", align = "center"),
  # headerPanel("Analysis and visualization of JUMP quantification results"),
  
  navbarPage(
    "JUMP Quantification",
    tabPanel(
      "Exploratory data analysis",
      ## Sidebar panel controlling an exploratory data analysis
      sidebarPanel(
        p("In this panel, you can perform explorative data analyses for your proteomics dataset."),
        p("1. PCA (principal component analysis) plot of samples"),
        p("2. Hierarchical clustering result of samples and peptides/proteins"),
        p("3. Data table of highly variant peptides/proteins"),
        br(),
        fileInput(
          "inputFile1", 
          label = HTML("Choose a file<br />
                      id_uni_pep_quan.xlsx or<br />
                      id_uni_prot_quan.xlsx")
        ),
        numericInput(
          "variant1", 
          label = HTML("Proportion of highly variant elements<br />
                      (e.g. if you choose 10, top 10% of highly variant peptides/proteins
                      will be used)"), value = 10
        ),
        selectInput(
          "metric1",
          label = "Select the measure of variation",
          choice = list("Coefficient of variation (CV)" = 1, 
                        "Median absolute deviation (MAD)" = 2),
                        selected = 1
        ),
        actionButton("submit1", "Submit")
      ),
      
      ## Main panel showing the results of the exploratory data analysis
      mainPanel(
        tabsetPanel(
          tabPanel(
            "Principal component analysis (PCA)",
            br(),
            plotOutput("pcaPlot", height = "600px")
          ),
          tabPanel(
            "Heatmap of the subset of peptides/proteins",
            br(),
            plotOutput("hclustPlot", height = "700px", width = "500px")
          ),
          tabPanel(
            "Data table",
            br(),
            downloadButton("download1", "Download"),
            br(), br(),
            DT::dataTableOutput("dataTable1"),
            br(),
            plotOutput("plotDataTable1")
          )
        )
      )
    ),
    
    ## Differential expression analysis
    tabPanel(
      "Differential expression",
      ## Sidebar panel controlling an exploratory data analysis
      sidebarPanel(
        p("This panel provides differential expression analysis results for your proteomics dataset."),
        p("1. For two group comparison, a volcano plot and heatmap will be provided"),
        p("2. For three or more group comparison, a heatmap will be provided"),
        p("3. Data table of differentially expressed peptides/proteins"),
        br(),
        fileInput(
          "inputFile2",
          label = HTML("Choose a file<br />
                       id_uni_pep_quan.xlsx or<br />
                       id_uni_prot_quan.xlsx")
        ),
        numericInput(
          "numGroups2",
          label = "Number of groups",
          value = 2
        ),
        uiOutput("groups2"),
        selectInput(
          "metric2",
          label = "Select the measure of significance",
          choice = list("p-value" = "p-value", "FDR" = "FDR"),
          selected = 1
        ),
        numericInput(
          "cutoff2", label = "Significance level",
          min = 0, max = 1, step = 0.01, value = 0.05
        ),
        numericInput("logfc2", label = "Log2-fold cutoff", value = 1),
        actionButton("submit2", "Submit")
      ),
      
      ## Main panel showing the results of the differential expression analysis
      mainPanel(
        tabsetPanel(
          tabPanel(
            "Volcano plot",
            br(),
            plotOutput("volcanoPlot", height = "600px")
          ),
          tabPanel(
            "Heatmap of differentially expressed peptides/proteins",
            br(),
            plotOutput("hclustDE", height = "700px", width = "500px")
          ),
          tabPanel(
            "Data table",
            br(),
            downloadButton("download2", "Download"),
            br(), br(),
            DT::dataTableOutput("dataTable2"),
            br(),
            plotOutput("plotDataTable2")
          )
        )
      )
    )
  )
)
