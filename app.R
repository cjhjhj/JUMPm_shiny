rm(list = ls())

library(shiny)
library(ggplot2)
library(gplots)
library(DT)

source("statTest.R")

server = function (input, output) {
    # Increase the maximum size of uploaded file (up to 30MB)
    options(shiny.maxRequestSize = 30 * (1024 ^ 2))
    
    ####################################################
    # Unsupervised analysis, i.e. explorative analysis #
    ####################################################
    # Load fully aligned feature table (xxx_fully_aligned.feature)
    data1 = reactive ({
        list(data = read.table(input$inputFile1$datapath, header = T, sep = "\t", check.names = F, comment.char = ""))
    })
    
    # Selection of a data subset (highly variable) for exploratory analysis
    subData1 = eventReactive(input$submit1, {
        rawData = data1()$data
        entry = paste0("feature", seq(1, dim(rawData)[1]))
        
        # CV or MAD calcuation is based on log2-transformed intensities, but output format is raw-intensity scale
        colInd = grep("intensity", tolower(colnames(rawData)))
        data = log(rawData[, colInd], 2)
        cv = apply(data, 1, sd) / rowMeans(data)
        mad = apply(abs(data - apply(data, 1, median)), 1, median)
        threshold = as.numeric(input$variant1)/ 100 ## Threshold percentage
        rowInd = NULL
        if (as.numeric(input$metric1) == 1) {
            rowInd = cv > quantile(cv, prob = 1 - threshold)
        } else if (as.numeric(input$metric1) == 2) {
            rowInd = mad > quantile(mad, prob = 1 - threshold)
        }
        
        # Return data for the following analyses
        # Column 1: entry (either protein accession or peptide sequence)
        # Column 2~ : log2-transformed intensity values for reporter ions
        entry = entry[rowInd]
        data = log(rawData[rowInd, colInd], 2)
        rownames(data) = entry
        return (data)
    })
    
    # PCA plot
    output$pcaPlot = renderPlot({
        data = subData1()
        colnames(data) = gsub("_Intensity", "", colnames(data))
        colnames(data) = gsub("_intensity", "", colnames(data))
            
        # Preparation of PCA result for visualization
        resPCA = prcomp(t(data), center = TRUE, scale = TRUE)
        eigs = resPCA$sdev ^ 2
        resPCA = data.frame(resPCA$x[, 1:2])
            
        # Parameter setup for ggplot
        xlabPCA = paste0("PC1 (", round((eigs[1] / sum(eigs)) * 100, 2),"%)")
        ylabPCA = paste0("PC2 (", round((eigs[2] / sum(eigs)) * 100, 2),"%)")
        ratioDisplay = 4/3
        ratioValue = (max(resPCA$PC1) - min(resPCA$PC1)) / (max(resPCA$PC2) - min(resPCA$PC2))
        g = ggplot(data = resPCA[, 1:2], aes(PC1, PC2)) +
            geom_jitter(size = 3) +
            geom_text(aes(label = rownames(resPCA)), vjust = "inward", hjust = "inward", size = 5) +
            labs(x = xlabPCA, y = ylabPCA) +
            coord_fixed(ratioValue / ratioDisplay) +
            theme(text = element_text(size = 12),
                  axis.text = element_text(size = 14),
                  axis.title = element_text(size = 14))
        plot(g)
    })
    
    # Heatmap and dendrogram
    output$hclustPlot = renderPlot({
        data = subData1()
        colnames(data) = gsub("_Intensity", "", colnames(data))
        colnames(data) = gsub("_intensity", "", colnames(data))
        mat = as.matrix(data)
        mat = t(scale(t(mat), center = T, scale = F)) # Only mean-centering
        limVal = round(min(abs(min(mat)), abs(max(mat))))
        myBreaks = seq(-limVal, limVal, length.out = 101)
        myColor <- colorRampPalette(c("blue", "white", "red"))(n = 100)
        par(oma = c(10, 3, 1, 3), mar = c(1, 1, 1, 1))
        h = heatmap.2(x = mat, density.info = "n", trace = "n", labRow = F, col = myColor,
                      hclust = function(x) hclust(x, method = "ward.D2"),
                      lhei = c(1, 6.5), lwid = c(2, 10), breaks = myBreaks,
                      key.par = list(mar= c(5, 0, 0, 0)), key.title = NA,
                      key.xlab = "scaled intensity")
    })
    
    # Data table
    output$dataTable1 = DT::renderDataTable({
        data = subData1()
        data = round(2 ** data, digits = 2)    # Transform back to the raw scale for visualization
    }, selection = 'single', options = list(scrollX = TRUE, pageLength = 5))
    
    # Plot of the selected rows from the data table
    output$plotDataTable1 = renderPlot({
        data = subData1()
        colnames(data) = gsub("_Intensity", "", colnames(data))
        colnames(data) = gsub("_intensity", "", colnames(data))
        data = round(2 ** data, digits = 2)    # Transform back to the raw scale for visualization
        rowInd = input$dataTable1_rows_selected
        if (length(rowInd) == 1) {
            x = as.numeric(data[rowInd, ])
            df = data.frame(samples = colnames(data), intensity = x)
            g = ggplot(df, aes(x = samples, y = intensity)) + 
                geom_bar(stat = "identity") + 
                theme(text = element_text(size = 15),
                      axis.text.x = element_text(angle = 90, hjust = 1)) + 
                scale_x_discrete(limits = colnames(data)) +
                coord_cartesian(ylim = c(0.8 * min(x), max(x)))
            plot(g)
        } else {
            rowInd = 1
            x = as.numeric(data[rowInd, ])
            df = data.frame(samples = colnames(data), intensity = x)
            g = ggplot(df, aes(x = samples, y = intensity)) + 
                geom_bar(stat = "identity") + 
                theme(text = element_text(size = 15), 
                      axis.text.x = element_text(angle = 90, hjust = 1)) + 
                scale_x_discrete(limits = colnames(data)) +
                coord_cartesian(ylim = c(0.8 * min(x), max(x)))
            plot(g)
        }
    })
    
    # Download the subset of data (exploratory analysis)
    output$download1= downloadHandler(
        filename = "exploratory_subset.txt",
        content = function(file) {
            write.table(subData1_for_download(), file, sep = "\t", row.names = FALSE)
        }
    )
    
    ##############################################################
    # Supervised analysis, i.e. differential expression analysis #
    ##############################################################
    # Load fully aligned feature table (xxx_fully_aligned.feature)
    data2 = reactive ({
        list(data = read.table(input$inputFile2$datapath, header = T, sep = "\t", check.names = F, comment.char = ""))
    })
    
    # Specificiation of groups of samples
    nGroups = reactive(as.integer(input$numGroups2))
    observeEvent(input$inputFile2, {
        output$groups2 = renderUI({
            data = data2()$data
            nGroups = nGroups()
            colSampleNames = grep('intensity', tolower(colnames(data)))
            sampleNames = colnames(data)[colSampleNames]
            lapply (1:nGroups, function(i) {
                checkboxGroupInput(inputId = paste0("Group", i), label = paste("Group", i),
                                   choiceNames = as.list(sampleNames), choiceValues = as.list(sampleNames))
            })
        })
    })
    
    # Differentially expressed peptides/proteins
    statRes = eventReactive(input$submit2, {
        data = data2()$data
        entry = paste0("feature", seq(1, dim(data)[1]))
        nGroups = nGroups()
        comparison = as.character()
        compSamples = as.character()
        for (g in 1:nGroups) {
            groupName = paste0("Group", g)
            comparison[g] = paste(input[[groupName]], collapse = ",")
        }
        groups = list()
        compSamples = NULL
        for (g in 1:nGroups) {
            groups[[g]] = unlist(strsplit(comparison[g], ","))
            compSamples = c(compSamples, groups[[g]])
        }
        statTest(data, entry, comparison)
    })
    
    # Prepare a set of differentiall expressed peptides/proteins from the "statRes" result
    subData2 = eventReactive(input$submit2, {
        nGroups = nGroups()
        statRes = statRes()
        data = statRes$data
        logFC = input$logfc2
        sigMetric = input$metric2
        sigCutoff = input$cutoff2
        resLogFC = statRes$res[, grep("Log2Fold", colnames(statRes$res))]
        if (nGroups > 2) {
            resLogFC = apply(cbind(abs(apply(resLogFC, 1, min)), abs(apply(resLogFC, 1, max))), 1, max)
        } else {
            resLogFC = abs(resLogFC)
        }
        rowInd = which(statRes$res[[sigMetric]] < sigCutoff & resLogFC >= logFC)
        data = data[rowInd, ]
    })
    
    subData2_for_download = eventReactive(input$submit2, {
        # For downloading, "data" should be a subset of raw data
        data = data2()$data
        
        # Row indices of differentially expressed peptides/proteins
        nGroups = nGroups()
        statRes = statRes()
        logFC = input$logfc2
        sigMetric = input$metric2
        sigCutoff = input$cutoff2
        resLogFC = statRes$res[, grep("Log2Fold", colnames(statRes$res))]
        if (nGroups > 2) {
            resLogFC = apply(cbind(abs(apply(resLogFC, 1, min)), abs(apply(resLogFC, 1, max))), 1, max)
        } else {
            resLogFC = abs(resLogFC)
        }
        rowInd = which(statRes$res[[sigMetric]] < sigCutoff & resLogFC >= logFC)
        
        # Re-organization of an output table
        colInd = max(grep('intensity', tolower(colnames(data))))
        data = cbind(data[rowInd, 1:colInd], statRes$res[rowInd, -1])
    })
    
    # Volcano plot of differential expression analysis - "statRes" is directly used
    output$volcanoPlot = renderPlot({
        res = statRes()$res
        logFC = input$logfc2
        sigMetric = input$metric2
        sigCutoff = input$cutoff2
        if (sigMetric == "p-value") {
            res = res[, 2:3]
            ylab = "-log10(p-value)"
        } else if (sigMetric == "FDR") {
            res = res[, c(2, 4)]
            ylab = "-log10(FDR)"
        }
        colnames(res) = c("logfc", "significance")
        res[, 2] = -log10(res[, 2])
        xlab = "log2(fold-change)"
        
        # Parameter setup for ggplot
        ratioDisplay = 4/3
        ratioValue = (max(res[, 1]) - min(res[, 1])) / (max(res[, 2]) - min(res[, 2]))
        g = ggplot(data = res, aes(logfc, significance)) +
            geom_point(alpha = 0.2, size = 2) + 
            geom_hline(aes(yintercept = -log10(sigCutoff))) + 
            geom_vline(aes(xintercept = -logFC)) + 
            geom_vline(aes(xintercept = logFC)) + 
            labs(x = xlab, y = ylab) +
            coord_fixed(ratioValue / ratioDisplay) + 
            theme(text = element_text(size = 20))
        plot(g)
    })
    
    # Heatmap of differentially expressed peptides/proteins
    output$hclustDE = renderPlot({
        data = subData2()
        mat = as.matrix(data)
        mat = t(scale(t(mat), center = T, scale = F)) # Only mean-centering
        limVal = round(min(abs(min(mat)), abs(max(mat))))
        myBreaks = seq(-limVal, limVal, length.out = 101)
        myColor <- colorRampPalette(c("blue", "white", "red"))(n = 100)
        par(oma = c(10, 3, 1, 3), mar = c(1, 1, 1, 1))
        h = heatmap.2(x = mat, density.info = "n", trace = "n", labRow = F, col = myColor,
                      hclust = function(x) hclust(x, method = "ward.D2"),
                      lhei = c(1, 6.5), lwid = c(2, 10), breaks = myBreaks,
                      key.par = list(mar= c(5, 0, 0, 0)), key.title = NA,
                      key.xlab = "scaled intensity")
    })
    
    # Data table of differentially expressed elements
    output$dataTable2 = DT::renderDataTable({
        # Since data is log2-transformed, it needs to be re-transformed to raw-scale intensity levels for visualization
        data = subData2()
        data = round(2 ** data, digits = 2)
    }, selection = 'single', options = list(scrollX = TRUE, pageLength = 5))
    
    # Plot of the selected rows from the data table
    output$plotDataTable2 = renderPlot({
        data = subData2()
        data = round(2 ** data, digits = 2)
        rowInd = input$dataTable2_rows_selected
        if (length(rowInd) == 1) {
            x = as.numeric(data[rowInd, ])
            df = data.frame(samples = colnames(data), intensity = x)
            g = ggplot(df, aes(x = samples, y = intensity)) + 
                geom_bar(stat = "identity") + 
                theme(text = element_text(size = 15),
                      axis.text.x = element_text(angle = 90, hjust = 1)) + 
                scale_x_discrete(limits = colnames(data)) +
                coord_cartesian(ylim = c(0.8 * min(x), max(x)))
            plot(g)
        } else {
            rowInd = 1
            x = as.numeric(data[rowInd, ])
            df = data.frame(samples = colnames(data), intensity = x)
            g = ggplot(df, aes(x = samples, y = intensity)) + 
                geom_bar(stat = "identity") + 
                theme(text = element_text(size = 15), 
                      axis.text.x = element_text(angle = 90, hjust = 1)) + 
                scale_x_discrete(limits = colnames(data)) +
                coord_cartesian(ylim = c(0.8 * min(x), max(x)))
            plot(g)
        }
    })
    
    # Download the subset of data
    output$download2= downloadHandler(
        filename = "differentially_expressed_subset.txt",
        content = function(file) {
            write.table(subData2_for_download(), file, sep = "\t", row.names = FALSE)
        }
    )
}


ui = fluidPage(
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

shinyApp(ui = ui, server = server)