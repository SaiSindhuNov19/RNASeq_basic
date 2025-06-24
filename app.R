# app.R
library(shiny)
library(shinythemes)
library(DT)
#library(edgeR)
library(DESeq2)
#library(limma)
library(pheatmap)
library(ggplot2)
library(plotly)
library(RColorBrewer)
library(rmarkdown)
library(htmltools)

source("R/deg_functions.R")
source("R/plot_functions.R")

# UI
ui <- fluidPage(
  theme = shinytheme("yeti"),
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "style.css"),
    tags$script(src = "https://cdn.jsdelivr.net/npm/chart.js")
  ),
  
  titlePanel(
  title = "RNA-Seq Differential Expression Analyzer",
  windowTitle = "RNA-Seq DEG Analyzer"
  ),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      fileInput("count_file", "Upload Count Matrix (CSV/RDS)", accept = c(".csv", ".rds")),
      fileInput("meta_file", "Upload Metadata (CSV)", accept = ".csv"),
      
      selectInput("method", "DEG Method", 
                  choices = c("DESeq2" = "deseq2")),
      
      uiOutput("condition_ui"),
      uiOutput("reference_ui"),
      uiOutput("comparison_ui"),
      
      numericInput("pval_cutoff", "Adjusted p-value cutoff", value = 0.05, min = 0, max = 1, step = 0.01),
      numericInput("lfc_cutoff", "Log2FC cutoff", value = 1, min = 0, step = 0.1),
      
      actionButton("run_analysis", "Run Analysis", class = "btn-primary"),
      br(), br(),
      downloadButton("download_report", "Download Report"),
      downloadButton("download_results", "Download Results")
    ),
    
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "main_tabs",
        tabPanel("Summary", 
                 verbatimTextOutput("summary_text"),
                 plotlyOutput("count_distribution")),
        
        tabPanel("DEG Results", 
                 DT::dataTableOutput("deg_table")),
        
        tabPanel("Volcano Plot", 
                 plotlyOutput("volcano_plot", height = "600px")),
        
        tabPanel("PCA", 
                 plotlyOutput("pca_plot", height = "600px")),
        
        tabPanel("MA Plot", 
                 plotlyOutput("ma_plot", height = "600px")),
        
        tabPanel("Package Info",
                 verbatimTextOutput("package_info"))
      )
    )
  )
)

# Server
server <- function(input, output, session) {
  # Reactive values
  rv <- reactiveValues(
    counts = NULL,
    metadata = NULL,
    dds = NULL,
    deg_results = NULL,
    norm_counts = NULL
  )
  
  # Load count data
  observeEvent(input$count_file, {
    req(input$count_file)
    
    if (tools::file_ext(input$count_file$name) == "rds") {
      rv$counts <- readRDS(input$count_file$datapath)
    } else {
      rv$counts <- as.matrix(read.csv(input$count_file$datapath, row.names = 1))
    }
  })
  
  # Load metadata
  observeEvent(input$meta_file, {
    req(input$meta_file)
    rv$metadata <- read.csv(input$meta_file$datapath, row.names = 1)
  })
  
  # Dynamic UI for condition selection
  output$condition_ui <- renderUI({
    req(rv$metadata)
    selectInput("condition", "Select Condition Column", choices = colnames(rv$metadata))
  })
  
  # Dynamic UI for reference level
  output$reference_ui <- renderUI({
    req(input$condition, rv$metadata)
    levels <- unique(rv$metadata[[input$condition]])
    selectInput("reference", "Reference Level", choices = levels)
  })
  
  # Dynamic UI for comparison level
  output$comparison_ui <- renderUI({
    req(input$condition, rv$metadata, input$reference)
    levels <- unique(rv$metadata[[input$condition]])
    levels <- setdiff(levels, input$reference)
    selectInput("comparison", "Comparison Level", choices = levels)
  })
  
  # Run DEG analysis
  observeEvent(input$run_analysis, {
    req(rv$counts, rv$metadata, input$condition, input$reference, input$comparison)
    
    withProgress(message = 'Running DEG analysis...', value = 0, {
      # Filter samples present in both counts and metadata
      common_samples <- intersect(colnames(rv$counts), rownames(rv$metadata))
      counts_filtered <- rv$counts[, common_samples]
      metadata_filtered <- rv$metadata[common_samples, , drop = FALSE]
      
      # Create design vector
      group <- factor(metadata_filtered[[input$condition]])
      group <- relevel(group, ref = input$reference)
      
      incProgress(0.2, detail = "Normalizing data...")
      
      # Run selected method
      if (input$method == "deseq2") {
        result <- run_deseq2(counts_filtered, group, input$comparison)
      }
      
      incProgress(0.8, detail = "Preparing results...")
      
      rv$deg_results <- result$deg_results
      rv$norm_counts <- result$norm_counts
      rv$dds <- result$dds
    })
  })
  
  # Summary text
  output$summary_text <- renderText({
    req(rv$deg_results)
    
    total_genes <- nrow(rv$deg_results)
    sig_genes <- sum(rv$deg_results$padj < input$pval_cutoff & 
                       abs(rv$deg_results$log2FoldChange) > input$lfc_cutoff, na.rm = TRUE)
    up_genes <- sum(rv$deg_results$padj < input$pval_cutoff & 
                      rv$deg_results$log2FoldChange > input$lfc_cutoff, na.rm = TRUE)
    down_genes <- sum(rv$deg_results$padj < input$pval_cutoff & 
                        rv$deg_results$log2FoldChange < -input$lfc_cutoff, na.rm = TRUE)
    
    paste(
      "Differential Expression Analysis Summary\n",
      "----------------------------------------\n",
      "Method: ", toupper(input$method), "\n",
      "Comparison: ", input$comparison, " vs ", input$reference, "\n",
      "Total genes: ", total_genes, "\n",
      "Significant genes (FDR < ", input$pval_cutoff, " & |log2FC| > ", input$lfc_cutoff, "): ", sig_genes, "\n",
      "Up-regulated: ", up_genes, "\n",
      "Down-regulated: ", down_genes, "\n",
      "----------------------------------------"
    )
  })
  
  # Count distribution plot
  output$count_distribution <- renderPlotly({
    req(rv$counts)
    
    counts_long <- reshape2::melt(log10(rv$counts + 1))
    colnames(counts_long) <- c("Gene", "Sample", "Value")
    
    p <- ggplot(counts_long, aes(x = Value)) +
      geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
      labs(title = "Log10 Count Distribution", x = "Log10(count + 1)", y = "Frequency") +
      theme_minimal()
    
    ggplotly(p)
  })
  
  # DEG table
  output$deg_table <- DT::renderDataTable({
    req(rv$deg_results)
    
    DT::datatable(
      rv$deg_results,
      extensions = 'Buttons',
      options = list(
        pageLength = 10,
        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel'),
        scrollX = TRUE
      ),
      rownames = TRUE,
      filter = 'top'
    ) %>% 
      DT::formatRound(columns = c("log2FoldChange", "pvalue", "padj"), digits = 4)
  })
  
  # Volcano plot
  output$volcano_plot <- renderPlotly({
    req(rv$deg_results)
    
    plot_volcano(rv$deg_results, 
                 pval_cutoff = input$pval_cutoff, 
                 lfc_cutoff = input$lfc_cutoff,
                 title = paste("Volcano Plot:", input$comparison, "vs", input$reference))
  })
  
  
  # PCA plot
  output$pca_plot <- renderPlotly({
    req(rv$norm_counts, rv$metadata, input$condition)
    
    plot_pca(rv$norm_counts, 
             metadata = rv$metadata,
             condition_col = input$condition,
             title = "PCA Plot")
  })
  
  # MA plot
  output$ma_plot <- renderPlotly({
    req(rv$deg_results)
    
    plot_ma(rv$deg_results, 
            pval_cutoff = input$pval_cutoff, 
            lfc_cutoff = input$lfc_cutoff,
            title = paste("MA Plot:", input$comparison, "vs", input$reference))
  })
  
  # Package info
  output$package_info <- renderPrint({
    pkgs <- c("shiny", "shinythemes", "DT", "edgeR", "DESeq2", "limma", 
              "pheatmap", "ggplot2", "plotly", "RColorBrewer", "rmarkdown")
    
    versions <- sapply(pkgs, function(pkg) {
      if (requireNamespace(pkg, quietly = TRUE)) {
        paste(pkg, packageVersion(pkg))
      } else {
        paste(pkg, "not installed")
      }
    })
    
    cat("Package Versions:\n")
    cat("----------------\n")
    cat(paste(versions, collapse = "\n"))
  })
  
  # Download report
  output$download_report <- downloadHandler(
    filename = function() {
      paste0("DEG_Report_", input$comparison, "_vs_", input$reference, "_", Sys.Date(), ".html")
    },
    content = function(file) {
      req(rv$deg_results, rv$norm_counts, rv$metadata)
      
      params <- list(
        deg_results = rv$deg_results,
        norm_counts = rv$norm_counts,
        metadata = rv$metadata,
        method = input$method,
        comparison = paste(input$comparison, "vs", input$reference),
        pval_cutoff = input$pval_cutoff,
        lfc_cutoff = input$lfc_cutoff,
        condition_col = input$condition
      )
      
      rmarkdown::render(
        input = "reports/template.Rmd",
        output_file = file,
        params = params,
        envir = new.env(parent = globalenv())
      )
    }
  )
  
  # Download results
  output$download_results <- downloadHandler(
    filename = function() {
      paste0("DEG_Results_", input$comparison, "_vs_", input$reference, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(rv$deg_results)
      write.csv(rv$deg_results, file)
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)
