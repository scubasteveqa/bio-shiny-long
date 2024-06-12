# Load necessary libraries
library(shiny)
library(BiocManager)
library(GenomicRanges)
library(SummarizedExperiment)

# Define UI for the application
ui <- fluidPage(
  titlePanel("Bioconductor Shiny App with GenomicRanges and SummarizedExperiment"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("file1", "Upload Genomic Ranges CSV File",
                accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
      fileInput("file2", "Upload Summarized Experiment CSV File",
                accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
      actionButton("loadData", "Load Data")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("GenomicRanges", tableOutput("grangesTable")),
        tabPanel("SummarizedExperiment", tableOutput("summarizedExpTable"))
      )
    )
  )
)

# Define server logic for the application
server <- function(input, output) {
  
  grangesData <- reactiveVal()
  summarizedExpData <- reactiveVal()
  
  observeEvent(input$loadData, {
    req(input$file1)
    req(input$file2)
    
    grangesFile <- input$file1$datapath
    summarizedExpFile <- input$file2$datapath
    
    granges_df <- read.csv(grangesFile)
    granges_obj <- GRanges(
      seqnames = Rle(granges_df$seqnames),
      ranges = IRanges(start = granges_df$start, end = granges_df$end),
      strand = Rle(strand(granges_df$strand)),
      score = granges_df$score,
      GC = granges_df$GC
    )
    grangesData(granges_obj)
    
    summarizedExp_df <- read.csv(summarizedExpFile)
    counts <- as.matrix(summarizedExp_df[, -1])
    colData <- DataFrame(row.names = colnames(counts))
    se_obj <- SummarizedExperiment(assays = list(counts = counts), colData = colData)
    summarizedExpData(se_obj)
  })
  
  output$grangesTable <- renderTable({
    req(grangesData())
    as.data.frame(grangesData())
  })
  
  output$summarizedExpTable <- renderTable({
    req(summarizedExpData())
    as.data.frame(assay(summarizedExpData()))
  })
}

# Run the application
shinyApp(ui = ui, server = server)
