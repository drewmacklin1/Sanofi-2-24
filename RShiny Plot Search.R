#shiny app to find any plot of a gene, celltype, disease combination (can find multiple combinations at once)
library(shiny, lib.loc = "/opt/R/site_lib/RLib_Common/4.0.0")

setwd("~/cloud-data/cloud-pipeline-tim-tri-culture-storage/Drew_Macklin/")
all_genes <- read.csv("results_target_credentialing/all_genes.csv", row.names = 1)
celltype_df <- read.csv("results_target_credentialing/celltype_df.csv")
celltype<- lapply(celltype_df, function(x) x[x!=""])
disease_study_df <-read.csv("results_target_credentialing/disease_study_df.csv")
disease_study <- lapply(disease_study_df, function(x) x[x!=""])

makePlotContainers <- function(n, ...) {
  lst <- lapply(seq.int(n), function(i){
    plotOutput(sprintf('%s_%g', "plot", i), height=400)
  })
  lst <- lapply(split(lst, (seq.int(n)-1)%/%2), function(x) column(12/2, x))
  do.call(tagList, lst)
}

renderPlots <- function(input, output) {
  
  #all celltype names
  checkedCells <- c(input$Astrocytes, input$Neurons, input$Macrophages, input$Endothelial, input$Ependymal, input$Fibroblasts, input$OPCs, input$Lymphocytes, input$Microglia, input$Oligodendrocytes, input$Pericytes)
  associatedCells<-NULL
  associatedNames <- NULL
  for (i in 1:length(checkedCells)){
    if(checkedCells[i]){
      associatedNames <- celltype[[i]]
    }
    associatedCells <- c(associatedCells, associatedNames)
  }
  associatedCells <- unique(associatedCells)
  
  #all study names
  checkedDiseases <- c(input$AD, input$ALS, input$MS, input$PD)
  associatedStudies <-NULL
  associatedStudy <- NULL
  for (i in 1:length(checkedDiseases)){
    if(checkedDiseases[i]){
      associatedStudy <- disease_study[[i]]
    }
    associatedStudies <- c(associatedStudies, associatedStudy)
  }
  associatedStudies <- unique(associatedStudies)
  
  #all plot locations
  plotLocations<-NULL
  for (x in 1:length(associatedStudies)){
    for (y in 1:length(associatedCells)){
        if(file.exists(paste0('datasets/', associatedStudies[x], '/bar_plots/', associatedCells[y], "/", associatedCells[y], '_', input$gene,'_barplot.png'))){
          plotLocations<-c(plotLocations, as.character(paste0('datasets/', associatedStudies[x], '/bar_plots/', associatedCells[y], "/", associatedCells[y], '_', input$gene,'_barplot.png'))) 
        }
    }
  }

  for (i in 1:length(plotLocations)) {
    local({
      ii <- i
      output[[sprintf('%s_%g', "plot", ii)]] <- renderImage({ list(src = plotLocations[ii], alt = "Alternate Text", width = 400, height = 400)}, deleteFile = FALSE)
    })
  }
}

ui <- shinyUI(
  fluidPage(
    sidebarLayout(
      sidebarPanel(
        selectizeInput("gene", "Gene:", all_genes, ""),
        checkboxInput("Astrocytes", "Astrocytes", FALSE),
        checkboxInput("Neurons", "Neurons", FALSE),
        checkboxInput("Macrophages", "Macrophages", FALSE),
        checkboxInput("Endothelial", "Endothelial", FALSE),
        checkboxInput("Ependymal", "Ependymal", FALSE),
        checkboxInput("Fibroblasts", "Fibroblasts", FALSE),
        checkboxInput("OPCs", "OPCs", FALSE),
        checkboxInput("Lymphocytes", "Lymphocytes", FALSE),
        checkboxInput("Microglia", "Microglia", FALSE),
        checkboxInput("Oligodendrocytes", "Oligodendrocytes", FALSE),
        checkboxInput("Pericytes", "Pericytes", FALSE),
        checkboxInput("AD", "AD", FALSE),
        checkboxInput("ALS", "ALS", FALSE),
        checkboxInput("MS", "MS", FALSE),
        checkboxInput("PD", "PD", FALSE)
      ),
      mainPanel(
        uiOutput('plots')
      )
    )
  )
)

server <- shinyServer(function(input, output) {
  output$plots <- renderUI({
    makePlotContainers(100)  #Make this variable?
  })
  observe({
    renderPlots(input, output)
    })
})


shinyApp(ui, server)