#shiny app that looks through the target credentialing folder for files of interest
library(shiny, lib.loc = "/opt/R/4.0.0/lib/R/site-library")
library(DT)

setwd("~/cloud-data/cloud-pipeline-tim-tri-culture-storage/Drew_Macklin")

program_list<-list.files(path="~/cloud-data/cloud-pipeline-tim-tri-culture-storage/Drew_Macklin/results_target_credentialing")

ui <- fluidPage(
  navbarPage("Single Cell Indication",
             tabPanel("Plot Search",
                      sidebarPanel(
                        selectizeInput("program", "Program:", program_list, ""),
                        selectizeInput("study", "Study:", NULL, ""),
                        selectizeInput("cell_type", "Cell Type:", NULL, ""),
                        selectizeInput("gene", "Gene:", NULL, "")
                      ),
                      mainPanel(
                        imageOutput("plot")
                      )
             ),
             tabPanel("Summary Search",
                      sidebarPanel(
                        selectizeInput("program2", "Program:", program_list, ""),
                        selectizeInput("study2", "Study:", NULL, ""),
                        selectizeInput("cell_type2", "Cell Type:", NULL, ""),
                        selectizeInput("gene2", "Gene:", NULL, "")
                      ),
                      mainPanel(
                        DTOutput("table")
                      )
             )
  )
)

server <- function(input, output, session) {
  observe({
    x <- as.character(input$program)
    
    updateSelectizeInput(session, "study",
                         choices = list.files(path=paste0('results_target_credentialing', "/", x)))
  })
  observe({
    x <- as.character(input$program)
    y <- as.character(input$study)
    
    updateSelectizeInput(session, "cell_type", 
                         choices = list.files(path=gsub(" ", "", paste( '~/cloud-data/cloud-pipeline-tim-tri-culture-storage/Drew_Macklin/results_target_credentialing', "/", x, "/", y))))
  })
  observe({
    x <- as.character(input$program)
    y <- as.character(input$study)
    z <- as.character(input$cell_type)
    
    updateSelectizeInput(session, "gene", 
                         choices = list.files(path=gsub(" ", "", paste( '~/cloud-data/cloud-pipeline-tim-tri-culture-storage/Drew_Macklin/results_target_credentialing', "/", x, "/", y, "/", z))))
  })
  observe({
    x <- as.character(input$program)
    y <- as.character(input$study)
    z <- as.character(input$cell_type)
    a <- as.character(input$gene)
    
    output$plot <- renderImage({
      list(src = gsub(" ", "", paste( '~/cloud-data/cloud-pipeline-tim-tri-culture-storage/Drew_Macklin/results_target_credentialing', "/", x, "/", y, "/", z, "/", a)),
           alt = "Alternate Text", width = 400, height = 400)
    }, deleteFile = FALSE)
  })
  
  observe({
    x <- as.character(input$program2)
    updateSelectizeInput(session, "study2",
                         choices = list.files(path=gsub(" ", "", paste( '~/cloud-data/cloud-pipeline-tim-tri-culture-storage/Drew_Macklin/results_target_credentialing', "/", x))))
  })
  observe({
    x <- as.character(input$program2)
    y <- as.character(input$study2)

      updateSelectizeInput(session, "cell_type2",
                           choices = list.files(path=gsub(" ", "", paste( '~/cloud-data/cloud-pipeline-tim-tri-culture-storage/Drew_Macklin/results_target_credentialing', "/", x, "/", y))))
  })
  observe({
    x <- as.character(input$program2)
    y <- as.character(input$study2)
    z <- as.character(input$cell_type2)

    updateSelectizeInput(session, "gene2",
                           choices = list.files(path=gsub(" ", "", paste( '~/cloud-data/cloud-pipeline-tim-tri-culture-storage/Drew_Macklin/results_target_credentialing', "/", x, "/", y, "/", z))))
  })
   observe({
    x <- as.character(input$program2)
    y <- as.character(input$study2)
    z <- as.character(input$cell_type2)
    a <- as.character(input$gene2)
  
    if(grepl(".csv", y)){
      output$table <- renderDT({as.data.frame(read.csv(paste0('results_target_credentialing/', x, "/", y)))
      })
    }else if(grepl(".csv", z)){
      output$table <- renderDT({as.data.frame(read.csv(paste0('results_target_credentialing/', x, "/", y, "/", z)))
      })
    }
    else{
      output$table <- renderDT({as.data.frame(read.csv(paste0('results_target_credentialing/', x, "/", y, "/", z, "/", a)))
      })
    }
  })
}

shinyApp(ui, server)