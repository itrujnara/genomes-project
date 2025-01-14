#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(dplyr)
library(tidyr)
library(vroom)
library(data.table)
library(ggplot2)
library(plotly)
library(umap)
library(Rtsne)

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Jellyfish k-mers"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            selectInput("k",
                        "k:",
                        list("5" = 5, "6" = 6)
            ),
            checkboxInput("normalize",
                          "Normalize genome size?",
                          TRUE),
            varSelectInput("colorcol",
                           "Color by:",
                           data.frame()),
            sliderInput("tsne_perplexity",
                        "tSNE perplexity",
                        min = 1,
                        max = 100,
                        step = 1,
                        value = 30),
            sliderInput("tsne_steps",
                        "tSNE steps",
                        min = 100,
                        max = 5000,
                        step = 100,
                        value = 2000),
            sliderInput("umap_neighbors",
                        "UMAP neighbors",
                        min = 1,
                        max = 30,
                        step = 1,
                        value = 15),
            sliderInput("umap_mindist",
                        "UMAP min distance",
                        min = 0.01,
                        max = 1.,
                        step = 0.01,
                        value = 0.1),
            sliderInput("umap_steps",
                        "UMAP steps",
                        min = 100,
                        max = 2000,
                        step = 100,
                        value = 1000)
        ),

        # Show a plot of the generated distribution
        mainPanel(
           textOutput("text"),
           h2("PCA"),
           plotlyOutput("pcaPlot"),
           h2("tSNE"),
           plotlyOutput("tsnePlot"),
           h2("UMAP"),
           plotlyOutput("umapPlot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    kmers <- reactive({
      folder <- "static/"
      filename <- paste0("all_", input$k, "mers.tsv")
      datapath <- paste0(folder, filename)
      vroom(datapath)
    })
    meta_sample <- vroom("static/p10k_sample_metadata.csv")
    meta_genome <- vroom("static/p10k_genome_metadata.csv")
    
    metadata <- meta_genome %>% left_join(meta_sample, by = c("SampleID" = "Sample ID"))
    
    kmers_sel <- reactive({
      kmers() %>% select(-c("kmer"))
    })
    
    kmers_T <- reactive({
      table1 <- kmers_sel() %>% transpose()
      
      table2 <- kmers_sel() %>% apply(2, function(row) {
        row / sum(row) * 10^6
      }) %>% data.frame() %>% transpose()
      
      table <- list(table1, table2)[[input$normalize + 1]]
      colnames(table) <- kmers()$kmer
      rownames(table) <- colnames(kmers_sel())
      
      table$id <- rownames(table)
      table$assembly <- sapply(strsplit(table$id, "_"), tail, n = 1)
      
      table
    })
    
    kmers_clean <- reactive({
      kmers_T() %>% select(-c("id", "assembly"))
    })
    
    pca <- reactive({
      prcomp(na.omit(kmers_clean()), center = TRUE, scale = TRUE)
    })
    
    pctable <- reactive({
      table <- bind_cols(na.omit(kmers_T())$assembly, pca()$x)
      colnames(table)[1] <- "AssemblyID"
      table <- table %>% left_join(metadata, by = "AssemblyID")
      table
    })
    
    umap_basic <- reactive({
      umap(na.omit(kmers_clean()),
           n_neighbors = input$umap_neighbors,
           min_dist = input$umap_mindist,
           n_epochs = input$umap_steps
      )
    })
    
    umap_table <- reactive({
      table <- bind_cols(na.omit(kmers_T())$assembly, umap_basic()$layout)
      colnames(table) <- c("AssemblyID", "x", "y")
      table <- table %>% left_join(metadata, by = "AssemblyID")
      table
    })
    
    tsne <- reactive ({
      Rtsne(na.omit(kmers_T()), dims = 2,
            perplexity = input$tsne_perplexity, max_iter = input$tsne_steps)
    })
    
    tsne_table <- reactive({
      table <- bind_cols(na.omit(kmers_T())$assembly, tsne()$Y[,c(1,2)])
      colnames(table) <- c("AssemblyID", "x", "y")
      table <- table %>% left_join(metadata, by = "AssemblyID")
      table
    })
    
    updateVarSelectInput(
      session,
      "colorcol",
      data = metadata,
      selected = "AnnotLevel"
    )

    output$pcaPlot <- renderPlotly({
        pca_plot <- plot_ly()
        
        pca_plot <- pca_plot %>%
          add_trace(
            type = "scatter",
            mode = "markers",
            x = pctable()$PC1,
            y = pctable()$PC2,
            color = pull(pctable()[,input$colorcol]),
            text = pctable()$Species.x,
            hovertemplate = paste('<b>%{text}</b>',
                                  '<br>X: %{x}<br>',
                                  'Y: %{y}')
          )
        
        imp <- summary(pca())$importance[2,]
        
        pca_plot <- pca_plot %>% layout(
          xaxis = list(title = paste0("PC1 (", imp[1]*100, "%)")),
          yaxis = list(title = paste0("PC2 (", imp[2]*100, "%)")),
          legend = list(title = list(text = as.character(input$colorcol)))
        )
        
        pca_plot
    })
    
    output$tsnePlot <- renderPlotly({
      tsne_plot <- plot_ly()
      
      tsne_plot <- tsne_plot %>%
        add_trace(
          type = "scatter",
          mode = "markers",
          x = tsne_table()$x,
          y = tsne_table()$y,
          color = pull(tsne_table()[,input$colorcol]),
          text = tsne_table()$Species.x,
          hovertemplate = paste('<b>%{text}</b>',
                                '<br>X: %{x}<br>',
                                'Y: %{y}')
        )
      
      tsne_plot <- tsne_plot %>% layout(
        xaxis = list(title = "tSNE 1"),
        yaxis = list(title = "tSNE 2"),
        legend = list(title = list(text = as.character(input$colorcol)))
      )
      
      tsne_plot
    })
    
    output$umapPlot <- renderPlotly({
      umap_plot <- plot_ly()
      
      umap_plot <- umap_plot %>%
        add_trace(
          type = "scatter",
          mode = "markers",
          x = umap_table()$x,
          y = umap_table()$y,
          color = pull(umap_table()[,input$colorcol]),
          text = umap_table()$Species.x,
          hovertemplate = paste('<b>%{text}</b>',
                                '<br>X: %{x}<br>',
                                'Y: %{y}')
        )
      
      umap_plot <- umap_plot %>% layout(
        xaxis = list(title = "UMAP1"),
        yaxis = list(title = "UMAP2"),
        legend = list(title = list(text = as.character(input$colorcol)))
      )
      
      umap_plot
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
