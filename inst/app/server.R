# Biome-shiny 0.9 - Server

library(shiny)
library(shinydashboard)
library(shinyBS)
library(microbiome)
if("speedyseq" %in% rownames(installed.packages()) == TRUE){
  library(speedyseq)
} else {
  paste0("Biome-Shiny can use the 'speedyseq' library to speed up some phyloseq functions. If you'd like, you can install it from GitHub: https://github.com/mikemc/speedyseq")
}
library(rmarkdown)
library(DT)
library(ggplot2)
library(plotly)
library(heatmaply)
library(knitr)
library(plyr)
library(dplyr)
library(ggpubr)
library(vegan)
library(hrbrthemes)


#Function to dynamically set plot width (and height) for plots
plot_width <- function(data, mult = 12, min.width = 1060, otu.or.tax = "otu"){
  if(mult <= 0){
    print("Error: Variable 'mult' requires a value higher than 0")
    return(NULL)
  }
  if(min.width <= 0){
    print("Error: Variable 'min.width' requires a value higher than 0")
    return(NULL)
  }
  if(otu.or.tax == "otu"){
    width <- ncol(otu_table(data))*mult
    if(width <= min.width){ #Value of width needs to be higher than minimum width, default 1060px
      width <- min.width
      return(width)
    } else {
      return(width)
    }
  }
  if(otu.or.tax == "tax"){
    width <- nrow(tax_table(data))*mult
    if(width <= min.width){
      width <- min.width
      return(width)
    } else {
      return(width)
    }
  }
}

# Modified the summarize package 
summarize_phyloseq_mod <- function(x){
  {
    ave <- minR <- maxR <- tR <- aR <- mR <- sR <- sR1 <- sR2 <- svar <- NULL
    sam_var <- zno <- comp <- NULL
    ave <- sum(sample_sums(x))/nsamples(x)
    comp <- length(which(colSums(abundances(x)) > 1))
    if (comp == 0) {
      a <- paste0("Compositional? YES")
    }
    else {
      a <- paste0("Compositional? NO")
    }
    minR <- paste0("Min. number of reads = ", min(sample_sums(x)))
    maxR <- paste0("Max. number of reads = ", max(sample_sums(x)))
    tR <- paste0("Total number of reads = ", sum(sample_sums(x)))
    aR <- paste0("Average number of reads = ", ave)
    mR <- paste0("Median number of reads = ", median(sample_sums(x)))
    if (any(taxa_sums(x) == 1) == TRUE) {
      sR <- paste0("Singletons? ", "YES")
    }
    else {
      sR <- paste0("Singletons? ", "NO")
    }
    zno <- paste0("Sparsity = ", length(which(abundances(x) == 
                                                0))/length(abundances(x)))
    sR1 <- paste0("Number of singletons = ", length(taxa_sums(x)[taxa_sums(x) == 
                                                                   1]))
    sR2 <- paste0("Percent of OTUs that are singletons = ", 
                  mean(taxa_sums(x) == 1) * 100)
    svar <- paste0("Number of sample variables: ", ncol(meta(x)))
    list(a,minR, maxR, tR, aR, mR, zno, sR, sR1, sR2, svar)
  }
}
#Function to fix the formatting on the sample variables
list_sample_variables <- function(x){
  a<-colnames(sample_data(x))
  as.list(a)
}

# Load sample datasets #
data("dietswap")
data("atlas1006")

# Server
server <- function(input, output, session) {
  
  has_data_initialized <- reactiveValues()
  has_data_initialized$data_initialized <- FALSE
  has_data_initialized$filters_initialized <-FALSE
  
  datasetChoice <- reactive({
    # Taxa Parse Function switch
    taxparse <- switch(
      input$taxParseFunction,
      "Default" = parse_taxonomy_default,
      "QIIME" = parse_taxonomy_qiime,
      "Greengenes" = parse_taxonomy_greengenes
    )
    
    # Sample dataset selector
    if (input$datasetChoice == "Sample dataset") {
      withProgress(message = 'Loading sample dataset...', style = "notification", value = 0.5, {
      has_data_initialized$data_initialized = FALSE
      sampleDataset <- switch(
        input$datasetSample,
        "dietswap" = dietswap,
        "atlas1006" = atlas1006
        )
        has_data_initialized$data_initialized = TRUE
        return(sampleDataset)
      })
      return(sampleDataset)
    }
    
    # Biom file upload
    if(input$datasetChoice == "Biom file") {
      req(input$dataset)
      withProgress(message = 'Loading biom file...', style = "notification", value = 0.5, {
      tryCatch({
        has_data_initialized$data_initialized = FALSE
        datapath <- input$dataset$datapath
        a <- import_biom(datapath, parseFunction = taxparse)
      }, error = function(e){
        simpleError("Error importing the .biom file.")
      })})
      if(input$datasetType == ".biom file with sample variables") {
        has_data_initialized$data_initialized = TRUE
        return(a)
      }
      if(input$datasetType == ".biom file with .csv metadata file"){
        withProgress(message = 'Loading metadata...', style = "notification", value = 0.5, {
        datapathMetadata <- input$datasetMetadata$datapath
        b <- as.data.frame(read.csv(datapathMetadata, skipNul = TRUE))
        setProgress(value = 0.7, message = "Merging metadata and biom file...")
        rownames(b) <- b[, 1]
        c <- sample_data(b)
        setProgress(value = 0.9, message = "Merging metadata and biom file...")
        biomfile <- merge_phyloseq(a,c)
        setProgress(value = 1, message = "Loading complete!")
        })
        has_data_initialized$data_initialized = TRUE
        return(biomfile)
      }
      if(input$datasetType == ".biom file without .csv metadata file"){
        withProgress(message = 'Generating metadata...', style = "notification", value = 0.5, {
        tryCatch({
          if(input$samplesAreColumns == TRUE){
            samples.out <- colnames(otu_table(a))
          }
          if(input$samplesAreColumns == FALSE){
            samples.out <- rownames(otu_table(a))
            otu_table(a) <- phyloseq::t(otu_table(a))
          }
          subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
          samdf <- data.frame(Subject=subject)
          rownames(samdf) <- samples.out
          b <- sample_data(samdf)
        }, error = function(e){
          simpleError("Error generating sample variables")
        })})
        tryCatch({
          biomfile <- merge_phyloseq(a, b)
          has_data_initialized$data_initialized = TRUE
          return(biomfile)
        }, error = function(e){
          simpleError("Error merging sample variables with .biom file")
        })
      }
    }
    
    # Phyloseq file upload
    if(input$datasetChoice == "Phyloseq files"){
      #tryCatch({
        otu <- input$phyloseqOTUTable$datapath
        tax <- input$phyloseqTaxTable$datapath
        otu.tab <- as.data.frame(read.csv(otu, skipNul = TRUE))
        tax.tab <- as.data.frame(read.csv(tax,skipNul = TRUE))
        pfile <- phyloseq(otu_table(otu.tab, taxa_are_rows = input$samplesAreColumnsPhyloseq), tax_table(tax.tab)) #Might get an error with taxa are rows
        if(input$phyloseqMetadataTable == is.null()){
          if(input$samplesAreColumnsPhyloseq == TRUE){
            samples.out <- colnames(otu_table(pfile))
          } else {
            samples.out <- rownames(otu_table(pfile))
          }
          subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
          samdf <- data.frame(Subject=subject)
          rownames(samdf) <- samples.out
          b <- sample_data(samdf)
          pfile <- merge_phyloseq(pfile, b)
          return(pfile)
        } else {
          meta.tab.path <- as.data.table(read.table(input$phyloseqMetadataTable$datapath), skipNul = TRUE)
          meta.tab <- sample_data(meta.tab.path)
          pfile <- merge(pfile, meta.tab)
          return(pfile)
        }
      #}, error = function(e){
       # simpleError()
      #})
    }

  })
  
  # uploadDataTableParams <- reactive(
  #   req(datasetInput),
  #   table <- as.data.frame(input$datasetInput$datapath),
  #   return(table)
  # )
  # 
  # output$uploadDataTable <- renderTable(
  #   uploadDataTableParams()
  # )
  # 
  # New DatasetInput function works as an intermediary that checks if the dataset has been altered
  datasetInput <- reactive({
    if(input$coreFilterDataset == TRUE ){ # Filters the dataset
      dataset <- filterData()
    }
    if(input$coreFilterDataset == FALSE ) { # Standard dataset input without filtering applied
      dataset <- datasetChoice()
    }
    return(dataset)
  })
  
  ## Dynamic Menu ##
  
  output$dynamicMenu <- renderMenu(
    
    if( has_data_initialized$data_initialized == FALSE){
      sidebarMenu(
        br(),
        strong(paste0("Waiting for dataset..."))
      )
    } else {
      sidebarMenu(
        menuItem("Filtering and Transformations (2/3)", tabName="dataprocessing"),
        menuItem("Phyloseq Summary  (3/3)", tabName="phyloseqsummary"),
        br(),
        paste("Microbiome analysis"),
        menuItem("Core microbiota", tabName = "coremicrobiota"),
        menuItem("Community composition", tabName = "communitycomposition"),
        menuItem("Alpha diversity", tabName = "alphadiversity"),
        menuItem("Beta diversity", tabName = "betadiversity"),
        br(),
        paste0("Statistical analysis"),
        menuItem("PERMANOVA", tabName = "permanova"),
        br(),
        paste("Outputs and Results"),
        menuItem("Results", tabName = "results"),
        br(),
        paste("Settings"),
        uiOutput("decimalSlider")
      )
    }
  )
  


  ## Dataset Filtering ##
  
  ## Populate SelectInput with taxonomic ranks ##
  observeEvent(input$datasetUpdate, {
    tryCatch({
      updateSelectInput(session, "subsetTaxaByRank",
                        choices = colnames(tax_table(datasetInput())))
    }, error = function(e) {
      simpleError(e)
    })
  }, ignoreNULL = FALSE)
  
  ## Update Checkbox Group based on the chosen taxonomic rank ##
  observeEvent(input$subsetTaxaByRank, {
    tryCatch({
      withProgress(message = 'Filtering and Transformations', detail = "Updating Taxonomic Ranks", style = "notification", min = 0, max = 1, value = 0.1, {
      updateCheckboxGroupInput(session, "subsetTaxaByRankTaxList",
                               choices = levels(data.frame(tax_table(datasetChoice()))[[input$subsetTaxaByRank]]),
                               selected = levels(data.frame(tax_table(datasetChoice()))[[input$subsetTaxaByRank]]),
                               inline = TRUE
      )})
    }, error = function(e) {
      simpleError(e)
    })
  })
  ## Update subsetSamples Checkbox Group ##
  observeEvent(input$datasetUpdate, {
    tryCatch({
      withProgress(message = 'Filtering and Transformations', detail = "Updating Samples", style = "notification", min = 0, max = 1, value = 0.1, {
      updateCheckboxGroupInput(session, "subsetSamples",
                               choices = colnames(otu_table(datasetChoice())),
                               selected = colnames(otu_table(datasetChoice())),
                               inline = TRUE
      )})
    }, error = function(e) {
      simpleError(e)
    })
  })
  
  #Check all samples
  observeEvent(input$subsetSamplesTickAll, {
    tryCatch({
      updateCheckboxGroupInput(session, "subsetSamples",
                               choices = colnames(otu_table(datasetChoice())),
                               selected = colnames(otu_table(datasetChoice())),
                               inline = TRUE
      )
    }, error = function(e) {
      simpleError(e)
    })
  })
  
  #Check all taxa
  observeEvent(input$subsetTaxaByRankTickAll, {
    tryCatch({
      updateCheckboxGroupInput(
        session, "subsetTaxaByRankTaxList",
        choices = levels(data.frame(tax_table(datasetChoice()))[[input$subsetTaxaByRank]]),
        selected = levels(data.frame(tax_table(datasetChoice()))[[input$subsetTaxaByRank]]),
        inline = TRUE
      )
    }, error = function(e) {
      simpleError(e)
    }
    )
  })
  
  
  #Uncheck all taxa
  observeEvent(input$subsetTaxaByRankUntickAll, {
    
  tryCatch({
        updateCheckboxGroupInput(
        session, "subsetTaxaByRankTaxList",
        choices = levels(data.frame(tax_table(datasetChoice()))[[input$subsetTaxaByRank]]),
        selected = NULL,
        inline = TRUE
      )
    }, error = function(e) {
      simpleError(e)
    }
    )
  })
  
  #Uncheck all samples
  observeEvent(input$subsetSamplesUntickAll, {
    tryCatch({
      updateCheckboxGroupInput(
        session, "subsetSamples",
        choices = colnames(otu_table(datasetChoice())),
        selected = NULL,
        inline = TRUE
      )
    }, error = function(e) {
      simpleError(e)
    }
    )
  })

  ## Function to apply filters to the dataset ##
  filterData <- reactive({
    withProgress(message = 'Applying filters to the dataset', detail = "Please wait...", style = "notification", min = 0, max = 1, value = 0.1, {
    physeq <- datasetChoice()
    # Rarefy data
    if(input$rarefactionCheck == TRUE) {
      physeqRare <- rarefy_even_depth(physeq, sample.size = min(sample_sums(physeq)),
                                      replace = input$rarefactionReplace,
                                      rngseed = input$rarefactionSeed )
      physeq <- physeqRare
    }
    # Subset data by taxonomic rank - commented out for now since I'm having issues implementing it
    if (input$subsetTaxaByRankCheck == TRUE){
      oldMA <- tax_table(physeq)
      oldDF <- data.frame(oldMA)
      newMA <- prune_taxa(oldDF[[input$subsetTaxaByRank]] %in% input$subsetTaxaByRankTaxList, oldMA)
      tax_table(physeq) <- tax_table(newMA)
    }
    #Filter out non-top taxa
    if (input$pruneTaxaCheck == TRUE){
      filterTaxa <- names(sort(taxa_sums(physeq), decreasing = TRUE)[1:input$pruneTaxa])
      physeq <- prune_taxa(filterTaxa, physeq)
    }
    #Filter out samples
    if (input$subsetSamplesCheck == TRUE){
      oldDF <- as(sample_data(physeq), "data.frame")
      newDF <- subset(oldDF, colnames(otu_table(physeq)) %in% input$subsetSamples)
      sample_data(physeq) <- sample_data(newDF)
    }})
    return(physeq)
  })

  output$corePhyloSummary <- renderPrint({ # Summary of corePhylo file
    summarize_phyloseq(filterData)
  })
  output$coreTaxa <- renderPrint({ # Reports the taxa in corePhylo
    taxa(filterData)
  })

  output$prevalenceAbsoluteOutput <- renderDT({
    datatable(prevalenceAbsolute())
  })

  output$downloadPrevalenceAbsolute <- downloadHandler(
    filename = function() {
      paste("PrevalenceAbsolute", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(prevalenceAbsolute(), file, row.names = TRUE)
    }
  )

  output$prevalenceRelativeOutput <- renderDT({
    datatable(prevalenceRelative())
  })

  output$downloadPrevalenceRelative <- downloadHandler(
    filename = function() {
      paste("PrevalenceRelative", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(prevalenceRelative(), file, row.names = TRUE)
    }
  )

  ## Core Microbiota ##

  coreHeatmapParams <- reactive({
    if ( input$samplesAreColumns == TRUE ) {
      if ( nrow(otu_table(datasetInput())) > 1000 ){
        simpleError("A maximum of 1000 OTUs are permitted. Please filter the dataset and try again.")
      } else {
        withProgress(message = 'Generating plot...', detail = "This may take a bit.", style = "notification", min = 0, max = 1, value = 0.1, {
          # for(i in 1:30){ 
          #     incProgress(1/30)
          #     Sys.sleep(1)
          #   }
        b <- heatmaply(otu_table(datasetInput()),
                       key.title = "Abundance", plot_method = "ggplot", height = "100%", width = "100%",
                       heatmap_layers = theme(
                         panel.background = element_rect(fill = "transparent"),
                         plot.background = element_rect(fill = "transparent"),
                         legend.background = element_rect(fill = "transparent")
                       )
        )})      }
    } else {
      if (ncol(otu_table(datasetInput())) > 1000){
        simpleError("A maximum of 1000 OTUs are permitted. Please filter the dataset and try again.")
      } else {
        withProgress(message = 'Generating plot...', detail = "This may take a bit.", style = "notification", min = 0, max = 1, value = 0.1, {
          b <- heatmaply(otu_table(datasetInput()),
                       key.title = "Abundance", plot_method = "ggplot",
                       heatmap_layers = theme(
                         panel.background = element_rect(fill = "transparent"),
                         plot.background = element_rect(fill = "transparent"),
                         legend.background = element_rect(fill = "transparent")
                    )
            )

        })
      }
    incProgress(amount = 0.1, message = "Plot is being generated.", detail = "This may take a bit.")
    return(b)
    }
  })
  output$coreHeatmap <- renderPlotly({
    ggplotly(coreHeatmapParams()) %>% config(editable = T, autosizable = T)
  })

  ## Community Composition ##


  # Updating SelectInputs when database changes #
  observeEvent(input$datasetUpdate, {
    tryCatch({
    updateSelectInput(session, "z1",
                      choices = colnames(meta(datasetInput())))
    updateSelectInput(session, "z2",
                      choices = colnames(meta(datasetInput())))
    updateSelectInput(session, "z3",
                      choices = colnames(meta(datasetInput())))
    updateSelectInput(session, "v4",
                      choices = colnames(tax_table(datasetInput())))
    updateSelectInput(session, "z1Average",
                      choices = colnames(meta(datasetInput())))
    updateSelectInput(session, "v4Plot",
                      choices = colnames(tax_table(datasetInput())))
    }, error = function(e) {
      simpleError(e)
    })
  }, ignoreNULL = FALSE)

  #Update metadata value selectInputs for CC Analysis
  observeEvent(input$z1, {
    updateSelectInput(session, "v1",
                      choices = (sample_data(datasetInput())[[input$z1]]))
  })
  observeEvent(input$z2, {
    updateSelectInput(session, "v2",
                      choices = (sample_data(datasetInput())[[input$z2]]))
  })
  observeEvent(input$z3, {
    updateSelectInput(session, "v3",
                      choices = (sample_data(datasetInput())[[input$z3]]))
  })

  # Abundance of taxa in sample variable by taxa
  communityPlotParams <- reactive ({
    withProgress(message = 'Generating plot...', detail = "This may take a bit.", style = "notification", min = 0, max = 1, value = 0.1, {
    taxglom <- tax_glom(datasetInput(), taxrank=input$v4)
    compositionplot <- plot_bar(taxglom, x=input$z1, y="Abundance", fill=input$v4, title=paste0("Abundance by ", input$v4, " in ", input$z1))  + geom_bar(stat="identity") + theme_pubr(base_size = 10, margin = TRUE, legend = "right", x.text.angle = 90) + rremove("xlab") + rremove("ylab")
    if(input$communityPlotFacetWrap == TRUE){
    compositionplot <- compositionplot + facet_grid(paste('~',input$z2), scales = "free", space = "free")
    }
    if(input$transparentCommunityPlot == TRUE){
      compositionplot <- compositionplot +
        theme(panel.background = element_rect(fill = "transparent", colour = NA), plot.background = element_rect(fill = "transparent", colour = NA), legend.background = element_rect(fill = "transparent", colour = NA), legend.box.background = element_rect(fill = "transparent", colour = NA))
    }
    p <- ggplotly(compositionplot, height = 500, width = plot_width(datasetInput())) %>% layout(xaxis = list(title = input$z1, automargin = TRUE), yaxis = list(title = "Abundance", automargin = TRUE))
    })
    return(p)
  })
  output$communityPlot <- renderPlotly({
    communityPlotParams() %>% config(editable = T, autosizable = T)
  })

  communityPlotGenusParams <- reactive({
    withProgress(message = 'Generating plot...', detail = "This may take a bit.", style = "notification", min = 0, max = 1, value = 0.1, {
    taxglom <- tax_glom(microbiome::transform(datasetInput(), "compositional"), taxrank=input$v4)
    compositionplot <- plot_bar(taxglom, x=input$z1, fill=input$v4, title=paste0("Relative abundance by ", input$v4, " in ", input$z1))  + geom_bar(stat="identity") +
      guides(fill = guide_legend(ncol = 1)) +
      scale_y_percent() +
      theme_pubr(base_size = 10, margin = TRUE, legend = "right", x.text.angle = 90) + rremove("xlab") + rremove("ylab")
    
    if(input$communityPlotFacetWrap == TRUE){
      compositionplot <- compositionplot + facet_grid(paste('~',input$z2),scales = "free", space = "free")
    }
    if(input$transparentCommunityPlot == TRUE){
      compositionplot <- compositionplot +
        theme(panel.background = element_rect(fill = "transparent", colour = NA), plot.background = element_rect(fill = "transparent", colour = NA), legend.background = element_rect(fill = "transparent", colour = NA), legend.box.background = element_rect(fill = "transparent", colour = NA))
    }
    p <- ggplotly(compositionplot, height = 500, width = plot_width(datasetInput())) %>% layout(xaxis = list(title = input$z1, automargin = TRUE), yaxis = list(title = "Abundance", automargin = TRUE))
    })
    return(p)
  })
  output$communityPlotGenus <- renderPlotly({
    communityPlotGenusParams()
  })

  # Taxa prevalence plot
  communityPrevalenceParams <- reactive({
    prevplot <- plot_taxa_prevalence(compositionalInput(), input$v4) + theme_pubr(base_size = 10, margin = TRUE, legend = "right", x.text.angle = 90) + #If OTUs > 25 it fails
    rremove("xlab") + rremove("ylab")
    if(input$transparentCommunityPlot == TRUE){
      prevplot <- prevplot +
        theme(panel.background = element_rect(fill = "transparent", colour = NA), plot.background = element_rect(fill = "transparent", colour = NA), legend.background = element_rect(fill = "transparent", colour = NA), legend.box.background = element_rect(fill = "transparent", colour = NA))
    }
    p <- ggplotly(prevplot, height = 500, width = 1000) %>%  layout(xaxis = list(title = "Average count abundance (log scale)", automargin = TRUE), yaxis = list(title = "Taxa prevalence", automargin = TRUE))
    return(p)
  })

  output$communityPrevalence <- renderPlotly({
    communityPrevalenceParams() 
  })


  # Phyloseq Summary #
  summaryParams <- reactive({
    req(datasetInput())
    withProgress(message = 'Generating summary...', style = "notification", min = 0, max = 1, value = 0.1, {
    as.character(summarize_phyloseq_mod(datasetInput()))})
  })
  
  output$summary <- renderPrint({
      summaryParams()
  })
  
  sampleVarsParams <- reactive({
      list_sample_variables(datasetInput())
  })
  
  output$sampleVars <- renderPrint({
    withProgress(message = 'Printing sample variables...', style = "notification", value = 0.1, {
      as.character(sampleVarsParams())
    })
  })
  
  ## Alpha Diversity ##

  # Slider for adjusting number of decimal cases

  output$decimalSlider <- renderUI({
    sliderInput(
      "decimalCases",
      "Tables: Number of decimal cases:",
      min = 1,
      max = 15,
      step = 1,
      value = 5,
      width = "80%"
    )
  })
  
  output$standardizationPick <- renderUI({
    checkboxInput("dataTransformation", strong("Dataset: Standardize data with Hellinger method"), value = TRUE, width = "80%")
  })
  
  decimalCalc <- eventReactive(input$decimalCases, ignoreNULL = TRUE, {
    input$decimalCases
  })

  
  #Abundance and Evenness tables#
  
  evennessParams <- reactive({
    withProgress(message = 'Generating table...', style = "notification", min = 0, max = 1, value = 0.1, {
    a <- evenness(datasetInput()) %>% round(digits = decimalCalc())
    colnames(a) <- c("Camargo","Pielou","Simpson","Smith and Wilson's Evar","Bulla")
    })
    datatable(a, options = list(scrollX = TRUE))
  })
  output$evennessTable <- renderDataTable({
    evennessParams()
  })

  output$downloadEvenness <- downloadHandler(
    filename = function() {
      paste("evenness", ".csv", sep = "")
    },
    content = function(file) {
      a <- evenness(datasetInput()) %>% round(digits = input$decimalCases)
      write.csv(a, file, row.names = TRUE)
    }
  )

  absoluteAbundanceParams <- reactive({
    withProgress(message = 'Generating table...', style = "notification", min = 0, max = 1, value = 0.1, {
    a <- abundances(datasetInput())
    })
    datatable(a, options = list(scrollX = TRUE))
  })
  output$absoluteAbundanceTable <- renderDataTable( server = FALSE, {
    absoluteAbundanceParams()
  })

  output$downloadAbundance <- downloadHandler(
    filename = function() {
      paste("evenness", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(abundances(datasetInput()), file, row.names = TRUE)
    }
  )

  relativeAbundanceParams <- reactive({
    withProgress(message = 'Generating table...', style = "notification", min = 0, max = 1, value = 0.1, {
    a <- abundances(datasetInput(), transform = "compositional") %>% round(digits = input$decimalCases)
    })
    datatable(a, options = list(scrollX = TRUE))
  })
  output$relativeAbundanceTable <- renderDataTable( server = FALSE, {
    relativeAbundanceParams()
  })
  output$downloadRelativeAbundance <- downloadHandler(
    filename = function() {
      paste("relativeAbundance", ".csv", sep = "")
    },
    content = function(file) {
      a <- abundances(datasetInput(), transform = "compositional") %>% round(digits = input$decimalCases)
      write.csv(a, file, row.names = TRUE)
    }
  )

  # Updating SelectInputs when database changes #
  observeEvent(input$datasetUpdate, {
    tryCatch({
    updateSelectInput(session, "x",
                      choices = colnames(meta(datasetInput())))
    updateSelectInput(session, "x2", choices = colnames(meta(datasetInput())))
    updateSelectInput(session, "x3", choices = colnames(meta(datasetInput())))
    updateSelectInput(session, "y",
                      choices = colnames(alpha(datasetInput())))
    }, error = function(e){
      simpleError(e)
    })
  }, ignoreNULL = FALSE)


  # Metadata table#
  metaTableParams <- reactive({
    withProgress(message = 'Generating table...', style = "notification", min = 0, max = 1, value = 0.1, {})
    datatable(sample_data(datasetInput()), options = list(scrollX = TRUE))
  })
  
  #Diversity measures table
  diversityMeasuresTableParams <- reactive({
    withProgress(message = 'Generating table...', style = "notification", min = 0, max = 1, value = 0.1, {
    a <- estimate_richness(datasetInput(), measures = c("Observed","Chao1","Ace","Shannon","Simpson","InvSimpson","Fisher")) %>% round(digits = input$decimalCases)
    colnames(a) <- c("Observed species", "Chao1", "Standard error (Chao1)", "ACE", "Standard error (ACE)", "Shannon's diversity index","Simpson's diversity index","Inverse Simpson","Fisher's alpha")
    })
    datatable(a, options = list(scrollX = TRUE))
  })

  output$metaTable <- DT::renderDataTable({
    metaTableParams()
  })

  output$downloadMetaTable <- downloadHandler(
    filename = function() {
      paste("MetadataTable", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(meta(datasetInput()), file, row.names = TRUE)
    }
  )

  output$diversityMeasuresTable <- DT::renderDataTable({
    diversityMeasuresTableParams()
  })
  
  output$downloadDiversityMeasuresTable <- downloadHandler(
    filename = function() {
      paste("DiversityMeasuresTable", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(meta(datasetInput()), file, row.names = TRUE)
    }
  )
  
  # Alpha Diversity Richness Plot #
  richnessPlotParams <- reactive({
    withProgress(message = 'Generating plot...', style = "notification", min = 0, max = 1, value = 0.1, {
    if(input$richnessPlotGridWrap == FALSE){
      richnessplot <- plot_richness(
        datasetInput(),
        x = input$x2,
        measures = input$richnessChoices,
        color = input$x3
      ) + theme_pubr(base_size = 10, margin = TRUE, legend = "right", x.text.angle = 90)
    } else {
      richnessplot <- plot_richness(
        datasetInput(),
        x = input$x2,
        measures = input$richnessChoices,
        color = input$x3 ) +
        facet_grid(paste('~',input$x),scales = "free", space = "free") + theme_pubr(base_size = 10, margin = TRUE, legend = "right", x.text.angle = 90)
    }
    if(input$transparentRichness == TRUE){
      richnessplot <- richnessplot +
      theme(panel.background = element_rect(fill = "transparent", colour = NA), plot.background = element_rect(fill = "transparent", colour = NA), legend.background = element_rect(fill = "transparent", colour = NA), legend.box.background = element_rect(fill = "transparent", colour = NA))
    }
    richnessplot <- richnessplot + rremove("xlab") + rremove("ylab")
    p <- ggplotly(richnessplot, height = 500, width = plot_width(datasetInput())) %>% layout(xaxis = list(title = input$x2, automargin = TRUE), yaxis = list(title = paste("Alpha Diversity Measure (", input$richnessChoices , ")"), automargin = TRUE))
    })
    print(p)
  })

  output$richnessPlot <- renderPlotly({
    richnessPlotParams() %>% config(editable = T, autosizable = T)
  })

  ## Beta Diversity ##

  # Updating SelectInputs when dataset changes#
  observeEvent(input$datasetUpdate, {
    tryCatch({
    updateSelectInput(session, "xb",
                      choices = colnames(meta(datasetInput())))
    updateSelectInput(session, "xb2",
                      choices = colnames(meta(datasetInput())))
    updateSelectInput(session, "xb3",
                      choices = colnames(meta(datasetInput())))
    updateSelectInput(session, "yb",
                      choices = colnames(meta(datasetInput())))
    updateSelectInput(session, "zb",
                      choices = colnames(tax_table(datasetInput())))
    updateSelectInput(session, "zbsplit",
                      choices = colnames(tax_table(datasetInput())))
    }, error = function(e){
      simpleError(e)
    })
  }, ignoreNULL = FALSE)

  compositionalInput <- reactive({
     a <- microbiome::transform(datasetInput(), transform = "compositional", target = "OTU")
    return(a)
  })
  
  betaDiversityInput <- reactive({
    if(input$dataTransformation == TRUE){
      a <- microbiome::transform(datasetInput(), transform = "hellinger", target = "OTU")
    } else {
      a <- compositionalInput()
    }
    return(a)
  })

  ordinateData <- reactive({
    ordinate(
      betaDiversityInput(),
      method = input$ordinate.method,
      distance = input$ordinate.distance
    )
  })

  ordinateDataTaxa <- reactive({
    ordinate(
      betaDiversityInput(),
      method = input$ordinate.method.taxa,
      distance = input$ordinate.distance.taxa
    )
  })

  ordinatePlotParams <- reactive({
    withProgress(message = 'Generating plot...', style = "notification", value = 0.5 + 0.1, min = 0, max = 1, {
    if (ncol(sample_data(datasetInput())) > 1){
      p <- phyloseq::plot_ordination(datasetInput(), ordinateData(), color = input$xb, label = NULL ) + geom_point(size = input$geom.size) + theme_pubr(base_size = 10, margin = TRUE, legend = "right")
    } else {
      a <- datasetInput()
      sample_data(a)[,2] <- sample_data(a)[,1]
      p <- phyloseq::plot_ordination(a, ordinateData(), color = input$xb, label = NULL ) + geom_point(size = input$geom.size) + theme_pubr(base_size = 10, margin = TRUE, legend = "right")
    }
    if(input$transparentOrdinatePlot){
      p <- p +
        theme(panel.background = element_rect(fill = "transparent", colour = NA), plot.background = element_rect(fill = "transparent", colour = NA), legend.background = element_rect(fill = "transparent", colour = NA), legend.box.background = element_rect(fill = "transparent", colour = NA))
    }
    })
    ggplotly(p, height = 500, width = 1050)
  })

  output$ordinatePlot <- renderPlotly({
    ordinatePlotParams() %>% config(editable = T, autosizable = T)
  })

  taxaOrdParams <- reactive({
    withProgress(message = 'Generating plot...', style = "notification", value = 0.5, {
    taxaOrdplot <-
      plot_ordination(
        datasetInput(),
        ordinateDataTaxa(),
        type = "taxa",
        color = input$zb,
        label = input$xb
      ) + geom_point(size = input$geom.size.taxa) + theme_pubr(base_size = 10, margin = TRUE, legend = "right")
    if(input$transparentTaxaOrd){
      taxaOrdplot <- taxaOrdplot +
        theme(panel.background = element_rect(fill = "transparent", colour = NA), plot.background = element_rect(fill = "transparent", colour = NA), legend.background = element_rect(fill = "transparent", colour = NA), legend.box.background = element_rect(fill = "transparent", colour = NA))
    }
    })
    ggplotly(taxaOrdplot, height = 500, width = 1050)
  })

  output$taxaOrd <- renderPlotly({
    taxaOrdParams() %>% config(editable = T, autosizable = T)
  })

  ###########################
  ## Statistical analysis ###
  ###########################

  ## PERMANOVA ##
  #Update metadata column when dataset changes
  observeEvent(input$datasetUpdate, {
    tryCatch({
    updateSelectInput(session, "permanovaColumn",
                      choices = colnames(meta(datasetInput())))
    updateSelectInput(session, "permanovaColumnP",
                      choices = colnames(meta(datasetInput())))
    updateSelectInput(session, "permanovaColumnFac",
                      choices = colnames(meta(datasetInput())))
    updateSelectInput(session, "permanovaMetadataNet",
                      choices = colnames(meta(datasetInput())))
    updateSelectInput(session, "permanovaMetaShapeNet",
                    choices = colnames(meta(datasetInput())))
    }, error = function(e){
      simpleError(e)
    })
  }, ignoreNULL = FALSE)

  permanova <- reactive({
    withProgress(message = 'Calculating PERMANOVA...', style = "notification", value = 0.5, {
    otu <- abundances(compositionalInput())
    meta <- meta(compositionalInput())
    permnumber <- input$permanovaPermutationsP
    metadata <- input$permanovaColumnP
    m <- meta[[metadata]]
    a <- adonis(t(otu) ~ m,
           data = meta, permutations = permnumber, method = input$permanovaDistanceMethodP, parallel = getOption("mc.cores")
    )
    b <- as.data.frame(a$aov.tab) %>% round(digits = input$decimalCases)
    names(b) <- c(metadata, "Df", "Sum Sq", "Mean Sq", "F value", "P value")
    })
    print(b)
  })

  output$pValue <- renderDataTable({
    permanova()
  })

  output$downloadPValue <- downloadHandler(
    filename = function() {
      paste("P-ValueTable", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(permanova(), file, row.names = TRUE)
    }
  )

  homogenietyParams <- reactive({
    otu <- abundances(compositionalInput())
    meta <- meta(compositionalInput())
    dist <- vegdist(t(otu))
    metadata <- input$permanovaColumnP
    anova(betadisper(dist, meta[[metadata]])) %>% round(digits = input$decimalCases)
  })

  output$homogeniety <- renderDataTable({
    homogenietyParams()
  })

  output$downloadHomogeniety <- downloadHandler(
    filename = function() {
      paste("HomogenietyTable", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(homogenietyParams(), file, row.names = TRUE)
    }
  )

  topFactorPlotParams <- reactive({
    withProgress(message = 'Generating plot...', style = "notification", value = 0.5, {
    otu <- abundances(compositionalInput())
    meta <- meta(compositionalInput())
    permnumber <- input$permanovaPermutationsFac
    metadata <- input$permanovaColumnFac
    column <- meta[[metadata]]
    permanova <- adonis(t(otu) ~ column,
                        data = meta, permutations = permnumber, method = input$permanovaDistanceMethodFac
    )
    coef <- coefficients(permanova)["column1",]
    top.coef <- coef[rev(order(abs(coef)))[1:20]] #top 20 coefficients
    par(mar = c(3, 14, 2, 1))
     p <- barplot(sort(top.coef), horiz = T, las = 1, main = "Top taxa")
    })
    print(p)
  })
  output$topFactorPlot <- renderPlot({
    topFactorPlotParams()
  })

  netPlotParams <- reactive({
      withProgress(message = 'Generating metadata...', style = "notification", value = 0.5, {
      n <- make_network(compositionalInput(), type = "samples", distance = input$permanovaDistanceMethodNet, max.dist = 0.2)
      p <- plot_network(n, compositionalInput(), type = "samples", shape = input$permanovaMetaShapeNet, color = input$permanovaMetadataNet, line_weight = 0.4)

    if(input$transparentPermanova == TRUE){
      p <- p + theme(panel.background = element_rect(fill = "transparent", colour = NA), plot.background = element_rect(fill = "transparent", colour = NA), legend.background = element_rect(fill = "transparent", colour = NA), legend.box.background = element_rect(fill = "transparent", colour = NA))
    }
    })
    ggplotly(p, height = 500, width = 1050)
  })
  output$netPlot <- renderPlotly({
    netPlotParams() %>% config(editable = T, autosizable = T)
  })


  ### RESULTS ### # Needs a serious rewrite all around


  output$downloadReportAlpha <- downloadHandler(
    filename = function() {
      paste('report', sep = '.', switch(
        input$format, PDF = 'pdf', HTML = 'html'
      ))
    },
    content = function(file) {
      src <- normalizePath('final_report.Rmd')

      # temporarily switch to the temp dir, in case you do not have write
      # permission to the current working directory
      owd <- setwd(tempdir())
      on.exit(setwd(owd))
      file.copy(src, 'final_report.Rmd', overwrite = TRUE)

      out <- rmarkdown::render('final_report.Rmd',
                               switch(input$format,
                                      PDF = pdf_document(),
                                      HTML = html_document()
                               ))
      file.rename(out, file)
    }
  )
}
