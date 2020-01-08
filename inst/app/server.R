# Biome-shiny 0.9 - Server

library(shiny)
library(shinydashboard)
library(shinyBS)
library(microbiome)
library(phyloseq)
library(rmarkdown)
library(DT)
library(ggplot2)
library(plotly)
library(heatmaply)
#library(ComplexHeatmap)
library(knitr)
library(plyr)
library(dplyr)
library(ggpubr)
library(hrbrthemes)
library(reshape2)
library(vegan)
library(biomformat)
library(ggplotify)
library(RColorBrewer)

#Plot_ordered_bar function | Created by pjames1 @ https://github.com/pjames1
plot_ordered_bar<-function (physeq, x = "Sample",
                            y = "Abundance",
                            fill = NULL,
                            leg_size = 0.5,
                            title = NULL) {
  require(ggplot2)
  require(phyloseq)
  require(plyr)
  require(grid)
  bb <- psmelt(physeq)


  samp_names <- aggregate(bb$Abundance, by=list(bb$Sample), FUN=sum)[,1]
  .e <- environment()
  bb[,fill]<- factor(bb[,fill], rev(sort(unique(bb[,fill])))) #fill to genus


  bb<- bb[order(bb[,fill]),] # genus to fill
  p = ggplot(bb, aes_string(x = x, y = y,
                            fill = fill),
             environment = .e, ordered = FALSE)


  p = p +geom_bar(stat = "identity",
                  position = "stack",
                  color = "black")

  p = p + theme(axis.text.x = element_text(angle = -90, hjust = 0))

  p = p + guides(fill = guide_legend(override.aes = list(colour = NULL), reverse=TRUE)) +
    theme(legend.key = element_rect(colour = "black"))

  p = p + theme(legend.key.size = unit(leg_size, "cm"))


  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}

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

# Functions to dynamically generate chunks for the final report
tidy_function_body <- function(fun) {
  paste(tidy_source(text = as.character(body(fun))[-1])$text.tidy, collapse="\n")
}

make_chunk_from_function_body <- function(fun, chunk.name="", chunk.options=list()) {
  opts <- paste(paste(names(chunk.options), chunk.options, sep="="), collapse=", ")
  header <- paste0("```{r ", chunk.name, " ", chunk.options, "}")
  paste(header, tidy_function_body(fun), "```", sep="\n")
}

report.source <- reactive({
  req(sessionData$import.params(),
      sessionData$filter.params())

  report <- readLines("sc_report_base.Rmd")

  insert.function <- function(report, tag, fun, chunk.name = "", chunk.options = list()) {
    w <- which(report == tag)
    report[w] <- make_chunk_from_function_body(fun, chunk.name = chunk.name, chunk.options = chunk.options)

    return(report)
  }

  # Import
  report <- insert.function(report, "<!-- import.fun -->", sessionData$import.fun(), chunk.name = "import")

  # Filter
  report <- insert.function(report, "<!-- filter.fun -->", sessionData$filter.fun(), chunk.name = "filter")


  return(report)
})

#Function to collect all tables and images

collect_content <- function(){
    orca(communityPlotParams(), file = "community_plot.png")
}


#New Microbiome update messed up the formatting on the Phyloseq summary.
summarize_phyloseq_mod <- function(x){
  {
    ave <- minR <- maxR <- tR <- aR <- mR <- sR <- sR1 <- sR2 <- svar <- NULL
    sam_var <- zno <- comp <- NULL
    ave <- sum(sample_sums(x))/nsamples(x)
    comp <- length(which(colSums(abundances(x)) > 1))
    if (comp == 0) {
      a <- paste0("Compositional = YES")
    }
    else {
      a <- paste0("Compositional = NO")
    }
    minR <- paste0("Min. number of reads = ", min(sample_sums(x)))
    maxR <- paste0("Max. number of reads = ", max(sample_sums(x)))
    tR <- paste0("Total number of reads = ", sum(sample_sums(x)))
    aR <- paste0("Average number of reads = ", ave)
    mR <- paste0("Median number of reads = ", median(sample_sums(x)))
    if (any(taxa_sums(x) <= 1) == TRUE) {
      sR <- paste0("Any OTU sum to 1 or less? ", "YES")
    }
    else {
      sR <- paste0("Any OTU sum to 1 or less? ", "NO")
    }
    zno <- paste0("Sparsity = ", length(which(abundances(x) == 
                                                0))/length(abundances(x)))
    sR1 <- paste0("Number of singletons = ", length(taxa_sums(x)[taxa_sums(x) <= 
                                                                   1]))
    sR2 <- paste0("Percent of OTUs that are singletons (i.e. exactly one read detected across all samples): ", 
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
data("peerj32")
peerj32 <- peerj32$phyloseq

# Server
server <- function(input, output, session) {
  datasetChoice <- reactive({
    if (input$datasetChoice == "Use sample dataset") {
      switch(
        input$datasetSample,
        "dietswap" = dietswap,
        "atlas1006" = atlas1006,
        "peerj32" = peerj32
      )
    } else {
      if(input$datasetType == ".biom file including sample variables") { #Simple .biom upload with sample_data() already set
          req(input$dataset)
          tryCatch({
            datapath <- input$dataset$datapath
            biomfile <- import_biom(datapath)
            return(biomfile)
          }, error = function(e){
            simpleError("Error importing the .biom file.")
          })
      }
      if(input$datasetType == ".biom file with .csv metadata file"){ #Loads a .csv along with the .biom
        # req(input$dataset2)
        #req(input$datapathMetadata)
        
        #tryCatch({
        datapath <- input$dataset2$datapath
        a <- import_biom(datapath)
        #}, error = function(e){
        simpleError("Error importing the .biom file.")
        #})
        
        #tryCatch({
        datapathMetadata <- input$datasetMetadata$datapath
        b <- as.data.frame(read.csv(datapathMetadata, skipNul = TRUE))
        rownames(b) <- b[, 1]
        c <- sample_data(b)
        #}, error = function(e){
        simpleError("Error importing the .csv metadata file.")
        #})
        
        # tryCatch({
        biomfile <- merge_phyloseq(a,c)
        return(biomfile)
        #}, error = function(e){
        # simpleError("Error in merging .biom file with .csv metadata file.")
        #})
      }

      if(input$datasetType == ".biom file without .csv metadata file"){ #Loads a .biom file and generates sample metadata
          req(input$dataset3)
          tryCatch({
            datapath <- input$dataset3$datapath
            a <- import_biom(datapath)
          }, error = function(e){
            simpleError("Error importing .biom file")
          })
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
          })
          tryCatch({
            biomfile <- merge_phyloseq(a, b)
            return(biomfile)
          }, error = function(e){
            simpleError("Error merging sample variables with .biom file")
          })
        }
      }
    }
  )

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
      updateCheckboxGroupInput(session, "subsetTaxaByRankTaxList",
                               choices = levels(data.frame(tax_table(datasetChoice()))[[input$subsetTaxaByRank]]),
                               selected = levels(data.frame(tax_table(datasetChoice()))[[input$subsetTaxaByRank]]),
                               inline = TRUE
      )
    }, error = function(e) {
      simpleError(e)
    })
  })
  ## Update subsetSamples Checkbox Group ##
  observeEvent(input$datasetUpdate, {
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
  
  #Chceck all samples
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
  
  ## Table generation functions ##
  prevalenceAbsolute <- reactive({
    a <- as.data.frame(prevalence(compositionalInput(), detection = input$detectionPrevalence2/100, sort = TRUE, count = TRUE))
    names(a) <- c("Prevalence (counts)")
    return(a)
  })
  prevalenceRelative <- reactive({
    a <- as.data.frame(prevalence(compositionalInput(), detection = input$detectionPrevalence2/100, sort = TRUE))
    names(a) <- c("Prevalence (relative)")
    return(a)
  })
  
  ## Function to apply filters to the dataset ##
  filterData <- reactive({
    physeq <- datasetChoice()
    # Subset data by taxonomic rank - commented out for now since I'm having issues implementing it
    if (input$subsetTaxaByRankCheck == TRUE){
      oldMA <- tax_table(physeq)
      oldDF <- data.frame(oldMA)
      newMA <- prune_taxa(oldDF[[input$subsetTaxaByRank]] %in% input$subsetTaxaByRankTaxList, oldMA)
      # newDF <- subset(oldDF, oldDF[[input$subsetTaxaByRank]] ==  input$subsetTaxaByRankTaxList )
      # newMA <- as(newDF, "matrix")
      # if (inherits(physeq, "taxonomyTable")) {
      #   return(tax_table(newMA))
      # }
      # else {
      tax_table(physeq) <- tax_table(newMA)
      # }
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
    }
    #physeq <- filter_taxa(physeq, prune = TRUE, flist = filterfun(kOverA(A = input$detectionPrevalence2, k = 1)) )
    physeq <- core(physeq, detection = 0, prevalence = input$prevalencePrevalence )
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
        b <- heatmaply(otu_table(datasetInput()),
                       key.title = "Abundance", plot_method = "ggplot", height = 800, width = 1200,
                       heatmap_layers = theme(
                         panel.background = element_rect(fill = "transparent"),
                         plot.background = element_rect(fill = "transparent"),
                         legend.background = element_rect(fill = "transparent")
                       )
        )      }
    } else {
      if (ncol(otu_table(datasetInput())) > 1000){
        simpleError("A maximum of 1000 OTUs are permitted. Please filter the dataset and try again.")
      } else {
        b <- heatmaply(otu_table(datasetInput()),
                       key.title = "Abundance", plot_method = "ggplot",
                       heatmap_layers = theme(
                         panel.background = element_rect(fill = "transparent"),
                         plot.background = element_rect(fill = "transparent"),
                         legend.background = element_rect(fill = "transparent")
                    )
            )
      }
    }
    return(b)
  })
  output$coreHeatmap <- renderPlotly({
    ggplotly(coreHeatmapParams(), height = 800, width = 1280)
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
    if(input$communityPlotFacetWrap == FALSE){
      compositionplot <- plot_ordered_bar(datasetInput(), x=input$z1, y="Abundance", fill=input$v4, title=paste0("Abundance by ", input$v4, " in ", input$z1))  + geom_bar(stat="identity") + theme_pubr(base_size = 10, margin = TRUE, legend = "right", x.text.angle = 90) + rremove("xlab") + rremove("ylab")
    } else {
      compositionplot <- plot_ordered_bar(datasetInput(), x=input$z1, y="Abundance", fill=input$v4, title=paste0("Abundance by ", input$v4, " in ", input$z1))  + geom_bar(stat="identity") + theme_pubr(base_size = 10, margin = TRUE, legend = "right", x.text.angle = 90) + facet_grid(paste('~',input$z2), scales = "free", space = "free") + rremove("xlab") + rremove("ylab")
    }
    if(input$transparentCommunityPlot == TRUE){
      compositionplot <- compositionplot +
        theme(panel.background = element_rect(fill = "transparent", colour = NA), plot.background = element_rect(fill = "transparent", colour = NA), legend.background = element_rect(fill = "transparent", colour = NA), legend.box.background = element_rect(fill = "transparent", colour = NA))
    }
    p <- ggplotly(compositionplot, height = 500, width = plot_width(datasetInput())) %>% layout(xaxis = list(title = input$z1, automargin = TRUE), yaxis = list(title = "Abundance", automargin = TRUE))
    return(p)
  })
  output$communityPlot <- renderPlotly({
    communityPlotParams()
  })

  communityPlotGenusParams <- reactive({
    if(input$communityPlotFacetWrap == FALSE){
      compositionplot <- plot_ordered_bar(compositionalInput(), x=input$z1, fill=input$v4, title=paste0("Relative abundance by ", input$v4, " in ", input$z1))  + geom_bar(stat="identity") +
        guides(fill = guide_legend(ncol = 1)) +
        scale_y_percent() +
        theme_pubr(base_size = 10, margin = TRUE, legend = "right", x.text.angle = 90) + rremove("xlab") + rremove("ylab")
    } else {
      compositionplot <- plot_ordered_bar(compositionalInput(), x="Sample",  fill=input$v4, title=paste0("Relative abundance by ", input$v4, " in ", input$z1))  + geom_bar(stat="identity") +
        guides(fill = guide_legend(ncol = 1)) +
        scale_y_percent() + facet_grid(paste('~',input$z2),scales = "free", space = "free") + rremove("xlab") + rremove("ylab")
    }
    if(input$transparentCommunityPlot == TRUE){
      compositionplot <- compositionplot +
        theme(panel.background = element_rect(fill = "transparent", colour = NA), plot.background = element_rect(fill = "transparent", colour = NA), legend.background = element_rect(fill = "transparent", colour = NA), legend.box.background = element_rect(fill = "transparent", colour = NA))
    }
    p <- ggplotly(compositionplot, height = 500, width = plot_width(datasetInput())) %>% layout(xaxis = list(title = "Sample", automargin = TRUE), yaxis = list(title = "Abundance", automargin = TRUE))
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
    as.character(summarize_phyloseq_mod(datasetInput()))
  })
  
  output$summary <- renderPrint({
    summaryParams()
  })
  
  sampleVarsParams <- reactive({
    list_sample_variables(datasetInput())
  })
  
  output$sampleVars <- renderPrint({
    as.character(sampleVarsParams())  
  })
  
  ## Alpha Diversity ##

  #Abundance and Evenness tables#

  evennessParams <- reactive({
    a <- evenness(datasetInput())
    colnames(a) <- c("Camargo","Pielou","Simpson","Smith and Wilson's Evar","Bulla")
    a <- datatable(a, options = list(scrollX = TRUE))
  })

  output$evennessTable <- renderDataTable({
    evennessParams()
  })

  output$downloadEvenness <- downloadHandler(
    filename = function() {
      paste("evenness", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(evenness(datasetInput()), file, row.names = TRUE)
    }
  )

  absoluteAbundanceParams <- reactive({
    datatable(abundances(datasetInput()), options = list(scrollX = TRUE))
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
    datatable(abundances(datasetInput(), transform = "compositional"), options = list(scrollX = TRUE))
  })
  output$relativeAbundanceTable <- renderDataTable( server = FALSE, {
    relativeAbundanceParams()
  })
  output$downloadRelativeAbundance <- downloadHandler(
    filename = function() {
      paste("relativeAbundance", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(abundances(datasetInput(), transform = "compositional"), file, row.names = TRUE)
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
    datatable(sample_data(datasetInput()), options = list(scrollX = TRUE))
  })
  
  #Diversity measures table
  diversityMeasuresTableParams <- reactive({
    a <- estimate_richness(datasetInput(), measures = c("Observed","Chao1","Ace","Shannon","Simpson","InvSimpson","Fisher"))
    colnames(a) <- c("Observed species", "Chao1", "Standard error (Chao1)", "ACE", "Standard error (ACE)", "Shannon's diversity index","Simpson's diversity index","Inverse Simpson","Fisher's alpha")
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
    print(p)
  })

  output$richnessPlot <- renderPlotly({
    richnessPlotParams()
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
    microbiome::transform(datasetInput(), "compositional")
  })

  ordinateData <- reactive({
    ordinate(
      compositionalInput(),
      method = input$ordinate.method,
      distance = input$ordinate.distance
    )
  })

  ordinateDataSplit <- reactive({
    ordinate(
      compositionalInput(),
      method = input$ordinate.method2,
      distance = input$ordinate.distance2
    )
  })

  ordinateDataTaxa <- reactive({
    ordinate(
      compositionalInput(),
      method = input$ordinate.method3,
      distance = input$ordinate.distance3
    )
  })

  ordinatePlotParams <- reactive({
    if (ncol(sample_data(datasetInput())) > 1){
      p <- phyloseq::plot_ordination(datasetInput(), ordinateData(), color = input$xb, label = input$yb ) + geom_point(size = input$geom.size) + theme_pubr(base_size = 10, margin = TRUE, legend = "right")
    } else {
      a <- datasetInput()
      sample_data(a)[,2] <- sample_data(a)[,1]
      p <- phyloseq::plot_ordination(a, ordinateData(), color = input$xb, label = input$yb ) + geom_point(size = input$geom.size) + theme_pubr(base_size = 10, margin = TRUE, legend = "right")
    }
    if(input$transparentOrdinatePlot){
      p <- p +
        theme(panel.background = element_rect(fill = "transparent", colour = NA), plot.background = element_rect(fill = "transparent", colour = NA), legend.background = element_rect(fill = "transparent", colour = NA), legend.box.background = element_rect(fill = "transparent", colour = NA))
    }
    ggplotly(p, height = 500, width = 1050)
  })

  output$ordinatePlot <- renderPlotly({
    ordinatePlotParams()
  })

  # Split plot - not happy with how it looks - commented out
  # splitOrdParams <- reactive({
  #   if (ncol(sample_data(datasetInput())) > 1){
  #     splitOrdplot <-
  #       plot_ordination(
  #         datasetInput(),
  #         ordinateDataSplit(),
  #         type = "split",
  #         shape = input$xb,
  #         #color = input$yb,
  #         color = input$zbsplit
  #       ) + geom_point(size = input$geom.size2) + theme_pubr(base_size = 10, margin = TRUE, legend = "right")
  #   } else {
  #     a <- datasetInput()
  #     sample_data(a)[,2] <- sample_data(a)[,1]
  #     splitOrdplot <-
  #       plot_ordination(
  #         a,
  #         ordinateDataSplit(),
  #         type = "split",
  #         shape = input$xb,
  #         #color = input$yb,
  #         color = input$zbsplit
  #       ) + geom_point(size = input$geom.size2) + theme_pubr(base_size = 10, margin = TRUE, legend = "right")
  #   }
  #   if(input$transparentSplitOrd){
  #     splitOrdplot <- splitOrdplot +
  #       theme(panel.background = element_rect(fill = "transparent", colour = NA), plot.background = element_rect(fill = "transparent", colour = NA), legend.background = element_rect(fill = "transparent", colour = NA), legend.box.background = element_rect(fill = "transparent", colour = NA))
  #   }
  #   ggplotly(splitOrdplot, height = 500, width = 1050)
  # })
  #
  # output$splitOrd <- renderPlotly({
  #   splitOrdParams()
  # })

  taxaOrdParams <- reactive({
    taxaOrdplot <-
      plot_ordination(
        datasetInput(),
        ordinateDataTaxa(),
        type = "taxa",
        color = input$zb,
        label = input$xb
      ) + geom_point(size = input$geom.size3) + theme_pubr(base_size = 10, margin = TRUE, legend = "right")
    if(input$transparentTaxaOrd){
      taxaOrdplot <- taxaOrdplot +
        theme(panel.background = element_rect(fill = "transparent", colour = NA), plot.background = element_rect(fill = "transparent", colour = NA), legend.background = element_rect(fill = "transparent", colour = NA), legend.box.background = element_rect(fill = "transparent", colour = NA))
    }
    ggplotly(taxaOrdplot, height = 500, width = 1050)
  })

  output$taxaOrd <- renderPlotly({
    taxaOrdParams()
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
    otu <- abundances(compositionalInput())
    meta <- meta(compositionalInput())
    permnumber <- input$permanovaPermutationsP
    metadata <- input$permanovaColumnP
    m <- meta[[metadata]]
    a <- adonis(t(otu) ~ m,
           data = meta, permutations = permnumber, method = input$permanovaDistanceMethodP, parallel = getOption("mc.cores")
    )
    b <- as.data.frame(a$aov.tab)
    names(b) <- c(metadata, "Df", "Sum Sq", "Mean Sq", "F value", "P value")
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
    anova(betadisper(dist, meta[[metadata]]))
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
    print(p)
  })
  output$topFactorPlot <- renderPlot({
    topFactorPlotParams()
  })

  netPlotParams <- reactive({
    #if(input$permanovaPlotTypeNet == "samples"){
      n <- make_network(compositionalInput(), type = "samples", distance = input$permanovaDistanceMethodNet, max.dist = 0.2)
      p <- plot_network(n, compositionalInput(), type = "samples", shape = input$permanovaMetaShapeNet, color = input$permanovaMetadataNet, line_weight = 0.4)
    #}
    # if(input$permanovaPlotTypeNet == "taxa"){ 
    #   n <- make_network(compositionalInput(), type = "taxa", distance = input$permanovaDistanceMethodNet)
    #   p <- plot_network(n, compositionalInput(), type = "taxa", color= ntaxa(otu_table(datasetInput())))
    #   data <- subset_samples(compositionalInput(), !is.na(colnames(otu_table(compositionalInput()))))
    #   p <- plot_net(data, color = input$permanovaMetadataNet, shape = input$permanovaMetaShapeNet )
    # }
    if(input$transparentPermanova == TRUE){
      p <- p + theme(panel.background = element_rect(fill = "transparent", colour = NA), plot.background = element_rect(fill = "transparent", colour = NA), legend.background = element_rect(fill = "transparent", colour = NA), legend.box.background = element_rect(fill = "transparent", colour = NA))
    }
    ggplotly(p, height = 500, width = 1050)
  })
  output$netPlot <- renderPlotly({
    netPlotParams()
  })


  ### RESULTS ###


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
