# Biome-shiny 0.9 - UI

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

####### FUNCTIONS #######

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


####### END OF FUNCTIONS #######


# Load sample datasets #
data("dietswap")
data("atlas1006")

# UI
ui <- dashboardPage(
  dashboardHeader(title = "Biome-Shiny v0.9"),
  dashboardSidebar(
    sidebarMenu(
      br(),
      paste0("Data Upload and Pre-Processing"),
      br(),
      menuItem(
        "Data Upload  (1/3)",
        tabName = "intro"
      ),

      menuItem("Filtering (2/3)", tabName="dataprocessing"),
      menuItem("Phyloseq Summary  (3/3)", tabName="phyloseqsummary"),
      br(),
      paste0("Microbiome analysis"),
      menuItem("Core microbiota", tabName = "coremicrobiota"),
      menuItem("Community composition", tabName = "communitycomposition"),
      menuItem("Alpha diversity", tabName = "alphadiversity"),
      menuItem("Beta diversity", tabName = "betadiversity"),
      br(),
      paste0("Statistical analysis"),
      menuItem("PERMANOVA", tabName = "permanova"),
      br(),
      paste0("Outputs and Results"),
      menuItem("Results", tabName = "results")
    ),

    br(),
    h5(
      "Made with ",
      img(src = "shiny.png", height = "35"),
      "by ",
      img(src = "biodata.png", height = "25")
    )
  ),

  dashboardBody(
    tabItems(
    #Introduction tab#
    tabItem(
      tabName = "intro",
      h1("Data Upload"),
      br(),
      box(
        width = "3",
        radioButtons(
          "datasetChoice",
          "Dataset Upload",
          c("Upload dataset", "Use sample dataset"),
          selected = "Upload dataset"
        ),
        conditionalPanel( condition = "input.datasetChoice == 'Upload dataset'",
                          radioButtons("datasetType", "Select dataset characteristics:", c(".biom file including sample variables",".biom file with .csv metadata file",".biom file without .csv metadata file")),
                          conditionalPanel(
                            condition = "input.datasetType == '.biom file including sample variables'",
                            fileInput(
                              "dataset",
                              "Dataset:",
                              multiple = FALSE,
                              accept = c(".biom"), placeholder="Phyloseq .biom files"
                            )
                          ),
                          conditionalPanel(
                            condition = "input.datasetType == '.biom file without .csv metadata file'",
                            fileInput(
                              "dataset3",
                              "Dataset:",
                              multiple = FALSE,
                              accept = c(".biom"), placeholder="Phyloseq .biom files"
                            ),
                            checkboxInput("samplesAreColumns","OTU Table: Samples are columns", value = TRUE)
                          ),
                          conditionalPanel(
                            condition = "input.datasetType == '.biom file with .csv metadata file'",
                            fileInput(
                              "dataset2",
                              "Dataset:",
                              multiple = FALSE,
                              accept = c(".biom"), placeholder="Phyloseq .biom files"
                            ),
                            fileInput("datasetMetadata", ".csv metadata file (sample variables):",
                                      multiple = FALSE,
                                      accept = c(".csv"), placeholder=".csv files"
                            )
                          ),
                          conditionalPanel( condition = "input.datasetChoice == 'Upload dataset'",
                            selectInput("taxParseFunction",label = "Taxonomy parsing function", choices = c("Default","QIIME","Greengenes"), selected = "Default")      
                          )
        ),
        conditionalPanel(
          condition = "input.datasetChoice == 'Use sample dataset'",
          selectInput(
            "datasetSample",
            "Choose a sample dataset:",
            choices = c("dietswap", "atlas1006")
          )
        ),
        actionButton("datasetUpdate", "Update Dataset"),
        bsTooltip("datasetUpdate", "Click to update metadata rows when changing datasets.", "bottom", options = list(container = "body"))

      ),
      box(
        paste0(
          "Biome-shiny is a microbiome analysis pipeline developed with the Shiny library for R, and based, primarily, on the \"microbiome\" and \"phyloseq\" libraries for analysis.\n\n\n\nThe application takes a .biom file, generated by programs such as QIIME, as an input. If necessary, it is possible to upload a .csv file containing the dataset's sample data. Finally, if the user does not happen to have any sample data, the application can generate sample data out of the sample headers. For more information on the .biom file format please visit the following link: http://biom-format.org/ "
        )
      )
    ),

    tabItem(tabName="dataprocessing",
            fluidRow(
              tabBox(
                title = "Data filtering", width = 10,
                tabPanel( "Top taxa",
                          checkboxInput("pruneTaxaCheck", "Remove non-top taxa"),
                          numericInput("pruneTaxa", label = "Number of top taxa:", value = "10", min = "1")
                ),
                tabPanel( "Subset taxa",
                          checkboxInput("subsetTaxaByRankCheck", "Subset taxa by taxonomy rank"),
                          actionButton("subsetTaxaByRankUntickAll", "Uncheck all taxa"),
                          actionButton("subsetTaxaByRankTickAll", "Check all taxa"),
                          selectInput("subsetTaxaByRank", label = "Taxa rank:", choices = ""),
                          checkboxGroupInput("subsetTaxaByRankTaxList", inline = TRUE, label = "Taxa:", choices = "")
                ),
                tabPanel(
                  title = "Remove samples",
                  actionButton("subsetSamplesUntickAll", "Uncheck all samples"),
                  actionButton("subsetSamplesTickAll", "Check all samples"),
                  checkboxInput("subsetSamplesCheck", label = "Remove unchecked samples"),
                  checkboxGroupInput("subsetSamples", inline = TRUE, label = "Samples:", choices = "")
                )),
              box(
                title = "Activate changes", collapsible = TRUE, width = 2,
                checkboxInput("coreFilterDataset", "Set as active dataset", value = FALSE)
              )
            )),


    tabItem(tabName = "phyloseqsummary",
            h1("Phyloseq Summary"),
            verbatimTextOutput("summary"),
            paste0("Sample variables:"),
            verbatimTextOutput("sampleVars")
    ),


    # Core microbiota #
    tabItem(
      tabName = "coremicrobiota",
          fixedRow(
                   plotlyOutput("coreHeatmap", height = 800, width = 1000)
            )
        ),

    # Community composition #
    tabItem(
      tabName = "communitycomposition",
      tabsetPanel(
        #tabPanel(title = "Abundance in samples by taxa", #Absolute abundance/counts
         #        tabsetPanel(
                   tabPanel("Variables",
                            box(
                              width = "2",
                              collapsible = TRUE,
                              selectInput(
                                "z1",
                                "Sample variable:",
                                choices = colnames("datasetMetadata"),
                                selected = "sample"
                              ),
                              bsTooltip("z1", "The variable that corresponds to the samples in the data.", "right", options = list(container = "body")),
                              checkboxInput("transparentCommunityPlot", "Transparent background", value = TRUE),
                              checkboxInput("communityPlotFacetWrap", "Group samples by metadata variable (facet_wrap)", value = FALSE),
                              conditionalPanel(condition = "input.communityPlotFacetWrap == 1",
                                               selectInput(
                                                 "z2",
                                                 "Metadata:",
                                                 choices = colnames("datasetMetadata"),
                                                 selected = "nationality"
                                               )
                              )
                            ),
                            box("Taxonomy", width = "2", collapsible = TRUE,
                                selectInput(
                                  "v4",
                                  "Taxonomy rank:",
                                  # Tax rank to analyze
                                  choices = c("Phylum", "Class", "Order", "Family", "Genus"),
                                  selected = "Phylum"
                                )
                            )
                   ),
                   tabPanel("Absolute Abundance Plot", div(style = 'overflow-x: scroll', plotlyOutput("communityPlot", height = "100%"))),
                   tabPanel("Relative Abundance Plot", div(style = 'overflow-x: scroll', plotlyOutput("communityPlotGenus", height = "100%")))
        #         )
        #),
        #tabPanel(title = "Taxonomy Prevalence Plot",
        #         plotlyOutput("communityPrevalence", height = "100%")
        #)
      )
    ),

    # Alpha Diversity #
    tabItem(
      tabName = "alphadiversity",
      br(),
      tabsetPanel(
        type = "tabs",
        id = "tabsetpanel",

        tabPanel(title = "Evenness Table",
                 dataTableOutput("evennessTable"), downloadButton("downloadEvenness")),
        tabPanel(title = "Abundance Table",
                 tabsetPanel(
                   tabPanel( title = "Abundance (Counts)",
                             dataTableOutput("absoluteAbundanceTable"), downloadButton("downloadAbundance")
                   ),
                   tabPanel( title = "Abundance (Relative)",
                             dataTableOutput("relativeAbundanceTable"), downloadButton("downloadRelativeAbundance")
                   )
                 )
        ),
        tabPanel(title = "Richness Table",
                 DT::dataTableOutput("diversityMeasuresTable"), downloadButton("downloadDiversityMeasuresTable")),
        #  Metadata table
        tabPanel(title = "Metadata Table",
                 DT::dataTableOutput("metaTable"), downloadButton("downloadMetaTable")),

        # A phyloseq richness plot
        tabPanel(title = "Richness Plot",
                 tabsetPanel(
                   tabPanel( title = "Variables",
                             box(
                               width = "2",
                               title = "Variables",
                               # X (the metadata) and Y (the diversity measure)
                               selectInput(
                                 "x2",
                                 "Samples:",
                                 choices = colnames("datasetMetadata")
                               ),
                               selectInput("x3", "Point color metadata:", choices = colnames("datasetMetadata"), selected = "subject"),
                               checkboxInput("transparentRichness", "Transparent background", value = TRUE),
                               checkboxInput("richnessPlotGridWrap", "Sort samples by metadata variable", value = FALSE),
                               conditionalPanel(condition = "input.richnessPlotGridWrap == 1",
                                                selectInput("x", "Sample sorting metadata:", choices = colnames("datasetMetadata"), selected = "nationality")
                               )
                             ),
                             box(radioButtons("richnessChoices", "Choose diversity measure:" ,choices = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"), selected = "Shannon"))),
                   tabPanel( title = "Plot",
                             div(style = 'overflow-x: scroll', plotlyOutput("richnessPlot", height = "100%")))))
      )
    ),

    #Beta diversity#
    tabItem(
      tabName = "betadiversity",
      tabsetPanel(
        tabPanel(title = "Ordination Plot",
                 tabsetPanel(
                   tabPanel(title = "Variables",
                            box(
                              title = "Variables",
                              width = "2",
                              collapsible = TRUE,
                              collapsed = FALSE,
                              selectInput(
                                "xb", "Metadata:", choices = colnames("datasetMetadata"), selected = "bmi_group"
                              ),
                              selectInput(
                                "ordinate.method",
                                "Ordination method:",
                                choices = c("DCA", "CCA", "RDA", "CAP", "DPCoA", "NMDS", "MDS", "PCoA"),
                                selected = "CCA"
                              ),
                              selectInput(
                                "ordinate.distance",
                                "Distance:",
                                choices = c("bray", "jaccard", "unifrac"),
                                selected = "bray"
                              ),
                              sliderInput(
                                "geom.size",
                                "Point size:",
                                min = 1,
                                max = 10,
                                step = 0.5,
                                value = "3"
                              ),
                              checkboxInput("transparentOrdinatePlot", "Transparent background", value = TRUE)
                            )),
                   tabPanel(title = "Plot",
                            plotlyOutput("ordinatePlot"),
                            textOutput("ordinatePrint"))
                 )),

        tabPanel(title = "Taxa Plot",
                 tabsetPanel(
                   tabPanel( title = "Variables",
                             box(
                               title = "Variables",
                               width = "2",
                               collapsible = TRUE,
                               collapsed = FALSE,
                               selectInput(
                                 "zb",
                                 "Taxonomy rank:",
                                 choices = c("Phylum", "Class", "Order", "Family", "Genus")
                               ),
                               selectInput(
                                 "ordinate.method3",
                                 "Ordination method:",
                                 choices = c("DCA", "CCA", "RDA", "CAP", "DPCoA", "NMDS", "MDS", "PCoA"),
                                 selected = "CCA"
                               ),
                               selectInput(
                                 "ordinate.distance3",
                                 "Distance:",
                                 choices = c("bray", "jaccard", "unifrac"),
                                 selected = "bray"
                               ),
                               sliderInput(
                                 "geom.size3",
                                 "Point size:",
                                 min = 1,
                                 max = 10,
                                 step = 0.5,
                                 value = "3"
                               ),
                               checkboxInput("transparentTaxaOrd", "Transparent background", value = TRUE)
                             )),
                   tabPanel(title = "Plot", plotlyOutput("taxaOrd"))
                 )
        )
      )),
    tabItem(
      tabName = "permanova",
      tabsetPanel(
        tabPanel( title = "P-Value",
                  tabsetPanel(
                    tabPanel(title = "Variables",
                             box( title = "Variables", width= "2", collapsible = TRUE,
                                  selectInput("permanovaDistanceMethodP","Dissimilarity index:", choices = c("manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao", "mahalanobis"), selected = "bray"),
                                  selectInput("permanovaColumnP","Sample variable:", choices = colnames("datasetMetadata")),
                                  numericInput("permanovaPermutationsP", "Number of permutations:", min = 1, step = 1, value = 99),
                                  bsTooltip("permanovaPermutationsP", "The more permutations, the slower the result.", "bottom", options = list(container = "body"))
                             )
                    ),
                    tabPanel(title = "Data Tables",
                             h2("P-Value table"),
                             dataTableOutput("pValue"), downloadButton("downloadPValue"),
                             h2("Homogeniety table"),
                             dataTableOutput("homogeniety"), downloadButton("downloadHomogeniety")
                    )
                  )
        ),
        tabPanel ( title = "Top Factors",
                   tabsetPanel(
                     tabPanel(title = "Variables",
                              box( title = "Variables", width= "2", collapsible = TRUE,
                                   selectInput("permanovaDistanceMethodFac","Dissimilarity index:", choices = c("manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao", "mahalanobis"), selected = "bray"),
                                   selectInput("permanovaColumnFac","Sample variable:", choices = colnames("datasetMetadata")),
                                   numericInput("permanovaPermutationsFac", "Number of permutations:", min = 1, step = 1, value = 99)#,
                              )
                     ),
                     tabPanel(title = "Plot",
                              plotOutput("topFactorPlot")
                     )
                   )
        ),
        
        tabPanel ( title = "Network Plot",
                   tabsetPanel(
                     tabPanel(title = "Variables",
                              box( title = "Variables", width= "2", collapsible = TRUE,
                                   selectInput("permanovaDistanceMethodNet","Distance method:", choices = c("bray","jaccard","unifrac"), selected = "bray"),
                                   checkboxInput("transparentPermanova", "Transparent background", value = TRUE),
                                      selectInput("permanovaMetadataNet", "Sample variable to cluster data samples:", c("Update")),
                                      selectInput("permanovaMetaShapeNet", "Sample variable to set different point shapes:", c("Update"))
                                  )
                              ),
                        tabPanel(title ="Plot",
                              plotlyOutput("netPlot")
                        )

                     )
                   )
        )
      ),
    tabItem( tabName = "results",
             fluidRow(
               box( radioButtons('format', 'Download', c('HTML'), inline = TRUE, selected = 'HTML'),
                    downloadButton('downloadReportAlpha', label = "Download report")
               )
             )
          )
    )
  )
)
