# Biome-shiny 0.9 - UI

library(shiny)
library(shinydashboard)
library(shinyBS)
library(microbiome)
library(speedyseq)
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

####### FUNCTIONS #######

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
ui <- dashboardPage( skin = "red",

  dashboardHeader(title = "Biome-Shiny v0.9"),
  dashboardSidebar(
    sidebarMenu(
      br(),
      paste0("Data Upload and Pre-Processing"),
      br(),
      menuItem( "Data Upload  (1/3)", tabName = "intro")
      ),
    #   sidebarMenu(
    # ),
    sidebarMenuOutput("dynamicMenu"),
    br(),
    h5(
      "Made with ",
      img(src = "shiny.png", height = "35"),
      "by ",
      img(src = "biodata.png", height = "25")
    )
  ),

  dashboardBody(
    tags$head(tags$style(
      HTML('
      /* body */
      .content-wrapper, .right-side { background-color: #ffffff; }' 
      ))),
    tabItems(
    #Introduction tab#
    tabItem(
      tabName = "intro",
      h2("Data Upload"),
      br(),
      box( solidHeader = TRUE,
           title = "Upload menu",
        width = "4",
        radioButtons(
          "datasetChoice",
          "Data Type",
          c("Biom file", 
            #"Phyloseq files", 
            "Sample dataset"),
          selected = "Biom file"
        ),
        conditionalPanel( condition = "input.datasetChoice == 'Biom file'",
                          radioButtons("datasetType", "Type of biom file:", c(".biom file with sample variables",".biom file with .csv metadata file",".biom file without .csv metadata file")),
                          fileInput(
                            "dataset",
                            "Dataset:",
                            multiple = FALSE,
                            accept = c(".biom"), placeholder="Phyloseq .biom files"
                          ),
                          conditionalPanel(
                            condition = "input.datasetType == '.biom file without .csv metadata file'",
                            checkboxInput("samplesAreColumns","OTU Table: Samples are columns", value = TRUE)
                          ),
                          conditionalPanel(
                            condition = "input.datasetType == '.biom file with .csv metadata file'",
                            fileInput("datasetMetadata", ".csv metadata file (sample variables):",
                                      multiple = FALSE,
                                      accept = c(".csv"), placeholder=".csv files"
                            )
                          ),
                          conditionalPanel( condition = "input.datasetChoice == 'Biom file'",
                            selectInput("taxParseFunction",label = "Taxonomy parsing function", choices = c("Default","QIIME","Greengenes"), selected = "Default")      
                          )
        ),
        
        # conditionalPanel(
        #   condition = "input.datasetChoice == 'Phyloseq files'",
        #   strong("*Required"),
        #   fileInput(
        #     "phyloseqOTUTable",
        #     "*OTU table file:",
        #     multiple = FALSE, 
        #     accept = c("text/csv"),
        #     placeholder = "File containing OTU table"
        #   ),
        #   fileInput(
        #     "phyloseqTaxTable",
        #     "*Taxonomy table file:",
        #     multiple = FALSE, 
        #     accept = c("text/csv"),
        #     placeholder = "File containing tax table"
        #   ),
        #   fileInput(
        #     "phyloseqMetadataTable",
        #     "Metadata table file:",
        #     multiple = FALSE, 
        #     accept = c("text/csv"),
        #     placeholder = "File containing metadata table"
        #   ),
        #   checkboxInput("samplesAreColumnsPhyloseq","OTU Table: Samples are columns", value = TRUE)
        # ),
        
        conditionalPanel(
          condition = "input.datasetChoice == 'Sample dataset'",
          selectInput(
            "datasetSample",
            "Choose a sample dataset:",
            choices = c("dietswap", "atlas1006")
          )
        ),
        actionButton("datasetUpdate", "Update Dataset"),
        bsTooltip("datasetUpdate", "Updates metadata rows when changing datasets", "bottom", options = list(container = "body"))
      )
      #tableOutput("uploadDataTable")
    ),

    tabItem(tabName="dataprocessing",
            fluidPage(
              column( width = 6,
              tabBox(
                title = "Data Filtering and Transformation", width = 10,
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
                ))
              ),
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
                            box( solidHeader = TRUE,
                              title = "Settings",
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
                            box(title = "Taxonomy", width = "2", collapsible = TRUE,
                                selectInput(
                                  "v4",
                                  "Taxonomic rank:",
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
        # A phyloseq richness plot
        tabPanel(title = "Richness Plot",
                 tabsetPanel(
                   tabPanel( title = "Plot Settings",
                    column( width = 6,
                        column( width = 5, 
                            box(
                               width = NULL,
                               title = "Variables",
                               # X (the metadata) and Y (the diversity measure)
                               selectInput(
                                 "x2",
                                 "Samples:",
                                 choices = colnames("datasetMetadata")
                               ),
                                selectInput("x3", "Point color:", choices = colnames("datasetMetadata"), selected = "subject"),
                                checkboxInput("transparentRichness", "Transparent background", value = TRUE),
                                checkboxInput("richnessPlotGridWrap", "Sort samples by metadata variable", value = FALSE),
                                conditionalPanel(condition = "input.richnessPlotGridWrap == 1",
                                                selectInput("x", "Sample sorting metadata:", choices = colnames("datasetMetadata"), selected = "nationality")
                               ),
                               radioButtons("richnessChoices", "Diversity measure:",
                                            choices = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"), selected = "Shannon")))
                        
                            )
                          ),
                        # column( width = 4,
                        #      box( width = NULL,
                        #        radioButtons("richnessChoices", "Diversity measure:",
                        #           choices = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"), selected = "Shannon")))
                        # )
                   #),
                   tabPanel( title = "Plot",
                             div(style = 'overflow-x: scroll', plotlyOutput("richnessPlot", height = "100%"))))
        ),
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
                 DT::dataTableOutput("metaTable"), downloadButton("downloadMetaTable"))
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
                              width = "4",
                              collapsible = TRUE,
                              collapsed = FALSE,
                              selectInput(
                                "xb", "Metadata:", choices = colnames("datasetMetadata"), selected = "bmi_group"
                              ),
                              selectInput(
                                "ordinate.method",
                                "Ordination method:",
                                choices = c("DCA", "CCA", "RDA", "CAP", "DPCoA", "NMDS", "MDS", "PCoA"),
                                selected = "NMDS"
                              ),
                              selectInput(
                                "ordinate.distance",
                                "Distance:",
                                choices = c("bray", "jaccard", "unifrac"),
                                selected = "bray"
                              ),
                              checkboxInput("dataTransformation", strong("Dataset:  Transform data with Hellinger method"), value = TRUE, width = "80%"),
                              sliderInput(
                                "geom.size",
                                "Point size:",
                                min = 1,
                                max = 10,
                                step = 0.5,
                                value = "2"
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
                                 "ordinate.method.taxa",
                                 "Ordination method:",
                                 choices = c("DCA", "CCA", "RDA", "CAP", "DPCoA", "NMDS", "MDS", "PCoA"),
                                 selected = "NMDS"
                               ),
                               selectInput(
                                 "ordinate.distance.taxa",
                                 "Distance:",
                                 choices = c("bray", "jaccard", "unifrac"),
                                 selected = "bray"
                               ),
                               sliderInput(
                                 "geom.size.taxa",
                                 "Point size:",
                                 min = 1,
                                 max = 10,
                                 step = 0.5,
                                 value = "2"
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
