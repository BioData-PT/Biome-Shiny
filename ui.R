# Biome-shiny 0.8 - UI

library(shiny)
library(shinydashboard)
library(shinyBS)
library(microbiome)
library(phyloseq)
library(rmarkdown)
library(DT)
library(ggplot2)
library(plotly)
library(knitr)
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

####### END OF FUNCTIONS #######


# Load sample datasets #
data("dietswap")
data("atlas1006")
data("peerj32")
peerj32 <- peerj32$phyloseq

# UI
ui <- dashboardPage(
  dashboardHeader(title = "Biome-shiny v0.8"),
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
      #menuItem("Density analysis", tabName = "landscaping"),
      #menuItem("DMM Clustering", tabName = "dirichlet"),
      br(),
      paste0("Statistical analysis"),
      menuItem("PERMANOVA", tabName = "permanova"),
      # menuItem("ANOSIM", tabName = "anosim"),
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
                          radioButtons("datasetType", "Select dataset characteristics:", c(".biom file including sample variables",".biom file with .csv metadata file",".biom file without .csv metadata file"))
        ),
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
        conditionalPanel(
          condition = "input.datasetChoice == 'Use sample dataset'",
          selectInput(
            "datasetSample",
            "Choose a sample dataset:",
            choices = c("dietswap", "atlas1006", "peerj32")
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
                     tabsetPanel(
                       tabPanel("Variables",
                                fluidRow(
                                box(
                                     title = "Fiter by abundance", collapsible = TRUE,
                                     #numericInput("detectionPrevalence2", "Minimum Abundance:", min = 0.00, max = 100, value = 0, step = 1),
                                     #bsTooltip("detectionPrevalence2", "Minimum abundance value of OTUs.", "right", options = list(container = "body")),
                                     numericInput("prevalencePrevalence","Prevalence [0-1]:", min = 0, max = 1, value = 0.0, step = 0.05),
                                     bsTooltip("prevalencePrevalence", "Ratio of OTU presence in samples. ", "right", options = list(container = "body")),
                                     checkboxInput("coreFilterDataset", "Set as active dataset", value = FALSE)
                                ),
                                box(
                                    title = "Prune and subset taxa",
                                    checkboxInput("pruneTaxaCheck", "Keep the top taxa"),
                                    conditionalPanel(
                                      condition = "input.pruneTaxaCheck == 1",
                                      numericInput("pruneTaxa", label = "Number of top taxa:", value = "10", min = "1")
                                    ),
                                    checkboxInput("subsetTaxaByRankCheck", "Subset taxa by taxonomy rank"),
                                    conditionalPanel(
                                      condition = "input.subsetTaxaByRankCheck == 1",
                                      selectInput("subsetTaxaByRank", label = "Taxa rank:", choices = ""),
                                      checkboxGroupInput("subsetTaxaByRankTaxList", label = "Taxa:", choices = "")
                                    )
                                ),
                                box(
                                    title = "Remove samples", collapsible = TRUE, collapsed = TRUE,
                                    checkboxInput("subsetSamplesCheck", label = "Remove unchecked samples"),
                                    checkboxGroupInput("subsetSamples", inline = TRUE, label = "Samples:", choices = "")
                                )
                       )),
                       tabPanel("Absolute prevalence", dataTableOutput("prevalenceAbsoluteOutput"), downloadButton("downloadPrevalenceAbsolute")),
                       tabPanel("Relative prevalence", dataTableOutput("prevalenceRelativeOutput"), downloadButton("downloadPrevalenceRelative"))#,
                       #tabPanel("Summary", verbatimTextOutput("corePhyloSummary")),
                       #tabPanel("Taxa", verbatimTextOutput("coreTaxa"))
                     )
    ),

    tabItem(tabName = "phyloseqsummary",
            h1("Phyloseq Summary"),
            verbatimTextOutput("summary")
    ),


    # Core microbiota #
    tabItem(
      tabName = "coremicrobiota",
          fixedRow(
              box( width = "2", collapsible = TRUE,
                   textInput("detectionMin", label = "Minimum detection threshold (Relative Abundance)", value = "0.0000001"),
                   bsTooltip("detectionMin", "Lowest detection value on the heatmap. Value must be higher than 0.", "left", options = list(container = "body")),
                   checkboxInput("transparentCoreHeatmap", "Transparent background", value = TRUE)
              ),
                   box(width = 10, div(style = 'overflow-x: scroll', plotlyOutput("coreHeatmap", width = "100%", height = "100%")))
            )
        ),

    # Community composition #
    tabItem(
      tabName = "communitycomposition",
      tabsetPanel(
        tabPanel(title = "Abundance in samples by taxa", #Absolute abundance/counts
                 tabsetPanel(
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
                 )
        ),
        tabPanel(title = "Taxonomy Prevalence Plot",
                 plotlyOutput("communityPrevalence", height = "100%")
        )
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

        # Alpha diversity measures with metadata
        tabPanel(title = "Metadata Table with diversity measures",
                 DT::dataTableOutput("view"), downloadButton("downloadMergedTable")),

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
                                "xb", "Sample variable:", choices = colnames("datasetMetadata"), selected = "bmi_group"
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

        # Split Ordination Plot tab - not happy with how it looks, commenting out for now
        # tabPanel(title = "Split Ordination Plot",
        #          tabsetPanel(
        #            tabPanel(title = "Variables",
        #                     box(
        #                       title = "Variables",
        #                       width = "2",
        #                       collapsible = TRUE,
        #                       collapsed = FALSE,
        #                       selectInput(
        #                         "xb2", "Sample variable:", choices = colnames("datasetMetadata"), selected = "bmi_group"
        #                       ),
        #
        #                       selectInput(
        #                         "zbsplit",
        #                         "Taxonomy rank:",
        #                         choices = c("Phylum", "Class", "Order", "Family", "Genus")
        #                       ),
        #
        #                       selectInput(
        #                         "ordinate.method2",
        #                         "Ordination method:",
        #                         choices = c("DCA", "CCA", "RDA", "CAP", "DPCoA", "NMDS", "MDS", "PCoA"),
        #                         selected = "CCA"
        #                       ),
        #
        #                       selectInput(
        #                         "ordinate.distance2",
        #                         "Distance:",
        #                         choices = c("bray", "jaccard", "unifrac"),
        #                         selected = "bray"
        #                       ),
        #
        #                       sliderInput(
        #                         "geom.size2",
        #                         "Point size:",
        #                         min = 1,
        #                         max = 10,
        #                         step = 0.5,
        #                         value = "3"
        #                       ),
        #                       checkboxInput("transparentSplitOrd", "Transparent background", value = TRUE)
        #                     )),
        #            tabPanel( title = "Plot", plotlyOutput("splitOrd")))),

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
      tabName = "landscaping",
      tabsetPanel(
        tabPanel("PCA",
                 plotlyOutput("landscapePCA")
        ),
        tabPanel("PCoA/MDS/NMDS",
                 tabsetPanel(
                   tabPanel("Variables",
                            box("Data Subsetting", width = "2", collapsible = TRUE,
                                numericInput("detectionLandscape", label="Input detection threshold for landscape analysis", min = 0, max = 100, step = 0.1, value = "0.1"),
                                numericInput("prevalenceLandscape", label="Input prevalence percentage of samples for landscape analysis", min = 0, max = 100, step = 0.1, value = "50")
                            ),
                            box("Data Subsetting", width = "2", collapsible = TRUE,
                                selectInput("metadataLandscape1", "Choose first metadata column to subset data:", choices=colnames("datasetMetadata")),
                                selectInput("metadataLandscapeValue1", "Choose first metadata value:", choices=sapply("datasetMetadata", levels)),
                                selectInput("metadataLandscape2", "Choose second metadata column to subset data:", choices=colnames("datasetMetadata")),
                                selectInput("metadataLandscapeValue2", "Choose second metadata value:", choices=sapply("datasetMetadata", levels))
                            ),
                            box(title = "Variables", width = "2", collapsible = TRUE,
                                selectInput("ordinateMethodLandscape", "Choose an ordination method:",
                                            choices=c("DCA", "CCA", "RDA", "CAP", "DPCoA", "NMDS", "MDS", "PCoA"), selected = "CCA"),
                                selectInput("ordinateDistanceLandscape", "Choose a distance method:",
                                            choices=c("bray","jaccard","unifrac"), selected = "unifrac")
                            )
                   ),
                   tabPanel("Plot", plotlyOutput("landscapeOrdination")
                            #,plotlyOutput("landscapeOrdinationSamplenames")
                   )
                 )
        ),
        tabPanel("t-SNE",
                 plotlyOutput("landscapeTSne")
        ),
        tabPanel("Abundance histogram",
                 tabsetPanel(
                   tabPanel("Taxa",
                            box(
                              selectInput("landscapeOTU", "Choose OTU:", choices = "")
                            )
                   ),
                   tabPanel("Absolute Abundance",
                            plotlyOutput("landscapeAbundanceHistAbs")
                   ),
                   tabPanel("Relative Abundance",
                            plotlyOutput("landscapeAbundanceHistRel")
                   )
                 )
        )
      )

    ),
    tabItem(
      tabName = "dirichlet",
      tabsetPanel(
        tabPanel("Variables",
                 box( title = "Variables", width = "2", collapsible = TRUE,
                      numericInput("detectionDMM", label="Input detection threshold for DMM", min = 0, max = 100, step = 0.1, value = "0.1"),
                      numericInput("prevalenceDMM", label="Input prevalence percentage of samples for DMM0", min = 0, max = 100, step = 0.1, value = "50"),
                      numericInput("maxModelsDMM", label="Input maximum number of community types for the DMM model.", value = "3")
                 )
        ),

        tabPanel("Model Verification Lineplot",
                 plotOutput("dmmModelCheck")
        ),

        tabPanel("Alpha and Theta Parameters",
                 dataTableOutput("dmmParameters"),
                 dataTableOutput("sampleAssignments")
        ),

        tabPanel("Taxa Contribution Per Community Type",
                 plotlyOutput("taxaContributionPerComponent")
        )
      )
    ),
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
                                   selectInput("permanovaDistanceMethodFac","Distance method:", choices = c("bray","jacard","unifrac"), selected = "bray"),
                                   selectInput("permanovaColumnFac","Sample variable:", choices = colnames("datasetMetadata")),
                                   numericInput("permanovaPermutationsFac", "Number of permutations:", min = 1, step = 1, value = 99)#,
                                   # checkboxInput("transparentPermanovaTopfactors", "Transparent background", value = TRUE)
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
                                   selectInput("permanovaPlotTypeNet", "Network Plot Type:", c("samples", "taxa"), selected = "samples"),
                                   selectInput("permanovaDistanceMethodNet","Distance method (not required by all ordination methods):", choices = c("bray","jacard","unifrac"), selected = "bray"),
                                   selectInput("permanovaMethodNet","Ordination method:",
                                               choices = c("DCA", "CCA", "RDA", "CAP", "DPCoA", "NMDS", "MDS", "PCoA"),
                                               selected = "CCA"),
                                   numericInput("permanovaPermutationsNet", "Number of permutations:", min = 1, step = 1, value = 99),
                                   checkboxInput("transparentPermanova", "Transparent background", value = TRUE),
                                   conditionalPanel(condition = "input.permanovaPlotTypeNet == 'samples'",
                                      selectInput("permanovaMetadataNet", "Sample variable to cluster data samples:", c("Update")),
                                      selectInput("permanovaMetaShapeNet", "Sample variable to set different point shapes:", c("Update"))
                                    )
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
             tabsetPanel(
               tabPanel("Core Microbiota",
                        box(
                          radioButtons("renderPrevTables", "Render Prevalence Tables", choices = c("Yes", "No")),
                          radioButtons("renderCoreTaxaSummary", "Render Core Taxa Summary", choices = c("Yes", "No")),
                          radioButtons("renderCoreHeatmap", "Render Core Taxa Heatmap", choices = c("Yes", "No"))
                        )
               ),
               tabPanel("Community Composition",
                        box(
                          radioButtons("renderSampleAbundancePlot", "Render Abundance in Samples Plot", choices = c("Yes", "No")),
                          radioButtons("renderRelativeAbundancePlot", "Render Relative Abundance in Samples Plot", choices = c("Yes", "No")),
                          radioButtons("renderPrevalencePlot", "Render Taxa Prevalence Plot", choices = c("Yes", "No"))
                        )
               ),
               tabPanel("Alpha Diversity",
                        box(
                          radioButtons("renderRichness", "Render Richness Plot", choices = c("Yes", "No")),
                          radioButtons("printTableAlpha", "Print Table Head", choices = c("Yes", "No"))
                        )
               ),
               tabPanel("Beta Diversity",
                        box(
                          radioButtons("renderOrdplot", "Render Ordination Plot", choices = c("Yes", "No")),
                          radioButtons("renderSplitOrdplot", "Render Split Ordination Plot", choices = c("Yes", "No")),
                          radioButtons("renderTaxaplot", "Render Taxa Plot", choices = c("Yes", "No"))
                        )
               ),

               tabPanel("PERMANOVA",
                        box(
                          radioButtons("renderDensityPlot", "Render density plot", choices = c("Yes","No")),
                          radioButtons("renderPValueTables", "Render P-Value tables", choices = c("Yes","No")),
                          radioButtons("renderFactorPlot", "Render top factors barplot", choices = c("Yes","No")),
                          radioButtons("renderNetworkMap", "Render network map", choices = c("Yes","No"))
                        )
               )#,
               # tabPanel("TEST",
               #          box(
               #            textInput("projectName", "Name your project", value = "Report1"),
               #            actionButton("fullDownloadTest", "Download Full Report")
               #          )
               #        )
             ),
             fluidRow(
               box( radioButtons('format', 'Document format (HTML only for now)', c('HTML'), inline = TRUE, selected = 'HTML'),
                    downloadButton('downloadReportAlpha', label = "Download report")
               )
             )
          )
    )
  )
)
