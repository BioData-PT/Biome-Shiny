## Biome-Shiny: A Shiny R app for microbiome visualization

### Introduction

Biome-Shiny is a graphical interface for visualizing microbiome data, primarily based around the "phyloseq" and "microbiome" libraries, with the "shiny" library being used to develop the interface. It provides a user-friendly interface to generate and visualize community composition and diversity through interactive plots.

### Installation

Biome-Shiny has some dependencies outside of R that need to be installed before downloading the package. In Ubuntu:

    sudo apt install libcairo2-dev

Biome-Shiny is an R package, and can installed, along with its dependencies using devtools.

    library("devtools")
    install_github("https://github.com/BioData-PT/Biome-Shiny")

### Launching the app

All you need to do to launch the app is open R and type in the R console:
        
        biomeshiny.package::launchApp()
    
Then type in your browser of choice:
    
        localhost:4800

    
### Features

Biome-Shiny has a variety of visualizations and tables that can be used to subset and visualize a dataset.

#### Uploading Data

Biome-Shiny makes use of phyloseq's .biom importing functions to upload data. It accepts a .biom file as input, and requires metadata (sample_data() in phyloseq format files).

    If the .biom file includes metadata, tick ".biom file including sample variables"
    
    If the metadata file is a separate .csv file, tick ".biom file with .csv metadata file"
    
    If the .biom file has no metadata nor a .csv metadata file, tick ".biom file without .csv metadata file". The app will attempt to generate metadata out of the OTU table's sample headers.
    
Then click "Update Dataset" to load the data.

Alternatively, you can use the sample data provided by microbiome for testing purposes, by clicking "Use sample dataset", followed by "Update Dataset"

#### Data Filtering

Biome-Shiny has several functions for subsetting and filtering datasets.

Once done with filtering, tick "Set as active dataset" to apply changes, and untick to return the dataset to default settings.

##### Subset by top taxa

In the "Top taxa" tab, the user can remove all but a number of taxa from the dataset, leaving only a small number of top taxa. To apply changes, the user must tick "Remove non-top taxa". Default value is top 10 taxa.

##### Remove taxa

The "Subset taxa" tab allows the user to remove all species or OTUs belonging to a certain taxa. The user must choose a taxonomic rank in the dropdown list, then untick all taxa that they wish to see removed. Once they're done, tick "Subset taxa by taxonomy rank" to apply changes.

##### Remove samples

The user can remove individual samples from the dataset in the "Remove samples" tab. Like the previous section, to remove samples just untick all samples that should be removed and tick "Removed unchecked samples".

##### Subset by prevalence

It's possible to subset the plot by a value of prevalence - the ratio of OTU/species presence in samples. For example, a prevalence of 0.5 will make it so that any species that doesn't appear in at least 50% of samples will be removed from the dataset. In the "Prevalence [0-1]:" box, input a number between 0 and 1. Inputting a value of 1 will likely cause an error, unless there is an OTU that's present in 100% of samples.

#### Phyloseq Summary

Displays a general summary of the entire dataset, based on the microbiome library's summarize_phyloseq() function. This includes smallest, largest and median reads, number of singletons, and whether the data is in compositional format or not.

#### Core Microbiota Heatmap

A heatmap that compares abundance values for each OTU on each sample.

#### Community Composition

Allows the user to view how abundant a given taxa is, generating abundance (absolute and relative) barplots which allow the samples to be sorted by metadata.

In the "Sample variable:" dropdown list, the user needs to pick the metadata variable which corresponds to the samples. If the uploaded dataset didn't have any metadata attached to it, the sample variable was automatically generated from the OTU table's sample headers.

By clicking "Group samples by metadata variable (facet_wrap)", the user can divide the plot bars, based on the levels of a given metadata variable (e.g. splitting samples between two different nationalities). The metadata which will be used to divide the samples can be picked in the "Metadata:" dropdown list.

The "Taxonomy Rank:" dropdown colors the bars based on the chosen taxonomic rank.

#### Alpha Diversity

Alpha diversity is the index that measures diversity within a sample.

Biome-Shiny displays several tables with diversity measures in the samples, such as an Evenness Table (measures how equal the OTU numbers are in each sample), Abundance Table (measures the number of each OTU in the table, whether as a simple count or in relation to the full number of species in the sample) and a table with assorted diversity measures (measure richness, among other things).

The Richness Plot can be used to visualize these diversity measures. The "Samples:" dropdown list, much like in Community Composition, asks for the metadata variable corresponding to the dataset's samples. "Point color metadata" will give each point on the plot (corresponding to a sample) a color based on the chosen metadata variable.

Finally, it's possible to group samples by a metadata sample in the same way as the Community Composition plots.


#### Beta Diversity

Beta diversity is an index that measures diversity between samples. Biome-Shiny measures Beta diversity by the samples' general composition. This is done by generating a distance matrix and representing it as a PCA plot, where each point corresponds to a sample.

The "Metadata:" dropdown list gives the points in the plot a color. There are a variety of ordination methods and distances to use, and not all ordination methods require distances. If in doubt, leave them on default settings.

#### PERMANOVA

The PERMANOVA test can be used to determine the importance of a sample variable in affecting the samples' composition. The test is applied over a user-input number of permutations, and with a dissimilarity index chosen by the user (default is "bray", though there are several more).

Visualization is done through a Top Factors plot, which applies the test in order to determine which OTUs are most important to differentiate samples in two different groups (i.e. given two different nationalities, which OTUs are more abundant in one nationality's sample group and which OTUs are more abundant in the other).

#### Saving Outputs

Individual tables and plots can be saved, either by pressing the Download button (for tables) or the <cameraiconhere> icon (for plots).

Alternatively, an HTML report can be downloaded by clicking "Download report" in the Results tab.
    

