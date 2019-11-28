# Biome-shiny
##  A Shiny R app for microbiome visualization, built around the "microbiome" package



#### 1. Features
Biome-shiny has, at its core, the Microbiome and Phyloseq libraries. As such, the test dataset used to demonstrate the app's features can also be found within the standalone Microbiome library.

![](https://i.gyazo.com/8916f6881ec52e16dc943949e1a72225.png)

*Dashboard menu with current features*

The app uses the plotly library to generatw interactive plots, to make it easier to visualize data and save results.

##### 2. Uploading data
Biome-shiny can take three kinds of inputs: The first, a .biom file with included metadata sample variables (in phyloseq, this is the sample_data() function). The second type of input is a .biom file without metadata and a matching .csv mapping file. Finally, the application can generate minimal metadata for a .biom file with no metadata nor mapping file, based on the sample headers of the OTU table.
In addition, microbiome's sample datasets can be used for testing purposes.
This guide, will use the *dietswap* sample dataset to showcase the app's features.

![](https://i.gyazo.com/31c1866b0c1e71d4d3e4cf31bc959de6.png, "Uploading a dataset")

*Dataset upload menu*

Be sure to click "Update Dataset" when selecting the first dataset and when switching from one dataset to another.

##### 3. Phyloseq Summary
The Microbiome library includes a simple dataset summarization function, summarize_phyloseq(). This overview can be examined in the "Phyloseq Summary" tab.

![](https://i.gyazo.com/cb6ad1787f7d273f2d8a4527da9f4924.png)

*Phyloseq summary window*

##### 4. Core Microbiota (Data filtering and visualization)
The "Core Microbiota" tab lets the user filter the active dataset by an input detection and prevalence value. It allows the user to view which species/OTUs remain in the filtered dataset, and get a Phyloseq summary of the filtered data. The application also outputs a prevalence table, allowing the user to know how the number of times a species shows up in samples.

![](https://i.gyazo.com/79c1c556e3bebf93b740095e9264e4bc.png)

*Data filtering variables. Ticking "set as active dataset" will set the filtered data as the working dataset*

The user can also visualize core microbiota through a heatmap, which requires the user to input a minimum detection threshold value.

![](https://i.gyazo.com/0cfbedd65f9ed931683c8fe2534a0417.png)

*Core OTU heatmap. Measures how often OTUs have a certain relative abundance value or higher*

##### 5. Community Composition
This section of the program allows the user to view the absolute and relative abundance of taxa in samples. The abundance plots can be separated by metadata values (i.e. dividing samples by nationality).

###### 5.1. Abundance Plots
![](https://i.gyazo.com/204a382ea6d4dbf414f2f9622a7329cd.png)

*Variables to generate abundance plots*


**Sample variable**, as the name indicates, refers to the metadata that refers to each individual sample in the dataset, although it's possible to use other variables, doing so might produce undesirable results.

Checking the **Group samples by metadata variable** box will prompt the user to pick a metadata variable to sort the samples.

The **Taxonomy** box prompts the user to pick the taxonomic level of abundance plots. 
Ticking **Transparent background** will make the plot's background transparent.

![](https://i.gyazo.com/1d02a7e96edb36f46505de60f0cf1a75.png)

*Abundance plot: Abundance of phylum in samples, sorted by nationality*

![](https://i.gyazo.com/ed4e654daab40b35458e130da9d85840.png)

*Abundance plot: Relative abundance of phylum in samples, sorted by nationality*

##### 6. Alpha Diversity
Alpha diversity refers to the diversity within a particular area, and is usually expressed by a species richness measure. Biome-shiny produces several tables and a richness plot to help illustrate species diversity in samples.
###### 6.1. Evenness Table
Evenness represents how evenly species are distributed in the sample. A high evenness score means species have similar representation in the ecosystem. This table shows various evenness measures per sample.

![](https://i.gyazo.com/909d3914a0f6bc8be54b6ddc5653746d.png)

*Evenness table*

###### 6.2. Abundance Table
Abundance is the number of each species present in the ecosystem. It can be represented as a count (absolute abundance) or a percentage (relative abundance). This section presents the abundance of each species/OTU as counts and percentages.

![](https://i.gyazo.com/76285f4c6154713db0ead40273778484.png)

*Counts table*

![](https://i.gyazo.com/1ddfd89b0e72d384ab7fd1a4084ccae9.png)

*Relative abundance table*

###### 6.3. Metadata Table (with diversity measures)
The metadata table shows the sample variables associated with each sample, with a series of measures of diversity for each sample.

![](https://i.gyazo.com/3600212dbaa096641e64be00cbe5a61e.png)

*Sample variable/Diversity measure table*

###### 6.4. Richness Plot
This plot represents alpha diversity in each sample, through user-selected diversity measures. Sample variables can be used to sort samples and gain more information from the plot.

![](https://i.gyazo.com/9b78a473611946a91b52c044fdcbe52f.png)

*Variables menu*

![](https://i.gyazo.com/a28f039c28e2703d5bb3dbe421fe5fc1.png)

*Richness plot: Simpson diversity index, samples sorted by nationality and colored by sex*

#### 7. Beta Diversity
Beta diversity refers to diversity between samples or ecosystems. It can be measured through the variation of species between ecosystems. Biome-shiny uses ordination plots to visualize beta diversity.

![](https://i.gyazo.com/736787b69c9e067d67e28597f0f3abf4.png)

*Variables menu for the ordination plot*

![](https://i.gyazo.com/a298b34da5e0a0937ee719d95e81010b.png)

*Ordination plot*

#### 8. PERMANOVA Test
The PERMANOVA test is a statistical test that can be used to, for example, find the importance of a certain sample variable for the genetic composition of a certain environmnet (for example, if two distinct diets will lead to a vastly different gut microbiome). In addition to performing the test and exporting results in a table, the application also creates a plot to verify which OTUs are the primary contributors to the differences between samples, and a sample network plot.

##### 8.1. PERMANOVA Test Table

![](https://i.gyazo.com/e117dc0a839fafdfa022a16dff3b58c5.png)

*PERMANOVA test variables*

**Dissimilarity index** refers to the index used for the PERMANOVA test.
**Sample variable** is the sample that will be tested by PERMANOVA.
**Number of permutations** is self-explanatory.

The output is a table with the P-Value, Pr(>F), which assumes a value between 0 and 1, and determines the significance of the tested variable to the genetic composition of the samples.

![](https://i.gyazo.com/e8c2104c6bd38d768cf0e2addd88ef31.png)

*P-Value table*

##### 8.2. Top Factors Plot
This plot represents the top taxa responsible for sample dissimilarity.

![](https://i.gyazo.com/b19512253933ce63b697779969e0ba7a.png)
*Variables menu*
![](https://i.gyazo.com/f2124b70ea959cd5204577e73ccdd50a.png)
*Top factors plot*
##### Network Plot
The network plot can be used to cluster samples or OTUs by similarity.

![](https://i.gyazo.com/e1b09cbcaa23e3258df490e303d5b920.png) ![](https://i.gyazo.com/9ebc5de20cf0a06b0d2f5f9675894e18.png)

*Variables menu. Sample-type network plots need a sample variable to cluster samples into, and optionally a shape to distinguish samples*

![](https://i.gyazo.com/12d0a9863fea71916e4548d770355e4a.png)

*Sample network plot*

#### 9. Exporting Results
Generated plots can be exported to an .html file, which will include a short description of the plot along with the interactive plot itself. This section is still in need of improvement and will be enhanced to include a code output (to generate plots without relying on Biome-shiny). Note that only plots that have been generated with the app will appear in the output file.

![](https://i.gyazo.com/0c6b9b435aaffd489d5cd07e75799c86.png)

*Output menu*

Right now, the only available output format is an .html file. The user can pick which plots and tables to generate in the final report.


#### 10. Citations

Winston Chang, Joe Cheng, JJ Allaire, Yihui Xie and Jonathan McPherson (2018). shiny: Web Application Framework for R. R package version 1.2.0. https://CRAN.R-project.org/package=shiny

Leo Lahti, Sudarshan Shetty et al. (2017). Tools for microbiome analysis in R. Version 1.5.28. URL: http://microbiome.github.com/microbiome.

Wickham H (2016). ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York. ISBN 978-3-319-24277-4, http://ggplot2.org.

Hadley Wickham (2007). Reshaping Data with the reshape Package. Journal of Statistical Software, 21(12), 1-20. URL http://www.jstatsoft.org/v21/i12/.

Hadley Wickham, Romain François, Lionel Henry and Kirill Müller (2018). dplyr: A Grammar of Data Manipulation. R package version package version 0.8.0.1. https://CRAN.R-project.org/package=dplyr

Alboukadel Kassambara (2018). ggpubr: 'ggplot2' Based Publication Ready Plots. R package version 0.2.999. http://www.sthda.com/english/rpkgs/ggpubr

Yihui Xie (2019). knitr: A General-Purpose Package for Dynamic Report Generation in R. R package version 1.22.

Bob Rudis (2019). hrbrthemes: Additional Themes, Theme Components and Utilities for 'ggplot2'. R package version 0.6.0. https://CRAN.R-project.org/package=hrbrthemes

pjames1 @ https://github.com/pjames1 for the plot_ordered_bar function
