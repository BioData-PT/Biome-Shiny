#Dependency installer script for Biome-shiny

#Functions that check if packages are installed and if not, install them
package_check_and_install = function(packages = c("package1","package2")){
  if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
    install.packages(setdiff(packages, rownames(installed.packages())))
  }
}
package_check_and_install_bioconductor = function(packages = c("package1","package2")){
  if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
    BiocManager::install(setdiff(packages, rownames(installed.packages())))
  }
}

#Packages installed through CRAN

package_check_and_install(c("shiny","shinydashboard","shinyBS","dplyr","knitr","rmarkdown","DT",
                            "ggplot2","ggpubr","hrbrthemes","ggplotify","RColorBrewer",
                            "plotly","heatmaply","reshape2","vegan","abundant","BiocManager"))

#BioConductor install
BiocManager::install()

#Packages installed by Bioconductor
package_check_and_install_bioconductor(c("phyloseq","microbiome","biomformat"))
