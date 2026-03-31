# Name: 00_setup.R 
# Description: Package Setup 

#Required Packages
required_pkgs <- c("Seurat", "ggplot2", "bnlearn", "tidyverse", "pcalg", "dplyr", "patchwork", "tidyr", "pheatmap", "RColorBrewer", "scales", "forcats", "SingleCellExperiment")

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% c("Seurat", "SingleCellExperiment")) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
  }
  
}
