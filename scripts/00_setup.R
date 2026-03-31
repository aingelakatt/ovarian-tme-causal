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
# The raw_object.rds contains cells that already passed:
#   - Min 200 genes expressed per cell
#   - Max 10% mitochondrial UMIs
#   - Doublet removal (cells expressing >1 major marker)
#   - NormalizeData (LogNormalize) + ScaleData
#
# The deposited object retained raw counts but not the
# normalized data slot. We re-run log-normalization
# (identical parameters) and save as raw_object_normalized.rds --> rerunning in 01_eda.R
#
# Source: Zheng et al. Nat Cancer 4, 1138–1156 (2023)