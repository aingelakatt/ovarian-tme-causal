
Sys.setenv("R_MAX_VSIZE" = "24Gb")

library(Seurat)
library(tidyverse)
library(scales)
library(pheatmap)

# --- Load normalized data ---
tme_data <- readRDS("~/Downloads/oc_tme_causal/raw_object_normalized_lite.rds")
plot_df <- tme_data@meta.data
cat("Loaded:", nrow(plot_df), "cells\n")

dir.create("tables", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)


