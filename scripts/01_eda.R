#-------Dataset Acquisition and Exploration-------

library(Seurat) 
library(tidyverse)
library(pcalg)
library(bnlearn)
library(ggplot2)


#Load Dataset 
tme_data <- readRDS("~/Downloads/oc_tme_causal/raw_object.rds")

#Metadata 
print(head(tme_data@meta.data))

#UMAP plot
library(ggplot2)

# 1. Pull the metadata into a standard data frame
plot_df <- tme_data@meta.data

# 2. UMAP plot of cell type annotation (finest grain)
umap_plot_celltypeannotation <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = Annotation)) +
  geom_point(size = 0.5, alpha = 0.7) +
  theme_minimal() +
  scale_color_viridis_d() +
  labs(
    title = "UMAP Projection",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  theme(legend.position = "bottom")

umap_plot_celltypeannotation

ggsave("figures/UMAP_finest_grain_celltype_annotation.png",
       umap_plot_celltypeannotation,
       width = 14,
       height = 10)

#Maintypes_3 visualization
umap_plot_m3 <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = maintypes_3)) +
  geom_point(size = 0.5, alpha = 0.7) +
  theme_minimal() +
  scale_color_viridis_d() + # Good for distinguishing many samples
  labs(title = "UMAP Projection",
       
       x = "UMAP 1",
       y = "UMAP 2")
#ggsave("figures/UMAP_maintypes3.png", plot=umap_plot_m3, width = 8, height=6, dpi=300)


