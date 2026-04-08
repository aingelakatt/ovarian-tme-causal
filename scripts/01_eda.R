#-------Dataset Acquisition and Exploration-------
Sys.setenv("R_MAX_VSIZE" = "32Gb") #adding more R
library(Seurat) 
library(tidyverse)
library(pcalg)
library(bnlearn)
library(ggplot2)
library(scales)
library(pheatmap)
# -------------------------------------------------
# Memory-safe loading + normalization
# -------------------------------------------------
#Run if starting the pipeline --> don't load the other libraries in only load Seurat --> the dataset that will be created is 15.7GB, need to save as much RAM as possible 
#library(Seurat)


#if (!file.exists("~/Downloads/oc_tme_causal/raw_object_normalized.rds")){
#  tme_data <- readRDS("~/Downloads/oc_tme_causal/raw_object.rds")
#  tme_data <- UpdateSeuratObject(tme_data)
#  gc()
  
  # Normalize in-place on sparse matrix — no extra memory
#  counts_mat <- tme_data@assays$RNA@counts
#  col_sums <- Matrix::colSums(counts_mat)
#  col_idx <- rep(seq_along(diff(counts_mat@p)), diff(counts_mat@p))
#  counts_mat@x <- log1p(counts_mat@x * (10000 / col_sums[col_idx]))
#  tme_data@assays$RNA@data <- counts_mat
#  rm(counts_mat, col_sums, col_idx); gc()
#  saveRDS(tme_data, "~/Downloads/oc_tme_causal/raw_object_normalized.rds")
#  cat("Done.\n")
#}else {
#  tme_data <- readRDS("~/Downloads/oc_tme_causal/raw_object_normalized.rds") #this is 15.4 GB don't run this lol 
#}
#tme_data@assays$RNA@counts <- new("dgCMatrix")
#gc()
#saveRDS(tme_data, "~/Downloads/oc_tme_causal/raw_object_normalized_lite.rds")

# -----------------------------------------------
# Metadata extraction
# -------------------------------------------------

plot_df <- tme_data@meta.data

cat("Ready:", nrow(plot_df), "cells\n")

print(head(plot_df))

#----------GLOBAL UMAP PLOTS --------------------


#1a. UMAP — finest grain cell type annotation
umap_plot_celltypeannotation <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = Annotation)) +
  geom_point(size = 0.5, alpha = 0.7) +
  theme_minimal() +
  scale_color_viridis_d() +
  labs(title = "UMAP — Fine Annotation", x = "UMAP 1", y = "UMAP 2") +
  theme(legend.position = "bottom")

umap_plot_celltypeannotation
ggsave("figures/fig01_umap_fine_annotation.png",
       umap_plot_celltypeannotation, width = 14, height = 10, dpi = 300)

# 1b. UMAP — broadest grouping (maintypes_3)
umap_plot_m3 <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = maintypes_3)) +
  geom_point(size = 0.5, alpha = 0.7) +
  theme_minimal() +
  scale_color_viridis_d() +
  labs(title = "UMAP — Broadest Grouping (maintypes_3)", x = "UMAP 1", y = "UMAP 2")
umap_plot_m3
ggsave("figures/fig02_umap_maintypes3.png", umap_plot_m3, width = 8, height = 6, dpi = 300)

# 1c. UMAP — broad cell types (maintypes_2)
umap_plot_m2 <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = maintypes_2)) +
  geom_point(size = 0.5, alpha = 0.7) +
  theme_minimal() +
  scale_color_viridis_d() +
  labs(title = "UMAP — Broad Cell Types (maintypes_2)", x = "UMAP 1", y = "UMAP 2")
umap_plot_m2
ggsave("figures/fig03_umap_maintypes2.png", umap_plot_m2, width = 10, height = 7, dpi = 300)

# 1d. UMAP — tissue site
umap_plot_tissue <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = Groups)) +
  geom_point(size = 0.3, alpha = 0.5) +
  theme_minimal() +
  scale_color_brewer(palette = "Set1") +
  labs(title = "UMAP — Tissue Site", x = "UMAP 1", y = "UMAP 2")
umap_plot_tissue
ggsave("figures/fig04_umap_tissue.png", umap_plot_tissue, width = 9, height = 7, dpi = 300)

# 1e. UMAP — patient
umap_plot_patient <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = Patients)) +
  geom_point(size = 0.3, alpha = 0.5) +
  theme_minimal() +
  scale_color_viridis_d() +
  labs(title = "UMAP — Patient", x = "UMAP 1", y = "UMAP 2")
umap_plot_patient
ggsave("figures/fig05_umap_patient.png", umap_plot_patient, width = 10, height = 7, dpi = 300)

# 1f. UMAP — cell types split by tissue site (small multiples)
umap_split_tissue <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = maintypes_2)) +
  geom_point(size = 0.2, alpha = 0.5) +
  facet_wrap(~Groups, ncol = 3) +
  theme_minimal(base_size = 10) +
  scale_color_viridis_d() +
  labs(title = "Cell types per tissue site", x = "UMAP 1", y = "UMAP 2") +
  theme(legend.position = "bottom")
umap_split_tissue
ggsave("figures/fig06_umap_split_by_tissue.png", umap_split_tissue, width = 16, height = 10, dpi = 300)

#---------------SUMMARY TABLES-----------------------

# The object is a Seurat v3 object loaded in v5, so use dim() on counts directly
n_cells <- ncol(tme_data@assays$RNA@counts)
n_genes <- nrow(tme_data@assays$RNA@counts)
cat("Cells:", n_cells, "\n")
cat("Genes:", n_genes, "\n")
cat("Patients:", length(unique(plot_df$Patients)), "\n")
cat("Samples:", length(unique(plot_df$Samples)), "\n")

# ----- Check actual factor levels -----
cat("\n=== Groups (tissue site) levels ===\n")
print(levels(plot_df$Groups))

cat("\n=== maintypes_2 levels ===\n")
print(levels(plot_df$maintypes_2))

cat("\n=== Patients levels ===\n")
print(levels(plot_df$Patients))



# Table 1: Dataset overview

tab01 <- data.frame(
  Metric = c("Total cells", "Total genes", "Patients", "Samples", "Tissue sites"),
  Value = c(n_cells, n_genes,
            length(unique(plot_df$Patients)),
            length(unique(plot_df$Samples)),
            length(unique(plot_df$Groups)))
)
write.csv(tab01, "tables/tab01_dataset_overview.csv", row.names = FALSE)
print(tab01)


# Table 2: Cells per patient
tab02 <- plot_df %>%
  count(Patients, name = "n_cells") %>%
  mutate(pct = round(100 * n_cells / sum(n_cells), 2)) %>%
  arrange(desc(n_cells))
print(tab02)
write.csv(tab02, "tables/tab02_cells_per_patient.csv", row.names = FALSE)

# Table 3: Cells per tissue site
tab03 <- plot_df %>%
  count(Groups, name = "n_cells") %>%
  mutate(pct = round(100 * n_cells / sum(n_cells), 2)) %>%
  arrange(desc(n_cells))
print(tab03)
write.csv(tab03, "tables/tab03_cells_per_tissue.csv", row.names = FALSE)

# Table 4: Cells per sample
tab04 <- plot_df %>%
  count(Samples, Patients, Groups, name = "n_cells") %>%
  mutate(pct = round(100 * n_cells / sum(n_cells), 2)) %>%
  arrange(Patients, Groups)
print(tab04)
write.csv(tab04, "tables/tab04_cells_per_sample.csv", row.names = FALSE)

# Table 5: Broad cell types (maintypes_2) — THE KEY TABLE
tab05 <- plot_df %>%
  count(maintypes_2, name = "n_cells") %>%
  mutate(pct = round(100 * n_cells / sum(n_cells), 2)) %>%
  arrange(desc(n_cells))
write.csv(tab05, "tables/tab05_broad_celltypes.csv", row.names = FALSE)
cat("\n=== Broad cell types (maintypes_2) ===\n")
print(tab05)

# Table 6: Broadest grouping (maintypes_3)
tab06 <- plot_df %>%
  count(maintypes_3, name = "n_cells") %>%
  mutate(pct = round(100 * n_cells / sum(n_cells), 2)) %>%
  arrange(desc(n_cells))
print(tab06)
write.csv(tab06, "tables/tab06_maintypes3.csv", row.names = FALSE)

# Table 7: All fine annotations
tab07 <- plot_df %>%
  count(Annotation, maintypes_2, maintypes_3, name = "n_cells") %>%
  mutate(pct = round(100 * n_cells / sum(n_cells), 2)) %>%
  arrange(desc(n_cells))
print(tab07)
write.csv(tab07, "tables/tab07_fine_annotations.csv", row.names = FALSE)

# Table 8: Patient x Tissue (cell counts)
tab08 <- plot_df %>%
  count(Patients, Groups) %>%
  pivot_wider(names_from = Groups, values_from = n, values_fill = 0)
write.csv(tab08, "tables/tab08_patient_by_tissue.csv", row.names = FALSE)

# Table 9: Cell type x Tissue (counts)
tab09 <- plot_df %>%
  count(maintypes_2, Groups) %>%
  pivot_wider(names_from = Groups, values_from = n, values_fill = 0)
print(tab09)
write.csv(tab09, "tables/tab09_celltype_by_tissue.csv", row.names = FALSE)

# Table 10: Cell type x Tissue (% within each tissue)
tab10 <- plot_df %>%
  count(Groups, maintypes_2) %>%
  group_by(Groups) %>%
  mutate(pct = round(100 * n / sum(n), 2)) %>%
  select(-n) %>%
  pivot_wider(names_from = Groups, values_from = pct, values_fill = 0)
print(tab10)
write.csv(tab10, "tables/tab10_celltype_by_tissue_pct.csv", row.names = FALSE)

# Table 11: Cell type x Patient (% within each patient)
tab11 <- plot_df %>%
  count(Patients, maintypes_2) %>%
  group_by(Patients) %>%
  mutate(pct = round(100 * n / sum(n), 2)) %>%
  select(-n) %>%
  pivot_wider(names_from = Patients, values_from = pct, values_fill = 0)
print(tab11)
write.csv(tab11, "tables/tab11_celltype_by_patient_pct.csv", row.names = FALSE)


#----------------Bar Plots--------------------------------


# 3a. Cell type proportions per patient
fig07 <- plot_df %>%
  count(Patients, maintypes_2) %>%
  group_by(Patients) %>%
  mutate(pct = n / sum(n)) %>%
  ggplot(aes(x = Patients, y = pct, fill = maintypes_2)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = percent) +
  scale_fill_viridis_d() +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Patient", y = "Proportion", fill = "Cell Type",
       title = "Cell type composition per patient")
fig07
ggsave("figures/fig07_composition_per_patient.png", fig07, width = 12, height = 7, dpi = 300)

# 3b. Cell type proportions per tissue site
fig08 <- plot_df %>%
  count(Groups, maintypes_2) %>%
  group_by(Groups) %>%
  mutate(pct = n / sum(n)) %>%
  ggplot(aes(x = Groups, y = pct, fill = maintypes_2)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = percent) +
  scale_fill_viridis_d() +
  theme_minimal(base_size = 12) + theme(axis.text.x = element_text(angle=45))+
  labs(x = "Tissue Site", y = "Proportion", fill = "Cell Type",
       title = "Cell type composition per tissue site")
fig08
ggsave("figures/fig08_composition_per_tissue.png", fig08, width = 9, height = 7, dpi = 300)

# 3c. Cell type proportions per sample (faceted by tissue)
fig09 <- plot_df %>%
  count(Samples, Groups, maintypes_2) %>%
  group_by(Samples) %>%
  mutate(pct = n / sum(n)) %>%
  ggplot(aes(x = Samples, y = pct, fill = maintypes_2)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = percent) +
  scale_fill_viridis_d() +
  facet_wrap(~Groups, scales = "free_x", nrow = 1) +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 70, hjust = 1, size = 6),
        legend.position = "bottom") +
  labs(x = "", y = "Proportion", fill = "Cell Type",
       title = "Cell type composition per sample")
fig09
ggsave("figures/fig09_composition_per_sample.png", fig09, width = 20, height = 8, dpi = 300)

# 3d. Total cells per sample
fig10 <- plot_df %>%
  count(Samples, Groups) %>%
  ggplot(aes(x = fct_reorder(Samples, -n), y = n, fill = Groups)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7)) +
  labs(x = "Sample", y = "Cell Count", fill = "Tissue",
       title = "Total cells per sample")
fig10
ggsave("figures/fig10_cells_per_sample.png", fig10, width = 14, height = 6, dpi = 300)

#--------------------- T cell exploration-------------------------------------------


# The 14 maintypes_2 levels include "CD4+ T", "CD8+ T",
# Filter from metadata
t_df <- plot_df %>%
  filter(maintypes_2 %in% c("CD4+ T", "CD8+ T"))
cat("CD4/CD8:", nrow(t_df), "\n")

# Table 12: T cell subtypes
tab12 <- t_df %>%
  count(Annotation, maintypes_2, name = "n_cells") %>%
  mutate(pct = round(100 * n_cells / sum(n_cells), 2)) %>%
  arrange(desc(n_cells))
write.csv(tab12, "tables/tab12_tcell_subtypes.csv", row.names = FALSE)
cat("\nT cell subtypes:\n")
print(tab12)


# --- Seurat subset for gene expression plots ---

t_cells_seurat <- subset(tme_data, cells = rownames(t_df))

# Add UMAP reduction to subset
umap_t <- as.matrix(t_df[, c("UMAP_1", "UMAP_2")])
rownames(umap_t) <- rownames(t_df)
t_cells_seurat[["umap"]] <- CreateDimReducObject(
  embeddings = umap_t, key = "UMAP_", assay = DefaultAssay(t_cells_seurat))

# Fig 11: T cell markers DotPlot — from Zheng et al. paper
t_markers <- c("CD4", "CD8A", "SELL", "CCR7", "ANXA1",
               "FOXP3", "IL2RA", "CTLA4", "CXCL13", "ANXA2",
               "GZMK", "CCL4", "CX3CR1", "GNLY", "FGFBP2",
               "GZMB", "TNFRSF9", "HAVCR2", "SLC4A10", "TRDV2")
t_markers_ok <- t_markers[t_markers %in% rownames(t_cells_seurat)]
cat("T markers found:", length(t_markers_ok), "/", length(t_markers), "\n")
cat("Missing:", paste(t_markers[!t_markers %in% rownames(t_cells_seurat)], collapse = ", "), "\n")

fig11 <- DotPlot(t_cells_seurat, features = t_markers_ok, group.by = "Annotation") +
  coord_flip() +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("T cell markers across subtypes (Zheng et al.)")
fig11
ggsave("figures/fig11_tcell_markers_dotplot.png", fig11, width = 14, height = 9, dpi = 300)


# --- Fig 12: Heatmap of T cell marker genes × clusters ---

# Define marker genes per cluster (from Zheng et al.)
# Organized by functional group
t_heatmap_genes <- c(
  # CD4 naive/memory
  "CCR7", "SELL", "LEF1",
  # CD4 effector memory
  "ANXA1", "S100A4", "IL7R",
  # Treg
  "FOXP3", "IL2RA", "CTLA4", "IKZF2",
  # CD4 cytotoxic
  "CX3CR1", "FGFBP2", "GNLY",
  # CD4 exhausted
  "PDCD1", "CXCL13", "TOX",
  # CD8 naive
  "TCF7", "CD8A", "CD8B",
  # CD8 effector memory
  "ANXA2", "GZMK", "CCL4", "CCL5",
  # CD8 cytotoxic
  "GZMB", "PRF1", "NKG7", "GNLY",
  # CD8 exhausted
  "HAVCR2", "LAG3", "TIGIT", "TNFRSF9",
  # MAIT
  "SLC4A10", "KLRB1", "ZBTB16",
  # gamma-delta
  "TRDV2", "TRDC"
)

# Remove duplicates and check availability
t_heatmap_genes <- unique(t_heatmap_genes)
genes_ok <- t_heatmap_genes[t_heatmap_genes %in% rownames(t_cells_seurat)]
cat("Heatmap genes found:", length(genes_ok), "/", length(t_heatmap_genes), "\n")
if (length(genes_ok) < length(t_heatmap_genes)) {
  cat("Missing:", paste(setdiff(t_heatmap_genes, genes_ok), collapse = ", "), "\n")
}

# Compute mean expression per cluster
# Pull normalized expression for these genes
expr_mat <- GetAssayData(t_cells_seurat, layer = "data")[genes_ok, ]
cluster_ids <- t_cells_seurat@meta.data$Annotation

# Average expression per cluster
clusters <- sort(unique(cluster_ids))
avg_expr <- matrix(0, nrow = length(genes_ok), ncol = length(clusters))
rownames(avg_expr) <- genes_ok
colnames(avg_expr) <- clusters

for (cl in clusters) {
  cells_in <- which(cluster_ids == cl)
  if (length(cells_in) > 0) {
    avg_expr[, cl] <- Matrix::rowMeans(expr_mat[, cells_in, drop = FALSE])
  }
}

# Scale rows (z-score per gene) for better visualization
avg_scaled <- t(scale(t(avg_expr)))


png("figures/fig12_tcell_marker_heatmap.png",
    width = 10, height = 12, units = "in", res = 300)
pheatmap(avg_scaled,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
         fontsize_row = 10,
         fontsize_col = 10,
         angle_col = 45,
         main = "T cell marker genes across clusters (row-scaled)",
         gaps_col = NULL,
         border_color = NA)
dev.off()
cat("Saved: fig12_tcell_marker_heatmap.png\n")

# Also display in RStudio
pheatmap(avg_scaled,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
         fontsize_row = 10,
         fontsize_col = 10,
         angle_col = 45,
         main = "T cell marker genes across clusters (row-scaled)",
         border_color = NA)


# Fig 13: Feature plots — cytotoxic vs exhaustion
feat_genes <- c("GZMB", "PRF1", "NKG7", "PDCD1", "HAVCR2", "TIGIT")
feat_ok <- feat_genes[feat_genes %in% rownames(t_cells_seurat)]

fig13 <- FeaturePlot(t_cells_seurat, features = feat_ok,
                     reduction = "umap", ncol = 3, order = TRUE) &
  theme_minimal(base_size = 9)
print(fig13)
ggsave("figures/fig13_tcell_cyto_exh_features.png", fig13, width = 14, height = 9, dpi = 300)


#------------Myeloid Cells-----------------

# Subset myeloid cells using M01-M15 annotation prefix
mac_df <- plot_df %>%
  filter(grepl("^M\\d", Annotation))
cat("Myeloid cells (M01-M15):", nrow(mac_df), "\n")

# Table 13: Myeloid cell subtypes
tab13 <- mac_df %>%
  count(Annotation, maintypes_2, name = "n_cells") %>%
  mutate(pct = round(100 * n_cells / sum(n_cells), 2)) %>%
  arrange(desc(n_cells))
write.csv(tab13, "tables/tab13_mac_cell_subtypes.csv", row.names = FALSE)
cat("\n Macrophage cell subtypes:\n")
print(tab13)








