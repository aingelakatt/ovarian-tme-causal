# ==============================================================================
# 02_feature_engineering.R
# Ovarian TME Causal Inference Pipeline — Phase 2
# 
# Purpose: Build sample-level feature matrix for causal DAG construction
#   - Cell type proportions (maintypes_2, fine Annotation)
#   - Tissue-specific features
#   - T cell / Myeloid subtype ratios
#
# Input:  Metadata from normalized Seurat object (or saved plot_df)
# Output: Feature matrix CSV for causal analysis
#
# Runs on: LOCAL (lightweight — only needs metadata, not full Seurat object)
# Heavy computation (MiloR): See 02a_miloR_crc.R for CRC
# ==============================================================================

library(tidyverse)
library(scales)

# --- Configuration ------------------------------------------------------------
DATA_DIR <- "~/Downloads/oc_tme_causal"
OUT_DIR <- "results"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# --- Load Data ----------------------------------------------------------------
# Option 1: If you saved plot_df from 01_eda.R
if (file.exists(file.path(DATA_DIR, "metadata.rds"))) {
  cat("Loading saved metadata...\n")
  plot_df <- readRDS(file.path(DATA_DIR, "metadata.rds"))
} else {
  # Option 2: Extract from Seurat object
  cat("Loading Seurat object to extract metadata...\n")
  library(Seurat)
  
  # Use the lite version if available (no counts matrix)
  if (file.exists(file.path(DATA_DIR, "raw_object_normalized_lite.rds"))) {
    tme_data <- readRDS(file.path(DATA_DIR, "raw_object_normalized_lite.rds"))
  } else {
    tme_data <- readRDS(file.path(DATA_DIR, "raw_object_normalized.rds"))
  }
  
  plot_df <- tme_data@meta.data
  
  # Save metadata separately for future lightweight runs
  saveRDS(plot_df, file.path(DATA_DIR, "metadata.rds"))
  cat("Saved metadata to metadata.rds for future runs\n")
  
  rm(tme_data); gc()
}

cat("Loaded metadata:", nrow(plot_df), "cells\n")
cat("Columns:", paste(colnames(plot_df), collapse = ", "), "\n\n")

# --- Verify required columns --------------------------------------------------
required_cols <- c("Patients", "Samples", "Groups", "maintypes_2", "Annotation")
missing <- setdiff(required_cols, colnames(plot_df))
if (length(missing) > 0) {
  stop("Missing required columns: ", paste(missing, collapse = ", "))
}

# ==============================================================================
# FEATURE 1: CELL TYPE PROPORTIONS (maintypes_2) PER SAMPLE
# ==============================================================================
cat("Computing cell type proportions per sample...\n")

# Proportions of broad cell types per sample
prop_maintypes2 <- plot_df %>%
  count(Samples, Patients, Groups, maintypes_2) %>%
  group_by(Samples) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  select(-n) %>%
  pivot_wider(
    names_from = maintypes_2,
    values_from = prop,
    values_fill = 0,
    names_prefix = "prop_"
  )

# Clean column names (remove spaces, special chars)
names(prop_maintypes2) <- gsub("[^[:alnum:]_]", "_", names(prop_maintypes2))
names(prop_maintypes2) <- gsub("__+", "_", names(prop_maintypes2))
names(prop_maintypes2) <- gsub("_$", "", names(prop_maintypes2))

cat("  Created", ncol(prop_maintypes2) - 3, "cell type proportion features\n")

# ==============================================================================
# FEATURE 2: T CELL SUBTYPE PROPORTIONS (within T cells)
# ==============================================================================
cat("Computing T cell subtype features...\n")

# T cell annotations (from your EDA: CD4+ T, CD8+ T in maintypes_2)
t_cell_df <- plot_df %>%
  filter(maintypes_2 %in% c("CD4+ T", "CD8+ T"))

# Proportions of T cell subtypes per sample (within T cells only)
prop_tcell <- t_cell_df %>%
  count(Samples, Annotation) %>%
  group_by(Samples) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  select(-n) %>%
  pivot_wider(
    names_from = Annotation,
    values_from = prop,
    values_fill = 0,
    names_prefix = "Tprop_"
  )

# Clean names
names(prop_tcell) <- gsub("[^[:alnum:]_]", "_", names(prop_tcell))
names(prop_tcell) <- gsub("__+", "_", names(prop_tcell))

cat("  Created", ncol(prop_tcell) - 1, "T cell subtype features\n")

# ==============================================================================
# FEATURE 3: MYELOID SUBTYPE PROPORTIONS (within myeloid)
# ==============================================================================
cat("Computing myeloid subtype features...\n")

# Myeloid annotations (M01-M15 pattern from your EDA)
myeloid_df <- plot_df %>%
  filter(grepl("^M\\d", Annotation))

if (nrow(myeloid_df) > 0) {
  prop_myeloid <- myeloid_df %>%
    count(Samples, Annotation) %>%
    group_by(Samples) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup() %>%
    select(-n) %>%
    pivot_wider(
      names_from = Annotation,
      values_from = prop,
      values_fill = 0,
      names_prefix = "Mprop_"
    )
  
  names(prop_myeloid) <- gsub("[^[:alnum:]_]", "_", names(prop_myeloid))
  cat("  Created", ncol(prop_myeloid) - 1, "myeloid subtype features\n")
} else {
  prop_myeloid <- NULL
  cat("  No myeloid cells found with M## pattern\n")
}

# ==============================================================================
# FEATURE 4: DERIVED RATIOS (biologically meaningful)
# ==============================================================================
cat("Computing derived ratio features...\n")

# Sample-level summary for ratios
sample_summary <- plot_df %>%
  group_by(Samples, Patients, Groups) %>%
  summarize(
    n_cells = n(),
    # T cell counts
    n_CD4 = sum(maintypes_2 == "CD4+ T"),
    n_CD8 = sum(maintypes_2 == "CD8+ T"),
    n_Tcells = n_CD4 + n_CD8,
    # Other immune
    n_Myeloid = sum(grepl("^M\\d", Annotation)),
    n_NK = sum(maintypes_2 == "NK"),
    n_B = sum(maintypes_2 == "B"),
    # Stromal
    n_Fibro = sum(maintypes_2 == "Fibroblasts"),
    n_Endo = sum(maintypes_2 == "Endothelial cells"),
    # Tumor
    n_Tumor = sum(maintypes_2 == "Epithelial cells"),
    .groups = "drop"
  ) %>%
  mutate(
    # Key ratios (add small constant to avoid div by zero)
    CD8_CD4_ratio = (n_CD8 + 1) / (n_CD4 + 1),
    Tcell_Myeloid_ratio = (n_Tcells + 1) / (n_Myeloid + 1),
    Immune_Tumor_ratio = (n_Tcells + n_Myeloid + n_NK + n_B + 1) / (n_Tumor + 1),
    Lymphoid_Myeloid_ratio = (n_Tcells + n_NK + n_B + 1) / (n_Myeloid + 1),
    # Proportions
    prop_immune = (n_Tcells + n_Myeloid + n_NK + n_B) / n_cells,
    prop_stromal = (n_Fibro + n_Endo) / n_cells
  )

cat("  Created ratio features: CD8/CD4, T/Myeloid, Immune/Tumor, etc.\n")

# ==============================================================================
# FEATURE 5: T CELL FUNCTIONAL STATES (based on Sun et al. 2023)
# ==============================================================================
cat("Computing T cell functional state features...\n")

# Detailed T cell states (7 categories)
t_cell_states <- t_cell_df %>%
  mutate(
    functional_state = case_when(
      # --- EXHAUSTED: HAVCR2 (TIM-3), CXCL13 ---
      Annotation %in% c("T10_CD8-HAVCR2", "T05_CD4-CXCL13") ~ "exhausted",
      
      # --- NAIVE: CCR7+ (lymph node homing, quiescent) ---
      Annotation %in% c("T01_CD4-CCR7", "T06_CD8-CCR7") ~ "naive",
      
      # --- CYTOTOXIC: CX3CR1+ (terminally differentiated killers) ---
      Annotation %in% c("T04_CD4-CX3CR1", "T09_CD8-CX3CR1") ~ "cytotoxic",
      
      # --- EFFECTOR MEMORY: ANXA1/2, GZMK (activated, tissue-resident capable) ---
      Annotation %in% c("T02_CD4-ANXA1", "T07_CD8-ANXA2", "T08_CD8-GZMK") ~ "effector_memory",
      
      # --- REGULATORY: FOXP3+ (immunosuppressive Tregs) ---
      Annotation == "T03_CD4-FOXP3" ~ "regulatory",
      
      # --- MAIT: SLC4A10+ (innate-like, mucosal) ---
      Annotation == "T11_CD8-SLC4A10" ~ "MAIT",
      
      # --- GAMMA-DELTA: TRDV2+ (unconventional T cells) ---
      Annotation == "T12_CD8-TRDV2" ~ "gamma_delta",
      
      TRUE ~ "other_T"
    )
  )

# Simplified T cell states (5 categories - better for causal analysis)
t_cell_states <- t_cell_states %>%
  mutate(
    functional_state_simple = case_when(
      functional_state == "exhausted" ~ "exhausted",
      functional_state == "cytotoxic" ~ "cytotoxic",
      functional_state == "regulatory" ~ "regulatory",
      functional_state == "naive" ~ "naive",
      TRUE ~ "effector_memory"  # Collapse ANXA, GZMK, MAIT, gamma-delta
    )
  )

# Detailed state proportions
prop_tcell_states <- t_cell_states %>%
  count(Samples, functional_state) %>%
  group_by(Samples) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  select(-n) %>%
  pivot_wider(
    names_from = functional_state,
    values_from = prop,
    values_fill = 0,
    names_prefix = "Tstate_"
  )

# Simplified state proportions (for causal DAG)
prop_tcell_states_simple <- t_cell_states %>%
  count(Samples, functional_state_simple) %>%
  group_by(Samples) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  select(-n) %>%
  pivot_wider(
    names_from = functional_state_simple,
    values_from = prop,
    values_fill = 0,
    names_prefix = "Tstate_simple_"
  )

cat("  Created", ncol(prop_tcell_states) - 1, "detailed T cell state features\n")
cat("  Created", ncol(prop_tcell_states_simple) - 1, "simplified T cell state features\n")

cat("  T cell states identified (detailed):\n")
print(table(t_cell_states$functional_state))
cat("  T cell states identified (simplified):\n")
print(table(t_cell_states$functional_state_simple))

# ==============================================================================
# FEATURE 6: MYELOID FUNCTIONAL STATES (based on Sun et al. 2023)
# ==============================================================================
cat("Computing myeloid functional state features...\n")

if (nrow(myeloid_df) > 0) {
  
  # Detailed myeloid states (10 categories)
  myeloid_states <- myeloid_df %>%
    mutate(
      functional_state = case_when(
        # --- cDC1: CLEC9A+ (cross-presenting, key for CD8 T cell activation) ---
        Annotation == "M02_DC-CLEC9A" ~ "cDC1",
        
        # --- cDC2: CD1C+ (activate CD4 T cells) ---
        Annotation == "M01_DC-CD1C" ~ "cDC2",
        
        # --- Mature DC: LAMP3+ (activated, migrating to lymph nodes) ---
        Annotation == "M03_DC-LAMP3" ~ "mature_DC",
        
        # --- Other DC ---
        Annotation == "M04_DC-LGALS2" ~ "DC_other",
        
        # --- Classical monocytes: CD14+ ---
        Annotation == "M05_Mono-CD14" ~ "classical_mono",
        
        # --- Non-classical monocytes: FCGR3A+ (CD16+) ---
        Annotation == "M06_Mono-FCGR3A" ~ "nonclassical_mono",
        
        # --- M2-like TAMs: FOLR2+ (immunosuppressive, pro-tumor) ---
        Annotation %in% c("M14_Macro-FOLR2", "M15_Proliferative-Macro-FOLR2") ~ "TAM_M2",
        
        # --- Complement-associated macrophages: C1QA, C3 ---
        Annotation %in% c("M10_Macro-C1QA", "M12_Macro-C3", "M13_Proliferative-Macro-C3") ~ "complement_mac",
        
        # --- Inflammatory/tissue remodeling macrophages ---
        Annotation %in% c("M07_Macro-EREG", "M08_Macro-FN1", "M11_Macro-VCAN") ~ "inflammatory_mac",
        
        # --- Lipid-associated macrophages: FABP5+ ---
        Annotation == "M09_Macro-FABP5" ~ "lipid_mac",
        
        TRUE ~ "other_myeloid"
      )
    )
  
  # Simplified myeloid states (4 categories - better for causal analysis)
  myeloid_states <- myeloid_states %>%
    mutate(
      functional_state_simple = case_when(
        functional_state %in% c("TAM_M2", "lipid_mac", "complement_mac") ~ "TAM_immunosuppressive",
        functional_state %in% c("cDC1", "cDC2", "mature_DC", "DC_other") ~ "dendritic_cell",
        functional_state %in% c("classical_mono", "nonclassical_mono") ~ "monocyte",
        TRUE ~ "macrophage_other"
      )
    )
  
  # Detailed state proportions
  prop_myeloid_states <- myeloid_states %>%
    count(Samples, functional_state) %>%
    group_by(Samples) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup() %>%
    select(-n) %>%
    pivot_wider(
      names_from = functional_state,
      values_from = prop,
      values_fill = 0,
      names_prefix = "Mstate_"
    )
  
  # Simplified state proportions (for causal DAG)
  prop_myeloid_states_simple <- myeloid_states %>%
    count(Samples, functional_state_simple) %>%
    group_by(Samples) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup() %>%
    select(-n) %>%
    pivot_wider(
      names_from = functional_state_simple,
      values_from = prop,
      values_fill = 0,
      names_prefix = "Mstate_simple_"
    )
  
  cat("  Created", ncol(prop_myeloid_states) - 1, "detailed myeloid state features\n")
  cat("  Created", ncol(prop_myeloid_states_simple) - 1, "simplified myeloid state features\n")
  
  cat("  Myeloid states identified (detailed):\n")
  print(table(myeloid_states$functional_state))
  cat("  Myeloid states identified (simplified):\n")
  print(table(myeloid_states$functional_state_simple))
  
} else {
  prop_myeloid_states <- NULL
  prop_myeloid_states_simple <- NULL
  cat("  No myeloid cells found\n")
}

# ==============================================================================
# MERGE ALL FEATURES INTO FINAL MATRIX
# ==============================================================================
cat("\nMerging all features into final matrix...\n")

# FULL feature matrix (all granularities - for exploration)
feature_matrix_full <- sample_summary %>%
  select(Samples, Patients, Groups, n_cells, 
         CD8_CD4_ratio, Tcell_Myeloid_ratio, Immune_Tumor_ratio,
         Lymphoid_Myeloid_ratio, prop_immune, prop_stromal) %>%
  left_join(prop_maintypes2 %>% select(-Patients, -Groups), by = "Samples") %>%
  left_join(prop_tcell, by = "Samples") %>%
  left_join(prop_tcell_states, by = "Samples") %>%
  left_join(prop_tcell_states_simple, by = "Samples")

if (!is.null(prop_myeloid)) {
  feature_matrix_full <- feature_matrix_full %>%
    left_join(prop_myeloid, by = "Samples")
}
if (!is.null(prop_myeloid_states)) {
  feature_matrix_full <- feature_matrix_full %>%
    left_join(prop_myeloid_states, by = "Samples")
}
if (!is.null(prop_myeloid_states_simple)) {
  feature_matrix_full <- feature_matrix_full %>%
    left_join(prop_myeloid_states_simple, by = "Samples")
}

# CAUSAL feature matrix (simplified - for DAG)
feature_matrix_causal <- sample_summary %>%
  select(Samples, Patients, Groups, n_cells, 
         CD8_CD4_ratio, Tcell_Myeloid_ratio, Immune_Tumor_ratio,
         Lymphoid_Myeloid_ratio, prop_immune, prop_stromal) %>%
  left_join(prop_maintypes2 %>% select(-Patients, -Groups), by = "Samples") %>%
  left_join(prop_tcell_states_simple, by = "Samples")

if (!is.null(prop_myeloid_states_simple)) {
  feature_matrix_causal <- feature_matrix_causal %>%
    left_join(prop_myeloid_states_simple, by = "Samples")
}

# Replace NAs with 0 for proportion columns
prop_cols_full <- grep("^prop_|^Tprop_|^Mprop_|^Tstate_|^Mstate_", names(feature_matrix_full))
feature_matrix_full[prop_cols_full] <- lapply(feature_matrix_full[prop_cols_full], function(x) ifelse(is.na(x), 0, x))

prop_cols_causal <- grep("^prop_|^Tstate_|^Mstate_", names(feature_matrix_causal))
feature_matrix_causal[prop_cols_causal] <- lapply(feature_matrix_causal[prop_cols_causal], function(x) ifelse(is.na(x), 0, x))

cat("Full feature matrix:", nrow(feature_matrix_full), "samples x", ncol(feature_matrix_full), "columns\n")
cat("Causal feature matrix:", nrow(feature_matrix_causal), "samples x", ncol(feature_matrix_causal), "columns\n")

# ==============================================================================
# SAVE OUTPUTS
# ==============================================================================

# Full feature matrix
write.csv(feature_matrix_full, file.path(OUT_DIR, "feature_matrix_full.csv"), row.names = FALSE)
saveRDS(feature_matrix_full, file.path(OUT_DIR, "feature_matrix_full.rds"))
cat("Saved: feature_matrix_full.csv\n")

# Causal feature matrix 
write.csv(feature_matrix_causal, file.path(OUT_DIR, "feature_matrix_causal.csv"), row.names = FALSE)
saveRDS(feature_matrix_causal, file.path(OUT_DIR, "feature_matrix_causal.rds"))
cat("Saved: feature_matrix_causal.csv\n")

# --- Select features for causal analysis -------------------------------------
# Rule of thumb: n_samples > 5 * n_variables
n_samples <- nrow(feature_matrix_causal)
max_vars <- floor(n_samples / 5)
cat("\n=== CAUSAL VARIABLE SELECTION ===\n")
cat("Samples:", n_samples, "\n")
cat("Recommended max variables:", max_vars, "\n\n")

# Identify numeric feature columns (exclude metadata)
numeric_features <- feature_matrix_causal %>%
  select(where(is.numeric)) %>%
  select(-n_cells)

# Compute variance of each feature
feature_var <- sapply(numeric_features, var, na.rm = TRUE) %>%
  sort(decreasing = TRUE)

cat("Top 15 features by variance:\n")
print(head(feature_var, 15))

# Create recommended subset for causal DAG
# Priority: key cell populations, ratios, functional states
# NOTE: Using ACTUAL column names from your data (simplified states use "Tstate_simple_" prefix)
recommended_vars <- c(
  # Key cell populations
  "prop_CD8_T", "prop_CD4_T", "prop_Epithelial_cells", "prop_Fibroblast",
  "prop_NK", "prop_B", "prop_Endothelial_cells",
  # Ratios
  "CD8_CD4_ratio", "Immune_Tumor_ratio", "Tcell_Myeloid_ratio",
  # Simplified T cell states
  "Tstate_simple_exhausted", "Tstate_simple_cytotoxic", "Tstate_simple_regulatory",
  # Simplified Myeloid states
  "Mstate_simple_TAM_immunosuppressive", "Mstate_simple_dendritic_cell"
)

# Keep only those that exist in our data
available_vars <- intersect(recommended_vars, names(numeric_features))
cat("\nRecommended vars found:", length(available_vars), "/", length(recommended_vars), "\n")

# If we need more, add high-variance features
if (length(available_vars) < max_vars && length(available_vars) < 15) {
  extra_vars <- setdiff(names(head(feature_var, max_vars + 5)), available_vars)
  n_extra <- min(max_vars - length(available_vars), length(extra_vars))
  if (n_extra > 0) {
    available_vars <- c(available_vars, extra_vars[1:n_extra])
  }
}

cat("\nSelected", length(available_vars), "variables for causal analysis:\n")
print(available_vars)

# Create final selected feature matrix
feature_matrix_selected <- feature_matrix_causal %>%
  select(Samples, Patients, Groups, all_of(available_vars))

write.csv(feature_matrix_selected, file.path(OUT_DIR, "feature_matrix_selected.csv"), row.names = FALSE)
saveRDS(feature_matrix_selected, file.path(OUT_DIR, "feature_matrix_selected.rds"))
cat("Saved: feature_matrix_selected.csv\n")

# ==============================================================================
# SUMMARY TABLE
# ==============================================================================

# Table: Feature summary statistics
feature_summary <- numeric_features %>%
  pivot_longer(everything(), names_to = "feature", values_to = "value") %>%
  group_by(feature) %>%
  summarize(
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    min = min(value, na.rm = TRUE),
    max = max(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(sd))

write.csv(feature_summary, "tables/tab_feature_summary.csv", row.names = FALSE)
cat("\nSaved: tables/tab_feature_summary.csv\n")

# ==============================================================================
# VISUALIZATION: FEATURE CORRELATIONS
# ==============================================================================
cat("\nGenerating feature correlation heatmap...\n")

library(pheatmap)

# Correlation matrix for causal variables
causal_numeric <- feature_matrix_selected %>% select(where(is.numeric))
if (ncol(causal_numeric) > 2) {
  cor_mat <- cor(causal_numeric, use = "pairwise.complete.obs")
  
  # Save heatmap
  png("figures/fig_feature_correlations.png", width = 10, height = 10, units = "in", res = 300)
  pheatmap(cor_mat,
           color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
           breaks = seq(-1, 1, length.out = 101),
           display_numbers = TRUE,
           number_format = "%.2f",
           fontsize_number = 7,
           main = "Feature Correlations (Causal Variables)")
  dev.off()
  
  cat("Saved: figures/fig_feature_correlations.png\n")
}

# ==============================================================================
# TISSUE-STRATIFIED ANALYSIS (for Groups comparison)
# ==============================================================================
cat("\nGenerating tissue-stratified feature comparison...\n")

# Compare features across tissue sites
tissue_comparison <- feature_matrix_selected %>%
  select(Groups, all_of(available_vars)) %>%
  pivot_longer(-Groups, names_to = "feature", values_to = "value") %>%
  group_by(Groups, feature) %>%
  summarize(
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = Groups, values_from = c(mean, sd))

write.csv(tissue_comparison, "tables/tab_features_by_tissue.csv", row.names = FALSE)
cat("Saved: tables/tab_features_by_tissue.csv\n")

# ==============================================================================
# NEXT STEPS
# ==============================================================================
cat("\n")
cat(strrep("=", 60), "\n")
cat("FEATURE ENGINEERING COMPLETE\n")
cat(strrep("=", 60), "\n")
cat("\nOutputs:\n")
cat("  - results/feature_matrix_full.csv       (all features)\n")
cat("  - results/feature_matrix_causal.csv     (simplified states)\n")
cat("  - results/feature_matrix_selected.csv   (top variables for DAG)\n")
cat("  - tables/tab_feature_summary.csv        (summary stats)\n")
cat("  - tables/tab_features_by_tissue.csv     (tissue comparison)\n")
cat("  - figures/fig_feature_correlations.png\n")
cat("\nNext step: Run 03_causal_dag.R on CRC\n")
cat("\n")
