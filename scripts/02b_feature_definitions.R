# ==============================================================================
# 02b_feature_definitions.R
# Display all numeric features and how they were constructed
# ==============================================================================

library(tidyverse)

# Load selected feature matrix
feature_matrix_selected <- read.csv("results/feature_matrix_selected.csv")

# Identify numeric features
numeric_features <- feature_matrix_selected %>%
  select(where(is.numeric)) %>%
  names()

# ------------------------------------------------------------------
# Feature definitions (based on how features were created in script)
# ------------------------------------------------------------------

feature_definitions <- tribble(
  ~feature, ~description, ~formula,
  
  "prop_CD8_T",
  "Proportion of CD8+ T cells among all cells in sample",
  "count(CD8+ T) / total_cells",
  
  "prop_CD4_T",
  "Proportion of CD4+ T cells among all cells in sample",
  "count(CD4+ T) / total_cells",
  
  "prop_Epithelial_cells",
  "Proportion of epithelial (tumor) cells",
  "count(Epithelial cells) / total_cells",
  
  "prop_Fibroblast",
  "Proportion of fibroblast stromal cells",
  "count(Fibroblasts) / total_cells",
  
  "prop_NK",
  "Proportion of natural killer cells",
  "count(NK) / total_cells",
  
  "prop_B",
  "Proportion of B cells",
  "count(B) / total_cells",
  
  "prop_Endothelial_cells",
  "Proportion of endothelial cells",
  "count(Endothelial cells) / total_cells",
  
  "CD8_CD4_ratio",
  "Ratio of CD8 T cells to CD4 T cells",
  "(n_CD8 + 1) / (n_CD4 + 1)",
  
  "Tcell_Myeloid_ratio",
  "Ratio of T cells to myeloid lineage cells",
  "(n_CD4 + n_CD8 + 1) / (n_Myeloid + 1)",
  
  "Immune_Tumor_ratio",
  "Ratio of immune cells to tumor epithelial cells",
  "(n_Tcells + n_Myeloid + n_NK + n_B + 1) / (n_Tumor + 1)",
  
  "Tstate_simple_exhausted",
  "Fraction of T cells classified as exhausted (TIM3/CXCL13 signature)",
  "count(T10_CD8-HAVCR2 + T05_CD4-CXCL13) / total_T_cells",
  
  "Tstate_simple_cytotoxic",
  "Fraction of cytotoxic T cells (CX3CR1+ terminal effectors)",
  "count(T04_CD4-CX3CR1 + T09_CD8-CX3CR1) / total_T_cells",
  
  "Tstate_simple_regulatory",
  "Fraction of regulatory T cells (FOXP3+ Tregs)",
  "count(T03_CD4-FOXP3) / total_T_cells",
  
  "Mstate_simple_TAM_immunosuppressive",
  "Fraction of tumor-associated macrophages with immunosuppressive phenotype",
  "count(FOLR2+, FABP5+, complement macrophages) / total_myeloid_cells",
  
  "Mstate_simple_dendritic_cell",
  "Fraction of dendritic cells (cDC1, cDC2, mature DC)",
  "count(CLEC9A+, CD1C+, LAMP3+ DC) / total_myeloid_cells"
)

# ------------------------------------------------------------------
# Match definitions to features present in matrix
# ------------------------------------------------------------------

feature_table <- tibble(feature = numeric_features) %>%
  left_join(feature_definitions, by = "feature") %>%
  mutate(
    description = ifelse(is.na(description), "Derived proportion feature", description),
    formula = ifelse(is.na(formula), "Computed as sample-level proportion", formula)
  )

# ------------------------------------------------------------------
# Print table
# ------------------------------------------------------------------

print(feature_table)

# ------------------------------------------------------------------
# Save table
# ------------------------------------------------------------------

dir.create("tables", showWarnings = FALSE)

write.csv(
  feature_table,
  "tables/tab_feature_definitions.csv",
  row.names = FALSE
)

cat("\nSaved: tables/tab_feature_definitions.csv\n")