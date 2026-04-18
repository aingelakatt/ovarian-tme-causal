# ==============================================================================
# 03_causal_dag.R
# Ovarian TME Causal Inference Pipeline — Phase 3
#
# Purpose: Construct causal DAGs from feature matrix
#   - Preliminary graphs (lenient, exploratory)
#   - Final graphs (stringent, bootstrap-validated)
#
# Input:  feature_matrix_selected.csv from 02_feature_engineering.R
# Output: DAG structures at multiple stringency levels, visualizations
# Running on select few groups for analysis 
# Runs on: local
# ==============================================================================
# --- Load Required Libraries --------------------------------------------------
library(tidyverse)
library(pcalg)
library(bnlearn)
library(parallel)
library(ggplot2)
library(pheatmap)

# --- Configuration ------------------------------------------------------------
# Detect environment: CRC (Slurm) vs Local
if (Sys.getenv("SLURM_JOB_ID") != "") {
  PROJECT_DIR <- Sys.getenv("PROJECT_DIR")
  if (PROJECT_DIR == "") {
    stop("ERROR: PROJECT_DIR environment variable not set.\n",
         "Add this to your sbatch script:\n",
         "  export PROJECT_DIR=/bgfs/your_lab/projects/ovarian-tme-causal")
  }
  N_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "16"))
  N_BOOTSTRAP <- 1000
  IS_CRC <- TRUE
} else {
  PROJECT_DIR <- "."
  N_CORES <- max(1, min(8, parallel::detectCores() - 1))
  N_BOOTSTRAP <- 1000
  IS_CRC <- FALSE
}

# File paths
INPUT_FILE <- file.path(PROJECT_DIR, "results/feature_matrix_causal.csv")
RUN_LABEL <- "drop_cytotoxic"
OUT_DIR <- file.path(PROJECT_DIR, "results", "sensitivity", RUN_LABEL)
FIG_DIR <- file.path(OUT_DIR, "figures")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# Algorithm parameters
ALPHA_LENIENT <- 0.10
ALPHA_STRINGENT <- 0.01
BOOT_THRESHOLD_LOW <- 0.30
BOOT_THRESHOLD_HIGH <- 0.50

SEED <- 42
set.seed(SEED)

# --- Logging ------------------------------------------------------------------
cat(strrep("=", 70), "\n")
cat("03_causal_dag.R — Causal Graph Construction -- Sensitivity analysis\n")
cat(strrep("=", 70), "\n")
cat("Environment:", ifelse(IS_CRC, "CRC (Slurm)", "Local"), "\n")
cat("Input:", INPUT_FILE, "\n")
cat("Output:", OUT_DIR, "\n")
cat("Figures:", FIG_DIR, "\n")
cat("Cores:", N_CORES, "\n")
cat("Bootstrap iterations:", N_BOOTSTRAP, "\n")
cat("Alpha (lenient):", ALPHA_LENIENT, "\n")
cat("Alpha (stringent):", ALPHA_STRINGENT, "\n")
cat("Bootstrap thresholds:", BOOT_THRESHOLD_LOW, "/", BOOT_THRESHOLD_HIGH, "\n\n")

# ==============================================================================
# LOAD AND PREPARE DATA
# ==============================================================================
cat("Loading feature matrix...\n")

if (!file.exists(INPUT_FILE)) {
  stop("Feature matrix not found: ", INPUT_FILE, "\n",
       "Run 02_feature_engineering.R first.")
}

feature_df <- read.csv(INPUT_FILE)
cat("Loaded:", nrow(feature_df), "samples x", ncol(feature_df), "columns\n")

# Keep Groups as an upstream context node
feature_df$Groups <- factor(
  feature_df$Groups,
  levels = c("Primary Tumor", "Metastatic Tumor", "Lymph Node", "Ascites", "PBMC")
)

# Revised biologically layered, lower-redundancy variable set
selected_vars <- c(
  "prop_Epithelial_cells",
  "prop_Fibroblast",
  "Mstate_simple_TAM_immunosuppressive", #LOO sensitivity 1
  "Mstate_simple_dendritic_cell",
  "Tstate_simple_regulatory",
  "Tstate_simple_exhausted",
  # LOO: sensitivity 2: "Tstate_simple_cytotoxic",
  "prop_NK"
)

missing_vars <- setdiff(selected_vars, names(feature_df))
if (length(missing_vars) > 0) {
  stop("Missing required variables: ", paste(missing_vars, collapse = ", "))
}

feature_df <- feature_df[, c("Samples", "Patients", selected_vars)]

cat("Using revised variable set:\n")
for (v in selected_vars) cat("  -", v, "\n")


meta_cols <- c("Samples", "Patients")
numeric_cols <- setdiff(names(feature_df), meta_cols)

causal_mat <- as.matrix(feature_df[, numeric_cols])
rownames(causal_mat) <- feature_df$Samples

# Remove constant columns
col_vars <- apply(causal_mat, 2, var, na.rm = TRUE)
constant_cols <- names(col_vars[col_vars < 1e-10])
if (length(constant_cols) > 0) {
  cat("\nRemoving constant columns (zero variance):\n")
  for (v in constant_cols) cat("  -", v, "\n")
  causal_mat <- causal_mat[, !colnames(causal_mat) %in% constant_cols, drop = FALSE]
}

# Impute missing values
if (any(is.na(causal_mat))) {
  n_missing <- sum(is.na(causal_mat))
  cat("\nImputing", n_missing, "missing values with column means\n")
  for (j in seq_len(ncol(causal_mat))) {
    na_idx <- is.na(causal_mat[, j])
    if (any(na_idx)) {
      causal_mat[na_idx, j] <- mean(causal_mat[, j], na.rm = TRUE)
    }
  }
}

# Standardize continuous variables
causal_mat_scaled <- scale(causal_mat)

n <- nrow(causal_mat_scaled)
p <- ncol(causal_mat_scaled)

cat("\nFinal data matrix:", n, "samples x", p, "variables\n")
cat("n/p ratio:", round(n / p, 2), "(recommend > 5 for stability)\n\n")

if (n / p < 3) {
  cat("\n*** WARNING: n/p < 3 — results will be highly unstable! ***\n")
  cat("Consider reducing variables or this analysis may not be reliable.\n\n")
} else if (n / p < 5) {
  cat("\n* Note: n/p < 5 — interpret results with caution.\n\n")
}

# Prepare sufficient statistics for pcalg
suffStat <- list(
  C = cor(causal_mat_scaled),
  n = n
)

# Prepare data frame for bnlearn
bn_data <- as.data.frame(causal_mat_scaled)
# ------------------------------------------------------------------------------
# Biological constraints
# ------------------------------------------------------------------------------

# Groups_num is upstream context and should not have incoming edges
blacklist <- NULL

# Optional biologically motivated whitelist
whitelist <- NULL

# Helper to extract PC edges
extract_pc_edges <- function(pc_obj, method_label) {
  amat <- as(pc_obj, "amat")
  edge_list <- vector("list", length = nrow(amat) * ncol(amat))
  idx <- 1
  
  for (i in seq_len(nrow(amat))) {
    for (j in seq_len(ncol(amat))) {
      if (amat[i, j] == 1) {
        if (amat[j, i] == 0) {
          edge_list[[idx]] <- data.frame(
            from = rownames(amat)[i],
            to = colnames(amat)[j],
            type = "directed",
            method = method_label,
            stringsAsFactors = FALSE
          )
          idx <- idx + 1
        } else if (i < j) {
          edge_list[[idx]] <- data.frame(
            from = rownames(amat)[i],
            to = colnames(amat)[j],
            type = "undirected",
            method = method_label,
            stringsAsFactors = FALSE
          )
          idx <- idx + 1
        }
      }
    }
  }
  
  edge_list <- edge_list[seq_len(idx - 1)]
  if (length(edge_list) == 0) {
    return(data.frame(from = character(), to = character(), type = character(),
                      method = character(), stringsAsFactors = FALSE))
  }
  bind_rows(edge_list)
}

# Helper to extract bnlearn arcs safely
extract_bn_edges <- function(bn_obj, method_label) {
  if (narcs(bn_obj) == 0) {
    return(data.frame(from = character(), to = character(), type = character(),
                      method = character(), stringsAsFactors = FALSE))
  }
  data.frame(
    from = arcs(bn_obj)[, 1],
    to = arcs(bn_obj)[, 2],
    type = "directed",
    method = method_label,
    stringsAsFactors = FALSE
  )
}

# ==============================================================================
# PRELIMINARY GRAPHS (LENIENT — EXPLORATORY)
# ==============================================================================
cat(strrep("=", 70), "\n")
cat("PHASE 1: PRELIMINARY GRAPHS (alpha =", ALPHA_LENIENT, ")\n")
cat(strrep("=", 70), "\n")
cat("\nThese graphs use lenient thresholds to reveal POTENTIAL causal\n")
cat("relationships. Many edges may be false positives.\n\n")

# --- PC Algorithm (Lenient) ---------------------------------------------------
cat(strrep("-", 50), "\n")
cat("PC Algorithm (constraint-based)\n")
cat(strrep("-", 50), "\n")

pc_lenient <- pc(
  suffStat = suffStat,
  indepTest = gaussCItest,
  alpha = ALPHA_LENIENT,
  labels = colnames(causal_mat_scaled),
  verbose = FALSE
)

pc_amat_lenient <- as(pc_lenient, "amat")
n_directed_len <- sum(pc_amat_lenient == 1 & t(pc_amat_lenient) == 0)
n_undirected_len <- sum(pc_amat_lenient == 1 & t(pc_amat_lenient) == 1) / 2

cat("  Directed edges:", n_directed_len, "\n")
cat("  Undirected edges:", n_undirected_len, "\n")
cat("  Total edges:", n_directed_len + n_undirected_len, "\n")

edges_pc_lenient <- extract_pc_edges(pc_lenient, "PC_lenient")

# --- Hill-Climbing (Lenient) --------------------------------------------------
cat(strrep("-", 50), "\n")
cat("Hill-Climbing (score-based)\n")
cat(strrep("-", 50), "\n")

hc_lenient <- hc(bn_data, score = "bic-g", blacklist = blacklist, whitelist = whitelist)

cat("  Edges:", narcs(hc_lenient), "\n")

edges_hc_lenient <- extract_bn_edges(hc_lenient, "HC_lenient")

# --- TABU Search (Lenient) ----------------------------------------------------
cat(strrep("-", 50), "\n")
cat("TABU Search (score-based, more thorough)\n")
cat(strrep("-", 50), "\n")

tabu_lenient <- tabu(bn_data, score = "bic-g", blacklist = blacklist, whitelist = whitelist)

cat("  Edges:", narcs(tabu_lenient), "\n")

edges_tabu_lenient <- extract_bn_edges(tabu_lenient, "TABU_lenient")

# --- FCI Algorithm (Lenient) --------------------------------------------------
cat(strrep("-", 50), "\n")
cat("FCI Algorithm (handles latent confounders)\n")
cat(strrep("-", 50), "\n")

fci_lenient <- fci(
  suffStat = suffStat,
  indepTest = gaussCItest,
  alpha = ALPHA_LENIENT,
  labels = colnames(causal_mat_scaled),
  verbose = FALSE
)

cat("  FCI complete (outputs PAG with special edge types)\n")

# --- Visualize Preliminary Graphs ---------------------------------------------
cat("\n")
cat(strrep("-", 50), "\n")
cat("Visualizing preliminary graphs\n")
cat(strrep("-", 50), "\n")

pdf(file.path(FIG_DIR, "fig_01_pc_preliminary.pdf"), width = 12, height = 10)
plot(pc_lenient, main = paste0("PC Algorithm (Preliminary, alpha = ", ALPHA_LENIENT, ")"))
dev.off()
cat("Saved: figures/fig_01_pc_preliminary.pdf\n")

if (requireNamespace("Rgraphviz", quietly = TRUE)) {
  pdf(file.path(FIG_DIR, "fig_02_hc_preliminary.pdf"), width = 12, height = 10)
  graphviz.plot(hc_lenient,
                main = paste0("Hill-Climbing (Preliminary, ", narcs(hc_lenient), " edges)"))
  dev.off()
  cat("Saved: figures/fig_02_hc_preliminary.pdf\n")
  
  pdf(file.path(FIG_DIR, "fig_03_tabu_preliminary.pdf"), width = 12, height = 10)
  graphviz.plot(tabu_lenient,
                main = paste0("TABU Search (Preliminary, ", narcs(tabu_lenient), " edges)"))
  dev.off()
  cat("Saved: figures/fig_03_tabu_preliminary.pdf\n")
} else {
  cat("Note: Install Rgraphviz for graph plots: BiocManager::install('Rgraphviz')\n")
  write.csv(edges_hc_lenient, file.path(OUT_DIR, "edges_hc_preliminary.csv"), row.names = FALSE)
  write.csv(edges_tabu_lenient, file.path(OUT_DIR, "edges_tabu_preliminary.csv"), row.names = FALSE)
  cat("Saved edge lists as CSV files.\n")
}

pdf(file.path(FIG_DIR, "fig_04_fci_preliminary.pdf"), width = 12, height = 10)
plot(fci_lenient, main = paste0("FCI Algorithm (Preliminary, alpha = ", ALPHA_LENIENT, ")"))
dev.off()
cat("Saved: figures/fig_04_fci_preliminary.pdf\n")

# --- Preliminary Consensus ----------------------------------------------------
cat("\n")
cat(strrep("-", 50), "\n")
cat("Preliminary edge consensus\n")
cat(strrep("-", 50), "\n")

all_edges_lenient <- bind_rows(
  edges_pc_lenient,
  edges_hc_lenient,
  edges_tabu_lenient
)

consensus_lenient <- all_edges_lenient %>%
  group_by(from, to) %>%
  summarize(
    n_methods = n_distinct(method),
    methods = paste(sort(unique(method)), collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(desc(n_methods), from, to)

write.csv(consensus_lenient, file.path(OUT_DIR, "preliminary_edge_consensus.csv"), row.names = FALSE)

cat("\nEdges found by multiple methods (preliminary):\n")
multi_method <- filter(consensus_lenient, n_methods >= 2)
if (nrow(multi_method) > 0) {
  print(multi_method)
} else {
  cat("  No edges found by 2+ methods.\n")
}

# --- Save Preliminary Results -------------------------------------------------
saveRDS(pc_lenient, file.path(OUT_DIR, "pc_lenient.rds"))
saveRDS(hc_lenient, file.path(OUT_DIR, "hc_lenient.rds"))
saveRDS(tabu_lenient, file.path(OUT_DIR, "tabu_lenient.rds"))
saveRDS(fci_lenient, file.path(OUT_DIR, "fci_lenient.rds"))

cat("\nPreliminary results saved.\n")

# ==============================================================================
# STRINGENT GRAPHS (CONSERVATIVE)
# ==============================================================================
cat("\n")
cat(strrep("=", 70), "\n")
cat("PHASE 2: STRINGENT GRAPHS (alpha =", ALPHA_STRINGENT, ")\n")
cat(strrep("=", 70), "\n")
cat("\nThese graphs use stringent thresholds. Fewer edges, but higher\n")
cat("confidence that remaining edges represent true causal relationships.\n\n")

# --- PC Algorithm (Stringent) -------------------------------------------------
cat(strrep("-", 50), "\n")
cat("PC Algorithm (stringent)\n")
cat(strrep("-", 50), "\n")

pc_stringent <- pc(
  suffStat = suffStat,
  indepTest = gaussCItest,
  alpha = ALPHA_STRINGENT,
  labels = colnames(causal_mat_scaled),
  verbose = FALSE
)

pc_amat_stringent <- as(pc_stringent, "amat")
n_directed_str <- sum(pc_amat_stringent == 1 & t(pc_amat_stringent) == 0)
n_undirected_str <- sum(pc_amat_stringent == 1 & t(pc_amat_stringent) == 1) / 2

cat("  Directed edges:", n_directed_str, "\n")
cat("  Undirected edges:", n_undirected_str, "\n")
cat("  Total edges:", n_directed_str + n_undirected_str, "\n")
cat("  (Lenient had", n_directed_len + n_undirected_len, "edges)\n")

edges_pc_stringent <- extract_pc_edges(pc_stringent, "PC_stringent")

# --- FCI Algorithm (Stringent) ------------------------------------------------
cat(strrep("-", 50), "\n")
cat("FCI Algorithm (stringent)\n")
cat(strrep("-", 50), "\n")

fci_stringent <- fci(
  suffStat = suffStat,
  indepTest = gaussCItest,
  alpha = ALPHA_STRINGENT,
  labels = colnames(causal_mat_scaled),
  verbose = FALSE
)

cat("  FCI complete\n")

# --- Visualize Stringent Graphs -----------------------------------------------
pdf(file.path(FIG_DIR, "fig_05_pc_stringent.pdf"), width = 12, height = 10)
plot(pc_stringent, main = paste0("PC Algorithm (Stringent, alpha = ", ALPHA_STRINGENT, ")"))
dev.off()
cat("Saved: figures/fig_05_pc_stringent.pdf\n")

pdf(file.path(FIG_DIR, "fig_06_fci_stringent.pdf"), width = 12, height = 10)
plot(fci_stringent, main = paste0("FCI Algorithm (Stringent, alpha = ", ALPHA_STRINGENT, ")"))
dev.off()
cat("Saved: figures/fig_06_fci_stringent.pdf\n")

# --- Save Stringent Results ---------------------------------------------------
saveRDS(pc_stringent, file.path(OUT_DIR, "pc_stringent.rds"))
saveRDS(fci_stringent, file.path(OUT_DIR, "fci_stringent.rds"))

cat("\nStringent results saved.\n")

# ==============================================================================
# BOOTSTRAP STABILITY ANALYSIS
# ==============================================================================
cat("\n")
cat(strrep("=", 70), "\n")
cat("PHASE 3: BOOTSTRAP STABILITY ANALYSIS (n =", N_BOOTSTRAP, ")\n")
cat(strrep("=", 70), "\n")
cat("\nBootstrap resamples the data", N_BOOTSTRAP, "times, fitting a DAG each time.\n")
cat("Edges that appear consistently across bootstraps are stable/robust.\n\n")

cat("Starting parallel cluster with", N_CORES, "cores...\n")
cl <- makeCluster(N_CORES)

cat("Running bootstrap (this may take a while)...\n")
boot_start <- Sys.time()

boot_strength <- boot.strength(
  bn_data,
  R = N_BOOTSTRAP,
  algorithm = "hc",
  algorithm.args = list(score = "bic-g", blacklist = blacklist, whitelist = whitelist),
  cpdag = FALSE,
  cluster = cl
)

boot_end <- Sys.time()
boot_time <- difftime(boot_end, boot_start, units = "mins")

stopCluster(cl)

cat("Bootstrap complete in", round(boot_time, 1), "minutes.\n")

saveRDS(boot_strength, file.path(OUT_DIR, "bootstrap_strength.rds"))

# --- Analyze Bootstrap at Multiple Thresholds ---------------------------------
cat("\n")
cat(strrep("-", 50), "\n")
cat("Bootstrap edge strength analysis\n")
cat(strrep("-", 50), "\n")

edges_boot_low <- boot_strength %>%
  filter(strength >= BOOT_THRESHOLD_LOW) %>%
  arrange(desc(strength))

cat("\nEdges with >=", BOOT_THRESHOLD_LOW * 100, "% bootstrap support:\n")
cat("  Count:", nrow(edges_boot_low), "\n")
if (nrow(edges_boot_low) > 0) {
  print(head(edges_boot_low, 15))
}

edges_boot_high <- boot_strength %>%
  filter(strength >= BOOT_THRESHOLD_HIGH) %>%
  arrange(desc(strength))

cat("\nEdges with >=", BOOT_THRESHOLD_HIGH * 100, "% bootstrap support (FINAL):\n")
cat("  Count:", nrow(edges_boot_high), "\n")
if (nrow(edges_boot_high) > 0) {
  print(edges_boot_high)
}

write.csv(edges_boot_low, file.path(OUT_DIR, "bootstrap_edges_30pct.csv"), row.names = FALSE)
write.csv(edges_boot_high, file.path(OUT_DIR, "bootstrap_edges_50pct.csv"), row.names = FALSE)

# --- Create Averaged Networks -------------------------------------------------
cat("\n")
cat(strrep("-", 50), "\n")
cat("Creating averaged networks\n")
cat(strrep("-", 50), "\n")

avg_net_lenient <- averaged.network(boot_strength, threshold = BOOT_THRESHOLD_LOW)
cat("Averaged network (", BOOT_THRESHOLD_LOW * 100, "% threshold):", narcs(avg_net_lenient), "edges\n")

avg_net_final <- averaged.network(boot_strength, threshold = BOOT_THRESHOLD_HIGH)
cat("Averaged network (", BOOT_THRESHOLD_HIGH * 100, "% threshold):", narcs(avg_net_final), "edges\n")

saveRDS(avg_net_lenient, file.path(OUT_DIR, "averaged_network_30pct.rds"))
saveRDS(avg_net_final, file.path(OUT_DIR, "averaged_network_50pct.rds"))

# --- Visualize Bootstrap Results ----------------------------------------------
cat("\n")
cat(strrep("-", 50), "\n")
cat("Generating bootstrap visualizations\n")
cat(strrep("-", 50), "\n")

strength_mat <- matrix(
  0,
  p,
  p,
  dimnames = list(colnames(causal_mat_scaled), colnames(causal_mat_scaled))
)

for (i in seq_len(nrow(boot_strength))) {
  from_var <- boot_strength$from[i]
  to_var <- boot_strength$to[i]
  if (from_var %in% rownames(strength_mat) && to_var %in% colnames(strength_mat)) {
    strength_mat[from_var, to_var] <- boot_strength$strength[i]
  }
}

png(file.path(FIG_DIR, "fig_07_bootstrap_strength_heatmap.png"),
    width = 12, height = 10, units = "in", res = 300)
pheatmap(
  strength_mat,
  color = colorRampPalette(c("white", "#FFF7BC", "#FEC44F", "#D95F0E", "#993404"))(100),
  main = paste0("Bootstrap Edge Strength (n = ", N_BOOTSTRAP, ")"),
  fontsize = 10,
  display_numbers = TRUE,
  number_format = "%.2f",
  fontsize_number = 8,
  cluster_rows = FALSE,
  cluster_cols = FALSE
)
dev.off()
cat("Saved: figures/fig_07_bootstrap_strength_heatmap.png\n")

if (nrow(edges_boot_low) > 0) {
  top_edges <- head(edges_boot_low, 20)
  
  p_stability <- ggplot(top_edges, aes(x = reorder(paste(from, "→", to), strength), y = strength)) +
    geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
    geom_hline(yintercept = BOOT_THRESHOLD_HIGH, linetype = "dashed", color = "red", linewidth = 1) +
    geom_hline(yintercept = BOOT_THRESHOLD_LOW, linetype = "dotted", color = "orange", linewidth = 0.8) +
    coord_flip() +
    labs(
      title = "Bootstrap Edge Stability",
      subtitle = paste0("n = ", N_BOOTSTRAP, " | Red = ", BOOT_THRESHOLD_HIGH * 100,
                        "% | Orange = ", BOOT_THRESHOLD_LOW * 100, "%"),
      x = "Edge",
      y = "Bootstrap Support"
    ) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
    theme_minimal(base_size = 12) +
    theme(axis.text.y = element_text(size = 10))
  
  ggsave(file.path(FIG_DIR, "fig_08_bootstrap_stability_barplot.png"),
         p_stability, width = 10, height = 8, dpi = 300)
  cat("Saved: figures/fig_08_bootstrap_stability_barplot.png\n")
}

if (requireNamespace("Rgraphviz", quietly = TRUE)) {
  if (narcs(avg_net_final) > 0) {
    pdf(file.path(FIG_DIR, "fig_09_final_causal_graph.pdf"), width = 14, height = 12)
    graphviz.plot(
      avg_net_final,
      main = paste0("Final Causal Graph (Bootstrap >= ",
                    BOOT_THRESHOLD_HIGH * 100, "%, n = ", N_BOOTSTRAP, ")"),
      shape = "ellipse",
      fontsize = 14
    )
    dev.off()
    cat("Saved: figures/fig_09_final_causal_graph.pdf\n")
  } else {
    cat("No edges in final graph at", BOOT_THRESHOLD_HIGH * 100, "% threshold.\n")
  }
  
  if (narcs(avg_net_lenient) > 0) {
    pdf(file.path(FIG_DIR, "fig_10_exploratory_causal_graph.pdf"), width = 14, height = 12)
    graphviz.plot(
      avg_net_lenient,
      main = paste0("Exploratory Causal Graph (Bootstrap >= ",
                    BOOT_THRESHOLD_LOW * 100, "%, n = ", N_BOOTSTRAP, ")"),
      shape = "ellipse",
      fontsize = 14
    )
    dev.off()
    cat("Saved: figures/fig_10_exploratory_causal_graph.pdf\n")
  }
} else {
  cat("Install Rgraphviz for network plots: BiocManager::install('Rgraphviz')\n")
}

# ==============================================================================
# FINAL CONSENSUS ACROSS ALL METHODS
# ==============================================================================
cat("\n")
cat(strrep("=", 70), "\n")
cat("FINAL CONSENSUS\n")
cat(strrep("=", 70), "\n")

edges_final <- bind_rows(
  edges_pc_stringent %>% mutate(method = "PC_stringent"),
  edges_boot_high %>%
    select(from, to) %>%
    mutate(type = "directed", method = "Bootstrap_50pct")
)

final_consensus <- edges_final %>%
  group_by(from, to) %>%
  summarize(
    n_methods = n_distinct(method),
    methods = paste(sort(unique(method)), collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(desc(n_methods), from, to)

write.csv(final_consensus, file.path(OUT_DIR, "final_edge_consensus.csv"), row.names = FALSE)

cat("\nHigh-confidence edges (found by stringent PC AND bootstrap):\n")
high_conf <- filter(final_consensus, n_methods >= 2)
if (nrow(high_conf) > 0) {
  for (i in seq_len(nrow(high_conf))) {
    cat("  ", high_conf$from[i], " → ", high_conf$to[i], "\n", sep = "")
  }
} else {
  cat("  No edges found by both methods.\n")
  cat("  This may indicate:\n")
  cat("    - Sample size too small for confident causal inference\n")
  cat("    - True causal relationships are weak\n")
  cat("    - Consider using lenient thresholds for hypothesis generation\n")
}

# ==============================================================================
# SUMMARY
# ==============================================================================
cat("\n")
cat(strrep("=", 70), "\n")
cat("CAUSAL DAG CONSTRUCTION COMPLETE\n")
cat(strrep("=", 70), "\n")

cat("\nSUMMARY OF RESULTS:\n")
cat(strrep("-", 40), "\n")
cat("Preliminary graphs (alpha =", ALPHA_LENIENT, "):\n")
cat("  PC:", n_directed_len + n_undirected_len, "edges\n")
cat("  Hill-climbing:", narcs(hc_lenient), "edges\n")
cat("  TABU:", narcs(tabu_lenient), "edges\n")
cat("\nStringent graphs (alpha =", ALPHA_STRINGENT, "):\n")
cat("  PC:", n_directed_str + n_undirected_str, "edges\n")
cat("\nBootstrap analysis (n =", N_BOOTSTRAP, "):\n")
cat("  Edges >=", BOOT_THRESHOLD_LOW * 100, "%:", nrow(edges_boot_low), "\n")
cat("  Edges >=", BOOT_THRESHOLD_HIGH * 100, "%:", nrow(edges_boot_high), "\n")
cat("  Final averaged network:", narcs(avg_net_final), "edges\n")

cat("\nOUTPUTS SAVED TO:", OUT_DIR, "\n")
cat(strrep("-", 40), "\n")
cat("Preliminary:\n")
cat("  - pc_lenient.rds, hc_lenient.rds, tabu_lenient.rds, fci_lenient.rds\n")
cat("  - fig_01-04_*_preliminary.pdf\n")
cat("  - preliminary_edge_consensus.csv\n")
cat("Stringent:\n")
cat("  - pc_stringent.rds, fci_stringent.rds\n")
cat("  - fig_05-06_*_stringent.pdf\n")
cat("Bootstrap:\n")
cat("  - bootstrap_strength.rds\n")
cat("  - averaged_network_30pct.rds, averaged_network_50pct.rds\n")
cat("  - bootstrap_edges_30pct.csv, bootstrap_edges_50pct.csv\n")
cat("  - fig_07_bootstrap_strength_heatmap.png\n")
cat("  - fig_08_bootstrap_stability_barplot.png\n")
cat("  - fig_09_final_causal_graph.pdf\n")
cat("  - fig_10_exploratory_causal_graph.pdf\n")
cat("Final:\n")
cat("  - final_edge_consensus.csv\n")

cat("\n")
cat("NEXT STEP: Run 04_ci_test.R to test specific causal hypotheses\n")
cat("\n")

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "session_info_03_causal_dag.txt"))
cat("Session info saved.\n")
cat("Done!\n")
