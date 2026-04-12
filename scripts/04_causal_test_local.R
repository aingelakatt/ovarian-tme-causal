# ==============================================================================
# 05_causal_test.R
# Local Test Run for Causal DAG Analysis
# 
# Purpose: Run 03_causal_dag.R logic on a subset of data to estimate
#          time and resources needed for full CRC run
#
# Strategy: 
#   - Use fewer variables (subset of features)
#   - Use fewer bootstrap iterations
#   - Track time and memory at each step
#   - Extrapolate to full run
#
# Runs on: Local RStudio (16GB RAM)
# Expected runtime: 3-10 minutes
# ==============================================================================

# --- Configuration ------------------------------------------------------------
PROJECT_DIR <- "."
INPUT_FILE <- file.path(PROJECT_DIR, "results/feature_matrix_selected.csv")
OUT_DIR <- file.path(PROJECT_DIR, "results/causal_test")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# TEST PARAMETERS
N_BOOTSTRAP_TEST <- 50      # 50 iterations (vs 1000 on CRC)
N_CORES_TEST <- 2           # 2 cores (vs 16 on CRC)
ALPHA <- 0.05
SEED <- 42

set.seed(SEED)

cat(strrep("=", 70), "\n")
cat("04_causal_test_local.R — Local Resource Estimation\n")
cat(strrep("=", 70), "\n\n")

# --- Load Libraries (timed) ---------------------------------------------------
cat("Loading libraries...\n")
t0 <- Sys.time()

library(tidyverse)
library(pcalg)
library(bnlearn)
library(parallel)

t_libs <- difftime(Sys.time(), t0, units = "secs")
cat("  Libraries loaded:", round(t_libs, 1), "sec\n\n")
# ==============================================================================
# LOAD AND SUBSET DATA
# ==============================================================================
cat("Loading data...\n")
t0 <- Sys.time()

feature_df <- read.csv(INPUT_FILE)
cat("  Full data:", nrow(feature_df), "samples x", ncol(feature_df), "columns\n")

# Extract numeric columns
meta_cols <- c("Samples", "Patients", "Groups")
numeric_cols <- setdiff(names(feature_df), meta_cols)
numeric_cols <- numeric_cols[sapply(feature_df[numeric_cols], is.numeric)]

cat("  Numeric features:", length(numeric_cols), "\n")
cat("  Variables:", paste(numeric_cols, collapse = ", "), "\n\n")

# --- Option: Subset variables for faster testing ------------------------------
# Uncomment to test with fewer variables (faster but less representative)
#
# SUBSET_VARS <- c(
#   "CD8_CD4_ratio", 
#   "Immune_Tumor_ratio",
#   "Tstate_simple_exhausted",
#   "Tstate_simple_cytotoxic",
#   "Mstate_simple_TAM_immunosuppressive"
# )
# numeric_cols <- intersect(SUBSET_VARS, numeric_cols)
# cat("  Using subset:", length(numeric_cols), "variables\n")

# Prepare matrix
causal_mat <- as.matrix(feature_df[, numeric_cols])
rownames(causal_mat) <- feature_df$Samples

# Remove constant columns
col_vars <- apply(causal_mat, 2, var, na.rm = TRUE)
constant_cols <- names(col_vars[col_vars < 1e-10])
if (length(constant_cols) > 0) {
  cat("  Removing constant columns:", paste(constant_cols, collapse = ", "), "\n")
  causal_mat <- causal_mat[, !colnames(causal_mat) %in% constant_cols]
}

# Impute missing
if (any(is.na(causal_mat))) {
  for (j in 1:ncol(causal_mat)) {
    causal_mat[is.na(causal_mat[,j]), j] <- mean(causal_mat[,j], na.rm = TRUE)
  }
  cat("  Imputed missing values\n")
}

# Scale
causal_mat_scaled <- scale(causal_mat)

n <- nrow(causal_mat_scaled)
p <- ncol(causal_mat_scaled)

t_load <- difftime(Sys.time(), t0, units = "secs")
cat("  Data prepared:", n, "samples x", p, "variables\n")
cat("  n/p ratio:", round(n/p, 2), "\n")
cat("  Load time:", round(t_load, 2), "sec\n\n")

# Prepare for algorithms
suffStat <- list(C = cor(causal_mat_scaled), n = n)
bn_data <- as.data.frame(causal_mat_scaled)

# ==============================================================================
# RUN ALGORITHMS WITH TIMING
# ==============================================================================
timing <- list()

cat(strrep("=", 70), "\n")
cat("ALGORITHM TIMING\n")
cat(strrep("=", 70), "\n\n")

# --- PC Algorithm -------------------------------------------------------------
cat("PC Algorithm...\n")
t0 <- Sys.time()

pc_fit <- pc(
  suffStat = suffStat,
  indepTest = gaussCItest,
  alpha = ALPHA,
  labels = colnames(causal_mat_scaled),
  verbose = FALSE
)

timing$pc <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
pc_amat <- as(pc_fit, "amat")
n_edges_pc <- sum(pc_amat == 1 & t(pc_amat) == 0) + sum(pc_amat == 1 & t(pc_amat) == 1) / 2

cat("  Time:", round(timing$pc, 2), "sec\n")
cat("  Edges:", n_edges_pc, "\n\n")

# --- FCI Algorithm ------------------------------------------------------------
cat("FCI Algorithm...\n")
t0 <- Sys.time()

fci_fit <- fci(
  suffStat = suffStat,
  indepTest = gaussCItest,
  alpha = ALPHA,
  labels = colnames(causal_mat_scaled),
  verbose = FALSE
)

timing$fci <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat("  Time:", round(timing$fci, 2), "sec\n\n")
# --- Hill-Climbing ------------------------------------------------------------
cat("Hill-Climbing...\n")
t0 <- Sys.time()

hc_fit <- hc(bn_data, score = "bic-g")

timing$hc <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat("  Time:", round(timing$hc, 2), "sec\n")
cat("  Edges:", narcs(hc_fit), "\n\n")
# --- TABU Search --------------------------------------------------------------
cat("TABU Search...\n")
t0 <- Sys.time()

tabu_fit <- tabu(bn_data, score = "bic-g")

timing$tabu <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat("  Time:", round(timing$tabu, 2), "sec\n")
cat("  Edges:", narcs(tabu_fit), "\n\n")
# --- Bootstrap ----------------------------------------------------------------
cat("Bootstrap (", N_BOOTSTRAP_TEST, " iterations, ", N_CORES_TEST, " cores)...\n", sep = "")
t0 <- Sys.time()

cl <- makeCluster(N_CORES_TEST)

boot_strength <- boot.strength(
  bn_data,
  R = N_BOOTSTRAP_TEST,
  algorithm = "hc",
  algorithm.args = list(score = "bic-g"),
  cpdag = FALSE,
  cluster = cl
)

stopCluster(cl)

timing$bootstrap <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
timing$per_bootstrap <- timing$bootstrap / N_BOOTSTRAP_TEST

cat("  Total time:", round(timing$bootstrap, 2), "sec\n")
cat("  Per iteration:", round(timing$per_bootstrap, 3), "sec\n")
cat("  Edges >=50%:", sum(boot_strength$strength >= 0.50), "\n\n")

# ==============================================================================
# RESOURCE ESTIMATION FOR CRC
# ==============================================================================
cat(strrep("=", 70), "\n")
cat("CRC RESOURCE ESTIMATES\n")
cat(strrep("=", 70), "\n\n")

# --- Time estimates -----------------------------------------------------------
# Non-bootstrap algorithms (scale linearly, roughly)
t_algorithms <- timing$pc + timing$fci + timing$hc + timing$tabu

# Bootstrap extrapolation
# With 16 cores vs 2 cores: ~8x speedup (not perfect due to overhead)
# Assume 75% parallel efficiency
parallel_efficiency <- 0.75
core_ratio <- 16 / N_CORES_TEST
speedup <- core_ratio * parallel_efficiency

t_boot_1000_16cores <- (timing$per_bootstrap * 1000) / speedup

# Total estimated time
t_total_estimated <- t_algorithms + t_boot_1000_16cores
t_total_with_buffer <- t_total_estimated * 1.5  # 50% safety buffer

cat("Time breakdown:\n")
cat("  PC + FCI + HC + TABU:", round(t_algorithms / 60, 1), "min\n")
cat("  Bootstrap (1000 iter, 16 cores):", round(t_boot_1000_16cores / 60, 1), "min\n")
cat("  Total estimated:", round(t_total_estimated / 60, 1), "min\n")
cat("  With 50% buffer:", round(t_total_with_buffer / 60, 1), "min\n\n")

# Round to SLURM-friendly values
slurm_hours <- ceiling(t_total_with_buffer / 3600)
if (slurm_hours < 1) slurm_hours <- 1
if (slurm_hours > 4) slurm_hours <- ceiling(slurm_hours / 2) * 2  # Round to even hours

cat("Recommended SLURM time:", slurm_hours, "hours\n\n")
# --- Memory estimates ---------------------------------------------------------
mem_data <- object.size(bn_data) / 1024^2  # MB
mem_results <- object.size(boot_strength) / 1024^2  # MB

# Each parallel worker needs a copy of data
# Plus overhead for R sessions
mem_per_worker <- mem_data * 1.5  # 50% overhead per worker
mem_16_workers <- mem_per_worker * 16
mem_total <- mem_16_workers + mem_results + 2000  # +2GB base R overhead

cat("Memory breakdown:\n")
cat("  Data size:", round(mem_data, 1), "MB\n")
cat("  Results size:", round(mem_results, 1), "MB\n")
cat("  Estimated for 16 workers:", round(mem_total / 1024, 1), "GB\n\n")

# Recommend SLURM memory
if (mem_total < 8000) {
  slurm_mem <- "16G"
} else if (mem_total < 16000) {
  slurm_mem <- "32G"
} else {
  slurm_mem <- "64G"
}

cat("Recommended SLURM memory:", slurm_mem, "\n\n")
# ==============================================================================
# SLURM SCRIPT RECOMMENDATION
# ==============================================================================
cat(strrep("=", 70), "\n")
cat("RECOMMENDED SLURM SETTINGS\n")
cat(strrep("=", 70), "\n\n")

cat("Copy this to slurm/run_causal_dag.sbatch:\n\n")

cat("```bash\n")
cat("#!/bin/bash\n")
cat("#SBATCH --job-name=causal_dag\n")
cat("#SBATCH --output=slurm/causal_dag_%j.out\n")
cat("#SBATCH --error=slurm/causal_dag_%j.err\n")
cat("#SBATCH --time=", sprintf("%02d", slurm_hours), ":00:00\n", sep = "")
cat("#SBATCH --nodes=1\n")
cat("#SBATCH --ntasks=1\n")
cat("#SBATCH --cpus-per-task=16\n")
cat("#SBATCH --mem=", slurm_mem, "\n", sep = "")
cat("#SBATCH --cluster=smp\n")
cat("#SBATCH --partition=smp\n")
cat("\n")
cat("# Set project directory\n")
cat("export PROJECT_DIR=/bgfs/your_lab/projects/ovarian-tme-causal\n")
cat("\n")
cat("# Load R module\n")
cat("module load r/4.3.0\n")
cat("\n")
cat("# Run script\n")
cat("cd $PROJECT_DIR\n")
cat("Rscript scripts/03_causal_dag.R\n")
cat("```\n\n")
# ==============================================================================
# PREVIEW RESULTS
# ==============================================================================
cat(strrep("=", 70), "\n")
cat("PREVIEW: TOP EDGES\n")
cat(strrep("=", 70), "\n\n")

# Top bootstrap edges
cat("Bootstrap edges (>=30% support in", N_BOOTSTRAP_TEST, "iterations):\n")
top_edges <- boot_strength %>%
  filter(strength >= 0.30) %>%
  arrange(desc(strength))

if (nrow(top_edges) > 0) {
  for (i in 1:min(nrow(top_edges), 10)) {
    cat("  ", top_edges$from[i], " → ", top_edges$to[i], 
        " (", round(top_edges$strength[i] * 100, 0), "%)\n", sep = "")
  }
} else {
  cat("  No edges with >=30% support\n")
}

cat("\n")
cat("Note: With only", N_BOOTSTRAP_TEST, "iterations, these estimates are rough.\n")
cat("Full 1000-iteration run on CRC will give stable estimates.\n\n")
# ==============================================================================
# SUMMARY
# ==============================================================================
cat(strrep("=", 70), "\n")
cat("SUMMARY\n")
cat(strrep("=", 70), "\n\n")

cat("Data:\n")
cat("  Samples:", n, "\n")
cat("  Variables:", p, "\n")
cat("  n/p ratio:", round(n/p, 2), "\n\n")

cat("Local test (", N_BOOTSTRAP_TEST, " bootstraps, ", N_CORES_TEST, " cores):\n", sep = "")
cat("  Total time:", round(sum(unlist(timing)), 1), "sec\n\n")

cat("CRC projection (1000 bootstraps, 16 cores):\n")
cat("  Estimated time:", round(t_total_estimated / 60, 1), "min\n")
cat("  SLURM allocation:", slurm_hours, "hours,", slurm_mem, "RAM\n\n")

cat("Preliminary findings:\n")
cat("  PC edges:", n_edges_pc, "\n")
cat("  HC edges:", narcs(hc_fit), "\n")
cat("  Bootstrap edges >=50%:", sum(boot_strength$strength >= 0.50), "\n\n")

# Save timing results
timing_df <- data.frame(
  step = c("PC", "FCI", "HC", "TABU", "Bootstrap_total", "Bootstrap_per_iter"),
  seconds = c(timing$pc, timing$fci, timing$hc, timing$tabu, 
              timing$bootstrap, timing$per_bootstrap),
  notes = c("", "", "", "", 
            paste0(N_BOOTSTRAP_TEST, " iterations"),
            "extrapolate to 1000")
)
write.csv(timing_df, file.path(OUT_DIR, "timing_results.csv"), row.names = FALSE)

cat("Timing saved to:", file.path(OUT_DIR, "timing_results.csv"), "\n")
cat("\n")
cat("Ready for CRC submission!\n")
cat(strrep("=", 70), "\n")






