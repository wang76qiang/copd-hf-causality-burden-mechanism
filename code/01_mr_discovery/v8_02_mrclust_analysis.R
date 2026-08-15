# ============================================================================
# Module 01 - MR discovery (LEGACY V8 script, kept for provenance)
# MR-Clust (Foley et al. 2021) clustering of Wald ratios for COPD -> HF.
# Input : data/derived/mr_discovery/discovery_harmonised_10snp_dat.xlsx
# Output: mrclust_results_* (archived in results/legacy/)
# ============================================================================

# --- repository paths -----------------------------------------------------------
# Run this script from the repository root, or set REPRO_ROOT to the repo root.
REPRO_ROOT <- Sys.getenv("REPRO_ROOT", unset = ".")
# GBD_RAW_DIR: folder with the raw IHME GBD 2021 csv downloads (see data/README.md)
GBD_RAW_DIR <- Sys.getenv("GBD_RAW_DIR", unset = file.path(REPRO_ROOT, "data", "raw", "gbd_2021"))
# V8_DIR: legacy V8 project tree (NOT shipped; only needed by v8_* legacy scripts)
V8_DIR <- Sys.getenv("V8_DIR", unset = "v8_legacy")
# 02_mrclust_analysis.R
# Run MR-Clust (Foley et al. 2021) on COPD -> Heart Failure summary statistics
# Requires: mrclust R package

# .libPaths(...)  # removed user-local library path; see r_packages.txt
library(mrclust)
library(TwoSampleMR)
library(tidyverse)
library(readxl)

data_file <- file.path(REPRO_ROOT, "data", "derived", "mr_discovery", "discovery_harmonised_10snp_dat.xlsx")
out_prefix <- file.path(V8_DIR, "R_scripts/results/mrclust_results")

# Read harmonized data
dat <- readxl::read_excel(data_file) %>%
  filter(mr_keep == TRUE)

# Prepare MR-Clust input (Wald ratios)
theta <- dat$beta.outcome / dat$beta.exposure
theta_se <- sqrt((dat$se.outcome^2 / dat$beta.exposure^2) +
                 (dat$beta.outcome^2 * dat$se.exposure^2 / dat$beta.exposure^4))

# Run MR-Clust
res <- mr_clust_em(
  theta = theta,
  theta_se = theta_se,
  bx = dat$beta.exposure,
  by = dat$beta.outcome,
  bxse = dat$se.exposure,
  byse = dat$se.outcome,
  obs_names = dat$SNP,
  plot_results = NULL  # disable plotting to avoid ggplot issues
)

# MR-Clust result extraction (summary() can fail on mixed data frames; extract manually)
clusters <- res$results$best$cluster
cluster_class <- res$results$best$cluster_class
cluster_mean <- res$results$best$cluster_mean
prob <- res$results$best$probability

cluster_summary <- data.frame(
  SNP = dat$SNP,
  cluster = clusters,
  cluster_class = cluster_class,
  cluster_mean = cluster_mean,
  probability = prob,
  wald_ratio = theta,
  stringsAsFactors = FALSE
)

# Tabulate clusters
cluster_table <- cluster_summary %>%
  group_by(cluster, cluster_class) %>%
  summarise(
    n = n(),
    mean_wald_ratio = mean(wald_ratio, na.rm = TRUE),
    snps = paste(SNP, collapse = ", "),
    .groups = "drop"
  )

print(cluster_table)

# Save results
saveRDS(res, file = paste0(out_prefix, "_mrclust_fit.rds"))
write_csv(cluster_summary, file = paste0(out_prefix, "_clusters.csv"))
write_csv(cluster_table, file = paste0(out_prefix, "_cluster_table.csv"))

# Create a compact JSON summary for the manuscript
summary_json <- list(
  n_snp = nrow(dat),
  n_clusters = length(unique(clusters)),
  bic_best = min(as.numeric(unlist(res$bic)), na.rm = TRUE),
  log_likelihood_best = max(as.numeric(unlist(res$log_likelihood)), na.rm = TRUE),
  clusters = jsonlite::fromJSON(jsonlite::toJSON(cluster_table))
)
jsonlite::write_json(summary_json, paste0(out_prefix, "_mrclust_summary.json"), pretty = TRUE, auto_unbox = TRUE)

cat("MR-Clust completed.\n")
