# ============================================================================
# Module 05 - GBD (LEGACY V8 script): frequentist/BAPC APC projections.
# Legacy BAPC projection run on the V8 age-stratified inputs; the frequentist
# APC comparator quoted in the manuscript (+73% to 6.25 M by 2050) originates
# from this pipeline. Inputs are NOT shipped (V8 data tree); kept for
# provenance. The primary projection is t1_bapc_run.R.
# ============================================================================

# --- repository paths -----------------------------------------------------------
# Run this script from the repository root, or set REPRO_ROOT to the repo root.
REPRO_ROOT <- Sys.getenv("REPRO_ROOT", unset = ".")
# GBD_RAW_DIR: folder with the raw IHME GBD 2021 csv downloads (see data/README.md)
GBD_RAW_DIR <- Sys.getenv("GBD_RAW_DIR", unset = file.path(REPRO_ROOT, "data", "raw", "gbd_2021"))
# V8_DIR: legacy V8 project tree (NOT shipped; only needed by v8_* legacy scripts)
V8_DIR <- Sys.getenv("V8_DIR", unset = "v8_legacy")
#!/usr/bin/env Rscript
# .libPaths(...)  # removed user-local library path; see r_packages.txt
# 08_bapc_projection.R
# Bayesian age-period-cohort projections of GBD burden measures using BAPC/INLA.

suppressPackageStartupMessages({
  library(BAPC)
  library(INLA)
})

base_dir <- V8_DIR
data_dir <- file.path(base_dir, "data/bapc")
out_dir  <- file.path(base_dir, "results/bapc")
plot_dir <- file.path(out_dir, "plots")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# Find all *_counts.csv files and match with *_population.csv
files <- list.files(data_dir, pattern = "_counts\\.csv$", full.names = TRUE)

run_one <- function(counts_path) {
  slug <- sub("_counts\\.csv$", "", basename(counts_path))
  pop_path <- file.path(data_dir, paste0(slug, "_population.csv"))
  if (!file.exists(pop_path)) {
    message("Population file missing for ", slug)
    return(NULL)
  }

  counts <- read.csv(counts_path, row.names = 1, check.names = FALSE)
  pop    <- read.csv(pop_path, row.names = 1, check.names = FALSE)

  # Ensure identical dimensions and order
  common_years <- intersect(rownames(counts), rownames(pop))
  common_ages  <- intersect(colnames(counts), colnames(pop))
  counts <- counts[common_years, common_ages, drop = FALSE]
  pop    <- pop[common_years, common_ages, drop = FALSE]

  # BAPC expects data.frames with periods as row names and age groups as columns
  counts <- as.data.frame(counts)
  pop    <- as.data.frame(pop)

  message("Running BAPC for ", slug, " (", nrow(counts), " years x ", ncol(counts), " ages)")

  apc <- APCList(epi = counts, pyrs = pop, gf = 5)

  # Fit with 10-year ahead projection and retrospective smoothing
  fit <- BAPC(apc,
              predict = list(npredict = 10, retro = TRUE),
              secondDiff = FALSE,
              verbose = FALSE)

  # Extract projected counts and rates
  proj_counts <- agespec.proj(fit)
  proj_rates  <- agespec.rate(fit)

  write.csv(proj_counts, file.path(out_dir, paste0(slug, "_projected_counts.csv")))
  write.csv(proj_rates,  file.path(out_dir, paste0(slug, "_projected_rates.csv")))

  # Summary of hyperparameters
  hyper <- summaryHyper(fit)
  write.csv(hyper, file.path(out_dir, paste0(slug, "_hyperparameters.csv")))

  # Plots
  pdf(file.path(plot_dir, paste0(slug, "_ageSpecProj.pdf")), width = 12, height = 8)
  plotBAPC(fit, type = "ageSpecProj", scale = 1)
  dev.off()

  pdf(file.path(plot_dir, paste0(slug, "_ageSpecRate.pdf")), width = 12, height = 8)
  plotBAPC(fit, type = "ageSpecRate", scale = 100000)
  dev.off()

  message("Finished ", slug)
  list(slug = slug, fit = fit)
}

results <- lapply(files, run_one)
names(results) <- sapply(results, function(x) x$slug)

message("BAPC projections complete. Outputs in ", out_dir)
