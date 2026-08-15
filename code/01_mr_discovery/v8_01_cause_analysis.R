# ============================================================================
# Module 01 - MR discovery (LEGACY V8 script, kept for provenance)
# CAUSE (Morrison et al. 2020) on the harmonised COPD -> HF dataset.
# With only 10 genome-wide significant instruments the full LOO-CV model
# comparison is underpowered; the causal-model posterior gamma is reported as
# a sensitivity measure (see script header caveat).
# Input : data/derived/mr_discovery/discovery_harmonised_10snp_dat.xlsx
# Output: cause_results_* (archived in results/legacy/)
# ============================================================================

# --- repository paths -----------------------------------------------------------
# Run this script from the repository root, or set REPRO_ROOT to the repo root.
REPRO_ROOT <- Sys.getenv("REPRO_ROOT", unset = ".")
# GBD_RAW_DIR: folder with the raw IHME GBD 2021 csv downloads (see data/README.md)
GBD_RAW_DIR <- Sys.getenv("GBD_RAW_DIR", unset = file.path(REPRO_ROOT, "data", "raw", "gbd_2021"))
# V8_DIR: legacy V8 project tree (NOT shipped; only needed by v8_* legacy scripts)
V8_DIR <- Sys.getenv("V8_DIR", unset = "v8_legacy")
# .libPaths(...)  # removed user-local library path; see r_packages.txt

# 01_cause_analysis.R
# Run CAUSE (Morrison et al. 2020) on COPD -> Heart Failure summary statistics
# Requires: cause R package, TwoSampleMR, tidyverse
# NOTE: CAUSE is designed for genome-wide summary statistics (>100k variants).
# The current local harmonised dataset contains only 10 genome-wide significant
# instruments. We therefore (i) attempt the full CAUSE fit, (ii) if the small-
# sample LOO-CV step fails, we fit the causal model separately and report the
# posterior causal effect (gamma) as a sensitivity measure, with an explicit
# caveat that full CAUSE model comparison requires genome-wide data.

library(cause)
library(TwoSampleMR)
library(tidyverse)
library(readxl)

# File paths
data_file <- file.path(REPRO_ROOT, "data", "derived", "mr_discovery", "discovery_harmonised_10snp_dat.xlsx")
out_prefix <- file.path(V8_DIR, "R_scripts/results/cause_results")

# Read harmonized data
dat <- readxl::read_excel(data_file) %>%
  filter(mr_keep == TRUE) %>%
  select(
    snp = SNP,
    beta_hat_1 = beta.exposure,
    beta_hat_2 = beta.outcome,
    seb1 = se.exposure,
    seb2 = se.outcome,
    A1 = effect_allele.exposure,
    A2 = other_allele.exposure
  ) %>%
  mutate(
    snp = as.character(snp),
    A1 = as.character(A1),
    A2 = as.character(A2)
  )

n_snp <- nrow(dat)
cat("Number of instruments:", n_snp, "\n")

# Estimate CAUSE parameters (mixture grid and rho)
params <- est_cause_params(dat, variants = dat$snp)

# Attempt full CAUSE fit
res <- tryCatch(
  cause(dat, params, variants = dat$snp),
  error = function(e) {
    cat("Full CAUSE fit failed:", conditionMessage(e), "\n")
    return(NULL)
  }
)

summary_out <- list()
summary_out$n_snp <- n_snp
summary_out$params_converged <- params$converged
summary_out$params_rho <- as.numeric(params$rho)

if (!is.null(res)) {
  summary_res <- summary(res)
  print(summary_res)
  elpd_diff <- as.numeric(summary_res$elpd["causal", ] - summary_res$elpd["sharing", ])
  causal_effect <- as.numeric(summary_res$z["causal"])
  summary_out$full_cause_success <- TRUE
  summary_out$elpd_diff_causal_minus_sharing <- elpd_diff
  summary_out$causal_effect_z <- causal_effect
  saveRDS(res, file = paste0(out_prefix, "_cause_fit.rds"))
} else {
  # Fallback: fit causal model separately and sample gamma posterior
  cat("Fitting causal model separately for sensitivity estimate...\n")
  cd <- new_cause_data(dat)
  sigma_g <- cause:::eta_gamma_prior(cd)
  fit_causal <- cause_grid_adapt(
    cd, params, max_post_per_bin = 0.001,
    params = c("gamma", "eta", "q"),
    priors = list(
      function(b) dnorm(b, 0, sigma_g),
      function(b) dnorm(b, 0, sigma_g),
      function(q) dbeta(q, 1, 10)
    ),
    n_start = c(21, 21, 11)
  )

  # Sample from the adapted grid using the normalized posterior
  jp <- fit_causal$joint_post
  probs <- jp$norm_log_post
  probs <- exp(probs - max(probs))
  probs <- probs / sum(probs)
  idx <- sample.int(nrow(jp), size = 10000, replace = TRUE, prob = probs)
  samps <- jp[idx, c("gamma", "eta", "q")]

  gamma_mean <- mean(samps$gamma)
  gamma_ci <- quantile(samps$gamma, c(0.025, 0.975))
  gamma_or <- exp(gamma_mean)
  gamma_or_ci <- exp(gamma_ci)

  cat("Causal model gamma (log-OR):", round(gamma_mean, 4), "\n")
  cat("95% CI:", paste(round(gamma_ci, 4), collapse = " to "), "\n")
  cat("OR:", round(gamma_or, 3), "(95% CI", paste(round(gamma_or_ci, 3), collapse = "-"), ")\n")

  summary_out$full_cause_success <- FALSE
  summary_out$fallback_causal_model <- TRUE
  summary_out$causal_gamma_log_or <- as.numeric(gamma_mean)
  summary_out$causal_gamma_95ci_low <- as.numeric(gamma_ci[1])
  summary_out$causal_gamma_95ci_high <- as.numeric(gamma_ci[2])
  summary_out$causal_or <- as.numeric(gamma_or)
  summary_out$causal_or_95ci_low <- as.numeric(gamma_or_ci[1])
  summary_out$causal_or_95ci_high <- as.numeric(gamma_or_ci[2])
  summary_out$caveat <- "Full CAUSE model comparison (causal vs sharing) could not be performed because fewer than 100,000 genome-wide variants were available. The causal-model gamma is reported as a sensitivity estimate only."

  saveRDS(list(params = params, fit_causal = fit_causal, samples = samps),
          file = paste0(out_prefix, "_cause_fit_fallback.rds"))
}

print(summary_out)
write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE, auto_unbox = TRUE)
write_json(summary_out, paste0(out_prefix, "_cause_summary.json"))
