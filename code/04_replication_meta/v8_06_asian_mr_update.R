# ============================================================================
# Module 04 - Replication: BBJ East Asian MR (LEGACY V8 script)
# Asian-ancestry COPD (bbj-a-103) -> congestive HF (bbj-a-109) two-sample MR
# via OpenGWAS (requires OPENGWAS_JWT). Result: IVW OR 0.933 (6 instruments).
# Output: 06_asian_* (archived in results/legacy/)
# ============================================================================

# --- repository paths -----------------------------------------------------------
# Run this script from the repository root, or set REPRO_ROOT to the repo root.
REPRO_ROOT <- Sys.getenv("REPRO_ROOT", unset = ".")
# GBD_RAW_DIR: folder with the raw IHME GBD 2021 csv downloads (see data/README.md)
GBD_RAW_DIR <- Sys.getenv("GBD_RAW_DIR", unset = file.path(REPRO_ROOT, "data", "raw", "gbd_2021"))
# V8_DIR: legacy V8 project tree (NOT shipped; only needed by v8_* legacy scripts)
V8_DIR <- Sys.getenv("V8_DIR", unset = "v8_legacy")
# 06_asian_mr_update.R
# Re-run Asian-ancestry COPD -> HF MR using BBJ Congestive heart failure (bbj-a-109)
# as the outcome, matched to BBJ COPD (bbj-a-103) exposure.

# .libPaths(...)  # removed user-local library path; see r_packages.txt
library(TwoSampleMR)
library(ieugwasr)
library(dplyr)
library(readr)

# NOTE: an OpenGWAS JWT is required for OpenGWAS/ieugwasr calls.
# Request your own token at https://api.opengwas.io and either
#   export OPENGWAS_JWT=<token>          (shell, picked up automatically), or
#   Sys.setenv(OPENGWAS_JWT = "<token>")  (in R, before running).
if (Sys.getenv("OPENGWAS_JWT") == "")
  message("WARNING: OPENGWAS_JWT is not set; OpenGWAS API calls will fail.")

# Asian exposure and outcome IDs
asian_copd_id <- "bbj-a-103"   # BBJ COPD, East Asian
asian_hf_id   <- "bbj-a-109"   # BBJ Congestive heart failure, East Asian

outdir <- "R_scripts/results"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

cat("Extracting Asian COPD instruments from", asian_copd_id, "\n")
exp_asia <- extract_instruments(
  outcomes = asian_copd_id,
  p1 = 5e-8,
  clump = TRUE,
  r2 = 0.001,
  kb = 10000
)

cat("N instruments:", nrow(exp_asia), "\n")
if (nrow(exp_asia) == 0) stop("No Asian COPD instruments found")

cat("Extracting outcome data from", asian_hf_id, "\n")
out_asia <- extract_outcome_data(
  snps = exp_asia$SNP,
  outcomes = asian_hf_id
)

dat_asia <- harmonise_data(exp_asia, out_asia, action = 2)

# Keep only harmonisable SNPs
dat_asia <- dat_asia %>% filter(mr_keep == TRUE)

cat("N harmonised SNPs:", nrow(dat_asia), "\n")

# MR
mr_methods <- c("mr_ivw", "mr_weighted_median", "mr_egger_regression")
if (nrow(dat_asia) == 1) mr_methods <- c("mr_wald_ratio")

res_asia <- mr(dat_asia, method_list = mr_methods)
res_asia <- generate_odds_ratios(res_asia)

print(res_asia)
write_csv(res_asia, file.path(outdir, "06_asian_mr_bbj_congestive_hf_result.csv"))

# Sensitivity
if (nrow(dat_asia) >= 3) {
  pleio <- mr_pleiotropy_test(dat_asia)
  het <- mr_heterogeneity(dat_asia)
  print(pleio)
  print(het)
  write_csv(pleio, file.path(outdir, "06_asian_mr_pleiotropy.csv"))
  write_csv(het, file.path(outdir, "06_asian_mr_heterogeneity.csv"))
}

# Single-SNP Wald ratios for forest plot
wald <- dat_asia %>%
  mutate(
    b = beta.outcome / beta.exposure,
    se = abs(se.outcome / beta.exposure),
    or = exp(b),
    or_lci95 = exp(b - 1.96 * se),
    or_uci95 = exp(b + 1.96 * se),
    pval = pnorm(abs(b / se), lower.tail = FALSE) * 2
  ) %>%
  select(SNP, exposure, outcome, b, se, or, or_lci95, or_uci95, pval)
write_csv(wald, file.path(outdir, "06_asian_mr_wald_ratios.csv"))

# Save harmonised data for potential plotting
write_csv(dat_asia, file.path(outdir, "06_asian_harmonised_data.csv"))

cat("Asian MR updated.\n")
