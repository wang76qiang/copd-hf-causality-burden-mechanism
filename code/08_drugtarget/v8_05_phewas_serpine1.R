# ============================================================================
# Module 08 - Drug target (LEGACY V8 script): SERPINE1 PheWAS
# PhenoScanner PheWAS of the SERPINE1 cis-pQTL lead SNP rs7860931.
# Output: phewas_serpine1_* (archived in results/legacy/)
# ============================================================================

# --- repository paths -----------------------------------------------------------
# Run this script from the repository root, or set REPRO_ROOT to the repo root.
REPRO_ROOT <- Sys.getenv("REPRO_ROOT", unset = ".")
# GBD_RAW_DIR: folder with the raw IHME GBD 2021 csv downloads (see data/README.md)
GBD_RAW_DIR <- Sys.getenv("GBD_RAW_DIR", unset = file.path(REPRO_ROOT, "data", "raw", "gbd_2021"))
# V8_DIR: legacy V8 project tree (NOT shipped; only needed by v8_* legacy scripts)
V8_DIR <- Sys.getenv("V8_DIR", unset = "v8_legacy")
# .libPaths(...)  # removed user-local library path; see r_packages.txt

# 05_phewas_serpine1.R
# PheWAS of SERPINE1 cis-pQTL instrument across UK Biobank and FinnGen phenotypes
# Requires: phenoscanner R package (internet access) or local PheWAS summary stats

library(phenoscanner)
library(jsonlite)

out_prefix <- file.path(V8_DIR, "R_scripts/results/phewas_serpine1")

# Top SERPINE1 cis-pQTL SNP (example; replace with actual lead SNP from pQTL data)
serpine1_snp <- "rs7860931"

results <- list(
  snp = serpine1_snp,
  status = "pending",
  note = "PheWAS requires a validated SERPINE1 cis-pQTL lead SNP and PhenoScanner/UKB-FinnGen access."
)

phewas_results <- tryCatch({
  res <- phenoscanner(snpquery = serpine1_snp, pvalue = 1e-5)
  df <- res$results
  cat("PhenoScanner returned", nrow(df), "associations.\n")
  write.csv(df, file = paste0(out_prefix, "_phenoscanner_results.csv"), row.names = FALSE)
  list(status = "completed", n_associations = nrow(df), top_results = head(df, 20))
}, error = function(e) {
  cat("PhenoScanner query failed:", conditionMessage(e), "\n")
  list(status = "failed", error = conditionMessage(e))
})

results$status <- phewas_results$status
results$phewas_results <- phewas_results

write_json(results, paste0(out_prefix, "_summary.json"), pretty = TRUE, auto_unbox = TRUE)
print(results)
