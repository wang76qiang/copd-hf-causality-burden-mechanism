# ============================================================================
# Module 02 - MVMR
# Build the model-1 (A1A: COPD+SMOK+BMI+SBP) and model-2 (B2:
# COPD+FEV1+FEV1/FVC+CPD) union SNP lists and per-dataset rsid list files.
# Input : data_downloads/snps10.rds, iv_final_<KEY>.tsv (data/derived/mvmr/)
# Output: rsidlist_A1A.txt, rsidlist_B2.txt, t1_iv_sets.rds (data/derived/mvmr/)
# ============================================================================

# --- repository paths -----------------------------------------------------------
# Run this script from the repository root, or set REPRO_ROOT to the repo root.
REPRO_ROOT <- Sys.getenv("REPRO_ROOT", unset = ".")
# GBD_RAW_DIR: folder with the raw IHME GBD 2021 csv downloads (see data/README.md)
GBD_RAW_DIR <- Sys.getenv("GBD_RAW_DIR", unset = file.path(REPRO_ROOT, "data", "raw", "gbd_2021"))
# V8_DIR: legacy V8 project tree (NOT shipped; only needed by v8_* legacy scripts)
V8_DIR <- Sys.getenv("V8_DIR", unset = "v8_legacy")
# .libPaths(...)  # removed user-local library path; see r_packages.txt
# Build A1A union SNP list (COPD task-10 + clumped covariate IVs) and write
# per-dataset rsid list files for shell extraction.
dl <- file.path(REPRO_ROOT, "data_downloads")
snps10 <- readRDS(paste0(dl, "/snps10.rds"))
mk <- function(key) {
  d <- read.delim(paste0(dl, "/iv_final_", key, ".tsv"), stringsAsFactors = FALSE)
  d$rsid
}
sets <- list(COPD = snps10, SMOK = mk("SMOK"), BMI = mk("BMI"), SBP = mk("SBP"))
union <- unique(unlist(sets))
cat("A1A union:", length(union), " (COPD", length(sets$COPD), "SMOK", length(sets$SMOK),
    "BMI", length(sets$BMI), "SBP", length(sets$SBP), ")\n")
writeLines(union, paste0(dl, "/rsidlist_A1A.txt"))
setsB2 <- list(COPD = snps10, FEV1 = mk("FEV1"), FEV1FVC = mk("FEV1FVC"), CPD = mk("CPD"))
unionB2 <- unique(unlist(setsB2))
cat("B2 union:", length(unionB2), "\n")
writeLines(unionB2, paste0(dl, "/rsidlist_B2.txt"))
saveRDS(list(A1A = sets, B2 = setsB2), paste0(dl, "/t1_iv_sets.rds"))
