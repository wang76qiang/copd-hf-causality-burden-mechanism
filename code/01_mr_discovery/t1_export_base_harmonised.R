# ============================================================================
# Module 01 - MR discovery (COPD -> HF, European discovery cohort)
# Export the local 10-SNP harmonised COPD (ebi-a-GCST90018807) -> HF
# (ebi-a-GCST009541) dataset to a plain CSV used as the exposure instrument
# basis for modules 02/03/04.
# Input : data/derived/mr_discovery/discovery_harmonised_10snp_dat.xlsx
# Output: results/t1_base_harmonised_10snp.csv
# ============================================================================

# --- repository paths -----------------------------------------------------------
# Run this script from the repository root, or set REPRO_ROOT to the repo root.
REPRO_ROOT <- Sys.getenv("REPRO_ROOT", unset = ".")
# GBD_RAW_DIR: folder with the raw IHME GBD 2021 csv downloads (see data/README.md)
GBD_RAW_DIR <- Sys.getenv("GBD_RAW_DIR", unset = file.path(REPRO_ROOT, "data", "raw", "gbd_2021"))
# V8_DIR: legacy V8 project tree (NOT shipped; only needed by v8_* legacy scripts)
V8_DIR <- Sys.getenv("V8_DIR", unset = "v8_legacy")
# .libPaths(...)  # removed user-local library path; see r_packages.txt
# Export the local 10-SNP harmonised COPD(ebi-a-GCST90018807) -> HF(ebi-a-GCST009541)
# dataset (source: data/derived/mr_discovery/discovery_harmonised_10snp_dat.xlsx) to a plain CSV
# used as the exposure instrument basis for tasks A1a / B1.
library(readxl)
d <- as.data.frame(read_excel(file.path(REPRO_ROOT, "data", "derived", "mr_discovery", "discovery_harmonised_10snp_dat.xlsx")))
d <- d[d$mr_keep == TRUE, ]
keep_cols <- c("SNP","chr.exposure","pos.exposure",
               "effect_allele.exposure","other_allele.exposure","eaf.exposure",
               "beta.exposure","se.exposure","pval.exposure","samplesize.exposure",
               "effect_allele.outcome","other_allele.outcome","eaf.outcome",
               "beta.outcome","se.outcome","pval.outcome","samplesize.outcome",
               "id.exposure","exposure","id.outcome","outcome")
out <- d[, keep_cols]
write.csv(out, file.path(REPRO_ROOT, "results/t1_base_harmonised_10snp.csv"), row.names = FALSE)
cat("Exported", nrow(out), "SNPs\n")
print(out[, c("SNP","effect_allele.exposure","other_allele.exposure","eaf.exposure","beta.exposure","se.exposure")])
