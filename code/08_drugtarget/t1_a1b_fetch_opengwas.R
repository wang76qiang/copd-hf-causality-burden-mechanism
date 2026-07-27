# ============================================================================
# Module 08 - Drug target: OpenGWAS retrieval (requires OPENGWAS_JWT)
# Fetch (1) HF betas for the 8 tissue lead cis SNPs, (2) HF full summary
# stats in the SERPINE1 region for coloc, (3) eQTLGen whole-blood SERPINE1
# cis region.
# Output: data_downloads/t1_a1b_hf_lead.rds, t1_a1b_hf_region.rds,
#         serpine1_eqtla_region.rds (shipped in data/derived/drugtarget|mediation)
# ============================================================================

# --- repository paths -----------------------------------------------------------
# Run this script from the repository root, or set REPRO_ROOT to the repo root.
REPRO_ROOT <- Sys.getenv("REPRO_ROOT", unset = ".")
# GBD_RAW_DIR: folder with the raw IHME GBD 2021 csv downloads (see data/README.md)
GBD_RAW_DIR <- Sys.getenv("GBD_RAW_DIR", unset = file.path(REPRO_ROOT, "data", "raw", "gbd_2021"))
# V8_DIR: legacy V8 project tree (NOT shipped; only needed by v8_* legacy scripts)
V8_DIR <- Sys.getenv("V8_DIR", unset = "v8_legacy")
# .libPaths(...)  # removed user-local library path; see r_packages.txt
# NOTE: an OpenGWAS JWT is required for OpenGWAS/ieugwasr calls.
# Request your own token at https://api.opengwas.io and either
#   export OPENGWAS_JWT=<token>          (shell, picked up automatically), or
#   Sys.setenv(OPENGWAS_JWT = "<token>")  (in R, before running).
if (Sys.getenv("OPENGWAS_JWT") == "")
  message("WARNING: OPENGWAS_JWT is not set; OpenGWAS API calls will fail.")
# A1b: fetch from OpenGWAS (3 requests):
#  1) HF betas for the 8 tissue lead cis SNPs (exploratory Wald ratios)
#  2) HF full summary stats in SERPINE1 region (b37 window) for coloc
#  3) eQTLGen whole-blood SERPINE1 cis region (eqtl-a-ENSG00000106366)
library(ieugwasr)
dl <- file.path(REPRO_ROOT, "data_downloads")
cand <- readRDS(paste0(dl, "/t1_a1b_candidates.rds"))
lead <- unique(sapply(names(cand), function(t) {
  g <- read.csv(paste0(dl, "/serpine1_gene/", t, ".csv"), stringsAsFactors = FALSE)
  g$rsid[which.min(g$pvalue)]
}))
cat("lead SNPs:", paste(lead, collapse = ", "), "\n")

f1 <- paste0(dl, "/t1_a1b_hf_lead.rds")
if (!file.exists(f1)) {
  a <- ieugwasr::associations(lead, "ebi-a-GCST009541", proxies = 0)
  saveRDS(a, f1); cat("HF lead assoc rows:", nrow(a), "\n")
}
Sys.sleep(3)
f2 <- paste0(dl, "/t1_a1b_hf_region.rds")
if (!file.exists(f2)) {
  r <- ieugwasr::associations("7:99770302-101782566", "ebi-a-GCST009541", proxies = 0)
  saveRDS(r, f2); cat("HF region rows:", nrow(r), "\n")
}
Sys.sleep(3)
f3 <- paste0(dl, "/t1_a1b_eqtlgen_region.rds")
if (!file.exists(f3)) {
  e <- tryCatch(ieugwasr::associations("7:100126167-101138431", "eqtl-a-ENSG00000106366", proxies = 0),
                error = function(e2) { cat("eqtlgen b38 window failed:", conditionMessage(e2), "\n"); NULL })
  if (is.null(e) || nrow(e) == 0) {
    e <- ieugwasr::associations("7:99770302-101782566", "eqtl-a-ENSG00000106366", proxies = 0)
    cat("eqtlgen b37 window rows:", nrow(e), "\n")
  }
  saveRDS(e, f3); cat("eQTLGen region rows:", nrow(e), "\n")
}
