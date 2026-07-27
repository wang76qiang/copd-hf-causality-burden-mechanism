# ============================================================================
# Module 03 - Mediation MR, step A (data retrieval, cached)
# Fetch (i) the eQTLGen whole-blood SERPINE1 (ENSG00000106366) cis region and
# (ii) TGFB1 (prot-a-2962) instruments from OpenGWAS. Requires OPENGWAS_JWT.
# Output: data_downloads/serpine1_eqtla_region.rds, iv_med_TGFB1.rds
#         (both shipped in data/derived/mediation/)
# ============================================================================

# --- repository paths -----------------------------------------------------------
# Run this script from the repository root, or set REPRO_ROOT to the repo root.
REPRO_ROOT <- Sys.getenv("REPRO_ROOT", unset = ".")
# GBD_RAW_DIR: folder with the raw IHME GBD 2021 csv downloads (see data/README.md)
GBD_RAW_DIR <- Sys.getenv("GBD_RAW_DIR", unset = file.path(REPRO_ROOT, "data", "raw", "gbd_2021"))
# V8_DIR: legacy V8 project tree (NOT shipped; only needed by v8_* legacy scripts)
V8_DIR <- Sys.getenv("V8_DIR", unset = "v8_legacy")
# .libPaths(...)  # removed user-local library path; see r_packages.txt
source(file.path(REPRO_ROOT, "code", "03_mediation", "t1_mediation_helpers.R"))
# ---- SERPINE1 (eqtl-a whole blood) region data: 1 API request, cached ----
reg_f <- paste0(DL, "/serpine1_eqtla_region.rds")
if (!file.exists(reg_f)) {
  reg <- ieugwasr::associations("7:99770302-101782566", "eqtl-a-ENSG00000106366", proxies = 0)
  saveRDS(reg, reg_f)
  cat("eqtl-a SERPINE1 region rows:", nrow(reg), "\n")
} else cat("eqtl-a region cached\n")
# ---- TGFB1 instruments (tophits, cached) ----
if (!file.exists(paste0(DL, "/iv_med_TGFB1.rds"))) {
  ex <- tryCatch(TwoSampleMR::extract_instruments("prot-a-2962", p1 = 5e-8, clump = TRUE, r2 = 0.01, kb = 10000),
                 error = function(e) NULL)
  if (is.null(ex) || nrow(ex) == 0)
    ex <- tryCatch(TwoSampleMR::extract_instruments("prot-a-2962", p1 = 5e-6, clump = TRUE, r2 = 0.01, kb = 10000),
                   error = function(e) NULL)
  if (is.null(ex)) ex <- data.frame()
  saveRDS(ex, paste0(DL, "/iv_med_TGFB1.rds"))
  cat("TGFB1 instruments:", nrow(ex), "\n")
} else cat("TGFB1 instruments cached\n")
cat("STEP A DONE\n")
