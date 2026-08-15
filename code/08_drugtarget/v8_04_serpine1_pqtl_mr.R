# ============================================================================
# Module 08 - Drug target (LEGACY V8 script): SERPINE1/PAI-1 cis-pQTL MR
# Protein-level MR of SERPINE1 (deCODE cis-pQTL, prot-c-2925_9_1) on HF via
# OpenGWAS (requires OPENGWAS_JWT or IEU_GWAS_TOKEN). Result: OR 1.01,
# P=0.69 (null).
# Output: serpine1_pqtl_mr_* (archived in results/legacy/)
# ============================================================================

# --- repository paths -----------------------------------------------------------
# Run this script from the repository root, or set REPRO_ROOT to the repo root.
REPRO_ROOT <- Sys.getenv("REPRO_ROOT", unset = ".")
# GBD_RAW_DIR: folder with the raw IHME GBD 2021 csv downloads (see data/README.md)
GBD_RAW_DIR <- Sys.getenv("GBD_RAW_DIR", unset = file.path(REPRO_ROOT, "data", "raw", "gbd_2021"))
# V8_DIR: legacy V8 project tree (NOT shipped; only needed by v8_* legacy scripts)
V8_DIR <- Sys.getenv("V8_DIR", unset = "v8_legacy")
# .libPaths(...)  # removed user-local library path; see r_packages.txt

# 04_serpine1_pqtl_mr.R
# Protein-level MR for SERPINE1 / PAI-1 on Heart Failure
# Requires: TwoSampleMR, ieugwasr (with valid JWT token), or local pQTL summary stats

library(TwoSampleMR)
library(ieugwasr)
library(jsonlite)

out_prefix <- file.path(V8_DIR, "R_scripts/results/serpine1_pqtl_mr")

# Try to obtain IEU token from environment or a local file
ieu_token <- Sys.getenv("IEU_GWAS_TOKEN")
if (ieu_token == "" && file.exists(file.path(V8_DIR, "ieu_token.txt"))) {
  ieu_token <- trimws(readLines(file.path(V8_DIR, "ieu_token.txt")))
}

results <- list(
  status = "pending_external_data",
  note = "SERPINE1 pQTL MR requires plasma pQTL summary statistics (INTERVAL/deCODE/UKB-PPP) and an IEU OpenGWAS token."
)

if (ieu_token != "") {
  cat("IEU token found; querying OpenGWAS for SERPINE1 pQTL instruments...\n")
  Sys.setenv(OPENGWAS_JWT = ieu_token)
  cat("Token set. Validating...\n")
  print(ieugwasr::check_access_token())

  # Attempt to locate SERPINE1 pQTL datasets in IEU OpenGWAS
  # Common IDs depend on the current IEU database; we try several documented examples.
  possible_pqtl_ids <- c("prot-c-2925_9_1", "prot-a-1295", "prot-b-1295")
  serpine1_pqtl <- NULL
  for (pid in possible_pqtl_ids) {
    tryCatch({
      inst <- extract_instruments(
        outcomes = pid,
        p1 = 5e-6,
        clump = TRUE,
        r2 = 0.001,
        kb = 10000
      )
      if (!is.null(inst) && nrow(inst) > 0) {
        serpine1_pqtl <- inst
        results$pqtl_dataset <- pid
        break
      }
    }, error = function(e) cat(pid, "failed:", conditionMessage(e), "\n"))
  }

  if (!is.null(serpine1_pqtl)) {
    cat("Found", nrow(serpine1_pqtl), "SERPINE1 cis-pQTL instruments.\n")
    hf_outcome <- extract_outcome_data(snps = serpine1_pqtl$SNP, outcomes = "ebi-a-GCST009541")
    harm <- harmonise_data(serpine1_pqtl, hf_outcome)
    res <- mr(harm, method_list = c("mr_ivw", "mr_weighted_median", "mr_egger_regression"))
    pleio <- mr_pleiotropy_test(harm)
    hetero <- mr_heterogeneity(harm)

    results$status <- "completed"
    results$note <- NULL
    results$n_instruments <- nrow(serpine1_pqtl)
    results$mr_results <- as.data.frame(res)
    results$pleiotropy <- as.data.frame(pleio)
    results$heterogeneity <- as.data.frame(hetero)

    write.csv(res, paste0(out_prefix, "_mr_results.csv"), row.names = FALSE)
    write.csv(pleio, paste0(out_prefix, "_pleiotropy.csv"), row.names = FALSE)
    write.csv(harm, paste0(out_prefix, "_harmonised.csv"), row.names = FALSE)
  } else {
    results$status <- "no_pqtl_dataset_found"
    cat("No SERPINE1 pQTL dataset found with tried IDs. Provide local pQTL file or correct IEU ID.\n")
  }
} else {
  cat("No IEU token provided. Saving placeholder result.\n")
  cat("To run this analysis, set IEU_GWAS_TOKEN environment variable or create {V8 legacy tree, not shipped}/ieu_token.txt\n")
}

write_json(results, paste0(out_prefix, "_summary.json"), pretty = TRUE, auto_unbox = TRUE)
print(results)
