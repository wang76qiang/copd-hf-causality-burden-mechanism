# ============================================================================
# Module 02 - MVMR
# LD clump the raw p<5e-8 candidate instrument sets (r2=0.001, kb=10000,
# OpenGWAS LD API, EUR). Large candidate sets (BMI/FEV1/FEV1FVC/SMOK) are
# greedily pre-thinned (1 Mb window) ONLY to respect the LD-API request size;
# final selection is always LD-based.
# Input : data_downloads/iv_raw_<KEY>.tsv  (produced by extract_sumstats.sh)
# Output: data_downloads/iv_final_<KEY>.tsv (shipped in data/derived/mvmr/)
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
# LD clump the raw p<5e-8 candidate sets (r2=0.001, kb=10000).
# For very large candidate sets (BMI/FEV1/FEV1FVC/SMOK) a greedy p-sorted
# 1 Mb window pre-thin is applied first ONLY to respect the LD-API request
# size; the final selection is always the LD-based clump.
library(ieugwasr)
dl <- file.path(REPRO_ROOT, "data_downloads")
keys <- c("CPD","SBP","SMOK","BMI","FEV1","FEV1FVC","COPD")

greedy_thin <- function(d, window = 1e6, max_n = 3000) {
  d <- d[order(d$p), ]
  keep <- logical(nrow(d))
  taken <- new.env()
  for (i in seq_len(nrow(d))) {
    chr <- as.character(d$chr[i]); pos <- d$pos[i]
    v <- taken[[chr]]; if (is.null(v)) v <- numeric(0)
    if (length(v) == 0 || all(abs(pos - v) > window)) {
      keep[i] <- TRUE; taken[[chr]] <- c(v, pos)
    }
    if (sum(keep) >= max_n) break
  }
  d[keep, ]
}

for (k in keys) {
  f_out <- paste0(dl, "/iv_final_", k, ".tsv")
  if (file.exists(f_out)) next
  d <- read.delim(paste0(dl, "/iv_raw_", k, ".tsv"), stringsAsFactors = FALSE)
  d <- d[!duplicated(d$rsid) & !is.na(d$p) & !is.na(d$beta) & d$se > 0, ]
  n0 <- nrow(d)
  if (n0 > 3000) { d <- greedy_thin(d); cat(k, "pre-thinned", n0, "->", nrow(d), "\n") }
  t0 <- Sys.time()
  cl <- tryCatch(ieugwasr::ld_clump(data.frame(rsid = d$rsid, pval = d$p),
                                    clump_r2 = 0.001, clump_kb = 10000, pop = "EUR"),
                 error = function(e) { cat(k, "clump ERROR:", conditionMessage(e), "\n"); NULL })
  if (is.null(cl)) next
  d2 <- d[d$rsid %in% cl$rsid, ]
  write.table(d2, f_out, sep = "\t", row.names = FALSE, quote = FALSE)
  cat(k, "clumped:", nrow(d2), " in", round(difftime(Sys.time(), t0, units="secs"),1), "s\n")
  Sys.sleep(1)
}
cat("ALL CLUMP DONE\n")
