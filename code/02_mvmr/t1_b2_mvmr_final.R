# ============================================================================
# Module 02 - MVMR model 2 (FINAL)
# MVMR COPD + FEV1 (GCST007432) + FEV1/FVC (GCST007431) + cigarettes-per-day
# (GCST009968) -> HF (GCST009541). 405 SNPs; COPD direct OR 1.165.
# Input : data_downloads/assoc_B2_<KEY>.tsv (shipped in data/derived/mvmr/)
# Output: results/t1_mvmr_lungfunction.csv, t1_b2_harmonised_mvmr.csv,
#         t1_b2_mvmr_qstat.csv
# ============================================================================

# --- repository paths -----------------------------------------------------------
# Run this script from the repository root, or set REPRO_ROOT to the repo root.
REPRO_ROOT <- Sys.getenv("REPRO_ROOT", unset = ".")
# GBD_RAW_DIR: folder with the raw IHME GBD 2021 csv downloads (see data/README.md)
GBD_RAW_DIR <- Sys.getenv("GBD_RAW_DIR", unset = file.path(REPRO_ROOT, "data", "raw", "gbd_2021"))
# V8_DIR: legacy V8 project tree (NOT shipped; only needed by v8_* legacy scripts)
V8_DIR <- Sys.getenv("V8_DIR", unset = "v8_legacy")
# .libPaths(...)  # removed user-local library path; see r_packages.txt
# B2 final: MVMR  COPD + FEV1 + FEV1/FVC + CPD -> HF  (all local EBI files)
library(MVMR)
dl <- file.path(REPRO_ROOT, "data_downloads")
OUT <- file.path(REPRO_ROOT, "results")
LOG <- file(file.path(REPRO_ROOT, "logs/t1_b2_mvmr_lungfunction.log"), open = "wt"); sink(LOG, split = TRUE)
cat("B2 MVMR run @", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Datasets: COPD=ebi-a-GCST90018807; FEV1=ebi-a-GCST007432 (Shrine 2019, N=321,047);\n")
cat("FEV1/FVC=ebi-a-GCST007431 (Shrine 2019, N=321,047);\n")
cat("CPD=ebi-a-GCST009968 (Buchwald 2020 cigarettes/day, N=4,772) -- SUBSTITUTE for\n")
cat("ieu-b-25 (GSCAN CPD, N=337k): GSCAN full sumstats unavailable (EBI: no;\n")
cat("conservancy.umn.edu: 403); GCST009968 is the largest EBI-file CPD GWAS but is\n")
cat("small -> only 2 instruments, weak. Interpret CPD-adjusted estimates with caution.\n")
cat("HF=ebi-a-GCST009541. Instruments: p<5e-8, clump r2<0.001 kb=10000 (LD API EUR).\n\n")

keys <- c("COPD","FEV1","FEV1FVC","CPD","HF")
D <- lapply(keys, function(k) {
  d <- read.delim(paste0(dl, "/assoc_B2_", k, ".tsv"), stringsAsFactors = FALSE)
  d <- d[!duplicated(d$rsid), ]; rownames(d) <- d$rsid
  d$ea <- toupper(d$ea); d$nea <- toupper(d$nea); d
})
names(D) <- keys
common <- Reduce(intersect, lapply(D, function(d) d$rsid))
cat("SNPs in all 5 datasets:", length(common), "\n")

comp <- function(a) chartr("ACGT", "TGCA", toupper(a))
align_one <- function(r1, r2, e1, e2) {
  if (e1 == r1 && e2 == r2) return(1); if (e1 == r2 && e2 == r1) return(-1)
  if (comp(e1) == r1 && comp(e2) == r2) return(1); if (comp(e1) == r2 && comp(e2) == r1) return(-1)
  0
}
BX <- SE <- matrix(NA_real_, length(common), 4, dimnames = list(common, keys[1:4]))
BY <- SEY <- setNames(rep(NA_real_, length(common)), common)
eaf_mat <- matrix(NA_real_, length(common), 5, dimnames = list(common, keys))
palin <- setNames(rep(FALSE, length(common)), common); incompat <- character(0)
for (s in common) {
  ref <- D$HF[s, ]
  palin[s] <- paste0(sort(c(ref$ea, ref$nea)), collapse = "") %in% c("AT","CG")
  ok <- TRUE
  for (k in keys) {
    d <- D[[k]][s, ]
    sgn <- align_one(ref$ea, ref$nea, d$ea, d$nea)
    if (sgn == 0) { ok <- FALSE; break }
    b <- sgn * d$beta
    if (k == "HF") { BY[s] <- b; SEY[s] <- d$se } else { BX[s, k] <- b; SE[s, k] <- d$se }
    e <- suppressWarnings(as.numeric(d$eaf))
    eaf_mat[s, k] <- if (is.na(e)) NA else if (sgn == 1) e else 1 - e
  }
  if (!ok) incompat <- c(incompat, s)
}
keep <- setdiff(common, incompat)
pal_drop <- character(0)
for (s in keep[palin[keep]]) {
  e <- eaf_mat[s, ]
  if (all(!is.na(e)) && max(e) - min(e) > 0.1) pal_drop <- c(pal_drop, s)
}
keep <- setdiff(keep, pal_drop)
cat("dropped incompatible:", length(incompat), "; palindromic:", length(pal_drop),
    "; final joint:", length(keep), "\n")

BX <- BX[keep, ]; SE <- SE[keep, ]; BYv <- BY[keep]; SEYv <- SEY[keep]
r_input <- MVMR::format_mvmr(BXGs = as.data.frame(BX), BYG = BYv,
                             seBXGs = as.data.frame(SE), seBYG = SEYv, RSID = keep)
sres <- MVMR::strength_mvmr(r_input, gencov = 0)
cat("\n=== conditional F ===\n"); print(sres)
mres <- MVMR::mvmr(r_input, gencov = 0)
cat("\n=== multivariable IVW ===\n"); print(mres$coef)

est <- as.data.frame(mres$coef); colexp <- keys[1:4]
mvmr_out <- data.frame(exposure = colexp, beta = est[, "Estimate"], se = est[, "Std. Error"],
                       pval = est[, grep("Pr", colnames(est), value = TRUE)[1]])
mvmr_out$or <- exp(mvmr_out$beta)
mvmr_out$or_lci95 <- exp(mvmr_out$beta - 1.96*mvmr_out$se)
mvmr_out$or_uci95 <- exp(mvmr_out$beta + 1.96*mvmr_out$se)
mvmr_out$condF <- as.numeric(as.data.frame(sres)[1, ])[match(mvmr_out$exposure, colexp)]
mvmr_out$nsnp_joint <- length(keep)
mvmr_out$dataset_id <- c(COPD = "ebi-a-GCST90018807", FEV1 = "ebi-a-GCST007432",
                         FEV1FVC = "ebi-a-GCST007431", CPD = "ebi-a-GCST009968")[mvmr_out$exposure]
write.csv(mvmr_out, file.path(OUT, "t1_mvmr_lungfunction.csv"), row.names = FALSE)
harm <- data.frame(SNP = keep, BX, SE, outcome_beta = BYv, outcome_se = SEYv,
                   palindromic = palin[keep], eaf_mat[keep, ])
write.csv(harm, file.path(OUT, "t1_b2_harmonised_mvmr.csv"), row.names = FALSE)
cat("\n=== FINAL ===\n"); print(mvmr_out)
sink(); close(LOG)
cat("B2 DONE\n")
