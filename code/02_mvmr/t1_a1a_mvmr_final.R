# ============================================================================
# Module 02 - MVMR model 1 (FINAL)
# Allele harmonisation across 5 local datasets + multivariable IVW (MVMR
# package) + Sanderson-Windmeijer conditional F.
# COPD + smoking (ebi-a-GCST90029014) + BMI (GCST90029007) + SBP
# (GCST90029011) -> HF (GCST009541). 709 SNPs; COPD direct OR 1.075.
# Input : data_downloads/assoc_A1A_<KEY>.tsv (shipped in data/derived/mvmr/)
# Output: results/t1_mvmr_results.csv, t1_a1a_harmonised_mvmr.csv,
#         t1_a1a_mvmr_qstat.csv
# ============================================================================

# --- repository paths -----------------------------------------------------------
# Run this script from the repository root, or set REPRO_ROOT to the repo root.
REPRO_ROOT <- Sys.getenv("REPRO_ROOT", unset = ".")
# GBD_RAW_DIR: folder with the raw IHME GBD 2021 csv downloads (see data/README.md)
GBD_RAW_DIR <- Sys.getenv("GBD_RAW_DIR", unset = file.path(REPRO_ROOT, "data", "raw", "gbd_2021"))
# V8_DIR: legacy V8 project tree (NOT shipped; only needed by v8_* legacy scripts)
V8_DIR <- Sys.getenv("V8_DIR", unset = "v8_legacy")
# .libPaths(...)  # removed user-local library path; see r_packages.txt
# A1A final: allele harmonisation across 5 local datasets + multivariable IVW
# (MVMR package) + Sanderson-Windmeijer conditional F.
library(MVMR)
dl <- file.path(REPRO_ROOT, "data_downloads")
OUT <- file.path(REPRO_ROOT, "results")
LOG <- file(file.path(REPRO_ROOT, "logs/t1_a1a_part1_mvmr.log"), open = "wt"); sink(LOG, split = TRUE)
cat("A1A MVMR final run @", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Datasets (all full summary statistics from EBI GWAS Catalog, harmonised files):\n")
cat("  COPD  = ebi-a-GCST90018807 (Sakaue 2021, N=468,475) [task-specified, 10-SNP IV]\n")
cat("  SMOK  = ebi-a-GCST90029014 (Smoking status, Loh PR 2018 UKB BOLT-LMM, N=468,170)\n")
cat("          SUBSTITUTE for ieu-b-4877 (GSCAN smoking initiation): GSCAN full sumstats\n")
cat("          not on EBI/OpenGWAS-files; OpenGWAS API too slow (~1.4 s/SNP).\n")
cat("  BMI   = ebi-a-GCST90029007 (Loh PR 2018, N=532,396); largest EBI-file BMI.\n")
cat("          (ieu-b-40 Yengo N=681k unavailable as file.)\n")
cat("  SBP   = ebi-a-GCST90029011 (Loh PR 2018, N=469,767)\n")
cat("          SUBSTITUTE for ebi-a-GCST90000066 (Surendran 2020): its EBI harmonised\n")
cat("          file contains only 241,748 variants (top-hits subset) -> unusable for MVMR.\n")
cat("  HF    = ebi-a-GCST009541 (Shah 2020 HERMES, N=977,323)\n")
cat("Instruments: p<5e-8, LD clump r2<0.001 kb=10000 (OpenGWAS LD API, EUR).\n")
cat("Large candidate sets pre-thinned (greedy p-sorted 1Mb window) ONLY to respect\n")
cat("LD-API request size; final selection always LD-based.\n\n")

keys <- c("COPD","SMOK","BMI","SBP","HF")
D <- lapply(keys, function(k) {
  d <- read.delim(paste0(dl, "/assoc_A1A_", k, ".tsv"), stringsAsFactors = FALSE)
  d <- d[!duplicated(d$rsid), ]
  rownames(d) <- d$rsid
  d
})
names(D) <- keys
common <- Reduce(intersect, lapply(D, function(d) d$rsid))
cat("SNPs present in all 5 datasets:", length(common), "\n")

comp <- function(a) chartr("ACGT", "TGCA", toupper(a))
align_one <- function(ref_ea, ref_nea, ea, nea) {
  ref_ea <- toupper(ref_ea); ref_nea <- toupper(ref_nea); ea <- toupper(ea); nea <- toupper(nea)
  # returns +1 (aligned), -1 (flipped), 0 (incompatible)
  if (ea == ref_ea && nea == ref_nea) return(1)
  if (ea == ref_nea && nea == ref_ea) return(-1)
  if (comp(ea) == ref_ea && comp(nea) == ref_nea) return(1)
  if (comp(ea) == ref_nea && comp(nea) == ref_ea) return(-1)
  0
}
snps10 <- readRDS(paste0(dl, "/snps10.rds"))
BX <- SE <- matrix(NA_real_, length(common), 4, dimnames = list(common, keys[1:4]))
BY <- SEY <- setNames(rep(NA_real_, length(common)), common)
eaf_mat <- matrix(NA_real_, length(common), 5, dimnames = list(common, keys))
palin <- setNames(rep(FALSE, length(common)), common)
incompat <- character(0)

for (s in common) {
  ref <- D$HF[s, ]                      # align everything to HF alleles
  refpair <- c(ref$ea, ref$nea)
  palin[s] <- paste0(sort(toupper(refpair)), collapse = "") %in% c("AT", "CG")
  ok <- TRUE
  for (k in keys) {
    d <- D[[k]][s, ]
    sgn <- align_one(ref$ea, ref$nea, d$ea, d$nea)
    if (sgn == 0) { ok <- FALSE; break }
    b <- sgn * d$beta
    if (k == "HF") { BY[s] <- b; SEY[s] <- d$se }
    else { BX[s, k] <- b; SE[s, k] <- d$se }
    e <- suppressWarnings(as.numeric(d$eaf))
    eaf_mat[s, k] <- if (is.na(e)) NA else if (sgn == 1) e else 1 - e
  }
  if (!ok) incompat <- c(incompat, s)
}
keep <- setdiff(common, incompat)
# palindromic SNPs: keep only if eaf consistent across datasets (|diff| < 0.1) or eaf missing
pal_drop <- character(0)
for (s in keep[palin[keep]]) {
  e <- eaf_mat[s, ]
  if (all(!is.na(e)) && max(e) - min(e) > 0.1) pal_drop <- c(pal_drop, s)
}
keep <- setdiff(keep, pal_drop)
cat("dropped allele-incompatible:", length(incompat), "; palindromic eaf-mismatch:", length(pal_drop),
    "; final joint SNPs:", length(keep), "\n")
cat("of the 10 task COPD SNPs in final set:", sum(snps10 %in% keep), "\n")

BX <- BX[keep, ]; SE <- SE[keep, ]; BYv <- BY[keep]; SEYv <- SEY[keep]

r_input <- MVMR::format_mvmr(BXGs = as.data.frame(BX), BYG = BYv,
                             seBXGs = as.data.frame(SE), seBYG = SEYv, RSID = keep)
sres <- MVMR::strength_mvmr(r_input, gencov = 0)
cat("\n=== conditional F (Sanderson-Windmeijer) ===\n"); print(sres)
mres <- MVMR::mvmr(r_input, gencov = 0)
cat("\n=== multivariable IVW (MVMR::mvmr) ===\n"); print(mres)

# MVMR::mvmr returns an MVMRIVW object; coefficients live in $coef
est <- as.data.frame(mres$coef)
colexp <- keys[1:4]
mvmr_out <- data.frame(
  exposure = colexp,
  beta = est[, "Estimate"], se = est[, "Std. Error"],
  pval = est[, grep("Pr", colnames(est), value = TRUE)[1]],
  stringsAsFactors = FALSE)
mvmr_out$or <- exp(mvmr_out$beta)
mvmr_out$or_lci95 <- exp(mvmr_out$beta - 1.96 * mvmr_out$se)
mvmr_out$or_uci95 <- exp(mvmr_out$beta + 1.96 * mvmr_out$se)
cf <- as.data.frame(sres); colnames(cf)[1] <- "condF"
mvmr_out$condF <- as.numeric(cf[1, ])[match(mvmr_out$exposure, colexp)]
mvmr_out$nsnp_joint <- length(keep)
mvmr_out$dataset_id <- c(COPD = "ebi-a-GCST90018807", SMOK = "ebi-a-GCST90029014",
                         BMI = "ebi-a-GCST90029007", SBP = "ebi-a-GCST90029011")[mvmr_out$exposure]
write.csv(mvmr_out, file.path(OUT, "t1_mvmr_results.csv"), row.names = FALSE)
harm <- data.frame(SNP = keep, BX, SE, outcome_beta = BYv, outcome_se = SEYv,
                   palindromic = palin[keep], eaf_mat[keep, ])
write.csv(harm, file.path(OUT, "t1_a1a_harmonised_mvmr.csv"), row.names = FALSE)
cat("\n=== FINAL TABLE ===\n"); print(mvmr_out)
sink(); close(LOG)
cat("A1A MVMR DONE\n")
