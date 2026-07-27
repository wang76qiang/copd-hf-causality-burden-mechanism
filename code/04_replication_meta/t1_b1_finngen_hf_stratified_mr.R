# ============================================================================
# Module 04 - Replication: FinnGen R12 strict HF
# COPD -> FinnGen R12 I9_HEARTFAIL ("Heart failure, strict"; 37,653 cases /
# 462,695 controls) two-sample MR with the 10-SNP discovery instrument set.
# FinnGen R12 has NO HFpEF/HFrEF-specific endpoints (manifest checked
# 2026-07-26); strict-HF is used as the stratified/external-replication
# outcome. Result: IVW OR 1.037.
# Inputs: results/t1_base_harmonised_10snp.csv;
#         data/derived/finngen/fg_header.txt, fg_10snps_raw.tsv
#         (re-extractable from finngen_R12_I9_HEARTFAIL.gz, see data/README.md)
# Output: results/t1_hfpef_hfref_*.csv
# ============================================================================

# --- repository paths -----------------------------------------------------------
# Run this script from the repository root, or set REPRO_ROOT to the repo root.
REPRO_ROOT <- Sys.getenv("REPRO_ROOT", unset = ".")
# GBD_RAW_DIR: folder with the raw IHME GBD 2021 csv downloads (see data/README.md)
GBD_RAW_DIR <- Sys.getenv("GBD_RAW_DIR", unset = file.path(REPRO_ROOT, "data", "raw", "gbd_2021"))
# V8_DIR: legacy V8 project tree (NOT shipped; only needed by v8_* legacy scripts)
V8_DIR <- Sys.getenv("V8_DIR", unset = "v8_legacy")
# .libPaths(...)  # removed user-local library path; see r_packages.txt
# ==========================================================================
# Task B1: COPD -> FinnGen R12 heart failure stratified MR
# FinnGen R12 has NO HFpEF/HFrEF-specific endpoints (manifest checked
# 2026-07-26: closest endpoints are I9_HEARTFAIL "Heart failure, strict",
# 37,653 cases / 462,695 controls, plus composites such as
# I9_HEARTFAIL_AND_CHD). Attempted endpoint names returning HTTP 404:
#   I9_HF, I9_HEARTF, I9_HF_STRICT, I9_HFPEF, I9_HFREF, I9_HEARTFAIL
#   (non-release path), I9_HF_ALL  -> all 404
# Working file: finngen-public-data-r12/summary_stats/release/
#   finngen_R12_I9_HEARTFAIL.gz  (HTTP 200)
# Therefore this analysis uses I9_HEARTFAIL (strict HF) as the stratified /
# external-replication outcome. Exposure instruments: local harmonised
# 10-SNP COPD (ebi-a-GCST90018807) set (t1_base_harmonised_10snp.csv).
# ==========================================================================

base_dir <- REPRO_ROOT
exp <- read.csv(file.path(base_dir, "results/t1_base_harmonised_10snp.csv"),
                stringsAsFactors = FALSE)

# ---- read FinnGen raw grep output (contains substring noise) ----
hdr <- readLines(file.path(base_dir, "data_downloads/fg_header.txt"))
cols <- strsplit(sub("^#", "", hdr), "\t")[[1]]
fg <- read.delim(file.path(base_dir, "data_downloads/fg_10snps_raw.tsv"),
                 header = FALSE, col.names = cols, stringsAsFactors = FALSE)

# exact rsid match (rsids column may be comma-separated)
keep <- sapply(strsplit(fg$rsids, ","), function(x) any(x %in% exp$SNP))
fg <- fg[keep, ]
fg$SNP <- sapply(strsplit(fg$rsids, ","), function(x) x[x %in% exp$SNP][1])
stopifnot(all(exp$SNP %in% fg$SNP))
fg <- fg[match(exp$SNP, fg$SNP), ]

comp <- function(a) chartr("ACGT", "TGCA", a)

res <- exp
res$fg_chrom <- fg$chrom; res$fg_pos_b38 <- fg$pos
res$fg_ref <- fg$ref;       res$fg_alt <- fg$alt
res$fg_af_alt <- fg$af_alt
res$fg_pval_raw <- fg$pval

res$beta.outcome <- NA_real_; res$se.outcome <- fg$sebeta
res$harmonise_action <- ""; res$palindromic <- FALSE; res$ambiguous <- FALSE

for (i in seq_len(nrow(res))) {
  ea <- res$effect_allele.exposure[i]; oa <- res$other_allele.exposure[i]
  ref <- fg$ref[i]; alt <- fg$alt[i]; b <- fg$beta[i]
  pal <- (paste0(ea, oa) %in% c("AT","TA","CG","GC"))
  res$palindromic[i] <- pal
  if (ea == alt && oa == ref)      { res$beta.outcome[i] <-  b; res$harmonise_action[i] <- "aligned_alt=EA" }
  else if (ea == ref && oa == alt) { res$beta.outcome[i] <- -b; res$harmonise_action[i] <- "flipped_alt=OA" }
  else if (comp(ea) == alt && comp(oa) == ref) { res$beta.outcome[i] <-  b; res$harmonise_action[i] <- "strand_flip_aligned" }
  else if (comp(ea) == ref && comp(oa) == alt) { res$beta.outcome[i] <- -b; res$harmonise_action[i] <- "strand_flip_flipped" }
  else { res$harmonise_action[i] <- "DROP_allele_mismatch"; next }
  # palindromic check: EA freq should match alt freq when aligned
  if (pal) {
    ea_freq_in_fg <- if (grepl("flip", res$harmonise_action[i]) && !grepl("strand", res$harmonise_action[i])) 1 - fg$af_alt[i] else fg$af_alt[i]
    if (abs(ea_freq_in_fg - res$eaf.exposure[i]) > 0.15) {
      res$ambiguous[i] <- TRUE
      res$harmonise_action[i] <- paste0(res$harmonise_action[i], ";AMBIGUOUS_palindromic_freq_mismatch")
    }
  }
}
res <- res[res$harmonise_action != "DROP_allele_mismatch" & !res$ambiguous, ]
cat("SNPs retained after harmonisation:", nrow(res), "of", nrow(exp), "\n")

# ---- Wald ratios + IVW ----
res$wald <- res$beta.outcome / res$beta.exposure
res$wald_se <- res$se.outcome / abs(res$beta.exposure)

ivw_fixed <- function(b, se) {
  w <- 1/se^2; est <- sum(w*b)/sum(w); se_est <- sqrt(1/sum(w))
  c(est = est, se = se_est, p = 2*pnorm(-abs(est/se_est)))
}
ivw_random <- function(b, se) {
  w <- 1/se^2; est_f <- sum(w*b)/sum(w)
  Q <- sum(w*(b-est_f)^2); df <- length(b)-1
  tau2 <- max(0, (Q-df)/(sum(w) - sum(w^2)/sum(w)))
  w2 <- 1/(se^2+tau2); est <- sum(w2*b)/sum(w2); se_est <- sqrt(1/sum(w2))
  c(est = est, se = se_est, p = 2*pnorm(-abs(est/se_est)), Q = Q, Q_df = df,
    Q_p = pchisq(Q, df, lower.tail = FALSE), I2 = max(0,(Q-df)/Q))
}
f <- ivw_fixed(res$wald, res$wald_se)
r <- ivw_random(res$wald, res$wald_se)

# MR-Egger (orientation: exposure betas positive)
sgn <- sign(res$beta.exposure)
bx <- abs(res$beta.exposure); by <- res$beta.outcome*sgn
egger <- summary(lm(by ~ bx, weights = 1/res$se.outcome^2))
eg_int <- coef(egger)[1, ]; eg_slo <- coef(egger)[2, ]

mk <- function(method, est, se, p, extra = list()) {
  data.frame(exposure = "COPD (ebi-a-GCST90018807, 10-SNP local instruments)",
             outcome = "FinnGen R12 I9_HEARTFAIL (Heart failure, strict; 37,653 cases / 462,695 controls)",
             note = "FinnGen R12 has no HFpEF/HFrEF-specific endpoint; strict-HF used as stratified outcome",
             method = method, nsnp = nrow(res), beta = est, se = se,
             or = exp(est), or_lci95 = exp(est-1.96*se), or_uci95 = exp(est+1.96*se),
             pval = p, stringsAsFactors = FALSE)
}
out <- rbind(
  mk("Wald-ratio IVW (fixed effect)", f["est"], f["se"], f["p"]),
  mk("Wald-ratio IVW (random effects, DL)", r["est"], r["se"], r["p"]),
  mk("MR-Egger slope", eg_slo["Estimate"], eg_slo["Std. Error"], eg_slo["Pr(>|t|)"])
)
het <- data.frame(statistic = c("Cochran_Q","Q_df","Q_pval","I2"),
                  value = c(r["Q"], r["Q_df"], r["Q_p"], r["I2"]))
egger_int <- data.frame(statistic = c("Egger_intercept","Egger_intercept_se","Egger_intercept_p"),
                        value = c(eg_int["Estimate"], eg_int["Std. Error"], eg_int["Pr(>|t|)"]))

write.csv(out, file.path(base_dir, "results/t1_hfpef_hfref_mr.csv"), row.names = FALSE)
write.csv(het, file.path(base_dir, "results/t1_hfpef_hfref_heterogeneity.csv"), row.names = FALSE)
write.csv(egger_int, file.path(base_dir, "results/t1_hfpef_hfref_egger_intercept.csv"), row.names = FALSE)
write.csv(res[, c("SNP","fg_chrom","fg_pos_b38","effect_allele.exposure","other_allele.exposure",
                  "fg_ref","fg_alt","eaf.exposure","fg_af_alt","beta.exposure","se.exposure",
                  "beta.outcome","se.outcome","wald","wald_se","palindromic","harmonise_action")],
          file.path(base_dir, "results/t1_hfpef_hfref_harmonised_snps.csv"), row.names = FALSE)

cat("\n===== RESULTS =====\n"); print(out)
cat("\nHeterogeneity:\n"); print(het)
cat("\nEgger intercept:\n"); print(egger_int)
cat("\nPer-SNP:\n"); print(res[, c("SNP","beta.exposure","beta.outcome","se.outcome","wald","harmonise_action")])
