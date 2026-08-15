# ============================================================================
# Module 08 - Drug target: 8-tissue cis-eQTL MR + colocalisation
# PRIMARY (prespecified): P<5e-6, r2<0.1 instruments -> 0 SNPs in ALL 8 GTEx
# v8 tissues, so IVW is not estimable. EXPLORATORY: per-tissue lead cis-SNP
# Wald ratio vs HF. COLOCALISATION: manual reimplementation of coloc.abf
# (Giambartolomei et al. 2014; W_eQTL=0.15, W_HF=0.2, p1=p2=1e-4, p12=1e-5).
# Input : data_downloads/serpine1_gene/*.csv, t1_a1b_hf_{lead,region}.rds
# Output: results/t1_drugtarget_mr_serpine1_7tissues.csv, t1_coloc_serpine1.csv
# ============================================================================

# --- repository paths -----------------------------------------------------------
# Run this script from the repository root, or set REPRO_ROOT to the repo root.
REPRO_ROOT <- Sys.getenv("REPRO_ROOT", unset = ".")
# GBD_RAW_DIR: folder with the raw IHME GBD 2021 csv downloads (see data/README.md)
GBD_RAW_DIR <- Sys.getenv("GBD_RAW_DIR", unset = file.path(REPRO_ROOT, "data", "raw", "gbd_2021"))
# V8_DIR: legacy V8 project tree (NOT shipped; only needed by v8_* legacy scripts)
V8_DIR <- Sys.getenv("V8_DIR", unset = "v8_legacy")
# .libPaths(...)  # removed user-local library path; see r_packages.txt
# A1b analysis (2026-07-26):
#  PRIMARY (prespecified): P<5e-6 r2<0.1 instruments -> 0 SNPs in ALL 8 GTEx v8
#    tissues (see logs/t1_a1b_serpine1_drugtarget.log) -> IVW not estimable.
#  EXPLORATORY (documented deviation): per-tissue lead cis-SNP Wald ratio vs HF.
#  COLOCALISATION: manual reimplementation of coloc.abf (Giambartolomei et al.
#    2014; coloc package unavailable: RcppArmadillo build failure on R 4.6.1).
#    W_eQTL=0.15 (quant), W_HF=0.2 (cc), p1=p2=1e-4, p12=1e-5.
dl <- file.path(REPRO_ROOT, "data_downloads")
OUT <- file.path(REPRO_ROOT, "results")
LOG <- file(file.path(REPRO_ROOT, "logs/t1_a1b_serpine1_drugtarget.log"), open = "wt")
sink(LOG, split = TRUE)
cat("A1b analysis @", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Data: eQTL Catalogue GTEx v8 (study QTS000015) per-tissue full cis summary\n")
cat("stats via remote tabix range-fetch (GRCh38 chr7:100,126,167-101,138,431).\n")
cat("Reason OpenGWAS not used for eQTL: catalogue snapshot contains NO GTEx\n")
cat("per-tissue eQTL datasets (only eqtl-a eQTLGen whole blood).\n\n")

tissues <- c("lung","blood","artery_aorta","artery_coronary","atrial_appendage",
             "left_ventricle","adipose_subcut","liver")
sizes <- c(lung=510, blood=670, artery_aorta=387, artery_coronary=213,
           atrial_appendage=372, left_ventricle=382, adipose_subcut=581, liver=208)

hf_lead <- as.data.frame(readRDS(paste0(dl, "/t1_a1b_hf_lead.rds")))
hf_reg  <- as.data.frame(readRDS(paste0(dl, "/t1_a1b_hf_region.rds")))

# ---------- exploratory lead-SNP Wald ratios ----------
res <- list()
for (t in tissues) {
  g <- read.csv(paste0(dl, "/serpine1_gene/", t, ".csv"), stringsAsFactors = FALSE)
  li <- which.min(g$pvalue)
  lead <- g[li, ]
  h <- hf_lead[hf_lead$rsid == lead$rsid, ]
  row <- data.frame(tissue = t, eqtl_dataset = paste0("eQTLCat GTExv8 ", t),
                    n_tissue = sizes[[t]], n_cis_snps = nrow(g),
                    min_p_eqtl = min(g$pvalue), n_instruments_p5e6 = sum(g$pvalue < 5e-6),
                    lead_snp = lead$rsid, lead_p_eqtl = lead$pvalue,
                    beta_eqtl = lead$beta, se_eqtl = lead$se)
  row$beta_hf <- NA_real_; row$se_hf <- NA_real_
  row$wald_beta <- NA_real_; row$wald_se <- NA_real_
  row$wald_or <- NA_real_; row$wald_lci <- NA_real_; row$wald_uci <- NA_real_
  row$wald_p <- NA_real_; row$note <- ""
  if (nrow(h) == 1) {
    # harmonise: eQTL effect allele = alt; HF effect allele = ea
    b_hf <- h$beta; se_hf <- h$se
    ok <- TRUE
    if (h$ea == lead$ref && h$nea == lead$alt) { b_hf <- -b_hf }
    else if (!(h$ea == lead$alt && h$nea == lead$ref)) {
      comp <- function(a) chartr("ACGT","TGCA",a)
      if (h$ea == comp(lead$alt) && h$nea == comp(lead$ref)) {}
      else if (h$ea == comp(lead$ref) && h$nea == comp(lead$alt)) { b_hf <- -b_hf }
      else ok <- FALSE
    }
    if (!ok) { row$note <- "allele_mismatch_dropped"; res[[t]] <- row; next }
    w <- b_hf / lead$beta; wse <- se_hf / abs(lead$beta)
    row$beta_hf <- b_hf; row$se_hf <- se_hf
    row$wald_beta <- w; row$wald_se <- wse
    row$wald_or <- exp(w); row$wald_lci <- exp(w - 1.96*wse); row$wald_uci <- exp(w + 1.96*wse)
    row$wald_p <- 2 * pnorm(-abs(w/wse))
  } else row$note <- "lead SNP absent from HF dataset"
  res[[t]] <- row
}
mr_df <- do.call(rbind, res)
write.csv(mr_df, file.path(OUT, "t1_drugtarget_mr_serpine1_7tissues.csv"), row.names = FALSE)
cat("=== exploratory lead-SNP Wald ratios (P<5e-6 IVW not estimable: 0 instruments) ===\n")
print(mr_df[, c("tissue","lead_snp","lead_p_eqtl","beta_eqtl","wald_or","wald_lci","wald_uci","wald_p")])

# ---------- manual coloc.abf ----------
lbf <- function(beta, se, W) {
  V <- se^2
  0.5 * log(V/(V+W)) + 0.5 * (beta^2/V) * (W/(V+W))
}
logsum <- function(x) { m <- max(x); m + log(sum(exp(x - m))) }
coloc_abf <- function(b1, se1, b2, se2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5) {
  l1 <- lbf(b1, se1, 0.15); l2 <- lbf(b2, se2, 0.20)
  S1 <- logsum(l1); S2 <- logsum(l2); S12 <- logsum(l1 + l2)
  lH0 <- 0
  lH1 <- log(p1) + S1
  lH2 <- log(p2) + S2
  # H3: sum over i!=j p1 p2 ABF1_i ABF2_j = p1p2 (S1*S2 - sum_i ABF1_i ABF2_i)
  lH3 <- log(p1) + log(p2) + log(exp(S1 + S2) - exp(S12))
  lH4 <- log(p12) + S12
  pp <- exp(c(H0 = lH0, H1 = lH1, H2 = lH2, H3 = lH3, H4 = lH4) -
            logsum(c(lH0, lH1, lH2, lH3, lH4)))
  pp
}
co_res <- list()
for (t in tissues) {
  g <- read.csv(paste0(dl, "/serpine1_gene/", t, ".csv"), stringsAsFactors = FALSE)
  m <- merge(g[, c("rsid","beta","se","maf")],
             hf_reg[, c("rsid","beta","se","eaf","ea","nea")],
             by = "rsid", suffixes = c(".eqtl",".hf"))
  if (nrow(m) < 100) { cat(t, ": too few overlapping SNPs (", nrow(m), ")\n"); next }
  co <- coloc_abf(m$beta.eqtl, m$se.eqtl, m$beta.hf, m$se.hf)
  co_res[[t]] <- data.frame(tissue = t, n_snps_overlap = nrow(m),
                            PP.H0 = co["H0"], PP.H1 = co["H1"], PP.H2 = co["H2"],
                            PP.H3 = co["H3"], PP.H4 = co["H4"],
                            note = "manual coloc.abf reimplementation; W=0.15/0.2; p1=p2=1e-4; p12=1e-5")
}
co_df <- do.call(rbind, co_res)
write.csv(co_df, file.path(OUT, "t1_coloc_serpine1.csv"), row.names = FALSE)
cat("\n=== coloc (manual ABF) ===\n"); print(co_df[, c("tissue","n_snps_overlap","PP.H3","PP.H4")])
sink(); close(LOG)
cat("\nA1b analysis DONE\n")
