# ============================================================================
# Module 03 - Mediation MR, step B (two-step MR + Sobel + Bonferroni)
# Two-step mediation for IL-6, CRP, TGFB1 and SERPINE1 between COPD and HF,
# plus the COPD -> IL-6 -> SERPINE1 -> HF chain.
# Mediators: IL6 ebi-a-GCST90012005; CRP ebi-a-GCST90029070; TGFB1 prot-a-2962;
# SERPINE1 eqtl-a-ENSG00000106366 (eQTLGen whole blood).
# Inputs: results/t1_base_harmonised_10snp.csv; data_downloads caches (shipped
#         in data/derived/mediation/); local GWAS files for IL6/CRP/HF
#         (GCST90012005/GCST90029070/GCST009541, see data/README.md).
# Output: results/t1_mediation_results.csv (CRP: 4.4% mediated),
#         t1_mediation_chained.csv, t1_a1a_harmonised_*_to_hf.csv,
#         t1_a1b_serpine1_eqtla_cis_instruments.csv
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
LOG <- file(file.path(REPRO_ROOT, "logs/t1_a1a_mediation.log"), open = "wt"); sink(LOG, split = TRUE)
cat("mediation stepB @", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Mediators: IL6=ebi-a-GCST90012005 (Folkersen 2020 SCALLOP N=21758, largest IL-6);\n")
cat("CRP=ebi-a-GCST90029070 (Said 2022 N=575531, largest European CRP);\n")
cat("TGFB1=prot-a-2962 (Sun 2018 INTERVAL TGF-beta-1 N=3301);\n")
cat("SERPINE1=eqtl-a-ENSG00000106366 (eQTLGen whole blood, Vosa 2018 N=31684)\n")
cat("  [substitute for prot-c-2925_9_1 pQTL: that study contains <3 of the 10 COPD\n")
cat("   instruments, so the a-path is unidentifiable there; eQTLGen whole-blood\n")
cat("   cis-eQTL used as SERPINE1 instrument source, per task fallback rule]\n")

snps10 <- readRDS(paste0(DL, "/snps10.rds"))
copd <- read.csv(file.path(REPRO_ROOT, "results/t1_base_harmonised_10snp.csv"),
                 stringsAsFactors = FALSE)
bx <- copd$beta.exposure

# ---- total effect (COPD -> HF, local) ----
hf10 <- read.delim(paste0(DL, "/ext_med_HF10.tsv"), stringsAsFactors = FALSE)
m <- match(copd$SNP, hf10$rsid)
by <- align_beta(copd$effect_allele.exposure, copd$other_allele.exposure,
                 hf10$ea[m], hf10$nea[m], hf10$beta[m])
ok <- !is.na(by) & !is.na(hf10$se[m])
tot <- ivw(bx[ok], by[ok], hf10$se[m][ok])
cat("total effect COPD->HF (local IVW, n=", sum(ok), "):", tot["b"], "se", tot["se"], "\n")

# ---- a-paths ----
med_files <- c(IL6 = "GCST90012005_IL6.h.tsv.gz", CRP = "GCST90029070_CRP.h.tsv.gz")
get_med <- function(snps, tag) {
  writeLines(snps, paste0(DL, "/tmp_iv_", tag, ".txt"))
  system2("bash", c(file.path(REPRO_ROOT, "code", "02_mvmr", "extract_sumstats.sh"),
                    shQuote(paste0(DL, "/", med_files[[tag]])), "rsids",
                    paste0(DL, "/tmp_iv_", tag, ".txt"), paste0(DL, "/tmp_med_", tag, ".tsv")))
  read.delim(paste0(DL, "/tmp_med_", tag, ".tsv"), stringsAsFactors = FALSE)
}
dist_clump_df <- function(d, window = 1e6) {
  d <- d[!duplicated(d$rsid) & !is.na(d$p) & !is.na(d$chr) & !is.na(d$pos), ]
  d <- d[!(d$chr == 6 & d$pos > 25e6 & d$pos < 34e6), ]
  d <- d[order(d$p), ]
  keep <- rep(TRUE, nrow(d))
  for (ch in unique(d$chr)) {
    idx <- which(d$chr == ch); selpos <- numeric(0)
    for (i in idx) { if (any(abs(selpos - d$pos[i]) < window)) keep[i] <- FALSE else selpos <- c(selpos, d$pos[i]) }
  }
  d[keep, ]
}
a_df <- data.frame()
for (mn in c("IL6","CRP")) {
  id <- c(IL6 = "ebi-a-GCST90012005", CRP = "ebi-a-GCST90029070")[[mn]]
  a <- get_med(snps10, mn)
  mm <- match(copd$SNP, a$rsid)
  bm <- align_beta(copd$effect_allele.exposure, copd$other_allele.exposure,
                   a$ea[mm], a$nea[mm], a$beta[mm])
  ok <- !is.na(bm) & !is.na(a$se[mm])
  ra <- ivw(bx[ok], bm[ok], a$se[mm][ok])
  a_df <- rbind(a_df, data.frame(mediator = mn, mediator_id = id, a = ra["b"], a_se = ra["se"],
                                 a_p = 2*pnorm(-abs(ra["b"]/ra["se"])), a_nsnp = sum(ok)))
  cat(mn, "a-path:", ra["b"], "se", ra["se"], "nsnp", sum(ok), "\n")
}
# TGFB1 a-path (prot-a-2962, 1 API request with retries)
{
  cf_ <- paste0(DL, "/apath_TGFB1.rds")
  if (file.exists(cf_)) { a <- readRDS(cf_) } else {
    a <- NULL
    for (att in 1:4) {
      a <- tryCatch(ieugwasr::associations(snps10, "prot-a-2962", proxies = 0),
                    error = function(e) { cat("TGFB1 attempt", att, "failed\n"); NULL })
      if (!is.null(a)) break
      Sys.sleep(15)
    }
    if (!is.null(a)) saveRDS(a, cf_)
  }
  if (!is.null(a)) {
    mm <- match(copd$SNP, a$rsid)
    bm <- align_beta(copd$effect_allele.exposure, copd$other_allele.exposure,
                     a$ea[mm], a$nea[mm], a$beta[mm])
    ok <- !is.na(bm) & !is.na(a$se[mm])
    ra <- ivw(bx[ok], bm[ok], a$se[mm][ok])
    a_df <- rbind(a_df, data.frame(mediator = "TGFB1", mediator_id = "prot-a-2962", a = ra["b"], a_se = ra["se"],
                                   a_p = 2*pnorm(-abs(ra["b"]/ra["se"])), a_nsnp = sum(ok)))
    cat("TGFB1 a-path:", ra["b"], "se", ra["se"], "nsnp", sum(ok), "\n")
  } else cat("TGFB1 a-path unavailable\n")
}
# SERPINE1 a-path: 10 COPD SNPs in the FULL eqtl-a dataset (1 API request)
{
  cf_ <- paste0(DL, "/apath_SERPINE1.rds")
  if (file.exists(cf_)) { a <- readRDS(cf_) } else {
    a <- NULL
    for (att in 1:4) {
      a <- tryCatch(ieugwasr::associations(snps10, "eqtl-a-ENSG00000106366", proxies = 0),
                    error = function(e) { cat("SERPINE1 a-path attempt", att, "failed\n"); NULL })
      if (!is.null(a)) break
      Sys.sleep(15)
    }
    if (!is.null(a)) saveRDS(a, cf_)
  }
  if (!is.null(a)) {
    mm <- match(copd$SNP, a$rsid)
    bm <- align_beta(copd$effect_allele.exposure, copd$other_allele.exposure,
                     a$ea[mm], a$nea[mm], a$beta[mm])
    ok <- !is.na(bm) & !is.na(a$se[mm])
    ra <- ivw(bx[ok], bm[ok], a$se[mm][ok])
    a_df <- rbind(a_df, data.frame(mediator = "SERPINE1", mediator_id = "eqtl-a-ENSG00000106366",
                                   a = ra["b"], a_se = ra["se"],
                                   a_p = 2*pnorm(-abs(ra["b"]/ra["se"])), a_nsnp = sum(ok)))
    cat("SERPINE1(eqtl-a) a-path:", ra["b"], "se", ra["se"], "nsnp", sum(ok), "\n")
  } else cat("SERPINE1 a-path unavailable\n")
}

# ---- b-paths: instruments -> local HF ----
reg <- readRDS(paste0(DL, "/serpine1_eqtla_region.rds"))   # eqtl-a SERPINE1 cis region (7,776 rows)
get_hf <- function(snps, tag) {
  writeLines(snps, paste0(DL, "/tmp_iv_", tag, ".txt"))
  system2("bash", c(file.path(REPRO_ROOT, "code", "02_mvmr", "extract_sumstats.sh"),
                    shQuote(paste0(DL, "/GCST009541_HF.h.tsv.gz")), "rsids",
                    paste0(DL, "/tmp_iv_", tag, ".txt"), paste0(DL, "/tmp_hf_", tag, ".tsv")))
  read.delim(paste0(DL, "/tmp_hf_", tag, ".tsv"), stringsAsFactors = FALSE)
}
b_df <- data.frame()
# IL6 / CRP instruments from LOCAL files (p<5e-8 + 1Mb distance clump)
for (mn in c("IL6","CRP")) {
  id <- c(IL6 = "ebi-a-GCST90012005", CRP = "ebi-a-GCST90029070")[[mn]]
  cf_ <- paste0(DL, "/iv_local_", mn, ".rds")
  if (file.exists(cf_)) { ex <- readRDS(cf_) } else {
    system2("bash", c(file.path(REPRO_ROOT, "code", "02_mvmr", "extract_sumstats.sh"),
                      shQuote(paste0(DL, "/", med_files[[mn]])), "p5e8", "-",
                      paste0(DL, "/tmp_p5e8_", mn, ".tsv")))
    raw <- read.delim(paste0(DL, "/tmp_p5e8_", mn, ".tsv"), stringsAsFactors = FALSE)
    dc <- dist_clump_df(raw)
    ex <- data.frame(SNP = dc$rsid, effect_allele.exposure = dc$ea, other_allele.exposure = dc$nea,
                     beta.exposure = dc$beta, se.exposure = dc$se, pval.exposure = dc$p,
                     eaf.exposure = dc$eaf, stringsAsFactors = FALSE)
    saveRDS(ex, cf_)
    cat(mn, "local instruments:", nrow(raw), "->", nrow(ex), "after clump\n")
  }
  hf <- get_hf(ex$SNP, mn)
  mm <- match(ex$SNP, hf$rsid); kp <- !is.na(mm)
  bby <- align_beta(ex$effect_allele.exposure[kp], ex$other_allele.exposure[kp],
                    hf$ea[mm[kp]], hf$nea[mm[kp]], hf$beta[mm[kp]])
  sey <- hf$se[mm[kp]]; bxx <- ex$beta.exposure[kp]
  ok <- !is.na(bby) & !is.na(sey)
  snps_ok <- ex$SNP[kp][ok]
  rb <- if (sum(ok) == 1) c(b = bby[ok]/bxx[ok], se = sey[ok]/abs(bxx[ok])) else ivw(bxx[ok], bby[ok], sey[ok])
  write.csv(data.frame(SNP = snps_ok, beta_exp = bxx[ok], beta_hf = bby[ok], se_hf = sey[ok]),
            paste0(OUT, "/t1_a1a_harmonised_", mn, "_to_hf.csv"), row.names = FALSE)
  b_df <- rbind(b_df, data.frame(mediator = mn, b_mediator_id = id, b = rb["b"], b_se = rb["se"],
                                 b_p = 2*pnorm(-abs(rb["b"]/rb["se"])), b_nsnp = sum(ok)))
  cat(mn, "b-path:", rb["b"], "se", rb["se"], "nsnp", sum(ok), "\n")
}
# TGFB1 b-path (cached 10 instruments)
ex <- readRDS(paste0(DL, "/iv_med_TGFB1.rds"))
hf <- get_hf(ex$SNP, "TGFB1")
mm <- match(ex$SNP, hf$rsid); kp <- !is.na(mm)
bby <- align_beta(ex$effect_allele.exposure[kp], ex$other_allele.exposure[kp],
                  hf$ea[mm[kp]], hf$nea[mm[kp]], hf$beta[mm[kp]])
sey <- hf$se[mm[kp]]; bxx <- ex$beta.exposure[kp]
ok <- !is.na(bby) & !is.na(sey)
rb <- if (sum(ok) == 1) c(b = bby[ok]/bxx[ok], se = sey[ok]/abs(bxx[ok])) else ivw(bxx[ok], bby[ok], sey[ok])
b_df <- rbind(b_df, data.frame(mediator = "TGFB1", b_mediator_id = "prot-a-2962", b = rb["b"], b_se = rb["se"],
                               b_p = 2*pnorm(-abs(rb["b"]/rb["se"])), b_nsnp = sum(ok)))
cat("TGFB1 b-path:", rb["b"], "se", rb["se"], "nsnp", sum(ok), "\n")
# SERPINE1 b-path: cis instruments from eqtl-a region (P<5e-6, 1Mb distance clump)
cis <- reg[!is.na(reg$p) & reg$p < 5e-6, ]
cis <- cis[order(cis$p), ]
selp <- numeric(0); keep <- rep(TRUE, nrow(cis))
for (i in seq_len(nrow(cis))) {
  if (any(abs(selp - cis$position[i]) < 1e6)) keep[i] <- FALSE else selp <- c(selp, cis$position[i])
}
cis <- cis[keep, ]
cat("SERPINE1 cis instruments (p<5e-6, 1Mb):", nrow(cis), "\n")
write.csv(cis, paste0(OUT, "/t1_a1b_serpine1_eqtla_cis_instruments.csv"), row.names = FALSE)
if (nrow(cis) >= 1) {
  hf <- get_hf(cis$rsid, "SERPINE1")
  mm <- match(cis$rsid, hf$rsid); kp <- !is.na(mm)
  bby <- align_beta(cis$ea[kp], cis$nea[kp], hf$ea[mm[kp]], hf$nea[mm[kp]], hf$beta[mm[kp]])
  sey <- hf$se[mm[kp]]; bxx <- cis$beta[kp]
  ok <- !is.na(bby) & !is.na(sey)
  rb <- if (sum(ok) == 1) c(b = bby[ok]/bxx[ok], se = sey[ok]/abs(bxx[ok])) else ivw(bxx[ok], bby[ok], sey[ok])
  b_df <- rbind(b_df, data.frame(mediator = "SERPINE1", b_mediator_id = "eqtl-a-ENSG00000106366",
                                 b = rb["b"], b_se = rb["se"],
                                 b_p = 2*pnorm(-abs(rb["b"]/rb["se"])), b_nsnp = sum(ok)))
  cat("SERPINE1(eqtl-a) b-path:", rb["b"], "se", rb["se"], "nsnp", sum(ok), "\n")
}

# ---- combine ----
med <- merge(a_df, b_df, by = "mediator", all = TRUE)
med$indirect <- med$a * med$b
med$indirect_se <- sqrt(med$b^2 * med$a_se^2 + med$a^2 * med$b_se^2)
med$indirect_z <- med$indirect / med$indirect_se
med$indirect_p <- 2 * pnorm(-abs(med$indirect_z))
med$bonferroni_sig_4tests <- med$indirect_p < 0.05/4
med$total_effect <- tot["b"]; med$total_se <- tot["se"]
med$prop_mediated <- med$indirect / tot["b"]
write.csv(med, paste0(OUT, "/t1_mediation_results.csv"), row.names = FALSE)
cat("\n=== mediation results ===\n"); print(med[, c("mediator","a","a_se","a_p","b","b_se","b_p","indirect","indirect_se","indirect_p","bonferroni_sig_4tests","prop_mediated")])

# ---- chain: COPD -> IL-6 -> SERPINE1(eqtl-a) -> HF ----
ex_il6 <- readRDS(paste0(DL, "/iv_local_IL6.rds"))
mm <- match(ex_il6$SNP, reg$rsid)
a2_beta <- align_beta(ex_il6$effect_allele.exposure, ex_il6$other_allele.exposure,
                      reg$ea[mm], reg$nea[mm], reg$beta[mm])
ok <- !is.na(a2_beta) & !is.na(reg$se[mm])
if (sum(ok) >= 1) {
  r2 <- if (sum(ok) == 1) c(b = a2_beta[ok]/ex_il6$beta.exposure[ok], se = reg$se[mm][ok]/abs(ex_il6$beta.exposure[ok])) else ivw(ex_il6$beta.exposure[ok], a2_beta[ok], reg$se[mm][ok])
  a1 <- med$a[med$mediator == "IL6"]; a1_se <- med$a_se[med$mediator == "IL6"]
  bb <- med$b[med$mediator == "SERPINE1"]; bb_se <- med$b_se[med$mediator == "SERPINE1"]
  a2 <- r2["b"]; a2_se <- r2["se"]
  ind <- a1 * a2 * bb
  se <- sqrt((a2*bb*a1_se)^2 + (a1*bb*a2_se)^2 + (a1*a2*bb_se)^2)
  z <- ind/se
  chain <- data.frame(chain = "COPD -> IL-6 -> SERPINE1(whole-blood eQTL) -> HF",
    a1_COPD_IL6 = a1, a1_se = a1_se, a2_IL6_SERPINE1 = a2, a2_se = a2_se, a2_nsnp = sum(ok),
    b_SERPINE1_HF = bb, b_se = bb_se, indirect = ind, indirect_se = se, indirect_z = z,
    indirect_p = 2*pnorm(-abs(z)), prop_of_total = ind/tot["b"])
  write.csv(chain, paste0(OUT, "/t1_mediation_chained.csv"), row.names = FALSE)
  cat("\n=== chain ===\n"); print(chain)
} else cat("chain: IL-6 instruments absent from eqtl-a region; chain not identifiable\n")

sink(); close(LOG)
cat("MEDIATION STEP B DONE\n")
