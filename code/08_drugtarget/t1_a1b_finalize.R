# ============================================================================
# Module 08 - Drug target: eQTLGen whole-blood finaliser
# eQTLGen (N=31,684) SERPINE1 cis instrument at P<5e-6 -> HF Wald ratio;
# appends the whole-blood row to the drug-target MR table.
# Input : data_downloads/serpine1_eqtla_region.rds, tmp_hf_SERPINE1.tsv
# Output: results/t1_drugtarget_mr_serpine1_7tissues.csv (updated)
# ============================================================================

# --- repository paths -----------------------------------------------------------
# Run this script from the repository root, or set REPRO_ROOT to the repo root.
REPRO_ROOT <- Sys.getenv("REPRO_ROOT", unset = ".")
# GBD_RAW_DIR: folder with the raw IHME GBD 2021 csv downloads (see data/README.md)
GBD_RAW_DIR <- Sys.getenv("GBD_RAW_DIR", unset = file.path(REPRO_ROOT, "data", "raw", "gbd_2021"))
# V8_DIR: legacy V8 project tree (NOT shipped; only needed by v8_* legacy scripts)
V8_DIR <- Sys.getenv("V8_DIR", unset = "v8_legacy")
# .libPaths(...)  # removed user-local library path; see r_packages.txt
DL <- file.path(REPRO_ROOT, "data_downloads")
OUT <- file.path(REPRO_ROOT, "results")
# eqtl-a (eQTLGen whole blood, N=31,684) SERPINE1 cis instrument at P<5e-6 -> HF (local)
reg <- readRDS(paste0(DL, "/serpine1_eqtla_region.rds"))
cis <- reg[!is.na(reg$p) & reg$p < 5e-6, ]
cis <- cis[order(cis$p), ]
sel <- numeric(0); keep <- rep(TRUE, nrow(cis))
for (i in seq_len(nrow(cis))) { if (any(abs(sel - cis$position[i]) < 1e6)) keep[i] <- FALSE else sel <- c(sel, cis$position[i]) }
cis <- cis[keep, ]
hf <- read.delim(paste0(DL, "/tmp_hf_SERPINE1.tsv"), stringsAsFactors = FALSE)
comp <- function(a) chartr("ACGT", "TGCA", a)
m <- match(cis$rsid, hf$rsid)
by <- rep(NA_real_, nrow(cis))
same <- toupper(hf$ea[m]) == toupper(cis$ea) & toupper(hf$nea[m]) == toupper(cis$nea)
flip <- toupper(hf$ea[m]) == toupper(cis$nea) & toupper(hf$nea[m]) == toupper(cis$ea)
by[same] <- hf$beta[m][same]; by[flip] <- -hf$beta[m][flip]
ok <- !is.na(by)
b <- by[ok]/cis$beta[ok]; se <- hf$se[m][ok]/abs(cis$beta[ok])
p <- 2*pnorm(-abs(b/se))
row <- data.frame(tissue = "whole_blood_eQTLGen", eqtl_dataset = "eqtl-a-ENSG00000106366",
  tissue_n = 31684, gene = "SERPINE1 (ENSG00000106366)",
  cis_window = "chr7:99,770,302-101,782,566 (GRCh37, gene+/-1Mb)",
  instrument_p_threshold = 5e-6,
  threshold_note = "pre-specified threshold (eQTLGen blood only tissue passing P<5e-6)",
  n_instruments = sum(ok), method = "Wald ratio", beta = b, se = se, or = exp(b),
  or_lci95 = exp(b-1.96*se), or_uci95 = exp(b+1.96*se), pval = p,
  lead_snp = cis$rsid[ok], lead_p_eqtl = cis$p[ok])
main <- read.csv(paste0(OUT, "/t1_drugtarget_mr_serpine1_7tissues.csv"), stringsAsFactors = FALSE)
main$lead_snp <- NA; main$lead_p_eqtl <- NA
final <- rbind(row, main)
# append the pre-specified-threshold null record
nullrows <- data.frame(
  tissue = c("whole_blood(GTEx)","lung","aorta","coronary_artery","atrial_appendage",
             "left_ventricle","subcutaneous_adipose","liver"),
  eqtl_dataset = c("QTD000356","QTD000271","QTD000131","QTD000136","QTD000251","QTD000256","QTD000116","QTD000266"),
  tissue_n = c(670,510,387,213,372,382,581,208),
  gene = "SERPINE1 (ENSG00000106366)",
  cis_window = "chr7:100,126,167-101,138,431 (GRCh38, gene+/-1Mb)",
  instrument_p_threshold = 5e-6,
  threshold_note = "PRIMARY pre-specified threshold: NO cis-eQTL instrument in ANY GTEx v8 tissue (min p: aorta 2.4e-5, blood 6.8e-4, lung 6.9e-3) - genuine null",
  n_instruments = 0, method = "none (no instrument)", beta = NA, se = NA, or = NA,
  or_lci95 = NA, or_uci95 = NA, pval = NA, lead_snp = NA, lead_p_eqtl = NA)
final <- rbind(nullrows, final)
write.csv(final, paste0(OUT, "/t1_drugtarget_mr_serpine1_7tissues.csv"), row.names = FALSE)
cat("final rows:", nrow(final), "\n")
print(final[, c("tissue","instrument_p_threshold","n_instruments","or","or_lci95","or_uci95","pval")])
