# ============================================================================
# Module 08 - Drug target: FinnGen R12 PheWAS of SERPINE1 lead SNP rs7860931
# rs7860931 (b37 chr9:72,211,325 C/T) = FinnGen R12 (GRCh38) chr9:69,596,409
# T/C. PheWeb API variant lookup across 2,470 endpoints; FinnGen effect
# allele = alt = C (the PAI-1-lowering allele, beta_pQTL = -0.2738).
# Input : data/derived/finngen/fg_r12_phewas_rs7860931.json (re-fetchable
#         from https://r12.finngen.fi/api/variant/9-69596409-T-C)
# Output: results/t1_c1_finngen_r12_phewas_rs7860931{,_cardiorespiratory}.csv
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
# Task C1: SERPINE1 lead SNP rs7860931 PheWAS extension via FinnGen R12
# rs7860931: b37 chr9:72,211,325 C/T -> FinnGen R12 (GRCh38) chr9:69,596,409 T/C
# (b38 coordinates confirmed by grepping the downloaded
#  finngen_R12_I9_HEARTFAIL.gz summary statistics).
# Source: PheWeb API https://r12.finngen.fi/api/variant/9-69596409-T-C (HTTP 200)
# NOTE: FinnGen effect allele = alt = C. The pQTL instrument orientation
# (prot-c-2925_9_1) had effect allele C with beta_pQTL = -0.2738, i.e.
# C allele lowers PAI-1. Sign convention below: raw FinnGen beta wrt alt(C).
# ==========================================================================
library(jsonlite)

j <- jsonlite::fromJSON(file.path(REPRO_ROOT, "data_downloads/fg_r12_phewas_rs7860931.json"),
                        simplifyVector = TRUE)
cat("top-level keys:", paste(names(j), collapse = ", "), "\n")

pw <- j$results
if (is.null(pw)) stop("no results element; inspect JSON keys")
cat("phewas rows:", nrow(pw), "; cols:", paste(colnames(pw), collapse = ", "), "\n")

df <- as.data.frame(pw, stringsAsFactors = FALSE)
# standardise column names across PheWeb flavours
cn <- colnames(df)
pick <- function(...) { for (x in c(...)) if (x %in% cn) return(x); NA }
df_out <- data.frame(
  phenocode = df[[pick("phenocode")]],
  phenotype = df[[pick("phenostring","phenotype")]],
  category  = df[[pick("category")]],
  pval      = as.numeric(df[[pick("pval")]]),
  beta_altC = as.numeric(df[[pick("beta")]]),
  sebeta    = as.numeric(df[[pick("sebeta")]]),
  mlogp     = as.numeric(df[[pick("mlogp")]]),
  maf       = as.numeric(df[[pick("maf")]]),
  n_case    = as.numeric(df[[pick("n_case")]]),
  n_control = as.numeric(df[[pick("n_control")]]),
  stringsAsFactors = FALSE
)
n_tot <- nrow(df_out)
bonf <- 0.05 / n_tot
df_out$bonferroni_sig <- df_out$pval < bonf
df_out$fdr <- p.adjust(df_out$pval, "BH")
df_out <- df_out[order(df_out$pval), ]
df_out$rank <- seq_len(nrow(df_out))

write.csv(df_out, file.path(REPRO_ROOT, "results/t1_c1_finngen_r12_phewas_rs7860931.csv"),
          row.names = FALSE)

cat("\nn endpoints tested:", n_tot, " Bonferroni threshold:", signif(bonf, 3), "\n")
cat("n Bonferroni-significant:", sum(df_out$bonferroni_sig, na.rm = TRUE), "\n\n")
cat("=== Top 25 associations (FinnGen R12, beta wrt alt allele C) ===\n")
print(head(df_out[, c("phenocode","phenotype","pval","beta_altC","bonferroni_sig")], 25),
      row.names = FALSE)

# highlight cardiometabolic / HF-relevant endpoints
hl <- grepl("HEART|CARDI|HF|COPD|J10|EMPHY|BRONCH|THROMB|EMBOL|DIAB|OBES|HYPERTEN",
            toupper(paste(df_out$phenocode, df_out$phenotype)))
sub <- df_out[which(hl & !is.na(df_out$pval) & df_out$pval < 1e-3), ]
cat("\n=== HF/cardiometabolic/respiratory-related endpoints with p<1e-3 ===\n")
print(sub[, c("phenocode","phenotype","pval","beta_altC","bonferroni_sig")], row.names = FALSE)
write.csv(sub, file.path(REPRO_ROOT, "results/t1_c1_finngen_r12_phewas_rs7860931_cardiorespiratory.csv"),
          row.names = FALSE)
