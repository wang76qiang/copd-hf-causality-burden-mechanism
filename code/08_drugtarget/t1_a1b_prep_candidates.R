# ============================================================================
# Module 08 - Drug target: SERPINE1 cis-eQTL candidate preparation
# Filter eQTL Catalogue GTEx v8 region files (8 tissues) to SERPINE1
# (ENSG00000106366), pick candidate instruments (P<5e-6), save per-tissue
# gene stats for coloc.
# Input : data_downloads/eqtl_<tissue>_serpine1.tsv (remote-tabix output; NOT
#         shipped - regenerate with t1_minitabix.py, see PIPELINE.md)
# Output: data_downloads/serpine1_gene/<tissue>.csv (shipped in
#         data/derived/drugtarget/), t1_a1b_candidates.rds
# ============================================================================

# --- repository paths -----------------------------------------------------------
# Run this script from the repository root, or set REPRO_ROOT to the repo root.
REPRO_ROOT <- Sys.getenv("REPRO_ROOT", unset = ".")
# GBD_RAW_DIR: folder with the raw IHME GBD 2021 csv downloads (see data/README.md)
GBD_RAW_DIR <- Sys.getenv("GBD_RAW_DIR", unset = file.path(REPRO_ROOT, "data", "raw", "gbd_2021"))
# V8_DIR: legacy V8 project tree (NOT shipped; only needed by v8_* legacy scripts)
V8_DIR <- Sys.getenv("V8_DIR", unset = "v8_legacy")
# .libPaths(...)  # removed user-local library path; see r_packages.txt
# A1b prep: filter eQTL Catalogue GTEx v8 region files to SERPINE1 (ENSG00000106366),
# pick candidate instruments (P<5e-6), save per-tissue gene stats for coloc.
dl <- file.path(REPRO_ROOT, "data_downloads")
tissues <- c(lung = "lung", blood = "blood", artery_aorta = "artery_aorta",
             artery_coronary = "artery_coronary", atrial_appendage = "atrial_appendage",
             left_ventricle = "left_ventricle", adipose_subcut = "adipose_subcut",
             liver = "liver")
dir.create(paste0(dl, "/serpine1_gene"), showWarnings = FALSE)
cand <- list()
for (t in names(tissues)) {
  f <- paste0(dl, "/eqtl_", t, "_serpine1.tsv")
  d <- read.delim(f, stringsAsFactors = FALSE)
  g <- d[d$gene_id == "ENSG00000106366" & d$type == "SNP" & !is.na(d$rsid) & d$rsid != "", ]
  g <- g[!duplicated(paste(g$variant, g$gene_id)), ]
  g$maf[g$maf == 0 | is.na(g$maf)] <- 0.001
  write.csv(g, paste0(dl, "/serpine1_gene/", t, ".csv"), row.names = FALSE)
  c5 <- g[g$pvalue < 5e-6, ]
  cand[[t]] <- c5
  cat(sprintf("%-18s N=%4d  cis SNPs=%5d  P<5e-6: %3d  minP=%.2e (lead %s)\n",
              t, c5$an[1], nrow(g), nrow(c5), min(g$pvalue),
              g$rsid[which.min(g$pvalue)]))
}
saveRDS(cand, paste0(dl, "/t1_a1b_candidates.rds"))
cat("\nunion candidate SNPs:", length(unique(unlist(lapply(cand, function(x) x$rsid)))), "\n")
