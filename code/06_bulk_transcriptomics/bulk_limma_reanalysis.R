# ============================================================================
# Module 06 - Bulk transcriptomics: independent GEOquery + limma re-analysis
#
# Re-downloads the two bulk cohorts from GEO and re-runs the differential
# expression analysis that originally produced the shared 41-gene
# inflammation/fibrosis program:
#   * GSE57148  COPD lung tissue        (95 COPD vs 89 control, GPL11154)
#   * GSE57338  failing myocardium      (177 HF vs 136 non-failing, GPL11532)
#
# For each cohort: log2-scale expression matrix -> limma eBayes two-group
# comparison -> FDR (BH). The shared program is the intersection of genes
# significant (FDR < 0.05, same direction) in BOTH cohorts, intersected with
# a prespecified inflammation/fibrosis candidate list
# (data/derived/bulk/gene_list_41_inflammation_fibrosis.csv).
#
# NOTE: this script is provided for full-pipeline reproducibility and was NOT
# executed as part of the shipped validation (network access to GEO is
# required and GEO SOFT matrices are ~100 MB each). The key-gene group
# statistics it reproduces are independently verified from the shipped
# investigator-curated expression workbooks by bulk_keygene_welch.py
# (see VALIDATION.md).
#
# Inputs : GEO (network); data/derived/bulk/gene_list_41_inflammation_fibrosis.csv
# Outputs: results/bulk_limma_gse57148_de.csv, results/bulk_limma_gse57338_de.csv,
#          results/bulk_limma_shared41_check.csv
# ============================================================================

# --- repository paths -----------------------------------------------------------
# Run this script from the repository root, or set REPRO_ROOT to the repo root.
REPRO_ROOT <- Sys.getenv("REPRO_ROOT", unset = ".")

library(GEOquery)
library(Biobase)
library(limma)

OUT <- file.path(REPRO_ROOT, "results")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

GENES41 <- read.csv(file.path(REPRO_ROOT, "data", "derived", "bulk",
                              "gene_list_41_inflammation_fibrosis.csv"),
                    header = FALSE)[, 1]

# --- helper: one-cohort limma two-group analysis --------------------------------
run_cohort <- function(gse_id, case_label, case_grep, platform_gpl) {
  message("== ", gse_id, " ==")
  gset <- getGEO(gse_id, GSEMatrix = TRUE, getGPL = TRUE)
  eset <- gset[[1]]
  # group labels from sample titles / source_name / characteristics
  meta <- pData(eset)
  grp_txt <- apply(meta[, c("title", "source_name_ch1"), drop = FALSE], 1,
                   paste, collapse = " ")
  group <- ifelse(grepl(case_grep, grp_txt, ignore.case = TRUE),
                  case_label, "control")
  group <- factor(group, levels = c("control", case_label))
  message("groups: ", paste(names(table(group)), "=", as.integer(table(group)),
                            collapse = ", "))
  expr <- exprs(eset)
  if (max(expr, na.rm = TRUE) > 100) expr <- log2(expr + 1)   # ensure log2 scale
  # gene symbols: prefer the platform "Gene Symbol" annotation column
  fd <- fData(eset)
  sym_col <- grep("gene.?symbol|^symbol$", colnames(fd), ignore.case = TRUE,
                  value = TRUE)[1]
  stopifnot(!is.na(sym_col))
  symbols <- sub(" ///.*", "", as.character(fd[[sym_col]]))
  keep <- !is.na(symbols) & symbols != ""
  expr <- expr[keep, , drop = FALSE]; symbols <- symbols[keep]
  # collapse duplicated symbols by max mean expression
  expr <- rowsum(expr, group = symbols, reorder = FALSE) /
          as.vector(table(symbols)[rownames(rowsum(expr, group = symbols))])
  design <- model.matrix(~ group)
  fit <- eBayes(lmFit(expr, design))
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  tt$gene <- rownames(tt)
  tt <- tt[, c("gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]
  write.csv(tt, file.path(OUT, paste0("bulk_limma_",
          tolower(gse_id), "_de.csv")), row.names = FALSE)
  tt
}

de_lung  <- run_cohort("GSE57148", "COPD", "COPD", "GPL11154")
de_heart <- run_cohort("GSE57338", "HF",   "fail|cardiomyopathy|HF",
                       "GPL11532")

# --- shared 41-gene check ---------------------------------------------------------
m <- merge(de_lung[, c("gene", "logFC", "adj.P.Val")],
           de_heart[, c("gene", "adj.P.Val")], by = "gene",
           suffixes = c("_lung", "_heart"))
chk <- m[m$gene %in% GENES41, ]
chk$concordant <- sign(chk$logFC_lung) == sign(chk$logFC_heart)
chk <- chk[order(chk$adj.P.Val_lung), ]
write.csv(chk, file.path(OUT, "bulk_limma_shared41_check.csv"),
          row.names = FALSE)
message("41-gene program recovered: ", sum(chk$concordant, na.rm = TRUE),
        " / ", nrow(chk), " concordant")
message("DONE")
