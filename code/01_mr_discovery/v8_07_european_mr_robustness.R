# ============================================================================
# Module 01 - MR discovery (LEGACY V8 script, kept for provenance)
# Direction-aware two-sample MR (COPD -> HF and reverse) with full sensitivity
# battery: IVW, weighted median, MR-Egger, MR-RAPS, radial IVW/Egger,
# MR-PRESSO, Steiger directionality, F-stat filtering, distance clumping.
# Inputs : full GWAS summary statistics (GCST90018807 COPD; Shah 2020 HERMES
#          HF) - NOT shipped, see data/README.md.
# Outputs: V8 results tree (discovery MR estimates used in Table_S2).
# ============================================================================

# --- repository paths -----------------------------------------------------------
# Run this script from the repository root, or set REPRO_ROOT to the repo root.
REPRO_ROOT <- Sys.getenv("REPRO_ROOT", unset = ".")
# GBD_RAW_DIR: folder with the raw IHME GBD 2021 csv downloads (see data/README.md)
GBD_RAW_DIR <- Sys.getenv("GBD_RAW_DIR", unset = file.path(REPRO_ROOT, "data", "raw", "gbd_2021"))
# V8_DIR: legacy V8 project tree (NOT shipped; only needed by v8_* legacy scripts)
V8_DIR <- Sys.getenv("V8_DIR", unset = "v8_legacy")
#!/usr/bin/env Rscript
# .libPaths(...)  # removed user-local library path; see r_packages.txt
# 07_european_mr_robustness.R
# Direction-aware two-sample MR (COPD -> HF and reverse) using summary statistics.
# Uses distance clumping (500 kb), F-stat filtering, allele harmonisation,
# and multiple sensitivity estimators: IVW, weighted median, MR-Egger,
# MR-RAPS, radial IVW/Egger, MR-PRESSO, Steiger directionality.

suppressPackageStartupMessages({
  library(data.table)
  library(MendelianRandomization)
  library(RadialMR)
  library(mr.raps)
  library(MRPRESSO)
  library(TwoSampleMR)
  library(jsonlite)
  library(ggplot2)
})

# -------------------------------------------------------------------------
# Paths and parameters
# -------------------------------------------------------------------------
base_dir <- V8_DIR
copd_file <- file.path(base_dir, "data/gwas/GCST90018807_buildGRCh37.tsv.gz")
hf_file   <- file.path(base_dir, "data/gwas/ShahS_31919418_HeartFailure.gz")
out_dir   <- file.path(base_dir, "results/mr_robustness")
plot_dir  <- file.path(out_dir, "plots")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

N_COPD <- 468475L
N_HF   <- 977323L
P_THRESHOLD <- 5e-8
CLUMP_KB    <- 500
MIN_F       <- 10

# -------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------
read_gwas <- function(path, n) {
  dt <- fread(path, header = TRUE, sep = "\t", data.table = TRUE, nrows = n)

  # Detect schema using lower-case column names, then rename from the original names
  cnames <- tolower(names(dt))
  orig <- names(dt)

  match_and_rename <- function(want_lower, new_names) {
    want_lower <- tolower(want_lower)
    idx <- match(want_lower, cnames)
    if (any(is.na(idx))) stop("Columns missing: ", paste(want_lower[is.na(idx)], collapse = ", "))
    setnames(dt, orig[idx], new_names)
  }

  if (all(c("chromosome", "base_pair_location", "effect_allele", "other_allele",
            "effect_allele_frequency", "beta", "standard_error", "p_value") %in% cnames)) {
    # EBI GCST format
    match_and_rename(
      c("chromosome", "base_pair_location", "effect_allele", "other_allele",
        "effect_allele_frequency", "beta", "standard_error", "p_value", "variant_id"),
      c("chr", "pos", "ea", "nea", "eaf", "beta", "se", "pval", "rsid"))
  } else if (all(c("snp", "chr", "bp", "a1", "a2", "freq", "b", "se", "p") %in% cnames)) {
    # Shah-style compressed format
    match_and_rename(
      c("snp", "chr", "bp", "a1", "a2", "freq", "b", "se", "p"),
      c("rsid", "chr", "pos", "ea", "nea", "eaf", "beta", "se", "pval"))
  } else {
    stop("Unrecognised GWAS column layout: ", paste(cnames, collapse = ", "))
  }

  dt[, chr := suppressWarnings(as.integer(chr))]
  dt[, pos := suppressWarnings(as.integer(pos))]
  dt <- dt[!is.na(chr) & !is.na(pos) & !is.na(beta) & !is.na(se) & !is.na(pval)]
  dt <- dt[chr %in% 1:22]
  dt[, ea := toupper(ea)]
  dt[, nea := toupper(nea)]
  dt[, id := paste(chr, pos, sep = ":")]
  dt
}

# Distance-based clumping: keep the SNP with the smallest p-value in each
# 500 kb window, iterating greedily genome-wide.
distance_clump <- function(dt, kb = CLUMP_KB) {
  dt <- dt[order(pval)]
  keep <- rep(FALSE, nrow(dt))
  kept_pos <- list()  # keyed by chromosome
  window_bp <- kb * 1000
  for (i in seq_len(nrow(dt))) {
    cr <- as.character(dt$chr[i])
    ps <- dt$pos[i]
    prev <- kept_pos[[cr]]
    if (is.null(prev) || all(abs(ps - prev) > window_bp)) {
      keep[i] <- TRUE
      kept_pos[[cr]] <- c(prev, ps)
    }
  }
  dt[keep, ]
}

harmonise <- function(exp_dt, out_dt) {
  m <- merge(exp_dt, out_dt, by = c("chr", "pos"), suffixes = c(".exp", ".out"), all = FALSE)
  m <- m[!is.na(eaf.exp) & !is.na(eaf.out)]
  if (nrow(m) == 0) return(NULL)

  same <- (m$ea.exp == m$ea.out) & (m$nea.exp == m$nea.out)
  flip <- (m$ea.exp == m$nea.out) & (m$nea.exp == m$ea.out)
  keep <- same | flip
  m <- m[keep, ]

  if (nrow(m) == 0) return(NULL)

  # recompute flip/same on the kept subset
  flip <- (m$ea.exp == m$nea.out) & (m$nea.exp == m$ea.out)
  m$beta.out[flip] <- -m$beta.out[flip]
  m$eaf.out[flip]  <- 1 - m$eaf.out[flip]
  # drop palindromes AT/GC that remain allele-frequency ambiguous
  pal <- ((m$ea.exp == "A" & m$nea.exp == "T") | (m$ea.exp == "T" & m$nea.exp == "A") |
          (m$ea.exp == "G" & m$nea.exp == "C") | (m$ea.exp == "C" & m$nea.exp == "G"))
  # keep palindrome only if EAF is not close to 0.5 (robustness)
  m <- m[!(pal & abs(m$eaf.exp - 0.5) > 0.38), ]

  # unique SNP id
  m[, SNP := id.exp]
  setnames(m,
           c("beta.exp", "se.exp", "eaf.exp", "pval.exp", "ea.exp", "nea.exp",
             "beta.out", "se.out", "eaf.out", "pval.out"),
           c("beta.exposure", "se.exposure", "eaf.exposure", "pval.exposure", "effect_allele.exposure", "other_allele.exposure",
             "beta.outcome", "se.outcome", "eaf.outcome", "pval.outcome"))
  m[, f_stat := beta.exposure^2 / se.exposure^2]
  m <- m[f_stat >= MIN_F]
  m
}

run_mr <- function(dat, exposure_name, outcome_name, n_exp, n_out, direction_label) {
  res <- list()

  if (nrow(dat) == 0) return(res)

  bx <- dat$beta.exposure
  by <- dat$beta.outcome
  bxse <- dat$se.exposure
  byse <- dat$se.outcome
  rs <- dat$SNP

  # IVW
  tryCatch({
    ivw <- mr_ivw(mr_input(bx, bxse, by, byse, snps = rs))
    res$ivw <- list(method = "IVW", b = ivw@Estimate, se = ivw@StdError,
                    p = ivw@Pvalue, lo = ivw@CILower, hi = ivw@CIUpper)
  }, error = function(e) res$ivw <- list(method = "IVW", b = NA, se = NA, p = NA, lo = NA, hi = NA, error = conditionMessage(e)))

  # Weighted median
  tryCatch({
    wm <- mr_median(mr_input(bx, bxse, by, byse, snps = rs))
    res$weighted_median <- list(method = "Weighted median", b = wm@Estimate, se = wm@StdError,
                                 p = wm@Pvalue, lo = wm@CILower, hi = wm@CIUpper)
  }, error = function(e) res$weighted_median <- list(method = "Weighted median", b = NA, se = NA, p = NA, lo = NA, hi = NA, error = conditionMessage(e)))

  # MR-Egger
  tryCatch({
    eg <- mr_egger(mr_input(bx, bxse, by, byse, snps = rs))
    res$egger <- list(method = "MR-Egger", b = eg@Estimate, se = eg@StdError.Est,
                      p = eg@Pvalue.Est, intercept = eg@Intercept,
                      intercept_p = eg@Pvalue.Int)
  }, error = function(e) res$egger <- list(method = "MR-Egger", b = NA, se = NA, p = NA, intercept = NA, intercept_p = NA, error = conditionMessage(e)))

  # MR-RAPS
  tryCatch({
    raps <- mr.raps(bx, by, bxse, byse, diagnostics = FALSE)
    res$raps <- list(method = "MR-RAPS", b = raps$beta.hat, se = raps$beta.se,
                     p = 2 * pnorm(-abs(raps$beta.hat / raps$beta.se)))
  }, error = function(e) res$raps <- list(method = "MR-RAPS", b = NA, se = NA, p = NA, error = conditionMessage(e)))

  # Radial IVW
  tryCatch({
    rad <- format_radial(bx, by, bxse, byse, rs)
    rivw <- ivw_radial(rad, weights = 3, alpha = 0.05)
    res$radial_ivw <- list(method = "Radial IVW", b = rivw$coef["Estimate", "B"],
                           se = rivw$coef["Estimate", "Std Error"],
                           p = rivw$coef["Estimate", "Pr(>|t|)"],
                           q = rivw$qstatistic, q_p = rivw$q_pval)
  }, error = function(e) res$radial_ivw <- list(method = "Radial IVW", b = NA, se = NA, p = NA, q = NA, q_p = NA, error = conditionMessage(e)))

  # MR-PRESSO
  tryCatch({
    presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
                        SdOutcome = "se.outcome", SdExposure = "se.exposure",
                        data = dat, OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
                        NbDistribution = 1000, SignifThreshold = 0.05)
    raw <- presso$`Main MR results`[1, ]
    out <- presso$`Main MR results`[2, ]
    res$presso_raw <- list(method = "MR-PRESSO raw", b = raw$`Causal Estimate`,
                            se = raw$Sd, p = raw$`P-value`)
    if (!is.null(out) && nrow(presso$`Main MR results`) > 1) {
      res$presso_outlier_corrected <- list(method = "MR-PRESSO outlier-corrected",
                                           b = out$`Causal Estimate`,
                                           se = out$Sd, p = out$`P-value`)
    }
    res$presso_distortion <- presso$`Distortion Test`$`Distortion Coefficient`
    res$presso_outliers <- presso$`MR-PRESSO results`$`Outlier Test`$SNP[presso$`MR-PRESSO results`$`Outlier Test`$Pvalue < 0.05]
  }, error = function(e) res$presso_error <- conditionMessage(e))

  # Steiger directionality
  tryCatch({
    steig <- steiger_filtering(dat, nsamp_exp = n_exp, nsamp_out = n_out)
    correct <- sum(steig$steiger_dir & steig$steiger_pval < 0.05, na.rm = TRUE)
    res$steiger <- list(correct_direction_snps = as.integer(correct), total_snps = nrow(steig))
  }, error = function(e) res$steiger <- list(correct_direction_snps = NA, total_snps = nrow(dat), error = conditionMessage(e)))

  res$n_snps <- nrow(dat)
  res$direction <- direction_label
  res$exposure <- exposure_name
  res$outcome <- outcome_name
  res
}

make_plots <- function(dat, prefix, label) {
  # scatter
  p1 <- ggplot(dat, aes(x = beta.exposure, y = beta.outcome)) +
    geom_errorbarh(aes(xmin = beta.exposure - 1.96 * se.exposure,
                       xmax = beta.exposure + 1.96 * se.exposure), alpha = 0.3) +
    geom_errorbar(aes(ymin = beta.outcome - 1.96 * se.outcome,
                      ymax = beta.outcome + 1.96 * se.outcome), alpha = 0.3) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, colour = "red") +
    labs(title = paste(label, "scatter"),
         x = paste0(dat$exposure[1], " beta"),
         y = paste0(dat$outcome[1], " beta")) +
    theme_bw()
  ggsave(file.path(plot_dir, paste0(prefix, "_scatter.pdf")), p1, width = 6, height = 5)

  # funnel
  dat$ivw_b <- sum(dat$beta.exposure * dat$beta.outcome / dat$se.outcome^2) / sum(dat$beta.exposure^2 / dat$se.outcome^2)
  dat$prec <- 1 / dat$se.outcome
  dat$effect <- dat$beta.outcome / dat$beta.exposure
  p2 <- ggplot(dat, aes(x = effect, y = prec)) +
    geom_vline(xintercept = dat$ivw_b[1], linetype = "dashed", colour = "red") +
    geom_point() +
    labs(title = paste(label, "funnel"), x = "SNP effect (outcome/exposure)", y = "Precision") +
    theme_bw()
  ggsave(file.path(plot_dir, paste0(prefix, "_funnel.pdf")), p2, width = 6, height = 5)
}

# -------------------------------------------------------------------------
# Main
# -------------------------------------------------------------------------
message("Reading COPD GWAS ...")
copd <- read_gwas(copd_file, n = Inf)
message("Reading HF GWAS ...")
hf   <- read_gwas(hf_file, n = Inf)

# Forward direction: COPD -> HF
message("Forward direction: COPD -> HF")
copd_ivs <- copd[pval < P_THRESHOLD]
copd_ivs <- distance_clump(copd_ivs)
forward <- harmonise(copd_ivs, hf)
if (!is.null(forward) && nrow(forward) > 0) {
  forward[, exposure := "COPD"]
  forward[, outcome := "Heart failure"]
  make_plots(forward, "forward_copd_hf", "COPD -> HF")
  fwrite(forward, file.path(out_dir, "forward_instruments.csv"))
  res_forward <- run_mr(forward, "COPD", "Heart failure", N_COPD, N_HF, "COPD -> HF")
} else {
  res_forward <- list(error = "No forward instruments after harmonisation")
}

# Reverse direction: HF -> COPD
message("Reverse direction: HF -> COPD")
hf_ivs <- hf[pval < P_THRESHOLD]
hf_ivs <- distance_clump(hf_ivs)
reverse <- harmonise(hf_ivs, copd)
if (!is.null(reverse) && nrow(reverse) > 0) {
  reverse[, exposure := "Heart failure"]
  reverse[, outcome := "COPD"]
  make_plots(reverse, "reverse_hf_copd", "HF -> COPD")
  fwrite(reverse, file.path(out_dir, "reverse_instruments.csv"))
  res_reverse <- run_mr(reverse, "Heart failure", "COPD", N_HF, N_COPD, "HF -> COPD")
} else {
  res_reverse <- list(error = "No reverse instruments after harmonisation")
}

# Save results
out <- list(forward = res_forward, reverse = res_reverse,
            params = list(p_threshold = P_THRESHOLD, clump_kb = CLUMP_KB, min_f = MIN_F,
                        N_COPD = N_COPD, N_HF = N_HF))
write_json(out, file.path(out_dir, "mr_robustness_results.json"), pretty = TRUE, auto_unbox = TRUE)

# Also produce a flat CSV summary
summary_rows <- lapply(list(res_forward, res_reverse), function(r) {
  if (is.null(r) || !is.null(r$error)) return(NULL)
  rows <- data.table()
  for (nm in names(r)) {
    x <- r[[nm]]
    if (is.list(x) && !is.null(x$method)) {
      rows <- rbind(rows, as.data.table(x), fill = TRUE)
    }
  }
  if (nrow(rows) > 0) rows[, direction := r$direction]
  rows
})
summary_dt <- rbindlist(summary_rows, fill = TRUE)
if (nrow(summary_dt) > 0) {
  fwrite(summary_dt, file.path(out_dir, "mr_robustness_summary.csv"))
}

message("Done. Results in ", out_dir)
