# ============================================================================
# Module 05 - GBD: BAPC projection to 2050
# Fully Bayesian age-period-cohort projection (BAPC; Poisson likelihood, RW2
# age/period/cohort effects, overdispersion, R-INLA) of global COPD-
# attributable HF cases to 2050, plus a rate-space sensitivity with the
# population frozen at the 2021 age structure.
# Input : results/t1_bapc_agepanel.csv
# Output: results/t1_bapc_projections.csv, t1_bapc_total_cases.csv,
#         t1_bapc_total_cases_frozenpop.csv, t1_bapc_summary.json
# ============================================================================

# --- repository paths -----------------------------------------------------------
# Run this script from the repository root, or set REPRO_ROOT to the repo root.
REPRO_ROOT <- Sys.getenv("REPRO_ROOT", unset = ".")
# GBD_RAW_DIR: folder with the raw IHME GBD 2021 csv downloads (see data/README.md)
GBD_RAW_DIR <- Sys.getenv("GBD_RAW_DIR", unset = file.path(REPRO_ROOT, "data", "raw", "gbd_2021"))
# V8_DIR: legacy V8 project tree (NOT shipped; only needed by v8_* legacy scripts)
V8_DIR <- Sys.getenv("V8_DIR", unset = "v8_legacy")
# T1: BAPC Bayesian APC projection of global COPD-attributed HF prevalence cases to 2050
# R 4.6.1 + INLA 24.12.11 (V8 R_libs_44, matching Windows binary) + BAPC 0.0.37
# Primary     : count-space BAPC (pyrs=1) -> projects case counts directly; implicitly
#               propagates historical demographic + epidemiologic trends (comparable to
#               the existing frequentist APC on counts).
# Sensitivity : rate-space BAPC with population frozen at the 2021 age structure
#               (isolates the pure epidemiologic APC signal).
# .libPaths(...)  # removed user-local library path; see r_packages.txt
suppressPackageStartupMessages(library(BAPC))

OUT <- REPRO_ROOT
panel <- read.csv(file.path(OUT, "results/t1_bapc_agepanel.csv"), check.names = FALSE)
years_obs <- panel$year                       # 1990-2021
case_cols <- grep("^cases_", names(panel), value = TRUE)
ages_all <- sub("^cases_", "", case_cols)
keep <- !ages_all %in% c("<5 years", "5-9 years", "10-14 years")  # structurally zero in GBD
ages <- ages_all[keep]; case_cols <- case_cols[keep]
pop_cols <- sub("^cases_", "pop_", case_cols)
YRS <- 1990:2050
NFUT <- 2050 - max(years_obs)                 # 29

mk_epi <- function() {
  e <- as.data.frame(round(as.matrix(panel[, case_cols])))  # GBD non-integer estimates -> round for Poisson
  names(e) <- ages
  fut <- as.data.frame(matrix(NA_integer_, NFUT, ncol(e))); names(fut) <- ages
  e <- rbind(e, fut); rownames(e) <- YRS; e
}
epi <- mk_epi()

run_bapc <- function(pyrs_df, tag) {
  apc <- APCList(epi, pyrs_df, gf = 5)
  res <- BAPC(apc, predict = list(npredict = NFUT, retro = FALSE), verbose = FALSE)
  proj <- agespec.proj(res)
  out <- do.call(rbind, lapply(ages, function(a) {
    m <- proj[[a]]
    data.frame(year = as.integer(rownames(m)), age = a,
               mean = m[, "mean"], sd = m[, "sd"], row.names = NULL)
  }))
  cat(tag, "done; periods", min(out$year), "-", max(out$year), "\n")
  out
}

# ---------- primary: count space ----------
one <- as.data.frame(matrix(1, nrow(epi), ncol(epi), dimnames = dimnames(epi)))
p <- run_bapc(one, "count-space")
p$lower <- pmax(0, p$mean - 1.96 * p$sd)
p$upper <- p$mean + 1.96 * p$sd
proj_df <- p[, c("year", "age", "mean", "lower", "upper")]
proj_df$model <- "BAPC_count_space"
write.csv(proj_df, file.path(OUT, "results/t1_bapc_projections.csv"), row.names = FALSE)

tot <- do.call(rbind, lapply(split(p, p$year), function(g) {
  m <- sum(g$mean); s <- sqrt(sum(g$sd^2))     # independent-age approximation
  data.frame(year = g$year[1], mean = m, lower = max(0, m - 1.96 * s), upper = m + 1.96 * s)
}))
write.csv(tot, file.path(OUT, "results/t1_bapc_total_cases.csv"), row.names = FALSE)

# ---------- sensitivity: rate space, population frozen at 2021 ----------
sens_tot <- tryCatch({
  pyrs <- as.data.frame(as.matrix(panel[, pop_cols])); names(pyrs) <- ages
  futp <- pyrs[rep(nrow(pyrs), NFUT), , drop = FALSE]         # freeze at 2021
  pyrs2 <- rbind(pyrs, futp); rownames(pyrs2) <- YRS
  s <- run_bapc(pyrs2, "rate-space")
  # NOTE: agespec.proj() in BAPC 0.0.37 returns expected COUNTS (rate x pyrs), not rates
  st <- do.call(rbind, lapply(split(s, s$year), function(g) {
    cm <- sum(g$mean); cv <- sum(g$sd^2)
    data.frame(year = g$year[1], cases_mean = cm,
               cases_lower = max(0, cm - 1.96 * sqrt(cv)), cases_upper = cm + 1.96 * sqrt(cv))
  }))
  write.csv(st, file.path(OUT, "results/t1_bapc_total_cases_frozenpop.csv"), row.names = FALSE)
  st
}, error = function(e) { cat("rate-space ERR:", conditionMessage(e), "\n"); NULL })

grab <- function(df, yr, col) df[df$year == yr, col][1]
base2021 <- grab(tot, 2021, "mean")
res_lines <- c("=== PRIMARY (count-space BAPC) ===")
for (yr in c(2030, 2035, 2040, 2050)) {
  res_lines <- c(res_lines, sprintf("%d: mean=%.0f  CrI=[%.0f, %.0f]  change vs 2021 fitted (%0.f) = %+.1f%%",
    yr, grab(tot, yr, "mean"), grab(tot, yr, "lower"), grab(tot, yr, "upper"),
    base2021, 100 * (grab(tot, yr, "mean") / base2021 - 1)))
}
if (!is.null(sens_tot)) {
  res_lines <- c(res_lines, "=== SENSITIVITY (rate-space, 2021 frozen population) ===")
  b21s <- grab(sens_tot, 2021, "cases_mean")
  for (yr in c(2035, 2050)) {
    res_lines <- c(res_lines, sprintf("%d: mean=%.0f  CrI=[%.0f, %.0f]  change vs 2021 = %+.1f%%",
      yr, grab(sens_tot, yr, "cases_mean"), grab(sens_tot, yr, "cases_lower"),
      grab(sens_tot, yr, "cases_upper"), 100 * (grab(sens_tot, yr, "cases_mean") / b21s - 1)))
  }
}
cat(paste(res_lines, collapse = "\n"), "\n")
writeLines(res_lines, file.path(OUT, "logs/t1_bapc_run.log"))
cat("DONE\n")
