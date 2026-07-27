# ============================================================================
# Module 03 - Mediation MR (shared helpers)
# Allele-alignment (align_beta) and IVW helper functions + shared paths.
# Sourced by t1_med_stepA.R and t1_med_stepB.R.
# ============================================================================

# --- repository paths -----------------------------------------------------------
# Run this script from the repository root, or set REPRO_ROOT to the repo root.
REPRO_ROOT <- Sys.getenv("REPRO_ROOT", unset = ".")
# GBD_RAW_DIR: folder with the raw IHME GBD 2021 csv downloads (see data/README.md)
GBD_RAW_DIR <- Sys.getenv("GBD_RAW_DIR", unset = file.path(REPRO_ROOT, "data", "raw", "gbd_2021"))
# V8_DIR: legacy V8 project tree (NOT shipped; only needed by v8_* legacy scripts)
V8_DIR <- Sys.getenv("V8_DIR", unset = "v8_legacy")
# .libPaths(...)  # removed user-local library path; see r_packages.txt
# shared helpers + objects for the mediation pipeline
# NOTE: an OpenGWAS JWT is required for OpenGWAS/ieugwasr calls.
# Request your own token at https://api.opengwas.io and either
#   export OPENGWAS_JWT=<token>          (shell, picked up automatically), or
#   Sys.setenv(OPENGWAS_JWT = "<token>")  (in R, before running).
if (Sys.getenv("OPENGWAS_JWT") == "")
  message("WARNING: OPENGWAS_JWT is not set; OpenGWAS API calls will fail.")
library(ieugwasr)
DL <- file.path(REPRO_ROOT, "data_downloads")
OUT <- file.path(REPRO_ROOT, "results")
comp <- function(a) chartr("ACGT", "TGCA", a)
align_beta <- function(ea_ref, nea_ref, ea, nea, beta) {
  ea <- toupper(ea); nea <- toupper(nea); ea_ref <- toupper(ea_ref); nea_ref <- toupper(nea_ref)
  ok0 <- !is.na(ea) & !is.na(nea) & !is.na(ea_ref) & !is.na(nea_ref)
  eq <- function(x, y) ok0 & !is.na(x) & !is.na(y) & x == y
  out <- rep(NA_real_, length(beta))
  same <- eq(ea, ea_ref) & eq(nea, nea_ref); out[same] <- beta[same]
  flip <- !same & eq(ea, nea_ref) & eq(nea, ea_ref); out[flip] <- -beta[flip]
  cs <- !same & !flip & eq(comp(ea), ea_ref) & eq(comp(nea), nea_ref); out[cs] <- beta[cs]
  cf <- !same & !flip & !cs & eq(comp(ea), nea_ref) & eq(comp(nea), ea_ref); out[cf] <- -beta[cf]
  out
}
ivw <- function(bx, by, sey) { b <- sum(bx*by/sey^2)/sum(bx^2/sey^2)
  se <- sqrt(1/sum(bx^2/sey^2)); c(b = b, se = se) }
