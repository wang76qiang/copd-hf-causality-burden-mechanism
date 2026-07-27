# ============================================================================
# Module 08 - Drug target: eQTLGen whole-blood MR + coloc addendum
# (requires OPENGWAS_JWT for any uncached calls)
# ============================================================================

# --- repository paths -----------------------------------------------------------
# Run this script from the repository root, or set REPRO_ROOT to the repo root.
REPRO_ROOT <- Sys.getenv("REPRO_ROOT", unset = ".")
# GBD_RAW_DIR: folder with the raw IHME GBD 2021 csv downloads (see data/README.md)
GBD_RAW_DIR <- Sys.getenv("GBD_RAW_DIR", unset = file.path(REPRO_ROOT, "data", "raw", "gbd_2021"))
# V8_DIR: legacy V8 project tree (NOT shipped; only needed by v8_* legacy scripts)
V8_DIR <- Sys.getenv("V8_DIR", unset = "v8_legacy")
# .libPaths(...)  # removed user-local library path; see r_packages.txt
# NOTE: an OpenGWAS JWT is required for OpenGWAS/ieugwasr calls.
# Request your own token at https://api.opengwas.io and either
#   export OPENGWAS_JWT=<token>          (shell, picked up automatically), or
#   Sys.setenv(OPENGWAS_JWT = "<token>")  (in R, before running).
if (Sys.getenv("OPENGWAS_JWT") == "")
  message("WARNING: OPENGWAS_JWT is not set; OpenGWAS API calls will fail.")
# A1b addendum: eQTLGen whole blood (eqtl-a-ENSG00000106366, N=31,684) MR + coloc
library(ieugwasr); library(jsonlite)
dl <- file.path(REPRO_ROOT, "data_downloads")
OUT <- file.path(REPRO_ROOT, "results")
LOG <- file(file.path(REPRO_ROOT, "logs/t1_a1b_serpine1_drugtarget.log"), open = "at"); sink(LOG, split = TRUE)
cat("\n== addendum: eQTLGen whole blood ==\n")

d <- as.data.frame(jsonlite::fromJSON(paste0(dl, "/tmp/eqtlgen_serpine1.json")))
d <- d[!is.na(d$beta) & d$se > 0 & d$rsid != "", ]
d <- d[order(d$p), ]; d <- d[!duplicated(d$rsid), ]
cat("eQTLGen SERPINE1 cis variants:", nrow(d), " min p:", signif(min(d$p), 3), "\n")

hf <- read.delim(paste0(dl, "/hf_serpine1_region.tsv"), stringsAsFactors = FALSE)
for (cc in c("beta","se","p","eaf")) hf[[cc]] <- suppressWarnings(as.numeric(hf[[cc]]))
hf <- hf[which(!duplicated(hf$rsid) & !is.na(hf$se) & hf$se > 0 & !is.na(hf$beta)), ]

comp <- function(a) chartr("ACGT", "TGCA", toupper(a))
align <- function(r1, r2, e1, e2) {
  r1 <- toupper(r1); r2 <- toupper(r2); e1 <- toupper(e1); e2 <- toupper(e2)
  if (e1 == r1 && e2 == r2) return(1); if (e1 == r2 && e2 == r1) return(-1)
  if (comp(e1) == r1 && comp(e2) == r2) return(1); if (comp(e1) == r2 && comp(e2) == r1) return(-1)
  NA_real_
}

# ---- MR ----
iv <- d[d$p < 5e-6, ]
cat("instruments p<5e-6:", nrow(iv), "\n")
cl <- ieugwasr::ld_clump(data.frame(rsid = iv$rsid, pval = iv$p), clump_r2 = 0.1,
                         clump_kb = 1000, pop = "EUR")
iv <- iv[iv$rsid %in% cl$rsid, ]
cat("after clump r2<0.1:", nrow(iv), "\n")
m <- merge(iv, hf, by = "rsid", suffixes = c(".eq",".hf"))
sgn <- mapply(align, m$ea.eq, m$nea.eq, m$ea.hf, m$nea.hf)
kix <- which(!is.na(sgn)); m <- m[kix, ]; sgn <- sgn[kix]
m$beta_hf <- m$beta.hf * sgn
m$ratio <- m$beta_hf / m$beta.eq; m$ratio_se <- m$se.hf / abs(m$beta.eq)
w <- 1/m$ratio_se^2; b <- sum(w*m$ratio)/sum(w); se <- sqrt(1/sum(w))
Q <- sum(w*(m$ratio-b)^2); Qp <- pchisq(Q, nrow(m)-1, lower.tail = FALSE)
mr <- data.frame(tissue = "blood_eQTLGen", dataset = "eqtl-a-ENSG00000106366",
                 method = "IVW", nsnp = nrow(m), beta = b, se = se,
                 or = exp(b), or_lci95 = exp(b-1.96*se), or_uci95 = exp(b+1.96*se),
                 pval = 2*pnorm(-abs(b/se)), Q = Q, Q_pval = Qp)
write.csv(m, file.path(OUT, "t1_a1b_harmonised_blood_eQTLGen.csv"), row.names = FALSE)
cat("=== eQTLGen MR ===\n"); print(mr)

# append to MR results file
mrfile <- file.path(OUT, "t1_drugtarget_mr_serpine1_7tissues.csv")
old <- if (file.exists(mrfile) && file.size(mrfile) > 0) {
  tryCatch(read.csv(mrfile, stringsAsFactors = FALSE), error = function(e) NULL)
} else NULL
old <- if (is.null(old) || nrow(old) == 0) mr else plyr::rbind.fill(old, mr)
write.csv(old, mrfile, row.names = FALSE)

# ---- coloc ----
m2 <- merge(d[, c("rsid","ea","nea","beta","se","eaf","p")],
            hf[, c("rsid","ea","nea","beta","se","eaf")],
            by = "rsid", suffixes = c(".eq",".hf"))
sgn2 <- mapply(align, m2$ea.eq, m2$nea.eq, m2$ea.hf, m2$nea.hf)
k2 <- which(!is.na(sgn2)); m2 <- m2[k2, ]; sgn2 <- sgn2[k2]
m2$beta_hf <- m2$beta.hf * sgn2
m2$maf <- suppressWarnings(as.numeric(m2$eaf.eq))
m2 <- m2[which(!is.na(m2$maf) & m2$maf > 0 & m2$se.eq > 0 & m2$se.hf > 0), ]
W <- 0.2^2
labf <- function(b, v) 0.5*log(v/(v+W)) + (b^2/(2*v))*(W/(v+W))
l1 <- labf(m2$beta.eq, m2$se.eq^2); l2 <- labf(m2$beta_hf, m2$se.hf^2)
s1 <- sum(exp(l1-max(l1)))*exp(max(l1)); s2 <- sum(exp(l2-max(l2)))*exp(max(l2))
ls12 <- l1+l2; s12 <- sum(exp(ls12-max(ls12)))*exp(max(ls12))
wts <- c(H0=1, H1=1e-4*s1, H2=1e-4*s2, H3=1e-8*s1*s2, H4=1e-5*s12)
pp <- wts/sum(wts)
co <- data.frame(tissue = "blood_eQTLGen", dataset = "eqtl-a-ENSG00000106366",
                 nsnps = nrow(m2), PP.H0 = pp["H0"], PP.H1 = pp["H1"], PP.H2 = pp["H2"],
                 PP.H3 = pp["H3"], PP.H4 = pp["H4"])
cat("=== eQTLGen coloc ===\n"); print(co)
co_old <- read.csv(file.path(OUT, "t1_coloc_serpine1.csv"), stringsAsFactors = FALSE)
co_all <- plyr::rbind.fill(co_old, co)
write.csv(co_all, file.path(OUT, "t1_coloc_serpine1.csv"), row.names = FALSE)
sink(); close(LOG)
cat("A1B ADDENDUM DONE\n")
