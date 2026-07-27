# ============================================================================
# Module 04 - Replication: three-cohort meta-analysis (COPD -> HF)
# Fixed-effect (inverse-variance) and DerSimonian-Laird random-effects
# meta-analysis of the three two-sample MR estimates:
#   1. EUR discovery  (COPD ebi-a-GCST90018807 -> HF ebi-a-GCST009541, HERMES)
#   2. FinnGen R12 I9_HEARTFAIL (strict HF replication)
#   3. BBJ bbj-a-109 (East Asian replication)
# Inputs : per-cohort log-OR and SE (below; identical to the values archived
#          in results/legacy/ and reproduced by the module-01/04 scripts)
# Outputs: results/t1_three_cohort_meta.json, results/t1_three_cohort_meta.csv
# Verified against the manuscript: FE OR 1.059 (1.022-1.098); Q = 14.5,
# df = 2, I2 = 86%; RE OR 1.042 (0.94-1.15) -- see VALIDATION.md.
# ============================================================================

REPRO_ROOT <- Sys.getenv("REPRO_ROOT", unset = ".")

cohorts <- data.frame(
  cohort = c("EUR discovery (GCST009541)", "FinnGen R12 I9_HEARTFAIL", "BBJ bbj-a-109"),
  beta   = c(0.1400, 0.0363, -0.0694),
  se     = c(0.0310, 0.0256,  0.0489),
  stringsAsFactors = FALSE
)

or_ci <- function(b, se) c(OR = exp(b), lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se))

# ---- per-cohort ----
per <- lapply(seq_len(nrow(cohorts)), function(i) {
  ci <- or_ci(cohorts$beta[i], cohorts$se[i])
  list(beta = cohorts$beta[i], se = cohorts$se[i],
       OR = unname(ci["OR"]), CI = c(unname(ci["lo"]), unname(ci["hi"])))
})
names(per) <- cohorts$cohort

# ---- fixed effect ----
w  <- 1 / cohorts$se^2
b_fe <- sum(w * cohorts$beta) / sum(w)
se_fe <- sqrt(1 / sum(w))
z_fe <- b_fe / se_fe
p_fe <- 2 * pnorm(-abs(z_fe))
ci_fe <- or_ci(b_fe, se_fe)

# ---- heterogeneity ----
Q  <- sum(w * (cohorts$beta - b_fe)^2)
df <- nrow(cohorts) - 1
p_Q <- pchisq(Q, df = df, lower.tail = FALSE)
I2 <- max(0, (Q - df) / Q) * 100
C <- sum(w) - sum(w^2) / sum(w)
tau2 <- max(0, (Q - df) / C)

# ---- DerSimonian-Laird random effect ----
wr <- 1 / (cohorts$se^2 + tau2)
b_re <- sum(wr * cohorts$beta) / sum(wr)
se_re <- sqrt(1 / sum(wr))
z_re <- b_re / se_re
p_re <- 2 * pnorm(-abs(z_re))
ci_re <- or_ci(b_re, se_re)

res <- list(
  per_cohort = per,
  fixed_effect = list(beta = b_fe, se = se_fe, OR = unname(ci_fe["OR"]),
                      CI = c(unname(ci_fe["lo"]), unname(ci_fe["hi"])), p = p_fe),
  heterogeneity = list(Q = Q, df = df, p = p_Q, I2_pct = I2, tau2 = tau2),
  random_effect = list(beta = b_re, se = se_re, OR = unname(ci_re["OR"]),
                       CI = c(unname(ci_re["lo"]), unname(ci_re["hi"])), p = p_re)
)

dir.create(file.path(REPRO_ROOT, "results"), showWarnings = FALSE)
jsonlite::write_json(res, file.path(REPRO_ROOT, "results", "t1_three_cohort_meta.json"),
                     auto_unbox = TRUE, digits = 6, pretty = TRUE)

tab <- data.frame(
  model = c("fixed_effect", "random_effect"),
  beta  = c(b_fe, b_re), se = c(se_fe, se_re),
  OR    = c(ci_fe["OR"], ci_re["OR"]),
  CI_lo = c(ci_fe["lo"], ci_re["lo"]), CI_hi = c(ci_fe["hi"], ci_re["hi"]),
  p     = c(p_fe, p_re),
  Q = Q, df = df, I2_pct = I2, tau2 = tau2
)
write.csv(tab, file.path(REPRO_ROOT, "results", "t1_three_cohort_meta.csv"), row.names = FALSE)

cat(sprintf("FE : beta %.4f se %.4f  OR %.3f (%.3f-%.3f)  p %.3g\n",
            b_fe, se_fe, ci_fe["OR"], ci_fe["lo"], ci_fe["hi"], p_fe))
cat(sprintf("Het: Q %.1f df %d p %.3g  I2 %.1f%%  tau2 %.4f\n", Q, df, p_Q, I2, tau2))
cat(sprintf("RE : beta %.4f se %.4f  OR %.3f (%.3f-%.3f)  p %.3g\n",
            b_re, se_re, ci_re["OR"], ci_re["lo"], ci_re["hi"], p_re))
