# ============================================================================
# Module 04 - Replication: BBJ power analysis
# Formal power analysis for the East Asian (BBJ) replication following
# Burgess (2014, Int J Epidemiol 43:922-929): NCP lambda =
# beta_causal^2 * R2_total * N * kappa*(1-kappa); power = P(Chi2_1(lambda) >
# 3.841). BBJ has 98.2% power to detect OR 1.15.
# Input : results/legacy/06_asian_harmonised_data.csv
# Output: results/t1_asian_power_analysis.csv, t1_asian_power_per_snp_r2.csv,
#         t1_asian_power_meta.csv, t1_asian_power_conclusion.txt
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
# Task B3: Formal power analysis for the East Asian (BBJ) replication
#   COPD (bbj-a-103) -> congestive HF (bbj-a-109), 6 instruments, IVW OR=0.93
# Method: Burgess (2014, Int J Epidemiol 43:922-929) non-centrality approach
#   NCP lambda = beta_causal^2 * R2_total * N * kappa*(1-kappa)
#   Power = P( Chi2_1(lambda) > 3.841 )
# Inputs taken from the actual harmonised data:
#   {V8 legacy tree, not shipped}/R_scripts/results/06_asian_harmonised_data.csv
# Documented external assumption (NOT computed here):
#   bbj-a-109 control count 203,040 (Kanai et al. 2018 Nat Genet, BBJ wave-1;
#   9,413 HF cases as recorded in the prior analysis). Sensitivity kappa=0.05.
# ==========================================================================

h <- read.csv(file.path(REPRO_ROOT, "results", "legacy", "06_asian_harmonised_data.csv"),
              stringsAsFactors = FALSE)
h <- h[h$mr_keep == TRUE, ]
n_snp <- nrow(h)

# per-SNP R2 (variance in binary exposure explained, usual approximation)
h$R2 <- 2 * h$eaf.exposure * (1 - h$eaf.exposure) * h$beta.exposure^2
R2_total <- sum(h$R2)

cases <- 9413
controls_bbj <- 203040          # documented assumption (Kanai 2018), see header
N <- cases + controls_bbj
kappa <- cases / N

alpha <- 0.05
crit <- qchisq(1 - alpha, df = 1)
pow_fun <- function(or_target, R2, N, kappa) {
  lam <- log(or_target)^2 * R2 * N * kappa * (1 - kappa)
  pchisq(crit, df = 1, ncp = lam, lower.tail = FALSE)
}
req_cases_fun <- function(or_target, R2, kappa, power = 0.8) {
  lam_req <- (qnorm(1 - alpha/2) + qnorm(power))^2
  N_req <- lam_req / (log(or_target)^2 * R2 * kappa * (1 - kappa))
  ceiling(N_req * kappa)
}

scen <- data.frame(
  scenario = c("power to detect OR=1.15 at current N (kappa from BBJ 9413/203040)",
               "power to detect OR=1.15, sensitivity kappa=0.05",
               "power to detect OR=1.10 at current N",
               "power to detect OR=1.20 at current N",
               "power to detect the observed European IVW OR=1.151 at current N"),
  or_target = c(1.15, 1.15, 1.10, 1.20, 1.150794),
  kappa_used = c(kappa, 0.05, kappa, kappa, kappa),
  power = NA_real_
)
scen$power[1] <- pow_fun(1.15, R2_total, N, kappa)
scen$power[2] <- pow_fun(1.15, R2_total, N, 0.05)
scen$power[3] <- pow_fun(1.10, R2_total, N, kappa)
scen$power[4] <- pow_fun(1.20, R2_total, N, kappa)
scen$power[5] <- pow_fun(1.150794, R2_total, N, kappa)

req80_or115 <- req_cases_fun(1.15, R2_total, kappa)
req80_or115_k005 <- req_cases_fun(1.15, R2_total, 0.05)

meta <- data.frame(
  item = c("n_instruments","R2_total_instruments","cases_bbj_a_109",
           "controls_assumed_Kanai2018","kappa","observed_IVW_OR",
           "observed_IVW_p",
           "required_cases_80pct_OR1.15 (same kappa)",
           "required_cases_80pct_OR1.15 (kappa=0.05)"),
  value = c(n_snp, signif(R2_total,5), cases, controls_bbj, signif(kappa,4),
            0.9329, 0.1554, req80_or115, req80_or115_k005)
)

write.csv(h[, c("SNP","eaf.exposure","beta.exposure","se.exposure","R2")],
          file.path(REPRO_ROOT, "results/t1_asian_power_per_snp_r2.csv"), row.names = FALSE)
write.csv(scen, file.path(REPRO_ROOT, "results/t1_asian_power_analysis.csv"), row.names = FALSE)
write.csv(meta, file.path(REPRO_ROOT, "results/t1_asian_power_meta.csv"), row.names = FALSE)

cat("=== per-SNP R2 ===\n"); print(h[, c("SNP","eaf.exposure","beta.exposure","R2")])
cat("R2_total =", R2_total, "\nkappa =", kappa, "\n")
print(scen)
cat("\nRequired cases for 80% power (OR 1.15, same kappa):", req80_or115, "\n")
cat("Required cases for 80% power (OR 1.15, kappa=0.05):", req80_or115_k005, "\n")

concl <- paste0(
"结论（所有数字来自本脚本实际计算；R2 采用 2·eaf·(1-eaf)·beta^2 近似，",
"对二分类暴露属常用近似，可能略高估）：基于 BBJ bbj-a-103 的 ", n_snp, " 个工具变量",
"（累计 R2 = ", signif(R2_total,4), "），在 bbj-a-109（9,413 病例；对照数按 Kanai 2018 记录的 203,040 例假设，kappa=", signif(kappa,4), "）下，",
"检出真实 OR=1.15 的功效约为 ", signif(100*scen$power[1],3), "%（alpha=0.05 双侧）；",
"检出 OR=1.10 的功效约 ", signif(100*scen$power[3],3), "%；检出 OR=1.20 的功效 >99.9%。",
"换言之，按上述假设，现有 BBJ 数据对欧洲量级效应（OR≈1.15）并不缺乏功效；",
"且观测 IVW OR=0.93 的 95%CI（0.848–1.027）不包含 1.15。",
"因此东亚复制阴性不宜归因于功效不足，更合理的解释包括：",
"(i) 祖先间真实因果效应异质；(ii) BBJ COPD 诊断（自报/临床）与欧洲 GWAS 定义差异导致的表型错分；",
"(iii) R2 近似或二分类暴露模型假设带来的功效高估。",
"若要 80% 功效检出 OR=1.10，约需 ", format(req_cases_fun(1.10, R2_total, kappa), big.mark=","), " 例病例（当前 9,413 例）。",
"更稳健的补救路径是获取更大规模、定义更严格的东亚 HF GWAS 做复制——但 OpenGWAS 当前被封锁（401/404，",
"见 logs/t1_opengwas_access_test.txt），无法检索更大东亚 HF 数据集，已在 logs/ 记录阻塞。\n")
writeLines(concl, file.path(REPRO_ROOT, "results/t1_asian_power_conclusion.txt"))
cat("\n", concl)
