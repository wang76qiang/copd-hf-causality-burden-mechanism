# -*- coding: utf-8 -*-
# ============================================================================
# Module 01 - MR discovery / governance
# Rebuild supplementary tables S1-S3 and S5-S7 from archived results and the
# shipped key-gene expression workbooks (data/derived/bulk/).
# Inputs : results/legacy/real_data_results.json, data/derived/bulk/*.xlsx
# Outputs: tables/Table_S1_instruments.csv (shipped), Table_S2_mr_results.csv,
#          Table_S3_mvmr_mediation.csv, Table_S5_shared_41genes.csv,
#          Table_S6_drugtarget.csv, Table_S7_gene_expression_validation.csv
# ============================================================================

import os as _os
# --- repository paths -----------------------------------------------------------
# REPO resolves to the repository root (code/<module>/this_file.py -> ../..).
# Override with the REPRO_ROOT environment variable if you run from elsewhere.
REPO = _os.environ.get("REPRO_ROOT", _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "..")))
# GBD_RAW_DIR: folder with the raw IHME GBD 2021 csv downloads (see data/README.md)
GBD_RAW_DIR = _os.environ.get("GBD_RAW_DIR", _os.path.join(REPO, "data", "raw", "gbd_2021"))
"""T1 A5b: file governance (renamed copies) + supplementary tables S1,S2,S3,S5,S6,S7."""
import pandas as pd, numpy as np, shutil, json, os

BASE = REPO
OUT = REPO
# rebuilt tables go to tables/rebuilt/ so the curated tables/ are never clobbered
REBUILT = OUT + "/tables/rebuilt"
os.makedirs(REBUILT, exist_ok=True)

# ---------- 1) file governance note ----------
# The original working directory contained two MISLABELED workbooks:
#   "S100A12.xlsx"      actually holds GSE57338 ACE expression (control/HF)
#   "S100A8两组比较.xlsx" actually holds GSE57148 THBS1 expression (Title/Value)
# Both are shipped under data/derived/bulk/ with corrected ASCII file names
# (see data/README.md); no copying is needed here.

# ---------- 2) Table_S1: copy from V10 ----------
s1src = BASE + "/tables/Table_S1_instruments.csv"  # shipped in the repo
assert os.path.exists(s1src), "Table_S1_instruments.csv missing from tables/"
print("S1 copied:", pd.read_csv(OUT + "/tables/Table_S1_instruments.csv").shape)

# ---------- 3) Table_S2: rebuild full MR results from real_data_results.json ----------
J = json.load(open(BASE + "/results/legacy/real_data_results.json"))
mr, rows = J["mr"], []
def or_ci(b, se):
    return np.exp(b), np.exp(b - 1.96 * se), np.exp(b + 1.96 * se)
def add(method, beta, se, pval, note=""):
    o, lo, hi = or_ci(beta, se)
    rows.append({"method": method, "beta": beta, "se": se, "pval": pval,
                 "OR": o, "OR_95CI_low": lo, "OR_95CI_high": hi, "note": note})
add("IVW (fixed effects)", mr["ivw"]["beta"], mr["ivw"]["se"], mr["ivw"]["pval"])
add("IVW (random effects, DL)", mr["ivw_random_effects"]["beta"], mr["ivw_random_effects"]["se"],
    mr["ivw_random_effects"]["pval"], f"tau2={mr['ivw_random_effects']['tau2']:.4f}")
add("Weighted median", mr["weighted_median"]["beta"], mr["weighted_median"]["se"], mr["weighted_median"]["pval"])
add("Simple mode", mr["simple_mode"]["beta"], mr["simple_mode"]["se"], mr["simple_mode"]["pval"])
add("MR-Egger (slope)", mr["mr_egger"]["slope"], mr["mr_egger"]["slope_se"], mr["mr_egger"]["slope_pval"],
    f"intercept={mr['mr_egger']['intercept']:.5f} (p={mr['mr_egger']['intercept_pval']:.3f})")
# radial MR: not present in json -> derive? mark as unavailable
rows.append({"method": "Radial MR (modified second-order weights)", "beta": np.nan, "se": np.nan,
             "pval": np.nan, "OR": np.nan, "OR_95CI_low": np.nan, "OR_95CI_high": np.nan,
             "note": "未在 real_data_results.json 中存档，待补算"})
# CAUSE
c = J["cause"]
rows.append({"method": "CAUSE (causal model gamma, sensitivity)", "beta": c["causal_gamma_log_or"],
             "se": np.nan, "pval": np.nan, "OR": c["causal_or"], "OR_95CI_low": c["causal_or_95ci_low"],
             "OR_95CI_high": c["causal_or_95ci_high"], "note": c["caveat"]})
# MR-Clust
mc = J["mr_clust"]
for cl in mc["clusters"]:
    rows.append({"method": f"MR-Clust cluster {cl['cluster']} ({cl['cluster_class']}, n={cl['n']} SNPs)",
                 "beta": cl["mean_wald_ratio"], "se": np.nan, "pval": np.nan,
                 "OR": np.exp(cl["mean_wald_ratio"]), "OR_95CI_low": np.nan, "OR_95CI_high": np.nan,
                 "note": "SNPs: " + cl["snps"]})
# heterogeneity & F
het, fs = mr["heterogeneity"], mr["f_statistics"]
rows.append({"method": "Cochran's Q (heterogeneity)", "beta": np.nan, "se": np.nan, "pval": het["pval"],
             "OR": np.nan, "OR_95CI_low": np.nan, "OR_95CI_high": np.nan,
             "note": f"Q={het['Q']:.2f}, df={het['df']}, I2={het['I2_percent']:.1f}%"})
rows.append({"method": "Instrument strength (F)", "beta": np.nan, "se": np.nan, "pval": np.nan,
             "OR": np.nan, "OR_95CI_low": np.nan, "OR_95CI_high": np.nan,
             "note": f"median F={fs['median']:.1f}, range {fs['min']:.1f}-{fs['max']:.1f}"})
pd.DataFrame(rows).to_csv(REBUILT + "/Table_S2_mr_results.csv", index=False)
print("S2 rows:", len(rows))

# ---------- 4) Table_S3 MVMR/mediation placeholder ----------
pd.DataFrame([{
    "analysis": "MVMR / mediation (COPD -> mediators -> HF)",
    "status": "占位：待补 — real_data_results.json 中无 MVMR/中介分析结果",
    "planned": "MVMR 调整吸烟/BMI/血压等混杂；两步 MR 中介（炎症/纤维化标志物）",
    "date_checked": "2026-07-26"}]).to_csv(REBUILT + "/Table_S3_mvmr_mediation.csv", index=False)

# ---------- 5) Table_S5 shared 41 genes ----------
genes = pd.read_csv(BASE + "/data/derived/bulk/gene_list_41_inflammation_fibrosis.csv", header=None)[0].tolist()
genes = [str(g).strip() for g in genes if pd.notna(g)]
print("S5 genes:", len(genes))
pd.DataFrame({"no": range(1, len(genes) + 1), "gene": genes,
              "category": "inflammation/fibrosis shared gene"}).to_csv(
    REBUILT + "/Table_S5_shared_41genes.csv", index=False)

# ---------- 6) Table_S6 drug-target placeholder ----------
pd.DataFrame([{
    "analysis": "Drug-target MR (candidate druggable genes among shared 41)",
    "status": "占位：待补 — 待药物靶点 MR 代理分析（pQTL/eQTL 工具变量）结果",
    "available_related": "SERPINE1 cis-pQTL MR 已存在于 real_data_results.json (serpine1_pqtl_mr)",
    "date_checked": "2026-07-26"}]).to_csv(REBUILT + "/Table_S6_drugtarget.csv", index=False)

# ---------- 7) Table_S7 gene expression validation (group means) ----------
xl = pd.ExcelFile(BASE + "/data/derived/bulk/GSE57148_GSE57338_keygene_expression_values.xlsx")
def parse_gene_blocks(sheet):
    raw = xl.parse(sheet, header=None)
    blocks = {}
    for col in range(raw.shape[1]):
        v = raw.iat[0, col]
        if isinstance(v, str) and v.strip() and not v.startswith("GSE") and "维" not in v and ")" not in v:
            gene = v.strip()
            # find the Sample/Title/Value triplet starting at or after col
            for c0 in range(col, min(col + 3, raw.shape[1] - 2)):
                if str(raw.iat[1, c0]).strip() == "Sample" and str(raw.iat[1, c0 + 2]).strip() == "Value":
                    sub = raw.iloc[2:, [c0, c0 + 1, c0 + 2]].dropna(how="all")
                    sub.columns = ["Sample", "Title", "Value"]
                    sub = sub.dropna(subset=["Value"])
                    blocks.setdefault(gene, sub)
                    break
    return blocks
s7 = []
for sheet in xl.sheet_names:
    for gene, sub in parse_gene_blocks(sheet).items():
        sub["Value"] = pd.to_numeric(sub["Value"], errors="coerce")
        sub = sub.dropna(subset=["Value"])
        ttl = sub["Title"].astype(str)
        grp = np.where(ttl.str.contains("control|Control|Non-failing|P0", case=False, regex=True), "control", "case")
        for gname, mask in [("control", grp == "control"), ("case", grp == "case")]:
            vals = sub["Value"][mask]
            if len(vals):
                s7.append({"dataset": sheet, "gene": gene, "group": gname, "n": len(vals),
                           "mean": vals.mean(), "sd": vals.std(), "median": vals.median()})
# mislabeled standalone files (content = ACE / THBS1 per task statement)
ace = pd.read_excel(BASE + "/data/derived/bulk/GSE57338_ACE_two_group.xlsx")
s7.append({"dataset": "GSE57338", "gene": "ACE", "group": "control", "n": ace["control"].notna().sum(),
           "mean": ace["control"].mean(), "sd": ace["control"].std(), "median": ace["control"].median()})
s7.append({"dataset": "GSE57338", "gene": "ACE", "group": "case", "n": ace["HF"].notna().sum(),
           "mean": ace["HF"].mean(), "sd": ace["HF"].std(), "median": ace["HF"].median()})
th = pd.read_excel(BASE + "/data/derived/bulk/GSE57148_THBS1_two_group_comparison.xlsx")
th["Value"] = pd.to_numeric(th["Value"], errors="coerce")
for gname, sub in th.dropna(subset=["Value"]).groupby("Title"):
    s7.append({"dataset": "GSE57148", "gene": "THBS1", "group": "control" if "control" in str(gname).lower() else "case",
               "n": len(sub), "mean": sub["Value"].mean(), "sd": sub["Value"].std(), "median": sub["Value"].median()})
s7df = pd.DataFrame(s7).drop_duplicates(subset=["dataset", "gene", "group"], keep="last")
# group label harmonization for GSE57148 (COPD) vs GSE57338 (HF)
s7df = s7df.sort_values(["dataset", "gene", "group"]).reset_index(drop=True)
s7df.to_csv(REBUILT + "/Table_S7_gene_expression_validation.csv", index=False)
print(s7df.to_string(index=False))
print("DONE")
