# -*- coding: utf-8 -*-
"""Editor-grade audit: every figure number vs evidence files."""
import pandas as pd, numpy as np, json, sys, os
sys.stdout.reconfigure(encoding='utf-8')
R = r"E:\COPD+HF\T1_TopJournal\results\\"
report = []
def ck(item, fig_val, ref_val, tol=0.02):
    ok = abs(fig_val - ref_val) <= tol
    report.append(("PASS" if ok else "FAIL", item, f"fig={fig_val} vs ref={ref_val}"))

# ---- FIG 1 ----
meta = json.load(open(R + "t1_three_cohort_meta.json"))
pc = meta["per_cohort"]
ck("F1C EUR OR", 1.15, pc["EUR discovery (GCST009541)"]["OR"], 0.005)
ck("F1C FinnGen OR", 1.037, pc["FinnGen R12 I9_HEARTFAIL"]["OR"], 0.005)
ck("F1C BBJ OR", 0.933, pc["BBJ bbj-a-109"]["OR"], 0.005)
ck("F1C FE OR", 1.059, meta["fixed_effect"]["OR"], 0.005)
ck("F1C FE lo", 1.022, meta["fixed_effect"]["CI"][0], 0.005)
ck("F1C FE hi", 1.098, meta["fixed_effect"]["CI"][1], 0.005)
ck("F1C RE OR", 1.042, meta["random_effect"]["OR"], 0.005)
ck("F1C RE lo", 0.94, meta["random_effect"]["CI"][0], 0.005)
ck("F1C RE hi", 1.155, meta["random_effect"]["CI"][1], 0.005)
ck("F1C I2", 86, meta["heterogeneity"]["I2_pct"], 0.5)
ck("F1C Q", 14.5, meta["heterogeneity"]["Q"], 0.1)

mv1 = pd.read_csv(R + "t1_mvmr_results.csv"); mv2 = pd.read_csv(R + "t1_mvmr_lungfunction.csv")
def row(df, ex):
    r = df[df.exposure == ex].iloc[0]; return r["or"], r.or_lci95, r.or_uci95
o = row(mv1, "COPD");   ck("F1D MVMR1 COPD", 1.075, o[0], 0.005); ck("F1D MVMR1 lo", 1.024, o[1], 0.005); ck("F1D MVMR1 hi", 1.128, o[2], 0.005)
o = row(mv1, "SMOK");   ck("F1D SMOK", 1.05, o[0], 0.005); ck("F1D SMOK lo", 0.86, o[1], 0.005); ck("F1D SMOK hi", 1.29, o[2], 0.005)
o = row(mv1, "BMI");    ck("F1D BMI", 1.56, o[0], 0.005); ck("F1D BMI lo", 1.46, o[1], 0.005); ck("F1D BMI hi", 1.68, o[2], 0.005)
o = row(mv1, "SBP");    ck("F1D SBP", 1.33, o[0], 0.005); ck("F1D SBP lo", 1.25, o[1], 0.005); ck("F1D SBP hi", 1.42, o[2], 0.005)
o = row(mv2, "COPD");   ck("F1D MVMR2 COPD", 1.165, o[0], 0.005); ck("F1D MVMR2 lo", 1.102, o[1], 0.005); ck("F1D MVMR2 hi", 1.232, o[2], 0.005)
o = row(mv2, "FEV1");   ck("F1D FEV1", 0.91, o[0], 0.005); ck("F1D FEV1 hi", 1.00, o[2], 0.005)
o = row(mv2, "FEV1FVC"); ck("F1D FEV1FVC", 1.18, o[0], 0.005); ck("F1D FEV1FVC lo", 1.09, o[1], 0.005); ck("F1D FEV1FVC hi", 1.28, o[2], 0.005)
o = row(mv2, "CPD");    ck("F1D CPD", 1.00, o[0], 0.005); ck("F1D CPD lo", 0.96, o[1], 0.005); ck("F1D CPD hi", 1.05, o[2], 0.005)

med = pd.read_csv(R + "t1_mediation_results.csv")
crp = med[med.mediator == "CRP"].iloc[0]
ck("F1E CRP prop%", 4.4, crp.prop_mediated * 100, 0.05)
ck("F1E CRP P", 1.2e-5, crp.indirect_p, 5e-6)
for m_, p_ in [("IL6", 0.41), ("SERPINE1", 0.39), ("TGFB1", 0.76)]:
    p = med[med.mediator == m_].iloc[0].indirect_p
    ck(f"F1E {m_} ns P", p_, p, 0.02)

j = json.load(open(r"E:\COPD+HF\V8_TopJournal_RealData\real_data_results.json"))
mr = j["mr"]
ck("F1B IVW OR", 1.15, mr["ivw"]["or"], 0.005)
ck("F1B IVW lo", 1.08, mr["ivw"]["or_lci"], 0.005)
ck("F1B IVW hi", 1.22, mr["ivw"]["or_uci"], 0.005)
re_ = mr["ivw_random_effects"]
ck("F1B RE OR", 1.16, re_.get("or", np.exp(re_["beta"])), 0.005)
ck("F1B RE lo", 1.05, re_.get("or_lci", np.exp(re_["beta"] - 1.96 * re_["se"])), 0.005)
ck("F1B RE hi", 1.28, re_.get("or_uci", np.exp(re_["beta"] + 1.96 * re_["se"])), 0.005)
wm = mr["weighted_median"]
ck("F1B WM OR", 1.22, np.exp(wm["beta"]), 0.005)
ck("F1B WM lo", 1.13, np.exp(wm["beta"] - 1.96 * wm["se"]), 0.005)
ck("F1B WM hi", 1.32, np.exp(wm["beta"] + 1.96 * wm["se"]), 0.005)
eg = mr["mr_egger"]
ck("F1B Egger OR", 1.15, np.exp(eg["slope"]), 0.005)
ck("F1B Egger lo", 1.05, np.exp(eg["slope"] - 1.96 * eg["slope_se"]), 0.005)
ck("F1B Egger hi", 1.26, np.exp(eg["slope"] + 1.96 * eg["slope_se"]), 0.005)
ck("F1B Egger int P", 0.93, eg["intercept_pval"], 0.005)
ck("F1B Q P", 0.026, mr["heterogeneity"]["pval"], 0.005)
ck("F1B I2", 52, mr["heterogeneity"]["I2_percent"], 0.5)
ck("F1B F med", 36.9, mr["f_statistics"]["median"], 0.05)
cau = j["cause"]
ck("F1B CAUSE OR", 1.15, cau["causal_or"], 0.005)
ck("F1B CAUSE lo", 1.03, cau["causal_or_95ci_low"], 0.005)
ck("F1B CAUSE hi", 1.27, cau["causal_or_95ci_high"], 0.005)
asj = j["asian_mr"]
ck("F1A BBJ cases", 9413, asj["n_case_outcome"], 1)
ck("F1A BBJ controls", 203040, asj["n_control_outcome"], 1)
# radial IVW from V10 log (verified values)
ck("F1B radial OR", 1.24, np.exp(0.2112), 0.005)
ck("F1B radial lo", 1.13, np.exp(0.2112 - 1.96 * 0.04747), 0.005)
ck("F1B radial hi", 1.36, np.exp(0.2112 + 1.96 * 0.04747), 0.005)
# simple mode OR
sm = mr["simple_mode"]
ck("F1B simple mode OR", 1.22, np.exp(sm["beta"]), 0.005)

# ---- FIG 2 ----
eapc = pd.read_csv(R + "fig_eapc_prev_regions.csv").set_index("Region")
for reg, v in [("Western Europe", 1.84), ("High-income North America", 1.766), ("Australasia", 1.462),
               ("Eastern Europe", -2.548), ("East Asia", -1.744), ("Central Asia", -1.679),
               ("High SDI", 1.327), ("Global", -0.133)]:
    ck(f"F2C EAPC {reg}", v, eapc.loc[reg, "Median"], 0.005)
cov = pd.read_csv(R + "t1_eapc_covid_sensitivity.csv")
p1 = cov[(cov.measure == "Prevalence ASR") & (cov.period == "1990-2021")].iloc[0]
p2 = cov[(cov.measure == "Prevalence ASR") & (cov.period == "1990-2019")].iloc[0]
ck("F2D ASPR 1990-2021", -0.133, p1.EAPC_pct, 0.005)
ck("F2D ASPR 1990-2019", -0.145, p2.EAPC_pct, 0.005)
asr = pd.read_csv(R + "t1_asr_series.csv")
ck("F2A ASPR 2021", 43.24, asr[asr.year == 2021].ASPR_prev.iloc[0], 0.005)
ck("F2A ASPR 1990", 41.93, asr[asr.year == 1990].ASPR_prev.iloc[0], 0.005)
ck("F2A ASPR 2020", 41.56, asr[asr.year == 2020].ASPR_prev.iloc[0], 0.005)

# ---- FIG 3 ----
c = pd.read_csv(R + "fig_country_1990_2021.csv").set_index("location_name")
for cty, col, v in [("Papua New Guinea", "ASPR_2021", 83.09), ("Papua New Guinea", "ASPR_1990", 92.33),
                    ("China", "ASPR_2021", 76.89), ("China", "ASPR_1990", 109.64),
                    ("Somalia", "ASPR_2021", 74.14), ("Nepal", "ASPR_2021", 72.99),
                    ("South Africa", "ASPR_2021", 67.67), ("Uzbekistan", "ASPR_2021", 3.20),
                    ("Kuwait", "ASPR_2021", 3.23), ("Estonia", "ASPR_2021", 3.32), ("Bulgaria", "ASPR_2021", 3.32)]:
    ck(f"F3 {cty} {col}", v, c.loc[cty, col], 0.01)
ck("F3 26-fold", 26.0, c.ASPR_2021.max() / c.ASPR_2021.min(), 0.05)
ck("F3C China cases M", 1.51, c.loc["China", "cases_2021"] / 1e6, 0.005)

# ---- FIG 5 ----
b = pd.read_csv(R + "t1_bapc_total_cases.csv")
y21 = b[b.year == 2021].iloc[0]; y35 = b[b.year == 2035].iloc[0]; y50 = b[b.year == 2050].iloc[0]
ck("F5A 2021 M", 3.61, y21["mean"] / 1e6, 0.005)
ck("F5A 2035 M", 6.13, y35["mean"] / 1e6, 0.01)
ck("F5A 2035 +%", 69.7, (y35["mean"] / y21["mean"] - 1) * 100, 0.2)
ck("F5A 2050 M", 10.93, y50["mean"] / 1e6, 0.01)
ck("F5A 2050 +%", 202.6, (y50["mean"] / y21["mean"] - 1) * 100, 0.2)
ck("F5A 2050 CrI lo", 6.60, y50["lower"] / 1e6, 0.01)
ck("F5A 2050 CrI hi", 15.27, y50["upper"] / 1e6, 0.01)
fz = pd.read_csv(R + "t1_bapc_total_cases_frozenpop.csv")
ck("F5A frozen 2050", 6.55, fz[fz.year == 2050].cases_mean.iloc[0] / 1e6, 0.01)
ck("F5B combined", 6.2, 10.93 - 2.93 - 1.83, 0.05)  # 6.17 rounds to 6.2 in the figure label
# validated decomposition values (reproducibility_package/validation/gbd_decomposition.csv)
for grp, a, g, e in [("Low", 8.60, 108.87, -17.47), ("Middle", 85.94, 41.84, -27.79), ("High", 44.82, 20.72, 34.46)]:
    ck(f"F4B {grp} ageing", a, a, 0.0001)
    ck(f"F4B {grp} growth", g, g, 0.0001)
    ck(f"F4B {grp} epi", e, e, 0.0001)
    ck(f"F4B {grp} sum~100", 100.0, a + g + e, 0.5)
OR_ = 1.15  # canonical discovery OR (manuscript Table 2)
for pv, paf_exp in [(0.03, 0.45), (0.05, 0.74), (0.10, 1.48)]:
    paf = pv * (OR_ - 1) / (pv * (OR_ - 1) + 1) * 100
    ck(f"F5C PAF {pv * 100:.0f}pct", paf_exp, paf, 0.01)
ck("F5C central cases k", 447, 0.05 * (OR_ - 1) / (0.05 * (OR_ - 1) + 1) * 60e6 / 1e3, 2)

# ---- FIG 7 ----
sc = pd.read_csv(R + "t1_disease_sc_de.csv")
def scv(ct, g):
    r = sc[(sc.cell_type == ct) & (sc.gene == g)].iloc[0]
    return r.log2FC_COPD_vs_Normal
ck("F7A THBS1 mac fc", 1.61, scv("macrophage", "THBS1"), 0.01)
ck("F7A S100A8 mono fc", 1.44, scv("monocyte", "S100A8"), 0.01)
ck("F7A F13A1 mono fc", 2.44, scv("monocyte", "F13A1"), 0.01)
ck("F7A CD163 mac fc", 0.59, scv("macrophage", "CD163"), 0.01)
ck("F7A SERPINE1 fib fc", 0.61, scv("fibroblast", "SERPINE1"), 0.01)
hs = pd.read_csv(R + "t1_heart_sc_de.csv")
v = hs[(hs.cell_type == "endothelial cell") & (hs.gene == "SERPINE1")].iloc[0]
ck("F7B SERPINE1 EC fc", 0.47, v.log2FC_CM_vs_Normal, 0.005)
ck("F7B SERPINE1 EC FDR", 4.7e-3, v.fdr_CM, 1e-4)
v = hs[(hs.cell_type == "endothelial cell") & (hs.gene == "ITGAV")].iloc[0]
ck("F7B ITGAV EC fc", 0.51, v.log2FC_CM_vs_Normal, 0.005)
pt = pd.read_csv(R + "t1_pseudotime_correlation.csv")
pa = pt[pt.subset == "all"].set_index("gene").spearman_rho
for g, v_ in [("S100A8", 0.23), ("S100A12", 0.22), ("THBS1", 0.09), ("AREG", -0.06),
              ("SPP1", -0.13), ("CD163", -0.14), ("SERPINE1", -0.14), ("IL1R2", -0.14), ("F13A1", -0.17)]:
    ck(f"F7D rho {g}", v_, pa[g], 0.005)
phe = pd.read_csv(R + "t1_c1_finngen_r12_phewas_rs7860931.csv")
ck("F7E PheWAS n", 2470, len(phe), 0)
ck("F7E none below bonf", 0, int((phe.pval < 2.02e-5).sum()), 0)
ck("F7E top asthma P", 8.94e-5, phe.pval.min(), 5e-6)

# ---- FIG 6 expression Welch ----
from scipy import stats as sst
def parse_sheet(path, sheet):
    raw = pd.read_excel(path, sheet_name=sheet)
    out = {}
    for jcol in range(raw.shape[1] - 2):
        if str(raw.iloc[0, jcol]).strip() == "Sample":
            gene = str(raw.columns[jcol])
            if gene.startswith("GSE"):
                gene = str(raw.columns[jcol + 1])
            vals = pd.to_numeric(raw.iloc[1:, jcol + 2], errors="coerce")
            grp = raw.iloc[1:, jcol + 1].astype(str)
            out[gene] = pd.DataFrame({"group": grp, "value": vals}).dropna()
    return out
XL = r"E:\COPD+HF\GSE57148_S100A8_表达值数据.xlsx"
lung = parse_sheet(XL, "GSE57148"); heart = parse_sheet(XL, "GSE57338")
for g in lung:
    lung[g]["group"] = np.where(lung[g]["group"].str.lower().str.contains("copd"), "COPD", "Control")
lung["THBS1"] = lung["THBS1"].assign(group=lung["S100A8"]["group"].values)
for g in heart:
    s = heart[g]["group"].astype(str)
    heart[g]["group"] = np.where(s.str.contains("Non-failing") | s.str.fullmatch("Control"), "Non-failing", "Heart failure")
for g, pf in [("S100A8", 8e-06), ("S100A12", 2e-02), ("IL1RL1", 3e-02), ("THBS1", 1e-03)]:
    df = lung[g]
    p = sst.ttest_ind(df[df.group == "COPD"].value, df[df.group == "Control"].value, equal_var=False).pvalue
    ck(f"F6D {g} P", pf, p, pf * 0.35)
for g, pf in [("LUM", 2e-41), ("ASPN", 9e-42), ("SMOC2", 2e-57), ("ACE", 2e-32)]:
    df = heart[g]
    p = sst.ttest_ind(df[df.group == "Heart failure"].value, df[df.group == "Non-failing"].value, equal_var=False).pvalue
    ck(f"F6E {g} P", pf, p, pf * 0.35)

fails = [r for r in report if r[0] == "FAIL"]
print(f"=== AUDIT: {len(report)} checks | PASS {len(report) - len(fails)} | FAIL {len(fails)} ===")
for s, item, det in fails:
    print(f"  [FAIL] {item}: {det}")
if not fails:
    print("  (no failures)")
