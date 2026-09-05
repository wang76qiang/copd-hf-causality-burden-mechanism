#!/usr/bin/env python3
"""Analysis 4: heterogeneity decomposition across all cohort/endpoint MR estimates.
- Assembles discovery (naive + winner's-curse-corrected), FinnGen strict + new endpoints, BBJ legacy + BBJ2025 strata
- Fixed/random meta-analysis overall and by ancestry and endpoint-breadth subgroups
- Illustrative meta-regression (k small): logOR ~ ancestry(EAS) + endpoint breadth
- Cohort-level funnel plot data + Egger asymmetry test (illustrative)
"""
import csv, json, math, os
import numpy as np
from scipy import stats

os.makedirs("results_new", exist_ok=True)

def or_ci(orv, lo, hi):
    b = math.log(orv); se = (math.log(hi) - math.log(lo)) / (2*1.959964)
    return b, se

cohorts = []
# discovery
b, se = or_ci(1.15, 1.08, 1.22)
cohorts.append(dict(id="Discovery (HERMES/UKB+MVP, EUR)", ancestry="EUR", breadth="composite", cases=47309, beta=b, se=se, note="published IVW-FE"))
# discovery winner's-curse corrected
d12 = json.load(open("results_new/analysis12_summary.json"))
wc = d12["winners_curse"]["ivw_corrected"]
cohorts.append(dict(id="Discovery, winner's-curse-corrected (EUR)", ancestry="EUR", breadth="composite", cases=47309,
                    beta=wc["ivw_fe"], se=wc["ivw_fe_se"], note="Zhong-Prentice conditional correction", pool=False))
# FinnGen strict (published)
b, se = or_ci(1.037, 0.986, 1.090)
cohorts.append(dict(id="FinnGen R12 strict HF", ancestry="EUR", breadth="strict", cases=37653, beta=b, se=se, note="published"))
# BBJ legacy congestive HF
b, se = or_ci(0.9329, 0.8477173501500259, 1.0267086605788982)
cohorts.append(dict(id="BBJ congestive HF (legacy 6-SNP)", ancestry="EAS", breadth="composite", cases=9413, beta=b, se=se, note="published"))
# FinnGen new endpoints
fg = list(csv.DictReader(open("results_new/finngen_extended_mr.csv")))
for r in fg:
    cohorts.append(dict(id=r["label"], ancestry="EUR", breadth=("broad" if "other" in r["label"].lower() else "composite"),
                        cases=int(r["ncases"]), beta=float(r["ivw_fe"]), se=float(r["ivw_fe_se"]),
                        note=f"I2={float(r['I2'])*100:.0f}%"))
# BBJ 2025 strata
bb = list(csv.DictReader(open("results_new/bbj2025_stratified_mr.csv")))
for r in bb:
    breadth = {"BBJ_allHF": "composite", "BBJ_HFpEF": "HFpEF", "BBJ_HFrEF": "HFrEF"}[r["tag"]]
    cohorts.append(dict(id=r["outcome"], ancestry="EAS", breadth=breadth,
                        cases={"BBJ_allHF":16251,"BBJ_HFpEF":7154,"BBJ_HFrEF":4254}[r["tag"]],
                        beta=float(r["ivw_fe"]), se=float(r["ivw_fe_se"]), note=f"I2={float(r['I2'])*100:.0f}%"))

def meta(estimates):
    b = np.array([e["beta"] for e in estimates]); s = np.array([e["se"] for e in estimates])
    w = 1/s**2
    est_f = float(np.sum(w*b)/np.sum(w)); se_f = math.sqrt(1/np.sum(w))
    Q = float(np.sum(w*(b-est_f)**2)); df = len(b)-1
    Q_p = float(stats.chi2.sf(Q, df)); I2 = max(0, (Q-df)/Q)
    tau2 = max(0.0, (Q-df)/(np.sum(w)-np.sum(w**2)/np.sum(w)))
    w2 = 1/(s**2+tau2); est_r = float(np.sum(w2*b)/np.sum(w2)); se_r = math.sqrt(1/np.sum(w2))
    return dict(k=len(b), fe_beta=est_f, fe_se=se_f, fe_p=float(2*stats.norm.sf(abs(est_f/se_f))),
                re_beta=est_r, re_se=se_r, re_p=float(2*stats.norm.sf(abs(est_r/se_r))),
                Q=Q, Q_df=df, Q_p=Q_p, I2=I2, tau2=tau2)

results = {"all": meta([c for c in cohorts if c.get("pool", True)])}
for anc in ["EUR", "EAS"]:
    sub = [c for c in cohorts if c["ancestry"] == anc and c.get("pool", True)]
    results[f"ancestry_{anc}"] = meta(sub)
for br in ["strict", "composite", "broad", "HFpEF", "HFrEF"]:
    sub = [c for c in cohorts if c["breadth"] == br and c.get("pool", True)]
    if len(sub) >= 2:
        results[f"breadth_{br}"] = meta(sub)

# meta-regression (illustrative): logOR ~ EAS + strict + HFpEF
pooled = [c for c in cohorts if c.get("pool", True)]
y = np.array([c["beta"] for c in pooled])
X = np.column_stack([
    np.ones(len(pooled)),
    [1 if c["ancestry"] == "EAS" else 0 for c in pooled],
    [1 if c["breadth"] == "strict" else 0 for c in pooled],
    [1 if c["breadth"] in ("HFpEF",) else 0 for c in pooled],
])
w = np.array([1/c["se"]**2 for c in pooled])
W = np.diag(w)
beta_hat = np.linalg.solve(X.T @ W @ X, X.T @ W @ y)
resid = y - X @ beta_hat
Q_res = float(np.sum(w * resid**2)); df_res = len(y) - X.shape[1]
cov = np.linalg.inv(X.T @ W @ X) * max(1.0, Q_res/df_res)
se_hat = np.sqrt(np.diag(cov))
p_hat = 2*stats.t.sf(np.abs(beta_hat/se_hat), df_res)
# QM test for moderators
beta0 = np.linalg.solve(X[:, :1].T @ W @ X[:, :1], X[:, :1].T @ W @ y)
resid0 = y - X[:, :1] @ beta0
Q_model = float(np.sum(w*resid0**2) - np.sum(w*resid**2))
QM_p = float(stats.chi2.sf(Q_model, 3))
R2 = 1 - np.sum(w*resid**2)/np.sum(w*resid0**2)
results["meta_regression"] = dict(
    coef_names=["intercept", "EastAsian", "strict_endpoint", "HFpEF_endpoint"],
    coefs=beta_hat.tolist(), ses=se_hat.tolist(), ps=p_hat.tolist(),
    QM=Q_model, QM_df=3, QM_p=QM_p, R2_analog=float(R2), k=len(y))

# funnel asymmetry (Egger, cohort-level, illustrative)
snd = y / np.array([c["se"] for c in pooled])
prec = 1 / np.array([c["se"] for c in pooled])
slope, intercept, r, p_eg, _ = stats.linregress(prec, snd)
results["funnel_egger"] = dict(intercept=float(intercept), slope=float(slope), p=float(p_eg), k=len(y))

# save
json.dump({"cohorts": cohorts, "meta": results}, open("results_new/heterogeneity_package.json", "w"), indent=1)
with open("results_new/heterogeneity_cohorts.csv", "w", newline="") as f:
    wcsv = csv.DictWriter(f, fieldnames=["id","ancestry","breadth","cases","beta","se","note","pool"], extrasaction="ignore")
    wcsv.writeheader(); wcsv.writerows(cohorts)

print(f"k = {len(cohorts)} cohort estimates")
for name, m in results.items():
    if name in ("meta_regression", "funnel_egger"): continue
    print(f"{name:20s} k={m['k']:2d} FE OR={math.exp(m['fe_beta']):.3f} (P={m['fe_p']:.3g})  RE OR={math.exp(m['re_beta']):.3f} (P={m['re_p']:.3g})  I2={m['I2']*100:.0f}%")
mr = results["meta_regression"]
print("\nMeta-regression (WLS, multiplicative dispersion):")
for n, c_, s_, p_ in zip(mr["coef_names"], mr["coefs"], mr["ses"], mr["ps"]):
    print(f"  {n:16s} {c_:+.4f} (SE {s_:.4f}) P={p_:.3g}")
print(f"  QM={mr['QM']:.2f} df=3 P={mr['QM_p']:.3g}; R2-analog={mr['R2_analog']*100:.0f}%")
print("Funnel Egger (cohort-level, illustrative): intercept P =", f"{results['funnel_egger']['p']:.3g}")
print("DONE")
