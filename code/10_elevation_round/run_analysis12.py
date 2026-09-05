#!/usr/bin/env python3
"""Analysis 1+2: Winner's curse correction (discovery) + BBJ 2025 HFpEF/HFrEF/allHF stratified MR.
Outputs CSVs into results_new/.
"""
import csv, json, math, os
from mr_engine import full_mr, winners_curse_conditional, harmonise_outcome

os.makedirs("results_new", exist_ok=True)
base = list(csv.DictReader(open("repo/results/t1_base_harmonised_10snp.csv")))
bx  = [float(r["beta.exposure"]) for r in base]
bse = [float(r["se.exposure"]) for r in base]
by0 = [float(r["beta.outcome"]) for r in base]
sy0 = [float(r["se.outcome"]) for r in base]
snps = [r["SNP"] for r in base]
eaf_x = [float(r["eaf.exposure"]) for r in base]
ea_x = [r["effect_allele.exposure"] for r in base]
oa_x = [r["other_allele.exposure"] for r in base]

# ---------- 2. Winner's curse ----------
wc_rows = []
bx_corr = []
for rs, b, se in zip(snps, bx, bse):
    z = b / se
    corr = winners_curse_conditional(b, se)
    bx_corr.append(corr)
    wc_rows.append({"SNP": rs, "beta_naive": b, "se": se, "Z": z,
                    "beta_corrected": corr, "shrinkage": corr / b})
res_wc_naive = full_mr(bx, by0, sy0)
res_wc_corr = full_mr(bx_corr, by0, sy0)
with open("results_new/winners_curse_persnp.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=wc_rows[0].keys()); w.writeheader(); w.writerows(wc_rows)
print("== Winner's curse (Zhong-Prentice conditional, threshold 5e-8) ==")
for r in wc_rows:
    print(f"  {r['SNP']}: {r['beta_naive']:.4f} -> {r['beta_corrected']:.4f} (x{r['shrinkage']:.3f})")
print(f"  IVW naive:    OR={math.exp(res_wc_naive['ivw_fe']):.4f} ({math.exp(res_wc_naive['ivw_fe']-1.96*res_wc_naive['ivw_fe_se']):.4f}-{math.exp(res_wc_naive['ivw_fe']+1.96*res_wc_naive['ivw_fe_se']):.4f})")
print(f"  IVW corrected:OR={math.exp(res_wc_corr['ivw_fe']):.4f} ({math.exp(res_wc_corr['ivw_fe']-1.96*res_wc_corr['ivw_fe_se']):.4f}-{math.exp(res_wc_corr['ivw_fe']+1.96*res_wc_corr['ivw_fe_se']):.4f})")

# ---------- 1. BBJ 2025 stratified MR ----------
def load_bbj(tag):
    rows = list(csv.DictReader(open(f"bbj2025_extracts/{tag}.tsv"), delimiter="\t"))
    d = {}
    for r in rows:
        d[r["rsid"]] = r
    return d

strat_results = []
harm_rows = []
for tag, label, ncases in [("BBJ_HFpEF", "BBJ 2025 HFpEF (7,154 cases)", 7154),
                            ("BBJ_HFrEF", "BBJ 2025 HFrEF (4,254 cases)", 4254),
                            ("BBJ_allHF", "BBJ 2025 all-cause HF (16,251 cases)", 16251)]:
    d = load_bbj(tag)
    bx2, by2, sy2, kept = [], [], [], []
    for rs, bx_i, ea, oa, eaf in zip(snps, bx, ea_x, oa_x, eaf_x):
        if rs not in d:
            harm_rows.append({"outcome": tag, "SNP": rs, "action": "missing_in_outcome"}); continue
        r = d[rs]
        af_y = float(r["effect_allele_frequency"])
        b, act = harmonise_outcome(ea, oa, eaf, r["effect_allele"], r["other_allele"],
                                   float(r["beta"]), af_y)
        harm_rows.append({"outcome": tag, "SNP": rs, "action": act,
                          "ea_x": ea, "oa_x": oa, "ea_y": r["effect_allele"], "oa_y": r["other_allele"],
                          "beta_y_raw": r["beta"], "beta_aligned": b})
        if b is None: continue
        bx2.append(bx_i); by2.append(b); sy2.append(float(r["standard_error"])); kept.append(rs)
    res = full_mr(bx2, by2, sy2)
    strat_results.append({"outcome": label, "tag": tag, "nsnp": res["n"],
                          **{k: v for k, v in res.items()}})
    print(f"\n== {label}: n={res['n']} SNPs ==")
    print(f"  IVW-FE OR={math.exp(res['ivw_fe']):.4f} ({math.exp(res['ivw_fe']-1.96*res['ivw_fe_se']):.4f}-{math.exp(res['ivw_fe']+1.96*res['ivw_fe_se']):.4f}) P={res['ivw_fe_p']:.3g}")
    print(f"  IVW-RE OR={math.exp(res['ivw_re']):.4f} ({math.exp(res['ivw_re']-1.96*res['ivw_re_se']):.4f}-{math.exp(res['ivw_re']+1.96*res['ivw_re_se']):.4f}) P={res['ivw_re_p']:.3g} I2={res['I2']*100:.0f}%")
    print(f"  WM    OR={math.exp(res['wm']):.4f} P={res['wm_p']:.3g}; Egger slope OR={math.exp(res['slope']):.4f} P={res['slope_p']:.3g}, int P={res['intercept_p']:.3g}")

fields = ["outcome","tag","nsnp","ivw_fe","ivw_fe_se","ivw_fe_p","ivw_re","ivw_re_se","ivw_re_p",
          "Q","Q_df","Q_p","I2","tau2","intercept","intercept_se","intercept_p","slope","slope_se","slope_p","wm","wm_se","wm_p"]
with open("results_new/bbj2025_stratified_mr.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fields); w.writeheader()
    for r in strat_results: w.writerow({k: r.get(k) for k in fields})
with open("results_new/bbj2025_harmonisation_log.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["outcome","SNP","action","ea_x","oa_x","ea_y","oa_y","beta_y_raw","beta_aligned"])
    w.writeheader(); w.writerows(harm_rows)

# summary JSON incl winner's curse
json.dump({
    "winners_curse": {"per_snp": wc_rows,
                      "ivw_naive": {k: res_wc_naive[k] for k in ["ivw_fe","ivw_fe_se","ivw_fe_p"]},
                      "ivw_corrected": {k: res_wc_corr[k] for k in ["ivw_fe","ivw_fe_se","ivw_fe_p"]}},
    "bbj2025_stratified": strat_results,
}, open("results_new/analysis12_summary.json", "w"), indent=1)
print("\nSaved results_new/analysis12_summary.json")
