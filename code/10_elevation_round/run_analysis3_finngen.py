#!/usr/bin/env python3
"""Analysis 3: MR for FinnGen R12 additional endpoints (HF composites) using the 10-SNP COPD set.
Also extracts FinnGen COPD-endpoint betas for the alternative-exposure (winner's curse) re-estimation.
"""
import csv, json, math, os, sys
from mr_engine import full_mr, harmonise_outcome

os.makedirs("results_new", exist_ok=True)
base = list(csv.DictReader(open("repo/results/t1_base_harmonised_10snp.csv")))
snps  = [r["SNP"] for r in base]
bx    = [float(r["beta.exposure"]) for r in base]
bse_x = [float(r["se.exposure"]) for r in base]
ea_x  = [r["effect_allele.exposure"] for r in base]
oa_x  = [r["other_allele.exposure"] for r in base]
eaf_x = [float(r["eaf.exposure"]) for r in base]
by_disc = [float(r["beta.outcome"]) for r in base]
sy_disc = [float(r["se.outcome"]) for r in base]

def load_fg(ep):
    rows = list(csv.DictReader(open(f"finngen_extracts/{ep}.tsv"), delimiter="\t"))
    out = {}
    for r in rows:
        rslist = r["rsids"].split(",")
        hit = [x for x in rslist if x in snps]
        if hit:
            out[hit[0]] = r
    return out

META = {
    "I9_HEARTFAIL_AND_CHD":       ("FinnGen R12 HF and CHD (17,532 cases)", "HF", "EUR-Finnish", 17532),
    "I9_HEARTFAIL_AND_OVERWEIGHT":("FinnGen R12 HF and BMI>=25 (22,965 cases)", "HF", "EUR-Finnish", 22965),
    "FG_OTHHEART":                ("FinnGen R12 other heart diseases (106,069 cases)", "HF-broad", "EUR-Finnish", 106069),
    "I9_OTHILLHEART":             ("FinnGen R12 other/ill-defined heart diseases (2,192 cases)", "HF-broad", "EUR-Finnish", 2192),
    "J10_COPD":                   ("FinnGen R12 COPD (24,138 cases)", "COPD", "EUR-Finnish", 24138),
    "COPD_LATER":                 ("FinnGen R12 later-onset COPD (13,535 cases)", "COPD", "EUR-Finnish", 13535),
    "COPD_EARLY":                 ("FinnGen R12 early-onset COPD (9,491 cases)", "COPD", "EUR-Finnish", 9491),
}

results = []
harm = []
fg_copd_betas = {}   # rs -> dict per COPD endpoint
for ep, (label, kind, ancestry, ncases) in META.items():
    d = load_fg(ep)
    print(f"\n== {ep} ({label}): {len(d)}/10 SNPs found")
    bx2, by2, sy2 = [], [], []
    for i, rs in enumerate(snps):
        if rs not in d:
            harm.append({"endpoint": ep, "SNP": rs, "action": "missing"}); continue
        r = d[rs]
        b, act = harmonise_outcome(ea_x[i], oa_x[i], eaf_x[i], r["alt"], r["ref"],
                                   float(r["beta"]), float(r["af_alt"]) if r["af_alt"] else None)
        harm.append({"endpoint": ep, "SNP": rs, "action": act, "beta_raw": r["beta"], "beta_aligned": b})
        if b is None: continue
        if kind == "COPD":
            fg_copd_betas.setdefault(rs, {})[ep] = {"beta": b, "se": float(r["sebeta"])}
        else:
            bx2.append(bx[i]); by2.append(b); sy2.append(float(r["sebeta"]))
    if kind != "COPD" and len(bx2) >= 4:
        res = full_mr(bx2, by2, sy2)
        results.append({"endpoint": ep, "label": label, "ancestry": ancestry, "ncases": ncases, **res})
        print(f"  IVW-FE OR={math.exp(res['ivw_fe']):.4f} ({math.exp(res['ivw_fe']-1.96*res['ivw_fe_se']):.4f}-{math.exp(res['ivw_fe']+1.96*res['ivw_fe_se']):.4f}) P={res['ivw_fe_p']:.3g} I2={res['I2']*100:.0f}%")

# FinnGen COPD as independent exposure -> discovery HF outcome (winner's curse alternative)
alt_rows = []
for ep in ["J10_COPD", "COPD_LATER", "COPD_EARLY"]:
    bx2, by2, sy2, used = [], [], [], []
    for i, rs in enumerate(snps):
        ent = fg_copd_betas.get(rs, {}).get(ep)
        if ent is None: continue
        bx2.append(ent["beta"]); by2.append(by_disc[i]); sy2.append(sy_disc[i]); used.append(rs)
    if len(bx2) >= 4:
        res = full_mr(bx2, by2, sy2)
        alt_rows.append({"exposure": f"FinnGen {ep}", "outcome": "Discovery HF (HERMES/ebi-a-GCST009541)", **res})
        print(f"\n== Alt-exposure {ep} -> discovery HF: n={res['n']}, IVW-FE OR={math.exp(res['ivw_fe']):.4f} ({math.exp(res['ivw_fe']-1.96*res['ivw_fe_se']):.4f}-{math.exp(res['ivw_fe']+1.96*res['ivw_fe_se']):.4f}) P={res['ivw_fe_p']:.3g}")

fields = ["endpoint","label","ancestry","ncases","n","ivw_fe","ivw_fe_se","ivw_fe_p","ivw_re","ivw_re_se","ivw_re_p",
          "Q","Q_df","Q_p","I2","tau2","intercept","intercept_se","intercept_p","slope","slope_se","slope_p","wm","wm_se","wm_p"]
with open("results_new/finngen_extended_mr.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fields); w.writeheader()
    for r in results: w.writerow({k: r.get(k) for k in fields})
with open("results_new/finngen_copd_alt_exposure_mr.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["exposure","outcome"]+fields[4:]); w.writeheader()
    for r in alt_rows: w.writerow({k: r.get(k) for k in ["exposure","outcome"]+fields[4:]})
with open("results_new/finngen_harmonisation_log.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["endpoint","SNP","action","beta_raw","beta_aligned"]); w.writeheader(); w.writerows(harm)
print("\nDONE")
