#!/usr/bin/env python3
"""M5: ancestry-stratified PAF sensitivity (Levin's formula) for COPD-attributable HF.
Discovery OR 1.15 (published) vs European subgroup 1.036 vs East Asian subgroup 0.922.
"""
import math, csv, os
os.makedirs("results_new", exist_ok=True)

def paf(rr, p):
    return p * (rr - 1) / (p * (rr - 1) + 1)

scenarios = [
    ("Discovery IVW (published; EUR discovery)", 1.15),
    ("European subgroup FE (k=6)", 1.036),
    ("East Asian subgroup FE (k=4)", 0.922),
    ("Overall FE (k=10)", 1.013),
]
HF_POOL = 60_000_000
rows = []
for name, orv in scenarios:
    for prev in [0.03, 0.05, 0.10]:
        pf = paf(orv, prev)
        rows.append({"scenario": name, "OR": orv, "copd_prevalence": prev,
                     "PAF": pf, "PAF_pct": pf * 100,
                     "implied_cases_of_60M_HF": round(pf * HF_POOL)})
with open("results_new/ancestry_paf_sensitivity.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
for r in rows:
    print(f"{r['scenario']:42s} prev={r['copd_prevalence']:.0%}  PAF={r['PAF_pct']:+.2f}%  cases={r['implied_cases_of_60M_HF']:,}")
print("\nKey readings:")
print("- East Asian OR<1 -> PAF negative: under the EAS estimate, COPD liability is NOT attributable-risk positive;")
print("  GBD East-Asia burden therefore reflects observational/GBD-framework attribution, not MR-implied causality.")
print("- European subgroup 1.036 -> PAF 0.11-0.36% at 3-10% prevalence (vs 0.45-1.48% at OR 1.15).")
print("DONE")
