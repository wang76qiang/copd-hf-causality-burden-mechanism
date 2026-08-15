# -*- coding: utf-8 -*-
# ============================================================================
# Module 05 - GBD: GBD-vs-MR bridge analysis (Supplementary Note 2).
# Compares the GBD comparative-risk attribution (3.61 M cases; implied RR and
# PAF against an assumed ~60 M global HF pool) with the MR-Levin PAF based on
# the discovery OR 1.15 at 3/5/10% assumed COPD prevalence.
# Output: results/t1_bridge_table.csv (+ manuscript numbers printed to stdout)
# ============================================================================

import os as _os
# --- repository paths -----------------------------------------------------------
# REPO resolves to the repository root (code/<module>/this_file.py -> ../..).
# Override with the REPRO_ROOT environment variable if you run from elsewhere.
REPO = _os.environ.get("REPRO_ROOT", _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "..")))
# GBD_RAW_DIR: folder with the raw IHME GBD 2021 csv downloads (see data/README.md)
GBD_RAW_DIR = _os.environ.get("GBD_RAW_DIR", _os.path.join(REPO, "data", "raw", "gbd_2021"))
"""T1 B4: GBD-MR bridge analysis -> results/t1_bridge_table.csv (+ prints md body numbers)."""
import pandas as pd, numpy as np, json

OUT = REPO
GBD_ATTR = 3613136.3631          # COPD-attributed HF cases 2021 (verified from raw csv)
HF_TOTAL = 60_000_000            # ASSUMPTION: global HF prevalence ~60 million
OR_MR = 1.1507944875580642       # IVW FE
P_COPD = {"low": 0.03, "central": 0.05, "high": 0.10}

rows = []
# --- GBD framework implied quantities ---
paf_gbd = GBD_ATTR / HF_TOTAL
for k, p in P_COPD.items():
    rr = 1 + paf_gbd / (p * (1 - paf_gbd))          # invert Levin: PAF = p(RR-1)/(1+p(RR-1))
    rows.append({"framework": "GBD comparative risk assessment",
                 "copd_prevalence_assumed": p, "scenario": k,
                 "effect_measure": "implied RR", "effect": rr,
                 "PAF": paf_gbd, "attributable_HF_cases_2021": GBD_ATTR,
                 "note": "implied PAF = GBD attributable / assumed 60M total HF"})
# --- MR-Levin framework ---
for k, p in P_COPD.items():
    paf = p * (OR_MR - 1) / (1 + p * (OR_MR - 1))
    rows.append({"framework": "MR + Levin PAF",
                 "copd_prevalence_assumed": p, "scenario": k,
                 "effect_measure": "OR (IVW FE)", "effect": OR_MR,
                 "PAF": paf, "attributable_HF_cases_2021": paf * HF_TOTAL,
                 "note": "PAF = p(OR-1)/(1+p(OR-1)); cases vs 60M total HF"})
df = pd.DataFrame(rows)
df.to_csv(OUT + "/results/t1_bridge_table.csv", index=False)
pd.set_option("display.width", 200)
print(df.to_string(index=False))
print("GBD implied PAF:", paf_gbd)
print("GBD/MR case ratio (central):", GBD_ATTR / (0.007483302340816935 * HF_TOTAL))
