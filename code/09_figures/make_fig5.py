# -*- coding: utf-8 -*-
# ============================================================================
# Module 09 - Figure 5: the 2050 horizon (BAPC projection, scenarios, PAF)
# NOTE: panel B scenario values are computed by
# code/05_gbd/gbd_intervention_scenarios.py; see VALIDATION.md for the
# numerical verification.
# Inputs: results/t1_bapc_total_cases{,_frozenpop}.csv, t1_bridge_table.csv
# Output: figures/Fig5_horizon_2050.{png,pdf}
# ============================================================================

import os as _os
# --- repository paths -----------------------------------------------------------
# REPO resolves to the repository root (code/<module>/this_file.py -> ../..).
# Override with the REPRO_ROOT environment variable if you run from elsewhere.
REPO = _os.environ.get("REPRO_ROOT", _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "..")))
# GBD_RAW_DIR: folder with the raw IHME GBD 2021 csv downloads (see data/README.md)
GBD_RAW_DIR = _os.environ.get("GBD_RAW_DIR", _os.path.join(REPO, "data", "raw", "gbd_2021"))
"""Figure 5 — The 2050 horizon: BAPC projection, scenarios, preventable fraction."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from figstyle import *
import pandas as pd
import numpy as np

R = os.path.join(REPO, "results") + os.sep
# MedComm revision: write into the working-copy folder, not the source repo
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")

bapc = pd.read_csv(R + "t1_bapc_total_cases.csv")
frozen = pd.read_csv(R + "t1_bapc_total_cases_frozenpop.csv")

fig = plt.figure(figsize=(13.5, 9.0))
gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 1], hspace=0.5, wspace=0.28,
                      left=0.08, right=0.96, top=0.87, bottom=0.09)

# ---------- A: fan chart ----------
axA = fig.add_subplot(gs[0, :])
obs = bapc[bapc.year <= 2021]; proj = bapc[bapc.year >= 2021]
axA.fill_between(proj.year, proj.lower / 1e6, proj.upper / 1e6,
                 color=CRIMSON, alpha=0.16, zorder=1, label="95% CrI")
axA.plot(obs.year, obs["mean"] / 1e6, color=INK, lw=2.6, zorder=4,
         label="Observed (GBD 2021)")
axA.plot(proj.year, proj["mean"] / 1e6, color=CRIMSON, lw=2.6, zorder=4,
         label="BAPC projection")
axA.plot(frozen.year, frozen.cases_mean / 1e6, color=SKY, lw=1.8, ls="--",
         zorder=3, label="Frozen-population sensitivity")
axA.axvline(2021, color=SLATE, ls=":", lw=1.1)
p2035 = bapc[bapc.year == 2035].iloc[0]; p2050 = bapc[bapc.year == 2050].iloc[0]
for p, col in [(p2035, AMBER), (p2050, CRIMSON)]:
    axC_ = axA.scatter([p.year], [p["mean"] / 1e6], s=90, color=col, zorder=5,
                       edgecolor="white", lw=1.2)
axA.annotate(f"2035: {p2035['mean']/1e6:.1f} M (+70%)", (2035, p2035["mean"] / 1e6),
             xytext=(2022, 9.4), fontsize=9, fontweight="bold", color=AMBER,
             arrowprops=dict(arrowstyle="-|>", color=AMBER, lw=1.2))
axA.annotate(f"2050: {p2050['mean']/1e6:.1f} M (+203%)\nCrI {p2050['lower']/1e6:.1f}–{p2050['upper']/1e6:.1f} M",
             (2050, p2050["mean"] / 1e6), xytext=(2037, 13.2), fontsize=9.5,
             fontweight="bold", color=CRIMSON,
             arrowprops=dict(arrowstyle="-|>", color=CRIMSON, lw=1.2))

axA.set_xlim(1990, 2051); axA.set_ylim(0, 16.5)
axA.set_xlabel("Year"); axA.set_ylabel("COPD-attributable HF cases (millions)")
axA.set_title("Bayesian age-period-cohort projection to 2050", loc="left")
axA.legend(loc="upper left", fontsize=8.5)
strip_spines(axA); grid(axA)
panel_tag(axA, "A", CRIMSON, x=-0.055)

# ---------- B: scenario waterfall ----------
axB = fig.add_subplot(gs[1, 0])
# scenario values from gbd_intervention_scenarios.py (validated): 2050 cases averted
# smoking cessation RR 0.40: 2.93 M; optimized pharmacotherapy RR 0.25: 1.83 M
base = 10.93; smoke = -2.93; treat = -1.83
labels = ["2050\nbaseline", "Smoking\ncessation", "Optimized COPD\npharmacotherapy", "Combined\npathway"]
comb = base + smoke + treat
starts = [0, base + smoke, base + smoke + treat, 0]
heights = [base, -smoke, -treat, comb]
cols = [INK, TEAL, TEAL, GREEN]
axB.bar(0, base, color=INK, width=0.62, zorder=3)
axB.bar(1, -smoke, bottom=base + smoke, color=TEAL, width=0.62, zorder=3)
axB.bar(2, -treat, bottom=base + smoke + treat, color=TEAL, width=0.62, zorder=3, alpha=0.75)
axB.bar(3, comb, color=GREEN, width=0.62, zorder=3)
# value labels above each floating bar (verified: 10.93 - 2.93 - 1.83 = 6.17 M combined)
for x, v, t in [(0, base, f"{base:.1f} M"), (1, base, "−2.9 M"),
                (2, base + smoke, "−1.8 M"), (3, comb, f"{comb:.1f} M")]:
    axB.text(x, v + 0.28, t, ha="center", fontsize=9.5, fontweight="bold",
             color="black")
axB.set_xticks(range(4)); axB.set_xticklabels(labels, fontsize=8)
axB.set_ylabel("Cases in 2050 (millions)")
axB.set_ylim(0, 12.5)
axB.set_title("What prevention can still buy", loc="left")
strip_spines(axB); grid(axB)
panel_tag(axB, "B", GREEN)

# ---------- C: PAF curve ----------
axC = fig.add_subplot(gs[1, 1])
OR = 1.15; HF_POOL = 60e6
p = np.linspace(0.01, 0.15, 200)
paf = p * (OR - 1) / (p * (OR - 1) + 1)
cases = paf * HF_POOL / 1e3
axC.plot(p * 100, cases, color=CRIMSON, lw=2.6, zorder=3)
axC.fill_between(p * 100, 0, cases, color=CRIMSON, alpha=0.10, zorder=1)
# MedComm revision: labels moved into the empty band between the PAF curve and
# the GBD attribution line, with leader lines, so no text touches curve/points.
offsets = {3: (0.9, 1250), 5: (2.6, 2100), 10: (6.0, 2950)}
for pv, lab in [(3, "3%"), (5, "5% central"), (10, "10%")]:
    pf = pv / 100 * (OR - 1) / (pv / 100 * (OR - 1) + 1)
    cs = pf * HF_POOL / 1e3
    axC.scatter([pv], [cs], s=70, color=INK, zorder=4, edgecolor="white")
    tx, ty = offsets[pv]
    axC.annotate(f"{lab}: PAF {pf*100:.2f}% ≈ {cs/1e3:.2f} M", (pv, cs),
                 xytext=(tx, ty), fontsize=8.5, color=INK, fontweight="bold",
                 arrowprops=dict(arrowstyle="-", color=SLATE, lw=0.8,
                                 shrinkA=2, shrinkB=3))
axC.axhline(3613, color=TEAL, ls="--", lw=1.4)
axC.text(1.2, 3850, "GBD attribution: 3.61 M",
         fontsize=8, color=TEAL, fontweight="bold")
axC.set_xlabel("Assumed COPD prevalence (%)")
axC.set_ylabel("Attributable HF cases (thousands)")
axC.set_xlim(0.5, 15); axC.set_ylim(0, 4600)
axC.set_title("Population-attributable fraction (MR-Levin)", loc="left")
strip_spines(axC); grid(axC)
panel_tag(axC, "C", CRIMSON)

supertitle(fig, "Figure 5  |  The 2050 horizon: a doubling tide and the fraction we can still prevent",
           "BAPC projection with credible intervals and frozen-population sensitivity (A), evidence-based intervention scenarios (B), and the MR-derived PAF reconciled against GBD attribution (C).")
source_note(fig, "Data: t1_bapc_total_cases.csv · t1_bapc_total_cases_frozenpop.csv · Levin PAF on OR 1.15 with 60 M HF pool · t1_bridge_table.csv")
save(fig, "Fig5_horizon_2050", OUT)
