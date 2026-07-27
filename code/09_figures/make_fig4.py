# -*- coding: utf-8 -*-
# ============================================================================
# Module 09 - Figure 4: engines of growth (decomposition, inequality, age)
# NOTE: panels A-C contain hard-coded decomposition and SII/CI values that
# are computed by code/05_gbd/gbd_decomposition.py and gbd_inequality.py;
# see VALIDATION.md for the numerical verification (and one documented
# sign discrepancy in the high-SDI panel B).
# Input : results/t1_bapc_agepanel.csv
# Output: figures/Fig4_engines_of_growth.{png,pdf}
# ============================================================================

import os as _os
# --- repository paths -----------------------------------------------------------
# REPO resolves to the repository root (code/<module>/this_file.py -> ../..).
# Override with the REPRO_ROOT environment variable if you run from elsewhere.
REPO = _os.environ.get("REPRO_ROOT", _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "..")))
# GBD_RAW_DIR: folder with the raw IHME GBD 2021 csv downloads (see data/README.md)
GBD_RAW_DIR = _os.environ.get("GBD_RAW_DIR", _os.path.join(REPO, "data", "raw", "gbd_2021"))
"""Figure 4 — Engines of growth: decomposition, inequality, age structure."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from figstyle import *
import pandas as pd
import numpy as np

R = os.path.join(REPO, "results") + os.sep
OUT = os.path.join(REPO, "figures")

age = pd.read_csv(R + "t1_bapc_agepanel.csv")

fig = plt.figure(figsize=(13.5, 9.4))
gs = fig.add_gridspec(2, 2, hspace=0.52, wspace=0.30,
                      left=0.09, right=0.96, top=0.87, bottom=0.09)

# ---------- A: decomposition waterfall (global) ----------
axA = fig.add_subplot(gs[0, 0])
drivers = ["Population\nageing", "Population\ngrowth", "Epidemiological\nchange"]
vals = [53.32, 44.74, 1.94]
cols = [CRIMSON, AMBER, TEAL]
cum = np.concatenate([[0], np.cumsum(vals)])
for i, (d, v, c) in enumerate(zip(drivers, vals, cols)):
    axA.bar(i, v, bottom=cum[i], color=c, width=0.62, zorder=3)
    if i == 2:  # epidemiological change: black label outside the small bar
        axA.text(i, cum[i] + v + 1.6, f"+{v:.1f}%", ha="center", va="bottom",
                 color="black", fontweight="bold", fontsize=10)
    else:
        axA.text(i, cum[i] + v / 2, f"+{v:.1f}%", ha="center", va="center",
                 color="white", fontweight="bold", fontsize=10)
axA.bar(3, cum[-1], color=INK, width=0.62, zorder=3)
axA.text(3, cum[-1] / 2, "100%", ha="center", va="center", color="white",
         fontweight="bold", fontsize=10)
axA.set_xticks(range(4)); axA.set_xticklabels(drivers + ["TOTAL"], fontsize=8.5)
axA.set_ylabel("Contribution to case growth (%)")
axA.set_ylim(0, 112)
axA.set_title("Global decomposition: demography rules", loc="left")
strip_spines(axA); grid(axA)
panel_tag(axA, "A", CRIMSON)

# ---------- B: decomposition by SDI ----------
axB = fig.add_subplot(gs[0, 1])
groups = ["Low SDI", "Middle SDI", "High SDI"]
# validated prevalence-based decomposition (gbd_decomposition.py, 2026-07-27):
aging = [8.60, 85.94, 44.82]
growth = [108.87, 41.84, 20.72]
epi = [-17.47, -27.79, 34.46]
x = np.arange(3); w = 0.26
for k, (arr, lab, c) in enumerate([(aging, "Ageing", CRIMSON),
                                   (growth, "Population growth", AMBER),
                                   (epi, "Epidemiological change", TEAL)]):
    b = axB.bar(x + (k - 1) * w, arr, width=w, color=c, label=lab, zorder=3)
    for xi, v in zip(x + (k - 1) * w, arr):
        axB.text(xi, v + (5 if v >= 0 else -9), f"{v:+.0f}", ha="center",
                 fontsize=7.5, color=c, fontweight="bold")
axB.axhline(0, color=INK, lw=0.9)
axB.set_xticks(x); axB.set_xticklabels(groups, fontsize=8.5)
axB.set_ylabel("Contribution (%)")
axB.set_ylim(-75, 135)
axB.set_title("Drivers invert across development", loc="left")
axB.legend(fontsize=7.5, loc="upper right", ncol=1)
strip_spines(axB); grid(axB)
panel_tag(axB, "B", AMBER)

# ---------- C: inequality slope chart ----------
axC = fig.add_subplot(gs[1, 0])
metrics = [("SII, ASPR", -40.05, -23.10), ("SII, ASYR", -3.52, -2.03)]
for i, (lab, v90, v21) in enumerate(metrics):
    axC.plot([0, 1], [v90, v21], color=CRIMSON if i == 0 else SKY, lw=2.4,
             marker="o", markersize=9, zorder=3)
    axC.text(-0.06, v90, f"{lab}  {v90:.2f}", ha="right", va="center", fontsize=8.5)
    axC.text(1.06, v21, f"{v21:.2f}", ha="left", va="center", fontsize=8.5,
             fontweight="bold", color=CRIMSON if i == 0 else SKY)
axC.set_xlim(-0.85, 1.5); axC.set_ylim(-46, 3)
axC.set_xticks([0, 1]); axC.set_xticklabels(["1990", "2021"])
axC.set_ylabel("Slope index of inequality")
axC.set_title("Inequality narrows but stays pro-poor\n(CI: −0.25 → −0.15)", loc="left")
strip_spines(axC); grid(axC)
panel_tag(axC, "C", SKY)

# ---------- D: age × year heatmap ----------
axD = fig.add_subplot(gs[1, 1])
casecols = [col for col in age.columns if col.startswith("cases_")]
agelab = [col.replace("cases_", "") for col in casecols]
M = age[casecols].to_numpy().T / 1000.0  # thousands
years = age.year.to_numpy()
keep = [2, 5, 8, 11, 14, 16, 18]  # 10-14 ... 85-89
M2 = M[keep]; labs = [agelab[i] for i in keep]
im = axD.imshow(M2, aspect="auto", cmap="inferno",
                extent=[years.min(), years.max(), -0.5, len(keep) - 0.5],
                origin="lower")
axD.set_yticks(range(len(keep))); axD.set_yticklabels(labs, fontsize=8)
cb = fig.colorbar(im, ax=axD, fraction=0.035, pad=0.02)
cb.set_label("Cases (thousands)", fontsize=8)
axD.set_xlabel("Year")
axD.set_title("Observed age structure of burden, 1990–2021\n(the ≥70 wave that BAPC extrapolates)", loc="left")
panel_tag(axD, "D", INK)

supertitle(fig, "Figure 4  |  Engines of growth: demography, inequality, and the ageing wave",
           "Decomposition of case growth 1990–2021 (A, B), slope-index inequality (C), and the age-resolved BAPC projection revealing the ≥70-year surge (D).")
source_note(fig, "Data: decomposition & SII/CI estimates (verified manuscript values) · t1_bapc_agepanel.csv (observed 1990-2021)")
save(fig, "Fig4_engines_of_growth", OUT)
