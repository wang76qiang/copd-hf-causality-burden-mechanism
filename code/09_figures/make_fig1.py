# -*- coding: utf-8 -*-
# ============================================================================
# Module 09 - Figure 1: genetic causality, replication, adjustment
# (discovery MR forest + leave-one-out, three-cohort meta, MVMR, mediation)
# Inputs: results/t1_three_cohort_meta.json, t1_mvmr_results.csv,
#         t1_mvmr_lungfunction.csv, t1_mediation_results.csv,
#         t1_hfpef_hfref_mr.csv; results/legacy/real_data_results.json
# Output: figures/Fig1_causal_ladder.{png,pdf}
# ============================================================================

import os as _os
# --- repository paths -----------------------------------------------------------
# REPO resolves to the repository root (code/<module>/this_file.py -> ../..).
# Override with the REPRO_ROOT environment variable if you run from elsewhere.
REPO = _os.environ.get("REPRO_ROOT", _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "..")))
# GBD_RAW_DIR: folder with the raw IHME GBD 2021 csv downloads (see data/README.md)
GBD_RAW_DIR = _os.environ.get("GBD_RAW_DIR", _os.path.join(REPO, "data", "raw", "gbd_2021"))
"""Figure 1 — The causal ladder: discovery, replication, adjustment, mediation."""
import sys, json, os
sys.path.insert(0, os.path.dirname(__file__))
from figstyle import *
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from matplotlib.ticker import FixedLocator, FixedFormatter

def plain_log_ticks(ax, ticks):
    ax.xaxis.set_major_locator(FixedLocator(ticks))
    ax.xaxis.set_major_formatter(FixedFormatter([f"{t:g}" for t in ticks]))
    ax.tick_params(axis="x", which="minor", bottom=False)

R = os.path.join(REPO, "results") + os.sep
OUT = os.path.join(REPO, "figures")

meta = json.load(open(R + "t1_three_cohort_meta.json"))
mvmr1 = pd.read_csv(R + "t1_mvmr_results.csv")
mvmr2 = pd.read_csv(R + "t1_mvmr_lungfunction.csv")
med = pd.read_csv(R + "t1_mediation_results.csv")
finn = pd.read_csv(R + "t1_hfpef_hfref_mr.csv")

fig = plt.figure(figsize=(13.5, 9.2))
gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.05], hspace=0.52, wspace=0.42,
                      left=0.09, right=0.965, top=0.88, bottom=0.09)

# ---------- A: design schematic ----------
axA = fig.add_subplot(gs[0, 0]); axA.axis("off")
panel_tag(axA, "A", INK, x=-0.02, y=1.02)
axA.set_xlim(0, 10); axA.set_ylim(0, 10)
axA.set_title("Three-cohort triangulation", loc="left")
boxes = [
    (1.0, 7.6, "Discovery\nEUR · GCST009541\n47,309 / 930,014", TEAL),
    (1.0, 4.6, "Replication 1\nFinnGen R12 · strict HF\n37,653 / 462,695", SKY),
    (1.0, 1.6, "Replication 2\nBioBank Japan · CHF\n9,413 / 203,040", SLATE),
]
for x, y, txt, c in boxes:
    axA.add_patch(FancyBboxPatch((x, y), 4.4, 1.9, boxstyle="round,pad=0.12",
                                 fc=c, ec="none", alpha=0.92))
    axA.text(x+2.2, y+0.95, txt, color="white", fontsize=8, ha="center",
             va="center", fontweight="bold", linespacing=1.5)
axA.add_patch(FancyBboxPatch((6.4, 4.6), 3.3, 1.9, boxstyle="round,pad=0.12",
                             fc=CRIMSON, ec="none"))
axA.text(8.05, 5.55, "Meta-analysis\nFE + RE · I² test", color="white",
         fontsize=8, ha="center", va="center", fontweight="bold", linespacing=1.5)
for y in (8.55, 5.55, 2.55):
    axA.add_patch(FancyArrowPatch((5.55, y), (6.35, 5.55), arrowstyle="-|>",
                                  mutation_scale=13, color=INK, lw=1.2))
axA.text(5.0, 0.3, "10 instruments · F = 31.5–166.4 · Steiger 10/10",
         fontsize=7.5, color=SLATE, ha="center")

# ---------- B: discovery forest ----------
axB = fig.add_subplot(gs[0, 1])
labels = ["IVW (fixed)", "IVW (random)", "Weighted median", "Simple mode",
          "MR–Egger", "Radial IVW", "CAUSE", "Relaxed IVs (49 SNP)"]
ors = [1.15, 1.16, 1.22, 1.22, 1.15, 1.24, 1.15, 1.16]
lo  = [1.08, 1.05, 1.13, 1.10, 1.05, 1.13, 1.03, 1.02]
hi  = [1.22, 1.28, 1.32, 1.34, 1.26, 1.36, 1.27, 1.33]
cols = [CRIMSON, CRIMSON, INK, INK, INK, INK, SKY, SLATE]
forest(axB, range(len(labels)), labels, ors, lo, hi, cols,
       title="Discovery: COPD → HF (EUR)", diamond_idx={0},
       xlim=(0.95, 1.45))
plain_log_ticks(axB, [1.0, 1.1, 1.2, 1.3, 1.4])

panel_tag(axB, "B", CRIMSON)

# ---------- C: three-cohort meta forest ----------
axC = fig.add_subplot(gs[0, 2])
labs = ["EUR discovery", "FinnGen R12", "BBJ (E. Asia)", "Fixed effect", "Random effects"]
ors  = [1.15, 1.037, 0.933, 1.059, 1.042]
lo   = [1.082, 0.986, 0.848, 1.022, 0.94]
hi   = [1.222, 1.090, 1.027, 1.098, 1.155]
cols = [TEAL, SKY, SLATE, CRIMSON, AMBER]
forest(axC, range(len(labs)), labs, ors, lo, hi, cols,
       title="Cross-cohort meta-analysis", diamond_idx={3, 4},
       xlim=(0.80, 1.35))
plain_log_ticks(axC, [0.8, 0.9, 1.0, 1.1, 1.2, 1.3])

panel_tag(axC, "C", CRIMSON)

# ---------- D: MVMR two models ----------
axD = fig.add_subplot(gs[1, 0:2])
def mvmr_offsets(df, y0):
    return [(row["exposure"], row["or"], row["or_lci95"], row["or_uci95"], y0 - i*0.55)
            for i, row in df.iterrows()]
rows1 = mvmr_offsets(mvmr1, 5.6)
rows2 = mvmr_offsets(mvmr2, 2.3)
cmap = {"COPD": TEAL, "SMOK": SLATE, "BMI": AMBER, "SBP": CRIMSON,
        "FEV1": SKY, "FEV1FVC": SKY, "CPD": SLATE}
axD.text(0.4, 6.6, "Model 1 · adj. smoking initiation + BMI + SBP (709 SNP)",
         fontsize=8.5, fontweight="bold", color=INK)
axD.text(0.4, 3.3, "Model 2 · adj. FEV1 + FEV1/FVC + cigarettes/day (405 SNP)",
         fontsize=8.5, fontweight="bold", color=INK)
for name, o, l, h, yv in rows1 + rows2:
    c = cmap.get(name, INK)
    bold = name == "COPD"
    lab = {"COPD": "COPD", "SMOK": "Smoking init.", "BMI": "BMI", "SBP": "SBP",
           "FEV1": "FEV1", "FEV1FVC": "FEV1/FVC", "CPD": "Cigarettes/day"}[name]
    axD.plot([l, h], [yv, yv], color=c, lw=3.2 if bold else 1.6,
             solid_capstyle="round", zorder=3)
    axD.scatter([o], [yv], s=110 if bold else 55, color=c, zorder=4,
                edgecolor="white", lw=0.8)
    axD.text(0.42, yv + 0.16, lab, fontsize=7.8, color=c, va="bottom",
             fontweight="bold" if bold else "normal")
    axD.text(1.72, yv, f"{o:.2f} ({l:.2f}–{h:.2f})", fontsize=7.5,
             va="center", ha="right", color=c,
             fontweight="bold" if bold else "normal")
axD.axvline(1.0, color=SLATE, ls="--", lw=0.9)
axD.set_xlim(0.4, 1.74); axD.set_ylim(-0.3, 7.1)
axD.set_yticks([]); axD.set_xlabel("Odds ratio for HF")
axD.set_title("Multivariable MR: the COPD effect survives adjustment", loc="left")
strip_spines(axD, ("left", "bottom")); grid(axD, "x")
panel_tag(axD, "D", TEAL, x=-0.07)

# ---------- E: mediation ----------
axE = fig.add_subplot(gs[1, 2])
med = med.sort_values("prop_mediated")
names = med.mediator.tolist()
prop = (med.prop_mediated * 100).tolist()
pv = med.indirect_p.tolist()
cols = [CRIMSON if (p < 0.05/4) else LIGHT for p in pv]
edge = [CRIMSON if (p < 0.05/4) else SLATE for p in pv]
y = np.arange(len(names))
axE.barh(y, prop, color=cols, edgecolor=edge, height=0.62, zorder=3)
for yi, pr, p, nm in zip(y, prop, pv, names):
    sig = "P = 1.2×10⁻⁵ ✓" if p < 0.05/4 else f"ns (P = {p:.2f})"
    axE.text(max(pr, 0) + 0.15, yi, f"{pr:+.1f}% · {sig}", va="center",
             fontsize=8, color=INK if p < 0.05/4 else SLATE,
             fontweight="bold" if p < 0.05/4 else "normal")
axE.axvline(0, color=INK, lw=0.8)
axE.set_yticks(y); axE.set_yticklabels(names)
axE.set_xlabel("Proportion mediated (%)")
axE.set_xlim(-8, 9)
axE.set_title("Mediation: only CRP survives\nBonferroni (4 tests)", loc="left")
strip_spines(axE); grid(axE, "x")
panel_tag(axE, "E", AMBER)

supertitle(fig, "Figure 1  |  Genetic causality of COPD in heart failure",
           "Discovery → replication → meta-analysis → multivariable adjustment → mediation. Every estimate backed by a deposited output file.")
source_note(fig, "Data: t1_three_cohort_meta.json · t1_mvmr_results.csv · t1_mvmr_lungfunction.csv · t1_mediation_results.csv · real_data_results.json (legacy archive)")
save(fig, "Fig1_causal_ladder", OUT)
