# -*- coding: utf-8 -*-
"""Figure 3 — The 26-fold divide: geography of COPD-attributable HF."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from figstyle import *
import pandas as pd
import numpy as np

R = r"E:\COPD+HF\T1_TopJournal\results\\"
OUT = r"E:\COPD+HF\T1_TopJournal\figures\final"

c = pd.read_csv(R + "fig_country_1990_2021.csv")

fig = plt.figure(figsize=(13.5, 9.6))
gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1], hspace=0.5, wspace=0.28,
                      left=0.12, right=0.97, top=0.87, bottom=0.09)

# ---------- A: dumbbell top-12 + bottom-4 ----------
axA = fig.add_subplot(gs[0, :])
top = c.nlargest(12, "ASPR_2021")
bot = c.nsmallest(4, "ASPR_2021")
d = pd.concat([top, bot]).iloc[::-1].reset_index(drop=True)
y = np.arange(len(d))
for i, row in d.iterrows():
    rise = row.ASPR_2021 >= row.ASPR_1990
    col = CRIMSON if rise else TEAL
    axA.plot([row.ASPR_1990, row.ASPR_2021], [i, i], color=col, lw=2.2, zorder=2, alpha=0.75)
    axA.scatter([row.ASPR_1990], [i], s=55, color="white", edgecolor=col, lw=1.8, zorder=3)
    axA.scatter([row.ASPR_2021], [i], s=90, color=col, edgecolor="white", lw=1.2, zorder=4)
    axA.text(max(row.ASPR_1990, row.ASPR_2021) + 2.0, i,
             f"{row.ASPR_2021:.1f}", va="center", fontsize=8, color=col, fontweight="bold")
axA.axvline(43.24, color=SLATE, ls=":", lw=1.1)
axA.text(43.7, -0.75, "global mean 43.2", fontsize=7.5, color=SLATE, va="top")
axA.set_yticks(y); axA.set_yticklabels(d.location_name, fontsize=8.5)
axA.set_xlabel("Age-standardized prevalence per 100,000  (○ 1990   ● 2021)")
axA.set_title("Dumbbell of extremes: 1990 → 2021 (top-12 and bottom-4 countries, 2021)", loc="left")
axA.set_xlim(0, 125)
strip_spines(axA); grid(axA, "x")
panel_tag(axA, "A", INK, x=-0.075)

# ---------- B: all-country ranked strip ----------
axB = fig.add_subplot(gs[1, 0])
s = c.sort_values("ASPR_2021").reset_index(drop=True)
x = np.arange(len(s))
axB.scatter(x, s.ASPR_2021, s=7, color=SLATE, alpha=0.55, zorder=2)
hi = s.nlargest(5, "ASPR_2021"); lo = s.nsmallest(5, "ASPR_2021")
axB.scatter(hi.index, hi.ASPR_2021, s=42, color=CRIMSON, zorder=4)
axB.scatter(lo.index, lo.ASPR_2021, s=42, color=TEAL, zorder=4)
for k, (_, r) in enumerate(hi.iterrows()):
    axB.annotate(r.location_name.replace("Papua New Guinea", "PNG"),
                 (r.name, r.ASPR_2021), textcoords="offset points",
                 xytext=(-18 + 8 * (k % 3), 7 + 9 * (k % 2)),
                 fontsize=7, ha="center", color=CRIMSON, fontweight="bold")
for k, (_, r) in enumerate(lo.iterrows()):
    axB.annotate(r.location_name, (r.name, r.ASPR_2021),
                 textcoords="offset points", xytext=(10 + 6 * (k % 2), 5 + 9 * k),
                 fontsize=7, ha="left", color=TEAL, fontweight="bold")

axB.set_xlabel("204 countries, ranked"); axB.set_ylabel("ASPR 2021 per 100,000")
axB.set_title("The full gradient", loc="left")
strip_spines(axB); grid(axB)
panel_tag(axB, "B", CRIMSON, x=-0.13)

# ---------- C: cases vs rate, China spotlight ----------
axC = fig.add_subplot(gs[1, 1])
axC.scatter(c.ASPR_2021, c.cases_2021 / 1e3, s=np.sqrt(c.cases_2021) / 12 + 6,
            color=SKY, alpha=0.45, edgecolor="white", lw=0.3, zorder=2)
cn = c[c.location_name == "China"].iloc[0]
axC.scatter([cn.ASPR_2021], [cn.cases_2021 / 1e3], s=220, color=CRIMSON,
            edgecolor="white", lw=1.4, zorder=4)
axC.annotate("China\nASPR 76.9 · 1.51 M cases\n(#2 rate · #1 caseload)",
             (cn.ASPR_2021, cn.cases_2021 / 1e3), xytext=(38, 1180),
             fontsize=8.5, fontweight="bold", color=CRIMSON,
             arrowprops=dict(arrowstyle="-|>", color=CRIMSON, lw=1.2,
                             connectionstyle="arc3,rad=0.15"))
for name, dx, dy in [("India", 10, -28), ("United States of America", 6, 20)]:
    r = c[c.location_name == name]
    if len(r):
        r = r.iloc[0]
        axC.annotate(name.replace("United States of America", "USA"),
                     (r.ASPR_2021, r.cases_2021 / 1e3), textcoords="offset points",
                     xytext=(dx, dy), fontsize=7.5, color=INK,
                     arrowprops=dict(arrowstyle="-", color=SLATE, lw=0.7))
axC.set_xlabel("ASPR 2021 per 100,000"); axC.set_ylabel("Cases 2021 (thousands)")
axC.set_title("Rate × volume: where the people are", loc="left")
strip_spines(axC); grid(axC)
panel_tag(axC, "C", CRIMSON, x=-0.13)

supertitle(fig, "Figure 3  |  The 26-fold divide: geography of COPD-attributable heart failure in 2021",
           "From Papua New Guinea (83.1/100k) to Uzbekistan (3.2/100k): national rates, trajectories since 1990, and the rate–volume duality that makes China the pivotal battleground.")
source_note(fig, "Data: fig_country_1990_2021.csv (raw IHME GBD 2021, 204 countries, both sexes)")
save(fig, "Fig3_26fold_divide", OUT)
