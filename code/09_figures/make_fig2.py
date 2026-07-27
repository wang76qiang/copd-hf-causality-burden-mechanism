# -*- coding: utf-8 -*-
# ============================================================================
# Module 09 - Figure 2: global trends 1990-2021 and the high-SDI resurgence
# Inputs: results/t1_asr_series.csv, fig_sdi_series_prev.csv,
#         fig_eapc_prev_regions.csv, t1_eapc_covid_sensitivity.csv
# Output: figures/Fig2_divergent_worlds.{png,pdf}
# ============================================================================

import os as _os
# --- repository paths -----------------------------------------------------------
# REPO resolves to the repository root (code/<module>/this_file.py -> ../..).
# Override with the REPRO_ROOT environment variable if you run from elsewhere.
REPO = _os.environ.get("REPRO_ROOT", _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "..")))
# GBD_RAW_DIR: folder with the raw IHME GBD 2021 csv downloads (see data/README.md)
GBD_RAW_DIR = _os.environ.get("GBD_RAW_DIR", _os.path.join(REPO, "data", "raw", "gbd_2021"))
"""Figure 2 — Divergent worlds: global trends and the high-SDI resurgence."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from figstyle import *
import pandas as pd

R = os.path.join(REPO, "results") + os.sep
OUT = os.path.join(REPO, "figures")

asr = pd.read_csv(R + "t1_asr_series.csv")
sdi = pd.read_csv(R + "fig_sdi_series_prev.csv", index_col=0)
eapc = pd.read_csv(R + "fig_eapc_prev_regions.csv")
covid = pd.read_csv(R + "t1_eapc_covid_sensitivity.csv")

fig = plt.figure(figsize=(13.5, 9.2))
gs = fig.add_gridspec(2, 2, hspace=0.55, wspace=0.30,
                      left=0.08, right=0.97, top=0.87, bottom=0.09)

# ---------- A: global trend with COVID band ----------
axA = fig.add_subplot(gs[0, 0])
axA.plot(asr.year, asr.ASPR_prev, color=TEAL, lw=2.4, label="ASPR (prevalence)")
axA.plot(asr.year, asr.ASYR_yld * 10, color=CRIMSON, lw=2.4, ls="-",
         label="ASYR (YLDs) ×10")
axA.axvspan(2020, 2021.9, color=AMBER, alpha=0.18, zorder=0)
axA.annotate("2021 rebound\n41.6 → 43.2", xy=(2020.9, 43.4), xytext=(2013.5, 50.2),
             fontsize=8, color=INK, fontweight="bold",
             arrowprops=dict(arrowstyle="-|>", color=INK, lw=1))
axA.annotate("COVID-19 era", xy=(2020.9, 34.6), fontsize=7.5, color=AMBER,
             fontweight="bold", ha="center")
axA.set_xlim(1989, 2022.5); axA.set_ylim(33, 52)
axA.set_xlabel("Year"); axA.set_ylabel("Age-standardized rate per 100,000")
axA.set_title("Global trajectory, 1990–2021", loc="left")
axA.legend(loc="upper left", fontsize=8)
strip_spines(axA); grid(axA)
panel_tag(axA, "A", TEAL)

# ---------- B: SDI quintile spaghetti ----------
axB = fig.add_subplot(gs[0, 1])
order = ["High SDI", "High-middle SDI", "Middle SDI", "Low-middle SDI", "Low SDI"]
cols = {"High SDI": CRIMSON, "High-middle SDI": AMBER, "Middle SDI": SLATE,
        "Low-middle SDI": SKY, "Low SDI": TEAL}
for loc in order:
    lw = 3.0 if loc == "High SDI" else 1.6
    axB.plot(sdi.index, sdi[loc], color=cols[loc], lw=lw, label=loc, zorder=3 if loc=="High SDI" else 2)
axB.annotate("High SDI · EAPC +1.33",
             xy=(2020.2, sdi["High SDI"].iloc[-2] + 0.4), xytext=(1990.5, 25.5),
             fontsize=8, color=CRIMSON, fontweight="bold",
             arrowprops=dict(arrowstyle="-|>", color=CRIMSON, lw=1.1,
                             connectionstyle="arc3,rad=-0.2"))
axB.set_xlim(1989, 2022); axB.set_xlabel("Year")
axB.set_ylabel("ASPR per 100,000")
axB.set_title("Five development strata diverge", loc="left")
axB.legend(loc="upper right", fontsize=7.5, ncol=1)
strip_spines(axB); grid(axB)
panel_tag(axB, "B", CRIMSON)

# ---------- C: regional EAPC diverging bars ----------
axC = fig.add_subplot(gs[1, 0])
e = eapc[~eapc.Region.isin(["Global", "Female", "Male"])].copy()
e = e.sort_values("Median")
y = np.arange(len(e))
colors = [CRIMSON if v > 0 else TEAL for v in e.Median]
axC.barh(y, e.Median, color=colors, height=0.66, zorder=3)
axC.errorbar(e.Median, y, xerr=[e.Median - e.Min, e.Max - e.Median],
             fmt="none", ecolor=INK, elinewidth=0.7, capsize=1.5, zorder=4)
axC.axvline(0, color=INK, lw=0.9)
axC.set_yticks(y); axC.set_yticklabels(e.Region, fontsize=7.5)
axC.set_xlabel("EAPC of ASPR, % per year (95% CI)")
axC.set_title("Regional EAPC: resurgence vs decline", loc="left")

strip_spines(axC); grid(axC, "x")
panel_tag(axC, "C", CRIMSON)

# ---------- D: COVID sensitivity ----------
axD = fig.add_subplot(gs[1, 1])
labels = {"Prevalence ASR": "ASPR", "YLD ASR": "ASYR"}
ypos = [3.2, 2.2, 1.2, 0.2]
for (i, row), yv in zip(covid.iterrows(), ypos):
    period = row["period"]; meas = labels[row["measure"]]
    c = INK if period == "1990-2021" else AMBER
    m = "s" if period == "1990-2021" else "o"
    axD.errorbar(row["EAPC_pct"], yv, xerr=[[row["EAPC_pct"] - row["CI_low"]], [row["CI_high"] - row["EAPC_pct"]]],
                 fmt=m, color=c, ecolor=c, elinewidth=2 if period == "1990-2021" else 1.2,
                 markersize=8 if period == "1990-2021" else 6, capsize=4, zorder=3)
    axD.text(-0.42, yv, f"{meas} · {period}", ha="right", va="center", fontsize=8.5,
             color=c, fontweight="bold" if period == "1990-2021" else "normal")
axD.axvline(0, color=SLATE, ls="--", lw=0.9)
axD.set_xlim(-0.44, 0.02); axD.set_ylim(-0.5, 4.0)
axD.set_yticks([])
axD.set_xlabel("EAPC, % per year (95% CI)")
axD.set_title("COVID-era sensitivity: trend is robust", loc="left")

strip_spines(axD, ("bottom",)); grid(axD, "x")
panel_tag(axD, "D", AMBER)

supertitle(fig, "Figure 2  |  Divergent worlds: decline almost everywhere, resurgence where least expected",
           "Age-standardized COPD-attributable HF rates, 1990–2021: global trends, SDI-stratified trajectories, regional EAPC, and COVID-era sensitivity.")
source_note(fig, "Data: t1_asr_series.csv · fig_sdi_series_prev.csv (raw IHME GBD 2021) · regional EAPC workbook · t1_eapc_covid_sensitivity.csv")
save(fig, "Fig2_divergent_worlds", OUT)
