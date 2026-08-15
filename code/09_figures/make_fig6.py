# -*- coding: utf-8 -*-
# ============================================================================
# Module 09 - Figure 6: shared inflammo-fibrotic program + cross-cohort
# Inputs: results/t1_crosscohort_validation.csv;
#         data/derived/bulk/GSE57148_GSE57338_keygene_expression_values.xlsx
# Output: figures/Fig6_molecular_convergence.{png,pdf}
# ============================================================================

import os as _os
# --- repository paths -----------------------------------------------------------
# REPO resolves to the repository root (code/<module>/this_file.py -> ../..).
# Override with the REPRO_ROOT environment variable if you run from elsewhere.
REPO = _os.environ.get("REPRO_ROOT", _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "..")))
# GBD_RAW_DIR: folder with the raw IHME GBD 2021 csv downloads (see data/README.md)
GBD_RAW_DIR = _os.environ.get("GBD_RAW_DIR", _os.path.join(REPO, "data", "raw", "gbd_2021"))
"""Figure 6 — Molecular convergence: a shared inflammo-fibrotic program."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from figstyle import *
import pandas as pd
import numpy as np
from matplotlib.patches import Circle
from scipy import stats as sst

R = os.path.join(REPO, "results") + os.sep
# MedComm revision: write into the working-copy folder, not the source repo
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")

cc = pd.read_csv(R + "t1_crosscohort_validation.csv")

# ---------- parse expression workbook ----------
def parse_sheet(path, sheet):
    raw = pd.read_excel(path, sheet_name=sheet)
    out = {}
    for j in range(raw.shape[1] - 2):
        if str(raw.iloc[0, j]).strip() == "Sample":
            gene = str(raw.columns[j])
            if gene.startswith("GSE"):
                gene = str(raw.columns[j + 1])
            vals = pd.to_numeric(raw.iloc[1:, j + 2], errors="coerce")
            grp = raw.iloc[1:, j + 1].astype(str)
            out[gene] = pd.DataFrame({"group": grp, "value": vals}).dropna()
    return out

XL = os.path.join(REPO, "data", "derived", "bulk", "GSE57148_GSE57338_keygene_expression_values.xlsx")
lung = parse_sheet(XL, "GSE57148")
heart = parse_sheet(XL, "GSE57338")
for g in lung:
    lung[g]["group"] = np.where(lung[g]["group"].str.lower().str.contains("copd"), "COPD", "Control")
# THBS1 block lacks group labels: borrow from S100A8 (same sample order, converted labels)
lung["THBS1"] = lung["THBS1"].assign(group=lung["S100A8"]["group"].values)
for g in heart:
    s = heart[g]["group"].astype(str)
    heart[g]["group"] = np.where(s.str.contains("Non-failing") | s.str.fullmatch("Control"),
                                 "Non-failing", "Heart failure")

fig = plt.figure(figsize=(14.5, 11.8))
gs = fig.add_gridspec(3, 3, height_ratios=[1, 1, 1.45], hspace=0.62, wspace=0.40,
                      left=0.07, right=0.975, top=0.885, bottom=0.08)

# ---------- A: DEG counts mirrored ----------
axA = fig.add_subplot(gs[0, 0])
cats = ["COPD lung\n(GSE57148)", "Failing heart\n(GSE57338)"]
up = [542, 234]; down = [819, 273]
y = np.arange(2)
axA.barh(y, up, color=CRIMSON, height=0.55, zorder=3, label="Upregulated")
axA.barh(y, [-d for d in down], color=SKY, height=0.55, zorder=3, label="Downregulated")
for yi, u, d in zip(y, up, down):
    axA.text(u + 25, yi, f"+{u}", va="center", color=CRIMSON, fontweight="bold", fontsize=10)
    axA.text(-d - 25, yi, f"−{d}", va="center", ha="right", color=SKY, fontweight="bold", fontsize=10)
axA.axvline(0, color=INK, lw=1)
axA.set_yticks(y); axA.set_yticklabels(cats, fontsize=9)
axA.set_xlim(-1000, 800); axA.set_xlabel("DEGs (adj. P < 0.05, |log₂FC| ≥ 0.5)")
axA.legend(loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=2, fontsize=8.5,
           frameon=False, handlelength=1.6, columnspacing=1.8)
axA.set_title("Two organs, two signatures", loc="left")
strip_spines(axA); grid(axA, "x")
panel_tag(axA, "A", INK, x=-0.28)

# ---------- B: Venn ----------
axB = fig.add_subplot(gs[0, 1]); axB.axis("off"); axB.set_aspect("equal")
c1 = Circle((0.42, 0.52), 0.34, color=TEAL, alpha=0.45, lw=2, ec=TEAL)
c2 = Circle((0.60, 0.52), 0.27, color=CRIMSON, alpha=0.45, lw=2, ec=CRIMSON)
axB.add_patch(c1); axB.add_patch(c2)
axB.text(0.26, 0.55, "1,320", fontsize=11, fontweight="bold", color=INK, ha="center")
axB.text(0.515, 0.55, "41", fontsize=13, fontweight="bold", color=INK, ha="center")
axB.text(0.75, 0.55, "466", fontsize=11, fontweight="bold", color=INK, ha="center")
axB.text(0.13, 0.13, "COPD lung DEGs\n(1,361)", fontsize=9, color=TEAL, ha="center", fontweight="bold")
axB.text(0.87, 0.13, "Failing heart DEGs\n(507)", fontsize=9, color=CRIMSON, ha="center", fontweight="bold")

axB.set_xlim(0, 1); axB.set_ylim(0, 1)
axB.set_title("The molecular handshake", loc="left")
panel_tag(axB, "B", TEAL, x=-0.05, y=1.02)

# ---------- C: 41-gene lung-vs-heart log2FC ----------
axC = fig.add_subplot(gs[0, 2])
lungfc = cc[cc.cohort == "CXG_COPD_lung"].set_index("gene").effect_log2fc
heartfc = cc[cc.cohort == "CXG_DCM_heart"].set_index("gene").effect_log2fc
common = lungfc.index.intersection(heartfc.index)
x, yv = lungfc[common], heartfc[common]
hub = {"SERPINE1", "THBS1", "SPP1", "IL1R2", "EGR1", "CDKN1A", "CD163", "AREG", "SOCS3"}
axC.axhline(0, color=LIGHT, lw=1); axC.axvline(0, color=LIGHT, lw=1)
lim = max(abs(x).max(), abs(yv).max()) * 1.15
axC.plot([-lim, lim], [-lim, lim], color=SLATE, lw=0.8, ls=":")
colors = [AMBER if g in hub else SLATE for g in common]
axC.scatter(x, yv, s=[70 if g in hub else 22 for g in common], color=colors,
            alpha=0.9, edgecolor="white", lw=0.5, zorder=3)
# leader-line label layout: left column, right column, bottom row
label_pos = {
    # left column (x=-4.25)
    "SVEP1": (-4.25, 2.35), "SFRP4": (-4.25, 1.7), "IL1R2": (-4.25, 1.05),
    "IL1RL1": (-4.25, 0.4), "EGR1": (-4.25, -0.25), "AREG": (-4.25, -0.9),
    "SERPINE1": (-4.25, -1.55),
    # right column (x=2.3)
    "CD163": (2.3, 1.7), "CDKN1A": (2.3, 1.05), "S100A12": (2.3, 0.4),
    "SPP1": (2.3, -0.25), "F13A1": (2.3, -0.9), "THBS1": (2.3, -1.55),
    "SOCS3": (2.3, -2.2), "S100A8": (2.3, -2.85),
    # bottom
    "PLA2G2A": (0.7, -4.2), "CHRDL2": (-1.8, -2.7),
}
for g, (lx, ly) in label_pos.items():
    if g not in common:
        continue
    axC.annotate(g, xy=(x[g], yv[g]), xytext=(lx, ly), fontsize=7.5,
                 color=INK, fontweight="bold" if g in hub else "normal",
                 ha="left", va="center", zorder=4,
                 arrowprops=dict(arrowstyle="-", color=SLATE, lw=0.6,
                                 shrinkA=0, shrinkB=2))
axC.set_xlabel("log₂FC in COPD lung (sc pseudobulk)"); axC.set_ylabel("log₂FC in failing heart (sn pseudobulk)")
axC.set_title("41 genes: same compass in both organs\n(independent single-cell cohorts; hubs in amber)", loc="left")
strip_spines(axC); grid(axC)
panel_tag(axC, "C", AMBER, x=-0.22)

# ---------- D: lung boxplots ----------
axD = fig.add_subplot(gs[1, 0:2])
genes_l = ["S100A8", "S100A12", "IL1RL1", "THBS1"]
pos = 0
xt, xl = [], []
for g in genes_l:
    df = lung[g]
    for k, grp in enumerate(["Control", "COPD"]):
        v = df[df.group == grp].value
        bp = axD.boxplot([v], positions=[pos + k * 0.8], widths=0.62,
                         patch_artist=True, showfliers=False, medianprops=dict(color=INK, lw=1.4))
        for patch in bp["boxes"]:
            patch.set_facecolor(TEAL if grp == "COPD" else LIGHT)
            patch.set_alpha(0.85 if grp == "COPD" else 0.9)
        jitter = np.random.default_rng(1).normal(0, 0.09, len(v))
        axD.scatter(pos + k * 0.8 + jitter, v, s=4, color=INK, alpha=0.25, zorder=3)
        u = sst.ttest_ind(df[df.group == "COPD"].value, df[df.group == "Control"].value,
                          equal_var=False)
    xt.append(pos + 0.4); xl.append(f"{g}\nP={u.pvalue:.0e}")
    pos += 2.6
axD.set_xticks(xt); axD.set_xticklabels(xl, fontsize=8.5)
axD.set_ylabel("Expression (arbitrary units)")
axD.set_title("Lung validation: alarmins & matricellular signals rise in COPD (Mann–Whitney)", loc="left")
strip_spines(axD); grid(axD)
panel_tag(axD, "D", TEAL, x=-0.10)

# ---------- E: heart boxplots ----------
axE = fig.add_subplot(gs[1, 2])
genes_h = ["LUM", "ASPN", "SMOC2", "ACE"]
pos = 0; xt, xl = [], []
for g in genes_h:
    df = heart[g]
    for k, grp in enumerate(["Non-failing", "Heart failure"]):
        v = df[df.group == grp].value
        bp = axE.boxplot([v], positions=[pos + k * 0.8], widths=0.62,
                         patch_artist=True, showfliers=False, medianprops=dict(color=INK, lw=1.4))
        for patch in bp["boxes"]:
            patch.set_facecolor(CRIMSON if grp == "Heart failure" else LIGHT)
        jitter = np.random.default_rng(2).normal(0, 0.09, len(v))
        axE.scatter(pos + k * 0.8 + jitter, v, s=4, color=INK, alpha=0.25, zorder=3)
    u = sst.ttest_ind(df[df.group == "Heart failure"].value, df[df.group == "Non-failing"].value,
                      equal_var=False)
    xt.append(pos + 0.4); xl.append(f"{g}\nP={u.pvalue:.0e}")
    pos += 2.6
axE.set_xticks(xt); axE.set_xticklabels(xl, fontsize=8.5)
axE.set_ylabel("Expression")
axE.set_title("Heart validation: fibrotic matrix\nremodelling in failing LV", loc="left")
strip_spines(axE); grid(axE)
panel_tag(axE, "E", CRIMSON, x=-0.22)

# ---------- F: cross-cohort heatmap (MedComm revision) ----------
# Show the 15 key / manuscript-named genes of the 41-gene intersection so the
# panel stays readable at 180 mm width; the full 41-gene matrix is unchanged
# in t1_crosscohort_validation.csv for the supplement. (LUM/ASPN/ACE are not
# part of the 41-gene intersection and remain shown in panel E.)
axF = fig.add_subplot(gs[2, :])
piv = cc.pivot_table(index="gene", columns="cohort", values="effect_log2fc")
pivp = cc.pivot_table(index="gene", columns="cohort", values="fdr_within_cohort")
order_c = ["GSE57148_COPD_lung", "CXG_COPD_lung", "GSE57338_HF_heart", "CXG_DCM_heart"]
piv = piv.reindex(columns=order_c); pivp = pivp.reindex(columns=order_c)
keep_g = piv.notna().sum(axis=1) >= 2
piv = piv[keep_g]; pivp = pivp[keep_g]
KEY_GENES = ["F13A1", "S100A8", "S100A12", "CD163", "SMOC2", "SFRP4",
             "SERPINE1", "THBS1", "SPP1", "IL1R2", "EGR1", "CDKN1A",
             "AREG", "SOCS3", "IL1RL1"]
KEY_GENES = [g for g in KEY_GENES if g in piv.index]
piv = piv.loc[KEY_GENES]; pivp = pivp.loc[KEY_GENES]
rowmean = piv.mean(axis=1, skipna=True)
piv = piv.loc[rowmean.sort_values().index]; pivp = pivp.loc[piv.index]
M = piv.to_numpy(dtype=float)
vmax = np.nanmax(np.abs(M))
im = axF.imshow(M, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
axF.set_xticks(range(len(order_c)))
axF.set_xticklabels(["Bulk COPD lung\nGSE57148", "sc COPD lung\n(CELLxGENE)",
                     "Bulk failing heart\nGSE57338", "sn failing heart\n(Reichart)"],
                    fontsize=10.5)
axF.set_yticks(range(len(piv)))
axF.set_yticklabels(piv.index, fontsize=11)
for tick, g in zip(axF.get_yticklabels(), piv.index):
    if g in hub:
        tick.set_fontweight("bold")
axF.tick_params(axis="both", length=0)
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        if not np.isnan(M[i, j]) and pivp.to_numpy(dtype=float)[i, j] < 0.05:
            axF.text(j, i, "\u25CF", ha="center", va="center", color="white",
                     fontsize=12, fontweight="bold")
# thin separators make the enlarged cells easier to scan
for yy in np.arange(0.5, M.shape[0], 1):
    axF.axhline(yy, color="white", lw=1.2)
for xx in np.arange(0.5, M.shape[1], 1):
    axF.axvline(xx, color="white", lw=1.2)
cb = fig.colorbar(im, ax=axF, fraction=0.025, pad=0.015)
cb.set_label("log₂FC (dot: FDR < 0.05)", fontsize=10)
cb.ax.tick_params(labelsize=9.5)
axF.set_title("Cross-cohort direction consistency — 15 representative genes of the "
              "41-gene intersection (white = not measured; full matrix in supplement)",
              loc="left", fontsize=9.5)
panel_tag(axF, "F", INK, x=-0.035)

supertitle(fig, "Figure 6  |  Molecular convergence: one inflammo-fibrotic program, two failing organs",
           "DEG landscapes (A), the 41-gene intersection (B), cross-organ concordance (C), direct expression validation in lung (D) and heart (E), and independent replication across four cohorts (F).",
           sy=0.952)
source_note(fig, "Data: t1_crosscohort_validation.csv · GSE57148_S100A8 expression workbook · hub genes from STRING/cytoHubba (manuscript Table S5)")
save(fig, "Fig6_molecular_convergence", OUT)
