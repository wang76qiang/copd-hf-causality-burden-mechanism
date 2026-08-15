# -*- coding: utf-8 -*-
# ============================================================================
# Module 09 - Figure 7: cellular executioners (single-cell + PheWAS)
# Inputs: results/t1_disease_sc_de.csv, t1_heart_sc_de.csv,
#         t1_disease_serpine1_lr.csv, t1_pseudotime_correlation.csv,
#         t1_c1_finngen_r12_phewas_rs7860931.csv;
#         data/derived/figure_assets/Fig7F.jpg
# Output: figures/Fig7_cellular_executioners.{png,pdf}
# ============================================================================

import os as _os
# --- repository paths -----------------------------------------------------------
# REPO resolves to the repository root (code/<module>/this_file.py -> ../..).
# Override with the REPRO_ROOT environment variable if you run from elsewhere.
REPO = _os.environ.get("REPRO_ROOT", _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "..")))
# GBD_RAW_DIR: folder with the raw IHME GBD 2021 csv downloads (see data/README.md)
GBD_RAW_DIR = _os.environ.get("GBD_RAW_DIR", _os.path.join(REPO, "data", "raw", "gbd_2021"))
"""Figure 7 (v2) — Cellular executioners + hero mechanism diagram."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from figstyle import *
import pandas as pd
import numpy as np
from matplotlib.patches import (FancyBboxPatch, FancyArrowPatch, Ellipse, Circle,
                                Polygon, Wedge, PathPatch, Rectangle, Arc)
from matplotlib.path import Path
import matplotlib.patheffects as pe

R = os.path.join(REPO, "results") + os.sep
# MedComm revision: write into the working-copy folder, not the source repo
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")

sc = pd.read_csv(R + "t1_disease_sc_de.csv")
hs = pd.read_csv(R + "t1_heart_sc_de.csv")
lr = pd.read_csv(R + "t1_disease_serpine1_lr.csv")
pt = pd.read_csv(R + "t1_pseudotime_correlation.csv")
phe = pd.read_csv(R + "t1_c1_finngen_r12_phewas_rs7860931.csv")

fig = plt.figure(figsize=(14.5, 13.2))
gs = fig.add_gridspec(3, 3, height_ratios=[1.0, 0.95, 1.5], hspace=0.55, wspace=0.40,
                      left=0.10, right=0.97, top=0.965, bottom=0.05)

# ================= A: COPD lung dotplot =================
axA = fig.add_subplot(gs[0, 0:2])
genes = ["S100A8", "S100A12", "THBS1", "F13A1", "CD163", "SERPINE1", "SPP1", "AREG", "IL1RL1", "IL1R2", "SOCS3", "EGR1", "CDKN1A"]
cts = ["monocyte", "macrophage", "neutrophil", "capillary endothelial cell",
       "fibroblast", "myofibroblast", "endothelial cell", "T cell", "B cell", "dendritic cell"]
d = sc[sc.gene.isin(genes) & sc.cell_type.isin(cts)].copy()
piv_fc = d.pivot_table(index="gene", columns="cell_type", values="log2FC_COPD_vs_Normal")
piv_pct = d.pivot_table(index="gene", columns="cell_type", values="pct_COPD")
piv_fdr = d.pivot_table(index="gene", columns="cell_type", values="fdr")
piv_fc = piv_fc.reindex(index=genes, columns=cts)
piv_pct = piv_pct.reindex(index=genes, columns=cts)
piv_fdr = piv_fdr.reindex(index=genes, columns=cts)
ok_cols = [c for c in cts if piv_fc[c].notna().any()]
cts = ok_cols
piv_fc = piv_fc[cts]; piv_pct = piv_pct[cts]; piv_fdr = piv_fdr[cts]
for i, g in enumerate(genes):
    for j, ct in enumerate(cts):
        fc = piv_fc.loc[g, ct]; pc = piv_pct.loc[g, ct]; fd = piv_fdr.loc[g, ct]
        if np.isnan(fc):
            continue
        axA.scatter(j, i, s=20 + 380 * (pc if not np.isnan(pc) else 0),
                    color=CRIMSON if fc > 0 else SKY,
                    alpha=min(1, abs(fc) / 1.5 + 0.15), edgecolor=INK if fd < 0.05 else "none",
                    lw=1.0, zorder=3)
axA.set_xticks(range(len(cts)))
axA.set_xticklabels([c.replace(" cell", "").replace("capillary endothelial", "capillary EC")
                     for c in cts], rotation=38, ha="right", fontsize=9)
axA.set_yticks(range(len(genes))); axA.set_yticklabels(genes, fontsize=9.5)
axA.text(0.01, 1.03, "color: log₂FC · size: % expressing · ring: FDR < 0.05",
         transform=axA.transAxes, fontsize=8, color=SLATE)
panel_tag(axA, "A", INK, x=-0.075)

# ================= B: heart dotplot =================
axB = fig.add_subplot(gs[0, 1]) if False else fig.add_subplot(gs[0, 2])
genes_h = ["SERPINE1", "THBS1", "SPP1", "ITGAV", "CD163", "S100A8"]
cts_h = ["endothelial cell", "fibroblast", "macrophage"]   # must match csv values
dh = hs[hs.gene.isin(genes_h) & hs.cell_type.isin(cts_h)]
for i, g in enumerate(genes_h):
    for j, ct in enumerate(cts_h):
        row = dh[(dh.gene == g) & (dh.cell_type == ct)]
        if not len(row):
            continue
        row = row.iloc[0]
        fc = row.log2FC_CM_vs_Normal
        axB.scatter(j, i, s=60 + 900 * row.pct_cardiomyopathy,
                    color=CRIMSON if fc > 0 else SKY,
                    alpha=min(1, abs(fc) / 0.5 + 0.2),
                    edgecolor=INK if row.fdr_CM < 0.05 else "none", lw=1.0, zorder=3)
axB.set_xticks(range(len(cts_h)))
axB.set_xticklabels(["  endothelial", "fibroblast", "macrophage"], fontsize=9)  # 2-char right shift
axB.set_yticks(range(len(genes_h))); axB.set_yticklabels(genes_h, fontsize=9.5)
panel_tag(axB, "B", INK, x=-0.22)

# ================= C: SERPINE1 LR edges =================
axC = fig.add_subplot(gs[1, 0])
top = lr[lr.receptor.isin(["PLAUR", "LRP1"])]
top = top[~top.sender.str.contains("_ALL_") & ~top.receiver.str.contains("_ALL_")]
top = top.groupby(["sender", "receiver", "receptor", "disease"], as_index=False).score.max()
edges = top[top.disease == "COPD"].nlargest(6, "score")[["sender", "receiver", "receptor"]]
rows = []
for _, e in edges.iterrows():
    for dis in ["Normal", "COPD"]:
        s = top[(top.sender == e.sender) & (top.receiver == e.receiver) &
                (top.receptor == e.receptor) & (top.disease == dis)].score
        rows.append((f"{e.sender[:14]}→{e.receiver[:12]}\n({e.receptor})", dis,
                     s.iloc[0] if len(s) else 0))
lab2idx = {}
for lab, dis, sc_ in rows:
    lab2idx.setdefault(lab, len(lab2idx))
for lab, dis, sc_ in rows:
    yv = lab2idx[lab] + (0.14 if dis == "COPD" else -0.14)
    axC.scatter(sc_, yv, s=90, color=CRIMSON if dis == "COPD" else SLATE, zorder=3,
                edgecolor="white")
axC.set_yticks(range(len(lab2idx))); axC.set_yticklabels(lab2idx.keys(), fontsize=8.2)
axC.set_xlabel("SERPINE1 ligand–receptor score")
axC.legend(handles=[plt.Line2D([], [], marker="o", ls="", color=CRIMSON, label="COPD"),
                    plt.Line2D([], [], marker="o", ls="", color=SLATE, label="Normal")],
           fontsize=8.2, loc="lower right")
strip_spines(axC); grid(axC, "x")
panel_tag(axC, "C", INK, x=-0.30)

# ================= D: pseudotime =================
axD = fig.add_subplot(gs[1, 1])
pa = pt[pt.subset == "all"].sort_values("spearman_rho")
y = np.arange(len(pa))
cols = [CRIMSON if v > 0 else SKY for v in pa.spearman_rho]
axD.barh(y, pa.spearman_rho, color=cols, height=0.62, zorder=3)
axD.axvline(0, color=INK, lw=0.9)
for yi, v in zip(y, pa.spearman_rho):
    axD.text(v + (0.008 if v > 0 else -0.008), yi, f"{v:+.2f}",
             va="center", ha="left" if v > 0 else "right", fontsize=8, fontweight="bold",
             color=CRIMSON if v > 0 else SKY)
axD.set_xlim(-0.24, 0.42)   # headroom keeps +value labels clear of panel E's wider y-labels
axD.set_yticks(y); axD.set_yticklabels(pa.gene, fontsize=9.2)
axD.set_xlabel("Spearman ρ with myeloid pseudotime")
strip_spines(axD); grid(axD, "x")
panel_tag(axD, "D", INK, x=-0.30)

# ================= E: PheWAS safety =================
axE = fig.add_subplot(gs[1, 2])
# MedComm revision: show the 8 smallest-P endpoints with FULL names (manual
# two-line breaks, no ellipsis truncation); "none significant" is unchanged.
_LABEL_BREAKS = {
    "Asthma (only as main-diagnosis) (more control exclusions)":
        "Asthma (only as main-diagnosis)\n(more control exclusions)",
    "Chronic rhinitis, nasopharyngitis and pharyngitis":
        "Chronic rhinitis, nasopharyngitis\nand pharyngitis",
}
top10 = phe.nsmallest(8, "pval").iloc[::-1]
y = np.arange(len(top10))
cols = [TEAL if b < 0 else CRIMSON for b in top10.beta_altC]
axE.scatter(top10.mlogp, y, s=65, color=cols, zorder=3, edgecolor="white")
for yi, (_, r) in zip(y, top10.iterrows()):
    axE.scatter([r.mlogp], [yi], s=65, color=TEAL if r.beta_altC < 0 else CRIMSON, zorder=3)
bonf = -np.log10(2.02e-5)
axE.axvline(bonf, color=CRIMSON, ls="--", lw=1.1)
axE.text(bonf - 0.06, len(top10) - 0.4, "Bonferroni\n2.0×10⁻⁵", fontsize=8,
         color=CRIMSON, va="top", ha="right")
axE.set_yticks(y)
axE.set_yticklabels([_LABEL_BREAKS.get(t, t) for t in top10.phenotype], fontsize=8)
axE.set_xlabel("−log₁₀ P  (rs7860931, SERPINE1 cis)")
axE.set_ylim(-0.7, len(top10) + 0.6)
strip_spines(axE); grid(axE, "x")
axE.text(0.02, 1.04, "2,470 endpoints · none significant — clean safety window\n(8 smallest-P endpoints shown)",
         transform=axE.transAxes, fontsize=8, color=SLATE, va="bottom")
panel_tag(axE, "E", INK, x=-0.32)

# ================= F: professional illustration (user-provided Fig7F.jpg) =================
axF = fig.add_subplot(gs[2, :])
axF.axis("off")
_img = plt.imread(os.path.join(REPO, "data", "derived", "figure_assets", "Fig7F.jpg"))
axF.imshow(_img, aspect="auto")
panel_tag(axF, "F", INK, x=-0.02, y=1.0)

save(fig, "Fig7_cellular_executioners", OUT)
