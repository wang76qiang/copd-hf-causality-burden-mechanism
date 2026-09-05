# -*- coding: utf-8 -*-
"""Figure 8 v2: replication architecture & heterogeneity forensics + interorgan crosstalk."""
import json, csv, math, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, "repo/code/09_figures")
from figstyle import INK, PAPER, TEAL, CRIMSON, AMBER, SLATE, SKY, LIGHT, GREEN, panel_tag, strip_spines, grid, save

os.makedirs("figures_new", exist_ok=True)
pkg = json.load(open("results_new/heterogeneity_package.json"))
cohorts = pkg["cohorts"]; meta = pkg["meta"]
pooled = [c for c in cohorts if c.get("pool", True)]

fig = plt.figure(figsize=(12.5, 9.2))
gs = fig.add_gridspec(2, 2, height_ratios=[1.45, 1.0], width_ratios=[1.55, 1.0],
                      hspace=0.42, wspace=0.34)

# ---------------- Panel A: forest ----------------
axA = fig.add_subplot(gs[0, 0])
order = sorted(range(len(cohorts)), key=lambda i: (cohorts[i]["ancestry"] != "EUR", -cohorts[i]["cases"]))
rows = [cohorts[i] for i in order]
labels = [c["id"] for c in rows]
ors = [math.exp(c["beta"]) for c in rows]
lo = [math.exp(c["beta"] - 1.96 * c["se"]) for c in rows]
hi = [math.exp(c["beta"] + 1.96 * c["se"]) for c in rows]
colors = [TEAL if c["ancestry"] == "EUR" else CRIMSON for c in rows]
y = np.arange(len(rows))[::-1].astype(float)
for yi, o, l, h, c, crow in zip(y, ors, lo, hi, colors, rows):
    axA.plot([l, h], [yi, yi], color=c, lw=1.6, zorder=2)
    axA.scatter([o], [yi], s=44, color=c, zorder=3, edgecolor="white", lw=0.6,
                marker="o" if crow.get("pool", True) else "s")
ann = [f"{o:.3f} ({l:.3f}–{h:.3f})" for o, l, h in zip(ors, lo, hi)]
extra_labels, extra_ann = [], []
extra_y = -1.2
for anc, col, name in [("EUR", TEAL, "European subgroup"), ("EAS", CRIMSON, "East Asian subgroup"), (None, INK, "Overall")]:
    m = meta["all"] if anc is None else meta[f"ancestry_{anc}"]
    o = math.exp(m["fe_beta"]); l = math.exp(m["fe_beta"] - 1.96*m["fe_se"]); h = math.exp(m["fe_beta"] + 1.96*m["fe_se"])
    axA.fill([l, o, h, o], [extra_y, extra_y + 0.30, extra_y, extra_y - 0.30], color=col, zorder=3)
    extra_labels.append(f"{name} (k={m['k']}, I²={m['I2']*100:.0f}%)")
    extra_ann.append(f"{o:.3f} ({l:.3f}–{h:.3f})")
    extra_y -= 1.2
labels_all = labels + extra_labels
ann_all = ann + extra_ann
yall = list(y) + [-1.2, -2.4, -3.6]
axA.axvline(1.0, color=SLATE, lw=0.9, ls="--", zorder=1)
axA.set_yticks(yall); axA.set_yticklabels(labels_all, fontsize=7.4)
axA.set_xscale("log"); axA.set_xlim(0.62, 2.9)
axA.set_xticks([0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4])
axA.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, p: f"{v:g}"))
axA.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
axA.set_xlabel("Odds ratio per log-unit higher COPD liability (log scale)")
for yi, txt in zip(yall, ann_all):
    axA.text(1.55, yi, txt, fontsize=7.2, va="center", color=INK)
axA.text(1.55, len(rows) + 0.1, "OR (95% CI)", fontsize=7.4, fontweight="bold", color=INK)
axA.scatter([], [], marker="s", color=SLATE)
axA.text(1.55, -4.6, "square = sensitivity (not pooled)", fontsize=6.8, color=SLATE, style="italic")
panel_tag(axA, "A")
strip_spines(axA); grid(axA, "x")

# ---------------- Panel B: funnel ----------------
axB = fig.add_subplot(gs[0, 1])
b = np.array([c["beta"] for c in pooled]); se = np.array([c["se"] for c in pooled])
prec = 1/se
colsB = [TEAL if c["ancestry"] == "EUR" else CRIMSON for c in pooled]
axB.scatter(b, prec, s=46, c=colsB, edgecolor="white", lw=0.6, zorder=3)
fe = meta["all"]["fe_beta"]
axB.axvline(fe, color=INK, lw=1.0)
xx = np.concatenate([np.linspace(-0.42, fe - 0.035, 80), np.linspace(fe + 0.035, 0.42, 80)])
axB.plot(xx, 1.96/np.abs(xx - fe), color=SLATE, lw=0.8, ls="--")
axB.set_ylim(0, prec.max() * 1.10); axB.set_xlim(-0.42, 0.42)
axB.set_xlabel("log odds ratio (cohort estimate)")
axB.set_ylabel("Precision (1/SE)")
axB.text(0.04, 0.04, f"Egger P = {meta['funnel_egger']['p']:.2f} (illustrative, k=10)",
         transform=axB.transAxes, fontsize=8, color=INK)
for c in pooled:
    if abs(c["beta"]) > 0.09 or c["se"] < 0.045:
        axB.annotate(c["id"].split(" (")[0].replace("FinnGen R12 ", "FG ").replace("BBJ 2025 ", "BBJ25 "),
                     (c["beta"], 1/c["se"]), fontsize=5.8, color=INK,
                     xytext=(4, 3), textcoords="offset points")
panel_tag(axB, "B")
strip_spines(axB); grid(axB, "y")

# ---------------- Panel C: winner's curse ----------------
axC = fig.add_subplot(gs[1, 0])
a12 = json.load(open("results_new/analysis12_summary.json"))
wc = a12["winners_curse"]["per_snp"]
xs = np.arange(len(wc))
naive = [r["beta_naive"] for r in wc]; corr = [r["beta_corrected"] for r in wc]
width = 0.38
axC.bar(xs - width/2, naive, width, color=SLATE, label="Naive (selected at P<5×10⁻⁸)")
axC.bar(xs + width/2, corr, width, color=AMBER, label="Winner's-curse-corrected")
axC.axhline(0, color=INK, lw=0.8)
axC.set_xticks(xs); axC.set_xticklabels([r["SNP"] for r in wc], rotation=45, ha="right", fontsize=6.8)
axC.set_ylabel("Exposure β (log-OR per allele)")
axC.set_ylim(-0.175, 0.285)
axC.legend(fontsize=7, loc="upper left", frameon=True, edgecolor=LIGHT)
axC.text(0.985, 0.04, f"IVW OR: {math.exp(a12['winners_curse']['ivw_naive']['ivw_fe']):.3f} → "
                     f"{math.exp(a12['winners_curse']['ivw_corrected']['ivw_fe']):.3f} after correction;\n"
                     f"independent FinnGen-exposure re-estimate: OR 1.19 (1.10–1.28)",
         transform=axC.transAxes, fontsize=7.6, color=INK, ha="right", va="bottom")
panel_tag(axC, "C")
strip_spines(axC); grid(axC, "y")

# ---------------- Panel D: crosstalk heatmap ----------------
axD = fig.add_subplot(gs[1, 1])
edges = pd.read_csv("results_new/b4_interorgan_edges_lung_to_heart.csv")
edges["pair"] = edges["ligand"] + "→" + edges["receptor"]
edges["max_score"] = edges[["score_normal", "score_disease"]].max(axis=1)
edges = edges[edges.groupby("pair")["max_score"].transform("max") >= 0.005]
pair_delta = edges.groupby("pair")["delta_log2"].apply(lambda s: s.iloc[s.abs().argmax()])
top_pairs = pair_delta.reindex(pair_delta.abs().sort_values(ascending=False).index).head(12).index.tolist()
senders_pref = ["macrophage", "monocyte", "capillary endothelial cell", "fibroblast"]
senders = [s for s in senders_pref if s in set(edges["sender"])]
mat = np.full((len(top_pairs), len(senders)), np.nan)
for i, p in enumerate(top_pairs):
    sub = edges[edges["pair"] == p]
    for j, s in enumerate(senders):
        v = sub[sub["sender"] == s]["delta_log2"]
        if len(v): mat[i, j] = v.mean()
keep = ~np.all(np.isnan(mat), axis=1)
mat = mat[keep]; top_pairs = [p for p, k in zip(top_pairs, keep) if k]
vmax = np.nanmax(np.abs(mat))
im = axD.imshow(mat, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
axD.set_xticks(range(len(senders))); axD.set_xticklabels([s.replace(" ", "\n") for s in senders], fontsize=7)
axD.set_yticks(range(len(top_pairs))); axD.set_yticklabels(top_pairs, fontsize=7)
for i in range(mat.shape[0]):
    for j in range(mat.shape[1]):
        if np.isfinite(mat[i, j]):
            axD.text(j, i, f"{mat[i,j]:+.2f}", ha="center", va="center", fontsize=6.5,
                     color="white" if abs(mat[i, j]) > 0.6*vmax else INK)
cb = fig.colorbar(im, ax=axD, shrink=0.8)
cb.set_label("Δ log₂(edge score)\nCOPD–failing vs normal–normal", fontsize=7)
cb.ax.tick_params(labelsize=6.5)
axD.set_xlabel("lung sender cell type (receiver-averaged)", fontsize=8)
panel_tag(axD, "D")
strip_spines(axD, keep=("left", "bottom"))

save(fig, "Figure_8_replication_heterogeneity_crosstalk", "figures_new")
