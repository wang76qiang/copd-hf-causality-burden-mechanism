# -*- coding: utf-8 -*-
# ============================================================================
# Module 07 - Single-cell: heart (Reichart 2022 DCM/ACM atlas, CELLxGENE)
# THBS1/SPP1 in fibroblasts, SERPINE1 in endothelial cells, macrophage panel;
# DCM vs normal Wilcoxon DE across three heart cell types.
# Input : data/t1_hf_{fibroblasts,endothelial,macrophages}.h5ad (NOT shipped)
# Output: results/t1_heart_sc_de.csv, figures/t1_heart_sc_dotplot.png
# ============================================================================

import os as _os
# --- repository paths -----------------------------------------------------------
# REPO resolves to the repository root (code/<module>/this_file.py -> ../..).
# Override with the REPRO_ROOT environment variable if you run from elsewhere.
REPO = _os.environ.get("REPRO_ROOT", _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "..")))
# GBD_RAW_DIR: folder with the raw IHME GBD 2021 csv downloads (see data/README.md)
GBD_RAW_DIR = _os.environ.get("GBD_RAW_DIR", _os.path.join(REPO, "data", "raw", "gbd_2021"))
"""T1 A2 step 5: heart-side validation in DCM/ACM heart cell atlas (Reichart 2022, CELLxGENE).
THBS1/SPP1 in fibroblasts, SERPINE1 in endothelial cells, macrophage panel;
DCM vs normal Wilcoxon.
Outputs: results/t1_heart_sc_de.csv, figures/t1_heart_sc_dotplot.png
"""
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import ranksums
from statsmodels.stats.multitest import multipletests

BASE = REPO
FILES = {"fibroblast": "t1_hf_fibroblasts.h5ad",
         "endothelial cell": "t1_hf_endothelial.h5ad",
         "macrophage": "t1_hf_macrophages.h5ad"}
GENES = ["SERPINE1", "THBS1", "SPP1", "S100A8", "S100A12", "IL1RL1", "F13A1",
         "CD163", "IL1R2", "AREG", "SOCS3", "EGR1", "CDKN1A"]
RECEPTORS = ["LRP1", "PLAUR", "ITGAV", "ITGB3", "VTN", "PLAU"]

def load_norm(path):
    a = sc.read_h5ad(path)
    if "feature_name" in a.var.columns:
        a.var_names = a.var["feature_name"].astype(str)
        a.var_names_make_unique()
    x0 = a.X[:100, :100]
    x0 = x0.toarray() if sp.issparse(x0) else np.asarray(x0)
    if not np.allclose(x0, np.round(x0)) and a.raw is not None:
        raw = a.raw.to_adata()
        a.X = raw.X
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    return a

def gvec(a, gene, mask):
    idx = a.var_names.get_loc(gene)
    col = a.X[mask, idx]
    return np.asarray(col.todense()).ravel()

rows = []
adatas = {}
for ct_label, f in FILES.items():
    print("loading", f)
    a = load_norm(BASE + "/data/" + f)
    a.obs["group"] = np.where(a.obs["disease"] == "normal", "Normal", "Cardiomyopathy")
    a.obs["group2"] = np.where(a.obs["disease"] == "dilated cardiomyopathy", "DCM", a.obs["group"])
    adatas[ct_label] = a
    print(ct_label, a.shape, a.obs["group"].value_counts().to_dict())
    genes_present = [g for g in GENES + RECEPTORS if g in a.var_names]
    print("genes present:", genes_present)
    m_d = (a.obs["group"] == "Cardiomyopathy").values
    m_n = (a.obs["group"] == "Normal").values
    m_dcm = (a.obs["disease"] == "dilated cardiomyopathy").values
    for g in genes_present:
        vd = gvec(a, g, m_d); vn = gvec(a, g, m_n)
        stat, p = ranksums(vd, vn)
        lfc = float(np.log2((np.expm1(vd.mean()) + 1e-9) / (np.expm1(vn.mean()) + 1e-9)))
        vdcm = gvec(a, g, m_dcm)
        stat2, p2 = ranksums(vdcm, vn)
        lfc2 = float(np.log2((np.expm1(vdcm.mean()) + 1e-9) / (np.expm1(vn.mean()) + 1e-9)))
        rows.append({"cell_type": ct_label, "gene": g,
                     "mean_cardiomyopathy": float(vd.mean()), "mean_normal": float(vn.mean()),
                     "pct_cardiomyopathy": float((vd > 0).mean()), "pct_normal": float((vn > 0).mean()),
                     "log2FC_CM_vs_Normal": lfc, "pvalue_CM": float(p),
                     "mean_DCM": float(vdcm.mean()), "log2FC_DCM_vs_Normal": lfc2, "pvalue_DCM": float(p2),
                     "n_CM": int(m_d.sum()), "n_normal": int(m_n.sum())})

de = pd.DataFrame(rows)
de["fdr_CM"] = multipletests(de["pvalue_CM"].values, method="fdr_bh")[1]
de["fdr_DCM"] = multipletests(de["pvalue_DCM"].values, method="fdr_bh")[1]
de.to_csv(BASE + "/results/t1_heart_sc_de.csv", index=False)
print(de.to_string(index=False))

# dotplot: genes x (cell type x group)
panels = []
for ct_label, a in adatas.items():
    for grp in ["Normal", "Cardiomyopathy"]:
        panels.append((ct_label, grp, a))
genes_plot = [g for g in GENES if any(g in a.var_names for a in adatas.values())]
fig, axes = plt.subplots(1, len(panels), figsize=(3.2 * len(panels), 0.5 + 0.45 * len(genes_plot) + 1.5), sharey=True)
for ax, (ct_label, grp, a) in zip(axes, panels):
    m = (a.obs["group"] == grp).values
    means = np.zeros(len(genes_plot)); pcts = np.zeros(len(genes_plot))
    for i, g in enumerate(genes_plot):
        if g in a.var_names:
            v = gvec(a, g, m)
            means[i] = v.mean(); pcts[i] = (v > 0).mean()
    for i in range(len(genes_plot)):
        ax.scatter(0, i, s=200 * pcts[i] + 5, c=means[i], cmap="Reds", vmin=0, vmax=2.5,
                   edgecolors="grey", linewidths=0.3)
    ax.set_xlim(-0.5, 0.5); ax.set_xticks([])
    ax.set_title(f"{ct_label}\n{grp}", fontsize=9)
    ax.set_yticks(range(len(genes_plot)))
    ax.set_yticklabels(genes_plot, fontsize=9)
    ax.invert_yaxis()
sm = plt.cm.ScalarMappable(cmap="Reds", norm=plt.Normalize(0, 2.5))
fig.colorbar(sm, ax=axes, label="mean log1p(CP10k)", fraction=0.02, pad=0.05)
fig.suptitle("DCM/ACM heart cell atlas (Reichart 2022): axis genes; dot size = fraction expressing", fontsize=11)
fig.tight_layout(rect=[0, 0, 0.94, 0.93])
fig.savefig(BASE + "/figures/t1_heart_sc_dotplot.png", dpi=300, bbox_inches="tight")
plt.close(fig)
print("HEART DONE")
