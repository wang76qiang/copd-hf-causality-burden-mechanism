# -*- coding: utf-8 -*-
# ============================================================================
# Module 07 - Single-cell: COPD lung (CELLxGENE dataset 8fbed309)
# Disease-state cell-level Wilcoxon DE, SERPINE1 ligand-receptor scoring and
# diffusion-pseudotime analysis in COPD vs control lung single cells.
# Input : data/t1_copd_lung.h5ad (NOT shipped, see data/README.md)
# Output: results/t1_disease_sc_de.csv, t1_disease_serpine1_lr.csv,
#         t1_disease_serpine1_lr_full.csv, t1_pseudotime_correlation.csv,
#         figures/t1_disease_sc_dotplot.png, t1_pseudotime_macrophage.png
# ============================================================================

import os as _os
# --- repository paths -----------------------------------------------------------
# REPO resolves to the repository root (code/<module>/this_file.py -> ../..).
# Override with the REPRO_ROOT environment variable if you run from elsewhere.
REPO = _os.environ.get("REPRO_ROOT", _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "..")))
# GBD_RAW_DIR: folder with the raw IHME GBD 2021 csv downloads (see data/README.md)
GBD_RAW_DIR = _os.environ.get("GBD_RAW_DIR", _os.path.join(REPO, "data", "raw", "gbd_2021"))
"""T1 A2: COPD disease-state single-cell validation (CELLxGENE dataset 8fbed309).
Outputs:
  results/t1_disease_sc_de.csv
  figures/t1_disease_sc_dotplot.png
  results/t1_disease_serpine1_lr.csv
  figures/t1_pseudotime_macrophage.png
  results/t1_pseudotime_correlation.csv
"""
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import ranksums, spearmanr

BASE = REPO
DATA = BASE + "/data/t1_copd_lung.h5ad"

GENES = ["SERPINE1", "THBS1", "SPP1", "S100A8", "S100A12", "IL1RL1", "F13A1",
         "CD163", "IL1R2", "AREG", "SOCS3", "EGR1", "CDKN1A"]
RECEPTORS = ["LRP1", "PLAUR", "ITGAV", "ITGB3", "VTN", "PLAU"]

print("Loading", DATA)
adata = sc.read_h5ad(DATA)
print("shape:", adata.shape)

# var_names are Ensembl IDs; use feature_name (gene symbol)
if "feature_name" in adata.var.columns:
    adata.var["gene_symbol"] = adata.var["feature_name"]
    adata.var_names = adata.var["gene_symbol"].astype(str)
    adata.var_names_make_unique()

# X should be raw counts (CELLxGENE primary data). Verify integer-ness
x0 = adata.X[:100, :100]
x0 = x0.toarray() if sp.issparse(x0) else np.asarray(x0)
is_counts = np.allclose(x0, np.round(x0))
print("X looks like raw counts:", is_counts)
if not is_counts and adata.raw is not None:
    print("using adata.raw as counts source")
    raw = adata.raw.to_adata()
    raw.var_names = adata.var_names
    adata.X = raw.X

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

adata.obs["disease_simple"] = adata.obs["disease"].map(
    {"chronic obstructive pulmonary disease": "COPD", "normal": "Normal"})
print(adata.obs["disease_simple"].value_counts().to_dict())

# restrict to cell types with >=100 cells
ct_counts = adata.obs["cell_type"].value_counts()
keep_cts = ct_counts[ct_counts >= 100].index.tolist()
print(f"{len(keep_cts)} cell types with >=100 cells:", keep_cts)

genes_present = [g for g in GENES if g in adata.var_names]
print("genes present:", genes_present)
recs = ["LRP1", "PLAUR", "ITGAV", "ITGB3", "VTN", "PLAU"]
recs_present = [g for g in recs if g in adata.var_names]
print("receptors present:", recs_present)

X = adata.X
if not sp.issparse(X):
    X = sp.csr_matrix(X)

def gvec(gene, mask):
    idx = adata.var_names.get_loc(gene)
    col = X[mask, idx]
    return np.asarray(col.todense()).ravel()

# ---------------------------------------------------------------------------
# 1) DE COPD vs Normal per cell type (Wilcoxon rank-sum on log1p CPM)
# ---------------------------------------------------------------------------
rows = []
obs = adata.obs
for ct in keep_cts:
    m_ct = (obs["cell_type"] == ct).values
    m_copd = m_ct & (obs["disease_simple"] == "COPD").values
    m_norm = m_ct & (obs["disease_simple"] == "Normal").values
    n_copd, n_norm = int(m_copd.sum()), int(m_norm.sum())
    if n_copd < 20 or n_norm < 20:
        continue
    for g in genes_present:
        a = gvec(g, m_copd)
        b = gvec(g, m_norm)
        stat, p = ranksums(a, b)
        lfc = float(np.log2((np.expm1(a.mean()) + 1e-9) / (np.expm1(b.mean()) + 1e-9)))
        rows.append({"cell_type": ct, "gene": g,
                     "mean_COPD": float(a.mean()), "mean_Normal": float(b.mean()),
                     "pct_COPD": float((a > 0).mean()), "pct_Normal": float((b > 0).mean()),
                     "log2FC_COPD_vs_Normal": lfc, "wilcoxon_stat": float(stat),
                     "pvalue": float(p), "n_COPD": n_copd, "n_Normal": n_norm})
de = pd.DataFrame(rows)
# BH FDR within all tests
from statsmodels.stats.multitest import multipletests
de["fdr"] = multipletests(de["pvalue"].values, method="fdr_bh")[1]
de = de.sort_values(["gene", "cell_type"])
de.to_csv(BASE + "/results/t1_disease_sc_de.csv", index=False)
print("DE rows:", len(de))

# key check: SERPINE1 in endothelial / fibroblast
key = de[de["gene"] == "SERPINE1"]
print("\nSERPINE1 COPD vs Normal by cell type:")
print(key[["cell_type", "mean_COPD", "mean_Normal", "log2FC_COPD_vs_Normal", "pvalue", "fdr"]].to_string(index=False))

# ---------------------------------------------------------------------------
# 2) Dotplot: genes x cell types, split Normal / COPD
# ---------------------------------------------------------------------------
top_cts = [ct for ct in keep_cts
           if ((obs["cell_type"] == ct) & (obs["disease_simple"] == "COPD")).sum() >= 20
           and ((obs["cell_type"] == ct) & (obs["disease_simple"] == "Normal")).sum() >= 20]
fig, axes = plt.subplots(1, 2, figsize=(4 + 1.1 * len(top_cts), 0.5 + 0.45 * len(genes_present) + 1.5),
                         sharey=True)
for ax, dis in zip(axes, ["Normal", "COPD"]):
    means = np.zeros((len(genes_present), len(top_cts)))
    pcts = np.zeros_like(means)
    for j, ct in enumerate(top_cts):
        m = ((obs["cell_type"] == ct) & (obs["disease_simple"] == dis)).values
        for i, g in enumerate(genes_present):
            v = gvec(g, m)
            means[i, j] = v.mean()
            pcts[i, j] = (v > 0).mean()
    # z-score per gene across cell types (for color comparability)
    mz = (means - means.mean(axis=1, keepdims=True)) / (means.std(axis=1, keepdims=True) + 1e-9)
    for i in range(len(genes_present)):
        for j in range(len(top_cts)):
            ax.scatter(j, i, s=120 * pcts[i, j] + 4, c=mz[i, j], cmap="RdBu_r",
                       vmin=-2, vmax=2, edgecolors="grey", linewidths=0.3)
    ax.set_xticks(range(len(top_cts)))
    ax.set_xticklabels([c.replace(" cell", "").replace("pulmonary ", "").replace("positive, alpha-beta ", "+")
                        for c in top_cts], rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(genes_present)))
    ax.set_yticklabels(genes_present, fontsize=9)
    ax.set_title(f"{dis} (n={int((obs['disease_simple']==dis).sum())} cells)", fontsize=10)
sm = plt.cm.ScalarMappable(cmap="RdBu_r", norm=plt.Normalize(-2, 2))
fig.colorbar(sm, ax=axes, label="z-scored mean log1p(CP10k)", fraction=0.02, pad=0.08)
fig.suptitle("COPD lung (CELLxGENE 8fbed309): 13-gene axis expression by cell type", fontsize=11)
fig.text(0.98, 0.5, "dot size = fraction expressing", rotation=90, va="center", fontsize=8)
fig.tight_layout(rect=[0, 0, 0.93, 0.95])
fig.savefig(BASE + "/figures/t1_disease_sc_dotplot.png", dpi=300, bbox_inches="tight")
plt.close(fig)
print("dotplot saved")

# ---------------------------------------------------------------------------
# 3) Custom SERPINE1 ligand-receptor scoring per disease
#    score = mean(ligand in sender ct) * mean(receptor in receiver ct)
# ---------------------------------------------------------------------------
lr_pairs = [("SERPINE1", r) for r in recs_present]
rows = []
for dis in ["Normal", "COPD"]:
    m_dis = (obs["disease_simple"] == dis).values
    cts_dis = [ct for ct in keep_cts if (m_dis & (obs["cell_type"] == ct).values).sum() >= 20]
    lig_mean = {}
    for ct in cts_dis:
        lig_mean[ct] = gvec("SERPINE1", m_dis & (obs["cell_type"] == ct).values).mean()
    rec_mean = {}
    for ct in cts_dis:
        m = m_dis & (obs["cell_type"] == ct).values
        for (l, r) in lr_pairs:
            rec_mean[(r, ct)] = gvec(r, m).mean()
    for s in cts_dis:
        for rc in cts_dis:
            for (l, r) in lr_pairs:
                rows.append({"disease": dis, "ligand": l, "receptor": r,
                             "sender": s, "receiver": rc,
                             "mean_ligand": float(lig_mean[s]),
                             "mean_receptor": float(rec_mean[(r, rc)]),
                             "score": float(lig_mean[s] * rec_mean[(r, rc)])})
lr = pd.DataFrame(rows)
lr.to_csv(BASE + "/results/t1_disease_serpine1_lr_full.csv", index=False)

# summarise key edges: sender in {fibroblast, endothelial}, receptor PLAUR/LRP1
summary_rows = []
for dis in ["Normal", "COPD"]:
    sub = lr[(lr["disease"] == dis) & lr["receptor"].isin(["PLAUR", "LRP1"])]
    tot = sub.groupby(["receptor"])["score"].sum()
    for r, v in tot.items():
        summary_rows.append({"disease": dis, "receptor": r, "total_score": v})
key_edges = lr[(lr["sender"].str.contains("fibroblast|endothelial|Endothelial|Fibroblast", case=False))
               & (lr["receptor"].isin(["PLAUR", "LRP1"]))]
piv = (key_edges.groupby(["sender", "receiver", "receptor", "disease"])["score"].sum()
       .reset_index())
summ = pd.DataFrame(summary_rows)
out = pd.concat([piv.assign(total_score=np.nan)[["disease", "sender", "receiver", "receptor", "score"]],
                 summ.rename(columns={"total_score": "score"})
                 .assign(sender="__ALL__", receiver="__ALL__")[["disease", "sender", "receiver", "receptor", "score"]]])
out.to_csv(BASE + "/results/t1_disease_serpine1_lr.csv", index=False)
print("\nSERPINE1 LR totals (sum over all sender->receiver):")
print(lr.groupby(["disease", "receptor"])["score"].sum().unstack().to_string())

# ---------------------------------------------------------------------------
# 4) C3: diffusion pseudotime on monocytes + macrophages
# ---------------------------------------------------------------------------
mono_cts = [ct for ct in adata.obs["cell_type"].unique()
            if any(k in ct.lower() for k in ["monocyte", "macrophage", "dendritic"])]
print("\npseudotime cell types:", mono_cts)
am = adata[adata.obs["cell_type"].isin(mono_cts)].copy()
print("myeloid cells:", am.shape)
sc.pp.highly_variable_genes(am, n_top_genes=2000, flavor="seurat")
am = am[:, am.var.highly_variable].copy()
sc.pp.scale(am, max_value=10)
sc.tl.pca(am, svd_solver="arpack")
sc.pp.neighbors(am, n_neighbors=15, n_pcs=30)
sc.tl.umap(am)
sc.tl.diffmap(am)
# root: classical monocyte-like cell = extreme of DC1 among monocytes
mono_mask = am.obs["cell_type"].str.lower().str.contains("monocyte").values
dc1 = am.obsm["X_diffmap"][:, 1]
cand = np.where(mono_mask)[0]
iroot = cand[np.argmax(dc1[cand])] if len(cand) else int(np.argmax(dc1))
am.uns["iroot"] = iroot
sc.tl.dpt(am)
am.obs["dpt"] = am.obs["dpt_pseudotime"]

# full-gene expression for correlation (use normalized values from adata)
corr_rows = []
for g in ["S100A8", "S100A12", "SPP1", "THBS1", "SERPINE1", "CD163", "F13A1", "IL1R2", "AREG"]:
    if g not in adata.var_names:
        continue
    idx = adata.var_names.get_loc(g)
    v = np.asarray(X[adata.obs_names.isin(am.obs_names), idx].todense()).ravel()
    for subset_name, mask in [("all", np.ones(am.n_obs, bool)),
                              ("COPD", (am.obs["disease_simple"] == "COPD").values),
                              ("Normal", (am.obs["disease_simple"] == "Normal").values)]:
        if mask.sum() < 50:
            continue
        rho, p = spearmanr(am.obs["dpt"].values[mask], v[mask])
        corr_rows.append({"gene": g, "subset": subset_name, "n": int(mask.sum()),
                          "spearman_rho": float(rho), "pvalue": float(p)})
corr = pd.DataFrame(corr_rows)
corr.to_csv(BASE + "/results/t1_pseudotime_correlation.csv", index=False)
print("\npseudotime correlations:")
print(corr.to_string(index=False))

# figure: UMAP colored by dpt / cell type / S100A8 / SPP1
am.var_names = am.var_names.astype(str)
expr_df = sc.get.obs_df(adata[am.obs_names], keys=["S100A8", "SPP1"] if "SPP1" in adata.var_names else ["S100A8"])
am.obs["S100A8_expr"] = expr_df["S100A8"].values
if "SPP1" in adata.var_names:
    am.obs["SPP1_expr"] = expr_df["SPP1"].values
colors = ["cell_type", "disease_simple", "dpt", "S100A8_expr"] + (["SPP1_expr"] if "SPP1" in adata.var_names else [])
fig, axes = plt.subplots(1, len(colors), figsize=(4.5 * len(colors), 4))
for ax, c in zip(np.atleast_1d(axes), colors):
    sc.pl.umap(am, color=c, ax=ax, show=False, s=6)
fig.tight_layout()
fig.savefig(BASE + "/figures/t1_pseudotime_macrophage.png", dpi=300, bbox_inches="tight")
plt.close(fig)
print("\nALL DONE")
