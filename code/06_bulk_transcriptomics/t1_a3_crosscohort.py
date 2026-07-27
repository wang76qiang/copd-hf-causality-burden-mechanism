# -*- coding: utf-8 -*-
# ============================================================================
# Module 06 - Bulk transcriptomics: 41-gene cross-cohort validation
# Welch t-tests + direction consistency for the 41 shared genes across four
# cohorts: GSE57148 (COPD lung bulk), GSE57338 (HF myocardium bulk),
# CELLxGENE 8fbed309 (COPD lung sc, donor pseudobulk), Reichart 2022 DCM/ACM
# heart atlas (donor pseudobulk).
# Inputs : data/derived/bulk/*.xlsx; data/*.h5ad (NOT shipped, see
#          data/README.md)
# Outputs: results/t1_crosscohort_validation.csv,
#          t1_crosscohort_direction_consistency.csv,
#          figures/t1_crosscohort_heatmap.png
# ============================================================================

import os as _os
# --- repository paths -----------------------------------------------------------
# REPO resolves to the repository root (code/<module>/this_file.py -> ../..).
# Override with the REPRO_ROOT environment variable if you run from elsewhere.
REPO = _os.environ.get("REPRO_ROOT", _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "..")))
# GBD_RAW_DIR: folder with the raw IHME GBD 2021 csv downloads (see data/README.md)
GBD_RAW_DIR = _os.environ.get("GBD_RAW_DIR", _os.path.join(REPO, "data", "raw", "gbd_2021"))
"""T1 A3 route (b): 41-gene signature cross-cohort validation.
Cohorts:
  1. GSE57148  (COPD lung bulk, values from local xlsx)
  2. GSE57338  (HF myocardium bulk, values from local xlsx)
  3. CXG_COPD  (CELLxGENE 8fbed309 COPD lung sc, donor-level pseudobulk)
  4. CXG_DCM   (DCM/ACM heart cell atlas sc, donor-level pseudobulk, DCM vs normal)
Outputs: results/t1_crosscohort_validation.csv, figures/t1_crosscohort_heatmap.png,
         logs/t1_a3_route_notes.md
"""
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind

BASE = REPO

# ---------------- 41-gene list ----------------
gl = pd.read_csv(os.path.join(REPO, "data", "derived", "bulk", "gene_list_41_inflammation_fibrosis.csv"), header=None)
GENES = gl[0].dropna().astype(str).str.strip().tolist()
GENES = [g for g in GENES if g != "SMOC2"] + ["SMOC2"]  # header cell is SMOC2
print(f"{len(GENES)} genes:", GENES)

records = []  # gene, cohort, mean_case, mean_control, effect(log2FC or diff), p, direction

def add_record(gene, cohort, case, ctrl, log_scale):
    case = np.asarray(case, float); ctrl = np.asarray(ctrl, float)
    case = case[~np.isnan(case)]; ctrl = ctrl[~np.isnan(ctrl)]
    if len(case) < 3 or len(ctrl) < 3:
        return
    t, p = ttest_ind(case, ctrl, equal_var=False)
    mc, mn = case.mean(), ctrl.mean()
    eff = (mc - mn) if log_scale else float(np.log2((mc + 1e-9) / (mn + 1e-9)))
    records.append({"gene": gene, "cohort": cohort, "n_case": len(case), "n_control": len(ctrl),
                    "mean_case": float(mc), "mean_control": float(mn),
                    "effect_log2fc": eff, "welch_t": float(t), "pvalue": float(p),
                    "direction": "up" if eff > 0 else ("down" if eff < 0 else "flat")})

# ---------------- 1) GSE57148 COPD lung (linear scale) ----------------
fp = os.path.join(REPO, "data", "derived", "bulk", "GSE57148_GSE57338_keygene_expression_values.xlsx")
df = pd.read_excel(fp, sheet_name="GSE57148", header=None)
blocks = {"S100A8": (1, 2, 3), "S100A12": (5, 6, 7), "IL1RL1": (9, 10, 11), "THBS1": (13, 14, 15)}
for gene, (cs, ct, cv) in blocks.items():
    sub = df.iloc[2:, [cs, ct, cv]].dropna(subset=[cv])
    sub.columns = ["sample", "title", "value"]
    val = pd.to_numeric(sub["value"], errors="coerce")
    case = val[sub["title"].astype(str).str.upper().str.contains("COPD")]
    ctrl = val[sub["title"].astype(str).str.lower().str.contains("control")]
    add_record(gene, "GSE57148_COPD_lung", case.values, ctrl.values, log_scale=False)
# GSE57148 THBS1: Title column empty in main sheet; S100A8两组比较.xlsx holds THBS1 values w/ labels
try:
    th = pd.read_excel(os.path.join(REPO, "data", "derived", "bulk", "GSE57148_THBS1_two_group_comparison.xlsx"))
    th.columns = ["title", "value"]
    tv = pd.to_numeric(th["value"], errors="coerce")
    case = tv[th["title"].astype(str).str.upper().str.contains("COPD")]
    ctrl = tv[th["title"].astype(str).str.lower().str.contains("control")]
    add_record("THBS1", "GSE57148_COPD_lung", case.values, ctrl.values, log_scale=False)
    print("GSE57148 THBS1 from S100A8两组比较.xlsx:", len(case), "vs", len(ctrl))
except Exception as e:
    print("THBS1 fallback failed:", e)

# ---------------- 2) GSE57338 HF myocardium (log scale) ----------------
df2 = pd.read_excel(fp, sheet_name="GSE57338", header=None)
blocks2 = {"CD163": (0, 1, 2), "IL1RL1": (9, 10, 11), "LUM": (13, 14, 15),
           "ASPN": (17, 18, 19), "SMOC2": (21, 22, 23), "ACE": (25, 26, 27)}
for gene, (cs, ct, cv) in blocks2.items():
    sub = df2.iloc[2:, [cs, ct, cv]].dropna(subset=[cv])
    sub.columns = ["sample", "title", "value"]
    val = pd.to_numeric(sub["value"], errors="coerce")
    tl = sub["title"].astype(str).str.lower()
    ctrl = val[tl.str.contains("non-failing|control")]
    case = val[~tl.str.contains("non-failing|control")]
    add_record(gene, "GSE57338_HF_heart", case.values, ctrl.values, log_scale=True)

print("bulk cohorts parsed:", len(records))

# ---------------- 3) CXG COPD pseudobulk ----------------
def pseudobulk_de(h5ad_paths, case_fn, cohort_name, genes):
    """Sum raw counts per donor across given files; CPM-log2; Welch t case vs control donors."""
    import anndata as ad
    donor_sum = {}
    for p in h5ad_paths:
        a = ad.read_h5ad(p)
        if "feature_name" in a.var.columns:
            names = a.var["feature_name"].astype(str).values
        else:
            names = a.var_names.astype(str)
        name2idx = {}
        for i, n in enumerate(names):
            name2idx.setdefault(n, i)
        gidx = {g: name2idx[g] for g in genes if g in name2idx}
        # raw counts
        x0 = a.X[:50, :50]
        x0 = x0.toarray() if sp.issparse(x0) else np.asarray(x0)
        if not np.allclose(x0, np.round(x0)) and a.raw is not None:
            Xc = a.raw.X
            if "feature_name" in a.raw.var.columns:
                rnames = a.raw.var["feature_name"].astype(str).values
            else:
                rnames = a.raw.var_names.astype(str)
            r2idx = {}
            for i, n in enumerate(rnames):
                r2idx.setdefault(n, i)
            gidx = {g: r2idx[g] for g in genes if g in r2idx}
        else:
            Xc = a.X
        gsorted = sorted(gidx)
        cols = [gidx[g] for g in gsorted]
        sub = Xc[:, cols]
        sub = sub.tocsr() if sp.issparse(sub) else sp.csr_matrix(sub)
        gpos = {g: genes.index(g) for g in gsorted}  # position in full list
        donors = a.obs["donor_id"].astype(str).values
        groups = a.obs["disease"].astype(str).map(case_fn).values
        libsize = np.asarray(Xc.sum(axis=1)).ravel()
        for grp in ["case", "control"]:
            for d in np.unique(donors[groups == grp]):
                m = donors == d
                s = np.asarray(sub[m].sum(axis=0)).ravel()
                lib = libsize[m].sum()
                full = np.zeros(len(genes))
                for k, g in enumerate(gsorted):
                    full[gpos[g]] = s[k]
                key = (d, grp)
                if key in donor_sum:
                    donor_sum[key][0] += full
                    donor_sum[key][1] += lib
                else:
                    donor_sum[key] = [full, lib]
        if hasattr(a, "file"):
            a.file.close()
        del a
    # build donor x gene matrix (log2 CPM)
    keys = sorted(donor_sum.keys())
    gsorted = genes
    M = np.zeros((len(keys), len(gsorted)))
    for i, k in enumerate(keys):
        s, lib = donor_sum[k]
        M[i] = np.log2((s / (lib + 1e-9)) * 1e6 + 1)
    grp_arr = np.array([k[1] for k in keys])
    print(cohort_name, "donors:", {g: int((grp_arr == g).sum()) for g in ["case", "control"]})
    for j, g in enumerate(gsorted):
        case = M[grp_arr == "case", j]
        ctrl = M[grp_arr == "control", j]
        add_record(g, cohort_name, case, ctrl, log_scale=True)

pseudobulk_de([BASE + "/data/t1_copd_lung.h5ad"],
              {"chronic obstructive pulmonary disease": "case", "normal": "control"},
              "CXG_COPD_lung", GENES)

pseudobulk_de([BASE + "/data/t1_hf_fibroblasts.h5ad",
               BASE + "/data/t1_hf_endothelial.h5ad",
               BASE + "/data/t1_hf_macrophages.h5ad"],
              {"dilated cardiomyopathy": "case", "normal": "control"},
              "CXG_DCM_heart", GENES)

res = pd.DataFrame(records)
from statsmodels.stats.multitest import multipletests
res["fdr_within_cohort"] = np.nan
for c in res["cohort"].unique():
    m = res["cohort"] == c
    pv = res.loc[m, "pvalue"].fillna(1.0).values
    res.loc[m, "fdr_within_cohort"] = multipletests(pv, method="fdr_bh")[1]
res = res.sort_values(["gene", "cohort"])
res.to_csv(BASE + "/results/t1_crosscohort_validation.csv", index=False)
print(res.to_string(index=False))

# ---------------- heatmap ----------------
cohorts = ["GSE57148_COPD_lung", "CXG_COPD_lung", "GSE57338_HF_heart", "CXG_DCM_heart"]
piv = res.pivot_table(index="gene", columns="cohort", values="effect_log2fc").reindex(GENES)[cohorts]
piv_p = res.pivot_table(index="gene", columns="cohort", values="fdr_within_cohort").reindex(GENES)[cohorts]
fig, ax = plt.subplots(figsize=(6.5, 0.42 * len(GENES) + 1.5))
vmax = np.nanmax(np.abs(piv.values)) if np.isfinite(piv.values).any() else 1
vmax = min(vmax, 3)
im = ax.imshow(piv.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
for i in range(piv.shape[0]):
    for j in range(piv.shape[1]):
        v = piv.values[i, j]
        if np.isnan(v):
            ax.text(j, i, "NA", ha="center", va="center", fontsize=6, color="grey")
        else:
            star = "*" if piv_p.values[i, j] < 0.05 else ""
            ax.text(j, i, f"{v:+.2f}{star}", ha="center", va="center", fontsize=6,
                    color="black" if abs(v) < 0.6 * vmax else "white")
ax.set_xticks(range(len(cohorts)))
ax.set_xticklabels(cohorts, rotation=20, ha="right", fontsize=8)
ax.set_yticks(range(len(GENES)))
ax.set_yticklabels(GENES, fontsize=7)
# mark the 13 highlighted genes
hl = ["SERPINE1", "THBS1", "SPP1", "S100A8", "S100A12", "IL1RL1", "F13A1",
      "CD163", "IL1R2", "AREG", "SOCS3", "EGR1", "CDKN1A"]
for tick, g in zip(ax.get_yticklabels(), GENES):
    if g in hl:
        tick.set_fontweight("bold")
fig.colorbar(im, ax=ax, label="effect (log2FC / log-scale diff), * FDR<0.05", fraction=0.03, pad=0.02)
ax.set_title("41-gene signature: cross-cohort direction validation", fontsize=11)
fig.tight_layout()
fig.savefig(BASE + "/figures/t1_crosscohort_heatmap.png", dpi=300, bbox_inches="tight")
plt.close(fig)

# direction consistency summary
sub = res[res["gene"].isin(hl)]
cons = sub.groupby("gene")["direction"].apply(lambda s: (s == "up").sum() / len(s)).sort_values(ascending=False)
print("\ndirection consistency (fraction 'up' across available cohorts, 13 key genes):")
print(cons.to_string())
cons.to_csv(BASE + "/results/t1_crosscohort_direction_consistency.csv")
print("A3 DONE")
