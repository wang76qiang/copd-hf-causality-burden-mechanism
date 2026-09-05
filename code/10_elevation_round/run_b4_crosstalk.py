# -*- coding: utf-8 -*-
"""B4 v2 (memory-safe, h5py streaming): Interorgan lung<->heart ligand-receptor crosstalk.
Streams CSR matrices in row-blocks; computes log1p-CP10k means per cell type x condition
and per donor, then curated interorgan edge scores (disease vs normal)."""
import numpy as np, pandas as pd, h5py, os, warnings, json
warnings.filterwarnings("ignore")
from scipy.stats import ttest_ind

OUT = "results_new"; os.makedirs(OUT, exist_ok=True)
BLOCK = 4000

LIGANDS = ["SERPINE1","THBS1","S100A8","S100A9","S100A12","SPP1","FN1","TGFB1","AREG","IL6","CCL2","CXCL8","IL33","PLAU"]
RECEPTORS = ["PLAUR","LRP1","VTN","CD36","CD47","TLR4","AGER","ITGAV","ITGB1","ITGB3","ITGA5","CD44",
             "TGFBR1","TGFBR2","EGFR","IL6R","IL6ST","CCR2","CXCR1","CXCR2","IL1RL1","IL1RAP"]
PAIRS = [("SERPINE1","PLAUR"),("SERPINE1","LRP1"),("SERPINE1","VTN"),("PLAU","PLAUR"),
         ("THBS1","CD36"),("THBS1","CD47"),("THBS1","LRP1"),("THBS1","ITGAV"),
         ("S100A8","TLR4"),("S100A8","AGER"),("S100A9","TLR4"),("S100A9","AGER"),
         ("S100A12","AGER"),("S100A12","TLR4"),
         ("SPP1","ITGAV"),("SPP1","CD44"),("SPP1","ITGB1"),
         ("FN1","ITGAV"),("FN1","ITGA5"),("FN1","ITGB1"),
         ("TGFB1","TGFBR1"),("TGFB1","TGFBR2"),
         ("AREG","EGFR"),("IL6","IL6R"),("IL6","IL6ST"),
         ("CCL2","CCR2"),("CXCL8","CXCR1"),("CXCL8","CXCR2"),
         ("IL33","IL1RL1"),("IL33","IL1RAP")]

def read_str_list(ds):
    return [x.decode() if isinstance(x, bytes) else str(x) for x in ds[:]]

def load_panel(path, genes, disease_map, disease_col="disease"):
    """Stream h5ad: returns obs df (donor, cond, ct) + dense log1pCP10k matrix for panel genes."""
    with h5py.File(path, "r") as f:
        # obs
        obs = {}
        for col in [disease_col, "donor_id", "cell_type"]:
            g = f["obs"][col]
            if isinstance(g, h5py.Group):  # categorical
                cats = read_str_list(g["categories"])
                codes = g["codes"][:]
                obs[col] = [cats[c] if c >= 0 else "" for c in codes]
            else:
                obs[col] = read_str_list(g)
        n = len(obs[disease_col])
        # var feature_name
        var = f["var"]
        if "feature_name" in var:
            fn = var["feature_name"]
            names = read_str_list(fn) if isinstance(fn, h5py.Dataset) else read_str_list(fn["categories"])
        else:
            raise RuntimeError("no feature_name")
        gidx = {}
        for g in genes:
            try: gidx[g] = names.index(g)
            except ValueError: pass
        # choose layer: raw if exists (counts), else X
        layer = "raw" if "raw" in f else "X"
        X = f[layer]["X"]
        data, indices, indptr = X["data"], X["indices"], X["indptr"]
        indptr = indptr[:]
        totals = np.zeros(n)
        panel = {g: np.zeros(n) for g in gidx}
        gi_set = {v: k for k, v in gidx.items()}
        for s in range(0, n, BLOCK):
            e = min(s + BLOCK, n)
            p0, p1 = int(indptr[s]), int(indptr[e])
            d = data[p0:p1]; ii = indices[p0:p1]
            # per-row pointer offsets
            row_ptr = indptr[s:e+1] - p0
            counts = np.diff(row_ptr)  # nnz per row, not sums; need sums -> accumulate
            sums = np.zeros(e - s)
            np.add.at(sums, np.repeat(np.arange(e - s), counts), d)
            totals[s:e] = sums
            # panel extraction
            keep_mask = np.isin(ii, list(gi_set.keys()))
            dk, ik = d[keep_mask], ii[keep_mask]
            rows_in_block = np.repeat(np.arange(e - s), counts)[keep_mask]
            for gcol, gname in gi_set.items():
                m = ik == gcol
                np.add.at(panel[gname], s + rows_in_block[m], dk[m])
        # log1p CP10k
        sf = 1e4 / np.maximum(totals, 1)
        M = np.column_stack([np.log1p(panel[g] * sf) for g in gidx])
        df = pd.DataFrame({"donor": obs["donor_id"], "cond": [disease_map.get(x) for x in obs[disease_col]],
                           "ct": obs["cell_type"]})
        return df, M, list(gidx.keys())

def summarize(df, M, genes):
    """means per (gene, ct, cond) and per (gene, donor, cond)."""
    res, donor = {}, {}
    mask_valid = df["cond"].notna().values
    df = df[mask_valid].reset_index(drop=True); M = M[mask_valid]
    for j, g in enumerate(genes):
        col = M[:, j]
        for ct in df["ct"].unique():
            for cond in ["COPD", "Normal"]:
                m = ((df["ct"] == ct) & (df["cond"] == cond)).values
                if m.sum() >= 20:
                    res[(g, ct, cond)] = float(col[m].mean())
        for dn in df["donor"].unique():
            for cond in ["COPD", "Normal"]:
                m = ((df["donor"] == dn) & (df["cond"] == cond)).values
                if m.sum() >= 10:
                    donor.setdefault(g, {}).setdefault(cond, []).append(float(col[m].mean()))
    return res, donor

print("=== lung ===")
lung_df, lung_M, lung_genes = load_panel("sc_data/t1_copd_lung.h5ad", LIGANDS + RECEPTORS,
    {"chronic obstructive pulmonary disease": "COPD", "normal": "Normal"})
print(lung_df.shape, lung_df["ct"].value_counts().to_dict())
lung_cts = [c for c in lung_df["ct"].value_counts().index if (lung_df["ct"] == c).sum() >= 100]
lung_mean, lung_donor = summarize(lung_df, lung_M, lung_genes)

heart_maps = {"dilated cardiomyopathy": "COPD", "arrhythmogenic right ventricular cardiomyopathy": "COPD",
              "non-compaction cardiomyopathy": "COPD", "normal": "Normal"}
heart_mean, heart_donor = {}, {}   # heart_donor[(tag, gene)][cond] = [values]
heart_genes_present = {}
heart_cts = []
for f, tag in [("sc_data/t1_hf_fibroblasts.h5ad", "cardiac fibroblast"),
               ("sc_data/t1_hf_endothelial.h5ad", "cardiac endothelial cell"),
               ("sc_data/t1_hf_macrophages.h5ad", "cardiac macrophage")]:
    print("===", tag, "===")
    hdf, hM, hgenes = load_panel(f, LIGANDS + RECEPTORS, heart_maps)
    hdf["ct"] = tag
    print(hdf.shape, hdf["cond"].value_counts().to_dict(), "| panel genes present:", len(hgenes))
    heart_genes_present[tag] = set(hgenes)
    m, dd = summarize(hdf, hM, hgenes)
    heart_mean.update(m)
    for g, cd in dd.items():
        for cond, vals in cd.items():
            heart_donor.setdefault((tag, g), {}).setdefault(cond, []).extend(vals)
    heart_cts.append(tag)

# ---- edges lung -> heart ----
lung_genes_present = set(lung_genes)
rows = []
for (lig, rec) in PAIRS:
    if lig not in lung_genes_present:
        continue
    for sct in lung_cts:
        ln = lung_mean.get((lig, sct, "Normal")); ld = lung_mean.get((lig, sct, "COPD"))
        if ln is None or ld is None: continue
        for rct in heart_cts:
            if rec not in heart_genes_present[rct]:
                continue
            rn = heart_mean.get((rec, rct, "Normal")); rd = heart_mean.get((rec, rct, "COPD"))
            if rn is None or rd is None: continue
            rows.append({"direction": "lung->heart", "ligand": lig, "receptor": rec,
                         "sender": sct, "receiver": rct,
                         "ligand_normal": ln, "ligand_disease": ld,
                         "receptor_normal": rn, "receptor_disease": rd,
                         "score_normal": ln * rn, "score_disease": ld * rd,
                         "delta": ld * rd - ln * rn,
                         "delta_log2": float(np.log2((ld * rd + 1e-6) / (ln * rn + 1e-6)))})
edges = pd.DataFrame(rows)
edges.to_csv(f"{OUT}/b4_interorgan_edges_lung_to_heart.csv", index=False)
print("edges:", len(edges))
top = edges.reindex(edges["delta"].abs().sort_values(ascending=False).index).head(25)
print(top[["ligand","receptor","sender","receiver","score_normal","score_disease","delta"]].to_string(index=False))

# ---- donor-level tests ----
donor_rows = []
for g in LIGANDS:
    d = lung_donor.get(g, {})
    a, b = d.get("COPD", []), d.get("Normal", [])
    if len(a) >= 4 and len(b) >= 4:
        t, p = ttest_ind(a, b, equal_var=False)
        donor_rows.append({"organ": "lung", "cell_type": "all lung (pooled)", "gene": g, "n_disease": len(a), "n_control": len(b),
                           "mean_disease": float(np.mean(a)), "mean_control": float(np.mean(b)), "welch_p": float(p)})
for g in RECEPTORS:
    for tag in heart_cts:
        d = heart_donor.get((tag, g), {})
        a, b = d.get("COPD", []), d.get("Normal", [])
        if len(a) >= 4 and len(b) >= 4:
            t, p = ttest_ind(a, b, equal_var=False)
            donor_rows.append({"organ": "heart", "cell_type": tag, "gene": g, "n_disease": len(a), "n_control": len(b),
                               "mean_disease": float(np.mean(a)), "mean_control": float(np.mean(b)), "welch_p": float(p)})
donor_df = pd.DataFrame(donor_rows)
from statsmodels.stats.multitest import multipletests
if "cell_type" not in donor_df.columns:
    donor_df["cell_type"] = ""
donor_df["welch_p"] = pd.to_numeric(donor_df["welch_p"], errors="coerce")
donor_df = donor_df.dropna(subset=["welch_p"])
donor_df["fdr"] = multipletests(donor_df["welch_p"].values.astype(float), method="fdr_bh")[1]
donor_df.to_csv(f"{OUT}/b4_donor_level_tests.csv", index=False)
print(donor_df.sort_values("welch_p").head(12).to_string(index=False))
print("DONE")
