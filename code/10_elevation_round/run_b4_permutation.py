# -*- coding: utf-8 -*-
"""M4: permutation test for interorgan edges.
Step 1: cache per-cell panel (log1p-CP10k) + obs to npz (one-time streaming).
Step 2: observed edges (as in B4) + 1000 donor-label permutations per organ -> empirical P per pair.
"""
import numpy as np, pandas as pd, h5py, os, warnings, sys
warnings.filterwarnings("ignore")
sys.path.insert(0, ".")
from scipy.stats import ttest_ind

OUT = "results_new"; os.makedirs(OUT, exist_ok=True)
BLOCK = 4000
rng = np.random.default_rng(20260904)
NPERM = 1000

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
    with h5py.File(path, "r") as f:
        obs = {}
        for col in [disease_col, "donor_id", "cell_type"]:
            g = f["obs"][col]
            if isinstance(g, h5py.Group):
                cats = read_str_list(g["categories"]); codes = g["codes"][:]
                obs[col] = [cats[c] if c >= 0 else "" for c in codes]
            else:
                obs[col] = read_str_list(g)
        n = len(obs[disease_col])
        var = f["var"]
        fn = var["feature_name"]
        names = read_str_list(fn) if isinstance(fn, h5py.Dataset) else read_str_list(fn["categories"])
        gidx = {g: names.index(g) for g in genes if g in names}
        layer = "raw" if "raw" in f else "X"
        X = f[layer]["X"]; data, indices = X["data"], X["indices"]
        indptr = X["indptr"][:]
        totals = np.zeros(n)
        panel = {g: np.zeros(n) for g in gidx}
        gi_set = {v: k for k, v in gidx.items()}
        gicol = np.array(list(gi_set.keys()))
        for s in range(0, n, BLOCK):
            e = min(s + BLOCK, n)
            p0, p1 = int(indptr[s]), int(indptr[e])
            d = data[p0:p1]; ii = indices[p0:p1]
            row_ptr = indptr[s:e+1] - p0
            counts = np.diff(row_ptr)
            sums = np.zeros(e - s); np.add.at(sums, np.repeat(np.arange(e - s), counts), d)
            totals[s:e] = sums
            keep = np.isin(ii, gicol)
            dk, ik = d[keep], ii[keep]
            rb = np.repeat(np.arange(e - s), counts)[keep]
            for gcol, gname in gi_set.items():
                m = ik == gcol
                np.add.at(panel[gname], s + rb[m], dk[m])
        sf = 1e4 / np.maximum(totals, 1)
        M = np.column_stack([np.log1p(panel[g] * sf) for g in gidx])
        df = pd.DataFrame({"donor": obs["donor_id"],
                           "cond": [disease_map.get(x) for x in obs[disease_col]],
                           "ct": obs["cell_type"]})
        return df, M, list(gidx.keys())

def cache_file(tag): return f"cache_{tag}.npz"

def get_panel(tag, path, genes, dmap):
    if os.path.exists(cache_file(tag)):
        z = np.load(cache_file(tag), allow_pickle=True)
        df = pd.DataFrame({"donor": z["donor"], "cond": z["cond"], "ct": z["ct"]})
        return df, z["M"], list(z["genes"])
    df, M, genes = load_panel(path, genes, dmap)
    np.savez_compressed(cache_file(tag), donor=df["donor"].values, cond=df["cond"].values,
                        ct=df["ct"].values, M=M, genes=np.array(genes))
    return df, M, genes

print("=== caching panels ===")
lung_df, lung_M, lung_genes = get_panel("lung", "sc_data/t1_copd_lung.h5ad", LIGANDS + RECEPTORS,
    {"chronic obstructive pulmonary disease": "COPD", "normal": "Normal"})
heart_maps = {"dilated cardiomyopathy": "COPD", "arrhythmogenic right ventricular cardiomyopathy": "COPD",
              "non-compaction cardiomyopathy": "COPD", "normal": "Normal"}
heart = {}
for tag, path in [("fibro", "sc_data/t1_hf_fibroblasts.h5ad"),
                  ("endo", "sc_data/t1_hf_endothelial.h5ad"),
                  ("macro", "sc_data/t1_hf_macrophages.h5ad")]:
    df, M, genes = get_panel(tag, path, LIGANDS + RECEPTORS, heart_maps)
    heart[tag] = (df, M, genes)
    print(tag, df.shape, len(genes), "genes")

# ---- group means from donor assignment ----
def group_means(df, M, genes, donor_cond):
    """donor_cond: dict donor->cond label. Returns dict[(gene, ct, cond)] = cell-weighted mean."""
    valid = df["donor"].map(donor_cond).notna().values
    d2 = df[valid]; M2 = M[valid]
    conds = d2["donor"].map(donor_cond).values
    out = {}
    df2 = pd.DataFrame({"ct": d2["ct"].values, "cond": conds, "donor": d2["donor"].values})
    for (ct, cond), idx in df2.groupby(["ct", "cond"]).groups.items():
        ii = df2.index.get_indexer(idx)
        out[(ct, cond)] = (M2[ii].mean(axis=0), len(ii))
    return out, genes

def pair_deltas(lung_assign, heart_assigns, lung_cts):
    """Returns dict pair -> mean delta_log2 across senders/receivers (edges with min score filter)."""
    lm, lg = group_means(lung_df, lung_M, lung_genes, lung_assign)
    deltas = {}
    hm = {}
    for tag, (df, M, genes) in heart.items():
        hm[tag] = (group_means(df, M, genes, heart_assigns[tag])[0], genes)
    for lig, rec in PAIRS:
        if lig not in lung_genes: continue
        li = lung_genes.index(lig)
        vals = []
        for sct in lung_cts:
            a = lm.get((sct, "COPD")); b = lm.get((sct, "Normal"))
            if not a or not b: continue
            ld, ln = a[0][li], b[0][li]
            for tag, (hmeans, hgenes) in hm.items():
                if rec not in hgenes: continue
                ri = hgenes.index(rec)
                ra = hmeans.get((tag + "_ct", "COPD")); rb_ = hmeans.get((tag + "_ct", "Normal"))
                # heart ct label: we set cell_type tag below
                if ra is None:
                    # heart files: ct column is original; we need unified tag - handle: use any ct key
                    keys_d = [k for k in hmeans if k[1] == "COPD"]; keys_n = [k for k in hmeans if k[1] == "Normal"]
                    if not keys_d or not keys_n: continue
                    rd = np.mean([hmeans[k][0][ri] for k in keys_d]); rn = np.mean([hmeans[k][0][ri] for k in keys_n])
                else:
                    rd, rn = ra[0][ri], rb_[0][ri]
                if max(ld*rd, ln*rn) < 0.005: continue
                vals.append(np.log2((ld*rd + 1e-6) / (ln*rn + 1e-6)))
        if vals:
            deltas[f"{lig}->{rec}"] = float(np.mean(vals))
    return deltas

# NOTE: for heart files we keep original cell_type column; simplify: treat each heart file as single receiver
# To match B4 convention, overwrite heart ct with tag before group means
for tag, (df, M, genes) in heart.items():
    df["ct"] = tag + "_ct"

lung_cts = [c for c in lung_df["ct"].value_counts().index if (lung_df["ct"] == c).sum() >= 100]

def donor_assign(df):
    donors = df.loc[df["cond"].notna(), "donor"].unique()
    cond_of = df.dropna(subset=["cond"]).groupby("donor")["cond"].first().to_dict()
    return cond_of

lung_cond_of = donor_assign(lung_df)
heart_cond_of = {tag: donor_assign(df) for tag, (df, M, g) in heart.items()}

print("=== observed ===")
obs = pair_deltas(lung_cond_of, heart_cond_of, lung_cts)
print(len(obs), "pairs with edges")

print("=== permutations ===")
pairs = list(obs.keys())
perm_stats = {p: [] for p in pairs}
for b in range(NPERM):
    la = {}
    donors = list(lung_cond_of.keys()); conds = [lung_cond_of[d] for d in donors]
    perm = rng.permutation(conds)
    la = dict(zip(donors, perm))
    ha = {}
    for tag, co in heart_cond_of.items():
        donors = list(co.keys()); conds = [co[d] for d in donors]
        ha[tag] = dict(zip(donors, rng.permutation(conds)))
    pd_ = pair_deltas(la, ha, lung_cts)
    for p in pairs:
        if p in pd_:
            perm_stats[p].append(pd_[p])
    if (b + 1) % 100 == 0:
        print(f"  perm {b+1}/{NPERM}", flush=True)

rows = []
for p in pairs:
    perms = np.array(perm_stats[p])
    obsv = obs[p]
    p_emp = (1 + np.sum(np.abs(perms) >= abs(obsv))) / (len(perms) + 1)
    rows.append({"pair": p, "observed_delta_log2": obsv, "n_perm": len(perms),
                 "perm_mean": float(perms.mean()), "perm_sd": float(perms.std()),
                 "empirical_p": p_emp})
res = pd.DataFrame(rows).sort_values("empirical_p")
from statsmodels.stats.multitest import multipletests
res["fdr"] = multipletests(res["empirical_p"].values, method="fdr_bh")[1]
res.to_csv(f"{OUT}/b4_edge_permutation_tests.csv", index=False)
print(res.head(15).to_string(index=False))
print("DONE")
