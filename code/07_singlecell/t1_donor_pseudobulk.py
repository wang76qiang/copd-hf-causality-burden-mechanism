# -*- coding: utf-8 -*-
# ============================================================================
# Module 07 - Single-cell: donor-level pseudobulk validation (= Table S8)
# Donor-level pseudobulk replication of the cell-level DE signals (raw h5py
# streaming, Ensembl -> symbol via feature_name; CSR/CSC safe).
# Input : data/*.h5ad (NOT shipped, see data/README.md)
# Output: results/t1_donor_pseudobulk_validation.csv (= tables/Table_S8)
# ============================================================================

import os as _os
# --- repository paths -----------------------------------------------------------
# REPO resolves to the repository root (code/<module>/this_file.py -> ../..).
# Override with the REPRO_ROOT environment variable if you run from elsewhere.
REPO = _os.environ.get("REPRO_ROOT", _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "..")))
# GBD_RAW_DIR: folder with the raw IHME GBD 2021 csv downloads (see data/README.md)
GBD_RAW_DIR = _os.environ.get("GBD_RAW_DIR", _os.path.join(REPO, "data", "raw", "gbd_2021"))
"""Donor-level pseudobulk validation v4 (raw h5py, Ensembl->symbol via feature_name, CSR/CSC safe)."""
import sys, numpy as np, pandas as pd, h5py, anndata as ad
from scipy import stats
sys.stdout.reconfigure(encoding='utf-8')

BASE = REPO
CHUNK = 200_000

class SparseReader:
    def __init__(self, path):
        self.f = h5py.File(path, 'r')
        X = self.f['X']
        enc = X.attrs.get('encoding-type', b'')
        self.enc = enc.decode() if isinstance(enc, bytes) else str(enc)
        self.data = X['data']; self.indices = X['indices']; self.indptr = X['indptr'][:]
        self.shape = tuple(int(x) for x in X.attrs['shape'])
        print(path.split('/')[-1], self.enc, self.shape, 'nnz:', self.data.shape[0])

    def _csr_row_segments(self):
        """Row boundaries so each segment spans <= CHUNK nonzeros."""
        ip = self.indptr
        nnz = int(ip[-1])
        bounds = [0]
        pos = 0
        while pos < nnz:
            pos = min(pos + CHUNK, nnz)
            bounds.append(int(np.searchsorted(ip, pos, side='left')))
        bounds = sorted(set(bounds + [self.shape[0]]))
        return [(bounds[k], bounds[k + 1]) for k in range(len(bounds) - 1)]

    def lib_sizes(self):
        n_obs, n_vars = self.shape
        lib = np.zeros(n_obs, dtype=np.float64)
        if self.enc == 'csc_matrix':
            N = self.data.shape[0]
            for s in range(0, N, CHUNK):
                e = min(N, s + CHUNK)
                np.add.at(lib, np.asarray(self.indices[s:e]), np.asarray(self.data[s:e], dtype=np.float64))
        else:  # csr: row chunks bounded by nnz
            ip = self.indptr
            for i0, i1 in self._csr_row_segments():
                s, e = int(ip[i0]), int(ip[i1])
                d = np.asarray(self.data[s:e], dtype=np.float64)
                counts = np.diff(ip[i0:i1 + 1]).astype(np.int64)
                keep = counts > 0
                starts = np.concatenate([[0], np.cumsum(counts)[:-1]])[keep]
                lib[i0:i1][keep] = np.add.reduceat(d, starts)
        return np.where(lib == 0, 1, lib)

    def gene_col(self, j):
        n_obs = self.shape[0]
        out = np.zeros(n_obs, dtype=np.float64)
        if self.enc == 'csc_matrix':
            s, e = int(self.indptr[j]), int(self.indptr[j + 1])
            out[np.asarray(self.indices[s:e])] = np.asarray(self.data[s:e], dtype=np.float64)
        else:
            ip = self.indptr
            for i0, i1 in self._csr_row_segments():
                s, e = int(ip[i0]), int(ip[i1])
                idx = np.asarray(self.indices[s:e])
                d = np.asarray(self.data[s:e], dtype=np.float64)
                counts = np.diff(ip[i0:i1 + 1]).astype(np.int64)
                rows = np.repeat(np.arange(i0, i1), counts)
                m = idx == j
                if m.any():
                    np.add.at(out, rows[m], d[m])
        return out

    def close(self):
        self.f.close()

def donor_test(obs, expr, cell_col, cell_type, case_label, ctrl_label, gene):
    mask = (obs[cell_col] == cell_type).values
    df = pd.DataFrame({'donor': obs['donor_id'].values[mask],
                       'cond': obs['disease'].values[mask],
                       'expr': expr[mask]})
    per = df.groupby(['donor', 'cond'], observed=True)['expr'].mean().reset_index()
    case = per.loc[per.cond == case_label, 'expr'].values
    ctrl = per.loc[per.cond == ctrl_label, 'expr'].values
    if len(case) < 2 or len(ctrl) < 2:
        print('  skip (few donors):', cell_type, gene, len(case), len(ctrl)); return None
    t, p = stats.ttest_ind(case, ctrl, equal_var=False)
    return dict(cell_type=cell_type, gene=gene,
                n_case_donors=len(case), n_ctrl_donors=len(ctrl),
                case_mean=round(float(case.mean()), 3), ctrl_mean=round(float(ctrl.mean()), 3),
                logFC_donor=round(float(case.mean() - ctrl.mean()), 3), welch_p=p)

rows = []

def run(path, tests, case_pred, tag):
    a = ad.read_h5ad(path, backed='r')
    obs = a.obs.copy()
    var = a.var.copy()
    a.file.close()
    if 'feature_name' in var.columns:
        sym2idx = {s: i for i, s in enumerate(var['feature_name'])}
    else:
        sym2idx = {s: i for i, s in enumerate(var.index)}
    dv = list(obs['disease'].unique())
    case_lab = [x for x in dv if case_pred(x)][0]
    ctrl_lab = [x for x in dv if x != case_lab][0]
    print(tag, '| case:', case_lab, '| ctrl:', ctrl_lab)
    sr = SparseReader(path)
    lib = sr.lib_sizes()
    for ct, g in tests:
        if g not in sym2idx:
            print('missing gene', g); continue
        expr = np.log1p(sr.gene_col(sym2idx[g]) / lib * 1e4)
        r = donor_test(obs, expr, 'cell_type', ct, case_lab, ctrl_lab, g)
        if r:
            r['dataset'] = tag; rows.append(r); print(r)
    sr.close()
    if rows:
        pd.DataFrame(rows).to_csv(f"{BASE}/results/t1_donor_pseudobulk_validation.csv", index=False)

run(f"{BASE}/data/t1_copd_lung.h5ad",
    [('macrophage', 'THBS1'), ('macrophage', 'CD163'), ('macrophage', 'SERPINE1'),
     ('monocyte', 'S100A8'), ('monocyte', 'S100A12'), ('monocyte', 'F13A1'),
     ('fibroblast', 'SERPINE1'), ('capillary endothelial cell', 'SERPINE1')],
    lambda x: 'chronic' in x.lower() or 'copd' in x.lower(), 'COPD lung')

run(f"{BASE}/data/t1_hf_fibroblasts.h5ad",
    [('fibroblast of cardiac tissue', 'SERPINE1'), ('fibroblast of cardiac tissue', 'ITGAV'), ('fibroblast of cardiac tissue', 'THBS1')],
    lambda x: 'cardiomyopathy' in x.lower(), 'failing heart (fibroblast)')

run(f"{BASE}/data/t1_hf_endothelial.h5ad",
    [('endothelial cell', 'SERPINE1'), ('endothelial cell', 'ITGAV'), ('endothelial cell', 'THBS1')],
    lambda x: 'cardiomyopathy' in x.lower(), 'failing heart (endothelial)')

out = pd.DataFrame(rows)
out['bonf_p'] = np.minimum(out.welch_p * len(out), 1.0)
out.to_csv(f"{BASE}/results/t1_donor_pseudobulk_validation.csv", index=False)
print('saved', len(out), 'tests')
