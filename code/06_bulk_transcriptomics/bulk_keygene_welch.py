# -*- coding: utf-8 -*-
# ============================================================================
# Module 06 - Bulk transcriptomics: key-gene Welch t-tests + normality checks
#
# Recomputes the two-group key-gene expression tests quoted in the manuscript
# from the investigator-curated expression workbooks shipped in
# data/derived/bulk/:
#   GSE57148 (COPD lung, 95 COPD vs 89 control): S100A8, S100A12, IL1RL1, THBS1
#   GSE57338 (HF myocardium, 177 failing vs 136 non-failing): LUM, ASPN, SMOC2,
#     ACE, IL1RL1
# Welch's unequal-variance t-test + Shapiro-Wilk normality per group.
#
# Input : data/derived/bulk/GSE57148_GSE57338_keygene_expression_values.xlsx
# Output: results/bulk_keygene_welch.csv
# Manuscript targets (see VALIDATION.md): S100A8 P=8e-6; S100A12 P=0.02;
# IL1RL1 P=0.03; THBS1 P=0.001; heart genes P <= ~2e-32.
# ============================================================================
import os
import numpy as np
import pandas as pd
from scipy import stats

REPO = os.environ.get("REPRO_ROOT", os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")))
XL = os.environ.get("KEYGENE_XLSX", os.path.join(
    REPO, "data", "derived", "bulk", "GSE57148_GSE57338_keygene_expression_values.xlsx"))
OUT_DIR = os.environ.get("OUT_DIR", os.path.join(REPO, "results"))

# gene -> first column of its Sample/Title/Value block (0-based) in each sheet
LUNG = {"sheet": "GSE57148", "case_label": "COPD", "control_label": "control",
        "genes": {"S100A8": 1, "S100A12": 5, "IL1RL1": 9, "THBS1": 13}}
HEART = {"sheet": "GSE57338",
         "genes": {"LUM": 13, "ASPN": 17, "SMOC2": 21, "ACE": 25, "IL1RL1": 9}}


def lung_groups(raw, c0):
    sub = raw.iloc[2:, [c0, c0 + 1, c0 + 2]].dropna(how="all")
    sub.columns = ["sample", "title", "val"]
    sub = sub.dropna(subset=["val"])
    if sub["title"].isna().all():
        # THBS1 block lost its Title column; recover groups via the S100A8 block
        smap = dict(zip(raw.iloc[2:, 1], raw.iloc[2:, 2]))
        sub["title"] = sub["sample"].map(smap)
    case = sub[sub.title == "COPD"]["val"].astype(float).values
    ctrl = sub[sub.title == "control"]["val"].astype(float).values
    return case, ctrl


def heart_groups(raw, c0):
    sub = raw.iloc[2:, [c0 + 1, c0 + 2]].dropna()
    sub.columns = ["title", "val"]
    t = sub["title"].astype(str)
    is_ctrl = t.str.contains("Non-failing") | (t.str.lower() == "control")
    ctrl = sub[is_ctrl]["val"].astype(float).values
    case = sub[~is_ctrl]["val"].astype(float).values
    return case, ctrl


def welch(dataset, gene, case, ctrl):
    t, p = stats.ttest_ind(case, ctrl, equal_var=False)
    sw_case = stats.shapiro(case)[1] if 3 <= len(case) <= 5000 else np.nan
    sw_ctrl = stats.shapiro(ctrl)[1] if 3 <= len(ctrl) <= 5000 else np.nan
    return {"dataset": dataset, "gene": gene,
            "n_case": len(case), "n_control": len(ctrl),
            "mean_case": case.mean(), "mean_control": ctrl.mean(),
            "log2FC_proxy": np.log2(case.mean() / ctrl.mean()),
            "welch_t": t, "welch_p": p,
            "shapiro_p_case": sw_case, "shapiro_p_control": sw_ctrl}


def main():
    rows = []
    raw = pd.read_excel(XL, sheet_name=LUNG["sheet"], header=None)
    for g, c0 in LUNG["genes"].items():
        case, ctrl = lung_groups(raw, c0)
        rows.append(welch("GSE57148", g, case, ctrl))
    raw = pd.read_excel(XL, sheet_name=HEART["sheet"], header=None)
    for g, c0 in HEART["genes"].items():
        case, ctrl = heart_groups(raw, c0)
        rows.append(welch("GSE57338", g, case, ctrl))
    out = pd.DataFrame(rows)
    os.makedirs(OUT_DIR, exist_ok=True)
    out.to_csv(os.path.join(OUT_DIR, "bulk_keygene_welch.csv"), index=False)
    print(out[["dataset", "gene", "n_case", "n_control", "mean_case",
               "mean_control", "welch_p"]].to_string(index=False))


if __name__ == "__main__":
    main()
