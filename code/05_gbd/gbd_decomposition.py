# -*- coding: utf-8 -*-
# ============================================================================
# Module 05 - GBD: Das Gupta decomposition of case growth, 1990 -> 2021
#
# Decomposes the change in COPD-attributable HF case numbers into three
# factors - population growth, population ageing (age-structure change) and
# epidemiological change (age-specific rates) - using Das Gupta's symmetric
# factor decomposition (averaging over all factor orderings), for the Global
# aggregate and the five SDI quintiles, for both Prevalence and YLDs.
#
# Population by age group is derived as cases/rate*1e5 from the raw IHME
# files; age groups with structurally zero cases (rate == 0, e.g. <15 years)
# are pooled into a single residual group whose population is the difference
# between the all-ages implied population and the sum of the case-carrying
# groups (this is required to capture the ageing effect of the declining
# child share).
#
# Input : GBD_RAW_DIR/IHME-GBD_2021_DATA-*.csv   (streamed, see data/README.md)
# Output: results/gbd_decomposition.csv
#
# Manuscript targets (Table 3 / Fig. 4a,b; verified in VALIDATION.md):
#   Global      Prevalence: ageing +53.3%, growth +44.7%, epi +1.9%
#   Low SDI     Prevalence: ageing +8.6%,  growth +108.9%, epi -17.5%
#   Middle SDI  YLDs      : ageing +85.6%, growth +41.9%,  epi -27.5%
#   High SDI    (Fig. 4b prints growth/ageing with flipped signs; this script
#               computes growth +20.7%, ageing +44.8%, epi +34.5% - the epi
#               share matches; see VALIDATION.md for the discrepancy note)
# ============================================================================
import os
import glob
import numpy as np
import pandas as pd

REPO = os.environ.get("REPRO_ROOT", os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")))
GBD_RAW_DIR = os.environ.get("GBD_RAW_DIR", os.path.join(REPO, "data", "raw", "gbd_2021"))
OUT_DIR = os.environ.get("OUT_DIR", os.path.join(REPO, "results"))

AGE5 = ['<5 years', '5-9 years', '10-14 years', '15-19 years', '20-24 years', '25-29 years',
        '30-34 years', '35-39 years', '40-44 years', '45-49 years', '50-54 years', '55-59 years',
        '60-64 years', '65-69 years', '70-74 years', '75-79 years', '80-84 years', '85-89 years',
        '90-94 years', '95+ years']
LOCS = ["Global", "High SDI", "High-middle SDI", "Middle SDI", "Low-middle SDI", "Low SDI"]
MEASURES = ["Prevalence", "YLDs (Years Lived with Disability)"]
USE = ["measure_name", "location_name", "sex_name", "age_name", "metric_name", "year", "val"]


def load():
    frames = []
    for f in sorted(glob.glob(os.path.join(GBD_RAW_DIR, "IHME-GBD_2021_DATA-*.csv"))):
        for ch in pd.read_csv(f, usecols=USE, chunksize=400_000):
            m = ch[(ch.sex_name == "Both") & (ch.location_name.isin(LOCS)) &
                   (ch.year.isin([1990, 2021]))]
            frames.append(m)
    return pd.concat(frames, ignore_index=True)


def get_series(sub, yr):
    """Population and rate vectors over case-carrying ages + residual (zero-rate) group."""
    aa_n = sub[(sub.year == yr) & (sub.age_name == "All ages") & (sub.metric_name == "Number")].val
    aa_r = sub[(sub.year == yr) & (sub.age_name == "All ages") & (sub.metric_name == "Rate")].val
    if len(aa_n) == 0 or len(aa_r) == 0:
        raise ValueError("all-ages rows missing")
    p_tot = aa_n.iloc[0] / aa_r.iloc[0] * 1e5
    n = (sub[(sub.year == yr) & (sub.age_name.isin(AGE5)) & (sub.metric_name == "Number")]
         .set_index("age_name").reindex(AGE5)["val"].values.astype(float))
    r = (sub[(sub.year == yr) & (sub.age_name.isin(AGE5)) & (sub.metric_name == "Rate")]
         .set_index("age_name").reindex(AGE5)["val"].values.astype(float))
    with np.errstate(divide="ignore", invalid="ignore"):
        p = np.where(r > 0, n / r * 1e5, 0.0)
    keep = r > 0
    # residual zero-rate population group (children etc.), rate = 0
    p_vec = np.concatenate([p[keep], [p_tot - p[keep].sum()]])
    r_vec = np.concatenate([r[keep], [0.0]])
    return p_vec, r_vec


def das_gupta(p1, r1, p2, r2):
    """Symmetric 3-factor Das Gupta decomposition of cases = P * sum(s_a * r_a)."""
    P1, P2 = p1.sum(), p2.sum()
    s1, s2 = p1 / P1, p2 / P2
    F11 = (s1 * r1).sum(); F22 = (s2 * r2).sum()
    F21 = (s2 * r1).sum(); F12 = (s1 * r2).sum()
    dcases = P2 * F22 - P1 * F11
    e_growth = (P2 - P1) / 6 * (2 * (F11 + F22) + F21 + F12)
    e_aging = ((s2 - s1) / 6 * (2 * (P1 * r1 + P2 * r2) + (P2 * r1 + P1 * r2))).sum()
    e_epi = ((r2 - r1) / 6 * (2 * (P1 * s1 + P2 * s2) + (P2 * s1 + P1 * s2))).sum()
    return dcases, e_growth, e_aging, e_epi


def main():
    d = load()
    rows = []
    for loc in LOCS:
        for meas in MEASURES:
            sub = d[(d.location_name == loc) & (d.measure_name == meas)]
            p1, r1 = get_series(sub, 1990)
            p2, r2 = get_series(sub, 2021)
            dc, g, a, e = das_gupta(p1, r1, p2, r2)
            rows.append({
                "location": loc, "measure": meas,
                "cases_1990": (p1 * r1).sum() / 1e5, "cases_2021": (p2 * r2).sum() / 1e5,
                "delta_cases": dc / 1e5,
                "pop_growth_pct": g / dc * 100,
                "ageing_pct": a / dc * 100,
                "epi_change_pct": e / dc * 100,
            })
    out = pd.DataFrame(rows)
    os.makedirs(OUT_DIR, exist_ok=True)
    out.to_csv(os.path.join(OUT_DIR, "gbd_decomposition.csv"), index=False)
    pd.set_option("display.width", 160)
    print(out.to_string(index=False,
                        float_format=lambda x: f"{x:,.2f}"))


if __name__ == "__main__":
    main()
