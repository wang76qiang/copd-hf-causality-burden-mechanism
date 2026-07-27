# -*- coding: utf-8 -*-
# ============================================================================
# Module 05 - GBD: cross-country health inequality (SII and CI)
#
# Slope index of inequality (SII) and concentration index (CI) of the
# age-standardized prevalence (ASPR) and YLD (ASYR) rates of COPD-attributable
# HF across 210 countries in 1990 and 2021, following WHO/GBD conventions:
#   * countries are ranked by SDI; each country's relative rank (ridit score)
#     is the midpoint of its cumulative population share;
#   * SII = slope of the population-weighted linear regression of the rate on
#     the ridit score;
#   * CI  = (2 / mu) * sum_i p_i * rate_i * R_i - 1, with p_i the population
#     share, R_i the ridit score and mu the population-weighted mean rate.
#
# Inputs:
#   1. Country rates + implied populations, recomputed here from the raw IHME
#      files (GBD_RAW_DIR/IHME-GBD_2021_DATA-*.csv; both-sex, all-ages number
#      and crude rate give the implied population; rates are the
#      age-standardized ones).
#   2. Country-year SDI values: data/derived/gbd2021_sdi.csv with columns
#      location_name, year, sdi  -- NOT SHIPPED.  Download "GBD 2021
#      Socio-Demographic Index (SDI) 1950-2021" from IHME GHDx
#      (https://ghdx.healthdata.org/record/global-burden-disease-study-2021-
#       gbd-2021-socio-demographic-index-sdi-1950%E2%80%932021; free account
#      required) and reformat to this long csv. GBD 2021 SDI values cannot be
#      redistributed and are not present in the raw burden downloads.
#
# Output: results/gbd_inequality.csv
# Manuscript targets (Table 3 / Fig. 4c-f): SII (ASPR) -40.05 (1990) ->
# -23.10 (2021); SII (ASYR) -3.52 -> -2.03; CI -0.25 -> -0.15.
# ============================================================================
import os
import glob
import numpy as np
import pandas as pd

REPO = os.environ.get("REPRO_ROOT", os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")))
GBD_RAW_DIR = os.environ.get("GBD_RAW_DIR", os.path.join(REPO, "data", "raw", "gbd_2021"))
OUT_DIR = os.environ.get("OUT_DIR", os.path.join(REPO, "results"))
SDI_CSV = os.environ.get("SDI_CSV", os.path.join(REPO, "data", "derived", "gbd2021_sdi.csv"))

PREV = "Prevalence"
YLD = "YLDs (Years Lived with Disability)"
import re as _re
AGG_PAT = _re.compile(r" - WB|Commonwealth|World Bank|SDI|WHO|OECD|G20|BRICS|Sub-Saharan|Latin America|"
                      r"Asia|Europe|Oceania|Caribbean|North Africa|Australasia|Global|High-income|income|"
                      r"Region|Countries|Bank|Middle East|Health System")

def is_country(name):
    return (not AGG_PAT.search(name)) and name not in ("Africa", "America")

USE = ["measure_name", "location_name", "sex_name", "age_name", "metric_name", "year", "val"]


def load_burden():
    frames = []
    for f in sorted(glob.glob(os.path.join(GBD_RAW_DIR, "IHME-GBD_2021_DATA-*.csv"))):
        for ch in pd.read_csv(f, usecols=USE, chunksize=400_000):
            m = ch[(ch.sex_name == "Both") & (ch.year.isin([1990, 2021])) &
                   (ch.age_name.isin(["Age-standardized", "All ages"]))]
            frames.append(m)
    d = pd.concat(frames, ignore_index=True)
    d = d[d.location_name.map(is_country)]
    asr = d[(d.age_name == "Age-standardized") & (d.metric_name == "Rate")]
    aa = d[(d.age_name == "All ages") & (d.metric_name.isin(["Number", "Rate"]))]
    # implied population from all-ages number / crude rate (prevalence)
    n = aa[(aa.measure_name == PREV) & (aa.metric_name == "Number")][
        ["location_name", "year", "val"]].rename(columns={"val": "num"})
    r = aa[(aa.measure_name == PREV) & (aa.metric_name == "Rate")][
        ["location_name", "year", "val"]].rename(columns={"val": "crude"})
    pop = n.merge(r, on=["location_name", "year"])
    pop["population"] = pop.num / pop.crude * 1e5
    out = asr.pivot_table(index=["location_name", "year"], columns="measure_name",
                          values="val").reset_index()
    out = out.rename(columns={PREV: "ASPR", YLD: "ASYR"})
    return out.merge(pop[["location_name", "year", "population"]],
                     on=["location_name", "year"])


def load_sdi():
    sdi = pd.read_csv(SDI_CSV)
    sdi.columns = [c.strip().lower().replace(" ", "_") for c in sdi.columns]
    ren = {}
    for c in sdi.columns:
        if c in ("location", "location_name", "country"):
            ren[c] = "location_name"
        elif c in ("year", "year_id"):
            ren[c] = "year"
        elif c in ("sdi", "sdi_value", "val", "mean_value"):
            ren[c] = "sdi"
    sdi = sdi.rename(columns=ren)
    return sdi[["location_name", "year", "sdi"]]


def sii_ci(df, rate_col):
    """Population-weighted SII and CI for one year-slice of country data."""
    df = df.dropna(subset=["sdi", rate_col, "population"]).sort_values("sdi").reset_index(drop=True)
    pop = df.population.values
    y = df[rate_col].values.astype(float)
    p = pop / pop.sum()
    cum = np.concatenate([[0.0], np.cumsum(p)])
    ridit = (cum[:-1] + cum[1:]) / 2.0
    # SII: weighted least squares of y on ridit
    X = np.column_stack([np.ones_like(ridit), ridit])
    W = p
    XtW = X.T * W
    slope = np.linalg.solve(XtW @ X, XtW @ y)[1]
    # CI
    mu = (p * y).sum()
    ci = 2.0 * (p * y * ridit).sum() / mu - 1.0
    return slope, ci


def main():
    if not os.path.exists(SDI_CSV):
        raise SystemExit(
            f"SDI file not found: {SDI_CSV}\n"
            "Download 'GBD 2021 Socio-Demographic Index (SDI) 1950-2021' from IHME GHDx\n"
            "(free account required; URL in this script's header and data/README.md),\n"
            "reformat to columns location_name,year,sdi and place it there.")
    burden = load_burden()
    sdi = load_sdi()
    d = burden.merge(sdi, on=["location_name", "year"], how="left")
    missing = d[d.sdi.isna()].location_name.unique()
    if len(missing):
        print("WARNING: no SDI for:", ", ".join(sorted(missing)[:10]),
              f"({len(missing)} locations dropped)")
    rows = []
    for yr in [1990, 2021]:
        for rate_col in ["ASPR", "ASYR"]:
            s, c = sii_ci(d[d.year == yr], rate_col)
            rows.append({"year": yr, "measure": rate_col, "SII": s, "CI": c})
            print(f"{yr} {rate_col}: SII {s:.2f}  CI {c:.3f}")
    out = pd.DataFrame(rows)
    os.makedirs(OUT_DIR, exist_ok=True)
    out.to_csv(os.path.join(OUT_DIR, "gbd_inequality.csv"), index=False)


if __name__ == "__main__":
    main()
