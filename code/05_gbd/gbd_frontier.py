# -*- coding: utf-8 -*-
# ============================================================================
# Module 05 - GBD: frontier analysis
#
# Identifies the lowest age-standardized prevalence (ASPR) of COPD-attributable
# HF observed at each level of socio-demographic development, using the full
# 1990-2021 country panel, then quantifies each country's 2021 "effective
# difference" (observed ASPR minus the frontier value at its SDI) - the
# unrealized improvement potential.
#
# Method (following GBD conventions):
#   * the frontier is the lower envelope of the ASPR-vs-SDI cloud across all
#     country-years 1990-2021: countries are sorted by SDI and the cumulative
#     minimum of a robust low quantile (here: the running minimum of the
#     observed ASPR, after removing the most extreme outliers with a 1st
#     percentile Winsorisation per SDI bin) defines the frontier line;
#   * effective difference = ASPR_2021 - frontier(SDI_2021);
#   * the 15 countries with the largest effective differences define the
#     frontier-lagging set reported in the manuscript (incl. China, India,
#     the Netherlands, Canada and Australia).
#
# Inputs:
#   1. Country-year ASPR, recomputed here from the raw IHME files
#      (GBD_RAW_DIR/IHME-GBD_2021_DATA-*.csv).
#   2. Country-year SDI values: data/derived/gbd2021_sdi.csv (NOT SHIPPED -
#      see gbd_inequality.py header and data/README.md for the GHDx source).
#
# Output: results/gbd_frontier.csv  (all countries, ranked by effective
#          difference, with a flag for the top 15)
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
import re as _re
AGG_PAT = _re.compile(r" - WB|Commonwealth|World Bank|SDI|WHO|OECD|G20|BRICS|Sub-Saharan|Latin America|"
                      r"Asia|Europe|Oceania|Caribbean|North Africa|Australasia|Global|High-income|income|"
                      r"Region|Countries|Bank|Middle East|Health System")

def is_country(name):
    return (not AGG_PAT.search(name)) and name not in ("Africa", "America")

USE = ["measure_name", "location_name", "sex_name", "age_name", "metric_name", "year", "val"]


def load_aspr():
    frames = []
    for f in sorted(glob.glob(os.path.join(GBD_RAW_DIR, "IHME-GBD_2021_DATA-*.csv"))):
        for ch in pd.read_csv(f, usecols=USE, chunksize=400_000):
            m = ch[(ch.sex_name == "Both") & (ch.measure_name == PREV) &
                   (ch.age_name == "Age-standardized") & (ch.metric_name == "Rate")]
            frames.append(m[["location_name", "year", "val"]])
    d = pd.concat(frames, ignore_index=True)
    d = d[d.location_name.map(is_country)]
    return d.rename(columns={"val": "ASPR"})


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
    return sdi.rename(columns=ren)[["location_name", "year", "sdi"]]


def frontier_line(sdi, aspr, grid):
    """Lower envelope: running minimum of the 1st-percentile-Winsorised ASPR."""
    order = np.argsort(sdi)
    s, y = sdi[order], aspr[order]
    # Winsorise within small SDI bins to remove downward outliers
    bins = np.floor(s * 50).astype(int)
    yw = y.copy()
    for b in np.unique(bins):
        m = bins == b
        if m.sum() >= 5:
            lo = np.percentile(y[m], 1)
            yw[m] = np.maximum(y[m], lo)
    env = np.minimum.accumulate(yw)
    return np.interp(grid, s, env)


def main():
    if not os.path.exists(SDI_CSV):
        raise SystemExit(
            f"SDI file not found: {SDI_CSV}\n"
            "See this script's header / data/README.md for the GHDx download source.")
    d = load_aspr().merge(load_sdi(), on=["location_name", "year"], how="inner")
    panel = d.dropna(subset=["sdi", "ASPR"])
    grid = np.linspace(panel.sdi.min(), panel.sdi.max(), 400)
    front = frontier_line(panel.sdi.values, panel.ASPR.values.astype(float), grid)

    cur = panel[panel.year == 2021].copy()
    cur["frontier"] = np.interp(cur.sdi, grid, front)
    cur["effective_difference"] = cur.ASPR - cur.frontier
    cur = cur.sort_values("effective_difference", ascending=False).reset_index(drop=True)
    cur["rank"] = cur.index + 1
    cur["frontier_lagging_top15"] = cur["rank"] <= 15
    os.makedirs(OUT_DIR, exist_ok=True)
    cur.to_csv(os.path.join(OUT_DIR, "gbd_frontier.csv"), index=False)
    print(cur.head(15)[["rank", "location_name", "sdi", "ASPR", "frontier",
                        "effective_difference"]].to_string(index=False))


if __name__ == "__main__":
    main()
