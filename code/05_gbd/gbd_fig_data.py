# -*- coding: utf-8 -*-
# ============================================================================
# Module 05 - GBD: regenerate the four figure-data csv files used by
# make_fig2.py / make_fig3.py, directly from the raw IHME GBD 2021 files:
#   results/fig_sdi_series_prev.csv    (year x SDI-quintile ASPR, 1990-2021)
#   results/fig_sdi_series_yld.csv     (year x SDI-quintile ASYR, 1990-2021)
#   results/fig_country_1990_2021.csv  (210 countries: ASPR + cases 1990/2021)
#   results/fig_eapc_prev_regions.csv  (EAPC of ASPR: Global/sexes/SDI/regions)
#
# Input : GBD_RAW_DIR/IHME-GBD_2021_DATA-*.csv   (streamed, see data/README.md)
# Output: the four csv files in OUT_DIR (default results/).
# Run with CHECK=1 to compare freshly computed values against the shipped
# results/ copies without overwriting them (used in VALIDATION.md).
# ============================================================================
import os
import glob
import numpy as np
import pandas as pd
from scipy import stats

REPO = os.environ.get("REPRO_ROOT", os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")))
GBD_RAW_DIR = os.environ.get("GBD_RAW_DIR", os.path.join(REPO, "data", "raw", "gbd_2021"))
OUT_DIR = os.environ.get("OUT_DIR", os.path.join(REPO, "results"))
CHECK = os.environ.get("CHECK", "") == "1"

PREV = "Prevalence"
YLD = "YLDs (Years Lived with Disability)"
SDI_LOCS = ["Global", "High SDI", "High-middle SDI", "Middle SDI", "Low-middle SDI", "Low SDI"]
REGIONS = ["Andean Latin America", "Australasia", "Caribbean", "Central Asia", "Central Europe",
           "Central Latin America", "Central Sub-Saharan Africa", "East Asia", "Eastern Europe",
           "Eastern Sub-Saharan Africa", "High-income Asia Pacific", "High-income North America",
           "North Africa and Middle East", "Oceania", "South Asia", "Southeast Asia",
           "Southern Latin America", "Southern Sub-Saharan Africa", "Tropical Latin America",
           "Western Europe", "Western Sub-Saharan Africa"]
# Aggregate (non-country) location filter, identical to t1_regen_s4_with_ui.py;
# yields exactly the 210 country-level locations used in the manuscript.
import re as _re
AGG_PAT = _re.compile(r" - WB|Commonwealth|World Bank|SDI|WHO|OECD|G20|BRICS|Sub-Saharan|Latin America|"
                      r"Asia|Europe|Oceania|Caribbean|North Africa|Australasia|Global|High-income|income|"
                      r"Region|Countries|Bank|Middle East|Health System")

def is_country(name):
    return (not AGG_PAT.search(name)) and name not in ("Africa", "America")

USE = ["measure_name", "location_name", "sex_name", "age_name", "metric_name", "year", "val"]


def load():
    frames = []
    for f in sorted(glob.glob(os.path.join(GBD_RAW_DIR, "IHME-GBD_2021_DATA-*.csv"))):
        for ch in pd.read_csv(f, usecols=USE, chunksize=400_000):
            m = ch[ch.age_name.isin(["Age-standardized", "All ages"])]
            frames.append(m)
    return pd.concat(frames, ignore_index=True)


def eapc(y, years):
    x = np.asarray(years, float)
    ly = np.log(np.asarray(y, float))
    b, a = np.polyfit(x, ly, 1)
    resid = ly - (a + b * x)
    n = len(x)
    s2 = (resid ** 2).sum() / (n - 2)
    se = np.sqrt(s2 / ((x - x.mean()) ** 2).sum())
    tc = stats.t.ppf(0.975, n - 2)
    return 100 * (np.exp(b) - 1), 100 * (np.exp(b - tc * se) - 1), 100 * (np.exp(b + tc * se) - 1)


def main():
    d = load()
    os.makedirs(OUT_DIR, exist_ok=True)
    asr = d[(d.age_name == "Age-standardized") & (d.metric_name == "Rate")]

    # ---- 1+2: SDI series (both sexes) ----
    for meas, short in [(PREV, "prev"), (YLD, "yld")]:
        s = asr[(asr.measure_name == meas) & (asr.sex_name == "Both") &
                (asr.location_name.isin(SDI_LOCS))]
        wide = s.pivot_table(index="year", columns="location_name", values="val").reset_index()
        wide = wide[["year", "Global", "High SDI", "High-middle SDI", "Low SDI",
                     "Low-middle SDI", "Middle SDI"]]
        wide.to_csv(os.path.join(OUT_DIR, f"fig_sdi_series_{short}.csv"), index=False)

    # ---- 3: country 1990/2021 ----
    nat = asr[(asr.measure_name == PREV) & (asr.sex_name == "Both") &
              (asr.location_name.map(is_country))]
    num = d[(d.measure_name == PREV) & (d.sex_name == "Both") &
            (d.age_name == "All ages") & (d.metric_name == "Number") &
            (d.location_name.map(is_country))]
    r90 = nat[nat.year == 1990][["location_name", "val"]].rename(columns={"val": "ASPR_1990"})
    r21 = nat[nat.year == 2021][["location_name", "val"]].rename(columns={"val": "ASPR_2021"})
    c90 = num[num.year == 1990][["location_name", "val"]].rename(columns={"val": "cases_1990"})
    c21 = num[num.year == 2021][["location_name", "val"]].rename(columns={"val": "cases_2021"})
    country = r90.merge(r21, on="location_name").merge(c90, on="location_name").merge(c21, on="location_name")
    country = country.sort_values("ASPR_2021", ascending=False).reset_index(drop=True)
    country.to_csv(os.path.join(OUT_DIR, "fig_country_1990_2021.csv"), index=False)
    print("countries:", len(country))

    # ---- 4: EAPC table ----
    rows = []
    def add_row(label, sub):
        sub = sub.sort_values("year")
        e, lo, hi = eapc(sub.val.values, sub.year.values)
        rows.append({"Region": label, "Min": lo, "Median": e, "Max": hi})
    for sex, lab in [("Both", "Global"), ("Female", "Female"), ("Male", "Male")]:
        add_row(lab, asr[(asr.measure_name == PREV) & (asr.location_name == "Global") &
                         (asr.sex_name == sex)])
    for loc in SDI_LOCS[1:] + REGIONS:
        add_row(loc, asr[(asr.measure_name == PREV) & (asr.location_name == loc) &
                         (asr.sex_name == "Both")])
    eapc_df = pd.DataFrame(rows)
    eapc_df.to_csv(os.path.join(OUT_DIR, "fig_eapc_prev_regions.csv"), index=False)

    # ---- optional check against shipped copies ----
    if CHECK:
        for fn in ["fig_sdi_series_prev.csv", "fig_sdi_series_yld.csv",
                   "fig_country_1990_2021.csv", "fig_eapc_prev_regions.csv"]:
            new = pd.read_csv(os.path.join(OUT_DIR, fn))
            old = pd.read_csv(os.path.join(REPO, "results", fn))
            numcols = new.select_dtypes("number").columns
            key = [c for c in old.columns if c in ("year", "location_name", "Region")]
            merged = old.merge(new, on=key, suffixes=("_old", "_new"),
                               how="outer", indicator=True)
            both = merged[merged._merge == "both"]
            maxdiff = 0.0
            for c in numcols:
                if c == "year":
                    continue
                dd = (both[f"{c}_old"] - both[f"{c}_new"]).abs().max()
                maxdiff = max(maxdiff, dd)
            print(f"CHECK {fn}: rows old={len(old)} new={len(new)} "
                  f"shared_keys={len(both)} max_abs_numeric_diff={maxdiff:.3g}")


if __name__ == "__main__":
    main()
