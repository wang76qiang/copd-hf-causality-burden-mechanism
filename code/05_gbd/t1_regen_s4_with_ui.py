# -*- coding: utf-8 -*-
# ============================================================================
# Module 05 - GBD: regenerate the 210-country 2021 ASPR table with 95% UIs.
# Input : GBD_RAW_DIR/IHME-GBD_2021_DATA-*.csv
# Output: results/t1_country_aspr_2021_strict.csv
# ============================================================================

import os as _os
# --- repository paths -----------------------------------------------------------
# REPO resolves to the repository root (code/<module>/this_file.py -> ../..).
# Override with the REPRO_ROOT environment variable if you run from elsewhere.
REPO = _os.environ.get("REPRO_ROOT", _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "..")))
# GBD_RAW_DIR: folder with the raw IHME GBD 2021 csv downloads (see data/README.md)
GBD_RAW_DIR = _os.environ.get("GBD_RAW_DIR", _os.path.join(REPO, "data", "raw", "gbd_2021"))
"""Regenerate S4 country table with 95% UI columns from raw IHME GBD 2021 files."""
import sys, re, glob, os, pandas as pd
sys.stdout.reconfigure(encoding='utf-8')

OUT = REPO
usecols = ['location_name', 'sex_name', 'age_name', 'measure_name', 'metric_name', 'year', 'val', 'upper', 'lower']
frames = []
for f in sorted(glob.glob(os.path.join(GBD_RAW_DIR, "IHME-GBD_2021_DATA-*.csv"))):
    df = pd.read_csv(f, usecols=usecols)
    df = df[(df.measure_name == 'Prevalence') & (df.year == 2021) & (df.sex_name == 'Both')]
    frames.append(df)
d = pd.concat(frames, ignore_index=True)

rate = d[(d.age_name == 'Age-standardized') & (d.metric_name == 'Rate')][
    ['location_name', 'val', 'lower', 'upper']].rename(
    columns={'val': 'ASPR_2021', 'lower': 'ASPR_lower', 'upper': 'ASPR_upper'})
num = d[(d.age_name == 'All ages') & (d.metric_name == 'Number')][
    ['location_name', 'val', 'lower', 'upper']].rename(
    columns={'val': 'cases_2021', 'lower': 'cases_lower', 'upper': 'cases_upper'})
m = rate.merge(num, on='location_name', how='inner')

aggpat = re.compile(r" - WB|Commonwealth|World Bank|SDI|WHO|OECD|G20|BRICS|Sub-Saharan|Latin America|"
                    r"Asia|Europe|Oceania|Caribbean|North Africa|Australasia|Global|High-income|income|"
                    r"Region|Countries|Bank|Middle East|Health System")
s4 = m[~m.location_name.str.contains(aggpat)].copy()
# exclude UN sub-regional aggregates that leak through the regex (2026-07-27 audit:
# these six rows made the count 210 instead of the true 204 GBD countries)
s4 = s4[~s4.location_name.isin(['Africa', 'America', 'Southern Africa', 'Eastern Africa',
                                'Western Africa', 'Central Africa', 'Northern Africa', 'North America'])]
s4 = s4.sort_values('ASPR_2021', ascending=False).reset_index(drop=True)
s4['rank'] = s4.index + 1
s4 = s4[['location_name', 'ASPR_2021', 'ASPR_lower', 'ASPR_upper', 'cases_2021', 'cases_lower', 'cases_upper', 'rank']]
s4.to_csv(os.path.join(OUT, "results", "t1_country_aspr_2021_strict.csv"), index=False)
print('rows:', len(s4))
for n in ['Papua New Guinea', 'China', 'Uzbekistan']:
    r = s4[s4.location_name == n].iloc[0]
    print(f"{n}: ASPR {r.ASPR_2021:.2f} ({r.ASPR_lower:.2f}–{r.ASPR_upper:.2f}), cases {r.cases_2021:,.0f}")
print('fold:', s4.ASPR_2021.max() / s4.ASPR_2021.min())
