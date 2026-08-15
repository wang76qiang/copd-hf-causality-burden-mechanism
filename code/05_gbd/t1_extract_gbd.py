# -*- coding: utf-8 -*-
# ============================================================================
# Module 05 - GBD 2021 burden (master extraction)
# Audit anchors, age panel, ASR series, EAPC + COVID sensitivity, country
# ranking - all computed from the raw IHME csv files (streamed).
# Input : GBD_RAW_DIR/IHME-GBD_2021_DATA-*.csv  (see data/README.md)
# Output: results/t1_extract_summary.json, t1_bapc_agepanel.csv,
#         t1_asr_series.csv, t1_eapc_covid_sensitivity.csv,
#         t1_country_aspr_2021_full.csv
# ============================================================================

import os as _os
# --- repository paths -----------------------------------------------------------
# REPO resolves to the repository root (code/<module>/this_file.py -> ../..).
# Override with the REPRO_ROOT environment variable if you run from elsewhere.
REPO = _os.environ.get("REPRO_ROOT", _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "..")))
# GBD_RAW_DIR: folder with the raw IHME GBD 2021 csv downloads (see data/README.md)
GBD_RAW_DIR = _os.environ.get("GBD_RAW_DIR", _os.path.join(REPO, "data", "raw", "gbd_2021"))
"""T1 master GBD extraction: audit anchors, V8-hypothesis tests, age panel, ASR series,
country ranking, China subnational check. All numbers computed from raw IHME csv files."""
import pandas as pd, numpy as np, glob, json, os

RAW = GBD_RAW_DIR
OUT = REPO
os.makedirs(OUT + "/results", exist_ok=True)

FILES = sorted(glob.glob(RAW + "/IHME-GBD_2021_DATA-*.csv"))
USECOLS = ["measure_name","location_name","sex_name","age_name","cause_name","rei_name",
           "metric_name","year","val","upper","lower"]

frames = []
for f in FILES:
    df = pd.read_csv(f, usecols=USECOLS)
    df = df[df.sex_name == "Both"]
    frames.append(df)
d = pd.concat(frames, ignore_index=True)
print("Both-sex rows:", len(d))

PREV = "Prevalence"; YLD = "YLDs (Years Lived with Disability)"

# ---------- A5a anchors ----------
def get(loc, age, metric, meas, yr):
    m = (d.location_name==loc)&(d.age_name==age)&(d.metric_name==metric)&(d.measure_name==meas)&(d.year==yr)
    r = d[m]
    return (r.val.iloc[0], r.lower.iloc[0], r.upper.iloc[0]) if len(r) else (np.nan,)*3

anchors = {}
for yr in [1990, 2021]:
    anchors[yr] = {
        "prev_number_all_ages": get("Global","All ages","Number",PREV,yr),
        "prev_rate_agestd":     get("Global","Age-standardized","Rate",PREV,yr),
        "yld_number_all_ages":  get("Global","All ages","Number",YLD,yr),
        "yld_rate_agestd":      get("Global","Age-standardized","Rate",YLD,yr),
    }
print(json.dumps(anchors, indent=1))

# ---------- V8 28.3 / 30M hypothesis tests ----------
g21 = d[(d.location_name=="Global")&(d.year==2021)&(d.measure_name==PREV)]
v8 = {}
# H1: sum Number over ALL sexes (we only kept Both; reload needed) -> approximate: Both*3 not valid; instead test on Global row sums
# H2: sum of Number over all age groups incl 'All ages' & 'Age-standardized' rows (double counting)
h2 = g21[(g21.metric_name=="Number")].val.sum()
v8["H2_sum_number_all_age_rows_global_2021"] = h2
# H3: sum over ALL locations (countries+regions+global, multi-level nesting) of All-ages Number
all_locs = d[(d.year==2021)&(d.measure_name==PREV)&(d.age_name=="All ages")&(d.metric_name=="Number")]
v8["H3_sum_number_all_locations_allages_2021"] = all_locs.val.sum()
v8["H3_n_locations"] = all_locs.location_name.nunique()
# H3b: sum over countries only (leaf locations) = should ~= global
v8["n_unique_locations_total"] = d.location_name.nunique()
# H4: mean of national ASPRs 2021 (unweighted) vs global ASPR
nat = d[(d.year==2021)&(d.measure_name==PREV)&(d.age_name=="Age-standardized")&(d.metric_name=="Rate")]
nat_only = nat[~nat.location_name.isin(["Global"]) ]
# identify "aggregate" locations heuristically: keep all, report mean of all non-global and of likely-country rows
v8["H4_mean_agerate_all_nonglobal_2021"] = nat_only.val.mean()
v8["H4_median_agerate_all_nonglobal_2021"] = nat_only.val.median()
# H5: V8 ratio test: 28.31/43.24 ~ 0.655; check YLD rate ratio too
v8["ratio_aspr_v8_over_anchor"] = 28.314137157943303/anchors[2021]["prev_rate_agestd"][0]
# H6: maybe V8 = mean of male+female rate or rate/1000? check Number-based implied pop
v8["implied_pop_v8_millions"] = 30025925.16711431/28.314137157943303*1e5/1e6
v8["implied_pop_anchor_millions"] = anchors[2021]["prev_number_all_ages"][0]/anchors[2021]["prev_rate_agestd"][0]*1e5/1e6
print(json.dumps(v8, indent=1))

# ---------- A4 age panel: Global Both-sex Prevalence, 5-yr age groups ----------
age5 = ['<5 years','5-9 years','10-14 years','15-19 years','20-24 years','25-29 years',
        '30-34 years','35-39 years','40-44 years','45-49 years','50-54 years','55-59 years',
        '60-64 years','65-69 years','70-74 years','75-79 years','80-84 years','85-89 years',
        '90-94 years','95+ years']
ap = d[(d.location_name=="Global")&(d.measure_name==PREV)&(d.age_name.isin(age5))]
num = ap[ap.metric_name=="Number"].pivot_table(index="year", columns="age_name", values="val")
rate = ap[ap.metric_name=="Rate"].pivot_table(index="year", columns="age_name", values="val")
num = num[age5].sort_index(); rate = rate[age5].sort_index()
panel = num.copy()
panel.columns = ["cases_"+c for c in panel.columns]
pop = num/rate*1e5  # derived mid-year population per age-year
pop.columns = ["pop_"+c for c in rate.columns]
panel = pd.concat([panel, pop], axis=1).reset_index()
panel.to_csv(OUT+"/results/t1_bapc_agepanel.csv", index=False)
print("age panel:", panel.shape)
# sanity: sum of 5yr cases vs All ages number
chk = pd.DataFrame({"sum5yr": num.sum(axis=1)})
chk["all_ages"] = [get("Global","All ages","Number",PREV,y)[0] for y in chk.index]
chk["ratio"] = chk.sum5yr/chk.all_ages
print(chk.loc[[1990,2021]])

# ---------- B5 ASR series ----------
asr = d[(d.location_name=="Global")&(d.age_name=="Age-standardized")&(d.metric_name=="Rate")]
series = asr.pivot_table(index="year", columns="measure_name", values="val").reset_index()
series.columns = ["year","ASPR_prev","ASYR_yld"]
series.to_csv(OUT+"/results/t1_asr_series.csv", index=False)

def eapc(y, years):
    x = np.asarray(years, float); ly = np.log(np.asarray(y, float))
    b, a = np.polyfit(x, ly, 1)
    resid = ly - (a + b*x); n = len(x)
    s2 = (resid**2).sum()/(n-2); se = np.sqrt(s2/((x-x.mean())**2).sum())
    from scipy import stats
    tc = stats.t.ppf(0.975, n-2)
    return 100*(np.exp(b)-1), 100*(np.exp(b-tc*se)-1), 100*(np.exp(b+tc*se)-1)

rows = []
for col, lab in [("ASPR_prev","Prevalence ASR"),("ASYR_yld","YLD ASR")]:
    for mask, per in [(series.year<=2021,"1990-2021"),(series.year<=2019,"1990-2019")]:
        sub = series[mask]
        e, lo, hi = eapc(sub[col], sub.year)
        rows.append({"measure":lab,"period":per,"n_years":len(sub),
                     "EAPC_pct":e,"CI_low":lo,"CI_high":hi})
eapc_df = pd.DataFrame(rows)
eapc_df.to_csv(OUT+"/results/t1_eapc_covid_sensitivity.csv", index=False)
print(eapc_df)

# ---------- Table S4: country ranking 2021 ASPR ----------
agg_names = [l for l in nat.location_name.unique() if l in
    ["Global"] ] 
# GBD aggregate locations list (standard)
AGG = set("""Global
High-income Asia Pacific|Central Asia|East Asia|South Asia|Southeast Asia|Australasia|Caribbean|Central Europe|Eastern Europe|Western Europe|Andean Latin America|Central Latin America|Southern Latin America|Tropical Latin America|North Africa and Middle East|Central Sub-Saharan Africa|Eastern Sub-Saharan Africa|Southern Sub-Saharan Africa|Western Sub-Saharan Africa|Oceania|High-income North America
Low SDI|Low-middle SDI|Middle SDI|High-middle SDI|High SDI
World Bank Low Income|World Bank Lower Middle Income|World Bank Upper Middle Income|World Bank High Income
Southeast Asia, East Asia, and Oceania|Central Europe, Eastern Europe, and Central Asia|High-income|Latin America and Caribbean|North Africa and Middle East|South Asia|Sub-Saharan Africa
Commonwealth Low Income|Commonwealth Lower Middle Income|Commonwealth Upper Middle Income|Commonwealth High Income
WHO African Region|WHO Eastern Mediterranean Region|WHO European Region|WHO Region of the Americas|WHO South-East Asia Region|WHO Western Pacific Region
G20|OECD Countries|BRICS|Asian Development Bank|Asia|Europe|Africa|America
Least Developed Countries|Small Island Developing States""".replace("|","\n").split("\n"))
cand = nat[~nat.location_name.isin(AGG)].copy()
# keep only locations present with 2021 both-sex age-standardized rate and that are countries:
# heuristic: country leaf = location whose Number all-ages row exists and not in AGG
countries = cand[["location_name","val"]].rename(columns={"val":"ASPR_2021"})
# attach cases
cn_num = d[(d.year==2021)&(d.measure_name==PREV)&(d.age_name=="All ages")&(d.metric_name=="Number")][["location_name","val"]].rename(columns={"val":"cases_2021"})
countries = countries.merge(cn_num, on="location_name")
countries = countries.sort_values("ASPR_2021", ascending=False).reset_index(drop=True)
countries["rank"] = countries.index+1
countries.to_csv(OUT+"/results/t1_country_aspr_2021_full.csv", index=False)
top10 = countries.head(10); bot10 = countries.tail(10)
print("n locations kept as countries:", len(countries))
print(top10[["rank","location_name","ASPR_2021"]]); print(bot10[["rank","location_name","ASPR_2021"]])

# ---------- C2 China subnational check ----------
all_locs_all = pd.concat([pd.read_csv(f, usecols=["location_name"]) for f in FILES]).location_name.unique()
prov = [l for l in all_locs_all if any(k in l for k in
    ["Beijing","Sichuan","Guangdong","Shanghai","Yunnan","Shandong","Henan","Hubei","Hunan",
     "Zhejiang","Jiangsu","Anhui","Fujian","Gansu","Guizhou","Heilongjiang","Jilin","Liaoning",
     "Shaanxi","Shanxi","Hebei","Jiangxi","Guangxi","Hainan","Ningxia","Qinghai","Tibet",
     "Xinjiang","Inner Mongolia","Chongqing","Tianjin","Hong Kong","Macau","Macao","Taiwan"])]
china = [l for l in all_locs_all if "China" in l]
print("China-related locations:", china)
print("Province-like locations:", prov)

json.dump({"anchors":anchors,"v8_hypotheses":v8,
           "china_locations":list(china),"province_locations":list(prov)},
          open(OUT+"/results/t1_extract_summary.json","w"), indent=1, default=float)
print("DONE")
