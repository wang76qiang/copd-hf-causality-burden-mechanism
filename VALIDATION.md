# VALIDATION

Date of validation run: 2026-07-27. Every check below was executed from this
repository against the raw external data or the shipped intermediates.
Comparison targets are the values printed in the manuscript.

## 1. Verified (re-computed, matches manuscript)

### 1.1 Three-cohort meta-analysis (`code/04_replication_meta/three_cohort_meta.R`)
| Quantity | Manuscript | Recomputed | Match |
|---|---|---|---|
| Fixed-effect OR (95% CI) | 1.059 (1.022-1.098) | 1.0593 (1.0220-1.0980) | YES |
| FE p | 0.0016 | 0.00164 | YES |
| Heterogeneity Q (df=2) | 14.5 | 14.503 | YES |
| I2 | 86% | 86.2% | YES |
| Random-effects OR | "random-effects null" | 1.042 (0.940-1.155), p=0.433 | YES (null) |
Output: `validation/t1_three_cohort_meta.{csv,json}`.

### 1.2 GBD 2021 master extraction (`code/05_gbd/t1_extract_gbd.py`)
| Quantity | Manuscript | Recomputed | Match |
|---|---|---|---|
| 2021 global cases | 3,613,136 (2.67-4.82 M UI) | 3,613,136.36 (2,666,978-4,801,973) | YES |
| 2021 ASPR /100k | 43.24 | 43.2373 | YES |
| 2021 ASYR /100k | 3.85 | 3.8464 | YES |
| 1990 global cases / ASPR / ASYR | 1,504,707 / 41.93 / 3.72 | identical | YES |
| EAPC ASPR 1990-2021 | -0.133 | -0.1327 (-0.215, -0.050) | YES |
| EAPC COVID sensitivity (1990-2019) | "unchanged when excluding COVID-era years" | -0.1454 (-0.236, -0.054) | YES |
| Country ranking (227 rows) | PNG top (83.09), Uzbekistan bottom | identical, exact | YES |
Outputs: `validation/extract/`. The V8 aggregation error (30.0 M / 28.31)
was also reproduced exactly (H3 hypothesis tests), confirming the audit.

### 1.3 GBD figure-data regeneration (`code/05_gbd/gbd_fig_data.py`, CHECK=1)
| File | Rows | Max abs diff vs shipped results/ |
|---|---|---|
| fig_sdi_series_prev.csv | 32 | 0 |
| fig_sdi_series_yld.csv | 32 | 0 |
| fig_country_1990_2021.csv | 210 | 0 |
| fig_eapc_prev_regions.csv | 29 | 0.00049 (EAPC rounding) |
The high-SDI EAPC of +1.33 quoted in the manuscript is reproduced inside
`fig_eapc_prev_regions.csv`.

### 1.4 Das Gupta decomposition (`code/05_gbd/gbd_decomposition.py`)
| Quantity | Manuscript | Recomputed | Match |
|---|---|---|---|
| Global ageing share | +53.3% | +53.32% | YES |
| Global population-growth share | +44.7% | +44.74% | YES |
| Global epidemiological share | +1.9% | +1.94% | YES |
| High-SDI epidemiological share | "+34.5%, only positive driver" | +34.46% | YES |
Output: `validation/gbd_decomposition.csv`.

### 1.5 Bulk key-gene Welch tests (`code/06_bulk_transcriptomics/bulk_keygene_welch.py`)
| Gene (cohort) | Manuscript | Recomputed | Match |
|---|---|---|---|
| S100A8 (GSE57148) | P = 8e-6 | P = 7.8e-6 | YES |
| S100A12 (GSE57148) | P = 0.02 | P = 0.0248 | YES |
| IL1RL1 (GSE57148) | P = 0.03 | P = 0.0328 | YES |
| THBS1 (GSE57148) | P = 0.001 | P = 9.8e-4 | YES |
| LUM (GSE57338) | P <= ~2e-32 | P = 2.0e-41 | YES |
| ASPN (GSE57338) | P <= ~2e-32 | P = 8.6e-42 | YES |
| SMOC2 (GSE57338) | P <= ~2e-32 | P = 2.5e-57 | YES |
| ACE (GSE57338) | P <= ~2e-32 | P = 2.0e-32 | YES |
| IL1RL1 (GSE57338) | strongly reduced in failing hearts | P = 1.1e-35 (down) | YES |
Shapiro-Wilk rejects normality for most groups, so the manuscript reports
Welch (unequal-variance) P values, which this script reproduces.
Output: `validation/bulk_keygene_welch.csv`.

### 1.6 End-to-end re-runs of the R pipeline from shipped caches
All of the following were re-executed from this repository after
`bash code/00_prepare_data_downloads.sh` (R 4.6.1); every output is
byte-identical to the shipped `results/` file:

| script | key output | result |
|---|---|---|
| `code/02_mvmr/t1_a1a_mvmr_final.R` | `t1_mvmr_results.csv` | COPD OR 1.075 (1.024-1.128), 709 SNP, condF 1.76 - identical |
| `code/02_mvmr/t1_b2_mvmr_final.R` | `t1_mvmr_lungfunction.csv` | COPD OR 1.165 (1.102-1.232), 405 SNP, condF 2.19 - identical |
| `code/04_replication_meta/t1_b1_finngen_hf_stratified_mr.R` | `t1_hfpef_hfref_mr.csv` | FinnGen R12 strict HF IVW OR 1.037 (0.986-1.090) - identical |
| `code/04_replication_meta/t1_b3_asian_power_analysis.R` | `t1_asian_power_analysis.csv` | BBJ power 98.2% for OR 1.15 - identical |
| `code/08_drugtarget/t1_a1b_analysis.R` | `t1_drugtarget_mr_serpine1_7tissues.csv`, `t1_coloc_serpine1.csv` | 8-tissue Wald ratios + coloc - identical |
| `code/08_drugtarget/t1_c1_finngen_phewas_rs7860931.R` | `t1_c1_finngen_r12_phewas_rs7860931.csv` | 2,470-endpoint PheWAS table - identical |
| `code/01_mr_discovery/t1_governance_tables.py` | `tables/rebuilt/Table_S2,S5,S7` | identical to curated `tables/` versions |
| `code/05_gbd/t1_bridge.py` | `t1_bridge_table.csv` | implied PAF 6.02%, GBD/MR ratio 8.0x |

### 1.7 Figure regeneration and audit
All seven `code/09_figures/make_fig*.py` scripts were executed from this
repository and regenerated `figures/Fig1-7`. `code/09_figures/audit_figures.py`
(131 automated checks of plotted numbers vs `results/`): **130 PASS, 1 FAIL**
- the single FAIL is `F5B combined: fig=6.6 vs ref=6.63`, a display-rounding
artifact in the figure annotation (10.93 - 2.5 - 1.8 = 6.63), not a data
error.

## 2. Discrepancies found during packaging (ALL RESOLVED 2026-07-27)

> **Resolution note (v1.0.1, 2026-07-27).** Both discrepancies below were
> corrected in the manuscript, tables and figures: Fig. 4b now displays the
> validated High-SDI decomposition (+44.8/+20.7/+34.5); Fig. 5b now uses the
> validated scenario values (10.93 − 2.93 − 1.83 = 6.17 M), and the text and
> Table 3 were updated accordingly (smoking cessation 2.9 M, optimized
> therapy 1.8 M; PAF labels 0.45/0.74/1.48%). The corrected figures are the
> ones shipped in this repository (`figures/Fig4_*`, `figures/Fig5_*`).
>
> **Resolution note (v1.0.2, 2026-07-27).** A third data issue was found and
> fixed: the strict country table (`results/t1_country_aspr_2021_strict.csv`)
> contained six UN sub-regional aggregate rows (Southern/Eastern/Western/
> Central/Northern Africa, North America) that had leaked through the
> aggregate filter, inflating the country count to 210. The filter in
> `code/05_gbd/t1_regen_s4_with_ui.py` and `gbd_fig_data.py` now excludes
> them explicitly; the corrected files contain exactly the **204** GBD 2021
> countries and territories (top: Papua New Guinea 83.09; bottom: Uzbekistan
> 3.20; 26-fold gradient and all rankings unchanged). The manuscript, Table 1,
> Table S4 and Fig. 3 were updated to the 204-country framing.

### 2.1 Intervention scenarios (`code/05_gbd/gbd_intervention_scenarios.py`) — RESOLVED
Manuscript/ Fig. 5b previously printed "~2.5 M (smoking cessation) and ~1.8 M
(optimized COPD pharmacotherapy) averted cases by 2050".
Recomputed under the documented convention (RR applied to the projected
growth above the 2021 caseload; baseline 2050 = 10.93 M):
- smoking cessation (RR 0.60): **2.93 M averted in 2050** (was ~2.5 M)
- pharmacotherapy (RR 0.75): **1.83 M averted in 2050** (matches ~1.8 M)
The ~2.5 M printed for smoking matched neither the applied-to-growth
convention (2.93 M) nor the naive whole-caseload convention (4.37 M). The
manuscript text, Table 3 and Fig. 5b were updated to the validated values on
2026-07-27 (2.9 M / 1.8 M; combined 6.2 M cases in 2050).
Output: `validation/gbd_intervention_scenarios.csv`.

### 2.2 Fig. 4b high-SDI decomposition sign labelling (`code/09_figures/make_fig4.py`) — RESOLVED
`gbd_decomposition.py` computes for High SDI: growth +20.7%, ageing +44.8%,
epidemiological +34.5%. Fig. 4b previously displayed growth = -20.8 and
ageing = -44.7 (flipped signs; the three shares summed to -32% instead of
+100%). The figure script was corrected on 2026-07-27 (prevalence-based
decomposition for all three SDI groups) and Fig. 4b was regenerated with the
validated values.

## 3. Not verified (inputs unavailable)

### 3.1 SII / concentration index (`code/05_gbd/gbd_inequality.py`)
Manuscript: SII (ASPR) -40.05 (1990) -> -23.10 (2021); CI -0.25 -> -0.15.
Requires country-year GBD 2021 SDI values (`data/derived/gbd2021_sdi.csv`),
which are **not shipped**: the SDI series is distributed only via IHME GHDx
behind a (free) login and no local copy existed at packaging time. The
script is complete and fail-fast; run it after downloading the SDI file
(see `data/README.md`, source 8). NOT VERIFIED.

### 3.2 Frontier analysis (`code/05_gbd/gbd_frontier.py`)
Manuscript: 15 frontier-lagging countries (incl. China, India, the
Netherlands, Canada, Australia). Same missing SDI input as 3.1.
NOT VERIFIED.

### 3.3 Components intentionally not re-executed
- `bulk_limma_reanalysis.R` (GEOquery + limma full re-run): runnable but not
  executed (GEO network + ~100 MB SOFT matrices); the same group statistics
  are independently verified in 1.5 from the shipped workbooks.
- BAPC/INLA projection (`t1_bapc_run.R`): requires R-INLA; the projection
  outputs are shipped and used by validated downstream steps (1.x, 2.1).
- Single-cell modules (07): h5ad files (~3.7 GB) not shipped; outputs are
  shipped and cross-checked by module 06/09.
- OpenGWAS-API-dependent fetch steps: require a personal JWT and are cached
  in `data/derived/`; the analysis steps reading those caches are runnable.


## MedComm v1.1.0 update (2026-08-15)

See `VALIDATION_ADDENDUM_MedComm.md`. The MedComm update removes the previously listed fourth author from the author list, updates the target journal after the earlier EClinicalMedicine submission formally ended, uses the corrected 204-country framing throughout, rebuilds Tables S1-S10 in English, and replaces Figures 5-7 with the corrected versions.
