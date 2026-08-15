# PIPELINE — module map, inputs, outputs, manuscript links

Run order: modules are largely independent once `data_downloads/` is
populated (`bash code/00_prepare_data_downloads.sh`). Scripts are run from
the repository root (or with `REPRO_ROOT` set). R scripts may require
`OPENGWAS_JWT` (see README § Authentication). `code/99_validation/run_validation.sh`
re-runs all numeric checks from VALIDATION.md in one command.

Legend: **[shipped]** = input included in `data/derived/` or `results/legacy/`;
**[download]** = large external input, see `data/README.md`;
**[API]** = needs an OpenGWAS token; **[cached]** = the API product is
shipped, so the step can be skipped.

## Module 01 — MR discovery (European discovery cohort)

| script | inputs | outputs | manuscript link |
|---|---|---|---|
| `code/01_mr_discovery/t1_export_base_harmonised.R` | `data/derived/mr_discovery/discovery_harmonised_10snp_dat.xlsx` **[shipped]** | `results/t1_base_harmonised_10snp.csv` | 10-SNP instrument set; Table S1 |
| `code/01_mr_discovery/v8_07_european_mr_robustness.R` *(legacy)* | full COPD + HF sumstats **[download]** | discovery IVW OR 1.15 (1.08–1.22), sensitivity battery | Results ¶1; Table S2; Fig. 1b |
| `code/01_mr_discovery/v8_01_cause_analysis.R` *(legacy)* | discovery xlsx **[shipped]** | CAUSE gamma (sensitivity) | Table S2 |
| `code/01_mr_discovery/v8_02_mrclust_analysis.R` *(legacy)* | discovery xlsx **[shipped]** | MR-Clust clusters | Table S2 |
| `code/01_mr_discovery/t1_governance_tables.py` | `results/legacy/real_data_results.json`, `data/derived/bulk/*` **[shipped]** | `tables/rebuilt/Table_S2, S3, S5, S6, S7` (rebuilt copies; curated finals ship in `tables/`) | Suppl. tables |

## Module 02 — MVMR

| step | script | inputs | outputs | manuscript link |
|---|---|---|---|---|
| 1 | `code/02_mvmr/extract_sumstats.sh` | `data_downloads/GCST*.h.tsv.gz` **[download]** | `iv_raw_*.tsv`, `assoc_A1A_*.tsv`, `assoc_B2_*.tsv` **[shipped in `data/derived/mvmr/`]** | — |
| 2 | `code/02_mvmr/t1_clump_all.R` | `iv_raw_*.tsv` + LD API **[API]** | `iv_final_*.tsv` **[shipped]** | Methods (LD r²<0.001, 10 Mb) |
| 3 | `code/02_mvmr/t1_make_unions2.R` | `iv_final_*.tsv`, `snps10.rds` **[shipped]** | `t1_iv_sets.rds` **[shipped]** | — |
| 4a | `code/02_mvmr/t1_a1a_mvmr_final.R` | `assoc_A1A_*.tsv` **[shipped]** | `results/t1_mvmr_results.csv` | Model 1: COPD+SMOK+BMI+SBP, 709 SNP, OR 1.075; Fig. 1d |
| 4b | `code/02_mvmr/t1_b2_mvmr_final.R` | `assoc_B2_*.tsv` **[shipped]** | `results/t1_mvmr_lungfunction.csv` | Model 2: COPD+FEV1+FEV1/FVC+CPD, 405 SNP, OR 1.165; Fig. 1d |

## Module 03 — Mediation MR

| step | script | inputs | outputs | manuscript link |
|---|---|---|---|---|
| 1 | `code/03_mediation/t1_med_stepA.R` | OpenGWAS **[API; cached]** | `serpine1_eqtla_region.rds`, `iv_med_TGFB1.rds` **[shipped in `data/derived/mediation/`]** | — |
| 2 | `code/03_mediation/t1_med_stepB.R` | `t1_base_harmonised_10snp.csv`, mediation caches **[shipped]**, IL6/CRP/HF sumstats **[download]** | `results/t1_mediation_results.csv` | CRP 4.4% mediated; IL-6/SERPINE1/TGFB1 null; Fig. 1e; Table S3 |

## Module 04 — Replication + meta-analysis + power

| script | inputs | outputs | manuscript link |
|---|---|---|---|
| `code/04_replication_meta/t1_b1_finngen_hf_stratified_mr.R` | `t1_base_harmonised_10snp.csv`, `data/derived/finngen/fg_*` **[shipped]** | `results/t1_hfpef_hfref_*.csv` | FinnGen R12 strict HF, OR 1.037; Fig. 1c |
| `code/04_replication_meta/v8_06_asian_mr_update.R` *(legacy)* | OpenGWAS BBJ **[API]** | `results/legacy/06_asian_*.csv` | BBJ OR 0.933 |
| `code/04_replication_meta/t1_b3_asian_power_analysis.R` | `results/legacy/06_asian_harmonised_data.csv` **[shipped]** | `results/t1_asian_power_*.csv/txt` | BBJ power 98.2% |
| `code/04_replication_meta/three_cohort_meta.R` | 3 cohort log-OR/SE (in script) | `results/t1_three_cohort_meta.{json,csv}` | FE OR 1.059 (1.022–1.098); Q=14.5, I²=86%; RE OR 1.042; Fig. 1c |

## Module 05 — GBD 2021 burden

| script | inputs | outputs | manuscript link |
|---|---|---|---|
| `code/05_gbd/t1_extract_gbd.py` | IHME raw csv **[download]** | `t1_extract_summary.json`, `t1_bapc_agepanel.csv`, `t1_asr_series.csv`, `t1_eapc_covid_sensitivity.csv`, `t1_country_aspr_2021_full.csv` | Table 3 anchors; EAPC −0.133; Fig. 2 |
| `code/05_gbd/t1_bapc_run.R` | `t1_bapc_agepanel.csv` | `t1_bapc_projections.csv`, `t1_bapc_total_cases{,_frozenpop}.csv` | 6.13 M by 2035; 10.93 M by 2050; Fig. 5a |
| `code/05_gbd/t1_bapc_figure_summary.py` | BAPC outputs | `t1_bapc_summary.json`, `tables/Table_S4`, `figures/t1_bapc_projection.png` | Table 3; Fig. 5a |
| `code/05_gbd/t1_regen_s4_with_ui.py` | IHME raw **[download]** | `t1_country_aspr_2021_strict.csv` | 26-fold gradient; Fig. 3 |
| `code/05_gbd/gbd_fig_data.py` *(new)* | IHME raw **[download]** | `fig_sdi_series_{prev,yld}.csv`, `fig_country_1990_2021.csv`, `fig_eapc_prev_regions.csv` | Figs. 2, 3; exact reproduction verified |
| `code/05_gbd/gbd_decomposition.py` *(new)* | IHME raw **[download]** | `gbd_decomposition.csv` | ageing +53.3% / growth +44.7% / epi +1.9%; Fig. 4a,b |
| `code/05_gbd/gbd_inequality.py` *(new)* | IHME raw + `gbd2021_sdi.csv` **[not shipped]** | `gbd_inequality.csv` | SII −40.05→−23.10; CI −0.25→−0.15; Fig. 4c–f |
| `code/05_gbd/gbd_frontier.py` *(new)* | IHME raw + `gbd2021_sdi.csv` **[not shipped]** | `gbd_frontier.csv` | 15 frontier-lagging countries; Fig. 4g,h |
| `code/05_gbd/gbd_intervention_scenarios.py` *(new)* | `t1_bapc_total_cases.csv` | `gbd_intervention_scenarios.csv` | ~2.5 M / ~1.8 M averted by 2050; Fig. 5c |
| `code/05_gbd/t1_bridge.py` | constants (verified anchors) | `t1_bridge_table.csv` | Suppl. Note 2; Fig. 5d |
| `code/05_gbd/v8_08_bapc_projection.R` *(legacy)* | V8 data tree **[not shipped]** | frequentist APC comparator (+73% → 6.25 M) | Results ¶ projection |

## Module 06 — Bulk transcriptomics

| script | inputs | outputs | manuscript link |
|---|---|---|---|
| `code/06_bulk_transcriptomics/bulk_limma_reanalysis.R` *(new)* | GEO GSE57148/GSE57338 **[download]** | limma DE tables; 41 shared genes | Fig. 6a–c; Table S5 |
| `code/06_bulk_transcriptomics/bulk_keygene_welch.py` *(new)* | `data/derived/bulk/*.xlsx` **[shipped]** | `bulk_keygene_welch.csv` | key-gene P values; Fig. 6f,g; Table S7 |
| `code/06_bulk_transcriptomics/t1_a3_crosscohort.py` | bulk xlsx **[shipped]** + h5ad **[download]** | `t1_crosscohort_validation.csv` | 41-gene 4-cohort consistency; Fig. 6h |
| — | STRING v12 + cytoHubba (web tools, not scripts) | PPI network; 9 hub genes: SERPINE1, THBS1, SPP1, IL1R2, EGR1, CDKN1A, CD163, AREG, SOCS3 | Fig. 6e; parameters: combined score > 0.4, cytoHubba MCC top-9 |

## Module 07 — Single-cell

| script | inputs | outputs | manuscript link |
|---|---|---|---|
| `code/07_singlecell/t1_a2_copd_analysis.py` | `data/t1_copd_lung.h5ad` **[download]** | `t1_disease_sc_de.csv`, `t1_disease_serpine1_lr*.csv`, `t1_pseudotime_correlation.csv` | COPD lung disease-state DE; SERPINE1 L–R; pseudotime; Fig. 7 |
| `code/07_singlecell/t1_a2_heart_analysis.py` | 3 heart h5ad **[download]** | `t1_heart_sc_de.csv` | THBS1/SPP1 fibroblasts, SERPINE1 endothelium; Fig. 7 |
| `code/07_singlecell/t1_donor_pseudobulk.py` | h5ad **[download]** | `t1_donor_pseudobulk_validation.csv` | donor-level replication = Table S8 |

## Module 08 — Drug target (SERPINE1)

| step | script | inputs | outputs | manuscript link |
|---|---|---|---|---|
| 1 | `code/08_drugtarget/t1_minitabix.py` | eQTL Catalogue FTP (remote tabix) | `eqtl_<tissue>_serpine1.tsv` **[not shipped]** | — |
| 2 | `code/08_drugtarget/t1_a1b_prep_candidates.R` | region tsv | `serpine1_gene/*.csv` **[shipped]** | — |
| 3 | `code/08_drugtarget/t1_a1b_fetch_opengwas.R` | OpenGWAS **[API; cached]** | `t1_a1b_hf_{lead,region}.rds` **[shipped]** | — |
| 4 | `code/08_drugtarget/t1_a1b_analysis.R` | shipped caches | `t1_drugtarget_mr_serpine1_7tissues.csv`, `t1_coloc_serpine1.csv` | 8-tissue cis-eQTL MR + coloc; Fig. 7; Table S6 |
| 5 | `code/08_drugtarget/t1_a1b_finalize.R` | shipped caches | updated MR table | eQTLGen whole blood |
| 6 | `code/08_drugtarget/t1_a1b_eqtlgen_addendum.R` | shipped caches + LD API **[API]** | `t1_a1b_harmonised_blood_eQTLGen.csv`, coloc update | — |
| 7 | `code/08_drugtarget/t1_c1_finngen_phewas_rs7860931.R` | `fg_r12_phewas_rs7860931.json` **[shipped]** | `t1_c1_finngen_r12_phewas_rs7860931*.csv` | PheWAS 2,470 endpoints; Fig. 7 |
| 8 | `code/08_drugtarget/v8_04_serpine1_pqtl_mr.R` *(legacy)* | OpenGWAS **[API]** | `results/legacy/serpine1_pqtl_mr_*` | deCODE cis-pQTL MR: OR 1.01, P=0.69 |
| 9 | `code/08_drugtarget/v8_05_phewas_serpine1.R` *(legacy)* | PhenoScanner | `results/phewas_serpine1_ieu_*` | IEU PheWAS |

## Module 09 — Figures

`code/09_figures/make_fig1.py` … `make_fig7.py` (+ shared `figstyle.py`)
rebuild `figures/Fig1–7 *.png/pdf` from `results/`. `audit_figures.py`
re-checks every plotted number against `results/` and prints PASS/FAIL
(latest run: 131 checks, 130 PASS, 1 display-rounding FAIL; see
VALIDATION.md). Fig. 4 panels A–C and Fig. 5 panel B contain hard-coded
values that are computed by the module-05 scripts marked *(new)* above; the
numerical verification (and one documented sign discrepancy) is in
`VALIDATION.md`.
