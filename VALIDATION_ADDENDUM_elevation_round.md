# Elevation-round addendum (JTM submission) — 2026-09-04

## Scope

Four pre-specified elevation analyses executed; two further analyses pre-registered but not executable with public data.

## Executed analyses and key outputs

1. **Winner's curse correction** (`code/10_elevation_round/run_analysis12.py`, pure-Python engine `mr_engine.py`)
   - Zhong–Prentice conditional-likelihood correction of the 10 discovery exposure betas (selection threshold P < 5×10⁻⁸).
   - Result: IVW OR 1.151 → **1.156** (95% CI 1.078–1.240); per-SNP shrinkage < 2% for 7/10 instruments.
   - Orthogonal check: exposure weights re-estimated in FinnGen R12 J10_COPD / COPD_LATER / COPD_EARLY → MR vs discovery HF: **OR 1.19 / 1.19 / 1.17** (all P ≤ 3.1×10⁻⁶).
   - Outputs: `results_elevation/winners_curse_persnp.csv`, `finngen_copd_alt_exposure_mr.csv`.

2. **Extended replication** (`run_analysis3_finngen.py`, `extract_gwascat.py`)
   - FinnGen R12: I9_HEARTFAIL_AND_CHD OR 1.048 (0.973–1.128); I9_HEARTFAIL_AND_OVERWEIGHT OR 1.034 (0.968–1.105); FG_OTHHEART OR 0.993; I9_OTHILLHEART OR 0.879.
   - BBJ 2025 suite (Enzan et al., Nat Commun 2025; GWAS Catalog GCST90668009/10/11, harmonised files fetched via HTTP-range tabix, `extract_gwascat.py`): all-cause HF OR 0.925 (0.859–0.996); **HFpEF OR 0.865 (0.776–0.964)**; HFrEF OR 1.002 (0.862–1.165).
   - 9/10 and 8/10 instruments available (rs114904431 and rs4488938 absent in some EAS files; MAF/imputation filters) — documented in harmonisation logs.
   - Outputs: `results_elevation/finngen_extended_mr.csv`, `bbj2025_stratified_mr.csv`, harmonisation logs.

3. **Heterogeneity forensics** (`run_analysis4_heterogeneity.py`)
   - k = 10 independent cohort–endpoint IVW estimates (winner's-curse-corrected discovery excluded from pooling, shown as sensitivity).
   - Overall FE OR 1.013 (0.991–1.035), I² = 77%; European k = 6: OR 1.036 (P = 0.004); East Asian k = 4: **OR 0.922 (P = 0.001), I² = 0%**.
   - Illustrative weighted meta-regression (ancestry, strictness, HFpEF): QM = 19.3, P = 2.4×10⁻⁴, R² analogue 50%.
   - Cohort funnel Egger P = 0.23 (illustrative).
   - Outputs: `results_elevation/heterogeneity_package.json`, `heterogeneity_cohorts.csv`.

4. **Interorgan ligand–receptor map** (`run_b4_crosstalk.py`, memory-safe h5py streaming)
   - 31 curated pairs × lung senders × heart receivers = 1,428 directed edges.
   - Attenuated in disease: capillary-endothelium SERPINE1 → PLAUR (−2.26 log₂) and SERPINE1 → VTN (−1.71); increased: fibroblast S100A9 → TLR4 (+1.03).
   - Donor-level robust (FDR < 0.05): cardiac fibroblast/endothelial/macrophage **CD36 ↑** (P = 1.6×10⁻⁴ / 2.4×10⁻³ / 2.4×10⁻³); macrophage **IL1RAP ↑** (P = 1.1×10⁻³); endothelial **ITGB3 ↑** (P = 3.7×10⁻³); macrophage **CD44 ↓** (P = 9.4×10⁻⁵).
   - Outputs: `results_elevation/b4_interorgan_edges_lung_to_heart.csv`, `b4_donor_level_tests.csv`.

5. **Extended SERPINE1 cis screen** (`run_a2_sqtl.py`)
   - eQTL Catalogue GTEx v8 exon/tx/txrev/leafcutter permuted files × 5 cardiorespiratory tissues (20 datasets).
   - No SERPINE1 cis-QTL in any molecular layer (all gene-level adjusted P > 0.06).
   - Output: `results_elevation/a2_sqtl_screen_summary_v2.csv`.

## Engine validation

- Pure-Python IVW/heterogeneity implementation reproduces the published discovery result exactly (OR 1.151, 95% CI 1.083–1.223; Q = 18.89; I² = 52.3%).
- FinnGen extraction validated against shipped `derived/finngen/fg_10snps_raw.tsv` column semantics.

## Pre-registered but not executable with public data

- European HFpEF/HFrEF MR: MVP GWAS restricted (dbGaP phs001672).
- Individual-level COPD-PRS → incident HF (UK Biobank access required).
- PM2.5 MR layer (OpenGWAS ukb-b-10817 requires a personal JWT; not available in this environment).

## New figure

- `figures/Figure_8_replication_heterogeneity_crosstalk.png/.pdf` (`code/10_elevation_round/make_fig8.py`).

## New data downloads (not shipped; re-fetchable)

- FinnGen R12 endpoint files (streamed extraction script `fetch_finngen.sh`; only 10-SNP rows kept).
- BBJ 2025 HF suite harmonised files (HTTP-range extraction, `extract_gwascat.py`).
- CELLxGENE COPD lung (`8fbed309...`, 849 MB) and Reichart heart fibroblast/endothelial/macrophage h5ad (`fetch_scdata.sh`).

---

## Round 2 (editorial-review-driven revisions, 2026-09-04)

6. **Instrument portability (Table S16)** — per-SNP exposure-effect-allele frequencies across EUR discovery, FinnGen R12, and BBJ 2025 files; availability flags per endpoint. Large EUR→EAS frequency shifts documented (e.g., rs11525583 0.175→0.628; rs34944514 0.467→0.858; rs114904431 absent in EAS files).
7. **Ancestry-stratified PAF sensitivity (Table S17; `run_m5_paf.py`)** — Levin's formula under discovery OR 1.15 (PAF 0.45–1.48%), European subgroup OR 1.036 (0.11–0.36%), East Asian subgroup OR 0.922 (negative PAF), overall OR 1.013 (0.04–0.13%).
8. **Interorgan edge permutation tests (Table S18; `run_b4_permutation.py`)** — 1,000 donor-label permutations per organ; pair-level empirical P. SERPINE1→VTN P = 0.012 (FDR 0.096); SERPINE1→PLAUR P = 0.097; others NS. Bug-fix note: an intermediate cached extraction had a gene-column mapping error (searchsorted on unsorted indices); caches were rebuilt and the fix validated against the original B4 donor tests.
9. **Formal SAP for the three pre-registered analyses** — `SAP_preregistered_analyses.md` (R1 MVP HFpEF/HFrEF MR; R2 UKB PRS→incident HF; R3 PM2.5 MR), with hypotheses, estimators, and success/failure criteria.

Manuscript text revisions in round 2: abstract ancestry qualification (M2) and lung-function softening (M3); instrument-availability and portability reporting (M1); ten-estimate provenance; ancestry-mechanism decomposition in Discussion; ancestry-stratified PAF paragraph in Results; donor-frame disclosure in Methods (M7); GBD-2021 version justification in Limitations (M6); Table 1 note disambiguation; Statistics software cross-reference; edge-permutation framework in Methods/Results/Fig. 8 legend (M4).

---

## v1.2.1-jtm update (2026-09-05)

Figure 8 in this archive was regenerated to the final manuscript revision: panel A now shows the European replication-only subgroup diamond (five FinnGen endpoints, discovery excluded; OR 1.016, 95% CI 0.990-1.043, P = 0.23, I2 = 20%) that the manuscript text and legend describe. Table S10 case-count uncertainty intervals were corrected to the audited extraction bounds (1990: 1.50 M, 95% UI 1.14-1.94 M; 2021: 3.61 M, 95% UI 2.67-4.80 M).
