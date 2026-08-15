# COPD and heart failure — reproducibility package

Code, derived data, results, tables and figures for the manuscript:

> **Chronic obstructive pulmonary disease and heart failure: genetic evidence
> of causality, global burden, and a tissue-confined inflammo-fibrotic
> mechanism** — Hongjuan Fang, Lijuan Bu, Qiang Wang.
> Prepared for submission to *MedComm* (2026).

The paper combines three evidence layers:

1. **Genetic causality** — two-sample Mendelian randomization (MR) of COPD on
   heart failure in a European discovery cohort (IVW OR 1.15, 95% CI
   1.08–1.22), replication in FinnGen R12 and BioBank Japan, three-cohort
   meta-analysis, multivariable MR (MVMR) and two-step mediation MR.
2. **Global burden** — GBD 2021 COPD-attributable HF in 204 countries and territories,
   1990–2021: EAPC trends, Das Gupta decomposition, SII/CI inequality,
   frontier analysis, BAPC projection to 2050, Levin PAF, GBD–MR bridge and
   intervention scenarios.
3. **Mechanism** — multi-scale transcriptomics: bulk limma re-analysis
   (GSE57148 lung, GSE57338 heart) yielding 41 shared inflammo-fibrotic
   genes; single-cell validation in COPD lung (CELLxGENE) and failing heart
   (Reichart 2022) with donor-level pseudobulk; SERPINE1 drug-target
   assessment (cis-eQTL MR, colocalization, PheWAS, cis-pQTL MR).

## Repository layout

```
├── README.md  LICENSE  CITATION.cff  requirements.txt  r_packages.txt
├── PIPELINE.md        # module map: script -> inputs -> outputs -> manuscript
├── VALIDATION.md      # executed numerical verifications vs the manuscript
├── code/
│   ├── 00_prepare_data_downloads.sh   # stage shipped intermediates -> data_downloads/
│   ├── 01_mr_discovery/      # discovery MR export + legacy discovery scripts
│   ├── 02_mvmr/              # sumstat extraction, LD clump, 2 MVMR models
│   ├── 03_mediation/         # two-step mediation MR (CRP/IL-6/SERPINE1/TGFB1)
│   ├── 04_replication_meta/  # FinnGen + BBJ replication, power, 3-cohort meta
│   ├── 05_gbd/               # GBD extraction, BAPC, decomposition, inequality,
│   │                         #   frontier, scenarios, bridge, fig data
│   ├── 06_bulk_transcriptomics/  # limma re-analysis, Welch key-gene tests,
│   │                           #   41-gene cross-cohort validation
│   ├── 07_singlecell/        # COPD lung + heart sc analyses, donor pseudobulk
│   ├── 08_drugtarget/        # SERPINE1 cis-eQTL MR/coloc, PheWAS, pQTL MR
│   ├── 09_figures/           # make_fig1–7 + audit
│   └── 99_validation/        # one-command re-run of the VALIDATION.md checks
├── data/
│   ├── README.md      # where to get every large dataset (required reading)
│   ├── derived/       # small shipped intermediates (see PIPELINE.md)
│   └── raw/           # (empty) place IHME GBD csv here or set GBD_RAW_DIR
├── results/           # all analysis outputs (csv/json/md), incl. legacy/
├── tables/            # Table_S1–S8 (curated; tables/rebuilt/ = regenerated)
├── figures/           # Fig1–7 (png + pdf) plus analysis-level figures
└── validation/        # outputs of the latest validation re-run
```

## Quick start

```bash
# 1) install dependencies (Python 3.10+, R 4.4+)
pip install -r requirements.txt
# R: install the packages listed in r_packages.txt (CRAN + Bioconductor +
# the INLA repo + jean997/cause + chr1swallace/mrclust)

# 2) stage the shipped intermediates into the pipeline working directory
bash code/00_prepare_data_downloads.sh

# 3) run modules (see PIPELINE.md for order and inputs), e.g.
Rscript code/04_replication_meta/three_cohort_meta.R
python code/06_bulk_transcriptomics/bulk_keygene_welch.py
GBD_RAW_DIR=/path/to/ihme_csv python code/05_gbd/gbd_decomposition.py

# 4) or re-run every numeric check in one go:
GBD_RAW_DIR=/path/to/ihme_csv bash code/99_validation/run_validation.sh
```

Scripts resolve the repository root automatically (or honour `REPRO_ROOT`).
Outputs are written to `results/` (and `logs/`, `data_downloads/` for working
files — both git-ignored). On Windows with Git Bash, `python` may be a broken
store alias — use `py -3` instead; the pipeline was validated that way.

## Authentication (OpenGWAS)

Several MR scripts call the OpenGWAS API and need a personal JWT (free):
register at <https://api.opengwas.io>, then `export OPENGWAS_JWT="<token>"`.
No tokens are stored in this repository (all previously hard-coded JWTs were
removed and replaced by environment-variable reads). Steps whose API products
are cached in `data/derived/` can be run **without** a token (see
PIPELINE.md).

## Large data — not shipped

IHME GBD 2021 raw downloads, full GWAS summary statistics, the GBD 2021 SDI
series, GEO series matrices and the single-cell `.h5ad` files are **not**
redistributed. `data/README.md` gives the exact source, query IDs and
expected local location for each. The GBD SDI csv is required only by
`code/05_gbd/gbd_inequality.py` and `code/05_gbd/gbd_frontier.py`.

## Verification

`VALIDATION.md` records which analyses were re-executed for this package and
compares the outputs number-by-number against the manuscript (meta-analysis,
GBD anchors, decomposition, figure data, intervention scenarios, key-gene
Welch tests, and end-to-end re-runs of the MVMR / replication / drug-target
R pipelines), including two honestly reported discrepancies and the two
components that could not be re-verified locally (SII/CI and frontier
analysis need the GHDx-restricted GBD 2021 SDI file).

## License / citation

Code and original documentation: MIT (see LICENSE). Third-party data remain
subject to their own terms (`data/README.md`). Please cite the manuscript and
this repository (see CITATION.cff).

Archived release: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21618246.svg)](https://doi.org/10.5281/zenodo.21618246)


## MedComm submission update

The earlier EClinicalMedicine submission has formally ended. This update removes the previously listed fourth author from the author list, uses the corrected 204-country GBD framing, and includes the revised Tables S1-S10 and Figures 5-7 prepared for MedComm.
