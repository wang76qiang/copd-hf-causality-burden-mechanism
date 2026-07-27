# Data

`data/derived/` ships all small intermediate inputs needed to run the
analysis steps without re-fetching. Large raw data are NOT shipped; obtain
them from the sources below and place them as indicated.

## Shipped (data/derived/)

| Path | Content |
|---|---|
| `mr_discovery/discovery_harmonised_10snp_dat.xlsx` | harmonised 10-SNP COPD(ebi-a-GCST90018807) -> HF(ebi-a-GCST009541) dataset used by modules 01-04 |
| `mvmr/` | clumped instrument sets (`iv_final_*.tsv`), union rsid lists, and per-dataset extracted associations (`assoc_A1A_*.tsv`, `assoc_B2_*.tsv`) for both MVMR models |
| `mediation/` | cached mediation inputs (a-path rds, `ext_med_HF10.tsv`, `tmp_med_*.tsv`, `tmp_hf_*.tsv`, `serpine1_eqtla_region.rds`, `iv_med_TGFB1.rds`, `iv_local_*.rds`) |
| `finngen/` | `fg_header.txt`, `fg_10snps_raw.tsv` (10-SNP rows of FinnGen R12 I9_HEARTFAIL), `fg_r12_phewas_rs7860931.json` (PheWeb API response) |
| `drugtarget/` | per-tissue SERPINE1 cis stats (`serpine1_gene/*.csv`), HF region/lead rds, `hf_serpine1_region.tsv`, `tmp_hf_SERPINE1.tsv`, `eqtlgen_serpine1.json` |
| `bulk/` | investigator-curated expression workbooks: `GSE57148_GSE57338_keygene_expression_values.xlsx` (sheets GSE57148 + GSE57338), `GSE57148_THBS1_two_group_comparison.xlsx`, `GSE57338_ACE_two_group.xlsx`; `gene_list_41_inflammation_fibrosis.csv` (41-gene program) |
| `figure_assets/Fig7F.jpg` | licensed microscopy panel used inside Figure 7 |

Run `bash code/00_prepare_data_downloads.sh` to stage these into the
`data_downloads/` working directory the R scripts expect.

## Not shipped - download instructions

1. **IHME GBD 2021, COPD-attributable heart failure** (module 05)
   GBD Results Tool queries **51d4f54c** and **68ad00ae**
   (cause = Chronic obstructive pulmonary disease; rei/risk = Heart failure;
   measure = Prevalence + YLDs; metric = Number + Rate; 1990-2021; Both sexes;
   all ages + age-standardized + 5-year age groups; all locations).
   Download the 5 csv parts (`IHME-GBD_2021_DATA-*.csv`, ~2.1 M rows) into
   `data/raw/gbd_2021/` or point `GBD_RAW_DIR` at them.
   https://vizhub.healthdata.org/gbd-results/

2. **GEO bulk cohorts** (module 06)
   - GSE57148 (COPD lung, 95 vs 89): https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE57148
   - GSE57338 (failing myocardium, 177 vs 136): https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE57338
   Downloaded automatically by `bulk_limma_reanalysis.R` via GEOquery.

3. **CELLxGENE COPD lung single-cell** (modules 06-07)
   Dataset `8fbed309-d3d4-441b-b3ff-e2dcbcec2d35` ("Anomalous Epithelial
   Variations and Ectopic Inflammatory Response in COPD", 57,918 cells,
   ~849 MB h5ad) -> save as `data/t1_copd_lung.h5ad`.
   https://cellxgene.cziscience.com/

4. **Reichart 2022 heart cell atlas (DCM/ACM)** (modules 06-07)
   CELLxGENE collection e75342a8, per-cell-type h5ad files:
   fibroblasts (ed2b673b), endothelial (1062c0f2), macrophages (9434b020)
   -> `data/t1_hf_fibroblasts.h5ad`, `data/t1_hf_endothelial.h5ad`,
   `data/t1_hf_macrophages.h5ad`.

5. **FinnGen R12 summary statistics** (module 04)
   `finngen_R12_I9_HEARTFAIL.gz` (~818 MB) from
   https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/
   Only the 10 instrument rows are re-used (shipped in `derived/finngen/`);
   the full file is needed only to re-run the extraction.

6. **GWAS summary statistics (EBI GWAS Catalog / OpenGWAS)** (modules 02-03)
   Harmonised full files: GCST90018807 (COPD, Sakaue 2021), GCST009541 (HF,
   Shah 2020 HERMES), GCST90029014 (smoking status), GCST90029007 (BMI),
   GCST90029011 (SBP), GCST007432 (FEV1), GCST007431 (FEV1/FVC),
   GCST009968 (cigarettes/day), GCST90012005 (IL-6), GCST90029070 (CRP).
   Via FTP `https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/` or
   OpenGWAS (ebi-a-* ids). Extracted subsets are shipped in `derived/mvmr/`
   and `derived/mediation/`; OpenGWAS-only datasets (prot-a-2962 TGFB1,
   eqtl-a-ENSG00000106366 SERPINE1 eQTLGen) need a personal `OPENGWAS_JWT`.

7. **eQTL Catalogue GTEx v8 cis summary statistics** (module 08)
   Study QTS000015 per-tissue `*.all.tsv.gz` (+ .tbi) at
   https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/
   `t1_minitabix.py` fetches only the SERPINE1 +/- 1 Mb region over HTTP
   Range; per-tissue extracts are already shipped in
   `derived/drugtarget/serpine1_gene/`.

8. **GBD 2021 Socio-Demographic Index (SDI) 1950-2021** (module 05,
   inequality + frontier only)
   GHDx record: https://ghdx.healthdata.org/record/global-burden-disease-study-2021-gbd-2021-socio-demographic-index-sdi-1950%E2%80%932021
   (free account required; redistribution not permitted). Reformat to a long
   csv `data/derived/gbd2021_sdi.csv` with columns `location_name,year,sdi`.
   Without it, `gbd_inequality.py` and `gbd_frontier.py` cannot run
   (see VALIDATION.md section 3).
