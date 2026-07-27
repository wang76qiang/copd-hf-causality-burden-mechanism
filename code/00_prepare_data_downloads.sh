#!/bin/bash
# ============================================================================
# Bootstrap: populate data_downloads/ (the working directory the R pipeline
# scripts read/write) with the small shipped intermediate files.
# Large inputs (full GWAS summary statistics, IHME GBD csv, h5ad) are NOT
# shipped - see data/README.md for download instructions.
# Run from the repository root:  bash code/00_prepare_data_downloads.sh
# ============================================================================
set -e
cd "$(dirname "$0")/.."
mkdir -p data_downloads logs

echo "== module 02 (MVMR) intermediates =="
cp -v data/derived/mvmr/* data_downloads/ 2>/dev/null || true

echo "== module 03 (mediation) intermediates =="
cp -v data/derived/mediation/* data_downloads/ 2>/dev/null || true

echo "== module 04 (FinnGen) intermediates =="
cp -v data/derived/finngen/fg_header.txt data/derived/finngen/fg_10snps_raw.tsv data_downloads/ 2>/dev/null || true

echo "== module 08 (drug target) intermediates =="
mkdir -p data_downloads/serpine1_gene
cp -v data/derived/drugtarget/serpine1_gene/*.csv data_downloads/serpine1_gene/ 2>/dev/null || true
cp -v data/derived/drugtarget/*.rds data/derived/drugtarget/*.tsv data_downloads/ 2>/dev/null || true
mkdir -p data_downloads/tmp
cp -v data/derived/finngen/fg_r12_phewas_rs7860931.json data_downloads/ 2>/dev/null || true
# eQTLGen SERPINE1 json used by t1_a1b_eqtlgen_addendum.R
[ -f data/derived/drugtarget/eqtlgen_serpine1.json ] && \
  cp -v data/derived/drugtarget/eqtlgen_serpine1.json data_downloads/tmp/ || true

echo "Done. data_downloads/ is ready for the pipeline steps in PIPELINE.md."
