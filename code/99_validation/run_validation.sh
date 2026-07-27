#!/bin/bash
# Re-run every check reported in VALIDATION.md.
# Prerequisites: raw IHME GBD 2021 csv files (set GBD_RAW_DIR; see data/README.md).
# Usage:  GBD_RAW_DIR=/path/to/raw bash code/99_validation/run_validation.sh
set -e
cd "$(dirname "$0")/../.."
export REPRO_ROOT="$PWD"
: "${GBD_RAW_DIR:=$PWD/data/raw/gbd_2021}"
OUT="$PWD/validation"; mkdir -p "$OUT/extract" "$OUT/figcheck"

echo "== 1. three-cohort meta-analysis =="
Rscript code/04_replication_meta/three_cohort_meta.R || \
  echo "Rscript unavailable - skipped (see VALIDATION.md 3.3)"

echo "== 2. GBD master extraction (anchors, EAPC, country ranking) =="
REPRO_ROOT="$OUT/extract" GBD_RAW_DIR="$GBD_RAW_DIR" py -3 code/05_gbd/t1_extract_gbd.py

echo "== 3. GBD figure data (CHECK mode vs shipped results/) =="
GBD_RAW_DIR="$GBD_RAW_DIR" OUT_DIR="$OUT/figcheck" CHECK=1 py -3 code/05_gbd/gbd_fig_data.py

echo "== 4. Das Gupta decomposition =="
GBD_RAW_DIR="$GBD_RAW_DIR" OUT_DIR="$OUT" py -3 code/05_gbd/gbd_decomposition.py

echo "== 5. intervention scenarios =="
OUT_DIR="$OUT" py -3 code/05_gbd/gbd_intervention_scenarios.py

echo "== 6. bulk key-gene Welch tests =="
OUT_DIR="$OUT" py -3 code/06_bulk_transcriptomics/bulk_keygene_welch.py

echo "== 7. SII/CI and frontier: require data/derived/gbd2021_sdi.csv (see data/README.md source 8) =="
[ -f data/derived/gbd2021_sdi.csv ] && {
  GBD_RAW_DIR="$GBD_RAW_DIR" OUT_DIR="$OUT" py -3 code/05_gbd/gbd_inequality.py
  GBD_RAW_DIR="$GBD_RAW_DIR" OUT_DIR="$OUT" py -3 code/05_gbd/gbd_frontier.py
} || echo "SDI csv missing - skipped (expected; see VALIDATION.md section 3)"

echo "Validation run complete. Compare $OUT against VALIDATION.md."
