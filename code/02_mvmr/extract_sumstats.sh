#!/bin/bash
# ============================================================================
# Module 02 - MVMR (shared utility)
# Extract subsets from a (harmonised) GWAS summary-statistics .tsv.gz file.
#   mode=p5e8  : all rows with p_value < 5e-8
#   mode=rsids : rows whose rsid is listed in a file
# Handles both harmonised (hm_*) and raw GWAS-Catalog column layouts.
# Usage: extract_sumstats.sh <file.tsv.gz> <mode> <arg> <out.tsv>
# ============================================================================
# extract_sumstats.sh <file.tsv.gz> <mode> <arg> <out.tsv>
#   mode=p5e8   : all rows with p_value < 5e-8 (arg unused, pass -)
#   mode=rsids  : rows whose rsid is in file <arg> (one rsid per line)
# Handles both harmonised (hm_*) and raw GWAS-Catalog column layouts.
# Output columns: rsid, chr, pos, ea, nea, eaf, beta, se, p, n
set -e
F="$1"; MODE="$2"; ARG="$3"; OUT="$4"
zcat "$F" | awk -v mode="$MODE" -v rsidfile="$ARG" '
BEGIN {
  FS = "\t"; OFS = "\t"
  if (mode == "rsids") { while ((getline line < rsidfile) > 0) want[line] = 1 }
}
NR == 1 {
  for (i = 1; i <= NF; i++) {
    h = $i
    if (h == "hm_rsid" || h == "rsid") c_rsid = i
    if (h == "hm_chrom" || h == "chromosome") c_chr = i
    if (h == "hm_pos" || h == "base_pair_location") c_pos = i
    if (h == "hm_effect_allele" || h == "effect_allele") c_ea = i
    if (h == "hm_other_allele" || h == "other_allele") c_nea = i
    if (h == "hm_effect_allele_frequency" || h == "effect_allele_frequency") c_eaf = i
    if (h == "hm_beta" || h == "beta") c_beta = i
    if (h == "standard_error") c_se = i
    if (h == "p_value") c_p = i
    if (h == "n") c_n = i
  }
  print "rsid","chr","pos","ea","nea","eaf","beta","se","p","n"; next
}
{
  rs = $c_rsid
  if (rs == "" || rs == ".") next
  if (mode == "p5e8") { if (($c_p + 0) >= 5e-8 || $c_p == "") next }
  else if (!(rs in want)) next
  b = $c_beta; if (b == "" || ($c_se+0) <= 0) next
  print rs, $c_chr, $c_pos, $c_ea, $c_nea, $c_eaf, b, $c_se, $c_p, (c_n ? $c_n : "")
}' > "$OUT"
wc -l "$OUT"
