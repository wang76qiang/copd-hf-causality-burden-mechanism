#!/bin/bash
# Stream-extract 10 SNP rows from FinnGen R12 endpoints (rsids = column 5, comma-separated)
cd "E:/投稿版 ECM1_submission/Medcomm_revised/_shenghua_work"
mkdir -p finngen_extracts
for ep in I9_HEARTFAIL_AND_CHD I9_HEARTFAIL_AND_OVERWEIGHT FG_OTHHEART I9_OTHILLHEART J10_COPD COPD_LATER COPD_EARLY; do
  echo "=== $ep ==="
  url="https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_${ep}.gz"
  curl -s --retry 3 "$url" | zcat 2>/dev/null | awk -F'\t' 'NR==1 {print; next} {n=split($5,a,","); for(i=1;i<=n;i++){ if (a[i] ~ /^rs(114904431|11525583|1246642|1512281|2273500|2395191|34944514|3822479|4488938|9788721)$/) {print; break} }}' > "finngen_extracts/${ep}.tsv"
  wc -l "finngen_extracts/${ep}.tsv"
done
echo DONE
