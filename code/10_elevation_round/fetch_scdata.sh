#!/bin/bash
cd "E:/投稿版 ECM1_submission/Medcomm_revised/_shenghua_work"
mkdir -p sc_data
cd sc_data
declare -A files=(
  ["t1_copd_lung.h5ad"]="https://datasets.cellxgene.cziscience.com/6eedf808-fb04-4c0c-a788-26e63bc51d92.h5ad"
  ["t1_hf_fibroblasts.h5ad"]="https://datasets.cellxgene.cziscience.com/9b7c7203-91cd-4e87-aff7-92ed572307dc.h5ad"
  ["t1_hf_endothelial.h5ad"]="https://datasets.cellxgene.cziscience.com/7206e56c-a6f9-45c9-a992-309ce1e7a1fa.h5ad"
  ["t1_hf_macrophages.h5ad"]="https://datasets.cellxgene.cziscience.com/48b2a43b-b04c-41c6-b0a7-69b44502da35.h5ad"
)
for f in "${!files[@]}"; do
  echo "=== $f ==="
  curl -sL -C - --retry 5 --retry-delay 10 -o "$f" "${files[$f]}" && echo "OK $f $(stat -c%s "$f")"
done
echo ALLDONE
