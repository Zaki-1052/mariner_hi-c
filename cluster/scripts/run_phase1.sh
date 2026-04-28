#!/usr/bin/env bash
# cluster/scripts/run_phase1.sh
# Phase 1 driver: data prep for Popay clustering pipeline (BAP1-KO late timepoint).
# Creates the output directory tree, builds the loop count file (Popay format)
# plus its differential-stats sidecar, and generates the mm10 promoter BED.
# Per project convention we do NOT use `set -euo pipefail` (would abort on
# benign warnings); the underlying Python and R scripts fail loudly themselves.

CLUSTER_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CLUSTER_PY="/opt/homebrew/anaconda3/envs/cluster/bin/python3"
SYS_RSCRIPT="/usr/local/bin/Rscript"

cd "${CLUSTER_ROOT}"

echo "============================================================"
echo "Phase 1: Data preparation"
echo "Cluster dir: ${CLUSTER_ROOT}"
echo "Started:     $(date)"
echo "============================================================"

echo ""
echo "[a] Ensuring output directory tree..."
mkdir -p \
  data \
  scripts \
  bap1_late/cluster3 \
  bap1_late/figures/annotation \
  bap1_late/figures/ChIP_intersect \
  bap1_late/figures/deeptools \
  bap1_late/figures/deeptools_input \
  bap1_late/chromHMM/binarized \
  bap1_late/chromHMM/learned_model \
  bap1_late/chromHMM/anchor_input \
  bap1_late/chromHMM/span_input \
  bap1_late/chromHMM/peak_beds \
  bap1_late/cooltools

echo ""
echo "[b] Step 1.1 — Building loop count file..."
"${CLUSTER_PY}" scripts/01_build_loop_count_file.py
py_status=$?
if [ ${py_status} -ne 0 ]; then
  echo "ERROR: 01_build_loop_count_file.py exited with status ${py_status}" >&2
  exit ${py_status}
fi

echo ""
echo "[c] Step 1.2 — Building mm10 promoter BED..."
"${SYS_RSCRIPT}" scripts/02_build_mm10_gene_annotation.R
r_status=$?
if [ ${r_status} -ne 0 ]; then
  echo "ERROR: 02_build_mm10_gene_annotation.R exited with status ${r_status}" >&2
  exit ${r_status}
fi

echo ""
echo "============================================================"
echo "Phase 1 outputs"
echo "============================================================"
for f in data/late_merged_loop_counts.txt \
         data/late_merged_loop_metadata.tsv \
         data/mm10_knownGene_pp.bed; do
  if [ -s "${f}" ]; then
    rows=$(wc -l < "${f}")
    cols=$(head -1 "${f}" | awk -F'\t' '{print NF}')
    printf "  %-50s  %8d rows  %3d cols\n" "${f}" "${rows}" "${cols}"
  else
    echo "  MISSING or EMPTY: ${f}" >&2
  fi
done

echo ""
echo "Phase 1 finished: $(date)"
