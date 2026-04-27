#!/usr/bin/env bash
# cluster/scripts/run_phase2.sh
# Phase 2 driver: ChromHMM segmentation for the BAP1-KO cerebellum (late, 250402).
# Wraps 03_chromhmm_segmentation.sh with banner logging mirroring run_phase1.sh.
# Per project convention we do NOT use `set -euo pipefail` (would abort on
# benign Java/awk warnings); the worker checks $? after each long step.

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
WORKER="${REPO_ROOT}/cluster/scripts/03_chromhmm_segmentation.sh"

cd "${REPO_ROOT}"

echo "============================================================"
echo "Phase 2 driver"
echo "Repo root: ${REPO_ROOT}"
echo "Started:   $(date)"
echo "============================================================"

bash "${WORKER}"
worker_status=$?
if [ ${worker_status} -ne 0 ]; then
  echo "ERROR: 03_chromhmm_segmentation.sh exited with status ${worker_status}" >&2
  exit ${worker_status}
fi

echo ""
echo "============================================================"
echo "Phase 2 outputs"
echo "============================================================"
for f in cluster/bap1_late/chromHMM/mm10_standard.txt \
         cluster/bap1_late/chromHMM/cellmarkfiletable.txt \
         cluster/bap1_late/chromHMM/learned_model/cerebellum_late_12_segments.bed \
         cluster/bap1_late/chromHMM/learned_model/emissions_12.txt; do
  if [ -s "${f}" ]; then
    rows=$(wc -l < "${f}")
    printf "  %-70s  %8d rows\n" "${f}" "${rows}"
  else
    echo "  MISSING or EMPTY: ${f}" >&2
  fi
done

echo ""
echo "Phase 2 finished: $(date)"
