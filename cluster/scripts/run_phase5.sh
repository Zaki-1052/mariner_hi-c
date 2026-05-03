#!/usr/bin/env bash
# cluster/scripts/run_phase5.sh
# Driver: Phase 5 -- deepTools metagene at loop anchors for BAP1-KO Hi-C clusters.
#
# Builds per-cluster anchor BEDs and runs deepTools bed_pileup with 8 BigWigs
# (4 marks x ctrl/mut). Single combined heatmap at cluster/bap1_late/figures/
# deeptools/histone_anchors/.
#
# Run from cluster/ (cd is automatic via $0).
#
# Env vars:
#   LOG  (default docs/phase5.txt, relative to cluster/)
#
# Usage:
#   bash cluster/scripts/run_phase5.sh
#   LOG=docs/phase5_test.txt bash cluster/scripts/run_phase5.sh

set -e

cd "$(dirname "$0")/.."   # cluster/
[ -n "${CLUSTER_CONF}" ] && source "${CLUSTER_CONF}"

# Cluster env's bin must be on PATH so the `computeMatrix` subprocess called
# by bed_pileup resolves. The cluster env has deepTools 3.5.5 installed,
# but is not the default-active env. PYTHON points to the same env explicitly.
_py="${CLUSTER_PYTHON:-$(command -v python3)}"
CLUSTER_ENV_BIN="${_py%/*}"
export PATH="${CLUSTER_ENV_BIN}:${PATH}"
# Force unbuffered Python so prints stream live through `tee` instead of
# block-buffering during the long subprocess.run('computeMatrix ...') call.
# Without this, the log appears to "hang" between the banner and the next
# Python print for 5-15 minutes while computeMatrix is genuinely running.
export PYTHONUNBUFFERED=1
PYTHON="${CLUSTER_ENV_BIN}/python3"
SCRIPT=scripts/06_deeptools_metagene.py
LOG=${LOG:-docs/phase5.txt}

mkdir -p "$(dirname "$LOG")"

{
  echo "============================================================"
  echo "Phase 5: deepTools metagene at loop anchors"
  echo "Cluster dir: $(pwd)"
  echo "Started:     $(date)"
  echo "log_path:    ${LOG}"
  echo "python:      ${PYTHON}"
  echo "computeMatrix: $(which computeMatrix 2>/dev/null || echo NOT_ON_PATH)"
  echo "plotHeatmap:   $(which plotHeatmap   2>/dev/null || echo NOT_ON_PATH)"
  echo "============================================================"
  echo
  echo "[1/1] Running 06_deeptools_metagene.py..."
  "$PYTHON" -u "$SCRIPT"
  echo
  echo "============================================================"
  echo "Phase 5 outputs"
  echo "============================================================"
  echo "--- per-cluster anchor BEDs ---"
  ls -lh bap1_late/figures/deeptools_input/clust*_anchors.bed 2>/dev/null
  echo
  echo "--- metagene heatmap ---"
  ls -lh bap1_late/figures/deeptools/histone_anchors/ 2>/dev/null
  echo
  echo "Phase 5 finished: $(date)"
} 2>&1 | tee "$LOG"
