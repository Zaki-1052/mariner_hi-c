#!/usr/bin/env bash
# cluster/scripts/run_oriented_metagene.sh
# Driver: Oriented anchor metagene -- tests K27me3 exterior/interior asymmetry
# at Hi-C loop anchor boundaries per cluster.
# Run from cluster/ (cd is automatic via $0).
#
# Env vars:
#   LOG  (default docs/oriented_metagene.txt, relative to cluster/)
#
# Usage:
#   bash cluster/scripts/run_oriented_metagene.sh
#   LOG=docs/oriented_metagene_test.txt bash cluster/scripts/run_oriented_metagene.sh

set -e

cd "$(dirname "$0")/.."   # cluster/
[ -n "${CLUSTER_CONF}" ] && source "${CLUSTER_CONF}"

_py="${CLUSTER_PYTHON:-$(command -v python3)}"
CLUSTER_ENV_BIN="${_py%/*}"
export PATH="${CLUSTER_ENV_BIN}:${PATH}"
export PYTHONUNBUFFERED=1
PYTHON="${CLUSTER_ENV_BIN}/python3"
SCRIPT=scripts/09_oriented_anchor_metagene.py
LOG=${LOG:-docs/oriented_metagene.txt}

mkdir -p "$(dirname "$LOG")"

{
  echo "============================================================"
  echo "Oriented anchor metagene: K27me3 exterior/interior asymmetry"
  echo "Cluster dir: $(pwd)"
  echo "Started:     $(date)"
  echo "log_path:    ${LOG}"
  echo "python:      ${PYTHON}"
  echo "computeMatrix: $(which computeMatrix 2>/dev/null || echo NOT_ON_PATH)"
  echo "============================================================"
  echo
  echo "[1/1] Running 09_oriented_anchor_metagene.py..."
  "$PYTHON" -u "$SCRIPT"
  echo
  echo "============================================================"
  echo "Oriented metagene outputs"
  echo "============================================================"
  echo "--- oriented anchor BEDs ---"
  ls -lh bap1_late/figures/deeptools_input/clust*_oriented_anchors.bed 2>/dev/null
  echo
  echo "--- oriented metagene heatmap ---"
  ls -lh bap1_late/figures/deeptools/oriented_anchors/ 2>/dev/null
  echo
  echo "Oriented metagene finished: $(date)"
} 2>&1 | tee "$LOG"
