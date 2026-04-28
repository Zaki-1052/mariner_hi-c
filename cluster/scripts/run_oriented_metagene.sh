#!/usr/bin/env bash
# cluster/scripts/run_oriented_metagene.sh
# Driver: Oriented anchor metagene -- tests K27me3 exterior/interior asymmetry
# at Hi-C loop anchor boundaries per cluster.
#
# Env vars:
#   LOG  (default cluster/oriented_metagene.txt)
#
# Usage:
#   bash cluster/scripts/run_oriented_metagene.sh
#   LOG=cluster/oriented_metagene_test.txt bash cluster/scripts/run_oriented_metagene.sh

set -e

cd "$(dirname "$0")/../.."   # repo root

CLUSTER_ENV_BIN=/opt/homebrew/anaconda3/envs/cluster/bin
export PATH="${CLUSTER_ENV_BIN}:${PATH}"
export PYTHONUNBUFFERED=1
PYTHON="${CLUSTER_ENV_BIN}/python3"
SCRIPT=cluster/scripts/09_oriented_anchor_metagene.py
LOG=${LOG:-cluster/oriented_metagene.txt}

{
  echo "============================================================"
  echo "Oriented anchor metagene: K27me3 exterior/interior asymmetry"
  echo "Repo root: $(pwd)"
  echo "Started:   $(date)"
  echo "log_path:  ${LOG}"
  echo "python:    ${PYTHON}"
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
  ls -lh cluster/bap1_late/figures/deeptools_input/clust*_oriented_anchors.bed 2>/dev/null
  echo
  echo "--- oriented metagene heatmap ---"
  ls -lh cluster/bap1_late/figures/deeptools/oriented_anchors/ 2>/dev/null
  echo
  echo "Oriented metagene finished: $(date)"
} 2>&1 | tee "$LOG"
