#!/usr/bin/env bash
# cluster/scripts/run_clust6_subgroups.sh
# Driver: Clust6 sub-group oriented anchor asymmetry (short <800kb vs long >=800kb).
# Runs computeMatrix for ~2,800 anchors only -- ~2-4 min vs ~96 min for full pipeline.
# Preserves existing oriented_anchors/ outputs untouched.
#
# Env vars:
#   LOG  (default docs/clust6_subgroups.txt, relative to cluster/)
#
# Usage:
#   bash cluster/scripts/run_clust6_subgroups.sh
#   LOG=docs/clust6_subgroups_test.txt bash cluster/scripts/run_clust6_subgroups.sh

set -e

cd "$(dirname "$0")/.."   # cluster/
[ -n "${CLUSTER_CONF}" ] && source "${CLUSTER_CONF}"

_py="${CLUSTER_PYTHON:-$(command -v python3)}"
CLUSTER_ENV_BIN="${_py%/*}"
export PATH="${CLUSTER_ENV_BIN}:${PATH}"
export PYTHONUNBUFFERED=1
PYTHON="${CLUSTER_ENV_BIN}/python3"
SCRIPT=scripts/10_clust6_subgroup_asymmetry.py
LOG=${LOG:-docs/clust6_subgroups.txt}

mkdir -p "$(dirname "$LOG")"

{
  echo "============================================================"
  echo "Clust6 sub-group asymmetry: short (<800kb) vs long (>=800kb)"
  echo "Cluster dir: $(pwd)"
  echo "Started:     $(date)"
  echo "log_path:    ${LOG}"
  echo "python:      ${PYTHON}"
  echo "computeMatrix: $(which computeMatrix 2>/dev/null || echo NOT_ON_PATH)"
  echo "============================================================"
  echo
  "$PYTHON" -u "$SCRIPT"
  echo
  echo "============================================================"
  echo "Output inventory"
  echo "============================================================"
  echo "--- oriented anchor BEDs ---"
  ls -lh outputs/bap1_late/figures/deeptools_input/clust6_{short,long}_oriented_anchors.bed 2>/dev/null
  echo
  echo "--- clust6 sub-group outputs ---"
  ls -lh outputs/bap1_late/figures/deeptools/clust6_subgroups/ 2>/dev/null
  echo
  echo "Finished: $(date)"
} 2>&1 | tee "$LOG"
