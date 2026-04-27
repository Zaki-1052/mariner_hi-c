#!/usr/bin/env bash
# cluster/scripts/run_phase3.sh
# Driver: Phase 3 — k-means clustering of nonredundant Hi-C loops via Cluster 3.0
#
# Positional args:
#   $1  K            (default 6)        — number of k-means clusters
#   $2  FILTER_PCT   (default 0.01)     — per-resolution percentile floor for ctrl_merge
#   $3  MIN_RATIO    (default unset)    — drop loops with mut/ctrl below this
#   $4  MAX_RATIO    (default unset)    — drop loops with mut/ctrl above this
#
# Env vars:
#   LOG  (default cluster/phase3.txt)   — full output log path
#
# Usage:
#   bash cluster/scripts/run_phase3.sh                                # k=6, 1%-tile filter, no ratio bounds
#   bash cluster/scripts/run_phase3.sh 6 0.01                          # explicit k+filter, no bounds
#   bash cluster/scripts/run_phase3.sh 6 0.01 0.333 3.0                 # add symmetric ratio bounds
#   LOG=cluster/phase3_v2.txt bash cluster/scripts/run_phase3.sh 6 0.01 0.333 3.0   # custom log path

set -e

cd "$(dirname "$0")/../.."   # repo root

PYTHON=/opt/homebrew/anaconda3/envs/cluster/bin/python3
SCRIPT=cluster/scripts/04_clustering.py
LOG=${LOG:-cluster/phase3.txt}

K=${1:-6}
FILTER_PCT=${2:-0.01}
MIN_RATIO=${3:-}
MAX_RATIO=${4:-}

# Build optional ratio-bound flags
RATIO_ARGS=""
[ -n "$MIN_RATIO" ] && RATIO_ARGS+=" --min-ratio $MIN_RATIO"
[ -n "$MAX_RATIO" ] && RATIO_ARGS+=" --max-ratio $MAX_RATIO"

{
  echo "============================================================"
  echo "Phase 3: K-means clustering"
  echo "Repo root: $(pwd)"
  echo "Started:   $(date)"
  echo "k=${K}  filter_pct=${FILTER_PCT}  min_ratio=${MIN_RATIO:-none}  max_ratio=${MAX_RATIO:-none}"
  echo "log_path:  ${LOG}"
  echo "============================================================"

  echo
  echo "[1/2] Generating elbow plot (k=1..14)..."
  "$PYTHON" "$SCRIPT" --elbow-only --filter-pct "$FILTER_PCT" $RATIO_ARGS

  echo
  echo "[2/2] Running clustering with k=${K}..."
  "$PYTHON" "$SCRIPT" --k "$K" --filter-pct "$FILTER_PCT" $RATIO_ARGS

  echo
  echo "============================================================"
  echo "Phase 3 outputs"
  echo "============================================================"
  ls -lh cluster/bap1_late/cluster3/elbow_plot/ 2>/dev/null || true
  ls -lh "cluster/bap1_late/cluster3/k-${K}/data/combined-clusters.txt" 2>/dev/null || true
  for sub in heatmap lineplot stripplot boxplot; do
    echo "--- ${sub}/ ---"
    ls -lh "cluster/bap1_late/cluster3/k-${K}/figures/${sub}/" 2>/dev/null || true
  done

  echo "Phase 3 finished: $(date)"
} 2>&1 | tee "$LOG"
