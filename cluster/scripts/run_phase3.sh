#!/usr/bin/env bash
# cluster/scripts/run_phase3.sh
# Driver: Phase 3 — k-means clustering of nonredundant Hi-C loops via Cluster 3.0
# Usage:
#   bash cluster/scripts/run_phase3.sh                # default: k=6, 1%-tile filter
#   bash cluster/scripts/run_phase3.sh 5              # override k
#   bash cluster/scripts/run_phase3.sh 6 0.05         # override k and filter percentile

set -e

cd "$(dirname "$0")/../.."   # repo root

PYTHON=/opt/homebrew/anaconda3/envs/cluster/bin/python3
SCRIPT=cluster/scripts/04_clustering.py
LOG=cluster/phase3.txt

K=${1:-6}
FILTER_PCT=${2:-0.01}

{
  echo "============================================================"
  echo "Phase 3: K-means clustering"
  echo "Repo root: $(pwd)"
  echo "Started:   $(date)"
  echo "k=${K}  filter_pct=${FILTER_PCT}"
  echo "============================================================"

  echo
  echo "[1/2] Generating elbow plot (k=1..14)..."
  "$PYTHON" "$SCRIPT" --elbow-only --filter-pct "$FILTER_PCT"

  echo
  echo "[2/2] Running clustering with k=${K}..."
  "$PYTHON" "$SCRIPT" --k "$K" --filter-pct "$FILTER_PCT"

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
