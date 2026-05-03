#!/usr/bin/env bash
# cluster/scripts/run_phase4.sh
# Driver: Phase 4 — grouped downstream analyses for BAP1-KO Hi-C loop clusters
#
# Positional args (all optional):
#   $1  ANALYSES  (default 'all')  — comma-separated subset of {4.1..4.8}, or 'all'
#
# Run from cluster/ (cd is automatic via $0).
#
# Env vars:
#   LOG  (default docs/phase4.txt, relative to cluster/)
#
# Usage:
#   bash cluster/scripts/run_phase4.sh                                   # all 8 sub-analyses
#   bash cluster/scripts/run_phase4.sh 4.4                               # KEY result only
#   LOG=docs/phase4_chip.txt bash cluster/scripts/run_phase4.sh 4.3,4.7   # subset + custom log

set -e

cd "$(dirname "$0")/.."   # cluster/
[ -n "${CLUSTER_CONF}" ] && source "${CLUSTER_CONF}"

PYTHON="${CLUSTER_PYTHON:-$(command -v python3)}"
SCRIPT=scripts/05_grouped_analyses.py
LOG=${LOG:-docs/phase4.txt}

mkdir -p "$(dirname "$LOG")"

ANALYSES=${1:-all}

{
  echo "============================================================"
  echo "Phase 4: Grouped downstream analyses"
  echo "Cluster dir: $(pwd)"
  echo "Started:     $(date)"
  echo "analyses:    ${ANALYSES}"
  echo "log_path:    ${LOG}"
  echo "python:      ${PYTHON}"
  echo "chromhmm:    $(which chromhmm 2>/dev/null || echo NOT_ON_PATH)"
  echo "bedtools:    $(which bedtools 2>/dev/null || echo NOT_ON_PATH)"
  echo "============================================================"

  echo
  echo "[1/1] Running 05_grouped_analyses.py..."
  "$PYTHON" "$SCRIPT" --analyses "$ANALYSES"

  echo
  echo "============================================================"
  echo "Phase 4 outputs"
  echo "============================================================"

  echo "--- 4.4 ChromHMM anchor vs span (KEY) ---"
  ls -lh bap1_late/chromHMM/anchor.txt bap1_late/chromHMM/span.txt 2>/dev/null || echo "  (not run)"
  ls -lh bap1_late/chromHMM/anchor.{png,pdf,svg,jpg} 2>/dev/null || true
  ls -lh bap1_late/chromHMM/span.{png,pdf,svg,jpg}   2>/dev/null || true

  echo "--- 4.1 loop size ---"
  ls -lh bap1_late/figures/loop_size/ 2>/dev/null || echo "  (not run)"

  echo "--- 4.8 cluster x differential ---"
  ls -lh bap1_late/figures/cluster_differential_status/ 2>/dev/null || echo "  (not run)"
  ls -lh bap1_late/figures/cluster_differential_status.stats.txt 2>/dev/null || true

  echo "--- 4.5 ChromHMM proportions ---"
  ls -lh bap1_late/figures/chromHMM_anchor/ 2>/dev/null || echo "  (not run)"

  echo "--- 4.6 gene annotation ---"
  ls -lh bap1_late/figures/annotation/clust*_annotation.txt 2>/dev/null | head -10 || echo "  (not run)"

  echo "--- 4.2 loop classification ---"
  ls -lh bap1_late/figures/loop_classification/ 2>/dev/null || echo "  (not run)"

  echo "--- 4.7 DiffBind ---"
  ls -d bap1_late/figures/ChIP_intersect/differential_binding_*  2>/dev/null || echo "  (not run)"
  ls -d bap1_late/figures/ChIP_intersect/ChIP_FC_*               2>/dev/null || true
  ls -lh bap1_late/figures/ChIP_intersect/diffbind_stats_*.txt    2>/dev/null || true

  echo "--- 4.3 anchor ChIP ---"
  ls -lh bap1_late/figures/ChIP_intersect/anchor_ChIP_box/ 2>/dev/null || echo "  (not run)"
  ls -lh bap1_late/figures/ChIP_intersect/anchor_ChIP.stats.txt 2>/dev/null || true

  echo
  echo "Phase 4 finished: $(date)"
} 2>&1 | tee "$LOG"
