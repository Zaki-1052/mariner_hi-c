#!/bin/bash
# cluster/scripts/run_summary.sh
# Driver: Phase 7 (Phase 8 in CLAUDE.md) -- lab-meeting summary figures.
# Run from cluster/ (cd is automatic via $0).
set -e
cd "$(dirname "$0")/.."   # cluster/

PYTHON=/opt/homebrew/anaconda3/envs/cluster/bin/python3
SCRIPT=scripts/08_summary_figures.py
LOG=${LOG:-docs/phase8_summary.txt}

mkdir -p "$(dirname "$LOG")"

{
  echo "Summary Figures: $(date)"
  echo "pwd: $(pwd)"
  echo "python: $PYTHON"
  echo ""

  "$PYTHON" "$SCRIPT"

  echo ""
  echo "--- Output inventory ---"
  for dir in dashboard mechanism heatmap; do
    echo "  ${dir}:"
    ls -lh "bap1_late/figures/summary_figures/${dir}/" 2>/dev/null || echo "    (missing)"
  done

  echo ""
  echo "Finished: $(date)"
} 2>&1 | tee "$LOG"
