#!/bin/bash
set -e
cd "$(dirname "$0")/../.."

PYTHON=/opt/homebrew/anaconda3/envs/cluster/bin/python3
SCRIPT=cluster/scripts/08_summary_figures.py
LOG=${LOG:-cluster/phase8_summary.txt}

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
    ls -lh "cluster/bap1_late/figures/summary_figures/${dir}/" 2>/dev/null || echo "    (missing)"
  done

  echo ""
  echo "Finished: $(date)"
} 2>&1 | tee "$LOG"
