#!/usr/bin/env bash
# loops/scripts/run_ctcf_analysis.sh
#
# Runs the CTCF Motif Analysis pipeline (2a, 2b, 2c)
# for both timepoints.
#
# Usage:
#   bash loops/scripts/run_ctcf_analysis.sh
#
# Must be run from the project root (mariner_hi-c/).

LOGFILE="loops/output/ctcf_motif_analysis/ctcf_analysis_run.txt"
mkdir -p "$(dirname "$LOGFILE")"

echo "===== CTCF Motif Analysis =====" | tee "$LOGFILE"
echo "Started: $(date)" | tee -a "$LOGFILE"
echo "Working directory: $(pwd)" | tee -a "$LOGFILE"
echo "" | tee -a "$LOGFILE"

for TP in late early; do
  echo "--- Running: $TP ---" | tee -a "$LOGFILE"
  Rscript loops/scripts/ctcf_motif_analysis.R --timepoint "$TP" 2>&1 | tee -a "$LOGFILE"
  echo "" | tee -a "$LOGFILE"
done

echo "===== Completed: $(date) =====" | tee -a "$LOGFILE"
echo "Log saved to: $LOGFILE"
