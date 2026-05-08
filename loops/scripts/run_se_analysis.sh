#!/usr/bin/env bash
# loops/scripts/run_se_analysis.sh
#
# Runs the full Superenhancer Hub Analysis pipeline (1a, 1b, 1c)
# for all timepoints and stripe set options.
#
# Usage:
#   bash loops/scripts/run_se_analysis.sh
#
# Must be run from the project root (mariner_hi-c/).

LOGFILE="loops/output/superenhancer_analysis/se_analysis_run.txt"
mkdir -p "$(dirname "$LOGFILE")"

echo "===== Superenhancer Hub Analysis =====" | tee "$LOGFILE"
echo "Started: $(date)" | tee -a "$LOGFILE"
echo "Working directory: $(pwd)" | tee -a "$LOGFILE"
echo "" | tee -a "$LOGFILE"

# --- 1b. Differential Stripes x Superenhancers ---

echo "=== 1b. SE-Stripe Overlap ===" | tee -a "$LOGFILE"

for TP in late early; do
  if [ "$TP" = "late" ]; then
    SETS="allsig highconf"
  else
    SETS="allsig concordant"
  fi
  for SS in $SETS; do
    echo "--- Running: $TP / $SS ---" | tee -a "$LOGFILE"
    Rscript loops/scripts/se_stripe_overlap.R --timepoint "$TP" --stripe-set "$SS" 2>&1 | tee -a "$LOGFILE"
    echo "" | tee -a "$LOGFILE"
  done
done

# --- 1a. SE-DEG Proximity ---

echo "=== 1a. SE-DEG Proximity ===" | tee -a "$LOGFILE"

for TP in late early; do
  echo "--- Running: $TP ---" | tee -a "$LOGFILE"
  Rscript loops/scripts/se_deg_proximity.R --timepoint "$TP" 2>&1 | tee -a "$LOGFILE"
  echo "" | tee -a "$LOGFILE"
done

# --- 1c. SE Hub Gene Classification ---

echo "=== 1c. SE Hub Classification ===" | tee -a "$LOGFILE"

for TP in late early; do
  echo "--- Running: $TP (all stripe sets) ---" | tee -a "$LOGFILE"
  Rscript loops/scripts/se_hub_classification.R --timepoint "$TP" 2>&1 | tee -a "$LOGFILE"
  echo "" | tee -a "$LOGFILE"
done

echo "===== Completed: $(date) =====" | tee -a "$LOGFILE"
echo "Log saved to: $LOGFILE"
