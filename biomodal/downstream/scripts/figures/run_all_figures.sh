#!/bin/bash
# biomodal/downstream/scripts/figures/run_all_figures.sh
# Run all consolidated publication figure scripts with logging to logs/figures/.
# Each figure sources _shared_config.R then _figure_config.R and produces a
# multi-panel patchwork figure in plots/figures/.
# Must be run from the downstream/ directory.
#
# NOTE: no 'set -euo pipefail' — benign Gviz/awk/Java warnings would abort the
# pipeline. Failures are tracked via explicit per-script exit-code checks.
#
# Usage:
#   cd downstream/
#   bash scripts/figures/run_all_figures.sh

SCRIPT_DIR="scripts/figures"
LOG_DIR="logs/figures"
mkdir -p "$LOG_DIR"

FIGURES=(
  "figure1_methylation_phenotype.R"
  "figure2_k119ub_chromatin_geography.R"
  "figure3_tet_mechanism.R"
  "figure4_mecp2_redistribution.R"
  "figure5_mecp2_chromatin_reader.R"
  "figure6_neuronal_specificity.R"
  "figure7_mechanism_summary.R"
  "figureS1_supplemental.R"
)

TOTAL=${#FIGURES[@]}
FAILED=0

for i in "${!FIGURES[@]}"; do
  script="${FIGURES[$i]}"
  name=$(echo "$script" | sed 's/\.R$//')
  logfile="${LOG_DIR}/${name}.txt"

  echo ""
  echo "========================================================================"
  echo "[$((i+1))/$TOTAL] Running $script"
  echo "  Log: $logfile"
  echo "========================================================================"
  echo $(date)
  Rscript "${SCRIPT_DIR}/${script}" 2>&1 | tee "$logfile"
  status=${PIPESTATUS[0]}

  if [ "$status" -ne 0 ]; then
    echo "  *** FAILED (exit $status) ***"
    FAILED=$((FAILED + 1))
  else
    echo "  Done."
  fi
done
echo $(date)
echo ""
echo "========================================================================"
echo "COMPLETE: $((TOTAL - FAILED))/$TOTAL figures succeeded"
if [ "$FAILED" -gt 0 ]; then
  echo "  $FAILED figure(s) failed — check logs in $LOG_DIR/"
fi
echo "========================================================================"
