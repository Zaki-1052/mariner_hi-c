# biomodal/downstream/scripts/viz_sections/run_sections_22_26.sh
# Run demethylation/TET/DNMT3A modeling sections (22-26) with logging.
# Must be run from the downstream/ directory.
#
# Usage:
#   cd downstream/
#   bash scripts/viz_sections/run_sections_22_26.sh

SCRIPT_DIR="scripts/viz_sections"
LOG_DIR="logs"
mkdir -p "$LOG_DIR"

SECTIONS=(
  "section_22_demethylation_ratio.R"
  "section_23_baseline_hmc_predictor.R"
  "section_24_dnmt3a_prediction.R"
  "section_25_delta_ratio_models.R"
  "section_26_tet_ko_comparison.R"
)

TOTAL=${#SECTIONS[@]}
FAILED=0

for i in "${!SECTIONS[@]}"; do
  script="${SECTIONS[$i]}"
  name="${script%.R}"
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
echo "COMPLETE: $((TOTAL - FAILED))/$TOTAL sections succeeded"
if [ "$FAILED" -gt 0 ]; then
  echo "  $FAILED section(s) failed — check logs in $LOG_DIR/"
fi
echo "========================================================================"
