# biomodal/downstream/scripts/viz_sections/run_sections_61_63.sh
# Run MeCP2 chromatin sensor / mechanism differentiation sections (61-63) with logging.
# Must be run from the downstream/ directory.
#
# Usage:
#   cd downstream/
#   bash scripts/viz_sections/run_sections_61_63.sh

SCRIPT_DIR="scripts/viz_sections"
LOG_DIR="logs"
mkdir -p "$LOG_DIR"

SECTIONS=(
  "section_61_stoichiometry_mechanism.R"
  "section_62_multifeature_chromatin_regression.R"
  "section_63_mecp2_master_heatmap.R"
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
