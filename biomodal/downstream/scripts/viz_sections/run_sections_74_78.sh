# biomodal/downstream/scripts/viz_sections/run_sections_74_77.sh
# Run sections 74-78 (Jai follow-up analyses) with logging to logs/a3/.
# Section 77 will skip gracefully if young-vs-adult DiffBind data is not available.
# Must be run from the downstream/ directory.
#
# Usage:
#   cd downstream/
#   bash scripts/viz_sections/run_sections_74_77.sh

SCRIPT_DIR="scripts/viz_sections"
LOG_DIR="logs/a3"
mkdir -p "$LOG_DIR"

SECTIONS=(
  "section_74_geneset_overlap_methylation.R"
  "section_75_mecp2_signal_reconciliation.R"
  "section_76_triple_overlap_synapse_chromatin.R"
  "section_77_mecp2_aging_trajectory.R"
  "section_78_stoichiometry_broad_neuronal.R"
)

TOTAL=${#SECTIONS[@]}
FAILED=0

for i in "${!SECTIONS[@]}"; do
  script="${SECTIONS[$i]}"
  # Extract section number for log filename
  secnum=$(echo "$script" | grep -o 'section_[0-9]*' | grep -o '[0-9]*')
  name=$(echo "$script" | sed 's/\.R$//')
  logfile="${LOG_DIR}/${secnum}_$(echo "$name" | sed 's/section_[0-9]*_//').txt"

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
