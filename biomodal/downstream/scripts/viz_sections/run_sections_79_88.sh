#!/bin/bash
# biomodal/downstream/scripts/viz_sections/run_sections_79_88.sh
# Run sections 79-88 (consolidated figure sections) with logging to logs/a4/.
# Must be run from the downstream/ directory.
#
# Usage:
#   cd downstream/
#   bash scripts/viz_sections/run_sections_79_88.sh

SCRIPT_DIR="scripts/viz_sections"
LOG_DIR="logs/a4"
mkdir -p "$LOG_DIR"

SECTIONS=(
  "section_79_methylation_phenotype.R"
  "section_80_k119ub_chromatin_geography.R"
  "section_81_tet_mechanism.R"
  "section_82_neuronal_specificity.R"
  "section_83_mecp2_redistribution.R"
  "section_84_mecp2_aging_methylation_overlay.R"
  "section_85_qc_pca_sex.R"
  "section_86_mecp2_chromatin_reader.R"
  "section_87_mechanism_summary.R"
  "section_88_supplemental_evidence.R"
)

TOTAL=${#SECTIONS[@]}
FAILED=0

for i in "${!SECTIONS[@]}"; do
  script="${SECTIONS[$i]}"
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
