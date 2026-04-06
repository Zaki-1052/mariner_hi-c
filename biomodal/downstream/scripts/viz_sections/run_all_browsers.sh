#!/bin/bash
# Submit parallel SLURM jobs for all KEY_GENES genome browser views.
# Each gene runs independently (~1-2 min on HPC).
# After all complete, run once without args for the composite figure.
#
# Usage:
#   bash scripts/viz_sections/run_all_browsers.sh

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SB="${SCRIPT_DIR}/run_browser.sb"

GENES=(Syt1 Zbtb20 Trpm3 Epha3 Mcu Cntnap2 Lpp Dlgap1 Arhgap26 Cdh8)

echo "Submitting ${#GENES[@]} genome browser jobs..."
JOBIDS=()

for gene in "${GENES[@]}"; do
  JOBID=$(sbatch --job-name="browser_${gene}" "${SB}" "${gene}" | awk '{print $4}')
  echo "  ${gene}: job ${JOBID}"
  JOBIDS+=("${JOBID}")
done

# Submit composite figure job that waits for all gene jobs to finish
DEPS=$(IFS=:; echo "${JOBIDS[*]}")
COMP_JOBID=$(sbatch --job-name="browser_composite" --dependency=afterok:${DEPS} "${SB}" | awk '{print $4}')
echo ""
echo "  Composite (after all genes): job ${COMP_JOBID}"
echo "  Depends on: ${DEPS}"
echo ""
echo "Done. Monitor with: squeue -u \$USER -n browser_"
