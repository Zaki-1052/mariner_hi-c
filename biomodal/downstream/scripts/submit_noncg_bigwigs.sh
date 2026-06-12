#!/bin/bash
# biomodal/downstream/scripts/submit_noncg_bigwigs.sh
# Submit individual SLURM jobs for each CHG/CHH bedGraph -> BigWig conversion
# Usage: bash scripts/submit_noncg_bigwigs.sh
# Run from: biomodal/downstream/

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SLURM_SCRIPT="${SCRIPT_DIR}/bedgraph_to_bigwig_single.sb"

INDIR="/expanse/lustre/projects/csd940/zalibhai/biomodal/modality/exports/bedgraph"

mkdir -p logs/noncg

COUNT=0

for CONTEXT in CHG CHH; do
    for BDG in "${INDIR}/${CONTEXT}_mc"/*.bdg.gz; do
        [[ -f "$BDG" ]] || continue

        SAMPLE=$(basename "$BDG" .bdg.gz)
        JOB_NAME="${CONTEXT}_${SAMPLE}"

        echo "Submitting: ${CONTEXT} ${SAMPLE}"
        sbatch --job-name="$JOB_NAME" "$SLURM_SCRIPT" "$BDG"

        COUNT=$((COUNT + 1))
    done
done

echo ""
echo "Submitted ${COUNT} jobs"
echo "Monitor: squeue -u \$USER"
echo "Logs: logs/noncg/"
