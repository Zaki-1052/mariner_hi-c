#!/bin/bash
# biomodal/non_cg_analysis/bigwigs/submit_remerge.sh
# Submit all 6 non-CG BigWig remerge jobs in parallel.
# Usage: bash submit_remerge.sh

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

for METHOD in wgbs emseq; do
    for GROUP in all ctrl mut; do
        echo "Submitting: ${METHOD} ${GROUP}"
        sbatch "${SCRIPT_DIR}/rerun_merge_bigwigs.sb" "$METHOD" "$GROUP"
    done
done

echo ""
echo "All 6 jobs submitted. Monitor with: squeue -u \$USER"
