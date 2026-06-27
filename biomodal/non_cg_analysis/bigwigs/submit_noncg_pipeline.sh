#!/bin/bash
# biomodal/non_cg_analysis/bigwigs/submit_noncg_pipeline.sh
#
# Submit the full non-CG BigWig pipeline with SLURM dependency chaining.
# WGBS and EM-seq paths run in parallel; within each path, steps are sequential.
#
# Usage: bash submit_noncg_pipeline.sh
#        (run from the bigwigs/ directory on Expanse)

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
LOGDIR="${SCRIPT_DIR}/logs"
OUTBASE="/expanse/lustre/projects/csd940/zalibhai/noncg_bigwigs"

mkdir -p "$LOGDIR" "${OUTBASE}/bedgraphs/wgbs" "${OUTBASE}/bedgraphs/emseq" \
         "${OUTBASE}/bigwigs/wgbs" "${OUTBASE}/bigwigs/emseq" "${OUTBASE}/tmp"

echo "========================================="
echo "Non-CG BigWig Pipeline"
echo "========================================="
echo "Script dir: $SCRIPT_DIR"
echo "Output dir: $OUTBASE"
echo "Log dir:    $LOGDIR"
echo ""

# Step 1: bedGraph creation (WGBS and EM-seq in parallel)
echo "--- Step 1: bedGraph creation ---"

WGBS_BDG_JID=$(sbatch --parsable --array=0-3 "${SCRIPT_DIR}/wgbs_noncg_to_bedgraph.sb")
echo "  WGBS bedGraph:  job $WGBS_BDG_JID (array 0-3)"

EMSEQ_BDG_JID=$(sbatch --parsable --array=0-3 "${SCRIPT_DIR}/emseq_noncg_to_bedgraph.sb")
echo "  EM-seq bedGraph: job $EMSEQ_BDG_JID (array 0-3)"

# Step 2: BigWig conversion (after respective bedGraph step)
echo ""
echo "--- Step 2: BigWig conversion ---"

WGBS_BW_JID=$(sbatch --parsable --dependency=afterok:${WGBS_BDG_JID} \
    "${SCRIPT_DIR}/noncg_bedgraph_to_bigwig.sb" wgbs)
echo "  WGBS BigWig:  job $WGBS_BW_JID (after $WGBS_BDG_JID)"

EMSEQ_BW_JID=$(sbatch --parsable --dependency=afterok:${EMSEQ_BDG_JID} \
    "${SCRIPT_DIR}/noncg_bedgraph_to_bigwig.sb" emseq)
echo "  EM-seq BigWig: job $EMSEQ_BW_JID (after $EMSEQ_BDG_JID)"

# Step 3: Merge BigWigs (after respective BigWig step)
echo ""
echo "--- Step 3: Merge BigWigs ---"

WGBS_MERGE_JID=$(sbatch --parsable --dependency=afterok:${WGBS_BW_JID} \
    "${SCRIPT_DIR}/noncg_merge_bigwigs.sb" wgbs)
echo "  WGBS merge:  job $WGBS_MERGE_JID (after $WGBS_BW_JID)"

EMSEQ_MERGE_JID=$(sbatch --parsable --dependency=afterok:${EMSEQ_BW_JID} \
    "${SCRIPT_DIR}/noncg_merge_bigwigs.sb" emseq)
echo "  EM-seq merge: job $EMSEQ_MERGE_JID (after $EMSEQ_BW_JID)"

echo ""
echo "========================================="
echo "All jobs submitted. Dependency chain:"
echo ""
echo "  WGBS:   $WGBS_BDG_JID -> $WGBS_BW_JID -> $WGBS_MERGE_JID"
echo "  EM-seq: $EMSEQ_BDG_JID -> $EMSEQ_BW_JID -> $EMSEQ_MERGE_JID"
echo ""
echo "Monitor: squeue -u \$USER"
echo "========================================="
