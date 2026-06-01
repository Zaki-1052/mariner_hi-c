#!/bin/bash
# biomodal/downstream/scripts/homer/submit_homer_pipeline.sh
# Master driver for the modular HOMER motif enrichment pipeline.
# Submits BED preparation, 13 parallel HOMER jobs, and a summary job.
#
# Usage (from login node):
#   cd /expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/biomodal/downstream
#   bash scripts/homer/submit_homer_pipeline.sh [--dry-run]

set -u

BASE_DIR="/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/biomodal/downstream"
SCRIPT_DIR="${BASE_DIR}/scripts/homer"
LOG_DIR="${BASE_DIR}/logs"

DRY_RUN=false
for arg in "$@"; do
    case "$arg" in
        --dry-run) DRY_RUN=true ;;
        *)
            echo "Unknown argument: $arg" >&2
            echo "Usage: bash submit_homer_pipeline.sh [--dry-run]" >&2
            exit 1
            ;;
    esac
done

mkdir -p "${LOG_DIR}"

echo "================================================================="
echo "HOMER Motif Enrichment Pipeline — Submission"
echo "================================================================="
echo "Start: $(date)"
echo "Scripts: ${SCRIPT_DIR}"
echo "Logs:    ${LOG_DIR}"
if ${DRY_RUN}; then
    echo "MODE:    DRY RUN (no jobs submitted)"
fi
echo ""

# ─────────────────────────────────────────────────────────────────────
# Stage 0: BED Preparation
# ─────────────────────────────────────────────────────────────────────

echo "--- Stage 0: BED Preparation ---"

PREP_CMD="sbatch --parsable --output=${LOG_DIR}/homer_prepare_beds_%j.out ${SCRIPT_DIR}/prepare_beds.sb"
echo "  ${PREP_CMD}"

if ${DRY_RUN}; then
    PREP_JID="DRY_PREP"
else
    PREP_JID=$(${PREP_CMD})
    if [ $? -ne 0 ]; then
        echo "ERROR: Failed to submit prepare_beds.sb"
        exit 1
    fi
fi
echo "  Job ID: ${PREP_JID}"
echo ""

# ─────────────────────────────────────────────────────────────────────
# Stage 1: HOMER Comparisons (all parallel, depend on Stage 0)
# ─────────────────────────────────────────────────────────────────────

echo "--- Stage 1: HOMER Comparisons (13 parallel jobs) ---"
echo ""

COMPARISONS=(
    "mc_hyper_homer.bed|GENOME|A1_mc_hyper_vs_genome"
    "mc_hypo_homer.bed|GENOME|A2_mc_hypo_vs_genome"
    "hmc_hyper_homer.bed|GENOME|A3_hmc_hyper_vs_genome"
    "hmc_hypo_homer.bed|GENOME|A4_hmc_hypo_vs_genome"
    "coordinated_q4.bed|discordant_q2.bed|A5_coordinated_vs_discordant"
    "k119ub_gained.bed|k119ub_lost.bed|B1_k119ub_gained_vs_lost"
    "k119ub_lost.bed|k119ub_gained.bed|B2_k119ub_lost_vs_gained"
    "diffbind_gained_fdr05.bed|diffbind_lost_fdr05.bed|B3_diffbind_gained_vs_lost"
    "diffbind_lost_fdr05.bed|diffbind_gained_fdr05.bed|B4_diffbind_lost_vs_gained"
    "k27ac_gained_fdr05.bed|k27ac_lost_fdr05.bed|C1_k27ac_gained_vs_lost"
    "k27ac_lost_fdr05.bed|k27ac_gained_fdr05.bed|C2_k27ac_lost_vs_gained"
    "k27ac_gained_fdr05.bed|GENOME|C3_k27ac_gained_vs_genome"
    "k27ac_lost_fdr05.bed|GENOME|C4_k27ac_lost_vs_genome"
)

HOMER_JIDS=""

for COMP in "${COMPARISONS[@]}"; do
    IFS='|' read -r FG BG NAME <<< "${COMP}"

    HOMER_CMD="sbatch --parsable \
        --dependency=afterok:${PREP_JID} \
        --job-name=homer_${NAME} \
        --output=${LOG_DIR}/homer_${NAME}_%j.out \
        ${SCRIPT_DIR}/homer_single_run.sb ${FG} ${BG} ${NAME}"

    echo "  [${NAME}]"
    echo "    ${HOMER_CMD}"

    if ${DRY_RUN}; then
        JID="DRY_${NAME}"
    else
        JID=$(eval ${HOMER_CMD})
        if [ $? -ne 0 ]; then
            echo "    ERROR: Failed to submit ${NAME}"
            continue
        fi
    fi

    echo "    Job ID: ${JID}"
    HOMER_JIDS="${HOMER_JIDS}:${JID}"
    echo ""
done

# ─────────────────────────────────────────────────────────────────────
# Stage 2: Results Summary (depends on ALL Stage 1 jobs)
# ─────────────────────────────────────────────────────────────────────

echo "--- Stage 2: Results Summary ---"

SUM_CMD="sbatch --parsable \
    --dependency=afterok${HOMER_JIDS} \
    --output=${LOG_DIR}/homer_summarize_%j.out \
    ${SCRIPT_DIR}/summarize_results.sb"

echo "  ${SUM_CMD}"

if ${DRY_RUN}; then
    SUM_JID="DRY_SUMMARY"
else
    SUM_JID=$(eval ${SUM_CMD})
    if [ $? -ne 0 ]; then
        echo "  ERROR: Failed to submit summarize_results.sb"
    fi
fi
echo "  Job ID: ${SUM_JID}"
echo ""

# ─────────────────────────────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────────────────────────────

echo "================================================================="
echo "Pipeline submitted:"
echo "  Stage 0 (BED prep):    ${PREP_JID}"
echo "  Stage 1 (13 HOMER):    ${HOMER_JIDS#:}"
echo "  Stage 2 (summary):     ${SUM_JID}"
echo ""
echo "Total jobs: 15 (1 prep + 13 HOMER + 1 summary)"
echo ""
if ! ${DRY_RUN}; then
    echo "Monitor with: squeue -u ${USER} --name=homer_%"
    echo "Cancel all:   scancel --name=homer_%"
fi
echo "================================================================="
