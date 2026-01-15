#!/bin/bash
# scripts/submit_pipeline.sh
# Submit TADCompare pipeline with SLURM job dependencies
# Jobs are submitted all at once; SLURM handles ordering. Exits immediately.
# Always re-runs (no output checking) - use step numbers to run specific steps.
#
# Usage: ./scripts/submit_pipeline.sh [option]
#   Options:
#     all          Full pipeline for both timepoints (default)
#     late         Full pipeline for late timepoint only
#     early        Full pipeline for early timepoint only
#     1-6          Specific step for both timepoints
#
# Examples:
#   ./scripts/submit_pipeline.sh all           # Full pipeline, both timepoints
#   ./scripts/submit_pipeline.sh late          # Full pipeline, late only
#   ./scripts/submit_pipeline.sh 3             # Just step 3 (consensus) for both

# =============================================================================
# Configuration
# =============================================================================

BASE_DIR="/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/tads"
cd "${BASE_DIR}" || exit 1

# Parse arguments
STEP="${1:-all}"

# Determine timepoints
case ${STEP} in
    late)
        TIMEPOINTS=("late")
        RUN_ALL_STEPS=true
        ;;
    early)
        TIMEPOINTS=("early")
        RUN_ALL_STEPS=true
        ;;
    all|1|2|3|4|5|6)
        TIMEPOINTS=("late" "early")
        RUN_ALL_STEPS=false
        if [[ "${STEP}" == "all" ]]; then
            RUN_ALL_STEPS=true
        fi
        ;;
    *)
        echo "Unknown option: ${STEP}"
        echo ""
        echo "Usage: ./scripts/submit_pipeline.sh [all|late|early|1-6]"
        echo ""
        echo "Steps:"
        echo "  1 - Matrix extraction (straw from .hic files)"
        echo "  2 - TADCompare differential analysis"
        echo "  3 - ConsensusTADs robustness assessment"
        echo "  4 - Post-processing (shift distances)"
        echo "  5 - Blacklist filtering"
        echo "  6 - Visualizations"
        echo ""
        echo "Note: Always re-runs jobs (no output checking)."
        exit 1
        ;;
esac

# Create logs directory
mkdir -p logs

echo "========================================="
echo "TADCompare Pipeline - SLURM Submission"
echo "========================================="
echo "Base directory: ${BASE_DIR}"
echo "Timepoints: ${TIMEPOINTS[*]}"
echo "Step: ${STEP}"
echo ""

# =============================================================================
# Helper Functions
# =============================================================================

# Submit job and return job ID
# Usage: submit_job script [timepoint] [dependency_jobid]
submit_job() {
    local script=$1
    local timepoint=$2
    local dep_jobid=$3

    local sbatch_args=""
    if [[ -n "${dep_jobid}" ]]; then
        sbatch_args="--dependency=afterok:${dep_jobid}"
    fi

    if [[ -n "${timepoint}" ]]; then
        sbatch --parsable ${sbatch_args} "${script}" "${timepoint}"
    else
        sbatch --parsable ${sbatch_args} "${script}"
    fi
}

# =============================================================================
# Full Pipeline Submission (with dependencies)
# =============================================================================

if ${RUN_ALL_STEPS}; then
    echo "Submitting full pipeline..."
    echo ""

    # Track final step 4 jobs for step 5 dependency
    declare -a STEP4_JOBS=()

    for tp in "${TIMEPOINTS[@]}"; do
        echo "--- ${tp} timepoint ---"

        # Step 1: Matrix extraction
        JOB1=$(submit_job "scripts/01_extract_matrices.sb" "${tp}")
        echo "  Step 1 (extract):     Job ${JOB1}"

        # Step 2: TADCompare (depends on step 1)
        JOB2=$(submit_job "scripts/02_run_tadcompare.sb" "${tp}" "${JOB1}")
        echo "  Step 2 (tadcompare):  Job ${JOB2}"

        # Step 3: ConsensusTADs (depends on step 2)
        JOB3=$(submit_job "scripts/03_run_consensus.sb" "${tp}" "${JOB2}")
        echo "  Step 3 (consensus):   Job ${JOB3}"

        # Step 4: Post-processing (depends on step 3)
        JOB4=$(submit_job "scripts/04_postprocess.sb" "${tp}" "${JOB3}")
        echo "  Step 4 (postprocess): Job ${JOB4}"

        STEP4_JOBS+=("${JOB4}")

        # Store job IDs for step 6
        eval "JOB4_${tp}=${JOB4}"

        echo ""
    done

    # Step 5: Blacklist filtering (depends on ALL step 4 jobs)
    DEP_STRING=$(IFS=:; echo "${STEP4_JOBS[*]}")
    JOB5=$(sbatch --parsable --dependency=afterok:${DEP_STRING} scripts/05_filter_blacklist.sb)
    echo "--- Blacklist filtering ---"
    echo "  Step 5 (filter):      Job ${JOB5} (depends on: ${DEP_STRING})"
    echo ""

    # Step 6: Visualizations for each timepoint (depends on step 5)
    echo "--- Visualizations ---"
    for tp in "${TIMEPOINTS[@]}"; do
        JOB6=$(submit_job "scripts/05_visualizations.sb" "${tp}" "${JOB5}")
        echo "  Step 6 (viz ${tp}):  Job ${JOB6}"
    done

# =============================================================================
# Single Step Submission
# =============================================================================

elif [[ "${STEP}" == "1" ]]; then
    echo "Submitting step 1 (extraction) only..."
    for tp in "${TIMEPOINTS[@]}"; do
        JOB=$(submit_job "scripts/01_extract_matrices.sb" "${tp}")
        echo "  ${tp}: Job ${JOB}"
    done

elif [[ "${STEP}" == "2" ]]; then
    echo "Submitting step 2 (TADCompare) only..."
    for tp in "${TIMEPOINTS[@]}"; do
        JOB=$(submit_job "scripts/02_run_tadcompare.sb" "${tp}")
        echo "  ${tp}: Job ${JOB}"
    done

elif [[ "${STEP}" == "3" ]]; then
    echo "Submitting step 3 (ConsensusTADs) only..."
    for tp in "${TIMEPOINTS[@]}"; do
        JOB=$(submit_job "scripts/03_run_consensus.sb" "${tp}")
        echo "  ${tp}: Job ${JOB}"
    done

elif [[ "${STEP}" == "4" ]]; then
    echo "Submitting step 4 (post-processing) only..."
    for tp in "${TIMEPOINTS[@]}"; do
        JOB=$(submit_job "scripts/04_postprocess.sb" "${tp}")
        echo "  ${tp}: Job ${JOB}"
    done

elif [[ "${STEP}" == "5" ]]; then
    echo "Submitting step 5 (blacklist filtering) only..."
    JOB=$(submit_job "scripts/05_filter_blacklist.sb")
    echo "  Job ${JOB}"

elif [[ "${STEP}" == "6" ]]; then
    echo "Submitting step 6 (visualizations) only..."
    for tp in "${TIMEPOINTS[@]}"; do
        JOB=$(submit_job "scripts/05_visualizations.sb" "${tp}")
        echo "  ${tp}: Job ${JOB}"
    done
fi

# =============================================================================
# Summary
# =============================================================================

echo ""
echo "========================================="
echo "All jobs submitted!"
echo "========================================="
echo ""
echo "Monitor with:"
echo "  squeue -u \$USER"
echo ""
echo "View logs:"
echo "  tail -f logs/*.out"
echo ""
echo "Cancel all pipeline jobs:"
echo "  scancel -u \$USER --name=extract_hic"
echo "  scancel -u \$USER --name=tadcompare"
echo "  # etc."
echo ""
