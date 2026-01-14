#!/bin/bash
# scripts/run_pipeline.sh
# Master script to run the complete TADCompare pipeline at 10kb resolution
# Includes matrix extraction, TADCompare, ConsensusTADs, post-processing,
# blacklist filtering, and visualizations
#
# Usage: ./scripts/run_pipeline.sh [options]
#   Options:
#     all          Run full pipeline for both timepoints (default)
#     late         Run pipeline for late timepoint only
#     early        Run pipeline for early timepoint only
#     extract      Run only matrix extraction for both timepoints
#     1-6          Run specific step for both timepoints
#     --force      Force re-run even if outputs exist
#
# Examples:
#   ./scripts/run_pipeline.sh all           # Full pipeline, both timepoints
#   ./scripts/run_pipeline.sh late          # Full pipeline, late only
#   ./scripts/run_pipeline.sh extract       # Just extraction
#   ./scripts/run_pipeline.sh 5             # Just blacklist filtering
#   ./scripts/run_pipeline.sh all --force   # Force re-run everything

#set -e  # Exit on error

# =============================================================================
# Configuration
# =============================================================================

BASE_DIR="/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/tads"
RESOLUTION=10000
RESOLUTION_KB="10kb"

# Parse arguments
STEP="${1:-all}"
FORCE=false
if [[ "$2" == "--force" ]] || [[ "$1" == "--force" ]]; then
    FORCE=true
    if [[ "$1" == "--force" ]]; then
        STEP="all"
    fi
fi

# Determine which timepoints to process
case ${STEP} in
    late)
        TIMEPOINTS=("late")
        RUN_ALL_STEPS=true
        ;;
    early)
        TIMEPOINTS=("early")
        RUN_ALL_STEPS=true
        ;;
    all|extract|1|2|3|4|5|6)
        TIMEPOINTS=("late" "early")
        RUN_ALL_STEPS=false
        if [[ "${STEP}" == "all" ]]; then
            RUN_ALL_STEPS=true
        fi
        ;;
    *)
        echo "Unknown option: ${STEP}"
        echo ""
        echo "Usage: ./scripts/run_pipeline.sh [all|late|early|extract|1-6] [--force]"
        echo ""
        echo "Steps:"
        echo "  1 - Matrix extraction (straw from .hic files)"
        echo "  2 - TADCompare differential analysis"
        echo "  3 - ConsensusTADs robustness assessment"
        echo "  4 - Post-processing (shift distances)"
        echo "  5 - Blacklist filtering"
        echo "  6 - Visualizations"
        exit 1
        ;;
esac

cd "${BASE_DIR}"

# =============================================================================
# Helper Functions
# =============================================================================

log_header() {
    echo ""
    echo "========================================="
    echo "$1"
    echo "========================================="
}

log_step() {
    echo ""
    echo "--- $1 ---"
}

check_slurm() {
    if command -v sbatch &> /dev/null; then
        return 0
    else
        return 1
    fi
}

# Submit job and wait for completion
submit_and_wait() {
    local script=$1
    local timepoint=$2
    local job_name=$3

    echo "  Submitting ${job_name} for ${timepoint}..."

    if check_slurm; then
        # SLURM submission
        JOB_ID=$(sbatch "${script}" "${timepoint}" | awk '{print $4}')
        echo "    Job ID: ${JOB_ID}"

        # Wait for job to complete
        echo "    Waiting for completion..."
        while squeue -j "${JOB_ID}" 2>/dev/null | grep -q "${JOB_ID}"; do
            sleep 30
        done

        # Check exit status from log
        LOG_FILE=$(ls -t logs/${job_name}_*.out 2>/dev/null | head -1)
        if [ -f "${LOG_FILE}" ] && grep -qE "(✅|Done!|completed successfully)" "${LOG_FILE}"; then
            echo "    ✅ ${job_name} completed successfully"
            return 0
        else
            echo "    ❌ ${job_name} failed - check ${LOG_FILE}"
            return 1
        fi
    else
        # Direct execution (no SLURM)
        echo "    Running directly (no SLURM)..."
        bash "${script}" "${timepoint}"
        return $?
    fi
}

# Run R script directly or via SLURM
run_r_script() {
    local script=$1
    local timepoint=$2
    local job_name=$3
    local slurm_wrapper=$4

    echo "  Running ${job_name} for ${timepoint}..."

    if check_slurm && [ -n "${slurm_wrapper}" ] && [ -f "${slurm_wrapper}" ]; then
        # Use SLURM wrapper if available
        submit_and_wait "${slurm_wrapper}" "${timepoint}" "${job_name}"
        return $?
    else
        # Direct R execution
        Rscript "${script}" "${timepoint}"
        if [ $? -eq 0 ]; then
            echo "    ✅ ${job_name} completed successfully"
            return 0
        else
            echo "    ❌ ${job_name} failed"
            return 1
        fi
    fi
}

# Check if step output exists
output_exists() {
    local timepoint=$1
    local step=$2

    case ${step} in
        1)
            # Check for extracted matrices
            local count=$(ls "data/${timepoint}/extracted/merged/"*.${RESOLUTION_KB}.txt 2>/dev/null | wc -l)
            [ "${count}" -ge 38 ]
            ;;
        2)
            [ -f "results/${timepoint}/tadcompare/tadcompare_all_boundaries.tsv" ]
            ;;
        3)
            [ -f "results/${timepoint}/consensus/tadcompare_with_robustness.tsv" ]
            ;;
        4)
            [ -f "results/${timepoint}/final/tadcompare_final_annotated.tsv" ]
            ;;
        5)
            [ -f "results/${timepoint}/final/tadcompare_final_filtered.tsv" ]
            ;;
        6)
            [ -d "results/visualizations/${timepoint}/overview" ]
            ;;
        *)
            return 1
            ;;
    esac
}

# =============================================================================
# Pipeline Steps
# =============================================================================

run_step_1() {
    local timepoint=$1
    log_step "Step 1: Matrix Extraction (${timepoint})"

    if ! ${FORCE} && output_exists "${timepoint}" 1; then
        echo "  Skipping - extraction already complete"
        return 0
    fi

    mkdir -p "data/${timepoint}/extracted/merged" "data/${timepoint}/extracted/replicates"
    submit_and_wait "scripts/01_extract_matrices.sb" "${timepoint}" "extract"
}

run_step_2() {
    local timepoint=$1
    log_step "Step 2: TADCompare Differential Analysis (${timepoint})"

    if ! ${FORCE} && output_exists "${timepoint}" 2; then
        echo "  Skipping - TADCompare already complete"
        return 0
    fi

    mkdir -p "results/${timepoint}/tadcompare"
    run_r_script "scripts/02_run_tadcompare.R" "${timepoint}" "tadcompare" "scripts/02_run_tadcompare.sb"
}

run_step_3() {
    local timepoint=$1
    log_step "Step 3: ConsensusTADs Robustness (${timepoint})"

    if ! ${FORCE} && output_exists "${timepoint}" 3; then
        echo "  Skipping - ConsensusTADs already complete"
        return 0
    fi

    mkdir -p "results/${timepoint}/consensus"
    run_r_script "scripts/03_run_consensus.R" "${timepoint}" "consensus" "scripts/03_run_consensus.sb"
}

run_step_4() {
    local timepoint=$1
    log_step "Step 4: Post-processing & Shift Distances (${timepoint})"

    if ! ${FORCE} && output_exists "${timepoint}" 4; then
        echo "  Skipping - post-processing already complete"
        return 0
    fi

    mkdir -p "results/${timepoint}/final"
    submit_and_wait "scripts/04_postprocess.sb" "${timepoint}" "postprocess"
}

run_step_5() {
    local timepoint=$1
    log_step "Step 5: Blacklist Filtering (${timepoint})"

    if ! ${FORCE} && output_exists "${timepoint}" 5; then
        echo "  Skipping - blacklist filtering already complete"
        return 0
    fi

    # Run blacklist filtering (R script, processes both timepoints at once)
    if [ "${timepoint}" == "${TIMEPOINTS[0]}" ]; then
        # Only run once for the first timepoint (script handles both)
        run_r_script "scripts/05_filter_blacklist.R" "" "blacklist_filter" "scripts/05_filter_blacklist.sb"
    fi
}

run_step_6() {
    local timepoint=$1
    log_step "Step 6: Visualizations (${timepoint})"

    if ! ${FORCE} && output_exists "${timepoint}" 6; then
        echo "  Skipping - visualizations already complete"
        return 0
    fi

    mkdir -p "results/visualizations/${timepoint}"
    run_r_script "scripts/tad_visualizations.R" "${timepoint}" "visualizations" "scripts/06_visualizations.sb"
}

# =============================================================================
# Main Execution
# =============================================================================

log_header "TADCompare Pipeline - 10kb Resolution"
echo "Base directory: ${BASE_DIR}"
echo "Resolution: ${RESOLUTION_KB}"
echo "Timepoints: ${TIMEPOINTS[*]}"
echo "Force re-run: ${FORCE}"
echo "Step: ${STEP}"

# Create directories
mkdir -p logs

# Determine what to run
if ${RUN_ALL_STEPS}; then
    # Run full pipeline for selected timepoints
    for tp in "${TIMEPOINTS[@]}"; do
        log_header "Processing ${tp} timepoint"

        run_step_1 "${tp}" || exit 1
        run_step_2 "${tp}" || exit 1
        run_step_3 "${tp}" || exit 1
        run_step_4 "${tp}" || exit 1
    done

    # Step 5 processes all timepoints at once
    log_step "Step 5: Blacklist Filtering (all timepoints)"
    run_r_script "scripts/05_filter_blacklist.R" "" "blacklist_filter" "scripts/05_filter_blacklist.sb" || exit 1

    # Step 6 for each timepoint
    for tp in "${TIMEPOINTS[@]}"; do
        run_step_6 "${tp}" || exit 1
    done

elif [[ "${STEP}" == "extract" ]] || [[ "${STEP}" == "1" ]]; then
    for tp in "${TIMEPOINTS[@]}"; do
        run_step_1 "${tp}" || exit 1
    done

elif [[ "${STEP}" == "2" ]]; then
    for tp in "${TIMEPOINTS[@]}"; do
        run_step_2 "${tp}" || exit 1
    done

elif [[ "${STEP}" == "3" ]]; then
    for tp in "${TIMEPOINTS[@]}"; do
        run_step_3 "${tp}" || exit 1
    done

elif [[ "${STEP}" == "4" ]]; then
    for tp in "${TIMEPOINTS[@]}"; do
        run_step_4 "${tp}" || exit 1
    done

elif [[ "${STEP}" == "5" ]]; then
    run_r_script "scripts/05_filter_blacklist.R" "" "blacklist_filter" "scripts/05_filter_blacklist.sb" || exit 1

elif [[ "${STEP}" == "6" ]]; then
    for tp in "${TIMEPOINTS[@]}"; do
        run_step_6 "${tp}" || exit 1
    done
fi

# =============================================================================
# Summary
# =============================================================================

log_header "Pipeline Complete!"

echo "Output locations:"
for tp in "${TIMEPOINTS[@]}"; do
    echo ""
    echo "  ${tp} timepoint:"
    if [ -f "results/${tp}/final/tadcompare_final_filtered.tsv" ]; then
        echo "    Final results: results/${tp}/final/tadcompare_final_filtered.tsv"
        wc -l "results/${tp}/final/tadcompare_final_filtered.tsv" 2>/dev/null | awk '{print "    Boundaries: " $1-1}'
    fi
    if [ -d "results/visualizations/${tp}" ]; then
        echo "    Visualizations: results/visualizations/${tp}/"
    fi
done

echo ""
echo "Done!"
