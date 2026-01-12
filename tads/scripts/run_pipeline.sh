#!/bin/bash
# run_pipeline.sh
# Master script to run the TADCompare pipeline with dependency tracking
# Usage: ./run_pipeline.sh [step]
#   step: 1, 2, 3, 4, or "all" (default: all)

BASE_DIR="/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/tads"
cd "${BASE_DIR}"

STEP=${1:-all}

submit_and_wait() {
    local script=$1
    local job_name=$2
    
    echo "Submitting ${job_name}..."
    JOB_ID=$(sbatch "${script}" | awk '{print $4}')
    echo "  Job ID: ${JOB_ID}"
    
    # Wait for job to complete
    echo "  Waiting for completion..."
    while squeue -j "${JOB_ID}" 2>/dev/null | grep -q "${JOB_ID}"; do
        sleep 30
    done
    
    # Check exit status from log
    LOG_FILE=$(ls -t logs/${job_name}_*.out 2>/dev/null | head -1)
    if [ -f "${LOG_FILE}" ] && grep -q "✅" "${LOG_FILE}"; then
        echo "  ✅ ${job_name} completed successfully"
        return 0
    else
        echo "  ❌ ${job_name} failed - check ${LOG_FILE}"
        return 1
    fi
}

echo "========================================="
echo "TADCompare Pipeline Runner"
echo "========================================="
echo "Base directory: ${BASE_DIR}"
echo "Step: ${STEP}"
echo ""

case ${STEP} in
    1)
        submit_and_wait scripts/01_extract_matrices.sh "extract"
        ;;
    2)
        submit_and_wait scripts/02_run_tadcompare.sh "tadcompare"
        ;;
    3)
        submit_and_wait scripts/03_run_consensus.sh "consensus"
        ;;
    4)
        submit_and_wait scripts/04_postprocess.sh "postprocess"
        ;;
    all)
        echo "Running full pipeline..."
        echo ""
        
        # Step 1: Extract (skip if already done)
        if [ -d "data/extracted/merged" ] && [ "$(ls data/extracted/merged/*.txt 2>/dev/null | wc -l)" -ge 38 ]; then
            echo "Step 1: Extraction already complete, skipping..."
        else
            submit_and_wait scripts/01_extract_matrices.sh "extract" || exit 1
        fi
        echo ""
        
        # Step 2: TADCompare (skip if already done)
        if [ -f "results/tadcompare/tadcompare_all_boundaries.tsv" ]; then
            echo "Step 2: TADCompare already complete, skipping..."
        else
            submit_and_wait scripts/02_run_tadcompare.sh "tadcompare" || exit 1
        fi
        echo ""
        
        # Step 3: ConsensusTADs
        submit_and_wait scripts/03_run_consensus.sh "consensus" || exit 1
        echo ""
        
        # Step 4: Post-processing
        submit_and_wait scripts/04_postprocess.sh "postprocess" || exit 1
        ;;
    *)
        echo "Unknown step: ${STEP}"
        echo "Usage: ./run_pipeline.sh [1|2|3|4|all]"
        exit 1
        ;;
esac

echo ""
echo "========================================="
echo "Pipeline complete!"
echo "========================================="
echo ""
echo "Final outputs in: ${BASE_DIR}/results/final/"
ls -lh "${BASE_DIR}/results/final/"/*.{tsv,bed,txt} 2>/dev/null
