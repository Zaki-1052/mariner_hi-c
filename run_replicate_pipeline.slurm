#!/bin/bash
#SBATCH --job-name=mariner_reps
#SBATCH --output=logs/mariner_reps_%j.out
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --account=csd940

# ==============================================================================
# Mariner Replicate-Aware Differential Loop Analysis Pipeline
# ==============================================================================
#
# Pipeline:
#   1. prep_loops.R      - Merge 6 BEDPE files → consensus positions
#   2. extract_counts.R  - Extract Hi-C from 6 .hic files
#   3. qc-val.R         - QC validation with replicate-aware analysis
#   4. edgeR.R          - Differential analysis with QL GLM (n=3 per group)
#
# Expected runtime: ~15 minutes
# Expected output: Hundreds to thousands of differential loops
# ==============================================================================

# Create logs directory
mkdir -p /expanse/lustre/projects/csd940/zalibhai/mariner/logs

echo "========================================="
echo "Mariner Replicate Pipeline"
echo "========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $HOSTNAME"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: 64GB"
echo "Start time: $(date)"
echo "========================================="
echo ""

# Activate conda environment
source ~/.bashrc
conda activate mariner_env

echo "Environment check:"
which R
R --version | head -n 1
echo ""

# Change to working directory
cd /expanse/lustre/projects/csd940/zalibhai/mariner

# ==============================================================================
# STEP 1: PREPARE LOOPS
# ==============================================================================

echo "========================================="
echo "STEP 1: Preparing loops (6 replicates)"
echo "========================================="
echo "Script: scripts/prep_loops.R"
echo ""

Rscript scripts/prep_loops.R 5000

if [ $? -eq 0 ]; then
    echo ""
    echo "✅ Loop preparation completed"
else
    echo ""
    echo "❌ Loop preparation failed with exit code $?"
    exit 1
fi

echo ""

# ==============================================================================
# STEP 2: EXTRACT COUNTS
# ==============================================================================

echo "========================================="
echo "STEP 2: Extracting Hi-C counts (6 files)"
echo "========================================="
echo "Script: scripts/extract_counts.R" 5000
echo ""

Rscript scripts/extract_counts.R

if [ $? -eq 0 ]; then
    echo ""
    echo "✅ Count extraction completed"
else
    echo ""
    echo "❌ Count extraction failed with exit code $?"
    exit 1
fi

echo ""

# ==============================================================================
# STEP 3: QC VALIDATION
# ==============================================================================

echo "========================================="
echo "STEP 3: QC Validation (replicate-aware)"
echo "========================================="
echo "Script: scripts/qc-val.R"
echo ""

Rscript scripts/qc-val.R 5000

if [ $? -eq 0 ]; then
    echo ""
    echo "✅ QC validation completed"
else
    echo ""
    echo "❌ QC validation failed with exit code $?"
    echo "Check outputs/qc_report/ for details"
    exit 1
fi

echo ""

# ==============================================================================
# STEP 4: DIFFERENTIAL ANALYSIS
# ==============================================================================

echo "========================================="
echo "STEP 4: edgeR differential analysis"
echo "========================================="
echo "Script: scripts/edgeR.R"
echo "Method: Quasi-likelihood GLM"
echo ""

Rscript scripts/edgeR.R 5000

if [ $? -eq 0 ]; then
    echo ""
    echo "✅ Differential analysis completed"
else
    echo ""
    echo "❌ Differential analysis failed with exit code $?"
    exit 1
fi

echo ""

# ==============================================================================
# PIPELINE COMPLETE
# ==============================================================================

echo "========================================="
echo "PIPELINE COMPLETE"
echo "========================================="
echo "End time: $(date)"
echo ""
echo "Output directories:"
echo "  - outputs/full/                   (intermediate files)"
echo "  - outputs/qc_report/              (QC plots & metrics)"
echo "  - outputs/edgeR_results/          (differential results)"
echo ""
echo "Key files to review:"
echo "  1. outputs/qc_report/qc_report_summary.rds"
echo "     → QC metrics"
echo ""
echo "  2. outputs/edgeR_results/plots/mds_plot.pdf"
echo "     → Sample relationships"
echo ""
echo "  3. outputs/edgeR_results/plots/bcv_plot.pdf"
echo "     → Dispersion estimates"
echo ""
echo "  4. outputs/edgeR_results/plots/ma_plot_primary.pdf"
echo "     → Differential loops visualization"
echo ""
echo "  5. outputs/edgeR_results/primary_analysis/all_results_primary.tsv"
echo "     → Full results table"
echo ""
echo "  6. outputs/edgeR_results/primary_analysis/summary_statistics.txt"
echo "     → Analysis summary"
echo ""
echo "Improvement over previous:"
echo "  - Previous (merged, n=1): 2 significant loops"
echo "  - Expected (replicates, n=3): 500-5,000+ loops"
echo "========================================="
