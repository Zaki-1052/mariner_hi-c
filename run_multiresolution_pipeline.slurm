#!/bin/bash
#SBATCH --job-name=mariner_multiRes
#SBATCH --output=logs/mariner_multiRes_%j.out
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --account=csd940

# ==============================================================================
# Mariner Multi-Resolution Differential Loop Analysis Pipeline
# ==============================================================================
#
# Pipeline:
#   For each resolution (5kb, 10kb, 25kb):
#     1. prep_loops.R      - Merge BEDPE files → consensus positions
#     2. extract_counts.R  - Extract Hi-C from 6 .hic files
#     3. qc-val.R         - QC validation with replicate-aware analysis
#     4. edgeR.R          - Differential analysis with QL GLM (n=3 per group)
#
# Expected runtime: ~6 hours (2 hours × 3 resolutions)
# Expected output: Separate analyses for each resolution
# ==============================================================================

# Create logs directory
mkdir -p /expanse/lustre/projects/csd940/zalibhai/mariner/logs

echo "========================================="
echo "Mariner Multi-Resolution Pipeline"
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
# MULTI-RESOLUTION ANALYSIS
# ==============================================================================

# Array of resolutions to process
RESOLUTIONS=(5000 10000 25000)

for RES in "${RESOLUTIONS[@]}"; do
  RES_KB=$((RES / 1000))

  echo ""
  echo "========================================="
  echo "PROCESSING RESOLUTION: ${RES} bp (${RES_KB} kb)"
  echo "========================================="
  echo "Start time: $(date)"
  echo ""

  # =============================================================================
  # STEP 1: PREPARE LOOPS
  # =============================================================================

  echo "========================================="
  echo "STEP 1: Preparing loops (6 replicates)"
  echo "========================================="
  echo "Script: scripts/prep_loops.R"
  echo "Resolution: ${RES} bp"
  echo ""

  Rscript scripts/prep_loops.R ${RES}

  if [ $? -eq 0 ]; then
    echo ""
    echo "✅ Loop preparation completed for ${RES_KB}kb"
  else
    echo ""
    echo "❌ Loop preparation failed for ${RES_KB}kb with exit code $?"
    echo "Stopping pipeline at resolution ${RES_KB}kb"
    exit 1
  fi

  echo ""

  # =============================================================================
  # STEP 2: EXTRACT COUNTS
  # =============================================================================

  echo "========================================="
  echo "STEP 2: Extracting Hi-C counts (6 files)"
  echo "========================================="
  echo "Script: scripts/extract_counts.R"
  echo "Resolution: ${RES} bp"
  echo ""

  Rscript scripts/extract_counts.R ${RES}

  if [ $? -eq 0 ]; then
    echo ""
    echo "✅ Count extraction completed for ${RES_KB}kb"
  else
    echo ""
    echo "❌ Count extraction failed for ${RES_KB}kb with exit code $?"
    echo "Stopping pipeline at resolution ${RES_KB}kb"
    exit 1
  fi

  echo ""

  # =============================================================================
  # STEP 3: QC VALIDATION
  # =============================================================================

  echo "========================================="
  echo "STEP 3: QC Validation (replicate-aware)"
  echo "========================================="
  echo "Script: scripts/qc-val.R"
  echo "Resolution: ${RES} bp"
  echo ""

  Rscript scripts/qc-val.R ${RES}

  if [ $? -eq 0 ]; then
    echo ""
    echo "✅ QC validation completed for ${RES_KB}kb"
  else
    echo ""
    echo "❌ QC validation failed for ${RES_KB}kb with exit code $?"
    echo "Check outputs/res_${RES_KB}kb/qc_report/ for details"
    echo "Stopping pipeline at resolution ${RES_KB}kb"
    exit 1
  fi

  echo ""

  # =============================================================================
  # STEP 4: DIFFERENTIAL ANALYSIS
  # =============================================================================

  echo "========================================="
  echo "STEP 4: edgeR differential analysis"
  echo "========================================="
  echo "Script: scripts/edgeR.R"
  echo "Resolution: ${RES} bp"
  echo "Method: Quasi-likelihood GLM"
  echo ""

  Rscript scripts/edgeR.R ${RES}

  if [ $? -eq 0 ]; then
    echo ""
    echo "✅ Differential analysis completed for ${RES_KB}kb"
  else
    echo ""
    echo "❌ Differential analysis failed for ${RES_KB}kb with exit code $?"
    echo "Stopping pipeline at resolution ${RES_KB}kb"
    exit 1
  fi

  echo ""
  echo "✓ COMPLETED RESOLUTION: ${RES_KB}kb"
  echo "End time: $(date)"
  echo ""

done

# ==============================================================================
# PIPELINE COMPLETE
# ==============================================================================

echo ""
echo "========================================="
echo "ALL RESOLUTIONS COMPLETE"
echo "========================================="
echo "End time: $(date)"
echo ""
echo "Output directories:"
for RES in "${RESOLUTIONS[@]}"; do
  RES_KB=$((RES / 1000))
  echo "  - outputs/res_${RES_KB}kb/                     (loops & counts)"
  echo "  - outputs/res_${RES_KB}kb/qc_report/           (QC metrics)"
  echo "  - outputs/edgeR_results_res_${RES_KB}kb/       (differential results)"
done
echo ""
echo "Key files to review per resolution:"
echo "  1. outputs/res_{RES}kb/qc_report/qc_report_summary.rds"
echo "     → QC metrics"
echo ""
echo "  2. outputs/edgeR_results_res_{RES}kb/plots/mds_plot.pdf"
echo "     → Sample relationships"
echo ""
echo "  3. outputs/edgeR_results_res_{RES}kb/plots/bcv_plot.pdf"
echo "     → Dispersion estimates"
echo ""
echo "  4. outputs/edgeR_results_res_{RES}kb/plots/ma_plot_primary.pdf"
echo "     → Differential loops visualization"
echo ""
echo "  5. outputs/edgeR_results_res_{RES}kb/primary_analysis/all_results_primary.tsv"
echo "     → Full results table"
echo ""
echo "  6. outputs/edgeR_results_res_{RES}kb/primary_analysis/summary_statistics.txt"
echo "     → Analysis summary"
echo ""
echo "Expected outcomes by resolution:"
echo "  - 5kb:  Most loops, short-range interactions"
echo "  - 10kb: Intermediate coverage"
echo "  - 25kb: Long-range loops, inter-TAD interactions"
echo ""
echo "Next step: Compare results across resolutions"
echo "  Rscript scripts/compare_resolutions.R"
echo "========================================="
