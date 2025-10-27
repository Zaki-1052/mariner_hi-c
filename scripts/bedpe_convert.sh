#!/bin/bash
# scripts/convert_all_to_bedpe.sh
# Batch convert all edgeR TSV results to BEDPE format for Juicebox
# (Standalone version - no SLURM)

set -e  # Exit on error

echo "======================================"
echo "BEDPE Conversion Batch Script"
echo "Started: $(date)"
echo "======================================"

# Create output directory
BEDPE_DIR="outputs/bedpe"
mkdir -p "$BEDPE_DIR"

# Python conversion script
SCRIPT="scripts/bedpe_convert.py"

# Resolutions to process
RESOLUTIONS=(5kb 10kb 25kb)

echo ""
echo "Processing loop files for all resolutions..."
echo ""

# Counter for tracking progress
total_conversions=0

# Loop through each resolution
for RES in "${RESOLUTIONS[@]}"; do
    echo "======================================"
    echo "Resolution: ${RES}"
    echo "======================================"

    RES_DIR="outputs/edgeR_results_res_${RES}"

    # Check if directory exists
    if [ ! -d "$RES_DIR/primary_analysis" ]; then
        echo "WARNING: Directory not found: $RES_DIR/primary_analysis"
        echo "Skipping ${RES}..."
        continue
    fi

    # =====================================
    # Process significant_loops_fdr05.tsv
    # =====================================
    SIG_FILE="${RES_DIR}/primary_analysis/significant_loops_fdr05.tsv"

    if [ -f "$SIG_FILE" ]; then
        echo ""
        echo "Processing: significant_loops_fdr05.tsv"
        echo "---------------------------------------"

        # All differential (up + down together)
        echo "[1/3] Creating ${RES}_sig_all_differential.bedpe..."
        python3 "$SCRIPT" \
            "$SIG_FILE" \
            "${BEDPE_DIR}/${RES}_sig_all_differential.bedpe" \
            --exclude-unchanged
        ((total_conversions++))

        # Up-regulated only
        echo "[2/3] Creating ${RES}_sig_up.bedpe..."
        python3 "$SCRIPT" \
            "$SIG_FILE" \
            "${BEDPE_DIR}/${RES}_sig_up.bedpe" \
            --direction up_in_mutant
        ((total_conversions++))

        # Down-regulated only
        echo "[3/3] Creating ${RES}_sig_down.bedpe..."
        python3 "$SCRIPT" \
            "$SIG_FILE" \
            "${BEDPE_DIR}/${RES}_sig_down.bedpe" \
            --direction down_in_mutant
        ((total_conversions++))

    else
        echo "WARNING: File not found: $SIG_FILE"
    fi

    # =====================================
    # Process all_results_primary.tsv
    # =====================================
    ALL_FILE="${RES_DIR}/primary_analysis/all_results_primary.tsv"

    if [ -f "$ALL_FILE" ]; then
        echo ""
        echo "Processing: all_results_primary.tsv"
        echo "---------------------------------------"

        # All differential (up + down, FDR<0.05)
        echo "[1/3] Creating ${RES}_all_differential.bedpe..."
        python3 "$SCRIPT" \
            "$ALL_FILE" \
            "${BEDPE_DIR}/${RES}_all_differential.bedpe" \
            --fdr-cutoff 0.03 \
            --exclude-unchanged
        ((total_conversions++))

        # Up-regulated only (FDR<0.05)
        echo "[2/3] Creating ${RES}_up.bedpe..."
        python3 "$SCRIPT" \
            "$ALL_FILE" \
            "${BEDPE_DIR}/${RES}_up.bedpe" \
            --fdr-cutoff 0.03 \
            --direction up_in_mutant
        ((total_conversions++))

        # Down-regulated only (FDR<0.05)
        echo "[3/3] Creating ${RES}_down.bedpe..."
        python3 "$SCRIPT" \
            "$ALL_FILE" \
            "${BEDPE_DIR}/${RES}_down.bedpe" \
            --fdr-cutoff 0.03 \
            --direction down_in_mutant
        ((total_conversions++))

    else
        echo "WARNING: File not found: $ALL_FILE"
    fi

    echo ""
done

echo ""
echo "======================================"
echo "Summary"
echo "======================================"
echo "Total conversions: $total_conversions"
echo ""
echo "Output directory: $BEDPE_DIR"
echo ""
echo "Generated files:"
ls -lh "$BEDPE_DIR"/*.bedpe 2>/dev/null || echo "No BEDPE files found"
echo ""
echo "======================================"
echo "Completed: $(date)"
echo "======================================"
