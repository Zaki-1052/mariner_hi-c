#!/bin/bash
# scripts/generate_final_results.sh
# Generate final_results.tsv from existing all_results_primary.tsv files
# Applies filters: |logFC| > 0.3, FDR < 0.03

set -e  # Exit on error

echo "======================================"
echo "Generate Final Results"
echo "======================================"
echo "Started: $(date)"
echo ""

# Check if Python script exists
SCRIPT="scripts/generate_final_results.py"
if [ ! -f "$SCRIPT" ]; then
    echo "ERROR: Python script not found: $SCRIPT"
    exit 1
fi

# Run Python script
python3 "$SCRIPT"

echo ""
echo "======================================"
echo "Completed: $(date)"
echo "======================================"
