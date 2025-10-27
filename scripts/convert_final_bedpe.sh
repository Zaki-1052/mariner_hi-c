#!/bin/bash
# scripts/convert_final_bedpe.sh
# Convert final_results.tsv files to BEDPE format
# Creates resolution-specific files plus a merged union

set -e  # Exit on error

echo "======================================"
echo "Final Results BEDPE Conversion"
echo "======================================"
echo "Started: $(date)"
echo ""

# Check if Python script exists
SCRIPT="scripts/convert_final_bedpe.py"
if [ ! -f "$SCRIPT" ]; then
    echo "ERROR: Python script not found: $SCRIPT"
    exit 1
fi

# Run Python conversion
python3 "$SCRIPT"

echo ""
echo "======================================"
echo "Completed: $(date)"
echo "======================================"
