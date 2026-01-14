#!/bin/bash
# peaks/scripts/extract_ctcf_motifs.sh
# Extract canonical CTCF motifs from HOMER pre-computed genome-wide motif file
# Run on HPC where the file is located

set -euo pipefail

# =============================================================================
# Configuration
# =============================================================================

# Input: HOMER pre-computed known motifs for mm10 (4GB compressed)
INPUT_FILE="/expanse/lustre/projects/csd940/ctea/homer/homer.KnownMotifs.mm10.191020.bed.gz"

# Output: Canonical CTCF motifs only
OUTPUT_DIR="/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/peaks"
OUTPUT_FILE="${OUTPUT_DIR}/ctcf_motifs_mm10.bed"

# =============================================================================
# Main
# =============================================================================

echo "=== Extracting Canonical CTCF Motifs ==="
echo "Date: $(date)"
echo ""

# Verify input exists
if [[ ! -f "${INPUT_FILE}" ]]; then
    echo "ERROR: Input file not found: ${INPUT_FILE}"
    exit 1
fi

echo "Input: ${INPUT_FILE}"
echo "Output: ${OUTPUT_FILE}"
echo ""

# Extract canonical CTCF(Zf) motifs, excluding satellite variants
# Format: chr  start  end  motif_name  score  strand
echo "Extracting CTCF(Zf) motifs (excluding Satellite variants)..."

zcat "${INPUT_FILE}" \
    | grep "CTCF(Zf)" \
    | grep -v "Satellite" \
    > "${OUTPUT_FILE}"

# Verify output
N_MOTIFS=$(wc -l < "${OUTPUT_FILE}")
echo ""
echo "Extracted ${N_MOTIFS} canonical CTCF motifs"

# Preview
echo ""
echo "Preview (first 5 lines):"
head -5 "${OUTPUT_FILE}"

# Summary stats
echo ""
echo "Chromosome distribution (top 5):"
cut -f1 "${OUTPUT_FILE}" | sort | uniq -c | sort -rn | head -5

echo ""
echo "=== DONE ==="
echo "Output file: ${OUTPUT_FILE}"
