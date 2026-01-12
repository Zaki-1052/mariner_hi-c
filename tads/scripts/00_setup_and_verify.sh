#!/bin/bash
# 00_setup_and_verify.sh
# Setup directory structure and verify input files for TADCompare pipeline
# Run this LOCALLY (not via SLURM) before submitting jobs

BASE_DIR="/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/tads"
HIC_DIR="/expanse/lustre/projects/csd940/zalibhai/stripes/StripeCaller/data/hic/250402"

echo "========================================="
echo "TADCompare Pipeline Setup"
echo "========================================="
echo ""

# === Create directory structure ===
echo "Creating directory structure..."

mkdir -p "${BASE_DIR}/scripts"
mkdir -p "${BASE_DIR}/data/extracted/merged"
mkdir -p "${BASE_DIR}/data/extracted/replicates"
mkdir -p "${BASE_DIR}/results/tadcompare"
mkdir -p "${BASE_DIR}/results/consensus"
mkdir -p "${BASE_DIR}/results/final"
mkdir -p "${BASE_DIR}/logs"

echo "  ✓ Directories created"
echo ""

# === Verify input .hic files ===
echo "Verifying input .hic files..."
echo "Source directory: ${HIC_DIR}"
echo ""

MISSING=0

# Merged files
for sample in ctrl_merged mut_merged; do
    if [ -f "${HIC_DIR}/${sample}.hic" ]; then
        size=$(ls -lh "${HIC_DIR}/${sample}.hic" | awk '{print $5}')
        echo "  ✓ ${sample}.hic (${size})"
    else
        echo "  ✗ ${sample}.hic - NOT FOUND"
        MISSING=$((MISSING + 1))
    fi
done

# Replicate files
for sample in ctrl_M1 ctrl_M2 ctrl_M3 mut_M1 mut_M2 mut_M3; do
    if [ -f "${HIC_DIR}/${sample}.hic" ]; then
        size=$(ls -lh "${HIC_DIR}/${sample}.hic" | awk '{print $5}')
        echo "  ✓ ${sample}.hic (${size})"
    else
        echo "  ✗ ${sample}.hic - NOT FOUND"
        MISSING=$((MISSING + 1))
    fi
done

echo ""

# === Check conda environment ===
echo "Checking conda environment..."
source ~/.bashrc

if conda activate mariner_env 2>/dev/null; then
    echo "  ✓ mariner_env activated"
    
    # Check straw
    if command -v straw &> /dev/null; then
        echo "  ✓ straw available"
    else
        echo "  ✗ straw not found in PATH"
        MISSING=$((MISSING + 1))
    fi
    
    # Check R and TADCompare
    if command -v R &> /dev/null; then
        echo "  ✓ R available ($(R --version 2>&1 | head -n1 | cut -d' ' -f3))"
        
        # Quick check for TADCompare
        if R -q -e "library(TADCompare)" &>/dev/null; then
            echo "  ✓ TADCompare package installed"
        else
            echo "  ✗ TADCompare package not installed"
            echo "    Install with: R -e 'BiocManager::install(\"TADCompare\")'"
            MISSING=$((MISSING + 1))
        fi
    else
        echo "  ✗ R not found"
        MISSING=$((MISSING + 1))
    fi
    
    # Check bedtools
    if command -v bedtools &> /dev/null; then
        echo "  ✓ bedtools available"
    else
        echo "  ⚠ bedtools not found (needed for Step 4)"
    fi
    
    conda deactivate
else
    echo "  ✗ Failed to activate mariner_env"
    MISSING=$((MISSING + 1))
fi

echo ""

# === Summary ===
echo "========================================="
if [ ${MISSING} -eq 0 ]; then
    echo "✅ Setup complete - ready to run pipeline"
    echo ""
    echo "Next steps:"
    echo "  1. cd ${BASE_DIR}"
    echo "  2. sbatch scripts/01_extract_matrices.sh"
    echo "  3. (wait for completion)"
    echo "  4. sbatch scripts/02_run_tadcompare.sh"
else
    echo "⚠️  ${MISSING} issues found - resolve before running"
fi
echo "========================================="

# === Print directory tree ===
echo ""
echo "Directory structure:"
find "${BASE_DIR}" -type d | head -20 | sed 's|'"${BASE_DIR}"'|.|g'
