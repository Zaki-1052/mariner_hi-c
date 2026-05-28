#!/bin/bash
# ML/cmpts/scripts/A0_setup_calder2_env.sh
# Stage A0: Create conda environment and install CALDER2.
#
# Run INTERACTIVELY on the Expanse login node (not via sbatch).
# Conda env creation ~5 min, R package installation ~15-20 min.
#
# Usage:
#   bash ML/cmpts/scripts/A0_setup_calder2_env.sh

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
CALDER2_SRC="${REPO_ROOT}/ML/cmpts/repos/CALDER2"
INSTALL_SCRIPT="${SCRIPT_DIR}/utils/install_calder2_deps.R"
ENV_NAME="calder2_env"

echo "==========================================="
echo "A0: CALDER2 Environment Setup"
echo "==========================================="
echo "Started:     $(date)"
echo "CALDER2 src: ${CALDER2_SRC}"
echo "==========================================="
echo ""

# ── 1. Verify source files ──

if [ ! -f "${CALDER2_SRC}/DESCRIPTION" ]; then
    echo "ERROR: CALDER2 source not found at ${CALDER2_SRC}" >&2
    exit 1
fi

if [ ! -f "${INSTALL_SCRIPT}" ]; then
    echo "ERROR: R install script not found at ${INSTALL_SCRIPT}" >&2
    exit 1
fi

CALDER_VERSION=$(grep '^Version:' "${CALDER2_SRC}/DESCRIPTION" | awk '{print $2}')
echo "[1/3] CALDER2 source found (v${CALDER_VERSION})"

# ── 2. Create conda environment ──

if conda info --envs 2>/dev/null | grep -q "^${ENV_NAME} "; then
    echo "[2/3] Conda env '${ENV_NAME}' already exists — skipping creation"
    echo "      To rebuild: conda env remove -n ${ENV_NAME} && re-run this script"
else
    echo "[2/3] Creating conda env '${ENV_NAME}'..."
    conda create -n "${ENV_NAME}" \
        r-base=4.2 \
        r-rcpp \
        r-rcpparmadillo \
        compilers \
        -c conda-forge -y
    echo "      Conda env created."
fi

# ── 3. Install R packages ──

echo "[3/3] Installing R packages via ${INSTALL_SCRIPT}..."

eval "$(conda shell.bash hook 2>/dev/null)"
conda activate "${ENV_NAME}"

Rscript "${INSTALL_SCRIPT}" "${CALDER2_SRC}"

echo ""
echo "==========================================="
echo "A0 COMPLETE"
echo "==========================================="
echo "Activate: conda activate ${ENV_NAME}"
echo "Finished: $(date)"
echo "==========================================="
