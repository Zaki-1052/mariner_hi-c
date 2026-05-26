#!/bin/bash
# ML/cmpts/scripts/A0_setup_calder2_env.sh
# Stage A0: Create conda environment for CALDER2 subcompartment calling.
#
# Run INTERACTIVELY on the Expanse login node (not via sbatch).
# Takes ~15-30 minutes for conda + R package installation.
#
# CALDER2 (library name: "CALDER", v2.0) is an R package with compiled
# C++ code (RcppArmadillo + OpenMP) and Bioconductor dependencies.
# It accepts .hic files directly via strawr, ships an mm10 reference,
# and is fully mm10-aware via genome='mm10'.
#
# Usage:
#   bash ML/cmpts/scripts/A0_setup_calder2_env.sh
#
# After completion, verify:
#   conda activate calder2_env
#   Rscript -e "library(CALDER); cat('CALDER loaded OK\n')"

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
CALDER2_SRC="${REPO_ROOT}/ML/cmpts/repos/CALDER2"

echo "==========================================="
echo "A0: CALDER2 Environment Setup"
echo "==========================================="
echo "Started:     $(date)"
echo "CALDER2 src: ${CALDER2_SRC}"
echo "==========================================="
echo ""

# ── 1. Verify CALDER2 source is available ──

if [ ! -f "${CALDER2_SRC}/DESCRIPTION" ]; then
    echo "ERROR: CALDER2 source not found at ${CALDER2_SRC}" >&2
    echo "       Expected cloned repo at ML/cmpts/repos/CALDER2/" >&2
    exit 1
fi

CALDER_VERSION=$(grep '^Version:' "${CALDER2_SRC}/DESCRIPTION" | awk '{print $2}')
echo "[1/4] CALDER2 source found (v${CALDER_VERSION})"

# ── 2. Create conda environment ──

ENV_NAME="calder2_env"

if conda info --envs 2>/dev/null | grep -q "^${ENV_NAME} "; then
    echo "[2/4] Conda env '${ENV_NAME}' already exists — skipping creation"
    echo "      To rebuild: conda env remove -n ${ENV_NAME} && re-run this script"
else
    echo "[2/4] Creating conda env '${ENV_NAME}' (R 4.2 + compilation toolchain)..."
    conda create -n "${ENV_NAME}" \
        r-base=4.2 \
        r-rcpp \
        r-rcpparmadillo \
        compilers \
        -c conda-forge -y

    echo "      Conda env created."
fi

# ── 3. Install R dependencies ──

echo "[3/4] Installing R packages (Bioconductor + CRAN + CALDER2 from source)..."
echo "      This takes ~10-20 minutes."

# Activate env for the Rscript call
eval "$(conda shell.bash hook 2>/dev/null)"
conda activate "${ENV_NAME}"

RSCRIPT_BIN="$(which Rscript)"
echo "      Using: ${RSCRIPT_BIN}"

"${RSCRIPT_BIN}" --vanilla -e "
cat('── Installing BiocManager + Bioconductor packages ──\n')
if (!requireNamespace('BiocManager', quietly=TRUE)) {
    install.packages('BiocManager', repos='https://cloud.r-project.org')
}
BiocManager::install(c('GenomicRanges', 'rhdf5'), update=FALSE, ask=FALSE)

cat('\n── Installing CRAN dependencies ──\n')
cran_pkgs <- c(
    'strawr',        # read .hic files (R interface to straw)
    'data.table',    # fast I/O
    'ape',           # phylogenetic tree operations
    'dendextend',    # dendrogram manipulation
    'fitdistrplus',  # distribution fitting
    'igraph',        # graph/tree operations
    'Matrix',        # sparse matrices
    'rARPACK',       # truncated SVD/PCA
    'factoextra',    # PCA utilities
    'fields',        # distance computation (rdist)
    'ggplot2',       # plotting
    'optparse',      # CLI argument parsing
    'R.utils',       # file utilities
    'doParallel',    # parallel foreach backend
    'foreach'        # foreach loop support
)
install.packages(cran_pkgs, repos='https://cloud.r-project.org')

cat('\n── Installing CALDER2 from local source ──\n')
install.packages('${CALDER2_SRC}', repos=NULL, type='source')

cat('\n── Verifying installation ──\n')
library(CALDER)
cat('CALDER package version:', as.character(packageVersion('CALDER')), '\n')
cat('All dependencies loaded successfully.\n')
"

STATUS=$?
if [ ${STATUS} -ne 0 ]; then
    echo "ERROR: R package installation failed (exit=${STATUS})." >&2
    echo "       Check the output above for missing system libraries." >&2
    echo "       Common fix: module load gcc (for C++11 + OpenMP compilation)" >&2
    exit 1
fi

# ── 4. Verification summary ──

echo ""
echo "[4/4] Verification..."

"${RSCRIPT_BIN}" --vanilla -e "
library(CALDER)

ref_file <- system.file('extdata', 'mm10_all_sub_compartments.bed', package='CALDER')
if (nchar(ref_file) == 0 || !file.exists(ref_file)) {
    stop('mm10 reference BED not found — CALDER2 installed incorrectly')
}
ref_lines <- length(readLines(ref_file))
cat('mm10 reference BED:', ref_lines, 'lines at', ref_file, '\n')

cat('\nKey function check: CALDER() exists:', exists('CALDER', mode='function'), '\n')

dep_pkgs <- c('strawr', 'GenomicRanges', 'rhdf5', 'data.table', 'ggplot2',
              'doParallel', 'igraph', 'dendextend', 'rARPACK', 'Matrix')
for (pkg in dep_pkgs) {
    v <- tryCatch(as.character(packageVersion(pkg)), error=function(e) 'MISSING')
    cat(sprintf('  %-15s %s\n', pkg, v))
}
"

VER_STATUS=$?
if [ ${VER_STATUS} -ne 0 ]; then
    echo "ERROR: Verification failed." >&2
    exit 1
fi

echo ""
echo "==========================================="
echo "A0 COMPLETE"
echo "==========================================="
echo "Environment: ${ENV_NAME}"
echo "Activate:    conda activate ${ENV_NAME}"
echo "Verify:      Rscript -e \"library(CALDER); cat('OK\n')\""
echo "Finished:    $(date)"
echo "==========================================="
