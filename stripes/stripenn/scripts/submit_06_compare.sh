#!/bin/bash
# stripes/stripenn/scripts/submit_06_compare.sh
# Wrapper: submits one SLURM job per timepoint for Stage 6
# cross-resolution comparison. 2 jobs total.

CODE_DIR="/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/stripes/stripenn"
DATA_DIR="/expanse/lustre/projects/csd940/zalibhai/stripes/stripenn"
SBATCH_SCRIPT="${CODE_DIR}/scripts/06_compare_resolutions.sb"
LOG_DIR="${DATA_DIR}/logs"

mkdir -p "${LOG_DIR}"

DEP_ARG=""
if [[ "${1:-}" == --dependency=* ]]; then
  DEP_ARG="$1"
fi

TIMEPOINTS=("250831" "250402")

for TP in "${TIMEPOINTS[@]}"; do
  LOG_OUT="${LOG_DIR}/06_crossres_${TP}_%j.out"
  JOBNAME="crossres_${TP}"

  JID=$(sbatch --parsable \
        ${DEP_ARG} \
        --job-name="${JOBNAME}" \
        --output="${LOG_OUT}" \
        "${SBATCH_SCRIPT}" "${TP}")
  echo "${JID}"
  echo "submitted ${JOBNAME}  jid=${JID}  log=${LOG_OUT}" >&2
done
