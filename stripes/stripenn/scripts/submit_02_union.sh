#!/bin/bash
# stripes/stripenn/scripts/submit_02_union.sh
# Wrapper: submits one SLURM job per (timepoint, resolution) for Stage 2
# union building. 4 jobs total (2tp x 2res).
# Prints each JobID (one per line) for downstream dependency chaining.

CODE_DIR="/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/stripes/stripenn"
DATA_DIR="/expanse/lustre/projects/csd940/zalibhai/stripes/stripenn"
SBATCH_SCRIPT="${CODE_DIR}/scripts/02_build_union.sb"
LOG_DIR="${DATA_DIR}/logs"

mkdir -p "${LOG_DIR}"

DEP_ARG=""
if [[ "${1:-}" == --dependency=* ]]; then
  DEP_ARG="$1"
fi

TIMEPOINTS=("250831" "250402")
RESOLUTIONS=("5000" "10000")

for TP in "${TIMEPOINTS[@]}"; do
  for RES in "${RESOLUTIONS[@]}"; do
    RES_KB=$((RES / 1000))
    LOG_OUT="${LOG_DIR}/02_union_${TP}_${RES_KB}kb_%j.out"
    JOBNAME="union_${TP}_${RES_KB}kb"

    JID=$(sbatch --parsable \
          ${DEP_ARG} \
          --job-name="${JOBNAME}" \
          --output="${LOG_OUT}" \
          "${SBATCH_SCRIPT}" "${TP}" "${RES}")
    echo "${JID}"
    echo "submitted ${JOBNAME}  jid=${JID}  log=${LOG_OUT}" >&2
  done
done
