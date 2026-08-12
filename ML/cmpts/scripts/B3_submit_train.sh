#!/bin/bash
# ML/cmpts/scripts/B3_submit_train.sh
# Wrapper: submits one SLURM job per timepoint for Stage B3.
# 2 jobs total: {250402, 250831}.
# Prints each JobID (one per line) so downstream stages can collect them
# via:   JIDS=$(bash B3_submit_train.sh | paste -sd:)
# Optional leading arg "--dependency=afterok:<jid[:<jid>...]>" is spliced
# into every sbatch call, for chaining after B2.

CODE_DIR="/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/ML/cmpts"
SBATCH_SCRIPT="${CODE_DIR}/scripts/B3_train_sniper.sb"
LOG_DIR="${CODE_DIR}/logs"

mkdir -p "${LOG_DIR}"

DEP_ARG=""
if [[ "${1:-}" == --dependency=* ]]; then
  DEP_ARG="$1"
fi

TIMEPOINTS=("250402" "250831")

for TP in "${TIMEPOINTS[@]}"; do
  LOG_OUT="${LOG_DIR}/b3_train_${TP}_%j.out"
  JOBNAME="b3_train_${TP}"

  JID=$(sbatch --parsable \
        ${DEP_ARG} \
        --job-name="${JOBNAME}" \
        --output="${LOG_OUT}" \
        "${SBATCH_SCRIPT}" "${TP}")

  if [ -z "${JID}" ]; then
    echo "ERROR: sbatch failed for TP=${TP}" >&2
    exit 1
  fi

  echo "${JID}"
  echo "submitted ${JOBNAME}  jid=${JID}  log=${LOG_OUT}" >&2
done
