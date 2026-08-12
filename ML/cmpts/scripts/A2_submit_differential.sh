#!/bin/bash
# ML/cmpts/scripts/A2_submit_differential.sh
# Wrapper: submits one SLURM job per timepoint for Stage A2.
# 2 jobs total: {250402, 250831}.
# Prints each JobID (one per line) so downstream stages can collect them
# via:   JIDS=$(bash A2_submit_differential.sh | paste -sd:)
# Optional leading arg "--dependency=afterok:<jid[:<jid>...]>" is spliced
# into every sbatch call, for chaining after A1.

CODE_DIR="/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/ML/cmpts"
SBATCH_SCRIPT="${CODE_DIR}/scripts/A2_run.sb"
LOG_DIR="${CODE_DIR}/logs"

mkdir -p "${LOG_DIR}"

DEP_ARG=""
if [[ "${1:-}" == --dependency=* ]]; then
  DEP_ARG="$1"
fi

TIMEPOINTS=("250402" "250831")

for TP in "${TIMEPOINTS[@]}"; do
  LOG_OUT="${LOG_DIR}/A2_diff_subcomp_${TP}_%j.out"
  JOBNAME="a2_diff_${TP}"

  JID=$(sbatch --parsable \
        ${DEP_ARG} \
        --job-name="${JOBNAME}" \
        --output="${LOG_OUT}" \
        "${SBATCH_SCRIPT}" "${TP}")
  echo "${JID}"
  echo "submitted ${JOBNAME}  jid=${JID}  log=${LOG_OUT}" >&2
done
