#!/bin/bash
# ML/cmpts/scripts/A1_submit_calder2.sh
# Wrapper: submits one SLURM job per (timepoint, sample) for Stage A1.
# 4 jobs total: {250402,250831} × {ctrl_merged,mut_merged}.
# Prints each JobID (one per line) so downstream stages can collect them
# via:   JIDS=$(bash A1_submit_calder2.sh | paste -sd:)
# Optional leading arg "--dependency=afterok:<jid[:<jid>...]>" is spliced
# into every sbatch call, for chaining into a prior stage.

CODE_DIR="/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/ML/cmpts"
SBATCH_SCRIPT="${CODE_DIR}/scripts/A1_run_calder2.sb"
LOG_DIR="${CODE_DIR}/logs"

mkdir -p "${LOG_DIR}"

DEP_ARG=""
if [[ "${1:-}" == --dependency=* ]]; then
  DEP_ARG="$1"
fi

TIMEPOINTS=("250402" "250831")
SAMPLES=("ctrl_merged" "mut_merged")

for TP in "${TIMEPOINTS[@]}"; do
  for SAMPLE in "${SAMPLES[@]}"; do
    LOG_OUT="${LOG_DIR}/A1_calder2_${TP}_${SAMPLE}_%j.out"
    JOBNAME="calder2_${TP}_${SAMPLE}"

    JID=$(sbatch --parsable \
          ${DEP_ARG} \
          --job-name="${JOBNAME}" \
          --output="${LOG_OUT}" \
          "${SBATCH_SCRIPT}" "${TP}" "${SAMPLE}")
    echo "${JID}"
    echo "submitted ${JOBNAME}  jid=${JID}  log=${LOG_OUT}" >&2
  done
done
