#!/bin/bash
# stripes/stripenn/scripts/submit_01_call.sh
# Wrapper: submits one SLURM job per (timepoint, resolution, merged sample)
# for Stage 1 stripe calling. 8 jobs total (2tp x 2res x 2 merged).
# Prints each JobID (one per line) so downstream stages can collect them
# via:   JIDS=$(bash submit_01_call.sh | paste -sd:)
# Optional trailing arg "--dependency=afterok:<jid[:<jid>...]>" is spliced
# into every sbatch call, for chaining into a prior stage.

CODE_DIR="/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/stripes/stripenn"
DATA_DIR="/expanse/lustre/projects/csd940/zalibhai/stripes/stripenn"
SBATCH_SCRIPT="${CODE_DIR}/scripts/01_call_stripes.sb"
LOG_DIR="${DATA_DIR}/logs"

mkdir -p "${LOG_DIR}"

DEP_ARG=""
if [[ "${1:-}" == --dependency=* ]]; then
  DEP_ARG="$1"
fi

TIMEPOINTS=("250831" "250402")
RESOLUTIONS=("5000" "10000")
SAMPLES=("ctrl_merged" "mut_merged")

for TP in "${TIMEPOINTS[@]}"; do
  for RES in "${RESOLUTIONS[@]}"; do
    RES_KB=$((RES / 1000))
    for S in "${SAMPLES[@]}"; do
      LOG_OUT="${LOG_DIR}/01_call_${TP}_${RES_KB}kb_${S}_%j.out"
      JOBNAME="stripenn_${TP}_${RES_KB}kb_${S}"

      JID=$(sbatch --parsable \
            ${DEP_ARG} \
            --job-name="${JOBNAME}" \
            --output="${LOG_OUT}" \
            "${SBATCH_SCRIPT}" "${TP}" "${RES}" "${S}")
      echo "${JID}"
      echo "submitted ${JOBNAME}  jid=${JID}  log=${LOG_OUT}" >&2
    done
  done
done
