#!/bin/bash
# stripes/stripenn/scripts/submit_03_score.sh
# Wrapper: submits one SLURM job per (timepoint, resolution, replicate)
# for Stage 3 scoring. 24 jobs total (2tp x 2res x 6 replicates).
# Prints each JobID (one per line) for downstream dependency chaining.

CODE_DIR="/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/stripes/stripenn"
DATA_DIR="/expanse/lustre/projects/csd940/zalibhai/stripes/stripenn"
SBATCH_SCRIPT="${CODE_DIR}/scripts/03_score_replicates.sb"
LOG_DIR="${DATA_DIR}/logs"

mkdir -p "${LOG_DIR}"

DEP_ARG=""
TP_FILTER=""
for arg in "$@"; do
  case "$arg" in
    --dependency=*) DEP_ARG="$arg" ;;
    --timepoint=*)  TP_FILTER="${arg#--timepoint=}" ;;
  esac
done

if [ -n "${TP_FILTER}" ]; then
  TIMEPOINTS=("${TP_FILTER}")
else
  TIMEPOINTS=("250831" "250402")
fi
RESOLUTIONS=("5000" "10000")
REPLICATES=("ctrl_M1" "ctrl_M2" "ctrl_M3" "mut_M1" "mut_M2" "mut_M3")

for TP in "${TIMEPOINTS[@]}"; do
  for RES in "${RESOLUTIONS[@]}"; do
    RES_KB=$((RES / 1000))
    for S in "${REPLICATES[@]}"; do
      LOG_OUT="${LOG_DIR}/03_score_${TP}_${RES_KB}kb_${S}_%j.out"
      JOBNAME="score_${TP}_${RES_KB}kb_${S}"

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
