#!/bin/bash
# stripes/stripenn/scripts/submit_00_convert.sh
# Wrapper: submits one SLURM job per (timepoint, sample) for Stage 0 conversion.
# Prints each JobID (one per line) so downstream stages can collect them
# via:   JIDS=$(bash submit_00_convert.sh | paste -sd:)
# Optional trailing arg "--dependency=afterok:<jid[:<jid>...]>" is spliced
# into every sbatch call, for chaining into a prior stage.

#set -e

CODE_DIR="/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/stripes/stripenn"
DATA_DIR="/expanse/lustre/projects/csd940/zalibhai/stripes/stripenn"
SBATCH_SCRIPT="${CODE_DIR}/scripts/00_convert_hic_to_cool.sb"
LOG_DIR="${DATA_DIR}/logs"

mkdir -p "${LOG_DIR}"

# Pass-through dependency flag (optional).
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
SAMPLES=("ctrl_M1" "ctrl_M2" "ctrl_M3" "ctrl_merged" \
         "mut_M1"  "mut_M2"  "mut_M3"  "mut_merged")

for TP in "${TIMEPOINTS[@]}"; do
  for S in "${SAMPLES[@]}"; do
    LOG_OUT="${LOG_DIR}/00_convert_${TP}_${S}_%j.out"
    JOBNAME="hic2cool_${TP}_${S}"

    # sbatch --parsable prints only the JobID to stdout.
    JID=$(sbatch --parsable \
          ${DEP_ARG} \
          --job-name="${JOBNAME}" \
          --output="${LOG_OUT}" \
          "${SBATCH_SCRIPT}" "${TP}" "${S}")
    echo "${JID}"
    echo "submitted ${JOBNAME}  jid=${JID}  log=${LOG_OUT}" >&2
  done
done
