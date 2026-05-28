#!/bin/bash
# stripes/stripenn/scripts/run_full_stripenn.sh
# Master orchestration driver for the full Stripenn differential stripe pipeline.
# Chains Stages 0-7 via SLURM --dependency=afterok. Run from login node.
#
# Usage: bash run_full_stripenn.sh [OPTIONS]
#   --skip-stage0       Skip .hic -> .mcool conversion (already done)
#   --timepoint=XXXX    Run only one timepoint (e.g., --timepoint=250831)
#
# Each stage's submit_*.sh captures job IDs and passes them as
# dependencies to the next stage. Within a stage, all jobs run in
# parallel; across stages, SLURM enforces ordering.

CODE_DIR="/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/stripes/stripenn"
DATA_DIR="/expanse/lustre/projects/csd940/zalibhai/stripes/stripenn"
cd "${CODE_DIR}/scripts"

SKIP_STAGE0=false
TP_ARG=""

for arg in "$@"; do
  case "$arg" in
    --skip-stage0)  SKIP_STAGE0=true ;;
    --timepoint=*)  TP_ARG="$arg" ;;
    *)
      echo "Unknown argument: $arg" >&2
      echo "Usage: bash run_full_stripenn.sh [--skip-stage0] [--timepoint=XXXX]" >&2
      exit 1
      ;;
  esac
done

echo "==========================================="
echo "Stripenn Full Pipeline"
echo "==========================================="
echo "Start: $(date)"
echo "Code:  ${CODE_DIR}"
if [ -n "${TP_ARG}" ]; then
  echo "Filter: ${TP_ARG}"
fi
echo "==========================================="
echo ""

if [ "${SKIP_STAGE0}" = true ]; then
  echo "NOTE: Skipping Stage 0 (conversion already done)"
  echo ""
fi

# Stage 0: hictk convert
if [ "${SKIP_STAGE0}" = false ]; then
  echo "=== Stage 0: hictk convert ==="
  JIDS_00=$(bash submit_00_convert.sh ${TP_ARG} 2>/dev/null | paste -sd:)
  if [ -z "${JIDS_00}" ]; then
    echo "ERROR: Stage 0 submission failed." >&2
    exit 1
  fi
  echo "  Job IDs: ${JIDS_00}"
  DEP_01="--dependency=afterok:${JIDS_00}"
else
  DEP_01=""
fi

# Stage 1: stripenn compute
echo ""
echo "=== Stage 1: stripenn compute ==="
JIDS_01=$(bash submit_01_call.sh ${DEP_01} ${TP_ARG} 2>/dev/null | paste -sd:)
if [ -z "${JIDS_01}" ]; then
  echo "ERROR: Stage 1 submission failed." >&2
  exit 1
fi
echo "  Job IDs: ${JIDS_01}"

# Stage 2: build union
echo ""
echo "=== Stage 2: build union ==="
JIDS_02=$(bash submit_02_union.sh --dependency=afterok:${JIDS_01} ${TP_ARG} 2>/dev/null | paste -sd:)
if [ -z "${JIDS_02}" ]; then
  echo "ERROR: Stage 2 submission failed." >&2
  exit 1
fi
echo "  Job IDs: ${JIDS_02}"

# Stage 3: score replicates
echo ""
echo "=== Stage 3: score replicates ==="
JIDS_03=$(bash submit_03_score.sh --dependency=afterok:${JIDS_02} ${TP_ARG} 2>/dev/null | paste -sd:)
if [ -z "${JIDS_03}" ]; then
  echo "ERROR: Stage 3 submission failed." >&2
  exit 1
fi
echo "  Job IDs: ${JIDS_03}"

# Stage 4: edgeR
echo ""
echo "=== Stage 4: edgeR differential ==="
JIDS_04=$(bash submit_04_edgeR.sh --dependency=afterok:${JIDS_03} ${TP_ARG} 2>/dev/null | paste -sd:)
if [ -z "${JIDS_04}" ]; then
  echo "ERROR: Stage 4 submission failed." >&2
  exit 1
fi
echo "  Job IDs: ${JIDS_04}"

# Stage 5: integration
echo ""
echo "=== Stage 5: integration ==="
JIDS_05=$(bash submit_05_integration.sh --dependency=afterok:${JIDS_04} ${TP_ARG} 2>/dev/null | paste -sd:)
if [ -z "${JIDS_05}" ]; then
  echo "ERROR: Stage 5 submission failed." >&2
  exit 1
fi
echo "  Job IDs: ${JIDS_05}"

# Stage 6: cross-resolution comparison
echo ""
echo "=== Stage 6: cross-resolution ==="
JIDS_06=$(bash submit_06_compare.sh --dependency=afterok:${JIDS_05} ${TP_ARG} 2>/dev/null | paste -sd:)
if [ -z "${JIDS_06}" ]; then
  echo "ERROR: Stage 6 submission failed." >&2
  exit 1
fi
echo "  Job IDs: ${JIDS_06}"

# Stage 7: visualization
echo ""
echo "=== Stage 7: visualization (1 job) ==="
JID_07=$(sbatch --parsable --dependency=afterok:${JIDS_06} \
  "${CODE_DIR}/scripts/stripenn_visualizations.sb" 2>/dev/null)
if [ -z "${JID_07}" ]; then
  echo "ERROR: Stage 7 submission failed." >&2
  exit 1
fi
echo "  Job ID: ${JID_07}"

ALL_JIDS="${JIDS_01}:${JIDS_02}:${JIDS_03}:${JIDS_04}:${JIDS_05}:${JIDS_06}:${JID_07}"
if [ "${SKIP_STAGE0}" = false ] && [ -n "${JIDS_00}" ]; then
  ALL_JIDS="${JIDS_00}:${ALL_JIDS}"
fi
N_JOBS=$(echo "${ALL_JIDS}" | tr ':' '\n' | wc -l | tr -d ' ')

echo ""
echo "==========================================="
echo "All stages submitted. ${N_JOBS} jobs total."
echo "Final stage job ID: ${JID_07}"
echo ""
echo "Monitor: squeue -u \$USER"
echo "Logs:    ${DATA_DIR}/logs/"
echo ""
echo "--- Rsync results to local (run after completion) ---"
TP_LABEL="${TP_ARG#--timepoint=}"
if [ -n "${TP_LABEL}" ]; then
  echo "rsync -avz --progress expanse:${DATA_DIR}/outputs/${TP_LABEL}/ outputs/${TP_LABEL}/"
  echo "rsync -avz --progress expanse:${DATA_DIR}/logs/*${TP_LABEL}* logs/"
else
  echo "rsync -avz --progress expanse:${DATA_DIR}/outputs/ outputs/"
  echo "rsync -avz --progress expanse:${DATA_DIR}/logs/ logs/"
fi
echo "==========================================="
