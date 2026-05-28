#!/bin/bash
# ML/cmpts/scripts/run_full_calder2.sh
# Master driver for Track A: CALDER2 subcompartment calling pipeline.
# Chains stages A1→A2→A3+A4 via SLURM --dependency=afterok.
# Run from login node: bash scripts/run_full_calder2.sh [options]
#
# Options:
#   --tp <250402|250831|late|early|both>  Timepoint(s) to process (default: both)
#   --start <A1|A2|A3|A4>  Start from this stage (default: A1).
#                           Prior stages must have completed successfully.
#
# Examples:
#   bash scripts/run_full_calder2.sh                        # full pipeline, both TPs
#   bash scripts/run_full_calder2.sh --tp early             # early only, A1-A4
#   bash scripts/run_full_calder2.sh --tp early --start A2  # early, skip A1
#
# Prerequisites: A0_setup_calder2_env.sh completed successfully.

CODE_DIR="/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/ML/cmpts"
DATA_DIR="/expanse/lustre/projects/csd940/zalibhai/sniper"
LOG_DIR="${CODE_DIR}/logs"

TP_FILTER="both"
START_STAGE="A1"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --tp)    TP_FILTER="$2"; shift 2 ;;
    --start) START_STAGE="$2"; shift 2 ;;
    -h|--help)
      sed -n '2,/^$/p' "$0" | grep '^#' | sed 's/^# \?//'
      exit 0 ;;
    *) echo "ERROR: Unknown option: $1" >&2; exit 1 ;;
  esac
done

case "${TP_FILTER}" in
  250402|late)  TIMEPOINTS=("250402") ;;
  250831|early) TIMEPOINTS=("250831") ;;
  both)         TIMEPOINTS=("250402" "250831") ;;
  *) echo "ERROR: Unknown timepoint: ${TP_FILTER}. Use 250402|250831|late|early|both." >&2; exit 1 ;;
esac

SAMPLES=("ctrl_merged" "mut_merged")

stage_ge() { [[ "A1 A2 A3 A4" == *"$1"*"$2"* ]] || [[ "$1" == "$2" ]]; }

mkdir -p "${LOG_DIR}"

echo "==========================================="
echo "CALDER2 Subcompartment Pipeline (Track A)"
echo "==========================================="
echo "Timepoints:  ${TIMEPOINTS[*]}"
echo "Start stage: ${START_STAGE}"
echo "Code dir:    ${CODE_DIR}"
echo "Started:     $(date)"
echo "==========================================="

JIDS_A1=()
JIDS_A2=()

# ── A1: CALDER2 calling ──────────────────────────────────────────────────────

if [[ "${START_STAGE}" == "A1" ]]; then
  echo ""
  echo "=== A1: CALDER2 subcompartment calling ==="

  for TP in "${TIMEPOINTS[@]}"; do
    for SAMPLE in "${SAMPLES[@]}"; do
      JID=$(sbatch --parsable \
        --job-name="calder2_${TP}_${SAMPLE}" \
        --output="${LOG_DIR}/A1_calder2_${TP}_${SAMPLE}_%j.out" \
        "${CODE_DIR}/scripts/A1_run_calder2.sb" "${TP}" "${SAMPLE}")

      if [ -z "${JID}" ]; then
        echo "ERROR: A1 submission failed for ${TP}/${SAMPLE}." >&2
        exit 1
      fi
      JIDS_A1+=("${JID}")
      echo "  ${TP}/${SAMPLE}  jid=${JID}"
    done
  done
  echo "  A1 total: ${#JIDS_A1[@]} jobs"
fi

# ── A2: Differential subcompartments ──────────────────────────────────────────

if stage_ge "${START_STAGE}" "A2"; then
  echo ""
  echo "=== A2: Differential subcompartment analysis ==="

  A1_DEP=""
  if [[ ${#JIDS_A1[@]} -gt 0 ]]; then
    A1_DEP="--dependency=afterok:$(IFS=:; echo "${JIDS_A1[*]}")"
  fi

  for TP in "${TIMEPOINTS[@]}"; do
    JID=$(sbatch --parsable ${A1_DEP} \
      --job-name="a2_diff_${TP}" \
      --output="${LOG_DIR}/a2_diff_${TP}_%j.out" \
      "${CODE_DIR}/scripts/A2_run.sb" "${TP}")

    if [ -z "${JID}" ]; then
      echo "ERROR: A2 submission failed for ${TP}." >&2
      exit 1
    fi
    JIDS_A2+=("${JID}")
    echo "  ${TP}  jid=${JID}"
  done
fi

# ── A3 + A4: Epigenomic validation + HOMER integration (parallel) ────────────

if stage_ge "${START_STAGE}" "A3"; then
  echo ""
  echo "=== A3+A4: Epigenomic validation + HOMER integration ==="

  A2_DEP=""
  if [[ ${#JIDS_A2[@]} -gt 0 ]]; then
    A2_DEP="--dependency=afterok:$(IFS=:; echo "${JIDS_A2[*]}")"
  fi

  A3_JID=$(sbatch --parsable ${A2_DEP} \
    --output="${LOG_DIR}/a3_epigenomic_%j.out" \
    "${CODE_DIR}/scripts/A3_run.sb")
  echo "  A3 (epigenomic)  jid=${A3_JID}"

  A4_JID=$(sbatch --parsable ${A2_DEP} \
    --output="${LOG_DIR}/a4_homer_integration_%j.out" \
    "${CODE_DIR}/scripts/A4_run.sb")
  echo "  A4 (HOMER)       jid=${A4_JID}"
fi

# ── Summary ───────────────────────────────────────────────────────────────────

N_TOTAL=$(( ${#JIDS_A1[@]} + ${#JIDS_A2[@]} ))
if stage_ge "${START_STAGE}" "A3"; then
  N_TOTAL=$(( N_TOTAL + 2 ))
fi

echo ""
echo "==========================================="
echo "Track A submitted. ${N_TOTAL} jobs total."
echo "Monitor: squeue -u \$USER"
echo "==========================================="
