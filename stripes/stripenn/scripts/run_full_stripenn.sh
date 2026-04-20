#!/bin/bash
# stripes/stripenn/scripts/run_full_stripenn.sh
# Master orchestration driver for the full Stripenn differential stripe pipeline.
# Chains Stages 0-6 via SLURM --dependency=afterok. Run from login node.
#
# Usage: bash run_full_stripenn.sh [--skip-stage0]
#   --skip-stage0 : Skip .hic -> .mcool conversion (already done)
#
# Each stage's submit_*.sh captures job IDs and passes them as
# dependencies to the next stage. Within a stage, all jobs run in
# parallel; across stages, SLURM enforces ordering.

CODE_DIR="/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/stripes/stripenn"
cd "${CODE_DIR}/scripts"

echo "==========================================="
echo "Stripenn Full Pipeline"
echo "==========================================="
echo "Start: $(date)"
echo "Code:  ${CODE_DIR}"
echo "==========================================="
echo ""

SKIP_STAGE0=false
if [[ "${1:-}" == "--skip-stage0" ]]; then
  SKIP_STAGE0=true
  echo "NOTE: Skipping Stage 0 (conversion already done)"
  echo ""
fi

# Stage 0: hictk convert (16 jobs)
if [ "${SKIP_STAGE0}" = false ]; then
  echo "=== Stage 0: hictk convert (16 jobs) ==="
  JIDS_00=$(bash submit_00_convert.sh 2>/dev/null | paste -sd:)
  if [ -z "${JIDS_00}" ]; then
    echo "ERROR: Stage 0 submission failed." >&2
    exit 1
  fi
  echo "  Job IDs: ${JIDS_00}"
  DEP_01="--dependency=afterok:${JIDS_00}"
else
  DEP_01=""
fi

# Stage 1: stripenn compute (8 jobs)
echo ""
echo "=== Stage 1: stripenn compute (8 jobs) ==="
JIDS_01=$(bash submit_01_call.sh ${DEP_01} 2>/dev/null | paste -sd:)
if [ -z "${JIDS_01}" ]; then
  echo "ERROR: Stage 1 submission failed." >&2
  exit 1
fi
echo "  Job IDs: ${JIDS_01}"

# Stage 2: build union (4 jobs)
echo ""
echo "=== Stage 2: build union (4 jobs) ==="
JIDS_02=$(bash submit_02_union.sh --dependency=afterok:${JIDS_01} 2>/dev/null | paste -sd:)
if [ -z "${JIDS_02}" ]; then
  echo "ERROR: Stage 2 submission failed." >&2
  exit 1
fi
echo "  Job IDs: ${JIDS_02}"

# Stage 3: score replicates (24 jobs)
echo ""
echo "=== Stage 3: score replicates (24 jobs) ==="
JIDS_03=$(bash submit_03_score.sh --dependency=afterok:${JIDS_02} 2>/dev/null | paste -sd:)
if [ -z "${JIDS_03}" ]; then
  echo "ERROR: Stage 3 submission failed." >&2
  exit 1
fi
echo "  Job IDs: ${JIDS_03}"

# Stage 4: edgeR (4 jobs)
echo ""
echo "=== Stage 4: edgeR differential (4 jobs) ==="
JIDS_04=$(bash submit_04_edgeR.sh --dependency=afterok:${JIDS_03} 2>/dev/null | paste -sd:)
if [ -z "${JIDS_04}" ]; then
  echo "ERROR: Stage 4 submission failed." >&2
  exit 1
fi
echo "  Job IDs: ${JIDS_04}"

# Stage 5: integration (4 jobs)
echo ""
echo "=== Stage 5: integration (4 jobs) ==="
JIDS_05=$(bash submit_05_integration.sh --dependency=afterok:${JIDS_04} 2>/dev/null | paste -sd:)
if [ -z "${JIDS_05}" ]; then
  echo "ERROR: Stage 5 submission failed." >&2
  exit 1
fi
echo "  Job IDs: ${JIDS_05}"

# Stage 6: cross-resolution comparison (2 jobs)
echo ""
echo "=== Stage 6: cross-resolution (2 jobs) ==="
JIDS_06=$(bash submit_06_compare.sh --dependency=afterok:${JIDS_05} 2>/dev/null | paste -sd:)
if [ -z "${JIDS_06}" ]; then
  echo "ERROR: Stage 6 submission failed." >&2
  exit 1
fi
echo "  Job IDs: ${JIDS_06}"

echo ""
echo "==========================================="
echo "All stages submitted."
echo "Total jobs: 0+8+4+24+4+4+2 = 46 (excluding Stage 0)"
echo "Final stage job IDs: ${JIDS_06}"
echo ""
echo "Monitor: squeue -u $USER"
echo "Logs:    /expanse/lustre/projects/csd940/zalibhai/stripes/stripenn/logs/"
echo "==========================================="
