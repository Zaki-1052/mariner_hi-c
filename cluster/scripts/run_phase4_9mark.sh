#!/usr/bin/env bash
# cluster/scripts/run_phase4_9mark.sh
# Driver: rerun Phase 4.4 (anchor-vs-span) and 4.5 (proportions) for the 9-mark
# ChromHMM model. Reads from chromHMM_9mark_{intersect,union}/ outputs.
# Does NOT touch existing 5-mark/12-state outputs.
#
# Usage:
#   bash cluster/scripts/run_phase4_9mark.sh                       # all 4 models
#   bash cluster/scripts/run_phase4_9mark.sh intersect              # intersect × 15+18
#   bash cluster/scripts/run_phase4_9mark.sh union 15               # union × 15 only
#   bash cluster/scripts/run_phase4_9mark.sh intersect 18           # intersect × 18 only
#
# Args:
#   $1  CONSENSUS  {intersect, union, both}  (default: both)
#   $2  STATES     {15, 18, both}            (default: both)
#
# Env vars:
#   LOG  (default: docs/phase4_9mark.txt, relative to cluster/)
#   PYTHON (default: cluster env interpreter)

cd "$(dirname "$0")/.."   # cluster/

PYTHON="${PYTHON:-/opt/homebrew/anaconda3/envs/cluster/bin/python3}"
SCRIPT=scripts/05_grouped_analyses.py
LOG="${LOG:-docs/phase4_9mark.txt}"
CONSENSUS="${1:-both}"
STATES="${2:-both}"

mkdir -p "$(dirname "$LOG")"

run_one() {
  local consensus="$1"
  local nstates="$2"
  echo ""
  echo "============================================================"
  echo "Phase 4.4+4.5: ${consensus} consensus, k=${nstates}"
  echo "============================================================"

  export CLUSTER_NSTATES="${nstates}"
  export CLUSTER_CHROMHMM_SUBDIR="chromHMM_9mark_${consensus}"
  export CLUSTER_MODEL_SUBDIR="learned_model_${nstates}"
  export CLUSTER_ENRICH_SUFFIX="_${nstates}"
  export CLUSTER_RENAME_SUFFIX="cerebellum"

  "$PYTHON" "$SCRIPT" --analyses 4.4,4.5
  local status=$?
  if [ ${status} -ne 0 ]; then
    echo "ERROR: Phase 4 failed for ${consensus}/k=${nstates} (exit ${status})" >&2
    return ${status}
  fi

  echo ""
  echo "--- Outputs for ${consensus}/k=${nstates} ---"
  local CHM="outputs/bap1_late/chromHMM_9mark_${consensus}"
  ls -lh "${CHM}/anchor_${nstates}.txt" "${CHM}/span_${nstates}.txt" 2>/dev/null || echo "  (enrichment tables missing)"
  ls -lh "${CHM}/anchor_${nstates}".{png,pdf} 2>/dev/null || true
  ls -lh "outputs/bap1_late/figures/chromHMM_anchor_${nstates}/" 2>/dev/null || echo "  (proportions missing)"
  return 0
}

{
  echo "============================================================"
  echo "Phase 4 (9-mark): anchor-vs-span + proportions"
  echo "Cluster dir: $(pwd)"
  echo "Started:     $(date)"
  echo "Consensus:   ${CONSENSUS}"
  echo "States:      ${STATES}"
  echo "Python:      ${PYTHON}"
  echo "chromhmm:    $(which chromhmm 2>/dev/null || echo NOT_ON_PATH)"
  echo "============================================================"

  # Build list of (consensus, nstates) pairs to run
  CONSENSUS_LIST=()
  if [ "${CONSENSUS}" = "both" ]; then
    CONSENSUS_LIST=(intersect union)
  else
    CONSENSUS_LIST=("${CONSENSUS}")
  fi

  STATES_LIST=()
  if [ "${STATES}" = "both" ]; then
    STATES_LIST=(15 18)
  else
    STATES_LIST=("${STATES}")
  fi

  FAIL=0
  for c in "${CONSENSUS_LIST[@]}"; do
    for s in "${STATES_LIST[@]}"; do
      run_one "${c}" "${s}"
      if [ $? -ne 0 ]; then FAIL=1; fi
    done
  done

  echo ""
  echo "============================================================"
  if [ ${FAIL} -ne 0 ]; then
    echo "Phase 4 (9-mark) finished WITH ERRORS: $(date)"
  else
    echo "Phase 4 (9-mark) finished: $(date)"
  fi
  echo "============================================================"
} 2>&1 | tee "$LOG"
