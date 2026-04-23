#!/usr/bin/env bash
# stripes/stripenn/scripts/run_stripe_integration.sh
#
# Runs all four Tier-1 stripe integration scripts across both timepoints
# and both resolutions. Designed to run from the stripes/stripenn/ directory.
#
# Usage:
#   cd stripes/stripenn
#   bash scripts/run_stripe_integration.sh [--n-permutations 10000]
#
# Default: 10000 permutations for T1.4 (loop crossref).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
N_PERM=10000

for arg in "$@"; do
  case "$arg" in
    --n-permutations=*) N_PERM="${arg#*=}" ;;
  esac
done

TIMEPOINTS=("late" "early")
RESOLUTIONS=("5000" "10000")

FAILED=0
TOTAL=0

run_script() {
  local label="$1"; shift
  TOTAL=$((TOTAL + 1))
  echo ""
  echo "================================================================"
  echo "  ${label}"
  echo "================================================================"
  if "$@"; then
    echo "  -> SUCCESS"
  else
    echo "  -> FAILED (exit $?)"
    FAILED=$((FAILED + 1))
  fi
}

echo "Stripe Integration Pipeline"
echo "Permutations: ${N_PERM}"
echo "Timepoints:   ${TIMEPOINTS[*]}"
echo "Resolutions:  ${RESOLUTIONS[*]}"
echo "Started:      $(date)"

for TP in "${TIMEPOINTS[@]}"; do
  for RES in "${RESOLUTIONS[@]}"; do
    TAG="${TP} res=${RES}"

    run_script "T1.1 Body-Gene Enrichment (${TAG})" \
      Rscript "${SCRIPT_DIR}/stripe_body_gene_enrichment.R" \
        --timepoint "$TP" --resolution "$RES"

    run_script "T1.2 Compartment Cross-Ref (${TAG})" \
      Rscript "${SCRIPT_DIR}/stripe_compartment_crossref.R" \
        --timepoint "$TP" --resolution "$RES"

    run_script "T1.3 DEG Stripe Violin (${TAG})" \
      Rscript "${SCRIPT_DIR}/deg_stripe_violin.R" \
        --timepoint "$TP" --resolution "$RES"

    run_script "T1.4 Loop x Stripe Cross-Ref (${TAG}, n_perm=${N_PERM})" \
      Rscript "${SCRIPT_DIR}/stripe_loop_crossref_extended.R" \
        --timepoint "$TP" --resolution "$RES" \
        --n-permutations "$N_PERM"
  done
done

echo ""
echo "================================================================"
echo "  COMPLETE: ${TOTAL} tasks, ${FAILED} failed"
echo "  Finished: $(date)"
echo "================================================================"

exit "${FAILED}"
