#!/usr/bin/env bash
# cluster/scripts/run_go_kegg.sh
# Driver: Phase 9 -- per-cluster GO/KEGG pathway enrichment.
# Run from cluster/ (cd is automatic via $0).
#
# Env vars:
#   LOG  (default docs/phase9_go_kegg.txt, relative to cluster/)
#
# Usage:
#   bash cluster/scripts/run_go_kegg.sh
#   LOG=docs/go_kegg_test.txt bash cluster/scripts/run_go_kegg.sh

set -e

cd "$(dirname "$0")/.."   # cluster/

SYS_RSCRIPT="/usr/local/bin/Rscript"
SCRIPT="scripts/12_go_kegg_enrichment.R"
LOG=${LOG:-docs/phase9_go_kegg.txt}

mkdir -p "$(dirname "$LOG")"

{
  echo "============================================================"
  echo "Phase 9: Per-cluster GO/KEGG enrichment"
  echo "Cluster dir: $(pwd)"
  echo "Started:     $(date)"
  echo "Log:         ${LOG}"
  echo "Rscript:     ${SYS_RSCRIPT}"
  echo "============================================================"
  echo

  echo "[1/1] Running 12_go_kegg_enrichment.R..."
  "${SYS_RSCRIPT}" "${SCRIPT}"
  r_status=$?
  if [ ${r_status} -ne 0 ]; then
    echo "ERROR: 12_go_kegg_enrichment.R exited with status ${r_status}" >&2
    exit ${r_status}
  fi

  echo
  echo "============================================================"
  echo "Phase 9 outputs (outputs/bap1_late/figures/go_enrichment/)"
  echo "============================================================"
  echo "TSV results:"
  ls -lh outputs/bap1_late/figures/go_enrichment/*.tsv 2>/dev/null || echo "  (no TSV outputs)"
  echo
  echo "Figure subfolders:"
  for d in outputs/bap1_late/figures/go_enrichment/*/; do
    if [ -d "$d" ]; then
      echo "  $(basename "$d")/"
      ls -1 "$d" 2>/dev/null | sed 's/^/    /'
    fi
  done

  echo
  echo "Phase 9 finished: $(date)"
} 2>&1 | tee "$LOG"
