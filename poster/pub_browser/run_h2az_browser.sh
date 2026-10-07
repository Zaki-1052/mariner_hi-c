#!/bin/bash
# poster/pub_browser/run_h2az_browser.sh
#
# Batch H2AZ/H2AZac browser views at 5 developmental x genotype cross-regulation loci.
#
# Each figure shows 3 marks (early ctrl, late ctrl, late mut) with
# H2AZ on top and H2AZac on bottom per mark — 6 tracks total per locus.
#
# Usage:
#   bash poster/pub_browser/run_h2az_browser.sh 2>&1 | tee poster/pub_browser/output/h2az/run_h2az.log

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PUB_BROWSER="${SCRIPT_DIR}/pub_browser.R"

configs=(
  "${SCRIPT_DIR}/configs/h2az_chr13_83M.yaml"
  "${SCRIPT_DIR}/configs/h2az_chr4_82M.yaml"
  "${SCRIPT_DIR}/configs/h2az_chr4_98M.yaml"
  "${SCRIPT_DIR}/configs/h2az_chr1_18M.yaml"
  "${SCRIPT_DIR}/configs/h2az_chr1_133M.yaml"
)

total=${#configs[@]}
succeeded=0
failed=0

for cfg in "${configs[@]}"; do
  name="$(basename "$cfg" .yaml)"
  echo "========================================"
  echo "[$(( succeeded + failed + 1 ))/${total}] ${name}"
  echo "========================================"
  if Rscript "$PUB_BROWSER" --config "$cfg"; then
    succeeded=$((succeeded + 1))
  else
    echo "  FAILED: ${name}"
    failed=$((failed + 1))
  fi
  echo ""
done

echo "========================================"
echo "Batch complete: ${succeeded}/${total} succeeded, ${failed} failed"
echo "Output: ${SCRIPT_DIR}/output/h2az/"
echo "========================================"
