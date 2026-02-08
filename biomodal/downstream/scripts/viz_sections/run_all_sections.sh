#\!/usr/bin/env bash
# biomodal/downstream/scripts/viz_sections/run_all_sections.sh
# Run all visualization section scripts sequentially
# Usage: cd downstream/ && bash scripts/viz_sections/run_all_sections.sh

echo "================================================================================"
echo "Running all biomodal visualization sections"
echo "================================================================================"
echo "Working directory: $(pwd)"
echo "Start time: $(date)"
echo ""

for section in scripts/viz_sections/section_*.R; do
  echo "--------------------------------------------------------------------------------"
  echo "Running ${section}..."
  echo "--------------------------------------------------------------------------------"
  Rscript "${section}"
  echo ""
done

echo "================================================================================"
echo "All sections complete: $(date)"
echo "================================================================================"
