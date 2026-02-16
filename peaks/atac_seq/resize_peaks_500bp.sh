#!/usr/bin/env bash
# peaks/atac_seq/resize_peaks_500bp.sh
# Resize consensus ATAC-seq peaks to 500bp centered on midpoint (standard for ABC)

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INPUT="${SCRIPT_DIR}/consensus_all.bed"
OUTPUT="${SCRIPT_DIR}/consensus_all_500bp.bed"
CHROMSIZES="${SCRIPT_DIR}/../../mm10.chrom.sizes"
WIDTH=500
HALF=$((WIDTH / 2))

if [ ! -f "${INPUT}" ]; then
  echo "ERROR: Input not found: ${INPUT}"
  exit 1
fi

if [ ! -f "${CHROMSIZES}" ]; then
  echo "ERROR: Chrom sizes not found: ${CHROMSIZES}"
  echo "  Download with: fetchChromSizes mm10 > ${CHROMSIZES}"
  exit 1
fi

N_BEFORE=$(wc -l < "${INPUT}")

awk -v h="${HALF}" 'BEGIN{OFS="\t"} {
  mid = int(($2 + $3) / 2)
  start = mid - h
  if (start < 0) start = 0
  print $1, start, mid + h
}' "${INPUT}" \
  | sort -k1,1 -k2,2n \
  | bedtools slop -i - -g "${CHROMSIZES}" -b 0 \
  | bedtools merge -i - \
  > "${OUTPUT}"

N_AFTER=$(wc -l < "${OUTPUT}")

echo "Resized ${N_BEFORE} peaks to ${WIDTH}bp → ${N_AFTER} peaks after merge"
echo "Output: ${OUTPUT}"
