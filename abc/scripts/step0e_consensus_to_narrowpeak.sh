#!/bin/bash
# Convert consensus_all.bed to narrowPeak format for ABC pipeline injection
# narrowPeak: chr start end name score strand signalValue pValue qValue peak
# 'peak' = summit offset from start (0-based). For consensus: midpoint - start.

CONSENSUS="/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/peaks/atac_seq/consensus_all.bed"
OUTPUT="/expanse/lustre/projects/csd940/zalibhai/abc/input/consensus_all.narrowPeak"

if [ ! -f "${CONSENSUS}" ]; then
  echo "ERROR: Consensus file not found: ${CONSENSUS}"
  exit 1
fi

awk 'BEGIN{OFS="\t"; i=0} {
  i++
  mid = int(($2 + $3) / 2)
  summit_offset = mid - $2
  print $1, $2, $3, "consensus_peak_"i, 0, ".", 0, -1, -1, summit_offset
}' "${CONSENSUS}" \
  | sort -k1,1 -k2,2n \
  > "${OUTPUT}"

N=$(wc -l < "${OUTPUT}")
echo "Generated narrowPeak with ${N} peaks (expect 75371)"

if [ "${N}" -ne 75371 ]; then
  echo "WARNING: Peak count mismatch. Expected 75371, got ${N}."
fi

echo "First 3 lines:"
head -3 "${OUTPUT}"

echo "Column count: $(head -1 "${OUTPUT}" | awk '{print NF}') (expect 10)"

echo "Summit offset stats:"
awk '{print $10}' "${OUTPUT}" | sort -n | awk '
  NR==1{min=$1} {sum+=$1; vals[NR]=$1}
  END{printf "  min=%d median=%d max=%d mean=%.0f\n", min, vals[int(NR/2)], vals[NR], sum/NR}'

echo "Output: ${OUTPUT}"
