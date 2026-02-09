#!/usr/bin/env bash
# peaks/atac_seq/generate_consensus.sh
# Generate consensus ATAC-seq peaks from raw replicate narrowPeak files
#
# Approach:
#   1. Filter non-standard chromosomes (chrUn_*, chr*_random, chrM)
#   2. Extract BED3 (chr, start, end), sort
#   3. Generate per-genotype consensus via bedtools multiIntersectBed
#   4. Filter to regions present in >=2 replicates, merge overlapping
#   5. Report stats for threshold sweep (>=1, >=2, >=3)
#
# Dual mode: batch-1-only (n=3/genotype, primary) AND all-samples (n=5/genotype)
#
# Usage: bash peaks/atac_seq/generate_consensus.sh
#   Run from the repository root (mariner_hi-c/)

# =============================================================================
# CONFIGURATION
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RAW_DIR="${SCRIPT_DIR}/raw"
OUTPUT_DIR="${SCRIPT_DIR}"
TMP_DIR="${SCRIPT_DIR}/.tmp_consensus"
THRESHOLD=2  # Minimum replicate overlap for consensus

# Standard chromosomes to retain
STANDARD_CHROMS="^chr[0-9XY]*\t"

echo "================================================================================"
echo "ATAC-seq Consensus Peak Generation"
echo "================================================================================"
echo "Date: $(date)"
echo "Raw peak directory: ${RAW_DIR}"
echo "Output directory: ${OUTPUT_DIR}"
echo "Replicate threshold: >= ${THRESHOLD}"
echo ""

# =============================================================================
# VALIDATE DEPENDENCIES
# =============================================================================

for tool in bedtools sort cut awk grep; do
  if ! command -v "${tool}" &>/dev/null; then
    echo "ERROR: Required tool '${tool}' not found in PATH"
    exit 1
  fi
done

# =============================================================================
# SAMPLE ASSIGNMENT
# =============================================================================

# Batch 1 (250310): 3 control + 3 mutant replicates
BATCH1_CTRL=()
BATCH1_MUT=()

# Batch 2 (250731): 2 control + 2 mutant replicates
BATCH2_CTRL=()
BATCH2_MUT=()

for f in "${RAW_DIR}"/*.bed; do
  fname="$(basename "$f")"
  case "${fname}" in
    250310*Ctr[123]_P100*)  BATCH1_CTRL+=("$f") ;;
    250310*[mM]ut[123]_P100*) BATCH1_MUT+=("$f") ;;
    250731*Ctr[12]_3[0-9]*)  BATCH2_CTRL+=("$f") ;;
    250731*[mM]ut[12]_3[0-9]*) BATCH2_MUT+=("$f") ;;
    *) echo "  WARNING: Unmatched file: ${fname}" ;;
  esac
done

echo "Sample assignment:"
echo "  Batch 1 Control (n=${#BATCH1_CTRL[@]}): $(printf '%s\n' "${BATCH1_CTRL[@]}" | xargs -I{} basename {} | tr '\n' ', ')"
echo "  Batch 1 Mutant  (n=${#BATCH1_MUT[@]}):  $(printf '%s\n' "${BATCH1_MUT[@]}" | xargs -I{} basename {} | tr '\n' ', ')"
echo "  Batch 2 Control (n=${#BATCH2_CTRL[@]}): $(printf '%s\n' "${BATCH2_CTRL[@]}" | xargs -I{} basename {} | tr '\n' ', ')"
echo "  Batch 2 Mutant  (n=${#BATCH2_MUT[@]}):  $(printf '%s\n' "${BATCH2_MUT[@]}" | xargs -I{} basename {} | tr '\n' ', ')"
echo ""

# Validate minimum sample counts
if [ "${#BATCH1_CTRL[@]}" -lt 2 ] || [ "${#BATCH1_MUT[@]}" -lt 2 ]; then
  echo "ERROR: Batch 1 requires at least 2 replicates per genotype"
  exit 1
fi

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# Filter to standard chromosomes, extract BED3, sort
filter_and_sort() {
  local infile="$1"
  local outfile="$2"
  grep -E "${STANDARD_CHROMS}" "${infile}" \
    | cut -f1-3 \
    | sort -k1,1 -k2,2n \
    > "${outfile}"
}

# Generate consensus from a set of BED files at a given threshold
# Args: threshold, output_bed, bed_file_1 bed_file_2 ...
generate_consensus() {
  local thresh="$1"
  local outbed="$2"
  shift 2
  local beds=("$@")

  bedtools multiinter -i "${beds[@]}" \
    | awk -v t="${thresh}" '$4 >= t' \
    | cut -f1-3 \
    | bedtools merge -i - \
    > "${outbed}"
}

# Report peak counts at multiple thresholds
report_threshold_sweep() {
  local label="$1"
  shift
  local beds=("$@")
  local n_files=${#beds[@]}

  echo "  ${label} (n=${n_files} replicates):"

  local total_raw=0
  for b in "${beds[@]}"; do
    local n
    n=$(wc -l < "$b")
    total_raw=$((total_raw + n))
  done
  local avg_raw=$((total_raw / n_files))
  echo "    Average peaks per replicate: ${avg_raw}"

  local max_thresh=${n_files}
  if [ "${max_thresh}" -gt 3 ]; then
    max_thresh=3
  fi

  for t in $(seq 1 "${max_thresh}"); do
    local count
    count=$(bedtools multiinter -i "${beds[@]}" \
              | awk -v t="${t}" '$4 >= t' \
              | cut -f1-3 \
              | bedtools merge -i - \
              | wc -l)
    echo "    Threshold >= ${t}: ${count} consensus peaks"

    # Sanity check at chosen threshold
    if [ "${t}" -eq "${THRESHOLD}" ]; then
      local pct_retained
      pct_retained=$(awk "BEGIN { printf \"%.1f\", 100 * ${count} / ${avg_raw} }")
      # Warn if <30% of average peaks retained (i.e. >70% lost)
      if [ $((count * 10)) -lt $((avg_raw * 3)) ]; then
        echo "    WARNING: >70% of average peaks lost at threshold >= ${t} (${pct_retained}% retained)"
      fi
    fi
  done
  echo ""
}

# =============================================================================
# PROCESS FILES
# =============================================================================

mkdir -p "${TMP_DIR}"

echo "Filtering chromosomes and extracting BED3..."

# Process all raw files
ALL_CTRL_CLEAN=()
ALL_MUT_CLEAN=()
B1_CTRL_CLEAN=()
B1_MUT_CLEAN=()

# Batch 1 control
for f in "${BATCH1_CTRL[@]}"; do
  out="${TMP_DIR}/$(basename "$f" .bed).clean.bed"
  filter_and_sort "$f" "$out"
  B1_CTRL_CLEAN+=("$out")
  ALL_CTRL_CLEAN+=("$out")
done

# Batch 1 mutant
for f in "${BATCH1_MUT[@]}"; do
  out="${TMP_DIR}/$(basename "$f" .bed).clean.bed"
  filter_and_sort "$f" "$out"
  B1_MUT_CLEAN+=("$out")
  ALL_MUT_CLEAN+=("$out")
done

# Batch 2 control
for f in "${BATCH2_CTRL[@]}"; do
  out="${TMP_DIR}/$(basename "$f" .bed).clean.bed"
  filter_and_sort "$f" "$out"
  ALL_CTRL_CLEAN+=("$out")
done

# Batch 2 mutant
for f in "${BATCH2_MUT[@]}"; do
  out="${TMP_DIR}/$(basename "$f" .bed).clean.bed"
  filter_and_sort "$f" "$out"
  ALL_MUT_CLEAN+=("$out")
done

echo "  Processed ${#ALL_CTRL_CLEAN[@]} control + ${#ALL_MUT_CLEAN[@]} mutant files"
echo ""

# =============================================================================
# MODE 1: BATCH 1 ONLY (PRIMARY - n=3 per genotype)
# =============================================================================

echo "================================================================================"
echo "MODE 1: BATCH 1 ONLY (n=3 per genotype) - PRIMARY"
echo "================================================================================"
echo ""

echo "Threshold sweep:"
report_threshold_sweep "Control" "${B1_CTRL_CLEAN[@]}"
report_threshold_sweep "Mutant" "${B1_MUT_CLEAN[@]}"

echo "Generating batch-1 consensus at threshold >= ${THRESHOLD}..."

generate_consensus "${THRESHOLD}" "${OUTPUT_DIR}/consensus_control.bed" "${B1_CTRL_CLEAN[@]}"
generate_consensus "${THRESHOLD}" "${OUTPUT_DIR}/consensus_mutant.bed" "${B1_MUT_CLEAN[@]}"

# Union: merge control + mutant consensus
cat "${OUTPUT_DIR}/consensus_control.bed" "${OUTPUT_DIR}/consensus_mutant.bed" \
  | sort -k1,1 -k2,2n \
  | bedtools merge -i - \
  > "${OUTPUT_DIR}/consensus_all.bed"

N_CTRL=$(wc -l < "${OUTPUT_DIR}/consensus_control.bed")
N_MUT=$(wc -l < "${OUTPUT_DIR}/consensus_mutant.bed")
N_ALL=$(wc -l < "${OUTPUT_DIR}/consensus_all.bed")

echo "  consensus_control.bed: ${N_CTRL} peaks"
echo "  consensus_mutant.bed:  ${N_MUT} peaks"
echo "  consensus_all.bed:     ${N_ALL} peaks (union)"
echo ""

# =============================================================================
# MODE 2: ALL SAMPLES (n=5 per genotype) - SECONDARY
# =============================================================================

echo "================================================================================"
echo "MODE 2: ALL SAMPLES (n=5 per genotype) - SECONDARY"
echo "  Note: Batch 2 identity uncertain (may be re-sequenced samples)"
echo "================================================================================"
echo ""

echo "Threshold sweep:"
report_threshold_sweep "Control" "${ALL_CTRL_CLEAN[@]}"
report_threshold_sweep "Mutant" "${ALL_MUT_CLEAN[@]}"

echo "Generating all-sample consensus at threshold >= ${THRESHOLD}..."

generate_consensus "${THRESHOLD}" "${OUTPUT_DIR}/consensus_control_all.bed" "${ALL_CTRL_CLEAN[@]}"
generate_consensus "${THRESHOLD}" "${OUTPUT_DIR}/consensus_mutant_all.bed" "${ALL_MUT_CLEAN[@]}"

cat "${OUTPUT_DIR}/consensus_control_all.bed" "${OUTPUT_DIR}/consensus_mutant_all.bed" \
  | sort -k1,1 -k2,2n \
  | bedtools merge -i - \
  > "${OUTPUT_DIR}/consensus_all_combined.bed"

N_CTRL_ALL=$(wc -l < "${OUTPUT_DIR}/consensus_control_all.bed")
N_MUT_ALL=$(wc -l < "${OUTPUT_DIR}/consensus_mutant_all.bed")
N_ALL_ALL=$(wc -l < "${OUTPUT_DIR}/consensus_all_combined.bed")

echo "  consensus_control_all.bed: ${N_CTRL_ALL} peaks"
echo "  consensus_mutant_all.bed:  ${N_MUT_ALL} peaks"
echo "  consensus_all_combined.bed: ${N_ALL_ALL} peaks (union)"
echo ""

# =============================================================================
# COMPARISON SUMMARY
# =============================================================================

echo "================================================================================"
echo "SUMMARY"
echo "================================================================================"
echo ""
echo "PRIMARY (batch-1 only, n=3):"
echo "  Control consensus: ${N_CTRL} | Mutant consensus: ${N_MUT} | Union: ${N_ALL}"
echo ""
echo "SECONDARY (all samples, n=5):"
echo "  Control consensus: ${N_CTRL_ALL} | Mutant consensus: ${N_MUT_ALL} | Union: ${N_ALL_ALL}"
echo ""
echo "Differential peaks (pre-existing):"
echo "  ATAC_up.bed:   $(wc -l < "${OUTPUT_DIR}/ATAC_up.bed") peaks"
echo "  ATAC_down.bed: $(wc -l < "${OUTPUT_DIR}/ATAC_down.bed") peaks"
echo ""

# =============================================================================
# CLEANUP
# =============================================================================

rm -rf "${TMP_DIR}"

echo "Temporary files cleaned up."
echo "Consensus generation complete: $(date)"
echo "================================================================================"
