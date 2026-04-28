#!/usr/bin/env bash
# cluster/scripts/03_chromhmm_segmentation.sh
# Phase 2 worker: ChromHMM 12-state segmentation from 5 cerebellum ChIP peak BEDs.
# Stages 3-col peak BEDs (filtered to chr1-19 + chrX + chrY) into the chromHMM
# subtree, builds the cellmarkfiletable, runs BinarizeBed, then LearnModel(k=12).
# Per project convention: NO `set -euo pipefail`; explicit $? checks instead.
#
# Inputs (verified by Phase 2 plan):
#   peaks/beds/H3K27acCerebellumLate2.bed      (15,105 peaks)
#   peaks/beds/H3K27me3CerebellumLate1.bed     (15,809 peaks)
#   peaks/beds/H3K4me1CerebellumLate1.bed      (113,781 peaks)
#   peaks/beds/H3K4me3CerebellumLate2.bed      (6,581 peaks)
#   peaks/CTCF.bed                             (32,487 peaks)
#
# Outputs (under cluster/bap1_late/chromHMM/):
#   peak_beds/{H3K27ac,H3K27me3,H3K4me1,H3K4me3,CTCF}.bed   filtered 3-col BEDs
#   mm10_standard.txt                                       21-row chromsizes
#   cellmarkfiletable.txt                                   5-row table
#   binarized/cerebellum_late_chr*_binary.txt               21 binary files
#   learned_model/cerebellum_late_12_segments.bed           segmentation
#   learned_model/emissions_12.txt                          12x5 prob matrix
#   learned_model/{transitions_12.txt,model_12.txt,webpage_12.html,...}
#
# Runtime: ~1-2 min for stages 1-4, ~30-90 min for LearnModel.

REPO_ROOT="${REPO_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
CLUSTER_OUT_DIR="${CLUSTER_OUT_DIR:-outputs/bap1_late}"
CLUSTER_CELL_NAME="${CLUSTER_CELL_NAME:-cerebellum_late}"
CHROMHMM="${REPO_ROOT}/cluster/ChromHMM/chromhmm"
CHM_DIR="${REPO_ROOT}/cluster/${CLUSTER_OUT_DIR}/chromHMM"
PEAK_DIR="${CHM_DIR}/peak_beds"
BIN_DIR="${CHM_DIR}/binarized"
MODEL_DIR="${CHM_DIR}/learned_model"
CMT="${CHM_DIR}/cellmarkfiletable.txt"
CHROMSIZES_SRC="${REPO_ROOT}/cluster/ChromHMM/CHROMSIZES/mm10.txt"
CHROMSIZES_DST="${CHM_DIR}/mm10_standard.txt"
CELL="${CLUSTER_CELL_NAME}"
NSTATES=12
ASSEMBLY="mm10"

cd "${REPO_ROOT}"

mkdir -p "${PEAK_DIR}" "${BIN_DIR}" "${MODEL_DIR}"

echo "============================================================"
echo "Phase 2: ChromHMM ${NSTATES}-state segmentation"
echo "Repo root: ${REPO_ROOT}"
echo "Cell:      ${CELL}"
echo "Started:   $(date)"
echo "============================================================"

# ---- [1/5] Stage peak BEDs ----------------------------------------------
echo ""
echo "[1/5] Staging peak BEDs (filter to chr1-19 + chrX + chrY, take cols 1-3)..."
# Wipe any stale staging from a previous run so chromsome filter is authoritative.
rm -f "${PEAK_DIR}"/*.bed

MARKS=(H3K27ac H3K27me3 H3K4me1 H3K4me3 CTCF)
SOURCES=(
  "${CLUSTER_PEAK_K27AC:-peaks/beds/H3K27acCerebellumLate2.bed}"
  "${CLUSTER_PEAK_K27ME3:-peaks/beds/H3K27me3CerebellumLate1.bed}"
  "${CLUSTER_PEAK_K4ME1:-peaks/beds/H3K4me1CerebellumLate1.bed}"
  "${CLUSTER_PEAK_K4ME3:-peaks/beds/H3K4me3CerebellumLate2.bed}"
  "${CLUSTER_PEAK_CTCF:-peaks/CTCF.bed}"
)
for i in "${!MARKS[@]}"; do
  MARK="${MARKS[$i]}"
  SRC="${SOURCES[$i]}"
  DST="${PEAK_DIR}/${MARK}.bed"
  if [ ! -s "${SRC}" ]; then
    echo "  ERROR: missing source BED ${SRC}" >&2
    exit 1
  fi
  IN_LINES=$(wc -l < "${SRC}")
  awk 'BEGIN{OFS="\t"} $1 ~ /^(chr[0-9]+|chrX|chrY)$/ {print $1,$2,$3}' "${SRC}" \
    | LC_ALL=C sort -k1,1 -k2,2n -u \
    > "${DST}"
  OUT_LINES=$(wc -l < "${DST}")
  printf "    %-9s  %s -> %s   (%7d -> %7d)\n" \
    "${MARK}" "${SRC}" "${DST}" "${IN_LINES}" "${OUT_LINES}"
done

# ---- [2/5] Filter chromsizes -------------------------------------------
echo ""
echo "[2/5] Filtering chromsizes to chr1-19 + chrX + chrY..."
awk '$1 ~ /^(chr[0-9]+|chrX|chrY)$/' "${CHROMSIZES_SRC}" > "${CHROMSIZES_DST}"
N_CHROMS=$(wc -l < "${CHROMSIZES_DST}")
SUM_BP=$(awk '{s+=$2} END{print s}' "${CHROMSIZES_DST}")
echo "    ${CHROMSIZES_DST}: ${N_CHROMS} chroms, total ${SUM_BP} bp"
if [ "${N_CHROMS}" -ne 21 ]; then
  echo "    ERROR: expected 21 standard chroms, got ${N_CHROMS}" >&2
  exit 1
fi

# ---- [3/5] Build cellmarkfiletable -------------------------------------
echo ""
echo "[3/5] Writing cellmarkfiletable..."
{
  printf '%s\tH3K27ac\tH3K27ac.bed\n'   "${CELL}"
  printf '%s\tH3K27me3\tH3K27me3.bed\n' "${CELL}"
  printf '%s\tH3K4me1\tH3K4me1.bed\n'   "${CELL}"
  printf '%s\tH3K4me3\tH3K4me3.bed\n'   "${CELL}"
  printf '%s\tCTCF\tCTCF.bed\n'         "${CELL}"
} > "${CMT}"
echo "    ${CMT}:"
sed 's/^/        /' "${CMT}"

# ---- [4/5] BinarizeBed --------------------------------------------------
echo ""
echo "[4/5] Running ChromHMM BinarizeBed..."
# Wipe stale binary files so chrom set is authoritative.
rm -f "${BIN_DIR}"/*.txt
"${CHROMHMM}" BinarizeBed -peaks \
    "${CHROMSIZES_DST}" \
    "${PEAK_DIR}" \
    "${CMT}" \
    "${BIN_DIR}"
status=$?
if [ ${status} -ne 0 ]; then
  echo "ERROR: BinarizeBed exited with status ${status}" >&2
  exit ${status}
fi
N_BIN=$(ls "${BIN_DIR}" | grep -c "_binary\.txt$")
echo "    binary files produced: ${N_BIN} (expect 21)"
if [ "${N_BIN}" -ne 21 ]; then
  echo "    WARNING: unexpected binary file count" >&2
fi

# ---- [5/5] LearnModel ---------------------------------------------------
echo ""
echo "[5/5] Running ChromHMM LearnModel (k=${NSTATES}, ${ASSEMBLY})..."
echo "    Started:  $(date)  -- expect ~30-90 min on Mac"
"${CHROMHMM}" LearnModel \
    -p 4 \
    -l "${CHROMSIZES_DST}" \
    "${BIN_DIR}" \
    "${MODEL_DIR}" \
    "${NSTATES}" \
    "${ASSEMBLY}"
status=$?
if [ ${status} -ne 0 ]; then
  echo "ERROR: LearnModel exited with status ${status}" >&2
  exit ${status}
fi
echo "    Finished: $(date)"

# ---- Verification -------------------------------------------------------
SEG_BED="${MODEL_DIR}/${CELL}_${NSTATES}_segments.bed"
EMI="${MODEL_DIR}/emissions_${NSTATES}.txt"

echo ""
echo "============================================================"
echo "Phase 2 verification"
echo "============================================================"
if [ -s "${SEG_BED}" ]; then
  echo "  OK  ${SEG_BED}: $(wc -l < "${SEG_BED}") segments"
else
  echo "  MISSING  ${SEG_BED}" >&2
  exit 1
fi
if [ -s "${EMI}" ]; then
  echo "  OK  ${EMI}"
  echo ""
  echo "  --- emissions_${NSTATES}.txt ---"
  cat "${EMI}"
  echo "  --- end emissions ---"
else
  echo "  MISSING  ${EMI}" >&2
  exit 1
fi

echo ""
echo "  Chrom coverage in segments BED (expect 21 standard chroms):"
cut -f1 "${SEG_BED}" | sort -u | tr '\n' ' '
echo ""

echo ""
echo "  State distribution (col 4 of segments.bed, sorted by count):"
awk '{print $4}' "${SEG_BED}" | sort | uniq -c | sort -rn | sed 's/^/      /'

TOTAL_BINS=$(awk '{s += int(($3-$2)/200)} END{print s}' "${SEG_BED}")
echo ""
echo "  Total 200bp bins covered: ${TOTAL_BINS}"

echo ""
echo "============================================================"
echo "NEXT (manual):"
echo "  Inspect ${EMI} above. Each row = one state's emission probabilities"
echo "  for {H3K27ac, H3K27me3, H3K4me1, H3K4me3, CTCF}."
echo "  Write rename file at:"
echo "    ${CHM_DIR}/12state_rename_cerebellum.txt"
echo "  Format: <state_id><TAB><biological_name>  (one per line, 12 lines)"
echo "  Use underscored names matching project taxonomy:"
echo "    Active_Promoter, Bivalent_Promoter, Weak_Promoter,"
echo "    Active_Enhancer, Strong_Enhancer, Poised_Enhancer,"
echo "    Polycomb, Insulator, Active_CTCF, Quiescent[_N]"
echo "============================================================"
echo ""
echo "Phase 2 finished: $(date)"
