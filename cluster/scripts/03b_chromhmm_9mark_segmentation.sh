#!/usr/bin/env bash
# cluster/scripts/03b_chromhmm_9mark_segmentation.sh
# Phase 2b worker: ChromHMM 15+18-state segmentation from 9 cerebellum marks.
# Mirrors 03_chromhmm_segmentation.sh structure. Same 5-stage pattern, explicit
# $? checks, no set -euo pipefail (benign Java/awk warnings would abort).
#
# 9 marks: H3K27ac, H3K27me3, H3K4me1, H3K4me3, CTCF (original 5)
#        + H2AK119ub, ATAC, H3K9ac, H3K9me3 (new 4)
#
# H3K9ac and H3K9me3 sources MUST be set via CLUSTER_PEAK_K9AC / CLUSTER_PEAK_K9ME3
# (no default — prevents accidentally using the wrong consensus).
#
# CLUSTER_CHROMHMM_SUBDIR controls the output subdirectory name (e.g.
# chromHMM_9mark_intersect or chromHMM_9mark_union).
#
# BinarizeBed runs once; LearnModel runs twice (k=15 and k=18).
#
# Timepoint caveat: H3K9me3 is early-only (no late peaks exist).
# H2AK119ub is P51 (~9 days from P60 late), treated as late.

REPO_ROOT="${REPO_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
CLUSTER_OUT_DIR="${CLUSTER_OUT_DIR:-outputs/bap1_late}"
CLUSTER_CELL_NAME="${CLUSTER_CELL_NAME:-cerebellum_late}"
CLUSTER_CHROMHMM_SUBDIR="${CLUSTER_CHROMHMM_SUBDIR:-chromHMM_9mark}"
CHROMHMM="${REPO_ROOT}/cluster/ChromHMM/chromhmm"
CHM_DIR="${REPO_ROOT}/cluster/${CLUSTER_OUT_DIR}/${CLUSTER_CHROMHMM_SUBDIR}"
PEAK_DIR="${CHM_DIR}/peak_beds"
BIN_DIR="${CHM_DIR}/binarized"
CMT="${CHM_DIR}/cellmarkfiletable.txt"
CHROMSIZES_SRC="${REPO_ROOT}/cluster/ChromHMM/CHROMSIZES/mm10.txt"
CHROMSIZES_DST="${CHM_DIR}/mm10_standard.txt"
CELL="${CLUSTER_CELL_NAME}"
NSTATES_LIST=(15 18)
ASSEMBLY="mm10"
NPROC="${SLURM_CPUS_PER_TASK:-4}"

if [ -z "${CLUSTER_PEAK_K9AC}" ] || [ -z "${CLUSTER_PEAK_K9ME3}" ]; then
  echo "ERROR: CLUSTER_PEAK_K9AC and CLUSTER_PEAK_K9ME3 must be set." >&2
  echo "  These point to the H3K9ac and H3K9me3 consensus BEDs (intersect or union)." >&2
  exit 1
fi

cd "${REPO_ROOT}"

mkdir -p "${PEAK_DIR}" "${BIN_DIR}"
for NS in "${NSTATES_LIST[@]}"; do
  mkdir -p "${CHM_DIR}/learned_model_${NS}"
done

echo "============================================================"
echo "Phase 2b: ChromHMM 9-mark segmentation (k=${NSTATES_LIST[*]})"
echo "Repo root:    ${REPO_ROOT}"
echo "Cell:         ${CELL}"
echo "Output dir:   ${CHM_DIR}"
echo "H3K9ac src:   ${CLUSTER_PEAK_K9AC}"
echo "H3K9me3 src:  ${CLUSTER_PEAK_K9ME3}"
echo "Parallelism:  ${NPROC} threads"
echo "Started:      $(date)"
echo "============================================================"

# ---- [1/5] Stage peak BEDs ------------------------------------------------
echo ""
echo "[1/5] Staging 9 peak BEDs (filter to chr1-19 + chrX + chrY, take cols 1-3)..."
rm -f "${PEAK_DIR}"/*.bed

MARKS=(H3K27ac H3K27me3 H3K4me1 H3K4me3 CTCF H2AK119ub ATAC H3K9ac H3K9me3)
SOURCES=(
  "${CLUSTER_PEAK_K27AC:-peaks/beds/H3K27acCerebellumLate2.bed}"
  "${CLUSTER_PEAK_K27ME3:-peaks/beds/H3K27me3CerebellumLate1.bed}"
  "${CLUSTER_PEAK_K4ME1:-peaks/beds/H3K4me1CerebellumLate1.bed}"
  "${CLUSTER_PEAK_K4ME3:-peaks/beds/H3K4me3CerebellumLate2.bed}"
  "${CLUSTER_PEAK_CTCF:-peaks/CTCF.bed}"
  "${CLUSTER_PEAK_K119UB:-peaks/intersect/P51_K119ub_ctrl_intersect.bed}"
  "${CLUSTER_PEAK_ATAC:-peaks/atac_seq/consensus_control.bed}"
  "${CLUSTER_PEAK_K9AC}"
  "${CLUSTER_PEAK_K9ME3}"
)
for i in "${!MARKS[@]}"; do
  MARK="${MARKS[$i]}"
  SRC="${SOURCES[$i]}"
  DST="${PEAK_DIR}/${MARK}.bed"
  if [ ! -s "${SRC}" ]; then
    echo "  ERROR: missing or empty source BED: ${SRC}" >&2
    exit 1
  fi
  IN_LINES=$(wc -l < "${SRC}")
  awk 'BEGIN{OFS="\t"} $1 ~ /^(chr[0-9]+|chrX|chrY)$/ {print $1,$2,$3}' "${SRC}" \
    | LC_ALL=C sort -k1,1 -k2,2n -u \
    > "${DST}"
  OUT_LINES=$(wc -l < "${DST}")
  printf "    %-12s  %-60s  (%7d -> %7d)\n" "${MARK}" "${SRC}" "${IN_LINES}" "${OUT_LINES}"
done

# ---- [2/5] Filter chromsizes ----------------------------------------------
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

# ---- [3/5] Build cellmarkfiletable ----------------------------------------
echo ""
echo "[3/5] Writing 9-mark cellmarkfiletable..."
{
  printf '%s\tH3K27ac\tH3K27ac.bed\n'   "${CELL}"
  printf '%s\tH3K27me3\tH3K27me3.bed\n' "${CELL}"
  printf '%s\tH3K4me1\tH3K4me1.bed\n'   "${CELL}"
  printf '%s\tH3K4me3\tH3K4me3.bed\n'   "${CELL}"
  printf '%s\tCTCF\tCTCF.bed\n'         "${CELL}"
  printf '%s\tH2AK119ub\tH2AK119ub.bed\n' "${CELL}"
  printf '%s\tATAC\tATAC.bed\n'          "${CELL}"
  printf '%s\tH3K9ac\tH3K9ac.bed\n'     "${CELL}"
  printf '%s\tH3K9me3\tH3K9me3.bed\n'   "${CELL}"
} > "${CMT}"
echo "    ${CMT}:"
sed 's/^/        /' "${CMT}"

# ---- [4/5] BinarizeBed ----------------------------------------------------
echo ""
echo "[4/5] Running ChromHMM BinarizeBed (9 marks)..."
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
HEADER_COLS=$(head -1 "${BIN_DIR}/${CELL}_chr1_binary.txt" | awk '{print NF}')
echo "    header columns in chr1 binary: ${HEADER_COLS} (expect 10 = cell + 9 marks)"

# ---- [5/5] LearnModel (loop over state counts) ----------------------------
for NSTATES in "${NSTATES_LIST[@]}"; do
  MODEL_DIR="${CHM_DIR}/learned_model_${NSTATES}"
  echo ""
  echo "[5/5] Running ChromHMM LearnModel (k=${NSTATES}, ${ASSEMBLY}, -p ${NPROC})..."
  echo "    Model dir: ${MODEL_DIR}"
  echo "    Started:   $(date)"
  "${CHROMHMM}" LearnModel \
      -p "${NPROC}" \
      -l "${CHROMSIZES_DST}" \
      "${BIN_DIR}" \
      "${MODEL_DIR}" \
      "${NSTATES}" \
      "${ASSEMBLY}"
  status=$?
  if [ ${status} -ne 0 ]; then
    echo "ERROR: LearnModel k=${NSTATES} exited with status ${status}" >&2
    exit ${status}
  fi
  echo "    Finished:  $(date)"

  SEG_BED="${MODEL_DIR}/${CELL}_${NSTATES}_segments.bed"
  EMI="${MODEL_DIR}/emissions_${NSTATES}.txt"

  echo ""
  echo "  --- Verification for k=${NSTATES} ---"
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
  echo "  Chrom coverage (expect 21):"
  cut -f1 "${SEG_BED}" | sort -u | tr '\n' ' '
  echo ""

  echo "  State distribution:"
  awk '{print $4}' "${SEG_BED}" | sort | uniq -c | sort -rn | sed 's/^/      /'
done

echo ""
echo "============================================================"
echo "NEXT (manual):"
for NSTATES in "${NSTATES_LIST[@]}"; do
  echo "  Inspect emissions_${NSTATES}.txt above (9 columns:"
  echo "    H3K27ac, H3K27me3, H3K4me1, H3K4me3, CTCF, H2AK119ub, ATAC, H3K9ac, H3K9me3)"
  echo "  Write rename file:"
  echo "    ${CHM_DIR}/${NSTATES}state_rename_cerebellum.txt"
  echo "  Format: E{N}<TAB>{Biological_Name}  (${NSTATES} lines)"
  echo ""
done
echo "============================================================"
echo ""
echo "Phase 2b finished: $(date)"
