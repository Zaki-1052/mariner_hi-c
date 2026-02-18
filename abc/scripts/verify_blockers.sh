#!/bin/bash
# abc/scripts/verify_blockers.sh
# Run from cluster. Captures all info needed before executing the plan.
# Usage: bash verify_blockers.sh 2>&1 | tee blocker_check_output.txt

echo "========================================================================"
echo "ABC BLOCKER VERIFICATION — $(date)"
echo "========================================================================"

# ---------- 1. AWK VERSION (Step 0A GTF parsing) ----------
echo ""
echo "=== 1. AWK VERSION ==="
awk --version 2>&1 | head -1
# If this says "mawk" or errors, the gawk match() syntax won't work.
# Need the Perl fallback in Step 0A.
which gawk 2>/dev/null && echo "gawk found at: $(which gawk)" || echo "gawk NOT found in PATH"
module avail 2>&1 | grep -i gawk || echo "No gawk module found"

# ---------- 2. LAB BLACKLIST FORMAT (Step 0B) ----------
echo ""
echo "=== 2. LAB BLACKLIST FORMAT ==="
LAB_BL="/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/tads/250123blacklist.bed"
ENCODE_BL="/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/tads/mm10-blacklist.v2.bed"

echo "--- Lab blacklist ---"
ls -la "${LAB_BL}" 2>/dev/null || echo "FILE NOT FOUND"
echo "Columns: $(head -1 "${LAB_BL}" | awk '{print NF}')"
echo "Lines: $(wc -l < "${LAB_BL}")"
echo "First 5 lines:"
head -5 "${LAB_BL}"
echo "Has header? (check if first line starts with chr):"
head -1 "${LAB_BL}" | grep -c "^chr"

echo ""
echo "--- ENCODE blacklist ---"
ls -la "${ENCODE_BL}" 2>/dev/null || echo "FILE NOT FOUND"
echo "Columns: $(head -1 "${ENCODE_BL}" | awk '{print NF}')"
echo "Lines: $(wc -l < "${ENCODE_BL}")"
echo "First 5 lines:"
head -5 "${ENCODE_BL}"

# ---------- 3. TAGALIGN STATUS (Step 0C) ----------
echo ""
echo "=== 3. TAGALIGN COMPRESSION STATUS ==="
INPUT_DIR="/expanse/lustre/projects/csd940/zalibhai/abc/input"

for f in ctrl_ATAC.tagAlign.gz mut_ATAC.tagAlign.gz; do
  FULL="${INPUT_DIR}/${f}"
  echo "--- ${f} ---"
  ls -lh "${FULL}" 2>/dev/null || echo "FILE NOT FOUND"
  # Check if .tbi already exists
  ls -la "${FULL}.tbi" 2>/dev/null || echo "  No .tbi index"
  # Check if .bgz version exists
  ls -la "${INPUT_DIR}/$(echo ${f} | sed 's/.gz/.bgz/')" 2>/dev/null || echo "  No .bgz version"
  # Check if bgzip-compatible (magic bytes)
  echo "  First 3 lines of content:"
  zcat "${FULL}" 2>/dev/null | head -3
  echo "  Column count: $(zcat "${FULL}" 2>/dev/null | head -1 | awk '{print NF}')"
  echo "  Total reads: $(zcat "${FULL}" 2>/dev/null | wc -l)"
done

# ---------- 4. SNAKEFILE STRUCTURE (Step 2 — consensus peak injection) ----------
echo ""
echo "=== 4. SNAKEFILE STRUCTURE ==="
ABC_REPO="/expanse/lustre/projects/csd940/zalibhai/abc/ABC-Enhancer-Gene-Prediction"

echo "--- Snakefile / rules location ---"
ls -la "${ABC_REPO}/workflow/Snakefile" 2>/dev/null || ls -la "${ABC_REPO}/Snakefile" 2>/dev/null || echo "Snakefile NOT FOUND"
ls "${ABC_REPO}/workflow/rules/" 2>/dev/null || echo "No rules/ directory"

echo ""
echo "--- MACS2 rule (how peaks are called) ---"
grep -rn "macs2\|MACS2\|call_peaks\|narrowPeak" "${ABC_REPO}/workflow/" 2>/dev/null | head -30

echo ""
echo "--- makeCandidateRegions references ---"
grep -rn "makeCandidateRegions\|candidate_regions\|CandidateRegions" "${ABC_REPO}/workflow/" 2>/dev/null | head -20

echo ""
echo "--- Expected MACS2 output path pattern ---"
grep -rn "narrowPeak" "${ABC_REPO}/workflow/" 2>/dev/null | grep -i "output\|rule" | head -10

# ---------- 5. EXISTING CONFIG.YAML (Step 1) ----------
echo ""
echo "=== 5. EXISTING CONFIG.YAML ==="
CONFIG="${ABC_REPO}/config/config.yaml"
if [ -f "${CONFIG}" ]; then
  echo "Full contents:"
  cat "${CONFIG}"
else
  echo "config.yaml NOT FOUND at ${CONFIG}"
fi

echo ""
echo "--- Existing biosample configs ---"
ls -la "${ABC_REPO}/config/"*.tsv 2>/dev/null

# ---------- 6. CONSENSUS PEAKS STATUS (Step 0E) ----------
echo ""
echo "=== 6. CONSENSUS PEAKS ==="
CONSENSUS="/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/peaks/atac_seq/consensus_all.bed"
echo "--- consensus_all.bed ---"
ls -la "${CONSENSUS}" 2>/dev/null || echo "FILE NOT FOUND"
echo "Lines: $(wc -l < "${CONSENSUS}" 2>/dev/null)"
echo "Columns: $(head -1 "${CONSENSUS}" | awk '{print NF}')"
echo "First 5 lines:"
head -5 "${CONSENSUS}"
echo "Chr prefix check: $(head -1 "${CONSENSUS}" | cut -f1)"

CONSENSUS_500="/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/peaks/atac_seq/consensus_all_500bp.bed"
echo ""
echo "--- consensus_all_500bp.bed ---"
ls -la "${CONSENSUS_500}" 2>/dev/null || echo "FILE NOT FOUND"
echo "Lines: $(wc -l < "${CONSENSUS_500}" 2>/dev/null)"

# ---------- 7. HI-C FILES (Step 4) ----------
echo ""
echo "=== 7. HI-C FILES ==="
for f in resorted_ctrl.hic resorted_mut.hic; do
  FULL="/expanse/lustre/projects/csd940/zalibhai/nf-hic/250402_Bap1_deepseq/juicerpre/merged/hic/${f}"
  echo "--- ${f} ---"
  ls -lh "${FULL}" 2>/dev/null || echo "FILE NOT FOUND"
done

# ---------- 8. RNA-SEQ FILE (Steps 6-7) ----------
echo ""
echo "=== 8. RNA-SEQ FILE ==="
RNASEQ="/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/tads/adult_timepoint_rna-seq-BAP1_WT_KO_v2_Results.xlsx"
ls -lh "${RNASEQ}" 2>/dev/null || echo "FILE NOT FOUND"

# ---------- 9. CONDA ENVIRONMENT ----------
echo ""
echo "=== 9. CONDA ENVIRONMENT ==="
conda env list 2>/dev/null | grep abc
conda activate abc-env 2>/dev/null && echo "abc-env activated OK" && conda list 2>/dev/null | grep -E "snakemake|bedtools|samtools|macs2|tabix|htslib|bgzip" || echo "abc-env activation failed or not found"

# ---------- 10. CHROM SIZES (sanity check) ----------
echo ""
echo "=== 10. CHROM SIZES ==="
CHROMSIZES="/expanse/lustre/projects/csd940/zalibhai/taiji-new/reference/mm10/mm10.chrom.sizes"
ls -la "${CHROMSIZES}" 2>/dev/null || echo "FILE NOT FOUND"
echo "First 5 lines:"
head -5 "${CHROMSIZES}"
echo "Has chr prefix: $(head -1 "${CHROMSIZES}" | grep -c "^chr")"
echo "Total chromosomes: $(wc -l < "${CHROMSIZES}")"

echo ""
echo "========================================================================"
echo "VERIFICATION COMPLETE — $(date)"
echo "========================================================================"
echo ""
echo "KEY DECISIONS NEEDED FROM OUTPUT:"
echo "  1. Is awk GNU awk? If not, Step 0A needs Perl fallback"
echo "  2. Lab blacklist format — does it need header stripping?"
echo "  3. TagAligns — need bgzip conversion? (almost certainly yes)"
echo "  4. Snakefile — what's the MACS2 output path for symlink injection?"
echo "  5. Config.yaml — what params exist and need changing?"
