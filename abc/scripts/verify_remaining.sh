#!/bin/bash
# abc/scripts/verify_remaining.sh
# Everything that was cut off or deferred from the first run.
# Usage: bash verify_remaining.sh 2>&1 | tee blocker_remaining_output.txt

ABC_REPO="/expanse/lustre/projects/csd940/zalibhai/abc/ABC-Enhancer-Gene-Prediction"

echo "========================================================================"
echo "ABC REMAINING VERIFICATION — $(date)"
echo "========================================================================"

# ==========================================================================
# A. SNAKEMAKE RULE FILES (full contents — need to understand the EXACT
#    input/output wiring before we can safely inject consensus peaks)
# ==========================================================================

echo ""
echo "=== A1. FULL macs2.smk ==="
cat "${ABC_REPO}/workflow/rules/macs2.smk"

echo ""
echo "=== A2. FULL candidate_regions.smk ==="
cat "${ABC_REPO}/workflow/rules/candidate_regions.smk"

echo ""
echo "=== A3. FULL neighborhoods.smk ==="
cat "${ABC_REPO}/workflow/rules/neighborhoods.smk"

echo ""
echo "=== A4. FULL Snakefile (main) ==="
cat "${ABC_REPO}/workflow/Snakefile"

echo ""
echo "=== A5. All rule files ==="
ls -la "${ABC_REPO}/workflow/rules/"

# ==========================================================================
# B. QUANTILE NORMALIZATION CONCERN
#    config.yaml has use_qnorm: True with K562 human reference.
#    Need to understand what this file looks like and whether it's
#    appropriate for mouse ATAC data.
# ==========================================================================

echo ""
echo "=== B1. QNorm reference file ==="
QNORM="${ABC_REPO}/reference/EnhancersQNormRef.K562.txt"
ls -la "${QNORM}" 2>/dev/null || echo "FILE NOT FOUND"
if [ -f "${QNORM}" ]; then
  echo "Lines: $(wc -l < "${QNORM}")"
  echo "First 10 lines:"
  head -10 "${QNORM}"
  echo "Last 5 lines:"
  tail -5 "${QNORM}"
fi

echo ""
echo "=== B2. ABC thresholds reference ==="
THRESH="${ABC_REPO}/reference/abc_thresholds.tsv"
ls -la "${THRESH}" 2>/dev/null || echo "FILE NOT FOUND"
if [ -f "${THRESH}" ]; then
  echo "Full contents:"
  cat "${THRESH}"
fi

# ==========================================================================
# C. EXISTING REFERENCE FILES IN THE REPO
#    Need to know what ships with ABC so we know exactly what to replace
# ==========================================================================

echo ""
echo "=== C. Reference directory contents ==="
ls -la "${ABC_REPO}/reference/" 2>/dev/null
ls -la "${ABC_REPO}/reference/hg38/" 2>/dev/null
echo ""
echo "Bundled gene annotation format (first 5 lines):"
head -5 "${ABC_REPO}/reference/hg38/CollapsedGeneBounds.hg38.bed" 2>/dev/null || echo "NOT FOUND"
echo ""
echo "Bundled TSS format (first 5 lines):"
head -5 "${ABC_REPO}/reference/hg38/CollapsedGeneBounds.hg38.TSS500bp.bed" 2>/dev/null || echo "NOT FOUND"
echo ""
echo "Bundled ubiquitous genes format (first 10 lines):"
head -10 "${ABC_REPO}/reference/UbiquitouslyExpressedGenes.txt" 2>/dev/null || echo "NOT FOUND"
echo ""
echo "Bundled blacklist format (first 5 lines):"
head -5 "${ABC_REPO}/reference/hg38/GRCh38_unified_blacklist.bed" 2>/dev/null || echo "NOT FOUND"

# ==========================================================================
# D. ITEMS CUT OFF FROM FIRST RUN
# ==========================================================================

echo ""
echo "=== D1. TAGALIGN READ COUNTS ==="
# The wc -l on 764MB gzipped file timed out last time. Use a faster method.
INPUT_DIR="/expanse/lustre/projects/csd940/zalibhai/abc/input"
for COND in ctrl mut; do
  F="${INPUT_DIR}/${COND}_ATAC.tagAlign.gz"
  if [ -f "$F" ]; then
    # Use pigz for parallel decompression if available, otherwise zcat
    echo "--- ${COND}_ATAC.tagAlign.gz ---"
    echo "  File size: $(ls -lh "$F" | awk '{print $5}')"
    # Approximate line count from file size (faster than full decompression)
    # But let's also get exact count — just be patient
    echo "  Counting reads (may take 1-2 min)..."
    COUNT=$(zcat "$F" | wc -l)
    echo "  Total reads: ${COUNT}"
    # Verify sort order (first 20 lines should be chr1 or chrN in order)
    echo "  Sort order check (first 5 chromosomes):"
    zcat "$F" | cut -f1 | uniq | head -5
  fi
done

echo ""
echo "=== D2. CONSENSUS PEAKS ==="
CONSENSUS="/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/peaks/atac_seq/consensus_all.bed"
CONSENSUS_500="/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/peaks/atac_seq/consensus_all_500bp.bed"

echo "--- consensus_all.bed ---"
ls -la "${CONSENSUS}" 2>/dev/null || echo "NOT FOUND"
echo "Lines: $(wc -l < "${CONSENSUS}" 2>/dev/null)"
echo "Columns: $(head -1 "${CONSENSUS}" | awk '{print NF}')"
echo "First 5:"
head -5 "${CONSENSUS}"
echo "Width stats (min/median/max):"
awk '{print $3-$2}' "${CONSENSUS}" | sort -n | awk '
  NR==1{min=$1} {vals[NR]=$1}
  END{printf "  min=%d median=%d max=%d n=%d\n", min, vals[int(NR/2)], vals[NR], NR}'

echo ""
echo "--- consensus_all_500bp.bed ---"
ls -la "${CONSENSUS_500}" 2>/dev/null || echo "NOT FOUND"
echo "Lines: $(wc -l < "${CONSENSUS_500}" 2>/dev/null)"
echo "Width stats:"
awk '{print $3-$2}' "${CONSENSUS_500}" | sort -n | awk '
  NR==1{min=$1} {vals[NR]=$1}
  END{printf "  min=%d median=%d max=%d n=%d\n", min, vals[int(NR/2)], vals[NR], NR}'

# Check chr prefix
echo "Chr prefix: $(head -1 "${CONSENSUS}" | cut -f1)"

# Check for non-standard chromosomes
echo "Chromosomes in consensus_all.bed:"
cut -f1 "${CONSENSUS}" | sort -u

echo ""
echo "=== D3. HI-C FILES ==="
for f in resorted_ctrl.hic resorted_mut.hic; do
  FULL="/expanse/lustre/projects/csd940/ctea/nf-hic/250402_Bap1_deepseq/juicerpre/merged/hic/${f}"
  echo "--- ${f} ---"
  ls -lh "${FULL}" 2>/dev/null || echo "FILE NOT FOUND"
done

echo ""
echo "=== D4. RNA-SEQ FILE ==="
RNASEQ="/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/tads/adult_timepoint_rna-seq-BAP1_WT_KO_v2_Results.xlsx"
ls -lh "${RNASEQ}" 2>/dev/null || echo "FILE NOT FOUND"
# Verify we can read it with python
python3 -c "
import pandas as pd
df = pd.read_excel('${RNASEQ}')
print(f'  Rows: {len(df)}')
print(f'  Columns: {list(df.columns)}')
print(f'  First gene: {df.iloc[0, 0]}')
print(f'  Gene column dtype: {df.iloc[:, 0].dtype}')
print(f'  Has baseMean: {\"baseMean\" in df.columns}')
print(f'  Has padj: {\"padj\" in df.columns}')
print(f'  Has log2FoldChange: {\"log2FoldChange\" in df.columns}')
" 2>&1

echo ""
echo "=== D5. CHROM SIZES ==="
CHROMSIZES="/expanse/lustre/projects/csd940/zalibhai/taiji-new/reference/mm10/mm10.chrom.sizes"
ls -la "${CHROMSIZES}" 2>/dev/null || echo "FILE NOT FOUND"
echo "First 5 lines:"
head -5 "${CHROMSIZES}"
echo "Has chr prefix: $(head -1 "${CHROMSIZES}" | grep -c "^chr")"
echo "Total entries: $(wc -l < "${CHROMSIZES}")"
echo "Standard chroms only:"
grep -c "^chr[0-9XY]*\b" "${CHROMSIZES}"

echo ""
echo "=== D6. CONDA ENVIRONMENT PACKAGES ==="
conda list 2>/dev/null | grep -E "^(snakemake|bedtools|samtools|macs2|htslib|tabix|bgzip|pysam|pybedtools|pandas|numpy|scipy|matplotlib|openpyxl|java|juicer)" | sort

echo ""
echo "=== D7. SNAKEMAKE VERSION + DRY-RUN TEST ==="
snakemake --version 2>&1
# Quick test: does snakemake even parse the config without errors?
cd "${ABC_REPO}"
snakemake -n --configfile config/config.yaml 2>&1 | tail -20

echo ""
echo "=== D8. JAVA (needed for juicertools / Hi-C extraction) ==="
java -version 2>&1
which java

echo ""
echo "=== D9. H3K27ac bigWig files (for future enhancement) ==="
ls -lh /expanse/lustre/projects/csd940/zalibhai/abc/input/H3K27ac_*.bw 2>/dev/null || echo "NOT FOUND"

echo ""
echo "=== D10. CHARACTERIZED LOOPS (for Step 9 cross-reference) ==="
LOOPS="/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/characterized_loops.tsv"
ls -la "${LOOPS}" 2>/dev/null || echo "FILE NOT FOUND"
if [ -f "${LOOPS}" ]; then
  echo "Columns:"
  head -1 "${LOOPS}" | tr '\t' '\n' | cat -n
  echo "Lines: $(wc -l < "${LOOPS}")"
  echo "First 3 data rows:"
  head -4 "${LOOPS}" | tail -3
fi

echo ""
echo "=== D11. H3K27ac PEAK BEDs (for Step 9 cross-reference) ==="
for f in H3K27acCerebellumLate2.bed H3K27acCerebellumEarly2.bed; do
  FULL="/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/peaks/beds/${f}"
  echo "--- ${f} ---"
  ls -la "${FULL}" 2>/dev/null || echo "NOT FOUND"
  if [ -f "${FULL}" ]; then
    echo "  Lines: $(wc -l < "${FULL}")"
    echo "  Columns: $(head -1 "${FULL}" | awk '{print NF}')"
    echo "  First 3 lines:"
    head -3 "${FULL}"
  fi
done

echo ""
echo "========================================================================"
echo "COMPLETE — $(date)"
echo "========================================================================"
echo ""
echo "DECISIONS PENDING FROM THIS OUTPUT:"
echo "  1. QNorm: Is K562 human qnorm appropriate for mouse ATAC?"
echo "     → If not, set use_qnorm: False and threshold: 0.02 manually"
echo "  2. Snakemake rules: Confirm symlink injection path from macs2.smk"
echo "  3. TagAlign sort order: Must match chrom.sizes chromosome order"
echo "  4. Java version: Must be compatible with juicertools for Hi-C"
echo "  5. characterized_loops.tsv column format: Needed for Step 9 implementation"
