# ABC Model Execution Plan — Detailed Implementation

> **Companion to:** `abc-prd-v2.md`  
> **Scope:** Steps with exact commands, validation, and scientific justification.  
> **Convention:** Each major step is a standalone SLURM script unless marked `[LOCAL]`.

---

## Table of Contents

1. [Step 0A: Gene Annotations from GENCODE vM25 GTF](#step-0a)
2. [Step 0B: Merge Blacklists](#step-0b)
3. [Step 0C: Fix TagAlign Compression + Index](#step-0c)
4. [Step 0D: Generate Ubiquitous Genes List](#step-0d)
5. [Step 0E: Format Consensus Peaks as narrowPeak](#step-0e)
6. [Step 1: Configure ABC Pipeline](#step-1)
7. [Step 2: Inject Consensus Peaks into Snakemake](#step-2)
8. [Step 3: Dry-Run Validation](#step-3)
9. [Step 4: Run ABC (WT and KO)](#step-4)
10. [Step 5: QC and Validation](#step-5)
11. [Step 6: Compute ΔABC with Unnormalized Scores](#step-6)
12. [Step 7: Integrate with RNA-seq](#step-7)
13. [Step 8: Gene-Level Aggregation](#step-8)
14. [Step 9: Cross-Reference with Loops and H3K27ac](#step-9)
15. [Validation Checklist](#validation-checklist)

---

## Resolved Decisions (from blocker verification)

| Decision | Resolution | Justification |
|----------|-----------|---------------|
| QNorm | **DISABLE** (`use_qnorm: False`) | Reference is human K562 DHS — cross-species AND cross-assay for mouse ATAC. Both conditions share same protocol, so no normalization needed for differential. |
| ABC threshold | **0.02 explicit** (not null/auto) | Auto would pick 0.021 (ATAC+no-H3K27ac+intact_hic). Disabling qnorm invalidates auto-calibration. 0.02 is close and matches PRD. |
| Gene annotation format | **8 columns with `#` header** | Must match bundled hg38 format: `#chr start end name score strand Ensembl_ID gene_type` |
| Gene biotype filter | **Exclude miRNA, pseudogenes, antisense** | Match bundled approach (keeps protein_coding + lincRNA + others) |
| Symlink injection path | `results/{biosample}/Peaks/macs2_peaks.narrowPeak` | Sort rule filters chroms against chrom.sizes and re-sorts — let it run |
| TagAlign sort order | Lexicographic is correct for tabix | Pipeline's sort rule re-orders by chrom.sizes internally |
| AWK version | GNU Awk 4.2.1 — `match()` with arrays works | No Perl fallback needed |

---

<a id="step-0a"></a>
## Step 0A: Gene Annotations from GENCODE vM25 GTF

### Purpose

ABC requires two reference files matching the bundled hg38 format exactly:
- **`mm10_genes.bed`**: 8 columns with header
- **`mm10_tss.bed`**: 500bp TSS regions, same 8-column format with header

### Format (must match exactly)

```
#chr    start   end     name    score   strand  Ensembl_ID      gene_type
chr1    3205901 3671498 Xkr4    0       -       ENSMUSG00000051951      protein_coding
```

### Script: `abc/scripts/step0a_gene_annotations.sh`

```bash
#!/bin/bash
#SBATCH --job-name=abc_gene_annot
#SBATCH --output=logs/abc_gene_annot_%j.out
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --account=csd940

mkdir -p /expanse/lustre/projects/csd940/zalibhai/abc/logs
mkdir -p /expanse/lustre/projects/csd940/zalibhai/abc/reference
source ~/.bashrc

REF_DIR="/expanse/lustre/projects/csd940/zalibhai/abc/reference"
GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz"
GTF_FILE="${REF_DIR}/gencode.vM25.annotation.gtf.gz"
CHROM_SIZES="/expanse/lustre/projects/csd940/zalibhai/taiji-new/reference/mm10/mm10.chrom.sizes"
GENES_OUT="${REF_DIR}/mm10_genes.bed"
TSS_OUT="${REF_DIR}/mm10_tss.bed"

# Biotypes to EXCLUDE (matching ABC bundled approach which filters these)
# ABC docs: "miRNAs, pseudogenes, antisense transcribed RNAs, bidirectional promoters"
EXCLUDE_BIOTYPES="miRNA|pseudogene|antisense|bidirectional_promoter"

# ---------- Download GTF ----------
if [ ! -f "${GTF_FILE}" ]; then
  echo "Downloading GENCODE vM25 GTF..."
  wget -O "${GTF_FILE}" "${GTF_URL}"
fi

echo "Verifying GTF integrity..."
GENE_COUNT=$(zcat "${GTF_FILE}" | awk '$3 == "gene"' | wc -l)
echo "  Total gene entries in GTF: ${GENE_COUNT}"
# Expected: ~55,000 for vM25

# ---------- Verify chr prefix ----------
echo "Checking chr prefix..."
zcat "${GTF_FILE}" | grep -v "^#" | head -1 | cut -f1
# GENCODE mouse GTF uses chr prefix — confirmed in GENCODE docs

# ---------- Extract gene annotation BED ----------
# Strategy:
#   - Filter to feature type "gene"
#   - Exclude problematic biotypes (miRNA, pseudogene, antisense, bidirectional)
#   - Extract 8 columns: chr, start, end, gene_name, 0, strand, gene_id, gene_type
#   - Strip Ensembl version suffix
#   - Add header matching bundled hg38 format

echo "Extracting gene annotations (excluding: ${EXCLUDE_BIOTYPES})..."

# Header (must match bundled format exactly)
echo -e "#chr\tstart\tend\tname\tscore\tstrand\tEnsembl_ID\tgene_type" > "${GENES_OUT}"

zcat "${GTF_FILE}" \
  | awk -F'\t' '$3 == "gene"' \
  | grep -v -E "gene_type \"(${EXCLUDE_BIOTYPES})" \
  | gawk -F'\t' 'BEGIN{OFS="\t"} {
      match($9, /gene_name "([^"]+)"/, gn)
      match($9, /gene_id "([^"]+)"/, gi)
      match($9, /gene_type "([^"]+)"/, gt)
      # Strip version suffix from Ensembl ID
      split(gi[1], eid, ".")
      print $1, $4-1, $5, gn[1], 0, $7, eid[1], gt[1]
    }' \
  | sort -k1,1 -k2,2n \
  >> "${GENES_OUT}"

# Count (excluding header)
N_GENES=$(tail -n +2 "${GENES_OUT}" | wc -l)
echo "  mm10_genes.bed: ${N_GENES} genes (after biotype filter)"

# ---------- Validate ----------
FIRST_CHR=$(tail -n +2 "${GENES_OUT}" | head -1 | cut -f1)
echo "  First chromosome entry: ${FIRST_CHR}"
if [[ ! "${FIRST_CHR}" =~ ^chr ]]; then
  echo "ERROR: Missing chr prefix in output. Aborting."
  exit 1
fi

# Report biotype distribution
echo "  Gene type distribution:"
tail -n +2 "${GENES_OUT}" | cut -f8 | sort | uniq -c | sort -rn | head -10

# Expect ~21,000 protein-coding + ~5,000-10,000 lincRNA + others
N_PC=$(tail -n +2 "${GENES_OUT}" | awk -F'\t' '$8 == "protein_coding"' | wc -l)
echo "  Protein-coding genes: ${N_PC}"
if [ "${N_PC}" -lt 20000 ] || [ "${N_PC}" -gt 25000 ]; then
  echo "WARNING: Unexpected protein-coding count (${N_PC}). Expected 20,000-25,000."
fi

# ---------- Check for duplicate gene symbols ----------
N_UNIQUE=$(tail -n +2 "${GENES_OUT}" | cut -f4 | sort -u | wc -l)
N_DUPS=$((N_GENES - N_UNIQUE))
echo "  Unique gene symbols: ${N_UNIQUE} (${N_DUPS} duplicates)"

if [ "${N_DUPS}" -gt 0 ]; then
  echo "  Resolving duplicates: keeping longest transcript per gene symbol..."
  # Preserve header, dedup body
  head -1 "${GENES_OUT}" > "${GENES_OUT}.dedup"
  tail -n +2 "${GENES_OUT}" | gawk -F'\t' 'BEGIN{OFS="\t"} {
    span = $3 - $2
    gene = $4
    if (!(gene in best) || span > best_span[gene]) {
      best[gene] = $0
      best_span[gene] = span
    }
  } END {
    for (g in best) print best[g]
  }' | sort -k1,1 -k2,2n >> "${GENES_OUT}.dedup"
  mv "${GENES_OUT}.dedup" "${GENES_OUT}"
  N_FINAL=$(tail -n +2 "${GENES_OUT}" | wc -l)
  echo "  After dedup: ${N_FINAL} genes"
fi

# ---------- Generate TSS file (500bp centered on TSS) ----------
echo "Generating 500bp TSS regions..."

# Header
echo -e "#chr\tstart\tend\tname\tscore\tstrand\tEnsembl_ID\tgene_type" > "${TSS_OUT}"

tail -n +2 "${GENES_OUT}" | awk -F'\t' 'BEGIN{OFS="\t"} {
  if ($6 == "+") {
    tss = $2
  } else {
    tss = $3
  }
  start = tss - 250
  if (start < 0) start = 0
  end = tss + 250
  print $1, start, end, $4, $5, $6, $7, $8
}' | sort -k1,1 -k2,2n >> "${TSS_OUT}"

N_TSS=$(tail -n +2 "${TSS_OUT}" | wc -l)
echo "  mm10_tss.bed: ${N_TSS} TSS regions"

# ---------- Cross-validate TSS within gene body ----------
echo "Validating TSS-gene body consistency..."
# Pair TSS and gene files (skip headers), check midpoint distance
paste <(tail -n +2 "${TSS_OUT}") <(tail -n +2 "${GENES_OUT}") \
  | awk -F'\t' '{
      tss_mid = int(($2 + $3) / 2)
      if ($6 == "+") {
        gene_tss = $9   # gene start (0-based)
      } else {
        gene_tss = $10  # gene end
      }
      diff = (tss_mid > gene_tss) ? tss_mid - gene_tss : gene_tss - tss_mid
      if (diff > 500) errors++
    } END {
      printf "  TSS > 500bp from gene boundary: %d / %d\n", errors+0, NR
    }'

# ---------- Verify column count ----------
NCOL_GENES=$(tail -n +2 "${GENES_OUT}" | head -1 | awk '{print NF}')
NCOL_TSS=$(tail -n +2 "${TSS_OUT}" | head -1 | awk '{print NF}')
echo "  Columns in genes file: ${NCOL_GENES} (expect 8)"
echo "  Columns in TSS file: ${NCOL_TSS} (expect 8)"

echo ""
echo "Step 0A complete."
echo "Outputs:"
echo "  ${GENES_OUT}"
echo "  ${TSS_OUT}"
```

### Validation Criteria
- Gene count: ~30,000-35,000 total (protein_coding + lincRNA + others after biotype filter)
- Protein-coding: 20,000-25,000
- All entries have `chr` prefix
- One entry per gene symbol (deduplicated)
- TSS midpoints within 500bp of gene boundaries
- 8 columns with `#` header line matching bundled format
- `Ensembl_ID` column has ENSMUSG IDs without version suffix

---

<a id="step-0b"></a>
## Step 0B: Merge Blacklists

### Purpose

Combine lab-specific (254 regions, BED3) and ENCODE (3,435 regions, BED6) blacklists.

### Verified formats
- Lab: BED3, no header, chr-prefixed
- ENCODE: BED6 (has annotation in cols 4-6), no header, chr-prefixed

### Commands `[LOCAL]`

```bash
LAB_BL="/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/tads/250123blacklist.bed"
ENCODE_BL="/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/tads/mm10-blacklist.v2.bed"
MERGED_BL="/expanse/lustre/projects/csd940/zalibhai/abc/reference/mm10_merged_blacklist.bed"

mkdir -p /expanse/lustre/projects/csd940/zalibhai/abc/reference

# Extract BED3 from both (ENCODE has 6 cols), cat, sort, merge
cat <(cut -f1-3 "${LAB_BL}") <(cut -f1-3 "${ENCODE_BL}") \
  | sort -k1,1 -k2,2n \
  | bedtools merge -i - \
  > "${MERGED_BL}"

echo "Merged blacklist: $(wc -l < "${MERGED_BL}") regions"
echo "  Lab-only:    $(bedtools subtract -a <(cut -f1-3 "${LAB_BL}") -b <(cut -f1-3 "${ENCODE_BL}") | wc -l) unique to lab"
echo "  ENCODE-only: $(bedtools subtract -a <(cut -f1-3 "${ENCODE_BL}") -b <(cut -f1-3 "${LAB_BL}") | wc -l) unique to ENCODE"
echo "  Overlapping: $(bedtools intersect -a <(cut -f1-3 "${LAB_BL}") -b <(cut -f1-3 "${ENCODE_BL}") -u | wc -l)"

# Check for chrM (common ATAC artifact region)
echo "  chrM entries: $(grep -c "^chrM" "${MERGED_BL}")"

# Verify chr prefix
echo "  First entry: $(head -1 "${MERGED_BL}")"
```

### Validation
- Merged count ≤ 254 + 3,435 = 3,689 (overlap reduces count)
- BED3 format, chr prefix, no header
- Matches bundled blacklist format (BED3, no header — confirmed from hg38 reference)

---

<a id="step-0c"></a>
## Step 0C: Fix TagAlign Compression + Index

### Verified status
- `ctrl_ATAC.tagAlign.gz`: 764MB, 246,139,345 reads, 6-column, lexicographic sort, gzip (not bgzip)
- `mut_ATAC.tagAlign.gz`: 749MB, 242,032,436 reads, same format
- No `.tbi` or `.bgz` files exist

### Script: `abc/scripts/step0c_fix_tagalign.sh`

```bash
#!/bin/bash
#SBATCH --job-name=abc_fix_tagalign
#SBATCH --output=logs/abc_fix_tagalign_%j.out
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --account=csd940

mkdir -p /expanse/lustre/projects/csd940/zalibhai/abc/logs
source ~/.bashrc
conda activate abc-env

INPUT_DIR="/expanse/lustre/projects/csd940/zalibhai/abc/input"

# Expected read counts from verification (for validation)
declare -A EXPECTED_COUNTS
EXPECTED_COUNTS[ctrl]=246139345
EXPECTED_COUNTS[mut]=242032436

for CONDITION in ctrl mut; do
  ORIG="${INPUT_DIR}/${CONDITION}_ATAC.tagAlign.gz"
  BGZIP="${INPUT_DIR}/${CONDITION}_ATAC.tagAlign.bgz"
  TBI="${BGZIP}.tbi"

  echo "=== Processing ${CONDITION} ==="

  if [ ! -f "${ORIG}" ]; then
    echo "ERROR: ${ORIG} not found"
    exit 1
  fi

  # Skip if already done
  if [ -f "${TBI}" ]; then
    echo "  ${BGZIP} and index already exist. Validating..."
    tabix -l "${BGZIP}" > /dev/null 2>&1
    if [ $? -eq 0 ]; then
      echo "  Valid. Skipping."
      continue
    else
      echo "  Index invalid. Regenerating..."
    fi
  fi

  # Sort lexicographically (required for tabix) and bgzip
  # TagAligns are already lexicographically sorted (verified: chr1, chr10, chr11...)
  # But re-sort to be safe — sort is idempotent
  echo "  Decompressing, sorting, and bgzipping..."
  echo "  (This will take several minutes for ~245M reads)"
  
  zcat "${ORIG}" \
    | sort -k1,1 -k2,2n -S 12G --parallel=4 \
    | bgzip -c -@ 4 \
    > "${BGZIP}"

  echo "  Creating tabix index..."
  tabix -p bed "${BGZIP}"

  if [ ! -f "${TBI}" ]; then
    echo "ERROR: Tabix index creation failed for ${BGZIP}"
    exit 1
  fi

  # Validate read count matches original
  echo "  Validating read count..."
  N_NEW=$(zcat "${BGZIP}" | wc -l)
  N_EXPECTED=${EXPECTED_COUNTS[${CONDITION}]}
  echo "  Original: ${N_EXPECTED} reads"
  echo "  New:      ${N_NEW} reads"
  
  if [ "${N_NEW}" -ne "${N_EXPECTED}" ]; then
    echo "ERROR: Read count mismatch! Expected ${N_EXPECTED}, got ${N_NEW}"
    exit 1
  fi

  # Quick tabix query validation
  N_CHR1=$(tabix "${BGZIP}" chr1:1-1000000 | wc -l)
  echo "  Tabix validation: ${N_CHR1} reads in chr1:1-1Mb"
  
  if [ "${N_CHR1}" -eq 0 ]; then
    echo "ERROR: Tabix query returned 0 reads. Index may be corrupt."
    exit 1
  fi

  echo "  Done."
  echo ""
done

echo "TagAlign conversion complete."
echo ""
echo "Files for ABC config (use these paths):"
echo "  WT ATAC: ${INPUT_DIR}/ctrl_ATAC.tagAlign.bgz"
echo "  KO ATAC: ${INPUT_DIR}/mut_ATAC.tagAlign.bgz"
echo ""
echo "IMPORTANT: config-biosamples-mm10.tsv must reference .bgz files, not .gz"
```

### Validation
- Read counts match exactly: ctrl=246,139,345 / mut=242,032,436
- `tabix -l file.bgz` returns chromosome list
- `tabix file.bgz chr1:1-1000000` returns non-zero reads

---

<a id="step-0d"></a>
## Step 0D: Generate Ubiquitous Genes List

### Verified format
Bundled file: one gene name per line, no header, uppercase human symbols (AARS, ABCF1, etc.)

### Script: `abc/scripts/step0d_ubiquitous_genes.py`

```python
#!/usr/bin/env python3
"""Generate ubiquitous genes list from DESeq2 baseMean.
Depends on Step 0A (needs mm10_genes.bed).

Bundled format: one gene name per line, no header.
"""

import pandas as pd
import sys

RNASEQ = "/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/tads/adult_timepoint_rna-seq-BAP1_WT_KO_v2_Results.xlsx"
GENES_BED = "/expanse/lustre/projects/csd940/zalibhai/abc/reference/mm10_genes.bed"
OUTPUT = "/expanse/lustre/projects/csd940/zalibhai/abc/reference/mm10_ubiquitous_genes.tsv"

# Read RNA-seq results (16,572 genes)
rna = pd.read_excel(RNASEQ)
print(f"RNA-seq entries: {len(rna)}")

# Column 'ensembl_gene_id' actually contains gene symbols (verified)
rna = rna.rename(columns={"ensembl_gene_id": "gene_symbol"})

# Read gene annotation (has header starting with #)
genes = pd.read_csv(GENES_BED, sep="\t", comment="#", header=None,
    names=["chr", "start", "end", "gene_symbol", "score", "strand", "ensembl_id", "gene_type"])
valid_genes = set(genes["gene_symbol"])
print(f"Genes in annotation: {len(valid_genes)}")

# Filter to genes present in our annotation
rna_matched = rna[rna["gene_symbol"].isin(valid_genes)].copy()
print(f"RNA-seq genes matching annotation: {len(rna_matched)}")

# Top 5% by baseMean = ubiquitously expressed
threshold = rna_matched["baseMean"].quantile(0.95)
ubiq = rna_matched[rna_matched["baseMean"] >= threshold]["gene_symbol"]
print(f"baseMean threshold (95th percentile): {threshold:.1f}")
print(f"Ubiquitous genes: {len(ubiq)}")

# Save: one gene per line, no header (matches bundled format)
ubiq.to_csv(OUTPUT, index=False, header=False)
print(f"Written to: {OUTPUT}")

# Spot-check known housekeeping genes
known_hk = ["Actb", "Gapdh", "Rpl13a", "Rps18", "Hsp90ab1"]
for g in known_hk:
    found = g in ubiq.values
    print(f"  {g}: {'FOUND' if found else 'MISSING'} in ubiquitous list")
```

### Validation
- ~800-1,100 genes (5% of matched genes)
- Contains known housekeeping genes (Actb, Gapdh, ribosomal)
- One gene per line, no header

---

<a id="step-0e"></a>
## Step 0E: Format Consensus Peaks as narrowPeak

### Verified status
- `consensus_all.bed`: 75,371 peaks, BED3, chr-prefixed, variable width (min=1, median=230, max=3457)
- Well under 150k → no top-N filtering needed
- 21 standard chromosomes only (chr1-19, chrX, chrY)

### Why use variable-width `consensus_all.bed` (not 500bp version)

The sort rule (`sort_narrowpeaks`) filters to chromosomes in chrom.sizes and re-sorts. Then `makeCandidateRegions.py` resizes peaks to 500bp centered on summit, removes blocklist regions, adds promoters, and merges overlaps. Feeding variable-width peaks with summit=midpoint lets the pipeline's own 500bp resizing work correctly.

### Script: `abc/scripts/step0e_consensus_to_narrowpeak.sh`

```bash
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
```

### Validation
- Exactly 75,371 peaks
- 10 columns per line
- Summit offsets positive and < peak width
- Median summit offset ~115 (half of median 230bp width)

---

<a id="step-1"></a>
## Step 1: Configure ABC Pipeline

### 1A. Biosample Table

Create `abc/ABC-Enhancer-Gene-Prediction/config/config-biosamples-mm10.tsv`:

**IMPORTANT:** This is a tab-separated file. All 10 columns must be present (including empty ones). Copy-paste from a text editor that preserves tabs, or create with the script below.

```bash
# Create biosample config with proper tabs
cat > /expanse/lustre/projects/csd940/zalibhai/abc/ABC-Enhancer-Gene-Prediction/config/config-biosamples-mm10.tsv << 'BIOSAMPLES'
biosample	DHS	ATAC	H3K27ac	default_accessibility_feature	HiC_file	HiC_type	HiC_resolution	alt_TSS	alt_genes
WT_cerebellum		/expanse/lustre/projects/csd940/zalibhai/abc/input/ctrl_ATAC.tagAlign.bgz		ATAC	/expanse/lustre/projects/csd940/zalibhai/nf-hic/250402_Bap1_deepseq/juicerpre/merged/hic/resorted_ctrl.hic	hic	5000		
KO_cerebellum		/expanse/lustre/projects/csd940/zalibhai/abc/input/mut_ATAC.tagAlign.bgz		ATAC	/expanse/lustre/projects/csd940/zalibhai/nf-hic/250402_Bap1_deepseq/juicerpre/merged/hic/resorted_mut.hic	hic	5000		
BIOSAMPLES

# Verify tab structure
echo "Column count per line:"
awk '{print NF}' /expanse/lustre/projects/csd940/zalibhai/abc/ABC-Enhancer-Gene-Prediction/config/config-biosamples-mm10.tsv
# All lines should show 10
```

### 1B. Config YAML

Overwrite `abc/ABC-Enhancer-Gene-Prediction/config/config.yaml`:

```bash
cat > /expanse/lustre/projects/csd940/zalibhai/abc/ABC-Enhancer-Gene-Prediction/config/config.yaml << 'CONFIGYAML'
# ABC config for mm10 BAP1-KO cerebellum analysis
# Modified from default hg38 config

### INPUT DATA
biosamplesTable: "config/config-biosamples-mm10.tsv"

### OUTPUT DATA
results_dir: "results/"

### REFERENCE FILES
ref:
        chrom_sizes: "/expanse/lustre/projects/csd940/zalibhai/taiji-new/reference/mm10/mm10.chrom.sizes"
        regions_blocklist: "/expanse/lustre/projects/csd940/zalibhai/abc/reference/mm10_merged_blacklist.bed"
        ubiquitous_genes: "/expanse/lustre/projects/csd940/zalibhai/abc/reference/mm10_ubiquitous_genes.tsv"
        genes: "/expanse/lustre/projects/csd940/zalibhai/abc/reference/mm10_genes.bed"
        genome_tss: "/expanse/lustre/projects/csd940/zalibhai/abc/reference/mm10_tss.bed"
        qnorm: "reference/EnhancersQNormRef.K562.txt"
        abc_thresholds: "reference/abc_thresholds.tsv"

### RULE SPECIFIC PARAMS
params_macs:
        pval: 0.1
        genome_size: mm   # CHANGED from hs — critical for MACS2 p-value calculation
        threads: 10

params_candidate:
        peakExtendFromSummit: 250    # produces 500bp regions (250 each side of summit)
        nStrongestPeaks: 150000      # our 75k peaks are under this, all retained

params_neighborhoods:
        # CHANGED to False: K562 human DHS qnorm reference is inappropriate
        # for mouse ATAC data (cross-species AND cross-assay)
        # Both conditions share same protocol, so differential is unaffected
        use_qnorm: False

params_predict:
        flags: "--scale_hic_using_powerlaw"
        hic_gamma: 1.024238616787792
        hic_scale: 5.9594510043736655
        hic_pseudocount_distance: 5000

params_filter_predictions:
        score_column: 'ABC.Score'
        threshold: 0.02              # CHANGED from null: explicit threshold per PRD
        include_self_promoter: True   # (auto-threshold invalid with qnorm disabled)
        only_expressed_genes: False

### INTERNAL USE ONLY
ABC_DIR_PATH: ""
CONFIGYAML

echo "Config written. Key changes from default:"
echo "  genome_size: mm (was hs)"
echo "  use_qnorm: False (was True)"
echo "  threshold: 0.02 (was null)"
echo "  All reference paths → mm10"
```

### ⚠️ Backup original config first

```bash
cp /expanse/lustre/projects/csd940/zalibhai/abc/ABC-Enhancer-Gene-Prediction/config/config.yaml \
   /expanse/lustre/projects/csd940/zalibhai/abc/ABC-Enhancer-Gene-Prediction/config/config.yaml.hg38_backup
```

---

<a id="step-2"></a>
## Step 2: Inject Consensus Peaks into Snakemake

### Verified pipeline chain

```
macs2.smk: call_macs_peaks → results/{biosample}/Peaks/macs2_peaks.narrowPeak
macs2.smk: sort_narrowpeaks → results/{biosample}/Peaks/macs2_peaks.narrowPeak.sorted
candidate_regions.smk: make_candidate_regions → ...sorted.candidateRegions.bed
```

### Strategy: Symlink at the unsorted narrowPeak level

We symlink our `consensus_all.narrowPeak` as `macs2_peaks.narrowPeak`. The `sort_narrowpeaks` rule then:
1. Filters to chromosomes present in chrom.sizes (keeps our 21 standard chroms, removes nothing since consensus is already clean)
2. Re-sorts by chrom.sizes order (size-descending: chr1, chr2, chrX, chr3...)

This is correct — the sort rule runs normally on our data.

### Commands

```bash
cd /expanse/lustre/projects/csd940/zalibhai/abc/ABC-Enhancer-Gene-Prediction

# Create output directories
mkdir -p results/WT_cerebellum/Peaks
mkdir -p results/KO_cerebellum/Peaks

# Symlink consensus narrowPeak as MACS2 output for BOTH conditions
# (same file — critical for consistent enhancer universe)
CONSENSUS_NP="/expanse/lustre/projects/csd940/zalibhai/abc/input/consensus_all.narrowPeak"

ln -sf "${CONSENSUS_NP}" results/WT_cerebellum/Peaks/macs2_peaks.narrowPeak
ln -sf "${CONSENSUS_NP}" results/KO_cerebellum/Peaks/macs2_peaks.narrowPeak

# Verify symlinks
ls -la results/WT_cerebellum/Peaks/macs2_peaks.narrowPeak
ls -la results/KO_cerebellum/Peaks/macs2_peaks.narrowPeak

# Tell Snakemake the MACS2 output already exists
# This prevents Snakemake from re-running MACS2 (which would overwrite our symlinks)
snakemake --touch results/WT_cerebellum/Peaks/macs2_peaks.narrowPeak \
                  results/KO_cerebellum/Peaks/macs2_peaks.narrowPeak
```

### ⚠️ Failure mode

If `snakemake --touch` doesn't prevent MACS2 from re-running (Snakemake may detect symlinks and re-trigger), use this fallback — copy instead of symlink:

```bash
cp "${CONSENSUS_NP}" results/WT_cerebellum/Peaks/macs2_peaks.narrowPeak
cp "${CONSENSUS_NP}" results/KO_cerebellum/Peaks/macs2_peaks.narrowPeak
```

---

<a id="step-3"></a>
## Step 3: Dry-Run Validation

### Commands

```bash
cd /expanse/lustre/projects/csd940/zalibhai/abc/ABC-Enhancer-Gene-Prediction
conda activate abc-env

# Full dry-run with printed commands
snakemake -n -p --configfile config/config.yaml 2>&1 | tee logs/dryrun_mm10.txt

# Summary table
snakemake -n --configfile config/config.yaml 2>&1 | tail -30
```

### What to Verify

1. **Rule count**: Expect ~17 rules (9 per biosample minus shared chrom_sizes_bed = ~17)
2. **`call_macs_peaks` should NOT appear** (satisfied by our narrowPeak)
3. **`sort_narrowpeaks` SHOULD appear** for both biosamples (sorts our consensus peaks)
4. **`make_candidate_regions` SHOULD appear** for both biosamples
5. **Reference paths** resolve correctly — grep for any remaining `hg38` references:
   ```bash
   grep -i "hg38\|GRCh38" logs/dryrun_mm10.txt
   # Should return NOTHING
   ```
6. **ATAC files** reference `.bgz` not `.gz`:
   ```bash
   grep "tagAlign" logs/dryrun_mm10.txt
   # Should show .bgz paths only
   ```
7. **H3K27ac is empty** in the neighborhoods command:
   ```bash
   grep "H3K27ac" logs/dryrun_mm10.txt
   # Should show --H3K27ac with empty value
   ```
8. **No --qnorm flag** in neighborhoods command (since use_qnorm: False)

### If MACS2 rule IS scheduled (unwanted)

```bash
# Check what triggered it
snakemake -n --configfile config/config.yaml --reason 2>&1 | grep macs

# If symlink timestamp issue, force-touch:
touch -r config/config.yaml results/*/Peaks/macs2_peaks.narrowPeak

# If still scheduled, copy the file instead of symlinking (see Step 2 fallback)
```

### Capture DAG

```bash
snakemake --dag --configfile config/config.yaml 2>/dev/null | dot -Tpdf > logs/abc_dag_mm10.pdf
# If dot is not available:
snakemake --dag --configfile config/config.yaml 2>/dev/null > logs/abc_dag_mm10.dot
```

---

<a id="step-4"></a>
## Step 4: Run ABC (WT and KO)

### Script: `abc/scripts/step4_run_abc.sh`

```bash
#!/bin/bash
#SBATCH --job-name=abc_run
#SBATCH --output=logs/abc_run_%j.out
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=40:00:00
#SBATCH --account=csd940

mkdir -p /expanse/lustre/projects/csd940/zalibhai/abc/logs
source ~/.bashrc
conda activate abc-env

cd /expanse/lustre/projects/csd940/zalibhai/abc/ABC-Enhancer-Gene-Prediction

echo "Starting ABC pipeline: $(date)"
echo "Using config: config/config.yaml"
echo "Snakemake version: $(snakemake --version)"
echo ""

snakemake \
  -j 16 \
  --configfile config/config.yaml \
  --rerun-incomplete \
  --printshellcmds \
  2>&1 | tee logs/abc_run_$(date +%Y%m%d_%H%M%S).log

EXIT_CODE=${PIPESTATUS[0]}

echo ""
echo "Pipeline finished with exit code: ${EXIT_CODE}"
echo "End time: $(date)"

# Quick output validation
for BIOSAMPLE in WT_cerebellum KO_cerebellum; do
  PRED_DIR="results/${BIOSAMPLE}/Predictions"
  echo ""
  echo "=== ${BIOSAMPLE} ==="
  
  if [ ! -d "${PRED_DIR}" ]; then
    echo "  ERROR: Predictions directory not found"
    continue
  fi

  # Thresholded predictions
  for f in "${PRED_DIR}"/EnhancerPredictions*.tsv; do
    if [ -f "$f" ]; then
      N=$(tail -n +2 "$f" | wc -l)
      echo "  $(basename $f): ${N} E-G pairs"
    fi
  done
  
  # Unthresholded
  ALLPUT="${PRED_DIR}/EnhancerPredictionsAllPutative.tsv.gz"
  if [ -f "${ALLPUT}" ]; then
    N=$(zcat "${ALLPUT}" | tail -n +2 | wc -l)
    echo "  AllPutative: ${N} E-G pairs (unthresholded)"
  else
    echo "  ERROR: AllPutative file not found"
  fi
  
  # QC metrics
  METRICS="results/${BIOSAMPLE}/Metrics"
  if [ -d "${METRICS}" ]; then
    echo "  QC metrics: $(ls "${METRICS}" | wc -l) files"
  fi
done

exit ${EXIT_CODE}
```

### Runtime Expectations
- 75k candidate elements, 2 biosamples, 16 cores
- Hi-C extraction: slowest step (reading 6.3GB + 5.9GB .hic files)
- Estimate: **4-12 hours**
- Memory: Hi-C extraction may spike to ~60-80 GB; 128 GB allocation provides headroom

### Monitoring

```bash
squeue -u zalibhai
tail -f /expanse/lustre/projects/csd940/zalibhai/abc/ABC-Enhancer-Gene-Prediction/logs/abc_run_*.log
grep -c "steps done\|Error\|error\|FAILED" logs/abc_run_*.log
```

---

<a id="step-5"></a>
## Step 5: QC and Validation

### Purpose

Before computing ΔABC, verify outputs are scientifically reasonable.

### Script: `abc/scripts/step5_qc.sh`

```bash
#!/bin/bash
#SBATCH --job-name=abc_qc
#SBATCH --output=logs/qc_%j.out
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --account=csd940

ABC_DIR="/expanse/lustre/projects/csd940/zalibhai/abc/ABC-Enhancer-Gene-Prediction"

# ===== CRITICAL: Inspect column headers =====
# This determines ALL column names for Steps 6-8
echo "=== COLUMN HEADERS (AllPutative) ==="
echo "--- WT ---"
zcat "${ABC_DIR}/results/WT_cerebellum/Predictions/EnhancerPredictionsAllPutative.tsv.gz" \
  | head -1 | tr '\t' '\n' | cat -n

echo ""
echo "--- KO (verify identical) ---"
zcat "${ABC_DIR}/results/KO_cerebellum/Predictions/EnhancerPredictionsAllPutative.tsv.gz" \
  | head -1 | tr '\t' '\n' | cat -n

echo ""
echo "=== SAMPLE DATA ROWS (first 3) ==="
zcat "${ABC_DIR}/results/WT_cerebellum/Predictions/EnhancerPredictionsAllPutative.tsv.gz" \
  | head -4

echo ""
echo "=== THRESHOLDED FILE HEADERS ==="
head -1 "${ABC_DIR}/results/WT_cerebellum/Predictions/EnhancerPredictions"*.tsv 2>/dev/null

for BIOSAMPLE in WT_cerebellum KO_cerebellum; do
  echo ""
  echo "========================================"
  echo "=== ${BIOSAMPLE} ==="
  echo "========================================"
  PRED="${ABC_DIR}/results/${BIOSAMPLE}/Predictions"
  
  # Find the thresholded file
  THRESH_FILE=$(ls "${PRED}"/EnhancerPredictions_threshold_*.tsv 2>/dev/null | head -1)
  if [ -z "${THRESH_FILE}" ]; then
    echo "ERROR: No thresholded prediction file found"
    continue
  fi
  
  echo "Thresholded file: $(basename "${THRESH_FILE}")"
  N_PAIRS=$(tail -n +2 "${THRESH_FILE}" | wc -l)
  echo "E-G pairs: ${N_PAIRS}"
  
  # ABC score distribution
  echo ""
  echo "--- ABC Score Distribution ---"
  # Need to find the ABC.Score column number from header
  ABC_COL=$(head -1 "${THRESH_FILE}" | tr '\t' '\n' | grep -n "ABC.Score" | head -1 | cut -d: -f1)
  echo "ABC.Score is column ${ABC_COL}"
  
  if [ -n "${ABC_COL}" ]; then
    tail -n +2 "${THRESH_FILE}" | cut -f"${ABC_COL}" | sort -g | awk '
      NR==1{min=$1} {sum+=$1; vals[NR]=$1}
      END{
        printf "  n=%d min=%.4f median=%.4f max=%.4f mean=%.4f\n",
          NR, min, vals[int(NR/2)], vals[NR], sum/NR
      }'
  fi

  # Enhancers per gene
  echo ""
  echo "--- Enhancers per Gene ---"
  GENE_COL=$(head -1 "${THRESH_FILE}" | tr '\t' '\n' | grep -n "TargetGene" | head -1 | cut -d: -f1)
  echo "TargetGene is column ${GENE_COL}"
  
  if [ -n "${GENE_COL}" ]; then
    N_GENES=$(tail -n +2 "${THRESH_FILE}" | cut -f"${GENE_COL}" | sort -u | wc -l)
    echo "  Unique genes with predictions: ${N_GENES}"
    echo "  Mean enhancers per gene: $(echo "scale=1; ${N_PAIRS} / ${N_GENES}" | bc)"
    
    echo "  Top 5 genes by enhancer count:"
    tail -n +2 "${THRESH_FILE}" | cut -f"${GENE_COL}" | sort | uniq -c | sort -rn | head -5
  fi
  
  # Per-chromosome distribution (check for missing chromosomes)
  echo ""
  echo "--- Per-Chromosome E-G Pairs ---"
  CHR_COL=$(head -1 "${THRESH_FILE}" | tr '\t' '\n' | grep -n "^chr$" | head -1 | cut -d: -f1)
  if [ -n "${CHR_COL}" ]; then
    tail -n +2 "${THRESH_FILE}" | cut -f"${CHR_COL}" | sort | uniq -c | sort -k2 -V
  else
    echo "  (chr column not found — check header)"
  fi
  
  # Check ABC scores sum to ~1 per gene (sample 5 genes)
  echo ""
  echo "--- ABC Score Sum per Gene (AllPutative, sample 5 genes) ---"
  ALLPUT="${PRED}/EnhancerPredictionsAllPutative.tsv.gz"
  if [ -f "${ALLPUT}" ]; then
    ALLPUT_GENE_COL=$(zcat "${ALLPUT}" | head -1 | tr '\t' '\n' | grep -n "TargetGene" | head -1 | cut -d: -f1)
    ALLPUT_ABC_COL=$(zcat "${ALLPUT}" | head -1 | tr '\t' '\n' | grep -n "ABC.Score" | head -1 | cut -d: -f1)
    
    # Get 5 random genes that have predictions
    SAMPLE_GENES=$(zcat "${ALLPUT}" | tail -n +2 | cut -f"${ALLPUT_GENE_COL}" | sort -u | shuf | head -5)
    for g in ${SAMPLE_GENES}; do
      SUM=$(zcat "${ALLPUT}" | awk -F'\t' -v col="${ALLPUT_ABC_COL}" -v gene="$g" -v gcol="${ALLPUT_GENE_COL}" \
        'NR>1 && $gcol == gene {sum += $col} END {printf "%.4f", sum}')
      N=$(zcat "${ALLPUT}" | awk -F'\t' -v gene="$g" -v gcol="${ALLPUT_GENE_COL}" \
        'NR>1 && $gcol == gene {n++} END {print n}')
      echo "  ${g}: sum(ABC)=${SUM} from ${N} enhancers"
    done
  fi
  
  # QC metrics directory
  echo ""
  echo "--- QC Metrics ---"
  ls "${ABC_DIR}/results/${BIOSAMPLE}/Metrics/" 2>/dev/null
done

echo ""
echo "========================================"
echo "CRITICAL: Copy the column header output above."
echo "Steps 6-8 column names MUST be updated to match."
echo "========================================"
```

### Expected Ranges
- ~30,000-60,000 E-G pairs per condition at threshold 0.02
- ~3 distal enhancers per expressed gene
- ABC scores sum to ~1 per gene in AllPutative
- All 21 standard chromosomes represented (chr1-19, chrX, chrY)

### ⚠️ STOP POINT

**Do not proceed to Step 6 until you have:**
1. Verified column headers from the output above
2. Updated column name variables in `step6_delta_abc.py` accordingly
3. Confirmed ABC scores sum to ~1 per gene

---

<a id="step-6"></a>
## Step 6: Compute ΔABC with Unnormalized Scores

### Purpose

For each enhancer-gene pair, compute:
- **ΔABC** = ABC_KO − ABC_WT (normalized — relative contribution change)
- **Δ(A×C)** = (Activity × Contact)_KO − (Activity × Contact)_WT (unnormalized — absolute regulatory change)

### Design Decisions

**Join key:** (enhancer_chr, enhancer_start, enhancer_end, target_gene). Coordinates match exactly since both conditions use identical consensus peaks.

**Use unthresholded predictions** (`AllPutative.tsv.gz`) for both conditions, then apply threshold post-join. This ensures pairs that drop below threshold in one condition are still captured.

### Script: `abc/scripts/step6_delta_abc.py`

```python
#!/usr/bin/env python3
"""
Compute ΔABC and Δ(Activity × Contact) between KO and WT conditions.

!! COLUMN NAMES BELOW ARE PLACEHOLDERS !!
!! UPDATE FROM STEP 5 OUTPUT BEFORE RUNNING !!

Usage: python step6_delta_abc.py
"""

import pandas as pd
import numpy as np
import sys
import os

# === CONFIGURATION ===
ABC_DIR = "/expanse/lustre/projects/csd940/zalibhai/abc/ABC-Enhancer-Gene-Prediction"
OUTPUT_DIR = "/expanse/lustre/projects/csd940/zalibhai/abc/results"
os.makedirs(OUTPUT_DIR, exist_ok=True)

WT_PRED = f"{ABC_DIR}/results/WT_cerebellum/Predictions/EnhancerPredictionsAllPutative.tsv.gz"
KO_PRED = f"{ABC_DIR}/results/KO_cerebellum/Predictions/EnhancerPredictionsAllPutative.tsv.gz"

# === LOAD DATA ===
print("Loading WT predictions...")
wt = pd.read_csv(WT_PRED, sep="\t")
print(f"  WT: {len(wt)} E-G pairs")
print(f"  Columns: {wt.columns.tolist()}")

print("Loading KO predictions...")
ko = pd.read_csv(KO_PRED, sep="\t")
print(f"  KO: {len(ko)} E-G pairs")

# === COLUMN NAME MAPPING ===
# !! VERIFY THESE FROM STEP 5 HEADER INSPECTION !!
# Typical ABC output columns include:
#   chr, start, end, name, class, TargetGene, TargetGeneEnsemblID,
#   TargetGeneTSS, TargetGeneExpression, TargetGenePromoterActivityQuantile,
#   TargetGeneIsExpressed, distance, isSelfPromoter, hic_contact,
#   powerlaw_contact, activity_base, TargetGenePromoterActivity,
#   ABC.Score, ABC.Score.Numerator, ...
#
# UPDATE THE FOLLOWING VARIABLES after inspecting Step 5 output:

ENHANCER_COLS = ["chr", "start", "end"]   # Enhancer coordinates
GENE_COL = "TargetGene"                    # Gene symbol column
ABC_SCORE_COL = "ABC.Score"                # Normalized ABC score
UNNORM_COL = "ABC.Score.Numerator"         # Activity × Contact (numerator)
# Fallback if ABC.Score.Numerator doesn't exist:
ACTIVITY_COL = "activity_base"             # Activity component
CONTACT_COL = "hic_contact"               # Contact component (Hi-C specific)

JOIN_COLS = ENHANCER_COLS + [GENE_COL]

# === VERIFY COLUMNS EXIST ===
for col in JOIN_COLS + [ABC_SCORE_COL]:
    if col not in wt.columns:
        print(f"ERROR: Column '{col}' not found in WT predictions.")
        print(f"  Available: {wt.columns.tolist()}")
        sys.exit(1)

# === COMPUTE UNNORMALIZED SCORE ===
for df, label in [(wt, "WT"), (ko, "KO")]:
    if UNNORM_COL in df.columns:
        print(f"  {label}: Using existing {UNNORM_COL} column")
        df["unnorm_score"] = df[UNNORM_COL]
    elif ACTIVITY_COL in df.columns and CONTACT_COL in df.columns:
        print(f"  {label}: Computing {ACTIVITY_COL} × {CONTACT_COL}")
        df["unnorm_score"] = df[ACTIVITY_COL] * df[CONTACT_COL]
    else:
        print(f"ERROR: Cannot compute unnormalized score for {label}.")
        print(f"  Need '{UNNORM_COL}' or both '{ACTIVITY_COL}' and '{CONTACT_COL}'")
        print(f"  Available: {df.columns.tolist()}")
        sys.exit(1)

# === OUTER JOIN ===
print("\nJoining WT and KO predictions on enhancer-gene pairs...")

# Select only needed columns to reduce memory
wt_slim = wt[JOIN_COLS + [ABC_SCORE_COL, "unnorm_score"]].copy()
ko_slim = ko[JOIN_COLS + [ABC_SCORE_COL, "unnorm_score"]].copy()

merged = pd.merge(
    wt_slim, ko_slim,
    on=JOIN_COLS,
    how="outer",
    suffixes=("_WT", "_KO"),
    indicator=True
)

merge_counts = merged["_merge"].value_counts()
print(f"  Both conditions: {merge_counts.get('both', 0)}")
print(f"  WT only:         {merge_counts.get('left_only', 0)}")
print(f"  KO only:         {merge_counts.get('right_only', 0)}")

total_pairs = len(merged)
both_count = merge_counts.get('both', 0)
pct_both = 100 * both_count / total_pairs if total_pairs > 0 else 0
print(f"  Overlap: {pct_both:.1f}% of pairs present in both conditions")

if pct_both < 90:
    print("WARNING: <90% overlap. Check if candidate regions differ between conditions.")

# Fill missing with 0
for col in [f"{ABC_SCORE_COL}_WT", f"{ABC_SCORE_COL}_KO", "unnorm_score_WT", "unnorm_score_KO"]:
    merged[col] = merged[col].fillna(0)

# === COMPUTE DELTAS ===
merged["delta_ABC"] = merged[f"{ABC_SCORE_COL}_KO"] - merged[f"{ABC_SCORE_COL}_WT"]
merged["delta_unnorm"] = merged["unnorm_score_KO"] - merged["unnorm_score_WT"]

# === APPLY THRESHOLD ===
# Keep pairs where EITHER condition has ABC ≥ 0.02
mask = (merged[f"{ABC_SCORE_COL}_WT"] >= 0.02) | (merged[f"{ABC_SCORE_COL}_KO"] >= 0.02)
filtered = merged[mask].copy()
print(f"\nTotal E-G pairs (union): {len(merged)}")
print(f"Pairs with ABC ≥ 0.02 in at least one condition: {len(filtered)}")

# Drop merge indicator
filtered = filtered.drop(columns=["_merge"])

# === SUMMARY STATISTICS ===
print("\n=== ΔABC Summary ===")
print(f"  Mean:   {filtered['delta_ABC'].mean():.5f}")
print(f"  Median: {filtered['delta_ABC'].median():.5f}")
print(f"  Std:    {filtered['delta_ABC'].std():.5f}")
print(f"  Gained (ΔABC > 0.01):    {(filtered['delta_ABC'] > 0.01).sum()}")
print(f"  Lost (ΔABC < -0.01):     {(filtered['delta_ABC'] < -0.01).sum()}")
print(f"  Unchanged (|Δ| ≤ 0.01): {(filtered['delta_ABC'].abs() <= 0.01).sum()}")

print("\n=== Δ(Activity × Contact) Summary ===")
print(f"  Mean:   {filtered['delta_unnorm'].mean():.5f}")
print(f"  Median: {filtered['delta_unnorm'].median():.5f}")

# Check for asymmetry (global gain/loss of enhancer activity in KO)
if abs(filtered['delta_ABC'].mean()) > 0.005:
    print(f"\n⚠️  ΔABC mean is notably non-zero ({filtered['delta_ABC'].mean():.5f}).")
    print("   This may indicate global activity shift in KO vs WT.")
    print("   The unnormalized Δ(A×C) may be more informative for directionality.")

# === SAVE ===
out_all = f"{OUTPUT_DIR}/delta_abc_all_pairs.tsv"
out_sig = f"{OUTPUT_DIR}/delta_abc_significant.tsv"

filtered.to_csv(out_all, sep="\t", index=False)
print(f"\nSaved all thresholded pairs: {out_all} ({len(filtered)} pairs)")

sig = filtered[filtered["delta_ABC"].abs() > 0.01]
sig.to_csv(out_sig, sep="\t", index=False)
print(f"Saved significant pairs (|ΔABC| > 0.01): {out_sig} ({len(sig)} pairs)")
```

### Validation
- No NaN values in delta columns
- >90% of E-G pairs present in both conditions
- ΔABC distribution roughly symmetric (if not, note this — it's biologically informative)
- Files written without error

---

<a id="step-7"></a>
## Step 7: Integrate with RNA-seq

### Verified RNA-seq format
- 16,572 rows, gene symbols in `ensembl_gene_id` column
- Has: `baseMean`, `log2FoldChange`, `padj`
- 7 ctrl + 7 mut sample columns

### Script: `abc/scripts/step7_rnaseq_integration.py`

```python
#!/usr/bin/env python3
"""
Integrate ΔABC with RNA-seq log2FC.
Produces scatter plots and concordance analysis.

Usage: python step7_rnaseq_integration.py
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend for cluster
import matplotlib.pyplot as plt
from scipy import stats
import os
import warnings
warnings.filterwarnings('ignore', category=FutureWarning)

# === CONFIGURATION ===
DELTA_ABC = "/expanse/lustre/projects/csd940/zalibhai/abc/results/delta_abc_all_pairs.tsv"
RNASEQ = "/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/tads/adult_timepoint_rna-seq-BAP1_WT_KO_v2_Results.xlsx"
OUTPUT_DIR = "/expanse/lustre/projects/csd940/zalibhai/abc/results"
GENE_COL = "TargetGene"  # UPDATE from Step 5 if different

os.makedirs(f"{OUTPUT_DIR}/figures", exist_ok=True)

# === LOAD DATA ===
delta = pd.read_csv(DELTA_ABC, sep="\t")
rna = pd.read_excel(RNASEQ)
rna = rna.rename(columns={"ensembl_gene_id": "gene_symbol"})

print(f"ΔABC pairs: {len(delta)}")
print(f"ΔABC unique genes: {delta[GENE_COL].nunique()}")
print(f"RNA-seq genes: {len(rna)}")

# === JOIN ===
merged = pd.merge(
    delta,
    rna[["gene_symbol", "log2FoldChange", "padj", "baseMean"]],
    left_on=GENE_COL,
    right_on="gene_symbol",
    how="inner"
)

n_genes_matched = merged["gene_symbol"].nunique()
n_genes_abc = delta[GENE_COL].nunique()
n_genes_rna = rna["gene_symbol"].nunique()
print(f"\nGenes in ABC predictions: {n_genes_abc}")
print(f"Genes in RNA-seq: {n_genes_rna}")
print(f"Genes matched: {n_genes_matched} ({100*n_genes_matched/n_genes_abc:.1f}% of ABC genes)")
print(f"E-G pairs with RNA-seq data: {len(merged)}")

# === CONCORDANCE ANALYSIS ===
ABC_THRESH = 0.01
RNA_PADJ = 0.05
RNA_LFC = 0.5

merged["abc_changed"] = merged["delta_ABC"].abs() > ABC_THRESH
merged["rna_changed"] = (merged["padj"] < RNA_PADJ) & (merged["log2FoldChange"].abs() > RNA_LFC)

both_changed = merged[merged["abc_changed"] & merged["rna_changed"]].copy()
if len(both_changed) > 0:
    both_changed["concordant"] = (
        np.sign(both_changed["delta_ABC"]) == np.sign(both_changed["log2FoldChange"])
    )
    
    n_conc = both_changed["concordant"].sum()
    n_disc = (~both_changed["concordant"]).sum()
    n_tot = len(both_changed)
    
    print(f"\n=== Concordance (both ABC and RNA-seq significantly changed) ===")
    print(f"  E-G pairs with both changed: {n_tot}")
    print(f"  Concordant (same direction): {n_conc} ({100*n_conc/n_tot:.1f}%)")
    print(f"  Discordant: {n_disc} ({100*n_disc/n_tot:.1f}%)")
    
    # Binomial test
    binom_p = stats.binomtest(n_conc, n_tot, 0.5).pvalue
    print(f"  Binomial test p-value (vs 50%): {binom_p:.2e}")
else:
    print("\nWARNING: No E-G pairs with both ABC and RNA-seq significantly changed.")
    print("  Try relaxing thresholds.")

# === GENE-LEVEL SUMMARY (strongest ΔABC per gene) ===
gene_summary = merged.groupby("gene_symbol").agg(
    max_delta_abc=("delta_ABC", lambda x: x.iloc[x.abs().argmax()]),
    max_delta_unnorm=("delta_unnorm", lambda x: x.iloc[x.abs().argmax()]),
    log2FC=("log2FoldChange", "first"),
    padj=("padj", "first"),
    baseMean=("baseMean", "first"),
    n_enhancers=("delta_ABC", "size")
).reset_index()

# === SCATTER PLOTS ===
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

sig_mask = gene_summary["padj"] < RNA_PADJ

# Panel A: Normalized ΔABC
ax = axes[0]
ax.scatter(gene_summary.loc[~sig_mask, "max_delta_abc"],
           gene_summary.loc[~sig_mask, "log2FC"],
           alpha=0.3, s=10, c="gray", label=f"padj ≥ {RNA_PADJ}", rasterized=True)
ax.scatter(gene_summary.loc[sig_mask, "max_delta_abc"],
           gene_summary.loc[sig_mask, "log2FC"],
           alpha=0.5, s=15, c="crimson", label=f"padj < {RNA_PADJ}", rasterized=True)
ax.axhline(0, color="black", linewidth=0.5, linestyle="--")
ax.axvline(0, color="black", linewidth=0.5, linestyle="--")
ax.set_xlabel("ΔABC (strongest enhancer per gene)")
ax.set_ylabel("log2FC (KO vs WT)")
ax.set_title("Normalized ABC Score")
ax.legend(fontsize=8)

de = gene_summary[sig_mask & gene_summary["max_delta_abc"].notna()]
if len(de) > 2:
    r, p = stats.pearsonr(de["max_delta_abc"], de["log2FC"])
    ax.text(0.05, 0.95, f"r = {r:.3f}, p = {p:.2e}\n(DE genes, n={len(de)})",
            transform=ax.transAxes, fontsize=8, va="top",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))

# Panel B: Unnormalized Δ(A×C)
ax = axes[1]
ax.scatter(gene_summary.loc[~sig_mask, "max_delta_unnorm"],
           gene_summary.loc[~sig_mask, "log2FC"],
           alpha=0.3, s=10, c="gray", label=f"padj ≥ {RNA_PADJ}", rasterized=True)
ax.scatter(gene_summary.loc[sig_mask, "max_delta_unnorm"],
           gene_summary.loc[sig_mask, "log2FC"],
           alpha=0.5, s=15, c="steelblue", label=f"padj < {RNA_PADJ}", rasterized=True)
ax.axhline(0, color="black", linewidth=0.5, linestyle="--")
ax.axvline(0, color="black", linewidth=0.5, linestyle="--")
ax.set_xlabel("Δ(Activity × Contact) (strongest enhancer per gene)")
ax.set_ylabel("log2FC (KO vs WT)")
ax.set_title("Unnormalized Score")
ax.legend(fontsize=8)

if len(de) > 2:
    r2, p2 = stats.pearsonr(de["max_delta_unnorm"], de["log2FC"])
    ax.text(0.05, 0.95, f"r = {r2:.3f}, p = {p2:.2e}\n(DE genes, n={len(de)})",
            transform=ax.transAxes, fontsize=8, va="top",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))

plt.tight_layout()
fig.savefig(f"{OUTPUT_DIR}/figures/delta_abc_vs_log2fc.pdf", dpi=300, bbox_inches="tight")
fig.savefig(f"{OUTPUT_DIR}/figures/delta_abc_vs_log2fc.png", dpi=150, bbox_inches="tight")
print(f"\nSaved: {OUTPUT_DIR}/figures/delta_abc_vs_log2fc.pdf")

# === SAVE ===
merged.to_csv(f"{OUTPUT_DIR}/delta_abc_with_rnaseq.tsv", sep="\t", index=False)
gene_summary.to_csv(f"{OUTPUT_DIR}/gene_level_summary.tsv", sep="\t", index=False)
print(f"Saved: {OUTPUT_DIR}/delta_abc_with_rnaseq.tsv ({len(merged)} pairs)")
print(f"Saved: {OUTPUT_DIR}/gene_level_summary.tsv ({len(gene_summary)} genes)")
```

---

<a id="step-8"></a>
## Step 8: Gene-Level Aggregation

### Script: `abc/scripts/step8_gene_aggregation.py`

```python
#!/usr/bin/env python3
"""
Gene-level aggregation of ABC predictions.
Extends the gene_summary from Step 7 with additional metrics.

Usage: python step8_gene_aggregation.py
"""

import pandas as pd
import numpy as np
import os

OUTPUT_DIR = "/expanse/lustre/projects/csd940/zalibhai/abc/results"
MERGED = f"{OUTPUT_DIR}/delta_abc_with_rnaseq.tsv"

# UPDATE these column names from Step 5
ABC_SCORE_WT = "ABC.Score_WT"
ABC_SCORE_KO = "ABC.Score_KO"

merged = pd.read_csv(MERGED, sep="\t")
print(f"Loaded: {len(merged)} E-G pairs")

gene_agg = merged.groupby("gene_symbol").agg(
    max_delta_abc=("delta_ABC", lambda x: x.iloc[x.abs().argmax()]),
    mean_delta_abc=("delta_ABC", "mean"),
    sum_abc_wt=(ABC_SCORE_WT, "sum"),
    sum_abc_ko=(ABC_SCORE_KO, "sum"),
    sum_delta_unnorm=("delta_unnorm", "sum"),
    n_enhancers_total=("delta_ABC", "size"),
    n_enhancers_gained=("delta_ABC", lambda x: (x > 0.01).sum()),
    n_enhancers_lost=("delta_ABC", lambda x: (x < -0.01).sum()),
    log2FC=("log2FoldChange", "first"),
    padj=("padj", "first"),
    baseMean=("baseMean", "first")
).reset_index()

# Dysregulated: significant ΔABC AND significant DE
gene_agg["dysregulated"] = (
    (gene_agg["max_delta_abc"].abs() > 0.01) &
    (gene_agg["padj"] < 0.05) &
    (gene_agg["log2FC"].abs() > 0.5)
)

n_dys = gene_agg["dysregulated"].sum()
print(f"\nDysregulated genes (|ΔABC|>0.01 AND padj<0.05 AND |log2FC|>0.5): {n_dys}")

# Breakdown
dys = gene_agg[gene_agg["dysregulated"]]
n_up = ((dys["max_delta_abc"] > 0) & (dys["log2FC"] > 0)).sum()
n_down = ((dys["max_delta_abc"] < 0) & (dys["log2FC"] < 0)).sum()
n_disc = n_dys - n_up - n_down
print(f"  Gained enhancer + upregulated: {n_up}")
print(f"  Lost enhancer + downregulated: {n_down}")
print(f"  Discordant: {n_disc}")

# Top dysregulated genes by |ΔABC|
print(f"\nTop 20 dysregulated genes by |ΔABC|:")
top = dys.nlargest(20, "max_delta_abc", key=abs)[
    ["gene_symbol", "max_delta_abc", "log2FC", "padj", "n_enhancers_total",
     "n_enhancers_gained", "n_enhancers_lost"]
]
print(top.to_string(index=False))

gene_agg.to_csv(f"{OUTPUT_DIR}/gene_level_aggregated.tsv", sep="\t", index=False)
print(f"\nSaved: {OUTPUT_DIR}/gene_level_aggregated.tsv ({len(gene_agg)} genes)")
```

---

<a id="step-9"></a>
## Step 9: Cross-Reference with Loops and H3K27ac

### 9A. Cross-Reference with Differential Loops

### Verified loop format
- 2,911 loops with 57 columns
- Key columns: `anchor1_chr/start/end`, `anchor2_chr/start/end`, `anchor1_nearest_gene`, `anchor2_nearest_gene`, `logFC`, `FDR`, `significant`, `direction`, `loop_type`

### Script: `abc/scripts/step9_cross_reference.py`

```python
#!/usr/bin/env python3
"""
Cross-reference ABC predictions with:
  A) Differential Hi-C loops (characterized_loops.tsv)
  B) H3K27ac peaks (H3K27acCerebellumLate2.bed)

Usage: python step9_cross_reference.py
"""

import pandas as pd
import numpy as np
import subprocess
import tempfile
import os

OUTPUT_DIR = "/expanse/lustre/projects/csd940/zalibhai/abc/results"
LOOPS = "/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/characterized_loops.tsv"
H3K27AC = "/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/peaks/beds/H3K27acCerebellumLate2.bed"
DELTA_ABC = f"{OUTPUT_DIR}/delta_abc_all_pairs.tsv"

# UPDATE column names from Step 5
ENHANCER_CHR = "chr"
ENHANCER_START = "start"
ENHANCER_END = "end"
GENE_COL = "TargetGene"

delta = pd.read_csv(DELTA_ABC, sep="\t")
print(f"ΔABC pairs: {len(delta)}")

# =====================================================================
# 9A: Overlap with differential loops
# =====================================================================
print("\n=== 9A: Differential Loop Overlap ===")

loops = pd.read_csv(LOOPS, sep="\t")
sig_loops = loops[loops["significant"] == True].copy()
print(f"Total loops: {len(loops)}")
print(f"Significant loops: {len(sig_loops)}")

# Strategy: For each ABC enhancer-gene pair, check if:
#   - The enhancer overlaps a loop anchor
#   - The loop's OTHER anchor is near the gene's promoter
# Use bedtools intersect via temp files

# Extract unique enhancers from ΔABC
enhancers = delta[[ENHANCER_CHR, ENHANCER_START, ENHANCER_END]].drop_duplicates()
print(f"Unique ABC enhancers: {len(enhancers)}")

# Write temp BED files
with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
    enh_bed = f.name
    enhancers.to_csv(f, sep="\t", header=False, index=False)

# Loop anchor1 and anchor2 as BED files
with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
    anchor1_bed = f.name
    sig_loops[["anchor1_chr", "anchor1_start", "anchor1_end"]].to_csv(
        f, sep="\t", header=False, index=False)

with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
    anchor2_bed = f.name
    sig_loops[["anchor2_chr", "anchor2_start", "anchor2_end"]].to_csv(
        f, sep="\t", header=False, index=False)

# Count ABC enhancers overlapping any significant loop anchor
for label, anchor_bed in [("anchor1", anchor1_bed), ("anchor2", anchor2_bed)]:
    result = subprocess.run(
        ["bedtools", "intersect", "-a", enh_bed, "-b", anchor_bed, "-u"],
        capture_output=True, text=True
    )
    n_overlap = len(result.stdout.strip().split("\n")) if result.stdout.strip() else 0
    print(f"  ABC enhancers overlapping {label} of sig loops: {n_overlap} / {len(enhancers)}")

# Gene-level: do ABC target genes appear in loop annotations?
abc_genes = set(delta[GENE_COL].unique())
loop_genes = set(sig_loops["anchor1_nearest_gene"].dropna()) | set(sig_loops["anchor2_nearest_gene"].dropna())
overlap_genes = abc_genes & loop_genes
print(f"  ABC target genes also in sig loop annotations: {len(overlap_genes)} / {len(abc_genes)}")

# Directional concordance: for overlapping genes, does ΔABC direction match loop direction?
loop_gene_direction = {}
for _, row in sig_loops.iterrows():
    for g in [row.get("anchor1_nearest_gene"), row.get("anchor2_nearest_gene")]:
        if pd.notna(g) and g in abc_genes:
            loop_gene_direction[g] = row["direction"]  # "up_in_mutant" or "down_in_mutant"

gene_summary = delta.groupby(GENE_COL).agg(
    max_delta_abc=("delta_ABC", lambda x: x.iloc[x.abs().argmax()])
).reset_index()

conc_count = 0
disc_count = 0
for g, direction in loop_gene_direction.items():
    row = gene_summary[gene_summary[GENE_COL] == g]
    if len(row) == 0:
        continue
    dabc = row.iloc[0]["max_delta_abc"]
    if direction == "up_in_mutant" and dabc > 0:
        conc_count += 1
    elif direction == "down_in_mutant" and dabc < 0:
        conc_count += 1
    else:
        disc_count += 1

total_dir = conc_count + disc_count
if total_dir > 0:
    print(f"  Directional concordance (loop direction vs ΔABC): {conc_count}/{total_dir} "
          f"({100*conc_count/total_dir:.1f}%)")

# Cleanup temp files
for f in [enh_bed, anchor1_bed, anchor2_bed]:
    os.unlink(f)

# =====================================================================
# 9B: H3K27ac peak overlap
# =====================================================================
print("\n=== 9B: H3K27ac Peak Overlap ===")

# Write enhancers BED
with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
    enh_bed = f.name
    enhancers.to_csv(f, sep="\t", header=False, index=False)

# Overlap with H3K27ac peaks (15,105 peaks, BED6)
result_u = subprocess.run(
    ["bedtools", "intersect", "-a", enh_bed, "-b", H3K27AC, "-u"],
    capture_output=True, text=True
)
result_v = subprocess.run(
    ["bedtools", "intersect", "-a", enh_bed, "-b", H3K27AC, "-v"],
    capture_output=True, text=True
)

n_with_k27ac = len(result_u.stdout.strip().split("\n")) if result_u.stdout.strip() else 0
n_without_k27ac = len(result_v.stdout.strip().split("\n")) if result_v.stdout.strip() else 0
print(f"  ABC enhancers overlapping H3K27ac peaks: {n_with_k27ac} / {len(enhancers)}")
print(f"  ABC enhancers WITHOUT H3K27ac: {n_without_k27ac} / {len(enhancers)}")

# Parse overlap results to annotate ΔABC table
overlap_enhancers = set()
if result_u.stdout.strip():
    for line in result_u.stdout.strip().split("\n"):
        parts = line.split("\t")
        overlap_enhancers.add((parts[0], int(parts[1]), int(parts[2])))

delta["has_H3K27ac"] = delta.apply(
    lambda r: (r[ENHANCER_CHR], int(r[ENHANCER_START]), int(r[ENHANCER_END])) in overlap_enhancers,
    axis=1
)

# Compare ABC scores for enhancers with vs without H3K27ac
with_k27 = delta[delta["has_H3K27ac"]]
without_k27 = delta[~delta["has_H3K27ac"]]
print(f"\n  Mean |ΔABC| for H3K27ac+ enhancers: {with_k27['delta_ABC'].abs().mean():.5f}")
print(f"  Mean |ΔABC| for H3K27ac- enhancers: {without_k27['delta_ABC'].abs().mean():.5f}")

# Save annotated table
delta.to_csv(f"{OUTPUT_DIR}/delta_abc_h3k27ac_annotated.tsv", sep="\t", index=False)
print(f"\nSaved: {OUTPUT_DIR}/delta_abc_h3k27ac_annotated.tsv")

os.unlink(enh_bed)
```

---

<a id="validation-checklist"></a>
## Validation Checklist

### Pre-Pipeline (Steps 0A-0E)

- [ ] `mm10_genes.bed`: 8 cols with `#` header, ~30k genes, chr prefix, one per symbol, protein_coding count 20-25k
- [ ] `mm10_tss.bed`: Same count, 500bp regions, TSS within gene body, 8 cols with header
- [ ] `mm10_merged_blacklist.bed`: ~3,500+ regions, BED3, chr prefix
- [ ] `ctrl_ATAC.tagAlign.bgz` + `.tbi`: 246,139,345 reads, tabix queries work
- [ ] `mut_ATAC.tagAlign.bgz` + `.tbi`: 242,032,436 reads, tabix queries work
- [ ] `mm10_ubiquitous_genes.tsv`: ~800-1,100 genes, one per line, includes Actb/Gapdh
- [ ] `consensus_all.narrowPeak`: 75,371 peaks, 10 columns

### Pipeline Configuration (Steps 1-3)

- [ ] `config-biosamples-mm10.tsv`: 10 columns per line, .bgz paths, ATAC-only
- [ ] `config.yaml`: `genome_size: mm`, `use_qnorm: False`, `threshold: 0.02`, all mm10 reference paths
- [ ] `config.yaml.hg38_backup` exists
- [ ] Dry-run: `call_macs_peaks` NOT scheduled, all other rules resolve
- [ ] No `hg38`/`GRCh38` references in dry-run output
- [ ] .bgz (not .gz) paths in dry-run output

### ABC Outputs (Steps 4-5)

- [ ] `EnhancerPredictionsAllPutative.tsv.gz` exists for both biosamples
- [ ] `EnhancerPredictions_threshold_0.02.tsv` exists for both
- [ ] ~30,000-60,000 E-G pairs per condition
- [ ] ~3 enhancers per expressed gene
- [ ] ABC scores sum to ~1 per gene
- [ ] All 21 chromosomes represented
- [ ] Column headers inspected and documented

### ΔABC and Integration (Steps 6-9)

- [ ] >90% E-G pairs present in both conditions
- [ ] ΔABC distribution inspected for symmetry
- [ ] Concordance with RNA-seq: reported (ideally >50%, binomial test)
- [ ] Gene-level table: no NaN in key fields
- [ ] Dysregulated gene count: biologically plausible (dozens to low hundreds)
- [ ] H3K27ac overlap annotated
- [ ] Loop cross-reference completed

---

## Execution Order

```
Step 0A: Gene annotations ──┐
Step 0B: Merge blacklists ──┤
Step 0C: Fix tagAligns ─────┼── Parallel (0D depends on 0A)
Step 0D: Ubiquitous genes ──┤
Step 0E: narrowPeak format ─┘
         │
Step 1:  Config files (depends on all Step 0)
Step 2:  Inject consensus peaks
Step 3:  Dry-run ← GATE: must pass before Step 4
         │
Step 4:  Run ABC pipeline ← COMPUTE: 4-12 hours, 128GB
         │
Step 5:  QC ← GATE: inspect headers, update Step 6-9 column names
         │
Step 6:  ΔABC computation
Step 7:  RNA-seq integration ─┐
Step 8:  Gene aggregation ────┼── Can combine into single script
Step 9:  Cross-reference ─────┘
```
