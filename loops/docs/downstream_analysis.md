# Downstream Analysis: Merged Loops & Characterization

## Script: `downstream_analysis.R`

## Overview

This script performs **post-processing integration** of differential loop results across multiple resolutions (5kb, 10kb, 25kb). It creates a unified, non-redundant set of differential chromatin loops and enriches them with **gene annotations** and **ChIP-seq-based functional classifications**.

## Purpose

After running differential analysis at three resolutions independently, we need to:
1. **Merge overlapping loops** across resolutions (avoid redundancy)
2. **Annotate loops with nearest genes** for biological interpretation
3. **Classify loop anchors** using ChIP-seq data (promoters vs. enhancers)
4. **Categorize loop types** (P-P, P-E, E-E, etc.)
5. **Export in multiple formats** for downstream visualization and analysis

## Input Data

### Required Files

**From edgeR Analysis:**
```
outputs/edgeR_results_res_5kb/primary_analysis/
├── all_results_primary.tsv    # All tested loops
└── final_results.tsv           # Stringent filters (|logFC|>0.3, FDR<0.03)

outputs/edgeR_results_res_10kb/primary_analysis/
├── all_results_primary.tsv
└── final_results.tsv

outputs/edgeR_results_res_25kb/primary_analysis/
├── all_results_primary.tsv
└── final_results.tsv
```

**Genomic Coordinates:**
```
outputs/res_5kb/03_binned.rds   # GInteractions with loop positions
outputs/res_10kb/03_binned.rds
outputs/res_25kb/03_binned.rds
```

**Genome Annotations:**
- **TxDb.Mmusculus.UCSC.mm10.knownGene** - Gene models
- **org.Mm.eg.db** - Gene symbol mappings

**ChIP-seq Peaks (from `peaks/beds/`, timepoint-specific):**
- `H3K27ac{Tissue}{Timepoint}{Rep}.bed` - H3K27ac (active enhancers/promoters)
- `H3K27me3{Tissue}{Timepoint}{Rep}.bed` - H3K27me3 (Polycomb repression)
- `H3K4me1{Tissue}{Timepoint}{Rep}.bed` - H3K4me1 (poised enhancers)
- `H3K4me3{Tissue}{Timepoint}{Rep}.bed` - H3K4me3 (active promoters)
- `Bivalent_{Tissue}_{Timepoint}.bed` - K4me3+K27me3 overlap (developmental poised)

### Data Structure

Each `final_results.tsv` contains:
```
loop_id    chr1  start1  end1  chr2  start2  end2  logFC  FDR  PValue  ...
loop_1     chr1  1000000 1005000 chr1 1500000 1505000 0.45 0.01 ...
loop_2     chr2  2000000 2005000 chr2 2800000 2805000 -0.52 0.02 ...
```

## Processing Workflow

### Section 1: Load Final Results from All Resolutions

```r
# Load stringent differential loops (|logFC| > 0.3, FDR < 0.03)
for (res in c(5000, 10000, 25000)) {
  final_file <- sprintf("outputs/edgeR_results_res_%dkb/primary_analysis/final_results.tsv", res/1000)
  coords_file <- sprintf("outputs/res_%dkb/03_binned.rds", res/1000)

  # Load TSV and GInteractions
  final_df <- read.table(final_file, sep="\t", header=TRUE)
  coords <- readRDS(coords_file)

  # Match coordinates to differential loops
  # Store filtered coordinates and data
}
```

**Output:** Lists of differential loops and coordinates per resolution

### Section 1B: Merge All Results (Complete Dataset)

In addition to high-confidence loops, also merge **all tested loops** (including non-significant) for:
- Creating complete volcano plots
- Background comparisons
- Validation of merging strategy

```r
# Load all_results_primary.tsv from all resolutions
# Match to coordinates
# Merge with overlap removal (10kb tolerance)
```

**Output:** `merged_all_results.tsv` - Complete dataset for visualizations

### Section 2: Non-Redundant Merging Across Resolutions

**Challenge:** Same biological loop may be called at multiple resolutions with slightly different coordinates.

**Solution:** Merge loops within **10kb tolerance** using custom matching algorithm.

#### Matching Algorithm

```r
match_loops <- function(coords1, coords2, tolerance_bp = 10000) {
  # For each loop in coords1:
  #   1. Find loops in coords2 on same chromosomes
  #   2. Calculate distance for both anchors
  #   3. If both anchors within tolerance → match
  #   4. Take closest match if multiple candidates

  # Returns: data.frame(index1, index2, distance)
}
```

#### Merging Strategy

```r
# Priority: Keep loop from finest resolution where detected
# 1. Start with 5kb loops
# 2. Add 10kb loops not matching 5kb (within 10kb tolerance)
# 3. Add 25kb loops not matching 5kb or 10kb

merged_coords <- c(
  coords_5kb,                          # All 5kb loops
  coords_10kb[!matched_to_5kb],        # Unique 10kb loops
  coords_25kb[!matched_to_10kb_or_5kb] # Unique 25kb loops
)
```

**Rationale:**
- **Finer resolution preferred** for precise localization
- **Union approach** captures resolution-specific features
- **Overlap removal** prevents double-counting

**Output:**
- `non_redundant_loops.tsv` - Merged loop set with metadata
- `overlap_report.tsv` - Statistics on resolution overlaps

### Section 3: Nearest Gene Annotation

For each loop anchor, find the **nearest transcription start site (TSS)**.

```r
# Extract TSS from TxDb
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
genes <- genes(txdb)
tss <- resize(genes, width=1, fix="start")

# For each anchor
anchor1 <- anchors(merged_coords, type="first")
anchor2 <- anchors(merged_coords, type="second")

# Find nearest gene
nearest1 <- distanceToNearest(anchor1, tss)
nearest2 <- distanceToNearest(anchor2, tss)

# Map Entrez ID to gene symbol
gene_symbols <- mapIds(org.Mm.eg.db,
                       keys=mcols(tss)$gene_id,
                       column="SYMBOL",
                       keytype="ENTREZID")
```

**Output columns added:**
- `anchor1_nearest_gene` - Gene symbol
- `anchor1_gene_distance` - Distance to TSS (bp)
- `anchor1_gene_id` - Entrez gene ID
- `anchor2_nearest_gene` (same for second anchor)
- `anchor2_gene_distance`
- `anchor2_gene_id`

**Interpretation:**
- **Distance ≤ 2kb:** Likely promoter interaction
- **Distance 2-50kb:** Proximal regulatory element
- **Distance > 50kb:** Distal enhancer or other regulatory element

### Section 4: ChIP-seq Based Anchor Classification

**Traditional approach (distance-only):**
- ≤ 2kb from TSS → "Promoter"
- \> 2kb from TSS → "Enhancer"

**Problem:** Misses true enhancers near genes and misclassifies non-functional regions

**Our approach (evidence-based):**

Uses actual ChIP-seq data to classify anchors:

#### Classification Logic

```r
classify_anchor <- function(anchor, tss_distance, h3k27ac_overlap, h3k4me1_overlap) {

  if (h3k27ac_overlap & tss_distance <= 2000) {
    return("Promoter")           # Active TSS-proximal
  }

  else if (h3k27ac_overlap & tss_distance > 2000) {
    return("Active_Enhancer")    # Active distal element
  }

  else if (h3k4me1_overlap & !h3k27ac_overlap & tss_distance > 2000) {
    return("Poised_Enhancer")    # Primed but not active
  }

  else {
    return("Other")              # No clear marks
  }
}
```

#### Chromatin State Definitions

| Classification | ChIP-seq Marks | Distance to TSS | Biological Function |
|----------------|---------------|-----------------|---------------------|
| **Promoter** | H3K27ac+ | ≤ 2kb | Active gene transcription |
| **Active_Enhancer** | H3K27ac+ | > 2kb | Active distal regulation |
| **Poised_Enhancer** | H3K4me1+ (no H3K27ac) | > 2kb | Primed for activation |
| **Other** | No marks | Any | Structural, silenced, or unknown |

**Biological significance:**
- **H3K27ac** - Mark of active regulatory elements
- **H3K4me1** - Mark of enhancer elements (active + poised)
- **Poised enhancers** - Important for developmental regulation
- **Distance thresholds** - Account for typical promoter width

**Output columns added:**
- `anchor1_classification` - Functional category
- `anchor2_classification`
- `anchor1_H3K27ac` - Boolean
- `anchor1_H3K4me1` - Boolean
- `anchor2_H3K27ac`
- `anchor2_H3K4me1`

### Section 5: Loop Type Classification

Combine anchor classifications to categorize **interaction types**.

#### 10 Loop Type Categories

```
1. Promoter-Promoter (P-P)
2. Promoter-Active_Enhancer (P-AE)
3. Promoter-Poised_Enhancer (P-PE)
4. Promoter-Other (P-O)
5. Active_Enhancer-Active_Enhancer (AE-AE)
6. Active_Enhancer-Poised_Enhancer (AE-PE)
7. Active_Enhancer-Other (AE-O)
8. Poised_Enhancer-Poised_Enhancer (PE-PE)
9. Poised_Enhancer-Other (PE-O)
10. Other-Other (O-O)
```

**Output column:** `loop_type`

#### Biological Interpretation by Loop Type

**Promoter-Promoter (P-P):**
- Gene regulatory hubs
- Co-regulation of nearby genes
- Often at TAD boundaries
- Example: Hox gene clusters

**Promoter-Active_Enhancer (P-AE):**
- Classic enhancer-promoter looping
- Direct transcriptional activation
- Most common type in differentiated cells
- Tissue-specific regulation

**Promoter-Poised_Enhancer (P-PE):**
- Developmental regulation
- Primed for rapid activation
- Important in stem cells
- Example: Lineage-specific genes

**Active_Enhancer-Active_Enhancer (AE-AE):**
- Super-enhancer formation
- Phase-separated condensates
- Strong transcriptional activation
- Often disease-associated

**Poised_Enhancer-Poised_Enhancer (PE-PE):**
- Regulatory potential
- Silent in current cell state
- May activate during differentiation

**Other types:**
- Structural roles
- Silencing mechanisms
- Technical artifacts
- Require validation

### Section 6: Export Results

**Format 1: TSV (Tab-Separated Values)**
```
outputs/merged_loops/characterized_loops.tsv
```

Full table with all annotations:
```
loop_id, chr1, start1, end1, chr2, start2, end2,
logFC, FDR, PValue,
source_resolution, source_resolutions_overlap,
anchor1_nearest_gene, anchor1_gene_distance,
anchor2_nearest_gene, anchor2_gene_distance,
anchor1_classification, anchor2_classification,
loop_type,
anchor1_H3K27ac, anchor1_H3K4me1,
anchor2_H3K27ac, anchor2_H3K4me1,
distance_bp, intrachromosomal
```

**Format 2: RDS (R Data Structure)**
```
outputs/merged_loops/non_redundant_loops.rds
```

GInteractions object with all metadata as mcols

**Format 3: BEDPE (for visualization)**
```
outputs/merged_loops/non_redundant_loops.bedpe
```

Standard BEDPE format for Juicebox and genome browsers:
```
chr1  start1  end1  chr2  start2  end2  name  score  strand1  strand2
```

Where:
- `name` = loop_id
- `score` = 1000 * abs(logFC) (for coloring)
- Color by direction: up (red), down (blue)

### Section 7: Summary Statistics & Plots

**Statistics Generated:**

1. **Resolution Contribution**
```
Loops from 5kb only:  1,234 (45%)
Loops from 10kb only: 567 (21%)
Loops from 25kb only: 234 (9%)
Shared 5kb + 10kb:    345 (12%)
Shared 5kb + 25kb:    123 (4%)
Shared 10kb + 25kb:   89 (3%)
All three resolutions: 156 (6%)
```

2. **Loop Type Distribution**
```
P-P:    234 loops (12%)
P-AE:   678 loops (34%)
P-PE:   123 loops (6%)
AE-AE:  456 loops (23%)
...
```

3. **Distance Statistics**
```
Median loop distance: 250 kb
Range: 10 kb - 5 Mb
Interchromosomal: 23 loops (1%)
```

**Plots Generated:** `outputs/merged_loops/plots/`

1. **Venn Diagram** - Resolution overlap
2. **Bar Plot** - Loop type distribution
3. **Histogram** - Distance distribution
4. **Scatter** - logFC vs. distance

## Output Files

### Primary Output

**`characterized_loops.tsv`** - **Use this for downstream analysis**

Complete table with all annotations. Each row is one differential loop.

### Supporting Files

| File | Description | Format |
|------|-------------|---------|
| `non_redundant_loops.tsv` | Minimal version | TSV |
| `non_redundant_loops.rds` | GInteractions object | RDS |
| `non_redundant_loops.bedpe` | Genome browser format | BEDPE |
| `merged_all_results.tsv` | All loops (for volcano) | TSV |
| `overlap_report.tsv` | Resolution overlap stats | TSV |
| `plots/*.pdf` | Summary visualizations | PDF |

## Usage Examples

### Run with Default Settings

```bash
Rscript scripts/downstream_analysis.R
```

Uses:
- TxDb.Mmusculus.UCSC.mm10.knownGene for genes
- 10kb tolerance for overlap detection
- ChIP-seq BED files (if available)

### Run with Custom GTF

```bash
Rscript scripts/downstream_analysis.R --gtf /path/to/genes.gtf
```

### Run with Different Tolerance

```bash
Rscript scripts/downstream_analysis.R --tolerance 25
```

Uses 25kb tolerance (more lenient merging)

## Interpretation Guide

### High-Priority Differential Loops

Look for loops that are:
1. **P-AE or P-PE type** (enhancer-promoter interactions)
2. **|logFC| > 1.0** (strong effect)
3. **FDR < 0.01** (very significant)
4. **Nearest gene is known regulator** (biological relevance)

### Validation Candidates

Prioritize:
- Loops near genes of interest
- Loops with extreme fold changes
- P-P loops at known regulatory loci
- Loops validated across multiple resolutions

### Red Flags

Be cautious of:
- **O-O loops with extreme logFC** - May be technical artifacts
- **Very short loops (< 50kb)** - Check if self-ligation
- **Interchromosomal loops** - Validate carefully
- **Loops with no ChIP-seq marks** - May be structural only

## Biological Insights

### Differential Loop Patterns

**Upregulated in BAP1-KO (logFC > 0):**
- Potential compensatory mechanisms
- Loss of loop repression
- Chromatin reorganization

**Downregulated in BAP1-KO (logFC < 0):**
- Direct BAP1-dependent loops
- Loss of regulatory interactions
- Gene silencing mechanisms

### Loop Type Enrichment

If differential loops are **enriched for P-AE:**
→ BAP1 regulates enhancer-promoter communication

If differential loops are **enriched for AE-AE:**
→ BAP1 affects super-enhancer organization

If differential loops are **enriched for P-PE:**
→ BAP1 affects developmental gene priming

## Performance

**Expected Runtime:** 10-15 minutes

**Memory Usage:** 4-8 GB RAM

**Disk Space:** ~100 MB output

## Dependencies

**R Packages:**
```r
library(mariner)
library(InteractionSet)
library(GenomicRanges)
library(rtracklayer)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(tidyverse)
library(ggplot2)
```

**External Data:**
- ChIP-seq BED files (optional)
- Gene annotations (automatic from Bioconductor)

## Troubleshooting

### ChIP-seq Files Not Found

**Error:** `H3K27ac bed file not found`

**Solution:**
- Script will run without ChIP-seq, using distance-only classification
- Anchor types will be limited to "Promoter_proximal", "Distal", "Other"
- For full classification, provide BED files in project root

### No Overlaps Between Resolutions

**Possible Causes:**
- Different filtering was applied
- Coordinate systems don't match
- Tolerance too stringent

**Solutions:**
- Verify all resolutions completed edgeR
- Check that final_results.tsv exists for all
- Increase tolerance (--tolerance 25)

### Gene Annotation Fails

**Error:** Cannot load TxDb

**Solution:**
```r
BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
BiocManager::install("org.Mm.eg.db")
```

## Quality Control

### Check Output Sizes

Expected loop counts:
- **5kb final_results:** 500-3,000 loops
- **10kb final_results:** 400-2,500 loops
- **25kb final_results:** 300-2,000 loops
- **Merged non-redundant:** 800-4,000 loops

If very different, investigate:
- Stringent filters too harsh/lenient
- Resolution bias in loop calling
- Poor merging efficiency

### Validate Merging

Check `overlap_report.tsv`:
- **High overlap (>80%):** Resolutions capturing same biology
- **Low overlap (<50%):** Resolution-specific effects (expected)
- **Zero overlap:** Problem with coordinate matching

### Verify Annotations

Sample random loops and check:
- Nearest genes make biological sense
- ChIP-seq classifications are reasonable
- Loop types match genomic context

## Next Steps

After running downstream_analysis.R:

1. **Visualize results:**
   ```bash
   Rscript scripts/visualizations.R
   ```

2. **Validate specific loops:**
   - Load BEDPE into Juicebox
   - Check Hi-C signal visually
   - Compare to published ChIA-PET data

3. **Functional analysis:**
   - GO enrichment on nearest genes
   - KEGG pathway analysis
   - Disease association studies

4. **Experimental validation:**
   - 3C/4C at selected loci
   - ChIP-qPCR for factor binding
   - Gene expression changes

---

**Last Updated:** November 10, 2025
