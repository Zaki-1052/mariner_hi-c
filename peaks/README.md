# Extended Loop Anchor Annotation with Chromatin State Categories

## Script: `annotate_loops_extended.R`

## Overview

This script annotates differential chromatin loop anchors with **8 chromatin state categories** based on histone modification ChIP-seq data. It extends the basic promoter/enhancer classification to include Polycomb-repressed regions, bivalent promoters, and CTCF-bound structural sites, providing a more nuanced view of loop anchor function.

## Purpose

After identifying differential loops between BAP1-KO and wildtype, we need to understand the **functional context** of these changes:

1. **Classify loop anchors** by chromatin state (not just distance to TSS)
2. **Distinguish active vs. repressed** regulatory elements
3. **Identify bivalent promoters** (developmental poised domains)
4. **Compare timepoints** (Early vs. Late cerebellum development)
5. **Compare peak sources** (single-replicate vs. consensus H3K4me3)

## Evolution & Development History

The script evolved through several iterations to address biological and technical challenges:

| Commit | Change | Rationale |
|--------|--------|-----------|
| `d903a1d` | Initial annotation script | Basic 6-category system with Polycomb |
| `4bdae62` | Loop annotations | Added loop type classification |
| `2442094` | Bivalent promoter | Split to 7 categories with H3K4me3-based Active_Promoter |
| `e335242` | Peak files - new | Standardized peak file naming (Cerebellum Early/Late) |
| `c5b3e92` | New peak results | Updated to use peaks/beds/ directory |
| `49f69b1` | New bivalent bed files | Created traceable bivalent generation script |
| `c235331` | Late consensus | Added late_consensus timepoint with 4-replicate H3K4me3 |

### Key Design Decisions

1. **H3K4me3 for Active_Promoter** (not H3K27ac): H3K4me3 is the canonical active promoter mark. Housekeeping genes can be H3K4me3+ without H3K27ac.

2. **Separate Repressed_Promoter category**: H3K27me3+ near TSS indicates Polycomb silencing at promoters, distinct from distal Polycomb.

3. **Pre-computed bivalent regions**: Using bedtools intersect of H3K4me3 + H3K27me3 ensures consistent bivalent calling.

4. **Traceable peak generation**: All bivalent files generated via `scripts/generate_bivalent_peaks.sb` for reproducibility.

## Chromatin State Categories

### 8 Anchor Types (Priority Order)

| Category | Definition | Biological Meaning |
|----------|------------|-------------------|
| **Active_Promoter** | H3K4me3+ AND NOT H3K27me3 AND ≤2kb from TSS | Actively transcribed gene promoter |
| **Repressed_Promoter** | H3K27me3+ AND NOT H3K27ac AND ≤2kb from TSS | Polycomb-silenced promoter |
| **Bivalent_Promoter** | H3K4me3+H3K27me3 overlap (pre-computed) | Developmental poised domain |
| **Polycomb** | H3K27me3+ AND >2kb from TSS | Distal repressive element |
| **Active_Enhancer** | H3K27ac+ AND >2kb from TSS | Active distal regulatory element |
| **Poised_Enhancer** | H3K4me1+ AND NOT H3K27ac AND NOT H3K27me3 AND >2kb | Primed but not active enhancer |
| **CTCF_Site** | CTCF ChIP+ OR CTCF motif+ AND not classified above | Structural/insulator loop anchor |
| **Other** | No ChIP-seq marks AND no CTCF evidence | Truly unmarked regions |

### Classification Logic

```
Priority 1: Active_Promoter
  └── H3K4me3+ AND NOT H3K27me3 AND distance_to_TSS ≤ 2kb

Priority 2: Repressed_Promoter
  └── H3K27me3+ AND NOT H3K27ac AND distance_to_TSS ≤ 2kb

Priority 3: Bivalent_Promoter
  └── Overlaps pre-computed K4me3+K27me3 intersection

Priority 4: Polycomb
  └── H3K27me3+ AND distance_to_TSS > 2kb

Priority 5: Active_Enhancer
  └── H3K27ac+ AND distance_to_TSS > 2kb

Priority 6: Poised_Enhancer
  └── H3K4me1+ AND NOT H3K27ac AND NOT H3K27me3 AND distance_to_TSS > 2kb

Priority 7: CTCF_Site
  └── (CTCF ChIP+ OR CTCF motif+) AND not classified above (structural/insulator)

Priority 8: Other
  └── Default (no marks, no CTCF evidence)
```

### 36 Loop Type Combinations

Loop types are named by combining anchor types in hierarchy order:
- `Active_Promoter-Active_Promoter`
- `Active_Promoter-Repressed_Promoter`
- `Active_Promoter-Bivalent_Promoter`
- `Active_Promoter-Polycomb`
- `Active_Promoter-CTCF_Site`
- ... etc.

## Peak File Configuration

### Directory Structure

```
peaks/
├── CTCF.bed                           # CTCF ChIP-seq peaks (32,487) - Bing Ren lab
├── ctcf_motifs_mm10.bed               # CTCF DNA motifs (114,081) - HOMER pre-computed
├── README.md                          # This file
├── beds/                              # Standardized ChIP-seq peak files
│   ├── H3K27acCerebellumEarly2.bed
│   ├── H3K27acCerebellumLate2.bed
│   ├── H3K27me3CerebellumEarly1.bed
│   ├── H3K27me3CerebellumLate1.bed
│   ├── H3K4me1CerebellumEarly1.bed
│   ├── H3K4me1CerebellumLate1.bed
│   ├── H3K4me3CerebellumEarly2.bed
│   ├── H3K4me3CerebellumLate2.bed    # 6,581 peaks
│   ├── Bivalent_Cerebellum_Early.bed  # 931 peaks
│   ├── Bivalent_Cerebellum_Late.bed   # 318 peaks
│   └── Bivalent_Consensus_Late.bed    # 688 peaks
├── scripts/                           # All analysis scripts
│   ├── annotate_loops_extended.R      # Main 8-category annotation
│   ├── extract_other_anchors.R        # Extract "Other" anchors to BED
│   ├── extract_other_up_loops.R       # Analyze "Other" up loops
│   ├── annotate_ctcf_motifs.R         # Integrate CTCF motif results
│   └── ctcf_motif_analysis.sb         # SLURM: Full CTCF motif pipeline
├── loop_annotation_extended/          # Output directories
│   ├── early/                         # Early timepoint results
│   ├── late/                          # Late timepoint results
│   └── extra/                         # Comparison analyses
│       ├── early_p12ctrl/
│       └── late_consensus/
└── old/peaks-v1/                      # Legacy peak files
    ├── consensus_H3K4me3_late_peaks.bed      # 9,651 peaks (4-replicate)
    ├── P12_ctrl_H3K27ac_early_peaks.bed      # 28,042 peaks (lenient)
    ├── P12_ctrl_H3K27me3_early_peaks.bed     # 23,491 peaks (lenient)
    └── 250224AddisonH3K4me3H3K27me3Early.bed # 933 peaks (bivalent)
```

### Peak File Naming Convention

Files follow the pattern: `{Mark}{Tissue}{Timepoint}{Replicate}.bed`

- **Mark**: H3K27ac, H3K27me3, H3K4me1, H3K4me3
- **Tissue**: Cerebellum
- **Timepoint**: Early, Late
- **Replicate**: 1, 2 (if multiple)

### Bivalent File Generation

Bivalent files are created using `scripts/generate_bivalent_peaks.sb`:

```bash
# Uses bedtools intersect -u (unique H3K4me3 peaks overlapping H3K27me3)
bedtools intersect -a H3K4me3.bed -b H3K27me3.bed -u > Bivalent.bed
```

The `-u` flag ensures each H3K4me3 peak is counted once, even if it overlaps multiple H3K27me3 peaks.

### Adding New Peak Files

To add new peak files for future timepoints or tissues:

1. **Add standardized BED files** to `peaks/beds/` following naming convention
2. **Generate bivalent file** using `generate_bivalent_peaks.sb` (add new section)
3. **Update PEAK_FILES** in `annotate_loops_extended.R`:

```r
PEAK_FILES <- list(
  ...
  new_timepoint = list(
    h3k27ac    = "beds/H3K27acNew.bed",
    h3k27me3   = "beds/H3K27me3New.bed",
    h3k4me1    = "beds/H3K4me1New.bed",
    h3k4me3    = "beds/H3K4me3New.bed",
    bivalent   = "beds/Bivalent_New.bed",
    ctcf       = "CTCF.bed",
    ctcf_motif = "ctcf_motifs_mm10.bed"
  )
)
```

4. **Update DEFAULT_INPUT_FILES** and **DEFAULT_OUTPUT_DIRS**
5. **Update timepoint validation** in `annotate_loops_extended()` function

## Usage

### Basic Usage

```bash
# Run from the peaks/ directory
cd peaks/

# Run both early and late timepoints (default)
Rscript scripts/annotate_loops_extended.R

# Run specific timepoint
Rscript scripts/annotate_loops_extended.R --timepoint early
Rscript scripts/annotate_loops_extended.R --timepoint late
Rscript scripts/annotate_loops_extended.R --timepoint late_consensus
Rscript scripts/annotate_loops_extended.R --timepoint early_p12ctrl

# Run all 4 timepoints for full comparison
Rscript scripts/annotate_loops_extended.R --all
```

### Command-Line Options

| Option | Description |
|--------|-------------|
| `--timepoint TP` | Run specific timepoint: `early`, `late`, `late_consensus`, `early_p12ctrl`, or `both` |
| `--all` | Run all 4 timepoints for full comparison |
| `--input FILE` | Override input file (single timepoint only) |
| `--output DIR` | Override output directory (single timepoint only) |
| `--help`, `-h` | Show help message |

### Timepoint Descriptions

| Timepoint | H3K4me3 Source | Bivalent Peaks | Purpose |
|-----------|---------------|----------------|---------|
| `early` | Cerebellum Early (10,370) | 931 | Early developmental stage |
| `late` | Cerebellum Late (6,581) | 318 | Late/adult stage |
| `late_consensus` | Consensus (9,651) | 688 | Compare with 4-replicate H3K4me3 |
| `early_p12ctrl` | Cerebellum Early (10,370) | 933 (Addison) | Compare with P12_ctrl peaks |

**Note:** `early_p12ctrl` uses P12_ctrl H3K27ac (28,042) and H3K27me3 (23,491) instead of Cerebellum (18,178 and 12,473). This demonstrates the ~54-88% more peaks from lenient peak calling.

## Input Data

### Required Files

**Differential Loops:**
```
{timepoint}_outputs/merged_loops/non_redundant_loops.tsv
```

Required columns:
- `chr1`, `start1`, `end1` - Anchor 1 coordinates
- `chr2`, `start2`, `end2` - Anchor 2 coordinates
- `logFC` - Log fold change
- `direction` - `up_in_mutant` or `down_in_mutant`

**ChIP-seq Peaks:**
- H3K27ac, H3K27me3, H3K4me1, H3K4me3, Bivalent BED files

**Gene Annotations:**
- `TxDb.Mmusculus.UCSC.mm10.knownGene` (loaded automatically)

## Output Files

### Per-Timepoint Output (`outputs/loop_annotation_extended/{timepoint}/`)

| File | Description |
|------|-------------|
| `extended_characterized_loops.tsv` | Full annotation table with all columns |
| `anchor_type_summary.tsv` | Anchor type counts by anchor |
| `loop_type_summary.tsv` | Loop type counts |
| `summary_statistics.txt` | Human-readable summary |
| `plots/` | Visualization PDFs |

### Visualization Plots

1. **`loop_type_piechart_comparison.pdf`** - Side-by-side pie charts for up vs. down loops
2. **`anchor_type_distribution.pdf`** - Stacked bar chart of anchor types by direction
3. **`loop_type_by_direction.pdf`** - Horizontal bar chart of loop types

### Output Columns

```
Original loop columns:
  chr1, start1, end1, chr2, start2, end2, logFC, FDR, direction, ...

ChIP-seq overlap columns:
  anchor1_H3K27ac_overlap, anchor1_H3K27me3_overlap, anchor1_H3K4me1_overlap,
  anchor1_H3K4me3_overlap, anchor1_Bivalent_Promoter_overlap, anchor1_CTCF_overlap,
  anchor1_CTCF_motif_overlap (same for anchor2)

TSS distance columns:
  anchor1_distance_to_tss_ext, anchor2_distance_to_tss_ext

Classification columns:
  anchor1_type_extended, anchor2_type_extended, loop_type_extended
```

## Comparing Late vs. Late_Consensus

The `late_consensus` timepoint allows comparison of annotation results using:
- **Cerebellum H3K4me3** (6,581 peaks) vs. **Consensus H3K4me3** (9,651 peaks)
- **318 bivalent peaks** vs. **688 bivalent peaks**

### Example Comparison (Late vs. Late_Consensus)

| Category | Late (Cerebellum) | Late_Consensus | Difference |
|----------|-------------------|----------------|------------|
| Active_Promoter (Anchor1) | 253 (8.7%) | 291 (10.0%) | +38 |
| Bivalent_Promoter (Anchor1) | 25 (0.9%) | 32 (1.1%) | +7 |
| Other (Anchor1) | 967 (33.2%) | 926 (31.8%) | -41 |

The consensus H3K4me3 captures more peaks, leading to:
- More Active_Promoter classifications (H3K4me3+ near TSS)
- More Bivalent_Promoter classifications (from 688 vs. 318 bivalent peaks)
- Fewer "Other" classifications

## Comparing Early vs. Early_P12ctrl (Peak Source Comparison)

The `early_p12ctrl` timepoint demonstrates the lack of a single source of truth for ChIP-seq annotations by comparing:
- **Cerebellum H3K27ac** (18,178 peaks) vs. **P12_ctrl H3K27ac** (28,042 peaks) - 54% more
- **Cerebellum H3K27me3** (12,473 peaks) vs. **P12_ctrl H3K27me3** (23,491 peaks) - 88% more

### Peak Count Comparison

| Mark | Cerebellum (early) | P12_ctrl (early_p12ctrl) | Difference |
|------|-------------------|-------------------------|------------|
| H3K27ac | 18,178 | 28,042 | +54% more |
| H3K27me3 | 12,473 | 23,491 | +88% more |
| H3K4me1 | 93,859 | 93,859 (same) | - |
| H3K4me3 | 10,370 | 10,370 (same) | - |
| Bivalent | 931 (generated) | 933 (Addison) | ~same |

**Key finding:** Cerebellum files are a ~97% subset of P12_ctrl (stricter filtering applied).

### Example Comparison (Early vs. Early_P12ctrl)

| Category | Early (Cerebellum) | Early_P12ctrl | Change |
|----------|-------------------|---------------|--------|
| Bivalent_Promoter (Anchor1) | 2 (1.2%) | 4 (2.4%) | +100% |
| Bivalent_Promoter (Anchor2) | 1 (0.6%) | 4 (2.4%) | +300% |
| Polycomb (Anchor2) | 21 (12.7%) | 24 (14.5%) | +14% |

Despite 54-88% more peaks in P12_ctrl, loop annotations are relatively stable. The biggest changes are in Bivalent and Polycomb categories.

## Biological Interpretation

### Anchor Type Patterns

**High Active_Promoter + Repressed_Promoter:**
- Loops connect promoters (gene regulatory hubs)
- May indicate promoter-promoter interactions

**High Polycomb:**
- Loops involve distal Polycomb-repressed regions
- Long-range repressive chromatin architecture

**High Active_Enhancer + Poised_Enhancer:**
- Classic enhancer-promoter or enhancer-enhancer interactions
- Regulatory loops for gene activation

**High Other:**
- Structural loops (CTCF-mediated)
- Regions without canonical histone marks

### Loop Type Patterns by Direction

**Upregulated in BAP1-KO (logFC > 0):**
- If enriched for Polycomb loops → Loss of BAP1 affects Polycomb targeting
- If enriched for Active loops → Compensatory activation

**Downregulated in BAP1-KO (logFC < 0):**
- If enriched for Active_Enhancer loops → BAP1 maintains enhancer contacts
- If enriched for Promoter loops → BAP1-dependent gene regulation

## Performance

**Runtime:** 2-5 minutes per timepoint

**Memory:** 4-8 GB RAM

**Dependencies:**
- R ≥ 4.3.0
- GenomicRanges, InteractionSet, rtracklayer
- TxDb.Mmusculus.UCSC.mm10.knownGene
- tidyverse, ggplot2, patchwork, RColorBrewer

## Troubleshooting

### Peak File Not Found

```
Error: H3K27ac bed file not found: peaks/beds/H3K27acCerebellumLate2.bed
```

**Solution:** Verify peak files exist in `peaks/beds/`. Check file naming matches PEAK_FILES configuration.

### Bivalent File Has Wrong Peak Count

**Cause:** Bivalent file may have been created with different bedtools parameters.

**Solution:** Re-run `bash scripts/generate_bivalent_peaks.sb` to regenerate with consistent `-u` flag.

### No Loops Annotated as Bivalent

**Possible causes:**
1. Bivalent peaks don't overlap loop anchors
2. Bivalent file is empty or wrong format
3. Priority order: Active_Promoter or Repressed_Promoter may take precedence

**Check:** Verify bivalent peak count in summary output.

### Memory Error

**Solution:** Run one timepoint at a time instead of `--all`:
```bash
Rscript scripts/annotate_loops_extended.R --timepoint early
Rscript scripts/annotate_loops_extended.R --timepoint late
```

## Related Scripts

All scripts are located in `peaks/scripts/`:

| Script | Purpose |
|--------|---------|
| `annotate_loops_extended.R` | Main 8-category chromatin state annotation |
| `extract_other_up_loops.R` | Extract and analyze "Other" category up loops |
| `extract_other_anchors.R` | Extract "Other" anchors to BED for motif analysis |
| `annotate_ctcf_motifs.R` | Integrate CTCF motif counts with loop annotations |
| `ctcf_motif_analysis.sb` | SLURM: Full CTCF motif validation pipeline |

---

## Extracting "Other" Category Up Loops

### Script: `extract_other_up_loops.R`

This script extracts and annotates loops with "Other" classification (truly unmarked anchors) that are upregulated in BAP1-KO, with gene annotations.

### Background

With CTCF motif integration, the early timepoint now shows comparable "Other" percentages to late:

| Anchor Type | Early (ChIP + Motif) | Late (ChIP + Motif) |
|-------------|----------------------|---------------------|
| CTCF_Site   | ~30% | ~24% |
| Other       | **~6-7%** | **~9%** |

**Key change:** CTCF_Site classification now uses BOTH adult CTCF ChIP-seq AND CTCF DNA motif data. This captures potential CTCF sites even without timepoint-specific ChIP-seq.

### Usage

```bash
cd peaks/
Rscript scripts/extract_other_up_loops.R
```

### Key Findings (Early Timepoint)

- **16 out of 73 up loops (22%)** have at least one "Other" anchor
  - 2 loops with both anchors "Other"
  - 14 loops with single anchor "Other"
- Only **17 unique "Other" anchors** remain (vs 87 before motif integration)
- These are truly unmarked regions - no histone marks, no CTCF ChIP, no CTCF motif

### Top Significant "Other" Up Loops

| Loop ID | Coordinates | logFC | FDR | Genes |
|---------|-------------|-------|-----|-------|
| loop_17936 | chr8:71.19Mb ↔ 71.81Mb | +1.01 | 2.56e-05 | Haus8 - Zfp709 |
| loop_4918 | chr8:49.80Mb ↔ 57.30Mb | +1.76 | 7.59e-04 | Gm2516 - Hand2 |
| loop_6571 | chr4:32.70Mb ↔ 34.63Mb | +0.51 | 4.70e-03 | Mdn1 - Rars2 |
| loop_4307 | chr8:71.14Mb ↔ 71.83Mb | +0.58 | 5.35e-03 | 1700026F02Rik - Zfp709 |

### Output Files

Located in `loop_annotation_extended/early/`:

| File | Description |
|------|-------------|
| `other_category_up_loops.tsv` | Main table with 16 loops, ranked by FDR |
| `other_up_loops_with_genes.tsv` | Full annotation with all original columns |
| `other_up_loops_summary.txt` | Summary statistics |

### Output Columns

```
loop_id, anchor1_coords, anchor2_coords
logFC, FDR, category (strong/moderate/weak_differential)
anchor1_type, anchor2_type, loop_type
other_category (both_other or single_other)
anchor1_nearest_gene, anchor1_gene_distance
anchor2_nearest_gene, anchor2_gene_distance
anchor1_overlaps_CTCF, anchor2_overlaps_CTCF
loop_distance, resolution_kb
```

### Notes

- With CTCF motif integration, "Other" now truly means unmarked regions (no histone marks, no CTCF evidence)
- The remaining 16 "Other" up loops represent potential novel structural elements or regions with marks below detection threshold
- Most former "Other" anchors are now classified as CTCF_Site due to motif overlap

---

## CTCF Motif Validation Pipeline

### Purpose

Since the early timepoint lacks CTCF ChIP-seq data, ~36.5% of anchors are classified as "Other" (vs 9.3% in late). To validate whether these "Other" anchors are likely CTCF-mediated, we check for CTCF DNA binding motifs using HOMER's pre-computed genome-wide motif scan.

### Method

Instead of running HOMER on-the-fly, we use a pre-computed file containing all known motif locations in mm10:
- **Source:** `homer.KnownMotifs.mm10.191020.bed.gz` (4GB, ~131k CTCF motifs)
- **Filter:** Canonical `CTCF(Zf)` motif only (exclude Satellite variants)
- **Tool:** `bedtools intersect` to count motifs per anchor

This approach is:
- Faster than running HOMER
- Reproducible (same pre-computed data)
- Biologically equivalent (same PWM-based motif detection)

### Usage

**On HPC (single command):**
```bash
cd /expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/peaks
sbatch scripts/ctcf_motif_analysis.sb
```

This master script runs 4 steps:
1. Verify input files exist
2. Extract canonical CTCF(Zf) motifs from pre-computed file
3. Run bedtools intersect to count motifs per anchor
4. Run R integration script to generate concordance matrix

### Prerequisites

Before running the pipeline, generate the "Other" anchors BED file:
```bash
cd peaks/
Rscript scripts/extract_other_anchors.R
```

This creates `loop_annotation_extended/early/other_anchors.bed` (87 unique anchors).

### Output Files

Located in `loop_annotation_extended/early/`:

| File | Description |
|------|-------------|
| `anchors_ctcf_motif_count.bed` | bedtools output with motif counts per anchor |
| `ctcf_motif_validation.txt` | Concordance matrix and interpretation |
| `other_category_up_loops.tsv` | Updated with `has_CTCF_motif` columns |

### Concordance Matrix

The analysis compares two approaches for identifying CTCF sites:

|                     | Late CTCF ChIP+ | Late CTCF ChIP- |
|---------------------|-----------------|-----------------|
| **CTCF Motif+**     | Strong CTCF     | Motif only      |
| **CTCF Motif-**     | Non-canonical   | True "Other"    |

**Interpretation:**
- **Motif+ ChIP+**: Strong evidence for CTCF site
- **Motif+ ChIP-**: Motif present but no late ChIP binding (may be early-specific)
- **Motif- ChIP+**: ChIP binding without canonical motif (indirect recruitment)
- **Motif- ChIP-**: True "Other" - structural loop without CTCF

### Expected Results

Based on CTCF biology:
- **50-70%** of "Other" anchors should contain ≥1 CTCF motif
- **High concordance** between motif+ and late ChIP+ (~70-80% agreement)
- Remaining **true "Other"** anchors may represent cohesin-only loops or novel structural elements

---

## Quality Control

### Expected Peak Counts (Late Timepoint)

| Mark | Cerebellum | Consensus |
|------|------------|-----------|
| H3K27ac | ~15,000 | - |
| H3K27me3 | ~15,800 | - |
| H3K4me1 | ~113,000 | - |
| H3K4me3 | ~6,500 | ~9,600 |
| Bivalent | ~320 | ~690 |

### Sanity Checks

1. **Bivalent count:** Should be a subset of H3K4me3 peaks
2. **Active_Promoter < H3K4me3:** Not all H3K4me3 peaks are near TSS
3. **Other not too high:** >50% "Other" may indicate missing peak files
4. **Anchor1 ≈ Anchor2:** Similar distributions expected

---

**Last Updated:** January 13, 2026

**Author:** Zakir Alibhai
