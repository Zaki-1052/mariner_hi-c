# Mariner Differential Loop Analysis - Context Document for edgeR

## Project Overview

**Research Question**: Identifying differential chromatin loops between BAP1 mutant and control mouse samples using Hi-C data.

**Current Status**: Completed mariner extraction pipeline and QC validation. Ready for differential analysis with edgeR.

**Key Challenge**: Hi-C loop calls can shift positions by 1-2 bins (5-10kb) between samples/replicates due to spatial imprecision. Mariner's buffer approach addresses this by extracting 5×5 pixel matrices (±10kb) around each loop.

---

## Current Pipeline State - What's Been Completed

### 1. Data Preparation & Merging ✅
**Script**: `prepare_loops.R`

**Input files used**:
```bash
Control: /expanse/lustre/projects/csd940/ctea/nf-hic/hiccups/hiccups_results/resorted_ctrl/merged_loops.bedpe
Mutant: /expanse/lustre/projects/csd940/ctea/nf-hic/hiccups/hiccups_results/resorted_mut/merged_loops.bedpe
```

**What was done**:
- Read BEDPE files (skipping 2 header lines)
- Filtered for 5kb resolution only: `resolution = x2 - x1 = 5000`
- Removed loops with `observed = 0`
- Starting counts: 17,568 control loops + 17,550 mutant loops
- **Merged with mariner**: Used 10kb radius, selected by highest `observed` count
- **Result**: 22,108 consensus loops (unified positions across both conditions)

**Key transformation**: Created unified genomic positions so same loops can be compared between conditions despite potential positional shifts.

### 2. Binning & Buffering ✅

**Binning** (`assignToBins`):
- Snapped loops to exact 5kb grid positions
- Used center of each range as reference
- Output: All 22,108 loops aligned to 5kb boundaries

**Buffering** (`pixelsToMatrices`):
- Created 5×5 bin regions around each loop
- Buffer = 2 bins = ±10kb around center position
- **Purpose**: Capture signal even if loop is shifted by 1-2 bins between conditions
- Output: 22,108 regions, each covering 25kb × 25kb (5×5 bins)

### 3. Hi-C Data Extraction ✅
**Script**: `05_extract_counts.R`

**.hic files**:
```bash
Control: /expanse/lustre/projects/csd940/ctea/nf-hic/hiccups/hicfiles/resorted_ctrl_merged.hic (6.3 GB)
Mutant: /expanse/lustre/projects/csd940/ctea/nf-hic/hiccups/hicfiles/resorted_mut_merged.hic (5.9 GB)
```

**Extraction parameters**:
- Resolution: 5kb
- Normalization: VC (Vanilla Coverage) - KR not available
- Matrix type: observed
- Storage: HDF5 format (on-disk for memory efficiency)

**Output dimensions**: 5 × 5 × 22,108 × 2
- 5×5: spatial matrix for each loop
- 22,108: number of loops
- 2: samples (ctrl, mut)

**Extraction time**: 32.4 seconds

### 4. Count Aggregation ✅

**Aggregation method**: Sum of all 25 pixels in each 5×5 matrix

**Alternative method available**: Weighted (center pixels weighted higher) - correlation with sum > 0.99

**Final count matrix**: 22,108 loops × 2 samples
- Control total counts: 10,789,577
- Mutant total counts: 10,377,560
- Library size ratio: 0.962 (very comparable)

### 5. Quality Control Validation ✅
**Script**: QC analysis (see outputs)

**QC Results Summary**:

✓ **Data Integrity**: 
- No NA values, no negatives, no zeros
- All extractions successful

✓ **Sample Quality**:
- Pearson correlation: 0.9015 (excellent technical quality)
- Spearman correlation: 0.9790
- Library sizes comparable (4% difference)

✓ **Biological Signal**:
- 54% of loops show positional shifts between conditions
  - **Meaning**: In 54% of loops, the peak signal in the 5×5 matrix is at different positions between ctrl and mut
  - **Per-loop shift status**: Saved as boolean vector for downstream enrichment analysis
  - **Validates buffer approach**: Without buffer, would miss shifted loops
- Mean center enrichment: 0.074 (ctrl), 0.075 (mut)
  - **Meaning**: Only 7-8% of total signal in the center pixel
  - **This is normal**: Signal naturally diffuse at 5kb resolution
  - **Why it matters**: Confirms need for buffer - 92% of signal is in surrounding pixels

⚠ **Warnings (expected and normal)**:
- Less than 50% of loops have perfectly centered peaks
- Low center enrichment (<0.2)
- Both expected for 5kb Hi-C resolution - signal is biologically diffuse

**Pre-differential fold changes** (raw, without statistics):
- Mean log2 FC: -0.066 ± 0.220
- |M| > 1 (>2-fold): 46 loops (0.2%)
- |M| > 2 (>4-fold): 2 loops (0.0%)

**Key takeaway**: Most loops show subtle changes - need proper statistics (edgeR) to identify significant differences.

---

## Output Files - Current Locations

**Base directory**: `/expanse/lustre/projects/csd940/zalibhai/mariner/outputs/full`

```
01_ginteractions.rds        # Initial loops (ctrl + mut separate)
02_merged.rds               # Merged consensus (22,108 loops)
03_binned.rds               # Loops snapped to 5kb grid
04_buffered.rds             # 5×5 buffered regions
05_extractedse.rds          # Full InteractionArray with 5×5 matrices
05_extractedassays.h5       # HDF5 storage for matrices
05_metadata.rds             # Metadata for extracted data
06_counts_matrix.rds        # **THIS IS WHAT YOU NEED FOR edgeR**
06_counts_matrix.tsv        # Same data, tab-separated format
06_edgeR_input.rds          # Pre-formatted for edgeR (if available)
06_all_strategies.rds       # Different aggregation strategies
```

**QC outputs**: `/expanse/lustre/projects/csd940/zalibhai/mariner/outputs/qc_report`
- 9 PDF plots validating data quality
- `qc_report_summary.rds` with all metrics
- `loop_shift_status.rds` - Boolean vector (length 22,108) of shift status per loop - `loop_shift_summary.tsv` - Human-readable table with peak positions and shift detection

---

## Understanding Key Technical Concepts

### The Buffer Approach - Why It's Necessary

**The Problem**: 
Hi-C loop calling is spatially imprecise. The same biological loop can be called at slightly different positions (±5-10kb) between:
- Different samples
- Different replicates
- Different conditions

**Why it happens**:
1. Biological: Chromatin loops aren't point contacts - they have spatial extent
2. Technical: 5kb bins average signal across large genomic windows
3. Statistical: Background noise spreads across neighboring bins
4. Imprecision: Loop callers can be off by 1-2 bins

**The Solution - 5×5 Buffer**:
```
Without buffer (center pixel only):
Control loop at position X: signal = 15
Mutant same loop shifted to X+5kb: signal at X = 8
→ Looks like it decreased (WRONG - just shifted!)

With buffer (5×5 sum):
Control: sum of 25 pixels = 200
Mutant: sum of 25 pixels = 198 (captures the shifted peak)
→ Correctly shows they're similar
```

### What "54% Show Positional Shifts" Means

In 54% of the 22,108 loops, the maximum signal in the 5×5 matrix falls at different pixel positions between control and mutant.

Example:
```
Control 5×5:          Mutant 5×5:
[ 2  3  4  3  2 ]    [ 3  5  7  4  2 ]
[ 3  5  8  5  3 ]    [ 5  9 12  8  4 ]
[ 4  8 [15] 8  4 ]    [ 7 12 18 11  5 ]  ← Peak shifted right
[ 3  5  8  5  3 ]    [ 4  8 11  7  3 ]
[ 2  3  4  3  2 ]    [ 2  4  5  3  2 ]

Peak at [3,3]        Peak at [3,4]
```

**This validates the buffer approach** - without it, you'd misinterpret shifted loops as decreased signal.

### What "Low Center Enrichment" Means

Center enrichment = 0.074 means only 7.4% of the total 5×5 signal is in the center pixel.

**This is normal and expected**:
- At 5kb resolution, signal naturally spreads across neighboring bins
- If you only used center pixel, you'd lose 92% of the signal
- Confirms the buffer is necessary to capture full loop signal

---

## Hiccups Differential Results - What You're Comparing Against

### Location & Structure

**Directory**: `/expanse/lustre/projects/csd940/ctea/nf-hic/juicerpre/merged/diffhiccups_results/resorted_ctrl_vs_resorted_mut/`

```
├── differential_loops1.bedpe  # Summary file
├── differential_loops2.bedpe  # Summary file
├── file1/                     # Loops enriched in CONTROL
│   └── postprocessed_pixels_5000.bedpe  (19,386 loops)
└── file2/                     # Loops enriched in MUTANT
    └── postprocessed_pixels_5000.bedpe  (19,207 loops)
```

**Confirmed identity**:
- file1 = control-enriched loops (total observed: 1,588,708)
- file2 = mutant-enriched loops (total observed: 1,509,188)
- Ratio 1.053 matches your data's ctrl/mut ratio of 1.040

### What Hiccups Found

**Total differential loops**: 38,593
- 19,386 significantly higher in control
- 19,207 significantly higher in mutant

**Key differences from mariner approach**:

| Aspect | Hiccups Differential | Mariner + edgeR |
|--------|---------------------|-----------------|
| Method | FDR-based, pixel-level comparison | Buffer aggregation + RNA-seq statistics |
| Positions | Tests exact bin positions | Robust to ±2 bin shifts |
| Sensitivity | Very sensitive (FDR=0.0 for many) | Depends on edgeR FDR cutoff |
| Output | Binary (significant/not) | Continuous fold changes + FDR |

**Why different loop counts?**:
- Your merged set: 22,108 consensus positions (merged within 10kb)
- Hiccups differential: 38,593 positions (different merging/calling strategy)
- These sets partially overlap but aren't identical

**Critical note**: Hiccups found many small but significant differences. Your pre-edgeR analysis showed only 46 loops with >2-fold change, but hiccups found 38k differential. This suggests:
1. Many loops have small but statistically significant changes
2. edgeR should identify these with proper variance modeling
3. Statistical significance ≠ large fold change

---

## What Needs to Be Done Next - edgeR Analysis

### Primary Task: Differential Loop Analysis

**Input file**: `06_counts_matrix.rds` (22,108 loops × 2 samples)
- class: matrix; columns: "ctrl", "mut"
- Rownames: "loop_1", "loop_2", ... "loop_22108" (generic IDs, not coordinates)
- Genomic coordinates are in: outputs/full/03_binned.rds (GInteractions object)
	- extract chr, start1, end1, start2, end2 for final results table

**edgeR workflow**:
1. Create DGEList object
2. Filter low-count loops (if needed)
3. Calculate normalization factors
4. Estimate dispersions
5. Fit model
6. Test for differential loops
7. Extract results with FDR < 0.05

**Output needed**:
- Table with: loop_id, chr, start1, end1, start2, end2, logFC, FDR, significance
- How many loops are significant?
- Distribution of fold changes
- Comparison statistics
- TSV with (loop_coords, logFC, FDR, significance), MA plot, volcano plot
- Differential edgeR analysis

### Secondary Task: Systematic Comparison with Hiccups

**Goal**: Understand method agreement and differences

**Analyses needed**:

1. **Overlap Analysis**:
   - How many of your 22,108 loops appear in hiccups differential files?
   - Of those, how many does edgeR also call differential?
   - Concordance in direction (both say up in control?)

2. **Method-Specific Hits**:
   - Loops significant in mariner but NOT hiccups → Buffer caught shifts
   - Loops significant in hiccups but NOT mariner → Why?
   - Are shifted loops more likely to be mariner-specific?

3. **Effect Size Comparison**:
   - For shared differential loops, do fold changes correlate?
   - Is mariner finding larger or smaller effects?

4. **Validation Questions**:
   - Of the 11,924 loops with positional shifts, how many are differential?
   - Does buffer approach specifically help for shifted loops?

### Expected Outcomes

**Likely scenario**: 
- edgeR will find hundreds to thousands of significant loops
- High agreement with hiccups for strong signals
- Mariner may catch some shifted loops that hiccups missed
- Both methods capture similar biology with different sensitivities

**Key questions to answer**:
- Is mariner+edgeR more or less sensitive than hiccups?
- Do the methods agree on biological conclusions?
- Are there condition-specific loop changes that make biological sense?

---

## Technical Notes for Implementation

### Count Matrix Format

**Structure**: 
```r
# 22,108 rows (loops) × 2 columns (samples)
# Column names: "ctrl", "mut"
# Row names: loop IDs or coordinates
```

**Statistics from QC**:
```
Control: Min=0.21, Median=412.89, Mean=488.04, Max=3322.95
Mutant:  Min=0.22, Median=397.20, Mean=469.40, Max=21740.08
```

**No filtering applied yet** - all 22,108 loops retained for edgeR

### Comparison Data Loading

To compare with hiccups, you'll need:

```r
# Load hiccups differential calls
hiccups_ctrl <- read.table(
  "/expanse/lustre/projects/csd940/ctea/nf-hic/juicerpre/merged/diffhiccups_results/resorted_ctrl_vs_resorted_mut/file1/postprocessed_pixels_5000.bedpe",
  header = FALSE, skip = 2, stringsAsFactors = FALSE
)

hiccups_mut <- read.table(
  "/expanse/lustre/projects/csd940/ctea/nf-hic/juicerpre/merged/diffhiccups_results/resorted_ctrl_vs_resorted_mut/file2/postprocessed_pixels_5000.bedpe",
  header = FALSE, skip = 2, stringsAsFactors = FALSE
)

# Assign column names (same as your BEDPE format)
# Need to extract coordinates and match to your 22,108 loops
```

### Key Questions to Keep in Mind

1. **How does library size normalization affect results?**
   - Control/mutant ratio is 1.04 - very balanced
   - edgeR's TMM normalization should handle this

1. **Should you filter low-count loops?** - yes, bolded
   - Minimum count is 0.21 - very low
   - **Consider filtering loops with mean count < 10 or similar**
   - Balance sensitivity vs false positives

3. **What FDR threshold to use?**
   - **Standard: FDR < 0.05** - yes
   - More stringent: FDR < 0.01
   - Compare to hiccups which used FDR = 0.0 threshold

**Important**: Use default edgeR `filterByExpr()`!
Avoid manual handling where possible.

4. **How to handle biological replicates?**
   - Your data is already merged/consensus
   - edgeR treating as n=1 per condition
   - Will have low power - rely on strong signals

---

## Relevant References from Notes

From ct-meeting1.md discussions:
- "hiccups best for loop calling"
- "stuck" on "loop diff analysis" 
- "multi-hi-compare" mentioned but mariner chosen instead
- "replicate aware" analysis was a goal
- Focus on "differential loop analysis" at 5kb resolution
- "merged replicates and control.mutant" for comparison

**Current approach**: Using mariner's buffer + edgeR statistics as replicate-aware alternative to hiccups differential.

---

## Environment Details

**System**: Expanse SDSC
**Conda environment**: `mariner_env`
**R version**: 4.4.2
**Key packages**: mariner, edgeR, limma, GenomicRanges, InteractionSet

**Cluster job that completed**:
- Job ID: 43522753  
- Runtime: ~54 minutes (QC analysis)
- Resources: Sufficient for full dataset

---

## Summary - Where We Are

✅ **Completed**:
- Loop preparation and merging
- Buffering strategy implementation
- Hi-C data extraction
- Count matrix aggregation
- Comprehensive QC validation

🔴 **Next immediate task**: Run edgeR differential analysis

📊 **Then**: Compare results with hiccups differential calls to validate approach

🎯 **End goal**: Identify biologically meaningful differential loops between BAP1 mutant and control with validated statistical framework