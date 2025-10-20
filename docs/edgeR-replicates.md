# edgeR Analysis Issue Documentation

## Current Analysis State

### What Was Actually Run

**Pipeline executed**: `01_run_edgeR_analysis.R`  
**Job ID**: 43532373  
**Date**: 2025-10-20  
**Runtime**: 28 seconds

**Input data used**:
```
Count matrix: outputs/full/06_counts_matrix.rds (22,108 loops × 2 samples)
Coordinates: outputs/full/03_binned.rds
Samples: "ctrl" and "mut"
```

**Sample structure**: n=1 per condition (merged/consensus samples)

### Results Obtained

**Primary analysis (BCV = 0.15)**:
- Loops tested: 22,104 (after filtering 4 loops)
- Significant at FDR < 0.05: **2 loops**
  - 1 up in mutant
  - 1 down in mutant
  - Both are strong differential (|logFC| > 1)

**Sensitivity analysis**:
- BCV = 0.10: 40 significant loops (31 up, 9 down)
- BCV = 0.15: 2 significant loops (1 up, 1 down)
- BCV = 0.20: 1 significant loop (1 up)

**Shifted loop analysis**: No enrichment detected (Fisher's p = 0.503)

### Why These Results Are Limited

**Statistical constraint**: Without biological replicates, edgeR cannot estimate variance from the data. Instead, it uses a fixed BCV value that assumes a certain level of biological variability.

**Effect size distribution from prior QC**:
- Mean log2FC: -0.066 ± 0.220
- Only 46 loops with |logFC| > 1
- Only 2 loops with |logFC| > 2

**Interpretation**: The biological changes between control and mutant are predominantly subtle (small fold changes). Without replicates to estimate variance, edgeR can only detect the most extreme outliers.

---

## Actual Data Structure Available

### Replicate Files Discovered

**Location**: `/expanse/lustre/projects/csd940/ctea/nf-hic/juicer_frompre/hiccups_results`

**Control replicates**:
- `ctrl_M1/` - contains postprocessed_pixels_5000.bedpe, merged_loops.bedpe
- `ctrl_M2/` - same structure
- `ctrl_M3/` - same structure

**Mutant replicates**:
- `mut_M1/` - contains postprocessed_pixels_5000.bedpe, merged_loops.bedpe
- `mut_M2/` - same structure
- `mut_M3/` - same structure

**Merged versions** (what was actually used in current analysis):
- `resorted_ctrl/` - consensus across ctrl_M1, M2, M3
- `resorted_mut/` - consensus across mut_M1, M2, M3

### Files in Each Replicate Directory

```
enriched_pixels_5000.bedpe
enriched_pixels_10000.bedpe
enriched_pixels_25000.bedpe
postprocessed_pixels_5000.bedpe
postprocessed_pixels_10000.bedpe
postprocessed_pixels_25000.bedpe
merged_loops.bedpe
fdr_thresholds_5000
fdr_thresholds_10000
fdr_thresholds_25000
```

---

## Current Pipeline Design vs Available Data

### What the Pipeline Was Built For

From `edgeR-prep.md`:
```
"Designed for merged replicate data (n=1 per condition)"
"Key Challenge: Hi-C loop calls can shift positions by 1-2 bins"
"Mariner's buffer approach addresses this"
```

From `prepare_loops.R` workflow:
1. Read merged BEDPE files (resorted_ctrl, resorted_mut)
2. Filter for 5kb resolution
3. Merge loops within 10kb radius to create consensus positions
4. Extract 5×5 matrices for all 22,108 consensus loops
5. Aggregate to counts (sum of 25 pixels)
6. Output: 22,108 × 2 matrix

### What Actually Exists

**True experimental design**: 3 biological replicates per condition
- Control: M1, M2, M3
- Mutant: M1, M2, M3

**Data that was used**: Merged/consensus versions that collapsed replicates

---

## Comparison Context

### Hiccups Differential Results

**Location**: `/expanse/lustre/projects/csd940/ctea/nf-hic/juicerpre/merged/diffhiccups_results/resorted_ctrl_vs_resorted_mut/`

**Found**:
- 19,386 loops enriched in control (file1/)
- 19,207 loops enriched in mutant (file2/)
- Total: 38,593 differential loops

**Method**: Hiccups' built-in differential analysis on merged samples

### Current edgeR Results vs Hiccups

- Hiccups: 38,593 differential loops
- edgeR (BCV=0.15): 2 differential loops
- edgeR (BCV=0.10): 40 differential loops

**Discrepancy factor**: ~1000x difference at primary BCV

---

## Technical Quality Indicators

From QC analysis (job 43522753):

**Sample correlation**:
- Pearson: 0.9015
- Spearman: 0.9790

**Library composition**:
- Control total: 10,789,577 counts
- Mutant total: 10,377,560 counts  
- Ratio: 1.040 (highly balanced)

**Data integrity**: All passed
- No NA values
- No negative values
- No zero-sum loops (except 4 filtered)

**Spatial properties**:
- 54% of loops show positional shifts between conditions
- Mean center enrichment: 0.074-0.075 (expected for 5kb resolution)

---

## Configuration Used

From `config/edgeR_config.yaml`:

```yaml
statistics:
  primary_bcv: 0.15
  sensitivity_bcv:
    lower: 0.10
    upper: 0.20
  
  filtering:
    min_count: 5
    min_total_count: 10
  
  fdr_primary: 0.05
  normalization_method: "TMM"
```

---

## Code That Ran

**Script**: `scripts/01_run_edgeR_analysis.R`

**Key sections executed**:
1. Data loading (counts, coordinates, QC)
2. DGEList creation with n=2 samples
3. filterByExpr (removed 4 loops)
4. TMM normalization
5. exactTest with fixed dispersion (no replicates)
6. Sensitivity analysis across 3 BCV values
7. Visualization and output generation

**No errors encountered** - analysis completed successfully with the merged data structure it was designed for.

---

## Goal Statement

**Objective**: Use the actual biological replicates (M1, M2, M3 per condition) instead of merged consensus samples to enable proper variance estimation in edgeR.

**Expected benefit**: With n=3 per group, edgeR can estimate dispersion from the data rather than assuming a fixed BCV, potentially increasing power to detect differential loops.

**Required modification**: Pipeline needs to be adapted to:
1. Work with individual replicate BEDPE/HiC files
2. Extract counts for all 6 samples (not 2)
3. Let edgeR estimate biological variance from replicates