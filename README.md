# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a Hi-C differential chromatin loop analysis pipeline using the **mariner** R package. The pipeline identifies differential chromatin loops between BAP1 mutant and control samples through a multi-stage workflow: loop preparation, Hi-C signal extraction, aggregation, quality control, and statistical testing with edgeR.

**Key biological concept:** The pipeline addresses a critical Hi-C analysis challenge: chromatin loops may be detected at slightly different genomic positions between samples (±1 bin or ±5kb) due to binning artifacts. The solution is to extract a "buffer zone" (5×5 pixel regions) around each loop position, then aggregate the signal for statistical testing.

## Pipeline Architecture

The pipeline follows a strictly sequential 6-step workflow, with multi-resolution support at each stage:

### Core Pipeline Stages

1. **prep_loops.R** - Loop preparation and consensus building
   - Reads BEDPE files from Hiccups at specified resolution
   - Merges loops across biological replicates (n=3 per condition) within 10kb radius
   - Bins merged loops to resolution-specific grid
   - Creates 5×5 pixel buffers (±2 bins) around each loop position
   - Input: `hiccups_results/*_M[1-3]/postprocessed_pixels_{RES}.bedpe`
   - Output: `outputs/res_{RES}kb/01-04_*.rds`

2. **extract_counts.R** - Hi-C signal extraction
   - Extracts contact matrices from .hic files at buffered positions
   - Processes 6 biological replicates (3 ctrl + 3 mut)
   - Uses HDF5 backend for memory efficiency with large datasets
   - Validates extraction with correlation analysis
   - Input: `.hic` files + buffered regions from stage 1
   - Output: `outputs/res_{RES}kb/05_extracted.h5` (HDF5-backed array)

3. **aggregate.R** - Buffer aggregation
   - Reduces 5×5 pixel matrices to single counts per loop
   - Implements two strategies: sum (primary) and center-weighted
   - Validates replicate correlations and library sizes
   - Creates edgeR-ready DGEList object with genomic annotations
   - Input: Extracted 5×5 matrices from stage 2
   - Output: `outputs/res_{RES}kb/06_counts_matrix.rds` + edgeR input

4. **qc-val.R** - Quality control and validation
   - Comprehensive QC with sample correlations, distributions, spatial analysis
   - Matrix-level validation: peak centering, center enrichment
   - **Critical:** Detects and saves per-loop shift status for downstream enrichment testing
   - Generates 9 diagnostic plots and QC report
   - Input: All outputs from stages 1-3
   - Output: `outputs/res_{RES}kb/qc_report/` directory with plots and `qc_report_summary.rds`

5. **edgeR.R** - Differential loop analysis (n=3 replicates)
   - Replicate-aware statistical testing with quasi-likelihood GLM
   - Data-driven dispersion estimation from biological replicates
   - Robust methods to handle outliers
   - Tests for enrichment of shifted loops in differential set
   - Input: Count matrix + coordinates + QC summary (for shift status)
   - Output: `outputs/edgeR_results_res_{RES}kb/` with results tables and plots

6. **compare_resolutions.R** - Multi-resolution comparison (optional)
   - Compares differential results across 5kb, 10kb, 25kb resolutions
   - Identifies high-confidence loops (significant at multiple resolutions)
   - Input: edgeR results from multiple resolutions
   - Output: `outputs/multiresolution_comparison/`

### Data Flow

```
BEDPE files (29k ctrl, 28k mut per replicate)
    ↓ prep_loops.R [filter, merge, bin, buffer]
Buffered GInteractions (150 loops for test, ~57k for full)
    ↓ extract_counts.R [pullHicMatrices from .hic files]
5×5×N_loops×6 HDF5 array (5×5 pixels, N loops, 6 replicates)
    ↓ aggregate.R [sum or weighted aggregation]
N_loops×6 count matrix (ready for edgeR)
    ↓ edgeR.R [QL-GLM testing with n=3 per group]
Differential loop results (FDR < 0.05)
```

## Multi-Resolution Support

**All pipeline scripts accept resolution as command-line argument:**
```bash
Rscript scripts/prep_loops.R 5000      # 5kb resolution
Rscript scripts/extract_counts.R 10000  # 10kb resolution
Rscript scripts/aggregate.R 25000       # 25kb resolution
```

**Resolution-specific paths:**
- Input: `hiccups_results/*/postprocessed_pixels_{RES}.bedpe`
- Output: `outputs/res_{RES}kb/`
- Results: `outputs/edgeR_results_res_{RES}kb/`

## Commands

### Run Full Pipeline (Single Resolution)
```bash
# Full pipeline at 5kb resolution
Rscript scripts/prep_loops.R 5000
Rscript scripts/extract_counts.R 5000
Rscript scripts/aggregate.R 5000
Rscript scripts/qc-val.R 5000          # Optional but recommended
Rscript scripts/edgeR.R 5000

# Change resolution by replacing 5000 with 10000 or 25000
```

### Run Multi-Resolution Pipeline
```bash
# Process all three resolutions in parallel (SLURM)
sbatch scripts/run_multiresolution_pipeline.sb

# Or run sequentially
for res in 5000 10000 25000; do
  Rscript scripts/prep_loops.R $res
  Rscript scripts/extract_counts.R $res
  Rscript scripts/aggregate.R $res
  Rscript scripts/qc-val.R $res
  Rscript scripts/edgeR.R $res
done

# Then compare results
Rscript scripts/compare_resolutions.R
```

## Critical Architecture Patterns

### 1. Resolution-Specific Directory Structure
All outputs are organized by resolution to support multi-resolution analysis:
```
outputs/
  res_5kb/
    01_ginteractions.rds
    02_merged.rds
    03_binned.rds
    04_buffered.rds
    05_extracted/  (HDF5 directory)
    05_metadata.rds
    06_counts_matrix.rds
    qc_report/
      qc_report_summary.rds  # Contains shift status
  res_10kb/
    [same structure]
  edgeR_results_res_5kb/
    primary_analysis/
      all_results_primary.tsv
      significant_loops_fdr05.tsv
    plots/
```

### 2. Biological Replicate Structure
All scripts handle **6 biological replicates** (3 ctrl + 3 mut):
- Sample order MUST be: `ctrl_M1, ctrl_M2, ctrl_M3, mut_M1, mut_M2, mut_M3`
- Count matrices are always N_loops × 6 (columns match sample order)
- Configuration file `config/edgeR_config.yaml` specifies sample metadata
- Replicate-aware correlation analysis in aggregation and QC stages

### 3. HDF5-Backed Arrays for Memory Efficiency
Extract counts uses HDF5 backend to handle large datasets:
```r
# In extract_counts.R
pixels <- pullHicMatrices(
  x = buffered,
  files = hicFiles,  # 6 .hic files
  h5File = h5_file_path,
  onDisk = TRUE,
  compressionLevel = 1
)

# Load in later scripts
pixels <- loadHDF5SummarizedExperiment(
  dir = input_dir,
  prefix = "05_extracted"
)
```

### 4. Buffer Aggregation Strategies
Two aggregation methods implemented:
- **Sum** (primary): Simple sum of all 25 pixels - robust to shifts
- **Weighted**: Gaussian-like weights favoring center pixel - balances precision and shift-tolerance
- Both strategies saved in `06_all_strategies.rds` for comparison

### 5. Shift Detection and Enrichment Testing
**Critical workflow:**
1. QC script (`qc-val.R`) analyzes 5×5 matrices and detects positional shifts
2. Saves boolean vector `loop_shift_status` in `qc_report_summary.rds`
3. edgeR script (`edgeR.R`) loads shift status and tests enrichment via Fisher's exact test
4. Tests hypothesis: "Are shifted loops more/less likely to be differential?"

### 6. Configuration-Driven edgeR Analysis
`config/edgeR_config.yaml` controls all statistical parameters:
- Sample metadata (names, groups)
- Filtering thresholds (`min_count`, `min_total_count`, `min_prop`)
- FDR thresholds (primary: 0.05, exploratory: 0.10)
- Fold change categories (strong: |logFC| > 1, moderate: |logFC| > 0.5)
- Visualization settings (colors, plot dimensions)

## Key Data Structures

### GInteractions Objects
Used throughout pipeline to represent chromatin loops:
```r
# Structure: paired genomic ranges (anchor1 ↔ anchor2)
GInteractions(
  anchor1: GRanges (chr, start, end),
  anchor2: GRanges (chr, start, end),
  metadata: observed, FDR values, etc.
)

# Created from BEDPE: as_ginteractions()
# Merged across replicates: mergePairs()
# Binned to grid: assignToBins()
# Buffered to regions: pixelsToMatrices()
```

### DGEList (edgeR Object)
Structure for differential testing:
```r
DGEList(
  $counts: N_loops × 6 matrix,
  $samples: data.frame with lib.size, norm.factors, group,
  $genes: data.frame with chr1, start1, end1, chr2, start2, end2
)
```

## Common Development Tasks

### Adding a New Aggregation Strategy
1. Modify `scripts/aggregate.R` section "AGGREGATION STRATEGY"
2. Add new strategy to `strategies` list
3. Update `all_strategies` object before saving
4. QC script will automatically compare correlations

### Changing Statistical Thresholds
Edit `config/edgeR_config.yaml`:
```yaml
statistics:
  fdr_primary: 0.05      # Primary FDR cutoff
  fdr_exploratory: 0.10  # Exploratory threshold
  fold_change_thresholds:
    strong: 1.0          # 2-fold change
    moderate: 0.5        # 1.4-fold change
```

### Running on Subset for Testing
Modify filtering in `prep_loops.R`:
```r
# Line ~38: Reduce loop count for testing
bedpe_subset <- bedpe_5kb[1:100, ]  # Test with 100 loops
```

### Validating New Resolution
1. Ensure Hiccups BEDPE files exist at that resolution
2. Check .hic files contain the resolution: `readHicBpResolutions()`
3. Run pipeline with new resolution value
4. Expected output directory: `outputs/res_{RES}kb/`

## Important Design Principles

### 1. Fail-Fast Validation
Scripts validate inputs immediately:
```r
# Check file exists
if (!file.exists(filepath)) {
  stop(sprintf("ERROR: File not found: %s", filepath))
}

# Check dimensions match
stopifnot(
  "Count matrix and coordinates dimension mismatch" =
    nrow(counts_matrix) == length(binned_gi)
)
```

### 2. Comprehensive Progress Logging
All scripts print detailed progress:
```r
cat("\n=== Stage Name ===\n")
cat(sprintf("Processing %d loops...\n", n_loops))
cat(sprintf("✓ Completed in %.1f seconds\n", elapsed_time))
```

### 3. Resolution-Aware Output Paths
Every script constructs resolution-specific paths:
```r
RESOLUTION <- if (length(args) > 0) as.numeric(args[1]) else 5000
input_dir <- sprintf("outputs/res_%dkb", RESOLUTION/1000)
output_dir <- sprintf("outputs/edgeR_results_res_%dkb", RESOLUTION/1000)
```

### 4. Data-Driven Statistical Testing
edgeR analysis with biological replicates:
- **No fixed BCV assumption** (unlike merged n=1 approach)
- Dispersion estimated from data: `estimateDisp(y, design, robust=TRUE)`
- Quasi-likelihood GLM for rigorous FDR control: `glmQLFit()` + `glmQLFTest()`
- Robust estimation to downweight outliers

## Troubleshooting Common Issues

### "File not found" errors
- Check resolution matches BEDPE filenames
- Verify .hic file paths in extract_counts.R (lines 37-44)
- Ensure previous pipeline stages completed successfully

### "Dimension mismatch" errors
- Likely from incomplete pipeline run
- Re-run from stage 1 with same resolution
- Check that all intermediate .rds files have matching loop counts

### "Very few significant loops"
Possible causes and solutions:
- **High BCV (>0.6)**: Check replicate quality in MDS plot
- **Weak biological differences**: Review sample correlations in QC report
- **Over-filtering**: Relax `min_count` in `config/edgeR_config.yaml`

### HDF5 loading issues
```r
# If HDF5 backend fails, check file existence:
list.files(file.path(input_dir, "temp_hdf5"))

# Reload with explicit path:
pixels <- loadHDF5SummarizedExperiment(
  dir = input_dir,
  prefix = "05_extracted"
)
```

## Performance Characteristics

### Computational Requirements
- **Memory**: 8-32 GB RAM (depends on loop count and resolution)
- **Storage**: ~1-5 GB per resolution (HDF5 arrays are largest)
- **Runtime**:
  - prep_loops: 2-5 minutes
  - extract_counts: 5-15 minutes (depends on .hic file size)
  - aggregate: 1-3 minutes
  - qc-val: 3-10 minutes (matrix analysis is slow)
  - edgeR: 5-10 minutes

### Scalability
- Test mode (100 loops): ~5 minutes total
- Full dataset (57k loops): ~30-45 minutes total
- Multi-resolution (3 resolutions): ~2-3 hours total

## Expected Results

### Pipeline Success Indicators
1. **QC report passes** with within-group correlation > 0.90
2. **BCV estimates** between 0.1-0.4 (typical for Hi-C)
3. **3,000-5,000 differential loops** at FDR < 0.05 (for n=3 replicates)
4. **Improvement over merged analysis**: ~1,500× more differential loops detected

### Typical Output Metrics
```
Loops tested:     18,000-22,000 (after filtering)
Significant:      3,000-5,000 (15-25% at FDR < 0.05)
Up in mutant:     1,500-2,500
Down in mutant:   1,200-2,000
Positional shifts: 10-30% of loops
BCV:              0.25-0.35 (data-driven estimate)
```

## Next Steps After Analysis

1. Review MDS plot to confirm replicate quality
2. Check BCV and QL dispersion plots for data behavior
3. Examine MA and volcano plots for biological patterns
4. Compare with Hiccups differential results if available
5. Annotate significant loops with nearby genes
6. Validate top hits with biological context or orthogonal methods

## Important Files Not to Modify

- `.hic` files: Read-only input from Juicer
- `*_ginteractions.rds`: Intermediate objects, regenerate if needed
- HDF5 arrays: Generated by mariner, do not edit manually

## Documentation References

- Pipeline documentation: `docs/mariner.md` (detailed stage-by-stage walkthrough)
- edgeR analysis: `docs/edgeR.md` (comprehensive statistical methods)
- Coding standards: `docs/guidelines.md`
- Configuration: `docs/config.md`
