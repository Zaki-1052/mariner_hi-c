# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**Mariner** is a differential chromatin loop analysis pipeline for Hi-C data comparing BAP1 mutant vs. control mouse samples. The pipeline quantifies Hi-C signal at loop positions across biological replicates (n=3 per condition) and uses edgeR's quasi-likelihood GLM framework for statistical testing.

**Key Innovation**: Treats loops as quantitative count data by extracting 5×5 pixel matrices around consensus positions, enabling replicate-aware statistical testing with proper variance modeling.

## Commands

### Run Full Pipeline

**Single resolution (5kb)**:
```bash
sbatch run_replicate_pipeline.sb
# Or locally:
Rscript scripts/prep_loops.R 5000
Rscript scripts/extract_counts.R 5000
Rscript scripts/aggregate.R 5000
Rscript scripts/edgeR.R 5000
```

**Multi-resolution (5kb, 10kb, 25kb)**:
```bash
sbatch run_multiresolution_pipeline.sb
# Runs all three resolutions sequentially (~6 hours)
```

**Resolution comparison**:
```bash
Rscript scripts/compare_resolutions.R
# Requires edgeR results from all three resolutions
```

### Quality Control

```bash
# Generate QC report
Rscript scripts/qc-val.R

# Check extraction quality
cat outputs/res_5kb/05_metadata.rds | grep "correlation"
```

### Format Conversion

```bash
# Convert to BEDPE for Juicebox visualization
python scripts/bedpe_convert.py \
  outputs/edgeR_results_res_5kb/primary_analysis/significant_loops_fdr05.tsv \
  outputs/bedpe/sig_5kb.bedpe \
  --fdr-cutoff 0.05
```

### Debugging

```bash
# Check SLURM job status
squeue -u $USER

# View logs
tail -f logs/mariner_multiRes_*.out

# Validate intermediate files
Rscript -e "x <- readRDS('outputs/res_5kb/03_binned.rds'); length(x)"
```

## Architecture

### Pipeline Flow (Sequential Stages)

```
1. prep_loops.R      → Merge 6 BEDPE files into consensus loops with 5×5 buffers
   Input:  loops_ctrl.bed, loops_mut.bed (6 replicates total)
   Output: outputs/res_{RES}kb/01-04_*.rds

2. extract_counts.R  → Pull Hi-C matrices from .hic files, store in HDF5
   Input:  04_buffered.rds + 6 .hic files
   Output: outputs/res_{RES}kb/05_extracted.h5 (HDF5-backed)

3. aggregate.R       → Aggregate 5×5 matrices to single count per loop
   Input:  05_extracted.h5
   Output: outputs/res_{RES}kb/06_counts_matrix.rds (loops × 6 samples)

4. edgeR.R          → Differential analysis via quasi-likelihood GLM
   Input:  06_counts_matrix.rds
   Output: edgeR_results_res_{RES}kb/ (TSV results + plots)

5. compare_resolutions.R → Meta-analysis across 5kb/10kb/25kb
   Input:  Results from all three resolutions
   Output: resolution_comparison/ (Venn diagrams, overlap matrices)
```

### Resolution-Aware Design Pattern

All scripts accept resolution as command-line argument and dynamically construct paths:

```r
args <- commandArgs(trailingOnly = TRUE)
RESOLUTION <- if (length(args) > 0) as.numeric(args[1]) else 5000

input_dir <- sprintf("outputs/res_%dkb", RESOLUTION/1000)
# Creates: outputs/res_5kb/, outputs/res_10kb/, outputs/res_25kb/
```

**Implication**: When modifying scripts, maintain this pattern. Don't hardcode resolution-specific paths.

### Data Structures

#### GInteractions (mariner)
- **Purpose**: Genomic interactions with paired anchor ranges
- **Created in**: `prep_loops.R` (read_bedpe → as_ginteractions)
- **Workflow**: Per-replicate list → merged consensus → binned → buffered
- **Key files**: 01_ginteractions.rds, 02_merged.rds, 03_binned.rds, 04_buffered.rds

#### HDF5Array (delayed evaluation)
- **Purpose**: Out-of-core storage for 5×5×22,108×6 matrix (~26GB if in-memory)
- **Created in**: `extract_counts.R` via `pullHicMatrices(..., onDisk=TRUE)`
- **Accessed in**: `aggregate.R` via `loadHDF5SummarizedExperiment()`
- **Pattern**: Lazy evaluation - materialize only 5×5 slices as needed
- **Critical**: Always clean HDF5 directory before fresh extraction to avoid dataset conflicts

#### DGEList (edgeR)
- **Purpose**: Count matrix + metadata for differential analysis
- **Structure**:
  - `counts`: 22,108 loops × 6 samples
  - `group`: factor(ctrl, ctrl, ctrl, mut, mut, mut)
  - `genes`: chr1, start1, end1, chr2, start2, end2, distance
- **Created in**: `aggregate.R` (DGEList construction)
- **Processed in**: `edgeR.R` (filtering, normalization, dispersion, testing)

### Configuration System

**Primary config**: `config/edgeR_config.yaml`
- Paths, sample metadata, statistical parameters
- Used by current replicate-aware pipeline

**Legacy config**: `config/edgeR_params.yaml`
- No-replicate analysis (fixed BCV=0.4)
- Kept for reference, not used in main pipeline

**Master framework**: `config/analysis_params.yaml`
- Comprehensive validation rules and QC thresholds
- Documents expected distributions and checks

**Loading pattern**:
```r
config <- yaml::read_yaml("config/edgeR_config.yaml")
# Command-line arg overrides paths dynamically:
input_dir <- sprintf("outputs/res_%dkb", RESOLUTION/1000)
```

### HDF5 Memory Management

**Problem**: Full dataset doesn't fit in memory
**Solution**: HDF5 + DelayedArray for lazy evaluation

**In extract_counts.R**:
```r
# ALWAYS clean before extraction
if (dir.exists(hdf5_dir)) {
  unlink(hdf5_dir, recursive = TRUE)  # Prevents dataset conflicts
}

pixels <- pullHicMatrices(
  ...,
  h5File = "temp_hdf5/extracted.h5",
  onDisk = TRUE,              # Store on disk
  compressionLevel = 1,       # Balance compression vs speed
  blockSize = 1e6             # 1MB chunks (tune for memory)
)
```

**In aggregate.R**:
```r
# Partial materialization - don't load full array
count_array <- counts(pixels)  # DelayedArray reference
for (i in 1:n_loops) {
  mat_5x5 <- as.matrix(count_array[, , i, j])  # Materialize 5×5 only
  counts_sum[i, j] <- sum(mat_5x5, na.rm = TRUE)
}
```

**Key principle**: Materialize minimally, only what's needed for current operation.

### Statistical Framework (edgeR)

**Current pipeline**: Replicate-aware (n=3 per condition)
```r
design <- model.matrix(~group)
y <- estimateDisp(y, design, robust=TRUE)  # Data-driven dispersion
fit <- glmQLFit(y, design, robust=TRUE)    # Quasi-likelihood fit
qlf <- glmQLFTest(fit, coef=2)             # Test mutant effect
```

**Key parameters**:
- `filterByExpr()`: min_count=5, min_total_count=15, min_prop=0.5
- Normalization: TMM (handles library size differences)
- Testing: Quasi-likelihood F-test (robust to outliers)
- FDR: 0.05 (primary), 0.01 (stringent)

**Advantages over no-replicate**:
- Data-driven dispersion (not fixed BCV=0.4)
- Residual df=4 for proper error rate control
- ~100× improvement in statistical power

## Critical Patterns and Constraints

### 1. Merge Conflict in extract_counts.R (Lines 79-87)
**Status**: Unresolved git conflict markers exist
**Impact**: Code still runs (one branch selected by interpreter)
**Action**: Clean up conflict markers before commits

### 2. Two-Stage Filtering Strategy
```r
# Stage 1: prep_loops.R - Remove junk loops
bedpe_filtered <- bedpe[bedpe$resolution == expected_resolution &
                        bedpe$observed > 0, ]

# Stage 2: edgeR.R - Remove low-signal loops
keep <- filterByExpr(y, group=y$samples$group)
y <- y[keep, , keep.lib.sizes=FALSE]
```
**Rationale**: BEDPE filtering removes technical artifacts, edgeR filtering optimizes statistical power.

### 3. Replicate Correlation Thresholds
```r
if (ctrl_within_mean > 0.95 && mut_within_mean > 0.95) {
  cat("✓ Excellent reproducibility\n")
} else if (ctrl_within_mean > 0.90 && mut_within_mean > 0.90) {
  cat("✓ Good reproducibility\n")
} else {
  cat("⚠ WARNING: Lower than expected correlation (<0.90)\n")
}
```
**Action on warning**: Check BEDPE quality, verify .hic file integrity, assess library quality.

### 4. HDF5 Dataset Conflicts
**Issue**: HDF5 datasets can't be overwritten
**Solution**: Always `unlink(hdf5_dir, recursive=TRUE)` before fresh extraction
**Code location**: `extract_counts.R:79-83`

### 5. Coordinate Matching Across Resolutions
**Issue**: Same loop has different anchor positions at 5kb vs 10kb vs 25kb
**Solution**: Use tolerance-based matching (±10kb window)
**Code location**: `compare_resolutions.R` uses fuzzy matching, not exact coordinates

### 6. Normalization Availability
**Issue**: KR normalization might not exist in .hic files
**Fallback**: Use VC (Vanilla Coverage) - always available
**Code location**: `extract_counts.R` checks `readHicNormTypes()` before extraction

## Quality Control Expectations

### Pre-Analysis Validation
- ✅ No NA values in count matrix
- ✅ No negative counts
- ✅ Expected loop count: ~22,108 at 5kb resolution
- ✅ All chromosomes present in anchors

### Post-Extraction Metrics
- **Sparsity**: Expect 30-60% zeros in 5×5 matrices
- **NA rate**: Should be <1%
- **Per-sample median**: Expect 50-200 counts per 5×5 matrix
- **Within-condition correlation**: >0.90 (ideally >0.95)
- **Between-condition correlation**: 0.70-0.90

### Post-Differential Results
- **P-value distribution**: Uniform under null (histogram should be flat for p>0.2)
- **LogFC symmetry**: Expect ~50/50 up/down regulated loops
- **Extreme logFC**: Flag if |logFC| > 5 (likely artifact)
- **Expected ranges**:
  - baseMean: 1-5000
  - logFC: -3 to +3
  - logCPM: 0-12

### Generated QC Plots
1. **MDS plot** - Sample clustering (groups should separate clearly)
2. **BCV plot** - Dispersion estimates (median BCV ~0.3-0.5 for Hi-C)
3. **QL dispersion** - Shrinkage toward prior
4. **MA plot** - Effect vs abundance (expect hourglass shape)
5. **Volcano plot** - Effect vs significance

## Common Development Scenarios

### Adding a New Aggregation Strategy

**File**: `scripts/aggregate.R`

```r
# Example: Add median aggregation
counts_median <- matrix(0, nrow=dims[3], ncol=dims[4])
for (i in 1:dims[3]) {
  for (j in 1:dims[4]) {
    mat_5x5 <- as.matrix(count_array[, , i, j])
    counts_median[i, j] <- median(mat_5x5, na.rm=TRUE)
  }
}

# Add to output
all_strategies$median <- counts_median
saveRDS(all_strategies, file.path(output_dir, "06_all_strategies.rds"))
```

**Then modify edgeR.R** to accept `--strategy median` argument.

### Testing a New Dispersion Value

**File**: `scripts/edgeR.R`

```r
# Current: Data-driven
y <- estimateDisp(y, design, robust=TRUE)

# Alternative: Fixed BCV for sensitivity analysis
bcv_fixed <- 0.3  # Test value
y$common.dispersion <- bcv_fixed^2
y$trended.dispersion <- rep(bcv_fixed^2, nrow(y))
y$tagwise.dispersion <- rep(bcv_fixed^2, nrow(y))
```

### Debugging Low Loop Retention

**Symptom**: "Warning: High filter rate (>50%)"

**Diagnosis**:
```r
# In edgeR.R, before filtering
zero_loops <- sum(rowSums(counts_matrix) == 0)
low_count_loops <- sum(rowSums(counts_matrix) < 10)
cat("Zero total: ", zero_loops, "\n")
cat("Low count (<10): ", low_count_loops, "\n")
```

**Solutions**:
1. Check BEDPE input quality (observed > threshold)
2. Verify .hic files match expected resolution
3. Lower `min_count` in filterByExpr (current: 5)
4. Check for systematic bias (one condition >> other)

### Adding a New Resolution

**Pattern**: Resolution argument propagates through all scripts

```bash
# Add 1kb resolution
Rscript scripts/prep_loops.R 1000
Rscript scripts/extract_counts.R 1000
Rscript scripts/aggregate.R 1000
Rscript scripts/edgeR.R 1000

# Output appears in: outputs/res_1kb/
```

**Caution**: 1kb requires 25× more loops than 5kb → memory/time increase

### Modifying Sample Metadata

**File**: `config/edgeR_config.yaml`

```yaml
samples:
  names: ["ctrl_M1", "ctrl_M2", "ctrl_M3", "mut_M1", "mut_M2", "mut_M3"]
  groups: ["ctrl", "ctrl", "ctrl", "mut", "mut", "mut"]
```

**If adding a 4th replicate**:
```yaml
samples:
  names: ["ctrl_M1", "ctrl_M2", "ctrl_M3", "ctrl_M4", "mut_M1", "mut_M2", "mut_M3", "mut_M4"]
  groups: ["ctrl", "ctrl", "ctrl", "ctrl", "mut", "mut", "mut", "mut"]
```

**Also update**:
- BEDPE file list in `prep_loops.R`
- .hic file paths in `extract_counts.R`

## Resource Requirements (HPC)

### SLURM Allocations

| Stage | CPUs | Memory | Time | Notes |
|-------|------|--------|------|-------|
| prep_loops.R | 4 | 16GB | 1 min | CPU rarely limiting |
| extract_counts.R | 16 | 64GB | 15 min | Memory-intensive, tune blockSize |
| aggregate.R | 4 | 8GB | 5 min | Light computation |
| edgeR.R | 4 | 16GB | 5 min | Statistical testing |

**Full pipeline (single resolution)**: 6:00:00 allocation (conservative, actual ~30 min)
**Multi-resolution pipeline**: 6:00:00 (runs 3 resolutions sequentially, actual ~2 hours)

### Scaling Considerations

- **Memory scales with**: Number of loops × number of samples × matrix size
- **Time scales with**: .hic file size and resolution (finer = more bins to extract)
- **Disk scales with**: HDF5 storage (~500MB per resolution)

**Current dataset**: 22,108 loops × 6 samples × 5×5 pixels = ~26GB uncompressed
**HDF5 compressed**: ~500MB on disk
**RAM usage**: <8GB due to lazy evaluation

## Known Limitations

1. **No inter-chromosomal loops**: Pipeline filters to intra-chromosomal only
2. **Fixed buffer size**: 5×5 matrices hardcoded (±10kb at 5kb resolution)
3. **No batch effect correction**: Assumes all samples processed identically
4. **Binary comparison only**: Two-group design (ctrl vs mut), no multi-factor
5. **Hiccups dependency**: Requires pre-called loops, not ab initio detection

## Documentation

Detailed stage-specific documentation in `docs/`:
- `mariner.md` - Pipeline overview
- `prep_loops.md` - Loop merging methodology
- `extract_counts.md` - Hi-C extraction details
- `edgeR.md` - Statistical framework
- `compare_resolutions.md` - Multi-resolution analysis
- `config.md` - Configuration reference
- `qc.md` - Quality control procedures

## File Naming Conventions

### Intermediate RDS Files
```
01_ginteractions.rds    # Per-replicate GInteractions list
02_merged.rds           # Consensus merged loops
03_binned.rds           # Binned to resolution grid
04_buffered.rds         # 5×5 buffered regions
05_extracted.h5         # HDF5-backed matrices
06_counts_matrix.rds    # Final count matrix (loops × samples)
```

**Pattern**: Sequential numbering indicates pipeline order. Don't skip stages.

### Output Results
```
edgeR_results_res_{RES}kb/
  primary_analysis/
    all_results_primary.tsv              # All loops
    significant_loops_fdr05.tsv          # FDR < 0.05
    significant_loops_fdr01.tsv          # FDR < 0.01
    top100_differential.tsv              # Top by |logFC|
  plots/
    ma_plot_primary.pdf
    volcano_plot_primary.pdf
    mds_plot.pdf
    bcv_plot.pdf
```

**Pattern**: Resolution-specific directories allow parallel analysis of 5kb/10kb/25kb without conflicts.
