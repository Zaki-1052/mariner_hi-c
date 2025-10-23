# Mariner: Differential Chromatin Loop Analysis Pipeline

A robust bioinformatics pipeline for identifying differential chromatin loops in Hi-C data using the edgeR statistical framework. This pipeline addresses bin-shifting artifacts in loop comparison through spatial buffering and provides multi-resolution analysis capabilities.

## Table of Contents

- [Overview](#overview)
- [Scientific Background](#scientific-background)
- [Key Features](#key-features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Pipeline Stages](#pipeline-stages)
- [Configuration](#configuration)
- [Usage Examples](#usage-examples)
- [Output Files](#output-files)
- [Quality Control](#quality-control)
- [HPC Deployment](#hpc-deployment)
- [Documentation](#documentation)
- [Citation](#citation)

## Overview

**Mariner** is a comprehensive pipeline for comparative Hi-C analysis that identifies statistically significant differences in chromatin loop structures between biological conditions. The current implementation compares BAP1 wild-type versus BAP1 mutant mouse samples across multiple resolutions (5kb, 10kb, 25kb).

### The Problem

When comparing chromatin loops between conditions, the same biological loop can be detected at slightly different genomic positions due to:
- Bin size discretization (±5kb shifts)
- Loop detection algorithm variations
- Stochastic peak position variation

### The Solution

Instead of requiring exact coordinate matches, Mariner extracts **5×5 bin matrices** (25kb × 25kb) centered on each loop position, creating a "buffer zone" that captures loops even if shifted by ±10kb. This enables:

- **Robust loop merging** across replicates and conditions
- **Sensitive detection** of differential loops
- **Quantification** via Hi-C contact frequencies within buffered regions
- **Statistical rigor** using edgeR's GLM framework

## Scientific Background

### Experimental Design

| Sample | Condition | Replicates | Description |
|--------|-----------|------------|-------------|
| `ctrl` | Control   | 3 (M1, M2, M3) | BAP1 wild-type |
| `mut`  | Mutant    | 3 (M1, M2, M3) | BAP1 mutant |

- **Total loops analyzed**: ~22,108 (at 5kb resolution)
- **Sample correlation**: Pearson 0.90, Spearman 0.98
- **Statistical approach**: Conservative dispersion estimates (BCV=0.4) for no-replicate design

### Why edgeR for Loop Data?

edgeR is traditionally used for RNA-seq differential expression, but adapts well to Hi-C loop analysis:
- Handles count-based data with appropriate statistical distributions
- Robust normalization (TMM) accounts for library size differences
- GLM framework enables complex experimental designs
- Conservative dispersion estimates prevent false positives with limited replicates

## Key Features

✅ **Multi-resolution analysis** - Supports 5kb, 10kb, and 25kb resolutions
✅ **Spatial tolerance** - Buffer zones handle ±10kb bin-shifting artifacts
✅ **Statistical rigor** - edgeR GLM + LRT with FDR correction
✅ **Flexible aggregation** - Multiple strategies (sum, max, center, Gaussian blur)
✅ **Comprehensive QC** - Extensive validation checks and diagnostic plots
✅ **HPC-ready** - SLURM job submission scripts for Expanse cluster
✅ **Visualization support** - BEDPE output for Juicebox integration
✅ **Hiccups comparison** - Benchmark against external loop caller

## Installation

### System Requirements

- **R** ≥ 4.1.0
- **Python** ≥ 3.8
- **HPC environment** (optional, for large-scale analysis)

### R Dependencies

Install required Bioconductor packages:

```r
# Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Core Bioconductor packages
BiocManager::install(c(
    "mariner",           # Hi-C loop analysis
    "edgeR",             # Differential expression
    "InteractionSet",    # Genomic interactions
    "GenomicRanges",     # Coordinate operations
    "rhdf5"              # HDF5 file handling
))

# CRAN packages
install.packages(c(
    "yaml",              # Configuration parsing
    "ggplot2",           # Plotting
    "dplyr",             # Data manipulation
    "tibble",            # Data frames
    "RColorBrewer",      # Color palettes
    "scales"             # Plot scaling
))
```

### Python Dependencies

```bash
pip install pandas numpy argparse
```

### Clone Repository

```bash
git clone https://github.com/yourusername/mariner-newest.git
cd mariner-newest
```

## Quick Start

### 1. Prepare Your Data

Organize input files:

```
mariner-newest/
├── loops_ctrl.bed       # Control condition loops (BEDPE format)
├── loops_mut.bed        # Mutant condition loops (BEDPE format)
├── data/
│   ├── ctrl_rep1.hic   # Control Hi-C contact matrices
│   ├── ctrl_rep2.hic
│   ├── ctrl_rep3.hic
│   ├── mut_rep1.hic    # Mutant Hi-C contact matrices
│   ├── mut_rep2.hic
│   └── mut_rep3.hic
```

### 2. Configure Analysis

Edit `config/analysis_params.yaml`:

```yaml
project:
  name: "BAP1_differential_loops"
  organism: "Mus musculus"
  resolution: 5000  # 5kb resolution

files:
  input_dir: "."
  output_dir: "outputs"
  hic_files:
    ctrl: "data/ctrl_merged.hic"
    mut: "data/mut_merged.hic"

edgeR:
  bcv: 0.4              # Biological coefficient of variation
  fdr_threshold: 0.05   # False discovery rate cutoff
```

### 3. Run Pipeline

**Single resolution (5kb):**

```bash
# Submit to SLURM
sbatch run.sb

# Or run locally
Rscript scripts/prep_loops.R 5000
Rscript scripts/extract_counts.R 5000
Rscript scripts/aggregate.R 5000
Rscript scripts/edgeR.R 5000
```

**Multi-resolution (5kb, 10kb, 25kb):**

```bash
sbatch run_multiresolution_pipeline.sb
```

### 4. Convert to BEDPE for Visualization

```bash
python scripts/bedpe_convert.py \
  outputs/edgeR_results_res_5kb/primary_analysis/all_results.tsv \
  outputs/bedpe/5kb_significant.bedpe \
  --fdr-cutoff 0.05 \
  --exclude-unchanged
```

## Pipeline Stages

### Stage 1: Loop Preparation (`prep_loops.R`)

**Purpose**: Merge loops from multiple replicates into consensus positions

**Input**: 6 BEDPE files (3 ctrl + 3 mut replicates)

**Process**:
1. Load and filter BEDPE files by resolution
2. Remove zero-count loops
3. Convert to GInteractions format
4. Merge loops within 10kb radius (select highest-signal representative)
5. Snap positions to 5kb grid boundaries
6. Create 5×5 bin buffer zones

**Output**: `outputs/res_5kb/04_buffered.rds`

**Key parameters**:
- `binSize`: 5000 (5kb bins)
- `radius`: 10000 (10kb merging tolerance)
- `buffer`: 2 bins (creates 5×5 matrices)

### Stage 2: Hi-C Data Extraction (`extract_counts.R`)

**Purpose**: Extract Hi-C contact matrices for each buffered loop

**Input**: Buffered loop regions + .hic files

**Process**:
1. Validate .hic files and resolution availability
2. Select VC (Vanilla Coverage) normalization
3. Extract 5×5 contact matrix for each loop from both samples
4. Store in HDF5 format for memory efficiency
5. Validate extraction quality (check NAs, sparsity)

**Output**: `outputs/res_5kb/05_extractedse.rds` (HDF5-backed)

**Key considerations**:
- Memory management via on-disk HDF5 storage
- Normalization method affects downstream comparisons
- Validation ensures data integrity before statistics

### Stage 3: Matrix Aggregation (`aggregate.R`)

**Purpose**: Convert 5×5 matrices to single count values per loop

**Input**: Extracted 5×5 Hi-C matrices

**Aggregation strategies**:
- **Sum** (default): Total contact frequency across 5×5 region
- **Max**: Maximum pixel value (peak detection)
- **Center**: Center 3×3 region only (reduce noise)
- **Gaussian blur**: Weighted sum with Gaussian kernel

**Output**: `outputs/res_5kb/06_counts_matrix.rds` (loops × samples)

**Recommendation**: Use **sum** for most analyses (captures total loop signal)

### Stage 4: Differential Analysis (`edgeR.R`)

**Purpose**: Identify statistically significant differential loops

**Input**: Count matrix + sample metadata

**Process**:
1. **Filtering**: Remove low-count loops (< 10 total counts)
2. **Normalization**: TMM (Trimmed Mean of M-values)
3. **Dispersion estimation**: Fixed BCV=0.4 for no-replicate design
4. **Statistical testing**: GLM + LRT (Likelihood Ratio Test)
5. **Multiple testing correction**: FDR (Benjamini-Hochberg)
6. **Sensitivity analysis**: Test across BCV 0.2-0.5 range

**Output**:
- `all_results_primary.tsv` - All loops with statistics
- `significant_loops_fdr05.tsv` - FDR < 0.05
- `significant_loops_fdr01.tsv` - FDR < 0.01
- MA plots, volcano plots, dispersion plots

**Statistical parameters**:
- **BCV**: 0.4 (conservative mammalian estimate)
- **Dispersion**: 0.16 (BCV²)
- **FDR threshold**: 0.05 (primary), 0.01 (stringent)
- **Effect size**: |logFC| > 1 (2-fold change)

### Stage 5: Format Conversion (`bedpe_convert.py`)

**Purpose**: Convert results to BEDPE for Juicebox visualization

**Input**: edgeR TSV results

**Output**: BEDPE files with statistical annotations

```bash
python scripts/bedpe_convert.py INPUT.tsv OUTPUT.bedpe \
  --fdr-cutoff 0.05 \
  --logfc-cutoff 1.0 \
  --exclude-unchanged
```

## Configuration

### Main Configuration: `config/analysis_params.yaml`

Comprehensive 377-line YAML file covering:

**Project metadata**:
```yaml
project:
  name: "BAP1_differential_loops"
  description: "Differential loop analysis comparing BAP1 mutant vs control"
  organism: "Mus musculus"
  resolution: 5000
```

**File paths**:
```yaml
files:
  input_dir: "."
  output_dir: "outputs"
  hic_files:
    ctrl: "path/to/ctrl.hic"
    mut: "path/to/mut.hic"
```

**edgeR parameters**:
```yaml
edgeR:
  filtering:
    method: "filterByExpr"
    min_count: 10
    min_total_count: 20
  normalization:
    method: "TMM"
  dispersion:
    bcv: 0.4
    dispersion: 0.16
  testing:
    method: "glmLRT"
    fdr_threshold: 0.05
```

**QC thresholds**:
```yaml
qc:
  pre_filter:
    max_na_allowed: 0
    max_negative_allowed: 0
  post_filter:
    min_loops_retained: 1000
    max_filter_rate: 0.5
```

See `docs/config.md` for complete documentation.

### edgeR-Specific Config: `config/edgeR_config.yaml`

Focused configuration for edgeR analysis:

```yaml
input:
  counts_matrix: "outputs/res_5kb/06_counts_matrix.rds"
  loop_coords: "outputs/res_5kb/03_binned.rds"

samples:
  ctrl:
    name: "Control"
    group: "ctrl"
  mut:
    name: "BAP1_mutant"
    group: "mut"

statistics:
  bcv: 0.4
  fdr_cutoffs: [0.01, 0.05, 0.1]
  logfc_cutoff: 1.0
```

## Usage Examples

### Example 1: Standard Analysis at 5kb

```bash
# Run complete pipeline
Rscript scripts/prep_loops.R 5000
Rscript scripts/extract_counts.R 5000
Rscript scripts/aggregate.R 5000
Rscript scripts/edgeR.R 5000

# Generate QC report
Rscript scripts/qc-val.R

# Convert significant loops to BEDPE
python scripts/bedpe_convert.py \
  outputs/edgeR_results_res_5kb/primary_analysis/significant_loops_fdr05.tsv \
  outputs/bedpe/sig_5kb.bedpe
```

### Example 2: Multi-Resolution Comparison

```bash
# Run pipeline at all resolutions
sbatch run_multiresolution_pipeline.sb

# Compare results across resolutions
Rscript scripts/compare_resolutions.R

# Output: outputs/resolution_comparison/overlap_analysis.tsv
```

### Example 3: Custom Dispersion Testing

```r
# Load custom script
source("scripts/edgeR.R")

# Test multiple BCV values
bcv_values <- c(0.2, 0.3, 0.4, 0.5)
for (bcv in bcv_values) {
  run_edgeR_analysis(
    counts_matrix = "outputs/res_5kb/06_counts_matrix.rds",
    bcv = bcv,
    output_dir = paste0("outputs/sensitivity_bcv", bcv)
  )
}
```

### Example 4: Hiccups Comparison

```r
# Compare with external Hiccups results
source("scripts/edgeR.R")

# Load Hiccups differential loops
hiccups_file <- "data/hiccups_differential.bedpe"

# Run comparison analysis (included in edgeR.R)
# Outputs Venn diagrams and overlap statistics
```

## Output Files

### Directory Structure

```
outputs/
├── res_5kb/                          # 5kb resolution intermediates
│   ├── 01_ginteractions.rds         # Initial GInteractions
│   ├── 02_merged.rds                # Merged consensus loops
│   ├── 03_binned.rds                # Snapped to 5kb grid
│   ├── 04_buffered.rds              # 5×5 buffered regions
│   ├── 05_extractedse.rds           # HDF5-backed matrices
│   └── 06_counts_matrix.rds         # Final count matrix
│
├── edgeR_results_res_5kb/            # 5kb differential results
│   ├── primary_analysis/
│   │   ├── all_results_primary.tsv  # All loops + statistics
│   │   ├── significant_loops_fdr05.tsv
│   │   ├── significant_loops_fdr01.tsv
│   │   ├── ma_plot.pdf
│   │   ├── volcano_plot.pdf
│   │   └── dispersion_sensitivity.pdf
│   └── hiccups_comparison/
│       ├── overlap_exact.tsv
│       ├── overlap_10kb.tsv
│       └── venn_diagrams.pdf
│
├── bedpe/                            # Juicebox-compatible BEDPE
│   ├── 5kb_all.bedpe
│   ├── 5kb_significant_fdr05.bedpe
│   └── 5kb_upregulated_fdr01.bedpe
│
├── resolution_comparison/            # Multi-resolution analysis
│   ├── overlap_matrix.tsv
│   └── resolution_sensitivity.pdf
│
└── qc_report/                        # Quality control
    ├── qc_summary.txt
    ├── correlation_plots.pdf
    └── validation_checks.pdf
```

### Key Output Formats

#### TSV Results Table

| Column | Description |
|--------|-------------|
| `seqnames1`, `start1`, `end1` | First anchor (chr, start, end) |
| `seqnames2`, `start2`, `end2` | Second anchor (chr, start, end) |
| `logFC` | Log2 fold-change (mut vs ctrl) |
| `logCPM` | Log2 counts per million |
| `PValue` | Raw p-value from LRT |
| `FDR` | Adjusted p-value (Benjamini-Hochberg) |
| `baseMean` | Mean normalized count across samples |

#### BEDPE Format (for Juicebox)

```
chr1  start1  end1  chr2  start2  end2  name  score  strand1  strand2  [custom_fields]
```

Custom fields include: `logFC`, `FDR`, `PValue`, `significance`

## Quality Control

### Automated QC Checks

The pipeline includes extensive validation (`scripts/qc-val.R`):

**Pre-filtering checks**:
- ✅ No NA values in count matrix
- ✅ No negative counts
- ✅ Expected column names present
- ✅ Expected loop count ranges

**Post-filtering validation**:
- ⚠️ Alert if >50% loops filtered
- ⚠️ Ensure minimum 1000 loops retained

**Post-differential validation**:
- 📊 P-value distribution (should be uniform under null)
- 📊 LogFC symmetry (expect ~50/50 up/down)
- ⚠️ Flag extreme logFC (|logFC| > 5)

**Expected distributions**:
- `baseMean`: 1-5000
- `logFC`: -3 to +3
- `logCPM`: 0-12

### Generate QC Report

```bash
Rscript scripts/qc-val.R
# Output: outputs/qc_report/qc_summary.txt
```

### Sample Correlation

Expected metrics for good data quality:
- **Pearson correlation**: > 0.85
- **Spearman correlation**: > 0.90
- **Library size ratio**: 0.8-1.2

## HPC Deployment

### SLURM Job Submission (SDSC Expanse)

**Single resolution**:

```bash
#!/bin/bash
#SBATCH --job-name=mariner_5kb
#SBATCH --account=csd940
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=6:00:00
#SBATCH --output=logs/mariner_%j.out
#SBATCH --error=logs/mariner_%j.err

module load r/4.2.1
module load python/3.9.0

Rscript scripts/prep_loops.R 5000
Rscript scripts/extract_counts.R 5000
Rscript scripts/aggregate.R 5000
Rscript scripts/edgeR.R 5000
```

**Multi-resolution**:

```bash
sbatch run_multiresolution_pipeline.sb
# Runs 5kb, 10kb, 25kb sequentially
# Estimated runtime: 6 hours
```

### Resource Requirements

| Stage | CPUs | Memory | Time |
|-------|------|--------|------|
| `prep_loops.R` | 4 | 16GB | 30min |
| `extract_counts.R` | 16 | 64GB | 2hr |
| `aggregate.R` | 4 | 8GB | 10min |
| `edgeR.R` | 4 | 16GB | 20min |

**Note**: Memory requirements scale with number of loops and Hi-C matrix size.

## Documentation

Comprehensive documentation is available in the `docs/` directory:

- **[mariner.md](docs/mariner.md)** - Pipeline overview and technical details
- **[prep_loops.md](docs/prep_loops.md)** - Loop preparation and merging
- **[extract_counts.md](docs/extract_counts.md)** - Hi-C data extraction
- **[edgeR.md](docs/edgeR.md)** - edgeR statistical analysis methodology
- **[compare_resolutions.md](docs/compare_resolutions.md)** - Multi-resolution analysis
- **[config.md](docs/config.md)** - Configuration file reference
- **[qc.md](docs/qc.md)** - Quality control procedures
- **[guidelines.md](docs/guidelines.md)** - Coding standards and best practices

## Troubleshooting

### Common Issues

**Issue**: "Error: .hic file not found"
**Solution**: Verify paths in `config/analysis_params.yaml` and ensure .hic files exist

**Issue**: "No loops retained after filtering"
**Solution**: Lower `min_count` threshold or check input loop quality

**Issue**: "HDF5 directory error"
**Solution**: Create `temp_hdf5/` directory or set `HDF5_DIR` environment variable

**Issue**: "Memory allocation error during extraction"
**Solution**: Increase SLURM `--mem` parameter or reduce number of loops

### Debug Mode

Enable verbose logging:

```r
options(verbose = TRUE)
options(warn = 2)  # Treat warnings as errors
```

## Citation

If you use this pipeline in your research, please cite:

**Mariner package**:
```
Kramer E, Davis E, Simeonov K, et al. (2024)
Mariner: Explore the Hi-Cs.
Bioconductor. doi:10.18129/B9.bioc.mariner
```

**edgeR**:
```
Robinson MD, McCarthy DJ, Smyth GK (2010)
edgeR: a Bioconductor package for differential expression analysis of digital gene expression data.
Bioinformatics 26(1):139-140. doi:10.1093/bioinformatics/btp616
```

## License

This pipeline is available under the MIT License. See LICENSE file for details.

## Contact

For questions or issues:
- **GitHub Issues**: [https://github.com/yourusername/mariner-newest/issues](https://github.com/yourusername/mariner-newest/issues)
- **Email**: your.email@institution.edu

## Acknowledgments

- **SDSC Expanse** cluster for computational resources
- **Bioconductor** community for Hi-C analysis tools
- **BAP1 project team** for biological insights and data generation
