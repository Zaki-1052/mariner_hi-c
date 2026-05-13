# Mariner Hi-C Differential Chromatin Loop Analysis Pipeline

## Overview

This pipeline performs **replicate-aware differential analysis of chromatin loops** from Hi-C data using the [mariner](https://bioconductor.org/packages/mariner/) R/Bioconductor package and [edgeR](https://bioconductor.org/packages/edgeR/) statistical framework. It compares chromatin looping between **BAP1-KO mutant** and **wildtype control** conditions in mouse (mm10 genome) with **biological replication** (n=3 per condition).

### Key Features

- **Multi-resolution analysis**: Processes loops at 5kb, 10kb, and 25kb resolutions
- **Biological replicates**: True replicate-aware statistics with n=3 per condition (6 samples total)
- **Data-driven dispersion**: Uses biological variation for robust differential testing
- **Comprehensive QC**: Extensive quality control with correlation analysis and shift detection
- **ChIP-seq integration**: Anchor classification using H3K27ac and H3K4me1 ChIP-seq data
- **Publication-ready outputs**: Volcano plots, MA plots, enrichment analyses, and APA heatmaps
- **Complete reproducibility**: Single-command execution with full logging

### Biological Context

**Chromatin loops** are 3D interactions between distant genomic regions that regulate gene expression by bringing enhancers into contact with promoters. This pipeline identifies loops that are differentially regulated between BAP1-KO and wildtype conditions, potentially revealing mechanisms of BAP1-mediated gene regulation.

**Key biological processes analyzed:**
- Loop formation and stability
- Enhancer-promoter interactions
- Chromatin architecture reorganization
- Gene regulatory network changes

---

## Pipeline Architecture

### Workflow Overview

```
Input: HiCCUPS loop calls (BEDPE) + Hi-C contact matrices (.hic)
  ↓
[1] prep_loops.R          → Merge 6 replicates, create consensus positions
  ↓
[2] extract_counts.R       → Extract 5×5 Hi-C matrices at each loop
  ↓
[3] aggregate.R            → Sum matrices → count matrix (loops × samples)
  ↓
[4] qc-val.R               → Quality control validation & diagnostics
  ↓
[5] edgeR.R                → Quasi-likelihood GLM differential analysis
  ↓
[6] compare_resolutions.R  → Cross-resolution comparison & overlap
  ↓
[7] generate_final_results.py → Apply stringent filters (|logFC|>0.3, FDR<0.03)
  ↓
[8] convert_final_bedpe.sh → Convert to BEDPE for Juicebox visualization
  ↓
[9] downstream_analysis.R  → Merge resolutions, annotate with genes & ChIP-seq
  ↓
[10] visualizations.R      → Publication plots, enrichment, volcano plots
  ↓
[11] apa_analysis.R        → Aggregate Peak Analysis (APA) heatmaps
  ↓
Output: Differential loops, visualizations, annotations, statistics
```

### Sample Structure

```
Control (wildtype):        Mutant (BAP1-KO):
├── ctrl_M1               ├── mut_M1
├── ctrl_M2               ├── mut_M2
└── ctrl_M3               └── mut_M3

Each sample has:
  - HiCCUPS loop calls: postprocessed_pixels_{resolution}.bedpe
  - Hi-C contact matrix: {sample}.hic
```

---

## Installation

### Requirements

**Software:**
- R ≥ 4.1.0
- Python ≥ 3.8
- Java (for Juicebox, optional)

**R Packages:**
```r
# Bioconductor core
BiocManager::install(c(
  "mariner",           # Hi-C loop analysis
  "edgeR",             # Differential analysis
  "InteractionSet",    # Genomic interactions
  "GenomicRanges",     # Genomic intervals
  "HDF5Array",         # Out-of-memory arrays
  "strawr",            # .hic file reading
  "ChIPseeker",        # Peak annotation
  "clusterProfiler",   # Enrichment analysis
  "TxDb.Mmusculus.UCSC.mm10.knownGene",  # Mouse gene annotations
  "org.Mm.eg.db"       # Mouse gene ontology
))

# CRAN packages
install.packages(c(
  "tidyverse",
  "ggplot2",
  "pheatmap",
  "viridis",
  "patchwork",
  "EnhancedVolcano",
  "yaml",
  "gtools"
))
```

**Python Packages:**
```bash
pip install pandas numpy
```

### Setup

```bash
# Clone repository
git clone <repository-url>
cd mariner_hi-c

# Install R packages (one-time setup)
Rscript install_R_packages.R

# Verify conda environment (if using SLURM on HPC)
conda activate mariner_env
```

---

## Quick Start

### Single-Command Execution (Recommended)

Run the complete pipeline across all three resolutions:

```bash
# Submit to SLURM (HPC environments)
sbatch scripts/run_full_pipeline.sb

# Or run directly (local/server)
bash scripts/run_full_pipeline.sb
```

**Expected runtime:** ~24 hours (includes time-intensive APA analysis)

**Output:** Complete differential analysis for 5kb, 10kb, and 25kb resolutions with all visualizations and statistics.

### Step-by-Step Execution

For manual control or debugging:

```bash
# Process each resolution individually
for RES in 5000 10000 25000; do
  echo "Processing ${RES} bp resolution..."

  # Step 1: Prepare loops (merge replicates)
  Rscript scripts/prep_loops.R ${RES}

  # Step 2: Extract Hi-C counts
  Rscript scripts/extract_counts.R ${RES}

  # Step 3: Aggregate matrices
  Rscript scripts/aggregate.R ${RES}

  # Step 4: Quality control
  Rscript scripts/qc-val.R ${RES}

  # Step 5: Differential analysis
  Rscript scripts/edgeR.R ${RES}
done

# Multi-resolution post-processing
Rscript scripts/compare_resolutions.R
python3 scripts/generate_final_results.py
bash scripts/convert_final_bedpe.sh
Rscript scripts/downstream_analysis.R
Rscript scripts/visualizations.R
Rscript scripts/apa_analysis.R
```

---

## Detailed Pipeline Steps

### Step 1: Loop Preparation (`prep_loops.R`)

**Purpose:** Merge HiCCUPS loop calls from 6 biological replicates into consensus positions.

**Input:**
- 6 BEDPE files (one per replicate) at specified resolution
- `/expanse/lustre/projects/csd940/ctea/nf-hic/juicer_frompre/hiccups_results/{sample}/postprocessed_pixels_{RES}.bedpe`

**Process:**
1. Load BEDPE files and filter for target resolution
2. Convert to GInteractions objects
3. Merge loops within 10kb radius using `mergePairs()`
4. Assign to genomic bins (resolution-specific)
5. Create 5×5 pixel buffers around each loop

**Output:** `outputs/res_{RES}kb/`
- `01_ginteractions.rds` - Individual replicate loops
- `02_merged.rds` - Consensus loop positions
- `03_binned.rds` - Resolution-aligned loops
- `04_buffered.rds` - 5×5 pixel regions (ready for extraction)

**Key biological insight:** The merging strategy creates a union set capturing loops from ANY replicate, maximizing discovery while handling technical variation in loop calling.

**Documentation:** See [docs/prep_loops.md](docs/prep_loops.md)

---

### Step 2: Hi-C Count Extraction (`extract_counts.R`)

**Purpose:** Extract 5×5 contact matrices at each consensus loop position from all 6 .hic files.

**Input:**
- Buffered loops from Step 1
- 6 .hic files: `/expanse/lustre/projects/csd940/ctea/nf-hic/juicerpre/{sample}.hic`

**Process:**
1. Validate .hic files and normalization methods
2. Extract 5×5 matrices using `pullHicMatrices()` with VC normalization
3. Store as HDF5-backed array for memory efficiency
4. Calculate per-sample and cross-sample correlations

**Output:** `outputs/res_{RES}kb/`
- `05_extracted/` - HDF5 SummarizedExperiment (dimensions: 5×5×loops×6)
- `05_metadata.rds` - Extraction parameters and QC metrics

**Key technical detail:** The 5×5 buffer handles positional uncertainty and bin-shift artifacts between replicates, ensuring robust signal capture.

**Documentation:** See [docs/extract_counts.md](docs/extract_counts.md)

---

### Step 3: Matrix Aggregation (`aggregate.R`)

**Purpose:** Convert 5×5 matrices to single count values per loop for differential analysis.

**Input:**
- HDF5 matrices from Step 2

**Process:**
1. **Sum aggregation** (primary): Sum all 25 pixels per loop
2. **Weighted aggregation** (alternative): Center-weighted sum with Gaussian-like weights
3. Calculate replicate correlations (within-group and between-group)
4. Create edgeR DGEList with proper group structure

**Output:** `outputs/res_{RES}kb/`
- `06_counts_matrix.tsv` - Count matrix (loops × 6 samples)
- `06_edgeR_input.rds` - DGEList ready for analysis
- `06_all_strategies.rds` - Both aggregation methods for comparison

**Expected quality metrics:**
- Within-group correlation: > 0.95 (excellent) or > 0.90 (good)
- Between-group correlation: 0.05-0.15 lower than within-group
- Library sizes: within 2-fold of each other

**Documentation:** See [docs/aggregate.md](docs/aggregate.md)

---

### Step 4: Quality Control Validation (`qc-val.R`)

**Purpose:** Comprehensive validation of data quality and technical reproducibility.

**Input:**
- All outputs from Steps 1-3

**QC Analyses:**
1. **Data Integrity**: NA values, zeros, negatives, outliers
2. **Sample Correlations**: Pearson/Spearman across all 6 replicates
3. **Library Sizes**: Total counts per sample
4. **MA Statistics**: Average expression vs. fold change
5. **Loop Variability**: Coefficient of variation across samples
6. **Spatial Analysis**: Distance decay, chromosome distribution
7. **Matrix Centering**: Peak position within 5×5 buffers
8. **Shift Detection**: Positional differences between conditions

**Output:** `outputs/res_{RES}kb/qc_report/`
- 9 diagnostic plots (PDF)
- `qc_report_summary.rds` - All QC metrics
- `loop_shift_status.rds` - Per-loop shift detection

**Pass Criteria:**
- ✓ Sample correlation > 0.85
- ✓ Library size ratio < 2-fold
- ✓ No negative values
- ✓ Clear biological signal (within > between correlations)

**Documentation:** See [docs/qc.md](docs/qc.md)

---

### Step 5: Differential Analysis (`edgeR.R`)

**Purpose:** Identify chromatin loops with significant differential contact frequency between conditions.

**Method:** Quasi-likelihood negative binomial GLM with robust dispersion estimation

**Input:**
- DGEList from Step 3
- QC metrics from Step 4

**Statistical Workflow:**
1. **Filtering**: `filterByExpr()` - Remove low-count loops
2. **Normalization**: TMM (Trimmed Mean of M-values)
3. **Dispersion Estimation**: Common, trended, and tagwise dispersions
4. **Model Fitting**: `glmQLFit()` with design `~group`
5. **Testing**: `glmQLFTest()` for mutant effect
6. **FDR Control**: Benjamini-Hochberg correction

**Output:** `outputs/edgeR_results_res_{RES}kb/`

**Primary Analysis:**
- `all_results_primary.tsv` - Full results table (all loops)
- `significant_loops_fdr05.tsv` - FDR < 0.05
- `final_results.tsv` - Stringent: |logFC| > 0.3 AND FDR < 0.03
- `top100_differential.tsv` - Top 100 by effect size
- `summary_statistics.txt` - Analysis summary

**Visualizations:**
- `mds_plot.pdf` - Sample relationships
- `bcv_plot.pdf` - Biological coefficient of variation
- `ql_dispersion_plot.pdf` - Quasi-likelihood dispersions
- `ma_plot_primary.pdf` - Average vs. fold change
- `volcano_plot_primary.pdf` - Fold change vs. significance
- `results_summary.pdf` - Bar plots of differential categories

**Expected Results:**
- **With replicates (n=3):** 500-5,000+ significant loops (FDR < 0.05)
- **Without replicates (n=1):** ~2 significant loops
- **Improvement:** 250-2,500× increase in statistical power

**Key Effect Size Categories:**
- **Strong**: |logFC| > 1.0 (2-fold change)
- **Moderate**: |logFC| > 0.5
- **Weak**: |logFC| ≤ 0.5 but FDR < 0.05

**Documentation:** See [docs/edgeR.md](docs/edgeR.md)

---

### Step 6: Resolution Comparison (`compare_resolutions.R`)

**Purpose:** Compare differential results across 5kb, 10kb, and 25kb resolutions.

**Input:**
- edgeR results from all three resolutions
- Genomic coordinates for each resolution

**Analyses:**
1. **Summary Statistics**: Loop counts, significance rates, effect sizes
2. **Coordinate Overlap**: Match loops across resolutions (10kb tolerance)
3. **Fold Change Correlation**: Agreement of effect sizes
4. **Venn Diagrams**: Overlap of differential loops
5. **Resolution-Specific Effects**: Unique to each resolution

**Output:** `outputs/resolution_comparison/`
- `summary_statistics_by_resolution.tsv`
- `loop_overlap_matrix.tsv`
- `differential_loops_by_resolution.pdf`
- `foldchange_correlation_5kb_vs_10kb.pdf` (and other pairs)
- `venn_diagram_differential_loops.pdf`

**Expected Patterns:**
- **5kb**: Most loops, short-range interactions enriched
- **10kb**: Intermediate coverage, balanced
- **25kb**: Fewest loops, long-range and inter-TAD interactions

**Documentation:** See [docs/compare_resolutions.md](docs/compare_resolutions.md)

---

### Step 7: Generate Final Results (`generate_final_results.py`)

**Purpose:** Apply stringent thresholds to identify high-confidence differential loops.

**Filters:**
- |logFC| > 0.3 (1.23-fold change)
- FDR < 0.03 (more stringent than standard 0.05)

**Input:** `all_results_primary.tsv` from each resolution

**Output:** `final_results.tsv` in each `edgeR_results_res_{RES}kb/primary_analysis/`

**Rationale:** These stringent thresholds reduce false positives and focus on loops with both statistical significance and meaningful biological effect sizes.

---

### Step 8: BEDPE Conversion (`convert_final_bedpe.sh`)

**Purpose:** Convert TSV results to BEDPE format for visualization in Juicebox.

**Output:** `outputs/bedpe_final/`
- `5kb_final.bedpe` - High-confidence differential loops at 5kb
- `10kb_final.bedpe` - High-confidence at 10kb
- `25kb_final.bedpe` - High-confidence at 25kb
- `merged_final.bedpe` - Union of all resolutions (non-redundant)
- `{RES}kb_all_loops.bedpe` - All tested loops (for context)
- `merged_all_loops.bedpe` - Union of all tested loops

**Juicebox Visualization:**
1. Open Juicebox
2. Load .hic file
3. Load .bedpe as 2D annotations
4. Color loops by logFC or significance

---

### Step 9: Downstream Analysis (`downstream_analysis.R`)

**Purpose:** Merge differential loops across resolutions and annotate with genomic features.

**Input:**
- `final_results.tsv` from all resolutions
- Mouse gene annotations (TxDb.Mmusculus.UCSC.mm10.knownGene)
- ChIP-seq peaks: H3K27ac (active marks), H3K4me1 (enhancer marks)

**Process:**
1. **Merge Resolutions**: Combine 5kb, 10kb, 25kb with overlap removal (10kb tolerance)
2. **Nearest Gene Annotation**: Find closest TSS for each anchor
3. **ChIP-seq Classification**:
   - **Promoter**: H3K27ac+ AND ≤2kb from TSS
   - **Active Enhancer**: H3K27ac+ AND >2kb from TSS
   - **Poised Enhancer**: H3K4me1+ (no H3K27ac) AND >2kb from TSS
   - **Other**: No ChIP-seq marks
4. **Loop Type Classification**: 10 categories based on anchor pairs
   - Promoter-Promoter (P-P)
   - Promoter-Active_Enhancer (P-E)
   - Promoter-Poised_Enhancer, etc.

**Output:** `outputs/merged_loops/`
- `non_redundant_loops.tsv` - Merged differential loops
- `characterized_loops.tsv` - With gene + ChIP-seq annotations
- `non_redundant_loops.bedpe` - BEDPE format
- `overlap_report.tsv` - Resolution overlap statistics
- `plots/` - Venn diagrams, distance distributions

**Key Insight:** This creates a unified, non-redundant set of differential loops with biological context for interpretation.

---

### Step 10: Visualization Analysis (`visualizations.R`)

**Purpose:** Generate publication-quality figures and perform functional enrichment.

**Input:**
- `characterized_loops.tsv` from Step 9
- `all_results_primary.tsv` from each resolution

**Visualizations Generated:**

**1. Volcano Plots** (`outputs/visualizations/volcano/`)
- EnhancedVolcano plots for each resolution
- Multi-resolution merged plot
- Color-coded by significance and direction

**2. Feature Distribution** (`outputs/visualizations/features/`)
- ChIPseeker annotation plots
- Distance to TSS distributions
- Genomic feature enrichment

**3. Enrichment Analysis** (`outputs/visualizations/enrichment/`)
- GO (Gene Ontology) enrichment
- KEGG pathway enrichment
- Dot plots and bar plots of top terms

**4. Loop Classification** (`outputs/visualizations/loop_classification/`)
- Loop type distribution (P-P, P-E, E-E, etc.)
- Differential by loop type
- Comparison across resolutions

**5. APA Previews** (if `apa_analysis.R` already run)
- Selected APA heatmaps
- Enrichment score comparisons

**Publication-Ready Formats:**
- High-resolution PDFs (300 dpi)
- Consistent color schemes
- Clear legends and annotations

---

### Step 11: Aggregate Peak Analysis (`apa_analysis.R`)

**Purpose:** Visualize aggregate Hi-C signal enrichment at differential loops.

**Method:** APA (Aggregate Peak Analysis) - average Hi-C contact matrices across many loops

**Input:**
- Differential loops (merged and resolution-specific)
- .hic files for all 6 samples

**Process:**
1. Extract larger matrices (e.g., 21×21 pixels for 50kb window)
2. Aggregate across loops by averaging
3. Calculate enrichment: P2LL = center / local background
4. Statistical comparison: ctrl vs. mut
5. Generate heatmaps

**Output:** `outputs/apa_results/res_{RES}kb/`

For each resolution and loop set:
- `apa_heatmap_ctrl.pdf` - Control average
- `apa_heatmap_mut.pdf` - Mutant average
- `apa_heatmap_difference.pdf` - Mut - Ctrl
- `enrichment_scores.tsv` - P2LL scores per loop
- `enrichment_comparison.pdf` - Violin plots ctrl vs. mut
- `statistical_test.txt` - Wilcoxon test results

**Expected Patterns:**
- **Upregulated loops**: Higher enrichment in mutant
- **Downregulated loops**: Higher enrichment in control
- **Visual confirmation** of differential contact frequency

**Runtime:** ~2 hours (time-intensive due to matrix extraction)

---

## Output Directory Structure

```
mariner_hi-c/
├── outputs/
│   ├── res_5kb/                          # 5kb resolution
│   │   ├── 01_ginteractions.rds
│   │   ├── 02_merged.rds
│   │   ├── 03_binned.rds
│   │   ├── 04_buffered.rds
│   │   ├── 05_extracted/                 # HDF5 matrices
│   │   ├── 05_metadata.rds
│   │   ├── 06_counts_matrix.tsv          # Count matrix
│   │   ├── 06_edgeR_input.rds
│   │   └── qc_report/                    # QC plots + metrics
│   │
│   ├── res_10kb/                         # 10kb resolution (same structure)
│   ├── res_25kb/                         # 25kb resolution (same structure)
│   │
│   ├── edgeR_results_res_5kb/
│   │   ├── primary_analysis/
│   │   │   ├── all_results_primary.tsv   # All loops
│   │   │   ├── final_results.tsv         # Stringent filters
│   │   │   ├── significant_loops_fdr05.tsv
│   │   │   ├── top100_differential.tsv
│   │   │   └── summary_statistics.txt
│   │   ├── plots/                        # MA, volcano, MDS, BCV plots
│   │   └── logs/
│   │
│   ├── edgeR_results_res_10kb/           # (same structure)
│   ├── edgeR_results_res_25kb/           # (same structure)
│   │
│   ├── resolution_comparison/
│   │   ├── summary_statistics_by_resolution.tsv
│   │   ├── loop_overlap_matrix.tsv
│   │   ├── venn_diagram_differential_loops.pdf
│   │   └── foldchange_correlation_*.pdf
│   │
│   ├── bedpe_final/
│   │   ├── 5kb_final.bedpe               # For Juicebox
│   │   ├── 10kb_final.bedpe
│   │   ├── 25kb_final.bedpe
│   │   ├── merged_final.bedpe            # Non-redundant union
│   │   └── *_all_loops.bedpe
│   │
│   ├── merged_loops/
│   │   ├── non_redundant_loops.tsv       # Merged across resolutions
│   │   ├── characterized_loops.tsv       # + gene + ChIP-seq annotations
│   │   ├── non_redundant_loops.bedpe
│   │   ├── overlap_report.tsv
│   │   └── plots/
│   │
│   ├── visualizations/
│   │   ├── volcano/                      # Volcano plots (all resolutions)
│   │   ├── features/                     # Genomic feature distribution
│   │   ├── enrichment/                   # GO/KEGG enrichment
│   │   └── loop_classification/          # P-P, P-E, E-E types
│   │
│   └── apa_results/
│       ├── res_5kb/
│       │   ├── merged_loops/             # APA for merged set
│       │   └── resolution_specific/      # APA for 5kb-only loops
│       ├── res_10kb/
│       └── res_25kb/
│
├── scripts/                              # All pipeline scripts
├── config/                               # Configuration files
├── docs/                                 # Detailed documentation
└── logs/                                 # Pipeline execution logs
```

---

## Key Results Files

### For Primary Analysis
- `outputs/edgeR_results_res_{RES}kb/primary_analysis/all_results_primary.tsv` - **Complete results**
- `outputs/edgeR_results_res_{RES}kb/primary_analysis/final_results.tsv` - **High-confidence differential**
- `outputs/merged_loops/characterized_loops.tsv` - **Annotated merged loops**

### For Visualization
- `outputs/visualizations/volcano/*.pdf` - **Volcano plots**
- `outputs/edgeR_results_res_{RES}kb/plots/ma_plot_primary.pdf` - **MA plots**
- `outputs/apa_results/**/*.pdf` - **APA heatmaps**

### For Juicebox
- `outputs/bedpe_final/merged_final.bedpe` - **Non-redundant differential loops**
- `outputs/bedpe_final/{RES}kb_final.bedpe` - **Resolution-specific**

---

## Configuration

### Main Configuration File

`config/edgeR_config.yaml` - Controls all analysis parameters

**Key settings:**
```yaml
statistics:
  fdr_primary: 0.05              # Primary FDR threshold
  fdr_exploratory: 0.10          # Exploratory threshold
  normalization_method: "TMM"    # Trimmed Mean of M-values

  filtering:
    min_count: 10                # Minimum counts per loop
    min_total_count: 15          # Minimum total across samples

  fold_change_thresholds:
    strong: 1.0                  # |logFC| > 1.0
    moderate: 0.5                # |logFC| > 0.5

samples:
  names: ["ctrl_M1", "ctrl_M2", "ctrl_M3", "mut_M1", "mut_M2", "mut_M3"]
  groups: ["ctrl", "ctrl", "ctrl", "mut", "mut", "mut"]
```

### File Paths Configuration

Edit `scripts/prep_loops.R`, `extract_counts.R`, and `apa_analysis.R` to update:
- HiCCUPS BEDPE file locations
- .hic file paths
- ChIP-seq BED file paths

---

## Quality Control Checkpoints

### After Step 1 (prep_loops.R)
✓ Check cluster size distribution - should have good overlap across replicates
✓ Verify loop counts are reasonable (~10,000-30,000 per resolution)
✓ Ensure similar contribution from each replicate

### After Step 3 (aggregate.R)
✓ **Within-group correlation > 0.90** (preferably > 0.95)
✓ **Between-group correlation 0.05-0.15 lower** than within
✓ Library sizes within 2-fold
✓ No systematic bias between conditions

### After Step 4 (qc-val.R)
✓ Review all QC plots in `qc_report/`
✓ Check `qc_report_summary.rds` for pass/fail status
✓ Verify MDS plot shows samples clustering by condition

### After Step 5 (edgeR.R)
✓ BCV (biological coefficient of variation) reasonable (~0.2-0.5)
✓ MDS plot shows clear condition separation
✓ Volcano plot has reasonable number of significant loops
✓ Summary statistics match expectations

---

## Troubleshooting

### Low Within-Group Correlation (< 0.90)
**Possible causes:**
- High biological variability
- Technical batch effects
- Poor quality replicate
- Different sequencing depths

**Solutions:**
- Check individual replicate quality
- Review raw .hic file statistics
- Consider removing outlier replicate
- Increase sequencing depth

### Very Few Significant Loops (< 100)
**Possible causes:**
- Weak biological effect
- Insufficient sequencing depth
- Overly stringent filtering
- Too few replicates

**Solutions:**
- Use exploratory FDR threshold (0.10)
- Reduce `min_count` in filtering
- Check BCV - high BCV reduces power
- Verify replicates are truly biological (not technical)

### Very High Significant Loops (> 10,000)
**Possible causes:**
- Strong biological effect (expected for major regulators)
- Batch effects (ctrl vs. mut from different batches)
- Different sequencing depths not properly normalized

**Solutions:**
- Verify this is biologically plausible
- Check MDS plot for batch effects
- Ensure TMM normalization is working
- Validate top hits with visualization

### Memory Issues
**Solutions:**
- Use HDF5 on-disk storage (already implemented)
- Process resolutions sequentially instead of parallel
- Reduce `blockSize` in `pullHicMatrices()`
- Use HPC with more RAM (current: 128GB)

### APA Analysis Fails
**Possible causes:**
- .hic files not accessible
- Insufficient disk space for HDF5 files
- Java memory issues (strawr backend)

**Solutions:**
- Verify .hic file paths exist
- Check disk space (need ~50GB free)
- Skip APA with `--skip-apa` flag
- Run APA separately with fewer loops

---

## Performance Metrics

### Expected Runtime (on HPC with 16 CPUs, 128GB RAM)

| Step | 5kb | 10kb | 25kb | Notes |
|------|-----|------|------|-------|
| 1. prep_loops | 2-3 min | 1-2 min | 30-60 sec | I/O bound |
| 2. extract_counts | 15-30 min | 10-20 min | 5-10 min | Depends on .hic size |
| 3. aggregate | 2-4 min | 1-3 min | 30-90 sec | HDF5 read speed |
| 4. qc-val | 5-10 min | 3-5 min | 2-3 min | Plot generation |
| 5. edgeR | 1-2 min | 1 min | 30 sec | Model fitting |
| 6. compare_resolutions | 5-10 min | - | - | Cross-resolution |
| 7-8. Python/bash | 1 min | - | - | Fast |
| 9. downstream_analysis | 10-15 min | - | - | Merging + annotation |
| 10. visualizations | 20-30 min | - | - | Enrichment analysis |
| 11. apa_analysis | ~2 hours | ~1.5 hours | ~1 hour | Matrix extraction |

**Total:** ~4-6 hours without APA, ~20-24 hours with APA

### Memory Requirements

- **prep_loops:** 2-4 GB
- **extract_counts:** 8-16 GB (increases with # replicates)
- **aggregate:** 4-8 GB
- **edgeR:** 2-4 GB
- **apa_analysis:** 16-32 GB

### Disk Space Requirements

Per resolution:
- Intermediate files: 1-5 GB
- HDF5 matrices: 500 MB - 2 GB
- Final outputs: 100-500 MB

**Total:** ~20-30 GB for all three resolutions

---

## Statistical Methodology

### Experimental Design

```
Condition:  Control (WT)       BAP1-KO (Mut)
Replicates: n = 3              n = 3
Model:      ~group
Test:       Mutant - Control (coefficient 2)
Method:     Quasi-likelihood GLM
```

### Why Quasi-Likelihood GLM?

1. **Accounts for uncertainty in dispersion estimates**
2. **More rigorous error rate control** than likelihood ratio tests
3. **Robust to outliers** with `robust=TRUE`
4. **Recommended for experiments with small sample sizes** (n=3)

### Dispersion Estimation

- **Common dispersion**: Shared across all loops
- **Trended dispersion**: Varies with mean expression
- **Tagwise dispersion**: Loop-specific, shrunk toward trend
- **QL dispersion**: Additional uncertainty in variance

### Multiple Testing Correction

- **Method**: Benjamini-Hochberg FDR
- **Primary threshold**: FDR < 0.05
- **Stringent threshold**: FDR < 0.03 AND |logFC| > 0.3

---

## Biological Interpretation

### Fold Change Magnitude

| logFC | Linear FC | Biological Interpretation |
|-------|-----------|--------------------------|
| 0.5 | 1.41× | Modest change, may be secondary effect |
| 1.0 | 2.0× | Strong change, likely direct regulation |
| 2.0 | 4.0× | Very strong, major reorganization |
| 3.0 | 8.0× | Extreme, validate carefully |

### Loop Types

**Promoter-Promoter (P-P):**
- Gene regulation hubs
- Often house clusters of co-regulated genes
- TAD boundary associations

**Promoter-Enhancer (P-E):**
- Classic enhancer-promoter interactions
- Direct gene activation
- Most common differential type

**Enhancer-Enhancer (E-E):**
- Super-enhancer formation
- Phase-separated condensates
- Often tissue-specific

**Poised Elements:**
- H3K4me1+ but H3K27ac-
- Ready for activation
- Developmental regulators

### Distance Patterns

- **Short-range (< 100kb):** Immediate neighbors, local regulation
- **Medium-range (100kb - 1Mb):** Within TADs, typical E-P distance
- **Long-range (> 1Mb):** Inter-TAD, often structural/architectural

---

## Citation

If you use this pipeline, please cite:

**mariner:**
- Kramer et al. (2022). *mariner: Explore the Hi-Cs*. Bioconductor. https://bioconductor.org/packages/mariner/

**edgeR:**
- Robinson et al. (2010). *edgeR: a Bioconductor package for differential expression analysis of digital gene expression data*. Bioinformatics 26(1):139-140.
- Chen et al. (2016). *From reads to genes to pathways: differential expression analysis of RNA-Seq experiments using Rsubread and the edgeR quasi-likelihood pipeline*. F1000Research 5:1438.

**HiCCUPS:**
- Rao et al. (2014). *A 3D map of the human genome at kilobase resolution reveals principles of chromatin looping*. Cell 159(7):1665-1680.

---

## Support and Documentation

### Detailed Documentation

Each script has comprehensive documentation in `/docs`:

- [prep_loops.md](docs/prep_loops.md) - Loop merging and preparation
- [extract_counts.md](docs/extract_counts.md) - Hi-C matrix extraction
- [aggregate.md](docs/aggregate.md) - Matrix aggregation strategies
- [qc.md](docs/qc.md) - Quality control and validation
- [edgeR.md](docs/edgeR.md) - Differential analysis methodology
- [compare_resolutions.md](docs/compare_resolutions.md) - Multi-resolution comparison

---

## TAD and Compartment Analysis

Standalone scripts for analyzing TAD and chromatin compartment (PC1) differential data from HOMER's `getDiffExpression.pl`. These analyses complement the loop analysis by examining higher-order chromatin organization.

### TAD Volcano Plot

Generates publication-quality volcano plots for TAD inclusion ratio differential analysis.

**Script:** `scripts/tad_volcano_plot.R`

**Input:** Tab-delimited file from `getDiffExpression.pl` with TAD scores
**Default input:** `tad_analysis/Bap1.diff.tad.txt`

```bash
# Basic usage (generates both relaxed and standard threshold plots)
Rscript scripts/tad_volcano_plot.R

# Custom input file
Rscript scripts/tad_volcano_plot.R tad_analysis/Bap1.diff.tad.txt

# Custom output directory
Rscript scripts/tad_volcano_plot.R --output outputs/custom_tad/

# With custom title
Rscript scripts/tad_volcano_plot.R --title "BAP1 TAD Analysis"

# Custom plot dimensions
Rscript scripts/tad_volcano_plot.R --width 12 --height 10
```

**Thresholds generated (automatically produces both):**

| Version | FDR Cutoff | Difference Cutoff | Purpose |
|---------|------------|-------------------|---------|
| Relaxed | < 0.15 | \|Diff\| > 0.15 | Exploratory analysis |
| Standard | < 0.05 | \|Diff\| > 0.30 | Publication quality |

**Output files (`outputs/tad_analysis/`):**
- `tad_volcano_relaxed.pdf` - Volcano plot with relaxed thresholds
- `tad_volcano_standard.pdf` - Volcano plot with standard thresholds
- `tad_significant_relaxed.tsv` - Significant TADs (relaxed criteria)
- `tad_significant_standard.tsv` - Significant TADs (standard criteria)
- `tad_volcano_summary.txt` - Summary statistics for both thresholds
- `tad_all_annotated.tsv` - Full dataset with annotations

### Compartment (PC1) Volcano Plot

Generates volcano plots for chromatin compartment switching analysis based on PC1 values.

**Script:** `scripts/compartment_volcano_plot.R`

**Input:** Tab-delimited file from `getDiffExpression.pl` with PC1 values + gene annotations
**Default input:** `tad_analysis/diffcompartments.txt`

```bash
# Basic usage (generates both threshold versions)
Rscript scripts/compartment_volcano_plot.R

# Custom input file
Rscript scripts/compartment_volcano_plot.R tad_analysis/diffcompartments.txt

# Custom output directory
Rscript scripts/compartment_volcano_plot.R --output outputs/custom_compartment/

# With gene labeling on plot (labels top significant genes)
Rscript scripts/compartment_volcano_plot.R --label-genes --n-labels 20

# Custom plot dimensions
Rscript scripts/compartment_volcano_plot.R --width 14 --height 12
```

**Biological interpretation:**

| Difference | Compartment Shift | Biological Meaning |
|------------|-------------------|-------------------|
| Positive | B → A | More active chromatin in mutant |
| Negative | A → B | More inactive chromatin in mutant |

**Output files (`outputs/compartment_analysis/`):**
- `compartment_volcano_relaxed.pdf` - Volcano plot with relaxed thresholds
- `compartment_volcano_standard.pdf` - Volcano plot with standard thresholds
- `compartment_significant_relaxed.tsv` - Significant regions (relaxed)
- `compartment_significant_standard.tsv` - Significant regions (standard)
- `compartment_volcano_summary.txt` - Summary with gene annotations
- `compartment_all_annotated.tsv` - Full dataset with all annotations

### TAD Analysis Data Files

Located in `tad_analysis/`:

| File | Description |
|------|-------------|
| `Bap1.diff.tad.txt` | TAD differential expression from HOMER getDiffExpression.pl |
| `diffcompartments.txt` | Compartment (PC1) differential expression with gene annotations |
| `BAP1.tad.scores.txt` | Raw TAD inclusion ratio scores |
| `all_PC1.txt` | Raw PC1 values per sample |

---

### Additional Scripts

Scripts without dedicated docs (helper/utility):
- `annotate_loop_anchors.R` - ChIP-seq based classification
- `generate_final_results.py` - Stringent filtering
- `convert_final_bedpe.sh` - BEDPE format conversion
- `install_R_packages.R` - Dependency installation
- `tad_volcano_plot.R` - TAD differential volcano plots
- `compartment_volcano_plot.R` - Compartment (PC1) differential volcano plots

### Getting Help

1. **Check documentation**: Read relevant `/docs/*.md` file
2. **Review QC reports**: `outputs/res_{RES}kb/qc_report/`
3. **Check logs**: `outputs/edgeR_results_res_{RES}kb/logs/`
4. **Examine example outputs**: Compare to expected structure above

---

## License

This pipeline is provided for academic research use. Component packages have their own licenses:
- mariner: Artistic-2.0
- edgeR: GPL-3
- Bioconductor packages: Various (mostly Artistic-2.0 or GPL)

---

## Authors

**Pipeline Development:** Zakir Alibhai

**Core Methods:**
- mariner: Eric S. Davis, Douglas H. Phanstiel
- edgeR: Mark D. Robinson, Davis J. McCarthy, Gordon K. Smyth, Yunshun Chen
- HiCCUPS: Suhas S.P. Rao, Erez Lieberman Aiden

---

## Changelog

### Version 2.0 (Current)
- ✨ Multi-resolution support (5kb, 10kb, 25kb)
- ✨ Biological replicates (n=3 per condition)
- ✨ Quasi-likelihood GLM for robust statistics
- ✨ ChIP-seq based anchor classification
- ✨ Comprehensive QC with shift detection
- ✨ Publication-ready visualizations
- ✨ Full APA analysis

### Version 1.0 (Legacy)
- Single resolution (5kb)
- Merged samples (n=1 per condition)
- Fixed BCV dispersion
- Limited statistical power

---

**Last Updated:** November 10, 2025
