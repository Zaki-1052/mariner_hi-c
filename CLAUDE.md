# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a genomics/epigenomics research project. Primary languages: R for analysis scripts, Python for pipelines and visualization (Streamlit + Plotly). SLURM is used for HPC job submission. Always use Unix line endings (LF, not CRLF) when writing shell/SLURM scripts.

**Mariner Hi-C Differential Chromatin Loop Analysis Pipeline** - A bioinformatics pipeline for replicate-aware differential analysis of chromatin loops from Hi-C data, comparing BAP1-KO mutant vs wildtype control conditions in mouse (mm10 genome) with biological replication (n=3 per condition).

**Key Technologies:**
- R/Bioconductor: mariner, edgeR, InteractionSet, GenomicRanges, strawr
- Python 3.8+: pandas, numpy
- Bash scripting with SLURM HPC integration
- HDF5 for out-of-memory array storage

## Running the Pipeline

### Single-Command Execution (Recommended)

```bash
# On HPC with SLURM
sbatch scripts/run_full_pipeline.sb

# On local/server (direct execution)
bash scripts/run_full_pipeline.sb
```

**Expected runtime:** ~24 hours (includes time-intensive APA analysis)

### Step-by-Step Execution

For manual control or debugging, process each resolution individually:

```bash
# For each resolution (5000, 10000, 25000)
for RES in 5000 10000 25000; do
  Rscript scripts/prep_loops.R ${RES}
  Rscript scripts/extract_counts.R ${RES}
  Rscript scripts/aggregate.R ${RES}
  Rscript scripts/qc-val.R ${RES}
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

### Package Installation

```bash
# One-time setup
Rscript install_R_packages.R

# Verify conda environment (if using SLURM on HPC)
conda activate mariner_env
```

## Pipeline Architecture

### 11-Step Workflow

The pipeline follows a sequential workflow across three resolutions (5kb, 10kb, 25kb):

**Per-Resolution Analysis (Steps 1-5):**
1. **prep_loops.R** - Merge 6 replicate BEDPE files → consensus loop positions
2. **extract_counts.R** - Extract 5×5 Hi-C matrices from .hic files (6 samples)
3. **aggregate.R** - Sum matrices → count matrix (loops × samples) → edgeR DGEList
4. **qc-val.R** - Quality control validation & diagnostics (9 QC plots)
5. **edgeR.R** - Quasi-likelihood GLM differential analysis (n=3 per group)

**Multi-Resolution Post-Processing (Steps 6-11):**
6. **compare_resolutions.R** - Cross-resolution comparison & overlap analysis
7. **generate_final_results.py** - Apply stringent filters (|logFC|>0.3, FDR<0.03)
8. **convert_final_bedpe.sh** - Convert to BEDPE format for Juicebox visualization
9. **downstream_analysis.R** - Merge resolutions, annotate with genes & ChIP-seq
10. **visualizations.R** - Publication plots, enrichment analyses, volcano plots
11. **apa_analysis.R** - Aggregate Peak Analysis (APA) heatmaps

### Key Data Structures

**GInteractions** - Bioconductor class for paired genomic ranges (chromatin loops)
**DGEList** - edgeR object containing count matrix + sample metadata
**HDF5Array** - Out-of-memory arrays for large 5×5 contact matrices
**SummarizedExperiment** - Container for matrices + genomic coordinates

### Critical File Paths

**Input data (configured in individual scripts):**
- HiCCUPS loop calls: `/expanse/lustre/projects/csd940/ctea/nf-hic/juicer_frompre/hiccups_results/{sample}/postprocessed_pixels_{RES}.bedpe`
- Hi-C matrices: `/expanse/lustre/projects/csd940/ctea/nf-hic/juicerpre/{sample}.hic`
- ChIP-seq peaks: `peaks/beds/` directory (H3K27ac, H3K27me3, H3K4me1, H3K4me3, Bivalent per timepoint)

**Output structure:**
- Per-resolution intermediate: `outputs/res_{RES}kb/`
- edgeR results: `outputs/edgeR_results_res_{RES}kb/`
- Merged/final: `outputs/merged_loops/`, `outputs/visualizations/`, `outputs/apa_results/`

## Configuration

### Main Config File

`config/edgeR_config.yaml` - Controls statistical parameters, sample metadata, visualization settings

**Key parameters:**
- `statistics.fdr_primary: 0.05` - Primary FDR threshold
- `statistics.normalization_method: "TMM"` - Trimmed Mean of M-values
- `statistics.filtering.min_count: 5` - Minimum counts per loop
- `samples.names/groups` - 6 sample IDs and their condition assignments
- `replicate_qc.min_within_group_correlation: 0.90` - QC threshold

### Resolution-Specific Settings

Resolution is passed as command-line argument: `Rscript script.R 5000`

Resolutions supported: 5000 (5kb), 10000 (10kb), 25000 (25kb)

## Development Guidelines

When pre-computed summary files (TSVs, CSVs, log outputs) exist, use them directly. Do not recalculate from raw data unless the user explicitly asks for recomputation.

### Script Modification Workflow

When modifying existing scripts:

1. **Test on single resolution first** - Use 25kb for faster iteration
2. **Check intermediate outputs** - Verify RDS files in `outputs/res_{RES}kb/`
3. **Review QC plots** - Always check `outputs/res_{RES}kb/qc_report/`
4. **Compare with expected results** - See README.md "Expected Results" sections

### Adding New Analysis Steps

1. **Follow naming convention**: `{action}_{subject}.R` (e.g., `extract_counts.R`)
2. **Accept resolution as first argument**: `args <- commandArgs(trailingOnly=TRUE); RES <- as.numeric(args[1])`
3. **Use consistent output directories**: `outputs/res_{RES_KB}kb/` for intermediates
4. **Implement progress logging**: Use `cat()` with timestamps for long operations
5. **Save both RDS and TSV**: RDS for R objects, TSV for human-readable results
6. **Add error handling**: Check file existence, validate data integrity, use `tryCatch()`

### Statistical Analysis Pattern

The pipeline uses edgeR's quasi-likelihood GLM framework:

```r
# Standard pattern in edgeR.R
dge <- DGEList(counts = count_matrix, group = sample_groups)
keep <- filterByExpr(dge, min.count = 5)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")
design <- model.matrix(~group)
dge <- estimateDisp(dge, design, robust = TRUE)
fit <- glmQLFit(dge, design, robust = TRUE)
qlf <- glmQLFTest(fit, coef = 2)  # Test mutant effect
results <- topTags(qlf, n = Inf)
```

**Key principles:**
- **Data-driven dispersion**: Uses biological replicates (no fixed BCV)
- **Robust estimation**: `robust=TRUE` protects against outliers
- **Proper filtering**: `filterByExpr()` with replicate awareness
- **TMM normalization**: Accounts for sequencing depth differences

When the user asks about a specific script or file, confirm which file they mean before answering. Do not assume based on similar names. If the user corrects a misidentification, acknowledge immediately and switch context.

### Quality Control Checkpoints

**After aggregate.R:**
- Within-group correlation > 0.90 (preferably > 0.95)
- Between-group correlation 0.05-0.15 lower than within
- Library sizes within 2-fold
- Check: `outputs/res_{RES}kb/06_edgeR_input.rds` correlation matrix

**After edgeR.R:**
- BCV (biological coefficient of variation) ~0.2-0.5
- MDS plot shows clear condition separation
- Volcano plot has reasonable number of significant loops (500-5,000+ at FDR<0.05)
- Check: `outputs/edgeR_results_res_{RES}kb}/plots/` for diagnostic plots

**Red flags:**
- Within-group correlation < 0.85 → Potential replicate quality issue
- Very few significant loops (< 100) → Check BCV, filtering, or biological effect
- Very high significant loops (> 10,000) → Verify not batch effects

## Code Architecture Patterns

### 1. Resolution-Parameterized Design

All per-resolution scripts follow this pattern:

```r
# Parse resolution from command line
args <- commandArgs(trailingOnly = TRUE)
RES <- as.numeric(args[1])  # e.g., 5000
RES_KB <- RES / 1000        # e.g., 5

# Define resolution-specific paths
output_dir <- sprintf("outputs/res_%dkb", RES_KB)
input_file <- sprintf("path/to/file_%d.bedpe", RES)
```

### 2. HDF5-Backed Out-of-Memory Arrays

Large matrices use HDF5 to avoid memory issues:

```r
# In extract_counts.R
library(HDF5Array)
matrices <- pullHicMatrices(
  loops, hic_files,
  half = buffer_size
)  # Returns HDF5-backed DelayedArray

# Access without loading into memory
counts <- colSums(matrices, dims = 2)  # Lazy evaluation
```

### 3. Bioconductor GInteractions Workflow

Standard pattern for loop manipulation:

```r
library(InteractionSet)
library(GenomicRanges)

# Load BEDPE → GInteractions
loops <- read.table(bedpe_file)
gi <- GInteractions(
  GRanges(loops$chr1, IRanges(loops$start1, loops$end2)),
  GRanges(loops$chr2, IRanges(loops$start2, loops$end2))
)

# Merge nearby loops
merged <- mergePairs(gi, radius = 10000)

# Bin to resolution
binned <- binPairs(merged, binSize = RES)

# Add buffer for 5×5 extraction
buffered <- pixelsToMatrices(binned, buffer = 2)
```

### 4. Multi-Resolution Aggregation

downstream_analysis.R demonstrates merging across resolutions:

```r
# Load all resolutions
loops_5kb <- read.table("outputs/edgeR_results_res_5kb/primary_analysis/final_results.tsv")
loops_10kb <- read.table("outputs/edgeR_results_res_10kb/primary_analysis/final_results.tsv")
loops_25kb <- read.table("outputs/edgeR_results_res_25kb/primary_analysis/final_results.tsv")

# Convert to GInteractions
gi_list <- lapply(list(loops_5kb, loops_10kb, loops_25kb), bedpe_to_gi)

# Remove duplicates with 10kb tolerance
merged <- reduce(c(gi_list[[1]], gi_list[[2]], gi_list[[3]]),
                 min.gapwidth = 10000)
```

### 5. ChIP-seq Based Anchor Classification (7 Categories)

`annotate_loops_extended.R` implements chromatin state annotation using 5 histone
marks from timepoint-specific peak files in `peaks/beds/`:

```r
# Peak files loaded per timepoint (early/late) from PEAK_FILES config
# H3K27ac, H3K27me3, H3K4me1, H3K4me3, Bivalent (K4me3+K27me3)

# 7-category priority-based classification:
#   1. Active_Promoter:    H3K4me3+ AND NOT H3K27me3 AND <=2kb from TSS
#   2. Repressed_Promoter: H3K27me3+ AND NOT H3K27ac AND <=2kb from TSS
#   3. Bivalent_Promoter:  K4me3+K27me3 overlap (pre-computed intersection)
#   4. Polycomb:           H3K27me3+ AND >2kb from TSS (distal repressive)
#   5. Active_Enhancer:    H3K27ac+ AND >2kb from TSS
#   6. Poised_Enhancer:    H3K4me1+ AND NOT H3K27ac AND NOT H3K27me3 AND >2kb
#   7. Other:              No ChIP-seq marks / structural elements

# Apply to both anchors → up to 28 loop type combinations
# (Active_Promoter-Active_Enhancer, Polycomb-Polycomb, etc.)
```

## Common Tasks

### Add New QC Plot

Edit `scripts/qc-val.R`:

```r
# Add after existing QC analyses
cat("\n=== New QC Analysis ===\n")

# Compute metric
new_metric <- compute_metric(count_matrix)

# Generate plot
pdf(file.path(qc_dir, "new_qc_plot.pdf"))
plot(new_metric, main = "New QC Metric")
dev.off()

# Add to summary
qc_summary$new_metric <- new_metric
saveRDS(qc_summary, file.path(qc_dir, "qc_report_summary.rds"))
```

### Modify Statistical Thresholds

Edit `config/edgeR_config.yaml`:

```yaml
statistics:
  fdr_primary: 0.03        # More stringent
  filtering:
    min_count: 10          # Higher filtering
```

Then re-run edgeR.R: `Rscript scripts/edgeR.R 5000`

### Extract Specific Loop Subsets

```r
# Load edgeR results
results <- read.table("outputs/edgeR_results_res_5kb/primary_analysis/all_results_primary.tsv",
                      header = TRUE)

# Filter by criteria
strong_up <- results[results$FDR < 0.05 & results$logFC > 1, ]
chr1_loops <- results[results$chr1 == "chr1" & results$chr2 == "chr1", ]

# Export
write.table(strong_up, "custom_subset.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
```

### Visualize Loops in Juicebox

1. Convert results to BEDPE (if not already done):
   ```bash
   bash scripts/convert_final_bedpe.sh
   ```

2. Load in Juicebox:
   - Open .hic file for sample
   - Load → 2D Annotations → Select `outputs/bedpe_final/5kb_final.bedpe`
   - Color by logFC or FDR column

### Debug Failed Pipeline Step

```bash
# Check exit code and logs
echo $?  # Non-zero = failure

# Review SLURM log
tail -100 logs/mariner_complete_*.out

# Check intermediate outputs exist
ls -lh outputs/res_5kb/
ls -lh outputs/edgeR_results_res_5kb/

# Re-run failed step manually with verbose output
Rscript scripts/extract_counts.R 5000 2>&1 | tee debug.log
```

## Important Files to Review

**Understanding results:**
- `outputs/edgeR_results_res_{RES}kb/primary_analysis/summary_statistics.txt` - Analysis summary
- `outputs/edgeR_results_res_{RES}kb/plots/volcano_plot_primary.pdf` - Visual overview
- `outputs/merged_loops/characterized_loops.tsv` - Annotated differential loops

**Quality control:**
- `outputs/res_{RES}kb/qc_report/` - All 9 QC diagnostic plots
- `outputs/edgeR_results_res_{RES}kb/plots/mds_plot.pdf` - Sample relationships
- `outputs/edgeR_results_res_{RES}kb/plots/bcv_plot.pdf` - Dispersion estimates

**Biological interpretation:**
- `outputs/visualizations/enrichment/` - GO/KEGG enrichment results
- `outputs/visualizations/loop_classification/` - Loop type distributions
- `outputs/apa_results/res_{RES}kb/` - Aggregate peak analysis heatmaps

## Key Dependencies & Versions

**Required R packages:**
- mariner (Hi-C loop analysis framework)
- edgeR ≥4.0 (differential analysis with quasi-likelihood GLM)
- InteractionSet (genomic interaction data structures)
- strawr (read .hic files via Straw API)
- HDF5Array (out-of-memory array support)
- ChIPseeker, clusterProfiler (functional enrichment)
- TxDb.Mmusculus.UCSC.mm10.knownGene (mouse gene annotations)

**System requirements:**
- R ≥ 4.3.0 (for Bioconductor 3.18+)
- Python ≥ 3.8
- HDF5 library (for HDF5Array support)
- 128GB RAM recommended (HPC), 32GB minimum (local with HDF5 backend)

See `R_PACKAGES.md` for complete list and installation instructions.

## Biological Context

**Chromatin loops** are 3D interactions between distant genomic regions that regulate gene expression by bringing enhancers into contact with promoters.

**This pipeline identifies:**
- Differential loops between BAP1-KO and wildtype
- Loop types based on 7 anchor categories (Active_Promoter, Repressed_Promoter, Bivalent_Promoter, Polycomb, Active_Enhancer, Poised_Enhancer, Other)
- Enrichment at gene regulatory elements (H3K27ac, H3K27me3, H3K4me1, H3K4me3, Bivalent marks)
- Functional enrichment (GO/KEGG pathways)

**Key effect sizes:**
- |logFC| > 1.0 (2-fold) = Strong change, likely direct regulation
- |logFC| > 0.5 (1.4-fold) = Moderate change
- |logFC| > 0.3 (1.23-fold) = Stringent threshold for final results

## Performance Notes

**Expected runtimes (16 CPUs, 128GB RAM):**
- prep_loops.R: 1-3 min per resolution
- extract_counts.R: 10-30 min per resolution (depends on .hic file size)
- aggregate.R: 1-4 min per resolution
- qc-val.R: 3-10 min per resolution
- edgeR.R: 1-2 min per resolution
- apa_analysis.R: ~2 hours (time-intensive, matrix extraction)

**Memory usage:**
- Most steps: 4-16 GB
- extract_counts.R: 8-32 GB (scales with # samples)
- apa_analysis.R: 16-64 GB (large matrix aggregation)

**Disk space per resolution:**
- Intermediate files: 1-5 GB
- HDF5 matrices: 500 MB - 2 GB
- Final outputs: 100-500 MB
- **Total for all 3 resolutions:** ~20-30 GB

## Troubleshooting Common Issues

**"Cannot open .hic file" error:**
- Verify paths in extract_counts.R and apa_analysis.R
- Check file permissions and network access (if on remote storage)
- Ensure strawr package installed correctly with Java backend

**"HDF5 library not found":**
```bash
# macOS
brew install hdf5
# Ubuntu
sudo apt-get install libhdf5-dev
# Then reinstall
BiocManager::install("HDF5Array", force = TRUE)
```

**Low within-group correlation (<0.90):**
- Check individual replicate quality in QC plots
- Review raw .hic file statistics
- Consider removing outlier replicate from analysis
- May indicate high biological variability

**Very few significant loops (<100):**
- Check BCV plot - high BCV (>0.5) reduces power
- Review filtering stringency in config
- Verify biological replicates (not technical replicates)
- May indicate weak biological effect

## TAD and Compartment Analysis Scripts

Standalone scripts for analyzing TAD and compartment (PC1) differential data from HOMER's `getDiffExpression.pl`.

### TAD Volcano Plot (`scripts/tad_volcano_plot.R`)

Generates volcano plots for TAD inclusion ratio differential analysis.

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
```

**Thresholds generated:**
- Relaxed: FDR < 0.15, |Difference| > 0.15 (exploratory)
- Standard: FDR < 0.05, |Difference| > 0.30 (publication)

**Output files (`outputs/tad_analysis/`):**
- `tad_volcano_relaxed.pdf` - Volcano plot with relaxed thresholds
- `tad_volcano_standard.pdf` - Volcano plot with standard thresholds
- `tad_significant_relaxed.tsv` - Significant TADs (relaxed)
- `tad_significant_standard.tsv` - Significant TADs (standard)
- `tad_volcano_summary.txt` - Summary statistics
- `tad_all_annotated.tsv` - Full dataset

### Compartment Volcano Plot (`scripts/compartment_volcano_plot.R`)

Generates volcano plots for PC1/compartment differential analysis.

**Input:** Tab-delimited file from `getDiffExpression.pl` with PC1 values
**Default input:** `tad_analysis/diffcompartments.txt`

```bash
# Basic usage (generates both threshold versions)
Rscript scripts/compartment_volcano_plot.R

# Custom input file
Rscript scripts/compartment_volcano_plot.R tad_analysis/diffcompartments.txt

# Custom output directory
Rscript scripts/compartment_volcano_plot.R --output outputs/custom_compartment/

# With gene labeling on plot
Rscript scripts/compartment_volcano_plot.R --label-genes --n-labels 20
```

**Biological interpretation:**
- Positive Difference = shift toward A compartment (more active) in mutant
- Negative Difference = shift toward B compartment (more inactive) in mutant

**Output files (`outputs/compartment_analysis/`):**
- `compartment_volcano_relaxed.pdf` - Volcano plot with relaxed thresholds
- `compartment_volcano_standard.pdf` - Volcano plot with standard thresholds
- `compartment_significant_relaxed.tsv` - Significant regions (relaxed)
- `compartment_significant_standard.tsv` - Significant regions (standard)
- `compartment_volcano_summary.txt` - Summary with gene annotations
- `compartment_all_annotated.tsv` - Full dataset

### TAD Analysis Data Files

Located in `tad_analysis/`:
- `Bap1.diff.tad.txt` - TAD differential expression from HOMER
- `diffcompartments.txt` - Compartment (PC1) differential expression
- `BAP1.tad.scores.txt` - Raw TAD scores
- `all_PC1.txt` - Raw PC1 values

## Additional Documentation

Detailed documentation for each major script in `/docs`:
- `prep_loops.md` - Loop merging methodology
- `extract_counts.md` - Hi-C matrix extraction details
- `aggregate.md` - Aggregation strategies (sum vs weighted)
- `qc.md` - All 9 QC analyses explained
- `edgeR.md` - Statistical framework and interpretation
- `compare_resolutions.md` - Multi-resolution comparison approach
- `apa_analysis.md` - Aggregate Peak Analysis methodology
- `downstream_analysis.md` - Merging and annotation workflow
- `visualizations.md` - Publication plot generation

## Citation

**Core methods:**
- mariner: Kramer et al. (2022) mariner: Explore the Hi-Cs. Bioconductor.
- edgeR: Robinson et al. (2010) Bioinformatics 26(1):139-140; Chen et al. (2016) F1000Research 5:1438
- HiCCUPS: Rao et al. (2014) Cell 159(7):1665-1680
