# Comprehensive Quality Control and Validation (qc-val.R)

## Biological Purpose

This script performs systematic quality control and validation of the entire Mariner pipeline, from loop preparation through aggregation. It generates comprehensive diagnostics across **6 biological replicates** (3 control + 3 mutant) to ensure data quality before differential analysis.

**Key QC areas:**
- **Data integrity**: Checking for corrupted or malformed data
- **Sample quality**: Assessing replicate consistency and technical reproducibility
- **Loop-level patterns**: Validating biological signal vs. noise
- **Spatial accuracy**: Detecting bin-shift artifacts and positional errors
- **Aggregation validation**: Comparing different counting strategies

## Multi-Resolution Support

The script accepts resolution as a command-line argument:
```bash
Rscript scripts/qc-val.R 5000   # 5kb resolution
Rscript scripts/qc-val.R 10000  # 10kb resolution
Rscript scripts/qc-val.R 25000  # 25kb resolution
```

## Input Data

### Pipeline Outputs (All Resolution-Specific)

**From `outputs/res_{RES}kb/`:**
- `01_ginteractions.rds` - Raw loop calls per replicate
- `02_merged.rds` - Consensus merged loops
- `03_binned.rds` - Resolution-aligned loops
- `04_buffered.rds` - 5×5 pixel buffers
- `05_extracted/` - HDF5-backed extraction matrices
- `05_metadata.rds` - Extraction metadata and correlations
- `06_counts_matrix.rds` - Aggregated count matrix (N_loops × 6)
- `06_edgeR_input.rds` - DGEList object
- `06_all_strategies.rds` - Alternative aggregation methods

## Output Directory

All QC outputs are saved to resolution-specific subdirectories:
```
outputs/res_5kb/qc_report/
outputs/res_10kb/qc_report/
outputs/res_25kb/qc_report/
```

## Section 1: Data Loading & Basic Integrity Checks

### File Validation

```r
# Verify all pipeline outputs exist and load successfully
tryCatch({
  gi_list <- readRDS(file.path(input_dir, "01_ginteractions.rds"))
  merged <- readRDS(file.path(input_dir, "02_merged.rds"))
  binned <- readRDS(file.path(input_dir, "03_binned.rds"))
  buffered <- readRDS(file.path(input_dir, "04_buffered.rds"))
  metadata <- readRDS(file.path(input_dir, "05_metadata.rds"))
  counts_matrix <- readRDS(file.path(input_dir, "06_counts_matrix.rds"))
  edger_obj <- readRDS(file.path(input_dir, "06_edgeR_input.rds"))
  pixels <- loadHDF5SummarizedExperiment(dir = input_dir,
                                         prefix = "05_extracted")
}, error = function(e) {
  stop("Cannot proceed with QC - data loading failed")
})
```

### Dimension Verification

```r
cat("Basic Dimensions:\n")
cat(sprintf("  Initial loops (ctrl_M1): %d\n", length(gi_list$ctrl_M1)))
cat(sprintf("  Initial loops (ctrl_M2): %d\n", length(gi_list$ctrl_M2)))
cat(sprintf("  Initial loops (ctrl_M3): %d\n", length(gi_list$ctrl_M3)))
cat(sprintf("  Initial loops (mut_M1):  %d\n", length(gi_list$mut_M1)))
cat(sprintf("  Initial loops (mut_M2):  %d\n", length(gi_list$mut_M2)))
cat(sprintf("  Initial loops (mut_M3):  %d\n", length(gi_list$mut_M3)))
cat(sprintf("  Merged loops:         %d\n", length(merged)))
cat(sprintf("  Count matrix:         %d loops × 6 samples\n",
            nrow(counts_matrix)))
```

**Expected patterns:**
- Individual replicate inputs: 15,000-25,000 loops each
- Merged consensus: 20,000-30,000 total positions
- Final count matrix: Same as merged (no loops lost in aggregation)

### Data Integrity Checks

```r
# NA values (should be minimal)
na_count <- sum(is.na(counts_matrix))
na_pct <- 100 * na_count / length(counts_matrix)
cat(sprintf("  NA values: %d (%.2f%%)\n", na_count, na_pct))

# Zero values (expected for sparse Hi-C data)
zero_count <- sum(counts_matrix == 0)
zero_pct <- 100 * zero_count / length(counts_matrix)
cat(sprintf("  Zero values: %d (%.2f%%)\n", zero_count, zero_pct))

# Negative values (MUST be zero!)
neg_count <- sum(counts_matrix < 0, na.rm = TRUE)
if (neg_count > 0) {
  cat(sprintf("  ✗ WARNING: %d negative values detected!\n", neg_count))
} else {
  cat("  ✓ No negative values detected\n")
}
```

**Quality indicators:**
- **NA < 5%**: Good data quality
- **Zeros 10-40%**: Normal Hi-C sparsity
- **Negatives = 0**: REQUIRED (data corruption if > 0)

## Section 2: Sample-Level Quality Control

### Replicate-Aware Statistics

```r
# Count statistics per replicate
sample_stats <- data.frame(
  Sample = colnames(counts_matrix),
  Min = apply(counts_matrix, 2, min, na.rm = TRUE),
  Q1 = apply(counts_matrix, 2, quantile, probs = 0.25, na.rm = TRUE),
  Median = apply(counts_matrix, 2, median, na.rm = TRUE),
  Mean = apply(counts_matrix, 2, mean, na.rm = TRUE),
  Q3 = apply(counts_matrix, 2, quantile, probs = 0.75, na.rm = TRUE),
  Max = apply(counts_matrix, 2, max, na.rm = TRUE),
  SD = apply(counts_matrix, 2, sd, na.rm = TRUE)
)
print(sample_stats)
```

**Example output (6 replicates):**
```
         Sample    Min    Q1 Median   Mean    Q3    Max     SD
1      ctrl_M1    0.0  12.3   45.2   67.8  89.1  567.3  78.4
2      ctrl_M2    0.0  11.8   43.7   65.2  86.3  543.1  76.1
3      ctrl_M3    0.0  12.7   46.9   69.1  91.4  589.2  80.3
4       mut_M1    0.0  13.1   48.5   71.3  94.2  612.7  82.6
5       mut_M2    0.0  12.5   46.8   68.9  90.7  587.4  79.8
6       mut_M3    0.0  13.4   49.2   72.1  95.8  625.3  84.1
```

### Library Size Analysis

```r
lib_sizes <- colSums(counts_matrix, na.rm = TRUE)
cat("Library Sizes (Total Counts):\n")
print(lib_sizes)

# For 6 replicates, analyze by group
ctrl_libs <- lib_sizes[1:3]
mut_libs <- lib_sizes[4:6]

cat(sprintf("\n  Control group:\n"))
cat(sprintf("    Mean: %.0f, SD: %.0f, CV: %.2f\n",
            mean(ctrl_libs), sd(ctrl_libs), sd(ctrl_libs)/mean(ctrl_libs)))
cat(sprintf("  Mutant group:\n"))
cat(sprintf("    Mean: %.0f, SD: %.0f, CV: %.2f\n",
            mean(mut_libs), sd(mut_libs), sd(mut_libs)/mean(mut_libs)))
cat(sprintf("  Ratio (mut/ctrl medians): %.3f\n",
            median(mut_libs)/median(ctrl_libs)))
```

**Quality indicators:**
- **CV < 0.20**: Excellent within-group consistency
- **CV 0.20-0.30**: Acceptable variability
- **Ratio 0.8-1.2**: Balanced between conditions
- **CV > 0.30**: Investigate outlier replicates

### Full 6×6 Correlation Matrix

```r
cor_matrix <- cor(counts_matrix, use = "complete.obs", method = "pearson")

cat("\nFull correlation matrix:\n")
print(round(cor_matrix, 3))

# Within-group correlations
ctrl_cors <- cor_matrix[1:3, 1:3][upper.tri(cor_matrix[1:3, 1:3])]
mut_cors <- cor_matrix[4:6, 4:6][upper.tri(cor_matrix[4:6, 4:6])]

cat(sprintf("  Control replicates: %.3f ± %.3f (range: %.3f-%.3f)\n",
            mean(ctrl_cors), sd(ctrl_cors), min(ctrl_cors), max(ctrl_cors)))
cat(sprintf("  Mutant replicates:  %.3f ± %.3f (range: %.3f-%.3f)\n",
            mean(mut_cors), sd(mut_cors), min(mut_cors), max(mut_cors)))

# Between-group correlations
between_cors <- as.vector(cor_matrix[1:3, 4:6])
cat(sprintf("  Between groups:     %.3f ± %.3f (range: %.3f-%.3f)\n",
            mean(between_cors), sd(between_cors),
            min(between_cors), max(between_cors)))
```

**Expected correlation patterns:**
```
           ctrl_M1  ctrl_M2  ctrl_M3  mut_M1  mut_M2  mut_M3
ctrl_M1      1.000    0.950    0.945   0.880   0.875   0.870
ctrl_M2      0.950    1.000    0.948   0.882   0.878   0.873
ctrl_M3      0.945    0.948    1.000   0.878   0.874   0.869
mut_M1       0.880    0.882    0.878   1.000   0.952   0.947
mut_M2       0.875    0.878    0.874   0.952   1.000   0.949
mut_M3       0.870    0.873    0.869   0.947   0.949   1.000
```

**Quality assessment:**
```r
if (mean(ctrl_cors) > 0.95 && mean(mut_cors) > 0.95) {
  cat("  ✓ Excellent within-group reproducibility (>0.95)\n")
} else if (mean(ctrl_cors) > 0.90 && mean(mut_cors) > 0.90) {
  cat("  ✓ Good within-group reproducibility (>0.90)\n")
} else {
  cat("  ⚠ Lower than expected within-group correlation\n")
}

if ((mean(ctrl_cors) - mean(between_cors)) > 0.05 ||
    (mean(mut_cors) - mean(between_cors)) > 0.05) {
  cat("  ✓ Clear biological signal (within > between)\n")
} else {
  cat("  ⚠ Weak biological differences\n")
}
```

### Visualization: Correlation Heatmap

**Output**: `01_sample_correlation.pdf`

- Displays 6×6 correlation matrix
- Color-coded by viridis palette
- Numerical values overlaid
- Clustering shows group structure

### Visualization: Distribution Comparison

**Output**: `03_count_distributions.pdf` (4-panel plot)

**Panel 1: Density plots**
- Overlapping density curves for all 6 samples
- Log10-transformed counts
- Shows distribution shape and overlap

**Panel 2: Box plots**
- Side-by-side comparison
- Identifies outliers and median differences
- Good for spotting replicate outliers

**Panel 3: Violin plots**
- Distribution shape + summary statistics
- Width shows density at each value
- Inner boxplot for quartiles

**Panel 4: Cumulative distribution**
- ECDF curves for each sample
- Shows what percentage of loops have ≤ X counts
- Useful for understanding data sparsity

## Section 3: Loop-Level Quality Assessment

### MA Plot Analysis

```r
# Calculate MA values using group averages
ctrl_avg <- rowMeans(counts_matrix[, 1:3])
mut_avg <- rowMeans(counts_matrix[, 4:6])

A <- 0.5 * (log2(ctrl_avg + 1) + log2(mut_avg + 1))
M <- log2(mut_avg + 1) - log2(ctrl_avg + 1)

cat("MA Statistics:\n")
cat(sprintf("  Average expression (A): %.2f ± %.2f\n",
            mean(A), sd(A)))
cat(sprintf("  Log2 fold change (M):   %.3f ± %.3f\n",
            mean(M), sd(M)))
cat(sprintf("  |M| > 1 (2-fold change): %d loops (%.1f%%)\n",
            sum(abs(M) > 1), 100 * sum(abs(M) > 1) / length(M)))
cat(sprintf("  |M| > 2 (4-fold change): %d loops (%.1f%%)\n",
            sum(abs(M) > 2), 100 * sum(abs(M) > 2) / length(M)))
```

**MA plot interpretation:**
- **X-axis (A)**: Average log2 expression (loop strength)
- **Y-axis (M)**: Log2 fold change between conditions
- **Centered at M=0**: No systematic bias
- **Symmetric distribution**: Balanced up/down regulation
- **|M| > 1 lines**: 2-fold change thresholds (blue dashed)
- **|M| > 2 lines**: 4-fold change thresholds (red dotted)

**Quality indicators:**
- **Mean M near 0**: No global shift
- **Symmetric scatter**: Equal up/down potential
- **2-10% with |M|>1**: Typical differential signal
- **0.5-5% with |M|>2**: Strong differential signal

### Visualization: MA Plot

**Output**: `04_ma_plot.pdf`

- Points colored by change magnitude
- Gray: |M| ≤ 1 (small change)
- Orange: |M| > 1 (moderate change, 2-fold)
- Red: |M| > 2 (large change, 4-fold)
- Reference lines at M = 0, ±1, ±2

### Loop Variability Analysis

```r
# Calculate coefficient of variation per loop
loop_means <- rowMeans(counts_matrix, na.rm = TRUE)
loop_sds <- apply(counts_matrix, 1, sd, na.rm = TRUE)
loop_cv <- loop_sds / loop_means

cat(sprintf("  Coefficient of Variation: %.3f ± %.3f\n",
            mean(loop_cv, na.rm = TRUE), sd(loop_cv, na.rm = TRUE)))
cat(sprintf("  High variability loops (CV > 1): %d (%.1f%%)\n",
            sum(loop_cv > 1, na.rm = TRUE),
            100 * sum(loop_cv > 1, na.rm = TRUE) / length(loop_cv)))
```

**CV interpretation:**
- **CV < 0.5**: Low variability, consistent signal
- **CV 0.5-1.0**: Moderate variability
- **CV > 1.0**: High variability (often low-count loops)
- **Very high CV**: May be filtered in edgeR

### Visualization: Variable Loops Heatmap

**Output**: `05_variable_loops_heatmap.pdf`

- Shows top 50 most variable loops
- Z-score scaled within each loop
- Hierarchical clustering of loops and samples
- Useful for identifying outliers or batch effects

## Section 4: Spatial & Distance Analysis

### Chromosome Distribution

```r
chr_summary <- loop_coords %>%
  group_by(chr) %>%
  summarise(
    n_loops = n(),
    mean_ctrl = mean(ctrl_count, na.rm = TRUE),
    mean_mut = mean(mut_count, na.rm = TRUE),
    mean_total = mean(total_count, na.rm = TRUE),
    mean_fc = mean(log2FC, na.rm = TRUE)
  ) %>%
  arrange(desc(n_loops))

print(head(chr_summary, 10))
```

**Expected patterns:**
- Larger chromosomes (chr1-5) have more loops
- Similar distribution between ctrl and mut
- Any chromosome-specific biases should be noted

### Genomic Distance Statistics

```r
cat("Genomic Distance Statistics (intrachromosomal loops):\n")
valid_distances <- loop_coords$distance[!is.na(loop_coords$distance)]
cat(sprintf("  Loops with distance: %d / %d (%.1f%%)\n",
            length(valid_distances), nrow(loop_coords),
            100 * length(valid_distances) / nrow(loop_coords)))
cat(sprintf("  Median distance: %s\n", comma(median(valid_distances))))
cat(sprintf("  Mean distance: %s\n", comma(mean(valid_distances))))
```

**Typical distance ranges:**
- **5kb resolution**: Median ~200-500 kb
- **10kb resolution**: Median ~300-600 kb
- **25kb resolution**: Median ~500 kb - 1 Mb

### Visualization: Distance Decay

**Output**: `07_distance_decay.pdf` (4-panel plot)

**Panel 1**: Total signal vs distance (ctrl + mut combined)
**Panel 2**: Control signal vs distance
**Panel 3**: Mutant signal vs distance
**Panel 4**: Log2 fold change vs distance

**Expected patterns:**
- Signal decreases with distance (power law decay)
- Similar decay curves for ctrl and mut
- FC variation independent of distance (centered at 0)
- Deviations may indicate distance-specific effects

## Section 5: Matrix-Level Buffer Validation

### 5×5 Matrix Inspection

```r
# Get HDF5-backed array
count_array <- counts(pixels)
array_dims <- dim(count_array)
cat(sprintf("  Array dimensions: %d x %d x %d x %d\n",
            array_dims[1], array_dims[2], array_dims[3], array_dims[4]))
# Expected: 5 x 5 x N_loops x 6
```

### Per-Loop Matrix Analysis

```r
analyze_loop_matrix <- function(loop_id, count_array, n_samples) {
  # For 6 replicates, average within groups
  ctrl_mat <- apply(count_array[, , loop_id, 1:3], c(1,2), mean)
  mut_mat <- apply(count_array[, , loop_id, 4:6], c(1,2), mean)

  # Detect peak position
  ctrl_peak_pos <- which(ctrl_mat == max(ctrl_mat), arr.ind = TRUE)[1, ]
  mut_peak_pos <- which(mut_mat == max(mut_mat), arr.ind = TRUE)[1, ]

  # Check if centered at (3,3)
  ctrl_centered <- all(ctrl_peak_pos == c(3, 3))
  mut_centered <- all(mut_peak_pos == c(3, 3))
  shift_detected <- !identical(ctrl_peak_pos, mut_peak_pos)

  # Center enrichment (fraction of total in center pixel)
  ctrl_center_enrichment <- ctrl_mat[3, 3] / sum(ctrl_mat)
  mut_center_enrichment <- mut_mat[3, 3] / sum(mut_mat)

  return(list(
    ctrl_centered = ctrl_centered,
    mut_centered = mut_centered,
    shift_detected = shift_detected,
    ctrl_center_enrichment = ctrl_center_enrichment,
    mut_center_enrichment = mut_center_enrichment,
    ctrl_peak_pos = ctrl_peak_pos,
    mut_peak_pos = mut_peak_pos
  ))
}
```

### Matrix Summary Statistics

```r
cat("Matrix-Level Summary:\n")
cat(sprintf("  Control peaks centered: %d / %d (%.1f%%)\n",
            centered_ctrl, n_loops, 100 * centered_ctrl / n_loops))
cat(sprintf("  Mutant peaks centered:  %d / %d (%.1f%%)\n",
            centered_mut, n_loops, 100 * centered_mut / n_loops))
cat(sprintf("  Positional shifts detected: %d / %d (%.1f%%)\n",
            shifts, n_loops, 100 * shifts / n_loops))
cat(sprintf("  Mean center enrichment (ctrl): %.3f ± %.3f\n",
            mean_ctrl_enrich, sd_ctrl_enrich))
cat(sprintf("  Mean center enrichment (mut):  %.3f ± %.3f\n",
            mean_mut_enrich, sd_mut_enrich))
```

**Quality indicators:**
- **60-80% centered**: Good loop calling precision
- **10-30% shifts detected**: Normal biological/technical variation
- **Center enrichment 0.20-0.40**: Appropriate signal concentration
- **< 50% centered**: May need larger buffer or different binning

### Shift Status Export for edgeR

**CRITICAL**: The script extracts per-loop shift status for downstream enrichment analysis:

```r
# Extract boolean shift status for each loop
loop_shift_status <- sapply(matrix_stats, function(x) x$shift_detected)

# Save for edgeR to use in enrichment testing
saveRDS(loop_shift_status, file.path(qc_dir, "loop_shift_status.rds"))

# Also save summary table
shift_summary_df <- data.frame(
  loop_id = paste0("loop_", 1:n_loops),
  shift_detected = loop_shift_status,
  ctrl_peak_row = sapply(matrix_stats, function(x) x$ctrl_peak_pos[1]),
  ctrl_peak_col = sapply(matrix_stats, function(x) x$ctrl_peak_pos[2]),
  mut_peak_row = sapply(matrix_stats, function(x) x$mut_peak_pos[1]),
  mut_peak_col = sapply(matrix_stats, function(x) x$mut_peak_pos[2])
)
write.table(shift_summary_df, "loop_shift_summary.tsv")
```

**Use in edgeR**:
- Loaded by edgeR script to test if shifted loops are enriched in significant set
- Fisher's exact test for shift enrichment
- Informs whether buffer approach is capturing real biological signal

### Visualization: Example Matrices

**Output**: `08_example_matrices.pdf`

- Shows top 6 loops by total signal
- Side-by-side ctrl vs mut for each loop
- 5×5 heatmaps with numerical values overlaid
- Viridis color scheme
- Validates extraction worked correctly

## Section 6: Aggregation Strategy Comparison

### Strategy Correlation Analysis

```r
if (!is.null(all_strategies)) {
  sum_counts <- all_strategies$counts$sum
  weighted_counts <- all_strategies$counts$weighted

  # Correlate strategies within each sample
  cat("Correlation between aggregation strategies (Control):\n")
  cor_strategies_ctrl <- cor(sum_counts[,1:3], weighted_counts[,1:3])
  print(round(cor_strategies_ctrl, 3))

  cat("Correlation between aggregation strategies (Mutant):\n")
  cor_strategies_mut <- cor(sum_counts[,4:6], weighted_counts[,4:6])
  print(round(cor_strategies_mut, 3))
}
```

**Expected correlation:**
- **r > 0.95**: Strategies highly concordant
- **r 0.90-0.95**: Some differences but generally similar
- **r < 0.90**: Strategies capture different aspects of data

### Strategy Interpretation

**Sum aggregation:**
- Simple sum of all 25 pixels
- Robust to positional shifts
- May include more background
- **Best for**: Shift-tolerant analysis

**Weighted aggregation:**
- Center pixels weighted higher
- More sensitive to well-positioned loops
- Reduces background noise
- **Best for**: Precision-focused analysis

### Visualization: Strategy Comparison

**Output**: `09_strategy_comparison.pdf` (4-panel plot)

**Panel 1**: Control - sum vs weighted scatter
**Panel 2**: Mutant - sum vs weighted scatter
**Panel 3**: Control - ratio distribution (weighted/sum)
**Panel 4**: Mutant - ratio distribution (weighted/sum)

**Interpretation:**
- Points near diagonal (y=x) → strategies agree
- Ratio near 1.0 → similar magnitude
- High correlation → choice has minimal impact on results

## Section 7: Final QC Assessment & Recommendations

### Automated Pass/Fail Criteria

```r
qc_pass <- TRUE
critical_issues <- c()
warnings <- c()

# Check 1: Sample correlation
if (cor_value < 0.70) {
  qc_pass <- FALSE
  critical_issues <- c(critical_issues, "Low sample correlation (<0.70)")
} else if (cor_value < 0.85) {
  warnings <- c(warnings, "Moderate sample correlation (0.70-0.85)")
}

# Check 2: Library size ratio
if (lib_ratio > 2.0) {
  warnings <- c(warnings, "Large library size difference (>2-fold)")
}

# Check 3: Data integrity
if (neg_count > 0) {
  qc_pass <- FALSE
  critical_issues <- c(critical_issues,
                      "Negative values detected in count matrix")
}

# Check 4: Matrix centering
if (centered_ctrl_pct < 50 || centered_mut_pct < 50) {
  warnings <- c(warnings, "Less than 50% of loops have centered peaks")
}

# Check 5: Center enrichment
if (mean_center_enrichment < 0.2) {
  warnings <- c(warnings,
              "Low center enrichment (<0.2) suggests diffuse signal")
}
```

### QC Status Reporting

**PASS - All checks passed:**
```
✓✓✓ PASS - All quality checks passed! ✓✓✓

Your data is high quality and ready for differential analysis.
```

**PASS with warnings:**
```
✓ PASS (with warnings)

Data quality is acceptable, but note the following:
  ⚠ Moderate sample correlation (0.70-0.85)
  ⚠ Less than 50% of loops have centered peaks
```

**FAIL - Critical issues:**
```
✗ FAIL - Critical issues detected

  ✗ Low sample correlation (<0.70)
  ✗ Negative values detected in count matrix

Additional warnings:
  ⚠ Large library size difference (>2-fold)
```

### Recommendations Based on Results

**If QC passes:**
```
RECOMMENDATIONS:
1. ✓ Proceed with edgeR differential analysis
2. ✓ Data quality supports biological interpretation
3. ⚠ 15.3% of loops show positional shifts - buffer approach validated

NEXT STEPS:
→ Run edgeR differential analysis script
→ Compare results with HiCCUPS differential calls
→ Validate top hits with biological context
```

**If QC fails:**
```
RECOMMENDATIONS:
1. ✗ Investigate critical issues before proceeding
2. ✗ Consider re-running pipeline with adjusted parameters
3. ✗ Check input .hic file quality

TROUBLESHOOTING:
→ Low correlation may indicate:
  - Batch effects between samples
  - Different sequencing depths
  - Sample mix-up
  Solution: Check raw .hic files and sequencing metrics

→ Negative values indicate data corruption
  Solution: Re-run extraction step
```

## Complete Output Files

All files saved to `outputs/res_{RES}kb/qc_report/`:

### Visualization PDFs
1. **`01_sample_correlation.pdf`**
   - 6×6 correlation heatmap
   - Numerical values overlaid
   - Validates replicate quality

2. **`02_sample_scatter.pdf`**
   - Sample-to-sample scatter plot
   - Perfect correlation line (y=x)
   - Linear fit with slope

3. **`03_count_distributions.pdf`**
   - 4-panel distribution comparison
   - Density, boxplot, violin, cumulative
   - Identifies distribution differences

4. **`04_ma_plot.pdf`**
   - MA plot with fold-change thresholds
   - Color-coded by magnitude
   - Assesses global bias

5. **`05_variable_loops_heatmap.pdf`**
   - Top 50 most variable loops
   - Hierarchical clustering
   - Z-score normalized

6. **`06_chromosome_distribution.pdf`**
   - Bar plot of loops per chromosome
   - Identifies chromosome-specific patterns

7. **`07_distance_decay.pdf`**
   - 4-panel distance analysis
   - Total, ctrl, mut, and FC vs distance
   - Validates expected Hi-C decay

8. **`08_example_matrices.pdf`**
   - Top 6 loops by signal
   - Side-by-side ctrl vs mut 5×5 matrices
   - Validates extraction quality

9. **`09_strategy_comparison.pdf`** *(if available)*
   - 4-panel aggregation comparison
   - Sum vs weighted strategies
   - Ratio distributions

### Data Files

10. **`qc_report_summary.rds`**
    - Complete QC metrics as R object
    - Can be loaded for further analysis
    - Contains all statistics and flags

11. **`loop_shift_status.rds`**
    - Boolean vector of shift detection
    - One TRUE/FALSE per loop
    - Used by edgeR for enrichment testing

12. **`loop_shift_summary.tsv`**
    - Human-readable shift table
    - Peak positions for ctrl and mut
    - Centered status per loop

## Typical QC Results

### Excellent Quality Example

```
QC ASSESSMENT:
✓✓✓ PASS - All quality checks passed! ✓✓✓

SUMMARY STATISTICS:
Dataset:             22,108 loops × 6 samples
Within-ctrl corr:    0.952 ± 0.004
Within-mut corr:     0.948 ± 0.006
Between-group corr:  0.883 ± 0.012
Library size CV:     0.08 (ctrl), 0.11 (mut)
Mean log2 FC:        0.023 ± 0.451
Loops with |M|>1:    2,156 (9.8%)
Positional shifts:   3,123 (14.1%)
Center enrichment:   0.312 (ctrl), 0.308 (mut)
```

### Good Quality with Warnings Example

```
QC ASSESSMENT:
✓ PASS (with warnings)

  ⚠ Moderate within-group correlation (0.88)
  ⚠ 23.4% positional shifts detected

SUMMARY STATISTICS:
Dataset:             18,543 loops × 6 samples
Within-ctrl corr:    0.882 ± 0.024
Within-mut corr:     0.891 ± 0.019
Between-group corr:  0.812 ± 0.031
Positional shifts:   4,339 (23.4%)

RECOMMENDATIONS:
→ Proceed with caution
→ One control replicate may be outlier (check correlation matrix)
→ Higher shift rate suggests more variability but buffer handles it
```

## Integration with Pipeline

### When to Run

**Option 1: After aggregation (before edgeR)**
```bash
Rscript scripts/prep_loops.R 5000
Rscript scripts/extract_counts.R 5000
Rscript scripts/aggregate.R 5000
Rscript scripts/qc-val.R 5000          # ← Run QC here
# Review QC report before proceeding
Rscript scripts/edgeR.R 5000
```

**Option 2: After full pipeline (retrospective)**
```bash
# Run full pipeline first
Rscript scripts/run_replicate_pipeline.sb
# Then review QC
Rscript scripts/qc-val.R 5000
```

### Multi-Resolution QC

```bash
# Run QC for all resolutions
for RES in 5000 10000 25000; do
  Rscript scripts/qc-val.R ${RES}
done

# Compare QC metrics across resolutions
# - Higher resolution → lower correlation (more noise)
# - Lower resolution → higher correlation (smoother signal)
# - Center enrichment may vary by resolution
```

## Interpreting QC for Biological Decisions

### High Within-Group Correlation (>0.95)
✓ **Good**: Excellent technical reproducibility
⚠ **Caution**: If between-group is also >0.95, biological signal may be weak

### Moderate Shift Rate (10-30%)
✓ **Expected**: Normal biological and technical variation
✓ **Buffer validated**: 5×5 buffer strategy is working as designed
⚠ **> 30%**: Consider larger buffer or different binning

### Low Center Enrichment (<0.20)
⚠ **Possible causes**:
- Loops are truly diffuse (biological)
- Loop calling imprecise (technical)
- Buffer too large for resolution

### Asymmetric Group Statistics
⚠ **One group has much lower correlation**:
- Check for outlier replicate
- Possible batch effect
- Different biological variability

### High MA Plot Scatter with Centered Mean
✓ **Good sign**: Biological variation with no systematic bias
✓ **Proceed** to edgeR which will model this variation

## Troubleshooting Common Issues

### Issue: "All data files not found"
**Cause**: Incomplete pipeline run or wrong directory
**Solution**: Verify all steps completed, check RESOLUTION parameter matches

### Issue: "Negative values detected"
**Cause**: Data corruption or incorrect aggregation
**Solution**: Re-run extraction and aggregation steps

### Issue: "Very low correlation (<0.5)"
**Cause**: Sample mix-up, severe batch effect, or data quality issues
**Solution**:
1. Check sample labels in input files
2. Review raw .hic file quality
3. Check for sequencing depth differences

### Issue: "Most loops not centered"
**Cause**: Loop calling imprecision or buffer size mismatch
**Solution**:
1. Check if HiCCUPS parameters appropriate for resolution
2. Consider buffer=3 instead of buffer=2
3. Compare with original HiCCUPS peak positions

## Next Steps

After reviewing QC report:

**If PASS**: Proceed to `edgeR.R` for differential analysis
**If PASS with warnings**: Note warnings in analysis interpretation
**If FAIL**: Address critical issues before continuing

The QC report provides foundation for interpreting all subsequent analysis results in proper technical and biological context.
