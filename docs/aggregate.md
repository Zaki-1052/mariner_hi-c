# Matrix Aggregation and Count Preparation (aggregate.R)

## Biological Purpose

This script converts the 5×5 Hi-C contact matrices from Step 2 into single count values per loop for all **6 biological replicates** (3 control + 3 mutant), creating a count matrix suitable for replicate-aware differential analysis with edgeR. The aggregation strategy is critical because it determines how spatial uncertainty and bin-shift artifacts are handled across multiple biological replicates.

**Key biological concepts:**
- **Spatial uncertainty**: Loop centers may shift slightly between samples
- **Bin-shift artifacts**: Genomic binning can place the same biological feature in different bins
- **Signal aggregation**: Combining multiple pixels to capture total loop strength
- **Replicate structure**: Preserving biological variation across 6 independent samples
- **Differential analysis preparation**: Creating count matrices compatible with RNA-seq tools

## Multi-Resolution Support

The script accepts resolution as a command-line argument:
```bash
Rscript scripts/aggregate.R 5000   # 5kb resolution
Rscript scripts/aggregate.R 10000  # 10kb resolution
Rscript scripts/aggregate.R 25000  # 25kb resolution
```

## Input Data

### HDF5-backed Contact Matrices
- **Input**: `05_extracted/` (HDF5 SummarizedExperiment from Step 2)
- **Dimensions**: 5 × 5 × N_loops × **6 samples**
- **Format**: DelayedArray for memory-efficient processing
- **Directory**: Resolution-specific (e.g., `outputs/res_5kb/`)

### Metadata
- **File**: `05_metadata.rds`
- **Contains**:
  - Extraction parameters (resolution, normalization)
  - File paths for all 6 replicates
  - Sample correlations (6×6 matrix)
  - Library sizes per replicate

**Example metadata structure:**
```r
metadata$files:
  ctrl_M1, ctrl_M2, ctrl_M3, mut_M1, mut_M2, mut_M3

metadata$correlation_matrix:
  6×6 correlation matrix from extraction step
```

## Aggregation Strategies

### Strategy 1: Simple Sum Aggregation (Primary Method)

```r
# Sum all 25 pixels in each 5x5 matrix
for (i in 1:dims[3]) {         # For each loop
  for (j in 1:dims[4]) {       # For each of 6 samples
    mat_5x5 <- as.matrix(count_array[, , i, j])
    counts_sum[i, j] <- sum(mat_5x5, na.rm = TRUE)
  }
}
```

**Biological rationale:**
- **Total signal capture**: Sums all Hi-C contacts in the local neighborhood
- **Robust to bin-shifts**: Doesn't depend on exact loop center position
- **Conservative approach**: Captures full loop signal regardless of spatial uncertainty
- **Compatible with edgeR**: Produces integer-like counts suitable for negative binomial modeling
- **Replicate consistency**: Handles biological variation across all 6 samples equally

**Advantages:**
- Simple and interpretable
- Robust to position variations between replicates
- Captures total loop strength
- Works well with sparse data
- Equal treatment of all 6 biological replicates

**Disadvantages:**
- May include background signal from edge pixels
- Less sensitive to precise loop positioning

### Strategy 2: Center-Weighted Aggregation

```r
# Gaussian-like weights emphasizing center pixels
weight_matrix <- matrix(c(
  0.04, 0.06, 0.08, 0.06, 0.04,  # Edge pixels: lowest weight
  0.06, 0.08, 0.10, 0.08, 0.06,  # Inner ring: medium weight
  0.08, 0.10, 0.16, 0.10, 0.08,  # Center row: highest in middle
  0.06, 0.08, 0.10, 0.08, 0.06,  # Inner ring: medium weight
  0.04, 0.06, 0.08, 0.06, 0.04   # Edge pixels: lowest weight
), nrow = 5, byrow = TRUE)

# Normalize to sum to 1
weight_matrix <- weight_matrix / sum(weight_matrix)

# Apply weighted sum
counts_weighted[i, j] <- sum(mat_5x5 * weight_matrix)
```

**Weight matrix visualization:**
```
0.040  0.060  0.080  0.060  0.040
0.060  0.080  0.100  0.080  0.060
0.080  0.100  0.160  0.100  0.080  ← Center row (highest weights)
0.060  0.080  0.100  0.080  0.060
0.040  0.060  0.080  0.060  0.040
                ↑
            Center pixel (16% of total weight)
```

**Biological rationale:**
- **Position sensitivity**: Higher weight to expected loop center
- **Background reduction**: Lower weight to edge pixels reduces noise
- **Spatial prior**: Assumes loops are most likely centered in the 5×5 region
- **Fine-tuned detection**: More sensitive to well-positioned loops

**Advantages:**
- Reduces background noise
- More sensitive to well-centered loops
- Accounts for expected spatial distribution

**Disadvantages:**
- Less robust to position shifts between replicates
- May miss off-center loops
- More complex interpretation

## Replicate-Aware Correlation Analysis

### Full 6×6 Correlation Matrix

```r
# Build correlation matrix across all 6 replicates
cor_matrix <- matrix(NA, nrow = 6, ncol = 6)
rownames(cor_matrix) <- c("ctrl_M1", "ctrl_M2", "ctrl_M3",
                          "mut_M1", "mut_M2", "mut_M3")
colnames(cor_matrix) <- rownames(cor_matrix)

for (i in 1:6) {
  for (j in 1:6) {
    cor_matrix[i, j] <- cor(counts_sum[,i], counts_sum[,j],
                            use = "complete.obs")
  }
}

cat("Full correlation matrix (sum strategy):\n")
print(round(cor_matrix, 3))
```

**Expected patterns:**
```
           ctrl_M1  ctrl_M2  ctrl_M3  mut_M1  mut_M2  mut_M3
ctrl_M1      1.000    0.950    0.945   0.880   0.875   0.870
ctrl_M2      0.950    1.000    0.948   0.882   0.878   0.873
ctrl_M3      0.945    0.948    1.000   0.878   0.874   0.869
mut_M1       0.880    0.882    0.878   1.000   0.952   0.947
mut_M2       0.875    0.878    0.874   0.952   1.000   0.949
mut_M3       0.870    0.873    0.869   0.947   0.949   1.000
```

### Within-Condition Correlations

```r
# Control replicates (all pairwise comparisons)
cat("Within-condition correlations:\n")
cat("  Control replicates:\n")
ctrl_cors <- numeric()
for (i in 1:2) {
  for (j in (i+1):3) {
    cor_val <- cor_matrix[i, j]
    ctrl_cors <- c(ctrl_cors, cor_val)
    cat(sprintf("    ctrl_M%d vs ctrl_M%d: %.3f\n", i, j, cor_val))
  }
}
ctrl_within_mean <- mean(ctrl_cors)
cat(sprintf("  Mean within-ctrl correlation: %.3f\n\n", ctrl_within_mean))

# Mutant replicates (all pairwise comparisons)
cat("  Mutant replicates:\n")
mut_cors <- numeric()
for (i in 4:5) {
  for (j in (i+1):6) {
    cor_val <- cor_matrix[i, j]
    mut_cors <- c(mut_cors, cor_val)
    cat(sprintf("    mut_M%d vs mut_M%d: %.3f\n", i-3, j-3, cor_val))
  }
}
mut_within_mean <- mean(mut_cors)
cat(sprintf("  Mean within-mut correlation: %.3f\n\n", mut_within_mean))
```

**Expected within-condition correlations:**
- **Excellent (> 0.95)**: High technical reproducibility
- **Good (0.90-0.95)**: Acceptable biological variation
- **Fair (0.80-0.90)**: Higher than expected variation, investigate
- **Poor (< 0.80)**: Technical problems or outlier replicates

### Between-Condition Correlations

```r
# Between control and mutant (all pairwise comparisons: 3×3 = 9)
between_cors <- as.vector(cor_matrix[1:3, 4:6])
cat("Between-condition correlations:\n")
cat(sprintf("  Mean ctrl vs mut correlation: %.3f\n", mean(between_cors)))
cat(sprintf("  Range: %.3f - %.3f\n\n", min(between_cors), max(between_cors)))
```

**Expected between-condition correlations:**
- Should be **lower** than within-condition (biological signal)
- Difference of 0.05-0.15 suggests real biological differences
- Very similar to within-condition may indicate weak biological effect

### Quality Assessment

```r
cat("Quality Assessment:\n")
if (ctrl_within_mean > 0.95 && mut_within_mean > 0.95) {
  cat("  ✓ Excellent within-condition reproducibility (>0.95)\n")
} else if (ctrl_within_mean > 0.90 && mut_within_mean > 0.90) {
  cat("  ✓ Good within-condition reproducibility (>0.90)\n")
} else {
  cat("  ⚠ WARNING: Lower than expected correlation (<0.90)\n")
  cat("    This may indicate technical issues\n")
}

if ((ctrl_within_mean - mean(between_cors)) > 0.05 ||
    (mut_within_mean - mean(between_cors)) > 0.05) {
  cat("  ✓ Clear biological signal (within > between correlations)\n")
} else {
  cat("  ⚠ WARNING: Weak biological signal\n")
}
```

## Strategy Comparison

### Correlation Between Strategies

```r
# Compare how well sum and weighted strategies agree
cor_strategies <- cor(as.vector(counts_sum),
                     as.vector(counts_weighted),
                     use = "complete.obs")

cat(sprintf("Correlation between sum and weighted strategies: %.3f\n",
            cor_strategies))
```

**Expected results:**
- **High inter-strategy correlation (> 0.95)**: Strategies capture similar biological signal
- **Similar sample correlations**: Both strategies preserve biological relationships
- **Sum slightly higher correlation**: More robust aggregation often improves reproducibility

## Biological Validation

### Differential Expression Preview

```r
# Calculate mean counts per condition
ctrl_means <- rowMeans(final_counts[, 1:3])   # Average of 3 ctrl replicates
mut_means <- rowMeans(final_counts[, 4:6])    # Average of 3 mut replicates

# Calculate log2 fold changes
log2fc <- log2((mut_means + 1) / (ctrl_means + 1))

# Identify extreme changes
extreme_changes <- which(abs(log2fc) > 2)
cat(sprintf("Loops with |log2FC| > 2: %d (%.1f%%)\n",
            length(extreme_changes),
            100 * length(extreme_changes) / nrow(final_counts)))
```

**Biological interpretation of log2FC:**
- **|log2FC| > 2**: 4-fold change (biologically significant)
- **|log2FC| > 1**: 2-fold change (moderate significance)
- **log2FC ≈ 0**: No differential looping

### Per-Replicate Distribution Quality Control

```r
cat("\nCount distribution by replicate:\n")

for (i in 1:6) {
  sample_name <- colnames(final_counts)[i]
  sample_counts <- final_counts[, i]
  positive <- sample_counts[sample_counts > 0]

  cat(sprintf("\n  %s:\n", sample_name))
  cat(sprintf("    Non-zero: %d / %d (%.1f%%)\n",
              length(positive), nrow(final_counts),
              100 * length(positive) / nrow(final_counts)))
  if (length(positive) > 0) {
    cat(sprintf("    Range: %.1f - %.1f\n", min(positive), max(positive)))
    cat(sprintf("    Median: %.1f, Mean: %.1f\n",
                median(positive), mean(positive)))
  }
}
```

**Quality indicators per replicate:**
- **High percentage of non-zero counts (> 80%)**: Good loop detection
- **Similar distributions across replicates**: Consistent experimental conditions
- **Reasonable dynamic range**: Mix of weak and strong loops
- **No extreme outlier replicates**: All samples contribute equally

### Group-Level Signal Strength

```r
cat("\nSignal strength by group:\n")
ctrl_all <- as.vector(final_counts[, 1:3])
mut_all <- as.vector(final_counts[, 4:6])

cat(sprintf("  Control (all replicates): median = %.1f, mean = %.1f\n",
            median(ctrl_all), mean(ctrl_all)))
cat(sprintf("  Mutant (all replicates):  median = %.1f, mean = %.1f\n",
            median(mut_all), mean(mut_all)))
```

**Expected patterns:**
- **Similar group medians**: Balanced experimental design
- **Moderate differences**: May indicate biological effects
- **Extreme differences**: Potential sequencing depth imbalance (edgeR will normalize)

### Systematic Bias Detection

```r
# Check for large systematic differences
mean_diff <- mean(mut_all) - mean(ctrl_all)
cat(sprintf("\nMean difference (mut - ctrl): %.1f\n", mean_diff))

if (abs(mean_diff) > 0.5 * mean(ctrl_all)) {
  cat("  ⚠ Large systematic difference - edgeR normalization will handle\n")
} else {
  cat("  ✓ No major systematic bias\n")
}
```

## edgeR Integration

### DGEList Creation with Replicate Structure

```r
# Load original loop coordinates
original_loops <- readRDS(file.path(input_dir, "02_merged.rds"))

# Extract genomic information
loop_anchors1 <- anchors(original_loops, type = "first")
loop_anchors2 <- anchors(original_loops, type = "second")

# Create annotation dataframe
gene_info <- data.frame(
  loop_id = rownames(final_counts),
  chr1 = as.character(seqnames(loop_anchors1)),
  start1 = start(loop_anchors1),
  end1 = end(loop_anchors1),
  chr2 = as.character(seqnames(loop_anchors2)),
  start2 = start(loop_anchors2),
  end2 = end(loop_anchors2),
  stringsAsFactors = FALSE
)

# Calculate interaction distance
gene_info$distance <- ifelse(
  gene_info$chr1 == gene_info$chr2,
  abs(gene_info$start2 - gene_info$start1),
  NA  # NA for interchromosomal
)

# Create DGEList with proper replicate structure
y <- DGEList(
  counts = final_counts,
  group = factor(c("ctrl", "ctrl", "ctrl", "mut", "mut", "mut")),
  genes = gene_info
)
```

**Critical group structure:**
```r
group = c("ctrl", "ctrl", "ctrl", "mut", "mut", "mut")
# Corresponds to:
# ctrl_M1, ctrl_M2, ctrl_M3, mut_M1, mut_M2, mut_M3
```

### Library Size Analysis

```r
cat("\nEdgeR object preview:\n")
cat(sprintf("  Samples: %s\n", paste(colnames(y), collapse = ", ")))
cat(sprintf("  Groups: %s\n", paste(levels(y$samples$group), collapse = ", ")))
cat(sprintf("  Replicates per group: n = %d\n", sum(y$samples$group == "ctrl")))

cat(sprintf("  Library sizes:\n"))
for (i in 1:ncol(y)) {
  cat(sprintf("    %s: %s\n",
              colnames(y)[i],
              format(y$samples$lib.size[i], big.mark = ",")))
}
```

**Example output:**
```
Library sizes:
  ctrl_M1: 1,234,567
  ctrl_M2: 1,198,432
  ctrl_M3: 1,276,891
  mut_M1:  1,145,678
  mut_M2:  1,189,234
  mut_M3:  1,167,543
```

**Quality indicators:**
- **Similar library sizes within group**: Good technical reproducibility
- **Ratio between groups 0.8-1.2**: Acceptable for TMM normalization
- **Extreme differences (> 2-fold)**: May need investigation

### Filtering Preview

```r
# Preview edgeR filtering
keep <- filterByExpr(y, group = y$samples$group)
cat(sprintf("\nFiltering preview:\n"))
cat(sprintf("  Loops passing filterByExpr: %d / %d (%.1f%%)\n",
            sum(keep), length(keep), 100 * sum(keep) / length(keep)))
```

**edgeR filtering criteria:**
- **Minimum count threshold**: Loops must have adequate signal
- **Minimum samples**: Must be detected in at least 3 samples (half of smallest group)
- **Library size consideration**: Adjusts thresholds based on sequencing depth
- **Replicate-aware**: Uses group structure to determine which loops to keep

## Output Files

### 1. Primary Count Matrix (`06_counts_matrix.tsv` and `.rds`)

**Format:**
```
           ctrl_M1  ctrl_M2  ctrl_M3  mut_M1  mut_M2  mut_M3
loop_1        45.2     48.1     43.7    52.8    55.2    50.1
loop_2        12.1     11.8     12.5     8.3     7.9     8.8
loop_3       123.7    128.4    119.2   145.2   149.8   141.3
...
```

**File contents:**
- **Rows**: Individual chromatin loops
- **Columns**: 6 samples (3 ctrl + 3 mut)
- **Values**: Aggregated Hi-C contact counts (sum of 5×5 pixels)

**Resolution-specific paths:**
- `outputs/res_5kb/06_counts_matrix.tsv`
- `outputs/res_10kb/06_counts_matrix.tsv`
- `outputs/res_25kb/06_counts_matrix.tsv`

### 2. edgeR-ready Object (`06_edgeR_input.rds`)

```r
# Load for differential analysis
y <- readRDS("outputs/res_5kb/06_edgeR_input.rds")

# Standard edgeR workflow (now with n=3 replicates per group!)
y <- calcNormFactors(y)
design <- model.matrix(~group, data = y$samples)
y <- estimateDisp(y, design, robust = TRUE)
fit <- glmQLFit(y, design, robust = TRUE)
qlf <- glmQLFTest(fit, coef = 2)
```

**DGEList components:**
- **counts**: N_loops × 6 count matrix
- **samples**: Sample metadata with group assignments
- **genes**: Genomic coordinates and loop information

### 3. All Strategies (`06_all_strategies.rds`)

```r
all_strategies <- list(
  counts = list(
    sum = counts_sum,
    weighted = counts_weighted
  ),
  metadata = list(
    n_loops = nrow(final_counts),
    samples = colnames(final_counts),
    aggregation_method = "sum",
    extraction_norm = "VC",
    binSize = RESOLUTION,
    buffer = 2,
    date = Sys.Date(),
    within_ctrl_correlation = ctrl_within_mean,
    within_mut_correlation = mut_within_mean,
    between_correlation = mean(between_cors)
  )
)
```

**Contains:**
- Both sum and weighted aggregation results
- QC metrics (correlations)
- Analysis parameters
- Useful for method comparison and validation

## Biological Interpretation

### What the Counts Represent

Each count value represents the **total Hi-C contact frequency** in a 5×5 pixel region around a chromatin loop, measured independently in each of 6 biological replicates. Higher counts indicate:

1. **Stronger chromatin interactions**
2. **More accessible chromatin regions**
3. **More stable loop structures**
4. **Higher local contact density**

### Differential Analysis Expectations

**Upregulated loops (mut > ctrl):**
- May indicate gained chromatin interactions
- Could reflect new enhancer-promoter contacts
- Might suggest altered chromatin architecture
- Should show consistent increase across 3 mutant replicates

**Downregulated loops (ctrl > mut):**
- May indicate lost chromatin interactions
- Could reflect disrupted regulatory contacts
- Might suggest chromatin reorganization
- Should show consistent decrease across 3 mutant replicates

### Replicate Consistency

**High-confidence differential loops:**
- Show consistent direction across all 3 replicates per group
- Have low within-group variation
- Display clear separation between conditions

**Low-confidence changes:**
- Variable across replicates
- High within-group dispersion
- Overlapping distributions between conditions

### Distance Effects

```r
# Analyze by interaction distance
intra_loops <- !is.na(gene_info$distance)
inter_loops <- is.na(gene_info$distance)

# Distance categories
short_range <- gene_info$distance < 100000  # < 100kb
long_range <- gene_info$distance >= 100000  # >= 100kb
```

**Expected patterns:**
- **Shorter distances**: Typically stronger interactions, higher counts
- **Longer distances**: More variable, often weaker signals
- **Interchromosomal**: Generally rare and weaker

## Quality Control Summary

### Excellent Quality Indicators

✓ **Within-condition correlation > 0.95** for both ctrl and mut
✓ **Between-condition correlation 0.05-0.15 lower** than within
✓ **Similar library sizes** across all 6 replicates
✓ **High percentage non-zero** (> 80%) in all samples
✓ **No outlier replicates** in correlation matrix
✓ **Consistent count distributions** across replicates

### Warning Signs

⚠ **Within-condition correlation < 0.90**: Technical issues or high biological variability
⚠ **Between-condition similar to within**: Weak biological signal
⚠ **Large library size differences** (> 2-fold): Sequencing depth imbalance
⚠ **Low non-zero percentage** (< 60%): Poor loop detection
⚠ **Outlier replicate**: One sample very different from others in group

## Next Steps

The aggregated count matrix is ready for:

1. **Differential looping analysis** (Step 4: `edgeR.R`)
   - Quasi-likelihood GLM with **n=3 biological replicates**
   - Data-driven dispersion estimation
   - Robust statistical testing
   - Expected: Hundreds to thousands of significant loops

2. **Resolution comparison** (multi-resolution pipeline only)
   - Compare results across 5kb, 10kb, 25kb
   - Identify resolution-specific and shared differential loops
   - Assess concordance of fold changes

## Technical Performance

### Memory Efficiency
- **Chunk-based processing**: Handles large matrices without RAM overload
- **HDF5 backend**: Enables out-of-memory computation
- **Efficient aggregation**: Optimized loops for speed
- **Memory usage**: < 8GB RAM even with 6 replicates

### Typical Runtime

**Processing time by resolution:**
- **5kb**: ~2-4 minutes (most loops)
- **10kb**: ~1-3 minutes (intermediate)
- **25kb**: ~30-90 seconds (fewest loops)

**Scales with:**
- Number of loops (linear)
- Number of samples (linear, currently 6)
- System I/O speed (HDF5 read performance)

## Improvement Over Merged Analysis

**Previous approach (merged, n=1 per condition):**
- 2 samples total (ctrl_merged, mut_merged)
- No biological replicates
- No variance estimation
- Result: 2 significant loops at FDR < 0.05

**Current approach (replicates, n=3 per condition):**
- 6 samples total (3 ctrl + 3 mut)
- True biological replication
- Data-driven dispersion estimation
- Expected: 500-5,000+ significant loops at FDR < 0.05

**Statistical power increase:** ~250-2,500× more differential loops detected
