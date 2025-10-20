# Matrix Aggregation and Count Preparation (aggregate.R)

## Biological Purpose

This script converts the 5x5 Hi-C contact matrices from Step 2 into single count values per loop, creating a count matrix suitable for differential expression analysis with edgeR. The aggregation strategy is critical because it determines how spatial uncertainty and bin-shift artifacts are handled.

**Key biological concepts:**
- **Spatial uncertainty**: Loop centers may shift slightly between samples
- **Bin-shift artifacts**: Genomic binning can place the same biological feature in different bins
- **Signal aggregation**: Combining multiple pixels to capture total loop strength
- **Differential analysis preparation**: Creating count matrices compatible with RNA-seq tools

## Input Data

### HDF5-backed Contact Matrices
- **Input**: `05_extracted/` (HDF5 SummarizedExperiment)
- **Dimensions**: 5 × 5 × N_loops × 2_samples
- **Format**: DelayedArray for memory-efficient processing

### Metadata
- **File**: `05_metadata.rds`
- **Contains**: Extraction parameters, file paths, sample correlations

## Aggregation Strategies

### Strategy 1: Simple Sum Aggregation
```r
# Sum all 25 pixels in each 5x5 matrix
for (i in 1:dims[3]) {
  for (j in 1:dims[4]) {
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

**Advantages:**
- Simple and interpretable
- Robust to position variations between samples
- Captures total loop strength
- Works well with sparse data

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
- **Spatial prior**: Assumes loops are most likely centered in the 5x5 region
- **Fine-tuned detection**: More sensitive to well-positioned loops

**Advantages:**
- Reduces background noise
- More sensitive to well-centered loops
- Accounts for expected spatial distribution

**Disadvantages:**
- Less robust to position shifts
- May miss off-center loops
- More complex interpretation

## Strategy Comparison

### Correlation Analysis
```r
# Compare sample correlations for each strategy
for (name in names(strategies)) {
  mat <- strategies[[name]]
  cor_val <- cor(mat[,1], mat[,2], use = "complete.obs")
  cat(sprintf("  %s: %.3f\n", name, cor_val))
}

# Correlation between strategies
cor_strategies <- cor(as.vector(counts_sum),
                     as.vector(counts_weighted),
                     use = "complete.obs")
```

**Expected results:**
- **High inter-strategy correlation (> 0.9)**: Strategies capture similar biological signal
- **Similar sample correlations**: Both strategies preserve biological relationships
- **Sum slightly higher correlation**: More robust aggregation often improves reproducibility

## Biological Validation

### Differential Expression Analysis Preparation
```r
# Calculate log2 fold changes
log2fc <- log2((final_counts[,2] + 1) / (final_counts[,1] + 1))

# Identify extreme changes
extreme_changes <- which(abs(log2fc) > 2)
```

**Biological interpretation of log2FC:**
- **|log2FC| > 2**: 4-fold change (biologically significant)
- **|log2FC| > 1**: 2-fold change (moderate significance)
- **log2FC ≈ 0**: No differential looping

### Distribution Quality Control
```r
# Check count distributions
ctrl_positive <- final_counts[final_counts[,1] > 0, 1]
mut_positive <- final_counts[final_counts[,2] > 0, 2]

# Distribution statistics
cat(sprintf("Non-zero: %d / %d (%.1f%%)\n",
            length(ctrl_positive), nrow(final_counts),
            100 * length(ctrl_positive) / nrow(final_counts)))
```

**Quality indicators:**
- **High percentage of non-zero counts (> 80%)**: Good loop detection
- **Similar distributions between samples**: Consistent experimental conditions
- **Reasonable dynamic range**: Mix of weak and strong loops

### Systematic Bias Detection
```r
mean_diff <- mean(final_counts[,2]) - mean(final_counts[,1])
if (abs(mean_diff) > 0.5 * mean(final_counts[,1])) {
  cat("Large systematic difference detected - may need normalization\n")
}
```

**Bias indicators:**
- **Large mean differences**: May indicate sequencing depth differences
- **Systematic shifts**: Could suggest technical artifacts
- **Solution**: edgeR normalization will correct for library size differences

## edgeR Integration

### DGEList Creation
```r
# Load original loop coordinates
original_loops <- readRDS("outputs/test/02_merged.rds")

# Extract genomic information
gene_info <- data.frame(
  loop_id = rownames(final_counts),
  chr1 = as.character(seqnames(loop_anchors1)),
  start1 = start(loop_anchors1),
  end1 = end(loop_anchors1),
  chr2 = as.character(seqnames(loop_anchors2)),
  start2 = start(loop_anchors2),
  end2 = end(loop_anchors2)
)

# Calculate interaction distances
gene_info$distance <- ifelse(
  gene_info$chr1 == gene_info$chr2,
  abs(gene_info$start2 - gene_info$start1),
  NA  # NA for interchromosomal
)

# Create edgeR object
y <- DGEList(
  counts = final_counts,
  group = factor(c("ctrl", "mut")),
  genes = gene_info
)
```

### Filtering Preview
```r
keep <- filterByExpr(y, group = y$samples$group)
cat(sprintf("Loops passing filterByExpr: %d / %d (%.1f%%)\n",
            sum(keep), length(keep), 100 * sum(keep) / length(keep)))
```

**edgeR filtering criteria:**
- **Minimum count threshold**: Loops must have adequate signal
- **Minimum samples**: Loops must be detected in sufficient replicates
- **Library size consideration**: Adjusts thresholds based on sequencing depth

## Output Files

### 1. Count Matrix (`06_counts_matrix.tsv` and `.rds`)
```
           ctrl    mut
loop_1     45.2   52.8
loop_2     12.1    8.3
loop_3    123.7  145.2
...
```

**Format description:**
- **Rows**: Individual chromatin loops
- **Columns**: Samples (control, mutant)
- **Values**: Aggregated Hi-C contact counts

### 2. edgeR-ready Object (`06_edgeR_input.rds`)
```r
# Load for differential analysis
y <- readRDS("outputs/full/06_edgeR_input.rds")

# Standard edgeR workflow
y <- calcNormFactors(y)
design <- model.matrix(~group, data = y$samples)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = 2)
```

### 3. All Strategies (`06_all_strategies.rds`)
- Contains both sum and weighted aggregation results
- Useful for method comparison and validation

## Biological Interpretation

### What the Counts Represent
Each count value represents the **total Hi-C contact frequency** in a 5×5 pixel region around a chromatin loop. Higher counts indicate:

1. **Stronger chromatin interactions**
2. **More accessible chromatin regions**
3. **More stable loop structures**
4. **Higher local contact density**

### Differential Analysis Expectations
**Upregulated loops (mut > ctrl):**
- May indicate gained chromatin interactions
- Could reflect new enhancer-promoter contacts
- Might suggest altered chromatin architecture

**Downregulated loops (ctrl > mut):**
- May indicate lost chromatin interactions
- Could reflect disrupted regulatory contacts
- Might suggest chromatin reorganization

### Distance Effects
```r
# Intrachromosomal vs interchromosomal
intra_loops <- !is.na(gene_info$distance)
inter_loops <- is.na(gene_info$distance)
```

**Expected patterns:**
- **Shorter distances**: Typically stronger interactions
- **Longer distances**: More variable, often weaker
- **Interchromosomal**: Generally rare and weaker

## Next Steps

The aggregated count matrix is ready for:
1. **Quality control visualization** (Step 4: `qc.R`)
2. **Exploratory data analysis** (Step 5: `exp.R`)
3. **Differential looping analysis** using edgeR
4. **Functional annotation** of significant loops

## Technical Considerations

### Memory Efficiency
- **Chunk-based processing**: Handles large matrices without RAM overload
- **HDF5 backend**: Enables out-of-memory computation
- **Efficient aggregation**: Optimized loops for speed

### Reproducibility
- **Fixed parameters**: Consistent aggregation across runs
- **Documented methods**: Clear strategy selection rationale
- **Version control**: Saved parameters in metadata