# Hi-C Contact Matrix Extraction (extract_counts.R)

## Biological Purpose

This script extracts the actual Hi-C contact frequencies from individual replicate `.hic` files at each loop position identified in Step 1. It creates 5×5 pixel contact matrices around each loop for **all 6 biological replicates** (3 control + 3 mutant), providing the quantitative data needed for replicate-aware differential analysis.

**Key biological concepts:**
- **Contact matrices**: 2D representations of chromatin interaction frequencies
- **Hi-C normalization**: Correction for systematic biases in Hi-C data
- **Biological replicates**: Independent measurements of the same biological system
- **Spatial resolution**: Resolution-specific bins provide appropriate interaction resolution
- **Local contact neighborhoods**: 5×5 regions capture interaction context around loop centers

## Multi-Resolution Support

The script accepts resolution as a command-line argument:
```bash
Rscript scripts/extract_counts.R 5000   # 5kb resolution
Rscript scripts/extract_counts.R 10000  # 10kb resolution
Rscript scripts/extract_counts.R 25000  # 25kb resolution
```

## Input Data

### .hic Files (Juicer format) - Individual Replicates

**Control replicates:**
```
/expanse/lustre/projects/csd940/ctea/nf-hic/juicerpre/ctrl_M1.hic
/expanse/lustre/projects/csd940/ctea/nf-hic/juicerpre/ctrl_M2.hic
/expanse/lustre/projects/csd940/ctea/nf-hic/juicerpre/ctrl_M3.hic
```

**Mutant replicates:**
```
/expanse/lustre/projects/csd940/ctea/nf-hic/juicerpre/mut_M1.hic
/expanse/lustre/projects/csd940/ctea/nf-hic/juicerpre/mut_M2.hic
/expanse/lustre/projects/csd940/ctea/nf-hic/juicerpre/mut_M3.hic
```

### Buffered Loop Positions
- **Input**: `04_buffered.rds` from Step 1
- **Format**: GInteractions with 5×5 pixel regions around loop centers
- **Directory**: Resolution-specific (e.g., `outputs/res_5kb/`)

## Technical Validation

### File Verification
```r
# Check file existence and sizes for all 6 replicates
for (name in names(hicFiles)) {
  if (!file.exists(hicFiles[name])) {
    stop(sprintf("ERROR: %s file not found", name))
  }
  size_gb <- file.info(hicFiles[name])$size / 1e9
  cat(sprintf("  %s: %.1f GB - OK\n", name, size_gb))
}
```

**Expected file sizes:**
- Individual replicates: ~3-8 GB each
- Total: ~36-48 GB across all 6 files

### Resolution and Normalization Checking
```r
# Verify resolution is available in all replicates
for (name in names(hicFiles)) {
  resolutions <- readHicBpResolutions(hicFiles[name])
  if (!RESOLUTION %in% resolutions) {
    stop(sprintf("ERROR: %dkb resolution not available in %s",
                 RESOLUTION/1000, name))
  }

  # Check available normalizations
  norms <- readHicNormTypes(hicFiles[name])
}
```

**Normalization types:**
- **VC (Vanilla Coverage)**: Corrects for coverage biases (used in this pipeline)
- **KR (Knight-Ruiz)**: More sophisticated matrix balancing (often unavailable)
- **NONE**: Raw counts without correction

## Matrix Extraction Process

### Core Extraction Function
```r
pixels <- pullHicMatrices(
  x = buffered,                 # Loop positions with 5×5 buffers
  files = hicFiles,            # All 6 replicate .hic files
  binSize = RESOLUTION,        # Resolution-specific binning
  h5File = h5_file_path,       # HDF5 output for memory efficiency
  norm = "VC",                 # Vanilla Coverage normalization
  matrix = "observed",         # Extract observed (not expected) counts
  blockSize = 1e6,            # Memory management parameter
  onDisk = TRUE,              # Store results on disk (not in RAM)
  compressionLevel = 1         # HDF5 compression
)
```

**Resolution-specific parameters:**
- All resolutions use same buffer=2 approach
- Extraction time scales with resolution (coarser = faster)
- HDF5 file size varies with loop count and resolution

### Output Dimensions
```
Dimensions: 5 × 5 × [N_loops] × 6
  = 5×5 pixels × [N] loops × 6 replicates
```

**Data structure interpretation:**
- **5×5 matrices**: Local contact neighborhood around each loop
- **N_loops**: Total number of consensus loops from Step 1
- **6 replicates**: 3 control + 3 mutant for replicate-aware analysis

**Example dimensions:**
- 5kb: 5 × 5 × 22,108 × 6 = ~6.6 million values
- 10kb: 5 × 5 × 18,000 × 6 = ~5.4 million values
- 25kb: 5 × 5 × 12,000 × 6 = ~3.6 million values

## Biological Data Validation

### Contact Matrix Quality Checks

#### 1. NA Value Assessment
```r
na_count <- sum(is.na(count_array))
na_percent <- 100 * na_count / length(count_array)

cat(sprintf("NA values: %d (%.2f%%)\n", na_count, na_percent))
```

**Expected patterns:**
- Low NA percentage (< 10%) indicates good data quality
- NAs typically from unmappable regions or chromosome boundaries
- Should be similar across all 6 replicates

#### 2. Value Distribution Analysis
```r
# Overall distribution statistics
non_na_values <- count_array[!is.na(count_array)]
summary(non_na_values)

# Sparsity calculation
zeros <- sum(non_na_values == 0)
sparsity <- 100 * zeros / length(non_na_values)
cat(sprintf("Sparsity: %.1f%% zeros\n", sparsity))
```

**Biological interpretation:**
- **High sparsity (70-90% zeros)**: Normal for Hi-C data due to chromatin accessibility
- **Extreme values**: May indicate technical artifacts or very strong interactions
- **Similar sparsity across replicates**: Indicates consistent data quality

#### 3. Per-Sample Statistics
```r
for (i in 1:dims[4]) {
  sample_name <- names(hicFiles)[i]
  sample_values <- count_array[,,,i][!is.na(count_array[,,,i])]

  cat(sprintf("  %s: median=%.2f, mean=%.2f, max=%.0f, total=%.0f\n",
              sample_name,
              median(sample_values),
              mean(sample_values),
              max(sample_values),
              sum(sample_values)))
}
```

**Key metrics:**
- **Median**: Typical contact strength (often low due to sparsity)
- **Mean**: Average contact strength (affected by high-value outliers)
- **Max**: Strongest interactions (should be reasonable, not extreme outliers)
- **Total**: Library size (should be comparable across replicates)

### Library Size Comparison
```r
# Compare library sizes by condition
ctrl_median <- median(sample_stats$total[1:3])
mut_median <- median(sample_stats$total[4:6])

cat(sprintf("Library size ratios:\n"))
cat(sprintf("  Median ctrl: %.0f\n", ctrl_median))
cat(sprintf("  Median mut: %.0f\n", mut_median))
cat(sprintf("  Ratio (ctrl/mut): %.3f\n", ctrl_median/mut_median))
```

**Quality indicators:**
- **Ratio near 1.0**: Balanced sequencing depth
- **Ratio 0.8-1.2**: Acceptable, edgeR normalization will handle
- **Ratio < 0.5 or > 2.0**: Large imbalance, investigate

## Comprehensive Correlation Analysis

### Full Correlation Matrix (6×6)
```r
# Build correlation matrix across all 6 replicates
cor_matrix <- matrix(NA, nrow = 6, ncol = 6)
rownames(cor_matrix) <- names(hicFiles)
colnames(cor_matrix) <- names(hicFiles)

for (i in 1:6) {
  for (j in 1:6) {
    values_i <- as.vector(count_array[,,,i])
    values_j <- as.vector(count_array[,,,j])
    complete_pairs <- complete.cases(values_i, values_j)
    cor_matrix[i, j] <- cor(values_i[complete_pairs],
                            values_j[complete_pairs])
  }
}

cat("Full correlation matrix:\n")
print(round(cor_matrix, 3))
```

### Within-Condition Correlations
```r
# Control replicates
cat("Within-condition correlations:\n")
cat("  Control replicates:\n")
ctrl_cors <- numeric()
for (i in 1:2) {
  for (j in (i+1):3) {
    cor_val <- cor_matrix[i, j]
    ctrl_cors <- c(ctrl_cors, cor_val)
    cat(sprintf("    %s vs %s: %.3f\n",
                names(hicFiles)[i], names(hicFiles)[j], cor_val))
  }
}
ctrl_within_mean <- mean(ctrl_cors)
cat(sprintf("  Mean within-ctrl correlation: %.3f\n\n", ctrl_within_mean))

# Mutant replicates
cat("  Mutant replicates:\n")
mut_cors <- numeric()
for (i in 4:5) {
  for (j in (i+1):6) {
    cor_val <- cor_matrix[i, j]
    mut_cors <- c(mut_cors, cor_val)
    cat(sprintf("    %s vs %s: %.3f\n",
                names(hicFiles)[i], names(hicFiles)[j], cor_val))
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
# Between control and mutant
between_cors <- as.vector(cor_matrix[1:3, 4:6])
cat(sprintf("Between-condition correlations:\n"))
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

## Example Matrix Visualization

### Individual Loop Matrices
```r
cat("Example extracted matrices:\n")
cat("Showing first 3 loops for control replicate 1:\n\n")

for (i in 1:min(3, dims[3])) {
  cat(sprintf("Loop %d:\n", i))

  # Extract and convert to regular matrix
  example_matrix <- count_array[, , i, 1]  # Loop i, ctrl_M1
  example_matrix <- as.matrix(example_matrix)

  # Display formatted matrix
  for (row in 1:5) {
    row_values <- round(example_matrix[row, ], 2)
    formatted_row <- format(row_values, width = 8, justify = "right")
    cat("  ", paste(formatted_row, collapse = " "), "\n")
  }

  # Summary stats
  matrix_vals <- as.vector(example_matrix)
  matrix_vals <- matrix_vals[!is.na(matrix_vals)]
  if (length(matrix_vals) > 0) {
    cat(sprintf("  Sum: %.0f, Max: %.0f, Non-zero: %d/25\n\n",
                sum(matrix_vals), max(matrix_vals), sum(matrix_vals > 0)))
  }
}
```

**Matrix interpretation example:**
```
Loop 1:
    12.5      8.2     15.4      6.1      3.2
     8.7     22.1     31.8     18.5      7.4
    15.9     31.2     45.6     29.3     14.1  ← Center row (loop focus)
     6.8     18.9     28.7     21.4      8.9
     3.1      7.6     13.8      9.2      4.5
                       ↑
                   Center pixel (loop center)
  Sum: 392, Max: 46, Non-zero: 25/25
```

**Expected patterns:**
- **Central enrichment**: Highest values near center (3,3) position
- **Symmetric decay**: Values decrease with distance from center
- **Biological variation**: Some asymmetry is normal due to chromatin context
- **Non-zero pixels**: Most pixels should have signal (not all zeros)

## Output Files

### 1. HDF5-backed InteractionArray (`05_extracted/`)
```r
saveHDF5SummarizedExperiment(
  x = pixels,
  dir = input_dir,              # Resolution-specific directory
  prefix = "05_extracted",
  replace = TRUE,
  verbose = TRUE
)
```

**Benefits of HDF5 format:**
- **Memory efficiency**: Handles large datasets (6 replicates) without loading into RAM
- **Chunk-based access**: Can read specific loops/samples without loading all data
- **Cross-platform compatibility**: Standard format for large genomic datasets
- **Compression**: Reduces disk space requirements

**File structure:**
```
outputs/res_5kb/05_extracted/
  ├── assays.h5          # HDF5 array with contact matrices
  ├── se.rds             # R object with metadata
  └── [other metadata]
```

### 2. Metadata File (`05_metadata.rds`)
```r
metadata <- list(
  resolution = RESOLUTION,
  n_loops = dims[3],
  n_replicates = dims[4],
  replicate_names = names(hicFiles),
  files = hicFiles,
  normalization = norm_to_use,
  extraction_time = extraction_time,
  dims = dims,
  sample_stats = sample_stats,
  correlation_matrix = cor_matrix,
  within_ctrl_correlation = ctrl_within_mean,
  within_mut_correlation = mut_within_mean,
  between_correlation = mean(between_cors)
)
```

**Metadata contents:**
- Extraction parameters (resolution, normalization)
- Sample information (names, file paths)
- Quality metrics (correlations, library sizes)
- Timing information (for performance tracking)

## Biological Significance

### Contact Matrix Interpretation

**What each 5×5 matrix represents:**
1. **Central pixel (3,3)**: Direct loop contact strength
2. **Adjacent pixels**: Local chromatin context
3. **Edge pixels**: Background interaction levels
4. **Asymmetry**: Directional chromatin organization effects

### Replicate-Specific Patterns

**Within-replicate consistency:**
- Similar matrices across replicates indicate robust loop calling
- Variation may reflect biological heterogeneity or technical noise
- Extreme outliers in single replicates warrant investigation

**Between-condition differences:**
- Systematic differences suggest biological regulation
- Condition-specific enrichment indicates differential looping
- Magnitude of differences informs statistical power

### Quality Control Metrics

**Good extraction indicators:**
- High within-condition correlation (> 0.90)
- Lower between-condition correlation (biological signal)
- Low percentage of NA values (< 15%)
- Central enrichment in example matrices
- Consistent library sizes across replicates

**Warning signs:**
- Very low within-condition correlation (< 0.80) - technical problems
- High NA percentage (> 30%) - sparse or problematic regions
- Uniform matrices - loss of spatial structure
- Extreme library size imbalances (> 2-fold)

## Next Steps

The extracted contact matrices are ready for aggregation (Step 3: `aggregate.R`), where the 5×5 matrices will be collapsed into single values per loop for all 6 replicates, creating a count matrix suitable for replicate-aware differential analysis with edgeR.

## Technical Performance

### Memory Management
- **HDF5 backend**: Enables analysis of 6 replicates simultaneously
- **Chunk processing**: Processes loops in batches
- **Compression**: Reduces storage requirements (~30-50% reduction)
- **On-disk computation**: Minimal RAM usage regardless of dataset size

### Typical Performance

**Extraction time by resolution:**
- **5kb**: ~5-10 minutes (most loops, largest file)
- **10kb**: ~3-7 minutes (intermediate)
- **25kb**: ~2-5 minutes (fewest loops, smallest file)

**Memory usage:** < 4GB RAM (due to on-disk HDF5 storage)

**File sizes (HDF5):**
- 5kb: ~500-800 MB
- 10kb: ~400-600 MB
- 25kb: ~300-400 MB

### Scalability

**Current pipeline (6 replicates):**
- Handles 20,000+ loops efficiently
- Extraction time linear with loop count
- Memory usage independent of replicate count (HDF5 backend)

**Potential scaling:**
- Could handle 10+ replicates with same approach
- Performance bottleneck is disk I/O, not RAM
- Parallel extraction possible for very large datasets
