# Hi-C Contact Matrix Extraction (extract_counts.R)

## Biological Purpose

This script extracts the actual Hi-C contact frequencies from the raw `.hic` files at each loop position identified in Step 1. It creates 5x5 pixel contact matrices around each loop, providing the quantitative data needed for differential analysis.

**Key biological concepts:**
- **Contact matrices**: 2D representations of chromatin interaction frequencies
- **Hi-C normalization**: Correction for systematic biases in Hi-C data
- **Spatial resolution**: 5kb bins provide gene-level interaction resolution
- **Local contact neighborhoods**: 5x5 regions capture interaction context around loop centers

## Input Data

### .hic Files (Juicer format)
- **Control**: `/expanse/lustre/projects/csd940/ctea/nf-hic/juicerpre/merged/hic/resorted_ctrl.hic`
- **Mutant**: `/expanse/lustre/projects/csd940/ctea/nf-hic/juicerpre/merged/hic/resorted_mut.hic`

### Buffered Loop Positions
- **Input**: `04_buffered.rds` from Step 1
- **Format**: GInteractions with 5x5 pixel regions around loop centers

## Technical Validation

### File Verification
```r
# Check file existence and sizes
for (name in names(hicFiles)) {
  size_gb <- file.info(filepath)$size / 1e9
  cat(sprintf("  %s: %.1f GB - OK\n", name, size_gb))
}
```

### Resolution and Normalization Checking
```r
# Verify 5kb resolution is available
resolutions <- readHicBpResolutions(hicFiles[name])
if (!5000 %in% resolutions) {
  stop("ERROR: 5kb resolution not available")
}

# Check available normalizations
norms <- readHicNormTypes(hicFiles[name])
# Uses VC normalization (Vanilla Coverage)
```

**Normalization types:**
- **VC (Vanilla Coverage)**: Corrects for coverage biases
- **KR (Knight-Ruiz)**: More sophisticated matrix balancing (often unavailable)
- **NONE**: Raw counts without correction

## Matrix Extraction Process

### Core Extraction Function
```r
pixels <- pullHicMatrices(
  x = buffered,                 # Loop positions with 5x5 buffers
  files = hicFiles,            # Control and mutant .hic files
  binSize = 5e3,               # 5kb resolution
  h5File = h5_file_path,       # HDF5 output for memory efficiency
  norm = "VC",                 # Vanilla Coverage normalization
  matrix = "observed",         # Extract observed (not expected) counts
  blockSize = 1e6,            # Memory management parameter
  onDisk = TRUE,              # Store results on disk (not in RAM)
  compressionLevel = 1         # HDF5 compression
)
```

### Output Dimensions
```
Dimensions: 5 x 5 x [N_loops] x 2
  = 5x5 pixels x [N] loops x 2 files (ctrl, mut)
```

**Data structure interpretation:**
- **5x5 matrices**: Local contact neighborhood around each loop
- **N_loops**: Total number of merged loops from Step 1
- **2 files**: Control and mutant samples for comparison

## Biological Data Validation

### Contact Matrix Quality Checks

#### 1. NA Value Assessment
```r
na_count <- sum(is.na(count_array))
na_percent <- 100 * na_count / length(count_array)
```
**Expected patterns:**
- Low NA percentage (< 10%) indicates good data quality
- High NAs may suggest sparse regions or technical issues

#### 2. Value Distribution Analysis
```r
# Overall distribution statistics
summary(non_na_values)

# Sparsity calculation
zeros <- sum(non_na_values == 0)
sparsity <- 100 * zeros / length(non_na_values)
```

**Biological interpretation:**
- **High sparsity (70-90% zeros)**: Normal for Hi-C data due to chromatin accessibility
- **Extreme values**: May indicate technical artifacts or very strong interactions

#### 3. Sample Correlation
```r
correlation <- cor(ctrl_values[complete_pairs],
                   mut_values[complete_pairs])
```

**Expected correlations:**
- **High correlation (> 0.7)**: Samples are biologically similar
- **Low correlation (< 0.5)**: May indicate significant biological differences or technical issues

### Per-Sample Statistics
```r
for (i in 1:dims[4]) {
  sample_values <- count_array[,,,i][!is.na(count_array[,,,i])]
  cat(sprintf("  %s: median=%.2f, mean=%.2f, max=%.0f\n",
              names(hicFiles)[i], median(sample_values),
              mean(sample_values), max(sample_values)))
}
```

**Key metrics:**
- **Median**: Typical contact strength (often low due to sparsity)
- **Mean**: Average contact strength (affected by high-value outliers)
- **Max**: Strongest interactions (should be reasonable, not extreme outliers)

## Example Matrix Visualization

### Individual Loop Matrices
```r
for (i in 1:min(3, dims[3])) {
  example_matrix <- count_array[, , i, 1]  # Loop i, control file

  # Display formatted matrix
  for (row in 1:5) {
    row_values <- round(example_matrix[row, ], 2)
    formatted_row <- format(row_values, width = 8, justify = "right")
    cat("  ", paste(formatted_row, collapse = " "), "\n")
  }
}
```

**Matrix interpretation:**
```
    12.5     8.2     15.4     6.1     3.2
     8.7    22.1     31.8    18.5     7.4
    15.9    31.2     45.6    29.3    14.1  ← Center row (loop focus)
     6.8    18.9     28.7    21.4     8.9
     3.1     7.6     13.8     9.2     4.5
                      ↑
                  Center pixel (loop center)
```

**Expected patterns:**
- **Central enrichment**: Highest values near center (2,2) position
- **Symmetric decay**: Values decrease with distance from center
- **Biological variation**: Some asymmetry is normal due to chromatin context

## Output Files

### 1. HDF5-backed InteractionArray (`05_extracted/`)
```r
saveHDF5SummarizedExperiment(
  x = pixels,
  dir = "outputs/full",
  prefix = "05_extracted",
  replace = TRUE
)
```

**Benefits of HDF5 format:**
- **Memory efficiency**: Handles large datasets without loading into RAM
- **Chunk-based access**: Can read specific loops/samples without loading all data
- **Cross-platform compatibility**: Standard format for large genomic datasets

### 2. Metadata File (`05_metadata.rds`)
```r
metadata <- list(
  n_loops = dims[3],
  n_files = dims[4],
  files = hicFiles,
  normalization = norm_to_use,
  extraction_time = extraction_time,
  dims = dims,
  correlation = correlation
)
```

## Biological Significance

### Contact Matrix Interpretation

**What each 5x5 matrix represents:**
1. **Central pixel (3,3)**: Direct loop contact strength
2. **Adjacent pixels**: Local chromatin context
3. **Edge pixels**: Background interaction levels
4. **Asymmetry**: Directional chromatin organization effects

### Quality Control Metrics

**Good extraction indicators:**
- Reasonable sample correlation (0.6-0.9)
- Low percentage of NA values (< 15%)
- Central enrichment in example matrices
- Consistent library sizes between samples

**Warning signs:**
- Very low correlation (< 0.4) - suggests technical problems
- High NA percentage (> 30%) - sparse or problematic regions
- Uniform matrices - loss of spatial structure

## Next Steps

The extracted contact matrices are ready for aggregation (Step 3: `aggregate.R`), where the 5x5 matrices will be collapsed into single values per loop for differential analysis.

## Technical Performance

### Memory Management
- **HDF5 backend**: Enables analysis of large datasets
- **Chunk processing**: Processes loops in batches
- **Compression**: Reduces storage requirements

### Typical Performance
- **Extraction time**: ~1-5 minutes for 150 loops
- **Memory usage**: < 1GB RAM due to on-disk storage
- **File size**: ~10-100MB depending on loop count and compression