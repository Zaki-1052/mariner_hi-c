# Mariner Usage in This Repository

## Overview

[Mariner](https://bioconductor.org/packages/mariner) is a Bioconductor package that provides a comprehensive framework for Hi-C chromatin loop analysis. This repository uses mariner extensively throughout the differential loop analysis pipeline for data structures, loop manipulation, and Hi-C matrix extraction.

**Core capabilities utilized:**
- Converting BEDPE loop calls to Bioconductor `GInteractions` objects
- Merging loop positions across biological replicates
- Binning loops to genomic resolution grids
- Creating buffered regions for robust quantification
- Extracting Hi-C contact matrices from `.hic` files
- Aggregate Peak Analysis (APA) for visualization

---

## 1. Data Structure Conversion

**Script**: `scripts/prep_loops.R` (lines 107-114)

### Function: `as_ginteractions()`

```r
gi <- as_ginteractions(
  bedpe_list[[name]],
  keep.extra.columns = TRUE,
  starts.in.df.are.0based = TRUE
)
```

**Purpose**: Convert HiCCUPS BEDPE files to Bioconductor `GInteractions` objects.

**Input**:
- BEDPE files from 6 replicates (3 control, 3 BAP1-KO mutant)
- Each file contains chromatin loop calls with metadata columns

**Output**:
- `GInteractions` object with paired genomic ranges
- Preserved metadata: observed counts, FDR values, centroid coordinates

**Key parameters**:
- `keep.extra.columns = TRUE`: Retains HiCCUPS statistics (observed, expectedBL, fdrBL, etc.)
- `starts.in.df.are.0based = TRUE`: Handles coordinate system conversion (BEDPE uses 0-based)

**What is GInteractions?**

A Bioconductor class for paired genomic ranges representing chromatin interactions:

```r
GInteractions object with 8500 interactions and 8 metadata columns:
        seqnames1   ranges1 seqnames2   ranges2 | observed fdrBL centroid1 centroid2
    [1]      chr1 1000-5000      chr1 2000-6000 |      125  1e-5    3000      4000
    [2]      chr1 3000-8000      chr2 1000-6000 |       89  2e-4    5500      3500
```

Each interaction has:
- **anchor1**: First genomic locus (e.g., promoter or enhancer)
- **anchor2**: Second genomic locus (e.g., enhancer or promoter)
- **metadata**: Observed counts, statistical significance, distance, etc.

---

## 2. Multi-Replicate Loop Merging

**Script**: `scripts/prep_loops.R` (lines 134-140)

### Function: `mergePairs()`

```r
merged <- mergePairs(
  x = gi_list,              # List of 6 GInteractions (one per replicate)
  radius = 10e3,            # 10kb merging radius
  column = "observed",      # Select based on Hi-C signal strength
  selectMax = TRUE,         # Keep loop with highest counts
  method = "manhattan"      # Distance calculation method
)
```

**Purpose**: Create consensus loop positions across 6 biological replicates.

**Strategy**: Union approach - includes loops detected in ANY replicate.

**Why merge?**
- Biological replicates often detect the same loop at slightly different coordinates
- HiCCUPS may call the same loop at bin positions that differ by 1-2 bins
- Merging creates a single consensus position for downstream analysis

**Parameters explained**:

| Parameter | Value | Purpose |
|-----------|-------|---------|
| `radius` | 10kb | Maximum distance between loop anchors to consider them the same loop |
| `column` | "observed" | Use Hi-C contact counts to select representative loop |
| `selectMax` | TRUE | Choose the loop with highest signal as the consensus position |
| `method` | "manhattan" | Calculate distance as sum of anchor1 distance + anchor2 distance |

**Example workflow**:

```
Replicate 1: chr1:1000000-1005000 <-> chr1:2000000-2005000  (observed=120)
Replicate 2: chr1:1005000-1010000 <-> chr1:2000000-2005000  (observed=150)
Replicate 3: chr1:1000000-1005000 <-> chr1:1995000-2000000  (observed=100)

→ Within 10kb radius, merged to single loop
→ Selects Rep2 position (highest observed=150)

Consensus: chr1:1005000-1010000 <-> chr1:2000000-2005000
```

**Typical reduction**:
- Input: ~15,000 loops across 6 replicates
- Output: ~8,500 consensus loop positions
- Reduction: ~40-45%

**Cluster metadata**:

The `mergePairs()` function attaches clustering information:

```r
clusters(merged)[[1]]
#   source       anchor1       anchor2
# 1 ctrl_M1 chr1:1000-5000 chr1:2000-6000
# 2 ctrl_M2 chr1:1005-5005 chr1:2000-6000
# 3 mut_M1  chr1:1000-5000 chr1:1995-5995
```

This shows which replicates contributed to each consensus loop.

---

## 3. Binning to Genomic Resolution

**Script**: `scripts/prep_loops.R` (lines 181-186)

### Function: `assignToBins()`

```r
binned <- assignToBins(
  x = merged,
  binSize = RESOLUTION,     # 5000, 10000, or 25000 bp
  pos1 = "center",          # Use loop centroid for anchor 1
  pos2 = "center"           # Use loop centroid for anchor 2
)
```

**Purpose**: Snap loop coordinates to a regular genomic grid matching `.hic` file bin structure.

**Why bin?**
- `.hic` files store Hi-C data in fixed-size bins (5kb, 10kb, 25kb, etc.)
- Loop coordinates must align with these bins for matrix extraction
- Ensures consistency between loop positions and Hi-C data

**How it works**:

```
Original merged loop: chr1:1,003,456-1,008,456 <-> chr1:2,012,789-2,017,789

At 5kb resolution:
  Anchor 1 center = 1,005,956
  → Binned to chr1:1,005,000-1,010,000 (bin covering center)

  Anchor 2 center = 2,015,289
  → Binned to chr1:2,015,000-2,020,000 (bin covering center)

Binned loop: chr1:1,005,000-1,010,000 <-> chr1:2,015,000-2,020,000
```

**pos1/pos2 options**:
- `"center"`: Use midpoint of loop region (recommended for merged loops with centroids)
- `"start"`: Use start coordinate
- `"end"`: Use end coordinate

---

## 4. Buffer Region Creation

**Script**: `scripts/prep_loops.R` (lines 198-201)

### Function: `pixelsToMatrices()`

```r
buffered <- pixelsToMatrices(
  x = binned,
  buffer = 2                # ±2 bins around center
)
```

**Purpose**: Expand single-pixel loop coordinates to rectangular regions (5×5 matrices at buffer=2).

**Why buffer?**

1. **Positional uncertainty**: Loops may shift by 1-2 bins between replicates
2. **Robustness**: Summing all pixels in a window is more stable than a single pixel
3. **Biological reality**: Chromatin loops are not point contacts but regions of enriched contact

**Buffer size calculation**:

```
buffer = 2
→ ±2 bins on each side of center
→ Total window = (2 + 1 + 2) = 5 bins

At 5kb resolution:
  5 bins × 5 kb/bin = 25 kb window per anchor
  Total region = 25kb × 25kb = 5×5 pixel matrix

At 10kb resolution:
  5 bins × 10 kb/bin = 50 kb window per anchor
  Total region = 50kb × 50kb = 5×5 pixel matrix
```

**Visual representation**:

```
     ←──────── 5 bins (25kb at 5kb res) ────────→
    ┌─────────────────────────────────────────────┐
    │  ·   ·   ·   ·   ·                          │  ↑
    │  ·   ·   ·   ·   ·   ← 2 bins above         │  │
    │  ·   ·   X   ·   ·   ← Expected loop pixel  │  5 bins
    │  ·   ·   ·   ·   ·   ← 2 bins below         │  │
    │  ·   ·   ·   ·   ·                          │  ↓
    └─────────────────────────────────────────────┘
         ↑       ↑       ↑
      2 bins  Center  2 bins
       left           right

  X = expected loop position (center pixel)
  · = flanking pixels (handle bin-shifts)
```

**Output coordinates**:

```r
# Before buffering (single pixel):
chr1:1,005,000-1,010,000 <-> chr1:2,015,000-2,020,000

# After buffering (buffer=2, 5kb resolution):
anchor1: chr1:995,000-1,020,000   (5 bins: 995-1000, 1000-1005, 1005-1010, 1010-1015, 1015-1020)
anchor2: chr1:2,005,000-2,030,000 (5 bins: 2005-2010, 2010-2015, 2015-2020, 2020-2025, 2025-2030)

→ Defines a 5×5 grid for matrix extraction
```

**Critical for aggregation strategy**:

The pipeline (in `scripts/aggregate.R`) sums all 25 pixels to get a robust count per loop:

```r
for each loop:
  Extract 5×5 matrix of Hi-C contacts
  Sum all 25 values → total_count
  Use total_count for edgeR differential analysis
```

This approach handles cases where the true loop signal is at pixel (3,2) instead of the expected center (3,3).

---

## 5. Hi-C Matrix Extraction

**Script**: `scripts/extract_counts.R` (lines 121-131)

### Function: `pullHicMatrices()`

```r
pixels <- pullHicMatrices(
  x = buffered,              # GInteractions with 5×5 regions
  files = hicFiles,          # Named vector of 6 .hic file paths
  binSize = RESOLUTION,      # 5000, 10000, or 25000
  h5File = h5_file_path,     # Path for HDF5 backend storage
  norm = "VC",               # Normalization method (VC or KR)
  matrix = "observed",       # Extract observed counts
  blockSize = 1e6,           # Process in 1Mb genomic chunks
  onDisk = TRUE,             # Store on disk, not in RAM
  compressionLevel = 1       # Light compression for speed
)
```

**Purpose**: Extract Hi-C contact matrices at each buffered loop position for all samples.

**What it extracts**:

For each loop in each sample, extract a 5×5 matrix of Hi-C contact values:

```
Loop: chr1:995,000-1,020,000 <-> chr1:2,005,000-2,030,000

Extract contacts between every pair of bins:

         Anchor 2 →
         2005-  2010-  2015-  2020-  2025-
         2010   2015   2020   2025   2030
        ┌─────────────────────────────────┐
 995-   │  12     8     15     7     4    │
 1000   │                                 │
        │                                 │
 1000-  │   9    14     18    11     6    │
 1005   │                                 │
        │                                 │
 1005-  │  11    16     25    14     8    │  ← Center row
 1010   │            ↑                    │
        │         Center                  │
 1010-  │   7    12     19    13     5    │
 1015   │                                 │
        │                                 │
 1015-  │   5     9     12     8     3    │
 1020   │                                 │
        └─────────────────────────────────┘

Each value = normalized Hi-C contact count
Center pixel (3,3) = 25 contacts (expected loop position)
```

**Output structure**: `InteractionArray` (4D HDF5-backed array)

```r
Dimensions: 5 × 5 × 8,500 × 6
  - Dimension 1: Anchor 1 bins (width)  = 5
  - Dimension 2: Anchor 2 bins (height) = 5
  - Dimension 3: Loop index             = 8,500 loops
  - Dimension 4: Sample index           = 6 samples

Total size: 5 × 5 × 8,500 × 6 = 1,275,000 individual matrices
```

**Memory management**: HDF5 backend

```r
# Without HDF5 (in-memory):
1,275,000 matrices × 8 bytes/value = ~10 GB RAM

# With HDF5 (on-disk):
Data stored in: outputs/res_5kb/temp_hdf5/extracted_matrices.h5
RAM usage: ~500 MB (only loads chunks as needed)
```

**Parameters explained**:

| Parameter | Value | Purpose |
|-----------|-------|---------|
| `binSize` | 5000 | Must match .hic file resolution |
| `norm` | "VC" | Vanilla Coverage normalization (adjusts for sequencing depth) |
| `matrix` | "observed" | Raw observed counts (not expected) |
| `blockSize` | 1e6 | Process 1Mb genomic chunks at a time (memory efficiency) |
| `onDisk` | TRUE | Use HDF5 file, don't load entire array into RAM |
| `compressionLevel` | 1 | Light compression (balance speed vs disk space) |

**Normalization options**:

- **VC (Vanilla Coverage)**: Adjusts for sequencing depth and mappability
- **KR (Knight-Ruiz)**: More sophisticated matrix balancing (not always available)
- **NONE**: Raw counts (not recommended)

This pipeline uses **VC** because KR normalization is not available in all `.hic` files.

**Example usage in pipeline**:

```r
# Step 1: Extract matrices (this step)
pixels <- pullHicMatrices(...)  # 5×5×8500×6 array

# Step 2: Aggregate to counts (aggregate.R)
for each loop i, sample j:
  mat_5x5 <- pixels[, , i, j]
  count[i, j] <- sum(mat_5x5)  # Sum all 25 pixels

# Result: 8500×6 count matrix
#         ctrl_M1  ctrl_M2  ctrl_M3  mut_M1  mut_M2  mut_M3
# loop_1    450      420      480      180     165     190
# loop_2    220      245      210      380     395     410
# ...

# Step 3: Differential analysis (edgeR.R)
dge <- DGEList(counts = count_matrix, ...)
```

---

## 6. Aggregate Peak Analysis (APA)

**Script**: `scripts/apa_analysis.R` (lines 168-180)

### Function: `pullHicMatrices()` (for APA)

```r
pixels <- pullHicMatrices(
  x = loops,
  files = hic_files,
  binSize = resolution,
  h5File = hdf5_file,
  half = "both",            # Extract both triangle halves (upper + lower)
  norm = norm,
  matrix = "observed",
  blockSize = 1e6,
  onDisk = TRUE,
  compressionLevel = 1
)
```

**Purpose**: Extract larger windows for Aggregate Peak Analysis (APA) visualization.

**Difference from count extraction**:

| Feature | Count Extraction | APA Extraction |
|---------|------------------|----------------|
| Window size | 5×5 (small, for summing) | 21×21 or larger (for visualization) |
| Buffer | 2 bins (±10-20kb) | 10+ bins (±50-100kb) |
| Purpose | Quantify loop signal for statistics | Visualize average loop structure |
| Aggregation | Sum across pixels | Mean/median across loops |

**APA workflow**:

```
1. Select significant differential loops (e.g., FDR < 0.05, |logFC| > 0.3)
2. Extract 21×21 matrices (50kb window at 5kb resolution) for each loop
3. Separate by direction: up-regulated vs down-regulated
4. Aggregate across loops: calculate mean matrix
5. Generate heatmaps comparing Control vs BAP1-KO
```

**Example APA matrix (aggregated across 500 up-regulated loops)**:

```
         ←────── 50kb window (21 bins at 5kb) ──────→
        -50kb                  0                  +50kb
         ┌───────────────────────────────────────────┐
         │   2   3   4   5   6   7   6   5   4   3  │
         │   3   4   5   7   9  11   9   7   5   4  │
         │   4   5   7   9  12  15  12   9   7   5  │
  -50kb  │   5   7   9  12  16  20  16  12   9   7  │
         │   6   9  12  16  22  28  22  16  12   9  │
         │   7  11  15  20  28  45  28  20  15  11  │ ← Peak at center
    0    │   6   9  12  16  22  28  22  16  12   9  │
         │   5   7   9  12  16  20  16  12   9   7  │
         │   4   5   7   9  12  15  12   9   7   5  │
  +50kb  │   3   4   5   7   9  11   9   7   5   4  │
         │   2   3   4   5   6   7   6   5   4   3  │
         └───────────────────────────────────────────┘

Center pixel (11,11) = 45 (mean contact enrichment)
Background (corners) = ~2-3 (local background)
P2LL enrichment = 45 / 2.5 = 18-fold
```

**APA output plots**:

1. **Control aggregate heatmap**: Average loop structure in wildtype
2. **BAP1-KO aggregate heatmap**: Average loop structure in mutant
3. **Difference heatmap**: Mutant - Control (shows changes)
4. **Enrichment boxplot**: Compares P2LL scores between conditions
5. **Per-replicate heatmaps**: Individual sample reproducibility check

**P2LL enrichment calculation**:

```r
# P2LL = Peak to Lower Left
center_value <- matrix[center_row, center_col]
corner_values <- c(matrix[1,1], matrix[1,21], matrix[21,1], matrix[21,21])
background <- mean(corner_values)

enrichment <- center_value / background
```

**Statistical comparison**:

```r
# Compare enrichment between conditions
ctrl_enrichment <- calculate_enrichment(ctrl_matrices)
mut_enrichment <- calculate_enrichment(mut_matrices)

wilcox_test <- wilcox.test(ctrl_enrichment, mut_enrichment)
# Tests: Do up-regulated loops have higher enrichment in BAP1-KO?
```

---

## Mariner Data Structures

### GInteractions

A Bioconductor class for paired genomic ranges representing chromatin interactions.

**Structure**:

```r
GInteractions object with 8500 interactions and 12 metadata columns:
        seqnames1    ranges1 seqnames2    ranges2 | observed expectedBL fdrBL distance
    [1]      chr1  1000-5000      chr1  2000-6000 |      125        45  1e-5   1000000
    [2]      chr1  3000-8000      chr2  1000-6000 |       89        62  2e-4        NA
    [3]      chr2 10000-15000     chr2 50000-55000 |      156        72  5e-6  40000000
```

**Access methods**:

```r
# Get first anchor
anchors(gi, type = "first")   # GRanges for anchor 1

# Get second anchor
anchors(gi, type = "second")  # GRanges for anchor 2

# Get metadata
mcols(gi)$observed            # Hi-C counts
mcols(gi)$fdrBL               # Statistical significance

# Subset
sig_loops <- gi[mcols(gi)$fdrBL < 0.05]
```

### InteractionArray (HDF5-backed)

4D array storing Hi-C contact matrices, backed by HDF5 for memory efficiency.

**Structure**:

```r
class: InteractionArray
dim: 5 5 8500 6
dimnames: NULL
assays(1): counts
backend: HDF5Array
seed: /path/to/extracted_matrices.h5

Dimensions explained:
  [1] = 5       (anchor 1 bins, width)
  [2] = 5       (anchor 2 bins, height)
  [3] = 8500    (number of loops)
  [4] = 6       (number of samples)
```

**Access methods**:

```r
# Extract counts array (4D DelayedArray)
count_array <- counts(pixels)

# Get specific matrix (loop i, sample j)
matrix_5x5 <- count_array[, , i, j]

# Aggregate across loops (mean matrix per sample)
sample1_mean <- apply(count_array[, , , 1], c(1, 2), mean)

# Load into memory (use sparingly!)
in_memory <- as.array(count_array)
```

**Why HDF5?**

```
Memory comparison for 8,500 loops × 6 samples:

In-memory array:
  5 × 5 × 8,500 × 6 = 1,275,000 values
  × 8 bytes (double precision)
  = 10.2 GB RAM

HDF5 on-disk:
  Same data stored in: extracted_matrices.h5
  RAM usage: ~500 MB (only active chunks)
  Disk usage: ~1-2 GB (with compression)

Pipeline can run on systems with 16GB RAM instead of requiring 64GB+
```

---

## Why Use Mariner?

### Advantages over manual implementation

| Task | Manual Approach | Mariner Approach |
|------|----------------|------------------|
| **BEDPE parsing** | Custom read.table + coordinate conversion | `as_ginteractions()` handles edge cases |
| **Loop merging** | Write clustering algorithm, handle edge cases | `mergePairs()` with customizable selection |
| **Binning** | Manual calculation, risk of off-by-one errors | `assignToBins()` consistent with .hic format |
| **Matrix extraction** | Direct strawr calls, manual chunking | `pullHicMatrices()` with HDF5 backend |
| **Memory management** | Manual HDF5 dataset creation | Automatic HDF5 integration |
| **Coordinate systems** | Track 0-based vs 1-based manually | Automatic conversion |

### Integration with Bioconductor ecosystem

Mariner objects work seamlessly with:

- **GenomicRanges**: Genomic coordinate manipulation
- **InteractionSet**: Interaction-specific operations
- **edgeR**: Differential analysis (via counts from matrices)
- **rtracklayer**: Import/export BED, BEDPE, bigWig files
- **strawr**: .hic file reading (mariner wrapper)
- **HDF5Array**: Out-of-memory array operations

### Replicate-aware design

```r
# Mariner handles replicates natively
mergePairs(
  x = list(rep1, rep2, rep3, rep4, rep5, rep6),  # Pass list
  ...
)
# Returns: Single consensus with cluster membership metadata

# Manual approach would require:
# 1. Nested loops comparing all pairs
# 2. Distance calculation
# 3. Greedy clustering algorithm
# 4. Selection of representative position
# 5. Bookkeeping of which replicates contributed
```

---

## Alternative Approaches (What if we didn't use mariner?)

### Direct strawr usage

```r
# Without mariner
library(strawr)

# For each loop, manually:
for (i in 1:nrow(loops)) {
  chr1 <- loops$chr1[i]
  start1 <- loops$start1[i] - 10000  # Manual buffer
  end1 <- loops$end1[i] + 10000

  chr2 <- loops$chr2[i]
  start2 <- loops$start2[i] - 10000
  end2 <- loops$end2[i] + 10000

  # Extract for each sample
  for (sample in hic_files) {
    matrix <- strawr::straw(
      norm = "VC",
      fname = sample,
      chr1loc = sprintf("%s:%d:%d", chr1, start1, end1),
      chr2loc = sprintf("%s:%d:%d", chr2, start2, end2),
      unit = "BP",
      binsize = resolution
    )

    # Convert sparse format to dense matrix
    # Handle missing bins
    # Store somewhere...
  }
}

# Problems:
# - 8,500 loops × 6 samples = 51,000 strawr calls
# - Manual memory management
# - No HDF5 integration
# - Error-prone coordinate arithmetic
```

### Manual HDF5 creation

```r
# Without mariner's HDF5 backend
library(rhdf5)

# Create file
h5createFile("matrices.h5")
h5createDataset("matrices.h5", "counts",
                dims = c(5, 5, 8500, 6),
                chunk = c(5, 5, 100, 1),  # Manual chunking
                compression = 1)

# Write data loop-by-loop
for (i in 1:8500) {
  for (j in 1:6) {
    # Extract matrix (strawr call)
    # Write to specific slice
    h5write(matrix, "matrices.h5", "counts",
            index = list(1:5, 1:5, i, j))
  }
}

# Problems:
# - Manual chunking strategy
# - No integration with Bioconductor objects
# - Have to track metadata separately
```

**Mariner handles all this automatically** with a single `pullHicMatrices()` call.

---

## Performance Characteristics

### Memory usage by pipeline step

| Step | Script | Peak RAM | Disk Output |
|------|--------|----------|-------------|
| BEDPE conversion | prep_loops.R | ~2 GB | 50 MB (RDS) |
| Loop merging | prep_loops.R | ~2 GB | 50 MB (RDS) |
| Matrix extraction | extract_counts.R | ~8 GB | 1.5 GB (HDF5) |
| Aggregation | aggregate.R | ~4 GB | 100 MB (TSV) |
| APA analysis | apa_analysis.R | ~16 GB | 500 MB (plots) |

### Runtime (16 CPUs, 128 GB RAM, HPC)

| Resolution | Extraction | Aggregation | Total per resolution |
|------------|-----------|-------------|---------------------|
| 5kb | 15-25 min | 2-4 min | ~30 min |
| 10kb | 10-20 min | 2-3 min | ~25 min |
| 25kb | 5-15 min | 1-2 min | ~20 min |

**Bottleneck**: Matrix extraction (I/O to .hic files)

**Optimization**: HDF5 onDisk storage prevents memory overflow at cost of disk I/O

---

## Common Troubleshooting

### "Cannot open .hic file"

**Cause**: File path incorrect or file not accessible

**Solution**:
```r
# Verify file exists
file.exists(hicFiles["ctrl_M1"])

# Check strawr can read it
strawr::readHicBpResolutions(hicFiles["ctrl_M1"])

# Common issues:
# - Relative path when absolute needed
# - Network mount not accessible
# - File permissions
```

### "Resolution not available in .hic file"

**Cause**: Requested resolution (e.g., 5000) not present in .hic file

**Solution**:
```r
# Check available resolutions
avail_res <- strawr::readHicBpResolutions(hicFiles["ctrl_M1"])
print(avail_res)
# [1]  2500  5000 10000 25000 50000 100000 250000 500000 1000000

# Use one of these resolutions
RESOLUTION <- 5000  # Must be in avail_res
```

### "HDF5 library not found"

**Cause**: HDF5Array package requires system HDF5 library

**Solution**:
```bash
# macOS
brew install hdf5

# Ubuntu/Debian
sudo apt-get install libhdf5-dev

# Then reinstall R package
R -e 'BiocManager::install("HDF5Array", force = TRUE)'
```

### "Out of memory during extraction"

**Cause**: Too many loops being processed at once, or blockSize too large

**Solution**:
```r
# In extract_counts.R, reduce blockSize
pixels <- pullHicMatrices(
  ...,
  blockSize = 5e5,  # Reduce from 1e6 (default 1Mb)
  ...
)

# Or process resolutions separately instead of in parallel
```

### "Dimensions do not match"

**Cause**: Mismatch between buffered regions and extracted matrices

**Solution**:
```r
# Ensure buffer matches extraction
buffer <- 2

# Create buffered regions
buffered <- pixelsToMatrices(x = binned, buffer = buffer)

# Extract with matching buffer
pixels <- pullHicMatrices(
  x = buffered,  # Use buffered, not binned
  ...
)
```

---

## Further Reading

**Mariner documentation**:
- Bioconductor: https://bioconductor.org/packages/mariner
- Vignettes: https://bioconductor.org/packages/release/bioc/vignettes/mariner/inst/doc/mariner.html
- GitHub: https://github.com/EricSDavis/mariner

**Related Bioconductor packages**:
- InteractionSet: https://bioconductor.org/packages/InteractionSet
- GenomicRanges: https://bioconductor.org/packages/GenomicRanges
- HDF5Array: https://bioconductor.org/packages/HDF5Array
- strawr: https://github.com/aidenlab/straw/tree/master/R

**Hi-C data formats**:
- .hic format specification: https://github.com/aidenlab/juicer/wiki/Data
- HiCCUPS loop calling: Rao et al. (2014) Cell 159:1665-1680

**Citation**:
- mariner: Kramer et al. (2022) mariner: Explore the Hi-Cs. Bioconductor.
