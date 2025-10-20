# Loop Preparation and Merging (prep_loops.R)

## Biological Purpose

This script is the first step in analyzing chromatin loops from Hi-C data. It processes raw loop calls from **HiCCUPS** (a peak-calling algorithm for Hi-C loops) and creates consensus loop positions across **6 biological replicates** (3 control + 3 mutant) for comparative analysis.

**Key biological concepts:**
- **Chromatin loops**: 3D interactions between distant genomic regions that regulate gene expression
- **Loop calling**: Computational identification of statistically significant chromatin contacts
- **HiCCUPS**: Algorithm that identifies enriched contact frequencies as potential loops
- **Biological replicates**: Independent experiments that capture biological variability
- **Multi-resolution analysis**: Testing loops at different genomic bin sizes (5kb, 10kb, 25kb)

## Multi-Resolution Support

The script accepts resolution as a command-line argument:
```bash
Rscript scripts/prep_loops.R 5000   # 5kb resolution
Rscript scripts/prep_loops.R 10000  # 10kb resolution
Rscript scripts/prep_loops.R 25000  # 25kb resolution
```

**Resolution effects:**
- **5kb**: Highest resolution, most loops, captures short-range interactions
- **10kb**: Intermediate, balances sensitivity and noise
- **25kb**: Coarsest resolution, captures long-range and inter-TAD interactions

## Input Data

### BEDPE Files from HiCCUPS (Individual Replicates)

**Control replicates:**
```
/expanse/lustre/projects/csd940/ctea/nf-hic/juicer_frompre/hiccups_results/ctrl_M1/postprocessed_pixels_{RES}.bedpe
/expanse/lustre/projects/csd940/ctea/nf-hic/juicer_frompre/hiccups_results/ctrl_M2/postprocessed_pixels_{RES}.bedpe
/expanse/lustre/projects/csd940/ctea/nf-hic/juicer_frompre/hiccups_results/ctrl_M3/postprocessed_pixels_{RES}.bedpe
```

**Mutant replicates:**
```
/expanse/lustre/projects/csd940/ctea/nf-hic/juicer_frompre/hiccups_results/mut_M1/postprocessed_pixels_{RES}.bedpe
/expanse/lustre/projects/csd940/ctea/nf-hic/juicer_frompre/hiccups_results/mut_M2/postprocessed_pixels_{RES}.bedpe
/expanse/lustre/projects/csd940/ctea/nf-hic/juicer_frompre/hiccups_results/mut_M3/postprocessed_pixels_{RES}.bedpe
```

### BEDPE Format Structure
```
chr1  x1    x2    chr2  y1    y2    name  score  strand1  strand2  color
observed  expectedBL  expectedDonut  expectedH  expectedV
fdrBL  fdrDonut  fdrH  fdrV  numCollapsed  centroid1  centroid2  radius
```

**Key columns:**
- `chr1, x1, x2`: First anchor coordinates (genomic position)
- `chr2, y1, y2`: Second anchor coordinates (interaction partner)
- `observed`: Raw Hi-C contact count
- `expected*`: Background models for statistical testing
- `fdr*`: False discovery rates for different background models

## Processing Steps

### 1. Data Loading and Filtering (`read_bedpe` function)
```r
read_bedpe <- function(filepath, sample_name, expected_resolution) {
  # Load BEDPE with proper column names
  bedpe <- read.table(filepath, header = FALSE, skip = 2,
                      stringsAsFactors = FALSE, sep = "\t")

  # Filter for specified resolution
  bedpe$resolution <- bedpe$x2 - bedpe$x1
  bedpe_filtered <- bedpe[bedpe$resolution == expected_resolution, ]

  # Remove zero-count interactions
  bedpe_subset <- bedpe_filtered[bedpe_filtered$observed > 0, ]

  return(bedpe_subset)
}
```

**Filtering rationale:**
- **Resolution-specific**: Only keeps loops at the requested resolution
- **observed > 0**: Removes technical artifacts and focuses on real interactions
- **Per-replicate processing**: Maintains biological replicate structure

### 2. Conversion to GInteractions Objects
```r
gi_list <- lapply(names(bedpe_list), function(name) {
  as_ginteractions(bedpe_list[[name]],
                   keep.extra.columns = TRUE,
                   starts.in.df.are.0based = TRUE)
})
names(gi_list) <- c("ctrl_M1", "ctrl_M2", "ctrl_M3",
                    "mut_M1", "mut_M2", "mut_M3")
```

**GInteractions benefits:**
- Standardized format for genomic interval pairs
- Efficient storage and manipulation
- Integration with Bioconductor ecosystem
- Maintains replicate identity

### 3. Loop Merging Across All 6 Replicates (`mergePairs`)
```r
merged <- mergePairs(
  x = gi_list,              # All 6 replicates
  radius = 10e3,            # 10kb merging radius
  column = "observed",      # Merge based on contact strength
  selectMax = TRUE,         # Keep strongest contact when merging
  method = "manhattan"      # Distance metric for merging
)
```

**Merging rationale:**
- **10kb radius**: Accounts for positional differences between replicates/conditions
- **Union approach**: Includes loops from ANY replicate (maximizes discovery)
- **Manhattan distance**: Appropriate for genomic coordinates
- **selectMax**: Preserves the strongest signal when multiple loops overlap

**Biological significance:**
- Creates unified loop set for differential analysis
- Handles technical variation in loop calling between replicates
- Enables direct comparison of same loop across all samples

### 4. Cluster Size Analysis
```r
cluster_sizes <- sapply(clusters(merged), nrow)

# Report cluster size distribution
for (size in sort(unique(cluster_sizes))) {
  count <- sum(cluster_sizes == size)
  pct <- 100 * count / length(cluster_sizes)
  cat(sprintf("  %d replicate%s: %d loops (%.1f%%)\n",
              size, ifelse(size > 1, "s", ""), count, pct))
}
```

**Cluster interpretation:**
- **1 replicate**: Loop found in only one sample (condition/replicate-specific)
- **2-5 replicates**: Loop with partial support across samples
- **6 replicates**: Loop consistently found in all samples (high confidence)

**Quality indicators:**
- High overlap suggests consistent loop calling
- Condition-specific loops may indicate biological differences
- Low overall overlap may suggest technical issues

### 5. Per-Replicate Contribution Analysis
```r
for (name in names(gi_list)) {
  has_replicate <- sapply(clusters(merged), function(cluster) {
    name %in% cluster$source
  })
  count <- sum(has_replicate)
  pct <- 100 * count / length(merged)
  cat(sprintf("  %s: %d/%d loops (%.1f%%)\n",
              name, count, length(merged), pct))
}
```

**Biological interpretation:**
- Similar contributions suggest balanced replicate quality
- Low contribution may indicate poor quality replicate
- Asymmetric contributions between conditions suggest differential looping

### 6. Binning for Matrix Extraction
```r
binned <- assignToBins(
  x = merged,
  binSize = RESOLUTION,    # Resolution-specific binning
  pos1 = "center",         # Center-align first anchor
  pos2 = "center"          # Center-align second anchor
)
```

**Resolution-aware binning:**
- 5kb: Fine-grained positioning
- 10kb: Moderate precision
- 25kb: Coarse positioning for long-range loops

### 7. Buffer Creation for Contact Matrices
```r
buffered <- pixelsToMatrices(
  x = binned,
  buffer = 2              # Creates 5x5 pixel regions
)
```

**Buffer significance:**
- **5x5 pixel regions**: Captures local contact neighborhood around each loop
- Accounts for spatial uncertainty in loop positions
- Handles bin-shift artifacts between replicates
- Enables robust aggregation strategies

**Buffer size calculation:**
```
buffer = 2 bins
5kb resolution: ±10kb around center (25kb × 25kb region)
10kb resolution: ±20kb around center (50kb × 50kb region)
25kb resolution: ±50kb around center (125kb × 125kb region)
```

## Output Files

### Resolution-Specific Directory Structure
```
outputs/res_5kb/    # 5kb resolution outputs
outputs/res_10kb/   # 10kb resolution outputs
outputs/res_25kb/   # 25kb resolution outputs
```

### 1. GInteractions Objects (`01_ginteractions.rds`)
- Raw loop calls converted to standardized format
- List of 6 GInteractions objects (one per replicate)
- Maintains original loop positions and metadata

### 2. Merged Loops (`02_merged.rds`)
- Combined loop set with positional merging across all 6 replicates
- Handles overlapping calls between samples
- **MergedGInteractions class** with cluster information

### 3. Binned Loops (`03_binned.rds`)
- Loops aligned to resolution-specific genomic bins
- Ready for matrix extraction

### 4. Buffered Loops (`04_buffered.rds`)
- Final loop set with 5x5 pixel buffers
- **Critical output**: Input for contact matrix extraction
- Dimensions: N_loops × 2 anchors × 5 pixels

## Biological Interpretation

### Loop Statistics Example
```
Resolution: 5000 bp (5 kb)
Input: 105,432 loops across 6 replicates
Output: 22,108 consensus loop positions

Cluster size distribution:
  1 replicate: 8,234 loops (37.2%)
  2 replicates: 4,567 loops (20.7%)
  3 replicates: 3,211 loops (14.5%)
  4 replicates: 2,456 loops (11.1%)
  5 replicates: 1,890 loops (8.5%)
  6 replicates: 1,750 loops (7.9%)

Per-replicate contribution:
  ctrl_M1: 18,456/22,108 loops (83.5%)
  ctrl_M2: 18,123/22,108 loops (82.0%)
  ctrl_M3: 17,890/22,108 loops (80.9%)
  mut_M1: 18,678/22,108 loops (84.5%)
  mut_M2: 18,234/22,108 loops (82.5%)
  mut_M3: 17,678/22,108 loops (80.0%)
```

**Quality indicators:**
- **Similar replicate contributions**: Good technical reproducibility
- **High-confidence loops (6 replicates)**: ~8% of total, suitable for validation
- **Condition-specific loops**: May reveal biological differences

### Resolution Comparison

**Expected loop counts by resolution:**
- **5kb**: Most loops (~20-25k), short-range enriched
- **10kb**: Intermediate (~15-20k), balanced coverage
- **25kb**: Fewest loops (~10-15k), long-range enriched

## Next Steps

The buffered loops are ready for Hi-C contact matrix extraction (Step 2: `extract_counts.R`), where 5×5 contact matrices will be extracted from the original .hic files at each loop position for all 6 replicates.

## Technical Notes

### Memory Efficiency
- Uses RDS format for fast loading
- Processes replicates iteratively to manage memory
- HDF5 backend in later steps for large datasets

### Reproducibility
- Fixed parameters ensure consistent merging
- Resolution parameter documented in outputs
- Cluster information preserved for validation

### Quality Control
- Removes zero-count and wrong-resolution loops
- Handles coordinate system conversions (0-based → 1-based)
- Validates file existence before processing

### Multi-Resolution Workflow
```bash
# Process all three resolutions
for RES in 5000 10000 25000; do
  Rscript scripts/prep_loops.R ${RES}
done

# Outputs stored in separate directories:
# outputs/res_5kb/
# outputs/res_10kb/
# outputs/res_25kb/
```

## Performance Metrics

**Typical runtime:**
- 5kb resolution: ~2-3 minutes
- 10kb resolution: ~1-2 minutes
- 25kb resolution: ~30-60 seconds

**Memory usage:** ~2-4 GB RAM

**Disk space:** ~100-500 MB per resolution
