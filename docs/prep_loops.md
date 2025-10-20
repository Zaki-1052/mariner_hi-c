# Loop Preparation and Merging (prep_loops.R)

## Biological Purpose

This script is the first step in analyzing chromatin loops from Hi-C data. It processes raw loop calls from **HiCCUPS** (a peak-calling algorithm for Hi-C loops) and prepares them for comparative analysis between control and mutant samples.

**Key biological concepts:**
- **Chromatin loops**: 3D interactions between distant genomic regions that regulate gene expression
- **Loop calling**: Computational identification of statistically significant chromatin contacts
- **HiCCUPS**: Algorithm that identifies enriched contact frequencies as potential loops
- **Resolution**: 5kb bins represent the spatial resolution of the Hi-C analysis

## Input Data

### BEDPE Files from HiCCUPS
- **Control**: `/expanse/lustre/projects/csd940/ctea/nf-hic/hiccups/hiccups_results/resorted_ctrl/merged_loops.bedpe`
- **Mutant**: `/expanse/lustre/projects/csd940/ctea/nf-hic/hiccups/hiccups_results/resorted_mut/merged_loops.bedpe`

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
read_bedpe <- function(filepath) {
  # Load BEDPE with proper column names
  bedpe <- read.table(filepath, header = FALSE, skip = 2,
                      stringsAsFactors = FALSE, sep = "\t")

  # Calculate resolution from coordinates
  bedpe$resolution <- bedpe$x2 - bedpe$x1

  # Filter for 5kb resolution only
  bedpe_5kb <- bedpe[bedpe$resolution == 5000, ]

  # Remove zero-count interactions
  bedpe_subset <- bedpe_subset[bedpe_subset$observed > 0, ]

  return(bedpe_subset)
}
```

**Filtering rationale:**
- **5kb resolution**: Standard resolution for loop analysis (balances sensitivity vs. noise)
- **observed > 0**: Removes technical artifacts and focuses on real interactions

### 2. Conversion to GInteractions Objects
```r
gi_ctrl <- as_ginteractions(ctrlLoops,
                           keep.extra.columns = TRUE,
                           starts.in.df.are.0based = TRUE)
```

**GInteractions benefits:**
- Standardized format for genomic interval pairs
- Efficient storage and manipulation
- Integration with Bioconductor ecosystem

### 3. Loop Merging (`mergePairs`)
```r
merged <- mergePairs(
  x = gi_list,
  radius = 10e3,           # 10kb merging radius
  column = "observed",     # Merge based on contact strength
  selectMax = TRUE,        # Keep strongest contact when merging
  method = "manhattan"     # Distance metric for merging
)
```

**Merging rationale:**
- **10kb radius**: Accounts for slight positional differences between samples
- **Manhattan distance**: Appropriate for genomic coordinates
- **selectMax**: Preserves the strongest signal when multiple loops overlap

### 4. Binning for Matrix Extraction
```r
binned <- assignToBins(
  x = merged,
  binSize = 5e3,          # 5kb bins
  pos1 = "center",        # Center-align first anchor
  pos2 = "center"         # Center-align second anchor
)
```

### 5. Buffer Creation for Contact Matrices
```r
buffered <- pixelsToMatrices(
  x = binned,
  buffer = 2              # Creates 5x5 pixel regions
)
```

**Buffer significance:**
- **5x5 pixel regions**: Captures local contact neighborhood around each loop
- Accounts for spatial uncertainty in loop positions
- Enables aggregation strategies to handle bin-shift artifacts

## Output Files

### 1. GInteractions Objects (`01_ginteractions.rds`)
- Raw loop calls converted to standardized format
- Separate objects for control and mutant samples

### 2. Merged Loops (`02_merged.rds`)
- Combined loop set with positional merging
- Handles overlapping calls between samples

### 3. Binned Loops (`03_binned.rds`)
- Loops aligned to 5kb genomic bins
- Ready for matrix extraction

### 4. Buffered Loops (`04_buffered.rds`)
- Final loop set with 5x5 pixel buffers
- **Critical output**: Input for contact matrix extraction

## Biological Interpretation

### Loop Statistics
```
ctrl: [N] loops, mut: [M] loops
merged from [N+M] input → [merged_count] final loops
clusters: [single] single, [paired] paired
```

**Cluster analysis:**
- **Single clusters**: Loops found in only one sample
- **Paired clusters**: Loops found in both samples (key for differential analysis)

### Quality Metrics
- **Merge efficiency**: High overlap suggests consistent loop calling
- **Resolution consistency**: All loops at 5kb resolution ensures comparable analysis
- **Contact strength**: Positive observed counts indicate real interactions

## Next Steps

The buffered loops are ready for Hi-C contact matrix extraction (Step 2: `extract_counts.R`), where 5x5 contact matrices will be extracted from the original .hic files at each loop position.

## Technical Notes

- **Memory efficiency**: Uses RDS format for fast loading
- **Reproducibility**: Fixed parameters ensure consistent merging
- **Quality control**: Removes zero-count and wrong-resolution loops
- **Coordinate system**: Handles 0-based vs 1-based coordinate conversions