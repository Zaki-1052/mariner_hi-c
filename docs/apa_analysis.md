# Aggregate Peak Analysis (APA): Hi-C Signal Visualization

## Script: `apa_analysis.R`

## Overview

**Aggregate Peak Analysis (APA)** visualizes Hi-C contact enrichment by averaging contact matrices across many loops. This provides visual confirmation of differential looping and quantifies signal strength differences between BAP1-KO and wildtype conditions.

## What is APA?

### Concept

Individual Hi-C contact matrices are noisy, but averaging across many loops reveals consistent patterns:

**Single Loop (Noisy):**
```
   2  1  3  1  2
   1  5  2  1  1
   3  2 12  4  2    ← Center shows peak but unclear
   1  1  4  3  1
   2  1  2  1  2
```

**Aggregate of 100 Loops (Clear Pattern):**
```
  15 20 25 20 15
  20 35 45 35 20
  25 45 85 45 25    ← Clear enrichment at center
  20 35 45 35 20
  15 20 25 20 15
```

### Biological Insight

**High center enrichment (strong APA signal):**
- Loops are real biological structures
- Strong Hi-C contact frequency
- Stable chromatin interactions

**Low center enrichment (weak APA signal):**
- Loops may be technical artifacts
- Weak or transient interactions
- Poor loop calling quality

**Differential APA:**
- **Mut > Ctrl:** Stronger loops in BAP1-KO
- **Ctrl > Mut:** Weaker loops in BAP1-KO
- **Similar:** Loop position changes but not strength

## Purpose

This script performs APA on:
1. **Merged differential loops** - Non-redundant set across all resolutions
2. **Resolution-specific loops** - Loops unique to 5kb, 10kb, or 25kb
3. **Control vs. Mutant comparison** - Test for differential enrichment
4. **Statistical validation** - Quantify significance of differences

---

## Input Data

### Loop Sets

**Primary Input:**
```
outputs/merged_loops/non_redundant_loops.rds    # GInteractions, merged across resolutions
```

**Resolution-Specific (optional):**
```
outputs/res_5kb/02_merged.rds     # All 5kb loops
outputs/res_10kb/02_merged.rds    # All 10kb loops
outputs/res_25kb/02_merged.rds    # All 25kb loops
```

### Hi-C Contact Matrices

**All 6 .hic files:**
```
/expanse/lustre/projects/csd940/ctea/nf-hic/juicerpre/ctrl_M1.hic
/expanse/lustre/projects/csd940/ctea/nf-hic/juicerpre/ctrl_M2.hic
/expanse/lustre/projects/csd940/ctea/nf-hic/juicerpre/ctrl_M3.hic
/expanse/lustre/projects/csd940/ctea/nf-hic/juicerpre/mut_M1.hic
/expanse/lustre/projects/csd940/ctea/nf-hic/juicerpre/mut_M2.hic
/expanse/lustre/projects/csd940/ctea/nf-hic/juicerpre/mut_M3.hic
```

Requirements:
- Must exist and be accessible
- Must contain target resolution
- Should have normalization available (KR or VC)

---

## Analysis Workflow

### Step 1: Configuration

**Choose which analyses to run:**

```bash
# All resolutions, all loop sets (default)
Rscript scripts/apa_analysis.R

# Single resolution only
Rscript scripts/apa_analysis.R --resolution 5000

# Merged loops only (faster)
Rscript scripts/apa_analysis.R --loops merged

# Resolution-specific loops only
Rscript scripts/apa_analysis.R --loops resolution_specific
```

**Parameters:**
- **Resolutions:** 5000, 10000, 25000 (bp)
- **Window size:** 50kb (±25kb from loop center)
- **Normalization:** KR (preferred) or VC (fallback)

### Step 2: Validate .hic Files

```r
# Check file existence
for (file in HIC_FILES) {
  if (!file.exists(file)) {
    stop("Missing .hic file")
  }
}

# Check available normalizations
norms <- readHicNormTypes(HIC_FILES[1])
# Prefer: "KR" > "VC" > "NONE"
```

**KR normalization:**
- Knight-Ruiz matrix balancing
- Best for removing systematic biases
- Gold standard for Hi-C

**VC normalization:**
- Vanilla Coverage
- Simpler normalization
- Used when KR not available

### Step 3: Extract Hi-C Matrices

For each loop, extract a larger window for visualization:

**Buffer Calculation:**
```r
# For 50kb window at each resolution
buffer_5kb  = 50,000 / 5,000  = 10 bins  → 21×21 matrix
buffer_10kb = 50,000 / 10,000 = 5 bins   → 11×11 matrix
buffer_25kb = 50,000 / 25,000 = 2 bins   → 5×5 matrix
```

**Extraction:**
```r
pixels <- pullHicMatrices(
  bedpe = loops,              # GInteractions with loop positions
  files = HIC_FILES,          # All 6 .hic files
  binSize = resolution,       # 5000, 10000, or 25000
  h5File = hdf5_file,         # On-disk storage
  norm = "KR",                # Normalization method
  matrix = "observed",        # Observed counts
  blockSize = 1e6,            # Memory management
  onDisk = TRUE               # Use HDF5 for large datasets
)
```

**Output:** InteractionArray with dimensions:
```
[bins × bins × loops × samples]
[21 × 21 × 456 × 6]  # Example for 5kb, 456 loops
```

### Step 4: Aggregate Across Loops

Average matrices across all loops for each sample:

**Method 1: Mean (default)**
```r
# For each sample
for (j in 1:n_samples) {
  agg_matrix[, , j] = mean(count_array[, , , j])
}
```

**Method 2: Median (robust to outliers)**
```r
agg_matrix[, , j] = median(count_array[, , , j])
```

**Output:** One aggregated matrix per sample
```
[21 × 21] matrix for ctrl_M1
[21 × 21] matrix for ctrl_M2
[21 × 21] matrix for ctrl_M3
[21 × 21] matrix for mut_M1
[21 × 21] matrix for mut_M2
[21 × 21] matrix for mut_M3
```

### Step 5: Average by Condition

Combine replicates within each condition:

```r
# Control average (n=3)
ctrl_avg = mean(agg_matrices[, , 1:3])

# Mutant average (n=3)
mut_avg = mean(agg_matrices[, , 4:6])

# Difference
difference = mut_avg - ctrl_avg
```

**Result:** Three heatmaps per loop set
- Control aggregate
- Mutant aggregate
- Difference (mut - ctrl)

### Step 6: Calculate Enrichment Scores

**P2LL (Peak to Lower Left) Metric:**

```
P2LL = Center_Pixel / Background_Mean
```

**Center pixel:**
```
matrix[center_row, center_col]
```

**Background:**
```
Four corners of matrix:
[1, 1], [1, n], [n, 1], [n, n]

background = mean(corners)
```

**Interpretation:**
- **P2LL > 2:** Strong enrichment
- **P2LL = 1:** No enrichment (background level)
- **P2LL < 1:** Depletion (artifact)

**Calculate per loop per sample:**
```r
for (loop in 1:n_loops) {
  for (sample in 1:6) {
    center_val = matrix[center, center]
    bg = mean(matrix[corners])
    enrichment[loop, sample] = center_val / bg
  }
}
```

**Output:** Data frame with enrichment scores
```
loop_id  sample    group  center_value  background_mean  enrichment
loop_1   ctrl_M1   ctrl   245.3         82.1             2.99
loop_1   ctrl_M2   ctrl   238.7         79.4             3.01
loop_1   ctrl_M3   ctrl   251.2         85.3             2.94
loop_1   mut_M1    mut    312.8         90.2             3.47
loop_1   mut_M2    mut    305.1         88.7             3.44
loop_1   mut_M3    mut    318.9         92.1             3.46
```

### Step 7: Statistical Testing

**Compare enrichment between conditions:**

**Wilcoxon Rank-Sum Test (per loop set):**
```r
ctrl_enrichments <- enrichment_df$enrichment[enrichment_df$group == "ctrl"]
mut_enrichments <- enrichment_df$enrichment[enrichment_df$group == "mut"]

wilcox_result <- wilcox.test(
  mut_enrichments,
  ctrl_enrichments,
  alternative = "two.sided"
)
```

**Output:**
```
Wilcoxon rank sum test
data:  mut_enrichments and ctrl_enrichments
W = 145678, p-value = 2.3e-12
alternative hypothesis: true location shift is not equal to 0
```

**Interpretation:**
- **p < 0.001:** Highly significant difference
- **Median enrichment:** Quantify direction and magnitude
- **Effect size:** Calculate Cohen's d or similar

### Step 8: Generate Visualizations

**Heatmaps:**

```r
pheatmap(
  agg_matrix,
  cluster_rows = FALSE,      # Preserve spatial structure
  cluster_cols = FALSE,
  color = viridis(100),      # Perceptually uniform colors
  main = "Control APA Heatmap",
  show_rownames = FALSE,
  show_colnames = FALSE
)
```

**Violin Plots:**

```r
ggplot(enrichment_df, aes(x = group, y = enrichment, fill = group)) +
  geom_violin() +
  geom_boxplot(width = 0.1, fill = "white") +
  stat_compare_means(method = "wilcox.test") +
  labs(title = "P2LL Enrichment by Condition",
       y = "P2LL Score",
       x = "Condition")
```

**Difference Heatmap:**

Shows where changes occur (center vs. periphery):

```r
difference_matrix = mut_avg - ctrl_avg

pheatmap(
  difference_matrix,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "Difference: Mutant - Control",
  breaks = seq(-max(abs(difference_matrix)),
               max(abs(difference_matrix)),
               length.out = 101)
)
```

---

## Output Structure

```
outputs/apa_results/
├── res_5kb/
│   ├── merged_loops/
│   │   ├── apa_heatmap_ctrl.pdf              # Control average
│   │   ├── apa_heatmap_mut.pdf               # Mutant average
│   │   ├── apa_heatmap_difference.pdf        # Mut - Ctrl
│   │   ├── enrichment_scores.tsv             # Per-loop P2LL scores
│   │   ├── enrichment_comparison.pdf         # Violin plot ctrl vs. mut
│   │   └── statistical_test.txt              # Wilcoxon test results
│   │
│   └── resolution_specific/
│       ├── apa_heatmap_ctrl.pdf
│       ├── apa_heatmap_mut.pdf
│       ├── apa_heatmap_difference.pdf
│       ├── enrichment_scores.tsv
│       ├── enrichment_comparison.pdf
│       └── statistical_test.txt
│
├── res_10kb/
│   ├── merged_loops/
│   └── resolution_specific/
│
└── res_25kb/
    ├── merged_loops/
    └── resolution_specific/
```

---

## Interpreting APA Results

### Heatmap Patterns

**Strong APA (Good Quality):**
```
     Low         High
      |           |
    [Dark]     [Bright]
      ↓           ↓
  ▓▓░░░░░░░░░░░░▓▓
  ▓▓░░░░░░░░░░░░▓▓
  ░░░░░░░░░░░░░░░░
  ░░░░░░▓▓▓▓░░░░░░    ← Bright center spot
  ░░░░░░▓▓▓▓░░░░░░
  ░░░░░░▓▓▓▓░░░░░░
  ░░░░░░░░░░░░░░░░
  ▓▓░░░░░░░░░░░░▓▓
  ▓▓░░░░░░░░░░░░▓▓
```

Features:
- Sharp peak at center
- Gradual decay to edges
- Symmetric pattern
- High P2LL (> 2.0)

**Weak APA (Poor Quality):**
```
  ░░░░░░░░░░░░░░░░
  ░░░░░░░░░░░░░░░░
  ░░░░░░░░░░░░░░░░
  ░░░░░░░░░░░░░░░░    ← Uniform, no peak
  ░░░░░░░░░░░░░░░░
  ░░░░░░░░░░░░░░░░
  ░░░░░░░░░░░░░░░░
  ░░░░░░░░░░░░░░░░
  ░░░░░░░░░░░░░░░░
```

Features:
- No clear center peak
- Uniform signal
- Low P2LL (< 1.5)
- Possible artifacts

### Difference Heatmap Interpretation

**Increased in Mutant:**
```
Blue   White   Red
(Ctrl) (Equal) (Mut)
  ↓      ↓      ↓
  ▓▓░░░░░░░░░░░░▓▓
  ▓▓░░░░░░░░░░░░▓▓
  ░░░░░░░░░░░░░░░░
  ░░░░░░▒▒▒▒░░░░░░    ← Red center = Higher in mut
  ░░░░░░▒▒▒▒░░░░░░
  ░░░░░░▒▒▒▒░░░░░░
  ░░░░░░░░░░░░░░░░
  ▓▓░░░░░░░░░░░░▓▓
  ▓▓░░░░░░░░░░░░▓▓
```

**Interpretation:** Loops have stronger Hi-C signal in BAP1-KO
→ Compensatory loop strengthening
→ Increased chromatin accessibility
→ Structural reorganization

**Decreased in Mutant:**
```
  ░░░░░░░░░░░░░░░░
  ░░░░░░░░░░░░░░░░
  ░░░░░░░░░░░░░░░░
  ░░░░░░▓▓▓▓░░░░░░    ← Blue center = Lower in mut
  ░░░░░░▓▓▓▓░░░░░░
  ░░░░░░▓▓▓▓░░░░░░
  ░░░░░░░░░░░░░░░░
  ░░░░░░░░░░░░░░░░
  ░░░░░░░░░░░░░░░░
```

**Interpretation:** Loops have weaker Hi-C signal in BAP1-KO
→ Loss of loop stability
→ Reduced chromatin interaction
→ Direct BAP1-dependent loops

### Enrichment Score Statistics

**High P2LL (> 3.0):**
- Very strong loops
- Core regulatory interactions
- High confidence

**Moderate P2LL (2.0 - 3.0):**
- Good quality loops
- Typical for enhancer-promoter
- Validate with other data

**Low P2LL (1.5 - 2.0):**
- Weak loops
- May be transient
- Require validation

**Very Low P2LL (< 1.5):**
- Likely artifacts
- Poor loop calling
- Consider filtering out

---

## Performance and Optimization

### Expected Runtime

**Per resolution:**
- 5kb: ~2 hours (most loops, finest resolution)
- 10kb: ~1.5 hours
- 25kb: ~1 hour

**Total (all resolutions, both loop sets):** ~8-10 hours

### Memory Requirements

**Peak memory:** 16-32 GB RAM

Factors:
- Number of loops
- Matrix size (resolution-dependent)
- HDF5 caching

### Disk Space

**Temporary HDF5 files:** 5-10 GB per resolution

**Final outputs:** 100-500 MB total

### Optimization Tips

**1. Process fewer loops:**
```bash
# Only merged loops (skip resolution-specific)
Rscript scripts/apa_analysis.R --loops merged
```

**2. Process single resolution:**
```bash
# Just 10kb (fastest)
Rscript scripts/apa_analysis.R --resolution 10000
```

**3. Reduce window size:**
```r
# In script, change:
buffer_bp <- 50000  # To: 25000 (smaller window, faster)
```

**4. Use SLURM:**
```bash
# Run on HPC with more resources
sbatch --cpus-per-task=8 --mem=64G apa_job.sb
```

---

## Troubleshooting

### .hic Files Not Found

**Error:**
```
ERROR: Missing .hic file at /path/to/file.hic
```

**Solution:**
- Verify file paths in script (lines 63-69)
- Check if files accessible from compute node
- Update paths if data moved

### Normalization Not Available

**Error:**
```
KR normalization not found, falling back to VC
```

**Solution:**
- This is a warning, not error
- VC normalization will be used instead
- Results still valid but may have more bias

### Out of Memory

**Error:**
```
Error: cannot allocate vector of size X GB
```

**Solution:**
- Increase SLURM memory request: `#SBATCH --mem=64G`
- Process fewer loops: `--loops merged`
- Process single resolution
- Reduce buffer size

### Very Low Enrichment Scores

**Observation:**
```
Mean P2LL < 1.5 for all loops
```

**Possible causes:**
- Poor loop calling quality
- Wrong coordinates (off-by-one error)
- Normalization issues
- Wrong resolution

**Solutions:**
- Validate loops visually in Juicebox
- Check coordinate system (0-based vs. 1-based)
- Try different normalization
- Verify resolution matches BEDPE

### Heatmaps Look Uniform

**Observation:**
```
No visible center enrichment in APA heatmaps
```

**Possible causes:**
- Loops are weak/transient
- Too many false positive loops
- Wrong window size
- Aggregating incompatible loop types

**Solutions:**
- Filter to high-confidence loops only
- Separate by loop type (P-E, E-E, etc.)
- Increase window size
- Check individual loop matrices

---

## Biological Validation

### Cross-Check with edgeR Results

**Expectation:** APA changes should correlate with edgeR logFC

**Upregulated loops (logFC > 0):**
→ Should show higher enrichment in mutant APA

**Downregulated loops (logFC < 0):**
→ Should show higher enrichment in control APA

**Validation:**
```r
# Merge enrichment scores with edgeR results
merged_data <- left_join(
  enrichment_df,
  edger_results,
  by = "loop_id"
)

# Correlation test
cor.test(merged_data$enrichment_diff, merged_data$logFC)
# Expected: positive correlation, r > 0.3
```

### Resolution-Specific vs. Merged

**Comparison:**
- Merged loops: Should have strong APA (validated across resolutions)
- Resolution-specific: May be weaker (unique features)

**Expected:**
```
Median P2LL:
  Merged loops:           2.8
  5kb-specific:           2.3
  10kb-specific:          2.5
  25kb-specific:          2.1
```

### Biological Coherence

**Check if APA matches biology:**

**BAP1 is a deubiquitinase involved in:**
- Chromatin remodeling
- Transcriptional regulation
- DNA damage response

**Expected differential APA:**
- Changes at developmental genes
- Effects on enhancer-promoter loops
- Compensatory mechanisms

---

## Advanced Analyses (Optional)

### Stratify by Loop Type

Run APA separately for each regulatory type:

```r
# Separate loops by type
p_e_loops <- loops[loops$loop_type == "Promoter-Active_Enhancer"]
e_e_loops <- loops[loops$loop_type == "Active_Enhancer-Active_Enhancer"]

# Run APA on each
apa_p_e <- run_apa(p_e_loops, ...)
apa_e_e <- run_apa(e_e_loops, ...)

# Compare enrichment patterns
```

**Biological insight:** Different loop types may show different APA patterns

### Distance Stratification

Separate short-range from long-range loops:

```r
short_loops <- loops[loops$distance < 100000]   # < 100kb
long_loops <- loops[loops$distance >= 100000]   # >= 100kb
```

**Expected:** Short-range loops often have higher APA enrichment

### Time-Series or Multiple Conditions

If you have additional conditions:

```r
# Compare multiple conditions
conditions <- c("WT", "BAP1-KO", "Rescue")

for (cond in conditions) {
  apa_results[[cond]] <- run_apa(loops, hic_files[[cond]])
}
```

---

## Citation

If you use APA methodology, cite:

**APA method:**
- Rao et al. (2014). *A 3D map of the human genome at kilobase resolution reveals principles of chromatin looping*. Cell 159(7):1665-1680.

**P2LL enrichment metric:**
- Nora et al. (2017). *Targeted degradation of CTCF decouples local insulation of chromosome domains from genomic compartmentalization*. Cell 169(5):930-944.

**mariner implementation:**
- Davis et al. (2023). *mariner: Explore the Hi-Cs*. https://bioconductor.org/packages/mariner/

---

**Last Updated:** November 10, 2025
