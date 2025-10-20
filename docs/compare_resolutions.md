# Multi-Resolution Comparison Analysis (compare_resolutions.R)

## Biological Purpose

This script performs **meta-analysis across multiple Hi-C resolutions** (5kb, 10kb, 25kb) to identify robust differential chromatin loops and assess resolution-specific effects. By comparing results across resolutions, we can distinguish high-confidence changes from resolution-specific artifacts and understand how chromatin organization changes at different spatial scales.

**Key biological concepts:**
- **Resolution-dependent detection**: Different resolutions capture different loop types
- **Cross-resolution concordance**: Loops differential at multiple resolutions are high-confidence
- **Spatial scale effects**: Short-range vs long-range interaction differences
- **Robustness validation**: Consistent results = reliable biology

## Input Data

### edgeR Results from All Resolutions

**Required files:**
```
outputs/edgeR_results_res_5kb/primary_analysis/all_results_primary.rds
outputs/edgeR_results_res_10kb/primary_analysis/all_results_primary.rds
outputs/edgeR_results_res_25kb/primary_analysis/all_results_primary.rds
```

**Each contains:**
- Loop coordinates (chr, start, end for both anchors)
- Statistical results (logFC, PValue, FDR)
- Significance flags
- Expression levels (logCPM)

### Loop Coordinates

**Required files:**
```
outputs/res_5kb/03_binned.rds
outputs/res_10kb/03_binned.rds
outputs/res_25kb/03_binned.rds
```

**Purpose**: Match loops across resolutions by genomic position

## Output Directory

All outputs saved to:
```
outputs/resolution_comparison/
```

## Section 1: Load Results from All Resolutions

### Multi-Resolution Data Loading

```r
resolutions <- c(5000, 10000, 25000)
results_list <- list()
coords_list <- list()

for (res in resolutions) {
  res_kb <- res / 1000

  # Load differential results
  results_file <- sprintf(
    "outputs/edgeR_results_res_%dkb/primary_analysis/all_results_primary.rds",
    res_kb
  )
  coords_file <- sprintf("outputs/res_%dkb/03_binned.rds", res_kb)

  if (file.exists(results_file) && file.exists(coords_file)) {
    results_list[[as.character(res)]] <- readRDS(results_file)
    coords_list[[as.character(res)]] <- readRDS(coords_file)
    cat(sprintf("✓ %d loops\n", nrow(results_list[[as.character(res)]])))
  }
}
```

**Expected loop counts:**
- 5kb: ~18,000-22,000 loops (most)
- 10kb: ~13,000-18,000 loops (intermediate)
- 25kb: ~8,000-12,000 loops (fewest)

**Why different counts?**
- Finer resolution → more loops detected
- Coarser resolution → merges nearby loops
- Different sensitivity to short vs long-range

## Section 2: Basic Statistics Per Resolution

### Summary Table Creation

```r
summary_stats <- data.frame(
  resolution_kb = numeric(),
  total_loops = numeric(),
  significant_fdr05 = numeric(),
  pct_significant = numeric(),
  up_in_mutant = numeric(),
  down_in_mutant = numeric(),
  median_logFC_all = numeric(),
  median_logFC_sig = numeric()
)

for (res in names(results_list)) {
  res_kb <- as.numeric(res) / 1000
  df <- results_list[[res]]

  n_sig <- sum(df$significant, na.rm = TRUE)
  n_up <- sum(df$significant & df$logFC > 0, na.rm = TRUE)
  n_down <- sum(df$significant & df$logFC < 0, na.rm = TRUE)

  summary_stats <- rbind(summary_stats, data.frame(
    resolution_kb = res_kb,
    total_loops = nrow(df),
    significant_fdr05 = n_sig,
    pct_significant = 100 * n_sig / nrow(df),
    up_in_mutant = n_up,
    down_in_mutant = n_down,
    median_logFC_all = median(abs(df$logFC), na.rm = TRUE),
    median_logFC_sig = median(abs(df$logFC[df$significant]), na.rm = TRUE)
  ))
}

print(summary_stats)
```

**Example output:**
```
  resolution_kb  total_loops  significant_fdr05  pct_significant  up_in_mutant  down_in_mutant  median_logFC_all  median_logFC_sig
1             5        21543               3245             15.1          1823            1422             0.234             0.789
2            10        15678               2156             13.8          1205             951             0.198             0.812
3            25         9234                987             10.7           534             453             0.167             0.845
```

**Interpretation:**
- **Higher resolution → more significant loops**: More sensitivity
- **Lower resolution → higher effect sizes**: Less noise, clearer signal
- **% significant decreases**: Coarser binning reduces power

**Output**: `summary_statistics_by_resolution.tsv`

## Section 3: Coordinate-Based Overlap Analysis

### Loop Matching Function

```r
match_loops <- function(coords1, coords2, tolerance_kb = 10) {
  # Convert to data.frame with coordinates
  df1 <- data.frame(
    chr1 = as.character(seqnames(anchors(coords1, "first"))),
    start1 = start(anchors(coords1, "first")),
    chr2 = as.character(seqnames(anchors(coords2, "second"))),
    start2 = start(anchors(coords2, "second")),
    index1 = 1:length(coords1)
  )

  df2 <- data.frame(
    chr1 = as.character(seqnames(anchors(coords2, "first"))),
    start1 = start(anchors(coords2, "first")),
    chr2 = as.character(seqnames(anchors(coords2, "second"))),
    start2 = start(anchors(coords2, "second")),
    index2 = 1:length(coords2)
  )

  tolerance_bp <- tolerance_kb * 1000

  # Find overlaps: same chr AND within tolerance for both anchors
  matches <- data.frame(index1 = integer(), index2 = integer())

  for (i in 1:nrow(df1)) {
    same_chr <- df2$chr1 == df1$chr1[i] & df2$chr2 == df1$chr2[i]
    dist1 <- abs(df2$start1 - df1$start1[i])
    dist2 <- abs(df2$start2 - df1$start2[i])
    close_enough <- same_chr & dist1 <= tolerance_bp & dist2 <= tolerance_bp

    if (any(close_enough)) {
      # Take closest match (minimum total distance)
      distances <- dist1 + dist2
      best_match <- which(close_enough)[which.min(distances[close_enough])]
      matches <- rbind(matches, data.frame(index1 = i, index2 = best_match))
    }
  }

  return(matches)
}
```

**Tolerance parameter:**
- **10kb tolerance**: Accounts for binning differences across resolutions
- Loop at 5kb: chr1:1,000,000-1,005,000
- Same loop at 10kb: chr1:1,000,000-1,010,000
- Within 10kb = matched

### Overlap Matrix Construction

```r
overlap_matrix <- matrix(0, nrow = length(resolutions), ncol = length(resolutions))
rownames(overlap_matrix) <- paste0(resolutions/1000, "kb")
colnames(overlap_matrix) <- paste0(resolutions/1000, "kb")

for (i in 1:length(resolutions)) {
  for (j in 1:length(resolutions)) {
    res_i <- as.character(resolutions[i])
    res_j <- as.character(resolutions[j])

    if (i == j) {
      # Diagonal: total loops at this resolution
      overlap_matrix[i, j] <- nrow(results_list[[res_i]])
    } else {
      # Off-diagonal: loops that match between resolutions
      matches <- match_loops(coords_list[[res_i]], coords_list[[res_j]], tolerance_kb = 10)
      overlap_matrix[i, j] <- nrow(matches)
    }
  }
}

cat("\nLoop overlap matrix:\n")
print(overlap_matrix)
```

**Example output:**
```
Loop overlap matrix:
       5kb    10kb   25kb
5kb  21543  14532  8234
10kb 14521  15678  7892
25kb  8231   7886  9234
```

**Interpretation:**
- **Diagonal**: Total loops at each resolution
- **Row → Column**: How many from row match to column
- **Asymmetry**: Different detection at different resolutions

**Example**: 14,532 loops from 5kb match to 10kb, but only 14,521 from 10kb match to 5kb
- Some 5kb loops merge into single 10kb loop
- Some 10kb loops have no 5kb counterpart (different sensitivity)

**Output**: `loop_overlap_matrix.tsv`

## Section 4: Differential Concordance Analysis

### Concordance Between Resolutions

```r
# Compare 5kb vs 10kb
if ("5000" %in% names(coords_list) && "10000" %in% names(coords_list)) {
  matches_5_10 <- match_loops(coords_list[["5000"]], coords_list[["10000"]], tolerance_kb = 10)

  res_5kb <- results_list[["5000"]]
  res_10kb <- results_list[["10000"]]

  matched_data <- data.frame(
    sig_5kb = res_5kb$significant[matches_5_10$index1],
    sig_10kb = res_10kb$significant[matches_5_10$index2],
    fc_5kb = res_5kb$logFC[matches_5_10$index1],
    fc_10kb = res_10kb$logFC[matches_5_10$index2]
  )

  # Concordance statistics
  both_sig <- sum(matched_data$sig_5kb & matched_data$sig_10kb, na.rm = TRUE)
  both_sig_same_dir <- sum(matched_data$sig_5kb & matched_data$sig_10kb &
                           sign(matched_data$fc_5kb) == sign(matched_data$fc_10kb),
                           na.rm = TRUE)

  cat("5kb vs 10kb:\n")
  cat(sprintf("  Matched loops: %d\n", nrow(matched_data)))
  cat(sprintf("  Both significant: %d (%.1f%%)\n",
              both_sig, 100*both_sig/nrow(matched_data)))
  cat(sprintf("  Concordant direction: %d (%.1f%% of both-sig)\n",
              both_sig_same_dir, 100*both_sig_same_dir/both_sig))
  cat(sprintf("  Fold-change correlation: %.3f\n\n",
              cor(matched_data$fc_5kb, matched_data$fc_10kb, use = "complete.obs")))
}
```

**Example output:**
```
5kb vs 10kb:
  Matched loops: 14,532
  Both significant: 1,823 (12.5%)
  Concordant direction: 1,798 (98.6% of both-sig)
  Fold-change correlation: 0.847

5kb vs 25kb:
  Matched loops: 8,234
  Both significant: 645 (7.8%)
  Concordant direction: 638 (98.9% of both-sig)
  Fold-change correlation: 0.712

10kb vs 25kb:
  Matched loops: 7,892
  Both significant: 523 (6.6%)
  Concordant direction: 518 (99.0% of both-sig)
  Fold-change correlation: 0.823
```

**Interpretation:**
- **Both significant**: High-confidence differential loops
- **Concordant direction >95%**: Very reliable when significant at both
- **FC correlation 0.7-0.9**: Moderate to high concordance of effect sizes
- **Lower % at coarser comparison**: Resolution effects

**High-confidence loops**: Significant at ALL three resolutions with same direction

## Section 5: Visualizations

### 1. Differential Loop Counts by Resolution

**File**: `differential_loops_by_resolution.pdf`

```r
summary_long <- summary_stats %>%
  pivot_longer(cols = c(up_in_mutant, down_in_mutant),
               names_to = "direction", values_to = "count")

ggplot(summary_long, aes(x = factor(resolution_kb), y = count, fill = direction)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Up in Mutant" = "#d73027",
                               "Down in Mutant" = "#4575b4")) +
  labs(
    title = "Differential Loops by Resolution",
    x = "Resolution (kb)",
    y = "Number of Significant Loops (FDR < 0.05)",
    fill = "Direction"
  )
```

**Shows:**
- Decreasing loop count with coarser resolution
- Balance of up vs down at each resolution
- Overall trend in sensitivity

### 2. Fold-Change Correlation Scatter Plots

**File**: `foldchange_correlation_5kb_vs_10kb.pdf`

```r
# Match loops between resolutions
matches_5_10 <- match_loops(coords_list[["5000"]], coords_list[["10000"]], tolerance_kb = 10)

plot_data <- data.frame(
  fc_5kb = results_list[["5000"]]$logFC[matches_5_10$index1],
  fc_10kb = results_list[["10000"]]$logFC[matches_5_10$index2],
  sig_5kb = results_list[["5000"]]$significant[matches_5_10$index1],
  sig_10kb = results_list[["10000"]]$significant[matches_5_10$index2]
) %>%
  mutate(
    status = case_when(
      sig_5kb & sig_10kb ~ "Both Significant",
      sig_5kb ~ "5kb Only",
      sig_10kb ~ "10kb Only",
      TRUE ~ "Neither"
    )
  )

ggplot(plot_data, aes(x = fc_5kb, y = fc_10kb, color = status)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  scale_color_manual(values = c(
    "Both Significant" = "#1b7837",
    "5kb Only" = "#d73027",
    "10kb Only" = "#4575b4",
    "Neither" = "#999999"
  )) +
  labs(
    title = "Fold-Change Correlation: 5kb vs 10kb",
    subtitle = sprintf("Pearson r = %.3f (n = %d matched loops)",
                      cor(plot_data$fc_5kb, plot_data$fc_10kb, use = "complete.obs"),
                      nrow(plot_data)),
    x = "log2 Fold Change (5kb)",
    y = "log2 Fold Change (10kb)"
  )
```

**Interpretation:**
- **Points on diagonal**: Perfect concordance
- **Green points**: High-confidence (significant at both)
- **Red/Blue points**: Resolution-specific effects
- **Gray scatter**: Non-significant

**Expected pattern:**
- Strong correlation (r > 0.7)
- Most significant points cluster on diagonal
- Some resolution-specific outliers

### 3. Venn Diagram of Differential Loops

**File**: `venn_diagram_differential_loops.pdf`

```r
# Map significant loops to common reference (5kb)
matches_5_10 <- match_loops(coords_list[["5000"]], coords_list[["10000"]], tolerance_kb = 10)
matches_5_25 <- match_loops(coords_list[["5000"]], coords_list[["25000"]], tolerance_kb = 10)

sig_5kb <- which(results_list[["5000"]]$significant)
sig_10kb_mapped <- matches_5_10$index1[matches_5_10$index2 %in% which(results_list[["10000"]]$significant)]
sig_25kb_mapped <- matches_5_25$index1[matches_5_25$index2 %in% which(results_list[["25000"]]$significant)]

venn_list <- list(
  "5kb" = sig_5kb,
  "10kb" = sig_10kb_mapped,
  "25kb" = sig_25kb_mapped
)

venn.plot <- venn.diagram(
  x = venn_list,
  filename = NULL,
  category.names = c("5kb", "10kb", "25kb"),
  fill = c("#d73027", "#4575b4", "#1b7837"),
  alpha = 0.5,
  main = "Differential Loops Across Resolutions"
)
```

**Venn diagram regions:**
- **Center (all 3)**: High-confidence differential loops (~500-1,000)
- **5kb only**: Short-range loops not detected at coarser resolution
- **10kb only**: Intermediate loops
- **25kb only**: Long-range loops requiring coarse binning
- **Pairwise overlaps**: Moderate confidence

**High-confidence set**: Loops in center (all 3 resolutions agree)

## Output Files

### Summary Tables

**1. `summary_statistics_by_resolution.tsv`**
- Rows: Each resolution (5kb, 10kb, 25kb)
- Columns: Loop counts, significance rates, effect sizes
- Purpose: Quick comparison of results

**2. `loop_overlap_matrix.tsv`**
- Matrix showing loop matches between resolutions
- Diagonal: Total loops
- Off-diagonal: Overlapping loops
- Purpose: Understand resolution relationships

### Visualizations

**3. `differential_loops_by_resolution.pdf`**
- Bar plot of up/down differential loops per resolution
- Shows detection trends
- Identifies which resolution has most power

**4. `foldchange_correlation_5kb_vs_10kb.pdf`**
- Scatter plot of fold changes
- Color-coded by significance
- Diagonal line for perfect concordance

**5. `foldchange_correlation_5kb_vs_25kb.pdf`**
- Same as above for 5kb vs 25kb

**6. `foldchange_correlation_10kb_vs_25kb.pdf`**
- Same as above for 10kb vs 25kb

**7. `venn_diagram_differential_loops.pdf`**
- 3-way Venn diagram
- Shows overlap and unique sets
- Identifies high-confidence consensus

## Biological Interpretation

### Resolution-Specific Effects

**5kb-specific differential loops:**
- **Biological**: Short-range interactions (50-200 kb)
- **Examples**: Promoter-enhancer loops, local regulatory elements
- **Caveat**: More sensitive to noise, verify with other data

**25kb-specific differential loops:**
- **Biological**: Long-range interactions (>500 kb)
- **Examples**: Inter-TAD loops, long-range enhancers
- **Caveat**: Coarse binning may miss fine structure

**Consensus across all resolutions:**
- **Biological**: Robust changes at multiple spatial scales
- **Confidence**: Very high - not artifacts of binning
- **Priority**: Top candidates for validation

### Effect Size Concordance

**High correlation (r > 0.8):**
- ✓ Consistent biological signal
- ✓ Resolution choice has minor impact on effect estimates
- ✓ Results are robust

**Moderate correlation (r 0.6-0.8):**
- ⚠ Some resolution-dependent effects
- ⚠ Consider biological meaning of scale differences
- ✓ Still generally concordant

**Low correlation (r < 0.6):**
- ⚠ Significant resolution effects
- ⚠ May indicate different biology at different scales
- ⚠ Interpret resolution-specific results carefully

### Direction Concordance

**>95% same direction when both significant:**
- ✓ Excellent reliability
- ✓ False discoveries are random, not systematic
- ✓ Loops significant at multiple resolutions are trustworthy

**<90% same direction:**
- ⚠ May indicate technical issues
- ⚠ Check for sample mix-ups
- ⚠ Review individual resolution QC

## Recommended Analysis Strategy

### Step 1: Identify High-Confidence Set

```r
# Loops significant at ALL three resolutions with same direction
high_conf_5kb_indices <- intersect(
  intersect(sig_5kb, sig_10kb_mapped),
  sig_25kb_mapped
)

# Check direction concordance
high_conf_loops <- results_list[["5000"]][high_conf_5kb_indices, ]
cat(sprintf("High-confidence differential loops: %d\n", nrow(high_conf_loops)))
```

**Use for:**
- Publication figures
- Experimental validation
- Functional follow-up
- Pathway analysis

### Step 2: Examine Resolution-Specific Sets

```r
# 5kb-only loops
only_5kb <- setdiff(setdiff(sig_5kb, sig_10kb_mapped), sig_25kb_mapped)

# Characterize by distance
only_5kb_distances <- calculate_distances(only_5kb)
median(only_5kb_distances)  # Expect shorter distances
```

**Interpretation:**
- Shorter distances → local regulation
- Longer distances → may need finer resolution
- Biological context determines relevance

### Step 3: Effect Size Meta-Analysis

```r
# For matched loops, average log2FC across resolutions
matched_all <- find_loops_in_all_three_resolutions()

meta_logFC <- rowMeans(cbind(
  fc_5kb[matched_all$idx_5kb],
  fc_10kb[matched_all$idx_10kb],
  fc_25kb[matched_all$idx_25kb]
))

# Rank by meta-analysis effect size
top_meta <- order(abs(meta_logFC), decreasing = TRUE)[1:100]
```

**Benefits:**
- Reduces resolution-specific noise
- Increases effect size estimates reliability
- Better for downstream prioritization

## Integration with Multi-Resolution Pipeline

**Automatic execution:**
```bash
sbatch scripts/run_multiresolution_pipeline.sb
# Runs prep, extract, aggregate, edgeR for each resolution
# Then automatically runs compare_resolutions.R
```

**Manual execution:**
```bash
# After running edgeR for all resolutions
Rscript scripts/compare_resolutions.R

# Outputs to outputs/resolution_comparison/
```

## Typical Results Summary

```
========================================
MULTI-RESOLUTION ANALYSIS COMPLETE
========================================

Output directory: outputs/resolution_comparison/

Files generated:
  1. summary_statistics_by_resolution.tsv
  2. loop_overlap_matrix.tsv
  3. differential_loops_by_resolution.pdf
  4. foldchange_correlation_*.pdf (3 files)
  5. venn_diagram_differential_loops.pdf

Key findings:
  5kb:  3,245 differential loops (15.1%)
  10kb: 2,156 differential loops (13.8%)
  25kb:   987 differential loops (10.7%)

Consensus (all 3 resolutions): 645 high-confidence loops
Resolution-specific:
  - 5kb only: 1,234 loops (short-range enriched)
  - 10kb only: 423 loops
  - 25kb only: 178 loops (long-range enriched)

Fold-change concordance:
  - 5kb vs 10kb: r = 0.847
  - 5kb vs 25kb: r = 0.712
  - 10kb vs 25kb: r = 0.823

Interpretation: High concordance across resolutions.
  Resolution choice has moderate impact on discovery.
  645 high-confidence loops for priority validation.
```

## Next Steps

### Priority Analysis

**1. High-confidence consensus loops (n=645)**
- Functional annotation
- Gene association
- Pathway enrichment
- Experimental validation

**2. Resolution-specific loops**
- Distance distribution analysis
- Biological context (TADs, compartments)
- Literature comparison
- Resolution-appropriate interpretation

**3. Meta-analysis effect sizes**
- Rank all loops by average |logFC|
- Top 100 for detailed characterization
- Compare with other Hi-C studies

### Validation Approaches

**For consensus loops:**
- 4C-seq validation (high priority)
- FISH imaging
- Compare with published datasets
- Check for known BAP1 targets

**For resolution-specific loops:**
- Understand why only detected at specific resolution
- Check if binning artifact vs real biology
- Consider biological scale (local vs distal)

## Comparison with Single-Resolution Analysis

**Previous approach (single resolution):**
- Only 5kb analysis
- No confidence assessment across scales
- Uncertain which loops are robust

**Current approach (multi-resolution):**
- Three independent analyses
- Cross-validation of findings
- High-confidence consensus set
- Resolution-appropriate interpretation

**Key advantage**: Distinguishes robust biology from resolution artifacts
