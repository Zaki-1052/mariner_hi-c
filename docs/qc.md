# Quality Control Visualization (qc.R)

## Biological Purpose

This script creates comprehensive visualizations to assess data quality and identify potential issues before proceeding with differential loop analysis. Quality control is critical in Hi-C analysis because technical artifacts can mimic biological signal.

**Key quality assessment areas:**
- **Count distribution patterns**: Are the data appropriately distributed for statistical modeling?
- **Sample comparability**: Do control and mutant samples show expected similarities?
- **Technical bias detection**: Are there systematic issues affecting interpretation?
- **Analysis readiness**: Is the data suitable for differential testing?

## Input Data

### Count Matrix
- **File**: `outputs/full/06_counts_matrix.rds`
- **Format**: N_loops × 2_samples matrix of aggregated Hi-C counts
- **Values**: Total contact frequencies per loop

## Visualization Components

### 1. Count Distribution Analysis

#### Density Plots
```r
# Convert to long format for ggplot
counts_long <- melt(counts)
colnames(counts_long) <- c("Loop", "Sample", "Count")

# Log-scale density plots
p1 <- ggplot(counts_long, aes(x = Count + 1, fill = Sample)) +
  geom_density(alpha = 0.5) +
  scale_x_log10() +
  theme_minimal() +
  labs(title = "Count Distribution by Sample",
       x = "Hi-C Contact Count + 1 (log10)",
       y = "Density")
```

**Biological interpretation:**
- **Similar distributions**: Indicates comparable chromatin architecture
- **Distribution shape**: Should be right-skewed (many weak, few strong loops)
- **Peak positions**: Similar peaks suggest consistent loop-calling sensitivity
- **Tail behavior**: Long tails may indicate technical artifacts or biological variation

**Expected patterns:**
- **Unimodal distributions**: Single peak for most loops
- **Log-normal shape**: Approximately normal on log scale
- **Sample overlap**: High overlap indicates biological consistency
- **No extreme outliers**: Distribution tails should be reasonable

#### Box Plot Comparison
```r
p2 <- ggplot(counts_long, aes(x = Sample, y = Count + 1, fill = Sample)) +
  geom_boxplot() +
  scale_y_log10() +
  theme_minimal() +
  labs(title = "Count Distribution Comparison",
       y = "Hi-C Contact Count + 1 (log10)")
```

**Quality indicators:**
- **Similar medians**: Indicates balanced data
- **Comparable IQRs**: Similar variability between samples
- **Outlier patterns**: Should be symmetric and reasonable
- **Whisker positions**: Should span expected biological range

### 2. MA Plot Analysis

#### MA Plot Construction
```r
# Calculate M (log fold change) and A (average expression)
A <- 0.5 * (log2(counts[,1] + 1) + log2(counts[,2] + 1))
M <- log2(counts[,2] + 1) - log2(counts[,1] + 1)

ma_data <- data.frame(A = A, M = M)

ggplot(ma_data, aes(x = A, y = M)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  geom_hline(yintercept = c(-1, 1), color = "blue", linetype = "dotted") +
  theme_minimal() +
  labs(title = "MA Plot: Mutant vs Control",
       x = "Average log2(Count + 1)",
       y = "log2 Fold Change (mut/ctrl)")
```

**MA Plot interpretation:**
- **X-axis (A)**: Average log2 expression (loop strength)
- **Y-axis (M)**: Log2 fold change between samples
- **Red line (y=0)**: No differential looping
- **Blue lines (y=±1)**: 2-fold change thresholds

**Quality assessment:**
- **Centered around y=0**: No systematic bias
- **Symmetric distribution**: Balanced up/down regulation
- **Intensity-dependent patterns**: May indicate normalization needs
- **Outlier identification**: Extreme changes for investigation

### 3. Sample Correlation Visualization

#### Scatter Plot Matrix
```r
# Scatter plot with correlation
correlation <- cor(counts[,1], counts[,2])

ggplot(data.frame(counts), aes(x = ctrl, y = mut)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", color = "red") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "blue") +
  scale_x_log10() + scale_y_log10() +
  labs(title = paste0("Sample Correlation (r = ", round(correlation, 3), ")"),
       x = "Control Counts (log10)",
       y = "Mutant Counts (log10)")
```

**Correlation interpretation:**
- **Points near diagonal**: High correlation between samples
- **Scatter pattern**: Reveals heteroscedasticity or systematic bias
- **Outliers**: Loops with unusual differential behavior
- **Linear fit**: Assesses proportional relationship

### 4. Statistical Distribution Testing

#### Diagnostic Plots for Normality
```r
# Q-Q plots for log-transformed data
par(mfrow = c(2, 2))

# Control sample Q-Q plot
log_ctrl <- log2(counts[,1] + 1)
qqnorm(log_ctrl, main = "Control Sample Q-Q Plot")
qqline(log_ctrl, col = "red")

# Mutant sample Q-Q plot
log_mut <- log2(counts[,2] + 1)
qqnorm(log_mut, main = "Mutant Sample Q-Q Plot")
qqline(log_mut, col = "red")

# Difference distribution
diff_counts <- log_mut - log_ctrl
qqnorm(diff_counts, main = "Log Fold Change Q-Q Plot")
qqline(diff_counts, col = "red")

# Mean-difference relationship
plot(A, M, main = "Mean vs Difference",
     xlab = "Average log2(Count + 1)",
     ylab = "log2 Fold Change")
abline(h = 0, col = "red", lty = 2)
```

**Normality assessment:**
- **Points on Q-Q line**: Indicates normal distribution
- **S-shaped curves**: Heavy tails (may need robust methods)
- **Systematic deviations**: Non-normal distributions

## Quality Control Metrics

### 1. Distribution Characteristics
```r
# Skewness and kurtosis
library(moments)

ctrl_skew <- skewness(log2(counts[,1] + 1))
mut_skew <- skewness(log2(counts[,2] + 1))
ctrl_kurt <- kurtosis(log2(counts[,1] + 1))
mut_kurt <- kurtosis(log2(counts[,2] + 1))

cat("Distribution characteristics:\n")
cat(sprintf("Control - Skewness: %.2f, Kurtosis: %.2f\n", ctrl_skew, ctrl_kurt))
cat(sprintf("Mutant - Skewness: %.2f, Kurtosis: %.2f\n", mut_skew, mut_kurt))
```

**Acceptable ranges:**
- **Skewness**: -1 to +1 (approximately symmetric)
- **Kurtosis**: 2-4 (normal to slightly heavy-tailed)

### 2. Outlier Detection
```r
# Identify outliers using IQR method
identify_outliers <- function(x) {
  Q1 <- quantile(x, 0.25)
  Q3 <- quantile(x, 0.75)
  IQR <- Q3 - Q1
  lower_bound <- Q1 - 1.5 * IQR
  upper_bound <- Q3 + 1.5 * IQR
  return(which(x < lower_bound | x > upper_bound))
}

ctrl_outliers <- identify_outliers(log2(counts[,1] + 1))
mut_outliers <- identify_outliers(log2(counts[,2] + 1))

cat(sprintf("Outliers - Control: %d (%.1f%%), Mutant: %d (%.1f%%)\n",
            length(ctrl_outliers), 100 * length(ctrl_outliers) / nrow(counts),
            length(mut_outliers), 100 * length(mut_outliers) / nrow(counts)))
```

### 3. Variance Homogeneity
```r
# Test for equal variances
var_test <- var.test(log2(counts[,1] + 1), log2(counts[,2] + 1))
cat(sprintf("Variance ratio (mut/ctrl): %.2f, p-value: %.3f\n",
            var_test$estimate, var_test$p.value))
```

**Variance interpretation:**
- **Ratio near 1.0**: Similar variability between samples
- **p-value > 0.05**: No significant difference in variances
- **Large ratios**: May indicate technical issues or biological differences

## Biological Quality Indicators

### 1. Dynamic Range Assessment
```r
# Calculate dynamic range
ctrl_range <- max(counts[,1]) / (min(counts[counts[,1] > 0, 1]))
mut_range <- max(counts[,2]) / (min(counts[counts[,2] > 0, 2]))

cat(sprintf("Dynamic range - Control: %.1f-fold, Mutant: %.1f-fold\n",
            ctrl_range, mut_range))
```

**Expected dynamic range:**
- **100-1000 fold**: Typical for Hi-C loop data
- **Similar between samples**: Indicates consistent sensitivity

### 2. Zero Inflation Assessment
```r
# Count zeros and low values
zeros_ctrl <- sum(counts[,1] == 0)
zeros_mut <- sum(counts[,2] == 0)
low_counts_ctrl <- sum(counts[,1] <= 5)
low_counts_mut <- sum(counts[,2] <= 5)

cat("Low count assessment:\n")
cat(sprintf("Zeros - Control: %d (%.1f%%), Mutant: %d (%.1f%%)\n",
            zeros_ctrl, 100 * zeros_ctrl / nrow(counts),
            zeros_mut, 100 * zeros_mut / nrow(counts)))
cat(sprintf("Low counts (≤5) - Control: %d (%.1f%%), Mutant: %d (%.1f%%)\n",
            low_counts_ctrl, 100 * low_counts_ctrl / nrow(counts),
            low_counts_mut, 100 * low_counts_mut / nrow(counts)))
```

**Acceptable levels:**
- **Zeros**: < 20% (higher suggests poor sensitivity)
- **Low counts**: < 40% (edgeR can handle some low counts)

## Technical Bias Detection

### 1. Systematic Bias Patterns
```r
# Check for intensity-dependent bias
smooth_fit <- loess(M ~ A, data = ma_data)
predicted_M <- predict(smooth_fit)
bias_range <- range(predicted_M)

cat(sprintf("Systematic bias range: %.3f to %.3f\n", bias_range[1], bias_range[2]))
if (max(abs(bias_range)) > 0.2) {
  cat("WARNING: Potential systematic bias detected\n")
}
```

### 2. Library Size Effects
```r
# Calculate effective library sizes
total_ctrl <- sum(counts[,1])
total_mut <- sum(counts[,2])
lib_ratio <- total_mut / total_ctrl

cat(sprintf("Total counts - Control: %d, Mutant: %d\n", total_ctrl, total_mut))
cat(sprintf("Library size ratio (mut/ctrl): %.2f\n", lib_ratio))

if (lib_ratio < 0.5 || lib_ratio > 2.0) {
  cat("WARNING: Large library size difference detected\n")
}
```

## Output Visualization Files

### Combined Plot Output
```r
# Combine plots using patchwork
library(patchwork)
combined_plots <- p1 + p2 + plot_layout(ncol = 2)
ggsave("qc_distributions.png", combined_plots, width = 12, height = 6)

# Save MA plot separately
ma_plot <- ggplot(ma_data, aes(x = A, y = M)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  geom_hline(yintercept = c(-1, 1), color = "blue", linetype = "dotted")
ggsave("ma_plot.png", ma_plot, width = 8, height = 6)
```

## Decision Making Framework

### Green Light Indicators (Proceed with Analysis)
- **High sample correlation (r > 0.7)**
- **Similar count distributions**
- **Low systematic bias in MA plot**
- **Reasonable outlier percentages (< 10%)**
- **Appropriate dynamic range (100-1000x)**

### Yellow Light Indicators (Proceed with Caution)
- **Moderate correlation (0.5 < r < 0.7)**
- **Some distribution differences**
- **Minor systematic bias**
- **Moderate outlier levels (10-20%)**

### Red Light Indicators (Address Issues First)
- **Low correlation (r < 0.5)**
- **Very different distributions**
- **Strong systematic bias**
- **High outlier levels (> 20%)**
- **Extreme library size differences**

## Troubleshooting Common Issues

### Issue 1: Low Sample Correlation
**Possible causes:**
- Genuine biological differences
- Technical batch effects
- Sequencing depth differences
- Different experimental conditions

**Solutions:**
- Check experimental metadata
- Consider batch correction
- Verify sample labeling
- Review library preparation protocols

### Issue 2: Systematic Bias in MA Plot
**Possible causes:**
- Library size differences
- GC content bias
- Fragment size distribution differences
- Normalization problems

**Solutions:**
- Apply TMM normalization
- Check Hi-C quality metrics
- Review processing pipeline
- Consider additional filtering

### Issue 3: Non-Normal Distributions
**Possible causes:**
- Outlier contamination
- Technical artifacts
- Inappropriate transformation
- Genuine biological heterogeneity

**Solutions:**
- Remove/investigate outliers
- Try alternative transformations
- Use robust statistical methods
- Increase sample size

## Biological Interpretation

The QC plots reveal important aspects of chromatin loop data:

1. **Data quality**: Technical adequacy for statistical analysis
2. **Biological signal**: Evidence of real vs. artifactual differences
3. **Analysis strategy**: Appropriate statistical approaches
4. **Interpretation confidence**: Reliability of downstream results

Good QC results increase confidence in biological conclusions, while poor QC results require investigation and potential re-analysis with modified parameters.

## Next Steps

Based on QC results:

**If QC passes**: Proceed with differential analysis using edgeR
**If QC shows issues**: Address problems before proceeding
**If results ambiguous**: Consider additional validation or alternative analysis methods

The QC visualization provides the foundation for interpreting all subsequent analysis results in the proper technical and biological context.