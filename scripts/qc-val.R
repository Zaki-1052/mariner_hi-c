# /expanse/lustre/projects/csd940/zalibhai/mariner/scripts/comprehensive_qc.R
# Comprehensive Quality Control and Validation of Mariner Pipeline Outputs
# Multi-resolution support: accepts resolution as command-line argument
# Purpose: Systematic validation of Hi-C loop extraction and aggregation
# Author: Generated for mariner pipeline validation
# Date: October 2025

# =============================================================================
# SETUP & CONFIGURATION
# =============================================================================

# Parse command-line arguments for resolution
args <- commandArgs(trailingOnly = TRUE)
RESOLUTION <- if (length(args) > 0) as.numeric(args[1]) else 5000

library(tidyverse)
library(mariner)
library(InteractionSet)
library(HDF5Array)
library(edgeR)
library(pheatmap)
library(viridis)
library(patchwork)
library(scales)

# Set working directory
setwd("/expanse/lustre/projects/csd940/zalibhai/mariner")

# Create resolution-specific QC output directory
input_dir <- sprintf("outputs/res_%dkb", RESOLUTION/1000)
qc_dir <- file.path(input_dir, "qc_report")
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

# Initialize QC report list
qc_report <- list()

cat("=============================================================================\n")
cat("MARINER PIPELINE QC VALIDATION\n")
cat(sprintf("RESOLUTION: %d bp (%d kb)\n", RESOLUTION, RESOLUTION/1000))
cat("=============================================================================\n\n")

# =============================================================================
# SECTION 1: DATA LOADING & BASIC INTEGRITY CHECKS
# =============================================================================

cat("SECTION 1: Data Loading & Integrity Validation\n")
cat("------------------------------------------------\n")

# Load all pipeline outputs
cat("Loading pipeline outputs...\n")

tryCatch({
  gi_list <- readRDS(file.path(input_dir, "01_ginteractions.rds"))
  merged <- readRDS(file.path(input_dir, "02_merged.rds"))
  binned <- readRDS(file.path(input_dir, "03_binned.rds"))
  buffered <- readRDS(file.path(input_dir, "04_buffered.rds"))
  metadata <- readRDS(file.path(input_dir, "05_metadata.rds"))
  counts_matrix <- readRDS(file.path(input_dir, "06_counts_matrix.rds"))
  edger_obj <- readRDS(file.path(input_dir, "06_edgeR_input.rds"))

  # Try to load HDF5-backed extracted data
  pixels <- loadHDF5SummarizedExperiment(
    dir = input_dir,
    prefix = "05_extracted"
  )

  # Try to load aggregation strategies comparison
  all_strategies <- tryCatch(
    readRDS(file.path(input_dir, "06_all_strategies.rds")),
    error = function(e) NULL
  )
  
  cat("✓ All data files loaded successfully\n\n")
  
}, error = function(e) {
  cat("✗ Error loading data files:\n")
  cat(paste("  ", e$message, "\n"))
  stop("Cannot proceed with QC - data loading failed")
})

# Basic dimensions check
cat("Basic Dimensions:\n")
cat(sprintf("  Initial loops (ctrl): %d\n", length(gi_list$ctrl)))
cat(sprintf("  Initial loops (mut):  %d\n", length(gi_list$mut)))
cat(sprintf("  Merged loops:         %d\n", length(merged)))
cat(sprintf("  Binned loops:         %d\n", length(binned)))
cat(sprintf("  Buffered regions:     %d\n", length(buffered)))
cat(sprintf("  Count matrix:         %d loops × %d samples\n", 
            nrow(counts_matrix), ncol(counts_matrix)))

qc_report$dimensions <- list(
  input_ctrl = length(gi_list$ctrl),
  input_mut = length(gi_list$mut),
  merged = length(merged),
  final_loops = nrow(counts_matrix),
  samples = ncol(counts_matrix)
)

# Check for data integrity issues
cat("\nData Integrity Checks:\n")

# NA values
na_count <- sum(is.na(counts_matrix))
na_pct <- 100 * na_count / length(counts_matrix)
cat(sprintf("  NA values: %d (%.2f%%)\n", na_count, na_pct))

# Zero values (expected for sparse Hi-C data)
zero_count <- sum(counts_matrix == 0)
zero_pct <- 100 * zero_count / length(counts_matrix)
cat(sprintf("  Zero values: %d (%.2f%%)\n", zero_count, zero_pct))

# Negative values (should never happen)
neg_count <- sum(counts_matrix < 0, na.rm = TRUE)
if (neg_count > 0) {
  cat(sprintf("  ✗ WARNING: %d negative values detected!\n", neg_count))
  qc_report$warnings <- c(qc_report$warnings, 
                          paste("Negative counts detected:", neg_count))
} else {
  cat("  ✓ No negative values detected\n")
}

qc_report$integrity <- list(
  na_count = na_count,
  na_percentage = na_pct,
  zero_count = zero_count,
  zero_percentage = zero_pct,
  negative_values = neg_count
)

cat("\n")

# =============================================================================
# SECTION 2: SAMPLE-LEVEL QUALITY CONTROL
# =============================================================================

cat("SECTION 2: Sample-Level Quality Assessment\n")
cat("-------------------------------------------\n")

# Basic count statistics per sample
cat("Count Statistics by Sample:\n")
sample_stats <- data.frame(
  Sample = colnames(counts_matrix),
  Min = apply(counts_matrix, 2, min, na.rm = TRUE),
  Q1 = apply(counts_matrix, 2, quantile, probs = 0.25, na.rm = TRUE),
  Median = apply(counts_matrix, 2, median, na.rm = TRUE),
  Mean = apply(counts_matrix, 2, mean, na.rm = TRUE),
  Q3 = apply(counts_matrix, 2, quantile, probs = 0.75, na.rm = TRUE),
  Max = apply(counts_matrix, 2, max, na.rm = TRUE),
  SD = apply(counts_matrix, 2, sd, na.rm = TRUE)
)
print(sample_stats)
cat("\n")

qc_report$sample_stats <- sample_stats

# Library size (total counts per sample)
lib_sizes <- colSums(counts_matrix, na.rm = TRUE)
cat("Library Sizes (Total Counts):\n")
print(lib_sizes)

# For replicates, compare within-group and between-group
n_samples <- ncol(counts_matrix)
if (n_samples == 6) {
  # Replicate-aware analysis
  ctrl_libs <- lib_sizes[1:3]
  mut_libs <- lib_sizes[4:6]
  cat(sprintf("\n  Control group:\n"))
  cat(sprintf("    Mean: %.0f, SD: %.0f, CV: %.2f\n",
              mean(ctrl_libs), sd(ctrl_libs), sd(ctrl_libs)/mean(ctrl_libs)))
  cat(sprintf("  Mutant group:\n"))
  cat(sprintf("    Mean: %.0f, SD: %.0f, CV: %.2f\n",
              mean(mut_libs), sd(mut_libs), sd(mut_libs)/mean(mut_libs)))
  cat(sprintf("  Ratio (mut/ctrl medians): %.3f\n", median(mut_libs)/median(ctrl_libs)))
} else {
  # Original 2-sample analysis
  cat(sprintf("  Size ratio (mut/ctrl): %.3f\n", lib_sizes[2]/lib_sizes[1]))
}
cat("\n")

# Check if library sizes are comparable (should be within 2-fold)
lib_ratio <- max(lib_sizes) / min(lib_sizes)
if (lib_ratio > 2.0) {
  cat("  ✗ WARNING: Large library size difference (>2-fold)\n")
  qc_report$warnings <- c(qc_report$warnings,
                          sprintf("Library size ratio: %.2f", lib_ratio))
} else {
  cat("  ✓ Library sizes are comparable\n")
}

qc_report$library_sizes <- lib_sizes

# Sample correlation
cat("Sample Correlation Analysis:\n")
cor_matrix <- cor(counts_matrix, use = "complete.obs", method = "pearson")

if (n_samples == 6) {
  # Replicate-aware correlation analysis
  cat("\nFull correlation matrix:\n")
  print(round(cor_matrix, 3))

  # Within-group correlations
  cat("\nWithin-group correlations:\n")
  ctrl_cors <- cor_matrix[1:3, 1:3][upper.tri(cor_matrix[1:3, 1:3])]
  mut_cors <- cor_matrix[4:6, 4:6][upper.tri(cor_matrix[4:6, 4:6])]
  cat(sprintf("  Control replicates: %.3f ± %.3f (range: %.3f-%.3f)\n",
              mean(ctrl_cors), sd(ctrl_cors), min(ctrl_cors), max(ctrl_cors)))
  cat(sprintf("  Mutant replicates:  %.3f ± %.3f (range: %.3f-%.3f)\n",
              mean(mut_cors), sd(mut_cors), min(mut_cors), max(mut_cors)))

  # Between-group correlations
  between_cors <- as.vector(cor_matrix[1:3, 4:6])
  cat(sprintf("  Between groups:     %.3f ± %.3f (range: %.3f-%.3f)\n",
              mean(between_cors), sd(between_cors), min(between_cors), max(between_cors)))

  # Quality assessment
  cat("\nQuality Assessment:\n")
  if (mean(ctrl_cors) > 0.95 && mean(mut_cors) > 0.95) {
    cat("  ✓ Excellent within-group reproducibility (>0.95)\n")
  } else if (mean(ctrl_cors) > 0.90 && mean(mut_cors) > 0.90) {
    cat("  ✓ Good within-group reproducibility (>0.90)\n")
  } else {
    cat("  ⚠ Lower than expected within-group correlation\n")
  }

  if ((mean(ctrl_cors) - mean(between_cors)) > 0.05 ||
      (mean(mut_cors) - mean(between_cors)) > 0.05) {
    cat("  ✓ Clear biological signal (within > between)\n")
  } else {
    cat("  ⚠ Weak biological differences\n")
  }

  qc_report$correlation <- list(
    within_ctrl = mean(ctrl_cors),
    within_mut = mean(mut_cors),
    between = mean(between_cors)
  )
  cor_value <- mean(between_cors)  # For downstream compatibility

} else {
  # Original 2-sample analysis
  cor_value <- cor_matrix[1, 2]
  cat(sprintf("  Pearson correlation: %.4f\n", cor_value))

  # Interpret correlation
  if (cor_value > 0.95) {
    cat("  ⚠ Very high correlation - may indicate limited biological difference\n")
  } else if (cor_value > 0.85) {
    cat("  ✓ High correlation - good technical quality\n")
  } else if (cor_value > 0.70) {
    cat("  ⚠ Moderate correlation - check for batch effects\n")
  } else {
    cat("  ✗ Low correlation - possible quality issues\n")
  }

  qc_report$correlation <- cor_value
}

# Spearman correlation (rank-based, robust to outliers)
spearman_cor <- cor(counts_matrix, use = "complete.obs", method = "spearman")
if (n_samples == 2) {
  cat(sprintf("  Spearman correlation: %.4f\n", spearman_cor[1,2]))
}
cat("\n")

# Visualization: Correlation heatmap
pdf(file.path(qc_dir, "01_sample_correlation.pdf"), width = 6, height = 5)
pheatmap(cor_matrix,
         display_numbers = TRUE,
         number_format = "%.3f",
         color = viridis(100),
         main = "Sample Correlation\n(Pearson on raw counts)",
         fontsize_number = 14,
         cellwidth = 80,
         cellheight = 80)
dev.off()

# Visualization: Scatter plot with marginal distributions
pdf(file.path(qc_dir, "02_sample_scatter.pdf"), width = 8, height = 8)
par(mar = c(5, 5, 4, 2))

# Log transform for visualization (add 1 to handle zeros)
log_counts <- log2(counts_matrix + 1)

# Main scatter plot
plot(log_counts[,1], log_counts[,2],
     xlab = "Control [log2(count + 1)]",
     ylab = "Mutant [log2(count + 1)]",
     main = sprintf("Sample Correlation\nr = %.3f", cor_value),
     pch = 16,
     col = rgb(0, 0, 0, 0.3),
     cex = 0.5)

# Add diagonal line (perfect correlation)
abline(a = 0, b = 1, col = "red", lwd = 2, lty = 2)

# Add regression line
fit <- lm(log_counts[,2] ~ log_counts[,1])
abline(fit, col = "blue", lwd = 2)

# Add legend
legend("topleft", 
       legend = c("y = x (perfect correlation)", 
                  sprintf("Linear fit (slope = %.3f)", coef(fit)[2])),
       col = c("red", "blue"),
       lty = c(2, 1),
       lwd = 2)

dev.off()

# Visualization: Distribution comparison
pdf(file.path(qc_dir, "03_count_distributions.pdf"), width = 12, height = 8)

counts_long <- counts_matrix %>%
  as.data.frame() %>%
  mutate(loop_id = 1:n()) %>%
  pivot_longer(cols = -loop_id, names_to = "Sample", values_to = "Count")

# Panel 1: Density plots (log scale)
p1 <- ggplot(counts_long, aes(x = Count + 1, fill = Sample)) +
  geom_density(alpha = 0.6) +
  scale_x_log10(labels = comma) +
  scale_fill_manual(values = c("ctrl" = "#377EB8", "mut" = "#E41A1C")) +
  theme_minimal(base_size = 12) +
  labs(title = "Count Distribution (Density)",
       x = "Hi-C Contact Count + 1 (log10)",
       y = "Density") +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5, face = "bold"))

# Panel 2: Box plots (log scale)
p2 <- ggplot(counts_long, aes(x = Sample, y = Count + 1, fill = Sample)) +
  geom_boxplot(outlier.alpha = 0.3, outlier.size = 0.5) +
  scale_y_log10(labels = comma) +
  scale_fill_manual(values = c("ctrl" = "#377EB8", "mut" = "#E41A1C")) +
  theme_minimal(base_size = 12) +
  labs(title = "Count Distribution (Boxplot)",
       x = "",
       y = "Hi-C Contact Count + 1 (log10)") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"))

# Panel 3: Violin plots
p3 <- ggplot(counts_long, aes(x = Sample, y = Count + 1, fill = Sample)) +
  geom_violin(alpha = 0.6) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.8) +
  scale_y_log10(labels = comma) +
  scale_fill_manual(values = c("ctrl" = "#377EB8", "mut" = "#E41A1C")) +
  theme_minimal(base_size = 12) +
  labs(title = "Count Distribution (Violin)",
       x = "",
       y = "Hi-C Contact Count + 1 (log10)") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"))

# Panel 4: Cumulative distribution
p4 <- ggplot(counts_long, aes(x = Count + 1, color = Sample)) +
  stat_ecdf(size = 1) +
  scale_x_log10(labels = comma) +
  scale_color_manual(values = c("ctrl" = "#377EB8", "mut" = "#E41A1C")) +
  theme_minimal(base_size = 12) +
  labs(title = "Cumulative Distribution",
       x = "Hi-C Contact Count + 1 (log10)",
       y = "Cumulative Probability") +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5, face = "bold"))

# Combine plots
(p1 | p2) / (p3 | p4)

dev.off()

cat("Sample-level QC plots saved to:", qc_dir, "\n\n")

# =============================================================================
# SECTION 3: LOOP-LEVEL QUALITY ASSESSMENT
# =============================================================================

cat("SECTION 3: Loop-Level Quality Assessment\n")
cat("-----------------------------------------\n")

# Calculate MA values (M = log fold change, A = average expression)
# For replicates, use group averages
if (n_samples == 6) {
  ctrl_avg <- rowMeans(counts_matrix[, 1:3])
  mut_avg <- rowMeans(counts_matrix[, 4:6])
  A <- 0.5 * (log2(ctrl_avg + 1) + log2(mut_avg + 1))
  M <- log2(mut_avg + 1) - log2(ctrl_avg + 1)
} else {
  A <- 0.5 * (log2(counts_matrix[,1] + 1) + log2(counts_matrix[,2] + 1))
  M <- log2(counts_matrix[,2] + 1) - log2(counts_matrix[,1] + 1)
}

cat("MA Statistics:\n")
cat(sprintf("  Average expression (A): %.2f ± %.2f\n", mean(A, na.rm = TRUE), sd(A, na.rm = TRUE)))
cat(sprintf("  Log2 fold change (M):   %.3f ± %.3f\n", mean(M, na.rm = TRUE), sd(M, na.rm = TRUE)))
cat(sprintf("  |M| > 1 (2-fold change): %d loops (%.1f%%)\n", 
            sum(abs(M) > 1, na.rm = TRUE), 
            100 * sum(abs(M) > 1, na.rm = TRUE) / length(M)))
cat(sprintf("  |M| > 2 (4-fold change): %d loops (%.1f%%)\n", 
            sum(abs(M) > 2, na.rm = TRUE), 
            100 * sum(abs(M) > 2, na.rm = TRUE) / length(M)))
cat("\n")

qc_report$ma_stats <- list(
  mean_A = mean(A, na.rm = TRUE),
  sd_A = sd(A, na.rm = TRUE),
  mean_M = mean(M, na.rm = TRUE),
  sd_M = sd(M, na.rm = TRUE),
  twofold = sum(abs(M) > 1, na.rm = TRUE),
  fourfold = sum(abs(M) > 2, na.rm = TRUE)
)

# Visualization: MA plot
pdf(file.path(qc_dir, "04_ma_plot.pdf"), width = 10, height = 8)

ma_data <- data.frame(A = A, M = M)

# Color by magnitude of change
ma_data$color <- case_when(
  abs(M) > 2 ~ "Large change (|M| > 2)",
  abs(M) > 1 ~ "Moderate change (|M| > 1)",
  TRUE ~ "Small change (|M| ≤ 1)"
)

ggplot(ma_data, aes(x = A, y = M, color = color)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_hline(yintercept = 0, color = "black", linetype = "solid", size = 0.5) +
  geom_hline(yintercept = c(-1, 1), color = "blue", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = c(-2, 2), color = "red", linetype = "dotted", size = 0.5) +
  scale_color_manual(values = c("Small change (|M| ≤ 1)" = "gray60",
                                "Moderate change (|M| > 1)" = "orange",
                                "Large change (|M| > 2)" = "red")) +
  theme_minimal(base_size = 12) +
  labs(title = "MA Plot: Mutant vs Control",
       subtitle = sprintf("Correlation = %.3f | Mean M = %.3f", cor_value, mean(M, na.rm = TRUE)),
       x = "Average Expression [A = 0.5 × (log2(ctrl+1) + log2(mut+1))]",
       y = "Log2 Fold Change [M = log2(mut+1) - log2(ctrl+1)]",
       color = "Change Magnitude") +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))

dev.off()

# Variance analysis - coefficient of variation
cat("Loop Variability Analysis:\n")

# Calculate CV for each loop across samples
loop_means <- rowMeans(counts_matrix, na.rm = TRUE)
loop_sds <- apply(counts_matrix, 1, sd, na.rm = TRUE)
loop_cv <- loop_sds / loop_means

# Remove infinite/NaN values (from loops with zero mean)
loop_cv <- loop_cv[is.finite(loop_cv)]

cat(sprintf("  Coefficient of Variation: %.3f ± %.3f\n", 
            mean(loop_cv, na.rm = TRUE), sd(loop_cv, na.rm = TRUE)))
cat(sprintf("  High variability loops (CV > 1): %d (%.1f%%)\n",
            sum(loop_cv > 1, na.rm = TRUE),
            100 * sum(loop_cv > 1, na.rm = TRUE) / length(loop_cv)))

qc_report$variability <- list(
  mean_cv = mean(loop_cv, na.rm = TRUE),
  high_var_count = sum(loop_cv > 1, na.rm = TRUE)
)

# Identify top variable loops for visualization
top_n <- 50
top_variable_idx <- order(loop_cv, decreasing = TRUE, na.last = TRUE)[1:top_n]

# Visualization: Heatmap of top variable loops
pdf(file.path(qc_dir, "05_variable_loops_heatmap.pdf"), width = 8, height = 10)

counts_scaled <- t(scale(t(counts_matrix[top_variable_idx, ])))

pheatmap(counts_scaled,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = sprintf("Top %d Most Variable Loops\n(Z-score scaled)", top_n),
         labels_col = c("Control", "Mutant"),
         show_rownames = FALSE,
         fontsize = 10)

dev.off()

cat("Loop-level QC plots saved\n\n")

# =============================================================================
# SECTION 4: SPATIAL & DISTANCE ANALYSIS
# =============================================================================

cat("SECTION 4: Spatial & Distance Analysis\n")
cat("---------------------------------------\n")

# Extract genomic coordinates from merged loops
if (n_samples == 6) {
  loop_coords <- data.frame(
    chr = as.character(seqnames(anchors(merged, type = "first"))),
    start1 = start(anchors(merged, type = "first")),
    end1 = end(anchors(merged, type = "first")),
    start2 = start(anchors(merged, type = "second")),
    end2 = end(anchors(merged, type = "second")),
    ctrl_count = rowMeans(counts_matrix[, 1:3]),
    mut_count = rowMeans(counts_matrix[, 4:6]),
    total_count = rowSums(counts_matrix),
    log2FC = M
  )
} else {
  loop_coords <- data.frame(
    chr = as.character(seqnames(anchors(merged, type = "first"))),
    start1 = start(anchors(merged, type = "first")),
    end1 = end(anchors(merged, type = "first")),
    start2 = start(anchors(merged, type = "second")),
    end2 = end(anchors(merged, type = "second")),
    ctrl_count = counts_matrix[, 1],
    mut_count = counts_matrix[, 2],
    total_count = rowSums(counts_matrix),
    log2FC = M
  )
}

# Calculate genomic distance (for intrachromosomal loops)
loop_coords$distance <- ifelse(
  loop_coords$start2 > loop_coords$start1,
  loop_coords$start2 - loop_coords$start1,
  NA
)

# Chromosome distribution
chr_summary <- loop_coords %>%
  group_by(chr) %>%
  summarise(
    n_loops = n(),
    mean_ctrl = mean(ctrl_count, na.rm = TRUE),
    mean_mut = mean(mut_count, na.rm = TRUE),
    mean_total = mean(total_count, na.rm = TRUE),
    mean_fc = mean(log2FC, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_loops))

cat("Loop Distribution by Chromosome:\n")
print(head(chr_summary, 10))
cat("\n")

qc_report$chr_distribution <- chr_summary

# Distance statistics
cat("Genomic Distance Statistics (intrachromosomal loops):\n")
valid_distances <- loop_coords$distance[!is.na(loop_coords$distance)]
cat(sprintf("  Loops with distance: %d / %d (%.1f%%)\n",
            length(valid_distances), nrow(loop_coords),
            100 * length(valid_distances) / nrow(loop_coords)))
cat(sprintf("  Min distance:  %s\n", comma(min(valid_distances))))
cat(sprintf("  Median distance: %s\n", comma(median(valid_distances))))
cat(sprintf("  Mean distance: %s\n", comma(mean(valid_distances))))
cat(sprintf("  Max distance:  %s\n", comma(max(valid_distances))))
cat("\n")

qc_report$distance_stats <- list(
  intrachromosomal_pct = 100 * length(valid_distances) / nrow(loop_coords),
  median_distance = median(valid_distances),
  mean_distance = mean(valid_distances)
)

# Visualization: Chromosome distribution
pdf(file.path(qc_dir, "06_chromosome_distribution.pdf"), width = 12, height = 6)

# Reorder chromosomes naturally
chr_order <- gtools::mixedsort(unique(chr_summary$chr))
chr_summary$chr <- factor(chr_summary$chr, levels = chr_order)

ggplot(chr_summary, aes(x = chr, y = n_loops)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
  geom_text(aes(label = n_loops), vjust = -0.5, size = 3) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(title = "Loop Distribution by Chromosome",
       x = "Chromosome",
       y = "Number of Loops")

dev.off()

# Visualization: Distance decay
pdf(file.path(qc_dir, "07_distance_decay.pdf"), width = 12, height = 8)

distance_data <- loop_coords %>%
  filter(!is.na(distance)) %>%
  mutate(distance_mb = distance / 1e6)

# Panel 1: Total signal vs distance
p1 <- ggplot(distance_data, aes(x = distance_mb, y = total_count)) +
  geom_point(alpha = 0.3, size = 1) +
  geom_smooth(method = "loess", color = "red", se = TRUE) +
  scale_x_log10(labels = comma) +
  scale_y_log10(labels = comma) +
  theme_minimal(base_size = 12) +
  labs(title = "Distance Decay: Total Signal",
       x = "Genomic Distance (Mb, log10)",
       y = "Total Count (ctrl + mut, log10)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Panel 2: Control signal vs distance
p2 <- ggplot(distance_data, aes(x = distance_mb, y = ctrl_count)) +
  geom_point(alpha = 0.3, size = 1, color = "#377EB8") +
  geom_smooth(method = "loess", color = "darkblue", se = TRUE) +
  scale_x_log10(labels = comma) +
  scale_y_log10(labels = comma) +
  theme_minimal(base_size = 12) +
  labs(title = "Distance Decay: Control",
       x = "Genomic Distance (Mb, log10)",
       y = "Control Count (log10)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Panel 3: Mutant signal vs distance
p3 <- ggplot(distance_data, aes(x = distance_mb, y = mut_count)) +
  geom_point(alpha = 0.3, size = 1, color = "#E41A1C") +
  geom_smooth(method = "loess", color = "darkred", se = TRUE) +
  scale_x_log10(labels = comma) +
  scale_y_log10(labels = comma) +
  theme_minimal(base_size = 12) +
  labs(title = "Distance Decay: Mutant",
       x = "Genomic Distance (Mb, log10)",
       y = "Mutant Count (log10)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Panel 4: Log2FC vs distance
p4 <- ggplot(distance_data, aes(x = distance_mb, y = log2FC)) +
  geom_point(alpha = 0.3, size = 1) +
  geom_smooth(method = "loess", color = "purple", se = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_log10(labels = comma) +
  theme_minimal(base_size = 12) +
  labs(title = "Fold Change vs Distance",
       x = "Genomic Distance (Mb, log10)",
       y = "Log2 Fold Change (mut/ctrl)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

(p1 | p2) / (p3 | p4)

dev.off()

cat("Spatial analysis plots saved\n\n")

# =============================================================================
# SECTION 5: MATRIX-LEVEL INSPECTION (5x5 BUFFERS)
# =============================================================================

cat("SECTION 5: Matrix-Level Buffer Validation\n")
cat("------------------------------------------\n")

cat("Inspecting 5x5 extraction matrices...\n")

# Get the count array from HDF5
count_array <- counts(pixels)
array_dims <- dim(count_array)
cat(sprintf("  Array dimensions: %d x %d x %d x %d\n", 
            array_dims[1], array_dims[2], array_dims[3], array_dims[4]))
cat("  (rows x cols x loops x samples)\n\n")

# Function to analyze a single loop's matrices
analyze_loop_matrix <- function(loop_id, count_array, n_samples) {

  if (n_samples == 6) {
    # For replicates, average across samples within each group
    ctrl_mat <- apply(count_array[, , loop_id, 1:3], c(1,2), mean, na.rm = TRUE)
    mut_mat <- apply(count_array[, , loop_id, 4:6], c(1,2), mean, na.rm = TRUE)
  } else {
    # Original 2-sample analysis
    ctrl_mat <- as.matrix(count_array[, , loop_id, 1])
    mut_mat <- as.matrix(count_array[, , loop_id, 2])
  }

  # Calculate statistics
  ctrl_sum <- sum(ctrl_mat, na.rm = TRUE)
  mut_sum <- sum(mut_mat, na.rm = TRUE)
  ctrl_center <- ctrl_mat[3, 3]
  mut_center <- mut_mat[3, 3]

  # Detect if peak is in center (expected) or shifted
  ctrl_peak_pos <- which(ctrl_mat == max(ctrl_mat, na.rm = TRUE), arr.ind = TRUE)[1, ]
  mut_peak_pos <- which(mut_mat == max(mut_mat, na.rm = TRUE), arr.ind = TRUE)[1, ]

  ctrl_centered <- all(ctrl_peak_pos == c(3, 3))
  mut_centered <- all(mut_peak_pos == c(3, 3))
  shift_detected <- !identical(ctrl_peak_pos, mut_peak_pos)

  # Calculate center enrichment (center pixel / total)
  ctrl_center_enrichment <- ctrl_center / ctrl_sum
  mut_center_enrichment <- mut_center / mut_sum

  return(list(
    ctrl_sum = ctrl_sum,
    mut_sum = mut_sum,
    ctrl_center = ctrl_center,
    mut_center = mut_center,
    ctrl_centered = ctrl_centered,
    mut_centered = mut_centered,
    shift_detected = shift_detected,
    ctrl_center_enrichment = ctrl_center_enrichment,
    mut_center_enrichment = mut_center_enrichment,
    ctrl_peak_pos = ctrl_peak_pos,
    mut_peak_pos = mut_peak_pos
  ))
}

# Analyze all loops
cat("Analyzing all loop matrices (this may take a moment)...\n")
n_loops <- array_dims[3]
matrix_stats <- lapply(1:n_loops, function(i) {
  if (i %% 5000 == 0) cat(sprintf("  Processed %d / %d loops\n", i, n_loops))
  analyze_loop_matrix(i, count_array, n_samples)
})

# Summarize matrix-level QC
cat("\nMatrix-Level Summary:\n")

centered_ctrl <- sum(sapply(matrix_stats, function(x) x$ctrl_centered))
centered_mut <- sum(sapply(matrix_stats, function(x) x$mut_centered))
shifts <- sum(sapply(matrix_stats, function(x) x$shift_detected))

cat(sprintf("  Control peaks centered: %d / %d (%.1f%%)\n",
            centered_ctrl, n_loops, 100 * centered_ctrl / n_loops))
cat(sprintf("  Mutant peaks centered:  %d / %d (%.1f%%)\n",
            centered_mut, n_loops, 100 * centered_mut / n_loops))
cat(sprintf("  Positional shifts detected: %d / %d (%.1f%%)\n",
            shifts, n_loops, 100 * shifts / n_loops))

# Center enrichment
center_enrichments <- data.frame(
  ctrl = sapply(matrix_stats, function(x) x$ctrl_center_enrichment),
  mut = sapply(matrix_stats, function(x) x$mut_center_enrichment)
)

cat(sprintf("  Mean center enrichment (ctrl): %.3f ± %.3f\n",
            mean(center_enrichments$ctrl, na.rm = TRUE),
            sd(center_enrichments$ctrl, na.rm = TRUE)))
cat(sprintf("  Mean center enrichment (mut):  %.3f ± %.3f\n",
            mean(center_enrichments$mut, na.rm = TRUE),
            sd(center_enrichments$mut, na.rm = TRUE)))
cat("\n")

qc_report$matrix_stats <- list(
  centered_ctrl_pct = 100 * centered_ctrl / n_loops,
  centered_mut_pct = 100 * centered_mut / n_loops,
  shifts_pct = 100 * shifts / n_loops,
  mean_center_enrichment_ctrl = mean(center_enrichments$ctrl, na.rm = TRUE),
  mean_center_enrichment_mut = mean(center_enrichments$mut, na.rm = TRUE)
)

# =============================================================================
# SECTION 5B: EXTRACT AND SAVE PER-LOOP SHIFT STATUS
# =============================================================================

cat("Extracting per-loop shift status for downstream analysis...\n")

# Extract shift status as boolean vector (TRUE = shift detected)
loop_shift_status <- sapply(matrix_stats, function(x) x$shift_detected)

# Validate dimensions
stopifnot(
	    "Shift status vector length mismatch" = length(loop_shift_status) == n_loops
	    )

# Validate it's actually boolean
stopifnot(
	    "Shift status must be logical type" = is.logical(loop_shift_status)
	    )

cat(sprintf("  ✓ Extracted shift status for %d loops\n", length(loop_shift_status)))
cat(sprintf("  - Shifts detected: %d (%.1f%%)\n", 
	                sum(loop_shift_status), 
			            100 * mean(loop_shift_status)))
cat(sprintf("  - Non-shifted: %d (%.1f%%)\n\n", 
	                sum(!loop_shift_status), 
			            100 * mean(!loop_shift_status)))

# Add to QC report with structured hierarchy
if (!"shift_analysis" %in% names(qc_report)) {
	  qc_report$shift_analysis <- list()
}

qc_report$shift_analysis$loop_shift_status <- loop_shift_status
qc_report$shift_analysis$n_shifted <- sum(loop_shift_status)
qc_report$shift_analysis$n_nonshifted <- sum(!loop_shift_status)
qc_report$shift_analysis$percent_shifted <- 100 * mean(loop_shift_status)

# Also save as standalone file for easy access
shift_status_file <- file.path(qc_dir, "loop_shift_status.rds")
saveRDS(loop_shift_status, shift_status_file)
cat(sprintf("  ✓ Shift status saved to: %s\n", shift_status_file))

# Create human-readable summary table
shift_summary_df <- data.frame(
			         loop_id = paste0("loop_", 1:n_loops),
				   shift_detected = loop_shift_status,
				   ctrl_peak_row = sapply(matrix_stats, function(x) x$ctrl_peak_pos[1]),
				     ctrl_peak_col = sapply(matrix_stats, function(x) x$ctrl_peak_pos[2]),
				     mut_peak_row = sapply(matrix_stats, function(x) x$mut_peak_pos[1]),
				       mut_peak_col = sapply(matrix_stats, function(x) x$mut_peak_pos[2]),
				       ctrl_centered = sapply(matrix_stats, function(x) x$ctrl_centered),
				         mut_centered = sapply(matrix_stats, function(x) x$mut_centered)
				       )

shift_summary_file <- file.path(qc_dir, "loop_shift_summary.tsv")
write.table(shift_summary_df, 
	                shift_summary_file,
			            sep = "\t", 
			            quote = FALSE, 
				                row.names = FALSE)
cat(sprintf("  ✓ Shift summary table saved to: %s\n\n", shift_summary_file))

cat("Per-loop shift status extraction complete\n\n")

# Visualize example matrices (top 6 by total signal)
pdf(file.path(qc_dir, "08_example_matrices.pdf"), width = 14, height = 10)

top_signal_idx <- order(rowSums(counts_matrix), decreasing = TRUE)[1:6]

par(mfrow = c(3, 4), mar = c(2, 2, 3, 2))

for (i in 1:6) {
  loop_id <- top_signal_idx[i]

  # Extract matrices (average across replicates if n=6)
  if (n_samples == 6) {
    ctrl_mat <- apply(count_array[, , loop_id, 1:3], c(1,2), mean, na.rm = TRUE)
    mut_mat <- apply(count_array[, , loop_id, 4:6], c(1,2), mean, na.rm = TRUE)
  } else {
    ctrl_mat <- as.matrix(count_array[, , loop_id, 1])
    mut_mat <- as.matrix(count_array[, , loop_id, 2])
  }
  
  # Plot control matrix
  image(1:5, 1:5, ctrl_mat,
        col = viridis(100),
        main = sprintf("Loop %d - Control\nSum=%.0f", loop_id, sum(ctrl_mat)),
        xlab = "", ylab = "", axes = FALSE)
  axis(1, at = 1:5, labels = 1:5)
  axis(2, at = 1:5, labels = 1:5)
  
  # Add text values
  for (row in 1:5) {
    for (col in 1:5) {
      text(col, row, round(ctrl_mat[row, col], 1), col = "white", cex = 0.8)
    }
  }
  
  # Plot mutant matrix
  image(1:5, 1:5, mut_mat,
        col = viridis(100),
        main = sprintf("Loop %d - Mutant\nSum=%.0f", loop_id, sum(mut_mat)),
        xlab = "", ylab = "", axes = FALSE)
  axis(1, at = 1:5, labels = 1:5)
  axis(2, at = 1:5, labels = 1:5)
  
  # Add text values
  for (row in 1:5) {
    for (col in 1:5) {
      text(col, row, round(mut_mat[row, col], 1), col = "white", cex = 0.8)
    }
  }
}

dev.off()

cat("Matrix inspection plots saved\n\n")

# =============================================================================
# SECTION 6: AGGREGATION STRATEGY COMPARISON
# =============================================================================

cat("SECTION 6: Aggregation Strategy Comparison\n")
cat("-------------------------------------------\n")

if (!is.null(all_strategies) && 
    !is.null(all_strategies$counts) &&
    length(all_strategies$counts) > 1) {
  
  cat("Comparing different aggregation methods...\n\n")
  
  # Extract the matrices
  sum_counts <- all_strategies$counts$sum
  weighted_counts <- all_strategies$counts$weighted
  
  # Get strategy names
  strategy_names <- names(all_strategies$counts)
  
  cat("Available aggregation strategies:\n")
  cat(paste("  -", strategy_names), sep = "\n")
  cat("\n")
  
  # Create data frames for correlation analysis
  # Each strategy as separate columns
  ctrl_strategies <- data.frame(
    sum = sum_counts[, "ctrl"],
    weighted = weighted_counts[, "ctrl"]
  )
  
  mut_strategies <- data.frame(
    sum = sum_counts[, "mut"],
    weighted = weighted_counts[, "mut"]
  )
  
  # Compare correlations between strategies
  cat("Correlation between aggregation strategies (Control):\n")
  cor_strategies_ctrl <- cor(ctrl_strategies, use = "complete.obs")
  print(round(cor_strategies_ctrl, 3))
  cat("\n")
  
  cat("Correlation between aggregation strategies (Mutant):\n")
  cor_strategies_mut <- cor(mut_strategies, use = "complete.obs")
  print(round(cor_strategies_mut, 3))
  cat("\n")
  
  # Compare actual values
  cat("Strategy comparison statistics:\n")
  cat(sprintf("  Control sum range:      %.1f - %.1f\n", 
              min(ctrl_strategies$sum, na.rm=TRUE),
              max(ctrl_strategies$sum, na.rm=TRUE)))
  cat(sprintf("  Control weighted range: %.1f - %.1f\n",
              min(ctrl_strategies$weighted, na.rm=TRUE),
              max(ctrl_strategies$weighted, na.rm=TRUE)))
  cat(sprintf("  Mutant sum range:       %.1f - %.1f\n",
              min(mut_strategies$sum, na.rm=TRUE),
              max(mut_strategies$sum, na.rm=TRUE)))
  cat(sprintf("  Mutant weighted range:  %.1f - %.1f\n",
              min(mut_strategies$weighted, na.rm=TRUE),
              max(mut_strategies$weighted, na.rm=TRUE)))
  cat("\n")
  
  # Check which matches the main count matrix
  sum_matches <- all.equal(sum_counts, counts_matrix, check.attributes = FALSE)
  weighted_matches <- all.equal(weighted_counts, counts_matrix, check.attributes = FALSE)
  
  if (isTRUE(sum_matches)) {
    used_strategy <- "sum"
    cat("✓ Main count matrix uses 'sum' aggregation\n")
  } else if (isTRUE(weighted_matches)) {
    used_strategy <- "weighted"
    cat("✓ Main count matrix uses 'weighted' aggregation\n")
  } else {
    used_strategy <- "unknown"
    cat("⚠ Could not determine which strategy was used for main analysis\n")
  }
  cat("\n")
  
  # Visualization: Strategy comparison
  pdf(file.path(qc_dir, "09_strategy_comparison.pdf"), width = 12, height = 10)
  
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))
  
  # Control: sum vs weighted
  plot(log2(ctrl_strategies$sum + 1), log2(ctrl_strategies$weighted + 1),
       xlab = "Sum Aggregation [log2(count+1)]",
       ylab = "Weighted Aggregation [log2(count+1)]",
       main = sprintf("Control: Sum vs Weighted\nr = %.3f", cor_strategies_ctrl[1,2]),
       pch = 16, col = rgb(0, 0, 1, 0.3), cex = 0.5)
  abline(a = 0, b = 1, col = "red", lwd = 2, lty = 2)
  fit <- lm(log2(ctrl_strategies$weighted + 1) ~ log2(ctrl_strategies$sum + 1))
  abline(fit, col = "blue", lwd = 2)
  legend("topleft", 
         legend = c("y = x", sprintf("Linear fit (slope=%.3f)", coef(fit)[2])),
         col = c("red", "blue"), lty = c(2, 1), lwd = 2, cex = 0.8)
  
  # Mutant: sum vs weighted
  plot(log2(mut_strategies$sum + 1), log2(mut_strategies$weighted + 1),
       xlab = "Sum Aggregation [log2(count+1)]",
       ylab = "Weighted Aggregation [log2(count+1)]",
       main = sprintf("Mutant: Sum vs Weighted\nr = %.3f", cor_strategies_mut[1,2]),
       pch = 16, col = rgb(1, 0, 0, 0.3), cex = 0.5)
  abline(a = 0, b = 1, col = "red", lwd = 2, lty = 2)
  fit <- lm(log2(mut_strategies$weighted + 1) ~ log2(mut_strategies$sum + 1))
  abline(fit, col = "blue", lwd = 2)
  legend("topleft", 
         legend = c("y = x", sprintf("Linear fit (slope=%.3f)", coef(fit)[2])),
         col = c("red", "blue"), lty = c(2, 1), lwd = 2, cex = 0.8)
  
  # Ratio plots: weighted/sum
  ctrl_ratio <- ctrl_strategies$weighted / ctrl_strategies$sum
  mut_ratio <- mut_strategies$weighted / mut_strategies$sum
  
  # Control ratio distribution
  hist(ctrl_ratio[is.finite(ctrl_ratio)], 
       breaks = 50,
       col = rgb(0, 0, 1, 0.5),
       xlab = "Weighted / Sum Ratio",
       main = "Control: Strategy Ratio Distribution",
       xlim = c(0, quantile(ctrl_ratio, 0.99, na.rm = TRUE)))
  abline(v = median(ctrl_ratio, na.rm = TRUE), col = "red", lwd = 2, lty = 2)
  legend("topright", 
         legend = sprintf("Median = %.3f", median(ctrl_ratio, na.rm = TRUE)),
         col = "red", lty = 2, lwd = 2, cex = 0.8)
  
  # Mutant ratio distribution
  hist(mut_ratio[is.finite(mut_ratio)], 
       breaks = 50,
       col = rgb(1, 0, 0, 0.5),
       xlab = "Weighted / Sum Ratio",
       main = "Mutant: Strategy Ratio Distribution",
       xlim = c(0, quantile(mut_ratio, 0.99, na.rm = TRUE)))
  abline(v = median(mut_ratio, na.rm = TRUE), col = "red", lwd = 2, lty = 2)
  legend("topright", 
         legend = sprintf("Median = %.3f", median(mut_ratio, na.rm = TRUE)),
         col = "red", lty = 2, lwd = 2, cex = 0.8)
  
  dev.off()
  
  # Interpretation
  cat("Strategy Interpretation:\n")
  cat("  - 'sum': Simple sum of all 25 pixels in 5×5 buffer\n")
  cat("           → Best for detecting ANY signal in buffer region\n")
  cat("           → Robust to positional shifts\n")
  cat("           → May include more background noise\n\n")
  
  cat("  - 'weighted': Distance-weighted sum (center pixels weighted higher)\n")
  cat("           → Gives more importance to center pixel\n")
  cat("           → Balance between precision and shift-tolerance\n")
  cat("           → Values typically lower than sum\n\n")
  
  high_cor <- min(cor_strategies_ctrl[1,2], cor_strategies_mut[1,2]) > 0.95
  
  if (high_cor) {
    cat("✓ High correlation between strategies (>0.95)\n")
    cat("  → Choice of aggregation method has minimal impact on results\n")
    cat("  → Both strategies capture similar biological signal\n\n")
  } else {
    cat("⚠ Moderate correlation between strategies\n")
    cat("  → Different strategies may lead to different biological conclusions\n")
    cat("  → Consider sensitivity analysis with both methods\n\n")
  }
  
  qc_report$aggregation_strategy <- used_strategy
  qc_report$strategy_correlation <- list(
    ctrl = cor_strategies_ctrl[1,2],
    mut = cor_strategies_mut[1,2]
  )
  
  cat("Strategy comparison plots saved\n\n")
  
} else {
  cat("⚠ Aggregation strategy comparison not available\n")
  cat("  File '06_all_strategies.rds' not found or has unexpected structure\n\n")
  cat("✓ Main analysis using default 'sum' aggregation\n\n")
  qc_report$aggregation_strategy <- "sum (no comparison available)"
}

# =============================================================================
# SECTION 7: COMPREHENSIVE REPORT & RECOMMENDATIONS
# =============================================================================

cat("SECTION 7: Final QC Assessment & Recommendations\n")
cat("=================================================\n\n")

# Define pass/fail criteria
qc_pass <- TRUE
critical_issues <- c()
warnings <- c()

# Check 1: Sample correlation
if (cor_value < 0.70) {
  qc_pass <- FALSE
  critical_issues <- c(critical_issues, "Low sample correlation (<0.70)")
} else if (cor_value < 0.85) {
  warnings <- c(warnings, "Moderate sample correlation (0.70-0.85)")
}

# Check 2: Library size ratio
if (lib_ratio > 2.0) {
  warnings <- c(warnings, "Large library size difference (>2-fold)")
}

# Check 3: Data integrity
if (neg_count > 0) {
  qc_pass <- FALSE
  critical_issues <- c(critical_issues, "Negative values detected in count matrix")
}

if (na_pct > 5) {
  warnings <- c(warnings, sprintf("High NA percentage (%.1f%%)", na_pct))
}

# Check 4: Matrix centering
if (qc_report$matrix_stats$centered_ctrl_pct < 50 || 
    qc_report$matrix_stats$centered_mut_pct < 50) {
  warnings <- c(warnings, "Less than 50% of loops have centered peaks")
}

# Check 5: Center enrichment
if (qc_report$matrix_stats$mean_center_enrichment_ctrl < 0.2 ||
    qc_report$matrix_stats$mean_center_enrichment_mut < 0.2) {
  warnings <- c(warnings, "Low center enrichment (<0.2) suggests diffuse signal")
}

# Print assessment
cat("QC ASSESSMENT:\n")
cat("==============\n")

if (qc_pass && length(warnings) == 0) {
  cat("✓✓✓ PASS - All quality checks passed! ✓✓✓\n\n")
  cat("Your data is high quality and ready for differential analysis.\n")
} else if (qc_pass) {
  cat("✓ PASS (with warnings)\n\n")
  cat("Data quality is acceptable, but note the following:\n")
  for (w in warnings) {
    cat(sprintf("  ⚠ %s\n", w))
  }
} else {
  cat("✗ FAIL - Critical issues detected\n\n")
  for (issue in critical_issues) {
    cat(sprintf("  ✗ %s\n", issue))
  }
  if (length(warnings) > 0) {
    cat("\nAdditional warnings:\n")
    for (w in warnings) {
      cat(sprintf("  ⚠ %s\n", w))
    }
  }
}

cat("\n")

# Summary statistics table
cat("SUMMARY STATISTICS:\n")
cat("===================\n")
cat(sprintf("Dataset:             %d loops × %d samples\n", nrow(counts_matrix), ncol(counts_matrix)))
cat(sprintf("Sample correlation:  %.3f (Pearson)\n", cor_value))
cat(sprintf("Library size ratio:  %.2f (mut/ctrl)\n", lib_ratio))
cat(sprintf("Mean log2 FC:        %.3f ± %.3f\n", mean(M, na.rm = TRUE), sd(M, na.rm = TRUE)))
cat(sprintf("Loops with |M|>1:    %d (%.1f%%)\n", 
            sum(abs(M) > 1, na.rm = TRUE),
            100 * sum(abs(M) > 1, na.rm = TRUE) / length(M)))
cat(sprintf("Positional shifts:   %d (%.1f%%)\n",
            qc_report$matrix_stats$shifts_pct * nrow(counts_matrix) / 100,
            qc_report$matrix_stats$shifts_pct))
cat(sprintf("Center enrichment:   %.3f (ctrl), %.3f (mut)\n",
            qc_report$matrix_stats$mean_center_enrichment_ctrl,
            qc_report$matrix_stats$mean_center_enrichment_mut))
cat("\n")

# Recommendations
cat("RECOMMENDATIONS:\n")
cat("================\n")

if (qc_pass) {
  cat("1. ✓ Proceed with edgeR differential analysis\n")
  cat("2. ✓ Data quality supports biological interpretation\n")
  
  if (qc_report$matrix_stats$shifts_pct > 10) {
    cat(sprintf("3. ⚠ %.1f%% of loops show positional shifts - buffer approach validated\n",
                qc_report$matrix_stats$shifts_pct))
  }
  
  if (cor_value > 0.95) {
    cat("4. ⚠ Very high correlation suggests limited biological differences\n")
    cat("   → May need larger sample size or deeper sequencing to detect subtle changes\n")
  }
  
  cat("\nNEXT STEPS:\n")
  cat("-----------\n")
  cat("→ Run edgeR differential analysis script\n")
  cat("→ Compare results with Hiccups differential calls\n")
  cat("→ Validate top hits with biological context\n")
  
} else {
  cat("1. ✗ Investigate critical issues before proceeding\n")
  cat("2. ✗ Consider re-running pipeline with adjusted parameters\n")
  cat("3. ✗ Check input .hic file quality\n")
  
  cat("\nTROUBLESHOOTING:\n")
  cat("----------------\n")
  if (cor_value < 0.70) {
    cat("→ Low correlation may indicate:\n")
    cat("  - Batch effects between samples\n")
    cat("  - Different sequencing depths\n")
    cat("  - Sample mix-up\n")
    cat("  Solution: Check raw .hic files and sequencing metrics\n\n")
  }
  
  if (neg_count > 0) {
    cat("→ Negative values indicate data corruption\n")
    cat("  Solution: Re-run extraction step\n\n")
  }
}

# Save QC report object
saveRDS(qc_report, file.path(qc_dir, "qc_report_summary.rds"))

cat("\n")
cat("=============================================================================\n")
cat("QC REPORT COMPLETE\n")
cat("=============================================================================\n")
cat(sprintf("All outputs saved to: %s\n", qc_dir))
cat("\nGenerated files:\n")
cat("  01_sample_correlation.pdf       - Sample similarity heatmap\n")
cat("  02_sample_scatter.pdf           - Sample-to-sample scatter plot\n")
cat("  03_count_distributions.pdf      - Distribution comparisons\n")
cat("  04_ma_plot.pdf                  - MA plot (fold change analysis)\n")
cat("  05_variable_loops_heatmap.pdf   - Top variable loops\n")
cat("  06_chromosome_distribution.pdf  - Loops per chromosome\n")
cat("  07_distance_decay.pdf           - Distance vs signal strength\n")
cat("  08_example_matrices.pdf         - Individual 5×5 matrix examples\n")
cat("  09_strategy_comparison.pdf      - Aggregation method comparison\n")
cat("  qc_report_summary.rds           - Complete QC metrics (R object)\n")
cat("\n")

if (qc_pass) {
  cat("✓ Ready for differential analysis!\n")
} else {
  cat("✗ Please address critical issues before proceeding\n")
}

cat("=============================================================================\n")
