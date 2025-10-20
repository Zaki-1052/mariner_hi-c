#!/usr/bin/env Rscript
# /expanse/lustre/projects/csd940/zalibhai/mariner/scripts/06_aggregate.R
# Purpose: Aggregate 5x5 Hi-C matrices to single counts per loop, handling bin-shifts
# Multi-resolution support: accepts resolution as command-line argument
# Input: 05_extracted HDF5 InteractionArray (5x5xN_loopsx6)
# Output: Count matrix (N_loops x 6 samples) ready for edgeR

# Parse command-line arguments for resolution
args <- commandArgs(trailingOnly = TRUE)
RESOLUTION <- if (length(args) > 0) as.numeric(args[1]) else 5000

cat("\n========================================\n")
cat("Step 6: Matrix Aggregation\n")
cat(sprintf("RESOLUTION: %d bp (%d kb)\n", RESOLUTION, RESOLUTION/1000))
cat("========================================\n\n")

# Load required libraries
suppressPackageStartupMessages({
  library(mariner)
  library(InteractionSet)
  library(HDF5Array)
  library(DelayedArray)
  library(Matrix)
  library(edgeR)  # Load here since we'll use it at the end
})

# Set up paths
base_dir <- "/expanse/lustre/projects/csd940/zalibhai/mariner"
setwd(base_dir)

# Resolution-specific directories
input_dir <- sprintf("outputs/res_%dkb", RESOLUTION/1000)
output_dir <- input_dir

# Load extracted matrices
cat("Loading extracted Hi-C matrices...\n")
pixels <- loadHDF5SummarizedExperiment(
  dir = input_dir,
  prefix = "05_extracted"
)

# Load metadata
metadata <- readRDS(file.path(input_dir, "05_metadata.rds"))
cat(sprintf("  Loaded: %d loops x %d samples\n", metadata$n_loops, metadata$n_replicates))
cat(sprintf("  Matrix dimensions per loop: 5x5 pixels\n"))
cat(sprintf("  Files: %s\n", paste(names(metadata$files), collapse=", ")))

# Extract the count array
count_array <- counts(pixels)
dims <- dim(count_array)
cat(sprintf("\nArray dimensions: %d x %d x %d x %d\n",
            dims[1], dims[2], dims[3], dims[4]))

# ============================================================================
# AGGREGATION STRATEGY 1: Simple Sum (Captures total signal)
# ============================================================================
cat("\n=== Strategy 1: Simple Sum Aggregation ===\n")
cat("Rationale: Captures total Hi-C signal regardless of position within 5x5 window\n")

# Initialize matrix for results
counts_sum <- matrix(NA, nrow = dims[3], ncol = dims[4])
colnames(counts_sum) <- names(metadata$files)
rownames(counts_sum) <- paste0("loop_", 1:dims[3])

# Aggregate each loop
cat("Aggregating...")
for (i in 1:dims[3]) {
  if (i %% 50 == 0) cat(sprintf(" %d", i))
  for (j in 1:dims[4]) {
    # Extract 5x5 matrix for this loop and sample
    mat_5x5 <- as.matrix(count_array[, , i, j])
    # Sum all pixels, removing NAs
    counts_sum[i, j] <- sum(mat_5x5, na.rm = TRUE)
  }
}
cat(" Done\n")

cat(sprintf("  Aggregated to: %d loops x %d samples\n",
            nrow(counts_sum), ncol(counts_sum)))

# Validation
sum_stats <- summary(as.vector(counts_sum))
cat("\nSum aggregation statistics:\n")
print(sum_stats)

# ============================================================================
# AGGREGATION STRATEGY 2: Center-Weighted (Emphasizes expected position)
# ============================================================================
cat("\n=== Strategy 2: Center-Weighted Aggregation ===\n")
cat("Rationale: Higher weight to center pixels, lower to edges\n")
cat("Weight matrix (normalized to sum=1):\n")

# Define weight matrix - Gaussian-like weights centered on middle pixel
weight_matrix <- matrix(c(
  0.04, 0.06, 0.08, 0.06, 0.04,  # Edge pixels: lowest weight
  0.06, 0.08, 0.10, 0.08, 0.06,  # Inner ring: medium weight
  0.08, 0.10, 0.16, 0.10, 0.08,  # Center row: highest in middle
  0.06, 0.08, 0.10, 0.08, 0.06,  # Inner ring: medium weight
  0.04, 0.06, 0.08, 0.06, 0.04   # Edge pixels: lowest weight
), nrow = 5, byrow = TRUE)

# Verify weights sum to 1
weight_matrix <- weight_matrix / sum(weight_matrix)

# Display weight matrix
for (row in 1:5) {
  cat(sprintf("  %.3f %.3f %.3f %.3f %.3f\n",
              weight_matrix[row,1], weight_matrix[row,2],
              weight_matrix[row,3], weight_matrix[row,4],
              weight_matrix[row,5]))
}

# Apply weighted aggregation
counts_weighted <- matrix(NA, nrow = dims[3], ncol = dims[4])
colnames(counts_weighted) <- names(metadata$files)
rownames(counts_weighted) <- paste0("loop_", 1:dims[3])

cat("\nApplying weighted aggregation...")
for (i in 1:dims[3]) {
  if (i %% 50 == 0) cat(sprintf(" %d", i))
  for (j in 1:dims[4]) {
    mat_5x5 <- as.matrix(count_array[, , i, j])
    # Replace NAs with 0 for multiplication
    mat_5x5[is.na(mat_5x5)] <- 0
    # Weighted sum
    counts_weighted[i, j] <- sum(mat_5x5 * weight_matrix)
  }
}
cat(" Done\n")

weighted_stats <- summary(as.vector(counts_weighted))
cat("\nWeighted aggregation statistics:\n")
print(weighted_stats)

# ============================================================================
# COMPARISON OF AGGREGATION STRATEGIES
# ============================================================================
cat("\n=================================\n")
cat("Strategy Comparison\n")
cat("=================================\n\n")

# Calculate correlations between strategies
strategies <- list(
  sum = counts_sum,
  weighted = counts_weighted
)

# Replicate-aware correlation analysis
cat("Replicate correlation analysis:\n\n")

# Build correlation matrix for sum strategy
n_samples <- dims[4]
cor_matrix <- matrix(NA, nrow = n_samples, ncol = n_samples)
rownames(cor_matrix) <- names(metadata$files)
colnames(cor_matrix) <- names(metadata$files)

for (i in 1:n_samples) {
  for (j in 1:n_samples) {
    cor_matrix[i, j] <- cor(counts_sum[,i], counts_sum[,j], use = "complete.obs")
  }
}

# Display correlation matrix
cat("Full correlation matrix (sum strategy):\n")
print(round(cor_matrix, 3))

# Within-condition correlations
cat("\nWithin-condition correlations:\n")
cat("  Control replicates:\n")
ctrl_cors <- numeric()
for (i in 1:2) {
  for (j in (i+1):3) {
    cor_val <- cor_matrix[i, j]
    ctrl_cors <- c(ctrl_cors, cor_val)
    cat(sprintf("    %s vs %s: %.3f\n",
                names(metadata$files)[i], names(metadata$files)[j], cor_val))
  }
}
ctrl_within_mean <- mean(ctrl_cors)
cat(sprintf("  Mean within-ctrl correlation: %.3f\n\n", ctrl_within_mean))

cat("  Mutant replicates:\n")
mut_cors <- numeric()
for (i in 4:5) {
  for (j in (i+1):6) {
    cor_val <- cor_matrix[i, j]
    mut_cors <- c(mut_cors, cor_val)
    cat(sprintf("    %s vs %s: %.3f\n",
                names(metadata$files)[i], names(metadata$files)[j], cor_val))
  }
}
mut_within_mean <- mean(mut_cors)
cat(sprintf("  Mean within-mut correlation: %.3f\n\n", mut_within_mean))

# Between-condition correlations
cat("Between-condition correlations:\n")
between_cors <- as.vector(cor_matrix[1:3, 4:6])
cat(sprintf("  Mean ctrl vs mut correlation: %.3f\n", mean(between_cors)))
cat(sprintf("  Range: %.3f - %.3f\n\n", min(between_cors), max(between_cors)))

# Quality assessment
cat("Quality Assessment:\n")
if (ctrl_within_mean > 0.95 && mut_within_mean > 0.95) {
  cat("  ✓ Excellent within-condition reproducibility (>0.95)\n")
} else if (ctrl_within_mean > 0.90 && mut_within_mean > 0.90) {
  cat("  ✓ Good within-condition reproducibility (>0.90)\n")
} else {
  cat("  ⚠ WARNING: Lower than expected within-condition correlation (<0.90)\n")
  cat("    This may indicate technical issues or high biological variability\n")
}

if ((ctrl_within_mean - mean(between_cors)) > 0.05 ||
    (mut_within_mean - mean(between_cors)) > 0.05) {
  cat("  ✓ Clear biological signal (within > between correlations)\n")
} else {
  cat("  ⚠ WARNING: Within-condition correlations not much higher than between\n")
  cat("    This may indicate weak biological differences\n")
}

# Correlation between strategies
cor_strategies <- cor(as.vector(counts_sum), as.vector(counts_weighted), use = "complete.obs")
cat(sprintf("\nCorrelation between sum and weighted strategies: %.3f\n", cor_strategies))

# ============================================================================
# BIOLOGICAL VALIDATION
# ============================================================================
cat("\n=================================\n")
cat("Biological Validation\n")
cat("=================================\n\n")

# Use sum aggregation as primary (most interpretable for edgeR)
final_counts <- counts_sum

# Calculate mean counts per condition
ctrl_means <- rowMeans(final_counts[, 1:3])
mut_means <- rowMeans(final_counts[, 4:6])

# Check for loops with extreme differences
log2fc <- log2((mut_means + 1) / (ctrl_means + 1))
extreme_changes <- which(abs(log2fc) > 2)
cat(sprintf("Loops with |log2FC| > 2: %d (%.1f%%)\n",
            length(extreme_changes),
            100 * length(extreme_changes) / nrow(final_counts)))

# Distribution check per replicate
cat("\nCount distribution by replicate:\n")

for (i in 1:dims[4]) {
  sample_name <- names(metadata$files)[i]
  sample_counts <- final_counts[, i]
  positive <- sample_counts[sample_counts > 0]

  cat(sprintf("\n  %s:\n", sample_name))
  cat(sprintf("    Non-zero: %d / %d (%.1f%%)\n",
              length(positive), nrow(final_counts),
              100 * length(positive) / nrow(final_counts)))
  if (length(positive) > 0) {
    cat(sprintf("    Range: %.1f - %.1f\n", min(positive), max(positive)))
    cat(sprintf("    Median: %.1f, Mean: %.1f\n",
                median(positive), mean(positive)))
  }
}

# Signal strength analysis by group
cat("\n\nSignal strength by group:\n")
ctrl_all <- as.vector(final_counts[, 1:3])
mut_all <- as.vector(final_counts[, 4:6])

cat(sprintf("  Control (all replicates): median = %.1f, mean = %.1f\n",
            median(ctrl_all), mean(ctrl_all)))
cat(sprintf("  Mutant (all replicates):  median = %.1f, mean = %.1f\n",
            median(mut_all), mean(mut_all)))

# Check for systematic bias
mean_diff <- mean(mut_all) - mean(ctrl_all)
cat(sprintf("\nMean difference (mut - ctrl): %.1f\n", mean_diff))
if (abs(mean_diff) > 0.5 * mean(ctrl_all)) {
  cat("  ⚠ Large systematic difference detected - edgeR normalization will handle this\n")
} else {
  cat("  ✓ No major systematic bias\n")
}

# ============================================================================
# SAVE OUTPUTS
# ============================================================================
cat("\n=================================\n")
cat("Saving Results\n")
cat("=================================\n\n")

# Save primary count matrix (sum aggregation)
write.table(final_counts,
            file = file.path(output_dir, "06_counts_matrix.tsv"),
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
cat(sprintf("  Count matrix saved: %s/06_counts_matrix.tsv\n", output_dir))

# Save as RDS for easy loading in R
saveRDS(final_counts, file.path(output_dir, "06_counts_matrix.rds"))
cat(sprintf("  RDS object saved: %s/06_counts_matrix.rds\n", output_dir))

# Save all strategies for comparison
all_strategies <- list(
  counts = strategies,
  metadata = list(
    n_loops = nrow(final_counts),
    samples = colnames(final_counts),
    aggregation_method = "sum",
    extraction_norm = metadata$normalization,
    binSize = RESOLUTION,
    buffer = 2,
    date = Sys.Date(),
    within_ctrl_correlation = ctrl_within_mean,
    within_mut_correlation = mut_within_mean,
    between_correlation = mean(between_cors)
  )
)
saveRDS(all_strategies, file.path(output_dir, "06_all_strategies.rds"))
cat(sprintf("  All strategies saved: %s/06_all_strategies.rds\n", output_dir))

# ============================================================================
# Create edgeR-ready object
# ============================================================================
cat("\n=== Preparing edgeR Input ===\n")

# Load the original merged loops to get genomic coordinates
original_loops <- readRDS(file.path(input_dir, "02_merged.rds"))

# Extract genomic information for the loops
loop_anchors1 <- anchors(original_loops, type = "first")
loop_anchors2 <- anchors(original_loops, type = "second")

# Create gene annotation dataframe for edgeR
gene_info <- data.frame(
  loop_id = rownames(final_counts),
  chr1 = as.character(seqnames(loop_anchors1)[1:nrow(final_counts)]),
  start1 = start(loop_anchors1)[1:nrow(final_counts)],
  end1 = end(loop_anchors1)[1:nrow(final_counts)],
  chr2 = as.character(seqnames(loop_anchors2)[1:nrow(final_counts)]),
  start2 = start(loop_anchors2)[1:nrow(final_counts)],
  end2 = end(loop_anchors2)[1:nrow(final_counts)],
  stringsAsFactors = FALSE
)

# Calculate interaction distance (intrachromosomal only)
gene_info$distance <- ifelse(
  gene_info$chr1 == gene_info$chr2,
  abs(gene_info$start2 - gene_info$start1),
  NA  # NA for interchromosomal
)

# Create DGEList with proper replicate structure
y <- DGEList(
  counts = final_counts,
  group = factor(c("ctrl", "ctrl", "ctrl", "mut", "mut", "mut")),
  genes = gene_info
)

# Save edgeR object
saveRDS(y, file.path(output_dir, "06_edgeR_input.rds"))
cat(sprintf("  edgeR DGEList saved: %s/06_edgeR_input.rds\n", output_dir))

# edgeR preview
cat("\nEdgeR object preview:\n")
cat(sprintf("  Samples: %s\n", paste(colnames(y), collapse = ", ")))
cat(sprintf("  Groups: %s\n", paste(levels(y$samples$group), collapse = ", ")))
cat(sprintf("  Replicates per group: n = %d\n", sum(y$samples$group == "ctrl")))
cat(sprintf("  Library sizes:\n"))
for (i in 1:ncol(y)) {
  cat(sprintf("    %s: %s\n", colnames(y)[i], format(y$samples$lib.size[i], big.mark = ",")))
}

# Filtering preview (what edgeR would typically do)
keep <- filterByExpr(y, group = y$samples$group)
cat(sprintf("\nFiltering preview:\n"))
cat(sprintf("  Loops passing filterByExpr: %d / %d (%.1f%%)\n",
            sum(keep), length(keep), 100 * sum(keep) / length(keep)))

if (sum(keep) > 0) {
  cat(sprintf("  Filtered library sizes would be:\n"))
  for (i in 1:ncol(y)) {
    cat(sprintf("    %s: %s\n",
                colnames(y)[i],
                format(sum(y$counts[keep, i]), big.mark = ",")))
  }
}

# Show a few example loops
cat("\nExample loop coordinates:\n")
print(head(gene_info, 3))

# ============================================================================
# FINAL SUMMARY
# ============================================================================
cat("\n=================================\n")
cat("Step 6: Aggregation Complete\n")
cat("=================================\n")

cat("\nPipeline Summary:\n")
cat(sprintf("  Resolution: %d bp (%d kb)\n", RESOLUTION, RESOLUTION/1000))
cat(sprintf("  Input: %d loops with 5x5 pixel buffers\n", dims[3]))
cat(sprintf("  Output: %d x %d count matrix\n",
            nrow(final_counts), ncol(final_counts)))
cat(sprintf("  Primary method: Sum aggregation\n"))
cat(sprintf("  Within-ctrl correlation: %.3f\n", ctrl_within_mean))
cat(sprintf("  Within-mut correlation: %.3f\n", mut_within_mean))
cat(sprintf("  Between-group correlation: %.3f\n", mean(between_cors)))
cat(sprintf("  Loops for analysis: %d (after filtering)\n", sum(keep)))

cat("\nOutput files:\n")
cat(sprintf("  1. %s/06_counts_matrix.tsv - Tab-separated counts\n", output_dir))
cat(sprintf("  2. %s/06_counts_matrix.rds - R object\n", output_dir))
cat(sprintf("  3. %s/06_edgeR_input.rds - Ready for differential analysis\n", output_dir))
cat(sprintf("  4. %s/06_all_strategies.rds - All aggregation methods\n", output_dir))

cat("\n" )
cat(paste(rep("=", 50), collapse = ""))
cat("\n")
cat("NEXT STEPS:\n")
cat("1. Review replicate correlations and quality metrics\n")
cat("2. Proceed with edgeR differential analysis using biological replicates\n")
cat(sprintf("3. Expected: Much higher statistical power than n=1 analysis\n"))
cat("\nTo start edgeR analysis:\n")
cat(sprintf("  Rscript scripts/edgeR.R %d\n", RESOLUTION))
cat("  # Uses quasi-likelihood GLM with n=3 replicates per group\n")
