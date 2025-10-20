#!/usr/bin/env Rscript
# /expanse/lustre/projects/csd940/zalibhai/mariner/scripts/06_aggregate.R
# Purpose: Aggregate 5x5 Hi-C matrices to single counts per loop, handling bin-shifts
# Input: 05_extracted HDF5 InteractionArray (5x5x150x2)
# Output: Count matrix (150 loops x 2 samples) ready for edgeR

cat("\n=================================\n")
cat("Step 6: Matrix Aggregation\n")
cat("=================================\n\n")

# Load required libraries (removed DelayedMatrixStats)
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

# Load extracted matrices
cat("Loading extracted Hi-C matrices...\n")
pixels <- loadHDF5SummarizedExperiment(
  dir = "outputs/full",
  prefix = "05_extracted"
)

# Load metadata
metadata <- readRDS("outputs/full/05_metadata.rds")
cat(sprintf("  Loaded: %d loops x %d samples\n", metadata$n_loops, metadata$n_files))
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

# Compare sample correlations for each strategy
cat("Control vs Mutant correlation by strategy:\n")
for (name in names(strategies)) {
  mat <- strategies[[name]]
  cor_val <- cor(mat[,1], mat[,2], use = "complete.obs")
  cat(sprintf("  %s: %.3f\n", name, cor_val))
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

# Check for loops with extreme differences
log2fc <- log2((final_counts[,2] + 1) / (final_counts[,1] + 1))
extreme_changes <- which(abs(log2fc) > 2)
cat(sprintf("Loops with |log2FC| > 2: %d (%.1f%%)\n", 
            length(extreme_changes), 
            100 * length(extreme_changes) / nrow(final_counts)))

# Distribution check
cat("\nCount distribution check:\n")
ctrl_positive <- final_counts[final_counts[,1] > 0, 1]
mut_positive <- final_counts[final_counts[,2] > 0, 2]

# Basic distribution statistics instead of Shapiro-Wilk test
cat("  Control distribution:\n")
cat(sprintf("    Non-zero: %d / %d (%.1f%%)\n", 
            length(ctrl_positive), nrow(final_counts),
            100 * length(ctrl_positive) / nrow(final_counts)))
if (length(ctrl_positive) > 0) {
  cat(sprintf("    Range: %.1f - %.1f\n", min(ctrl_positive), max(ctrl_positive)))
  cat(sprintf("    Median: %.1f, Mean: %.1f\n", 
              median(ctrl_positive), mean(ctrl_positive)))
}

cat("  Mutant distribution:\n")
cat(sprintf("    Non-zero: %d / %d (%.1f%%)\n", 
            length(mut_positive), nrow(final_counts),
            100 * length(mut_positive) / nrow(final_counts)))
if (length(mut_positive) > 0) {
  cat(sprintf("    Range: %.1f - %.1f\n", min(mut_positive), max(mut_positive)))
  cat(sprintf("    Median: %.1f, Mean: %.1f\n", 
              median(mut_positive), mean(mut_positive)))
}

# Signal strength analysis
cat("\nSignal strength by sample:\n")
cat(sprintf("  Control: median = %.1f, mean = %.1f\n",
            median(final_counts[,1]), mean(final_counts[,1])))
cat(sprintf("  Mutant:  median = %.1f, mean = %.1f\n",
            median(final_counts[,2]), mean(final_counts[,2])))

# Check for systematic bias
mean_diff <- mean(final_counts[,2]) - mean(final_counts[,1])
cat(sprintf("\nMean difference (mut - ctrl): %.1f\n", mean_diff))
if (abs(mean_diff) > 0.5 * mean(final_counts[,1])) {
  cat(" Large systematic difference detected - may need normalization\n")
} else {
  cat(" No major systematic bias\n")
}

# Sample correlation
final_cor <- cor(final_counts[,1], final_counts[,2], use = "complete.obs")
cat(sprintf("\nFinal sample correlation: %.3f\n", final_cor))

# ============================================================================
# SAVE OUTPUTS
# ============================================================================
cat("\n=================================\n")
cat("Saving Results\n")
cat("=================================\n\n")

output_dir <- "outputs/full"

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
    binSize = 5000,
    buffer = 2,
    date = Sys.Date(),
    sample_correlation = final_cor
  )
)
saveRDS(all_strategies, file.path(output_dir, "06_all_strategies.rds"))
cat(sprintf("  All strategies saved: %s/06_all_strategies.rds\n", output_dir))

# ============================================================================
# Create edgeR-ready object
# ============================================================================
cat("\n=== Preparing edgeR Input ===\n")

# Load the original merged loops to get genomic coordinates
original_loops <- readRDS("outputs/test/02_merged.rds")

# Extract genomic information for the loops
# Note: We need to handle the InteractionSet structure properly
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

# Create DGEList
y <- DGEList(
  counts = final_counts,
  group = factor(c("ctrl", "mut")),
  genes = gene_info
)

# Save edgeR object
saveRDS(y, file.path(output_dir, "06_edgeR_input.rds"))
cat("  edgeR DGEList saved: outputs/test/06_edgeR_input.rds\n")

# edgeR preview
cat("\nEdgeR object preview:\n")
cat(sprintf("  Samples: %s\n", paste(colnames(y), collapse = ", ")))
cat(sprintf("  Groups: %s\n", paste(levels(y$samples$group), collapse = ", ")))
cat(sprintf("  Library sizes: %s\n", 
            paste(format(y$samples$lib.size, big.mark = ","), collapse = ", ")))

# Filtering preview (what edgeR would typically do)
keep <- filterByExpr(y, group = y$samples$group)
cat(sprintf("\nFiltering preview:\n"))
cat(sprintf("  Loops passing filterByExpr: %d / %d (%.1f%%)\n",
            sum(keep), length(keep), 100 * sum(keep) / length(keep)))

if (sum(keep) > 0) {
  cat(sprintf("  Filtered library sizes would be: %s\n",
              paste(format(colSums(y$counts[keep,]), big.mark = ","), collapse = ", ")))
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
cat(sprintf("  Input: %d loops with 5x5 pixel buffers\n", dims[3]))
cat(sprintf("  Output: %d x %d count matrix\n", 
            nrow(final_counts), ncol(final_counts)))
cat(sprintf("  Primary method: Sum aggregation\n"))
cat(sprintf("  Sample correlation: %.3f\n", final_cor))
cat(sprintf("  Loops for analysis: %d (after filtering)\n", sum(keep)))

cat("\nOutput files:\n")
cat("  1. outputs/test/06_counts_matrix.tsv - Tab-separated counts\n")
cat("  2. outputs/test/06_counts_matrix.rds - R object\n")
cat("  3. outputs/test/06_edgeR_input.rds - Ready for differential analysis\n")
cat("  4. outputs/test/06_all_strategies.rds - All aggregation methods\n")

cat("\n" )
cat(paste(rep("=", 50), collapse = ""))
cat("\n")
cat("NEXT STEPS:\n")
cat("1. Review count distributions and sample correlation\n")
cat("2. If satisfied, proceed with edgeR differential analysis\n")
cat("3. For full dataset: Update scripts to process all loops (not just 100)\n")
cat("\nTo start edgeR analysis:\n")
cat("  y <- readRDS('outputs/test/06_edgeR_input.rds')\n")
cat("  # Then proceed with normalization, dispersion, and testing\n")
