#!/usr/bin/env Rscript
# stripes/scripts/phase2_quantification.R
# Phase 2: Hi-C Matrix Extraction and Count Aggregation for Stripe Quantification
#
# Converts stripes to GInteractions (anchor + 100kb sampling point), extracts
# 5x5 Hi-C matrices at sampling positions, and aggregates to count matrix.
#
# Usage: Rscript phase2_quantification.R <timepoint>
#   timepoint: "early" (250831) or "late" (250402)

suppressPackageStartupMessages({
  library(mariner)
  library(InteractionSet)
  library(GenomicRanges)
  library(IRanges)
  library(strawr)
  library(HDF5Array)
  library(DelayedArray)
  library(yaml)
  library(dplyr)
})

# ==============================================================================
# PARSE ARGUMENTS
# ==============================================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript phase2_quantification.R <timepoint>\n",
       "  timepoint: 'early' or 'late'")
}
TIMEPOINT <- args[1]

if (!TIMEPOINT %in% c("early", "late")) {
  stop("Timepoint must be 'early' or 'late'")
}

cat("\n========================================\n")
cat("Phase 2: Hi-C Quantification\n")
cat(sprintf("Timepoint: %s\n", TIMEPOINT))
cat("========================================\n\n")

# ==============================================================================
# LOAD CONFIGURATION
# ==============================================================================
cat("Loading configuration...\n")

config_file <- "config/stripe_config.yaml"
if (!file.exists(config_file)) {
  stop(sprintf("Config file not found: %s\n  Run from project base directory",
               config_file))
}

config <- yaml::read_yaml(config_file)
cat(sprintf("  Config: %s\n", config_file))
cat(sprintf("  Project base: %s\n", config$project$base_dir))

# ==============================================================================
# HELPER: Build Paths from Config
# ==============================================================================
get_output_dir <- function(timepoint, resolution_kb, config) {
  file.path(
    config$outputs$base_dir,
    timepoint,
    paste0("res_", resolution_kb, "kb")
  )
}

get_hic_file_path <- function(timepoint, sample, config) {
  timepoint_dir <- config$stripe_data$timepoints[[timepoint]]
  file.path(
    config$hic_files$base_path,
    timepoint_dir,
    paste0(sample, ".hic")
  )
}

# ==============================================================================
# CONFIGURATION PARAMETERS
# ==============================================================================
timepoint_dir <- config$stripe_data$timepoints[[TIMEPOINT]]
resolution <- config$stripe_data$resolutions$primary  # 5000 bp
resolution_kb <- resolution / 1000
sampling_distance <- config$quantification$sampling_distance_bp  # 100000 bp
buffer_bins <- config$quantification$buffer_bins  # 2
norm_preferred <- config$quantification$normalization  # "KR"
norm_fallback <- config$quantification$fallback_norm  # "VC"

cat(sprintf("\nParameters:\n"))
cat(sprintf("  Timepoint dir: %s\n", timepoint_dir))
cat(sprintf("  Resolution: %d bp (%d kb)\n", resolution, resolution_kb))
cat(sprintf("  Sampling distance: %d bp\n", sampling_distance))
cat(sprintf("  Buffer: %d bins (= %dx%d matrix)\n",
            buffer_bins, 2*buffer_bins + 1, 2*buffer_bins + 1))
cat(sprintf("  Normalization: %s (fallback: %s)\n", norm_preferred, norm_fallback))

# ==============================================================================
# LOAD PHASE 1 OUTPUT
# ==============================================================================
cat("\n=== Loading Phase 1 Unified Stripes ===\n")

input_dir <- get_output_dir(TIMEPOINT, resolution_kb, config)
stripe_file <- file.path(input_dir, "01_unified_stripes.tsv")

if (!file.exists(stripe_file)) {
  stop(sprintf("Phase 1 output not found: %s\n  Run phase1_detection.R first",
               stripe_file))
}

stripes <- read.table(stripe_file, header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE)
cat(sprintf("  Loaded: %d stripes from %s\n", nrow(stripes), stripe_file))

# Validate required columns
required_cols <- c("stripe_id", "chr", "anchor_x1", "anchor_x2",
                   "span_y1", "span_y2", "direction_type", "source",
                   "anchor_center")
missing <- setdiff(required_cols, colnames(stripes))
if (length(missing) > 0) {
  stop(sprintf("Missing required columns: %s", paste(missing, collapse = ", ")))
}

cat(sprintf("  Stripe counts by source:\n"))
print(table(stripes$source))

# ==============================================================================
# STRIPE TO GINTERACTIONS CONVERSION
# ==============================================================================
cat("\n=== Converting Stripes to GInteractions ===\n")

# Calculate sampling point based on direction
# 3'prime stripes: sample downstream (+ distance)
# 5'prime stripes: sample upstream (- distance)
stripes$sampling_point <- ifelse(
  stripes$direction_type == "3prime",
  stripes$anchor_center + sampling_distance,
  stripes$anchor_center - sampling_distance
)

# Validate sampling points are within stripe span
stripes$in_span <- ifelse(
  stripes$direction_type == "3prime",
  stripes$sampling_point >= stripes$span_y1 & stripes$sampling_point <= stripes$span_y2,
  stripes$sampling_point >= stripes$span_y1 & stripes$sampling_point <= stripes$span_y2
)

cat(sprintf("  Sampling points within span: %d / %d (%.1f%%)\n",
            sum(stripes$in_span), nrow(stripes),
            100 * mean(stripes$in_span)))

# Create GRanges for anchor (region) and sampling point (single position)
anchor1_gr <- GRanges(
  seqnames = stripes$chr,
  ranges = IRanges(start = stripes$anchor_x1, end = stripes$anchor_x2)
)

# For anchor2, use a small region at the sampling point
anchor2_gr <- GRanges(
  seqnames = stripes$chr,
  ranges = IRanges(start = stripes$sampling_point,
                   end = stripes$sampling_point + 1)
)

# Create GInteractions
gi <- GInteractions(anchor1_gr, anchor2_gr)

# Add stripe metadata
mcols(gi)$stripe_id <- stripes$stripe_id
mcols(gi)$source <- stripes$source
mcols(gi)$direction_type <- stripes$direction_type
mcols(gi)$confidence <- stripes$confidence

cat(sprintf("  Created GInteractions: %d interactions\n", length(gi)))

# ==============================================================================
# BINNING AND BUFFERING
# ==============================================================================
cat("\n=== Binning to Resolution Grid ===\n")

# Assign to resolution bins
binned <- assignToBins(gi, binSize = resolution, pos1 = "center", pos2 = "center")
cat(sprintf("  Binned to %d bp resolution\n", resolution))

# Add buffer for 5x5 matrix extraction
buffered <- pixelsToMatrices(binned, buffer = buffer_bins)
cat(sprintf("  Added buffer: %dx%d pixel regions\n",
            2*buffer_bins + 1, 2*buffer_bins + 1))

# ==============================================================================
# LOAD HI-C FILES
# ==============================================================================
cat("\n=== Loading Hi-C Files ===\n")

sample_names <- config$samples$names
hicFiles <- c()

for (sample in sample_names) {
  filepath <- get_hic_file_path(TIMEPOINT, sample, config)
  hicFiles[sample] <- filepath
}

# Verify files exist
cat("Verifying Hi-C files:\n")
for (name in names(hicFiles)) {
  filepath <- hicFiles[name]
  if (!file.exists(filepath)) {
    stop(sprintf("ERROR: %s file not found at %s", name, filepath))
  }
  size_gb <- file.info(filepath)$size / 1e9
  cat(sprintf("  %s: %.1f GB - OK\n", name, size_gb))
}

# Check available normalizations
cat("\nChecking available normalizations:\n")
first_file <- hicFiles[1]
available_norms <- readHicNormTypes(first_file)
cat(sprintf("  Available: %s\n", paste(available_norms, collapse = ", ")))

# Determine which normalization to use
if (norm_preferred %in% available_norms) {
  norm_to_use <- norm_preferred
  cat(sprintf("  Using: %s (preferred)\n", norm_to_use))
} else if (norm_fallback %in% available_norms) {
  norm_to_use <- norm_fallback
  cat(sprintf("  Using: %s (fallback - %s not available)\n",
              norm_to_use, norm_preferred))
} else {
  stop(sprintf("Neither %s nor %s normalization available",
               norm_preferred, norm_fallback))
}

# Check resolution availability
available_res <- readHicBpResolutions(first_file)
if (!resolution %in% available_res) {
  stop(sprintf("%d bp resolution not available. Available: %s",
               resolution, paste(head(available_res, 5), collapse = ", ")))
}

# ==============================================================================
# EXTRACT HI-C MATRICES
# ==============================================================================
cat("\n=== Extracting Hi-C Matrices ===\n")

# Set up HDF5 file path
hdf5_dir <- file.path(input_dir, "temp_hdf5")

# Clean up existing temp HDF5 directory
if (dir.exists(hdf5_dir)) {
  cat("Cleaning up existing temp HDF5 directory...\n")
  unlink(hdf5_dir, recursive = TRUE)
}
dir.create(hdf5_dir, recursive = TRUE)

h5_file_path <- file.path(hdf5_dir, "stripe_matrices.h5")
cat(sprintf("  HDF5 output: %s\n", h5_file_path))

cat(sprintf("\nExtracting %d stripes x %d samples x %d pixels...\n",
            length(buffered), length(hicFiles),
            (2*buffer_bins + 1)^2))

start_time <- Sys.time()

pixels <- pullHicMatrices(
  x = buffered,
  files = hicFiles,
  binSize = resolution,
  h5File = h5_file_path,
  norm = norm_to_use,
  matrix = "observed",
  blockSize = config$quantification$block_size,
  onDisk = config$quantification$on_disk,
  compressionLevel = config$quantification$compression_level
)

end_time <- Sys.time()
extraction_time <- difftime(end_time, start_time, units = "secs")
cat(sprintf("\nExtraction completed in %.1f seconds\n", extraction_time))

# ==============================================================================
# VALIDATION OF EXTRACTED DATA
# ==============================================================================
cat("\n=== Validation of Extracted Data ===\n")

count_array <- counts(pixels)
dims <- dim(count_array)
cat(sprintf("Dimensions: %d x %d x %d x %d\n", dims[1], dims[2], dims[3], dims[4]))
cat(sprintf("  = %dx%d pixels x %d stripes x %d samples\n",
            dims[1], dims[2], dims[3], dims[4]))

# NA check
na_count <- sum(is.na(count_array))
na_percent <- 100 * na_count / length(count_array)
cat(sprintf("\nNA values: %d (%.2f%%)\n", na_count, na_percent))

if (na_percent > 5) {
  cat("  WARNING: High NA percentage - may indicate extraction issues\n")
}

# Value distribution
non_na_values <- count_array[!is.na(count_array)]
if (length(non_na_values) > 0) {
  cat("\nValue distribution (excluding NAs):\n")
  print(summary(non_na_values))

  zeros <- sum(non_na_values == 0)
  sparsity <- 100 * zeros / length(non_na_values)
  cat(sprintf("\nSparsity: %.1f%% zeros\n", sparsity))
}

# ==============================================================================
# AGGREGATE TO COUNT MATRIX
# ==============================================================================
cat("\n=== Aggregating to Count Matrix ===\n")
cat("Strategy: Sum all pixels in 5x5 region\n")

counts_matrix <- matrix(NA, nrow = dims[3], ncol = dims[4])
colnames(counts_matrix) <- names(hicFiles)
rownames(counts_matrix) <- stripes$stripe_id

cat("Aggregating...")
for (i in 1:dims[3]) {
  if (i %% 50 == 0) cat(sprintf(" %d", i))
  for (j in 1:dims[4]) {
    mat_5x5 <- as.matrix(count_array[, , i, j])
    counts_matrix[i, j] <- sum(mat_5x5, na.rm = TRUE)
  }
}
cat(" Done\n")

cat(sprintf("  Count matrix: %d stripes x %d samples\n",
            nrow(counts_matrix), ncol(counts_matrix)))

# Summary statistics
cat("\nCount summary:\n")
print(summary(as.vector(counts_matrix)))

# ==============================================================================
# CORRELATION ANALYSIS
# ==============================================================================
cat("\n=== Correlation Analysis ===\n")

n_samples <- dims[4]
cor_matrix <- matrix(NA, nrow = n_samples, ncol = n_samples)
rownames(cor_matrix) <- names(hicFiles)
colnames(cor_matrix) <- names(hicFiles)

for (i in 1:n_samples) {
  for (j in 1:n_samples) {
    cor_matrix[i, j] <- cor(counts_matrix[, i], counts_matrix[, j],
                            use = "complete.obs")
  }
}

cat("\nCorrelation matrix:\n")
print(round(cor_matrix, 3))

# Within-condition correlations
cat("\nWithin-condition correlations:\n")
ctrl_indices <- 1:3
mut_indices <- 4:6

ctrl_cors <- cor_matrix[ctrl_indices, ctrl_indices]
ctrl_within_mean <- mean(ctrl_cors[upper.tri(ctrl_cors)])
cat(sprintf("  Control (within): %.3f\n", ctrl_within_mean))

mut_cors <- cor_matrix[mut_indices, mut_indices]
mut_within_mean <- mean(mut_cors[upper.tri(mut_cors)])
cat(sprintf("  Mutant (within): %.3f\n", mut_within_mean))

between_cors <- as.vector(cor_matrix[ctrl_indices, mut_indices])
between_mean <- mean(between_cors)
cat(sprintf("  Between groups: %.3f\n", between_mean))

# Quality assessment
cat("\nQuality Assessment:\n")
if (ctrl_within_mean > 0.95 && mut_within_mean > 0.95) {
  cat("  Excellent within-condition reproducibility (>0.95)\n")
} else if (ctrl_within_mean > 0.90 && mut_within_mean > 0.90) {
  cat("  Good within-condition reproducibility (>0.90)\n")
} else {
  cat("  WARNING: Lower than expected within-condition correlation (<0.90)\n")
}

if ((ctrl_within_mean - between_mean) > 0.03 ||
    (mut_within_mean - between_mean) > 0.03) {
  cat("  Clear biological signal (within > between correlations)\n")
} else {
  cat("  WARNING: Within-condition correlations not much higher than between\n")
}

# ==============================================================================
# SAVE OUTPUTS
# ==============================================================================
cat("\n=== Saving Outputs ===\n")

output_dir <- input_dir  # Same as input directory

# Save count matrix
counts_tsv <- file.path(output_dir, "02_stripe_counts.tsv")
counts_rds <- file.path(output_dir, "02_stripe_counts.rds")

write.table(counts_matrix, counts_tsv, sep = "\t", quote = FALSE,
            row.names = TRUE, col.names = NA)
saveRDS(counts_matrix, counts_rds)

cat(sprintf("  Count matrix: %s\n", counts_tsv))
cat(sprintf("  Count matrix RDS: %s\n", counts_rds))

# Save extraction metadata
metadata <- list(
  timepoint = TIMEPOINT,
  timepoint_dir = timepoint_dir,
  resolution = resolution,
  sampling_distance = sampling_distance,
  buffer_bins = buffer_bins,
  normalization = norm_to_use,
  n_stripes = nrow(stripes),
  n_samples = length(hicFiles),
  sample_names = names(hicFiles),
  hic_files = hicFiles,
  dims = dims,
  extraction_time = as.numeric(extraction_time),
  na_percent = na_percent,
  correlation = list(
    matrix = cor_matrix,
    within_ctrl = ctrl_within_mean,
    within_mut = mut_within_mean,
    between = between_mean
  ),
  date = Sys.Date()
)

metadata_file <- file.path(output_dir, "02_extraction_metadata.rds")
saveRDS(metadata, metadata_file)
cat(sprintf("  Metadata: %s\n", metadata_file))

# Save extracted pixels (HDF5)
cat("\nSaving HDF5 data...\n")
saveHDF5SummarizedExperiment(
  x = pixels,
  dir = output_dir,
  prefix = "02_extracted",
  replace = TRUE,
  verbose = TRUE
)

# ==============================================================================
# SUMMARY
# ==============================================================================
cat("\n=================================\n")
cat("Phase 2 Complete\n")
cat("=================================\n\n")

cat(sprintf("Timepoint: %s (%s)\n", TIMEPOINT, timepoint_dir))
cat(sprintf("Resolution: %d bp\n", resolution))
cat(sprintf("Stripes quantified: %d\n", nrow(stripes)))
cat(sprintf("Samples: %d\n", length(hicFiles)))
cat(sprintf("Normalization: %s\n", norm_to_use))

cat("\nQuality Metrics:\n")
cat(sprintf("  Within-ctrl correlation: %.3f\n", ctrl_within_mean))
cat(sprintf("  Within-mut correlation: %.3f\n", mut_within_mean))
cat(sprintf("  Between-group correlation: %.3f\n", between_mean))
cat(sprintf("  NA percent: %.2f%%\n", na_percent))
cat(sprintf("  Extraction time: %.1f sec\n", extraction_time))

cat(sprintf("\nOutput: %s\n", output_dir))
cat(sprintf("\nNext: Rscript scripts/phase3_edgeR.R %s\n\n", TIMEPOINT))
