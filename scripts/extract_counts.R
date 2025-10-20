#!/usr/bin/env Rscript
# /expanse/lustre/projects/csd940/zalibhai/mariner/scripts/05_extract_counts.R
# Purpose: Extract Hi-C contact matrices at buffered loop positions
# Input: 04_buffered.rds (GInteractions with 5x5 pixel regions)
# Output: 05_extracted.h5 (HDF5-backed InteractionArray)

cat("\n=================================\n")
cat("Step 5: Hi-C Contact Extraction\n")
cat("=================================\n\n")

# Load required libraries
suppressPackageStartupMessages({
  library(mariner)
  library(InteractionSet)
  library(strawr)
  library(DelayedArray)
  library(HDF5Array)
})

# Set up paths
base_dir <- "/expanse/lustre/projects/csd940/zalibhai/mariner"
setwd(base_dir)

# Load buffered loops
cat("Loading buffered loops...\n")
buffered <- readRDS("outputs/full/04_buffered.rds")
cat(sprintf("  Loaded %d loops with 5x5 pixel regions\n", length(buffered)))

# Define .hic files
hicFiles <- c(
  ctrl = "/expanse/lustre/projects/csd940/ctea/nf-hic/juicerpre/merged/hic/resorted_ctrl.hic",
  mut = "/expanse/lustre/projects/csd940/ctea/nf-hic/juicerpre/merged/hic/resorted_mut.hic"
)

# Verify files exist
cat("\nVerifying .hic files...\n")
for (name in names(hicFiles)) {
  filepath <- hicFiles[name]
  if (!file.exists(filepath)) {
    stop(sprintf("ERROR: %s file not found at %s", name, filepath))
  }
  size_gb <- file.info(filepath)$size / 1e9
  cat(sprintf("  %s: %.1f GB - OK\n", name, size_gb))
}

# Check available resolutions and normalizations
cat("\nChecking .hic file properties...\n")
for (name in names(hicFiles)) {
  cat(sprintf("\n%s file:\n", name))
  resolutions <- readHicBpResolutions(hicFiles[name])
  cat("  Resolutions (bp):", paste(head(resolutions, 5), collapse=", "), "...\n")
  if (!5000 %in% resolutions) {
    stop(sprintf("ERROR: 5kb resolution not available in %s", name))
  }
  norms <- readHicNormTypes(hicFiles[name])
  cat("  Normalizations:", paste(norms, collapse=", "), "\n")
}

# Use VC normalization (KR not available)
norm_to_use <- "VC"
cat(sprintf("\nUsing normalization: %s\n", norm_to_use))

# Set up HDF5 file path
hdf5_dir <- file.path(base_dir, "temp_hdf5")
if (!dir.exists(hdf5_dir)) {
  dir.create(hdf5_dir, recursive = TRUE)
}
h5_file_path <- file.path(hdf5_dir, "extracted_matrices.h5")

# Pre-extraction validation
cat("\nValidating loop positions...\n")
cat(sprintf("  Chromosomes in loops: %d\n", length(seqlevels(buffered))))
cat(sprintf("  Example chromosomes: %s\n", 
            paste(head(seqlevels(buffered), 5), collapse=", ")))

anchor1_width <- max(width(anchors(buffered, type = "first")))
anchor2_width <- max(width(anchors(buffered, type = "second")))
min_blockSize <- 2 * max(anchor1_width, anchor2_width)
cat(sprintf("  Maximum anchor width: %d bp\n", max(anchor1_width, anchor2_width)))
cat(sprintf("  Minimum required blockSize: %d bp\n", min_blockSize))

# Main extraction
cat("\n=== Starting extraction ===\n")
cat(sprintf("Processing %d loops x %d files x 25 pixels = %d total extractions\n",
            length(buffered), length(hicFiles), length(buffered) * length(hicFiles) * 25))

start_time <- Sys.time()

pixels <- pullHicMatrices(
  x = buffered,
  files = hicFiles,
  binSize = 5e3,
  h5File = h5_file_path,
  norm = norm_to_use,
  matrix = "observed",
  blockSize = 1e6,
  onDisk = TRUE,
  compressionLevel = 1
)

end_time <- Sys.time()
extraction_time <- difftime(end_time, start_time, units = "secs")
cat(sprintf("\nExtraction completed in %.1f seconds\n", extraction_time))

# Validation
cat("\n=== Validation of Extracted Data ===\n")

# NOTE: Dimensions are width x height x loops x files (5 x 5 x 150 x 2)
dims <- dim(counts(pixels))
cat(sprintf("Dimensions: %d x %d x %d x %d\n", dims[1], dims[2], dims[3], dims[4]))
cat(sprintf("  = %dx%d pixels x %d loops x %d files\n", 
            dims[1], dims[2], dims[3], dims[4]))

# Extract counts array
count_array <- counts(pixels)

# NA check
na_count <- sum(is.na(count_array))
na_percent <- 100 * na_count / length(count_array)
cat(sprintf("\nNA values: %d (%.2f%%)\n", na_count, na_percent))

# Value distribution
non_na_values <- count_array[!is.na(count_array)]
if (length(non_na_values) > 0) {
  cat("\nValue distribution (excluding NAs):\n")
  print(summary(non_na_values))
  
  zeros <- sum(non_na_values == 0)
  sparsity <- 100 * zeros / length(non_na_values)
  cat(sprintf("\nSparsity: %.1f%% zeros\n", sparsity))
}

# Per-sample statistics (files are in dimension 4)
cat("\nPer-sample statistics:\n")
for (i in 1:dims[4]) {
  sample_values <- count_array[,,,i][!is.na(count_array[,,,i])]
  if (length(sample_values) > 0) {
    cat(sprintf("  %s: median=%.2f, mean=%.2f, max=%.0f\n",
                names(hicFiles)[i],
                median(sample_values),
                mean(sample_values),
                max(sample_values)))
  }
}

# Correlation between samples
cat("\nSample correlation:\n")
ctrl_values <- as.vector(count_array[,,,1])
mut_values <- as.vector(count_array[,,,2])
complete_pairs <- complete.cases(ctrl_values, mut_values)
if (sum(complete_pairs) > 0) {
  correlation <- cor(ctrl_values[complete_pairs], mut_values[complete_pairs])
  cat(sprintf("  Correlation (ctrl vs mut): %.3f\n", correlation))
}

# Show example matrices (FIXED)
cat("\n=== Example Extracted Matrices ===\n")
cat("Showing first 3 loops for control sample:\n\n")

for (i in 1:min(3, dims[3])) {
  cat(sprintf("Loop %d:\n", i))
  
  # Extract and convert to regular matrix
  example_matrix <- count_array[, , i, 1]  # Loop i, control file
  example_matrix <- as.matrix(example_matrix)  # Convert DelayedMatrix to regular matrix
  
  # Display matrix
  for (row in 1:5) {
    row_values <- round(example_matrix[row, ], 2)
    formatted_row <- format(row_values, width = 8, justify = "right")
    cat("  ", paste(formatted_row, collapse = " "), "\n")
  }
  
  # Summary stats
  matrix_vals <- as.vector(example_matrix)
  matrix_vals <- matrix_vals[!is.na(matrix_vals)]
  if (length(matrix_vals) > 0) {
    cat(sprintf("  Sum: %.0f, Max: %.0f, Non-zero: %d/25\n\n", 
                sum(matrix_vals), max(matrix_vals), sum(matrix_vals > 0)))
  }
}

# Save using HDF5Array method
cat("\n=== Saving Extracted Data ===\n")
output_dir <- "outputs/full"
output_prefix <- "05_extracted"

# Save as HDF5
cat("Saving as HDF5 (required for out-of-memory objects)...\n")

# Use the correct saving function for HDF5-backed objects
saveHDF5SummarizedExperiment(
  x = pixels,
  dir = output_dir,
  prefix = output_prefix,
  replace = TRUE,
  verbose = TRUE
)

cat(sprintf("  Saved to: %s/\n", file.path(output_dir, output_prefix)))

# Also save metadata for easy loading
metadata <- list(
  n_loops = dims[3],
  n_files = dims[4],
  files = hicFiles,
  normalization = norm_to_use,
  extraction_time = extraction_time,
  dims = dims,
  correlation = correlation
)
saveRDS(metadata, file.path(output_dir, "05_metadata.rds"))
cat("  Metadata saved to: outputs/test/05_metadata.rds\n")

# Final summary
cat("\n=================================\n")
cat("✅ Step 5: Extraction Complete\n")
cat("=================================\n")
cat(sprintf("- Processed: %d loops\n", dims[3]))
cat(sprintf("- Samples: %s\n", paste(names(hicFiles), collapse = ", ")))
cat(sprintf("- Sample correlation: %.3f\n", correlation))
cat(sprintf("- Normalization: %s\n", norm_to_use))
cat(sprintf("- Time: %.1f seconds\n", extraction_time))
cat(sprintf("- Output: %s/\n", file.path(output_dir, output_prefix)))
cat(sprintf("- Ready for aggregation: %d matrices\n", dims[3]))
cat("\nNext: Run aggregation script (06_aggregate.R)\n")
cat("\nTo load in next script:\n")
cat("  pixels <- loadHDF5SummarizedExperiment(dir='outputs/test', prefix='05_extracted')\n")
