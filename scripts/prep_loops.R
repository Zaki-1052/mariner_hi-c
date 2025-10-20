#!/usr/bin/env Rscript
# /expanse/lustre/projects/csd940/zalibhai/mariner/scripts/prepare_loops.R
# Modified for biological replicates (n=3 per condition)
# Multi-resolution support: accepts resolution as command-line argument

library(mariner)
library(InteractionSet)
library(GenomicRanges)

# Parse command-line arguments for resolution
args <- commandArgs(trailingOnly = TRUE)
RESOLUTION <- if (length(args) > 0) as.numeric(args[1]) else 5000

cat("\n========================================\n")
cat(sprintf("RESOLUTION: %d bp (%d kb)\n", RESOLUTION, RESOLUTION/1000))
cat("========================================\n\n")

# Define paths to individual replicate BEDPE files (resolution-aware)
base_path <- "/expanse/lustre/projects/csd940/ctea/nf-hic/juicer_frompre/hiccups_results"
bedpeFiles <- c(
  ctrl_M1 = sprintf("%s/ctrl_M1/postprocessed_pixels_%d.bedpe", base_path, RESOLUTION),
  ctrl_M2 = sprintf("%s/ctrl_M2/postprocessed_pixels_%d.bedpe", base_path, RESOLUTION),
  ctrl_M3 = sprintf("%s/ctrl_M3/postprocessed_pixels_%d.bedpe", base_path, RESOLUTION),
  mut_M1 = sprintf("%s/mut_M1/postprocessed_pixels_%d.bedpe", base_path, RESOLUTION),
  mut_M2 = sprintf("%s/mut_M2/postprocessed_pixels_%d.bedpe", base_path, RESOLUTION),
  mut_M3 = sprintf("%s/mut_M3/postprocessed_pixels_%d.bedpe", base_path, RESOLUTION)
)

bedpe_colnames <- c(
		      "chr1", "x1", "x2", "chr2", "y1", "y2",
		        "name", "score", "strand1", "strand2", "color",
		        "observed", "expectedBL", "expectedDonut", "expectedH", "expectedV",
			  "fdrBL", "fdrDonut", "fdrH", "fdrV",
			  "numCollapsed", "centroid1", "centroid2", "radius"
			  )

read_bedpe <- function(filepath, sample_name, expected_resolution) {
  cat(sprintf("Reading %s... ", sample_name))

  # Check file exists
  if (!file.exists(filepath)) {
    stop(sprintf("ERROR: File not found: %s", filepath))
  }

  bedpe <- read.table(filepath, header = FALSE, skip = 2,
                      stringsAsFactors = FALSE, sep = "\t",
                      comment.char = "")
  colnames(bedpe) <- bedpe_colnames

  # Filter for specified resolution
  bedpe$resolution <- bedpe$x2 - bedpe$x1
  bedpe_filtered <- bedpe[bedpe$resolution == expected_resolution, ]

  # Remove zero observed counts
  bedpe_subset <- bedpe_filtered[bedpe_filtered$observed > 0, ]

  cat(sprintf("%d loops (%dkb, observed>0)\n", nrow(bedpe_subset), expected_resolution/1000))
  return(bedpe_subset)
}

# Read all 6 BEDPE files
cat("\n=== Reading BEDPE files for all replicates ===\n")
bedpe_list <- lapply(names(bedpeFiles), function(name) {
  read_bedpe(bedpeFiles[name], name, RESOLUTION)
})
names(bedpe_list) <- names(bedpeFiles)

# Summary statistics
cat("\nPer-replicate loop counts:\n")
for (name in names(bedpe_list)) {
  cat(sprintf("  %s: %d loops\n", name, nrow(bedpe_list[[name]])))
}
cat(sprintf("\nTotal loops across all replicates: %d\n",
            sum(sapply(bedpe_list, nrow))))

# Convert all to GInteractions
cat("\n=== Converting to GInteractions ===\n")
gi_list <- lapply(names(bedpe_list), function(name) {
  cat(sprintf("Converting %s... ", name))
  gi <- as_ginteractions(bedpe_list[[name]],
                         keep.extra.columns = TRUE,
                         starts.in.df.are.0based = TRUE)
  cat(sprintf("%d interactions\n", length(gi)))
  return(gi)
})
names(gi_list) <- names(bedpe_list)

# Create resolution-specific output directory
output_dir <- sprintf("outputs/res_%dkb", RESOLUTION/1000)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat(sprintf("Created output directory: %s\n", output_dir))
}

# Save individual GInteractions
saveRDS(gi_list, file.path(output_dir, "01_ginteractions.rds"))
cat(sprintf("✓ Saved to: %s/01_ginteractions.rds\n", output_dir))

# Merge loops across all 6 replicates (union approach)
cat("\n=== Merging loops across all replicates ===\n")
cat("Strategy: Union (include loops from ANY replicate)\n")
cat("Radius: 10 kb\n")
cat("Selection: Maximum observed count\n\n")

merged <- mergePairs(
  x = gi_list,
  radius = 10e3,
  column = "observed",
  selectMax = TRUE,
  method = "manhattan"
)

# Detailed clustering statistics
cluster_sizes <- sapply(clusters(merged), nrow)
total_input <- sum(sapply(gi_list, length))

cat(sprintf("✓ Merged %d consensus loops from %d input loops\n",
            length(merged), total_input))
cat(sprintf("  Reduction: %.1f%%\n\n",
            100 * (1 - length(merged)/total_input)))

# Cluster size distribution
cat("Cluster size distribution:\n")
size_table <- table(cluster_sizes)
for (size in sort(unique(cluster_sizes))) {
  count <- sum(cluster_sizes == size)
  pct <- 100 * count / length(cluster_sizes)
  cat(sprintf("  %d replicate%s: %d loops (%.1f%%)\n",
              size, ifelse(size > 1, "s", ""), count, pct))
}

# Per-replicate contribution to merged set
cat("\nPer-replicate contribution to consensus:\n")
for (name in names(gi_list)) {
  # Count how many merged loops have this replicate in their cluster
  has_replicate <- sapply(clusters(merged), function(cluster) {
    name %in% cluster$source
  })
  count <- sum(has_replicate)
  pct <- 100 * count / length(merged)
  cat(sprintf("  %s: %d/%d loops (%.1f%%)\n",
              name, count, length(merged), pct))
}


cat("\n")
saveRDS(merged, file.path(output_dir, "02_merged.rds"))
cat(sprintf("✓ Saved to: %s/02_merged.rds\n", output_dir))

# Binning to resolution-specific grid
cat(sprintf("\n=== Binning to %dkb grid ===\n", RESOLUTION/1000))
binned <- assignToBins(
  x = merged,
  binSize = RESOLUTION,
  pos1 = "center",
  pos2 = "center"
)

cat(sprintf("✓ Binned %d loops to %dkb resolution\n", length(binned), RESOLUTION/1000))
saveRDS(binned, file.path(output_dir, "03_binned.rds"))
cat(sprintf("✓ Saved to: %s/03_binned.rds\n", output_dir))

# Create buffered regions (5x5 pixels)
cat("\n=== Creating buffered regions ===\n")
buffer_kb <- 2 * RESOLUTION / 1000
cat(sprintf("Buffer: ±2 bins (±%dkb) around center\n", buffer_kb))
cat("Purpose: Handle positional shifts between replicates\n\n")

buffered <- pixelsToMatrices(
  x = binned,
  buffer = 2
)

anchor1_width <- unique(width(anchors(buffered, type = "first")))
anchor2_width <- unique(width(anchors(buffered, type = "second")))

buffer_size_kb <- (anchor1_width[1] - 1) / 1000

cat(sprintf("✓ Created %d buffered regions\n", length(buffered)))
cat(sprintf("  Region size: %.0f × %.0f pixels (%dkb × %dkb)\n",
            (anchor1_width-1)/RESOLUTION, (anchor2_width-1)/RESOLUTION,
            buffer_size_kb, buffer_size_kb))

saveRDS(buffered, file.path(output_dir, "04_buffered.rds"))
cat(sprintf("✓ Saved to: %s/04_buffered.rds\n", output_dir))

# Final summary
cat("\n=================================\n")
cat("LOOP PREPARATION COMPLETE\n")
cat("=================================\n")
cat(sprintf("Resolution: %d bp (%d kb)\n", RESOLUTION, RESOLUTION/1000))
cat(sprintf("Input: %d loops across 6 replicates\n", total_input))
cat(sprintf("Output: %d consensus loop positions\n", length(buffered)))
cat(sprintf("Ready for extraction: %d loops × 6 replicates × 25 pixels\n",
            length(buffered)))
cat(sprintf("                    = %d total extractions\n\n",
            length(buffered) * 6 * 25))
cat(sprintf("Output directory: %s\n", output_dir))
cat(sprintf("Next step: Rscript scripts/extract_counts.R %d\n", RESOLUTION))
