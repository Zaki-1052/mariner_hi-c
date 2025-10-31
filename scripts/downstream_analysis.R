#!/usr/bin/env Rscript
# scripts/downstream_analysis.R
# Downstream Analysis: Merged Loops & Anchor Characterization
# Author: Zakir Alibhai
# Date: 2025-10-30
#
# Purpose:
#   1. Merge final results from 5kb, 10kb, 25kb with overlap removal
#   2. Validate all loops files for volcano plots
#   3. Characterize loop anchors with genomic features and gene annotations
#   4. Export to multiple formats (TSV, RDS, BEDPE)
#
# Usage:
#   Rscript scripts/downstream_analysis.R [--gtf PATH] [--tolerance KB]
#
# Arguments:
#   --gtf PATH       Path to GTF file for gene annotations (optional)
#   --tolerance KB   Overlap tolerance in kb (default: 10)

# =============================================================================
# SETUP
# =============================================================================

cat("\n")
cat("========================================\n")
cat("Downstream Analysis: Merged Loops\n")
cat("========================================\n\n")

suppressPackageStartupMessages({
  library(mariner)
  library(InteractionSet)
  library(GenomicRanges)
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
})

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
gtf_file <- NULL
tolerance_kb <- 10

if (length(args) > 0) {
  for (i in seq_along(args)) {
    if (args[i] == "--gtf" && i < length(args)) {
      gtf_file <- args[i + 1]
    } else if (args[i] == "--tolerance" && i < length(args)) {
      tolerance_kb <- as.numeric(args[i + 1])
    }
  }
}

cat(sprintf("Configuration:\n"))
cat(sprintf("  Overlap tolerance: %d kb\n", tolerance_kb))
if (!is.null(gtf_file)) {
  cat(sprintf("  GTF file: %s\n", gtf_file))
} else {
  cat("  GTF file: Using TxDb.Mmusculus.UCSC.mm10.knownGene\n")
}
cat("\n")

# Set working directory (adjust if needed)
if (dir.exists("/expanse/lustre/projects/csd940/zalibhai/mariner")) {
  setwd("/expanse/lustre/projects/csd940/zalibhai/mariner")
} else if (basename(getwd()) == "mariner-final") {
  # Already in project directory
} else {
  warning("Working directory may not be correct. Current: ", getwd())
}

# Create output directory
output_dir <- "outputs/merged_loops"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "plots"), recursive = TRUE, showWarnings = FALSE)

cat(sprintf("Output directory: %s\n\n", output_dir))

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' Match loops across GInteractions with tolerance
#'
#' @param coords1 GInteractions object
#' @param coords2 GInteractions object
#' @param tolerance_bp Tolerance in base pairs
#' @return data.frame with index1, index2, distance
match_loops <- function(coords1, coords2, tolerance_bp) {
  # Convert to data.frame with coordinates
  df1 <- data.frame(
    chr1 = as.character(seqnames(anchors(coords1, "first"))),
    start1 = start(anchors(coords1, "first")),
    end1 = end(anchors(coords1, "first")),
    chr2 = as.character(seqnames(anchors(coords1, "second"))),
    start2 = start(anchors(coords1, "second")),
    end2 = end(anchors(coords1, "second")),
    index1 = 1:length(coords1)
  )

  df2 <- data.frame(
    chr1 = as.character(seqnames(anchors(coords2, "first"))),
    start1 = start(anchors(coords2, "first")),
    end1 = end(anchors(coords2, "first")),
    chr2 = as.character(seqnames(anchors(coords2, "second"))),
    start2 = start(anchors(coords2, "second")),
    end2 = end(anchors(coords2, "second")),
    index2 = 1:length(coords2)
  )

  # Find overlaps
  matches <- data.frame(index1 = integer(), index2 = integer(), distance = numeric())

  for (i in 1:nrow(df1)) {
    # Same chromosome check
    same_chr <- df2$chr1 == df1$chr1[i] & df2$chr2 == df1$chr2[i]

    # Distance check (within tolerance for both anchors)
    dist1 <- abs(df2$start1 - df1$start1[i])
    dist2 <- abs(df2$start2 - df1$start2[i])

    close_enough <- same_chr & dist1 <= tolerance_bp & dist2 <= tolerance_bp

    if (any(close_enough)) {
      # Take closest match
      distances <- dist1 + dist2
      best_match <- which(close_enough)[which.min(distances[close_enough])]
      matches <- rbind(matches, data.frame(
        index1 = i,
        index2 = best_match,
        distance = distances[best_match]
      ))
    }
  }

  return(matches)
}

# =============================================================================
# SECTION 1: LOAD FINAL RESULTS FROM ALL RESOLUTIONS
# =============================================================================

cat("="*70, "\n", sep="")
cat("SECTION 1: Loading Final Results\n")
cat("="*70, "\n\n", sep="")

resolutions <- c(5000, 10000, 25000)
final_results_list <- list()
coords_list <- list()

for (res in resolutions) {
  res_kb <- res / 1000
  cat(sprintf("Loading %dkb results...\n", res_kb))

  # Load final results TSV
  final_file <- sprintf("outputs/edgeR_results_res_%dkb/primary_analysis/final_results.tsv", res_kb)
  coords_file <- sprintf("outputs/res_%dkb/03_binned.rds", res_kb)

  if (!file.exists(final_file)) {
    warning(sprintf("  Final results file not found: %s", final_file))
    next
  }

  if (!file.exists(coords_file)) {
    warning(sprintf("  Coordinates file not found: %s", coords_file))
    next
  }

  # Load data
  final_df <- read.table(final_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  coords <- readRDS(coords_file)

  cat(sprintf("  Final results: %d loops\n", nrow(final_df)))
  cat(sprintf("  Coordinates: %d interactions\n", length(coords)))

  # Match final results to coordinates
  coords_df <- data.frame(
    chr1 = as.character(seqnames(anchors(coords, "first"))),
    start1 = start(anchors(coords, "first")),
    end1 = end(anchors(coords, "first")),
    chr2 = as.character(seqnames(anchors(coords, "second"))),
    start2 = start(anchors(coords, "second")),
    end2 = end(anchors(coords, "second"))
  )

  # Find matching indices
  matches <- numeric(nrow(final_df))
  for (i in 1:nrow(final_df)) {
    idx <- which(
      coords_df$chr1 == final_df$chr1[i] &
      coords_df$start1 == final_df$start1[i] &
      coords_df$end1 == final_df$end1[i] &
      coords_df$chr2 == final_df$chr2[i] &
      coords_df$start2 == final_df$start2[i] &
      coords_df$end2 == final_df$end2[i]
    )
    if (length(idx) > 0) {
      matches[i] <- idx[1]
    }
  }

  # Filter to matched loops
  matches <- matches[matches > 0]

  if (length(matches) == 0) {
    warning(sprintf("  No matches found between final results and coordinates"))
    next
  }

  # Store filtered coordinates and results
  coords_filtered <- coords[matches]
  final_df_filtered <- final_df[matches > 0, ]

  # Add resolution column
  final_df_filtered$resolution <- res
  final_df_filtered$resolution_kb <- res_kb

  final_results_list[[as.character(res)]] <- final_df_filtered
  coords_list[[as.character(res)]] <- coords_filtered

  cat(sprintf("  Matched: %d loops\n\n", length(matches)))
}

if (length(final_results_list) == 0) {
  stop("ERROR: No final results loaded. Run generate_final_results.py first.")
}

cat(sprintf("✓ Loaded final results from %d resolutions\n", length(final_results_list)))
cat(sprintf("  Total loops before merging: %d\n\n", sum(sapply(final_results_list, nrow))))

# =============================================================================
# SECTION 2: MERGE LOOPS WITH OVERLAP REMOVAL
# =============================================================================

cat("="*70, "\n", sep="")
cat("SECTION 2: Merging Loops (Overlap Removal)\n")
cat("="*70, "\n\n", sep="")

cat(sprintf("Overlap tolerance: %d kb\n", tolerance_kb))
cat("Priority: Lowest FDR (most significant)\n\n")

# Combine all loops into a single data structure
all_loops_df <- bind_rows(final_results_list)
all_loops_coords <- do.call(c, coords_list)

cat(sprintf("Total loops to process: %d\n", nrow(all_loops_df)))
cat(sprintf("  5kb: %d loops\n", sum(all_loops_df$resolution == 5000)))
cat(sprintf("  10kb: %d loops\n", sum(all_loops_df$resolution == 10000)))
cat(sprintf("  25kb: %d loops\n\n", sum(all_loops_df$resolution == 25000)))

# Build overlap graph
cat("Detecting overlaps...\n")

tolerance_bp <- tolerance_kb * 1000
n_loops <- length(all_loops_coords)

# Create adjacency matrix for overlaps
overlap_graph <- matrix(FALSE, nrow = n_loops, ncol = n_loops)

for (i in 1:(n_loops - 1)) {
  if (i %% 100 == 0) {
    cat(sprintf("  Processed %d/%d loops\r", i, n_loops))
  }

  # Get coordinates for loop i
  chr1_i <- as.character(seqnames(anchors(all_loops_coords[i], "first")))
  start1_i <- start(anchors(all_loops_coords[i], "first"))
  chr2_i <- as.character(seqnames(anchors(all_loops_coords[i], "second")))
  start2_i <- start(anchors(all_loops_coords[i], "second"))

  # Check all subsequent loops
  for (j in (i + 1):n_loops) {
    chr1_j <- as.character(seqnames(anchors(all_loops_coords[j], "first")))
    start1_j <- start(anchors(all_loops_coords[j], "first"))
    chr2_j <- as.character(seqnames(anchors(all_loops_coords[j], "second")))
    start2_j <- start(anchors(all_loops_coords[j], "second"))

    # Check if same chromosome
    if (chr1_i == chr1_j && chr2_i == chr2_j) {
      # Check distance within tolerance
      dist1 <- abs(start1_i - start1_j)
      dist2 <- abs(start2_i - start2_j)

      if (dist1 <= tolerance_bp && dist2 <= tolerance_bp) {
        overlap_graph[i, j] <- TRUE
        overlap_graph[j, i] <- TRUE
      }
    }
  }
}

cat(sprintf("  Processed %d/%d loops\n", n_loops, n_loops))

# Find connected components (clusters of overlapping loops)
cat("\nFinding overlap clusters...\n")

visited <- rep(FALSE, n_loops)
clusters <- list()
cluster_id <- 0

for (i in 1:n_loops) {
  if (!visited[i]) {
    # Start new cluster
    cluster_id <- cluster_id + 1
    current_cluster <- c(i)
    visited[i] <- TRUE

    # Find all connected loops (BFS)
    queue <- c(i)

    while (length(queue) > 0) {
      current <- queue[1]
      queue <- queue[-1]

      # Find neighbors
      neighbors <- which(overlap_graph[current, ] & !visited)

      if (length(neighbors) > 0) {
        visited[neighbors] <- TRUE
        current_cluster <- c(current_cluster, neighbors)
        queue <- c(queue, neighbors)
      }
    }

    clusters[[cluster_id]] <- current_cluster
  }
}

n_singleton <- sum(sapply(clusters, length) == 1)
n_merged <- sum(sapply(clusters, length) > 1)

cat(sprintf("  Found %d clusters\n", length(clusters)))
cat(sprintf("    Singleton clusters (no overlaps): %d\n", n_singleton))
cat(sprintf("    Multi-loop clusters (will merge): %d\n", n_merged))
cat(sprintf("    Total loops to merge: %d\n\n", sum(sapply(clusters[sapply(clusters, length) > 1], length))))

# Select representative loop from each cluster (lowest FDR)
cat("Selecting representative loops (lowest FDR)...\n")

keep_indices <- numeric(length(clusters))
overlap_report <- data.frame(
  cluster_id = integer(),
  n_overlaps = integer(),
  kept_index = integer(),
  kept_resolution = character(),
  kept_fdr = numeric(),
  merged_indices = character(),
  merged_resolutions = character(),
  stringsAsFactors = FALSE
)

for (i in seq_along(clusters)) {
  cluster_indices <- clusters[[i]]
  cluster_fdr <- all_loops_df$FDR[cluster_indices]

  # Select loop with lowest FDR
  best_idx_in_cluster <- which.min(cluster_fdr)
  best_idx_global <- cluster_indices[best_idx_in_cluster]
  keep_indices[i] <- best_idx_global

  # Record in overlap report
  if (length(cluster_indices) > 1) {
    overlap_report <- rbind(overlap_report, data.frame(
      cluster_id = i,
      n_overlaps = length(cluster_indices),
      kept_index = best_idx_global,
      kept_resolution = sprintf("%dkb", all_loops_df$resolution_kb[best_idx_global]),
      kept_fdr = all_loops_df$FDR[best_idx_global],
      merged_indices = paste(cluster_indices, collapse = ","),
      merged_resolutions = paste(sprintf("%dkb", all_loops_df$resolution_kb[cluster_indices]), collapse = ","),
      stringsAsFactors = FALSE
    ))
  }
}

cat(sprintf("  Selected %d representative loops\n", length(keep_indices)))
cat(sprintf("  Removed %d overlapping loops\n\n", n_loops - length(keep_indices)))

# Create merged dataset
merged_loops_df <- all_loops_df[keep_indices, ]
merged_loops_coords <- all_loops_coords[keep_indices]

# Add overlap annotation
merged_loops_df$n_overlaps <- sapply(clusters, length)
merged_loops_df$is_merged <- merged_loops_df$n_overlaps > 1

cat(sprintf("✓ Non-redundant loop set: %d loops\n", nrow(merged_loops_df)))
cat(sprintf("  From 5kb: %d loops\n", sum(merged_loops_df$resolution == 5000)))
cat(sprintf("  From 10kb: %d loops\n", sum(merged_loops_df$resolution == 10000)))
cat(sprintf("  From 25kb: %d loops\n\n", sum(merged_loops_df$resolution == 25000)))

# Save overlap report
if (nrow(overlap_report) > 0) {
  write.table(overlap_report,
              file.path(output_dir, "overlap_report.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("✓ Overlap report saved: %s\n\n", file.path(output_dir, "overlap_report.tsv")))
}

# =============================================================================
# SECTION 3: VALIDATE ALL LOOPS FILES FOR VOLCANO PLOTS
# =============================================================================

cat("="*70, "\n", sep="")
cat("SECTION 3: Validating All Loops Files\n")
cat("="*70, "\n\n", sep="")

volcano_status <- data.frame(
  resolution = character(),
  file_exists = logical(),
  n_loops = integer(),
  has_logFC = logical(),
  has_FDR = logical(),
  has_PValue = logical(),
  stringsAsFactors = FALSE
)

for (res_kb in c(5, 10, 25)) {
  all_results_file <- sprintf("outputs/edgeR_results_res_%dkb/primary_analysis/all_results_primary.tsv", res_kb)

  file_exists <- file.exists(all_results_file)

  if (file_exists) {
    df <- read.table(all_results_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE, nrows = 1)

    volcano_status <- rbind(volcano_status, data.frame(
      resolution = sprintf("%dkb", res_kb),
      file_exists = TRUE,
      n_loops = NA,  # Will count later
      has_logFC = "logFC" %in% colnames(df),
      has_FDR = "FDR" %in% colnames(df),
      has_PValue = "PValue" %in% colnames(df),
      stringsAsFactors = FALSE
    ))

    # Count actual loops
    df_full <- read.table(all_results_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    volcano_status$n_loops[nrow(volcano_status)] <- nrow(df_full)

    cat(sprintf("%s: ✓ Valid (%d loops)\n", sprintf("%dkb", res_kb), nrow(df_full)))
  } else {
    volcano_status <- rbind(volcano_status, data.frame(
      resolution = sprintf("%dkb", res_kb),
      file_exists = FALSE,
      n_loops = NA,
      has_logFC = FALSE,
      has_FDR = FALSE,
      has_PValue = FALSE,
      stringsAsFactors = FALSE
    ))

    cat(sprintf("%s: ✗ File not found\n", sprintf("%dkb", res_kb)))
  }
}

cat("\n")

# Save validation report
write.table(volcano_status,
            file.path(output_dir, "volcano_files_validation.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("✓ Validation complete\n"))
cat(sprintf("  Files saved in: outputs/edgeR_results_res_{RES}kb/primary_analysis/all_results_primary.tsv\n\n"))

# =============================================================================
# SECTION 4: LOOP ANCHOR CHARACTERIZATIONS
# =============================================================================

cat("="*70, "\n", sep="")
cat("SECTION 4: Loop Anchor Characterizations\n")
cat("="*70, "\n\n", sep="")

# --- 4A. Basic Genomic Features ---

cat("4A. Computing basic genomic features...\n")

# Extract anchor coordinates
anchor1_gr <- anchors(merged_loops_coords, "first")
anchor2_gr <- anchors(merged_loops_coords, "second")

# Calculate midpoints
anchor1_mid <- start(anchor1_gr) + (end(anchor1_gr) - start(anchor1_gr)) / 2
anchor2_mid <- start(anchor2_gr) + (end(anchor2_gr) - start(anchor2_gr)) / 2

# Basic features
merged_loops_df$anchor1_chr <- as.character(seqnames(anchor1_gr))
merged_loops_df$anchor1_start <- start(anchor1_gr)
merged_loops_df$anchor1_end <- end(anchor1_gr)
merged_loops_df$anchor1_midpoint <- anchor1_mid
merged_loops_df$anchor1_size <- end(anchor1_gr) - start(anchor1_gr)

merged_loops_df$anchor2_chr <- as.character(seqnames(anchor2_gr))
merged_loops_df$anchor2_start <- start(anchor2_gr)
merged_loops_df$anchor2_end <- end(anchor2_gr)
merged_loops_df$anchor2_midpoint <- anchor2_mid
merged_loops_df$anchor2_size <- end(anchor2_gr) - start(anchor2_gr)

# Loop distance
merged_loops_df$loop_distance <- anchor2_mid - anchor1_mid

# Distance categories
merged_loops_df$distance_category <- cut(
  merged_loops_df$loop_distance,
  breaks = c(0, 100000, 500000, 1000000, Inf),
  labels = c("<100kb", "100-500kb", "500kb-1Mb", ">1Mb"),
  include.lowest = TRUE
)

# Chromosome type
merged_loops_df$is_intrachromosomal <- merged_loops_df$anchor1_chr == merged_loops_df$anchor2_chr

cat(sprintf("  Computed features for %d loops\n", nrow(merged_loops_df)))
cat(sprintf("  Distance range: %d bp - %d bp\n",
            min(merged_loops_df$loop_distance),
            max(merged_loops_df$loop_distance)))
cat(sprintf("  All loops are intrachromosomal: %s\n\n",
            ifelse(all(merged_loops_df$is_intrachromosomal), "Yes", "No")))

# --- 4B. Gene Annotations ---

cat("4B. Annotating with gene information...\n")

# Load gene annotations
if (!is.null(gtf_file) && file.exists(gtf_file)) {
  cat(sprintf("  Loading GTF: %s\n", gtf_file))
  suppressPackageStartupMessages({
    library(GenomicFeatures)
  })
  txdb <- makeTxDbFromGFF(gtf_file, format = "gtf")
} else {
  cat("  Loading TxDb.Mmusculus.UCSC.mm10.knownGene\n")
  suppressPackageStartupMessages({
    library(TxDb.Mmusculus.UCSC.mm10.knownGene)
    library(org.Mm.eg.db)
  })
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
}

# Extract genes
genes <- genes(txdb)

cat(sprintf("  Loaded %d genes\n", length(genes)))

# Get TSS (transcription start sites)
tss <- resize(genes, width = 1, fix = "start")

# Find nearest gene for each anchor
cat("  Finding nearest genes for anchor1...\n")
nearest1 <- distanceToNearest(anchor1_gr, tss)
merged_loops_df$anchor1_nearest_gene_idx <- NA
merged_loops_df$anchor1_distance_to_tss <- NA

if (length(nearest1) > 0) {
  merged_loops_df$anchor1_nearest_gene_idx[queryHits(nearest1)] <- subjectHits(nearest1)
  merged_loops_df$anchor1_distance_to_tss[queryHits(nearest1)] <- mcols(nearest1)$distance
}

cat("  Finding nearest genes for anchor2...\n")
nearest2 <- distanceToNearest(anchor2_gr, tss)
merged_loops_df$anchor2_nearest_gene_idx <- NA
merged_loops_df$anchor2_distance_to_tss <- NA

if (length(nearest2) > 0) {
  merged_loops_df$anchor2_nearest_gene_idx[queryHits(nearest2)] <- subjectHits(nearest2)
  merged_loops_df$anchor2_distance_to_tss[queryHits(nearest2)] <- mcols(nearest2)$distance
}

# Add gene IDs
merged_loops_df$anchor1_nearest_gene <- NA
merged_loops_df$anchor2_nearest_gene <- NA

for (i in 1:nrow(merged_loops_df)) {
  if (!is.na(merged_loops_df$anchor1_nearest_gene_idx[i])) {
    gene_idx <- merged_loops_df$anchor1_nearest_gene_idx[i]
    merged_loops_df$anchor1_nearest_gene[i] <- names(genes)[gene_idx]
  }

  if (!is.na(merged_loops_df$anchor2_nearest_gene_idx[i])) {
    gene_idx <- merged_loops_df$anchor2_nearest_gene_idx[i]
    merged_loops_df$anchor2_nearest_gene[i] <- names(genes)[gene_idx]
  }
}

# Classify anchors by proximity to TSS
# Promoter: within 2kb of TSS
merged_loops_df$anchor1_is_promoter <- !is.na(merged_loops_df$anchor1_distance_to_tss) &
                                        merged_loops_df$anchor1_distance_to_tss <= 2000
merged_loops_df$anchor2_is_promoter <- !is.na(merged_loops_df$anchor2_distance_to_tss) &
                                        merged_loops_df$anchor2_distance_to_tss <= 2000

# Classify loop type
merged_loops_df$loop_type <- case_when(
  merged_loops_df$anchor1_is_promoter & merged_loops_df$anchor2_is_promoter ~ "Promoter-Promoter",
  merged_loops_df$anchor1_is_promoter | merged_loops_df$anchor2_is_promoter ~ "Promoter-Enhancer",
  TRUE ~ "Enhancer-Enhancer"
)

cat(sprintf("\n  Gene annotation complete\n"))
cat(sprintf("  Loops with genes near both anchors: %d (%.1f%%)\n",
            sum(!is.na(merged_loops_df$anchor1_nearest_gene) & !is.na(merged_loops_df$anchor2_nearest_gene)),
            100 * mean(!is.na(merged_loops_df$anchor1_nearest_gene) & !is.na(merged_loops_df$anchor2_nearest_gene))))
cat(sprintf("\n  Loop type classification:\n"))
cat(sprintf("    Promoter-Promoter: %d (%.1f%%)\n",
            sum(merged_loops_df$loop_type == "Promoter-Promoter"),
            100 * mean(merged_loops_df$loop_type == "Promoter-Promoter")))
cat(sprintf("    Promoter-Enhancer: %d (%.1f%%)\n",
            sum(merged_loops_df$loop_type == "Promoter-Enhancer"),
            100 * mean(merged_loops_df$loop_type == "Promoter-Enhancer")))
cat(sprintf("    Enhancer-Enhancer: %d (%.1f%%)\n\n",
            sum(merged_loops_df$loop_type == "Enhancer-Enhancer"),
            100 * mean(merged_loops_df$loop_type == "Enhancer-Enhancer")))

# --- 4C. Characterization Summaries ---

cat("4C. Generating characterization summaries and plots...\n")

# Summary statistics
summary_text <- c(
  "="*70,
  "LOOP ANCHOR CHARACTERIZATION SUMMARY",
  "="*70,
  "",
  sprintf("Total non-redundant loops: %d", nrow(merged_loops_df)),
  sprintf("  From 5kb resolution: %d", sum(merged_loops_df$resolution == 5000)),
  sprintf("  From 10kb resolution: %d", sum(merged_loops_df$resolution == 10000)),
  sprintf("  From 25kb resolution: %d", sum(merged_loops_df$resolution == 25000)),
  "",
  "Loop Distance Distribution:",
  sprintf("  <100kb: %d (%.1f%%)",
          sum(merged_loops_df$distance_category == "<100kb", na.rm = TRUE),
          100 * mean(merged_loops_df$distance_category == "<100kb", na.rm = TRUE)),
  sprintf("  100-500kb: %d (%.1f%%)",
          sum(merged_loops_df$distance_category == "100-500kb", na.rm = TRUE),
          100 * mean(merged_loops_df$distance_category == "100-500kb", na.rm = TRUE)),
  sprintf("  500kb-1Mb: %d (%.1f%%)",
          sum(merged_loops_df$distance_category == "500kb-1Mb", na.rm = TRUE),
          100 * mean(merged_loops_df$distance_category == "500kb-1Mb", na.rm = TRUE)),
  sprintf("  >1Mb: %d (%.1f%%)",
          sum(merged_loops_df$distance_category == ">1Mb", na.rm = TRUE),
          100 * mean(merged_loops_df$distance_category == ">1Mb", na.rm = TRUE)),
  "",
  "Gene Association:",
  sprintf("  Anchor1 with nearby gene (<50kb): %d (%.1f%%)",
          sum(!is.na(merged_loops_df$anchor1_distance_to_tss) &
              merged_loops_df$anchor1_distance_to_tss < 50000),
          100 * mean(!is.na(merged_loops_df$anchor1_distance_to_tss) &
                     merged_loops_df$anchor1_distance_to_tss < 50000)),
  sprintf("  Anchor2 with nearby gene (<50kb): %d (%.1f%%)",
          sum(!is.na(merged_loops_df$anchor2_distance_to_tss) &
              merged_loops_df$anchor2_distance_to_tss < 50000),
          100 * mean(!is.na(merged_loops_df$anchor2_distance_to_tss) &
                     merged_loops_df$anchor2_distance_to_tss < 50000)),
  "",
  "Loop Type Classification:",
  sprintf("  Promoter-Promoter: %d (%.1f%%)",
          sum(merged_loops_df$loop_type == "Promoter-Promoter"),
          100 * mean(merged_loops_df$loop_type == "Promoter-Promoter")),
  sprintf("  Promoter-Enhancer: %d (%.1f%%)",
          sum(merged_loops_df$loop_type == "Promoter-Enhancer"),
          100 * mean(merged_loops_df$loop_type == "Promoter-Enhancer")),
  sprintf("  Enhancer-Enhancer: %d (%.1f%%)",
          sum(merged_loops_df$loop_type == "Enhancer-Enhancer"),
          100 * mean(merged_loops_df$loop_type == "Enhancer-Enhancer")),
  "",
  "Differential Loops:",
  sprintf("  Up in mutant: %d (%.1f%%)",
          sum(merged_loops_df$logFC > 0),
          100 * mean(merged_loops_df$logFC > 0)),
  sprintf("  Down in mutant: %d (%.1f%%)",
          sum(merged_loops_df$logFC < 0),
          100 * mean(merged_loops_df$logFC < 0)),
  "",
  "="*70
)

writeLines(summary_text, file.path(output_dir, "characterization_summary.txt"))

cat("  ✓ Summary saved\n")

# Generate plots

# Plot 1: Distance distribution
p1 <- ggplot(merged_loops_df, aes(x = loop_distance / 1e6)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "black", alpha = 0.7) +
  scale_x_continuous(trans = "log10") +
  labs(
    title = "Distribution of Loop Distances",
    x = "Loop Distance (Mb, log10)",
    y = "Count"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(output_dir, "plots", "distance_distribution.pdf"), p1, width = 8, height = 6)
cat("  ✓ Distance distribution plot saved\n")

# Plot 2: Chromosome distribution
chr_counts <- merged_loops_df %>%
  group_by(anchor1_chr) %>%
  summarise(n = n()) %>%
  arrange(desc(n))

p2 <- ggplot(chr_counts, aes(x = reorder(anchor1_chr, -n), y = n)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "black", alpha = 0.7) +
  labs(
    title = "Loop Distribution by Chromosome",
    x = "Chromosome",
    y = "Number of Loops"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(file.path(output_dir, "plots", "chromosome_distribution.pdf"), p2, width = 10, height = 6)
cat("  ✓ Chromosome distribution plot saved\n")

# Plot 3: Gene proximity distribution
gene_proximity_df <- data.frame(
  anchor = rep(c("Anchor1", "Anchor2"), each = nrow(merged_loops_df)),
  distance_to_tss = c(merged_loops_df$anchor1_distance_to_tss,
                      merged_loops_df$anchor2_distance_to_tss)
) %>%
  filter(!is.na(distance_to_tss))

p3 <- ggplot(gene_proximity_df, aes(x = distance_to_tss / 1000, fill = anchor)) +
  geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
  scale_x_log10() +
  labs(
    title = "Distance to Nearest Gene TSS",
    x = "Distance to TSS (kb, log10)",
    y = "Count",
    fill = "Anchor"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(output_dir, "plots", "gene_proximity_distribution.pdf"), p3, width = 8, height = 6)
cat("  ✓ Gene proximity distribution plot saved\n")

# Plot 4: Loop type classification
loop_type_counts <- merged_loops_df %>%
  group_by(loop_type, direction) %>%
  summarise(n = n(), .groups = "drop")

p4 <- ggplot(loop_type_counts, aes(x = loop_type, y = n, fill = direction)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", alpha = 0.7) +
  scale_fill_manual(values = c("up_in_mutant" = "#d73027", "down_in_mutant" = "#4575b4")) +
  labs(
    title = "Loop Type Classification by Direction",
    x = "Loop Type",
    y = "Number of Loops",
    fill = "Direction"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(file.path(output_dir, "plots", "loop_type_classification.pdf"), p4, width = 8, height = 6)
cat("  ✓ Loop type classification plot saved\n\n")

# =============================================================================
# SECTION 5: EXPORT TO MULTIPLE FORMATS
# =============================================================================

cat("="*70, "\n", sep="")
cat("SECTION 5: Exporting Results\n")
cat("="*70, "\n\n", sep="")

# --- 5A. TSV Export ---

cat("5A. Exporting TSV files...\n")

# Main results file with all features
write.table(merged_loops_df,
            file.path(output_dir, "characterized_loops.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  ✓ %s (%d loops)\n",
            "characterized_loops.tsv", nrow(merged_loops_df)))

# Non-redundant loops (basic info)
non_redundant_df <- merged_loops_df %>%
  select(loop_id, chr1, start1, end1, chr2, start2, end2,
         logFC, logCPM, FDR, PValue, significant,
         resolution_kb, n_overlaps, loop_distance, loop_type,
         direction, category)

write.table(non_redundant_df,
            file.path(output_dir, "non_redundant_loops.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  ✓ %s (%d loops)\n",
            "non_redundant_loops.tsv", nrow(non_redundant_df)))

# Anchor annotations (per-anchor detail)
anchor1_df <- merged_loops_df %>%
  select(loop_id, anchor = anchor1_chr,
         start = anchor1_start, end = anchor1_end, midpoint = anchor1_midpoint,
         nearest_gene = anchor1_nearest_gene,
         distance_to_tss = anchor1_distance_to_tss,
         is_promoter = anchor1_is_promoter) %>%
  mutate(anchor_num = 1)

anchor2_df <- merged_loops_df %>%
  select(loop_id, anchor = anchor2_chr,
         start = anchor2_start, end = anchor2_end, midpoint = anchor2_midpoint,
         nearest_gene = anchor2_nearest_gene,
         distance_to_tss = anchor2_distance_to_tss,
         is_promoter = anchor2_is_promoter) %>%
  mutate(anchor_num = 2)

anchor_annotations <- bind_rows(anchor1_df, anchor2_df)

write.table(anchor_annotations,
            file.path(output_dir, "anchor_annotations.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  ✓ %s (%d anchors)\n\n",
            "anchor_annotations.tsv", nrow(anchor_annotations)))

# --- 5B. RDS Export ---

cat("5B. Exporting RDS files...\n")

# Add metadata to GInteractions
mcols(merged_loops_coords) <- merged_loops_df

saveRDS(merged_loops_coords,
        file.path(output_dir, "non_redundant_loops.rds"))
cat(sprintf("  ✓ %s (GInteractions object)\n\n",
            "non_redundant_loops.rds"))

# --- 5C. BEDPE Export ---

cat("5C. Exporting BEDPE files...\n")

# Create BEDPE format
bedpe_df <- data.frame(
  chr1 = merged_loops_df$chr1,
  start1 = merged_loops_df$start1,
  end1 = merged_loops_df$end1,
  chr2 = merged_loops_df$chr2,
  start2 = merged_loops_df$start2,
  end2 = merged_loops_df$end2,
  name = merged_loops_df$loop_id,
  score = -log10(merged_loops_df$FDR),  # For visualization
  strand1 = ".",
  strand2 = ".",
  color = ifelse(merged_loops_df$logFC > 0, "255,0,0", "0,0,255")  # Red=up, Blue=down
)

write.table(bedpe_df,
            file.path(output_dir, "non_redundant_loops.bedpe"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
cat(sprintf("  ✓ %s (for Juicebox/IGV)\n\n",
            "non_redundant_loops.bedpe"))

# =============================================================================
# FINAL SUMMARY
# =============================================================================

cat("="*70, "\n", sep="")
cat("DOWNSTREAM ANALYSIS COMPLETE\n")
cat("="*70, "\n\n", sep="")

cat("Output directory: outputs/merged_loops/\n\n")

cat("Files generated:\n\n")

cat("Main results:\n")
cat("  1. non_redundant_loops.tsv\n")
cat("     → Merged loop list with basic info (%d loops)\n", nrow(non_redundant_df))
cat("\n")
cat("  2. characterized_loops.tsv\n")
cat("     → Complete characterization with all features (%d loops)\n", nrow(merged_loops_df))
cat("\n")
cat("  3. anchor_annotations.tsv\n")
cat("     → Per-anchor gene annotations (%d anchors)\n", nrow(anchor_annotations))
cat("\n")

if (nrow(overlap_report) > 0) {
  cat("  4. overlap_report.tsv\n")
  cat(sprintf("     → Details of merged overlapping loops (%d clusters)\n", nrow(overlap_report)))
  cat("\n")
}

cat("  5. non_redundant_loops.rds\n")
cat("     → GInteractions object for R analysis\n")
cat("\n")

cat("  6. non_redundant_loops.bedpe\n")
cat("     → BEDPE format for genome browser visualization\n")
cat("\n")

cat("Quality control:\n")
cat("  7. volcano_files_validation.tsv\n")
cat("     → Validation of all_results files for volcano plots\n")
cat("\n")

cat("  8. characterization_summary.txt\n")
cat("     → Summary statistics report\n")
cat("\n")

cat("Plots (outputs/merged_loops/plots/):\n")
cat("  9. distance_distribution.pdf\n")
cat("  10. chromosome_distribution.pdf\n")
cat("  11. gene_proximity_distribution.pdf\n")
cat("  12. loop_type_classification.pdf\n")
cat("\n")

cat("Summary statistics:\n")
cat(sprintf("  Input loops: %d (from 3 resolutions)\n", n_loops))
cat(sprintf("  Non-redundant loops: %d\n", nrow(merged_loops_df)))
cat(sprintf("  Loops removed (overlaps): %d\n", n_loops - nrow(merged_loops_df)))
cat(sprintf("  Promoter-involved loops: %d (%.1f%%)\n",
            sum(merged_loops_df$loop_type != "Enhancer-Enhancer"),
            100 * mean(merged_loops_df$loop_type != "Enhancer-Enhancer")))
cat("\n")

cat("="*70, "\n", sep="")
cat("\n")
