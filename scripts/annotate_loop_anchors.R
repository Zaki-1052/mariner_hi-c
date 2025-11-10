#!/usr/bin/env Rscript
# scripts/annotate_loop_anchors.R
# Proper ChIP-seq-Based Loop Anchor Annotation
# Author: Zakir Alibhai
# Date: 2025-11-10
#
# Purpose:
#   Accurately classify loop anchors using actual ChIP-seq data
#   instead of naive distance-to-TSS thresholds.
#
# Methodology:
#   - Promoter: H3K27ac+ AND ≤2kb from TSS
#   - Active Enhancer: H3K27ac+ AND >2kb from TSS
#   - Poised Enhancer: H3K4me1+ (no H3K27ac) AND >2kb from TSS
#   - Other/Unclassified: No ChIP-seq marks
#
# Loop Types (10 categories):
#   Promoter-Promoter, Promoter-Active_Enhancer, Promoter-Poised_Enhancer,
#   Promoter-Other, Active_Enhancer-Active_Enhancer, Active_Enhancer-Poised_Enhancer,
#   Active_Enhancer-Other, Poised_Enhancer-Poised_Enhancer, Poised_Enhancer-Other,
#   Other-Other
#
# Usage:
#   Can be sourced from other scripts or run standalone:
#   source("scripts/annotate_loop_anchors.R")
#   result <- annotate_loop_anchors(ginteractions, txdb, tss_threshold=2000)

# =============================================================================
# LIBRARY LOADING
# =============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(InteractionSet)
  library(rtracklayer)  # For importing BED files
  library(dplyr)
})

# =============================================================================
# BED FILE LOADING FUNCTIONS
# =============================================================================

#' Load H3K27ac ChIP-seq peaks (active enhancer/promoter mark)
#'
#' @param bed_path Path to H3K27ac bed file
#' @return GRanges object with H3K27ac peaks
load_h3k27ac_peaks <- function(bed_path = "220310index25H3K27acLatePeakRegions.bed") {
  if (!file.exists(bed_path)) {
    stop(sprintf("H3K27ac bed file not found: %s", bed_path))
  }

  cat(sprintf("  Loading H3K27ac peaks from: %s\n", bed_path))

  # Import bed file
  k27ac_gr <- import(bed_path, format = "bed")

  cat(sprintf("    Loaded %d H3K27ac peaks\n", length(k27ac_gr)))

  return(k27ac_gr)
}

#' Load H3K4me1 ChIP-seq peaks (enhancer mark - active + poised)
#'
#' @param bed_path Path to H3K4me1 bed file
#' @return GRanges object with H3K4me1 peaks
load_h3k4me1_peaks <- function(bed_path = "K4me1_aligned_reads_peaks.broadPeak-filtered.bed") {
  if (!file.exists(bed_path)) {
    stop(sprintf("H3K4me1 bed file not found: %s", bed_path))
  }

  cat(sprintf("  Loading H3K4me1 peaks from: %s\n", bed_path))

  # Import bed file
  k4me1_gr <- import(bed_path, format = "bed")

  cat(sprintf("    Loaded %d H3K4me1 peaks\n", length(k4me1_gr)))

  return(k4me1_gr)
}

# =============================================================================
# ANCHOR ANNOTATION FUNCTIONS
# =============================================================================

#' Annotate anchors with ChIP-seq overlaps
#'
#' @param anchor_gr GRanges object with anchor coordinates
#' @param k27ac_gr GRanges object with H3K27ac peaks
#' @param k4me1_gr GRanges object with H3K4me1 peaks
#' @return data.frame with H3K27ac_overlap, H3K4me1_overlap columns
annotate_chip_overlaps <- function(anchor_gr, k27ac_gr, k4me1_gr) {
  # Find overlaps with H3K27ac
  k27ac_overlaps <- countOverlaps(anchor_gr, k27ac_gr) > 0

  # Find overlaps with H3K4me1
  k4me1_overlaps <- countOverlaps(anchor_gr, k4me1_gr) > 0

  # Return as data.frame
  return(data.frame(
    H3K27ac_overlap = k27ac_overlaps,
    H3K4me1_overlap = k4me1_overlaps
  ))
}

#' Classify anchor type based on ChIP-seq marks and TSS distance
#'
#' @param h3k27ac_overlap Logical - overlaps H3K27ac peak
#' @param h3k4me1_overlap Logical - overlaps H3K4me1 peak
#' @param distance_to_tss Numeric - distance to nearest TSS (NA if no gene)
#' @param tss_threshold Numeric - promoter distance threshold (default 2000bp)
#' @return Character vector with anchor type
classify_anchor_type <- function(h3k27ac_overlap, h3k4me1_overlap,
                                 distance_to_tss, tss_threshold = 2000) {
  # Initialize as "Other"
  anchor_type <- rep("Other", length(h3k27ac_overlap))

  # Promoter: H3K27ac+ AND ≤2kb from TSS
  is_promoter <- h3k27ac_overlap &
                 !is.na(distance_to_tss) &
                 distance_to_tss <= tss_threshold
  anchor_type[is_promoter] <- "Promoter"

  # Active Enhancer: H3K27ac+ AND >2kb from TSS (or no nearby TSS)
  is_active_enhancer <- h3k27ac_overlap &
                        (is.na(distance_to_tss) | distance_to_tss > tss_threshold)
  anchor_type[is_active_enhancer] <- "Active_Enhancer"

  # Poised Enhancer: H3K4me1+ (no H3K27ac) AND >2kb from TSS (or no nearby TSS)
  is_poised_enhancer <- !h3k27ac_overlap &
                        h3k4me1_overlap &
                        (is.na(distance_to_tss) | distance_to_tss > tss_threshold)
  anchor_type[is_poised_enhancer] <- "Poised_Enhancer"

  return(anchor_type)
}

# =============================================================================
# LOOP TYPE CLASSIFICATION
# =============================================================================

#' Classify loop type based on anchor types
#'
#' @param anchor1_type Character vector - anchor1 types
#' @param anchor2_type Character vector - anchor2 types
#' @return Character vector with loop types (10 categories)
classify_loop_type <- function(anchor1_type, anchor2_type) {
  # Initialize
  loop_type <- rep(NA_character_, length(anchor1_type))

  # Define type hierarchy for consistent ordering
  type_order <- c("Promoter", "Active_Enhancer", "Poised_Enhancer", "Other")

  for (i in seq_along(anchor1_type)) {
    t1 <- anchor1_type[i]
    t2 <- anchor2_type[i]

    # Ensure consistent ordering (alphabetical for ties)
    types_sorted <- sort(c(t1, t2))

    # Classify based on sorted pair
    if (t1 == t2) {
      # Homotypic loops
      loop_type[i] <- sprintf("%s-%s", t1, t2)
    } else {
      # Heterotypic loops - maintain consistent ordering
      if (match(t1, type_order) < match(t2, type_order)) {
        loop_type[i] <- sprintf("%s-%s", t1, t2)
      } else {
        loop_type[i] <- sprintf("%s-%s", t2, t1)
      }
    }
  }

  return(loop_type)
}

# =============================================================================
# MAIN ANNOTATION FUNCTION
# =============================================================================

#' Annotate loop anchors with ChIP-seq-based classifications
#'
#' @param gi GInteractions object with loops
#' @param txdb TxDb object with gene annotations (for TSS)
#' @param k27ac_path Path to H3K27ac bed file (default in project root)
#' @param k4me1_path Path to H3K4me1 bed file (default in project root)
#' @param tss_threshold Promoter distance threshold in bp (default 2000)
#' @param existing_tss_distances Optional - pre-computed TSS distances (data.frame)
#' @return List with annotated GInteractions and summary statistics
annotate_loop_anchors <- function(gi, txdb,
                                   k27ac_path = "220310index25H3K27acLatePeakRegions.bed",
                                   k4me1_path = "K4me1_aligned_reads_peaks.broadPeak-filtered.bed",
                                   tss_threshold = 2000,
                                   existing_tss_distances = NULL) {

  cat("\n========================================\n")
  cat("ChIP-seq-Based Loop Anchor Annotation\n")
  cat("========================================\n\n")

  cat(sprintf("Parameters:\n"))
  cat(sprintf("  TSS threshold: %d bp\n", tss_threshold))
  cat(sprintf("  H3K27ac bed: %s\n", k27ac_path))
  cat(sprintf("  H3K4me1 bed: %s\n\n", k4me1_path))

  # --- Step 1: Load ChIP-seq data ---
  cat("Step 1: Loading ChIP-seq peak files...\n")
  k27ac_gr <- load_h3k27ac_peaks(k27ac_path)
  k4me1_gr <- load_h3k4me1_peaks(k4me1_path)
  cat("\n")

  # --- Step 2: Extract anchors ---
  cat("Step 2: Extracting loop anchors...\n")
  anchor1_gr <- anchors(gi, "first")
  anchor2_gr <- anchors(gi, "second")
  cat(sprintf("  Extracted %d loop anchors (first)\n", length(anchor1_gr)))
  cat(sprintf("  Extracted %d loop anchors (second)\n\n", length(anchor2_gr)))

  # --- Step 3: Compute or use existing TSS distances ---
  if (!is.null(existing_tss_distances)) {
    cat("Step 3: Using pre-computed TSS distances...\n")
    anchor1_distance_to_tss <- existing_tss_distances$anchor1_distance_to_tss
    anchor2_distance_to_tss <- existing_tss_distances$anchor2_distance_to_tss
  } else {
    cat("Step 3: Computing TSS distances...\n")

    # Extract genes and TSS
    genes <- genes(txdb)
    tss <- resize(genes, width = 1, fix = "start")
    cat(sprintf("  Loaded %d genes\n", length(genes)))

    # Find nearest TSS for each anchor
    cat("  Finding nearest TSS for anchor1...\n")
    nearest1 <- distanceToNearest(anchor1_gr, tss)
    anchor1_distance_to_tss <- rep(NA_real_, length(anchor1_gr))
    if (length(nearest1) > 0) {
      anchor1_distance_to_tss[queryHits(nearest1)] <- mcols(nearest1)$distance
    }

    cat("  Finding nearest TSS for anchor2...\n")
    nearest2 <- distanceToNearest(anchor2_gr, tss)
    anchor2_distance_to_tss <- rep(NA_real_, length(anchor2_gr))
    if (length(nearest2) > 0) {
      anchor2_distance_to_tss[queryHits(nearest2)] <- mcols(nearest2)$distance
    }
  }
  cat("\n")

  # --- Step 4: Annotate ChIP-seq overlaps ---
  cat("Step 4: Annotating ChIP-seq overlaps...\n")
  cat("  Checking anchor1 overlaps...\n")
  anchor1_chip <- annotate_chip_overlaps(anchor1_gr, k27ac_gr, k4me1_gr)

  cat("  Checking anchor2 overlaps...\n")
  anchor2_chip <- annotate_chip_overlaps(anchor2_gr, k27ac_gr, k4me1_gr)

  cat(sprintf("  Anchor1: %d H3K27ac+, %d H3K4me1+\n",
              sum(anchor1_chip$H3K27ac_overlap),
              sum(anchor1_chip$H3K4me1_overlap)))
  cat(sprintf("  Anchor2: %d H3K27ac+, %d H3K4me1+\n\n",
              sum(anchor2_chip$H3K27ac_overlap),
              sum(anchor2_chip$H3K4me1_overlap)))

  # --- Step 5: Classify anchor types ---
  cat("Step 5: Classifying anchor types...\n")
  anchor1_type <- classify_anchor_type(
    anchor1_chip$H3K27ac_overlap,
    anchor1_chip$H3K4me1_overlap,
    anchor1_distance_to_tss,
    tss_threshold
  )

  anchor2_type <- classify_anchor_type(
    anchor2_chip$H3K27ac_overlap,
    anchor2_chip$H3K4me1_overlap,
    anchor2_distance_to_tss,
    tss_threshold
  )

  # Print anchor type distribution
  cat("\n  Anchor1 type distribution:\n")
  for (type in c("Promoter", "Active_Enhancer", "Poised_Enhancer", "Other")) {
    count <- sum(anchor1_type == type)
    pct <- 100 * mean(anchor1_type == type)
    cat(sprintf("    %-20s: %5d (%.1f%%)\n", type, count, pct))
  }

  cat("\n  Anchor2 type distribution:\n")
  for (type in c("Promoter", "Active_Enhancer", "Poised_Enhancer", "Other")) {
    count <- sum(anchor2_type == type)
    pct <- 100 * mean(anchor2_type == type)
    cat(sprintf("    %-20s: %5d (%.1f%%)\n", type, count, pct))
  }
  cat("\n")

  # --- Step 6: Classify loop types ---
  cat("Step 6: Classifying loop types...\n")
  loop_type <- classify_loop_type(anchor1_type, anchor2_type)

  # Print loop type distribution
  loop_type_counts <- table(loop_type)
  loop_type_counts <- loop_type_counts[order(loop_type_counts, decreasing = TRUE)]

  cat("\n  Loop type distribution:\n")
  for (i in seq_along(loop_type_counts)) {
    type_name <- names(loop_type_counts)[i]
    count <- loop_type_counts[i]
    pct <- 100 * count / length(loop_type)
    cat(sprintf("    %-35s: %5d (%.1f%%)\n", type_name, count, pct))
  }
  cat("\n")

  # --- Step 7: Add annotations to GInteractions ---
  cat("Step 7: Adding annotations to GInteractions object...\n")

  # Add to mcols
  mcols(gi)$anchor1_H3K27ac_overlap <- anchor1_chip$H3K27ac_overlap
  mcols(gi)$anchor1_H3K4me1_overlap <- anchor1_chip$H3K4me1_overlap
  mcols(gi)$anchor1_type <- anchor1_type
  mcols(gi)$anchor1_is_promoter <- (anchor1_type == "Promoter")  # Backward compatibility

  mcols(gi)$anchor2_H3K27ac_overlap <- anchor2_chip$H3K27ac_overlap
  mcols(gi)$anchor2_H3K4me1_overlap <- anchor2_chip$H3K4me1_overlap
  mcols(gi)$anchor2_type <- anchor2_type
  mcols(gi)$anchor2_is_promoter <- (anchor2_type == "Promoter")  # Backward compatibility

  mcols(gi)$loop_type <- loop_type

  cat("  Added columns to GInteractions mcols:\n")
  cat("    - anchor1_H3K27ac_overlap, anchor2_H3K27ac_overlap\n")
  cat("    - anchor1_H3K4me1_overlap, anchor2_H3K4me1_overlap\n")
  cat("    - anchor1_type, anchor2_type\n")
  cat("    - anchor1_is_promoter, anchor2_is_promoter (backward compat)\n")
  cat("    - loop_type\n\n")

  # --- Step 8: Prepare summary statistics ---
  summary_stats <- list(
    total_loops = length(gi),
    tss_threshold_bp = tss_threshold,
    anchor1_summary = list(
      h3k27ac_positive = sum(anchor1_chip$H3K27ac_overlap),
      h3k4me1_positive = sum(anchor1_chip$H3K4me1_overlap),
      promoter = sum(anchor1_type == "Promoter"),
      active_enhancer = sum(anchor1_type == "Active_Enhancer"),
      poised_enhancer = sum(anchor1_type == "Poised_Enhancer"),
      other = sum(anchor1_type == "Other")
    ),
    anchor2_summary = list(
      h3k27ac_positive = sum(anchor2_chip$H3K27ac_overlap),
      h3k4me1_positive = sum(anchor2_chip$H3K4me1_overlap),
      promoter = sum(anchor2_type == "Promoter"),
      active_enhancer = sum(anchor2_type == "Active_Enhancer"),
      poised_enhancer = sum(anchor2_type == "Poised_Enhancer"),
      other = sum(anchor2_type == "Other")
    ),
    loop_type_summary = as.list(loop_type_counts)
  )

  cat("========================================\n")
  cat("Annotation Complete\n")
  cat("========================================\n\n")

  return(list(
    ginteractions = gi,
    summary = summary_stats
  ))
}

# =============================================================================
# CONVENIENCE WRAPPER FOR DATA.FRAME
# =============================================================================

#' Add ChIP-seq annotations to existing data.frame with loop information
#'
#' @param loops_df data.frame with loop coordinates and TSS distances
#' @param k27ac_path Path to H3K27ac bed file
#' @param k4me1_path Path to H3K4me1 bed file
#' @param tss_threshold Promoter distance threshold in bp
#' @return data.frame with added annotation columns
annotate_loops_dataframe <- function(loops_df,
                                     k27ac_path = "220310index25H3K27acLatePeakRegions.bed",
                                     k4me1_path = "K4me1_aligned_reads_peaks.broadPeak-filtered.bed",
                                     tss_threshold = 2000) {

  cat("\nAnnotating loops data.frame with ChIP-seq data...\n")

  # Load ChIP-seq data
  k27ac_gr <- load_h3k27ac_peaks(k27ac_path)
  k4me1_gr <- load_h3k4me1_peaks(k4me1_path)

  # Create GRanges for anchors
  anchor1_gr <- GRanges(
    seqnames = loops_df$anchor1_chr,
    ranges = IRanges(start = loops_df$anchor1_start, end = loops_df$anchor1_end)
  )

  anchor2_gr <- GRanges(
    seqnames = loops_df$anchor2_chr,
    ranges = IRanges(start = loops_df$anchor2_start, end = loops_df$anchor2_end)
  )

  # Annotate ChIP overlaps
  anchor1_chip <- annotate_chip_overlaps(anchor1_gr, k27ac_gr, k4me1_gr)
  anchor2_chip <- annotate_chip_overlaps(anchor2_gr, k27ac_gr, k4me1_gr)

  # Classify anchor types
  anchor1_type <- classify_anchor_type(
    anchor1_chip$H3K27ac_overlap,
    anchor1_chip$H3K4me1_overlap,
    loops_df$anchor1_distance_to_tss,
    tss_threshold
  )

  anchor2_type <- classify_anchor_type(
    anchor2_chip$H3K27ac_overlap,
    anchor2_chip$H3K4me1_overlap,
    loops_df$anchor2_distance_to_tss,
    tss_threshold
  )

  # Classify loop types
  loop_type <- classify_loop_type(anchor1_type, anchor2_type)

  # Add to data.frame
  loops_df$anchor1_H3K27ac_overlap <- anchor1_chip$H3K27ac_overlap
  loops_df$anchor1_H3K4me1_overlap <- anchor1_chip$H3K4me1_overlap
  loops_df$anchor1_type <- anchor1_type
  loops_df$anchor1_is_promoter <- (anchor1_type == "Promoter")

  loops_df$anchor2_H3K27ac_overlap <- anchor2_chip$H3K27ac_overlap
  loops_df$anchor2_H3K4me1_overlap <- anchor2_chip$H3K4me1_overlap
  loops_df$anchor2_type <- anchor2_type
  loops_df$anchor2_is_promoter <- (anchor2_type == "Promoter")

  loops_df$loop_type <- loop_type

  cat("✓ Annotation complete\n\n")

  return(loops_df)
}

# =============================================================================
# STANDALONE EXECUTION
# =============================================================================

if (!interactive() && sys.nframe() == 0) {
  cat("\n")
  cat("========================================\n")
  cat("Standalone Execution Not Recommended\n")
  cat("========================================\n\n")
  cat("This script is designed to be sourced from other scripts.\n")
  cat("See scripts/downstream_analysis.R for usage examples.\n\n")
  cat("To use in your own script:\n")
  cat("  source('scripts/annotate_loop_anchors.R')\n")
  cat("  result <- annotate_loop_anchors(gi, txdb)\n\n")
}
