#!/usr/bin/env Rscript
# scripts/annotate_loops_extended.R
# Extended ChIP-seq-Based Loop Anchor Annotation with Chromatin State Categories
# Author: Zakir Alibhai
# Date: 2025-01-05 (Updated: 2026-01-12)
#
# Purpose:
#   Annotate loop anchors with chromatin state categories using ChIP-seq data.
#   Categories are biologically defined based on histone modifications and TSS proximity.
#
# Anchor Categories (7 types, priority order):
#   1. Active_Promoter:    H3K4me3+ AND NOT H3K27me3 AND ≤2kb from TSS
#   2. Repressed_Promoter: H3K27me3+ AND NOT H3K27ac AND ≤2kb from TSS
#   3. Bivalent:           K4me3+K27me3 overlap (Addison file, early timepoint)
#   4. Polycomb:           H3K27me3+ AND >2kb from TSS (distal repressive)
#   5. Active_Enhancer:    H3K27ac+ AND >2kb from TSS
#   6. Poised_Enhancer:    H3K4me1+ AND NOT H3K27ac AND NOT H3K27me3 AND >2kb from TSS
#   7. Other:              No ChIP-seq marks / structural elements
#
# Loop Types (28 combinations):
#   AP-AP, AP-RP, AP-B, AP-Pc, AP-AE, AP-PE, AP-O
#   RP-RP, RP-B, RP-Pc, RP-AE, RP-PE, RP-O
#   B-B, B-Pc, B-AE, B-PE, B-O
#   Pc-Pc, Pc-AE, Pc-PE, Pc-O
#   AE-AE, AE-PE, AE-O
#   PE-PE, PE-O
#   O-O
#
# Usage:
#   Rscript scripts/annotate_loops_extended.R [--input FILE] [--output DIR]
#
#   Default input:  25042-late_outputs/merged_loops/non_redundant_loops.tsv
#   Default output: outputs/loop_annotation_extended/

# =============================================================================
# LIBRARY LOADING
# =============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(InteractionSet)
  library(rtracklayer)
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(RColorBrewer)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
})

# =============================================================================
# CONFIGURATION
# =============================================================================

# Default ChIP-seq peak files
# Note: These paths can be easily updated when standardized files are provided
DEFAULT_H3K27AC_PATH <- "peaks/220310index25H3K27acLatePeakRegions.bed"
DEFAULT_H3K27ME3_PATH <- "peaks/220310index29H3K27me3LatePeakRegions.bed"
DEFAULT_H3K4ME1_PATH <- "peaks/K4me1_aligned_reads_peaks.broadPeak-filtered.bed"
DEFAULT_H3K4ME3_PATH <- "consensus_H3K4me3_late_peaks.bed"  # Adult H3K4me3 consensus

# Bivalent regions: K4me3+K27me3 overlap from early timepoint (Addison file)
DEFAULT_BIVALENT_PATH <- "250224AddisonH3K4me3H3K27me3Early.bed"

# Default input/output paths
DEFAULT_INPUT_FILE <- "25042-late_outputs/merged_loops/non_redundant_loops.tsv"
DEFAULT_OUTPUT_DIR <- "outputs/loop_annotation_extended"

# Anchor type hierarchy (for consistent loop type ordering)
# 7 categories reflecting chromatin states
ANCHOR_TYPE_ORDER <- c("Active_Promoter", "Repressed_Promoter", "Bivalent",
                       "Polycomb", "Active_Enhancer", "Poised_Enhancer", "Other")

# Color scheme for visualizations
ANCHOR_COLORS <- c(
  "Active_Promoter" = "#e41a1c",     # Red - active transcription
  "Repressed_Promoter" = "#756bb1", # Purple - Polycomb-silenced promoter
  "Bivalent" = "#984ea3",            # Magenta - developmental poised
  "Polycomb" = "#4daf4a",            # Green - distal repressive
  "Active_Enhancer" = "#377eb8",     # Blue - active enhancer
  "Poised_Enhancer" = "#ff7f00",     # Orange - primed enhancer
  "Other" = "#999999"                # Gray - structural/unmarked
)

# =============================================================================
# ChIP-seq PEAK LOADING FUNCTIONS
# =============================================================================

#' Load H3K27ac ChIP-seq peaks
#'
#' @param bed_path Path to H3K27ac bed file
#' @return GRanges object
load_h3k27ac_peaks <- function(bed_path = DEFAULT_H3K27AC_PATH) {
  if (!file.exists(bed_path)) {
    stop(sprintf("H3K27ac bed file not found: %s", bed_path))
  }
  cat(sprintf("  Loading H3K27ac peaks from: %s\n", bed_path))

  # Read as table first (handles non-standard BED formats)
  df <- read.table(bed_path, sep = "\t", header = FALSE,
                   stringsAsFactors = FALSE)
  gr <- GRanges(
    seqnames = df$V1,
    ranges = IRanges(start = df$V2, end = df$V3)
  )

  cat(sprintf("    Loaded %d peaks\n", length(gr)))
  return(gr)
}

#' Load H3K27me3 ChIP-seq peaks
#'
#' @param bed_path Path to H3K27me3 bed file
#' @return GRanges object
load_h3k27me3_peaks <- function(bed_path = DEFAULT_H3K27ME3_PATH) {
  if (!file.exists(bed_path)) {
    stop(sprintf("H3K27me3 bed file not found: %s", bed_path))
  }
  cat(sprintf("  Loading H3K27me3 peaks from: %s\n", bed_path))

  # Read as table first (handles non-standard BED formats)
  df <- read.table(bed_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  gr <- GRanges(
    seqnames = df$V1,
    ranges = IRanges(start = df$V2, end = df$V3)
  )

  cat(sprintf("    Loaded %d peaks\n", length(gr)))
  return(gr)
}

#' Load H3K4me1 ChIP-seq peaks
#'
#' @param bed_path Path to H3K4me1 bed file
#' @return GRanges object
load_h3k4me1_peaks <- function(bed_path = DEFAULT_H3K4ME1_PATH) {
  if (!file.exists(bed_path)) {
    stop(sprintf("H3K4me1 bed file not found: %s", bed_path))
  }
  cat(sprintf("  Loading H3K4me1 peaks from: %s\n", bed_path))

  # Read as table first (handles non-standard BED formats)
  df <- read.table(bed_path, sep = "\t", header = FALSE,
                   stringsAsFactors = FALSE)
  gr <- GRanges(
    seqnames = df$V1,
    ranges = IRanges(start = df$V2, end = df$V3)
  )

  cat(sprintf("    Loaded %d peaks\n", length(gr)))
  return(gr)
}

#' Load H3K4me3 ChIP-seq peaks (active promoter mark)
#'
#' @param bed_path Path to H3K4me3 bed file (consensus peaks)
#' @return GRanges object
load_h3k4me3_peaks <- function(bed_path = DEFAULT_H3K4ME3_PATH) {
  if (!file.exists(bed_path)) {
    stop(sprintf("H3K4me3 bed file not found: %s", bed_path))
  }
  cat(sprintf("  Loading H3K4me3 peaks from: %s\n", bed_path))

  # Read as table first (handles non-standard BED formats)
  df <- read.table(bed_path, sep = "\t", header = FALSE,
                   stringsAsFactors = FALSE)
  gr <- GRanges(
    seqnames = df$V1,
    ranges = IRanges(start = df$V2, end = df$V3)
  )

  cat(sprintf("    Loaded %d peaks\n", length(gr)))
  return(gr)
}

#' Load bivalent (K4me3+K27me3) regions
#'
#' Bivalent domains contain both active (H3K4me3) and repressive (H3K27me3) marks.
#' These regions are characteristic of developmental poised states.
#'
#' @param bed_path Path to K4me3+K27me3 overlap bed file (Addison file)
#' @return GRanges object
load_bivalent_peaks <- function(bed_path = DEFAULT_BIVALENT_PATH) {
  if (!file.exists(bed_path)) {
    stop(sprintf("Bivalent bed file not found: %s", bed_path))
  }
  cat(sprintf("  Loading bivalent (K4me3+K27me3) regions from: %s\n", bed_path))

  # Read as table first (handles non-standard BED formats)
  df <- read.table(bed_path, sep = "\t", header = FALSE,
                   stringsAsFactors = FALSE)
  gr <- GRanges(
    seqnames = df$V1,
    ranges = IRanges(start = df$V2, end = df$V3)
  )

  cat(sprintf("    Loaded %d regions\n", length(gr)))
  return(gr)
}

# =============================================================================
# ANCHOR ANNOTATION FUNCTIONS
# =============================================================================

#' Annotate anchors with all ChIP-seq overlaps
#'
#' @param anchor_gr GRanges with anchor coordinates
#' @param k27ac_gr H3K27ac peaks
#' @param k27me3_gr H3K27me3 peaks
#' @param k4me1_gr H3K4me1 peaks
#' @param k4me3_gr H3K4me3 peaks (active promoter mark)
#' @param bivalent_gr Bivalent regions (K4me3+K27me3 overlap)
#' @return data.frame with overlap columns
annotate_chip_overlaps_extended <- function(anchor_gr, k27ac_gr, k27me3_gr,
                                            k4me1_gr, k4me3_gr, bivalent_gr) {
  data.frame(
    H3K27ac_overlap = countOverlaps(anchor_gr, k27ac_gr) > 0,
    H3K27me3_overlap = countOverlaps(anchor_gr, k27me3_gr) > 0,
    H3K4me1_overlap = countOverlaps(anchor_gr, k4me1_gr) > 0,
    H3K4me3_overlap = countOverlaps(anchor_gr, k4me3_gr) > 0,
    Bivalent_overlap = countOverlaps(anchor_gr, bivalent_gr) > 0
  )
}

#' Classify anchor type with chromatin state categories (7 types)
#'
#' Priority order:
#'   Active_Promoter > Repressed_Promoter > Bivalent > Polycomb >
#'   Active_Enhancer > Poised_Enhancer > Other
#'
#' Biological rationale:
#'   - Active promoters: H3K4me3 is the canonical active promoter mark
#'     (housekeeping genes can have H3K4me3 without H3K27ac)
#'   - Repressed promoters: H3K27me3 at TSS indicates Polycomb silencing
#'   - Bivalent: K4me3+K27me3 overlap marks developmental poised domains
#'   - Polycomb: Distal H3K27me3 regions (long-range repressive loops)
#'   - Poised enhancers: H3K4me1 without repressive or active marks
#'
#' @param h3k27ac_overlap Logical - overlaps H3K27ac peak
#' @param h3k27me3_overlap Logical - overlaps H3K27me3 peak
#' @param h3k4me1_overlap Logical - overlaps H3K4me1 peak
#' @param h3k4me3_overlap Logical - overlaps H3K4me3 peak (active promoter)
#' @param bivalent_overlap Logical - overlaps K4me3+K27me3 region (Addison file)
#' @param distance_to_tss Numeric - distance to nearest TSS
#' @param tss_threshold Numeric - promoter distance threshold (default 2000bp)
#' @return Character vector with anchor types
classify_anchor_type_extended <- function(h3k27ac_overlap, h3k27me3_overlap,
                                          h3k4me1_overlap, h3k4me3_overlap,
                                          bivalent_overlap, distance_to_tss,
                                          tss_threshold = 2000) {
  n <- length(h3k27ac_overlap)
  anchor_type <- rep("Other", n)

  # 1. Active_Promoter: H3K4me3+ AND NOT H3K27me3 AND ≤2kb from TSS
  # H3K4me3 is the canonical active promoter mark; H3K27ac not required
  # (housekeeping genes can be H3K4me3+ without K27ac)
  is_active_promoter <- h3k4me3_overlap & !h3k27me3_overlap &
                        !is.na(distance_to_tss) &
                        distance_to_tss <= tss_threshold
  anchor_type[is_active_promoter] <- "Active_Promoter"

  # 2. Repressed_Promoter: H3K27me3+ AND NOT H3K27ac AND ≤2kb from TSS
  # Polycomb-silenced promoter near TSS
  is_repressed_promoter <- !is_active_promoter &
                           h3k27me3_overlap & !h3k27ac_overlap &
                           !is.na(distance_to_tss) &
                           distance_to_tss <= tss_threshold
  anchor_type[is_repressed_promoter] <- "Repressed_Promoter"

  # 3. Bivalent: K4me3+K27me3 overlap (not already classified as promoter)
  # Developmental poised domains from early timepoint (Addison file)
  is_bivalent <- !is_active_promoter & !is_repressed_promoter & bivalent_overlap
  anchor_type[is_bivalent] <- "Bivalent"

  # 4. Polycomb: H3K27me3+ AND >2kb from TSS (distal repressive)
  # Polycomb loops tend to be long-range
  is_polycomb <- !is_active_promoter & !is_repressed_promoter & !is_bivalent &
                 h3k27me3_overlap &
                 (is.na(distance_to_tss) | distance_to_tss > tss_threshold)
  anchor_type[is_polycomb] <- "Polycomb"

  # 5. Active_Enhancer: H3K27ac+ AND >2kb from TSS
  is_active_enhancer <- !is_active_promoter & !is_repressed_promoter &
                        !is_bivalent & !is_polycomb &
                        h3k27ac_overlap &
                        (is.na(distance_to_tss) | distance_to_tss > tss_threshold)
  anchor_type[is_active_enhancer] <- "Active_Enhancer"

  # 6. Poised_Enhancer: H3K4me1+ AND NOT H3K27ac AND NOT H3K27me3 AND >2kb
  # Primed enhancer without active or repressive marks
  is_poised_enhancer <- !is_active_promoter & !is_repressed_promoter &
                        !is_bivalent & !is_polycomb & !is_active_enhancer &
                        h3k4me1_overlap & !h3k27ac_overlap & !h3k27me3_overlap &
                        (is.na(distance_to_tss) | distance_to_tss > tss_threshold)
  anchor_type[is_poised_enhancer] <- "Poised_Enhancer"

  # 7. Other: default (no ChIP-seq marks - structural elements, CTCF sites, etc.)
  return(anchor_type)
}

# =============================================================================
# LOOP TYPE CLASSIFICATION
# =============================================================================

#' Classify loop type based on anchor types (21 combinations)
#'
#' @param anchor1_type Character vector - anchor1 types
#' @param anchor2_type Character vector - anchor2 types
#' @return Character vector with loop types
classify_loop_type_extended <- function(anchor1_type, anchor2_type) {
  n <- length(anchor1_type)
  loop_type <- rep(NA_character_, n)

  for (i in seq_len(n)) {
    t1 <- anchor1_type[i]
    t2 <- anchor2_type[i]

    if (t1 == t2) {
      # Homotypic loops
      loop_type[i] <- sprintf("%s-%s", t1, t2)
    } else {
      # Heterotypic loops - order by hierarchy
      pos1 <- match(t1, ANCHOR_TYPE_ORDER)
      pos2 <- match(t2, ANCHOR_TYPE_ORDER)

      if (pos1 < pos2) {
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

#' Annotate loops with chromatin state categories
#'
#' @param input_file Path to input TSV (non_redundant_loops.tsv or similar)
#' @param output_dir Output directory
#' @param h3k27ac_path Path to H3K27ac peaks
#' @param h3k27me3_path Path to H3K27me3 peaks
#' @param h3k4me1_path Path to H3K4me1 peaks
#' @param h3k4me3_path Path to H3K4me3 peaks (active promoter mark)
#' @param bivalent_path Path to bivalent regions (K4me3+K27me3 overlap)
#' @param tss_threshold Promoter distance threshold (bp)
#' @return data.frame with annotated loops
annotate_loops_extended <- function(
  input_file = DEFAULT_INPUT_FILE,
  output_dir = DEFAULT_OUTPUT_DIR,
  h3k27ac_path = DEFAULT_H3K27AC_PATH,
  h3k27me3_path = DEFAULT_H3K27ME3_PATH,
  h3k4me1_path = DEFAULT_H3K4ME1_PATH,
  h3k4me3_path = DEFAULT_H3K4ME3_PATH,
  bivalent_path = DEFAULT_BIVALENT_PATH,
  tss_threshold = 2000
) {

  cat("\n")
  cat("========================================\n")
  cat("Extended Loop Anchor Annotation\n")
  cat("========================================\n\n")

  cat("Configuration:\n")
  cat(sprintf("  Input file:     %s\n", input_file))
  cat(sprintf("  Output dir:     %s\n", output_dir))
  cat(sprintf("  TSS threshold:  %d bp\n\n", tss_threshold))

  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "plots"), recursive = TRUE, showWarnings = FALSE)

  # --- Step 1: Load input data ---
  cat("Step 1: Loading input data...\n")

  if (!file.exists(input_file)) {
    stop(sprintf("Input file not found: %s", input_file))
  }

  loops_df <- read.table(input_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  cat(sprintf("  Loaded %d loops\n", nrow(loops_df)))

  # Check for required columns
  required_cols <- c("chr1", "start1", "end1", "chr2", "start2", "end2", "logFC", "direction")
  missing_cols <- setdiff(required_cols, colnames(loops_df))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }

  cat(sprintf("  Up in mutant: %d\n", sum(loops_df$direction == "up_in_mutant")))
  cat(sprintf("  Down in mutant: %d\n\n", sum(loops_df$direction == "down_in_mutant")))

  # --- Step 2: Load ChIP-seq peaks ---
  cat("Step 2: Loading ChIP-seq peak files...\n")
  k27ac_gr <- load_h3k27ac_peaks(h3k27ac_path)
  k27me3_gr <- load_h3k27me3_peaks(h3k27me3_path)
  k4me1_gr <- load_h3k4me1_peaks(h3k4me1_path)
  k4me3_gr <- load_h3k4me3_peaks(h3k4me3_path)
  bivalent_gr <- load_bivalent_peaks(bivalent_path)
  cat("\n")

  # --- Step 3: Create anchor GRanges ---
  cat("Step 3: Creating anchor GRanges...\n")

  anchor1_gr <- GRanges(
    seqnames = loops_df$chr1,
    ranges = IRanges(start = loops_df$start1, end = loops_df$end1)
  )

  anchor2_gr <- GRanges(
    seqnames = loops_df$chr2,
    ranges = IRanges(start = loops_df$start2, end = loops_df$end2)
  )

  cat(sprintf("  Anchor1: %d regions\n", length(anchor1_gr)))
  cat(sprintf("  Anchor2: %d regions\n\n", length(anchor2_gr)))

  # --- Step 4: Compute TSS distances ---
  cat("Step 4: Computing TSS distances...\n")

  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  genes <- genes(txdb)
  tss <- resize(genes, width = 1, fix = "start")
  cat(sprintf("  Loaded %d genes\n", length(genes)))

  # Anchor1 TSS distance
  nearest1 <- distanceToNearest(anchor1_gr, tss)
  anchor1_distance_to_tss <- rep(NA_real_, length(anchor1_gr))
  if (length(nearest1) > 0) {
    anchor1_distance_to_tss[queryHits(nearest1)] <- mcols(nearest1)$distance
  }

  # Anchor2 TSS distance
  nearest2 <- distanceToNearest(anchor2_gr, tss)
  anchor2_distance_to_tss <- rep(NA_real_, length(anchor2_gr))
  if (length(nearest2) > 0) {
    anchor2_distance_to_tss[queryHits(nearest2)] <- mcols(nearest2)$distance
  }

  cat("\n")

  # --- Step 5: Annotate ChIP-seq overlaps ---
  cat("Step 5: Annotating ChIP-seq overlaps...\n")

  anchor1_chip <- annotate_chip_overlaps_extended(anchor1_gr, k27ac_gr, k27me3_gr,
                                                  k4me1_gr, k4me3_gr, bivalent_gr)
  anchor2_chip <- annotate_chip_overlaps_extended(anchor2_gr, k27ac_gr, k27me3_gr,
                                                  k4me1_gr, k4me3_gr, bivalent_gr)

  cat("  Anchor1 overlaps:\n")
  cat(sprintf("    H3K27ac+:  %d (%.1f%%)\n",
              sum(anchor1_chip$H3K27ac_overlap),
              100 * mean(anchor1_chip$H3K27ac_overlap)))
  cat(sprintf("    H3K27me3+: %d (%.1f%%)\n",
              sum(anchor1_chip$H3K27me3_overlap),
              100 * mean(anchor1_chip$H3K27me3_overlap)))
  cat(sprintf("    H3K4me1+:  %d (%.1f%%)\n",
              sum(anchor1_chip$H3K4me1_overlap),
              100 * mean(anchor1_chip$H3K4me1_overlap)))
  cat(sprintf("    H3K4me3+:  %d (%.1f%%)\n",
              sum(anchor1_chip$H3K4me3_overlap),
              100 * mean(anchor1_chip$H3K4me3_overlap)))
  cat(sprintf("    Bivalent:  %d (%.1f%%)\n",
              sum(anchor1_chip$Bivalent_overlap),
              100 * mean(anchor1_chip$Bivalent_overlap)))

  cat("  Anchor2 overlaps:\n")
  cat(sprintf("    H3K27ac+:  %d (%.1f%%)\n",
              sum(anchor2_chip$H3K27ac_overlap),
              100 * mean(anchor2_chip$H3K27ac_overlap)))
  cat(sprintf("    H3K27me3+: %d (%.1f%%)\n",
              sum(anchor2_chip$H3K27me3_overlap),
              100 * mean(anchor2_chip$H3K27me3_overlap)))
  cat(sprintf("    H3K4me1+:  %d (%.1f%%)\n",
              sum(anchor2_chip$H3K4me1_overlap),
              100 * mean(anchor2_chip$H3K4me1_overlap)))
  cat(sprintf("    H3K4me3+:  %d (%.1f%%)\n",
              sum(anchor2_chip$H3K4me3_overlap),
              100 * mean(anchor2_chip$H3K4me3_overlap)))
  cat(sprintf("    Bivalent:  %d (%.1f%%)\n\n",
              sum(anchor2_chip$Bivalent_overlap),
              100 * mean(anchor2_chip$Bivalent_overlap)))

  # --- Step 6: Classify anchor types ---
  cat("Step 6: Classifying anchor types (7 categories)...\n")

  anchor1_type <- classify_anchor_type_extended(
    anchor1_chip$H3K27ac_overlap,
    anchor1_chip$H3K27me3_overlap,
    anchor1_chip$H3K4me1_overlap,
    anchor1_chip$H3K4me3_overlap,
    anchor1_chip$Bivalent_overlap,
    anchor1_distance_to_tss,
    tss_threshold
  )

  anchor2_type <- classify_anchor_type_extended(
    anchor2_chip$H3K27ac_overlap,
    anchor2_chip$H3K27me3_overlap,
    anchor2_chip$H3K4me1_overlap,
    anchor2_chip$H3K4me3_overlap,
    anchor2_chip$Bivalent_overlap,
    anchor2_distance_to_tss,
    tss_threshold
  )

  cat("\n  Anchor1 type distribution:\n")
  for (type in ANCHOR_TYPE_ORDER) {
    count <- sum(anchor1_type == type)
    pct <- 100 * mean(anchor1_type == type)
    cat(sprintf("    %-20s: %5d (%.1f%%)\n", type, count, pct))
  }

  cat("\n  Anchor2 type distribution:\n")
  for (type in ANCHOR_TYPE_ORDER) {
    count <- sum(anchor2_type == type)
    pct <- 100 * mean(anchor2_type == type)
    cat(sprintf("    %-20s: %5d (%.1f%%)\n", type, count, pct))
  }
  cat("\n")

  # --- Step 7: Classify loop types ---
  cat("Step 7: Classifying loop types (28 combinations)...\n")

  loop_type <- classify_loop_type_extended(anchor1_type, anchor2_type)

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

  # --- Step 8: Add annotations to data.frame ---
  cat("Step 8: Adding annotations to data.frame...\n")

  # ChIP-seq overlaps
  loops_df$anchor1_H3K27ac_overlap <- anchor1_chip$H3K27ac_overlap
  loops_df$anchor1_H3K27me3_overlap <- anchor1_chip$H3K27me3_overlap
  loops_df$anchor1_H3K4me1_overlap <- anchor1_chip$H3K4me1_overlap
  loops_df$anchor1_H3K4me3_overlap <- anchor1_chip$H3K4me3_overlap
  loops_df$anchor1_Bivalent_overlap <- anchor1_chip$Bivalent_overlap

  loops_df$anchor2_H3K27ac_overlap <- anchor2_chip$H3K27ac_overlap
  loops_df$anchor2_H3K27me3_overlap <- anchor2_chip$H3K27me3_overlap
  loops_df$anchor2_H3K4me1_overlap <- anchor2_chip$H3K4me1_overlap
  loops_df$anchor2_H3K4me3_overlap <- anchor2_chip$H3K4me3_overlap
  loops_df$anchor2_Bivalent_overlap <- anchor2_chip$Bivalent_overlap

  # TSS distances
  loops_df$anchor1_distance_to_tss_ext <- anchor1_distance_to_tss
  loops_df$anchor2_distance_to_tss_ext <- anchor2_distance_to_tss

  # Anchor types (extended)
  loops_df$anchor1_type_extended <- anchor1_type
  loops_df$anchor2_type_extended <- anchor2_type

  # Loop type (extended)
  loops_df$loop_type_extended <- loop_type

  cat("  Added columns:\n")
  cat("    - anchor1/2_H3K27ac_overlap\n")
  cat("    - anchor1/2_H3K27me3_overlap\n")
  cat("    - anchor1/2_H3K4me1_overlap\n")
  cat("    - anchor1/2_H3K4me3_overlap\n")
  cat("    - anchor1/2_Bivalent_overlap\n")
  cat("    - anchor1/2_distance_to_tss_ext\n")
  cat("    - anchor1/2_type_extended\n")
  cat("    - loop_type_extended\n\n")

  # --- Step 9: Save annotated data ---
  cat("Step 9: Saving annotated data...\n")

  output_file <- file.path(output_dir, "extended_characterized_loops.tsv")
  write.table(loops_df, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  Saved: %s\n", output_file))

  # Anchor type summary
  anchor_summary <- data.frame(
    anchor = rep(c("Anchor1", "Anchor2"), each = length(ANCHOR_TYPE_ORDER)),
    type = rep(ANCHOR_TYPE_ORDER, 2),
    count = c(
      sapply(ANCHOR_TYPE_ORDER, function(t) sum(anchor1_type == t)),
      sapply(ANCHOR_TYPE_ORDER, function(t) sum(anchor2_type == t))
    )
  )
  anchor_summary$percentage <- 100 * anchor_summary$count / nrow(loops_df)

  anchor_summary_file <- file.path(output_dir, "anchor_type_summary.tsv")
  write.table(anchor_summary, anchor_summary_file, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  Saved: %s\n", anchor_summary_file))

  # Loop type summary
  loop_type_summary <- data.frame(
    loop_type = names(loop_type_counts),
    count = as.numeric(loop_type_counts),
    percentage = 100 * as.numeric(loop_type_counts) / nrow(loops_df)
  )

  loop_type_summary_file <- file.path(output_dir, "loop_type_summary.tsv")
  write.table(loop_type_summary, loop_type_summary_file, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  Saved: %s\n\n", loop_type_summary_file))

  # --- Step 10: Generate visualizations ---
  cat("Step 10: Generating visualizations...\n")

  # 10a. Side-by-side pie charts
  create_loop_type_piechart_comparison(loops_df, output_dir)

  # 10b. Anchor type bar chart
  create_anchor_type_barplot(loops_df, output_dir)

  # 10c. Loop type by direction bar chart
  create_loop_type_direction_barplot(loops_df, output_dir)

  # --- Step 11: Save summary statistics ---
  cat("\nStep 11: Saving summary statistics...\n")

  summary_text <- c(
    "========================================",
    "Extended Loop Anchor Classification Summary",
    "========================================",
    "",
    sprintf("Input file: %s", input_file),
    sprintf("Output dir: %s", output_dir),
    sprintf("TSS threshold: %d bp", tss_threshold),
    "",
    sprintf("Total loops: %d", nrow(loops_df)),
    sprintf("  Up in mutant: %d (%.1f%%)",
            sum(loops_df$direction == "up_in_mutant"),
            100 * mean(loops_df$direction == "up_in_mutant")),
    sprintf("  Down in mutant: %d (%.1f%%)",
            sum(loops_df$direction == "down_in_mutant"),
            100 * mean(loops_df$direction == "down_in_mutant")),
    "",
    "ChIP-seq Peak Files:",
    sprintf("  H3K27ac:  %s (%d peaks)", h3k27ac_path, length(k27ac_gr)),
    sprintf("  H3K27me3: %s (%d peaks)", h3k27me3_path, length(k27me3_gr)),
    sprintf("  H3K4me1:  %s (%d peaks)", h3k4me1_path, length(k4me1_gr)),
    sprintf("  H3K4me3:  %s (%d peaks)", h3k4me3_path, length(k4me3_gr)),
    sprintf("  Bivalent: %s (%d regions)", bivalent_path, length(bivalent_gr)),
    "",
    "Anchor Type Distribution:",
    "  Anchor1:"
  )

  for (type in ANCHOR_TYPE_ORDER) {
    count <- sum(anchor1_type == type)
    pct <- 100 * mean(anchor1_type == type)
    summary_text <- c(summary_text, sprintf("    %-20s: %5d (%.1f%%)", type, count, pct))
  }

  summary_text <- c(summary_text, "  Anchor2:")
  for (type in ANCHOR_TYPE_ORDER) {
    count <- sum(anchor2_type == type)
    pct <- 100 * mean(anchor2_type == type)
    summary_text <- c(summary_text, sprintf("    %-20s: %5d (%.1f%%)", type, count, pct))
  }

  summary_text <- c(summary_text, "", "Loop Type Distribution:")
  for (i in seq_along(loop_type_counts)) {
    type_name <- names(loop_type_counts)[i]
    count <- loop_type_counts[i]
    pct <- 100 * count / nrow(loops_df)
    summary_text <- c(summary_text, sprintf("  %-35s: %5d (%.1f%%)", type_name, count, pct))
  }

  summary_text <- c(summary_text, "", "========================================")

  summary_file <- file.path(output_dir, "summary_statistics.txt")
  writeLines(summary_text, summary_file)
  cat(sprintf("  Saved: %s\n", summary_file))

  cat("\n========================================\n")
  cat("Annotation Complete\n")
  cat("========================================\n\n")

  return(loops_df)
}

# =============================================================================
# VISUALIZATION FUNCTIONS
# =============================================================================

#' Create side-by-side pie charts for up vs down loops
#'
#' @param loops_df Annotated loops data.frame
#' @param output_dir Output directory
create_loop_type_piechart_comparison <- function(loops_df, output_dir) {
  cat("  Creating side-by-side pie charts...\n")

  # Generate colors for loop types
  loop_types <- unique(loops_df$loop_type_extended)
  n_types <- length(loop_types)

  # Use a colorful palette with enough colors
  if (n_types <= 12) {
    loop_type_colors <- RColorBrewer::brewer.pal(max(3, n_types), "Set3")
  } else {
    loop_type_colors <- c(
      RColorBrewer::brewer.pal(12, "Set3"),
      RColorBrewer::brewer.pal(min(n_types - 12, 8), "Pastel1")
    )
  }
  names(loop_type_colors) <- sort(loop_types)[1:length(loop_type_colors)]

  # Summarize by direction
  up_summary <- loops_df %>%
    filter(direction == "up_in_mutant") %>%
    group_by(loop_type_extended) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(percentage = 100 * count / sum(count))

  down_summary <- loops_df %>%
    filter(direction == "down_in_mutant") %>%
    group_by(loop_type_extended) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(percentage = 100 * count / sum(count))

  n_up <- sum(loops_df$direction == "up_in_mutant")
  n_down <- sum(loops_df$direction == "down_in_mutant")

  # Up-regulated pie chart
  p_up <- ggplot(up_summary, aes(x = "", y = percentage, fill = loop_type_extended)) +
    geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.3) +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = loop_type_colors, name = "Loop Type") +
    labs(
      title = "Up-regulated Loops",
      subtitle = sprintf("n = %d", n_up)
    ) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      legend.text = element_text(size = 8),
      legend.key.size = unit(0.4, "cm"),
      legend.position = "none"
    ) +
    geom_text(aes(label = ifelse(percentage > 4, sprintf("%.1f%%", percentage), "")),
              position = position_stack(vjust = 0.5), size = 2.5)

  # Down-regulated pie chart
  p_down <- ggplot(down_summary, aes(x = "", y = percentage, fill = loop_type_extended)) +
    geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.3) +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = loop_type_colors, name = "Loop Type") +
    labs(
      title = "Down-regulated Loops",
      subtitle = sprintf("n = %d", n_down)
    ) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      legend.text = element_text(size = 8),
      legend.key.size = unit(0.4, "cm"),
      legend.position = "none"
    ) +
    geom_text(aes(label = ifelse(percentage > 4, sprintf("%.1f%%", percentage), "")),
              position = position_stack(vjust = 0.5), size = 2.5)

  # Create legend plot
  legend_data <- data.frame(
    loop_type = sort(loop_types),
    y = seq_along(loop_types)
  )

  p_legend <- ggplot(legend_data, aes(x = 1, y = y, fill = loop_type)) +
    geom_tile(width = 0.5, height = 0.8) +
    geom_text(aes(label = loop_type), hjust = 0, x = 1.4, size = 3) +
    scale_fill_manual(values = loop_type_colors) +
    theme_void() +
    theme(legend.position = "none") +
    xlim(0.5, 5) +
    coord_fixed(ratio = 0.3)

  # Combine plots
  p_combined <- (p_up | p_down) / p_legend +
    plot_layout(heights = c(3, 2)) +
    plot_annotation(
      title = "Loop Type Classification by Direction (Extended Categories)",
      subtitle = "BAP1-KO vs Wildtype Differential Loops",
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        plot.subtitle = element_text(hjust = 0.5, size = 12)
      )
    )

  output_file <- file.path(output_dir, "plots", "loop_type_piechart_comparison.pdf")
  ggsave(output_file, p_combined, width = 14, height = 12)
  cat(sprintf("    Saved: %s\n", output_file))
}

#' Create anchor type distribution bar chart
#'
#' @param loops_df Annotated loops data.frame
#' @param output_dir Output directory
create_anchor_type_barplot <- function(loops_df, output_dir) {
  cat("  Creating anchor type bar chart...\n")

  # Reshape data
  anchor_data <- bind_rows(
    loops_df %>%
      dplyr::select(direction, anchor_type = anchor1_type_extended) %>%
      dplyr::mutate(anchor = "Anchor1"),
    loops_df %>%
      dplyr::select(direction, anchor_type = anchor2_type_extended) %>%
      dplyr::mutate(anchor = "Anchor2")
  )

  anchor_summary <- anchor_data %>%
    group_by(direction, anchor, anchor_type) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(direction, anchor) %>%
    mutate(percentage = 100 * count / sum(count))

  # Set factor levels for ordering
  anchor_summary$anchor_type <- factor(anchor_summary$anchor_type, levels = ANCHOR_TYPE_ORDER)
  anchor_summary$direction <- factor(anchor_summary$direction,
                                     levels = c("up_in_mutant", "down_in_mutant"),
                                     labels = c("Up in Mutant", "Down in Mutant"))

  p <- ggplot(anchor_summary, aes(x = anchor, y = percentage, fill = anchor_type)) +
    geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.3) +
    facet_wrap(~direction) +
    scale_fill_manual(values = ANCHOR_COLORS, name = "Anchor Type") +
    labs(
      title = "Anchor Type Distribution by Direction",
      subtitle = "Extended classification with Polycomb and Bivalent categories",
      x = "Anchor",
      y = "Percentage (%)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      legend.position = "right"
    )

  output_file <- file.path(output_dir, "plots", "anchor_type_distribution.pdf")
  ggsave(output_file, p, width = 10, height = 6)
  cat(sprintf("    Saved: %s\n", output_file))
}

#' Create loop type by direction bar chart
#'
#' @param loops_df Annotated loops data.frame
#' @param output_dir Output directory
create_loop_type_direction_barplot <- function(loops_df, output_dir) {
  cat("  Creating loop type by direction bar chart...\n")

  loop_summary <- loops_df %>%
    group_by(direction, loop_type_extended) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(direction) %>%
    mutate(percentage = 100 * count / sum(count))

  loop_summary$direction <- factor(loop_summary$direction,
                                   levels = c("up_in_mutant", "down_in_mutant"),
                                   labels = c("Up in Mutant", "Down in Mutant"))

  # Sort loop types by total count
  loop_type_order <- loops_df %>%
    group_by(loop_type_extended) %>%
    summarise(total = n()) %>%
    arrange(desc(total)) %>%
    pull(loop_type_extended)

  loop_summary$loop_type_extended <- factor(loop_summary$loop_type_extended,
                                            levels = rev(loop_type_order))

  p <- ggplot(loop_summary, aes(x = loop_type_extended, y = count, fill = direction)) +
    geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.3) +
    scale_fill_manual(
      values = c("Up in Mutant" = "#d73027", "Down in Mutant" = "#4575b4"),
      name = "Direction"
    ) +
    labs(
      title = "Loop Type Distribution by Direction",
      subtitle = "Extended classification with Polycomb and Bivalent categories",
      x = "Loop Type",
      y = "Count"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      legend.position = "top"
    ) +
    coord_flip()

  output_file <- file.path(output_dir, "plots", "loop_type_by_direction.pdf")
  ggsave(output_file, p, width = 10, height = 10)
  cat(sprintf("    Saved: %s\n", output_file))
}

# =============================================================================
# COMMAND-LINE ARGUMENT PARSING
# =============================================================================

parse_arguments <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  # Defaults
  input_file <- DEFAULT_INPUT_FILE
  output_dir <- DEFAULT_OUTPUT_DIR

  # Parse arguments
  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--input" && i < length(args)) {
      input_file <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--output" && i < length(args)) {
      output_dir <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--help" || args[i] == "-h") {
      cat("\n")
      cat("Extended Loop Anchor Annotation\n")
      cat("================================\n\n")
      cat("Usage: Rscript scripts/annotate_loops_extended.R [OPTIONS]\n\n")
      cat("Options:\n")
      cat("  --input FILE    Input TSV file (default: 25042-late_outputs/merged_loops/non_redundant_loops.tsv)\n")
      cat("  --output DIR    Output directory (default: outputs/loop_annotation_extended)\n")
      cat("  --help, -h      Show this help message\n\n")
      cat("Description:\n")
      cat("  Annotates differential loops with chromatin state categories:\n")
      cat("    - Active_Promoter, Repressed_Promoter, Bivalent, Polycomb,\n")
      cat("      Active_Enhancer, Poised_Enhancer, Other\n\n")
      cat("Output files:\n")
      cat("  - extended_characterized_loops.tsv  Full annotation table\n")
      cat("  - anchor_type_summary.tsv           Per-anchor statistics\n")
      cat("  - loop_type_summary.tsv             Loop type counts\n")
      cat("  - plots/                            Visualization PDFs\n")
      cat("  - summary_statistics.txt            Text summary\n\n")
      quit(status = 0)
    } else {
      i <- i + 1
    }
  }

  return(list(input_file = input_file, output_dir = output_dir))
}

# =============================================================================
# STANDALONE EXECUTION
# =============================================================================

if (!interactive()) {
  args <- parse_arguments()

  result <- annotate_loops_extended(
    input_file = args$input_file,
    output_dir = args$output_dir
  )

  cat(sprintf("\nOutput saved to: %s\n", args$output_dir))
  cat("Files generated:\n")
  cat("  - extended_characterized_loops.tsv\n")
  cat("  - anchor_type_summary.tsv\n")
  cat("  - loop_type_summary.tsv\n")
  cat("  - plots/loop_type_piechart_comparison.pdf\n")
  cat("  - plots/anchor_type_distribution.pdf\n")
  cat("  - plots/loop_type_by_direction.pdf\n")
  cat("  - summary_statistics.txt\n")
}
