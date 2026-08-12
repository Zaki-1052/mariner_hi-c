#!/usr/bin/env Rscript
# tads/scripts/tad_chip_classification.R
# TAD Boundary ChIP-seq Classification with 8-Category Chromatin State System
# Author: Zakir Alibhai
# Date: 2026-01-19
#
# Purpose:
#   Classify TAD boundaries from TADCompare differential analysis using the
#   8-category chromatin state system (identical to annotate_loops_extended.R).
#   Generates visualizations showing chromatin state distribution by:
#   - Overall boundary distribution
#   - Differential vs Non-Differential status
#   - TAD boundary type (Shifted, Merge, Split, etc.)
#   - Enrichment direction (Control vs Mutant)
#
# Usage:
#   Rscript data/scripts/figure_1_tads_boundaries_compartments/tad_chip_classification.R                    # Default: both
#   Rscript data/scripts/figure_1_tads_boundaries_compartments/tad_chip_classification.R --timepoint late   # Late only
#   Rscript data/scripts/figure_1_tads_boundaries_compartments/tad_chip_classification.R --timepoint early  # Early only
#   Rscript data/scripts/figure_1_tads_boundaries_compartments/tad_chip_classification.R --timepoint both   # Both timepoints
#   Original usage: cd tads/ && Rscript scripts/tad_chip_classification.R
#
# Output:
#   data/plots/figure_3_epigenetic_integration/tad_chip_{early,late}/  (plots)
#   data/tsvs/figure_1_tads_boundaries_compartments/                   (1F TSVs)
#   data/tsvs/figure_3_epigenetic_integration/                         (3B TSVs)
#   Original: results/visualizations/chip/{early,late}/

# =============================================================================
# SECTION 1: SETUP AND CONFIGURATION
# =============================================================================

cat("=== TAD Boundary ChIP-seq Classification ===\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# Load required packages
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(RColorBrewer)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
})

# Load multi-format output utility
# Script runs from tads/ directory
source("data/scripts/_shared/multi_format_output.R")  # Original: source("../scripts/utils/multi_format_output.R")

# =============================================================================
# CONFIGURATION
# =============================================================================

# Input files by timepoint
INPUT_FILES <- list(
  early = "data/tsvs/figure_1_tads_boundaries_compartments/tadcompare_final_annotated_early.tsv",  # Original: results/early/final/tadcompare_final_annotated.tsv
  late  = "data/tsvs/figure_1_tads_boundaries_compartments/tadcompare_final_annotated_late.tsv"    # Original: results/late/final/tadcompare_final_annotated.tsv
)

# ChIP-seq peak files by timepoint
CHIP_PEAK_FILES <- list(
  early = list(
    H3K27ac  = "peaks/beds/H3K27acCerebellumEarly2.bed",       # Note: repo-relative path, not bundled in data/  # Original: ../peaks/beds/H3K27acCerebellumEarly2.bed
    H3K27me3 = "peaks/beds/H3K27me3CerebellumEarly1.bed",      # Note: repo-relative path, not bundled in data/  # Original: ../peaks/beds/H3K27me3CerebellumEarly1.bed
    H3K4me1  = "peaks/beds/H3K4me1CerebellumEarly1.bed",       # Note: repo-relative path, not bundled in data/  # Original: ../peaks/beds/H3K4me1CerebellumEarly1.bed
    H3K4me3  = "peaks/beds/H3K4me3CerebellumEarly2.bed",       # Note: repo-relative path, not bundled in data/  # Original: ../peaks/beds/H3K4me3CerebellumEarly2.bed
    bivalent = "peaks/beds/Bivalent_Cerebellum_Early.bed",      # Note: repo-relative path, not bundled in data/  # Original: ../peaks/beds/Bivalent_Cerebellum_Early.bed
    ctcf     = "peaks/CTCF.bed",                                # Note: repo-relative path, not bundled in data/  # Original: ../peaks/CTCF.bed
    ctcf_motif = "peaks/ctcf_motifs_mm10.bed"                   # Note: repo-relative path, not bundled in data/  # Original: ../peaks/ctcf_motifs_mm10.bed
  ),
  late = list(
    H3K27ac  = "peaks/beds/H3K27acCerebellumLate2.bed",        # Note: repo-relative path, not bundled in data/  # Original: ../peaks/beds/H3K27acCerebellumLate2.bed
    H3K27me3 = "peaks/beds/H3K27me3CerebellumLate1.bed",       # Note: repo-relative path, not bundled in data/  # Original: ../peaks/beds/H3K27me3CerebellumLate1.bed
    H3K4me1  = "peaks/beds/H3K4me1CerebellumLate1.bed",        # Note: repo-relative path, not bundled in data/  # Original: ../peaks/beds/H3K4me1CerebellumLate1.bed
    H3K4me3  = "peaks/beds/H3K4me3CerebellumLate2.bed",        # Note: repo-relative path, not bundled in data/  # Original: ../peaks/beds/H3K4me3CerebellumLate2.bed
    bivalent = "peaks/beds/Bivalent_Cerebellum_Late.bed",       # Note: repo-relative path, not bundled in data/  # Original: ../peaks/beds/Bivalent_Cerebellum_Late.bed
    ctcf     = "peaks/CTCF.bed",                                # Note: repo-relative path, not bundled in data/  # Original: ../peaks/CTCF.bed
    ctcf_motif = "peaks/ctcf_motifs_mm10.bed"                   # Note: repo-relative path, not bundled in data/  # Original: ../peaks/ctcf_motifs_mm10.bed
  )
)

# Output directories by timepoint (split plots vs TSVs)
# Original: results/visualizations/chip/{early,late}
OUTPUT_PLOT_DIRS <- list(
  early = "data/plots/figure_3_epigenetic_integration/tad_chip_early",
  late  = "data/plots/figure_3_epigenetic_integration/tad_chip_late"
)
OUTPUT_TSV_DIRS <- list(
  early = "data/tsvs/figure_1_tads_boundaries_compartments",
  late  = "data/tsvs/figure_1_tads_boundaries_compartments"
)
OUTPUT_TSV_3B_DIRS <- list(
  early = "data/tsvs/figure_3_epigenetic_integration",
  late  = "data/tsvs/figure_3_epigenetic_integration"
)

# 8-category chromatin state system (priority order)
CHROMATIN_STATE_ORDER <- c(
  "Active_Promoter", "Repressed_Promoter", "Bivalent_Promoter",
  "Polycomb", "Active_Enhancer", "Poised_Enhancer",
  "CTCF_Site", "Other"
)

# Color scheme (consistent with annotate_loops_extended.R)
CHROMATIN_STATE_COLORS <- c(
  "Active_Promoter" = "#e41a1c",     # Red - active transcription
  "Repressed_Promoter" = "#756bb1",  # Purple - Polycomb-silenced promoter
  "Bivalent_Promoter" = "#984ea3",   # Magenta - developmental poised
  "Polycomb" = "#4daf4a",            # Green - distal repressive
  "Active_Enhancer" = "#377eb8",     # Blue - active enhancer
  "Poised_Enhancer" = "#ff7f00",     # Orange - primed enhancer
  "CTCF_Site" = "#a65628",           # Brown - structural/insulator
  "Other" = "#999999"                # Gray - unmarked
)

# TAD boundary type colors (from tad_visualizations.R)
BOUNDARY_TYPE_COLORS <- c(
  "Non-Differential" = "#999999",  # Gray
  "Strength Change" = "#ff7f00",   # Orange
  "Complex" = "#984ea3",           # Purple
  "Shifted" = "#e41a1c",           # Red
  "Merge" = "#377eb8",             # Blue
  "Split" = "#4daf4a"              # Green
)

# Enrichment direction colors
ENRICHMENT_COLORS <- c(
  "Matrix 1" = "#4575b4",  # Control = Blue
  "Matrix 2" = "#d73027"   # Mutant = Red
)

# Window size around TAD boundary for ChIP overlap (25kb total)
BOUNDARY_WINDOW <- 12500  # 12.5kb on each side

# TSS threshold for promoter classification
TSS_THRESHOLD <- 2000  # 2kb

# =============================================================================
# SECTION 2: HELPER FUNCTIONS
# =============================================================================

#' Load ChIP-seq peaks from BED file
#'
#' @param bed_path Path to BED file
#' @param mark_name Name of the mark (for logging)
#' @return GRanges object
load_chip_peaks <- function(bed_path, mark_name = "ChIP") {
  if (!file.exists(bed_path)) {
    stop(sprintf("%s BED file not found: %s", mark_name, bed_path))
  }

  df <- read.table(bed_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  gr <- GRanges(
    seqnames = df$V1,
    ranges = IRanges(start = df$V2, end = df$V3)
  )

  cat(sprintf("  %s: %d peaks from %s\n", mark_name, length(gr), basename(bed_path)))
  return(gr)
}

#' Compute ChIP-seq overlaps for boundary GRanges
#'
#' @param boundary_gr GRanges with boundary coordinates (with window)
#' @param chip_peaks_list List of GRanges for each ChIP mark
#' @param ctcf_motif_gr GRanges for CTCF motifs (optional)
#' @return data.frame with overlap columns
compute_chip_overlaps <- function(boundary_gr, chip_peaks_list, ctcf_motif_gr = NULL) {
  overlaps <- data.frame(
    H3K27ac_overlap = countOverlaps(boundary_gr, chip_peaks_list$H3K27ac) > 0,
    H3K27me3_overlap = countOverlaps(boundary_gr, chip_peaks_list$H3K27me3) > 0,
    H3K4me1_overlap = countOverlaps(boundary_gr, chip_peaks_list$H3K4me1) > 0,
    H3K4me3_overlap = countOverlaps(boundary_gr, chip_peaks_list$H3K4me3) > 0,
    bivalent_overlap = countOverlaps(boundary_gr, chip_peaks_list$bivalent) > 0,
    ctcf_overlap = countOverlaps(boundary_gr, chip_peaks_list$ctcf) > 0
  )

  # Add CTCF motif overlap if provided
  if (!is.null(ctcf_motif_gr)) {
    overlaps$ctcf_motif_overlap <- countOverlaps(boundary_gr, ctcf_motif_gr) > 0
  } else {
    overlaps$ctcf_motif_overlap <- FALSE
  }

  return(overlaps)
}

#' Classify boundary chromatin state using 8-category priority system
#'
#' Priority order (identical to annotate_loops_extended.R):
#'   1. Active_Promoter:    H3K4me3+ AND NOT H3K27me3 AND ≤2kb from TSS
#'   2. Repressed_Promoter: H3K27me3+ AND NOT H3K27ac AND ≤2kb from TSS
#'   3. Bivalent_Promoter:  Overlaps pre-computed K4me3+K27me3 intersection
#'   4. Polycomb:           H3K27me3+ AND >2kb from TSS
#'   5. Active_Enhancer:    H3K27ac+ AND >2kb from TSS
#'   6. Poised_Enhancer:    H3K4me1+ AND NOT H3K27ac AND NOT H3K27me3 AND >2kb
#'   7. CTCF_Site:          CTCF motif+ (DNA motifs for all timepoints)
#'   8. Other:              Default (no marks, no CTCF motif)
#'
#' @param overlaps data.frame with overlap columns
#' @param distance_to_tss Numeric vector of distances to nearest TSS
#' @param tss_threshold Distance threshold for promoter (default 2000bp)
#' @param use_motif_for_ctcf Use CTCF motif instead of ChIP (for early timepoint)
#' @return Character vector with chromatin state classifications
classify_boundary_chromatin_state <- function(overlaps, distance_to_tss,
                                               tss_threshold = 2000,
                                               use_motif_for_ctcf = FALSE) {
  n <- nrow(overlaps)
  chromatin_state <- rep("Other", n)

  # Extract overlap columns
  h3k27ac <- overlaps$H3K27ac_overlap
  h3k27me3 <- overlaps$H3K27me3_overlap
  h3k4me1 <- overlaps$H3K4me1_overlap
  h3k4me3 <- overlaps$H3K4me3_overlap
  bivalent <- overlaps$bivalent_overlap
  ctcf <- overlaps$ctcf_overlap
  ctcf_motif <- overlaps$ctcf_motif_overlap

  # 1. Active_Promoter: H3K4me3+ AND NOT H3K27me3 AND ≤2kb from TSS
  is_active_promoter <- h3k4me3 & !h3k27me3 &
    !is.na(distance_to_tss) & distance_to_tss <= tss_threshold
  chromatin_state[is_active_promoter] <- "Active_Promoter"

  # 2. Repressed_Promoter: H3K27me3+ AND NOT H3K27ac AND ≤2kb from TSS
  is_repressed_promoter <- !is_active_promoter &
    h3k27me3 & !h3k27ac &
    !is.na(distance_to_tss) & distance_to_tss <= tss_threshold
  chromatin_state[is_repressed_promoter] <- "Repressed_Promoter"

  # 3. Bivalent_Promoter: K4me3+K27me3 overlap (not already classified)
  is_bivalent <- !is_active_promoter & !is_repressed_promoter & bivalent
  chromatin_state[is_bivalent] <- "Bivalent_Promoter"

  # 4. Polycomb: H3K27me3+ AND >2kb from TSS
  is_polycomb <- !is_active_promoter & !is_repressed_promoter & !is_bivalent &
    h3k27me3 & (is.na(distance_to_tss) | distance_to_tss > tss_threshold)
  chromatin_state[is_polycomb] <- "Polycomb"

  # 5. Active_Enhancer: H3K27ac+ AND >2kb from TSS
  is_active_enhancer <- !is_active_promoter & !is_repressed_promoter &
    !is_bivalent & !is_polycomb &
    h3k27ac & (is.na(distance_to_tss) | distance_to_tss > tss_threshold)
  chromatin_state[is_active_enhancer] <- "Active_Enhancer"

  # 6. Poised_Enhancer: H3K4me1+ AND NOT H3K27ac AND NOT H3K27me3 AND >2kb
  is_poised_enhancer <- !is_active_promoter & !is_repressed_promoter &
    !is_bivalent & !is_polycomb & !is_active_enhancer &
    h3k4me1 & !h3k27ac & !h3k27me3 &
    (is.na(distance_to_tss) | distance_to_tss > tss_threshold)
  chromatin_state[is_poised_enhancer] <- "Poised_Enhancer"

  # 7. CTCF_Site: Use motif for early, ChIP for late
  if (use_motif_for_ctcf) {
    is_ctcf_site <- !is_active_promoter & !is_repressed_promoter &
      !is_bivalent & !is_polycomb & !is_active_enhancer &
      !is_poised_enhancer & ctcf_motif
  } else {
    is_ctcf_site <- !is_active_promoter & !is_repressed_promoter &
      !is_bivalent & !is_polycomb & !is_active_enhancer &
      !is_poised_enhancer & ctcf
  }
  chromatin_state[is_ctcf_site] <- "CTCF_Site"

  # 8. Other: Default (no marks)
  return(chromatin_state)
}

# =============================================================================
# SECTION 3: MAIN ANALYSIS FUNCTION
# =============================================================================

#' Run TAD boundary ChIP-seq classification for a single timepoint
#'
#' @param timepoint "early" or "late"
#' @return List with analysis results
run_tad_chip_classification <- function(timepoint) {
  cat("\n")
  cat("============================================================\n")
  cat(sprintf("Processing %s timepoint\n", toupper(timepoint)))
  cat("============================================================\n\n")

  # Set paths for this timepoint
  input_file <- INPUT_FILES[[timepoint]]
  chip_files <- CHIP_PEAK_FILES[[timepoint]]
  PLOT_DIR <- OUTPUT_PLOT_DIRS[[timepoint]]   # Original: OUTPUT_DIRS[[timepoint]]
  TSV_DIR  <- OUTPUT_TSV_DIRS[[timepoint]]    # Original: OUTPUT_DIRS[[timepoint]]
  TSV_3B_DIR <- OUTPUT_TSV_3B_DIRS[[timepoint]]  # Original: OUTPUT_DIRS[[timepoint]]
  dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)   # Original: dir.create(OUTPUT_DIR, ...)
  dir.create(TSV_DIR, showWarnings = FALSE, recursive = TRUE)    # Original: (was same dir)
  dir.create(TSV_3B_DIR, showWarnings = FALSE, recursive = TRUE) # Original: (was same dir)

  cat("Input file:", input_file, "\n")
  cat("Output plot directory:", PLOT_DIR, "\n")
  cat("Output TSV directory (1F):", TSV_DIR, "\n")
  cat("Output TSV directory (3B):", TSV_3B_DIR, "\n\n")

  # Validate input file exists
  if (!file.exists(input_file)) {
    stop(sprintf("ERROR: Input file not found: %s", input_file))
  }

  # ==========================================================================
  # Step 1: Load TAD Boundaries
  # ==========================================================================

  cat("=== Step 1: Loading TAD Boundary Data ===\n")

  tad_df <- read.table(input_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  cat("Loaded", nrow(tad_df), "TAD boundaries\n")

  # Summarize by Differential status
  cat("  Differential:", sum(tad_df$Differential == "Differential", na.rm = TRUE), "\n")
  cat("  Non-Differential:", sum(tad_df$Differential == "Non-Differential", na.rm = TRUE), "\n")

  # Summarize by Type
  cat("\n  Boundary types:\n")
  type_counts <- table(tad_df$Type)
  for (t in names(type_counts)) {
    cat(sprintf("    %s: %d (%.1f%%)\n", t, type_counts[t], 100 * type_counts[t] / nrow(tad_df)))
  }
  cat("\n")

  # ==========================================================================
  # Step 2: Create GRanges for Boundaries
  # ==========================================================================

  cat("=== Step 2: Creating Boundary GRanges ===\n")

  # Use 25kb window (12.5kb each side) around boundary position
  boundary_gr <- GRanges(
    seqnames = tad_df$chr,
    ranges = IRanges(
      start = pmax(1, tad_df$Boundary - BOUNDARY_WINDOW),
      end = tad_df$Boundary + BOUNDARY_WINDOW
    )
  )

  cat(sprintf("  Created GRanges for %d boundaries\n", length(boundary_gr)))
  cat(sprintf("  Window size: %d bp on each side (%d kb total)\n\n",
              BOUNDARY_WINDOW, BOUNDARY_WINDOW * 2 / 1000))

  # ==========================================================================
  # Step 3: Load ChIP-seq Peaks
  # ==========================================================================

  cat("=== Step 3: Loading ChIP-seq Peak Files ===\n")

  chip_peaks <- list(
    H3K27ac = load_chip_peaks(chip_files$H3K27ac, "H3K27ac"),
    H3K27me3 = load_chip_peaks(chip_files$H3K27me3, "H3K27me3"),
    H3K4me1 = load_chip_peaks(chip_files$H3K4me1, "H3K4me1"),
    H3K4me3 = load_chip_peaks(chip_files$H3K4me3, "H3K4me3"),
    bivalent = load_chip_peaks(chip_files$bivalent, "Bivalent"),
    ctcf = load_chip_peaks(chip_files$ctcf, "CTCF")
  )

  # Load CTCF motifs
  ctcf_motif_gr <- NULL
  if (!is.null(chip_files$ctcf_motif) && file.exists(chip_files$ctcf_motif)) {
    ctcf_motif_gr <- load_chip_peaks(chip_files$ctcf_motif, "CTCF_motif")
  }
  cat("\n")

  # ==========================================================================
  # Step 4: Compute TSS Distances
  # ==========================================================================

  cat("=== Step 4: Computing TSS Distances ===\n")

  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  genes <- genes(txdb)
  tss_gr <- resize(genes, width = 1, fix = "start")
  cat(sprintf("  Loaded %d genes for TSS distance computation\n", length(genes)))

  nearest <- distanceToNearest(boundary_gr, tss_gr)
  distance_to_tss <- rep(NA_real_, length(boundary_gr))
  if (length(nearest) > 0) {
    distance_to_tss[queryHits(nearest)] <- mcols(nearest)$distance
  }

  # Summary of TSS distances
  cat(sprintf("  Boundaries <= 2kb from TSS: %d (%.1f%%)\n",
              sum(distance_to_tss <= TSS_THRESHOLD, na.rm = TRUE),
              100 * mean(distance_to_tss <= TSS_THRESHOLD, na.rm = TRUE)))
  cat(sprintf("  Boundaries > 2kb from TSS: %d (%.1f%%)\n\n",
              sum(distance_to_tss > TSS_THRESHOLD, na.rm = TRUE),
              100 * mean(distance_to_tss > TSS_THRESHOLD, na.rm = TRUE)))

  # ==========================================================================
  # Step 5: Compute ChIP-seq Overlaps
  # ==========================================================================

  cat("=== Step 5: Computing ChIP-seq Overlaps ===\n")

  overlaps <- compute_chip_overlaps(boundary_gr, chip_peaks, ctcf_motif_gr)

  # Print overlap summary
  cat("  Mark overlap summary:\n")
  cat(sprintf("    H3K27ac:  %d boundaries (%.1f%%)\n",
              sum(overlaps$H3K27ac_overlap),
              100 * mean(overlaps$H3K27ac_overlap)))
  cat(sprintf("    H3K27me3: %d boundaries (%.1f%%)\n",
              sum(overlaps$H3K27me3_overlap),
              100 * mean(overlaps$H3K27me3_overlap)))
  cat(sprintf("    H3K4me1:  %d boundaries (%.1f%%)\n",
              sum(overlaps$H3K4me1_overlap),
              100 * mean(overlaps$H3K4me1_overlap)))
  cat(sprintf("    H3K4me3:  %d boundaries (%.1f%%)\n",
              sum(overlaps$H3K4me3_overlap),
              100 * mean(overlaps$H3K4me3_overlap)))
  cat(sprintf("    Bivalent: %d boundaries (%.1f%%)\n",
              sum(overlaps$bivalent_overlap),
              100 * mean(overlaps$bivalent_overlap)))
  cat(sprintf("    CTCF:     %d boundaries (%.1f%%)\n",
              sum(overlaps$ctcf_overlap),
              100 * mean(overlaps$ctcf_overlap)))
  cat(sprintf("    CTCF_motif: %d boundaries (%.1f%%)\n\n",
              sum(overlaps$ctcf_motif_overlap),
              100 * mean(overlaps$ctcf_motif_overlap)))

  # ==========================================================================
  # Step 6: Classify Chromatin States
  # ==========================================================================

  cat("=== Step 6: Classifying Chromatin States (8 categories) ===\n")

  # Use CTCF DNA motifs for ALL timepoints for methodological consistency
  # DNA motifs are genome-wide sequence features (not timepoint-specific)
  # ChIP-seq overlap columns are still computed for reference/validation
  use_motif_for_ctcf <- TRUE
  cat("  CTCF classification: Using DNA motifs (consistent across timepoints)\n")

  chromatin_state <- classify_boundary_chromatin_state(
    overlaps, distance_to_tss,
    tss_threshold = TSS_THRESHOLD,
    use_motif_for_ctcf = use_motif_for_ctcf
  )

  # Print chromatin state distribution
  cat("\n  Chromatin state distribution:\n")
  for (state in CHROMATIN_STATE_ORDER) {
    count <- sum(chromatin_state == state)
    pct <- 100 * count / length(chromatin_state)
    cat(sprintf("    %-20s: %5d (%.1f%%)\n", state, count, pct))
  }
  cat("\n")

  # ==========================================================================
  # Step 7: Add Classifications to Data Frame
  # ==========================================================================

  cat("=== Step 7: Adding Classifications to Data Frame ===\n")

  # Add all overlap columns
  tad_df$H3K27ac_overlap <- overlaps$H3K27ac_overlap
  tad_df$H3K27me3_overlap <- overlaps$H3K27me3_overlap
  tad_df$H3K4me1_overlap <- overlaps$H3K4me1_overlap
  tad_df$H3K4me3_overlap <- overlaps$H3K4me3_overlap
  tad_df$bivalent_overlap <- overlaps$bivalent_overlap
  tad_df$CTCF_overlap <- overlaps$ctcf_overlap
  tad_df$CTCF_motif_overlap <- overlaps$ctcf_motif_overlap

  # Add TSS distance and chromatin state
  tad_df$distance_to_tss <- distance_to_tss
  tad_df$chromatin_state <- chromatin_state

  # Set factor levels for proper ordering
  tad_df$chromatin_state <- factor(tad_df$chromatin_state, levels = CHROMATIN_STATE_ORDER)

  cat("  Added columns: H3K27ac_overlap, H3K27me3_overlap, H3K4me1_overlap,\n")
  cat("                 H3K4me3_overlap, bivalent_overlap, CTCF_overlap,\n")
  cat("                 CTCF_motif_overlap, distance_to_tss, chromatin_state\n\n")

  # ==========================================================================
  # Step 8: Generate Visualizations
  # ==========================================================================

  cat("=== Step 8: Generating Visualizations ===\n")

  # Figure 1: Overall Chromatin State Distribution
  cat("  Creating Figure 1: Overall distribution...\n")
  create_chromatin_state_distribution(tad_df, PLOT_DIR, timepoint)  # Original: OUTPUT_DIR

  # Figure 2: By Differential Status
  cat("  Creating Figure 2: By differential status...\n")
  create_chromatin_by_differential(tad_df, PLOT_DIR, timepoint)  # Original: OUTPUT_DIR

  # Figure 3: By TAD Boundary Type
  cat("  Creating Figure 3: By boundary type...\n")
  create_chromatin_by_boundary_type(tad_df, PLOT_DIR, timepoint)  # Original: OUTPUT_DIR

  # Figure 4: By Enrichment Direction
  cat("  Creating Figure 4: By enrichment direction...\n")
  create_chromatin_by_enrichment(tad_df, PLOT_DIR, timepoint)  # Original: OUTPUT_DIR

  # Figure 5: ChIP Mark Overlap Heatmap
  cat("  Creating Figure 5: ChIP mark heatmap...\n")
  create_chip_overlap_heatmap(tad_df, PLOT_DIR, timepoint)  # Original: OUTPUT_DIR

  # Figure 6: Combined Summary
  cat("  Creating Figure 6: Combined summary...\n")
  create_combined_summary(tad_df, PLOT_DIR, timepoint)  # Original: OUTPUT_DIR

  cat("\n")

  # ==========================================================================
  # Step 9: Export Results
  # ==========================================================================

  cat("=== Step 9: Exporting Results ===\n")

  # Save full annotated data (1F TSV)
  annotated_file <- file.path(TSV_DIR, sprintf("boundaries_with_chromatin_state_%s.tsv", timepoint))  # Original: file.path(OUTPUT_DIR, "boundaries_with_chromatin_state.tsv")
  write.table(tad_df, annotated_file, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  Saved: %s\n", annotated_file))

  # Save chromatin state summary (3B TSV)
  state_summary <- tad_df %>%
    group_by(chromatin_state) %>%
    summarise(
      count = n(),
      percentage = 100 * n() / nrow(tad_df),
      n_differential = sum(Differential == "Differential", na.rm = TRUE),
      n_shifted = sum(Type == "Shifted", na.rm = TRUE),
      .groups = "drop"
    )

  summary_file <- file.path(TSV_3B_DIR, sprintf("boundary_chromatin_state_summary_%s.tsv", timepoint))  # Original: file.path(OUTPUT_DIR, "boundary_chromatin_state_summary.tsv")
  write.table(state_summary, summary_file, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  Saved: %s\n", summary_file))

  # Save statistics report (3B TSV)
  stats_file <- file.path(TSV_3B_DIR, sprintf("chip_classification_statistics_%s.txt", timepoint))  # Original: file.path(OUTPUT_DIR, "chip_classification_statistics.txt")
  create_statistics_report(tad_df, timepoint, chip_files, stats_file)
  cat(sprintf("  Saved: %s\n", stats_file))

  # ==========================================================================
  # Return Results
  # ==========================================================================

  cat("\n=== Analysis Complete for", toupper(timepoint), "===\n")
  cat("Plot directory:", PLOT_DIR, "\n")
  cat("TSV directory (1F):", TSV_DIR, "\n")
  cat("TSV directory (3B):", TSV_3B_DIR, "\n")
  cat("Files generated:\n")
  cat("  Plots:\n")
  cat("  - 01_boundary_chromatin_state_distribution.{pdf,svg,jpg}\n")
  cat("  - 02_chromatin_state_by_differential.{pdf,svg,jpg}\n")
  cat("  - 03_chromatin_state_by_boundary_type.{pdf,svg,jpg}\n")
  cat("  - 04_chromatin_state_by_enrichment.{pdf,svg,jpg}\n")
  cat("  - 05_chip_mark_overlap_heatmap.{pdf,svg,jpg}\n")
  cat("  - 06_summary_comparison.{pdf,svg,jpg}\n")
  cat("  TSVs:\n")
  cat(sprintf("  - boundaries_with_chromatin_state_%s.tsv (1F)\n", timepoint))
  cat(sprintf("  - boundary_chromatin_state_summary_%s.tsv (3B)\n", timepoint))
  cat(sprintf("  - chip_classification_statistics_%s.txt (3B)\n\n", timepoint))

  return(list(
    timepoint = timepoint,
    n_boundaries = nrow(tad_df),
    state_distribution = table(chromatin_state),
    tad_df = tad_df
  ))
}

# =============================================================================
# SECTION 4: VISUALIZATION FUNCTIONS
# =============================================================================

#' Create overall chromatin state distribution plot
create_chromatin_state_distribution <- function(tad_df, output_dir, timepoint) {
  # Summary data
  state_summary <- tad_df %>%
    group_by(chromatin_state) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(percentage = 100 * count / sum(count))

  # Bar chart
  p_bar <- ggplot(state_summary, aes(x = chromatin_state, y = percentage, fill = chromatin_state)) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.1f%%", percentage)), vjust = -0.3, size = 3) +
    scale_fill_manual(values = CHROMATIN_STATE_COLORS) +
    scale_y_continuous(limits = c(0, max(state_summary$percentage) * 1.1), expand = c(0, 0)) +
    labs(
      title = "TAD Boundary Chromatin State Distribution",
      subtitle = sprintf("%s timepoint | n = %d boundaries", toupper(timepoint), nrow(tad_df)),
      x = "Chromatin State",
      y = "Percentage (%)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )

  # Pie chart
  p_pie <- ggplot(state_summary, aes(x = "", y = percentage, fill = chromatin_state)) +
    geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.3) +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = CHROMATIN_STATE_COLORS, name = "Chromatin State") +
    labs(title = "Distribution") +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      legend.position = "right"
    ) +
    geom_text(aes(label = ifelse(percentage > 4, sprintf("%.1f%%", percentage), "")),
              position = position_stack(vjust = 0.5), size = 2.5)

  # Combine
  p_combined <- p_bar + p_pie +
    plot_layout(widths = c(2, 1)) +
    plot_annotation(
      title = "TAD Boundary Chromatin State Classification",
      subtitle = sprintf("%s timepoint | 8-category system", toupper(timepoint)),
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey40")
      )
    )

  save_multiformat_ggplot(p_combined,
                          file.path(output_dir, "01_boundary_chromatin_state_distribution"),
                          width = 14, height = 6)
}

#' Create chromatin state by differential status plot
create_chromatin_by_differential <- function(tad_df, output_dir, timepoint) {
  # Summary by differential status
  diff_summary <- tad_df %>%
    group_by(Differential, chromatin_state) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(Differential) %>%
    mutate(percentage = 100 * count / sum(count)) %>%
    ungroup()

  n_diff <- sum(tad_df$Differential == "Differential", na.rm = TRUE)
  n_nondiff <- sum(tad_df$Differential == "Non-Differential", na.rm = TRUE)

  # Stacked bar chart
  p_stacked <- ggplot(diff_summary, aes(x = Differential, y = percentage, fill = chromatin_state)) +
    geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.2) +
    scale_fill_manual(values = CHROMATIN_STATE_COLORS, name = "Chromatin State") +
    scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
    labs(
      title = "Chromatin State by Differential Status",
      subtitle = sprintf("Differential: n=%d | Non-Differential: n=%d", n_diff, n_nondiff),
      x = "Differential Status",
      y = "Percentage (%)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"),
      legend.position = "right"
    )

  # Grouped bar chart
  p_grouped <- ggplot(diff_summary, aes(x = chromatin_state, y = percentage, fill = Differential)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8),
             width = 0.7, color = "black", linewidth = 0.2) +
    scale_fill_manual(values = c("Differential" = "#d73027", "Non-Differential" = "#4575b4"),
                      name = "Status") +
    scale_y_continuous(limits = c(0, max(diff_summary$percentage) * 1.1), expand = c(0, 0)) +
    labs(
      title = "Chromatin State Comparison",
      subtitle = "Differential vs Non-Differential boundaries",
      x = "Chromatin State",
      y = "Percentage (%)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey40"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top"
    )

  # Combine
  p_combined <- p_stacked / p_grouped +
    plot_annotation(
      title = "Chromatin State Distribution by Differential Status",
      subtitle = sprintf("%s timepoint", toupper(timepoint)),
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey40")
      )
    )

  save_multiformat_ggplot(p_combined,
                          file.path(output_dir, "02_chromatin_state_by_differential"),
                          width = 12, height = 10)
}

#' Create chromatin state by boundary type plot
create_chromatin_by_boundary_type <- function(tad_df, output_dir, timepoint) {
  # Filter to get boundary types (exclude Non-Differential for cleaner viz)
  type_summary <- tad_df %>%
    group_by(Type, chromatin_state) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(Type) %>%
    mutate(percentage = 100 * count / sum(count)) %>%
    ungroup()

  # Get type counts for labels
  type_counts <- tad_df %>%
    group_by(Type) %>%
    summarise(n = n(), .groups = "drop")

  # Stacked bar chart - all types
  p_all <- ggplot(type_summary, aes(x = Type, y = percentage, fill = chromatin_state)) +
    geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.2) +
    scale_fill_manual(values = CHROMATIN_STATE_COLORS, name = "Chromatin State") +
    scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
    labs(
      title = "Chromatin State by TAD Boundary Type",
      subtitle = "All boundary types",
      x = "Boundary Type",
      y = "Percentage (%)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )

  # Focus on differential types only
  diff_types <- c("Shifted", "Merge", "Split", "Strength Change", "Complex")
  type_summary_diff <- type_summary %>%
    filter(Type %in% diff_types)

  p_diff <- ggplot(type_summary_diff, aes(x = Type, y = percentage, fill = chromatin_state)) +
    geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.2) +
    scale_fill_manual(values = CHROMATIN_STATE_COLORS, name = "Chromatin State") +
    scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
    labs(
      title = "Differential Boundary Types Only",
      subtitle = "Excluding Non-Differential",
      x = "Boundary Type",
      y = "Percentage (%)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey40"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )

  # Combine
  p_combined <- p_all + p_diff +
    plot_layout(widths = c(1.5, 1)) +
    plot_annotation(
      title = "Chromatin State Classification by TAD Boundary Type",
      subtitle = sprintf("%s timepoint", toupper(timepoint)),
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey40")
      )
    )

  save_multiformat_ggplot(p_combined,
                          file.path(output_dir, "03_chromatin_state_by_boundary_type"),
                          width = 14, height = 7)
}

#' Create chromatin state by enrichment direction plot
create_chromatin_by_enrichment <- function(tad_df, output_dir, timepoint) {
  # Filter to differential boundaries only
  diff_df <- tad_df %>%
    filter(Differential == "Differential" | Type != "Non-Differential") %>%
    filter(!is.na(Enriched_In))

  if (nrow(diff_df) == 0) {
    cat("    No differential boundaries with enrichment info. Skipping.\n")
    return(NULL)
  }

  # Summary by enrichment
  enrich_summary <- diff_df %>%
    group_by(Enriched_In, chromatin_state) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(Enriched_In) %>%
    mutate(percentage = 100 * count / sum(count)) %>%
    ungroup()

  # Add labels
  enrich_summary$Enriched_Label <- case_when(
    enrich_summary$Enriched_In == "Matrix 1" ~ "Control-enriched",
    enrich_summary$Enriched_In == "Matrix 2" ~ "Mutant-enriched",
    TRUE ~ enrich_summary$Enriched_In
  )

  n_ctrl <- sum(diff_df$Enriched_In == "Matrix 1")
  n_mut <- sum(diff_df$Enriched_In == "Matrix 2")

  # Stacked bar chart
  p_stacked <- ggplot(enrich_summary, aes(x = Enriched_Label, y = percentage, fill = chromatin_state)) +
    geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.2) +
    scale_fill_manual(values = CHROMATIN_STATE_COLORS, name = "Chromatin State") +
    scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
    labs(
      title = "Chromatin State by Enrichment Direction",
      subtitle = sprintf("Control-enriched: n=%d | Mutant-enriched: n=%d", n_ctrl, n_mut),
      x = "Enrichment Direction",
      y = "Percentage (%)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"),
      legend.position = "right"
    )

  # Grouped comparison
  p_grouped <- ggplot(enrich_summary, aes(x = chromatin_state, y = percentage, fill = Enriched_Label)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8),
             width = 0.7, color = "black", linewidth = 0.2) +
    scale_fill_manual(values = c("Control-enriched" = "#4575b4", "Mutant-enriched" = "#d73027"),
                      name = "Direction") +
    scale_y_continuous(limits = c(0, max(enrich_summary$percentage) * 1.1), expand = c(0, 0)) +
    labs(
      title = "Chromatin State: Control vs Mutant",
      subtitle = "Differential boundaries only",
      x = "Chromatin State",
      y = "Percentage (%)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey40"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top"
    )

  # Combine
  p_combined <- p_stacked / p_grouped +
    plot_annotation(
      title = "Chromatin State by Enrichment Direction",
      subtitle = sprintf("%s timepoint | Differential boundaries only", toupper(timepoint)),
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey40")
      )
    )

  save_multiformat_ggplot(p_combined,
                          file.path(output_dir, "04_chromatin_state_by_enrichment"),
                          width = 11, height = 10)
}

#' Create ChIP mark overlap heatmap
create_chip_overlap_heatmap <- function(tad_df, output_dir, timepoint) {
  # Calculate overlap percentages by boundary type
  mark_cols <- c("H3K27ac_overlap", "H3K27me3_overlap", "H3K4me1_overlap",
                 "H3K4me3_overlap", "CTCF_overlap")

  type_overlap <- tad_df %>%
    group_by(Type) %>%
    summarise(
      H3K27ac = 100 * mean(H3K27ac_overlap, na.rm = TRUE),
      H3K27me3 = 100 * mean(H3K27me3_overlap, na.rm = TRUE),
      H3K4me1 = 100 * mean(H3K4me1_overlap, na.rm = TRUE),
      H3K4me3 = 100 * mean(H3K4me3_overlap, na.rm = TRUE),
      CTCF = 100 * mean(CTCF_overlap, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(cols = -Type, names_to = "Mark", values_to = "Percentage")

  # By differential status
  diff_overlap <- tad_df %>%
    group_by(Differential) %>%
    summarise(
      H3K27ac = 100 * mean(H3K27ac_overlap, na.rm = TRUE),
      H3K27me3 = 100 * mean(H3K27me3_overlap, na.rm = TRUE),
      H3K4me1 = 100 * mean(H3K4me1_overlap, na.rm = TRUE),
      H3K4me3 = 100 * mean(H3K4me3_overlap, na.rm = TRUE),
      CTCF = 100 * mean(CTCF_overlap, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(cols = -Differential, names_to = "Mark", values_to = "Percentage")

  # Heatmap by boundary type
  p_type_heatmap <- ggplot(type_overlap, aes(x = Mark, y = Type, fill = Percentage)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.1f%%", Percentage)), size = 3) +
    scale_fill_gradient2(low = "white", mid = "#fee090", high = "#d73027",
                         midpoint = 50, name = "% Overlap", limits = c(0, 100)) +
    labs(
      title = "ChIP Mark Overlap by Boundary Type",
      x = "Histone Mark / CTCF",
      y = "Boundary Type"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      axis.text.x = element_text(size = 10, face = "bold"),
      axis.text.y = element_text(size = 10),
      panel.grid = element_blank()
    )

  # Heatmap by differential status
  p_diff_heatmap <- ggplot(diff_overlap, aes(x = Mark, y = Differential, fill = Percentage)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.1f%%", Percentage)), size = 3.5) +
    scale_fill_gradient2(low = "white", mid = "#fee090", high = "#d73027",
                         midpoint = 50, name = "% Overlap", limits = c(0, 100)) +
    labs(
      title = "ChIP Mark Overlap by Differential Status",
      x = "Histone Mark / CTCF",
      y = "Status"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      axis.text.x = element_text(size = 10, face = "bold"),
      axis.text.y = element_text(size = 10),
      panel.grid = element_blank()
    )

  # Combine
  p_combined <- p_type_heatmap / p_diff_heatmap +
    plot_layout(heights = c(2, 1)) +
    plot_annotation(
      title = "ChIP-seq Mark Overlap at TAD Boundaries",
      subtitle = sprintf("%s timepoint | %% of boundaries overlapping each mark", toupper(timepoint)),
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey40")
      )
    )

  save_multiformat_ggplot(p_combined,
                          file.path(output_dir, "05_chip_mark_overlap_heatmap"),
                          width = 10, height = 10)
}

#' Create combined summary figure
create_combined_summary <- function(tad_df, output_dir, timepoint) {
  # Panel A: Overall distribution bar
  state_summary <- tad_df %>%
    group_by(chromatin_state) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(percentage = 100 * count / sum(count))

  p_a <- ggplot(state_summary, aes(x = chromatin_state, y = percentage, fill = chromatin_state)) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.2) +
    scale_fill_manual(values = CHROMATIN_STATE_COLORS) +
    scale_y_continuous(limits = c(0, max(state_summary$percentage) * 1.1), expand = c(0, 0)) +
    labs(title = "A. Overall Distribution", x = "", y = "% Boundaries") +
    theme_minimal(base_size = 9) +
    theme(
      plot.title = element_text(face = "bold", size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      legend.position = "none"
    )

  # Panel B: Differential vs Non-Differential
  diff_summary <- tad_df %>%
    group_by(Differential, chromatin_state) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(Differential) %>%
    mutate(percentage = 100 * count / sum(count))

  p_b <- ggplot(diff_summary, aes(x = Differential, y = percentage, fill = chromatin_state)) +
    geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.2) +
    scale_fill_manual(values = CHROMATIN_STATE_COLORS) +
    scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
    labs(title = "B. By Differential Status", x = "", y = "% Boundaries") +
    theme_minimal(base_size = 9) +
    theme(
      plot.title = element_text(face = "bold", size = 10),
      axis.text.x = element_text(size = 8),
      legend.position = "none"
    )

  # Panel C: By enrichment (differential only)
  diff_df <- tad_df %>%
    filter(Differential == "Differential" | Type != "Non-Differential") %>%
    filter(!is.na(Enriched_In))

  if (nrow(diff_df) > 0) {
    enrich_summary <- diff_df %>%
      group_by(Enriched_In, chromatin_state) %>%
      summarise(count = n(), .groups = "drop") %>%
      group_by(Enriched_In) %>%
      mutate(percentage = 100 * count / sum(count)) %>%
      mutate(Enriched_Label = ifelse(Enriched_In == "Matrix 1", "Control", "Mutant"))

    p_c <- ggplot(enrich_summary, aes(x = Enriched_Label, y = percentage, fill = chromatin_state)) +
      geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.2) +
      scale_fill_manual(values = CHROMATIN_STATE_COLORS) +
      scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
      labs(title = "C. By Enrichment", x = "", y = "% Boundaries") +
      theme_minimal(base_size = 9) +
      theme(
        plot.title = element_text(face = "bold", size = 10),
        axis.text.x = element_text(size = 8),
        legend.position = "none"
      )
  } else {
    p_c <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "No data", size = 4) +
      labs(title = "C. By Enrichment") +
      theme_void() +
      theme(plot.title = element_text(face = "bold", size = 10))
  }

  # Panel D: Legend
  legend_data <- data.frame(
    state = factor(CHROMATIN_STATE_ORDER, levels = CHROMATIN_STATE_ORDER),
    y = seq_along(CHROMATIN_STATE_ORDER)
  )

  p_legend <- ggplot(legend_data, aes(x = 1, y = y, fill = state)) +
    geom_tile(width = 0.3, height = 0.8) +
    geom_text(aes(label = state), hjust = 0, x = 1.25, size = 2.5) +
    scale_fill_manual(values = CHROMATIN_STATE_COLORS) +
    theme_void() +
    theme(legend.position = "none") +
    labs(title = "D. Legend") +
    theme(plot.title = element_text(face = "bold", size = 10, hjust = 0)) +
    xlim(0.5, 3) +
    coord_fixed(ratio = 0.3)

  # Combine panels
  p_combined <- (p_a | p_b | p_c | p_legend) +
    plot_annotation(
      title = "TAD Boundary Chromatin State Classification Summary",
      subtitle = sprintf("%s timepoint | n = %d boundaries | 8-category system",
                        toupper(timepoint), nrow(tad_df)),
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey40")
      )
    )

  save_multiformat_ggplot(p_combined,
                          file.path(output_dir, "06_summary_comparison"),
                          width = 14, height = 5)
}

#' Create statistics report
create_statistics_report <- function(tad_df, timepoint, chip_files, output_file) {
  sink(output_file)

  cat("===========================================\n")
  cat("TAD BOUNDARY ChIP-seq CLASSIFICATION REPORT\n")
  cat("===========================================\n")
  cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Timepoint:", toupper(timepoint), "\n\n")

  cat("--- DATA OVERVIEW ---\n")
  cat("Total boundaries:", nrow(tad_df), "\n")
  cat("  Differential:", sum(tad_df$Differential == "Differential", na.rm = TRUE), "\n")
  cat("  Non-Differential:", sum(tad_df$Differential == "Non-Differential", na.rm = TRUE), "\n\n")

  cat("--- BOUNDARY TYPES ---\n")
  type_tbl <- table(tad_df$Type)
  for (t in names(type_tbl)) {
    cat(sprintf("  %s: %d (%.1f%%)\n", t, type_tbl[t], 100 * type_tbl[t] / nrow(tad_df)))
  }
  cat("\n")

  cat("--- ChIP-seq PEAK FILES ---\n")
  cat("  H3K27ac:", chip_files$H3K27ac, "\n")
  cat("  H3K27me3:", chip_files$H3K27me3, "\n")
  cat("  H3K4me1:", chip_files$H3K4me1, "\n")
  cat("  H3K4me3:", chip_files$H3K4me3, "\n")
  cat("  Bivalent:", chip_files$bivalent, "\n")
  cat("  CTCF:", chip_files$ctcf, "\n")
  cat("  CTCF_motif:", chip_files$ctcf_motif, "\n\n")

  cat("--- ChIP MARK OVERLAP ---\n")
  cat(sprintf("  H3K27ac:  %.1f%% of boundaries\n", 100 * mean(tad_df$H3K27ac_overlap)))
  cat(sprintf("  H3K27me3: %.1f%% of boundaries\n", 100 * mean(tad_df$H3K27me3_overlap)))
  cat(sprintf("  H3K4me1:  %.1f%% of boundaries\n", 100 * mean(tad_df$H3K4me1_overlap)))
  cat(sprintf("  H3K4me3:  %.1f%% of boundaries\n", 100 * mean(tad_df$H3K4me3_overlap)))
  cat(sprintf("  Bivalent: %.1f%% of boundaries\n", 100 * mean(tad_df$bivalent_overlap)))
  cat(sprintf("  CTCF:     %.1f%% of boundaries\n", 100 * mean(tad_df$CTCF_overlap)))
  cat(sprintf("  CTCF_motif: %.1f%% of boundaries\n\n", 100 * mean(tad_df$CTCF_motif_overlap)))

  cat("--- CHROMATIN STATE DISTRIBUTION ---\n")
  state_tbl <- table(tad_df$chromatin_state)
  for (s in CHROMATIN_STATE_ORDER) {
    count <- ifelse(s %in% names(state_tbl), state_tbl[s], 0)
    cat(sprintf("  %-20s: %5d (%.1f%%)\n", s, count, 100 * count / nrow(tad_df)))
  }
  cat("\n")

  cat("--- TSS DISTANCE SUMMARY ---\n")
  tss_dist <- tad_df$distance_to_tss
  cat(sprintf("  Median: %.1f bp\n", median(tss_dist, na.rm = TRUE)))
  cat(sprintf("  Mean: %.1f bp\n", mean(tss_dist, na.rm = TRUE)))
  cat(sprintf("  <= 2kb from TSS: %d (%.1f%%)\n",
              sum(tss_dist <= 2000, na.rm = TRUE),
              100 * mean(tss_dist <= 2000, na.rm = TRUE)))
  cat(sprintf("  > 2kb from TSS: %d (%.1f%%)\n\n",
              sum(tss_dist > 2000, na.rm = TRUE),
              100 * mean(tss_dist > 2000, na.rm = TRUE)))

  cat("--- CTCF CLASSIFICATION NOTE ---\n")
  cat("  Used CTCF DNA MOTIFS for classification (consistent across all timepoints)\n")
  cat("  DNA motifs are genome-wide sequence features, not timepoint-specific\n")
  cat("  CTCF ChIP-seq overlap columns are saved for reference only\n")
  cat("\n")

  cat("===========================================\n")

  sink()
}

# =============================================================================
# SECTION 5: ARGUMENT PARSING
# =============================================================================

parse_arguments <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  timepoint <- "both"  # Default

  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--timepoint" && i < length(args)) {
      timepoint <- args[i + 1]
      i <- i + 2
    } else if (args[i] %in% c("--help", "-h")) {
      cat("Usage: Rscript scripts/tad_chip_classification.R [OPTIONS]\n\n")
      cat("Options:\n")
      cat("  --timepoint TP  Timepoint: 'early', 'late', or 'both' (default: both)\n")
      cat("  --help, -h      Show this help message\n\n")
      cat("Output:\n")
      cat("  data/plots/figure_3_epigenetic_integration/tad_chip_{early,late}/\n")  # Original: results/visualizations/chip/{early,late}/
      cat("  data/tsvs/figure_1_tads_boundaries_compartments/\n")
      cat("  data/tsvs/figure_3_epigenetic_integration/\n\n")
      cat("8-Category Chromatin State System:\n")
      cat("  1. Active_Promoter    - H3K4me3+ AND NOT H3K27me3 AND <=2kb TSS\n")
      cat("  2. Repressed_Promoter - H3K27me3+ AND NOT H3K27ac AND <=2kb TSS\n")
      cat("  3. Bivalent_Promoter  - K4me3+K27me3 overlap\n")
      cat("  4. Polycomb           - H3K27me3+ AND >2kb TSS\n")
      cat("  5. Active_Enhancer    - H3K27ac+ AND >2kb TSS\n")
      cat("  6. Poised_Enhancer    - H3K4me1+ only AND >2kb TSS\n")
      cat("  7. CTCF_Site          - CTCF motif+ (DNA motifs for all timepoints)\n")
      cat("  8. Other              - No marks, no CTCF motif\n")
      quit(save = "no", status = 0)
    } else {
      i <- i + 1
    }
  }

  if (!timepoint %in% c("early", "late", "both")) {
    stop("ERROR: timepoint must be 'early', 'late', or 'both'")
  }

  return(timepoint)
}

# =============================================================================
# SECTION 6: MAIN EXECUTION
# =============================================================================

TIMEPOINT_ARG <- parse_arguments()

if (TIMEPOINT_ARG == "both") {
  TIMEPOINTS_TO_RUN <- c("early", "late")
} else {
  TIMEPOINTS_TO_RUN <- TIMEPOINT_ARG
}

cat("Timepoint(s) to process:", paste(TIMEPOINTS_TO_RUN, collapse = ", "), "\n\n")

results <- list()

for (tp in TIMEPOINTS_TO_RUN) {
  tryCatch({
    results[[tp]] <- run_tad_chip_classification(tp)
  }, error = function(e) {
    cat(sprintf("\nERROR processing %s timepoint: %s\n", tp, e$message))
  })
}

# Print summary across timepoints
cat("\n=== All Timepoints Complete ===\n")
cat("Output directories:\n")
for (tp in TIMEPOINTS_TO_RUN) {
  cat(sprintf("  - Plots: %s\n", OUTPUT_PLOT_DIRS[[tp]]))   # Original: OUTPUT_DIRS[[tp]]
  cat(sprintf("  - TSVs (1F): %s\n", OUTPUT_TSV_DIRS[[tp]]))
  cat(sprintf("  - TSVs (3B): %s\n", OUTPUT_TSV_3B_DIRS[[tp]]))
}

if (length(results) > 0) {
  cat("\nChromatin State Summary:\n")
  for (tp in names(results)) {
    r <- results[[tp]]
    cat(sprintf("\n%s timepoint (%d boundaries):\n", toupper(tp), r$n_boundaries))
    for (state in CHROMATIN_STATE_ORDER) {
      count <- ifelse(state %in% names(r$state_distribution),
                      r$state_distribution[state], 0)
      cat(sprintf("  %-20s: %5d (%.1f%%)\n", state, count,
                  100 * count / r$n_boundaries))
    }
  }
}

cat("\nEnd time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
