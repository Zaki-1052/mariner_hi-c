# scripts/chip_distance_analysis.R
# ChIP-seq Mark Distance Analysis: Testing H3K27ac vs H3K27me3 Trends
#
# Hypothesis: As loop distance increases, H3K27ac decreases and H3K27me3 increases
# This extends Figure 7 from loop_distance_analysis.R by adding H3K27me3 and using
# timepoint-specific ChIP-seq peaks.
#
# Usage:
#   Rscript scripts/chip_distance_analysis.R                    # Default: late
#   Rscript scripts/chip_distance_analysis.R --timepoint late   # Late timepoint
#   Rscript scripts/chip_distance_analysis.R --timepoint early  # Early timepoint
#   Rscript scripts/chip_distance_analysis.R --timepoint both   # Run both
#
# Output: output/chip_distance/{early,late}/

# ==============================================================================
# SECTION 1: SETUP AND CONFIGURATION
# ==============================================================================

cat("=== ChIP-seq Mark Distance Analysis ===\n")
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
source("scripts/utils/multi_format_output.R")

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Input files by timepoint
INPUT_FILES <- list(
  late = "25042-late_outputs/merged_loops/non_redundant_loops.tsv",
  early = "250831-early_outputs/merged_loops/non_redundant_loops.tsv"
)

# ChIP-seq peak files by timepoint (relative paths from project root)
CHIP_PEAK_FILES <- list(
  early = list(
    H3K27ac  = "peaks/beds/H3K27acCerebellumEarly2.bed",
    H3K27me3 = "peaks/beds/H3K27me3CerebellumEarly1.bed",
    H3K4me1  = "peaks/beds/H3K4me1CerebellumEarly1.bed",
    H3K4me3  = "peaks/beds/H3K4me3CerebellumEarly2.bed",
    bivalent = "peaks/beds/Bivalent_Cerebellum_Early.bed",
    ctcf     = "peaks/CTCF.bed",
    ctcf_motif = "peaks/ctcf_motifs_mm10.bed"
  ),
  late = list(
    H3K27ac  = "peaks/beds/H3K27acCerebellumLate2.bed",
    H3K27me3 = "peaks/beds/H3K27me3CerebellumLate1.bed",
    H3K4me1  = "peaks/beds/H3K4me1CerebellumLate1.bed",
    H3K4me3  = "peaks/beds/H3K4me3CerebellumLate2.bed",
    bivalent = "peaks/beds/Bivalent_Cerebellum_Late.bed",
    ctcf     = "peaks/CTCF.bed",
    ctcf_motif = "peaks/ctcf_motifs_mm10.bed"
  )
)

# Anchor type order and colors (from annotate_loops_extended.R)
ANCHOR_TYPE_ORDER <- c("Active_Promoter", "Repressed_Promoter", "Bivalent_Promoter",
                       "Polycomb", "Active_Enhancer", "Poised_Enhancer",
                       "CTCF_Site", "Other")

ANCHOR_COLORS <- c(
  "Active_Promoter" = "#e41a1c",
  "Repressed_Promoter" = "#756bb1",
  "Bivalent_Promoter" = "#984ea3",
  "Polycomb" = "#4daf4a",
  "Active_Enhancer" = "#377eb8",
  "Poised_Enhancer" = "#ff7f00",
  "CTCF_Site" = "#a65628",
  "Other" = "#999999"
)

# Output directories by timepoint
OUTPUT_DIRS <- list(
  late = "output/chip_distance/late",
  early = "output/chip_distance/early"
)

# Distance category order and thresholds
DISTANCE_ORDER <- c("<100kb", "100-500kb", "500kb-1Mb", ">1Mb")
DISTANCE_BREAKS <- c(0, 100000, 500000, 1000000, Inf)

# Color palettes (consistent with existing pipeline)
COLORS <- list(
  down = "#d73027",      # Red for down/lost in mutant
  up = "#4575b4",        # Blue for up/gained in mutant
  neutral = "#999999"
)

# Direction labels
DIRECTION_LABELS <- c(
  "down_in_mutant" = "Lost in BAP1-KO",
  "up_in_mutant" = "Gained in BAP1-KO"
)

# ==============================================================================
# ARGUMENT PARSING
# ==============================================================================

parse_arguments <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  timepoint <- "late"  # Default

  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--timepoint" && i < length(args)) {
      timepoint <- args[i + 1]
      i <- i + 2
    } else if (args[i] %in% c("--help", "-h")) {
      cat("Usage: Rscript scripts/chip_distance_analysis.R [OPTIONS]\n\n")
      cat("Options:\n")
      cat("  --timepoint TP  Timepoint: 'early', 'late', or 'both' (default: late)\n")
      cat("  --help, -h      Show this help message\n\n")
      cat("Hypothesis: As loop distance increases, H3K27ac decreases and H3K27me3 increases\n")
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

TIMEPOINT_ARG <- parse_arguments()

if (TIMEPOINT_ARG == "both") {
  TIMEPOINTS_TO_RUN <- c("late", "early")
} else {
  TIMEPOINTS_TO_RUN <- TIMEPOINT_ARG
}

cat("Timepoint(s) to process:", paste(TIMEPOINTS_TO_RUN, collapse = ", "), "\n\n")

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

#' Load ChIP-seq peaks from BED file
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

#' Assign distance category based on loop distance
assign_distance_category <- function(distance) {
  cut(distance,
      breaks = DISTANCE_BREAKS,
      labels = DISTANCE_ORDER,
      include.lowest = TRUE)
}

#' Compute ChIP-seq overlaps for anchor GRanges (extended for 8-category system)
compute_chip_overlaps <- function(anchor_gr, chip_peaks_list, ctcf_motif_gr = NULL) {
  overlaps <- data.frame(
    H3K27ac_overlap = countOverlaps(anchor_gr, chip_peaks_list$H3K27ac) > 0,
    H3K27me3_overlap = countOverlaps(anchor_gr, chip_peaks_list$H3K27me3) > 0,
    H3K4me1_overlap = countOverlaps(anchor_gr, chip_peaks_list$H3K4me1) > 0,
    H3K4me3_overlap = countOverlaps(anchor_gr, chip_peaks_list$H3K4me3) > 0,
    bivalent_overlap = countOverlaps(anchor_gr, chip_peaks_list$bivalent) > 0,
    ctcf_overlap = countOverlaps(anchor_gr, chip_peaks_list$ctcf) > 0
  )
  # Add CTCF motif overlap if provided
  if (!is.null(ctcf_motif_gr)) {
    overlaps$ctcf_motif_overlap <- countOverlaps(anchor_gr, ctcf_motif_gr) > 0
  } else {
    overlaps$ctcf_motif_overlap <- FALSE
  }
  return(overlaps)
}

#' Classify anchor type with 8-category chromatin state system
#' (Priority order from annotate_loops_extended.R)
#'
#' Priority order:
#'   1. Active_Promoter:    H3K4me3+ AND NOT H3K27me3 AND ≤2kb from TSS
#'   2. Repressed_Promoter: H3K27me3+ AND NOT H3K27ac AND ≤2kb from TSS
#'   3. Bivalent_Promoter:  Overlaps pre-computed K4me3+K27me3 intersection
#'   4. Polycomb:           H3K27me3+ AND >2kb from TSS
#'   5. Active_Enhancer:    H3K27ac+ AND >2kb from TSS
#'   6. Poised_Enhancer:    H3K4me1+ AND NOT H3K27ac AND NOT H3K27me3 AND >2kb
#'   7. CTCF_Site:          CTCF+ (ChIP for late, motif for early)
#'   8. Other:              Default (no marks)
#'
#' @param overlaps data.frame with overlap columns
#' @param distance_to_tss Numeric vector of distances to nearest TSS
#' @param tss_threshold Distance threshold for promoter (default 2000bp)
#' @param use_motif_for_ctcf Use CTCF motif instead of ChIP (for early timepoint)
#' @return Character vector with anchor types
classify_anchor_type_extended <- function(overlaps, distance_to_tss,
                                          tss_threshold = 2000,
                                          use_motif_for_ctcf = FALSE) {
  n <- nrow(overlaps)
  anchor_type <- rep("Other", n)

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
  anchor_type[is_active_promoter] <- "Active_Promoter"

  # 2. Repressed_Promoter: H3K27me3+ AND NOT H3K27ac AND ≤2kb from TSS
  is_repressed_promoter <- !is_active_promoter &
    h3k27me3 & !h3k27ac &
    !is.na(distance_to_tss) & distance_to_tss <= tss_threshold
  anchor_type[is_repressed_promoter] <- "Repressed_Promoter"

  # 3. Bivalent_Promoter: K4me3+K27me3 overlap (not already classified)
  is_bivalent <- !is_active_promoter & !is_repressed_promoter & bivalent
  anchor_type[is_bivalent] <- "Bivalent_Promoter"

  # 4. Polycomb: H3K27me3+ AND >2kb from TSS
  is_polycomb <- !is_active_promoter & !is_repressed_promoter & !is_bivalent &
    h3k27me3 & (is.na(distance_to_tss) | distance_to_tss > tss_threshold)
  anchor_type[is_polycomb] <- "Polycomb"

  # 5. Active_Enhancer: H3K27ac+ AND >2kb from TSS
  is_active_enhancer <- !is_active_promoter & !is_repressed_promoter &
    !is_bivalent & !is_polycomb &
    h3k27ac & (is.na(distance_to_tss) | distance_to_tss > tss_threshold)
  anchor_type[is_active_enhancer] <- "Active_Enhancer"

  # 6. Poised_Enhancer: H3K4me1+ AND NOT H3K27ac AND NOT H3K27me3 AND >2kb
  is_poised_enhancer <- !is_active_promoter & !is_repressed_promoter &
    !is_bivalent & !is_polycomb & !is_active_enhancer &
    h3k4me1 & !h3k27ac & !h3k27me3 &
    (is.na(distance_to_tss) | distance_to_tss > tss_threshold)
  anchor_type[is_poised_enhancer] <- "Poised_Enhancer"

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
  anchor_type[is_ctcf_site] <- "CTCF_Site"

  # 8. Other: Default (no marks)
  return(anchor_type)
}

# ==============================================================================
# MAIN ANALYSIS FUNCTION
# ==============================================================================

run_chip_distance_analysis <- function(timepoint) {
  cat("\n")
  cat("============================================================\n")
  cat(sprintf("Processing %s timepoint\n", toupper(timepoint)))
  cat("============================================================\n\n")

  # Set paths for this timepoint
  input_file <- INPUT_FILES[[timepoint]]
  chip_files <- CHIP_PEAK_FILES[[timepoint]]
  OUTPUT_DIR <- OUTPUT_DIRS[[timepoint]]
  dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

  cat("Input file:", input_file, "\n")
  cat("Output directory:", OUTPUT_DIR, "\n\n")

  # Validate input file exists
  if (!file.exists(input_file)) {
    stop(sprintf("ERROR: Input file not found: %s", input_file))
  }

  # ==========================================================================
  # SECTION 2: DATA LOADING
  # ==========================================================================

  cat("=== Step 1: Loading Loop Data ===\n")

  loops <- read_tsv(input_file, show_col_types = FALSE)
  cat("Loaded", nrow(loops), "differential loops\n")

  # Check required columns
  required_cols <- c("chr1", "start1", "end1", "chr2", "start2", "end2", "logFC", "direction")
  missing_cols <- setdiff(required_cols, colnames(loops))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }

  # Calculate loop distance if not present
  if (!"loop_distance" %in% colnames(loops)) {
    loops <- loops %>%
      mutate(
        mid1 = (start1 + end1) / 2,
        mid2 = (start2 + end2) / 2,
        loop_distance = abs(mid2 - mid1)
      )
  }

  # Add distance category and direction labels
  loops <- loops %>%
    mutate(
      distance_category = factor(assign_distance_category(loop_distance), levels = DISTANCE_ORDER),
      direction_label = factor(
        case_when(
          direction == "down_in_mutant" ~ "Lost in BAP1-KO",
          direction == "up_in_mutant" ~ "Gained in BAP1-KO",
          TRUE ~ "Other"
        ),
        levels = c("Lost in BAP1-KO", "Gained in BAP1-KO")
      )
    )

  # Filter to directional loops only
  loops_directional <- loops %>%
    filter(direction %in% c("up_in_mutant", "down_in_mutant"))

  # Separate by direction for chi-square tests
  lost_loops <- loops_directional %>% filter(direction == "down_in_mutant")
  gained_loops <- loops_directional %>% filter(direction == "up_in_mutant")

  cat("Directional loops:", nrow(loops_directional), "\n")
  cat("  - Lost (down_in_mutant):", nrow(lost_loops), "\n")
  cat("  - Gained (up_in_mutant):", nrow(gained_loops), "\n\n")

  # ==========================================================================
  # SECTION 3: LOAD ChIP-seq PEAKS
  # ==========================================================================

  cat("=== Step 2: Loading ChIP-seq Peak Files ===\n")

  chip_peaks <- list(
    H3K27ac = load_chip_peaks(chip_files$H3K27ac, "H3K27ac"),
    H3K27me3 = load_chip_peaks(chip_files$H3K27me3, "H3K27me3"),
    H3K4me1 = load_chip_peaks(chip_files$H3K4me1, "H3K4me1"),
    H3K4me3 = load_chip_peaks(chip_files$H3K4me3, "H3K4me3"),
    bivalent = load_chip_peaks(chip_files$bivalent, "Bivalent (K4me3+K27me3)"),
    ctcf = load_chip_peaks(chip_files$ctcf, "CTCF")
  )

  # Load CTCF motifs if available (for early timepoint validation)
  ctcf_motif_gr <- NULL
  if (!is.null(chip_files$ctcf_motif) && file.exists(chip_files$ctcf_motif)) {
    ctcf_motif_gr <- load_chip_peaks(chip_files$ctcf_motif, "CTCF_motif")
  }
  cat("\n")

  # ==========================================================================
  # SECTION 4: COMPUTE ChIP-seq OVERLAPS
  # ==========================================================================

  cat("=== Step 3: Computing ChIP-seq Overlaps ===\n")

  # Create GRanges for anchors
  anchor1_gr <- GRanges(
    seqnames = loops_directional$chr1,
    ranges = IRanges(start = loops_directional$start1, end = loops_directional$end1)
  )

  anchor2_gr <- GRanges(
    seqnames = loops_directional$chr2,
    ranges = IRanges(start = loops_directional$start2, end = loops_directional$end2)
  )

  # Compute overlaps for both anchors (extended for 8-category system)
  anchor1_overlaps <- compute_chip_overlaps(anchor1_gr, chip_peaks, ctcf_motif_gr)
  anchor2_overlaps <- compute_chip_overlaps(anchor2_gr, chip_peaks, ctcf_motif_gr)

  # Add overlap columns to dataframe
  loops_directional <- loops_directional %>%
    mutate(
      # Anchor 1 overlaps
      anchor1_H3K27ac_overlap = anchor1_overlaps$H3K27ac_overlap,
      anchor1_H3K27me3_overlap = anchor1_overlaps$H3K27me3_overlap,
      anchor1_H3K4me1_overlap = anchor1_overlaps$H3K4me1_overlap,
      anchor1_H3K4me3_overlap = anchor1_overlaps$H3K4me3_overlap,
      # Anchor 2 overlaps
      anchor2_H3K27ac_overlap = anchor2_overlaps$H3K27ac_overlap,
      anchor2_H3K27me3_overlap = anchor2_overlaps$H3K27me3_overlap,
      anchor2_H3K4me1_overlap = anchor2_overlaps$H3K4me1_overlap,
      anchor2_H3K4me3_overlap = anchor2_overlaps$H3K4me3_overlap,
      # Loop-level: at least one anchor has the mark
      has_H3K27ac = anchor1_H3K27ac_overlap | anchor2_H3K27ac_overlap,
      has_H3K27me3 = anchor1_H3K27me3_overlap | anchor2_H3K27me3_overlap,
      has_H3K4me1 = anchor1_H3K4me1_overlap | anchor2_H3K4me1_overlap,
      has_H3K4me3 = anchor1_H3K4me3_overlap | anchor2_H3K4me3_overlap
    )

  # Print overlap summary
  cat("Loop-level mark overlap summary:\n")
  cat(sprintf("  H3K27ac:  %d loops (%.1f%%)\n",
              sum(loops_directional$has_H3K27ac),
              100 * mean(loops_directional$has_H3K27ac)))
  cat(sprintf("  H3K27me3: %d loops (%.1f%%)\n",
              sum(loops_directional$has_H3K27me3),
              100 * mean(loops_directional$has_H3K27me3)))
  cat(sprintf("  H3K4me1:  %d loops (%.1f%%)\n",
              sum(loops_directional$has_H3K4me1),
              100 * mean(loops_directional$has_H3K4me1)))
  cat(sprintf("  H3K4me3:  %d loops (%.1f%%)\n\n",
              sum(loops_directional$has_H3K4me3),
              100 * mean(loops_directional$has_H3K4me3)))

  # ==========================================================================
  # SECTION 4b: TSS DISTANCE AND ANCHOR TYPE CLASSIFICATION
  # ==========================================================================

  cat("=== Step 3b: Computing TSS Distances and Anchor Types ===\n")

  # Load TxDb for TSS annotations
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  genes <- genes(txdb)
  tss_gr <- resize(genes, width = 1, fix = "start")
  cat(sprintf("  Loaded %d genes for TSS distance computation\n", length(genes)))

  # Compute TSS distance for each anchor
  nearest1 <- distanceToNearest(anchor1_gr, tss_gr)
  anchor1_distance_to_tss <- rep(NA_real_, length(anchor1_gr))
  if (length(nearest1) > 0) {
    anchor1_distance_to_tss[queryHits(nearest1)] <- mcols(nearest1)$distance
  }

  nearest2 <- distanceToNearest(anchor2_gr, tss_gr)
  anchor2_distance_to_tss <- rep(NA_real_, length(anchor2_gr))
  if (length(nearest2) > 0) {
    anchor2_distance_to_tss[queryHits(nearest2)] <- mcols(nearest2)$distance
  }

  # Classify anchor types using 8-category system
  # Early timepoints use CTCF motif, late timepoints use CTCF ChIP
  use_motif_for_ctcf <- timepoint == "early"
  if (use_motif_for_ctcf) {
    cat("  CTCF classification: Using motif (early timepoint)\n")
  } else {
    cat("  CTCF classification: Using ChIP-seq (late timepoint)\n")
  }

  anchor1_type <- classify_anchor_type_extended(
    anchor1_overlaps, anchor1_distance_to_tss,
    tss_threshold = 2000, use_motif_for_ctcf = use_motif_for_ctcf
  )

  anchor2_type <- classify_anchor_type_extended(
    anchor2_overlaps, anchor2_distance_to_tss,
    tss_threshold = 2000, use_motif_for_ctcf = use_motif_for_ctcf
  )

  # Add anchor types to dataframe
  loops_directional <- loops_directional %>%
    mutate(
      anchor1_distance_to_tss = anchor1_distance_to_tss,
      anchor2_distance_to_tss = anchor2_distance_to_tss,
      anchor1_type = anchor1_type,
      anchor2_type = anchor2_type
    )

  # Print anchor type distribution
  cat("\n  Anchor type distribution (8-category system):\n")
  for (type in ANCHOR_TYPE_ORDER) {
    n1 <- sum(anchor1_type == type)
    n2 <- sum(anchor2_type == type)
    pct1 <- 100 * n1 / length(anchor1_type)
    pct2 <- 100 * n2 / length(anchor2_type)
    cat(sprintf("    %-20s: Anchor1 %5d (%.1f%%), Anchor2 %5d (%.1f%%)\n",
                type, n1, pct1, n2, pct2))
  }
  cat("\n")

  # ==========================================================================
  # SECTION 5: STATISTICAL TESTS
  # ==========================================================================

  cat("=== Step 4: Statistical Tests ===\n")

  # Chi-square tests for mark vs distance (by direction)
  # For lost loops
  chisq_lost_h3k27ac <- chisq.test(table(lost_loops$distance_category,
                                          loops_directional$has_H3K27ac[loops_directional$direction == "down_in_mutant"]))
  chisq_lost_h3k27me3 <- chisq.test(table(lost_loops$distance_category,
                                           loops_directional$has_H3K27me3[loops_directional$direction == "down_in_mutant"]))

  # For gained loops
  chisq_gained_h3k27ac <- chisq.test(table(gained_loops$distance_category,
                                            loops_directional$has_H3K27ac[loops_directional$direction == "up_in_mutant"]))
  chisq_gained_h3k27me3 <- chisq.test(table(gained_loops$distance_category,
                                             loops_directional$has_H3K27me3[loops_directional$direction == "up_in_mutant"]))

  cat("Chi-square tests (mark vs distance):\n")
  cat(sprintf("  Lost loops - H3K27ac:  p = %.2e\n", chisq_lost_h3k27ac$p.value))
  cat(sprintf("  Lost loops - H3K27me3: p = %.2e\n", chisq_lost_h3k27me3$p.value))
  cat(sprintf("  Gained loops - H3K27ac:  p = %.2e\n", chisq_gained_h3k27ac$p.value))
  cat(sprintf("  Gained loops - H3K27me3: p = %.2e\n\n", chisq_gained_h3k27me3$p.value))

  # Trend tests using logistic regression
  cat("Trend tests (logistic regression):\n")

  glm_h3k27ac <- glm(has_H3K27ac ~ loop_distance, data = loops_directional, family = binomial)
  glm_h3k27me3 <- glm(has_H3K27me3 ~ loop_distance, data = loops_directional, family = binomial)

  coef_h3k27ac <- summary(glm_h3k27ac)$coefficients["loop_distance", ]
  coef_h3k27me3 <- summary(glm_h3k27me3)$coefficients["loop_distance", ]

  cat(sprintf("  H3K27ac:  beta = %.2e, p = %.2e\n",
              coef_h3k27ac["Estimate"], coef_h3k27ac["Pr(>|z|)"]))
  cat(sprintf("  H3K27me3: beta = %.2e, p = %.2e\n\n",
              coef_h3k27me3["Estimate"], coef_h3k27me3["Pr(>|z|)"]))

  # ==========================================================================
  # FIGURE 1: ChIP-seq Mark x Distance x Direction Analysis
  # (Same style as Figure 7 from loop_distance_analysis.R)
  # ==========================================================================

  cat("=== Step 5: Generating Figures ===\n")
  cat("Creating Figure 1: ChIP-seq x Distance Analysis (matching Fig 7 style)...\n")

  # Summary of ChIP-seq overlaps by distance and direction
  # Now including H3K27me3!
  chipseq_summary <- loops_directional %>%
    group_by(direction_label, distance_category) %>%
    summarise(
      pct_H3K27ac = mean(has_H3K27ac) * 100,
      pct_H3K27me3 = mean(has_H3K27me3) * 100,
      pct_H3K4me1 = mean(has_H3K4me1) * 100,
      pct_H3K4me3 = mean(has_H3K4me3) * 100,
      .groups = "drop"
    ) %>%
    pivot_longer(cols = c(pct_H3K27ac, pct_H3K27me3, pct_H3K4me1, pct_H3K4me3),
                 names_to = "mark", values_to = "percentage") %>%
    mutate(mark = gsub("pct_", "", mark))

  # Panel A: Focus on H3K27ac vs H3K27me3 (primary hypothesis)
  chipseq_focus <- chipseq_summary %>%
    filter(mark %in% c("H3K27ac", "H3K27me3"))

  p1a <- ggplot(chipseq_focus,
                aes(x = distance_category, y = percentage, fill = mark)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8),
             width = 0.7, color = "black", linewidth = 0.2) +
    facet_wrap(~direction_label) +
    scale_fill_manual(values = c("H3K27ac" = "#e41a1c", "H3K27me3" = "#377eb8"),
                      name = "Histone Mark") +
    scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
    labs(
      title = "H3K27ac vs H3K27me3 Overlap by Distance",
      x = "Distance Category",
      y = "% Loops with Mark at Anchor"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )

  # Panel B: All four marks for complete picture
  p1b <- ggplot(chipseq_summary,
                aes(x = distance_category, y = percentage, fill = mark)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8),
             width = 0.7, color = "black", linewidth = 0.2) +
    facet_wrap(~direction_label) +
    scale_fill_manual(values = c("H3K27ac" = "#e41a1c", "H3K27me3" = "#377eb8",
                                 "H3K4me1" = "#4daf4a", "H3K4me3" = "#984ea3"),
                      name = "Histone Mark") +
    scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
    labs(
      title = "All ChIP-seq Marks Overlap by Distance",
      x = "Distance Category",
      y = "% Loops with Mark at Anchor"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )

  # Combine panels
  p1_combined <- p1a / p1b +
    plot_annotation(
      title = "ChIP-seq Features by Loop Distance and Direction",
      subtitle = sprintf("H3K27ac trend: p=%.2e, H3K27me3 trend: p=%.2e (%s timepoint)",
                        coef_h3k27ac["Pr(>|z|)"], coef_h3k27me3["Pr(>|z|)"], timepoint),
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey40")
      )
    )

  save_multiformat_ggplot(p1_combined, file.path(OUTPUT_DIR, "01_chipseq_distance_analysis"),
                          width = 11, height = 10)

  # ==========================================================================
  # FIGURE 2: Line Plot - Mark Trends by Distance (Key Hypothesis Visualization)
  # ==========================================================================

  cat("Creating Figure 2: Mark Trend Line Plot...\n")

  # Aggregate across directions for overall trend
  distance_summary <- loops_directional %>%
    group_by(distance_category) %>%
    summarise(
      n_loops = n(),
      pct_H3K27ac = mean(has_H3K27ac) * 100,
      pct_H3K27me3 = mean(has_H3K27me3) * 100,
      pct_H3K4me1 = mean(has_H3K4me1) * 100,
      pct_H3K4me3 = mean(has_H3K4me3) * 100,
      .groups = "drop"
    )

  distance_long <- distance_summary %>%
    pivot_longer(cols = starts_with("pct_"),
                 names_to = "mark", values_to = "percentage") %>%
    mutate(mark = gsub("pct_", "", mark))

  # Focus on H3K27ac vs H3K27me3
  distance_focus <- distance_long %>%
    filter(mark %in% c("H3K27ac", "H3K27me3"))

  p2_line <- ggplot(distance_focus,
                    aes(x = distance_category, y = percentage, color = mark, group = mark)) +
    geom_line(linewidth = 1.5) +
    geom_point(size = 4) +
    scale_color_manual(
      values = c("H3K27ac" = "#e41a1c", "H3K27me3" = "#377eb8"),
      name = "Histone Mark",
      labels = c("H3K27ac (Active)", "H3K27me3 (Repressive)")
    ) +
    scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
    labs(
      title = "ChIP-seq Mark Enrichment Trend by Loop Distance",
      subtitle = sprintf("Hypothesis: H3K27ac decreases, H3K27me3 increases with distance\n(%s timepoint, n = %d loops)",
                        timepoint, nrow(loops_directional)),
      x = "Loop Distance Category",
      y = "% Loops with Mark at Anchor"
    ) +
    annotate("text", x = 3.5, y = 90,
             label = sprintf("H3K27ac trend: p = %.2e\nH3K27me3 trend: p = %.2e",
                            coef_h3k27ac["Pr(>|z|)"], coef_h3k27me3["Pr(>|z|)"]),
             hjust = 0.5, size = 3.5, fontface = "italic") +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"),
      legend.position = "top",
      axis.text.x = element_text(size = 11),
      panel.grid.minor = element_blank()
    )

  save_multiformat_ggplot(p2_line, file.path(OUTPUT_DIR, "02_mark_trend_lineplot"),
                          width = 9, height = 7)

  # ==========================================================================
  # FIGURE 3: Heatmap of Mark Enrichment by Distance
  # ==========================================================================

  cat("Creating Figure 3: Mark Distance Heatmap...\n")

  heatmap_data <- distance_long %>%
    mutate(
      mark = factor(mark, levels = c("H3K27ac", "H3K27me3", "H3K4me1", "H3K4me3")),
      distance_category = factor(distance_category, levels = DISTANCE_ORDER)
    )

  p3_heatmap <- ggplot(heatmap_data,
                       aes(x = mark, y = distance_category, fill = percentage)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.1f%%", percentage)), size = 4) +
    scale_fill_gradient2(
      low = "white", mid = "#fee090", high = "#d73027",
      midpoint = 50, name = "% Loops",
      limits = c(0, 100)
    ) +
    labs(
      title = "ChIP-seq Mark Enrichment Heatmap",
      subtitle = sprintf("%s timepoint: %% of loops with mark at one or more anchors", timepoint),
      x = "Histone Mark",
      y = "Loop Distance Category"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"),
      axis.text.x = element_text(size = 11, face = "bold"),
      axis.text.y = element_text(size = 11),
      panel.grid = element_blank()
    )

  save_multiformat_ggplot(p3_heatmap, file.path(OUTPUT_DIR, "03_mark_distance_heatmap"),
                          width = 8, height = 6)

  # ==========================================================================
  # FIGURE 4: Combined Summary Figure
  # ==========================================================================

  cat("Creating Figure 4: Combined Summary Figure...\n")

  # Panel A: Line plot (simplified)
  p4a <- ggplot(distance_focus,
                aes(x = distance_category, y = percentage, color = mark, group = mark)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    scale_color_manual(
      values = c("H3K27ac" = "#e41a1c", "H3K27me3" = "#377eb8"),
      name = ""
    ) +
    scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
    labs(title = "A. Mark Trend by Distance", x = "", y = "% Loops") +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      legend.position = "top",
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
    )

  # Panel B: Barplot for Lost loops
  chipseq_lost <- chipseq_focus %>% filter(direction_label == "Lost in BAP1-KO")

  p4b <- ggplot(chipseq_lost,
                aes(x = distance_category, y = percentage, fill = mark)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8),
             width = 0.7, color = "black", linewidth = 0.2) +
    scale_fill_manual(values = c("H3K27ac" = "#e41a1c", "H3K27me3" = "#377eb8"), name = "") +
    scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
    labs(title = "B. Lost Loops", x = "", y = "% Loops") +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
    )

  # Panel C: Barplot for Gained loops
  chipseq_gained <- chipseq_focus %>% filter(direction_label == "Gained in BAP1-KO")

  p4c <- ggplot(chipseq_gained,
                aes(x = distance_category, y = percentage, fill = mark)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8),
             width = 0.7, color = "black", linewidth = 0.2) +
    scale_fill_manual(values = c("H3K27ac" = "#e41a1c", "H3K27me3" = "#377eb8"), name = "") +
    scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
    labs(title = "C. Gained Loops", x = "", y = "% Loops") +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
    )

  # Panel D: Statistics summary
  h3k27ac_trend <- ifelse(coef_h3k27ac["Estimate"] < 0, "DECREASES", "INCREASES")
  h3k27me3_trend <- ifelse(coef_h3k27me3["Estimate"] > 0, "INCREASES", "DECREASES")
  hypothesis_supported <- coef_h3k27ac["Estimate"] < 0 && coef_h3k27me3["Estimate"] > 0

  stats_text <- sprintf(
    "Statistical Summary:\n\n
Trend Tests (logistic regression):\n
  H3K27ac:  beta = %.2e\n
            p = %.2e\n
  H3K27me3: beta = %.2e\n
            p = %.2e\n\n
Results:\n
  H3K27ac %s with distance\n
  H3K27me3 %s with distance\n\n
Conclusion:\n
  Hypothesis %s",
    coef_h3k27ac["Estimate"], coef_h3k27ac["Pr(>|z|)"],
    coef_h3k27me3["Estimate"], coef_h3k27me3["Pr(>|z|)"],
    h3k27ac_trend, h3k27me3_trend,
    ifelse(hypothesis_supported, "SUPPORTED", "NOT SUPPORTED")
  )

  p4d <- ggplot() +
    annotate("text", x = 0, y = 0.5, label = stats_text,
             hjust = 0, vjust = 0.5, size = 2.8, family = "mono") +
    xlim(-0.1, 1) + ylim(0, 1) +
    labs(title = "D. Statistics") +
    theme_void() +
    theme(plot.title = element_text(face = "bold", size = 11))

  # Combine
  p4_combined <- (p4a | p4b | p4c) / p4d +
    plot_layout(heights = c(2, 1)) +
    plot_annotation(
      title = "ChIP-seq Mark vs Loop Distance Analysis",
      subtitle = sprintf("%s timepoint: Testing H3K27ac/H3K27me3 distance relationship", toupper(timepoint)),
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey40")
      )
    )

  save_multiformat_ggplot(p4_combined, file.path(OUTPUT_DIR, "04_mark_comparison_summary"),
                          width = 12, height = 9)

  # ==========================================================================
  # FIGURE 5: Anchor Type Distribution by Distance (8-Category System)
  # ==========================================================================

  cat("Creating Figure 5: Anchor Type Distribution by Distance...\n")

  # Prepare anchor type data by combining anchor1 and anchor2
  anchor_type_data <- bind_rows(
    loops_directional %>%
      dplyr::select(direction_label, distance_category, anchor_type = anchor1_type) %>%
      dplyr::mutate(anchor = "Anchor 1"),
    loops_directional %>%
      dplyr::select(direction_label, distance_category, anchor_type = anchor2_type) %>%
      dplyr::mutate(anchor = "Anchor 2")
  )

  # Summarize anchor types by distance category and direction
  anchor_type_summary <- anchor_type_data %>%
    group_by(direction_label, distance_category, anchor_type) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(direction_label, distance_category) %>%
    mutate(percentage = 100 * count / sum(count)) %>%
    ungroup()

  # Set factor levels for proper ordering
  anchor_type_summary$anchor_type <- factor(
    anchor_type_summary$anchor_type, levels = ANCHOR_TYPE_ORDER
  )

  # Panel A: Stacked bar chart by distance and direction
  p5a <- ggplot(anchor_type_summary,
                aes(x = distance_category, y = percentage, fill = anchor_type)) +
    geom_bar(stat = "identity", position = "stack",
             color = "white", linewidth = 0.2) +
    facet_wrap(~direction_label) +
    scale_fill_manual(values = ANCHOR_COLORS, name = "Anchor Type") +
    scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
    labs(
      title = "Anchor Type Distribution by Loop Distance",
      subtitle = sprintf("%s timepoint: 8-category chromatin state classification",
                        toupper(timepoint)),
      x = "Loop Distance Category",
      y = "% of Anchors"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )

  # Panel B: Grouped bar for specific anchor types by distance
  # Focus on key regulatory categories
  key_types <- c("Active_Promoter", "Active_Enhancer", "Polycomb", "CTCF_Site")
  anchor_type_key <- anchor_type_summary %>%
    filter(anchor_type %in% key_types)

  p5b <- ggplot(anchor_type_key,
                aes(x = distance_category, y = percentage, fill = anchor_type)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8),
             width = 0.7, color = "black", linewidth = 0.2) +
    facet_wrap(~direction_label) +
    scale_fill_manual(values = ANCHOR_COLORS[key_types], name = "Anchor Type") +
    scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
    labs(
      title = "Key Regulatory Anchor Types by Distance",
      subtitle = "Focus on Active Promoter, Active Enhancer, Polycomb, CTCF",
      x = "Loop Distance Category",
      y = "% of Anchors"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey40"),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )

  # Combine panels
  p5_combined <- p5a / p5b +
    plot_annotation(
      title = "Anchor Chromatin State Classification by Loop Distance",
      subtitle = sprintf("%s timepoint: Using proper 8-category system from annotate_loops_extended.R",
                        toupper(timepoint)),
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey40")
      )
    )

  save_multiformat_ggplot(p5_combined, file.path(OUTPUT_DIR, "05_anchor_type_by_distance"),
                          width = 12, height = 12)

  # Also save anchor type summary table
  anchor_type_summary_wide <- anchor_type_summary %>%
    pivot_wider(
      names_from = anchor_type,
      values_from = c(count, percentage),
      values_fill = 0
    )
  write_tsv(anchor_type_summary_wide,
            file.path(OUTPUT_DIR, "anchor_type_distance_summary.tsv"))
  cat("Saved: anchor_type_distance_summary.tsv\n")

  # ==========================================================================
  # SECTION 6: EXPORT RESULTS
  # ==========================================================================

  cat("\n=== Step 6: Exporting Results ===\n")

  # Save distance summary
  write_tsv(distance_summary, file.path(OUTPUT_DIR, "chip_distance_summary.tsv"))
  cat("Saved: chip_distance_summary.tsv\n")

  # Save detailed chipseq summary by direction
  chipseq_wide <- chipseq_summary %>%
    pivot_wider(names_from = mark, values_from = percentage)
  write_tsv(chipseq_wide, file.path(OUTPUT_DIR, "chip_distance_direction_summary.tsv"))
  cat("Saved: chip_distance_direction_summary.tsv\n")

  # Save annotated loops with ChIP-seq overlaps
  write_tsv(loops_directional, file.path(OUTPUT_DIR, "loops_with_chip_overlaps.tsv"))
  cat("Saved: loops_with_chip_overlaps.tsv\n")

  # Save statistics report
  stats_file <- file.path(OUTPUT_DIR, "chip_distance_statistics.txt")
  sink(stats_file)
  cat("=== ChIP-seq Mark vs Loop Distance: Statistical Report ===\n")
  cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Timepoint:", timepoint, "\n\n")

  cat("--- DATA OVERVIEW ---\n")
  cat("Total loops analyzed:", nrow(loops_directional), "\n")
  cat("  Lost in BAP1-KO:", nrow(lost_loops), "\n")
  cat("  Gained in BAP1-KO:", nrow(gained_loops), "\n\n")

  cat("--- ChIP-seq PEAK FILES ---\n")
  cat("  H3K27ac:", chip_files$H3K27ac, sprintf("(%d peaks)\n", length(chip_peaks$H3K27ac)))
  cat("  H3K27me3:", chip_files$H3K27me3, sprintf("(%d peaks)\n", length(chip_peaks$H3K27me3)))
  cat("  H3K4me1:", chip_files$H3K4me1, sprintf("(%d peaks)\n", length(chip_peaks$H3K4me1)))
  cat("  H3K4me3:", chip_files$H3K4me3, sprintf("(%d peaks)\n\n", length(chip_peaks$H3K4me3)))

  cat("--- OVERALL MARK PRESENCE ---\n")
  cat(sprintf("  H3K27ac:  %.1f%% of loops\n", 100 * mean(loops_directional$has_H3K27ac)))
  cat(sprintf("  H3K27me3: %.1f%% of loops\n", 100 * mean(loops_directional$has_H3K27me3)))
  cat(sprintf("  H3K4me1:  %.1f%% of loops\n", 100 * mean(loops_directional$has_H3K4me1)))
  cat(sprintf("  H3K4me3:  %.1f%% of loops\n\n", 100 * mean(loops_directional$has_H3K4me3)))

  cat("--- DISTANCE CATEGORY SUMMARY ---\n")
  print(as.data.frame(distance_summary), row.names = FALSE)
  cat("\n")

  cat("--- STATISTICAL TESTS ---\n\n")
  cat("Chi-square tests (mark vs distance by direction):\n")
  cat(sprintf("  Lost loops - H3K27ac:  p = %.2e\n", chisq_lost_h3k27ac$p.value))
  cat(sprintf("  Lost loops - H3K27me3: p = %.2e\n", chisq_lost_h3k27me3$p.value))
  cat(sprintf("  Gained loops - H3K27ac:  p = %.2e\n", chisq_gained_h3k27ac$p.value))
  cat(sprintf("  Gained loops - H3K27me3: p = %.2e\n\n", chisq_gained_h3k27me3$p.value))

  cat("Logistic regression trend test (all loops):\n")
  cat("  H3K27ac:\n")
  print(summary(glm_h3k27ac)$coefficients)
  cat("\n  H3K27me3:\n")
  print(summary(glm_h3k27me3)$coefficients)

  cat("\n--- HYPOTHESIS EVALUATION ---\n")
  cat("Hypothesis: As loop distance increases, H3K27ac decreases and H3K27me3 increases\n\n")

  cat(sprintf("Results:\n"))
  cat(sprintf("  - H3K27ac %s with distance (beta = %.2e, p = %.2e)\n",
              h3k27ac_trend, coef_h3k27ac["Estimate"], coef_h3k27ac["Pr(>|z|)"]))
  cat(sprintf("  - H3K27me3 %s with distance (beta = %.2e, p = %.2e)\n",
              h3k27me3_trend, coef_h3k27me3["Estimate"], coef_h3k27me3["Pr(>|z|)"]))

  cat(sprintf("\nConclusion: Hypothesis %s\n",
              ifelse(hypothesis_supported, "SUPPORTED", "NOT SUPPORTED in expected direction")))

  cat("\n=================================================\n")
  sink()
  cat("Saved: chip_distance_statistics.txt\n")

  # ==========================================================================
  # COMPLETION
  # ==========================================================================

  cat("\n=== Analysis Complete for", toupper(timepoint), "===\n")
  cat("Output directory:", OUTPUT_DIR, "\n")
  cat("Files generated:\n")
  cat("  - 01_chipseq_distance_analysis.{pdf,svg,jpg}\n")
  cat("  - 02_mark_trend_lineplot.{pdf,svg,jpg}\n")
  cat("  - 03_mark_distance_heatmap.{pdf,svg,jpg}\n")
  cat("  - 04_mark_comparison_summary.{pdf,svg,jpg}\n")
  cat("  - 05_anchor_type_by_distance.{pdf,svg,jpg}\n")
  cat("  - chip_distance_summary.tsv\n")
  cat("  - chip_distance_direction_summary.tsv\n")
  cat("  - anchor_type_distance_summary.tsv\n")
  cat("  - loops_with_chip_overlaps.tsv\n")
  cat("  - chip_distance_statistics.txt\n\n")

  # Return key results
  return(list(
    timepoint = timepoint,
    n_loops = nrow(loops_directional),
    distance_summary = distance_summary,
    h3k27ac_trend = coef_h3k27ac,
    h3k27me3_trend = coef_h3k27me3,
    hypothesis_supported = hypothesis_supported
  ))
}

# ==============================================================================
# MAIN EXECUTION
# ==============================================================================

results <- list()

for (tp in TIMEPOINTS_TO_RUN) {
  tryCatch({
    results[[tp]] <- run_chip_distance_analysis(tp)
  }, error = function(e) {
    cat(sprintf("\nERROR processing %s timepoint: %s\n", tp, e$message))
  })
}

# Print summary across timepoints
cat("\n=== All Timepoints Complete ===\n")
cat("Output directories:\n")
for (tp in TIMEPOINTS_TO_RUN) {
  cat(sprintf("  - %s\n", OUTPUT_DIRS[[tp]]))
}

if (length(results) > 0) {
  cat("\nHypothesis Summary:\n")
  for (tp in names(results)) {
    r <- results[[tp]]
    cat(sprintf("  %s: %s (H3K27ac beta=%.2e, H3K27me3 beta=%.2e)\n",
                toupper(tp),
                ifelse(r$hypothesis_supported, "SUPPORTED", "NOT SUPPORTED"),
                r$h3k27ac_trend["Estimate"],
                r$h3k27me3_trend["Estimate"]))
  }
}

cat("\nEnd time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
