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
    H3K4me3  = "peaks/beds/H3K4me3CerebellumEarly2.bed"
  ),
  late = list(
    H3K27ac  = "peaks/beds/H3K27acCerebellumLate2.bed",
    H3K27me3 = "peaks/beds/H3K27me3CerebellumLate1.bed",
    H3K4me1  = "peaks/beds/H3K4me1CerebellumLate1.bed",
    H3K4me3  = "peaks/beds/H3K4me3CerebellumLate2.bed"
  )
)

# Output directories by timepoint
OUTPUT_DIRS <- list(
  late = "output/chip_distance/late",
  early = "output/chip_distance/early"
)

# Distance category order and thresholds
DISTANCE_ORDER <- c("<100kb", "100-500kb", "500kb-1Mb", ">1Mb")
DISTANCE_BREAKS <- c(0, 100000, 500000, 1000000, Inf)

# Color palettes
COLORS <- list(
  down = "#d73027",      # Red for down/lost in mutant
  up = "#4575b4",        # Blue for up/gained in mutant
  neutral = "#999999"
)

# ChIP-seq mark colors (distinct, colorblind-friendly)
MARK_COLORS <- c(
  "H3K27ac" = "#e41a1c",   # Red - active mark

  "H3K27me3" = "#377eb8",  # Blue - repressive mark
  "H3K4me1" = "#4daf4a",   # Green - enhancer mark
  "H3K4me3" = "#984ea3"    # Purple - promoter mark
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

#' Compute ChIP-seq overlaps for anchor GRanges
compute_chip_overlaps <- function(anchor_gr, chip_peaks_list) {
  overlaps <- data.frame(
    H3K27ac_overlap = countOverlaps(anchor_gr, chip_peaks_list$H3K27ac) > 0,
    H3K27me3_overlap = countOverlaps(anchor_gr, chip_peaks_list$H3K27me3) > 0,
    H3K4me1_overlap = countOverlaps(anchor_gr, chip_peaks_list$H3K4me1) > 0,
    H3K4me3_overlap = countOverlaps(anchor_gr, chip_peaks_list$H3K4me3) > 0
  )
  return(overlaps)
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

  # Add distance category
  loops <- loops %>%
    mutate(
      distance_category = assign_distance_category(loop_distance),
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

  cat("Directional loops:", nrow(loops_directional), "\n")
  cat("  - Lost (down_in_mutant):", sum(loops_directional$direction == "down_in_mutant"), "\n")
  cat("  - Gained (up_in_mutant):", sum(loops_directional$direction == "up_in_mutant"), "\n\n")

  # ==========================================================================
  # SECTION 3: LOAD ChIP-seq PEAKS
  # ==========================================================================

  cat("=== Step 2: Loading ChIP-seq Peak Files ===\n")

  chip_peaks <- list(
    H3K27ac = load_chip_peaks(chip_files$H3K27ac, "H3K27ac"),
    H3K27me3 = load_chip_peaks(chip_files$H3K27me3, "H3K27me3"),
    H3K4me1 = load_chip_peaks(chip_files$H3K4me1, "H3K4me1"),
    H3K4me3 = load_chip_peaks(chip_files$H3K4me3, "H3K4me3")
  )
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

  # Compute overlaps for both anchors
  anchor1_overlaps <- compute_chip_overlaps(anchor1_gr, chip_peaks)
  anchor2_overlaps <- compute_chip_overlaps(anchor2_gr, chip_peaks)

  # Add overlap columns to dataframe
  loops_directional <- loops_directional %>%
    mutate(
      # Anchor 1 overlaps
      anchor1_H3K27ac = anchor1_overlaps$H3K27ac_overlap,
      anchor1_H3K27me3 = anchor1_overlaps$H3K27me3_overlap,
      anchor1_H3K4me1 = anchor1_overlaps$H3K4me1_overlap,
      anchor1_H3K4me3 = anchor1_overlaps$H3K4me3_overlap,
      # Anchor 2 overlaps
      anchor2_H3K27ac = anchor2_overlaps$H3K27ac_overlap,
      anchor2_H3K27me3 = anchor2_overlaps$H3K27me3_overlap,
      anchor2_H3K4me1 = anchor2_overlaps$H3K4me1_overlap,
      anchor2_H3K4me3 = anchor2_overlaps$H3K4me3_overlap,
      # Loop-level: at least one anchor has the mark
      has_H3K27ac = anchor1_H3K27ac | anchor2_H3K27ac,
      has_H3K27me3 = anchor1_H3K27me3 | anchor2_H3K27me3,
      has_H3K4me1 = anchor1_H3K4me1 | anchor2_H3K4me1,
      has_H3K4me3 = anchor1_H3K4me3 | anchor2_H3K4me3
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
  # SECTION 5: DISTANCE-STRATIFIED ANALYSIS
  # ==========================================================================

  cat("=== Step 4: Distance-Stratified Analysis ===\n")

  # Summarize by distance category
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

  # Summarize by distance and direction
  distance_direction_summary <- loops_directional %>%
    group_by(distance_category, direction_label) %>%
    summarise(
      n_loops = n(),
      pct_H3K27ac = mean(has_H3K27ac) * 100,
      pct_H3K27me3 = mean(has_H3K27me3) * 100,
      pct_H3K4me1 = mean(has_H3K4me1) * 100,
      pct_H3K4me3 = mean(has_H3K4me3) * 100,
      .groups = "drop"
    )

  cat("Distance category summary (all loops):\n")
  print(as.data.frame(distance_summary), row.names = FALSE)
  cat("\n")

  # ==========================================================================
  # SECTION 6: STATISTICAL TESTS
  # ==========================================================================

  cat("=== Step 5: Statistical Tests ===\n")

  # Chi-square test: H3K27ac presence independent of distance?
  h3k27ac_table <- table(loops_directional$distance_category, loops_directional$has_H3K27ac)
  chisq_h3k27ac <- chisq.test(h3k27ac_table)
  cat("Chi-square test (H3K27ac vs distance):\n")
  cat(sprintf("  X-squared = %.2f, df = %d, p = %.2e\n",
              chisq_h3k27ac$statistic, chisq_h3k27ac$parameter, chisq_h3k27ac$p.value))

  # Chi-square test: H3K27me3 presence independent of distance?
  h3k27me3_table <- table(loops_directional$distance_category, loops_directional$has_H3K27me3)
  chisq_h3k27me3 <- chisq.test(h3k27me3_table)
  cat("Chi-square test (H3K27me3 vs distance):\n")
  cat(sprintf("  X-squared = %.2f, df = %d, p = %.2e\n",
              chisq_h3k27me3$statistic, chisq_h3k27me3$parameter, chisq_h3k27me3$p.value))

  # Spearman correlation: distance vs mark presence
  # Assign numeric distance (use midpoint of category)
  loops_directional <- loops_directional %>%
    mutate(
      distance_numeric = case_when(
        distance_category == "<100kb" ~ 50000,
        distance_category == "100-500kb" ~ 300000,
        distance_category == "500kb-1Mb" ~ 750000,
        distance_category == ">1Mb" ~ 1500000,
        TRUE ~ NA_real_
      )
    )

  cor_h3k27ac <- cor.test(loops_directional$loop_distance, as.numeric(loops_directional$has_H3K27ac),
                          method = "spearman")
  cor_h3k27me3 <- cor.test(loops_directional$loop_distance, as.numeric(loops_directional$has_H3K27me3),
                           method = "spearman")

  cat("\nSpearman correlation (distance vs mark presence):\n")
  cat(sprintf("  H3K27ac:  rho = %.3f, p = %.2e\n", cor_h3k27ac$estimate, cor_h3k27ac$p.value))
  cat(sprintf("  H3K27me3: rho = %.3f, p = %.2e\n\n", cor_h3k27me3$estimate, cor_h3k27me3$p.value))

  # Cochran-Armitage trend test approximation using logistic regression
  cat("Trend test (logistic regression):\n")

  # H3K27ac trend
  glm_h3k27ac <- glm(has_H3K27ac ~ loop_distance, data = loops_directional, family = binomial)
  coef_h3k27ac <- summary(glm_h3k27ac)$coefficients["loop_distance", ]
  cat(sprintf("  H3K27ac:  beta = %.2e, z = %.2f, p = %.2e\n",
              coef_h3k27ac["Estimate"], coef_h3k27ac["z value"], coef_h3k27ac["Pr(>|z|)"]))

  # H3K27me3 trend
  glm_h3k27me3 <- glm(has_H3K27me3 ~ loop_distance, data = loops_directional, family = binomial)
  coef_h3k27me3 <- summary(glm_h3k27me3)$coefficients["loop_distance", ]
  cat(sprintf("  H3K27me3: beta = %.2e, z = %.2f, p = %.2e\n\n",
              coef_h3k27me3["Estimate"], coef_h3k27me3["z value"], coef_h3k27me3["Pr(>|z|)"]))

  # ==========================================================================
  # FIGURE 1: Line Plot - Mark Presence vs Distance Category (KEY FIGURE)
  # ==========================================================================

  cat("=== Step 6: Generating Figures ===\n")
  cat("Creating Figure 1: Mark vs Distance Line Plot (key hypothesis test)...\n")

  # Reshape for plotting
  distance_long <- distance_summary %>%
    pivot_longer(
      cols = starts_with("pct_"),
      names_to = "mark",
      values_to = "percentage"
    ) %>%
    mutate(
      mark = gsub("pct_", "", mark),
      mark = factor(mark, levels = c("H3K27ac", "H3K27me3", "H3K4me1", "H3K4me3"))
    )

  # Focus on H3K27ac vs H3K27me3 for primary hypothesis
  distance_focus <- distance_long %>%
    filter(mark %in% c("H3K27ac", "H3K27me3"))

  p1_line <- ggplot(distance_focus,
                    aes(x = distance_category, y = percentage, color = mark, group = mark)) +
    geom_line(linewidth = 1.5) +
    geom_point(size = 4) +
    scale_color_manual(
      values = c("H3K27ac" = MARK_COLORS["H3K27ac"], "H3K27me3" = MARK_COLORS["H3K27me3"]),
      name = "Histone Mark",
      labels = c("H3K27ac (Active)", "H3K27me3 (Repressive)")
    ) +
    labs(
      title = "ChIP-seq Mark Enrichment by Loop Distance",
      subtitle = sprintf("Hypothesis: H3K27ac decreases, H3K27me3 increases with distance\n(n = %d loops, %s timepoint)",
                        nrow(loops_directional), timepoint),
      x = "Loop Distance Category",
      y = "% Loops with Mark at Anchor"
    ) +
    annotate("text", x = 3.5, y = max(distance_focus$percentage) * 0.95,
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

  save_multiformat_ggplot(p1_line, file.path(OUTPUT_DIR, "01_mark_vs_distance_lineplot"),
                          width = 9, height = 7)

  # ==========================================================================
  # FIGURE 2: All Four Marks Line Plot
  # ==========================================================================

  cat("Creating Figure 2: All Four Marks Line Plot...\n")

  p2_all_marks <- ggplot(distance_long,
                         aes(x = distance_category, y = percentage, color = mark, group = mark)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    scale_color_manual(values = MARK_COLORS, name = "Histone Mark") +
    labs(
      title = "All ChIP-seq Marks by Loop Distance",
      subtitle = sprintf("%s timepoint (n = %d loops)", timepoint, nrow(loops_directional)),
      x = "Loop Distance Category",
      y = "% Loops with Mark at Anchor"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"),
      legend.position = "right",
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
      panel.grid.minor = element_blank()
    )

  save_multiformat_ggplot(p2_all_marks, file.path(OUTPUT_DIR, "02_all_marks_lineplot"),
                          width = 10, height = 7)

  # ==========================================================================
  # FIGURE 3: Grouped Bar Chart by Distance and Direction
  # ==========================================================================

  cat("Creating Figure 3: Mark Distance Barplot by Direction...\n")

  # Reshape for grouped bar
  distance_direction_long <- distance_direction_summary %>%
    pivot_longer(
      cols = starts_with("pct_"),
      names_to = "mark",
      values_to = "percentage"
    ) %>%
    mutate(
      mark = gsub("pct_", "", mark),
      mark = factor(mark, levels = c("H3K27ac", "H3K27me3", "H3K4me1", "H3K4me3"))
    ) %>%
    filter(mark %in% c("H3K27ac", "H3K27me3"))

  p3_bar <- ggplot(distance_direction_long,
                   aes(x = distance_category, y = percentage, fill = mark)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8),
             width = 0.7, color = "black", linewidth = 0.3) +
    facet_wrap(~direction_label) +
    scale_fill_manual(
      values = c("H3K27ac" = MARK_COLORS["H3K27ac"], "H3K27me3" = MARK_COLORS["H3K27me3"]),
      name = "Histone Mark",
      labels = c("H3K27ac (Active)", "H3K27me3 (Repressive)")
    ) +
    labs(
      title = "ChIP-seq Mark Enrichment by Distance and Direction",
      subtitle = sprintf("%s timepoint", timepoint),
      x = "Loop Distance Category",
      y = "% Loops with Mark at Anchor"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"),
      legend.position = "top",
      strip.text = element_text(face = "bold", size = 11),
      axis.text.x = element_text(size = 9, angle = 45, hjust = 1),
      panel.grid.minor = element_blank()
    )

  save_multiformat_ggplot(p3_bar, file.path(OUTPUT_DIR, "02_mark_distance_barplot"),
                          width = 11, height = 7)

  # ==========================================================================
  # FIGURE 4: Heatmap of Mark Enrichment by Distance
  # ==========================================================================

  cat("Creating Figure 4: Mark Distance Heatmap...\n")

  # Create matrix for heatmap
  heatmap_data <- distance_summary %>%
    dplyr::select(distance_category, pct_H3K27ac, pct_H3K27me3, pct_H3K4me1, pct_H3K4me3) %>%
    pivot_longer(cols = -distance_category, names_to = "mark", values_to = "percentage") %>%
    mutate(
      mark = gsub("pct_", "", mark),
      mark = factor(mark, levels = c("H3K27ac", "H3K27me3", "H3K4me1", "H3K4me3")),
      distance_category = factor(distance_category, levels = DISTANCE_ORDER)
    )

  p4_heatmap <- ggplot(heatmap_data,
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
      subtitle = sprintf("%s timepoint: % of loops with mark at one or more anchors", timepoint),
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

  save_multiformat_ggplot(p4_heatmap, file.path(OUTPUT_DIR, "03_mark_distance_heatmap"),
                          width = 8, height = 6)

  # ==========================================================================
  # FIGURE 5: Combined Summary Figure
  # ==========================================================================

  cat("Creating Figure 5: Combined Summary Figure...\n")

  # Panel A: Line plot (simplified)
  p5a <- ggplot(distance_focus,
                aes(x = distance_category, y = percentage, color = mark, group = mark)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    scale_color_manual(
      values = c("H3K27ac" = MARK_COLORS["H3K27ac"], "H3K27me3" = MARK_COLORS["H3K27me3"]),
      name = ""
    ) +
    labs(title = "A. Mark Trend by Distance", x = "", y = "% Loops") +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
    )

  # Panel B: Bar by direction (simplified)
  p5b <- ggplot(distance_direction_long %>% filter(direction_label == "Lost in BAP1-KO"),
                aes(x = distance_category, y = percentage, fill = mark)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    scale_fill_manual(
      values = c("H3K27ac" = MARK_COLORS["H3K27ac"], "H3K27me3" = MARK_COLORS["H3K27me3"]),
      name = ""
    ) +
    labs(title = "B. Lost Loops", x = "", y = "% Loops") +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
    )

  p5c <- ggplot(distance_direction_long %>% filter(direction_label == "Gained in BAP1-KO"),
                aes(x = distance_category, y = percentage, fill = mark)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    scale_fill_manual(
      values = c("H3K27ac" = MARK_COLORS["H3K27ac"], "H3K27me3" = MARK_COLORS["H3K27me3"]),
      name = ""
    ) +
    labs(title = "C. Gained Loops", x = "", y = "% Loops") +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
    )

  # Panel D: Statistics summary
  stats_text <- sprintf(
    "Statistical Summary:\n\n
Trend Tests (logistic regression):\n
  H3K27ac:  beta = %.2e (p = %.2e)\n
  H3K27me3: beta = %.2e (p = %.2e)\n\n
Spearman Correlation:\n
  H3K27ac vs distance:  rho = %.3f\n
  H3K27me3 vs distance: rho = %.3f\n\n
Interpretation:\n
  %s",
    coef_h3k27ac["Estimate"], coef_h3k27ac["Pr(>|z|)"],
    coef_h3k27me3["Estimate"], coef_h3k27me3["Pr(>|z|)"],
    cor_h3k27ac$estimate, cor_h3k27me3$estimate,
    ifelse(coef_h3k27ac["Estimate"] < 0 && coef_h3k27me3["Estimate"] > 0,
           "H3K27ac decreases and\nH3K27me3 increases with\ndistance (supports hypothesis)",
           "Trend pattern differs\nfrom hypothesis")
  )

  p5d <- ggplot() +
    annotate("text", x = 0, y = 0.5, label = stats_text,
             hjust = 0, vjust = 0.5, size = 2.8, family = "mono") +
    xlim(-0.1, 1) + ylim(0, 1) +
    labs(title = "D. Statistics") +
    theme_void() +
    theme(plot.title = element_text(face = "bold", size = 11))

  # Legend plot
  legend_data <- data.frame(
    x = c(1, 2),
    y = c(1, 1),
    label = c("H3K27ac (Active)", "H3K27me3 (Repressive)"),
    color = c(MARK_COLORS["H3K27ac"], MARK_COLORS["H3K27me3"])
  )

  p_legend <- ggplot(legend_data, aes(x = x, y = y, color = label)) +
    geom_point(size = 4) +
    geom_text(aes(label = label), hjust = -0.15, size = 3.5) +
    scale_color_manual(values = c("H3K27ac (Active)" = MARK_COLORS["H3K27ac"],
                                  "H3K27me3 (Repressive)" = MARK_COLORS["H3K27me3"])) +
    xlim(0.5, 5) +
    theme_void() +
    theme(legend.position = "none")

  # Combine
  p5_combined <- (p5a | p5b | p5c) / (p5d | p_legend) +
    plot_layout(heights = c(2, 1)) +
    plot_annotation(
      title = "ChIP-seq Mark vs Loop Distance Analysis",
      subtitle = sprintf("%s timepoint: Testing H3K27ac/H3K27me3 distance relationship", toupper(timepoint)),
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey40")
      )
    )

  save_multiformat_ggplot(p5_combined, file.path(OUTPUT_DIR, "04_mark_comparison_summary"),
                          width = 12, height = 9)

  # ==========================================================================
  # SECTION 7: EXPORT RESULTS
  # ==========================================================================

  cat("\n=== Step 7: Exporting Results ===\n")

  # Save distance summary
  write_tsv(distance_summary, file.path(OUTPUT_DIR, "chip_distance_summary.tsv"))
  cat("Saved: chip_distance_summary.tsv\n")

  # Save detailed distance x direction summary
  write_tsv(distance_direction_summary, file.path(OUTPUT_DIR, "chip_distance_direction_summary.tsv"))
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
  cat("  Lost in BAP1-KO:", sum(loops_directional$direction == "down_in_mutant"), "\n")
  cat("  Gained in BAP1-KO:", sum(loops_directional$direction == "up_in_mutant"), "\n\n")

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
  cat("Chi-square test (H3K27ac vs distance):\n")
  print(chisq_h3k27ac)
  cat("\nChi-square test (H3K27me3 vs distance):\n")
  print(chisq_h3k27me3)
  cat("\nSpearman correlation (H3K27ac vs distance):\n")
  print(cor_h3k27ac)
  cat("\nSpearman correlation (H3K27me3 vs distance):\n")
  print(cor_h3k27me3)
  cat("\nLogistic regression trend test:\n")
  cat("  H3K27ac:\n")
  print(summary(glm_h3k27ac)$coefficients)
  cat("\n  H3K27me3:\n")
  print(summary(glm_h3k27me3)$coefficients)

  cat("\n--- HYPOTHESIS EVALUATION ---\n")
  cat("Hypothesis: As loop distance increases, H3K27ac decreases and H3K27me3 increases\n\n")

  h3k27ac_trend <- ifelse(coef_h3k27ac["Estimate"] < 0, "DECREASES", "INCREASES")
  h3k27me3_trend <- ifelse(coef_h3k27me3["Estimate"] > 0, "INCREASES", "DECREASES")

  cat(sprintf("Results:\n"))
  cat(sprintf("  - H3K27ac %s with distance (beta = %.2e, p = %.2e)\n",
              h3k27ac_trend, coef_h3k27ac["Estimate"], coef_h3k27ac["Pr(>|z|)"]))
  cat(sprintf("  - H3K27me3 %s with distance (beta = %.2e, p = %.2e)\n",
              h3k27me3_trend, coef_h3k27me3["Estimate"], coef_h3k27me3["Pr(>|z|)"]))

  hypothesis_supported <- coef_h3k27ac["Estimate"] < 0 && coef_h3k27me3["Estimate"] > 0
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
  cat("  - 01_mark_vs_distance_lineplot.{pdf,svg,jpg}\n")
  cat("  - 02_all_marks_lineplot.{pdf,svg,jpg}\n")
  cat("  - 02_mark_distance_barplot.{pdf,svg,jpg}\n")
  cat("  - 03_mark_distance_heatmap.{pdf,svg,jpg}\n")
  cat("  - 04_mark_comparison_summary.{pdf,svg,jpg}\n")
  cat("  - chip_distance_summary.tsv\n")
  cat("  - chip_distance_direction_summary.tsv\n")
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
