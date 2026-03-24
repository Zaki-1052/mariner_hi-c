# scripts/loop_distance_k27me3_filtered.R
# H3K27me3-Filtered Loop Distance Analysis
# Generates CDF and density plots for loops filtered by H3K27me3 and bivalent overlaps
#
# Purpose:
#   Create distance distribution plots (CDF + density) for three H3K27me3-related subsets:
#   1. K27me3-anchored: at least one anchor overlaps H3K27me3
#   2. K27me3-K27me3: BOTH anchors overlap H3K27me3
#   3. Bivalent: at least one anchor overlaps bivalent (K4me3+K27me3) peaks
#
# Usage:
#   Rscript scripts/loop_distance_k27me3_filtered.R                    # Default: late
#   Rscript scripts/loop_distance_k27me3_filtered.R --timepoint late   # Late timepoint
#   Rscript scripts/loop_distance_k27me3_filtered.R --timepoint early  # Early timepoint
#   Rscript scripts/loop_distance_k27me3_filtered.R --timepoint both   # Run both
#
# Output:
#   output/loops_k27me3_filtered/{early,late}/
#   - 01_cdf_k27me3_anchored.{pdf,svg,jpg}
#   - 01_cdf_k27me3_both.{pdf,svg,jpg}
#   - 01_cdf_bivalent.{pdf,svg,jpg}
#   - 03_density_k27me3_anchored.{pdf,svg,jpg}
#   - 03_density_k27me3_both.{pdf,svg,jpg}
#   - 03_density_bivalent.{pdf,svg,jpg}
#   - filter_summary.txt

# ==============================================================================
# SECTION 1: SETUP
# ==============================================================================

cat("=== H3K27me3-Filtered Loop Distance Analysis ===\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(scales)
  library(GenomicRanges)
})

source("data/scripts/_shared/multi_format_output.R")  # Original: source("scripts/utils/multi_format_output.R")

# ==============================================================================
# SECTION 2: CONFIGURATION
# ==============================================================================

# Input files by timepoint (characterized_loops.tsv from downstream_analysis.R)
INPUT_FILES <- list(
  late = "data/upstream/loop_calls/late_characterized_loops.tsv",  # Original: "25042-late_outputs/merged_loops/characterized_loops.tsv"
  early = "250831-early_outputs/merged_loops/characterized_loops.tsv"  # TODO: not in data/
)

# ChIP-seq peak files by timepoint
PEAK_FILES <- list(
  late = list(
    h3k27me3 = "peaks/beds/H3K27me3CerebellumLate1.bed",  # TODO: not in data/
    bivalent = "peaks/beds/Bivalent_Cerebellum_Late.bed"  # TODO: not in data/
  ),
  early = list(
    h3k27me3 = "peaks/beds/H3K27me3CerebellumEarly1.bed",  # TODO: not in data/
    bivalent = "peaks/beds/Bivalent_Cerebellum_Early.bed"  # TODO: not in data/
  )
)

# Output directories by timepoint
# Original: OUTPUT_DIRS <- list(late = "output/loops_k27me3_filtered/late", early = "output/loops_k27me3_filtered/early")
PLOT_DIR <- "data/plots/figure_2_loop_rewiring"
TSV_DIR  <- "data/tsvs/figure_2_loop_rewiring"

# Color scheme (consistent with loop_distance_analysis.R)
COLORS <- list(
  down = "#d73027",      # Red for down/lost in mutant
  up = "#4575b4",        # Blue for up/gained in mutant
  down_light = "#f4a582",
  up_light = "#92c5de"
)

# Direction labels
DIRECTION_LABELS <- c(
  "down_in_mutant" = "Lost in BAP1-KO",
  "up_in_mutant" = "Gained in BAP1-KO"
)

# ==============================================================================
# SECTION 3: ARGUMENT PARSING
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
      cat("Usage: Rscript scripts/loop_distance_k27me3_filtered.R [OPTIONS]\n\n")
      cat("Options:\n")
      cat("  --timepoint TP  Timepoint: 'early', 'late', or 'both' (default: late)\n")
      cat("  --help, -h      Show this help message\n\n")
      cat("Output:\n")
      cat("  output/loops_k27me3_filtered/{early,late}/\n")
      cat("    - 01_cdf_k27me3_anchored.{pdf,svg,jpg}\n")
      cat("    - 01_cdf_k27me3_both.{pdf,svg,jpg}\n")
      cat("    - 01_cdf_bivalent.{pdf,svg,jpg}\n")
      cat("    - 03_density_k27me3_anchored.{pdf,svg,jpg}\n")
      cat("    - 03_density_k27me3_both.{pdf,svg,jpg}\n")
      cat("    - 03_density_bivalent.{pdf,svg,jpg}\n")
      cat("    - filter_summary.txt\n\n")
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
# SECTION 4: HELPER FUNCTIONS
# ==============================================================================

#' Load ChIP-seq peaks from BED file
#'
#' @param bed_path Path to BED file
#' @param mark_name Name for logging
#' @return GRanges object
load_chip_peaks <- function(bed_path, mark_name = "ChIP") {
  if (!file.exists(bed_path)) {
    stop(sprintf("%s bed file not found: %s", mark_name, bed_path))
  }
  cat(sprintf("  Loading %s peaks from: %s\n", mark_name, bed_path))

  df <- read.table(bed_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  gr <- GRanges(
    seqnames = df$V1,
    ranges = IRanges(start = df$V2, end = df$V3)
  )

  cat(sprintf("    Loaded %d peaks\n", length(gr)))
  gr
}

#' Generate CDF plot for a filtered subset
#'
#' @param loops_df Data frame with filtered loops
#' @param subset_name Name of the subset for title
#' @param output_path Output path (without extension)
#' @param colors Color list
#' @return ggplot object (invisibly)
generate_cdf_plot <- function(loops_df, subset_name, output_path, colors) {
  # Separate by direction
  lost_loops <- loops_df %>% filter(direction == "down_in_mutant")
  gained_loops <- loops_df %>% filter(direction == "up_in_mutant")

  # Skip if no data
  if (nrow(lost_loops) < 5 || nrow(gained_loops) < 5) {
    cat(sprintf("    Skipping CDF for %s: insufficient data (Lost: %d, Gained: %d)\n",
                subset_name, nrow(lost_loops), nrow(gained_loops)))
    return(invisible(NULL))
  }

  # Calculate medians
  median_lost <- median(lost_loops$loop_distance)
  median_gained <- median(gained_loops$loop_distance)

  # KS test
  ks_test <- ks.test(lost_loops$loop_distance, gained_loops$loop_distance)

  # Create plot
  p <- ggplot(loops_df, aes(x = loop_distance_kb, color = direction_label)) +
    stat_ecdf(geom = "step", linewidth = 1.2) +
    scale_color_manual(
      values = c("Lost in BAP1-KO" = colors$down, "Gained in BAP1-KO" = colors$up),
      name = "Direction"
    ) +
    scale_x_log10(
      labels = comma,
      breaks = c(10, 100, 1000, 10000),
      limits = c(10, 100000)
    ) +
    geom_vline(xintercept = median_lost / 1000, color = colors$down,
               linetype = "dashed", alpha = 0.7) +
    geom_vline(xintercept = median_gained / 1000, color = colors$up,
               linetype = "dashed", alpha = 0.7) +
    annotate("text", x = 50000, y = 0.15,
             label = sprintf("KS test p = %.2e", ks_test$p.value),
             hjust = 1, size = 4) +
    annotate("text", x = 50000, y = 0.08,
             label = sprintf("Median: Lost = %d kb, Gained = %d kb",
                            round(median_lost / 1000), round(median_gained / 1000)),
             hjust = 1, size = 3.5) +
    annotate("text", x = 50000, y = 0.22,
             label = sprintf("n = %d loops (%d lost, %d gained)",
                            nrow(loops_df), nrow(lost_loops), nrow(gained_loops)),
             hjust = 1, size = 3.5, fontface = "italic") +
    labs(
      title = sprintf("Loop Distance CDF: %s", subset_name),
      subtitle = "BAP1-KO preferentially loses long-range loops",
      x = "Loop Distance (kb, log scale)",
      y = "Cumulative Proportion"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey40"),
      legend.position = c(0.15, 0.85),
      legend.background = element_rect(fill = "white", color = NA),
      panel.grid.minor = element_blank()
    )

  save_multiformat_ggplot(p, output_path, width = 8, height = 6)
  invisible(p)
}

#' Generate density plot for a filtered subset
#'
#' @param loops_df Data frame with filtered loops
#' @param subset_name Name of the subset for title
#' @param output_path Output path (without extension)
#' @param colors Color list
#' @return ggplot object (invisibly)
generate_density_plot <- function(loops_df, subset_name, output_path, colors) {
  # Separate by direction
  lost_loops <- loops_df %>% filter(direction == "down_in_mutant")
  gained_loops <- loops_df %>% filter(direction == "up_in_mutant")

  # Skip if no data
  if (nrow(lost_loops) < 5 || nrow(gained_loops) < 5) {
    cat(sprintf("    Skipping density for %s: insufficient data (Lost: %d, Gained: %d)\n",
                subset_name, nrow(lost_loops), nrow(gained_loops)))
    return(invisible(NULL))
  }

  # Calculate medians
  median_lost <- median(lost_loops$loop_distance)
  median_gained <- median(gained_loops$loop_distance)

  # Wilcoxon test
  wilcox_test <- wilcox.test(lost_loops$loop_distance, gained_loops$loop_distance)

  # Create plot
  p <- ggplot(loops_df, aes(x = loop_distance_kb, fill = direction_label)) +
    geom_density(alpha = 0.5, color = "black", linewidth = 0.5) +
    geom_vline(xintercept = median_lost / 1000, color = colors$down,
               linetype = "dashed", linewidth = 1) +
    geom_vline(xintercept = median_gained / 1000, color = colors$up,
               linetype = "dashed", linewidth = 1) +
    scale_fill_manual(
      values = c("Lost in BAP1-KO" = colors$down, "Gained in BAP1-KO" = colors$up),
      name = ""
    ) +
    scale_x_log10(
      labels = comma,
      breaks = c(10, 100, 1000, 10000),
      limits = c(10, 100000)
    ) +
    annotate("text", x = median_lost / 1000 * 1.3, y = Inf,
             label = sprintf("Median\n%d kb", round(median_lost / 1000)),
             color = colors$down, vjust = 1.5, hjust = 0, size = 3.5, fontface = "bold") +
    annotate("text", x = median_gained / 1000 * 0.7, y = Inf,
             label = sprintf("Median\n%d kb", round(median_gained / 1000)),
             color = colors$up, vjust = 1.5, hjust = 1, size = 3.5, fontface = "bold") +
    annotate("text", x = 20, y = 0,
             label = sprintf("Wilcoxon p = %.2e | n = %d", wilcox_test$p.value, nrow(loops_df)),
             hjust = 0, vjust = -0.5, size = 4) +
    labs(
      title = sprintf("Loop Distance Density: %s", subset_name),
      subtitle = "Lost loops are systematically longer than gained loops",
      x = "Loop Distance (kb, log scale)",
      y = "Density"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey40"),
      legend.position = "top",
      panel.grid.minor = element_blank()
    )

  save_multiformat_ggplot(p, output_path, width = 8, height = 6)
  invisible(p)
}

# ==============================================================================
# SECTION 5: MAIN ANALYSIS FUNCTION
# ==============================================================================

run_k27me3_filtered_analysis <- function(timepoint) {
  cat("\n")
  cat("============================================================\n")
  cat(sprintf("Processing %s timepoint\n", toupper(timepoint)))
  cat("============================================================\n\n")

  # Set paths
  input_file <- INPUT_FILES[[timepoint]]
  peak_files <- PEAK_FILES[[timepoint]]
  dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)  # Original: dir.create(OUTPUT_DIRS[[timepoint]], ...)
  dir.create(TSV_DIR, showWarnings = FALSE, recursive = TRUE)  # Original: dir.create(OUTPUT_DIRS[[timepoint]], ...)

  cat("Input file:", input_file, "\n")
  cat("Plot directory:", PLOT_DIR, "\n")
  cat("TSV directory:", TSV_DIR, "\n\n")

  # Validate input
  if (!file.exists(input_file)) {
    stop(sprintf("ERROR: Input file not found: %s\nRun downstream_analysis.R first.", input_file))
  }

  # ============================================================================
  # Step 1: Load characterized loops
  # ============================================================================

  cat("=== Step 1: Loading Data ===\n")
  loops <- read_tsv(input_file, show_col_types = FALSE)
  cat("Loaded", nrow(loops), "differential loops\n")

  # Prepare data
  loops <- loops %>%
    mutate(
      direction_label = factor(
        case_when(
          direction == "down_in_mutant" ~ "Lost in BAP1-KO",
          direction == "up_in_mutant" ~ "Gained in BAP1-KO",
          TRUE ~ "Other"
        ),
        levels = c("Lost in BAP1-KO", "Gained in BAP1-KO")
      ),
      loop_distance_kb = loop_distance / 1000
    )

  # Filter to directional loops only
  loops_directional <- loops %>%
    filter(direction %in% c("up_in_mutant", "down_in_mutant"))

  cat("Directional loops:", nrow(loops_directional), "\n")
  cat("  - Lost (down_in_mutant):", sum(loops_directional$direction == "down_in_mutant"), "\n")
  cat("  - Gained (up_in_mutant):", sum(loops_directional$direction == "up_in_mutant"), "\n\n")

  # ============================================================================
  # Step 2: Load ChIP-seq peaks
  # ============================================================================

  cat("=== Step 2: Loading ChIP-seq Peaks ===\n")
  k27me3_gr <- load_chip_peaks(peak_files$h3k27me3, "H3K27me3")
  bivalent_gr <- load_chip_peaks(peak_files$bivalent, "Bivalent")
  cat("\n")

  # ============================================================================
  # Step 3: Create anchor GRanges and compute overlaps
  # ============================================================================

  cat("=== Step 3: Computing Overlaps ===\n")

  # Create GRanges for anchors
  anchor1_gr <- GRanges(
    seqnames = loops_directional$anchor1_chr,
    ranges = IRanges(start = loops_directional$anchor1_start,
                     end = loops_directional$anchor1_end)
  )

  anchor2_gr <- GRanges(
    seqnames = loops_directional$anchor2_chr,
    ranges = IRanges(start = loops_directional$anchor2_start,
                     end = loops_directional$anchor2_end)
  )

  # Compute H3K27me3 overlaps
  loops_directional$anchor1_K27me3 <- countOverlaps(anchor1_gr, k27me3_gr) > 0
  loops_directional$anchor2_K27me3 <- countOverlaps(anchor2_gr, k27me3_gr) > 0

  # Compute Bivalent overlaps
  loops_directional$anchor1_Bivalent <- countOverlaps(anchor1_gr, bivalent_gr) > 0
  loops_directional$anchor2_Bivalent <- countOverlaps(anchor2_gr, bivalent_gr) > 0

  cat("  H3K27me3 overlaps:\n")
  cat(sprintf("    Anchor1: %d loops (%.1f%%)\n",
              sum(loops_directional$anchor1_K27me3),
              100 * mean(loops_directional$anchor1_K27me3)))
  cat(sprintf("    Anchor2: %d loops (%.1f%%)\n",
              sum(loops_directional$anchor2_K27me3),
              100 * mean(loops_directional$anchor2_K27me3)))

  cat("  Bivalent overlaps:\n")
  cat(sprintf("    Anchor1: %d loops (%.1f%%)\n",
              sum(loops_directional$anchor1_Bivalent),
              100 * mean(loops_directional$anchor1_Bivalent)))
  cat(sprintf("    Anchor2: %d loops (%.1f%%)\n\n",
              sum(loops_directional$anchor2_Bivalent),
              100 * mean(loops_directional$anchor2_Bivalent)))

  # ============================================================================
  # Step 4: Create filtered subsets
  # ============================================================================

  cat("=== Step 4: Creating Filtered Subsets ===\n")

  # Subset 1: K27me3-anchored (at least one anchor overlaps H3K27me3)
  loops_k27me3_anchored <- loops_directional %>%
    filter(anchor1_K27me3 | anchor2_K27me3)

  # Subset 2: K27me3-K27me3 (BOTH anchors overlap H3K27me3)
  loops_k27me3_both <- loops_directional %>%
    filter(anchor1_K27me3 & anchor2_K27me3)

  # Subset 3: Bivalent (at least one anchor overlaps bivalent)
  loops_bivalent <- loops_directional %>%
    filter(anchor1_Bivalent | anchor2_Bivalent)

  cat(sprintf("  K27me3-anchored (at least one anchor): %d loops\n", nrow(loops_k27me3_anchored)))
  cat(sprintf("    - Lost: %d, Gained: %d\n",
              sum(loops_k27me3_anchored$direction == "down_in_mutant"),
              sum(loops_k27me3_anchored$direction == "up_in_mutant")))

  cat(sprintf("  K27me3-K27me3 (both anchors): %d loops\n", nrow(loops_k27me3_both)))
  cat(sprintf("    - Lost: %d, Gained: %d\n",
              sum(loops_k27me3_both$direction == "down_in_mutant"),
              sum(loops_k27me3_both$direction == "up_in_mutant")))

  cat(sprintf("  Bivalent (at least one anchor): %d loops\n", nrow(loops_bivalent)))
  cat(sprintf("    - Lost: %d, Gained: %d\n\n",
              sum(loops_bivalent$direction == "down_in_mutant"),
              sum(loops_bivalent$direction == "up_in_mutant")))

  # ============================================================================
  # Step 5: Generate plots for each subset
  # ============================================================================

  cat("=== Step 5: Generating Plots ===\n")

  # K27me3-anchored plots
  cat("  Generating K27me3-anchored plots...\n")
  generate_cdf_plot(
    loops_k27me3_anchored,
    "K27me3-Anchored Loops",
    file.path(PLOT_DIR, "2E_cdf_k27me3_anchored"),  # Original: file.path(OUTPUT_DIR, "01_cdf_k27me3_anchored")
    COLORS
  )
  generate_density_plot(
    loops_k27me3_anchored,
    "K27me3-Anchored Loops",
    file.path(PLOT_DIR, "2E_density_k27me3_anchored"),  # Original: file.path(OUTPUT_DIR, "03_density_k27me3_anchored")
    COLORS
  )

  # K27me3-K27me3 plots
  cat("  Generating K27me3-K27me3 (both anchors) plots...\n")
  generate_cdf_plot(
    loops_k27me3_both,
    "K27me3-K27me3 Loops (Both Anchors)",
    file.path(PLOT_DIR, "2E_cdf_k27me3_both"),  # Original: file.path(OUTPUT_DIR, "01_cdf_k27me3_both")
    COLORS
  )
  generate_density_plot(
    loops_k27me3_both,
    "K27me3-K27me3 Loops (Both Anchors)",
    file.path(PLOT_DIR, "2E_density_k27me3_both"),  # Original: file.path(OUTPUT_DIR, "03_density_k27me3_both")
    COLORS
  )

  # Bivalent plots
  cat("  Generating Bivalent plots...\n")
  generate_cdf_plot(
    loops_bivalent,
    "Bivalent-Anchored Loops",
    file.path(PLOT_DIR, "2E_cdf_bivalent"),  # Original: file.path(OUTPUT_DIR, "01_cdf_bivalent")
    COLORS
  )
  generate_density_plot(
    loops_bivalent,
    "Bivalent-Anchored Loops",
    file.path(PLOT_DIR, "2E_density_bivalent"),  # Original: file.path(OUTPUT_DIR, "03_density_bivalent")
    COLORS
  )

  cat("\n")

  # ============================================================================
  # Step 6: Save summary statistics
  # ============================================================================

  cat("=== Step 6: Saving Summary ===\n")

  # Calculate statistics for each subset
  calc_subset_stats <- function(df, name) {
    if (nrow(df) < 10) {
      return(list(
        name = name,
        total = nrow(df),
        n_lost = sum(df$direction == "down_in_mutant"),
        n_gained = sum(df$direction == "up_in_mutant"),
        median_lost = NA,
        median_gained = NA,
        ks_pvalue = NA,
        wilcox_pvalue = NA
      ))
    }

    lost <- df %>% filter(direction == "down_in_mutant")
    gained <- df %>% filter(direction == "up_in_mutant")

    ks_p <- if (nrow(lost) >= 2 && nrow(gained) >= 2) {
      ks.test(lost$loop_distance, gained$loop_distance)$p.value
    } else NA

    wilcox_p <- if (nrow(lost) >= 2 && nrow(gained) >= 2) {
      wilcox.test(lost$loop_distance, gained$loop_distance)$p.value
    } else NA

    list(
      name = name,
      total = nrow(df),
      n_lost = nrow(lost),
      n_gained = nrow(gained),
      median_lost = median(lost$loop_distance),
      median_gained = median(gained$loop_distance),
      ks_pvalue = ks_p,
      wilcox_pvalue = wilcox_p
    )
  }

  stats_all <- calc_subset_stats(loops_directional, "All Loops")
  stats_k27me3_anchored <- calc_subset_stats(loops_k27me3_anchored, "K27me3-Anchored")
  stats_k27me3_both <- calc_subset_stats(loops_k27me3_both, "K27me3-K27me3 (Both)")
  stats_bivalent <- calc_subset_stats(loops_bivalent, "Bivalent")

  # Write summary file
  summary_file <- file.path(TSV_DIR, "2E_k27me3_filter_summary.txt")  # Original: file.path(OUTPUT_DIR, "filter_summary.txt")
  sink(summary_file)
  cat("=== H3K27me3-Filtered Loop Distance Analysis Summary ===\n")
  cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Timepoint:", timepoint, "\n\n")

  cat("Input file:", input_file, "\n")
  cat("ChIP-seq peaks:\n")
  cat("  H3K27me3:", peak_files$h3k27me3, "\n")
  cat("  Bivalent:", peak_files$bivalent, "\n\n")

  cat("=== Subset Statistics ===\n\n")

  for (stats in list(stats_all, stats_k27me3_anchored, stats_k27me3_both, stats_bivalent)) {
    cat(sprintf("--- %s ---\n", stats$name))
    cat(sprintf("  Total loops: %d\n", stats$total))
    cat(sprintf("  Lost in BAP1-KO: %d (%.1f%%)\n",
                stats$n_lost, 100 * stats$n_lost / max(stats$total, 1)))
    cat(sprintf("  Gained in BAP1-KO: %d (%.1f%%)\n",
                stats$n_gained, 100 * stats$n_gained / max(stats$total, 1)))
    if (!is.na(stats$median_lost)) {
      cat(sprintf("  Median distance (Lost): %.1f kb\n", stats$median_lost / 1000))
      cat(sprintf("  Median distance (Gained): %.1f kb\n", stats$median_gained / 1000))
      cat(sprintf("  Median ratio (Lost/Gained): %.2f\n", stats$median_lost / stats$median_gained))
      cat(sprintf("  KS test p-value: %.2e\n", stats$ks_pvalue))
      cat(sprintf("  Wilcoxon p-value: %.2e\n", stats$wilcox_pvalue))
    }
    cat("\n")
  }

  cat("=== Overlap Breakdown ===\n\n")
  cat(sprintf("Anchor1 H3K27me3+: %d / %d (%.1f%%)\n",
              sum(loops_directional$anchor1_K27me3),
              nrow(loops_directional),
              100 * mean(loops_directional$anchor1_K27me3)))
  cat(sprintf("Anchor2 H3K27me3+: %d / %d (%.1f%%)\n",
              sum(loops_directional$anchor2_K27me3),
              nrow(loops_directional),
              100 * mean(loops_directional$anchor2_K27me3)))
  cat(sprintf("Either anchor H3K27me3+: %d / %d (%.1f%%)\n",
              nrow(loops_k27me3_anchored),
              nrow(loops_directional),
              100 * nrow(loops_k27me3_anchored) / nrow(loops_directional)))
  cat(sprintf("Both anchors H3K27me3+: %d / %d (%.1f%%)\n\n",
              nrow(loops_k27me3_both),
              nrow(loops_directional),
              100 * nrow(loops_k27me3_both) / nrow(loops_directional)))

  cat(sprintf("Anchor1 Bivalent+: %d / %d (%.1f%%)\n",
              sum(loops_directional$anchor1_Bivalent),
              nrow(loops_directional),
              100 * mean(loops_directional$anchor1_Bivalent)))
  cat(sprintf("Anchor2 Bivalent+: %d / %d (%.1f%%)\n",
              sum(loops_directional$anchor2_Bivalent),
              nrow(loops_directional),
              100 * mean(loops_directional$anchor2_Bivalent)))
  cat(sprintf("Either anchor Bivalent+: %d / %d (%.1f%%)\n",
              nrow(loops_bivalent),
              nrow(loops_directional),
              100 * nrow(loops_bivalent) / nrow(loops_directional)))

  sink()

  cat("Saved:", summary_file, "\n")

  # ============================================================================
  # Completion
  # ============================================================================

  cat("\n=== Analysis Complete for", toupper(timepoint), "===\n")
  cat("Plot directory:", PLOT_DIR, "\n")
  cat("TSV directory:", TSV_DIR, "\n")
  cat("Files generated:\n")
  cat("  - 2E_cdf_k27me3_anchored.{pdf,svg,jpg}\n")
  cat("  - 2E_cdf_k27me3_both.{pdf,svg,jpg}\n")
  cat("  - 2E_cdf_bivalent.{pdf,svg,jpg}\n")
  cat("  - 2E_density_k27me3_anchored.{pdf,svg,jpg}\n")
  cat("  - 2E_density_k27me3_both.{pdf,svg,jpg}\n")
  cat("  - 2E_density_bivalent.{pdf,svg,jpg}\n")
  cat("  - 2E_k27me3_filter_summary.txt\n")
}

# ==============================================================================
# SECTION 6: MAIN EXECUTION
# ==============================================================================

for (tp in TIMEPOINTS_TO_RUN) {
  tryCatch({
    run_k27me3_filtered_analysis(tp)
  }, error = function(e) {
    cat(sprintf("\nERROR processing %s timepoint: %s\n", tp, e$message))
  })
}

cat("\n=== All timepoints complete ===\n")
cat("Output directories:\n")
cat(sprintf("  Plots: %s\n", PLOT_DIR))
cat(sprintf("  TSVs:  %s\n", TSV_DIR))
cat("\nEnd time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
