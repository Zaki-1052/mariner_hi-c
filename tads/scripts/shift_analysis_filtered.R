#!/usr/bin/env Rscript
# tads/scripts/shift_analysis_filtered.R
# Shift distance analysis visualizations using filtered shifted TAD boundaries
# Only includes shifts >= 10kb (1 bin at 10kb resolution)
#
# Usage:
#   Rscript scripts/shift_analysis_filtered.R [--timepoint TIMEPOINT]

cat("\n")
cat("========================================\n")
cat("Filtered Shift Analysis Visualizations\n")
cat("(>= 10kb shifts only)\n")
cat("========================================\n\n")

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
timepoint <- "late"  # Default

if (length(args) > 0) {
  for (i in seq_along(args)) {
    if (args[i] == "--timepoint" && i < length(args)) {
      timepoint <- args[i + 1]
    }
  }
}

cat(sprintf("Timepoint: %s\n\n", timepoint))

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(patchwork)
})

# Load multi-format output utility
tads_dir <- getwd()
base_dir <- dirname(tads_dir)
source(file.path(base_dir, "scripts/utils/multi_format_output.R"))

# Define paths
input_file <- file.path("results", timepoint, "final/shifted_boundaries_min1bin.bed")
output_dir <- file.path("results/visualizations", timepoint, "shift_analysis")

# Validate input
if (!file.exists(input_file)) {
  stop(sprintf("ERROR: Filtered shift file not found: %s\nRun the filtering step first.", input_file))
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Color schemes
enrichment_colors <- c(
  "Control" = "#4575b4",
  "Mutant" = "#d73027"
)

# Common theme
theme_tad <- function() {
  theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      legend.position = "right",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      plot.margin = margin(10, 10, 10, 10)
    )
}

# =============================================================================
# LOAD DATA
# =============================================================================

cat("Loading filtered shift data...\n")

shifted_df <- read.table(input_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Rename Enriched_In to cleaner labels
shifted_df$Direction <- ifelse(shifted_df$Enriched_In == "Control", "Control", "Mutant")

n_shifted <- nrow(shifted_df)
cat(sprintf("  Loaded %d shifted boundaries (>= 10kb shift)\n", n_shifted))

# Summary statistics
shift_stats <- shifted_df %>%
  summarise(
    n = n(),
    median_kb = median(shift_distance_kb, na.rm = TRUE),
    mean_kb = mean(shift_distance_kb, na.rm = TRUE),
    sd_kb = sd(shift_distance_kb, na.rm = TRUE),
    min_kb = min(shift_distance_kb, na.rm = TRUE),
    max_kb = max(shift_distance_kb, na.rm = TRUE)
  )

cat(sprintf("  Median shift: %.1f kb\n", shift_stats$median_kb))
cat(sprintf("  Mean shift: %.1f kb (SD: %.1f kb)\n", shift_stats$mean_kb, shift_stats$sd_kb))
cat(sprintf("  Range: %.1f - %.1f kb\n", shift_stats$min_kb, shift_stats$max_kb))

# By direction
shift_by_dir <- shifted_df %>%
  group_by(Direction) %>%
  summarise(
    n = n(),
    median_kb = median(shift_distance_kb, na.rm = TRUE),
    mean_kb = mean(shift_distance_kb, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nBy enrichment direction:\n")
for (i in 1:nrow(shift_by_dir)) {
  cat(sprintf("  %s: n=%d, median=%.1f kb, mean=%.1f kb\n",
              shift_by_dir$Direction[i], shift_by_dir$n[i],
              shift_by_dir$median_kb[i], shift_by_dir$mean_kb[i]))
}

# Binned distribution (round to nearest 10 to handle 10.001, 20.001 values)
shifted_df$shift_bin <- cut(round(shifted_df$shift_distance_kb),
                            breaks = c(0, 10, 20, 30, 40, Inf),
                            labels = c("10kb", "20kb", "30kb", "40kb", ">40kb"),
                            include.lowest = TRUE)
shift_bins <- shifted_df %>%
  count(shift_bin)

cat("\nShift distance bins:\n")
for (i in 1:nrow(shift_bins)) {
  cat(sprintf("  %s: %d (%.1f%%)\n", shift_bins$shift_bin[i], shift_bins$n[i],
              100 * shift_bins$n[i] / n_shifted))
}

# =============================================================================
# VISUALIZATIONS
# =============================================================================

cat("\n========================================\n")
cat("Generating Visualizations\n")
cat("========================================\n\n")

# 1. Shift Distance Histogram
cat("Creating shift distance histogram...\n")

p_hist <- ggplot(shifted_df, aes(x = shift_distance_kb, fill = Direction)) +
  geom_histogram(binwidth = 5, alpha = 0.7, position = "identity", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = enrichment_colors, name = "Enriched In") +
  labs(
    title = "Shift Distance Distribution (Filtered)",
    subtitle = sprintf("Shifted boundaries >= 10kb (n = %d) | Median = %.1f kb | %s",
                       n_shifted, shift_stats$median_kb, timepoint),
    x = "Shift Distance (kb)",
    y = "Count"
  ) +
  theme_tad()

save_multiformat_ggplot(p_hist,
                        file.path(output_dir, "shift_distance_histogram_filtered"),
                        width = 10, height = 6)

# 2. Shift Distance Violin Plot by Direction
cat("Creating shift distance violin plot...\n")

p_violin <- ggplot(shifted_df, aes(x = Direction, y = shift_distance_kb, fill = Direction)) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  geom_jitter(alpha = 0.4, width = 0.1, size = 1.5, color = "gray30") +
  scale_fill_manual(values = enrichment_colors) +
  labs(
    title = "Shift Distance by Enrichment Direction (Filtered)",
    subtitle = sprintf("n = %d boundaries | %s", n_shifted, timepoint),
    x = "Enriched In",
    y = "Shift Distance (kb)"
  ) +
  theme_tad() +
  theme(legend.position = "none")

save_multiformat_ggplot(p_violin,
                        file.path(output_dir, "shift_distance_violin_filtered"),
                        width = 8, height = 6)

# 3. Shift Distance vs Gap Score
cat("Creating shift vs Gap Score scatter...\n")

p_scatter <- ggplot(shifted_df, aes(x = shift_distance_kb, y = Gap_Score, color = Direction)) +
  geom_point(alpha = 0.7, size = 2.5) +
  geom_smooth(method = "lm", se = TRUE, linetype = "dashed", linewidth = 0.8) +
  scale_color_manual(values = enrichment_colors, name = "Enriched In") +
  labs(
    title = "Shift Distance vs Gap Score (Filtered)",
    subtitle = sprintf("Do larger shifts correlate with stronger differential signal? | %s", timepoint),
    x = "Shift Distance (kb)",
    y = "Gap Score"
  ) +
  theme_tad()

save_multiformat_ggplot(p_scatter,
                        file.path(output_dir, "shift_vs_gap_score_filtered"),
                        width = 10, height = 7)

# 4. Shift Bin Distribution Bar Plot
cat("Creating shift bin distribution plot...\n")

# Use pre-calculated shift_bin column
shift_bin_dir <- shifted_df %>%
  count(shift_bin, Direction)

p_bins <- ggplot(shift_bin_dir, aes(x = shift_bin, y = n, fill = Direction)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.3) +
  scale_fill_manual(values = enrichment_colors, name = "Enriched In") +
  labs(
    title = "Shift Distance Bins by Enrichment Direction",
    subtitle = sprintf("Filtered to >= 10kb shifts | n = %d | %s", n_shifted, timepoint),
    x = "Shift Distance Bin",
    y = "Count"
  ) +
  theme_tad()

save_multiformat_ggplot(p_bins,
                        file.path(output_dir, "shift_bins_by_direction_filtered"),
                        width = 10, height = 6)

# 5. Shift Distance per Chromosome
cat("Creating per-chromosome shift analysis...\n")

chr_shifts <- shifted_df %>%
  group_by(chr) %>%
  summarise(
    n = n(),
    median_kb = median(shift_distance_kb, na.rm = TRUE),
    mean_kb = mean(shift_distance_kb, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(chr = factor(chr, levels = paste0("chr", c(1:19, "X"))))

p_chr <- ggplot(chr_shifts, aes(x = chr, y = n)) +
  geom_bar(stat = "identity", fill = "#4575b4", color = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.0fkb", median_kb)), vjust = -0.3, size = 2.5, color = "gray30") +
  labs(
    title = "Shifted Boundaries per Chromosome (Filtered)",
    subtitle = sprintf("Numbers above bars = median shift distance | n = %d | %s", n_shifted, timepoint),
    x = "Chromosome",
    y = "Number of Shifted Boundaries"
  ) +
  theme_tad() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_multiformat_ggplot(p_chr,
                        file.path(output_dir, "shift_per_chromosome_filtered"),
                        width = 14, height = 6)

# 6. Combined Summary Panel
cat("Creating combined summary panel...\n")

p_combined <- (p_hist | p_violin) / (p_scatter | p_bins) +
  plot_annotation(
    title = sprintf("Filtered Shift Analysis Summary (%s)", timepoint),
    subtitle = sprintf("Shifted TAD boundaries with >= 10kb (1 bin) shift distance | n = %d", n_shifted),
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5)
    )
  )

save_multiformat_ggplot(p_combined,
                        file.path(output_dir, "shift_analysis_summary_filtered"),
                        width = 16, height = 12)

# =============================================================================
# SUMMARY STATISTICS FILE
# =============================================================================

cat("Writing summary statistics...\n")

summary_lines <- c(
  "==========================================",
  sprintf("FILTERED SHIFT ANALYSIS SUMMARY - %s", toupper(timepoint)),
  "==========================================",
  sprintf("Generated: %s", Sys.time()),
  sprintf("Filter: shift_distance >= 10kb (1 bin minimum)"),
  "",
  "--- OVERALL STATISTICS ---",
  sprintf("Total shifted boundaries (filtered): %d", n_shifted),
  sprintf("Median shift distance: %.1f kb", shift_stats$median_kb),
  sprintf("Mean shift distance: %.1f kb", shift_stats$mean_kb),
  sprintf("SD: %.1f kb", shift_stats$sd_kb),
  sprintf("Range: %.1f - %.1f kb", shift_stats$min_kb, shift_stats$max_kb),
  "",
  "--- BY ENRICHMENT DIRECTION ---"
)

for (i in 1:nrow(shift_by_dir)) {
  summary_lines <- c(summary_lines,
                     sprintf("  %s: n=%d (%.1f%%), median=%.1f kb, mean=%.1f kb",
                             shift_by_dir$Direction[i], shift_by_dir$n[i],
                             100 * shift_by_dir$n[i] / n_shifted,
                             shift_by_dir$median_kb[i], shift_by_dir$mean_kb[i]))
}

summary_lines <- c(summary_lines,
                   "",
                   "--- SHIFT DISTANCE BINS ---")
for (i in 1:nrow(shift_bins)) {
  summary_lines <- c(summary_lines,
                     sprintf("  %s: %d (%.1f%%)", shift_bins$shift_bin[i], shift_bins$n[i],
                             100 * shift_bins$n[i] / n_shifted))
}

summary_lines <- c(summary_lines,
                   "",
                   "--- COMPARISON TO UNFILTERED ---",
                   sprintf("Unfiltered shifted boundaries: see original file"),
                   sprintf("Filtered (>= 10kb): %d", n_shifted),
                   sprintf("Removed (< 10kb, likely noise): see original file"),
                   "",
                   "==========================================",
                   "OUTPUT FILES",
                   "==========================================",
                   "  - shift_distance_histogram_filtered.pdf",
                   "  - shift_distance_violin_filtered.pdf",
                   "  - shift_vs_gap_score_filtered.pdf",
                   "  - shift_bins_by_direction_filtered.pdf",
                   "  - shift_per_chromosome_filtered.pdf",
                   "  - shift_analysis_summary_filtered.pdf",
                   "  - shift_analysis_filtered_summary.txt",
                   "==========================================")

writeLines(summary_lines, file.path(output_dir, "shift_analysis_filtered_summary.txt"))

cat("\n")
cat(paste(summary_lines, collapse = "\n"))
cat("\n\n")

cat("========================================\n")
cat("FILTERED SHIFT ANALYSIS COMPLETE\n")
cat("========================================\n\n")
cat(sprintf("Output directory: %s\n", output_dir))
cat("All visualizations generated successfully!\n\n")
