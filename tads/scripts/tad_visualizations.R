#!/usr/bin/env Rscript
# tads/scripts/tad_visualizations.R
# TADCompare Visualization Analysis for Differential TAD Boundaries
# Author: Zakir Alibhai
# Date: 2025-01
#
# Purpose:
#   Generate publication-quality visualizations from TADCompare differential
#   TAD boundary analysis:
#   1. Overview plots (Gap Score, TAD Score scatter)
#   2. Boundary type classification
#   3. Shift distance analysis
#   4. Robustness analysis
#   5. Chromosome-level analysis
#   6. ChIP-seq integration (H3K27ac, H3K27me3, H3K4me1)
#   7. GO/KEGG enrichment analysis
#   8. syt1/nav3 locus focus (chr10)
#   9. Summary statistics
#
# Usage:
#   Rscript tads/scripts/tad_visualizations.R [--timepoint TIMEPOINT]
#
# Arguments:
#   --timepoint   Timepoint to analyze: "late" (default) or "early"
#
# Modular design: Functions accept timepoint parameter for easy extension

# =============================================================================
# SETUP
# =============================================================================

cat("\n")
cat("========================================\n")
cat("TADCompare Visualization Analysis\n")
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

cat(sprintf("Configuration:\n"))
cat(sprintf("  Timepoint: %s\n\n", timepoint))

# Load required libraries
cat("Loading required packages...\n")
suppressPackageStartupMessages({
  # Core visualization
  library(ggplot2)
  library(patchwork)
  library(RColorBrewer)

  # Data wrangling
  library(tidyverse)
  library(scales)

  # Genomic annotation
  library(GenomicRanges)
  library(rtracklayer)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(org.Mm.eg.db)

  # Enrichment analysis
  library(clusterProfiler)
  library(enrichplot)
})

cat("Packages loaded\n\n")

# Load multi-format output utility
# Note: path is relative to base_dir (mariner_hi-c/)
source(file.path(dirname(getwd()), "scripts/utils/multi_format_output.R"))

# Define base directories
# Script is run from tads/ directory via: cd tads && Rscript scripts/tad_visualizations.R
tads_dir <- getwd()  # Should be tads/
base_dir <- dirname(tads_dir)  # Parent = mariner_hi-c/

cat(sprintf("  Working directory: %s\n", tads_dir))
cat(sprintf("  Base directory: %s\n", base_dir))

# Input/output paths (relative to tads/)
# Use timepoint-specific input directory
input_file <- file.path("results", timepoint, "final/tadcompare_final_annotated.tsv")
output_base <- file.path("results/visualizations", timepoint)

# ChIP-seq peak paths (in parent directory peaks/beds/)
# Timepoint-specific Cerebellum peaks (standardized)
peaks_dir <- file.path(base_dir, "peaks", "beds")

if (timepoint == "late") {
  h3k27ac_path <- file.path(peaks_dir, "H3K27acCerebellumLate2.bed")
  h3k27me3_path <- file.path(peaks_dir, "H3K27me3CerebellumLate1.bed")
  h3k4me1_path <- file.path(peaks_dir, "H3K4me1CerebellumLate1.bed")
} else if (timepoint == "early") {
  h3k27ac_path <- file.path(peaks_dir, "H3K27acCerebellumEarly2.bed")
  h3k27me3_path <- file.path(peaks_dir, "H3K27me3CerebellumEarly1.bed")
  h3k4me1_path <- file.path(peaks_dir, "H3K4me1CerebellumEarly1.bed")
} else {
  stop(sprintf("Unknown timepoint: %s. Use 'late' or 'early'.", timepoint))
}

cat(sprintf("ChIP-seq files:\n"))
cat(sprintf("  H3K27ac:  %s\n", basename(h3k27ac_path)))
cat(sprintf("  H3K27me3: %s\n", basename(h3k27me3_path)))
cat(sprintf("  H3K4me1:  %s\n", if (!is.null(h3k4me1_path)) basename(h3k4me1_path) else "NOT AVAILABLE"))

# Create output directories
subdirs <- c("overview", "classification", "shift_analysis", "robustness",
             "chromosome", "chipseq", "enrichment", "syt1_nav3_focus", "summary")
for (subdir in subdirs) {
  dir.create(file.path(output_base, subdir), recursive = TRUE, showWarnings = FALSE)
}

cat(sprintf("Output directory: %s\n\n", output_base))

# =============================================================================
# COLOR SCHEMES
# =============================================================================

# Direction colors (Control vs Mutant enriched)
enrichment_colors <- c(
  "Matrix 1" = "#4575b4",      # Control = Blue
  "Matrix 2" = "#d73027"       # Mutant = Red
)

# Boundary type colors (6 categories)
type_colors <- c(
  "Non-Differential" = "#999999",  # Gray
  "Strength Change" = "#ff7f00",   # Orange
  "Complex" = "#984ea3",           # Purple
  "Shifted" = "#e41a1c",           # Red
  "Merge" = "#377eb8",             # Blue
  "Split" = "#4daf4a"              # Green
)

# Robustness colors
robustness_colors <- c(
  "Both" = "#2166ac",             # Dark blue
  "Control_only" = "#92c5de",     # Light blue
  "Mutant_only" = "#f4a582",      # Light red
  "Neither" = "#d9d9d9"           # Light gray
)

# Anchor classification colors
anchor_colors <- c(
  "Promoter" = "#e41a1c",
  "Active_Enhancer" = "#ff7f00",
  "Polycomb_Domain" = "#984ea3",
  "Poised_Enhancer" = "#4daf4a",
  "Other" = "#999999"
)

# Common theme for all plots
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

cat("========================================\n")
cat("Loading Input Data\n")
cat("========================================\n\n")

# Validate input file exists
if (!file.exists(input_file)) {
  stop(sprintf("ERROR: Input file not found: %s", input_file))
}

# Load TADCompare results
tad_df <- read.table(input_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
cat(sprintf("Loaded TADCompare results: %d boundaries\n", nrow(tad_df)))

# Create derived columns
tad_df <- tad_df %>%
  mutate(
    TAD_Score_diff = TAD_Score2 - TAD_Score1,
    abs_Gap_Score = abs(Gap_Score),
    # Ensure consistent factor levels
    Type = factor(Type, levels = names(type_colors)),
    Differential = factor(Differential, levels = c("Non-Differential", "Differential", "Shifted")),
    robustness = factor(robustness, levels = c("Both", "Control_only", "Mutant_only", "Neither"))
  )

# Summary statistics
n_total <- nrow(tad_df)
n_diff <- sum(tad_df$Differential %in% c("Differential", "Shifted"), na.rm = TRUE)
n_shifted <- sum(tad_df$Type == "Shifted", na.rm = TRUE)

cat(sprintf("  Total boundaries: %d\n", n_total))
cat(sprintf("  Differential: %d (%.1f%%)\n", n_diff, 100 * n_diff / n_total))
cat(sprintf("  Shifted: %d (%.1f%%)\n", n_shifted, 100 * n_shifted / n_total))
cat("\n")

# =============================================================================
# SECTION 1: OVERVIEW PLOTS
# =============================================================================

cat("========================================\n")
cat("SECTION 1: Overview Plots\n")
cat("========================================\n\n")

# 1.1 Gap Score Distribution by Differential Status
cat("Creating Gap Score distribution plot...\n")

# Create simplified differential category for plotting
tad_df$diff_category <- case_when(
  tad_df$Type == "Shifted" ~ "Shifted",
  tad_df$Differential == "Differential" ~ "Differential",
  TRUE ~ "Non-Differential"
)
tad_df$diff_category <- factor(tad_df$diff_category,
                                levels = c("Non-Differential", "Differential", "Shifted"))

p_gap_dist <- ggplot(tad_df, aes(x = diff_category, y = Gap_Score, fill = diff_category)) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  scale_fill_manual(values = c("Non-Differential" = "#999999",
                               "Differential" = "#ff7f00",
                               "Shifted" = "#e41a1c")) +
  labs(
    title = "Gap Score Distribution by Differential Status",
    subtitle = sprintf("n = %d boundaries", n_total),
    x = "Differential Status",
    y = "Gap Score (Z-score)"
  ) +
  theme_tad() +
  theme(legend.position = "none")

save_multiformat_ggplot(p_gap_dist,
                        file.path(output_base, "overview", "gap_score_distribution"),
                        width = 8, height = 6)

# 1.2 TAD Score Scatter Plot (MA-plot style)
cat("Creating TAD Score scatter plot...\n")

p_scatter <- ggplot(tad_df, aes(x = TAD_Score1, y = TAD_Score2, color = Type)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  scale_color_manual(values = type_colors, name = "Boundary Type") +
  labs(
    title = "TAD Boundary Strength: Control vs Mutant",
    subtitle = sprintf("n = %d boundaries | %d differential | %d shifted", n_total, n_diff, n_shifted),
    x = "TAD Score (Control)",
    y = "TAD Score (Mutant)"
  ) +
  coord_fixed(ratio = 1) +
  theme_tad()

save_multiformat_ggplot(p_scatter,
                        file.path(output_base, "overview", "tad_score_scatter"),
                        width = 10, height = 8)

# 1.3 Gap Score vs TAD Score Difference (Volcano-like)
cat("Creating differential landscape plot...\n")

p_landscape <- ggplot(tad_df %>% filter(!is.na(Enriched_In)),
                      aes(x = TAD_Score_diff, y = abs_Gap_Score, color = Enriched_In)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  scale_color_manual(values = enrichment_colors, name = "Enriched In",
                     labels = c("Matrix 1" = "Control", "Matrix 2" = "Mutant")) +
  labs(
    title = "Differential TAD Boundary Landscape",
    subtitle = "Boundary strength difference vs differential score magnitude",
    x = "TAD Score Difference (Mutant - Control)",
    y = "|Gap Score|"
  ) +
  theme_tad()

save_multiformat_ggplot(p_landscape,
                        file.path(output_base, "overview", "differential_landscape"),
                        width = 10, height = 7)

cat("Section 1 complete: Overview plots generated\n\n")

# =============================================================================
# SECTION 2: BOUNDARY TYPE CLASSIFICATION
# =============================================================================

cat("========================================\n")
cat("SECTION 2: Boundary Type Classification\n")
cat("========================================\n\n")

# 2.1 Boundary Type Pie Charts
cat("Creating boundary type pie charts...\n")

# Summary for all boundaries
type_summary_all <- tad_df %>%
  count(Type) %>%
  mutate(
    percentage = 100 * n / sum(n),
    label = sprintf("%s\n%.1f%%", Type, percentage)
  )

# Summary for differential only
type_summary_diff <- tad_df %>%
  filter(Differential %in% c("Differential", "Shifted") | Type == "Shifted") %>%
  count(Type) %>%
  mutate(
    percentage = 100 * n / sum(n),
    label = sprintf("%s\n%.1f%%", Type, percentage)
  )

# Pie chart - All boundaries
p_pie_all <- ggplot(type_summary_all, aes(x = "", y = n, fill = Type)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = type_colors) +
  labs(title = "All Boundaries", subtitle = sprintf("n = %d", n_total)) +
  theme_void() +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    legend.position = "none"
  ) +
  geom_text(aes(label = ifelse(percentage > 3, sprintf("%.1f%%", percentage), "")),
            position = position_stack(vjust = 0.5), size = 3)

# Pie chart - Differential only
p_pie_diff <- ggplot(type_summary_diff, aes(x = "", y = n, fill = Type)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = type_colors) +
  labs(title = "Differential Boundaries", subtitle = sprintf("n = %d", sum(type_summary_diff$n))) +
  theme_void() +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    legend.position = "right",
    legend.title = element_blank()
  ) +
  geom_text(aes(label = ifelse(percentage > 5, sprintf("%.1f%%", percentage), "")),
            position = position_stack(vjust = 0.5), size = 3)

p_pie_combined <- p_pie_all | p_pie_diff
save_multiformat_ggplot(p_pie_combined,
                        file.path(output_base, "classification", "boundary_type_pie"),
                        width = 12, height = 5)

# 2.2 Boundary Type by Chromosome
cat("Creating boundary type by chromosome plot...\n")

chr_type_summary <- tad_df %>%
  count(chr, Type) %>%
  group_by(chr) %>%
  mutate(percentage = 100 * n / sum(n)) %>%
  ungroup() %>%
  mutate(chr = factor(chr, levels = paste0("chr", c(1:19, "X"))))

p_chr_type <- ggplot(chr_type_summary, aes(x = chr, y = n, fill = Type)) +
  geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.2) +
  scale_fill_manual(values = type_colors, name = "Boundary Type") +
  labs(
    title = "Boundary Type Distribution by Chromosome",
    x = "Chromosome",
    y = "Number of Boundaries"
  ) +
  theme_tad() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_multiformat_ggplot(p_chr_type,
                        file.path(output_base, "classification", "boundary_type_by_chromosome"),
                        width = 14, height = 6)

# 2.3 Enrichment Direction by Type
cat("Creating enrichment direction by type plot...\n")

enrichment_by_type <- tad_df %>%
  filter(!is.na(Enriched_In) & Type != "Non-Differential") %>%
  count(Type, Enriched_In) %>%
  mutate(Enriched_In = factor(Enriched_In, levels = c("Matrix 1", "Matrix 2")))

p_enrich_type <- ggplot(enrichment_by_type, aes(x = Type, y = n, fill = Enriched_In)) +
  geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.3) +
  scale_fill_manual(values = enrichment_colors, name = "Enriched In",
                    labels = c("Matrix 1" = "Control", "Matrix 2" = "Mutant")) +
  labs(
    title = "Enrichment Direction by Boundary Type",
    subtitle = "Differential boundaries only",
    x = "Boundary Type",
    y = "Count"
  ) +
  theme_tad() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_multiformat_ggplot(p_enrich_type,
                        file.path(output_base, "classification", "enrichment_direction_by_type"),
                        width = 10, height = 6)

cat("Section 2 complete: Boundary type classification plots generated\n\n")

# =============================================================================
# SECTION 3: SHIFT DISTANCE ANALYSIS
# =============================================================================

cat("========================================\n")
cat("SECTION 3: Shift Distance Analysis\n")
cat("========================================\n\n")

# Filter to shifted boundaries only
shifted_df <- tad_df %>%
  filter(Type == "Shifted" & !is.na(shift_distance_kb))

n_shifted_valid <- nrow(shifted_df)
cat(sprintf("Shifted boundaries with distance data: %d\n", n_shifted_valid))

if (n_shifted_valid > 0) {
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
  cat(sprintf("  Mean shift: %.1f kb\n", shift_stats$mean_kb))
  cat(sprintf("  Range: %.1f - %.1f kb\n", shift_stats$min_kb, shift_stats$max_kb))

  # 3.1 Shift Distance Histogram
  cat("Creating shift distance histogram...\n")

  p_shift_hist <- ggplot(shifted_df, aes(x = shift_distance_kb, fill = Enriched_In)) +
    geom_histogram(bins = 30, alpha = 0.7, position = "identity", color = "black", linewidth = 0.2) +
    scale_fill_manual(values = enrichment_colors, name = "Enriched In",
                      labels = c("Matrix 1" = "Control", "Matrix 2" = "Mutant")) +
    labs(
      title = "Shift Distance Distribution",
      subtitle = sprintf("Shifted boundaries (n = %d) | Median = %.1f kb", n_shifted_valid, shift_stats$median_kb),
      x = "Shift Distance (kb)",
      y = "Count"
    ) +
    theme_tad()

  save_multiformat_ggplot(p_shift_hist,
                          file.path(output_base, "shift_analysis", "shift_distance_histogram"),
                          width = 10, height = 6)

  # 3.2 Shift Distance Violin Plot
  cat("Creating shift distance violin plot...\n")

  p_shift_violin <- ggplot(shifted_df %>% filter(!is.na(Enriched_In)),
                           aes(x = Enriched_In, y = shift_distance_kb, fill = Enriched_In)) +
    geom_violin(alpha = 0.7, trim = FALSE) +
    geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
    geom_jitter(alpha = 0.4, width = 0.1, size = 1.5) +
    scale_fill_manual(values = enrichment_colors) +
    scale_x_discrete(labels = c("Matrix 1" = "Control", "Matrix 2" = "Mutant")) +
    labs(
      title = "Shift Distance by Enrichment Direction",
      x = "Enriched In",
      y = "Shift Distance (kb)"
    ) +
    theme_tad() +
    theme(legend.position = "none")

  save_multiformat_ggplot(p_shift_violin,
                          file.path(output_base, "shift_analysis", "shift_distance_violin"),
                          width = 8, height = 6)

  # 3.3 Shift Distance vs Gap Score
  cat("Creating shift vs Gap Score scatter...\n")

  p_shift_gap <- ggplot(shifted_df, aes(x = shift_distance_kb, y = Gap_Score, color = Enriched_In)) +
    geom_point(alpha = 0.7, size = 2.5) +
    geom_smooth(method = "lm", se = TRUE, linetype = "dashed", linewidth = 0.8) +
    scale_color_manual(values = enrichment_colors, name = "Enriched In",
                       labels = c("Matrix 1" = "Control", "Matrix 2" = "Mutant")) +
    labs(
      title = "Shift Distance vs Gap Score",
      subtitle = "Do larger shifts correlate with stronger differential signal?",
      x = "Shift Distance (kb)",
      y = "Gap Score"
    ) +
    theme_tad()

  save_multiformat_ggplot(p_shift_gap,
                          file.path(output_base, "shift_analysis", "shift_vs_gap_score"),
                          width = 10, height = 7)

} else {
  cat("  No shifted boundaries with valid distance data. Skipping shift analysis plots.\n")
}

cat("\nSection 3 complete: Shift distance analysis plots generated\n\n")

# =============================================================================
# SECTION 4: ROBUSTNESS ANALYSIS
# =============================================================================

cat("========================================\n")
cat("SECTION 4: Robustness Analysis\n")
cat("========================================\n\n")

# 4.1 Robustness by Differential Status Heatmap
cat("Creating robustness x differential heatmap...\n")

robust_diff_summary <- tad_df %>%
  count(diff_category, robustness) %>%
  complete(diff_category, robustness, fill = list(n = 0))

p_robust_heatmap <- ggplot(robust_diff_summary, aes(x = robustness, y = diff_category, fill = n)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = n), size = 5, fontface = "bold") +
  scale_fill_gradient(low = "white", high = "#2166ac", name = "Count") +
  labs(
    title = "Robustness vs Differential Status",
    subtitle = "Cross-tabulation of boundary classifications",
    x = "Robustness Category",
    y = "Differential Status"
  ) +
  theme_tad() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

save_multiformat_ggplot(p_robust_heatmap,
                        file.path(output_base, "robustness", "robustness_differential_heatmap"),
                        width = 10, height = 6)

# 4.2 Robustness by Boundary Type
cat("Creating robustness by type plot...\n")

robust_type_summary <- tad_df %>%
  filter(Type != "Non-Differential") %>%
  count(Type, robustness) %>%
  group_by(Type) %>%
  mutate(percentage = 100 * n / sum(n)) %>%
  ungroup()

p_robust_type <- ggplot(robust_type_summary, aes(x = Type, y = percentage, fill = robustness)) +
  geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.3) +
  scale_fill_manual(values = robustness_colors, name = "Robustness") +
  labs(
    title = "Robustness by Boundary Type",
    subtitle = "Differential boundaries only",
    x = "Boundary Type",
    y = "Percentage (%)"
  ) +
  theme_tad() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_multiformat_ggplot(p_robust_type,
                        file.path(output_base, "robustness", "robustness_by_type"),
                        width = 10, height = 6)

cat("Section 4 complete: Robustness analysis plots generated\n\n")

# =============================================================================
# SECTION 5: CHROMOSOME-LEVEL ANALYSIS
# =============================================================================

cat("========================================\n")
cat("SECTION 5: Chromosome-Level Analysis\n")
cat("========================================\n\n")

# 5.1 Per-Chromosome Differential Percentage
cat("Creating per-chromosome differential percentage plot...\n")

chr_diff_summary <- tad_df %>%
  group_by(chr) %>%
  summarise(
    total = n(),
    n_diff = sum(Differential %in% c("Differential", "Shifted") | Type == "Shifted", na.rm = TRUE),
    pct_diff = 100 * n_diff / total
  ) %>%
  mutate(chr = factor(chr, levels = paste0("chr", c(1:19, "X"))))

overall_pct_diff <- 100 * n_diff / n_total

p_chr_diff <- ggplot(chr_diff_summary, aes(x = chr, y = pct_diff)) +
  geom_bar(stat = "identity", fill = "#4575b4", color = "black", linewidth = 0.3) +
  geom_hline(yintercept = overall_pct_diff, linetype = "dashed", color = "#d73027", linewidth = 1) +
  annotate("text", x = 1, y = overall_pct_diff + 1.5,
           label = sprintf("Overall: %.1f%%", overall_pct_diff),
           hjust = 0, color = "#d73027", fontface = "bold", size = 3.5) +
  labs(
    title = "Differential Boundary Percentage by Chromosome",
    subtitle = sprintf("Red line = overall differential rate (%.1f%%)", overall_pct_diff),
    x = "Chromosome",
    y = "% Differential Boundaries"
  ) +
  theme_tad() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_multiformat_ggplot(p_chr_diff,
                        file.path(output_base, "chromosome", "per_chromosome_differential"),
                        width = 12, height = 6)

cat("Section 5 complete: Chromosome-level analysis plots generated\n\n")

# =============================================================================
# SECTION 6: ChIP-seq INTEGRATION
# =============================================================================

cat("========================================\n")
cat("SECTION 6: ChIP-seq Integration\n")
cat("========================================\n\n")

# Load TxDb for TSS annotation
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
tss <- resize(genes(txdb), width = 1, fix = "start")
cat(sprintf("Loaded %d TSS positions from mm10\n", length(tss)))

# Function to load ChIP-seq peaks
load_peaks <- function(peak_file, peak_type) {
  if (!file.exists(peak_file)) {
    warning(sprintf("Peak file not found: %s", peak_file))
    return(NULL)
  }

  cat(sprintf("  Loading %s: %s\n", peak_type, basename(peak_file)))

  peak_df <- tryCatch({
    read.table(peak_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  }, error = function(e) {
    warning(sprintf("Could not read peak file: %s", e$message))
    return(NULL)
  })

  if (is.null(peak_df) || nrow(peak_df) == 0) {
    return(NULL)
  }

  # Create GRanges (BED format: 0-based start)
  peaks_gr <- GRanges(
    seqnames = peak_df$V1,
    ranges = IRanges(start = peak_df$V2 + 1, end = peak_df$V3)
  )

  cat(sprintf("    Loaded %d peaks\n", length(peaks_gr)))
  return(peaks_gr)
}

# Load ChIP-seq peaks
h3k27ac_peaks <- load_peaks(h3k27ac_path, "H3K27ac")
h3k27me3_peaks <- load_peaks(h3k27me3_path, "H3K27me3")
# H3K4me1 not available for early timepoint
if (!is.null(h3k4me1_path)) {
  h3k4me1_peaks <- load_peaks(h3k4me1_path, "H3K4me1")
} else {
  cat("  H3K4me1: Not available for this timepoint\n")
  h3k4me1_peaks <- NULL
}

# Create GRanges for TAD boundaries (use 25kb window around boundary)
boundary_gr <- GRanges(
  seqnames = tad_df$chr,
  ranges = IRanges(start = pmax(1, tad_df$Boundary - 12500),
                   end = tad_df$Boundary + 12500)
)

# Compute distance to nearest TSS
cat("Computing distance to nearest TSS...\n")
nearest_tss <- distanceToNearest(boundary_gr, tss)

tad_df$distance_to_tss <- NA
tad_df$nearest_gene_id <- NA

if (length(nearest_tss) > 0) {
  tad_df$distance_to_tss[queryHits(nearest_tss)] <- mcols(nearest_tss)$distance
  tad_df$nearest_gene_id[queryHits(nearest_tss)] <- names(tss)[subjectHits(nearest_tss)]
}

# Compute ChIP-seq overlaps
cat("Computing ChIP-seq overlaps...\n")

tad_df$h3k27ac_overlap <- FALSE
tad_df$h3k27me3_overlap <- FALSE
tad_df$h3k4me1_overlap <- FALSE

if (!is.null(h3k27ac_peaks)) {
  tad_df$h3k27ac_overlap <- countOverlaps(boundary_gr, h3k27ac_peaks) > 0
  cat(sprintf("  H3K27ac overlaps: %d (%.1f%%)\n",
              sum(tad_df$h3k27ac_overlap), 100 * mean(tad_df$h3k27ac_overlap)))
}

if (!is.null(h3k27me3_peaks)) {
  tad_df$h3k27me3_overlap <- countOverlaps(boundary_gr, h3k27me3_peaks) > 0
  cat(sprintf("  H3K27me3 overlaps: %d (%.1f%%)\n",
              sum(tad_df$h3k27me3_overlap), 100 * mean(tad_df$h3k27me3_overlap)))
}

if (!is.null(h3k4me1_peaks)) {
  tad_df$h3k4me1_overlap <- countOverlaps(boundary_gr, h3k4me1_peaks) > 0
  cat(sprintf("  H3K4me1 overlaps: %d (%.1f%%)\n",
              sum(tad_df$h3k4me1_overlap), 100 * mean(tad_df$h3k4me1_overlap)))
}

# Classify boundary anchors
cat("Classifying boundary types...\n")
if (is.null(h3k4me1_peaks)) {
  cat("  Note: Poised_Enhancer classification unavailable (no H3K4me1 data)\n")
}
tss_threshold <- 2000  # 2kb

tad_df$anchor_type <- case_when(
  # Promoter: H3K27ac+ AND within 2kb of TSS
  tad_df$h3k27ac_overlap & tad_df$distance_to_tss <= tss_threshold ~ "Promoter",
  # Active Enhancer: H3K27ac+ AND > 2kb from TSS
  tad_df$h3k27ac_overlap & tad_df$distance_to_tss > tss_threshold ~ "Active_Enhancer",
  # Polycomb Domain: H3K27me3+
  tad_df$h3k27me3_overlap ~ "Polycomb_Domain",
  # Poised Enhancer: H3K4me1+ (no H3K27ac) AND > 2kb from TSS
  tad_df$h3k4me1_overlap & !tad_df$h3k27ac_overlap &
    tad_df$distance_to_tss > tss_threshold ~ "Poised_Enhancer",
  # Other: No marks
  TRUE ~ "Other"
)

# Print summary
anchor_summary <- table(tad_df$anchor_type)
cat("\nAnchor type distribution:\n")
for (atype in names(anchor_summary)) {
  cat(sprintf("  %s: %d (%.1f%%)\n", atype, anchor_summary[atype],
              100 * anchor_summary[atype] / nrow(tad_df)))
}

# 6.1 Anchor Classification Plot
cat("\nCreating anchor classification plot...\n")

anchor_by_diff <- tad_df %>%
  count(diff_category, anchor_type) %>%
  group_by(diff_category) %>%
  mutate(percentage = 100 * n / sum(n)) %>%
  ungroup()

p_anchor_class <- ggplot(anchor_by_diff, aes(x = diff_category, y = percentage, fill = anchor_type)) +
  geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.3) +
  scale_fill_manual(values = anchor_colors, name = "Anchor Type") +
  labs(
    title = "Anchor Classification by Differential Status",
    subtitle = "Based on ChIP-seq overlaps and TSS proximity",
    x = "Differential Status",
    y = "Percentage (%)"
  ) +
  theme_tad()

save_multiformat_ggplot(p_anchor_class,
                        file.path(output_base, "chipseq", "anchor_classification"),
                        width = 10, height = 6)

# 6.2 ChIP-seq Overlap Heatmap by Boundary Type
cat("Creating ChIP-seq overlap heatmap...\n")

# Build mark list based on available data
chip_marks <- c("H3K27ac", "H3K27me3")
if (!is.null(h3k4me1_peaks)) chip_marks <- c(chip_marks, "H3K4me1")

chip_overlap_summary <- tad_df %>%
  filter(Type != "Non-Differential") %>%
  group_by(Type) %>%
  summarise(
    H3K27ac = 100 * mean(h3k27ac_overlap, na.rm = TRUE),
    H3K27me3 = 100 * mean(h3k27me3_overlap, na.rm = TRUE),
    H3K4me1 = 100 * mean(h3k4me1_overlap, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = all_of(chip_marks),
               names_to = "Mark", values_to = "Percentage")

chip_subtitle <- if (is.null(h3k4me1_peaks)) {
  sprintf("Differential boundaries only | %s timepoint (no H3K4me1)", timepoint)
} else {
  "Differential boundaries only"
}

p_chip_heatmap <- ggplot(chip_overlap_summary, aes(x = Mark, y = Type, fill = Percentage)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%.1f%%", Percentage)), size = 4, fontface = "bold") +
  scale_fill_gradient(low = "white", high = "#2166ac", name = "% Overlap") +
  labs(
    title = "ChIP-seq Mark Overlap by Boundary Type",
    subtitle = chip_subtitle,
    x = "Histone Mark",
    y = "Boundary Type"
  ) +
  theme_tad() +
  theme(panel.grid = element_blank())

save_multiformat_ggplot(p_chip_heatmap,
                        file.path(output_base, "chipseq", "chipseq_overlap_heatmap"),
                        width = 8, height = 6)

# 6.3 ChIP-seq Overlap by Enrichment Direction
cat("Creating ChIP-seq by enrichment direction plot...\n")

chip_by_enrich <- tad_df %>%
  filter(!is.na(Enriched_In) & Type != "Non-Differential") %>%
  group_by(Enriched_In) %>%
  summarise(
    H3K27ac = 100 * mean(h3k27ac_overlap, na.rm = TRUE),
    H3K27me3 = 100 * mean(h3k27me3_overlap, na.rm = TRUE),
    H3K4me1 = 100 * mean(h3k4me1_overlap, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = all_of(chip_marks),
               names_to = "Mark", values_to = "Percentage")

p_chip_enrich <- ggplot(chip_by_enrich, aes(x = Mark, y = Percentage, fill = Enriched_In)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.3) +
  scale_fill_manual(values = enrichment_colors, name = "Enriched In",
                    labels = c("Matrix 1" = "Control", "Matrix 2" = "Mutant")) +
  labs(
    title = "ChIP-seq Overlap by Enrichment Direction",
    subtitle = chip_subtitle,
    x = "Histone Mark",
    y = "% of Boundaries with Overlap"
  ) +
  theme_tad()

save_multiformat_ggplot(p_chip_enrich,
                        file.path(output_base, "chipseq", "chipseq_by_enrichment_direction"),
                        width = 10, height = 6)

cat("Section 6 complete: ChIP-seq integration plots generated\n\n")

# =============================================================================
# SECTION 7: GENE ENRICHMENT ANALYSIS
# =============================================================================

cat("========================================\n")
cat("SECTION 7: Gene Enrichment Analysis\n")
cat("========================================\n\n")

# Function to get genes near boundaries
get_nearby_genes <- function(df_subset, max_dist = 10000) {
  if (nrow(df_subset) == 0) return(character(0))

  boundary_gr <- GRanges(
    seqnames = df_subset$chr,
    ranges = IRanges(start = pmax(1, df_subset$Boundary - max_dist),
                     end = df_subset$Boundary + max_dist)
  )

  genes_txdb <- genes(txdb)
  overlaps <- findOverlaps(boundary_gr, genes_txdb, maxgap = 0)

  gene_ids <- names(genes_txdb)[subjectHits(overlaps)]
  return(unique(gene_ids[!is.na(gene_ids)]))
}

# Get genes for control-enriched and mutant-enriched differential boundaries
diff_df <- tad_df %>%
  filter(Differential %in% c("Differential", "Shifted") | Type == "Shifted")

ctrl_enriched <- diff_df %>% filter(Enriched_In == "Matrix 1")
mut_enriched <- diff_df %>% filter(Enriched_In == "Matrix 2")

cat(sprintf("Control-enriched differential boundaries: %d\n", nrow(ctrl_enriched)))
cat(sprintf("Mutant-enriched differential boundaries: %d\n", nrow(mut_enriched)))

ctrl_genes <- get_nearby_genes(ctrl_enriched)
mut_genes <- get_nearby_genes(mut_enriched)

cat(sprintf("  Genes near control-enriched boundaries: %d\n", length(ctrl_genes)))
cat(sprintf("  Genes near mutant-enriched boundaries: %d\n", length(mut_genes)))

# Save gene lists
gene_df <- data.frame(
  gene_id = c(ctrl_genes, mut_genes),
  enriched_in = c(rep("Control", length(ctrl_genes)), rep("Mutant", length(mut_genes))),
  stringsAsFactors = FALSE
)

# Add gene symbols
if (nrow(gene_df) > 0) {
  symbol_map <- tryCatch({
    AnnotationDbi::select(org.Mm.eg.db,
                          keys = unique(gene_df$gene_id),
                          columns = "SYMBOL",
                          keytype = "ENTREZID")
  }, error = function(e) {
    data.frame(ENTREZID = character(), SYMBOL = character())
  })
  gene_df$symbol <- symbol_map$SYMBOL[match(gene_df$gene_id, symbol_map$ENTREZID)]

  write.table(gene_df, file.path(output_base, "enrichment", "boundary_genes.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("  Saved: boundary_genes.tsv\n")
}

# Create gene list for enrichment
gene_list <- list()
if (length(ctrl_genes) >= 10) gene_list$Control_enriched <- ctrl_genes
if (length(mut_genes) >= 10) gene_list$Mutant_enriched <- mut_genes

if (length(gene_list) > 0) {

  # GO Biological Process
  cat("\nRunning GO Biological Process enrichment...\n")
  go_bp <- tryCatch({
    compareCluster(
      geneCluster = gene_list,
      fun = "enrichGO",
      OrgDb = org.Mm.eg.db,
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05
    )
  }, error = function(e) {
    cat(sprintf("  Error: %s\n", e$message))
    return(NULL)
  })

  if (!is.null(go_bp) && nrow(go_bp@compareClusterResult) > 0) {
    p_go_bp <- dotplot(go_bp, showCategory = 15) +
      labs(title = "GO Biological Process Enrichment") +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    save_multiformat_ggplot(p_go_bp,
                            file.path(output_base, "enrichment", "go_bp_dotplot"),
                            width = 12, height = 10)
  } else {
    cat("  No significant GO BP terms found\n")
  }

  # GO Cellular Component
  cat("Running GO Cellular Component enrichment...\n")
  go_cc <- tryCatch({
    compareCluster(
      geneCluster = gene_list,
      fun = "enrichGO",
      OrgDb = org.Mm.eg.db,
      ont = "CC",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05
    )
  }, error = function(e) {
    cat(sprintf("  Error: %s\n", e$message))
    return(NULL)
  })

  if (!is.null(go_cc) && nrow(go_cc@compareClusterResult) > 0) {
    p_go_cc <- dotplot(go_cc, showCategory = 15) +
      labs(title = "GO Cellular Component Enrichment") +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    save_multiformat_ggplot(p_go_cc,
                            file.path(output_base, "enrichment", "go_cc_dotplot"),
                            width = 10, height = 8)
  } else {
    cat("  No significant GO CC terms found\n")
  }

  # GO Molecular Function
  cat("Running GO Molecular Function enrichment...\n")
  go_mf <- tryCatch({
    compareCluster(
      geneCluster = gene_list,
      fun = "enrichGO",
      OrgDb = org.Mm.eg.db,
      ont = "MF",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05
    )
  }, error = function(e) {
    cat(sprintf("  Error: %s\n", e$message))
    return(NULL)
  })

  if (!is.null(go_mf) && nrow(go_mf@compareClusterResult) > 0) {
    p_go_mf <- dotplot(go_mf, showCategory = 15) +
      labs(title = "GO Molecular Function Enrichment") +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    save_multiformat_ggplot(p_go_mf,
                            file.path(output_base, "enrichment", "go_mf_dotplot"),
                            width = 10, height = 8)
  } else {
    cat("  No significant GO MF terms found\n")
  }

  # KEGG pathways
  cat("Running KEGG pathway enrichment...\n")
  kegg <- tryCatch({
    compareCluster(
      geneCluster = gene_list,
      fun = "enrichKEGG",
      organism = "mmu",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05
    )
  }, error = function(e) {
    cat(sprintf("  Error: %s\n", e$message))
    return(NULL)
  })

  if (!is.null(kegg) && nrow(kegg@compareClusterResult) > 0) {
    p_kegg <- dotplot(kegg, showCategory = 15) +
      labs(title = "KEGG Pathway Enrichment") +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    save_multiformat_ggplot(p_kegg,
                            file.path(output_base, "enrichment", "kegg_dotplot"),
                            width = 12, height = 10)
  } else {
    cat("  No significant KEGG pathways found\n")
  }

} else {
  cat("Insufficient genes for enrichment analysis (need >= 10 per group)\n")
}

cat("\nSection 7 complete: Gene enrichment analysis done\n\n")

# =============================================================================
# SECTION 8: syt1/nav3 LOCUS FOCUS (chr10)
# =============================================================================

cat("========================================\n")
cat("SECTION 8: syt1/nav3 Locus Focus (chr10)\n")
cat("========================================\n\n")

# Define syt1/nav3 region on chr10
# Syt1 is around chr10:108,500,000-108,600,000
# Nav3 is around chr10:109,500,000-109,900,000
# We'll look at a broader region covering both

syt1_nav3_region <- list(
  chr = "chr10",
  start = 107000000,  # ~1.5Mb upstream of syt1
  end = 111000000     # ~1Mb downstream of nav3
)

# Filter boundaries to region
regional_df <- tad_df %>%
  filter(chr == syt1_nav3_region$chr &
           Boundary >= syt1_nav3_region$start &
           Boundary <= syt1_nav3_region$end)

cat(sprintf("Boundaries in syt1/nav3 region (chr10:%d-%d): %d\n",
            syt1_nav3_region$start, syt1_nav3_region$end, nrow(regional_df)))

if (nrow(regional_df) > 0) {

  # Regional statistics
  n_regional_diff <- sum(regional_df$Differential %in% c("Differential", "Shifted") |
                           regional_df$Type == "Shifted", na.rm = TRUE)
  n_regional_shifted <- sum(regional_df$Type == "Shifted", na.rm = TRUE)

  cat(sprintf("  Differential: %d (%.1f%%)\n", n_regional_diff,
              100 * n_regional_diff / nrow(regional_df)))
  cat(sprintf("  Shifted: %d (%.1f%%)\n", n_regional_shifted,
              100 * n_regional_shifted / nrow(regional_df)))

  # 8.1 Regional Overview Plot
  cat("Creating syt1/nav3 regional overview...\n")

  # Create position in Mb for cleaner axis
  regional_df$Position_Mb <- regional_df$Boundary / 1e6

  p_regional <- ggplot(regional_df, aes(x = Position_Mb, y = Gap_Score, color = Type)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = c(108.5, 109.7), linetype = "dotted", color = "gray30", linewidth = 0.8) +
    annotate("text", x = 108.5, y = max(regional_df$Gap_Score) * 0.95,
             label = "Syt1", hjust = 0.5, fontface = "italic", size = 3.5) +
    annotate("text", x = 109.7, y = max(regional_df$Gap_Score) * 0.95,
             label = "Nav3", hjust = 0.5, fontface = "italic", size = 3.5) +
    scale_color_manual(values = type_colors, name = "Boundary Type") +
    labs(
      title = "TAD Boundaries in syt1/nav3 Region",
      subtitle = sprintf("chr10:%d-%d | %d boundaries | %d differential",
                         syt1_nav3_region$start, syt1_nav3_region$end,
                         nrow(regional_df), n_regional_diff),
      x = "Position (Mb)",
      y = "Gap Score"
    ) +
    theme_tad()

  save_multiformat_ggplot(p_regional,
                          file.path(output_base, "syt1_nav3_focus", "syt1_nav3_regional_overview"),
                          width = 12, height = 6)

  # 8.2 Regional Statistics Table
  regional_stats <- regional_df %>%
    select(chr, Boundary, Gap_Score, TAD_Score1, TAD_Score2,
           Differential, Enriched_In, Type, robustness,
           shift_distance_kb, anchor_type)

  write.table(regional_stats,
              file.path(output_base, "syt1_nav3_focus", "syt1_nav3_statistics.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("  Saved: syt1_nav3_statistics.tsv\n")

} else {
  cat("  No boundaries found in specified region.\n")
}

cat("\nSection 8 complete: syt1/nav3 locus focus done\n\n")

# =============================================================================
# SECTION 9: FINAL SUMMARY
# =============================================================================

cat("========================================\n")
cat("SECTION 9: Final Summary\n")
cat("========================================\n\n")

# Generate comprehensive summary
summary_lines <- c(
  "==========================================",
  "TADCompare VISUALIZATION ANALYSIS SUMMARY",
  "==========================================",
  sprintf("Generated: %s", Sys.time()),
  sprintf("Timepoint: %s", timepoint),
  sprintf("Output directory: %s", output_base),
  "",
  "--- DATA OVERVIEW ---",
  sprintf("Total TAD boundaries: %d", n_total),
  sprintf("Differential boundaries: %d (%.1f%%)", n_diff, 100 * n_diff / n_total),
  sprintf("Shifted boundaries: %d (%.1f%%)", n_shifted, 100 * n_shifted / n_total),
  "",
  "--- BOUNDARY TYPE DISTRIBUTION ---"
)

type_table <- table(tad_df$Type)
for (t in names(type_table)) {
  summary_lines <- c(summary_lines,
                     sprintf("  %s: %d (%.1f%%)", t, type_table[t], 100 * type_table[t] / n_total))
}

summary_lines <- c(summary_lines,
                   "",
                   "--- ROBUSTNESS DISTRIBUTION ---")
robust_table <- table(tad_df$robustness)
for (r in names(robust_table)) {
  summary_lines <- c(summary_lines,
                     sprintf("  %s: %d (%.1f%%)", r, robust_table[r], 100 * robust_table[r] / n_total))
}

summary_lines <- c(summary_lines,
                   "",
                   "--- ENRICHMENT DIRECTION (Differential Only) ---")
enrich_table <- table(tad_df$Enriched_In[tad_df$Type != "Non-Differential"])
for (e in names(enrich_table)) {
  label <- ifelse(e == "Matrix 1", "Control-enriched", "Mutant-enriched")
  summary_lines <- c(summary_lines,
                     sprintf("  %s: %d (%.1f%%)", label, enrich_table[e],
                             100 * enrich_table[e] / sum(enrich_table)))
}

if (n_shifted_valid > 0) {
  summary_lines <- c(summary_lines,
                     "",
                     "--- SHIFT DISTANCE STATISTICS ---",
                     sprintf("  Median: %.1f kb", shift_stats$median_kb),
                     sprintf("  Mean: %.1f kb", shift_stats$mean_kb),
                     sprintf("  SD: %.1f kb", shift_stats$sd_kb),
                     sprintf("  Range: %.1f - %.1f kb", shift_stats$min_kb, shift_stats$max_kb))
}

summary_lines <- c(summary_lines,
                   "",
                   "--- ANCHOR TYPE DISTRIBUTION ---")
for (a in names(anchor_summary)) {
  summary_lines <- c(summary_lines,
                     sprintf("  %s: %d (%.1f%%)", a, anchor_summary[a],
                             100 * anchor_summary[a] / n_total))
}

summary_lines <- c(summary_lines,
                   "",
                   "==========================================",
                   "OUTPUT FILES GENERATED",
                   "==========================================",
                   "",
                   "Overview:",
                   "  - gap_score_distribution.pdf",
                   "  - tad_score_scatter.pdf",
                   "  - differential_landscape.pdf",
                   "",
                   "Classification:",
                   "  - boundary_type_pie.pdf",
                   "  - boundary_type_by_chromosome.pdf",
                   "  - enrichment_direction_by_type.pdf",
                   "",
                   "Shift Analysis:",
                   "  - shift_distance_histogram.pdf",
                   "  - shift_distance_violin.pdf",
                   "  - shift_vs_gap_score.pdf",
                   "",
                   "Robustness:",
                   "  - robustness_differential_heatmap.pdf",
                   "  - robustness_by_type.pdf",
                   "",
                   "Chromosome:",
                   "  - per_chromosome_differential.pdf",
                   "",
                   "ChIP-seq:",
                   "  - anchor_classification.pdf",
                   "  - chipseq_overlap_heatmap.pdf",
                   "  - chipseq_by_enrichment_direction.pdf",
                   "",
                   "Enrichment:",
                   "  - go_bp_dotplot.pdf",
                   "  - go_cc_dotplot.pdf",
                   "  - go_mf_dotplot.pdf",
                   "  - kegg_dotplot.pdf",
                   "  - boundary_genes.tsv",
                   "",
                   "syt1/nav3 Focus:",
                   "  - syt1_nav3_regional_overview.pdf",
                   "  - syt1_nav3_statistics.tsv",
                   "",
                   "Summary:",
                   "  - visualization_summary.txt",
                   "",
                   "==========================================")

# Write summary
summary_file <- file.path(output_base, "summary", "visualization_summary.txt")
writeLines(summary_lines, summary_file)
cat(sprintf("Summary saved: %s\n\n", summary_file))

# Print summary to console
cat(paste(summary_lines, collapse = "\n"))
cat("\n\n")

cat("========================================\n")
cat("TADCompare VISUALIZATION ANALYSIS COMPLETE\n")
cat("========================================\n\n")

cat(sprintf("Output directory: %s\n\n", output_base))
cat("All visualizations have been generated successfully!\n\n")

