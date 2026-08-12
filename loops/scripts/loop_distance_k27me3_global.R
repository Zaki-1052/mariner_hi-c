# scripts/loop_distance_k27me3_global.R
# H3K27me3-Filtered Loop Distance Analysis - GLOBAL (All Loops)
# Compares Ctrl-biased vs Mut-biased loops across ALL loops (not just differential)
#
# Purpose:
#   Test whether the distance shift pattern is visible at a global scale by comparing
#   loops with negative logFC (higher in ctrl) vs positive logFC (higher in mut).
#   Filtered by H3K27me3 and bivalent overlaps.
#
# Subsets:
#   1. K27me3-anchored: at least one anchor overlaps H3K27me3
#   2. K27me3-K27me3: BOTH anchors overlap H3K27me3
#   3. Bivalent: at least one anchor overlaps bivalent (K4me3+K27me3) peaks
#
# Usage:
#   Rscript scripts/loop_distance_k27me3_global.R                    # Default: late
#   Rscript scripts/loop_distance_k27me3_global.R --timepoint late   # Late timepoint
#   Rscript scripts/loop_distance_k27me3_global.R --timepoint early  # Early timepoint
#   Rscript scripts/loop_distance_k27me3_global.R --timepoint both   # Run both
#
# Output:
#   output/loops_k27me3_global/{early,late}/

# ==============================================================================
# SECTION 1: SETUP
# ==============================================================================

cat("=== H3K27me3-Filtered Loop Distance Analysis (GLOBAL) ===\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(scales)
  library(GenomicRanges)
})

source("scripts/utils/multi_format_output.R")

# ==============================================================================
# SECTION 2: CONFIGURATION
# ==============================================================================

# Input files - ALL loops (not just differential)
INPUT_FILES <- list(
  late = "25042-late_outputs/merged_loops/merged_all_results.tsv",
  early = "250831-early_outputs/merged_loops/merged_all_results.tsv"
)

# ChIP-seq peak files by timepoint
PEAK_FILES <- list(
  late = list(
    h3k27me3 = "peaks/beds/H3K27me3CerebellumLate1.bed",
    bivalent = "peaks/beds/Bivalent_Cerebellum_Late.bed"
  ),
  early = list(
    h3k27me3 = "peaks/beds/H3K27me3CerebellumEarly1.bed",
    bivalent = "peaks/beds/Bivalent_Cerebellum_Early.bed"
  )
)

# Output directories
OUTPUT_DIRS <- list(
  late = "output/loops_k27me3_global/late",
  early = "output/loops_k27me3_global/early"
)

# Color scheme - using different colors to distinguish from differential analysis
COLORS <- list(
  ctrl = "#2166ac",    # Blue for ctrl-biased (logFC < 0)
  mut = "#b2182b",     # Red for mut-biased (logFC > 0)
  ctrl_light = "#92c5de",
  mut_light = "#f4a582"
)

# Direction labels for global comparison
DIRECTION_LABELS <- c(
  "ctrl_biased" = "Ctrl-biased (logFC < 0)",
  "mut_biased" = "Mut-biased (logFC > 0)"
)

# ==============================================================================
# SECTION 3: ARGUMENT PARSING
# ==============================================================================

parse_arguments <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  timepoint <- "late"

  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--timepoint" && i < length(args)) {
      timepoint <- args[i + 1]
      i <- i + 2
    } else if (args[i] %in% c("--help", "-h")) {
      cat("Usage: Rscript scripts/loop_distance_k27me3_global.R [OPTIONS]\n\n")
      cat("Options:\n")
      cat("  --timepoint TP  Timepoint: 'early', 'late', or 'both' (default: late)\n")
      cat("  --help, -h      Show this help message\n\n")
      cat("Description:\n")
      cat("  Compares loop distance distributions for ALL loops (not just differential)\n")
      cat("  Ctrl-biased: logFC < 0 (higher signal in control)\n")
      cat("  Mut-biased:  logFC > 0 (higher signal in mutant)\n\n")
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

generate_cdf_plot <- function(loops_df, subset_name, output_path, colors) {
  ctrl_loops <- loops_df %>% filter(global_direction == "ctrl_biased")
  mut_loops <- loops_df %>% filter(global_direction == "mut_biased")

  if (nrow(ctrl_loops) < 10 || nrow(mut_loops) < 10) {
    cat(sprintf("    Skipping CDF for %s: insufficient data (Ctrl: %d, Mut: %d)\n",
                subset_name, nrow(ctrl_loops), nrow(mut_loops)))
    return(invisible(NULL))
  }

  median_ctrl <- median(ctrl_loops$loop_distance)
  median_mut <- median(mut_loops$loop_distance)

  ks_test <- ks.test(ctrl_loops$loop_distance, mut_loops$loop_distance)

  p <- ggplot(loops_df, aes(x = loop_distance_kb, color = direction_label)) +
    stat_ecdf(geom = "step", linewidth = 1.2) +
    scale_color_manual(
      values = c("Ctrl-biased (logFC < 0)" = colors$ctrl,
                 "Mut-biased (logFC > 0)" = colors$mut),
      name = "Direction"
    ) +
    scale_x_log10(
      labels = comma,
      breaks = c(10, 100, 1000, 10000),
      limits = c(10, 100000)
    ) +
    geom_vline(xintercept = median_ctrl / 1000, color = colors$ctrl,
               linetype = "dashed", alpha = 0.7) +
    geom_vline(xintercept = median_mut / 1000, color = colors$mut,
               linetype = "dashed", alpha = 0.7) +
    annotate("text", x = 50000, y = 0.15,
             label = sprintf("KS test p = %.2e", ks_test$p.value),
             hjust = 1, size = 4) +
    annotate("text", x = 50000, y = 0.08,
             label = sprintf("Median: Ctrl = %d kb, Mut = %d kb",
                            round(median_ctrl / 1000), round(median_mut / 1000)),
             hjust = 1, size = 3.5) +
    annotate("text", x = 50000, y = 0.22,
             label = sprintf("n = %d loops (%d ctrl, %d mut)",
                            nrow(loops_df), nrow(ctrl_loops), nrow(mut_loops)),
             hjust = 1, size = 3.5, fontface = "italic") +
    labs(
      title = sprintf("Global Loop Distance CDF: %s", subset_name),
      subtitle = "All loops (not just differential) - Ctrl-biased vs Mut-biased",
      x = "Loop Distance (kb, log scale)",
      y = "Cumulative Proportion"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey40"),
      legend.position = c(0.20, 0.85),
      legend.background = element_rect(fill = "white", color = NA),
      panel.grid.minor = element_blank()
    )

  save_multiformat_ggplot(p, output_path, width = 8, height = 6)
  invisible(p)
}

generate_density_plot <- function(loops_df, subset_name, output_path, colors) {
  ctrl_loops <- loops_df %>% filter(global_direction == "ctrl_biased")
  mut_loops <- loops_df %>% filter(global_direction == "mut_biased")

  if (nrow(ctrl_loops) < 10 || nrow(mut_loops) < 10) {
    cat(sprintf("    Skipping density for %s: insufficient data (Ctrl: %d, Mut: %d)\n",
                subset_name, nrow(ctrl_loops), nrow(mut_loops)))
    return(invisible(NULL))
  }

  median_ctrl <- median(ctrl_loops$loop_distance)
  median_mut <- median(mut_loops$loop_distance)

  wilcox_test <- wilcox.test(ctrl_loops$loop_distance, mut_loops$loop_distance)

  p <- ggplot(loops_df, aes(x = loop_distance_kb, fill = direction_label)) +
    geom_density(alpha = 0.5, color = "black", linewidth = 0.5) +
    geom_vline(xintercept = median_ctrl / 1000, color = colors$ctrl,
               linetype = "dashed", linewidth = 1) +
    geom_vline(xintercept = median_mut / 1000, color = colors$mut,
               linetype = "dashed", linewidth = 1) +
    scale_fill_manual(
      values = c("Ctrl-biased (logFC < 0)" = colors$ctrl,
                 "Mut-biased (logFC > 0)" = colors$mut),
      name = ""
    ) +
    scale_x_log10(
      labels = comma,
      breaks = c(10, 100, 1000, 10000),
      limits = c(10, 100000)
    ) +
    annotate("text", x = median_ctrl / 1000 * 1.3, y = Inf,
             label = sprintf("Median\n%d kb", round(median_ctrl / 1000)),
             color = colors$ctrl, vjust = 1.5, hjust = 0, size = 3.5, fontface = "bold") +
    annotate("text", x = median_mut / 1000 * 0.7, y = Inf,
             label = sprintf("Median\n%d kb", round(median_mut / 1000)),
             color = colors$mut, vjust = 1.5, hjust = 1, size = 3.5, fontface = "bold") +
    annotate("text", x = 20, y = 0,
             label = sprintf("Wilcoxon p = %.2e | n = %d", wilcox_test$p.value, nrow(loops_df)),
             hjust = 0, vjust = -0.5, size = 4) +
    labs(
      title = sprintf("Global Loop Distance Density: %s", subset_name),
      subtitle = "All loops - Ctrl-biased (logFC<0) tend to be longer",
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

run_k27me3_global_analysis <- function(timepoint) {
  cat("\n")
  cat("============================================================\n")
  cat(sprintf("Processing %s timepoint (GLOBAL ANALYSIS)\n", toupper(timepoint)))
  cat("============================================================\n\n")

  input_file <- INPUT_FILES[[timepoint]]
  peak_files <- PEAK_FILES[[timepoint]]
  OUTPUT_DIR <- OUTPUT_DIRS[[timepoint]]
  dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

  cat("Input file:", input_file, "\n")
  cat("Output directory:", OUTPUT_DIR, "\n\n")

  if (!file.exists(input_file)) {
    stop(sprintf("ERROR: Input file not found: %s", input_file))
  }

  # ============================================================================
  # Step 1: Load ALL loops
  # ============================================================================

  cat("=== Step 1: Loading ALL Loops ===\n")
  loops <- read_tsv(input_file, show_col_types = FALSE)
  cat("Loaded", nrow(loops), "total loops\n")

  # Calculate loop distance if not present
  if (!"loop_distance" %in% colnames(loops)) {
    loops <- loops %>%
      mutate(
        anchor1_midpoint = (start1 + end1) / 2,
        anchor2_midpoint = (start2 + end2) / 2,
        loop_distance = abs(anchor2_midpoint - anchor1_midpoint)
      )
  }

  # Create global direction based on logFC sign
  loops <- loops %>%
    mutate(
      global_direction = case_when(
        logFC < 0 ~ "ctrl_biased",
        logFC > 0 ~ "mut_biased",
        TRUE ~ "neutral"
      ),
      direction_label = factor(
        case_when(
          global_direction == "ctrl_biased" ~ "Ctrl-biased (logFC < 0)",
          global_direction == "mut_biased" ~ "Mut-biased (logFC > 0)",
          TRUE ~ "Neutral"
        ),
        levels = c("Ctrl-biased (logFC < 0)", "Mut-biased (logFC > 0)")
      ),
      loop_distance_kb = loop_distance / 1000
    )

  # Filter to non-neutral loops
  loops_directional <- loops %>%
    filter(global_direction %in% c("ctrl_biased", "mut_biased"))

  cat("Directional loops:", nrow(loops_directional), "\n")
  cat("  - Ctrl-biased (logFC < 0):", sum(loops_directional$global_direction == "ctrl_biased"), "\n")
  cat("  - Mut-biased (logFC > 0):", sum(loops_directional$global_direction == "mut_biased"), "\n\n")

  # ============================================================================
  # Step 2: Load ChIP-seq peaks
  # ============================================================================

  cat("=== Step 2: Loading ChIP-seq Peaks ===\n")
  k27me3_gr <- load_chip_peaks(peak_files$h3k27me3, "H3K27me3")
  bivalent_gr <- load_chip_peaks(peak_files$bivalent, "Bivalent")
  cat("\n")

  # ============================================================================
  # Step 3: Compute overlaps
  # ============================================================================

  cat("=== Step 3: Computing Overlaps ===\n")

  # Create GRanges for anchors
  anchor1_gr <- GRanges(
    seqnames = loops_directional$chr1,
    ranges = IRanges(start = loops_directional$start1,
                     end = loops_directional$end1)
  )

  anchor2_gr <- GRanges(
    seqnames = loops_directional$chr2,
    ranges = IRanges(start = loops_directional$start2,
                     end = loops_directional$end2)
  )

  loops_directional$anchor1_K27me3 <- countOverlaps(anchor1_gr, k27me3_gr) > 0
  loops_directional$anchor2_K27me3 <- countOverlaps(anchor2_gr, k27me3_gr) > 0
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

  loops_k27me3_anchored <- loops_directional %>%
    filter(anchor1_K27me3 | anchor2_K27me3)

  loops_k27me3_both <- loops_directional %>%
    filter(anchor1_K27me3 & anchor2_K27me3)

  loops_bivalent <- loops_directional %>%
    filter(anchor1_Bivalent | anchor2_Bivalent)

  cat(sprintf("  All loops: %d\n", nrow(loops_directional)))
  cat(sprintf("    - Ctrl-biased: %d, Mut-biased: %d\n",
              sum(loops_directional$global_direction == "ctrl_biased"),
              sum(loops_directional$global_direction == "mut_biased")))

  cat(sprintf("  K27me3-anchored: %d loops\n", nrow(loops_k27me3_anchored)))
  cat(sprintf("    - Ctrl-biased: %d, Mut-biased: %d\n",
              sum(loops_k27me3_anchored$global_direction == "ctrl_biased"),
              sum(loops_k27me3_anchored$global_direction == "mut_biased")))

  cat(sprintf("  K27me3-K27me3 (both): %d loops\n", nrow(loops_k27me3_both)))
  cat(sprintf("    - Ctrl-biased: %d, Mut-biased: %d\n",
              sum(loops_k27me3_both$global_direction == "ctrl_biased"),
              sum(loops_k27me3_both$global_direction == "mut_biased")))

  cat(sprintf("  Bivalent: %d loops\n", nrow(loops_bivalent)))
  cat(sprintf("    - Ctrl-biased: %d, Mut-biased: %d\n\n",
              sum(loops_bivalent$global_direction == "ctrl_biased"),
              sum(loops_bivalent$global_direction == "mut_biased")))

  # ============================================================================
  # Step 5: Generate plots
  # ============================================================================

  cat("=== Step 5: Generating Plots ===\n")

  # All loops
  cat("  Generating All Loops plots...\n")
  generate_cdf_plot(loops_directional, "All Loops",
                    file.path(OUTPUT_DIR, "01_cdf_all_loops"), COLORS)
  generate_density_plot(loops_directional, "All Loops",
                        file.path(OUTPUT_DIR, "03_density_all_loops"), COLORS)

  # K27me3-anchored
  cat("  Generating K27me3-anchored plots...\n")
  generate_cdf_plot(loops_k27me3_anchored, "K27me3-Anchored Loops",
                    file.path(OUTPUT_DIR, "01_cdf_k27me3_anchored"), COLORS)
  generate_density_plot(loops_k27me3_anchored, "K27me3-Anchored Loops",
                        file.path(OUTPUT_DIR, "03_density_k27me3_anchored"), COLORS)

  # K27me3-K27me3
  cat("  Generating K27me3-K27me3 plots...\n")
  generate_cdf_plot(loops_k27me3_both, "K27me3-K27me3 Loops (Both Anchors)",
                    file.path(OUTPUT_DIR, "01_cdf_k27me3_both"), COLORS)
  generate_density_plot(loops_k27me3_both, "K27me3-K27me3 Loops (Both Anchors)",
                        file.path(OUTPUT_DIR, "03_density_k27me3_both"), COLORS)

  # Bivalent
  cat("  Generating Bivalent plots...\n")
  generate_cdf_plot(loops_bivalent, "Bivalent-Anchored Loops",
                    file.path(OUTPUT_DIR, "01_cdf_bivalent"), COLORS)
  generate_density_plot(loops_bivalent, "Bivalent-Anchored Loops",
                        file.path(OUTPUT_DIR, "03_density_bivalent"), COLORS)

  cat("\n")

  # ============================================================================
  # Step 6: Save summary
  # ============================================================================

  cat("=== Step 6: Saving Summary ===\n")

  calc_subset_stats <- function(df, name) {
    if (nrow(df) < 10) {
      return(list(name = name, total = nrow(df),
                  n_ctrl = sum(df$global_direction == "ctrl_biased"),
                  n_mut = sum(df$global_direction == "mut_biased"),
                  median_ctrl = NA, median_mut = NA,
                  ks_pvalue = NA, wilcox_pvalue = NA))
    }

    ctrl <- df %>% filter(global_direction == "ctrl_biased")
    mut <- df %>% filter(global_direction == "mut_biased")

    ks_p <- if (nrow(ctrl) >= 2 && nrow(mut) >= 2) {
      ks.test(ctrl$loop_distance, mut$loop_distance)$p.value
    } else NA

    wilcox_p <- if (nrow(ctrl) >= 2 && nrow(mut) >= 2) {
      wilcox.test(ctrl$loop_distance, mut$loop_distance)$p.value
    } else NA

    list(
      name = name,
      total = nrow(df),
      n_ctrl = nrow(ctrl),
      n_mut = nrow(mut),
      median_ctrl = median(ctrl$loop_distance),
      median_mut = median(mut$loop_distance),
      ks_pvalue = ks_p,
      wilcox_pvalue = wilcox_p
    )
  }

  stats_all <- calc_subset_stats(loops_directional, "All Loops")
  stats_k27me3_anchored <- calc_subset_stats(loops_k27me3_anchored, "K27me3-Anchored")
  stats_k27me3_both <- calc_subset_stats(loops_k27me3_both, "K27me3-K27me3 (Both)")
  stats_bivalent <- calc_subset_stats(loops_bivalent, "Bivalent")

  summary_file <- file.path(OUTPUT_DIR, "filter_summary.txt")
  sink(summary_file)
  cat("=== H3K27me3-Filtered Loop Distance Analysis (GLOBAL) Summary ===\n")
  cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Timepoint:", timepoint, "\n\n")

  cat("GLOBAL ANALYSIS: Comparing ALL loops by logFC sign\n")
  cat("  Ctrl-biased: logFC < 0 (higher signal in control)\n")
  cat("  Mut-biased:  logFC > 0 (higher signal in mutant)\n\n")

  cat("Input file:", input_file, "\n")
  cat("ChIP-seq peaks:\n")
  cat("  H3K27me3:", peak_files$h3k27me3, "\n")
  cat("  Bivalent:", peak_files$bivalent, "\n\n")

  cat("=== Subset Statistics ===\n\n")

  for (stats in list(stats_all, stats_k27me3_anchored, stats_k27me3_both, stats_bivalent)) {
    cat(sprintf("--- %s ---\n", stats$name))
    cat(sprintf("  Total loops: %d\n", stats$total))
    cat(sprintf("  Ctrl-biased: %d (%.1f%%)\n",
                stats$n_ctrl, 100 * stats$n_ctrl / max(stats$total, 1)))
    cat(sprintf("  Mut-biased: %d (%.1f%%)\n",
                stats$n_mut, 100 * stats$n_mut / max(stats$total, 1)))
    if (!is.na(stats$median_ctrl)) {
      cat(sprintf("  Median distance (Ctrl-biased): %.1f kb\n", stats$median_ctrl / 1000))
      cat(sprintf("  Median distance (Mut-biased): %.1f kb\n", stats$median_mut / 1000))
      cat(sprintf("  Median ratio (Ctrl/Mut): %.2f\n", stats$median_ctrl / stats$median_mut))
      cat(sprintf("  KS test p-value: %.2e\n", stats$ks_pvalue))
      cat(sprintf("  Wilcoxon p-value: %.2e\n", stats$wilcox_pvalue))
    }
    cat("\n")
  }

  cat("=== Interpretation ===\n\n")
  cat("If Ctrl-biased loops (logFC < 0) have longer median distances than\n")
  cat("Mut-biased loops (logFC > 0), this suggests the distance shift pattern\n")
  cat("is visible even at a global scale across ALL loops, not just the\n")
  cat("significantly differential ones.\n\n")
  cat("This supports the hypothesis that BAP1-KO globally affects long-range\n")
  cat("chromatin architecture, particularly at H3K27me3-marked regions.\n")

  sink()

  cat("Saved:", summary_file, "\n")

  cat("\n=== Analysis Complete for", toupper(timepoint), "(GLOBAL) ===\n")
  cat("Output directory:", OUTPUT_DIR, "\n")
}

# ==============================================================================
# SECTION 6: MAIN EXECUTION
# ==============================================================================

for (tp in TIMEPOINTS_TO_RUN) {
  tryCatch({
    run_k27me3_global_analysis(tp)
  }, error = function(e) {
    cat(sprintf("\nERROR processing %s timepoint: %s\n", tp, e$message))
  })
}

cat("\n=== All timepoints complete ===\n")
cat("Output directories:\n")
for (tp in TIMEPOINTS_TO_RUN) {
  cat(sprintf("  - %s\n", OUTPUT_DIRS[[tp]]))
}
cat("\nEnd time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
