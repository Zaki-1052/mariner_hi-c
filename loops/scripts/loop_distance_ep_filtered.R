# scripts/loop_distance_ep_filtered.R
# Enhancer-Promoter Loop Distance Analysis
# Generates CDF and density plots for loops where one anchor is enhancer (H3K27ac+)
# and the other is promoter (H3K4me3+)
#
# Purpose:
#   Create distance distribution plots (CDF + density) for E-P loop subsets:
#   1. Lenient E-P: (a1_K27ac & a2_K4me3) OR (a1_K4me3 & a2_K27ac)
#   2. Classified E-P: loop_type_extended contains both "Active_Enhancer" AND "Active_Promoter"
#
# Biological rationale:
#   - Lenient: Captures any enhancer-promoter combination based on histone marks
#   - Classified: Uses 8-category anchor classification (Active_Enhancer requires K27ac+
#     AND >2kb from TSS; Active_Promoter requires K4me3+ AND NOT K27me3 AND <=2kb from TSS)
#
# Usage:
#   Rscript scripts/loop_distance_ep_filtered.R                    # Default: late
#   Rscript scripts/loop_distance_ep_filtered.R --timepoint late   # Late timepoint
#   Rscript scripts/loop_distance_ep_filtered.R --timepoint early  # Early timepoint
#   Rscript scripts/loop_distance_ep_filtered.R --timepoint both   # Run both
#
# Output:
#   output/loops_ep_filtered/{early,late}/
#   - 01_cdf_ep_lenient.{pdf,svg,jpg}
#   - 01_cdf_ep_classified.{pdf,svg,jpg}
#   - 03_density_ep_lenient.{pdf,svg,jpg}
#   - 03_density_ep_classified.{pdf,svg,jpg}
#   - filter_summary.txt

# ==============================================================================
# SECTION 1: SETUP
# ==============================================================================

cat("=== Enhancer-Promoter Loop Distance Analysis ===\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(scales)
})

source("scripts/utils/multi_format_output.R")

# ==============================================================================
# SECTION 2: CONFIGURATION
# ==============================================================================

# Input files by timepoint (extended_characterized_loops.tsv from annotate_loops_extended.R)
INPUT_FILES <- list(
  late = "peaks/loop_annotation_extended/late/extended_characterized_loops.tsv",
  early = "peaks/loop_annotation_extended/early/extended_characterized_loops.tsv"
)

# Output directories by timepoint
OUTPUT_DIRS <- list(
  late = "output/loops_ep_filtered/late",
  early = "output/loops_ep_filtered/early"
)

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

# Minimum samples required for statistical tests
MIN_SAMPLES <- 5

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
      cat("Usage: Rscript scripts/loop_distance_ep_filtered.R [OPTIONS]\n\n")
      cat("Options:\n")
      cat("  --timepoint TP  Timepoint: 'early', 'late', or 'both' (default: late)\n")
      cat("  --help, -h      Show this help message\n\n")
      cat("Output:\n")
      cat("  output/loops_ep_filtered/{early,late}/\n")
      cat("    - 01_cdf_ep_lenient.{pdf,svg,jpg}\n")
      cat("    - 01_cdf_ep_classified.{pdf,svg,jpg}\n")
      cat("    - 03_density_ep_lenient.{pdf,svg,jpg}\n")
      cat("    - 03_density_ep_classified.{pdf,svg,jpg}\n")
      cat("    - filter_summary.txt\n\n")
      cat("Filter modes:\n")
      cat("  - Lenient: (K27ac_a1 & K4me3_a2) OR (K4me3_a1 & K27ac_a2)\n")
      cat("  - Classified: loop_type_extended contains 'Active_Enhancer' AND 'Active_Promoter'\n\n")
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

#' Filter loops using lenient E-P criteria
#'
#' Lenient filter: either (anchor1 is K27ac+ AND anchor2 is K4me3+)
#' OR (anchor1 is K4me3+ AND anchor2 is K27ac+)
#'
#' @param df Data frame with extended characterized loops
#' @return Filtered data frame
filter_ep_lenient <- function(df) {
  df %>% filter(
    (anchor1_H3K27ac_overlap & anchor2_H3K4me3_overlap) |
    (anchor1_H3K4me3_overlap & anchor2_H3K27ac_overlap)
  )
}

#' Filter loops using classification-based E-P criteria
#'
#' Classified filter: loop_type_extended contains both "Active_Enhancer" AND "Active_Promoter"
#' This is more stringent because Active_Enhancer requires K27ac+ AND >2kb from TSS,
#' and Active_Promoter requires K4me3+ AND NOT K27me3 AND <=2kb from TSS
#'
#' @param df Data frame with extended characterized loops
#' @return Filtered data frame
filter_ep_classified <- function(df) {
  df %>% filter(
    grepl("Active_Enhancer", loop_type_extended) &
    grepl("Active_Promoter", loop_type_extended)
  )
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

  # Skip if insufficient data
  if (nrow(lost_loops) < MIN_SAMPLES || nrow(gained_loops) < MIN_SAMPLES) {
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
      subtitle = "Enhancer-Promoter loop distance distributions",
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

  # Skip if insufficient data
  if (nrow(lost_loops) < MIN_SAMPLES || nrow(gained_loops) < MIN_SAMPLES) {
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
      subtitle = "Enhancer-Promoter loop distance distributions",
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

#' Calculate statistics for a filtered subset
#'
#' @param df Data frame with filtered loops
#' @param name Name of the subset for reporting
#' @return List with statistics
calc_subset_stats <- function(df, name) {
  if (nrow(df) < 10) {
    return(list(
      name = name,
      total = nrow(df),
      n_lost = sum(df$direction == "down_in_mutant"),
      n_gained = sum(df$direction == "up_in_mutant"),
      median_lost = NA,
      median_gained = NA,
      median_ratio = NA,
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

  median_lost <- if (nrow(lost) > 0) median(lost$loop_distance) else NA
  median_gained <- if (nrow(gained) > 0) median(gained$loop_distance) else NA
  median_ratio <- if (!is.na(median_lost) && !is.na(median_gained) && median_gained > 0) {
    median_lost / median_gained
  } else NA

  list(
    name = name,
    total = nrow(df),
    n_lost = nrow(lost),
    n_gained = nrow(gained),
    median_lost = median_lost,
    median_gained = median_gained,
    median_ratio = median_ratio,
    ks_pvalue = ks_p,
    wilcox_pvalue = wilcox_p
  )
}

# ==============================================================================
# SECTION 5: MAIN ANALYSIS FUNCTION
# ==============================================================================

run_ep_filtered_analysis <- function(timepoint) {
  cat("\n")
  cat("============================================================\n")
  cat(sprintf("Processing %s timepoint\n", toupper(timepoint)))
  cat("============================================================\n\n")

  # Set paths
  input_file <- INPUT_FILES[[timepoint]]
  OUTPUT_DIR <- OUTPUT_DIRS[[timepoint]]
  dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

  cat("Input file:", input_file, "\n")
  cat("Output directory:", OUTPUT_DIR, "\n\n")

  # Validate input
  if (!file.exists(input_file)) {
    stop(sprintf("ERROR: Input file not found: %s\nRun annotate_loops_extended.R first.", input_file))
  }

  # ============================================================================
  # Step 1: Load extended characterized loops
  # ============================================================================

  cat("=== Step 1: Loading Data ===\n")
  loops <- read_tsv(input_file, show_col_types = FALSE)
  cat("Loaded", nrow(loops), "differential loops\n")

  # Validate required columns
  required_cols <- c("anchor1_H3K27ac_overlap", "anchor2_H3K27ac_overlap",
                     "anchor1_H3K4me3_overlap", "anchor2_H3K4me3_overlap",
                     "loop_type_extended", "direction", "loop_distance")

  missing_cols <- setdiff(required_cols, names(loops))
  if (length(missing_cols) > 0) {
    stop(sprintf("ERROR: Missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }

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
  # Step 2: Report ChIP-seq mark overlaps
  # ============================================================================

  cat("=== Step 2: ChIP-seq Mark Overlaps ===\n")

  cat("  H3K27ac overlaps (enhancer mark):\n")
  cat(sprintf("    Anchor1: %d loops (%.1f%%)\n",
              sum(loops_directional$anchor1_H3K27ac_overlap),
              100 * mean(loops_directional$anchor1_H3K27ac_overlap)))
  cat(sprintf("    Anchor2: %d loops (%.1f%%)\n",
              sum(loops_directional$anchor2_H3K27ac_overlap),
              100 * mean(loops_directional$anchor2_H3K27ac_overlap)))

  cat("  H3K4me3 overlaps (promoter mark):\n")
  cat(sprintf("    Anchor1: %d loops (%.1f%%)\n",
              sum(loops_directional$anchor1_H3K4me3_overlap),
              100 * mean(loops_directional$anchor1_H3K4me3_overlap)))
  cat(sprintf("    Anchor2: %d loops (%.1f%%)\n\n",
              sum(loops_directional$anchor2_H3K4me3_overlap),
              100 * mean(loops_directional$anchor2_H3K4me3_overlap)))

  # ============================================================================
  # Step 3: Create filtered subsets
  # ============================================================================

  cat("=== Step 3: Creating E-P Filtered Subsets ===\n")

  # Subset 1: Lenient E-P filter
  loops_ep_lenient <- filter_ep_lenient(loops_directional)

  # Subset 2: Classified E-P filter
  loops_ep_classified <- filter_ep_classified(loops_directional)

  cat(sprintf("  Lenient E-P (mark-based): %d loops\n", nrow(loops_ep_lenient)))
  cat(sprintf("    Filter: (a1_K27ac & a2_K4me3) OR (a1_K4me3 & a2_K27ac)\n"))
  cat(sprintf("    - Lost: %d, Gained: %d\n",
              sum(loops_ep_lenient$direction == "down_in_mutant"),
              sum(loops_ep_lenient$direction == "up_in_mutant")))

  cat(sprintf("  Classified E-P (Active_Enhancer + Active_Promoter): %d loops\n",
              nrow(loops_ep_classified)))
  cat(sprintf("    Filter: loop_type_extended contains 'Active_Enhancer' AND 'Active_Promoter'\n"))
  cat(sprintf("    - Lost: %d, Gained: %d\n\n",
              sum(loops_ep_classified$direction == "down_in_mutant"),
              sum(loops_ep_classified$direction == "up_in_mutant")))

  # ============================================================================
  # Step 4: Generate plots for each subset
  # ============================================================================

  cat("=== Step 4: Generating Plots ===\n")

  # Lenient E-P plots
  cat("  Generating Lenient E-P plots...\n")
  generate_cdf_plot(
    loops_ep_lenient,
    "Enhancer-Promoter Loops (Lenient)",
    file.path(OUTPUT_DIR, "01_cdf_ep_lenient"),
    COLORS
  )
  generate_density_plot(
    loops_ep_lenient,
    "Enhancer-Promoter Loops (Lenient)",
    file.path(OUTPUT_DIR, "03_density_ep_lenient"),
    COLORS
  )

  # Classified E-P plots
  cat("  Generating Classified E-P plots...\n")
  generate_cdf_plot(
    loops_ep_classified,
    "Enhancer-Promoter Loops (Classified)",
    file.path(OUTPUT_DIR, "01_cdf_ep_classified"),
    COLORS
  )
  generate_density_plot(
    loops_ep_classified,
    "Enhancer-Promoter Loops (Classified)",
    file.path(OUTPUT_DIR, "03_density_ep_classified"),
    COLORS
  )

  cat("\n")

  # ============================================================================
  # Step 5: Save summary statistics
  # ============================================================================

  cat("=== Step 5: Saving Summary ===\n")

  # Calculate statistics for each subset
  stats_all <- calc_subset_stats(loops_directional, "All Loops")
  stats_ep_lenient <- calc_subset_stats(loops_ep_lenient, "Lenient E-P")
  stats_ep_classified <- calc_subset_stats(loops_ep_classified, "Classified E-P")

  # Write summary file
  summary_file <- file.path(OUTPUT_DIR, "filter_summary.txt")
  sink(summary_file)
  cat("=== Enhancer-Promoter Loop Distance Analysis Summary ===\n")
  cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Timepoint:", timepoint, "\n\n")

  cat("Input file:", input_file, "\n\n")

  cat("=== Filter Definitions ===\n\n")
  cat("Lenient E-P filter:\n")
  cat("  (anchor1_H3K27ac+ AND anchor2_H3K4me3+) OR\n")
  cat("  (anchor1_H3K4me3+ AND anchor2_H3K27ac+)\n\n")

  cat("Classified E-P filter:\n")
  cat("  loop_type_extended contains 'Active_Enhancer' AND 'Active_Promoter'\n")
  cat("  (Active_Enhancer = K27ac+ AND >2kb from TSS)\n")
  cat("  (Active_Promoter = K4me3+ AND NOT K27me3 AND <=2kb from TSS)\n\n")

  cat("=== Subset Statistics ===\n\n")

  for (stats in list(stats_all, stats_ep_lenient, stats_ep_classified)) {
    cat(sprintf("--- %s ---\n", stats$name))
    cat(sprintf("  Total loops: %d\n", stats$total))
    cat(sprintf("  Lost in BAP1-KO: %d (%.1f%%)\n",
                stats$n_lost, 100 * stats$n_lost / max(stats$total, 1)))
    cat(sprintf("  Gained in BAP1-KO: %d (%.1f%%)\n",
                stats$n_gained, 100 * stats$n_gained / max(stats$total, 1)))
    if (!is.na(stats$median_lost) && !is.na(stats$median_gained)) {
      cat(sprintf("  Median distance (Lost): %.1f kb\n", stats$median_lost / 1000))
      cat(sprintf("  Median distance (Gained): %.1f kb\n", stats$median_gained / 1000))
      if (!is.na(stats$median_ratio)) {
        cat(sprintf("  Median ratio (Lost/Gained): %.2f\n", stats$median_ratio))
      }
      if (!is.na(stats$ks_pvalue)) {
        cat(sprintf("  KS test p-value: %.2e\n", stats$ks_pvalue))
      }
      if (!is.na(stats$wilcox_pvalue)) {
        cat(sprintf("  Wilcoxon p-value: %.2e\n", stats$wilcox_pvalue))
      }
    } else {
      cat("  (Insufficient data for statistical comparison)\n")
    }
    cat("\n")
  }

  cat("=== ChIP-seq Mark Overlap Breakdown ===\n\n")
  cat(sprintf("Anchor1 H3K27ac+: %d / %d (%.1f%%)\n",
              sum(loops_directional$anchor1_H3K27ac_overlap),
              nrow(loops_directional),
              100 * mean(loops_directional$anchor1_H3K27ac_overlap)))
  cat(sprintf("Anchor2 H3K27ac+: %d / %d (%.1f%%)\n",
              sum(loops_directional$anchor2_H3K27ac_overlap),
              nrow(loops_directional),
              100 * mean(loops_directional$anchor2_H3K27ac_overlap)))
  cat(sprintf("Anchor1 H3K4me3+: %d / %d (%.1f%%)\n",
              sum(loops_directional$anchor1_H3K4me3_overlap),
              nrow(loops_directional),
              100 * mean(loops_directional$anchor1_H3K4me3_overlap)))
  cat(sprintf("Anchor2 H3K4me3+: %d / %d (%.1f%%)\n\n",
              sum(loops_directional$anchor2_H3K4me3_overlap),
              nrow(loops_directional),
              100 * mean(loops_directional$anchor2_H3K4me3_overlap)))

  cat("=== Loop Type Breakdown (Classified E-P) ===\n\n")
  if (nrow(loops_ep_classified) > 0) {
    type_counts <- loops_ep_classified %>%
      count(loop_type_extended, direction) %>%
      pivot_wider(names_from = direction, values_from = n, values_fill = 0) %>%
      arrange(desc(down_in_mutant + up_in_mutant))

    print(as.data.frame(type_counts))
  } else {
    cat("No loops matched the classified E-P filter\n")
  }

  sink()

  cat("Saved:", summary_file, "\n")

  # ============================================================================
  # Completion
  # ============================================================================

  cat("\n=== Analysis Complete for", toupper(timepoint), "===\n")
  cat("Output directory:", OUTPUT_DIR, "\n")
  cat("Files generated:\n")
  cat("  - 01_cdf_ep_lenient/{pdf,svg,jpg}\n")
  cat("  - 01_cdf_ep_classified/{pdf,svg,jpg}\n")
  cat("  - 03_density_ep_lenient/{pdf,svg,jpg}\n")
  cat("  - 03_density_ep_classified/{pdf,svg,jpg}\n")
  cat("  - filter_summary.txt\n")

  # Return summary for programmatic access
  invisible(list(
    timepoint = timepoint,
    n_total = nrow(loops_directional),
    n_ep_lenient = nrow(loops_ep_lenient),
    n_ep_classified = nrow(loops_ep_classified),
    stats_lenient = stats_ep_lenient,
    stats_classified = stats_ep_classified
  ))
}

# ==============================================================================
# SECTION 6: MAIN EXECUTION
# ==============================================================================

results <- list()

for (tp in TIMEPOINTS_TO_RUN) {
  tryCatch({
    results[[tp]] <- run_ep_filtered_analysis(tp)
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
