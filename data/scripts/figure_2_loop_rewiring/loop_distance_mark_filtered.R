# scripts/loop_distance_mark_filtered.R
# Generalized ChIP-seq Mark-Filtered Loop Distance Analysis
# Generates CDF and density plots for loops filtered by anchor ChIP-seq overlaps
#
# Purpose:
#   Create distance distribution plots (CDF + density) for mark-specific subsets:
#   - One-anchor: at least one anchor overlaps the mark
#   - Both-anchors: both anchors overlap the mark (e.g., super-enhancer for H3K27ac)
#
# Supported marks: H3K27ac, H3K27me3, H3K4me1, H3K4me3, Bivalent, CTCF
#
# Usage:
#   Rscript scripts/loop_distance_mark_filtered.R                          # All marks, late
#   Rscript scripts/loop_distance_mark_filtered.R --timepoint both         # All marks, both
#   Rscript scripts/loop_distance_mark_filtered.R --marks H3K27ac          # Single mark
#   Rscript scripts/loop_distance_mark_filtered.R --marks H3K27ac,H3K4me3  # Multiple marks
#
# Output:
#   output/loops_mark_filtered/{early,late}/{mark}/
#   - 01_cdf_one_anchor/{pdf,svg,jpg}
#   - 01_cdf_both_anchors/{pdf,svg,jpg}
#   - 03_density_one_anchor/{pdf,svg,jpg}
#   - 03_density_both_anchors/{pdf,svg,jpg}
#   - filter_summary.txt

# ==============================================================================
# SECTION 1: SETUP
# ==============================================================================

cat("=== Generalized Mark-Filtered Loop Distance Analysis ===\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(scales)
})

source("data/scripts/_shared/multi_format_output.R")  # Original: source("scripts/utils/multi_format_output.R")

# ==============================================================================
# SECTION 2: CONFIGURATION
# ==============================================================================

# Mark configuration: column names and display info for each mark
# Uses pre-computed overlap columns from extended_characterized_loops.tsv
MARK_CONFIG <- list(
  H3K27ac = list(
    col1 = "anchor1_H3K27ac_overlap",
    col2 = "anchor2_H3K27ac_overlap",
    display_name = "H3K27ac",
    biological_role = "Active Enhancer/Super-Enhancer",
    dir_name = "h3k27ac"
  ),
  H3K27me3 = list(
    col1 = "anchor1_H3K27me3_overlap",
    col2 = "anchor2_H3K27me3_overlap",
    display_name = "H3K27me3",
    biological_role = "Polycomb Repression",
    dir_name = "h3k27me3"
  ),
  H3K4me1 = list(
    col1 = "anchor1_H3K4me1_overlap",
    col2 = "anchor2_H3K4me1_overlap",
    display_name = "H3K4me1",
    biological_role = "Poised Enhancer",
    dir_name = "h3k4me1"
  ),
  H3K4me3 = list(
    col1 = "anchor1_H3K4me3_overlap",
    col2 = "anchor2_H3K4me3_overlap",
    display_name = "H3K4me3",
    biological_role = "Active Promoter",
    dir_name = "h3k4me3"
  ),
  Bivalent = list(
    col1 = "anchor1_Bivalent_Promoter_overlap",
    col2 = "anchor2_Bivalent_Promoter_overlap",
    display_name = "Bivalent",
    biological_role = "Bivalent Promoter (K4me3+K27me3)",
    dir_name = "bivalent"
  ),
  CTCF = list(
    col1 = "anchor1_CTCF_overlap",
    col2 = "anchor2_CTCF_overlap",
    display_name = "CTCF",
    biological_role = "CTCF/Cohesin Anchor (Loop Extrusion)",
    dir_name = "ctcf"
  )
)

# Input files by timepoint (extended_characterized_loops.tsv with pre-computed overlaps)
INPUT_FILES <- list(
  late = "data/upstream/loop_calls/late_characterized_loops.tsv",  # Original: "peaks/loop_annotation_extended/late/extended_characterized_loops.tsv"
  early = "peaks/loop_annotation_extended/early/extended_characterized_loops.tsv"  # Note: repo-relative path, not bundled in data/
)

# Output directories
# Original: OUTPUT_BASE <- "output/loops_mark_filtered"
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

# Minimum sample size for reliable statistics
MIN_SAMPLES <- 10

# ==============================================================================
# SECTION 3: ARGUMENT PARSING
# ==============================================================================

parse_arguments <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  timepoint <- "late"  # Default
  marks <- "all"       # Default: run all marks

  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--timepoint" && i < length(args)) {
      timepoint <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--marks" && i < length(args)) {
      marks <- args[i + 1]
      i <- i + 2
    } else if (args[i] %in% c("--help", "-h")) {
      cat("Usage: Rscript scripts/loop_distance_mark_filtered.R [OPTIONS]\n\n")
      cat("Options:\n")
      cat("  --timepoint TP   Timepoint: 'early', 'late', or 'both' (default: late)\n")
      cat("  --marks MARKS    Marks to analyze: 'all' or comma-separated list\n")
      cat("                   Available: H3K27ac,H3K27me3,H3K4me1,H3K4me3,Bivalent,CTCF\n")
      cat("                   (default: all)\n")
      cat("  --help, -h       Show this help message\n\n")
      cat("Output:\n")
      cat("  output/loops_mark_filtered/{early,late}/{mark}/\n")
      cat("    - 01_cdf_one_anchor/{pdf,svg,jpg}\n")
      cat("    - 01_cdf_both_anchors/{pdf,svg,jpg}\n")
      cat("    - 03_density_one_anchor/{pdf,svg,jpg}\n")
      cat("    - 03_density_both_anchors/{pdf,svg,jpg}\n")
      cat("    - filter_summary.txt\n\n")
      cat("Examples:\n")
      cat("  Rscript scripts/loop_distance_mark_filtered.R\n")
      cat("  Rscript scripts/loop_distance_mark_filtered.R --timepoint both\n")
      cat("  Rscript scripts/loop_distance_mark_filtered.R --marks H3K27ac --timepoint late\n")
      cat("  Rscript scripts/loop_distance_mark_filtered.R --marks H3K27ac,H3K4me3 --timepoint both\n")
      quit(save = "no", status = 0)
    } else {
      i <- i + 1
    }
  }

  if (!timepoint %in% c("early", "late", "both")) {
    stop("ERROR: timepoint must be 'early', 'late', or 'both'")
  }

  return(list(timepoint = timepoint, marks = marks))
}

args_parsed <- parse_arguments()

# Resolve timepoints to run
if (args_parsed$timepoint == "both") {
  TIMEPOINTS_TO_RUN <- c("late", "early")
} else {
  TIMEPOINTS_TO_RUN <- args_parsed$timepoint
}

# Resolve marks to run
if (args_parsed$marks == "all") {
  MARKS_TO_RUN <- names(MARK_CONFIG)
} else {
  MARKS_TO_RUN <- strsplit(args_parsed$marks, ",")[[1]]
  invalid_marks <- setdiff(MARKS_TO_RUN, names(MARK_CONFIG))
  if (length(invalid_marks) > 0) {
    stop(sprintf("ERROR: Invalid marks: %s\nAvailable: %s",
                 paste(invalid_marks, collapse = ", "),
                 paste(names(MARK_CONFIG), collapse = ", ")))
  }
}

cat("Timepoint(s) to process:", paste(TIMEPOINTS_TO_RUN, collapse = ", "), "\n")
cat("Mark(s) to analyze:", paste(MARKS_TO_RUN, collapse = ", "), "\n\n")

# ==============================================================================
# SECTION 4: HELPER FUNCTIONS
# ==============================================================================

#' Apply mark filter to loops data
#'
#' @param df Data frame with overlap columns
#' @param col1 Column name for anchor1 overlap
#' @param col2 Column name for anchor2 overlap
#' @param mode Filter mode: "one" (at least one anchor) or "both" (both anchors)
#' @return Filtered data frame
create_mark_filter <- function(df, col1, col2, mode) {
  if (mode == "one") {
    df %>% filter(.data[[col1]] | .data[[col2]])
  } else if (mode == "both") {
    df %>% filter(.data[[col1]] & .data[[col2]])
  } else {
    stop("Invalid mode: must be 'one' or 'both'")
  }
}

#' Calculate statistics for a filtered subset
#'
#' @param df Data frame with filtered loops
#' @param name Name of the subset
#' @return List with statistics
calc_subset_stats <- function(df, name) {
  if (nrow(df) < MIN_SAMPLES) {
    return(list(
      name = name,
      total = nrow(df),
      n_lost = sum(df$direction == "down_in_mutant"),
      n_gained = sum(df$direction == "up_in_mutant"),
      median_lost = NA,
      median_gained = NA,
      ks_pvalue = NA,
      wilcox_pvalue = NA,
      sufficient_data = FALSE
    ))
  }

  lost <- df %>% filter(direction == "down_in_mutant")
  gained <- df %>% filter(direction == "up_in_mutant")

  ks_p <- if (nrow(lost) >= 2 && nrow(gained) >= 2) {
    tryCatch(ks.test(lost$loop_distance, gained$loop_distance)$p.value, error = function(e) NA)
  } else NA

  wilcox_p <- if (nrow(lost) >= 2 && nrow(gained) >= 2) {
    tryCatch(wilcox.test(lost$loop_distance, gained$loop_distance)$p.value, error = function(e) NA)
  } else NA

  list(
    name = name,
    total = nrow(df),
    n_lost = nrow(lost),
    n_gained = nrow(gained),
    median_lost = if (nrow(lost) > 0) median(lost$loop_distance) else NA,
    median_gained = if (nrow(gained) > 0) median(gained$loop_distance) else NA,
    ks_pvalue = ks_p,
    wilcox_pvalue = wilcox_p,
    sufficient_data = TRUE
  )
}

#' Generate CDF plot for a filtered subset
#'
#' @param loops_df Data frame with filtered loops
#' @param subset_name Name of the subset for title
#' @param output_path Output path (without extension)
#' @param colors Color list
#' @param mark_config Mark configuration for subtitle
#' @return ggplot object (invisibly)
generate_cdf_plot <- function(loops_df, subset_name, output_path, colors, mark_config = NULL) {
  # Separate by direction
  lost_loops <- loops_df %>% filter(direction == "down_in_mutant")
  gained_loops <- loops_df %>% filter(direction == "up_in_mutant")

  # Skip if insufficient data
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

  # Determine subtitle based on median comparison
  if (median_lost > median_gained) {
    subtitle_text <- "BAP1-KO preferentially loses long-range loops"
  } else if (median_lost < median_gained) {
    subtitle_text <- "BAP1-KO preferentially loses short-range loops"
  } else {
    subtitle_text <- "Similar distance distributions"
  }

  if (!is.null(mark_config)) {
    subtitle_text <- paste0(subtitle_text, " | ", mark_config$biological_role)
  }

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
      subtitle = subtitle_text,
      x = "Loop Distance (kb, log scale)",
      y = "Cumulative Proportion"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"),
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
#' @param mark_config Mark configuration for subtitle
#' @return ggplot object (invisibly)
generate_density_plot <- function(loops_df, subset_name, output_path, colors, mark_config = NULL) {
  # Separate by direction
  lost_loops <- loops_df %>% filter(direction == "down_in_mutant")
  gained_loops <- loops_df %>% filter(direction == "up_in_mutant")

  # Skip if insufficient data
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

  # Determine subtitle based on median comparison
  if (median_lost > median_gained) {
    subtitle_text <- "Lost loops are systematically longer than gained loops"
  } else if (median_lost < median_gained) {
    subtitle_text <- "Lost loops are systematically shorter than gained loops"
  } else {
    subtitle_text <- "Similar distance distributions"
  }

  if (!is.null(mark_config)) {
    subtitle_text <- paste0(subtitle_text, " | ", mark_config$biological_role)
  }

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
      subtitle = subtitle_text,
      x = "Loop Distance (kb, log scale)",
      y = "Density"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"),
      legend.position = "top",
      panel.grid.minor = element_blank()
    )

  save_multiformat_ggplot(p, output_path, width = 8, height = 6)
  invisible(p)
}

# ==============================================================================
# SECTION 5: MAIN ANALYSIS FUNCTION
# ==============================================================================

#' Run analysis for a single mark and timepoint
#'
#' @param timepoint Character: "early" or "late"
#' @param mark_name Character: name of mark (e.g., "H3K27ac")
#' @param mark_config List: configuration for this mark
#' @param loops_directional Data frame: pre-loaded directional loops
#' @return List with statistics for both filter modes
run_mark_analysis <- function(timepoint, mark_name, mark_config, loops_directional) {
  cat(sprintf("\n--- %s (%s) ---\n", mark_config$display_name, mark_config$biological_role))

  # Original: mark_dir <- file.path(output_base, mark_config$dir_name); dir.create(mark_dir, ...)
  dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)
  dir.create(TSV_DIR, showWarnings = FALSE, recursive = TRUE)

  # Check that required columns exist
  if (!mark_config$col1 %in% names(loops_directional) ||
      !mark_config$col2 %in% names(loops_directional)) {
    cat(sprintf("  WARNING: Columns %s or %s not found. Skipping %s.\n",
                mark_config$col1, mark_config$col2, mark_name))
    return(NULL)
  }

  # ============================================================================
  # Filter 1: One-anchor (at least one anchor has mark)
  # ============================================================================

  cat("  Filter 1: One+ anchor overlaps mark\n")
  loops_one <- create_mark_filter(loops_directional, mark_config$col1, mark_config$col2, "one")
  stats_one <- calc_subset_stats(loops_one, paste0(mark_config$display_name, " (One+ Anchor)"))

  cat(sprintf("    Total: %d loops (Lost: %d, Gained: %d)\n",
              stats_one$total, stats_one$n_lost, stats_one$n_gained))

  if (stats_one$sufficient_data) {
    subset_name_one <- sprintf("%s-Anchored Loops (One+)", mark_config$display_name)

    generate_cdf_plot(
      loops_one,
      subset_name_one,
      file.path(PLOT_DIR, sprintf("2E_%s_cdf_one_anchor", mark_config$dir_name)),  # Original: file.path(mark_dir, "01_cdf_one_anchor")
      COLORS,
      mark_config
    )
    generate_density_plot(
      loops_one,
      subset_name_one,
      file.path(PLOT_DIR, sprintf("2I_%s_density_one_anchor", mark_config$dir_name)),  # Original: file.path(mark_dir, "03_density_one_anchor")
      COLORS,
      mark_config
    )
  } else {
    cat("    Insufficient data for one-anchor plots\n")
  }

  # ============================================================================
  # Filter 2: Both-anchors (both anchors have mark)
  # ============================================================================

  cat("  Filter 2: Both anchors overlap mark\n")
  loops_both <- create_mark_filter(loops_directional, mark_config$col1, mark_config$col2, "both")
  stats_both <- calc_subset_stats(loops_both, paste0(mark_config$display_name, " (Both Anchors)"))

  cat(sprintf("    Total: %d loops (Lost: %d, Gained: %d)\n",
              stats_both$total, stats_both$n_lost, stats_both$n_gained))

  if (stats_both$sufficient_data) {
    # Determine biological label for both-anchor
    both_label <- switch(mark_name,
                         "H3K27ac" = "Super-Enhancer",
                         "H3K27me3" = "Polycomb-Polycomb",
                         "H3K4me1" = "Enhancer-Enhancer",
                         "H3K4me3" = "Promoter-Promoter",
                         "Bivalent" = "Bivalent-Bivalent",
                         "CTCF" = "CTCF-CTCF",
                         "Both Anchors")

    subset_name_both <- sprintf("%s Loops (Both Anchors = %s)",
                                mark_config$display_name, both_label)

    generate_cdf_plot(
      loops_both,
      subset_name_both,
      file.path(PLOT_DIR, sprintf("2E_%s_cdf_both_anchors", mark_config$dir_name)),  # Original: file.path(mark_dir, "01_cdf_both_anchors")
      COLORS,
      mark_config
    )
    generate_density_plot(
      loops_both,
      subset_name_both,
      file.path(PLOT_DIR, sprintf("2I_%s_density_both_anchors", mark_config$dir_name)),  # Original: file.path(mark_dir, "03_density_both_anchors")
      COLORS,
      mark_config
    )
  } else {
    cat(sprintf("    Insufficient data for both-anchor plots (n=%d < %d minimum)\n",
                stats_both$total, MIN_SAMPLES))
  }

  # Return statistics for summary
  list(one_anchor = stats_one, both_anchors = stats_both)
}

#' Run analysis for a single timepoint across all marks
#'
#' @param timepoint Character: "early" or "late"
#' @return NULL (invisibly)
run_timepoint_analysis <- function(timepoint) {
  cat("\n")
  cat("============================================================\n")
  cat(sprintf("Processing %s timepoint\n", toupper(timepoint)))
  cat("============================================================\n\n")

  # Set paths
  input_file <- INPUT_FILES[[timepoint]]
  dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)  # Original: output_base <- file.path(OUTPUT_BASE, timepoint); dir.create(output_base, ...)
  dir.create(TSV_DIR, showWarnings = FALSE, recursive = TRUE)

  cat("Input file:", input_file, "\n")
  cat("Plot directory:", PLOT_DIR, "\n")
  cat("TSV directory:", TSV_DIR, "\n\n")

  # Validate input
  if (!file.exists(input_file)) {
    stop(sprintf("ERROR: Input file not found: %s\nRun peaks/annotate_loops_extended.R first.",
                 input_file))
  }

  # ==========================================================================
  # Step 1: Load data
  # ==========================================================================

  cat("=== Step 1: Loading Data ===\n")
  loops <- read_tsv(input_file, show_col_types = FALSE)
  cat("Loaded", nrow(loops), "annotated loops\n")

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

  # ==========================================================================
  # Step 2: Run analysis for each mark
  # ==========================================================================

  cat("=== Step 2: Analyzing Marks ===\n")

  all_stats <- list()
  for (mark_name in MARKS_TO_RUN) {
    mark_config <- MARK_CONFIG[[mark_name]]
    tryCatch({
      stats <- run_mark_analysis(timepoint, mark_name, mark_config, loops_directional)
      if (!is.null(stats)) {
        all_stats[[mark_name]] <- stats
      }
    }, error = function(e) {
      cat(sprintf("  ERROR processing %s: %s\n", mark_name, e$message))
    })
  }

  # ==========================================================================
  # Step 3: Write summary file
  # ==========================================================================

  cat("\n=== Step 3: Writing Summary ===\n")

  summary_file <- file.path(TSV_DIR, "2E_mark_filter_summary.txt")  # Original: file.path(output_base, "filter_summary.txt")
  sink(summary_file)
  cat("=== Mark-Filtered Loop Distance Analysis Summary ===\n")
  cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Timepoint:", timepoint, "\n")
  cat("Input file:", input_file, "\n\n")

  cat("=== Overall Loop Counts ===\n")
  cat(sprintf("Total directional loops: %d\n", nrow(loops_directional)))
  cat(sprintf("  Lost in BAP1-KO: %d (%.1f%%)\n",
              sum(loops_directional$direction == "down_in_mutant"),
              100 * mean(loops_directional$direction == "down_in_mutant")))
  cat(sprintf("  Gained in BAP1-KO: %d (%.1f%%)\n\n",
              sum(loops_directional$direction == "up_in_mutant"),
              100 * mean(loops_directional$direction == "up_in_mutant")))

  cat("=== Per-Mark Statistics ===\n\n")

  for (mark_name in names(all_stats)) {
    mark_config <- MARK_CONFIG[[mark_name]]
    stats <- all_stats[[mark_name]]

    cat(sprintf("--- %s (%s) ---\n", mark_config$display_name, mark_config$biological_role))

    # One-anchor stats
    s <- stats$one_anchor
    cat(sprintf("\n[One+ Anchor]\n"))
    cat(sprintf("  Total loops: %d\n", s$total))
    cat(sprintf("  Lost: %d (%.1f%%), Gained: %d (%.1f%%)\n",
                s$n_lost, 100 * s$n_lost / max(s$total, 1),
                s$n_gained, 100 * s$n_gained / max(s$total, 1)))
    if (s$sufficient_data && !is.na(s$median_lost)) {
      cat(sprintf("  Median distance (Lost): %.1f kb\n", s$median_lost / 1000))
      cat(sprintf("  Median distance (Gained): %.1f kb\n", s$median_gained / 1000))
      cat(sprintf("  Median ratio (Lost/Gained): %.2f\n", s$median_lost / s$median_gained))
      cat(sprintf("  KS test p-value: %.2e\n", s$ks_pvalue))
      cat(sprintf("  Wilcoxon p-value: %.2e\n", s$wilcox_pvalue))
    }

    # Both-anchor stats
    s <- stats$both_anchors
    cat(sprintf("\n[Both Anchors]\n"))
    cat(sprintf("  Total loops: %d\n", s$total))
    cat(sprintf("  Lost: %d (%.1f%%), Gained: %d (%.1f%%)\n",
                s$n_lost, 100 * s$n_lost / max(s$total, 1),
                s$n_gained, 100 * s$n_gained / max(s$total, 1)))
    if (s$sufficient_data && !is.na(s$median_lost)) {
      cat(sprintf("  Median distance (Lost): %.1f kb\n", s$median_lost / 1000))
      cat(sprintf("  Median distance (Gained): %.1f kb\n", s$median_gained / 1000))
      cat(sprintf("  Median ratio (Lost/Gained): %.2f\n", s$median_lost / s$median_gained))
      cat(sprintf("  KS test p-value: %.2e\n", s$ks_pvalue))
      cat(sprintf("  Wilcoxon p-value: %.2e\n", s$wilcox_pvalue))
    } else if (!s$sufficient_data) {
      cat(sprintf("  INSUFFICIENT DATA (n=%d < %d minimum)\n", s$total, MIN_SAMPLES))
    }

    cat("\n")
  }

  cat("=== Anchor Overlap Summary ===\n\n")
  for (mark_name in names(all_stats)) {
    mark_config <- MARK_CONFIG[[mark_name]]
    col1 <- mark_config$col1
    col2 <- mark_config$col2

    if (col1 %in% names(loops_directional) && col2 %in% names(loops_directional)) {
      n_anchor1 <- sum(loops_directional[[col1]])
      n_anchor2 <- sum(loops_directional[[col2]])
      n_either <- sum(loops_directional[[col1]] | loops_directional[[col2]])
      n_both <- sum(loops_directional[[col1]] & loops_directional[[col2]])

      cat(sprintf("%s:\n", mark_config$display_name))
      cat(sprintf("  Anchor1+: %d (%.1f%%)\n", n_anchor1, 100 * n_anchor1 / nrow(loops_directional)))
      cat(sprintf("  Anchor2+: %d (%.1f%%)\n", n_anchor2, 100 * n_anchor2 / nrow(loops_directional)))
      cat(sprintf("  Either anchor+: %d (%.1f%%)\n", n_either, 100 * n_either / nrow(loops_directional)))
      cat(sprintf("  Both anchors+: %d (%.1f%%)\n\n", n_both, 100 * n_both / nrow(loops_directional)))
    }
  }

  sink()
  cat("Saved:", summary_file, "\n")

  # ==========================================================================
  # Completion
  # ==========================================================================

  cat("\n=== Analysis Complete for", toupper(timepoint), "===\n")
  cat("Plot directory:", PLOT_DIR, "\n")
  cat("TSV directory:", TSV_DIR, "\n")
  cat("Marks analyzed:", paste(names(all_stats), collapse = ", "), "\n")

  invisible(NULL)
}

# ==============================================================================
# SECTION 6: MAIN EXECUTION
# ==============================================================================

for (tp in TIMEPOINTS_TO_RUN) {
  tryCatch({
    run_timepoint_analysis(tp)
  }, error = function(e) {
    cat(sprintf("\nERROR processing %s timepoint: %s\n", tp, e$message))
  })
}

cat("\n=== All analyses complete ===\n")
cat("Output directories:\n")
cat(sprintf("  Plots: %s\n", PLOT_DIR))
cat(sprintf("  TSVs:  %s\n", TSV_DIR))
cat("\nMarks processed:", paste(MARKS_TO_RUN, collapse = ", "), "\n")
cat("\nEnd time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
