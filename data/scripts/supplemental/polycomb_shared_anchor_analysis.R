# scripts/polycomb_shared_anchor_analysis.R
#
# Task 3c: Polycomb-Specific Shared Anchor Analysis
#
# Purpose:
#   Filter shared anchors to those classified as "Polycomb" or "Repressed_Promoter"
#   and demonstrate that these anchors specifically drive the long-to-short loop
#   switching phenomenon observed in BAP1-KO.
#
# Biological Question:
#   Is the distance difference (lost loops longer than gained loops) LARGER at
#   Polycomb-anchored loops compared to non-Polycomb-anchored loops?
#
# Usage:
#   Rscript scripts/polycomb_shared_anchor_analysis.R                     # Late (default)
#   Rscript scripts/polycomb_shared_anchor_analysis.R --timepoint early   # Early
#   Rscript scripts/polycomb_shared_anchor_analysis.R --timepoint both    # Both
#   Rscript scripts/polycomb_shared_anchor_analysis.R --include-bivalent  # Include Bivalent_Promoter
#
# Output:
#   output/shared_anchor_analysis/{early,late}/polycomb_specific/
#     tables/    - TSV files with filtered loops and statistics
#     plots/     - Multi-format plots (PDF, SVG, JPEG)
#     polycomb_analysis_summary.txt

# ==============================================================================
# 1. SETUP AND PACKAGES
# ==============================================================================

cat("=== Polycomb-Specific Shared Anchor Analysis ===\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

source("data/scripts/_shared/multi_format_output.R")  # Original: source("scripts/utils/multi_format_output.R")

# ==============================================================================
# 2. CONFIGURATION
# ==============================================================================

# Input files by timepoint (from shared_anchor_analysis.R output)
SHARED_LOOP_FILES <- list(
  late = "data/tsvs/supplemental/shared_anchor_loops.tsv",  # Original: output/shared_anchor_analysis/late/tables/shared_anchor_loops.tsv
  early = "output/shared_anchor_analysis/early/tables/shared_anchor_loops.tsv"  # TODO: not in data/
)

# Background files for enrichment testing (all characterized loops)
BACKGROUND_FILES <- list(
  late = "data/upstream/loop_calls/late_characterized_loops.tsv",  # Original: 25042-late_outputs/merged_loops/characterized_loops.tsv
  early = "250831-early_outputs/merged_loops/characterized_loops.tsv"  # TODO: not in data/
)

# Output base directories
OUTPUT_BASE <- list(
  late = list(
    tsvs = "data/tsvs/supplemental",    # Original: output/shared_anchor_analysis/late/polycomb_specific (tables)
    plots = "data/plots/supplemental"    # Original: output/shared_anchor_analysis/late/polycomb_specific (plots)
  ),
  early = list(
    tsvs = "data/tsvs/supplemental",    # Original: output/shared_anchor_analysis/early/polycomb_specific (tables)
    plots = "data/plots/supplemental"    # Original: output/shared_anchor_analysis/early/polycomb_specific (plots)
  )
)

# Polycomb anchor types (can be extended with --include-bivalent)
POLYCOMB_TYPES_BASE <- c("Polycomb", "Repressed_Promoter")

# Color scheme
COLORS <- list(
  lost = "#d73027",       # Red for lost in mutant
  gained = "#4575b4",     # Blue for gained in mutant
  polycomb = "#4DAF4A",   # Green for Polycomb
  non_polycomb = "#A65628" # Brown for non-Polycomb
)

# Direction labels
DIRECTION_LABELS <- c(
  "down_in_mutant" = "Lost",
  "up_in_mutant" = "Gained"
)

# ==============================================================================
# 3. ARGUMENT PARSING
# ==============================================================================

parse_arguments <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  params <- list(
    timepoint = "late",
    include_bivalent = FALSE
  )

  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--timepoint" && i < length(args)) {
      params$timepoint <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--include-bivalent") {
      params$include_bivalent <- TRUE
      i <- i + 1
    } else if (args[i] %in% c("--help", "-h")) {
      cat("Usage: Rscript scripts/polycomb_shared_anchor_analysis.R [OPTIONS]\n\n")
      cat("Options:\n")
      cat("  --timepoint TP     Timepoint: 'early', 'late', or 'both' (default: late)\n")
      cat("  --include-bivalent Include Bivalent_Promoter in Polycomb category\n")
      cat("  --help, -h         Show this help message\n\n")
      cat("Output:\n")
      cat("  output/shared_anchor_analysis/{early,late}/polycomb_specific/\n")
      cat("    tables/  - Filtered loops and statistics\n")
      cat("    plots/   - CDF, density, violin, enrichment plots\n")
      quit(save = "no", status = 0)
    } else {
      i <- i + 1
    }
  }

  if (!params$timepoint %in% c("early", "late", "both")) {
    stop("ERROR: timepoint must be 'early', 'late', or 'both'")
  }

  return(params)
}

PARAMS <- parse_arguments()

if (PARAMS$timepoint == "both") {
  TIMEPOINTS_TO_RUN <- c("late", "early")
} else {
  TIMEPOINTS_TO_RUN <- PARAMS$timepoint
}

# Define Polycomb types based on flag
POLYCOMB_TYPES <- if (PARAMS$include_bivalent) {
  c("Polycomb", "Repressed_Promoter", "Bivalent_Promoter")
} else {
  POLYCOMB_TYPES_BASE
}

cat("Configuration:\n")
cat("  Timepoint(s):", paste(TIMEPOINTS_TO_RUN, collapse = ", "), "\n")
cat("  Polycomb types:", paste(POLYCOMB_TYPES, collapse = ", "), "\n\n")

# ==============================================================================
# 4. HELPER FUNCTIONS
# ==============================================================================

#' Classify loops by Polycomb anchor status
#'
#' @param loops_df Data frame with anchor1_type and anchor2_type columns
#' @param polycomb_types Character vector of types considered "Polycomb"
#' @return Data frame with additional classification columns
classify_polycomb_anchors <- function(loops_df, polycomb_types = POLYCOMB_TYPES) {
  loops_df %>%
    mutate(
      anchor1_is_polycomb = anchor1_type %in% polycomb_types,
      anchor2_is_polycomb = anchor2_type %in% polycomb_types,
      polycomb_anchor_count = as.integer(anchor1_is_polycomb) + as.integer(anchor2_is_polycomb),
      polycomb_status = case_when(
        polycomb_anchor_count == 2 ~ "Both_Polycomb",
        polycomb_anchor_count == 1 ~ "One_Polycomb",
        TRUE ~ "Neither_Polycomb"
      ),
      is_polycomb_anchored = polycomb_anchor_count >= 1
    )
}

#' Compute distance statistics for a subset of loops
#'
#' @param loops_df Data frame with direction and loop_distance columns
#' @param group_name Name for the group (for output)
#' @return List with statistics
compute_distance_stats <- function(loops_df, group_name = "All") {
  lost <- loops_df %>% filter(direction == "down_in_mutant")
  gained <- loops_df %>% filter(direction == "up_in_mutant")

  n_lost <- nrow(lost)
  n_gained <- nrow(gained)

  # Check for sufficient data
  if (n_lost < 2 || n_gained < 2) {
    return(list(
      group = group_name,
      n_total = nrow(loops_df),
      n_lost = n_lost,
      n_gained = n_gained,
      median_lost = if (n_lost > 0) median(lost$loop_distance) else NA,
      median_gained = if (n_gained > 0) median(gained$loop_distance) else NA,
      median_ratio = NA,
      wilcox_p = NA,
      ks_p = NA
    ))
  }

  # Wilcoxon rank-sum test (one-sided: lost > gained)
  wilcox_test <- wilcox.test(
    lost$loop_distance,
    gained$loop_distance,
    alternative = "greater"
  )

  # Kolmogorov-Smirnov test (two-sided)
  ks_test <- ks.test(
    lost$loop_distance,
    gained$loop_distance
  )

  list(
    group = group_name,
    n_total = nrow(loops_df),
    n_lost = n_lost,
    n_gained = n_gained,
    median_lost = median(lost$loop_distance),
    median_gained = median(gained$loop_distance),
    mean_lost = mean(lost$loop_distance),
    mean_gained = mean(gained$loop_distance),
    median_ratio = median(lost$loop_distance) / median(gained$loop_distance),
    wilcox_p = wilcox_test$p.value,
    ks_p = ks_test$p.value,
    ks_d = ks_test$statistic
  )
}

#' Perform interaction test (is distance difference larger at Polycomb?)
#'
#' @param loops_df Data frame with polycomb_status and direction columns
#' @return List with interaction test results
test_polycomb_direction_interaction <- function(loops_df) {
  # Filter to directional loops
  df <- loops_df %>%
    filter(direction %in% c("down_in_mutant", "up_in_mutant")) %>%
    mutate(
      is_lost = direction == "down_in_mutant",
      log_distance = log10(loop_distance)
    )

  if (nrow(df) < 20) {
    return(list(
      interaction_p = NA,
      method = "insufficient data",
      interpretation = "Cannot test interaction with < 20 loops"
    ))
  }

  # Two-way ANOVA on log-transformed distance
  # If interaction is significant, the effect of direction differs by Polycomb status
  model <- tryCatch({
    aov(log_distance ~ is_polycomb_anchored * is_lost, data = df)
  }, error = function(e) NULL)

  if (is.null(model)) {
    return(list(
      interaction_p = NA,
      method = "ANOVA failed",
      interpretation = "Model fitting failed"
    ))
  }

  anova_table <- summary(model)[[1]]

  # Extract interaction p-value
  interaction_p <- anova_table["is_polycomb_anchored:is_lost", "Pr(>F)"]
  main_polycomb_p <- anova_table["is_polycomb_anchored", "Pr(>F)"]
  main_direction_p <- anova_table["is_lost", "Pr(>F)"]

  # Also do a rank-based approach (more robust)
  df$distance_rank <- rank(df$loop_distance)
  rank_model <- tryCatch({
    aov(distance_rank ~ is_polycomb_anchored * is_lost, data = df)
  }, error = function(e) NULL)

  rank_interaction_p <- NA
  if (!is.null(rank_model)) {
    rank_anova <- summary(rank_model)[[1]]
    rank_interaction_p <- rank_anova["is_polycomb_anchored:is_lost", "Pr(>F)"]
  }

  # Interpretation
  interpretation <- if (!is.na(interaction_p) && interaction_p < 0.05) {
    "Significant interaction: Distance difference between lost/gained is larger at Polycomb anchors"
  } else if (!is.na(interaction_p) && interaction_p < 0.10) {
    "Marginal interaction trend (p < 0.10)"
  } else {
    "No significant interaction: Distance difference is similar at Polycomb vs non-Polycomb anchors"
  }

  list(
    method = "Two-way ANOVA (log10 distance)",
    interaction_p = interaction_p,
    rank_interaction_p = rank_interaction_p,
    main_polycomb_p = main_polycomb_p,
    main_direction_p = main_direction_p,
    interpretation = interpretation,
    anova_table = anova_table
  )
}

#' Fisher's exact test for Polycomb enrichment at shared anchors
#'
#' @param shared_loops Loops at shared anchors
#' @param all_loops All characterized loops (background)
#' @return List with enrichment statistics
test_polycomb_enrichment <- function(shared_loops, all_loops) {
  # Classify both datasets
  shared_classified <- classify_polycomb_anchors(shared_loops)
  all_classified <- classify_polycomb_anchors(all_loops)

  # Count Polycomb-anchored loops in each set
  shared_polycomb <- sum(shared_classified$is_polycomb_anchored)
  shared_other <- sum(!shared_classified$is_polycomb_anchored)
  all_polycomb <- sum(all_classified$is_polycomb_anchored)
  all_other <- sum(!all_classified$is_polycomb_anchored)

  # Background = all loops NOT in shared
  bg_polycomb <- all_polycomb - shared_polycomb
  bg_other <- all_other - shared_other

  # Ensure non-negative (in case shared has loops not in all)
  bg_polycomb <- max(bg_polycomb, 0)
  bg_other <- max(bg_other, 0)

  # 2x2 contingency table
  #                    Shared    Background
  # Polycomb           a         b
  # Non-Polycomb       c         d

  mat <- matrix(c(shared_polycomb, bg_polycomb,
                  shared_other, bg_other),
                nrow = 2, byrow = TRUE,
                dimnames = list(c("Polycomb", "Non-Polycomb"),
                               c("Shared", "Background")))

  # Fisher's exact test
  fisher_result <- fisher.test(mat)

  list(
    contingency_table = mat,
    odds_ratio = fisher_result$estimate,
    conf_int = fisher_result$conf.int,
    p_value = fisher_result$p.value,
    shared_polycomb_pct = 100 * shared_polycomb / (shared_polycomb + shared_other),
    background_polycomb_pct = 100 * bg_polycomb / (bg_polycomb + bg_other)
  )
}

# ==============================================================================
# 5. VISUALIZATION FUNCTIONS
# ==============================================================================

#' Generate CDF plot for Polycomb-anchored loops
#'
#' @param loops_df Data frame with classified loops
#' @param polycomb_filter "polycomb" or "non_polycomb"
#' @param stats_list Statistics from compute_distance_stats
#' @return ggplot object
create_cdf_plot <- function(loops_df, polycomb_filter, stats_list) {
  if (polycomb_filter == "polycomb") {
    df <- loops_df %>% filter(is_polycomb_anchored)
    title_suffix <- "Polycomb-Anchored"
    subtitle <- "Loops with at least one Polycomb/Repressed_Promoter anchor"
  } else {
    df <- loops_df %>% filter(!is_polycomb_anchored)
    title_suffix <- "Non-Polycomb"
    subtitle <- "Loops without Polycomb/Repressed_Promoter anchors"
  }

  # Filter to directional only
  df <- df %>%
    filter(direction %in% c("down_in_mutant", "up_in_mutant")) %>%
    mutate(
      direction_label = ifelse(direction == "down_in_mutant", "Lost", "Gained"),
      loop_distance_kb = loop_distance / 1000
    )

  if (nrow(df) < 10) {
    return(NULL)
  }

  # Get stats
  stats <- stats_list

  p <- ggplot(df, aes(x = loop_distance_kb, color = direction_label)) +
    stat_ecdf(geom = "step", linewidth = 1.2) +
    scale_color_manual(
      values = c("Lost" = COLORS$lost, "Gained" = COLORS$gained),
      name = "Direction"
    ) +
    scale_x_log10(
      labels = comma,
      breaks = c(10, 100, 1000, 10000),
      limits = c(10, 100000)
    ) +
    geom_vline(xintercept = stats$median_lost / 1000, color = COLORS$lost,
               linetype = "dashed", alpha = 0.7) +
    geom_vline(xintercept = stats$median_gained / 1000, color = COLORS$gained,
               linetype = "dashed", alpha = 0.7) +
    annotate("text", x = 50000, y = 0.12,
             label = sprintf("KS test p = %.2e", stats$ks_p),
             hjust = 1, size = 4) +
    annotate("text", x = 50000, y = 0.05,
             label = sprintf("Median: Lost = %.0f kb, Gained = %.0f kb\nRatio = %.1fx",
                            stats$median_lost / 1000, stats$median_gained / 1000,
                            stats$median_ratio),
             hjust = 1, size = 3.5) +
    annotate("text", x = 50000, y = 0.20,
             label = sprintf("n = %d loops (%d lost, %d gained)",
                            stats$n_total, stats$n_lost, stats$n_gained),
             hjust = 1, size = 3.5, fontface = "italic") +
    labs(
      title = sprintf("Loop Distance CDF: %s Loops at Shared Anchors", title_suffix),
      subtitle = subtitle,
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

  return(p)
}

#' Generate faceted CDF comparison plot
#'
#' @param loops_df Data frame with classified loops
#' @param polycomb_stats Stats for Polycomb subset
#' @param non_polycomb_stats Stats for non-Polycomb subset
#' @return ggplot object
create_faceted_cdf_plot <- function(loops_df, polycomb_stats, non_polycomb_stats) {
  # Prepare data with proper faceting
  df <- loops_df %>%
    filter(direction %in% c("down_in_mutant", "up_in_mutant")) %>%
    mutate(
      direction_label = ifelse(direction == "down_in_mutant", "Lost", "Gained"),
      loop_distance_kb = loop_distance / 1000,
      polycomb_label = ifelse(is_polycomb_anchored,
                              "Polycomb-Anchored",
                              "Non-Polycomb")
    )

  if (nrow(df) < 20) {
    return(NULL)
  }

  # Create annotation data frame
  annot_df <- tibble(
    polycomb_label = c("Polycomb-Anchored", "Non-Polycomb"),
    median_lost = c(polycomb_stats$median_lost, non_polycomb_stats$median_lost),
    median_gained = c(polycomb_stats$median_gained, non_polycomb_stats$median_gained),
    ratio = c(polycomb_stats$median_ratio, non_polycomb_stats$median_ratio),
    n = c(polycomb_stats$n_total, non_polycomb_stats$n_total),
    ks_p = c(polycomb_stats$ks_p, non_polycomb_stats$ks_p)
  )

  p <- ggplot(df, aes(x = loop_distance_kb, color = direction_label)) +
    stat_ecdf(geom = "step", linewidth = 1.0) +
    scale_color_manual(
      values = c("Lost" = COLORS$lost, "Gained" = COLORS$gained),
      name = ""
    ) +
    scale_x_log10(
      labels = comma,
      breaks = c(10, 100, 1000, 10000),
      limits = c(10, 100000)
    ) +
    facet_wrap(~ polycomb_label, ncol = 2) +
    geom_text(
      data = annot_df,
      aes(x = 30000, y = 0.15, label = sprintf("Ratio = %.1fx\np = %.2e", ratio, ks_p)),
      inherit.aes = FALSE, hjust = 1, size = 3.5
    ) +
    labs(
      title = "Loop Distance CDF: Polycomb vs Non-Polycomb Anchors",
      subtitle = "Distance difference is amplified at Polycomb-anchored loops",
      x = "Loop Distance (kb, log scale)",
      y = "Cumulative Proportion"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey40"),
      legend.position = "top",
      strip.text = element_text(face = "bold", size = 11),
      panel.grid.minor = element_blank()
    )

  return(p)
}

#' Generate density plot for a subset
#'
#' @param loops_df Data frame with classified loops
#' @param polycomb_filter "polycomb" or "non_polycomb"
#' @param stats_list Statistics from compute_distance_stats
#' @return ggplot object
create_density_plot <- function(loops_df, polycomb_filter, stats_list) {
  if (polycomb_filter == "polycomb") {
    df <- loops_df %>% filter(is_polycomb_anchored)
    title_suffix <- "Polycomb-Anchored"
  } else {
    df <- loops_df %>% filter(!is_polycomb_anchored)
    title_suffix <- "Non-Polycomb"
  }

  df <- df %>%
    filter(direction %in% c("down_in_mutant", "up_in_mutant")) %>%
    mutate(
      direction_label = ifelse(direction == "down_in_mutant", "Lost", "Gained"),
      loop_distance_kb = loop_distance / 1000
    )

  if (nrow(df) < 10) {
    return(NULL)
  }

  stats <- stats_list

  p <- ggplot(df, aes(x = loop_distance_kb, fill = direction_label)) +
    geom_density(alpha = 0.5, color = "black", linewidth = 0.5) +
    geom_vline(xintercept = stats$median_lost / 1000, color = COLORS$lost,
               linetype = "dashed", linewidth = 1) +
    geom_vline(xintercept = stats$median_gained / 1000, color = COLORS$gained,
               linetype = "dashed", linewidth = 1) +
    scale_fill_manual(
      values = c("Lost" = COLORS$lost, "Gained" = COLORS$gained),
      name = ""
    ) +
    scale_x_log10(
      labels = comma,
      breaks = c(10, 100, 1000, 10000),
      limits = c(10, 100000)
    ) +
    annotate("text", x = 15, y = Inf,
             label = sprintf("Wilcoxon p = %.2e\nn = %d", stats$wilcox_p, stats$n_total),
             hjust = 0, vjust = 1.5, size = 4) +
    labs(
      title = sprintf("Loop Distance Density: %s", title_suffix),
      subtitle = sprintf("Median: Lost = %.0f kb, Gained = %.0f kb (Ratio = %.1fx)",
                        stats$median_lost / 1000, stats$median_gained / 1000,
                        stats$median_ratio),
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

  return(p)
}

#' Generate 2x2 violin interaction plot
#'
#' @param loops_df Data frame with classified loops
#' @param interaction_result Result from test_polycomb_direction_interaction
#' @return ggplot object
create_interaction_violin_plot <- function(loops_df, interaction_result) {
  df <- loops_df %>%
    filter(direction %in% c("down_in_mutant", "up_in_mutant")) %>%
    mutate(
      direction_label = factor(
        ifelse(direction == "down_in_mutant", "Lost", "Gained"),
        levels = c("Lost", "Gained")
      ),
      polycomb_label = factor(
        ifelse(is_polycomb_anchored, "Polycomb-Anchored", "Non-Polycomb"),
        levels = c("Polycomb-Anchored", "Non-Polycomb")
      ),
      loop_distance_mb = loop_distance / 1e6
    )

  if (nrow(df) < 20) {
    return(NULL)
  }

  # Calculate group stats for annotation
  group_stats <- df %>%
    group_by(polycomb_label, direction_label) %>%
    summarise(
      median_dist = median(loop_distance_mb),
      n = n(),
      .groups = "drop"
    )

  p <- ggplot(df, aes(x = direction_label, y = loop_distance_mb, fill = direction_label)) +
    geom_violin(alpha = 0.7, trim = FALSE, color = "black", linewidth = 0.5) +
    geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA,
                 color = "black", linewidth = 0.5) +
    facet_wrap(~ polycomb_label) +
    scale_fill_manual(
      values = c("Lost" = COLORS$lost, "Gained" = COLORS$gained),
      guide = "none"
    ) +
    scale_y_log10(labels = comma) +
    geom_text(
      data = group_stats,
      aes(y = median_dist, label = sprintf("%.0f kb\nn=%d", median_dist * 1000, n)),
      vjust = -1, size = 3, fontface = "bold"
    ) +
    labs(
      title = "Loop Distance by Polycomb Status and Direction",
      subtitle = sprintf("Interaction test: p = %.2e\n%s",
                        interaction_result$interaction_p,
                        interaction_result$interpretation),
      x = NULL,
      y = "Loop Distance (Mb, log scale)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"),
      strip.text = element_text(face = "bold", size = 11),
      panel.grid.minor = element_blank()
    )

  return(p)
}

#' Generate anchor type distribution plot
#'
#' @param loops_df Data frame with anchor type columns
#' @return ggplot object
create_anchor_type_distribution_plot <- function(loops_df) {
  df <- loops_df %>%
    filter(direction %in% c("down_in_mutant", "up_in_mutant")) %>%
    mutate(
      direction_label = ifelse(direction == "down_in_mutant", "Lost", "Gained")
    )

  # Gather anchor types (both anchors)
  type_data <- bind_rows(
    df %>% select(direction_label, anchor_type = anchor1_type),
    df %>% select(direction_label, anchor_type = anchor2_type)
  )

  # Calculate proportions
  type_summary <- type_data %>%
    count(direction_label, anchor_type) %>%
    group_by(direction_label) %>%
    mutate(
      prop = n / sum(n),
      is_polycomb = anchor_type %in% POLYCOMB_TYPES
    ) %>%
    ungroup()

  # Order anchor types by total frequency
  type_order <- type_summary %>%
    group_by(anchor_type) %>%
    summarise(total = sum(n)) %>%
    arrange(desc(total)) %>%
    pull(anchor_type)

  type_summary$anchor_type <- factor(type_summary$anchor_type, levels = rev(type_order))

  p <- ggplot(type_summary, aes(x = direction_label, y = prop, fill = anchor_type)) +
    geom_col(color = "black", linewidth = 0.3) +
    scale_fill_brewer(palette = "Set2", name = "Anchor Type") +
    scale_y_continuous(labels = percent) +
    labs(
      title = "Anchor Type Distribution at Shared Anchors",
      subtitle = "Lost loops have higher Polycomb/Repressed_Promoter proportion",
      x = NULL,
      y = "Proportion of Anchors"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey40"),
      legend.position = "right"
    )

  return(p)
}

#' Generate Polycomb enrichment dotplot
#'
#' @param enrichment_result Result from test_polycomb_enrichment
#' @return ggplot object
create_enrichment_dotplot <- function(enrichment_result) {
  df <- tibble(
    comparison = "Shared vs Background",
    log2OR = log2(enrichment_result$odds_ratio),
    log2_ci_low = log2(enrichment_result$conf_int[1]),
    log2_ci_high = log2(enrichment_result$conf_int[2]),
    p_value = enrichment_result$p_value,
    sig = ifelse(enrichment_result$p_value < 0.05, "p < 0.05", "NS")
  )

  p <- ggplot(df, aes(x = log2OR, y = comparison, color = sig)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = log2_ci_low, xmax = log2_ci_high), height = 0.2, linewidth = 1) +
    geom_point(size = 5) +
    scale_color_manual(values = c("p < 0.05" = "#E41A1C", "NS" = "gray50"), name = "") +
    annotate("text", x = df$log2OR + 0.1, y = 1.2,
             label = sprintf("OR = %.2f\np = %.2e",
                            enrichment_result$odds_ratio,
                            enrichment_result$p_value),
             hjust = 0, size = 4) +
    labs(
      title = "Polycomb Enrichment at Shared Anchors",
      subtitle = sprintf("Shared: %.1f%% Polycomb, Background: %.1f%% Polycomb",
                        enrichment_result$shared_polycomb_pct,
                        enrichment_result$background_polycomb_pct),
      x = "log2(Odds Ratio)",
      y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey40"),
      legend.position = "bottom",
      axis.text.y = element_text(face = "bold")
    ) +
    coord_cartesian(xlim = c(-1, max(df$log2_ci_high, 1) + 0.5))

  return(p)
}

#' Create summary panel figure
#'
#' @param cdf_polycomb CDF plot for Polycomb subset
#' @param cdf_non_polycomb CDF plot for non-Polycomb subset
#' @param violin_plot Interaction violin plot
#' @param enrichment_plot Enrichment dotplot
#' @return Combined patchwork plot
create_summary_panel <- function(cdf_polycomb, cdf_non_polycomb, violin_plot, enrichment_plot) {
  # Handle NULL plots
  plots <- list(cdf_polycomb, cdf_non_polycomb, violin_plot, enrichment_plot)
  valid_plots <- plots[!sapply(plots, is.null)]

  if (length(valid_plots) == 0) {
    return(NULL)
  }

  if (length(valid_plots) == 4) {
    panel <- (cdf_polycomb | cdf_non_polycomb) / (violin_plot | enrichment_plot) +
      plot_annotation(
        title = "Polycomb-Specific Loop Switching at Shared Anchors",
        subtitle = "Long-range loop loss and short-range loop gain occur preferentially at Polycomb/Repressed_Promoter anchors",
        theme = theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
          plot.subtitle = element_text(hjust = 0.5, size = 12, color = "grey40")
        )
      ) +
      plot_layout(guides = "collect") &
      theme(legend.position = "bottom")
  } else {
    # Combine whatever is available
    panel <- wrap_plots(valid_plots, ncol = 2) +
      plot_annotation(
        title = "Polycomb-Specific Loop Switching at Shared Anchors",
        theme = theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
        )
      )
  }

  return(panel)
}

# ==============================================================================
# 6. MAIN ANALYSIS FUNCTION
# ==============================================================================

run_polycomb_analysis <- function(timepoint) {
  cat("\n")
  cat("============================================================\n")
  cat(sprintf("Processing %s timepoint\n", toupper(timepoint)))
  cat("============================================================\n\n")

  # Set paths
  shared_file <- SHARED_LOOP_FILES[[timepoint]]
  background_file <- BACKGROUND_FILES[[timepoint]]
  output_base <- OUTPUT_BASE[[timepoint]]

  # Create output directories
  tables_dir <- output_base$tsvs   # Original: file.path(output_dir, "tables")
  plots_dir <- output_base$plots   # Original: file.path(output_dir, "plots")
  output_dir <- output_base$tsvs   # For summary report path compatibility
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

  cat("Input file:", shared_file, "\n")
  cat("Background file:", background_file, "\n")
  cat("Output TSVs:", tables_dir, "\n")
  cat("Output plots:", plots_dir, "\n\n")

  # ============================================================================
  # Step 1: Load and validate data
  # ============================================================================

  cat("=== Step 1: Loading Data ===\n")

  if (!file.exists(shared_file)) {
    stop(sprintf("Shared loops file not found: %s\nRun shared_anchor_analysis.R first.", shared_file))
  }

  shared_loops <- read_tsv(shared_file, show_col_types = FALSE)
  cat("Loaded", nrow(shared_loops), "loops at shared anchors\n")

  # Check for required columns
  required_cols <- c("anchor1_type", "anchor2_type", "direction", "loop_distance")
  missing_cols <- setdiff(required_cols, colnames(shared_loops))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }

  # Load background if available
  all_loops <- NULL
  if (file.exists(background_file)) {
    all_loops <- read_tsv(background_file, show_col_types = FALSE)
    cat("Loaded", nrow(all_loops), "background loops\n")
  } else {
    cat("WARNING: Background file not found, enrichment test will be skipped\n")
  }

  # Filter to directional loops
  shared_directional <- shared_loops %>%
    filter(direction %in% c("down_in_mutant", "up_in_mutant"))

  cat("\nDirectional loops at shared anchors:", nrow(shared_directional), "\n")
  cat("  - Lost (down_in_mutant):", sum(shared_directional$direction == "down_in_mutant"), "\n")
  cat("  - Gained (up_in_mutant):", sum(shared_directional$direction == "up_in_mutant"), "\n\n")

  # Check for sufficient data
  if (nrow(shared_directional) < 20) {
    cat("\nWARNING: Insufficient data (< 20 directional loops) for robust analysis\n")
    cat("Results should be interpreted with caution.\n\n")
  }

  # ============================================================================
  # Step 2: Classify by Polycomb status
  # ============================================================================

  cat("=== Step 2: Classifying by Polycomb Status ===\n")

  shared_classified <- classify_polycomb_anchors(shared_directional)

  # Summarize
  polycomb_loops <- shared_classified %>% filter(is_polycomb_anchored)
  non_polycomb_loops <- shared_classified %>% filter(!is_polycomb_anchored)
  both_polycomb <- shared_classified %>% filter(polycomb_anchor_count == 2)

  cat("Polycomb-anchored (at least one anchor):", nrow(polycomb_loops), "\n")
  cat("  - Lost:", sum(polycomb_loops$direction == "down_in_mutant"), "\n")
  cat("  - Gained:", sum(polycomb_loops$direction == "up_in_mutant"), "\n")
  cat("Both anchors Polycomb:", nrow(both_polycomb), "\n")
  cat("Non-Polycomb:", nrow(non_polycomb_loops), "\n")
  cat("  - Lost:", sum(non_polycomb_loops$direction == "down_in_mutant"), "\n")
  cat("  - Gained:", sum(non_polycomb_loops$direction == "up_in_mutant"), "\n\n")

  # ============================================================================
  # Step 3: Compute distance statistics
  # ============================================================================

  cat("=== Step 3: Computing Distance Statistics ===\n")

  stats_all <- compute_distance_stats(shared_classified, "All Shared")
  stats_polycomb <- compute_distance_stats(polycomb_loops, "Polycomb-Anchored")
  stats_non_polycomb <- compute_distance_stats(non_polycomb_loops, "Non-Polycomb")
  stats_both_polycomb <- compute_distance_stats(both_polycomb, "Both-Polycomb")

  cat("\n--- All Shared Loops ---\n")
  cat(sprintf("  n = %d (Lost: %d, Gained: %d)\n", stats_all$n_total, stats_all$n_lost, stats_all$n_gained))
  cat(sprintf("  Median Lost: %.0f kb, Median Gained: %.0f kb\n",
              stats_all$median_lost / 1000, stats_all$median_gained / 1000))
  cat(sprintf("  Ratio: %.2fx\n", stats_all$median_ratio))
  cat(sprintf("  Wilcoxon p (lost > gained): %.2e\n", stats_all$wilcox_p))

  cat("\n--- Polycomb-Anchored ---\n")
  cat(sprintf("  n = %d (Lost: %d, Gained: %d)\n", stats_polycomb$n_total, stats_polycomb$n_lost, stats_polycomb$n_gained))
  if (!is.na(stats_polycomb$median_lost)) {
    cat(sprintf("  Median Lost: %.0f kb, Median Gained: %.0f kb\n",
                stats_polycomb$median_lost / 1000, stats_polycomb$median_gained / 1000))
    cat(sprintf("  Ratio: %.2fx\n", stats_polycomb$median_ratio))
    cat(sprintf("  Wilcoxon p: %.2e\n", stats_polycomb$wilcox_p))
  }

  cat("\n--- Non-Polycomb ---\n")
  cat(sprintf("  n = %d (Lost: %d, Gained: %d)\n", stats_non_polycomb$n_total, stats_non_polycomb$n_lost, stats_non_polycomb$n_gained))
  if (!is.na(stats_non_polycomb$median_lost)) {
    cat(sprintf("  Median Lost: %.0f kb, Median Gained: %.0f kb\n",
                stats_non_polycomb$median_lost / 1000, stats_non_polycomb$median_gained / 1000))
    cat(sprintf("  Ratio: %.2fx\n", stats_non_polycomb$median_ratio))
    cat(sprintf("  Wilcoxon p: %.2e\n", stats_non_polycomb$wilcox_p))
  }
  cat("\n")

  # ============================================================================
  # Step 4: Interaction test
  # ============================================================================

  cat("=== Step 4: Testing Polycomb x Direction Interaction ===\n")

  interaction_result <- test_polycomb_direction_interaction(shared_classified)

  cat(sprintf("Method: %s\n", interaction_result$method))
  cat(sprintf("Interaction p-value: %.4f\n", interaction_result$interaction_p))
  if (!is.na(interaction_result$rank_interaction_p)) {
    cat(sprintf("Rank-based interaction p: %.4f\n", interaction_result$rank_interaction_p))
  }
  cat(sprintf("Interpretation: %s\n\n", interaction_result$interpretation))

  # ============================================================================
  # Step 5: Enrichment test (if background available)
  # ============================================================================

  enrichment_result <- NULL
  if (!is.null(all_loops)) {
    cat("=== Step 5: Testing Polycomb Enrichment ===\n")

    enrichment_result <- test_polycomb_enrichment(shared_loops, all_loops)

    cat("Contingency table:\n")
    print(enrichment_result$contingency_table)
    cat(sprintf("\nOdds Ratio: %.2f (95%% CI: %.2f - %.2f)\n",
                enrichment_result$odds_ratio,
                enrichment_result$conf_int[1],
                enrichment_result$conf_int[2]))
    cat(sprintf("Fisher p-value: %.2e\n", enrichment_result$p_value))
    cat(sprintf("Shared anchors: %.1f%% Polycomb\n", enrichment_result$shared_polycomb_pct))
    cat(sprintf("Background: %.1f%% Polycomb\n\n", enrichment_result$background_polycomb_pct))
  }

  # ============================================================================
  # Step 6: Save tables
  # ============================================================================

  cat("=== Step 6: Saving Tables ===\n")

  # Polycomb-anchored loops
  write_tsv(polycomb_loops, file.path(tables_dir, "polycomb_shared_loops.tsv"))  # Original: tables/polycomb_shared_loops.tsv
  cat("  Saved: polycomb_shared_loops.tsv\n")

  # Both-Polycomb loops
  if (nrow(both_polycomb) > 0) {
    write_tsv(both_polycomb, file.path(tables_dir, "polycomb_both_shared_loops.tsv"))  # Original: tables/both_polycomb_shared_loops.tsv
    cat("  Saved: polycomb_both_shared_loops.tsv\n")
  }

  # Distance statistics
  stats_df <- bind_rows(
    as_tibble(stats_all),
    as_tibble(stats_polycomb),
    as_tibble(stats_non_polycomb),
    as_tibble(stats_both_polycomb)
  )
  write_tsv(stats_df, file.path(tables_dir, "polycomb_distance_statistics.tsv"))  # Original: tables/distance_statistics.tsv
  cat("  Saved: polycomb_distance_statistics.tsv\n")

  # Interaction test results
  interaction_df <- tibble(
    method = interaction_result$method,
    interaction_p = interaction_result$interaction_p,
    rank_interaction_p = interaction_result$rank_interaction_p,
    main_polycomb_p = interaction_result$main_polycomb_p,
    main_direction_p = interaction_result$main_direction_p,
    interpretation = interaction_result$interpretation
  )
  write_tsv(interaction_df, file.path(tables_dir, "polycomb_interaction_test.tsv"))  # Original: tables/interaction_test.tsv
  cat("  Saved: polycomb_interaction_test.tsv\n")

  # Enrichment statistics
  if (!is.null(enrichment_result)) {
    enrichment_df <- tibble(
      odds_ratio = enrichment_result$odds_ratio,
      ci_low = enrichment_result$conf_int[1],
      ci_high = enrichment_result$conf_int[2],
      p_value = enrichment_result$p_value,
      shared_polycomb_pct = enrichment_result$shared_polycomb_pct,
      background_polycomb_pct = enrichment_result$background_polycomb_pct
    )
    write_tsv(enrichment_df, file.path(tables_dir, "polycomb_enrichment_statistics.tsv"))  # Original: tables/enrichment_statistics.tsv
    cat("  Saved: polycomb_enrichment_statistics.tsv\n")
  }
  cat("\n")

  # ============================================================================
  # Step 7: Generate plots
  # ============================================================================

  cat("=== Step 7: Generating Plots ===\n")

  # CDF plots
  cdf_polycomb <- create_cdf_plot(shared_classified, "polycomb", stats_polycomb)
  if (!is.null(cdf_polycomb)) {
    save_multiformat_ggplot(cdf_polycomb, file.path(plots_dir, "01_cdf_polycomb_anchored"),
                            width = 8, height = 6)
  }

  cdf_non_polycomb <- create_cdf_plot(shared_classified, "non_polycomb", stats_non_polycomb)
  if (!is.null(cdf_non_polycomb)) {
    save_multiformat_ggplot(cdf_non_polycomb, file.path(plots_dir, "01_cdf_non_polycomb_anchored"),
                            width = 8, height = 6)
  }

  # Faceted CDF comparison
  cdf_faceted <- create_faceted_cdf_plot(shared_classified, stats_polycomb, stats_non_polycomb)
  if (!is.null(cdf_faceted)) {
    save_multiformat_ggplot(cdf_faceted, file.path(plots_dir, "01_cdf_comparison_faceted"),
                            width = 12, height = 5)
  }

  # Density plots
  density_polycomb <- create_density_plot(shared_classified, "polycomb", stats_polycomb)
  if (!is.null(density_polycomb)) {
    save_multiformat_ggplot(density_polycomb, file.path(plots_dir, "02_density_polycomb_anchored"),
                            width = 8, height = 6)
  }

  density_non_polycomb <- create_density_plot(shared_classified, "non_polycomb", stats_non_polycomb)
  if (!is.null(density_non_polycomb)) {
    save_multiformat_ggplot(density_non_polycomb, file.path(plots_dir, "02_density_non_polycomb_anchored"),
                            width = 8, height = 6)
  }

  # Anchor type distribution
  type_plot <- create_anchor_type_distribution_plot(shared_classified)
  if (!is.null(type_plot)) {
    save_multiformat_ggplot(type_plot, file.path(plots_dir, "03_anchor_type_distribution"),
                            width = 8, height = 6)
  }

  # Enrichment dotplot
  enrichment_plot <- NULL
  if (!is.null(enrichment_result)) {
    enrichment_plot <- create_enrichment_dotplot(enrichment_result)
    save_multiformat_ggplot(enrichment_plot, file.path(plots_dir, "03_polycomb_enrichment_dotplot"),
                            width = 7, height = 4)
  }

  # Interaction violin plot
  violin_plot <- create_interaction_violin_plot(shared_classified, interaction_result)
  if (!is.null(violin_plot)) {
    save_multiformat_ggplot(violin_plot, file.path(plots_dir, "04_distance_violin_interaction"),
                            width = 10, height = 6)
  }

  # Summary panel
  summary_panel <- create_summary_panel(cdf_polycomb, cdf_non_polycomb, violin_plot, enrichment_plot)
  if (!is.null(summary_panel)) {
    save_multiformat_ggplot(summary_panel, file.path(plots_dir, "05_summary_panel"),
                            width = 14, height = 10)
  }

  cat("\n")

  # ============================================================================
  # Step 8: Write summary report
  # ============================================================================

  cat("=== Step 8: Writing Summary Report ===\n")

  summary_file <- file.path(output_dir, "polycomb_analysis_summary.txt")

  sink(summary_file)
  cat("===============================================================\n")
  cat("POLYCOMB-SPECIFIC SHARED ANCHOR ANALYSIS\n")
  cat("===============================================================\n\n")
  cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Timepoint:", timepoint, "\n")
  cat("Polycomb types:", paste(POLYCOMB_TYPES, collapse = ", "), "\n\n")

  cat("--- INPUT FILES ---\n")
  cat("Shared loops:", shared_file, "\n")
  if (!is.null(all_loops)) {
    cat("Background:", background_file, "\n")
  }
  cat("\n")

  cat("--- LOOP COUNTS ---\n")
  cat(sprintf("Total directional loops at shared anchors: %d\n", nrow(shared_classified)))
  cat(sprintf("  - Lost (down_in_mutant): %d\n", sum(shared_classified$direction == "down_in_mutant")))
  cat(sprintf("  - Gained (up_in_mutant): %d\n", sum(shared_classified$direction == "up_in_mutant")))
  cat("\n")
  cat(sprintf("Polycomb-anchored: %d (%.1f%%)\n",
              nrow(polycomb_loops),
              100 * nrow(polycomb_loops) / nrow(shared_classified)))
  cat(sprintf("  - Lost: %d\n", sum(polycomb_loops$direction == "down_in_mutant")))
  cat(sprintf("  - Gained: %d\n", sum(polycomb_loops$direction == "up_in_mutant")))
  cat(sprintf("Both anchors Polycomb: %d\n", nrow(both_polycomb)))
  cat("\n")
  cat(sprintf("Non-Polycomb: %d (%.1f%%)\n",
              nrow(non_polycomb_loops),
              100 * nrow(non_polycomb_loops) / nrow(shared_classified)))
  cat(sprintf("  - Lost: %d\n", sum(non_polycomb_loops$direction == "down_in_mutant")))
  cat(sprintf("  - Gained: %d\n", sum(non_polycomb_loops$direction == "up_in_mutant")))
  cat("\n")

  cat("--- DISTANCE STATISTICS ---\n\n")

  for (stats in list(stats_all, stats_polycomb, stats_non_polycomb, stats_both_polycomb)) {
    cat(sprintf("--- %s ---\n", stats$group))
    cat(sprintf("  n = %d (Lost: %d, Gained: %d)\n", stats$n_total, stats$n_lost, stats$n_gained))
    if (!is.na(stats$median_lost)) {
      cat(sprintf("  Median distance (Lost): %.0f kb (%.2f Mb)\n",
                  stats$median_lost / 1000, stats$median_lost / 1e6))
      cat(sprintf("  Median distance (Gained): %.0f kb (%.2f Mb)\n",
                  stats$median_gained / 1000, stats$median_gained / 1e6))
      cat(sprintf("  Median ratio (Lost/Gained): %.2fx\n", stats$median_ratio))
      cat(sprintf("  Wilcoxon p (lost > gained): %.2e\n", stats$wilcox_p))
      cat(sprintf("  KS test p: %.2e (D = %.3f)\n", stats$ks_p, stats$ks_d))
    }
    cat("\n")
  }

  cat("--- INTERACTION TEST ---\n")
  cat(sprintf("Method: %s\n", interaction_result$method))
  cat(sprintf("Interaction p-value: %.4f\n", interaction_result$interaction_p))
  if (!is.na(interaction_result$rank_interaction_p)) {
    cat(sprintf("Rank-based interaction p: %.4f\n", interaction_result$rank_interaction_p))
  }
  if (!is.na(interaction_result$main_polycomb_p)) {
    cat(sprintf("Main effect (Polycomb) p: %.4f\n", interaction_result$main_polycomb_p))
    cat(sprintf("Main effect (Direction) p: %.4f\n", interaction_result$main_direction_p))
  }
  cat(sprintf("Interpretation: %s\n\n", interaction_result$interpretation))

  if (!is.null(enrichment_result)) {
    cat("--- POLYCOMB ENRICHMENT ---\n")
    cat(sprintf("Odds Ratio: %.2f (95%% CI: %.2f - %.2f)\n",
                enrichment_result$odds_ratio,
                enrichment_result$conf_int[1],
                enrichment_result$conf_int[2]))
    cat(sprintf("Fisher p-value: %.2e\n", enrichment_result$p_value))
    cat(sprintf("Shared anchors: %.1f%% Polycomb-anchored\n", enrichment_result$shared_polycomb_pct))
    cat(sprintf("Background: %.1f%% Polycomb-anchored\n\n", enrichment_result$background_polycomb_pct))
  }

  cat("--- KEY FINDINGS ---\n")

  if (!is.na(stats_polycomb$median_ratio) && !is.na(stats_non_polycomb$median_ratio)) {
    cat(sprintf("1. Distance ratio at Polycomb anchors: %.2fx (Lost/Gained)\n", stats_polycomb$median_ratio))
    cat(sprintf("2. Distance ratio at non-Polycomb anchors: %.2fx (Lost/Gained)\n", stats_non_polycomb$median_ratio))
    cat(sprintf("3. Ratio amplification: %.2fx larger at Polycomb\n",
                stats_polycomb$median_ratio / stats_non_polycomb$median_ratio))
  }

  if (!is.na(interaction_result$interaction_p) && interaction_result$interaction_p < 0.05) {
    cat("4. Significant interaction confirms Polycomb anchors drive the loop switching pattern\n")
  }

  cat("\n")
  cat("--- OUTPUT FILES ---\n")
  cat("Tables: ", tables_dir, "\n")
  cat("  - polycomb_shared_loops.tsv\n")
  cat("  - both_polycomb_shared_loops.tsv\n")
  cat("  - distance_statistics.tsv\n")
  cat("  - interaction_test.tsv\n")
  if (!is.null(enrichment_result)) {
    cat("  - enrichment_statistics.tsv\n")
  }
  cat("\nPlots: ", plots_dir, "\n")
  cat("  - 01_cdf_polycomb_anchored\n")
  cat("  - 01_cdf_non_polycomb_anchored\n")
  cat("  - 01_cdf_comparison_faceted\n")
  cat("  - 02_density_polycomb_anchored\n")
  cat("  - 02_density_non_polycomb_anchored\n")
  cat("  - 03_anchor_type_distribution\n")
  cat("  - 03_polycomb_enrichment_dotplot\n")
  cat("  - 04_distance_violin_interaction\n")
  cat("  - 05_summary_panel\n")

  cat("\n")
  cat("===============================================================\n")
  cat("END REPORT\n")
  cat("===============================================================\n")
  sink()

  cat("  Saved: polycomb_analysis_summary.txt\n")

  # ============================================================================
  # Return results
  # ============================================================================

  cat("\n=== Analysis Complete for", toupper(timepoint), "===\n")
  cat("Output directory:", output_dir, "\n")

  return(list(
    shared_classified = shared_classified,
    polycomb_loops = polycomb_loops,
    non_polycomb_loops = non_polycomb_loops,
    stats_polycomb = stats_polycomb,
    stats_non_polycomb = stats_non_polycomb,
    interaction_result = interaction_result,
    enrichment_result = enrichment_result
  ))
}

# ==============================================================================
# 7. MAIN EXECUTION
# ==============================================================================

results <- list()

for (tp in TIMEPOINTS_TO_RUN) {
  tryCatch({
    results[[tp]] <- run_polycomb_analysis(tp)
  }, error = function(e) {
    cat(sprintf("\nERROR processing %s timepoint: %s\n", tp, e$message))
    cat("Stack trace:\n")
    traceback()
  })
}

cat("\n")
cat("============================================================\n")
cat("ALL TIMEPOINTS COMPLETE\n")
cat("============================================================\n")

for (tp in TIMEPOINTS_TO_RUN) {
  cat(sprintf("\n%s: %s\n", toupper(tp), OUTPUT_BASE[[tp]]))
}

cat("\nEnd time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
