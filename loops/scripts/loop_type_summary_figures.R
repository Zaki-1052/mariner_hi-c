# loops/scripts/loop_type_summary_figures.R
#
# Loop Type Summary Figures for Paper
#
# Produces pie charts of anchor type composition (8-category) and
# CRE vs structural classification (4-category) for gained vs lost loops.
# Reads pre-computed extended_characterized_loops.tsv — no re-annotation.
#
# Usage:
#   Rscript loops/scripts/loop_type_summary_figures.R --timepoint late
#   Rscript loops/scripts/loop_type_summary_figures.R --timepoint early
#   Rscript loops/scripts/loop_type_summary_figures.R --timepoint both
#
# Output:
#   loops/output/loop_type_summary/{timepoint}/

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(RColorBrewer)
})

# =============================================================================
# 1. CONFIGURATION
# =============================================================================

BASE_DIR <- getwd()

source(file.path(BASE_DIR, "loops/scripts/utils/multi_format_output.R"))

ANCHOR_TYPE_ORDER <- c("Active_Promoter", "Repressed_Promoter", "Bivalent_Promoter",
                       "Polycomb", "Active_Enhancer", "Poised_Enhancer",
                       "CTCF_Site", "Other")

ANCHOR_COLORS <- c(
  "Active_Promoter"    = "#e41a1c",
  "Repressed_Promoter" = "#756bb1",
  "Bivalent_Promoter"  = "#984ea3",
  "Polycomb"           = "#4daf4a",
  "Active_Enhancer"    = "#377eb8",
  "Poised_Enhancer"    = "#ff7f00",
  "CTCF_Site"          = "#a65628",
  "Other"              = "#999999"
)

CLASSIFICATION_ORDER <- c("structural", "CRE", "mixed", "unclassified")

CLASSIFICATION_COLORS <- c(
  "structural"   = "#7f7f7f",
  "CRE"          = "#2ca02c",
  "mixed"        = "#ff7f0e",
  "unclassified" = "#c7c7c7"
)

CLASSIFICATION_LABELS <- c(
  "structural"   = "Structural\n(CTCF-CTCF)",
  "CRE"          = "CRE\n(Enh/Prom)",
  "mixed"        = "Mixed\n(CTCF + CRE)",
  "unclassified" = "Unclassified"
)

EOR_P_TYPES <- c("Active_Promoter", "Repressed_Promoter", "Bivalent_Promoter",
                 "Active_Enhancer", "Poised_Enhancer")

INPUT_FILES <- list(
  late  = file.path(BASE_DIR, "peaks/loop_annotation_extended/late/extended_characterized_loops.tsv"),
  early = file.path(BASE_DIR, "peaks/loop_annotation_extended/early/extended_characterized_loops.tsv")
)

OUTPUT_BASE <- file.path(BASE_DIR, "loops/output/loop_type_summary")

PIE_LABEL_THRESHOLD <- 4.0

DIRECTION_LABELS <- c(
  "up_in_mutant"   = "Strengthened\n(Gained in BAP1-KO)",
  "down_in_mutant" = "Weakened\n(Lost in BAP1-KO)"
)

# =============================================================================
# 2. ARGUMENT PARSING
# =============================================================================

parse_args <- function(args) {
  timepoint <- "both"
  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--timepoint" && i < length(args)) {
      timepoint <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--help" || args[i] == "-h") {
      cat("Usage: Rscript loops/scripts/loop_type_summary_figures.R [--timepoint late|early|both]\n")
      quit(save = "no", status = 0)
    } else {
      i <- i + 1
    }
  }
  list(timepoint = timepoint)
}

# =============================================================================
# 3. CORE FUNCTIONS
# =============================================================================

load_loops <- function(timepoint) {
  input_file <- INPUT_FILES[[timepoint]]
  if (!file.exists(input_file)) stop(sprintf("Input file not found: %s", input_file))

  loops_df <- readr::read_tsv(input_file, show_col_types = FALSE)

  required_cols <- c("direction", "anchor1_type_extended", "anchor2_type_extended",
                     "anchor1_CTCF_motif_overlap", "anchor2_CTCF_motif_overlap")
  missing <- setdiff(required_cols, colnames(loops_df))
  if (length(missing) > 0) stop(sprintf("Missing columns: %s", paste(missing, collapse = ", ")))

  n_up   <- sum(loops_df$direction == "up_in_mutant")
  n_down <- sum(loops_df$direction == "down_in_mutant")
  cat(sprintf("  Loaded %d loops (%d up, %d down) for %s\n", nrow(loops_df), n_up, n_down, timepoint))
  if (nrow(loops_df) < 200) {
    cat(sprintf("  NOTE: Small sample size (n=%d); statistical tests may be unreliable\n", nrow(loops_df)))
  }
  loops_df
}

classify_cre_structural <- function(loops_df) {
  anchor1_eorp <- loops_df$anchor1_type_extended %in% EOR_P_TYPES
  anchor2_eorp <- loops_df$anchor2_type_extended %in% EOR_P_TYPES
  eorp_count   <- as.integer(anchor1_eorp) + as.integer(anchor2_eorp)

  ctcf_count <- as.integer(loops_df$anchor1_CTCF_motif_overlap) +
                as.integer(loops_df$anchor2_CTCF_motif_overlap)

  loops_df$classification <- factor(
    dplyr::case_when(
      ctcf_count == 2 & eorp_count <  2 ~ "structural",
      ctcf_count <  2 & eorp_count == 2 ~ "CRE",
      ctcf_count == 2 & eorp_count == 2 ~ "mixed",
      TRUE                               ~ "unclassified"
    ),
    levels = CLASSIFICATION_ORDER
  )
  loops_df
}

# =============================================================================
# 4. FIGURE 1: 8-CATEGORY ANCHOR TYPE PIE CHARTS
# =============================================================================

build_anchor_pie <- function(direction_data, direction_label, n_loops) {
  pie_data <- direction_data %>%
    dplyr::count(anchor_type, .drop = FALSE) %>%
    dplyr::mutate(
      percentage = 100 * n / sum(n),
      label = ifelse(percentage >= PIE_LABEL_THRESHOLD,
                     sprintf("%.1f%%", percentage), "")
    )

  ggplot(pie_data, aes(x = "", y = percentage, fill = anchor_type)) +
    geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.4) +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = ANCHOR_COLORS, name = "Anchor Type", drop = FALSE) +
    geom_text(aes(label = label),
              position = position_stack(vjust = 0.5),
              size = 3.2, fontface = "bold") +
    labs(title = direction_label,
         subtitle = sprintf("n = %d loops (%d anchors pooled)", n_loops, 2L * n_loops)) +
    theme_void(base_size = 12) +
    theme(
      plot.title    = element_text(hjust = 0.5, face = "bold", size = 13),
      plot.subtitle = element_text(hjust = 0.5, size = 9, color = "gray40"),
      legend.position = "right",
      legend.text     = element_text(size = 9),
      legend.key.size = unit(0.4, "cm")
    )
}

create_anchor_pie_charts <- function(loops_df, timepoint, output_dir) {
  anchor_pooled <- dplyr::bind_rows(
    loops_df %>% dplyr::transmute(direction, anchor_type = anchor1_type_extended),
    loops_df %>% dplyr::transmute(direction, anchor_type = anchor2_type_extended)
  ) %>%
    dplyr::filter(direction %in% c("up_in_mutant", "down_in_mutant")) %>%
    dplyr::mutate(anchor_type = factor(anchor_type, levels = ANCHOR_TYPE_ORDER))

  up_data   <- anchor_pooled %>% dplyr::filter(direction == "up_in_mutant")
  down_data <- anchor_pooled %>% dplyr::filter(direction == "down_in_mutant")
  n_up   <- sum(loops_df$direction == "up_in_mutant")
  n_down <- sum(loops_df$direction == "down_in_mutant")

  p_up   <- build_anchor_pie(up_data,   DIRECTION_LABELS[["up_in_mutant"]],   n_up)
  p_down <- build_anchor_pie(down_data, DIRECTION_LABELS[["down_in_mutant"]], n_down)

  p_combined <- (p_up | p_down) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = sprintf("Chromatin State Composition of Differential Loop Anchors (%s)",
                      toupper(timepoint)),
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15))
    )

  save_multiformat_ggplot(p_combined,
                          file.path(output_dir, "plots", "anchor_type_pie_charts"),
                          width = 12, height = 7)
  cat("  Saved: anchor_type_pie_charts/{pdf,svg,jpg}\n")

  invisible(list(up = p_up, down = p_down))
}

# =============================================================================
# 5. FIGURE 2: CRE VS STRUCTURAL CLASSIFICATION
# =============================================================================

build_classification_pie <- function(direction_data, direction_label, n_loops) {
  pie_data <- direction_data %>%
    dplyr::count(classification, .drop = FALSE) %>%
    dplyr::mutate(
      percentage = 100 * n / sum(n),
      label = ifelse(percentage >= PIE_LABEL_THRESHOLD,
                     sprintf("%.1f%%\n(%d)", percentage, n), "")
    )

  ggplot(pie_data, aes(x = "", y = percentage, fill = classification)) +
    geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.5) +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = CLASSIFICATION_COLORS, name = "Classification",
                      labels = CLASSIFICATION_LABELS, drop = FALSE) +
    geom_text(aes(label = label),
              position = position_stack(vjust = 0.5),
              size = 3.5, fontface = "bold") +
    labs(title = direction_label,
         subtitle = sprintf("n = %d loops", n_loops)) +
    theme_void(base_size = 12) +
    theme(
      plot.title    = element_text(hjust = 0.5, face = "bold", size = 13),
      plot.subtitle = element_text(hjust = 0.5, size = 9, color = "gray40"),
      legend.position = "right",
      legend.text     = element_text(size = 9),
      legend.key.size = unit(0.45, "cm")
    )
}

create_classification_figures <- function(loops_df, timepoint, output_dir) {
  diff_loops <- loops_df %>%
    dplyr::filter(direction %in% c("up_in_mutant", "down_in_mutant"))

  up_data   <- diff_loops %>% dplyr::filter(direction == "up_in_mutant")
  down_data <- diff_loops %>% dplyr::filter(direction == "down_in_mutant")
  n_up   <- nrow(up_data)
  n_down <- nrow(down_data)

  # --- Pie charts ---
  p_up   <- build_classification_pie(up_data,   DIRECTION_LABELS[["up_in_mutant"]],   n_up)
  p_down <- build_classification_pie(down_data, DIRECTION_LABELS[["down_in_mutant"]], n_down)

  contingency <- diff_loops %>%
    dplyr::count(direction, classification) %>%
    tidyr::pivot_wider(names_from = classification, values_from = n, values_fill = 0L) %>%
    tibble::column_to_rownames("direction") %>%
    as.matrix()

  chi2_result  <- NULL
  fisher_result <- NULL

  chi2_result <- tryCatch(
    chisq.test(contingency),
    warning = function(w) {
      msg <- conditionMessage(w)
      cat(sprintf("  WARNING: %s\n", msg))
      suppressWarnings(chisq.test(contingency))
    }
  )

  expected <- chi2_result$expected
  if (any(expected < 5)) {
    cat("  Running Fisher's exact test (expected cell count < 5)\n")
    fisher_result <- fisher.test(contingency, simulate.p.value = TRUE, B = 10000)
  }

  p_val <- if (!is.null(fisher_result)) fisher_result$p.value else chi2_result$p.value
  test_label <- if (!is.null(fisher_result)) {
    sprintf("Fisher p = %.2g", p_val)
  } else {
    sprintf("chi2 = %.1f, p = %.2g", chi2_result$statistic, p_val)
  }

  p_pie_combined <- (p_up | p_down) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = sprintf("CRE vs Structural Loop Classification (%s)", toupper(timepoint)),
      subtitle = test_label,
      theme = theme(
        plot.title    = element_text(hjust = 0.5, face = "bold", size = 15),
        plot.subtitle = element_text(hjust = 0.5, size = 11, face = "italic")
      )
    )

  save_multiformat_ggplot(p_pie_combined,
                          file.path(output_dir, "plots", "cre_structural_pie_charts"),
                          width = 12, height = 7)
  cat("  Saved: cre_structural_pie_charts/{pdf,svg,jpg}\n")

  # --- Stacked bar variant ---
  bar_summary <- diff_loops %>%
    dplyr::count(direction, classification, .drop = FALSE) %>%
    dplyr::group_by(direction) %>%
    dplyr::mutate(percentage = 100 * n / sum(n)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      direction = factor(direction,
                         levels = c("up_in_mutant", "down_in_mutant"),
                         labels = c("Strengthened\n(Up in Mutant)",
                                    "Weakened\n(Down in Mutant)"))
    )

  p_bar <- ggplot(bar_summary, aes(x = direction, y = percentage, fill = classification)) +
    geom_bar(stat = "identity", position = "stack", width = 0.6,
             color = "white", linewidth = 0.4) +
    scale_fill_manual(values = CLASSIFICATION_COLORS, name = "Classification",
                      labels = CLASSIFICATION_LABELS, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    labs(x = NULL, y = "Percentage of loops (%)",
         title = sprintf("CRE vs Structural Classification (%s)", toupper(timepoint)),
         subtitle = test_label) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title    = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 11, face = "italic"),
      axis.text.x   = element_text(size = 11),
      legend.position = "right",
      panel.grid.major.x = element_blank()
    )

  save_multiformat_ggplot(p_bar,
                          file.path(output_dir, "plots", "cre_structural_bar"),
                          width = 8, height = 6)
  cat("  Saved: cre_structural_bar/{pdf,svg,jpg}\n")

  invisible(list(up = p_up, down = p_down, chi2 = chi2_result, fisher = fisher_result))
}

# =============================================================================
# 6. FIGURE 3: COMBINED PANEL
# =============================================================================

create_combined_panel <- function(anchor_pies, class_pies, timepoint, output_dir) {
  p_combined <- (anchor_pies$up + anchor_pies$down) /
                (class_pies$up + class_pies$down) +
    plot_annotation(
      title = sprintf("Loop Anchor Classification Summary (%s)", toupper(timepoint)),
      subtitle = "Top: 8-category chromatin state | Bottom: CRE vs structural",
      tag_levels = "A",
      theme = theme(
        plot.title    = element_text(hjust = 0.5, face = "bold", size = 16),
        plot.subtitle = element_text(hjust = 0.5, size = 11, face = "italic")
      )
    )

  save_multiformat_ggplot(p_combined,
                          file.path(output_dir, "plots", "combined_summary_panel"),
                          width = 14, height = 13)
  cat("  Saved: combined_summary_panel/{pdf,svg,jpg}\n")
}

# =============================================================================
# 7. STATISTICS OUTPUT
# =============================================================================

write_summary_statistics <- function(loops_df, timepoint, output_dir,
                                     chi2_result, fisher_result) {
  n_total <- nrow(loops_df)
  n_up    <- sum(loops_df$direction == "up_in_mutant")
  n_down  <- sum(loops_df$direction == "down_in_mutant")

  lines <- c(
    "========================================",
    sprintf("Loop Type Summary Statistics (%s)", toupper(timepoint)),
    "========================================",
    "",
    "Script: loops/scripts/loop_type_summary_figures.R",
    sprintf("Input: %s", INPUT_FILES[[timepoint]]),
    sprintf("Output: %s", output_dir),
    sprintf("Generated: %s", Sys.Date()),
    "",
    sprintf("Total loops: %d", n_total),
    sprintf("  Up in mutant (strengthened):   %d (%.1f%%)", n_up, 100 * n_up / n_total),
    sprintf("  Down in mutant (weakened):     %d (%.1f%%)", n_down, 100 * n_down / n_total)
  )

  if (n_total < 200) {
    lines <- c(lines, "",
      sprintf("NOTE: Small sample size (n=%d). Statistical tests may be unreliable.", n_total))
  }

  # --- Figure 1: 8-category ---
  lines <- c(lines, "",
    "--- FIGURE 1: 8-Category Anchor Type Composition ---", "")

  anchor_pooled <- dplyr::bind_rows(
    loops_df %>% dplyr::transmute(direction, anchor_type = anchor1_type_extended),
    loops_df %>% dplyr::transmute(direction, anchor_type = anchor2_type_extended)
  ) %>%
    dplyr::filter(direction %in% c("up_in_mutant", "down_in_mutant"))

  for (dir in c("up_in_mutant", "down_in_mutant")) {
    dir_data <- anchor_pooled %>% dplyr::filter(direction == dir)
    n_dir <- nrow(dir_data)
    dir_label <- if (dir == "up_in_mutant") "Strengthened (up_in_mutant)" else "Weakened (down_in_mutant)"
    n_loops_dir <- if (dir == "up_in_mutant") n_up else n_down
    lines <- c(lines, sprintf("  %s, n_loops=%d, n_anchors=%d:", dir_label, n_loops_dir, n_dir))
    type_counts <- dir_data %>%
      dplyr::count(anchor_type) %>%
      dplyr::arrange(match(anchor_type, ANCHOR_TYPE_ORDER))
    for (i in seq_len(nrow(type_counts))) {
      lines <- c(lines, sprintf("    %-22s %5d (%.1f%%)",
                                type_counts$anchor_type[i],
                                type_counts$n[i],
                                100 * type_counts$n[i] / n_dir))
    }
    lines <- c(lines, "")
  }

  # --- Figure 2: CRE/structural ---
  lines <- c(lines,
    "--- FIGURE 2: CRE vs Structural Classification ---", "",
    "Classification logic:",
    "  CTCF: anchor1_CTCF_motif_overlap + anchor2_CTCF_motif_overlap",
    "  EorP: anchor type in {Active_Promoter, Repressed_Promoter, Bivalent_Promoter,",
    "                        Active_Enhancer, Poised_Enhancer}",
    "  structural:   CTCF==2 AND EorP<2",
    "  CRE:          CTCF<2  AND EorP==2",
    "  mixed:        CTCF==2 AND EorP==2",
    "  unclassified: all other combinations", "")

  diff_loops <- loops_df %>%
    dplyr::filter(direction %in% c("up_in_mutant", "down_in_mutant"))

  for (dir in c("up_in_mutant", "down_in_mutant")) {
    dir_data <- diff_loops %>% dplyr::filter(direction == dir)
    n_dir <- nrow(dir_data)
    dir_label <- if (dir == "up_in_mutant") "Strengthened (up_in_mutant)" else "Weakened (down_in_mutant)"
    lines <- c(lines, sprintf("  %s (n=%d):", dir_label, n_dir))
    class_counts <- dir_data %>%
      dplyr::count(classification) %>%
      dplyr::arrange(match(classification, CLASSIFICATION_ORDER))
    for (i in seq_len(nrow(class_counts))) {
      lines <- c(lines, sprintf("    %-16s %5d (%.1f%%)",
                                class_counts$classification[i],
                                class_counts$n[i],
                                100 * class_counts$n[i] / n_dir))
    }
    lines <- c(lines, "")
  }

  # Contingency table
  contingency <- diff_loops %>%
    dplyr::count(direction, classification) %>%
    tidyr::pivot_wider(names_from = classification, values_from = n, values_fill = 0L)

  lines <- c(lines, "Contingency table (counts):")
  lines <- c(lines, sprintf("  %-18s %s",
                             "",
                             paste(sprintf("%14s", CLASSIFICATION_ORDER), collapse = "")))
  for (i in seq_len(nrow(contingency))) {
    dir_short <- if (contingency$direction[i] == "up_in_mutant") "up" else "down"
    vals <- sapply(CLASSIFICATION_ORDER, function(cl) {
      v <- contingency[[cl]][i]
      if (is.null(v) || is.na(v)) 0L else v
    })
    lines <- c(lines, sprintf("  %-18s %s", paste0(dir_short, ":"),
                               paste(sprintf("%14d", vals), collapse = "")))
  }
  lines <- c(lines, "")

  if (!is.null(chi2_result)) {
    lines <- c(lines, sprintf("Chi-squared test: chi2 = %.2f, df = %d, p = %.4g",
                               chi2_result$statistic, chi2_result$parameter,
                               chi2_result$p.value))
  }
  if (!is.null(fisher_result)) {
    lines <- c(lines, sprintf("Fisher's exact test (simulated): p = %.4g",
                               fisher_result$p.value))
  }

  lines <- c(lines, "", "========================================")

  stats_file <- file.path(output_dir, "summary_statistics.txt")
  writeLines(lines, stats_file)
  cat(sprintf("  Statistics written to: %s\n", stats_file))

  # --- TSV outputs ---
  anchor_type_tsv <- anchor_pooled %>%
    dplyr::count(direction, anchor_type) %>%
    dplyr::group_by(direction) %>%
    dplyr::mutate(percentage = round(100 * n / sum(n), 2)) %>%
    dplyr::ungroup() %>%
    dplyr::rename(count = n)

  readr::write_tsv(anchor_type_tsv,
                   file.path(output_dir, "anchor_type_by_direction.tsv"))

  class_tsv <- diff_loops %>%
    dplyr::count(direction, classification) %>%
    dplyr::group_by(direction) %>%
    dplyr::mutate(percentage = round(100 * n / sum(n), 2)) %>%
    dplyr::ungroup() %>%
    dplyr::rename(count = n)

  readr::write_tsv(class_tsv,
                   file.path(output_dir, "classification_by_direction.tsv"))

  cat("  Tables written to output directory\n")
}

# =============================================================================
# 8. MAIN ORCHESTRATOR
# =============================================================================

run_timepoint <- function(timepoint) {
  cat(sprintf("\n%s\n", paste(rep("=", 50), collapse = "")))
  cat(sprintf("Loop Type Summary Figures (%s)\n", toupper(timepoint)))
  cat(sprintf("%s\n\n", paste(rep("=", 50), collapse = "")))

  output_dir <- file.path(OUTPUT_BASE, timepoint)
  dir.create(file.path(output_dir, "plots"), recursive = TRUE, showWarnings = FALSE)

  cat("Step 1: Loading annotated loops...\n")
  loops_df <- load_loops(timepoint)

  cat("\nStep 2: Classifying loops (CRE/structural)...\n")
  loops_df <- classify_cre_structural(loops_df)
  class_tab <- table(loops_df$classification[loops_df$direction %in% c("up_in_mutant", "down_in_mutant")])
  cat(sprintf("  %s\n", paste(sprintf("%s=%d", names(class_tab), class_tab), collapse = " | ")))

  cat("\nStep 3: Creating anchor type pie charts...\n")
  anchor_pies <- create_anchor_pie_charts(loops_df, timepoint, output_dir)

  cat("\nStep 4: Creating CRE vs structural figures...\n")
  class_result <- create_classification_figures(loops_df, timepoint, output_dir)

  cat("\nStep 5: Creating combined panel...\n")
  create_combined_panel(anchor_pies, class_result, timepoint, output_dir)

  cat("\nStep 6: Writing statistics...\n")
  write_summary_statistics(loops_df, timepoint, output_dir,
                           class_result$chi2, class_result$fisher)

  cat(sprintf("\nCompleted: %s\n", timepoint))
  cat(sprintf("Output: %s\n", output_dir))
}

# =============================================================================
# 9. DISPATCH
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
parsed <- parse_args(args)

if (parsed$timepoint == "both") {
  timepoints <- c("late", "early")
} else {
  timepoints <- parsed$timepoint
}

for (tp in timepoints) {
  run_timepoint(tp)
}

cat("\nDone.\n")
