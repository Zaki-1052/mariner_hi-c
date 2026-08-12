# scripts/se_stripe_overlap.R
#
# 1b. Differential Stripes x Superenhancers
#
# Overlaps differential stripes (stripenn) with superenhancer BEDs.
# Tests: are gained stripes enriched for SE overlap? Does K27ac change
# at SE-overlapping stripe anchors match stripe direction?
#
# Usage:
#   Rscript scripts/se_stripe_overlap.R --timepoint late
#   Rscript scripts/se_stripe_overlap.R --timepoint early
#   Rscript scripts/se_stripe_overlap.R --timepoint both
#   Rscript scripts/se_stripe_overlap.R --timepoint late --stripe-set highconf
#
# Output:
#   output/superenhancer_analysis/1b_se_stripe_overlap/{timepoint}/{stripe_set}/

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
})

# =============================================================================
# 1. CONFIGURATION
# =============================================================================

BASE_DIR <- getwd()

source(file.path(BASE_DIR, "loops/scripts/utils/multi_format_output.R"))
source(file.path(BASE_DIR, "loops/scripts/utils/se_utils.R"))

SE_FILES <- list(
  P60    = file.path(BASE_DIR, "peaks/Superenhancers_P60.bed"),
  ENCODE = file.path(BASE_DIR, "peaks/Superenhancers_encode.bed")
)

DIFFBIND_K27AC <- file.path(BASE_DIR, "peaks/diffbind/K27ac_diffbind_results_summit_appended_ap.txt")

STRIPE_FILES <- list(
  late = list(
    allsig   = file.path(BASE_DIR, "stripes/stripenn/outputs/250402/250402_stripes_allsig.bedpe"),
    highconf = file.path(BASE_DIR, "stripes/stripenn/outputs/250402/250402_stripes_highconf.bedpe")
  ),
  early = list(
    allsig     = file.path(BASE_DIR, "stripes/stripenn/outputs/250831/250831_stripes_allsig.bedpe"),
    concordant = file.path(BASE_DIR, "stripes/stripenn/outputs/250831/250831_stripes_concordant.bedpe")
  )
)

DEFAULT_STRIPE_SET <- list(late = "allsig", early = "concordant")

# =============================================================================
# 2. ARGUMENT PARSING
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  timepoint <- "late"
  stripe_set <- NULL

  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--timepoint" && i < length(args)) {
      timepoint <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--stripe-set" && i < length(args)) {
      stripe_set <- args[i + 1]
      i <- i + 2
    } else {
      i <- i + 1
    }
  }

  return(list(timepoint = timepoint, stripe_set = stripe_set))
}

parsed <- parse_args(args)

# =============================================================================
# 3. ANALYSIS FUNCTIONS
# =============================================================================

load_stripes <- function(stripe_file) {
  if (!file.exists(stripe_file)) stop(sprintf("Stripe file not found: %s", stripe_file))

  stripes <- read_tsv(stripe_file, show_col_types = FALSE)

  required_cols <- c("chr1", "x1", "x2", "direction", "logFC", "FDR", "h3k27ac")
  missing <- setdiff(required_cols, colnames(stripes))
  if (length(missing) > 0) {
    stop(sprintf("Missing stripe columns: %s", paste(missing, collapse = ", ")))
  }

  stripes <- stripes %>%
    dplyr::filter(chr1 %in% VALID_CHROMS) %>%
    dplyr::mutate(stripe_idx = row_number())

  cat(sprintf("  Loaded %d stripes from %s\n", nrow(stripes), basename(stripe_file)))
  cat(sprintf("    - gained: %d, lost: %d\n",
              sum(stripes$direction == "gained"),
              sum(stripes$direction == "lost")))

  return(stripes)
}

overlap_stripes_with_ses <- function(stripes, se_gr, se_label) {
  anchor_gr <- GRanges(
    seqnames = stripes$chr1,
    ranges = IRanges(start = stripes$x1, end = stripes$x2)
  )

  hits <- findOverlaps(anchor_gr, se_gr, ignore.strand = TRUE)

  overlap_col <- paste0("se_", se_label, "_overlap")
  stripes[[overlap_col]] <- FALSE
  stripes[[overlap_col]][unique(queryHits(hits))] <- TRUE

  n_overlap <- sum(stripes[[overlap_col]])
  cat(sprintf("  %s SE overlap: %d / %d stripes (%.1f%%)\n",
              se_label, n_overlap, nrow(stripes),
              100 * n_overlap / nrow(stripes)))

  return(stripes)
}

run_fisher_test <- function(stripes, se_label) {
  overlap_col <- paste0("se_", se_label, "_overlap")

  diff_stripes <- stripes %>%
    dplyr::filter(direction %in% c("gained", "lost"))

  tbl <- table(
    direction = diff_stripes$direction,
    se_overlap = diff_stripes[[overlap_col]]
  )

  if (ncol(tbl) < 2 || nrow(tbl) < 2) {
    cat(sprintf("  Fisher test skipped for %s — insufficient cells\n", se_label))
    return(NULL)
  }

  ft <- fisher.test(tbl)

  result <- tibble(
    se_set = se_label,
    gained_with_se = tbl["gained", "TRUE"],
    gained_without_se = tbl["gained", "FALSE"],
    lost_with_se = tbl["lost", "TRUE"],
    lost_without_se = tbl["lost", "FALSE"],
    odds_ratio = ft$estimate,
    ci_lower = ft$conf.int[1],
    ci_upper = ft$conf.int[2],
    p_value = ft$p.value
  )

  cat(sprintf("  Fisher test (%s): OR = %.2f [%.2f-%.2f], p = %.2e\n",
              se_label, ft$estimate, ft$conf.int[1], ft$conf.int[2], ft$p.value))

  return(result)
}

classify_stripe_k27ac <- function(stripes, se_gr, diffbind_gr, se_label) {
  overlap_col <- paste0("se_", se_label, "_overlap")
  se_stripes <- stripes %>% dplyr::filter(.data[[overlap_col]] == TRUE)

  if (nrow(se_stripes) == 0) {
    cat(sprintf("  No SE-overlapping stripes for %s — skipping K27ac classification\n", se_label))
    return(NULL)
  }

  anchor_gr <- GRanges(
    seqnames = se_stripes$chr1,
    ranges = IRanges(start = se_stripes$x1, end = se_stripes$x2)
  )

  k27ac_result <- classify_k27ac_change(anchor_gr, diffbind_gr)
  se_stripes$k27ac_class <- k27ac_result$k27ac_class

  concordance <- se_stripes %>%
    dplyr::filter(direction %in% c("gained", "lost")) %>%
    dplyr::mutate(
      concordant = (direction == "gained" & k27ac_class == "gained_k27ac") |
                   (direction == "lost"   & k27ac_class == "lost_k27ac")
    )

  cat(sprintf("  K27ac concordance (%s): %d / %d SE-stripes concordant (%.1f%%)\n",
              se_label, sum(concordance$concordant), nrow(concordance),
              100 * sum(concordance$concordant) / max(1, nrow(concordance))))

  return(concordance)
}

# =============================================================================
# 4. PLOTTING FUNCTIONS
# =============================================================================

plot_se_overlap_bar <- function(stripes, se_label, output_dir) {
  overlap_col <- paste0("se_", se_label, "_overlap")

  plot_df <- stripes %>%
    dplyr::filter(direction %in% c("gained", "lost")) %>%
    dplyr::group_by(direction) %>%
    dplyr::summarise(
      with_se = sum(.data[[overlap_col]]),
      without_se = sum(!.data[[overlap_col]]),
      total = n(),
      pct_se = 100 * with_se / total,
      .groups = "drop"
    ) %>%
    tidyr::pivot_longer(cols = c(with_se, without_se),
                        names_to = "overlap_status", values_to = "count") %>%
    dplyr::mutate(overlap_status = factor(overlap_status,
                                          levels = c("with_se", "without_se"),
                                          labels = c("SE overlap", "No SE overlap")))

  p <- ggplot(plot_df, aes(x = direction, y = count, fill = overlap_status)) +
    geom_bar(stat = "identity", position = "stack", width = 0.7) +
    scale_fill_manual(values = c("SE overlap" = "#e0a730", "No SE overlap" = "#cccccc")) +
    scale_x_discrete(labels = c("gained" = "Gained", "lost" = "Lost")) +
    labs(
      title = sprintf("Differential Stripes x %s Superenhancers", se_label),
      x = "Stripe Direction", y = "Number of Stripes", fill = ""
    ) +
    theme_classic(base_size = 14) +
    theme(legend.position = "top")

  save_multiformat_ggplot(p, file.path(output_dir, sprintf("se_stripe_bar_%s", tolower(se_label))),
                          width = 5, height = 6)
  return(p)
}

plot_logfc_by_se <- function(stripes, se_label, output_dir) {
  overlap_col <- paste0("se_", se_label, "_overlap")

  plot_df <- stripes %>%
    dplyr::filter(direction %in% c("gained", "lost")) %>%
    dplyr::mutate(se_status = ifelse(.data[[overlap_col]], "SE overlap", "No SE overlap"))

  if (length(unique(plot_df$se_status)) < 2) {
    cat(sprintf("  Skipping logFC violin for %s — only one overlap class\n", se_label))
    return(NULL)
  }

  p <- ggplot(plot_df, aes(x = se_status, y = logFC, fill = se_status)) +
    geom_violin(trim = TRUE, alpha = 0.7) +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.5) +
    scale_fill_manual(values = c("SE overlap" = "#e0a730", "No SE overlap" = "#cccccc")) +
    stat_compare_means(method = "wilcox.test", label = "p.format",
                       label.x.npc = "center", label.y.npc = "top") +
    labs(
      title = sprintf("Stripe Effect Size by %s SE Overlap", se_label),
      x = "", y = "logFC (gained > 0, lost < 0)", fill = ""
    ) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none")

  save_multiformat_ggplot(p, file.path(output_dir, sprintf("se_stripe_logfc_violin_%s", tolower(se_label))),
                          width = 5, height = 6)
  return(p)
}

plot_k27ac_concordance <- function(concordance_df, se_label, output_dir) {
  if (is.null(concordance_df) || nrow(concordance_df) == 0) return(NULL)

  plot_df <- concordance_df %>%
    dplyr::filter(k27ac_class != "no_diffbind_peak") %>%
    dplyr::count(direction, k27ac_class)

  if (nrow(plot_df) == 0) return(NULL)

  p <- ggplot(plot_df, aes(x = direction, y = n, fill = k27ac_class)) +
    geom_bar(stat = "identity", position = "fill", width = 0.7) +
    scale_fill_manual(values = K27AC_CLASS_COLORS,
                      labels = c("gained_k27ac" = "Gained K27ac",
                                 "lost_k27ac" = "Lost K27ac",
                                 "stable_k27ac" = "Stable K27ac")) +
    scale_x_discrete(labels = c("gained" = "Gained stripes", "lost" = "Lost stripes")) +
    scale_y_continuous(labels = scales::percent) +
    labs(
      title = sprintf("K27ac Change at %s SE-Overlapping Stripe Anchors", se_label),
      x = "", y = "Proportion of SE-Overlapping Stripes", fill = "K27ac Status"
    ) +
    theme_classic(base_size = 14) +
    theme(legend.position = "right")

  save_multiformat_ggplot(p, file.path(output_dir, sprintf("k27ac_concordance_%s", tolower(se_label))),
                          width = 7, height = 6)
  return(p)
}

plot_p60_vs_encode <- function(stripes, output_dir) {
  diff_stripes <- stripes %>%
    dplyr::filter(direction %in% c("gained", "lost"))

  if (!"se_P60_overlap" %in% colnames(diff_stripes) ||
      !"se_ENCODE_overlap" %in% colnames(diff_stripes)) {
    return(NULL)
  }

  summary_df <- diff_stripes %>%
    dplyr::mutate(
      overlap_class = dplyr::case_when(
        se_P60_overlap & se_ENCODE_overlap ~ "Both",
        se_P60_overlap & !se_ENCODE_overlap ~ "P60 only",
        !se_P60_overlap & se_ENCODE_overlap ~ "ENCODE only",
        TRUE ~ "Neither"
      )
    ) %>%
    dplyr::count(direction, overlap_class) %>%
    dplyr::mutate(overlap_class = factor(overlap_class,
                                         levels = c("Both", "P60 only", "ENCODE only", "Neither")))

  p <- ggplot(summary_df, aes(x = direction, y = n, fill = overlap_class)) +
    geom_bar(stat = "identity", position = "stack", width = 0.7) +
    scale_fill_manual(values = c("Both" = "#d73027", "P60 only" = "#e0a730",
                                 "ENCODE only" = "#4575b4", "Neither" = "#cccccc")) +
    scale_x_discrete(labels = c("gained" = "Gained", "lost" = "Lost")) +
    labs(
      title = "SE Set Comparison: P60 vs ENCODE",
      x = "Stripe Direction", y = "Number of Stripes", fill = "SE Overlap"
    ) +
    theme_classic(base_size = 14) +
    theme(legend.position = "right")

  save_multiformat_ggplot(p, file.path(output_dir, "se_p60_vs_encode_comparison"),
                          width = 7, height = 6)
  return(p)
}

# =============================================================================
# 5. MAIN ANALYSIS
# =============================================================================

run_analysis <- function(timepoint, stripe_set) {
  cat(sprintf("\n========== 1b. SE-Stripe Overlap: %s (%s) ==========\n\n",
              toupper(timepoint), stripe_set))

  stripe_file <- STRIPE_FILES[[timepoint]][[stripe_set]]
  if (is.null(stripe_file) || !file.exists(stripe_file)) {
    stop(sprintf("Stripe file not available for %s / %s", timepoint, stripe_set))
  }

  output_dir <- file.path(BASE_DIR, "loops/output/superenhancer_analysis/1b_se_stripe_overlap",
                           timepoint, stripe_set)
  tables_dir <- file.path(output_dir, "tables")
  stats_dir <- file.path(output_dir, "statistics")
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(stats_dir, recursive = TRUE, showWarnings = FALSE)

  # --- Load data ---
  cat("Loading stripes...\n")
  stripes <- load_stripes(stripe_file)

  cat("\nLoading superenhancers...\n")
  p60_gr <- load_se_bed(SE_FILES$P60, "P60")
  encode_gr <- load_se_bed(SE_FILES$ENCODE, "ENCODE")

  cat("\nLoading K27ac DiffBind...\n")
  diffbind_gr <- load_k27ac_diffbind(DIFFBIND_K27AC)

  # --- Overlap analysis ---
  cat("\n--- Overlap Analysis ---\n")
  stripes <- overlap_stripes_with_ses(stripes, p60_gr, "P60")
  stripes <- overlap_stripes_with_ses(stripes, encode_gr, "ENCODE")

  # --- Fisher tests ---
  cat("\n--- Contingency Tests ---\n")
  fisher_p60 <- run_fisher_test(stripes, "P60")
  fisher_encode <- run_fisher_test(stripes, "ENCODE")
  fisher_results <- bind_rows(fisher_p60, fisher_encode)

  # --- K27ac concordance ---
  cat("\n--- K27ac Concordance ---\n")
  concordance_p60 <- classify_stripe_k27ac(stripes, p60_gr, diffbind_gr, "P60")
  concordance_encode <- classify_stripe_k27ac(stripes, encode_gr, diffbind_gr, "ENCODE")

  # --- Save tables ---
  cat("\n--- Saving Tables ---\n")
  write_tsv(stripes, file.path(tables_dir, "stripe_se_overlap.tsv"))
  cat(sprintf("  Wrote stripe_se_overlap.tsv (%d rows)\n", nrow(stripes)))

  if (nrow(fisher_results) > 0) {
    write_tsv(fisher_results, file.path(tables_dir, "fisher_test_results.tsv"))
    cat(sprintf("  Wrote fisher_test_results.tsv\n"))
  }

  if (!is.null(concordance_p60)) {
    write_tsv(concordance_p60, file.path(tables_dir, "k27ac_concordance_P60.tsv"))
  }
  if (!is.null(concordance_encode)) {
    write_tsv(concordance_encode, file.path(tables_dir, "k27ac_concordance_ENCODE.tsv"))
  }

  # --- Summary statistics ---
  stats_file <- file.path(stats_dir, "se_stripe_statistics.txt")
  sink(stats_file)
  cat(sprintf("SE-Stripe Overlap Analysis\n"))
  cat(sprintf("Timepoint: %s\n", timepoint))
  cat(sprintf("Stripe set: %s\n", stripe_set))
  cat(sprintf("Date: %s\n\n", Sys.time()))

  cat(sprintf("Total stripes: %d\n", nrow(stripes)))
  cat(sprintf("  gained: %d\n", sum(stripes$direction == "gained")))
  cat(sprintf("  lost: %d\n\n", sum(stripes$direction == "lost")))

  cat(sprintf("P60 SE overlaps: %d / %d (%.1f%%)\n",
              sum(stripes$se_P60_overlap), nrow(stripes),
              100 * sum(stripes$se_P60_overlap) / nrow(stripes)))
  cat(sprintf("ENCODE SE overlaps: %d / %d (%.1f%%)\n\n",
              sum(stripes$se_ENCODE_overlap), nrow(stripes),
              100 * sum(stripes$se_ENCODE_overlap) / nrow(stripes)))

  if (nrow(fisher_results) > 0) {
    cat("Fisher Exact Tests (direction x SE overlap):\n")
    for (i in seq_len(nrow(fisher_results))) {
      r <- fisher_results[i, ]
      cat(sprintf("  %s: OR = %.3f [%.3f - %.3f], p = %.2e\n",
                  r$se_set, r$odds_ratio, r$ci_lower, r$ci_upper, r$p_value))
    }
    cat("\n")
  }

  if (!is.null(concordance_p60)) {
    cat("K27ac concordance (P60):\n")
    print(table(concordance_p60$direction, concordance_p60$k27ac_class))
    cat("\n")
  }
  if (!is.null(concordance_encode)) {
    cat("K27ac concordance (ENCODE):\n")
    print(table(concordance_encode$direction, concordance_encode$k27ac_class))
    cat("\n")
  }
  sink()
  cat(sprintf("  Wrote %s\n", stats_file))

  # --- Plots ---
  cat("\n--- Generating Plots ---\n")
  plot_se_overlap_bar(stripes, "P60", output_dir)
  plot_se_overlap_bar(stripes, "ENCODE", output_dir)
  plot_logfc_by_se(stripes, "P60", output_dir)
  plot_logfc_by_se(stripes, "ENCODE", output_dir)
  plot_k27ac_concordance(concordance_p60, "P60", output_dir)
  plot_k27ac_concordance(concordance_encode, "ENCODE", output_dir)
  plot_p60_vs_encode(stripes, output_dir)

  cat(sprintf("\n========== Done: %s / %s ==========\n", timepoint, stripe_set))
  return(invisible(stripes))
}

# =============================================================================
# 6. DISPATCH
# =============================================================================

if (parsed$timepoint == "both") {
  timepoints <- c("late", "early")
} else {
  timepoints <- parsed$timepoint
}

for (tp in timepoints) {
  available_sets <- names(STRIPE_FILES[[tp]])
  available_sets <- available_sets[!sapply(STRIPE_FILES[[tp]], is.null)]

  if (!is.null(parsed$stripe_set)) {
    if (!parsed$stripe_set %in% available_sets) {
      stop(sprintf("Stripe set '%s' not available for %s. Available: %s",
                   parsed$stripe_set, tp, paste(available_sets, collapse = ", ")))
    }
    sets_to_run <- parsed$stripe_set
  } else {
    sets_to_run <- available_sets
  }

  for (ss in sets_to_run) {
    run_analysis(tp, ss)
  }
}
