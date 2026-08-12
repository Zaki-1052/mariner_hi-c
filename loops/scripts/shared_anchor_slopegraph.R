# scripts/shared_anchor_slopegraph.R
#
# Slopegraph (paired dot-and-line plot) showing that at 212 shared anchor
# hubs, the lost loop is systematically longer than the gained loop.
# Each line connects one anchor's median lost distance to its median
# gained distance.  Poster figure for Beat 1 (distance shift).
#
# Reads pre-computed paired distance table; does not recompute from raw loops.

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(scales)
})

source("scripts/utils/multi_format_output.R")

# ---- Configuration -----------------------------------------------------------

BASE      <- getwd()
DATA_TSV  <- file.path(BASE, "output/shared_anchor_analysis/late/tables/paired_distance_stats.tsv")
LOOPS_TSV <- file.path(BASE, "output/shared_anchor_analysis/late/tables/shared_anchor_loops.tsv")
OUT_DIR   <- file.path(BASE, "output/shared_anchor_analysis/late/plots")

COL_LOST   <- "#d73027"
COL_GAINED <- "#4575b4"

# ---- Load & classify ---------------------------------------------------------

df <- read_tsv(DATA_TSV, show_col_types = FALSE)
loops_df <- read_tsv(LOOPS_TSV, show_col_types = FALSE)

n_lost_loops   <- sum(loops_df$direction == "down_in_mutant")
n_gained_loops <- sum(loops_df$direction == "up_in_mutant")

cat(sprintf("Loaded %d shared anchors\n", nrow(df)))
cat(sprintf("  Unique loops — lost: %d, gained: %d\n", n_lost_loops, n_gained_loops))

df <- df %>%
  mutate(
    direction = ifelse(median_lost_distance > median_gained_distance,
                       "lost_longer", "gained_longer")
  ) %>%
  arrange(direction == "lost_longer")

n_total    <- nrow(df)
n_lost_gt  <- sum(df$direction == "lost_longer")
pct_lost   <- round(n_lost_gt / n_total * 100)

med_lost   <- median(df$median_lost_distance)
med_gained <- median(df$median_gained_distance)

wilcox_p <- wilcox.test(
  df$median_lost_distance,
  df$median_gained_distance,
  paired = TRUE, alternative = "greater"
)$p.value

cat(sprintf("  %d/%d (%d%%) lost > gained\n", n_lost_gt, n_total, pct_lost))
cat(sprintf("  Median lost: %.0f kb, gained: %.0f kb\n", med_lost / 1e3, med_gained / 1e3))
cat(sprintf("  Paired Wilcoxon p = %.2e\n", wilcox_p))

# ---- Median markers ----------------------------------------------------------

median_df <- tibble(
  x     = c(1, 2),
  y     = c(med_lost, med_gained),
  label = c(sprintf("%.0f kb", med_lost / 1e3),
            sprintf("%.0f kb", med_gained / 1e3)),
  hjust = c(1.15, -0.15)
)

# ---- Build slopegraph --------------------------------------------------------

build_slopegraph <- function(df, median_df) {

  ggplot(df) +
    geom_segment(
      aes(x = 1, xend = 2,
          y = median_lost_distance, yend = median_gained_distance,
          color = direction, alpha = direction),
      linewidth = 0.4
    ) +
    geom_point(
      aes(x = 1, y = median_lost_distance, color = direction),
      size = 1.0, alpha = 0.6
    ) +
    geom_point(
      aes(x = 2, y = median_gained_distance, color = direction),
      size = 1.0, alpha = 0.6
    ) +
    geom_point(
      data = median_df, aes(x = x, y = y),
      shape = 18, size = 5, color = "black"
    ) +
    geom_text(
      data = median_df, aes(x = x, y = y, label = label, hjust = hjust),
      size = 3.8, fontface = "bold", vjust = 0.5
    ) +
    scale_color_manual(
      values = c("lost_longer" = COL_LOST, "gained_longer" = COL_GAINED),
      guide  = "none"
    ) +
    scale_alpha_manual(
      values = c("lost_longer" = 0.35, "gained_longer" = 0.65),
      guide  = "none"
    ) +
    scale_y_log10(
      breaks = c(1e4, 2e4, 5e4, 1e5, 2e5, 5e5, 1e6, 2e6, 5e6, 1e7, 2e7, 5e7),
      minor_breaks = NULL,
      labels = label_number(scale = 1e-3, suffix = " kb", big.mark = ","),
      expand = expansion(mult = c(0.05, 0.12))
    ) +
    scale_x_continuous(
      breaks = c(1, 2),
      labels = c(sprintf("Lost\n(n = %d)", n_lost_loops),
                 sprintf("Gained\n(n = %d)", n_gained_loops)),
      limits = c(0.55, 2.45),
      expand = expansion(0)
    ) +
    annotate(
      "text", x = 1.5,
      y = max(df$median_lost_distance, df$median_gained_distance) * 2.0,
      label = sprintf("%d%% of anchors: lost loop is longer\nPaired Wilcoxon p = %.2e",
                      pct_lost, wilcox_p),
      size = 3.5, hjust = 0.5, vjust = 1,
      fontface = "italic", color = "gray30", lineheight = 0.95
    ) +
    labs(
      title    = "Loop distance at shared anchors",
      subtitle = sprintf("Each line = one anchor hub (n = %d)", n_total),
      x = NULL,
      y = "Median loop distance"
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title    = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, color = "gray40", size = 11),
      axis.text.x   = element_text(size = 14, face = "bold"),
      axis.text.y   = element_text(size = 11),
      axis.title.y  = element_text(size = 12),
      panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
      axis.ticks.x  = element_blank(),
      axis.line.x   = element_blank()
    )
}

# ---- Poster variant ----------------------------------------------------------

build_slopegraph_poster <- function(df, median_df) {

  ggplot(df) +
    geom_segment(
      aes(x = 1, xend = 2,
          y = median_lost_distance, yend = median_gained_distance,
          color = direction, alpha = direction),
      linewidth = 0.4
    ) +
    geom_point(
      aes(x = 1, y = median_lost_distance, color = direction),
      size = 1.0, alpha = 0.3
    ) +
    geom_point(
      aes(x = 2, y = median_gained_distance, color = direction),
      size = 1.0, alpha = 0.3
    ) +
    geom_segment(
      data = tibble(x = 1, xend = 2, y = med_lost, yend = med_gained),
      aes(x = x, xend = xend, y = y, yend = yend),
      linewidth = 2.5, color = "black", lineend = "round"
    ) +
    geom_point(
      data = median_df, aes(x = x, y = y),
      shape = 18, size = 8, color = "black"
    ) +
    geom_text(
      data = median_df, aes(x = x, y = y, label = label, hjust = hjust),
      size = 6.0, fontface = "bold", vjust = 0.5
    ) +
    scale_color_manual(
      values = c("lost_longer" = COL_LOST, "gained_longer" = COL_GAINED),
      guide  = "none"
    ) +
    scale_alpha_manual(
      values = c("lost_longer" = 0.12, "gained_longer" = 0.20),
      guide  = "none"
    ) +
    scale_y_log10(
      breaks = c(1e4, 2e4, 5e4, 1e5, 2e5, 5e5, 1e6, 2e6, 5e6, 1e7, 2e7, 5e7),
      minor_breaks = NULL,
      labels = label_number(scale = 1e-3, suffix = " kb", big.mark = ","),
      expand = expansion(mult = c(0.05, 0.05))
    ) +
    scale_x_continuous(
      breaks = c(1, 2),
      labels = c(sprintf("Lost\n(n = %d)", n_lost_loops),
                 sprintf("Gained\n(n = %d)", n_gained_loops)),
      limits = c(0.3, 2.7),
      expand = expansion(0)
    ) +
    annotate(
      "text", x = 1.5,
      y = max(df$median_lost_distance, df$median_gained_distance) * 2.0,
      label = sprintf("%d%% of anchors: lost loop is longer\nPaired Wilcoxon p = %.2e",
                      pct_lost, wilcox_p),
      size = 5.5, hjust = 0.5, vjust = 1,
      fontface = "italic", color = "gray30", lineheight = 0.95
    ) +
    labs(
      title    = "Loop distance at shared anchors",
      subtitle = sprintf("Each line = one anchor hub (n = %d)", n_total),
      x = NULL,
      y = "Median loop distance"
    ) +
    theme_classic(base_size = 20) +
    theme(
      plot.title    = element_text(hjust = 0.5, face = "bold", size = 22),
      plot.subtitle = element_text(hjust = 0.5, color = "gray40", size = 14),
      axis.text.x   = element_text(size = 18, face = "bold"),
      axis.text.y   = element_text(size = 14),
      axis.title.y  = element_text(size = 16),
      panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
      axis.ticks.x  = element_blank(),
      axis.line.x   = element_blank(),
      plot.margin   = margin(10, 15, 10, 15)
    )
}

# ---- Save --------------------------------------------------------------------

p <- build_slopegraph(df, median_df)

save_multiformat_ggplot(
  p,
  file.path(OUT_DIR, "shared_anchor_slopegraph"),
  width = 4.5, height = 6, dpi = 300
)

p_poster <- build_slopegraph_poster(df, median_df)

REPO_ROOT  <- normalizePath(file.path(BASE, ".."))
POSTER_OUT <- file.path(REPO_ROOT, "poster/figures/shared_anchor_slopegraph_poster")
save_multiformat_ggplot(
  p_poster,
  POSTER_OUT,
  width = 6.5, height = 7, dpi = 300
)

cat("Done.\n")
