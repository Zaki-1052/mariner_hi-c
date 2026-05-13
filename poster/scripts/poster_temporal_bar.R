# poster/scripts/poster_temporal_bar.R
# Beat 3 poster figure: paired stacked bars showing 18x temporal amplification
# Run from repo root: Rscript poster/scripts/poster_temporal_bar.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(scales)
})

BASE_DIR <- getwd()
source(file.path(BASE_DIR, "data/scripts/_shared/multi_format_output.R"))

# ---- Configuration ----------------------------------------------------------

EARLY_TSV <- file.path(BASE_DIR, "data/upstream/loop_calls/early_characterized_loops.tsv")
LATE_TSV  <- file.path(BASE_DIR, "data/upstream/loop_calls/late_characterized_loops.tsv")
OUT_BASE  <- file.path(BASE_DIR, "poster/figures/temporal_loop_bar")

stopifnot(file.exists(EARLY_TSV), file.exists(LATE_TSV))

COL_LOST   <- "#d73027"
COL_GAINED <- "#4575b4"

# ---- Data loading -----------------------------------------------------------

cat("Loading loop counts...\n")
early <- read_tsv(EARLY_TSV, col_select = "direction", show_col_types = FALSE)
late  <- read_tsv(LATE_TSV,  col_select = "direction", show_col_types = FALSE)

n_early_lost   <- sum(early$direction == "down_in_mutant")
n_early_gained <- sum(early$direction == "up_in_mutant")
n_late_lost    <- sum(late$direction == "down_in_mutant")
n_late_gained  <- sum(late$direction == "up_in_mutant")

total_early <- n_early_lost + n_early_gained
total_late  <- n_late_lost  + n_late_gained
fold_change <- round(total_late / total_early)

cat(sprintf("  P12:   %d lost + %d gained = %d\n", n_early_lost, n_early_gained, total_early))
cat(sprintf("  Adult: %d lost + %d gained = %d\n", n_late_lost, n_late_gained, total_late))
cat(sprintf("  Fold change: %d×\n", fold_change))

stopifnot(total_early == 165, total_late == 2910)

# ---- Theme ------------------------------------------------------------------

theme_poster <- function() {
  theme_minimal(base_size = 20) +
    theme(
      plot.title    = element_text(face = "bold", size = 22, hjust = 0.5),
      plot.subtitle = element_text(size = 16, hjust = 0.5, color = "gray40"),
      axis.title    = element_text(face = "bold", size = 18),
      axis.text     = element_text(size = 16),
      legend.text   = element_text(size = 16),
      legend.title  = element_text(size = 13, color = "gray40", hjust = 0.5),
      legend.title.position = "bottom",
      legend.position  = "top",
      panel.grid.minor = element_blank()
    )
}

# ---- Build plot dataframe ---------------------------------------------------

plot_df <- tibble(
  timepoint = factor(c("P12", "P12", "Adult", "Adult"),
                     levels = c("P12", "Adult")),
  direction = factor(c("Lost", "Gained", "Lost", "Gained"),
                     levels = c("Lost", "Gained")),
  count     = c(n_early_lost, n_early_gained, n_late_lost, n_late_gained)
)

totals_df <- tibble(
  timepoint = factor(c("P12", "Adult"), levels = c("P12", "Adult")),
  total     = c(total_early, total_late)
)

# ---- Plot -------------------------------------------------------------------

cat("Building figure...\n")

arrow_y <- total_late * 1.16
label_y <- total_late * 1.28

p <- ggplot(plot_df, aes(x = timepoint, y = count, fill = direction)) +
  geom_col(width = 0.45, color = "white", linewidth = 0.6) +
  scale_fill_manual(
    values = c("Lost" = COL_LOST, "Gained" = COL_GAINED),
    name   = NULL
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.38)),
    labels = comma
  ) +
  geom_text(
    data      = totals_df,
    aes(x = timepoint, y = total, label = comma(total)),
    inherit.aes = FALSE,
    vjust     = -0.4,
    size      = 5.5,
    fontface  = "bold",
    color     = "grey20"
  ) +
  annotate(
    "segment",
    x = 1, xend = 2,
    y = arrow_y, yend = arrow_y,
    color     = "grey30",
    linewidth = 0.7
  ) +
  annotate(
    "segment",
    x = 2, xend = 2,
    y = arrow_y, yend = total_late * 1.09,
    arrow     = arrow(length = unit(0.12, "inches"), type = "closed"),
    color     = "grey30",
    linewidth = 0.7
  ) +
  annotate(
    "text",
    x     = 1.5,
    y     = label_y,
    label = paste0(fold_change, "×"),
    size  = 7,
    fontface = "bold",
    color = "grey20",
    hjust = 0.5
  ) +
  labs(x = NULL, y = "Differential loops") +
  theme_poster() +
  theme(panel.grid.major.x = element_blank())

# ---- Save -------------------------------------------------------------------

save_multiformat_ggplot(p, OUT_BASE, width = 4, height = 5, dpi = 300)
cat("Done.\n")
