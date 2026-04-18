# scripts/shared_anchor_violin_variants.R
#
# Standalone regeneration of the "Loop Distance at Shared Anchors" violin
# for the LATE timepoint in four layout/color variants so the user can
# pick the most readable version.
#
# Variants (order | Lost color / Gained color):
#   v1_current        Lost   | Gained    blue   / red     (current)
#   v2_swap_order     Gained | Lost      blue   / red
#   v3_swap_colors    Lost   | Gained    red    / blue
#   v4_swap_both      Gained | Lost      red    / blue
#
# Reads pre-computed tables; does not recompute statistics on raw loops
# beyond the Mann-Whitney test already used in the main script.

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
})

source("scripts/utils/multi_format_output.R")

# ---- Inputs ------------------------------------------------------------------
BASE <- getwd()
LOOPS_TSV   <- file.path(BASE, "output/shared_anchor_analysis/late/tables/shared_anchor_loops.tsv")
ANCHORS_TSV <- file.path(BASE, "output/shared_anchor_analysis/late/tables/shared_anchors.tsv")
OUT_DIR     <- file.path(BASE, "output/shared_anchor_analysis/late/plots/distance_violin_shared_variants")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- Load & prep -------------------------------------------------------------
loops <- read_tsv(LOOPS_TSV, show_col_types = FALSE) %>%
  filter(anchor_status == "shared") %>%
  mutate(direction_label = case_when(
    direction == "down_in_mutant" ~ "Lost",
    direction == "up_in_mutant"   ~ "Gained",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(direction_label))

n_shared <- nrow(read_tsv(ANCHORS_TSV, show_col_types = FALSE))
n_lost   <- sum(loops$direction_label == "Lost")
n_gained <- sum(loops$direction_label == "Gained")

wt <- wilcox.test(
  loops$loop_distance[loops$direction_label == "Lost"],
  loops$loop_distance[loops$direction_label == "Gained"],
  alternative = "greater"
)

# ---- Plot builder ------------------------------------------------------------
# order_vec: c("Lost","Gained") or c("Gained","Lost")
# color_map: named vector mapping "Lost"/"Gained" -> hex
build_plot <- function(order_vec, color_map, title_suffix) {
  labels_named <- setNames(
    sprintf("%s\n(n = %d)", order_vec,
            c(Lost = n_lost, Gained = n_gained)[order_vec]),
    order_vec
  )
  plot_data <- loops %>%
    mutate(direction_label = factor(
      direction_label, levels = order_vec, labels = labels_named[order_vec]
    ))
  fill_vals <- setNames(
    unname(color_map[order_vec]),
    labels_named[order_vec]
  )

  ggplot(plot_data, aes(x = direction_label, y = loop_distance / 1e6,
                        fill = direction_label)) +
    geom_violin(alpha = 0.7, trim = FALSE, color = "black", linewidth = 0.5) +
    geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA,
                 color = "black", linewidth = 0.5) +
    scale_y_log10() +
    scale_fill_manual(values = fill_vals) +
    labs(
      title = "Loop Distance at Shared Anchors",
      subtitle = sprintf(
        "%d shared anchor regions  |  Mann-Whitney p = %.2e (lost > gained)  |  %s",
        n_shared, wt$p.value, title_suffix
      ),
      x = NULL, y = "Loop Distance (Mb, log scale)"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
      legend.position = "none"
    )
}

BLUE <- "#4575b4"; RED <- "#d73027"

variants <- list(
  list(name = "v1_current",
       order = c("Lost", "Gained"),
       colors = c(Lost = BLUE, Gained = RED),
       suffix = "current (Lost blue left, Gained red right)"),
  list(name = "v2_swap_order",
       order = c("Gained", "Lost"),
       colors = c(Lost = BLUE, Gained = RED),
       suffix = "swap order (Gained red left, Lost blue right)"),
  list(name = "v3_swap_colors",
       order = c("Lost", "Gained"),
       colors = c(Lost = RED, Gained = BLUE),
       suffix = "swap colors (Lost red left, Gained blue right)"),
  list(name = "v4_swap_both",
       order = c("Gained", "Lost"),
       colors = c(Lost = RED, Gained = BLUE),
       suffix = "swap both (Gained blue left, Lost red right)")
)

for (v in variants) {
  p <- build_plot(v$order, v$colors, v$suffix)
  save_multiformat_ggplot(
    p,
    file.path(OUT_DIR, v$name),
    width = 5, height = 5, dpi = 300
  )
  cat(sprintf("Saved: %s/{pdf,svg,jpg}\n", v$name))
}

cat(sprintf("\nAll variants written to: %s\n", OUT_DIR))
cat(sprintf("n_shared = %d | n_lost loops = %d | n_gained loops = %d | p = %.3e\n",
            n_shared, n_lost, n_gained, wt$p.value))
