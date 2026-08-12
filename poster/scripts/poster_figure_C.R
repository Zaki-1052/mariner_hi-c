# poster/scripts/poster_figure_C.R
# Chromatin composition of gained vs lost differential loop anchors
# Version A: Grouped bar chart  |  Version B: Paired donut chart

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(patchwork)
})

BASE_DIR <- getwd()
source(file.path(BASE_DIR, "data/scripts/_shared/multi_format_output.R"))

# --- Shared palette and category definitions ---

CAT_COLORS <- c(
  "Polycomb"         = "#8E44AD",
  "Active"           = "#E74C3C",
  "Poised Enhancer"  = "#F39C12",
  "CTCF/Structural"  = "#2E86AB",
  "Other"            = "#95A5A6"
)
CAT_ORDER <- c("Polycomb", "Active", "Poised Enhancer", "CTCF/Structural", "Other")

DIRECTION_LABELS <- c(
  up_in_mutant   = "Gained in BAP1-KO",
  down_in_mutant = "Lost in BAP1-KO"
)

# --- 1. Read anchor_type_by_direction.tsv ---

input_path <- file.path(
  BASE_DIR,
  "loops/output/loop_type_summary/late/anchor_type_by_direction.tsv"
)
stopifnot(file.exists(input_path))

raw <- read.delim(input_path, stringsAsFactors = FALSE)
stopifnot(all(c("direction", "anchor_type", "count", "percentage") %in% colnames(raw)))

# --- 2. Collapse 8 anchor types into 5 categories ---

collapse_category <- function(anchor_type) {
  case_when(
    anchor_type %in% c("Polycomb", "Repressed_Promoter", "Bivalent_Promoter") ~ "Polycomb",
    anchor_type %in% c("Active_Promoter", "Active_Enhancer")                  ~ "Active",
    anchor_type == "Poised_Enhancer"                                           ~ "Poised Enhancer",
    anchor_type == "CTCF_Site"                                                 ~ "CTCF/Structural",
    TRUE                                                                       ~ "Other"
  )
}

collapsed <- raw %>%
  mutate(category = collapse_category(anchor_type)) %>%
  group_by(direction, category) %>%
  summarise(count = sum(count), .groups = "drop")

# --- 3. Recompute percentages within each direction ---

collapsed <- collapsed %>%
  group_by(direction) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup()

# Set display order
collapsed$category <- factor(collapsed$category, levels = CAT_ORDER)
collapsed$direction_label <- factor(
  DIRECTION_LABELS[collapsed$direction],
  levels = c("Gained in BAP1-KO", "Lost in BAP1-KO")
)

# Compute per-direction totals for donut center labels
direction_totals <- collapsed %>%
  group_by(direction) %>%
  summarise(total = sum(count), .groups = "drop")

# ============================================================
# VERSION A: Grouped Bar Chart
# ============================================================

p_bar <- ggplot(collapsed, aes(x = category, y = percentage, fill = category)) +
  geom_col(width = 0.7) +
  geom_text(
    aes(label = sprintf("%.1f%%", percentage)),
    vjust = -0.5,
    size  = 4.5,
    fontface = "bold"
  ) +
  facet_wrap(~direction_label, ncol = 2) +
  scale_fill_manual(values = CAT_COLORS) +
  scale_y_continuous(
    limits = c(0, 45),
    expand = c(0, 0)
  ) +
  labs(y = "% of loop anchors", x = NULL) +
  theme_minimal(base_size = 16) +
  theme(
    legend.position  = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.x = element_text(angle = 27, hjust = 1, size = 13),
    strip.text  = element_text(size = 16, face = "bold")
  )

dir.create(file.path(BASE_DIR, "poster/figures"), recursive = TRUE, showWarnings = FALSE)

save_multiformat_ggplot(
  p_bar,
  file.path(BASE_DIR, "poster/figures/chromatin_composition_bar"),
  width  = 12,
  height = 7
)

# ============================================================
# VERSION B: Paired Donut Chart
# ============================================================

# Helper: build a single donut plot for one direction
build_donut <- function(df, dir_key, totals_df) {
  total_n <- totals_df$total[totals_df$direction == dir_key]
  label   <- DIRECTION_LABELS[dir_key]

  # Order slices by CAT_ORDER and compute cumulative positions
  ring <- df %>%
    filter(direction == dir_key) %>%
    arrange(category) %>%
    mutate(
      ymax = cumsum(percentage),
      ymin = lag(ymax, default = 0),
      ymid = (ymin + ymax) / 2,
      slice_label = ifelse(percentage >= 5, sprintf("%.1f%%", percentage), "")
    )

  ggplot(ring) +
    geom_rect(
      aes(xmin = 1.2, xmax = 2.0, ymin = ymin, ymax = ymax, fill = category),
      color = "white",
      linewidth = 0.8
    ) +
    geom_text(
      aes(x = 1.6, y = ymid, label = slice_label),
      size     = 4,
      fontface = "bold"
    ) +
    annotate(
      "text",
      x     = 0.3,
      y     = 0,
      label = paste0("n=", format(total_n, big.mark = ",")),
      size  = 5.5,
      fontface = "bold"
    ) +
    coord_polar(theta = "y") +
    xlim(0.3, 2.5) +
    scale_fill_manual(values = CAT_COLORS, name = NULL) +
    labs(title = label) +
    theme_void(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
    )
}

p_gained <- build_donut(collapsed, "up_in_mutant",   direction_totals)
p_lost   <- build_donut(collapsed, "down_in_mutant", direction_totals)

p_donut <- (p_gained + p_lost) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

save_multiformat_ggplot(
  p_donut,
  file.path(BASE_DIR, "poster/figures/chromatin_composition_donut"),
  width  = 12,
  height = 6
)

cat("Figure C complete (both versions).\n")
