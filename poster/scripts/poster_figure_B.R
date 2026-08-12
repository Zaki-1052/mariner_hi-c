# poster/scripts/poster_figure_B.R
# Poster-ready anchor vs span ChromHMM enrichment figure with collapsed state groups
# Two-panel horizontal bar chart: clust5 (gained) vs clust6 (lost)
# 18 ChromHMM states collapsed into 4 biologically meaningful groups via
# genome-coverage-weighted averaging

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(ggpattern)
})

BASE_DIR <- getwd()
source(file.path(BASE_DIR, "data/scripts/_shared/multi_format_output.R"))

# ---------------------------------------------------------------------------
# 1. Read input files
# ---------------------------------------------------------------------------

anchor_path <- file.path(
  BASE_DIR,
  "cluster/outputs/bap1_late/chromHMM_9mark_intersect/anchor_18.txt"
)
span_path <- file.path(
  BASE_DIR,
  "cluster/outputs/bap1_late/chromHMM_9mark_intersect/span_18.txt"
)
overlap_path <- file.path(
  BASE_DIR,
  "cluster/outputs/bap1_late/chromHMM_9mark_intersect/learned_model_18/cerebellum_late_18_overlap.txt"
)

stopifnot(file.exists(anchor_path))
stopifnot(file.exists(span_path))
stopifnot(file.exists(overlap_path))

read_enrichment <- function(path) {
  df <- read.delim(path, header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE, check.names = FALSE)
  # Drop the "Base" row
  df <- df[df[[1]] != "Base", , drop = FALSE]
  # Extract integer state ID from "1_K119ub_Only" -> 1
  df$state_id <- as.integer(sub("_.*", "", df[[1]]))
  df
}

anchor_df <- read_enrichment(anchor_path)
span_df   <- read_enrichment(span_path)

# Read overlap/genome-percent file
overlap_raw <- read.delim(overlap_path, header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE, check.names = FALSE)
overlap_raw <- overlap_raw[overlap_raw[[1]] != "Base", , drop = FALSE]

# Build named vector: state_id -> genome percent
genome_pct <- setNames(
  as.numeric(overlap_raw[["Genome %"]]),
  as.integer(overlap_raw[[1]])
)

cat("Input files loaded.\n")
cat(sprintf("  anchor_18: %d states\n", nrow(anchor_df)))
cat(sprintf("  span_18:   %d states\n", nrow(span_df)))
cat(sprintf("  genome %%:  %d states, sum = %.2f%%\n",
            length(genome_pct), sum(genome_pct)))

# ---------------------------------------------------------------------------
# 2. Define state groupings (drop Quiescent = state 13)
# ---------------------------------------------------------------------------

state_groups <- list(
  "Polycomb/repressive" = c(1, 2, 3, 17, 18),
  "Active regulatory"   = c(6, 7, 8, 9, 11, 12),
  "Poised"              = c(4, 5, 10, 15),
  "Structural"          = c(14, 16)
)

GROUP_COLORS <- c(
  "Polycomb/repressive" = "#7b2d8e",
  "Active regulatory"   = "#d62728",
  "Poised"              = "#ff7f0e",
  "Structural"          = "#2077b4"
)

# Order for y-axis (bottom to top in the plot)
group_order <- c("Structural", "Poised", "Active regulatory", "Polycomb/repressive")

# ---------------------------------------------------------------------------
# 3. Compute genome-coverage-weighted averages per group
# ---------------------------------------------------------------------------

compute_weighted_avg <- function(enrichment_df, col_name, group_ids, genome_pct_vec) {
  # Subset to states in this group
  idx <- enrichment_df$state_id %in% group_ids
  sub_df <- enrichment_df[idx, , drop = FALSE]

  enrichments <- as.numeric(sub_df[[col_name]])
  weights     <- genome_pct_vec[as.character(sub_df$state_id)]

  sum(enrichments * weights) / sum(weights)
}

# Clusters of interest
clusters <- c("clust5.bed", "clust6.bed")
cluster_labels <- c(
  "clust5.bed" = "Gained loops (n=667)",
  "clust6.bed" = "Lost loops (n=2,359)"
)

# Build the long data frame
records <- list()
for (grp_name in names(state_groups)) {
  grp_ids <- state_groups[[grp_name]]
  for (clust_col in clusters) {
    anchor_val <- compute_weighted_avg(anchor_df, clust_col, grp_ids, genome_pct)
    span_val   <- compute_weighted_avg(span_df,   clust_col, grp_ids, genome_pct)

    records[[length(records) + 1]] <- data.frame(
      group       = grp_name,
      cluster     = cluster_labels[[clust_col]],
      region_type = "Anchor",
      enrichment  = anchor_val,
      stringsAsFactors = FALSE
    )
    records[[length(records) + 1]] <- data.frame(
      group       = grp_name,
      cluster     = cluster_labels[[clust_col]],
      region_type = "Span",
      enrichment  = span_val,
      stringsAsFactors = FALSE
    )
  }
}

plot_data <- do.call(rbind, records)

# Set factor levels for ordering
plot_data$group <- factor(plot_data$group, levels = group_order)
plot_data$region_type <- factor(plot_data$region_type, levels = c("Span", "Anchor"))
plot_data$cluster <- factor(
  plot_data$cluster,
  levels = c("Gained loops (n=667)", "Lost loops (n=2,359)")
)

# ---------------------------------------------------------------------------
# 4. Print computed values for verification
# ---------------------------------------------------------------------------

cat("\n--- Weighted-average enrichment values ---\n")
for (i in seq_len(nrow(plot_data))) {
  row <- plot_data[i, ]
  cat(sprintf("  %-22s | %-24s | %-6s | %.3f\n",
              as.character(row$group),
              as.character(row$cluster),
              as.character(row$region_type),
              row$enrichment))
}

# ---------------------------------------------------------------------------
# 5. Build ggplot: dark+solid (anchor) vs light+striped (span)
# ---------------------------------------------------------------------------

lighten <- function(hex, amount = 0.45) {
  rgb_vals <- col2rgb(hex) / 255
  light <- rgb_vals + (1 - rgb_vals) * amount
  rgb(light[1], light[2], light[3])
}

# Build fill color per group x region_type
plot_data$fill_key <- paste(plot_data$group, plot_data$region_type, sep = ".")

fill_colors <- c()
for (grp in names(GROUP_COLORS)) {
  fill_colors[[paste0(grp, ".Anchor")]] <- GROUP_COLORS[[grp]]
  fill_colors[[paste0(grp, ".Span")]]   <- lighten(GROUP_COLORS[[grp]])
}

plot_data$fill_key <- factor(plot_data$fill_key, levels = names(fill_colors))

p <- ggplot(plot_data,
            aes(x = enrichment, y = group,
                fill = fill_key, pattern = region_type)) +
  geom_col_pattern(
    position        = position_dodge(width = 0.75),
    width           = 0.65,
    color           = "grey30",
    linewidth       = 0.3,
    pattern_fill    = "white",
    pattern_colour  = "grey50",
    pattern_spacing = 0.03,
    pattern_density = 0.3,
    pattern_angle   = 45
  ) +
  geom_vline(xintercept = 1.0, linetype = "dashed", color = "grey40",
             linewidth = 0.6) +
  facet_wrap(~ cluster, ncol = 2) +
  scale_fill_manual(values = fill_colors, guide = "none") +
  scale_pattern_manual(
    values = c("Anchor" = "none", "Span" = "stripe"),
    name   = NULL
  ) +
  labs(
    x = "Fold enrichment (vs genome)",
    y = NULL
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y        = element_text(size = 14, face = "bold"),
    axis.text.x        = element_text(size = 14, face = "bold"),
    axis.title.x       = element_text(size = 14, face = "bold"),
    strip.text         = element_text(size = 16, face = "bold"),
    legend.text        = element_text(size = 14),
    legend.position    = "bottom",
    legend.key.size    = unit(1.2, "lines"),
    plot.margin        = margin(t = 10, r = 15, b = 40, l = 10, unit = "pt")
  ) +
  guides(
    pattern = guide_legend(
      override.aes = list(
        fill            = c("grey40", "grey80"),
        pattern_fill    = "white",
        pattern_spacing = 0.02
      )
    )
  )

# ---------------------------------------------------------------------------
# 6. Add mechanism labels below each panel
# ---------------------------------------------------------------------------

# Use annotation_custom via a grob approach for panel-specific labels.
# The simpler approach: use a secondary caption constructed from facet positions.
# ggplot2's facet annotation is tricky, so we use grid::textGrob injected via
# the gtable for precise panel placement.

# Build the plot first, then modify the gtable
library(grid)
library(gtable)

gt <- ggplotGrob(p)

# Find the panel positions in the gtable layout
panel_positions <- gt$layout[grepl("panel", gt$layout$name), ]

# Mechanism labels (in facet order: Gained = left, Lost = right)
mechanism_labels <- c(
  "Polycomb domain compaction",
  "Anchor disruption"
)

# Add a new row below the panels for the labels
# Find the bottom-most panel row
panel_bottom <- max(panel_positions$b)

# Add spacing row + label row below panels (before the axis)
gt <- gtable_add_rows(gt, heights = unit(0.6, "lines"), pos = panel_bottom)
gt <- gtable_add_rows(gt, heights = unit(1.2, "lines"), pos = panel_bottom + 1)

for (i in seq_len(nrow(panel_positions))) {
  label_grob <- textGrob(
    mechanism_labels[i],
    gp = gpar(fontsize = 13, fontface = "italic", col = "grey30")
  )
  gt <- gtable_add_grob(
    gt, label_grob,
    t = panel_bottom + 2,
    l = panel_positions$l[i],
    r = panel_positions$r[i],
    name = paste0("mechanism_label_", i)
  )
}

# ---------------------------------------------------------------------------
# 7. Save output
# ---------------------------------------------------------------------------

dir.create(file.path(BASE_DIR, "poster/figures"), recursive = TRUE, showWarnings = FALSE)

# save_multiformat_ggplot expects a ggplot or grob-wrappable object.
# Since we modified the gtable, use ggsave directly on the gtable via a
# thin ggplot wrapper, or call ggsave on the grob.

output_base <- file.path(BASE_DIR, "poster/figures/anchor_vs_span_enrichment")
output_dir  <- file.path(dirname(output_base), basename(output_base))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

fig_width  <- 12
fig_height <- 6

# PDF
pdf_path <- file.path(output_dir, paste0(basename(output_base), ".pdf"))
pdf(pdf_path, width = fig_width, height = fig_height)
grid.draw(gt)
dev.off()

# SVG
svg_path <- file.path(output_dir, paste0(basename(output_base), ".svg"))
svglite::svglite(svg_path, width = fig_width, height = fig_height)
grid.draw(gt)
dev.off()

# JPEG
jpg_path <- file.path(output_dir, paste0(basename(output_base), ".jpg"))
jpeg(jpg_path, width = fig_width * 300, height = fig_height * 300,
     res = 300, quality = 95)
grid.draw(gt)
dev.off()

cat(sprintf("  Saved: anchor_vs_span_enrichment/{pdf,svg,jpg}\n"))

cat("Figure B complete.\n")
