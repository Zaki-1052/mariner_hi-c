# poster/scripts/poster_figure_A.R
# Poster-ready 2-row cluster feature heatmap: Gained loops (clust5) vs Lost loops (clust6)
# Z-score normalized tiles with raw value annotations

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
})

BASE_DIR <- getwd()
source(file.path(BASE_DIR, "data/scripts/_shared/multi_format_output.R"))

# --- 1. Read data (skip multi-line quoted headers, parse by position) ---

input_path <- file.path(
  BASE_DIR,
  "cluster/outputs/bap1_late/figures/summary_figures/heatmap/feature_values.tsv"
)
stopifnot(file.exists(input_path))

lines <- readLines(input_path)
data_lines <- lines[grep("^clust", lines)]
raw <- read.table(text = data_lines, sep = "\t", stringsAsFactors = FALSE)
# V1  = cluster name
# V2  = Median logFC
# V3  = % Gained
# V4  = % Lost
# V5  = Median size (kb)
# V6  = % Structural
# V7  = % CRE
# V8  = Polycomb anchor enrichment
# V9  = Polycomb span enrichment
# V10 = Bivalent anchor enrichment
# V11 = % Polycomb (anchors)
# V12 = K119ub ctrl signal
# V13 = K119ub mut signal

rownames(raw) <- raw$V1

# --- 2. Extract clust5 and clust6 rows ---

stopifnot(all(c("clust5", "clust6") %in% raw$V1))

# --- 3. Build full 6-cluster x 7-feature matrix (z-score across all 6) ---

feature_names <- c(
  "n", "Median\nlogFC", "Median\nsize (kb)",
  "% CRE", "Polycomb\nanchor", "Polycomb\nspan", "K119ub\n(mutant)"
)

# Hardcoded sample sizes per cluster (not in the input file)
cluster_n <- c(
  clust1 = 12298, clust2 = 10970, clust3 = 8738,
  clust4 = 3916,  clust5 = 667,   clust6 = 2359
)

all_clusters <- raw$V1
mat_all <- data.frame(
  cluster         = all_clusters,
  median_logfc    = raw$V2,
  n               = cluster_n[all_clusters],
  median_size_kb  = raw$V5,
  pct_cre         = raw$V7,
  polycomb_anchor = raw$V8,
  polycomb_span   = raw$V9,
  k119ub_mut      = raw$V13,
  stringsAsFactors = FALSE
)

# --- 4. Z-score normalize each feature across all 6 clusters, then subset ---

value_cols <- c(
  "n", "median_logfc", "median_size_kb",
  "pct_cre", "polycomb_anchor", "polycomb_span", "k119ub_mut"
)

zscore_all <- mat_all
for (col in value_cols) {
  zscore_all[[col]] <- as.numeric(scale(mat_all[[col]]))
}

# Subset to clust5 and clust6 for display
show_clusters <- c("clust5", "clust6")
mat <- mat_all[mat_all$cluster %in% show_clusters, ]
zscore_mat <- zscore_all[zscore_all$cluster %in% show_clusters, ]

# --- 5. Format raw annotations per feature ---

format_raw <- function(col_name, value) {
  switch(col_name,
    median_logfc    = sprintf("%+.2f", value),
    n               = format(as.integer(value), big.mark = ","),
    median_size_kb  = sprintf("%.0f", value),
    pct_cre         = sprintf("%.1f%%", value),
    polycomb_anchor = sprintf("%.2fx", value),
    polycomb_span   = sprintf("%.2fx", value),
    k119ub_mut      = sprintf("%.2f", value),
    as.character(value)
  )
}

# --- 6. Pivot to long format ---

# Z-scores long
zscore_long <- zscore_mat %>%
  pivot_longer(
    cols      = all_of(value_cols),
    names_to  = "feature_id",
    values_to = "zscore"
  )

# Raw values long
raw_long <- mat %>%
  pivot_longer(
    cols      = all_of(value_cols),
    names_to  = "feature_id",
    values_to = "raw_value"
  )

# Merge and add formatted labels
plot_df <- left_join(zscore_long, raw_long, by = c("cluster", "feature_id"))

plot_df$label <- mapply(format_raw, plot_df$feature_id, plot_df$raw_value)

# n is sample size, not a biological feature — use neutral grey instead of z-score
plot_df$zscore[plot_df$feature_id == "n"] <- NA

# Map internal column names to display names
feature_display <- setNames(feature_names, value_cols)
plot_df$feature <- feature_display[plot_df$feature_id]

# Set factor levels for axis ordering
plot_df$feature <- factor(plot_df$feature, levels = feature_names)

# Row labels (simple — n already shown in its own column)
row_labels <- c(clust5 = "Gained", clust6 = "Lost")
plot_df$row_label <- row_labels[plot_df$cluster]

# Y-axis order: Lost on top, Gained on bottom
plot_df$row_label <- factor(
  plot_df$row_label,
  levels = c("Gained", "Lost")
)

# --- 7. Build ggplot ---

p <- ggplot(plot_df, aes(x = feature, y = row_label, fill = zscore)) +
  geom_tile(color = "white", linewidth = 1.5) +
  geom_text(
    aes(label = label),
    fontface = "bold",
    size     = 5
  ) +
  scale_fill_gradient2(
    low      = "#2166ac",
    mid      = "white",
    high     = "#b2182b",
    midpoint = 0,
    name     = "Z-score",
    na.value = "grey85"
  ) +
  labs(title = "Differential loop cluster features") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid      = element_blank(),
    axis.title      = element_blank(),
    axis.text.x     = element_text(size = 12, face = "bold", lineheight = 0.9),
    axis.text.y     = element_text(size = 16, face = "bold"),
    plot.title      = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.position = "right"
  )

# --- 8. Save via multi-format output ---

dir.create(file.path(BASE_DIR, "poster/figures"), recursive = TRUE, showWarnings = FALSE)
save_multiformat_ggplot(
  p,
  file.path(BASE_DIR, "poster/figures/cluster_feature_heatmap"),
  width  = 12,
  height = 3.5
)

cat("Figure A complete.\n")
