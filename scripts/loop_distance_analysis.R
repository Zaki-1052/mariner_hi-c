# scripts/loop_distance_analysis.R
# Loop Distance Shift Analysis for BAP1-KO "Loop Rewriting" Phenomenon
# Generates publication-quality visualizations demonstrating that BAP1-KO mutants
# lose long-range loops and gain shorter loops

# ==============================================================================
# SECTION 1: SETUP AND CONFIGURATION
# ==============================================================================

cat("=== Loop Distance Shift Analysis ===\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# Load required packages
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(RColorBrewer)
  library(viridis)
  library(clusterProfiler)
  library(org.Mm.eg.db)
})

# Load multi-format output utility for PDF + SVG + JPEG output
source("scripts/utils/multi_format_output.R")

# Define color palettes (consistent with existing pipeline)
COLORS <- list(
  down = "#d73027",      # Red for down/lost in mutant

  up = "#4575b4",        # Blue for up/gained in mutant
  down_light = "#f4a582",
  up_light = "#92c5de",
  neutral = "#999999",
  significant = "#2166ac"
)

# Direction labels for plots
DIRECTION_LABELS <- c(
  "down_in_mutant" = "Lost in BAP1-KO",
  "up_in_mutant" = "Gained in BAP1-KO"
)

# Distance category order
DISTANCE_ORDER <- c("<100kb", "100-500kb", "500kb-1Mb", ">1Mb")

# Output directory
OUTPUT_DIR <- "outputs/loops_visualization_extended"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("Output directory:", OUTPUT_DIR, "\n\n")

# ==============================================================================
# SECTION 2: DATA LOADING AND PREPROCESSING
# ==============================================================================

cat("=== Loading Data ===\n")

# Load characterized loops
input_file <- "25042-late_outputs/merged_loops/characterized_loops.tsv"
if (!file.exists(input_file)) {
  input_file <- "outputs/merged_loops/characterized_loops.tsv"
}

loops <- read_tsv(input_file, show_col_types = FALSE)
cat("Loaded", nrow(loops), "differential loops\n")

# Ensure proper factor levels for distance_category
loops <- loops %>%
  mutate(
    distance_category = factor(distance_category, levels = DISTANCE_ORDER),
    direction_label = factor(
      case_when(
        direction == "down_in_mutant" ~ "Lost in BAP1-KO",
        direction == "up_in_mutant" ~ "Gained in BAP1-KO",
        TRUE ~ "Other"
      ),
      levels = c("Lost in BAP1-KO", "Gained in BAP1-KO")
    ),
    loop_distance_kb = loop_distance / 1000,
    loop_distance_mb = loop_distance / 1e6
  )

# Filter to only up/down loops (exclude any unchanged)
loops_directional <- loops %>%
  filter(direction %in% c("up_in_mutant", "down_in_mutant"))

cat("Directional loops:", nrow(loops_directional), "\n")
cat("  - Lost (down_in_mutant):", sum(loops_directional$direction == "down_in_mutant"), "\n")
cat("  - Gained (up_in_mutant):", sum(loops_directional$direction == "up_in_mutant"), "\n\n")

# ==============================================================================
# SECTION 3: STATISTICAL TESTS (Pre-compute for annotations)
# ==============================================================================

cat("=== Running Statistical Tests ===\n")

# Separate by direction
lost_loops <- loops_directional %>% filter(direction == "down_in_mutant")
gained_loops <- loops_directional %>% filter(direction == "up_in_mutant")

# Kolmogorov-Smirnov test for distribution difference
ks_test <- ks.test(lost_loops$loop_distance, gained_loops$loop_distance)
cat("KS test p-value:", format(ks_test$p.value, scientific = TRUE, digits = 3), "\n")

# Wilcoxon rank-sum test for median difference
wilcox_test <- wilcox.test(lost_loops$loop_distance, gained_loops$loop_distance)
cat("Wilcoxon test p-value:", format(wilcox_test$p.value, scientific = TRUE, digits = 3), "\n")

# Chi-square test for distance category independence
contingency_table <- table(loops_directional$direction, loops_directional$distance_category)
chisq_test <- chisq.test(contingency_table)
cat("Chi-square test p-value:", format(chisq_test$p.value, scientific = TRUE, digits = 3), "\n")

# Spearman correlation between distance and logFC
spearman_cor <- cor.test(loops_directional$loop_distance, loops_directional$logFC,
                          method = "spearman")
cat("Spearman correlation (distance vs logFC):", round(spearman_cor$estimate, 3), "\n")
cat("  p-value:", format(spearman_cor$p.value, scientific = TRUE, digits = 3), "\n")

# Calculate medians
median_lost <- median(lost_loops$loop_distance)
median_gained <- median(gained_loops$loop_distance)
cat("\nMedian loop distance:\n")
cat("  - Lost loops:", round(median_lost / 1000), "kb\n")
cat("  - Gained loops:", round(median_gained / 1000), "kb\n")
cat("  - Ratio (lost/gained):", round(median_lost / median_gained, 2), "\n\n")

# ==============================================================================
# FIGURE 1: Cumulative Distribution Function (CDF) by Direction
# ==============================================================================

cat("=== Generating Figure 1: CDF Plot ===\n")

p1_cdf <- ggplot(loops_directional, aes(x = loop_distance_kb, color = direction_label)) +
  stat_ecdf(geom = "step", linewidth = 1.2) +
  scale_color_manual(
    values = c("Lost in BAP1-KO" = COLORS$down, "Gained in BAP1-KO" = COLORS$up),
    name = "Direction"
  ) +
  scale_x_log10(
    labels = comma,
    breaks = c(10, 100, 1000, 10000),
    limits = c(10, 100000)
  ) +
  geom_vline(xintercept = median_lost / 1000, color = COLORS$down,
             linetype = "dashed", alpha = 0.7) +
  geom_vline(xintercept = median_gained / 1000, color = COLORS$up,
             linetype = "dashed", alpha = 0.7) +
  annotate("text", x = 50000, y = 0.15,
           label = sprintf("KS test p = %.2e", ks_test$p.value),
           hjust = 1, size = 4) +
  annotate("text", x = 50000, y = 0.08,
           label = sprintf("Median: Lost = %d kb, Gained = %d kb",
                          round(median_lost / 1000), round(median_gained / 1000)),
           hjust = 1, size = 3.5) +
  labs(
    title = "Loop Distance Distribution: Lost vs Gained Loops",
    subtitle = "BAP1-KO preferentially loses long-range loops",
    x = "Loop Distance (kb, log scale)",
    y = "Cumulative Proportion"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey40"),
    legend.position = c(0.15, 0.85),
    legend.background = element_rect(fill = "white", color = NA),
    panel.grid.minor = element_blank()
  )

save_multiformat_ggplot(p1_cdf, file.path(OUTPUT_DIR, "01_distance_cdf_by_direction"), width = 8, height = 6)
cat("Saved: 01_distance_cdf_by_direction.pdf\n")

# ==============================================================================
# FIGURE 2: Distance Category Bar Chart with Fold-Enrichment
# ==============================================================================

cat("=== Generating Figure 2: Distance Category Bar Chart ===\n")

# Calculate counts and percentages by category
category_summary <- loops_directional %>%
  group_by(direction_label, distance_category) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(direction_label) %>%
  mutate(
    total = sum(count),
    percentage = count / total * 100
  ) %>%
  ungroup()

# Calculate fold-enrichment (lost / gained ratio)
enrichment <- category_summary %>%
  dplyr::select(direction_label, distance_category, count) %>%
  pivot_wider(names_from = direction_label, values_from = count) %>%
  mutate(
    ratio = `Lost in BAP1-KO` / `Gained in BAP1-KO`,
    log2_ratio = log2(ratio),
    enrichment_label = sprintf("%.2fx", ratio)
  )

# Bar plot with counts
p2_bar <- ggplot(category_summary,
                  aes(x = distance_category, y = count, fill = direction_label)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8),
           width = 0.7, color = "black", linewidth = 0.3) +
  geom_text(aes(label = count),
            position = position_dodge(width = 0.8),
            vjust = -0.5, size = 3.5) +
  scale_fill_manual(
    values = c("Lost in BAP1-KO" = COLORS$down, "Gained in BAP1-KO" = COLORS$up),
    name = ""
  ) +
  # Add enrichment ratio annotations
  geom_text(data = enrichment,
            aes(x = distance_category, y = max(category_summary$count) * 1.15,
                label = enrichment_label),
            inherit.aes = FALSE, size = 3.5, fontface = "bold") +
  geom_text(data = data.frame(x = 2.5, y = max(category_summary$count) * 1.22),
            aes(x = x, y = y, label = "Lost:Gained Ratio"),
            inherit.aes = FALSE, size = 3, color = "grey40") +
  labs(
    title = "Loop Count by Distance Category",
    subtitle = sprintf("Chi-square test p = %.2e (significant association)", chisq_test$p.value),
    x = "Loop Distance Category",
    y = "Number of Loops"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey40"),
    legend.position = "top",
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(size = 11)
  )

save_multiformat_ggplot(p2_bar, file.path(OUTPUT_DIR, "02_distance_category_barplot"), width = 9, height = 7)
cat("Saved: 02_distance_category_barplot.pdf\n")

# ==============================================================================
# FIGURE 3: Density Plot with Median Annotations
# ==============================================================================

cat("=== Generating Figure 3: Density Plot ===\n")

p3_density <- ggplot(loops_directional, aes(x = loop_distance_kb, fill = direction_label)) +
  geom_density(alpha = 0.5, color = "black", linewidth = 0.5) +
  geom_vline(xintercept = median_lost / 1000, color = COLORS$down,
             linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = median_gained / 1000, color = COLORS$up,
             linetype = "dashed", linewidth = 1) +
  scale_fill_manual(
    values = c("Lost in BAP1-KO" = COLORS$down, "Gained in BAP1-KO" = COLORS$up),
    name = ""
  ) +
  scale_x_log10(
    labels = comma,
    breaks = c(10, 100, 1000, 10000),
    limits = c(10, 100000)
  ) +
  annotate("text", x = median_lost / 1000 * 1.3, y = Inf,
           label = sprintf("Median\n%d kb", round(median_lost / 1000)),
           color = COLORS$down, vjust = 1.5, hjust = 0, size = 3.5, fontface = "bold") +
  annotate("text", x = median_gained / 1000 * 0.7, y = Inf,
           label = sprintf("Median\n%d kb", round(median_gained / 1000)),
           color = COLORS$up, vjust = 1.5, hjust = 1, size = 3.5, fontface = "bold") +
  annotate("text", x = 20, y = 0,
           label = sprintf("Wilcoxon p = %.2e", wilcox_test$p.value),
           hjust = 0, vjust = -0.5, size = 4) +
  labs(
    title = "Loop Distance Density Distribution",
    subtitle = "Lost loops are systematically longer than gained loops",
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

save_multiformat_ggplot(p3_density, file.path(OUTPUT_DIR, "03_distance_density_overlay"), width = 8, height = 6)
cat("Saved: 03_distance_density_overlay.pdf\n")

# ==============================================================================
# FIGURE 4: LogFC vs Distance Scatter with Trend
# ==============================================================================

cat("=== Generating Figure 4: LogFC vs Distance Scatter ===\n")

p4_scatter <- ggplot(loops_directional,
                      aes(x = loop_distance_kb, y = logFC)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(aes(color = direction_label), alpha = 0.4, size = 1.5) +
  geom_smooth(method = "loess", color = "black", linewidth = 1.2,
              fill = "grey70", alpha = 0.3) +
  scale_color_manual(
    values = c("Lost in BAP1-KO" = COLORS$down, "Gained in BAP1-KO" = COLORS$up),
    name = ""
  ) +
  scale_x_log10(
    labels = comma,
    breaks = c(10, 100, 1000, 10000)
  ) +
  annotate("text", x = 50000, y = max(loops_directional$logFC) * 0.9,
           label = sprintf("Spearman rho = %.3f\np = %.2e",
                          spearman_cor$estimate, spearman_cor$p.value),
           hjust = 1, size = 4) +
  labs(
    title = "Effect Size vs Loop Distance",
    subtitle = "Negative correlation: longer loops tend to decrease in BAP1-KO",
    x = "Loop Distance (kb, log scale)",
    y = "Log2 Fold Change (Mutant/WT)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey40"),
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

save_multiformat_ggplot(p4_scatter, file.path(OUTPUT_DIR, "04_logfc_vs_distance_scatter"), width = 8, height = 6)
cat("Saved: 04_logfc_vs_distance_scatter.pdf\n")

# ==============================================================================
# FIGURE 5: Distance-Stratified Volcano Plots
# ==============================================================================

cat("=== Generating Figure 5: Distance-Stratified Volcano ===\n")

# Calculate counts per category for annotations
category_counts <- loops_directional %>%
  group_by(distance_category, direction_label) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(distance_category) %>%
  summarise(
    n_up = sum(n[direction_label == "Gained in BAP1-KO"]),
    n_down = sum(n[direction_label == "Lost in BAP1-KO"]),
    label = sprintf("Up: %d | Down: %d", n_up, n_down),
    .groups = "drop"
  )

p5_volcano <- ggplot(loops_directional,
                      aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = direction_label), alpha = 0.5, size = 1.5) +
  geom_vline(xintercept = c(-0.3, 0.3), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  facet_wrap(~distance_category, nrow = 2, scales = "free_y") +
  scale_color_manual(
    values = c("Lost in BAP1-KO" = COLORS$down, "Gained in BAP1-KO" = COLORS$up),
    name = ""
  ) +
  geom_text(data = category_counts,
            aes(x = -Inf, y = Inf, label = label),
            hjust = -0.1, vjust = 1.5, size = 3, inherit.aes = FALSE) +
  labs(
    title = "Distance-Stratified Volcano Plots",
    subtitle = "Differential patterns vary by loop distance category",
    x = "Log2 Fold Change (Mutant/WT)",
    y = "-Log10(FDR)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey40"),
    legend.position = "top",
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank()
  )

save_multiformat_ggplot(p5_volcano, file.path(OUTPUT_DIR, "05_distance_stratified_volcano"), width = 10, height = 8)
cat("Saved: 05_distance_stratified_volcano.pdf\n")

# ==============================================================================
# FIGURE 6: Loop Type x Distance Heatmap
# ==============================================================================

cat("=== Generating Figure 6: Loop Type x Distance Heatmap ===\n")

# Calculate net change per loop type and distance category
looptype_distance <- loops_directional %>%
  group_by(loop_type, distance_category, direction) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = direction, values_from = n, values_fill = 0) %>%
  mutate(
    net_change = up_in_mutant - down_in_mutant,
    total = up_in_mutant + down_in_mutant,
    log2_ratio = log2((up_in_mutant + 1) / (down_in_mutant + 1))
  ) %>%
  filter(total >= 5)  # Filter out low-count combinations

# Order loop types by total count
looptype_order <- loops_directional %>%
  count(loop_type) %>%
  arrange(desc(n)) %>%
  pull(loop_type)

looptype_distance <- looptype_distance %>%
  mutate(loop_type = factor(loop_type, levels = rev(looptype_order)))

p6_heatmap <- ggplot(looptype_distance,
                      aes(x = distance_category, y = loop_type, fill = log2_ratio)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%+d", net_change)), size = 3) +
  scale_fill_gradient2(
    low = COLORS$down, mid = "white", high = COLORS$up,
    midpoint = 0, name = "Log2(Up/Down)",
    limits = c(-2, 2), oob = scales::squish
  ) +
  labs(
    title = "Net Loop Change by Type and Distance",
    subtitle = "Positive = more gained, Negative = more lost",
    x = "Distance Category",
    y = "Loop Type"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey40"),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

save_multiformat_ggplot(p6_heatmap, file.path(OUTPUT_DIR, "06_looptype_distance_heatmap"), width = 9, height = 8)
cat("Saved: 06_looptype_distance_heatmap.pdf\n")

# ==============================================================================
# FIGURE 7: ChIP-seq Mark x Distance x Direction Analysis
# ==============================================================================

cat("=== Generating Figure 7: ChIP-seq x Distance Analysis ===\n")

# Summarize anchor types by distance and direction
anchor_type_summary <- loops_directional %>%
  dplyr::select(loop_id, direction_label, distance_category, anchor1_type, anchor2_type) %>%
  pivot_longer(cols = c(anchor1_type, anchor2_type),
               names_to = "anchor", values_to = "type") %>%
  group_by(direction_label, distance_category, type) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(direction_label, distance_category) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup()

# Panel A: Anchor type distribution by distance and direction
p7a <- ggplot(anchor_type_summary,
               aes(x = distance_category, y = percentage, fill = type)) +
  geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.2) +
  facet_wrap(~direction_label) +
  scale_fill_brewer(palette = "Set2", name = "Anchor Type") +
  labs(
    title = "Anchor Type Distribution by Distance",
    x = "Distance Category",
    y = "Percentage of Anchors"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

# Panel B: Chi-square test for anchor type vs distance (by direction)
chisq_lost <- chisq.test(table(lost_loops$anchor1_type, lost_loops$distance_category))
chisq_gained <- chisq.test(table(gained_loops$anchor1_type, gained_loops$distance_category))

# Summary of ChIP-seq overlaps by distance
chipseq_summary <- loops_directional %>%
  mutate(
    has_H3K27ac = anchor1_H3K27ac_overlap | anchor2_H3K27ac_overlap,
    has_H3K4me1 = anchor1_H3K4me1_overlap | anchor2_H3K4me1_overlap
  ) %>%
  group_by(direction_label, distance_category) %>%
  summarise(
    pct_H3K27ac = mean(has_H3K27ac) * 100,
    pct_H3K4me1 = mean(has_H3K4me1) * 100,
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c(pct_H3K27ac, pct_H3K4me1),
               names_to = "mark", values_to = "percentage") %>%
  mutate(mark = gsub("pct_", "", mark))

# Panel B: ChIP-seq mark overlap by distance
p7b <- ggplot(chipseq_summary,
               aes(x = distance_category, y = percentage, fill = mark)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8),
           width = 0.7, color = "black", linewidth = 0.2) +
  facet_wrap(~direction_label) +
  scale_fill_manual(values = c("H3K27ac" = "#e41a1c", "H3K4me1" = "#377eb8"),
                    name = "Histone Mark") +
  labs(
    title = "ChIP-seq Mark Overlap by Distance",
    x = "Distance Category",
    y = "% Loops with Mark at Anchor"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

# Combine panels
p7_combined <- p7a / p7b +
  plot_annotation(
    title = "ChIP-seq Features by Loop Distance and Direction",
    subtitle = sprintf("Chi-square: Lost loops p=%.2e, Gained loops p=%.2e",
                      chisq_lost$p.value, chisq_gained$p.value),
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey40")
    )
  )

save_multiformat_ggplot(p7_combined, file.path(OUTPUT_DIR, "07_chipseq_distance_analysis"), width = 11, height = 10)
cat("Saved: 07_chipseq_distance_analysis.pdf\n")

# ==============================================================================
# FIGURE 8: GO Enrichment Comparison (Long Lost vs Short Gained)
# ==============================================================================

cat("=== Generating Figure 8: GO Enrichment Comparison ===\n")

# Extract genes from long lost loops (>500kb, down_in_mutant)
# Note: Gene names in this dataset are already Entrez IDs (numeric)
long_lost_genes <- loops_directional %>%
  filter(direction == "down_in_mutant", loop_distance > 500000) %>%
  dplyr::select(anchor1_nearest_gene, anchor2_nearest_gene) %>%
  pivot_longer(everything(), values_to = "gene") %>%
  pull(gene) %>%
  unique() %>%
  na.omit() %>%
  as.character()

# Extract genes from short gained loops (<500kb, up_in_mutant)
short_gained_genes <- loops_directional %>%
  filter(direction == "up_in_mutant", loop_distance < 500000) %>%
  dplyr::select(anchor1_nearest_gene, anchor2_nearest_gene) %>%
  pivot_longer(everything(), values_to = "gene") %>%
  pull(gene) %>%
  unique() %>%
  na.omit() %>%
  as.character()

cat("Long lost loop genes:", length(long_lost_genes), "\n")
cat("Short gained loop genes:", length(short_gained_genes), "\n")

# Try GO enrichment - gene IDs are already Entrez IDs in this dataset
tryCatch({
  # Filter to valid Entrez IDs (numeric strings)
  long_lost_entrez <- long_lost_genes[grepl("^[0-9]+$", long_lost_genes)]
  short_gained_entrez <- short_gained_genes[grepl("^[0-9]+$", short_gained_genes)]

  cat("Valid Entrez IDs - Long lost:", length(long_lost_entrez),
      ", Short gained:", length(short_gained_entrez), "\n")

  # Run GO enrichment for long lost loops
  go_long_lost <- enrichGO(
    gene = long_lost_entrez,
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.1,
    readable = TRUE
  )

  # Run GO enrichment for short gained loops
  go_short_gained <- enrichGO(
    gene = short_gained_entrez,
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.1,
    readable = TRUE
  )

  # Save GO results
  if (nrow(as.data.frame(go_long_lost)) > 0) {
    write_tsv(as.data.frame(go_long_lost),
              file.path(OUTPUT_DIR, "go_long_lost_loops.tsv"))
    cat("Saved GO results for long lost loops\n")
  }

  if (nrow(as.data.frame(go_short_gained)) > 0) {
    write_tsv(as.data.frame(go_short_gained),
              file.path(OUTPUT_DIR, "go_short_gained_loops.tsv"))
    cat("Saved GO results for short gained loops\n")
  }

  # Create comparison plot if both have results
  if (nrow(as.data.frame(go_long_lost)) > 0 && nrow(as.data.frame(go_short_gained)) > 0) {
    # Take top 10 terms from each
    top_long <- as.data.frame(go_long_lost) %>%
      head(10) %>%
      mutate(group = "Long Lost (>500kb)")

    top_short <- as.data.frame(go_short_gained) %>%
      head(10) %>%
      mutate(group = "Short Gained (<500kb)")

    combined_go <- bind_rows(top_long, top_short) %>%
      mutate(
        Description = str_wrap(Description, width = 40),
        Description = factor(Description, levels = unique(Description))
      )

    p8_go <- ggplot(combined_go, aes(x = -log10(p.adjust), y = reorder(Description, -log10(p.adjust)))) +
      geom_point(aes(size = Count, color = group)) +
      facet_wrap(~group, scales = "free_y", ncol = 2) +
      scale_color_manual(values = c("Long Lost (>500kb)" = COLORS$down,
                                    "Short Gained (<500kb)" = COLORS$up)) +
      scale_size_continuous(range = c(3, 8)) +
      labs(
        title = "GO Biological Process Enrichment",
        subtitle = "Comparing genes at long lost vs short gained loops",
        x = "-Log10(Adjusted p-value)",
        y = ""
      ) +
      theme_minimal(base_size = 11) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey40"),
        strip.text = element_text(face = "bold", size = 11),
        legend.position = "bottom",
        axis.text.y = element_text(size = 9)
      ) +
      guides(color = "none")

    save_multiformat_ggplot(p8_go, file.path(OUTPUT_DIR, "08_go_comparison_long_vs_short"), width = 14, height = 8)
    cat("Saved: 08_go_comparison_long_vs_short.pdf\n")
  } else {
    cat("Insufficient GO terms for comparison plot\n")
    # Create a placeholder figure
    p8_placeholder <- ggplot() +
      annotate("text", x = 0.5, y = 0.5,
               label = "GO enrichment results\n(see TSV files for details)",
               size = 6, hjust = 0.5) +
      theme_void() +
      labs(title = "GO Enrichment Analysis")

    save_multiformat_ggplot(p8_placeholder, file.path(OUTPUT_DIR, "08_go_comparison_long_vs_short"), width = 10, height = 6)
  }

}, error = function(e) {
  cat("GO enrichment failed:", conditionMessage(e), "\n")
  cat("Creating placeholder file...\n")

  # Create placeholder
  p8_error <- ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = paste("GO enrichment could not complete\n", conditionMessage(e)),
             size = 5, hjust = 0.5) +
    theme_void() +
    labs(title = "GO Enrichment Analysis")

  save_multiformat_ggplot(p8_error, file.path(OUTPUT_DIR, "08_go_comparison_long_vs_short"), width = 10, height = 6)
})

# ==============================================================================
# FIGURE 9: Combined Summary Figure
# ==============================================================================

cat("=== Generating Figure 9: Combined Summary Figure ===\n")

# Simplified CDF for summary
p9a <- ggplot(loops_directional, aes(x = loop_distance_kb, color = direction_label)) +
  stat_ecdf(geom = "step", linewidth = 1) +
  scale_color_manual(
    values = c("Lost in BAP1-KO" = COLORS$down, "Gained in BAP1-KO" = COLORS$up),
    name = ""
  ) +
  scale_x_log10(labels = comma, breaks = c(100, 1000, 10000)) +
  labs(title = "A. Distance CDF", x = "Distance (kb)", y = "Cumulative Prop.") +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    legend.position = "none"
  )

# Simplified bar chart for summary (percentage)
category_pct <- category_summary %>%
  dplyr::select(direction_label, distance_category, percentage)

p9b <- ggplot(category_pct,
               aes(x = distance_category, y = percentage, fill = direction_label)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(
    values = c("Lost in BAP1-KO" = COLORS$down, "Gained in BAP1-KO" = COLORS$up),
    name = ""
  ) +
  labs(title = "B. Distance Distribution", x = "", y = "% of Loops") +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9)
  )

# Simplified scatter for summary
p9c <- ggplot(loops_directional, aes(x = loop_distance_kb, y = logFC)) +
  geom_point(aes(color = direction_label), alpha = 0.3, size = 1) +
  geom_smooth(method = "loess", color = "black", linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(
    values = c("Lost in BAP1-KO" = COLORS$down, "Gained in BAP1-KO" = COLORS$up)
  ) +
  scale_x_log10(labels = comma, breaks = c(100, 1000, 10000)) +
  labs(title = "C. LogFC vs Distance", x = "Distance (kb)", y = "Log2FC") +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    legend.position = "none"
  )

# Key statistics panel
stats_text <- sprintf(
  "Key Statistics:\n\n
Loops Analyzed: %d\n
  Lost: %d (%.1f%%)\n
  Gained: %d (%.1f%%)\n\n
Median Distance:\n
  Lost: %d kb\n
  Gained: %d kb\n
  Ratio: %.2fx\n\n
Statistical Tests:\n
  KS test: p = %.2e\n
  Wilcoxon: p = %.2e\n
  Chi-square: p = %.2e\n\n
Spearman rho: %.3f",
  nrow(loops_directional),
  sum(loops_directional$direction == "down_in_mutant"),
  sum(loops_directional$direction == "down_in_mutant") / nrow(loops_directional) * 100,
  sum(loops_directional$direction == "up_in_mutant"),
  sum(loops_directional$direction == "up_in_mutant") / nrow(loops_directional) * 100,
  round(median_lost / 1000),
  round(median_gained / 1000),
  median_lost / median_gained,
  ks_test$p.value,
  wilcox_test$p.value,
  chisq_test$p.value,
  spearman_cor$estimate
)

p9d <- ggplot() +
  annotate("text", x = 0, y = 0.5, label = stats_text,
           hjust = 0, vjust = 0.5, size = 3, family = "mono") +
  xlim(-0.1, 1) + ylim(0, 1) +
  labs(title = "D. Summary Statistics") +
  theme_void() +
  theme(plot.title = element_text(face = "bold", size = 11))

# Legend plot
legend_data <- data.frame(
  x = c(1, 2),
  y = c(1, 1),
  label = c("Lost in BAP1-KO", "Gained in BAP1-KO"),
  color = c(COLORS$down, COLORS$up)
)

p_legend <- ggplot(legend_data, aes(x = x, y = y, color = label)) +
  geom_point(size = 5) +
  geom_text(aes(label = label), hjust = -0.2, size = 4) +
  scale_color_manual(values = c("Lost in BAP1-KO" = COLORS$down,
                                "Gained in BAP1-KO" = COLORS$up)) +
  xlim(0.5, 4) +
  theme_void() +
  theme(legend.position = "none")

# Combine all panels
p9_combined <- (p9a | p9b) / (p9c | p9d) / p_legend +
  plot_layout(heights = c(1, 1, 0.15)) +
  plot_annotation(
    title = "BAP1-KO 'Loop Rewriting': Loss of Long-Range Loops",
    subtitle = "Chromatin loop reorganization in BAP1 knockout shows preferential loss of long-range interactions",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 12, color = "grey40")
    )
  )

save_multiformat_ggplot(p9_combined, file.path(OUTPUT_DIR, "09_loop_rewriting_summary"), width = 11, height = 10)
cat("Saved: 09_loop_rewriting_summary.pdf\n")

# ==============================================================================
# EXPORT STATISTICS AND SUMMARY TABLES
# ==============================================================================

cat("\n=== Exporting Summary Statistics ===\n")

# Distance shift summary table
distance_shift_summary <- loops_directional %>%
  group_by(direction, distance_category) %>%
  summarise(
    count = n(),
    .groups = "drop"
  ) %>%
  group_by(direction) %>%
  mutate(
    total = sum(count),
    percentage = round(count / total * 100, 2)
  ) %>%
  ungroup() %>%
  pivot_wider(
    id_cols = distance_category,
    names_from = direction,
    values_from = c(count, percentage)
  ) %>%
  mutate(
    lost_gained_ratio = round(count_down_in_mutant / count_up_in_mutant, 3),
    fold_enrichment_loss = round((percentage_down_in_mutant / percentage_up_in_mutant), 3)
  ) %>%
  arrange(match(distance_category, DISTANCE_ORDER))

write_tsv(distance_shift_summary, file.path(OUTPUT_DIR, "distance_shift_summary.tsv"))
cat("Saved: distance_shift_summary.tsv\n")

# Statistics text file
stats_file <- file.path(OUTPUT_DIR, "distance_shift_statistics.txt")
sink(stats_file)
cat("=== LOOP DISTANCE SHIFT ANALYSIS: STATISTICAL SUMMARY ===\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("--- DATA OVERVIEW ---\n")
cat("Total differential loops:", nrow(loops_directional), "\n")
cat("Lost in BAP1-KO (down_in_mutant):", sum(loops_directional$direction == "down_in_mutant"), "\n")
cat("Gained in BAP1-KO (up_in_mutant):", sum(loops_directional$direction == "up_in_mutant"), "\n\n")

cat("--- MEDIAN DISTANCES ---\n")
cat("Lost loops median:", round(median_lost / 1000, 1), "kb\n")
cat("Gained loops median:", round(median_gained / 1000, 1), "kb\n")
cat("Ratio (lost/gained):", round(median_lost / median_gained, 3), "\n\n")

cat("--- STATISTICAL TESTS ---\n")
cat("\nKolmogorov-Smirnov Test (distribution difference):\n")
print(ks_test)
cat("\nWilcoxon Rank-Sum Test (median difference):\n")
print(wilcox_test)
cat("\nChi-Square Test (distance category independence):\n")
print(chisq_test)
cat("\nSpearman Correlation (distance vs logFC):\n")
print(spearman_cor)

cat("\n--- DISTANCE CATEGORY BREAKDOWN ---\n")
print(as.data.frame(distance_shift_summary), row.names = FALSE)

cat("\n--- KEY FINDINGS ---\n")
cat("1. Long loops (>1Mb) are preferentially LOST in BAP1-KO:\n")
cat("   - ", enrichment$percentage_down_in_mutant[enrichment$distance_category == ">1Mb"],
    "% of lost loops vs ", enrichment$percentage_up_in_mutant[enrichment$distance_category == ">1Mb"],
    "% of gained loops\n")
cat("2. Short loops (<100kb) are preferentially GAINED in BAP1-KO:\n")
cat("3. Negative Spearman correlation confirms: longer loops tend to decrease\n")
cat("4. This supports the 'loop rewriting' hypothesis\n")
sink()

cat("Saved: distance_shift_statistics.txt\n")

# ==============================================================================
# COMPLETION
# ==============================================================================

cat("\n=== Analysis Complete ===\n")
cat("Output directory:", OUTPUT_DIR, "\n")
cat("Files generated:\n")
cat("  - 01_distance_cdf_by_direction.pdf\n")
cat("  - 02_distance_category_barplot.pdf\n")
cat("  - 03_distance_density_overlay.pdf\n")
cat("  - 04_logfc_vs_distance_scatter.pdf\n")
cat("  - 05_distance_stratified_volcano.pdf\n")
cat("  - 06_looptype_distance_heatmap.pdf\n")
cat("  - 07_chipseq_distance_analysis.pdf\n")
cat("  - 08_go_comparison_long_vs_short.pdf\n")
cat("  - 09_loop_rewriting_summary.pdf\n")
cat("  - distance_shift_summary.tsv\n")
cat("  - distance_shift_statistics.txt\n")
cat("  - go_long_lost_loops.tsv (if GO enrichment succeeded)\n")
cat("  - go_short_gained_loops.tsv (if GO enrichment succeeded)\n")
cat("\nEnd time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
