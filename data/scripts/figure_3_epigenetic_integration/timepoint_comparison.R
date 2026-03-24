# tads/scripts/timepoint_comparison.R
# Cross-Timepoint Comparison of TADCompare Results
# Visualizes the enrichment direction flip between early and late timepoints

# =============================================================================
# SETUP
# =============================================================================

cat("\n")
cat("========================================\n")
cat("TADCompare Timepoint Comparison Analysis\n")
cat("========================================\n\n")

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(scales)
})

# Load multi-format output utility
source("data/scripts/_shared/multi_format_output.R") # Original: source(file.path(base_dir, "scripts/utils/multi_format_output.R"))

# Define paths
early_file <- "results/early/final/tadcompare_final_filtered.tsv"  # TODO: not in data/
late_file <- "results/late/final/tadcompare_final_filtered.tsv"    # TODO: not in data/
TSV_DIR  <- "data/tsvs/figure_3_epigenetic_integration"   # Original: output_dir <- "results/visualizations/comparison"
PLOT_DIR <- "data/plots/figure_3_epigenetic_integration"   # Original: (same)
output_dir <- PLOT_DIR  # kept for save_multiformat_ggplot calls
dir.create(TSV_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# LOAD DATA
# =============================================================================

cat("Loading data from both timepoints...\n")

early_df <- read.table(early_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
late_df <- read.table(late_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

early_df$timepoint <- "Early (P13)"
late_df$timepoint <- "Late (P60)"

cat(sprintf("  Early: %d boundaries\n", nrow(early_df)))
cat(sprintf("  Late:  %d boundaries\n", nrow(late_df)))

# =============================================================================
# COMPUTE SUMMARY STATISTICS
# =============================================================================

compute_stats <- function(df, tp_name) {
  # Total counts
  n_total <- nrow(df)
  n_diff <- sum(df$Differential == "Differential", na.rm = TRUE)
  n_shifted <- sum(df$Differential == "Shifted", na.rm = TRUE)
  n_all_diff <- n_diff + n_shifted

  # Enrichment direction for differential boundaries
  diff_only <- df %>% filter(Differential %in% c("Differential", "Shifted"))
  n_ctrl <- sum(diff_only$Enriched_In == "Matrix 1", na.rm = TRUE)
  n_mut <- sum(diff_only$Enriched_In == "Matrix 2", na.rm = TRUE)

  # Type breakdown
  type_counts <- table(df$Type)

  data.frame(
    timepoint = tp_name,
    total = n_total,
    differential = n_diff,
    shifted = n_shifted,
    all_differential = n_all_diff,
    pct_differential = 100 * n_all_diff / n_total,
    ctrl_enriched = n_ctrl,
    mut_enriched = n_mut,
    pct_ctrl = 100 * n_ctrl / (n_ctrl + n_mut),
    pct_mut = 100 * n_mut / (n_ctrl + n_mut),
    n_merge = as.numeric(type_counts["Merge"]),
    n_split = as.numeric(type_counts["Split"]),
    n_complex = as.numeric(type_counts["Complex"]),
    n_strength = as.numeric(type_counts["Strength Change"])
  )
}

early_stats <- compute_stats(early_df, "Early (P13)")
late_stats <- compute_stats(late_df, "Late (P60)")
stats_df <- bind_rows(early_stats, late_stats)

cat("\n--- SUMMARY STATISTICS ---\n")
cat(sprintf("\nEarly (P13):\n"))
cat(sprintf("  Total boundaries: %d\n", early_stats$total))
cat(sprintf("  Differential: %d (%.1f%%)\n", early_stats$all_differential, early_stats$pct_differential))
cat(sprintf("  Control-enriched: %d (%.1f%%)\n", early_stats$ctrl_enriched, early_stats$pct_ctrl))
cat(sprintf("  Mutant-enriched:  %d (%.1f%%)\n", early_stats$mut_enriched, early_stats$pct_mut))

cat(sprintf("\nLate (P60):\n"))
cat(sprintf("  Total boundaries: %d\n", late_stats$total))
cat(sprintf("  Differential: %d (%.1f%%)\n", late_stats$all_differential, late_stats$pct_differential))
cat(sprintf("  Control-enriched: %d (%.1f%%)\n", late_stats$ctrl_enriched, late_stats$pct_ctrl))
cat(sprintf("  Mutant-enriched:  %d (%.1f%%)\n", late_stats$mut_enriched, late_stats$pct_mut))

# =============================================================================
# COLOR SCHEMES
# =============================================================================

enrichment_colors <- c(
  "Control" = "#4575b4",      # Blue
  "Mutant" = "#d73027"        # Red
)

timepoint_colors <- c(
  "Early (P13)" = "#7570b3",  # Purple
  "Late (P60)" = "#1b9e77"    # Teal
)

# =============================================================================
# PLOT 1: ENRICHMENT DIRECTION FLIP (Main Finding)
# =============================================================================

cat("\n\nCreating enrichment direction comparison plots...\n")

# Prepare data for stacked bar
enrich_data <- data.frame(
  timepoint = factor(c("Early (P13)", "Early (P13)", "Late (P60)", "Late (P60)"),
                     levels = c("Early (P13)", "Late (P60)")),
  direction = c("Control", "Mutant", "Control", "Mutant"),
  count = c(early_stats$ctrl_enriched, early_stats$mut_enriched,
            late_stats$ctrl_enriched, late_stats$mut_enriched),
  percentage = c(early_stats$pct_ctrl, early_stats$pct_mut,
                 late_stats$pct_ctrl, late_stats$pct_mut)
)

# 1A: Stacked percentage bar chart
p1a <- ggplot(enrich_data, aes(x = timepoint, y = percentage, fill = direction)) +
  geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.5, width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", percentage, count)),
            position = position_stack(vjust = 0.5), size = 4, fontface = "bold") +
  scale_fill_manual(values = enrichment_colors, name = "Enriched In") +
  labs(
    title = "TAD Boundary Enrichment Direction: Early vs Late",
    subtitle = "Direction of differential boundary strength change",
    x = "Timepoint",
    y = "Percentage of Differential Boundaries (%)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 11),
    legend.position = "right",
    legend.title = element_text(size = 11, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 105))

save_multiformat_ggplot(p1a, file.path(output_dir, "3D_enrichment_direction_flip"),
                        width = 8, height = 6)

# 1B: Side-by-side grouped bar chart (raw counts)
p1b <- ggplot(enrich_data, aes(x = timepoint, y = count, fill = direction)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8),
           color = "black", linewidth = 0.5, width = 0.7) +
  geom_text(aes(label = count, group = direction),
            position = position_dodge(width = 0.8), vjust = -0.5, size = 4, fontface = "bold") +
  scale_fill_manual(values = enrichment_colors, name = "Enriched In") +
  labs(
    title = "Differential Boundary Counts by Enrichment Direction",
    subtitle = "Number of boundaries stronger in Control vs Mutant",
    x = "Timepoint",
    y = "Number of Differential Boundaries"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 11),
    legend.position = "right",
    legend.title = element_text(size = 11, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

save_multiformat_ggplot(p1b, file.path(output_dir, "3D_enrichment_direction_counts"),
                        width = 8, height = 6)

# =============================================================================
# PLOT 2: DIVERGING BAR (Net Direction Change)
# =============================================================================

# Net change: positive = more mutant-enriched, negative = more control-enriched
net_data <- data.frame(
  timepoint = factor(c("Early (P13)", "Late (P60)"), levels = c("Early (P13)", "Late (P60)")),
  net_pct = c(early_stats$pct_mut - 50, late_stats$pct_mut - 50),  # Deviation from 50%
  interpretation = c("Boundaries WEAKER\nin mutant", "Boundaries STRONGER\nin mutant")
)

p2 <- ggplot(net_data, aes(x = timepoint, y = net_pct, fill = net_pct > 0)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.5, width = 0.6) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 1) +
  geom_text(aes(label = interpretation, y = net_pct/2),
            size = 4, fontface = "bold", color = "white") +
  geom_text(aes(label = sprintf("%.1f%% mutant-enriched", 50 + net_pct)),
            vjust = ifelse(net_data$net_pct > 0, -0.5, 1.5), size = 3.5) +
  scale_fill_manual(values = c("TRUE" = "#d73027", "FALSE" = "#4575b4"), guide = "none") +
  labs(
    title = "Net Direction of TAD Boundary Changes",
    subtitle = "Deviation from 50% (balanced) | Positive = mutant-enriched majority",
    x = "Timepoint",
    y = "Net Direction (% above/below 50%)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 11),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  ) +
  scale_y_continuous(limits = c(-20, 20), breaks = seq(-20, 20, 10))

save_multiformat_ggplot(p2, file.path(output_dir, "3D_net_direction_diverging"),
                        width = 8, height = 6)

# =============================================================================
# PLOT 3: BOUNDARY TYPE COMPARISON
# =============================================================================

cat("Creating boundary type comparison plots...\n")

# Combine for comparison
type_data <- bind_rows(
  early_df %>% mutate(timepoint = "Early (P13)"),
  late_df %>% mutate(timepoint = "Late (P60)")
) %>%
  filter(Differential %in% c("Differential", "Shifted")) %>%
  count(timepoint, Type) %>%
  group_by(timepoint) %>%
  mutate(percentage = 100 * n / sum(n)) %>%
  ungroup() %>%
  mutate(timepoint = factor(timepoint, levels = c("Early (P13)", "Late (P60)")))

# Filter to key types
key_types <- c("Merge", "Split", "Strength Change", "Complex", "Shifted")
type_data_filtered <- type_data %>% filter(Type %in% key_types)

type_colors <- c(
  "Strength Change" = "#ff7f00",   # Orange
  "Complex" = "#984ea3",           # Purple
  "Shifted" = "#e41a1c",           # Red
  "Merge" = "#377eb8",             # Blue
  "Split" = "#4daf4a"              # Green
)

p3 <- ggplot(type_data_filtered, aes(x = Type, y = percentage, fill = timepoint)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8),
           color = "black", linewidth = 0.3, width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", percentage), group = timepoint),
            position = position_dodge(width = 0.8), vjust = -0.5, size = 3) +
  scale_fill_manual(values = timepoint_colors, name = "Timepoint") +
  labs(
    title = "Differential Boundary Types: Early vs Late",
    subtitle = "Distribution of boundary change mechanisms",
    x = "Boundary Type",
    y = "Percentage of Differential Boundaries (%)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

save_multiformat_ggplot(p3, file.path(output_dir, "3D_boundary_type_comparison"),
                        width = 10, height = 6)

# =============================================================================
# PLOT 4: COMBINED SUMMARY FIGURE (Multi-panel)
# =============================================================================

cat("Creating combined summary figure...\n")

# Panel A: Overview comparison
overview_data <- data.frame(
  timepoint = factor(c("Early (P13)", "Late (P60)"), levels = c("Early (P13)", "Late (P60)")),
  total = c(early_stats$total, late_stats$total),
  differential = c(early_stats$all_differential, late_stats$all_differential),
  pct_diff = c(early_stats$pct_differential, late_stats$pct_differential)
)

p4a <- ggplot(overview_data, aes(x = timepoint, y = pct_diff, fill = timepoint)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.5, width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%\n(%d/%d)", pct_diff, differential, total)),
            vjust = -0.3, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = timepoint_colors, guide = "none") +
  labs(
    title = "A. % Differential Boundaries",
    x = NULL, y = "Percent (%)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2)), limits = c(0, 25))

# Panel B: Enrichment direction (simplified)
p4b <- ggplot(enrich_data, aes(x = timepoint, y = percentage, fill = direction)) +
  geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.3, width = 0.6) +
  geom_text(aes(label = sprintf("%.0f%%", percentage)),
            position = position_stack(vjust = 0.5), size = 3.5, fontface = "bold", color = "white") +
  scale_fill_manual(values = enrichment_colors, name = "Enriched In") +
  labs(
    title = "B. Enrichment Direction",
    x = NULL, y = "Percent (%)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9),
    legend.position = "bottom",
    legend.title = element_text(size = 9),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  ) +
  scale_y_continuous(expand = c(0, 0))

# Panel C: Merge/Split comparison (dramatic changes)
merge_split_data <- data.frame(
  timepoint = factor(rep(c("Early (P13)", "Late (P60)"), each = 2),
                     levels = c("Early (P13)", "Late (P60)")),
  type = rep(c("Merge", "Split"), 2),
  count = c(early_stats$n_merge, early_stats$n_split,
            late_stats$n_merge, late_stats$n_split)
)

p4c <- ggplot(merge_split_data, aes(x = timepoint, y = count, fill = type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7),
           color = "black", linewidth = 0.3, width = 0.6) +
  geom_text(aes(label = count, group = type),
            position = position_dodge(width = 0.7), vjust = -0.3, size = 3) +
  scale_fill_manual(values = c("Merge" = "#377eb8", "Split" = "#4daf4a"), name = "Type") +
  labs(
    title = "C. Merge vs Split Events",
    x = NULL, y = "Count"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9),
    legend.position = "bottom",
    legend.title = element_text(size = 9),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2)))

# Combine panels
combined <- (p4a | p4b | p4c) +
  plot_annotation(
    title = "TADCompare Cross-Timepoint Comparison: BAP1-KO vs Control",
    subtitle = sprintf("Early: %.1f%% control-enriched | Late: %.1f%% mutant-enriched",
                       early_stats$pct_ctrl, late_stats$pct_mut),
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray30")
    )
  )

save_multiformat_ggplot(combined, file.path(output_dir, "3D_combined_comparison_summary"),
                        width = 14, height = 6)

# =============================================================================
# PLOT 5: ENRICHMENT DIRECTION COMPARISON (Corrected)
# =============================================================================

cat("Creating enrichment direction comparison figure...\n")

# Data with exact values (no rounding)
direction_data <- data.frame(
  timepoint = factor(rep(c("Early (P13)", "Late (P60)"), each = 2),
                     levels = c("Early (P13)", "Late (P60)")),
  direction = factor(rep(c("Control-enriched\n(stronger in ctrl)",
                           "Mutant-enriched\n(stronger in mut)"), 2),
                     levels = c("Control-enriched\n(stronger in ctrl)",
                                "Mutant-enriched\n(stronger in mut)")),
  count = c(early_stats$ctrl_enriched, early_stats$mut_enriched,
            late_stats$ctrl_enriched, late_stats$mut_enriched),
  pct = c(early_stats$pct_ctrl, early_stats$pct_mut,
          late_stats$pct_ctrl, late_stats$pct_mut)
)

p5 <- ggplot(direction_data, aes(x = timepoint, y = pct, fill = direction)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8),
           color = "black", linewidth = 0.5, width = 0.7) +
  # Exact percentages with counts
geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", pct, count), group = direction),
            position = position_dodge(width = 0.8), vjust = -0.3, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = c("Control-enriched\n(stronger in ctrl)" = "#4575b4",
                               "Mutant-enriched\n(stronger in mut)" = "#d73027"),
                    name = "Enrichment Direction") +
  geom_hline(yintercept = 50, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  annotate("text", x = 0.55, y = 52, label = "50%", size = 3, color = "gray40") +
  labs(
    title = "TAD Boundary Strength: Control vs Mutant Enrichment",
    subtitle = sprintf("Early (P13): n=%d differential | Late (P60): n=%d differential",
                       early_stats$all_differential, late_stats$all_differential),
    x = "Timepoint",
    y = "Percentage of Differential Boundaries (%)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 11),
    legend.position = "bottom",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  ) +
  scale_y_continuous(limits = c(0, 75), expand = c(0, 0))

save_multiformat_ggplot(p5, file.path(output_dir, "3D_enrichment_direction_comparison"),
                        width = 10, height = 7)

# =============================================================================
# SAVE SUMMARY STATISTICS
# =============================================================================

cat("\nSaving summary statistics...\n")

# Write detailed comparison table
comparison_table <- data.frame(
  Metric = c(
    "Total boundaries",
    "Differential boundaries",
    "% Differential",
    "Shifted boundaries",
    "Control-enriched",
    "Mutant-enriched",
    "% Control-enriched",
    "% Mutant-enriched",
    "Merge events",
    "Split events",
    "Strength Change",
    "Complex"
  ),
  Early_P13 = c(
    early_stats$total,
    early_stats$all_differential,
    sprintf("%.1f%%", early_stats$pct_differential),
    early_stats$shifted,
    early_stats$ctrl_enriched,
    early_stats$mut_enriched,
    sprintf("%.1f%%", early_stats$pct_ctrl),
    sprintf("%.1f%%", early_stats$pct_mut),
    early_stats$n_merge,
    early_stats$n_split,
    early_stats$n_strength,
    early_stats$n_complex
  ),
  Late_P60 = c(
    late_stats$total,
    late_stats$all_differential,
    sprintf("%.1f%%", late_stats$pct_differential),
    late_stats$shifted,
    late_stats$ctrl_enriched,
    late_stats$mut_enriched,
    sprintf("%.1f%%", late_stats$pct_ctrl),
    sprintf("%.1f%%", late_stats$pct_mut),
    late_stats$n_merge,
    late_stats$n_split,
    late_stats$n_strength,
    late_stats$n_complex
  )
)

write.table(comparison_table, file.path(TSV_DIR, "3D_timepoint_comparison_stats.tsv"),  # Original: file.path(output_dir, "timepoint_comparison_stats.tsv")
            sep = "\t", quote = FALSE, row.names = FALSE)

# =============================================================================
# FINAL SUMMARY
# =============================================================================

cat("\n")
cat("========================================\n")
cat("TIMEPOINT COMPARISON COMPLETE\n")
cat("========================================\n\n")

cat("KEY FINDING: Enrichment Direction Flip\n")
cat("----------------------------------------\n")
cat(sprintf("Early (P13): %.1f%% Control-enriched -> Boundaries WEAKER in mutant\n", early_stats$pct_ctrl))
cat(sprintf("Late (P60):  %.1f%% Mutant-enriched  -> Boundaries STRONGER in mutant\n", late_stats$pct_mut))
cat("\n")
cat("INTERPRETATION:\n")
cat("  - Early: Initial destabilization of TAD boundaries in BAP1-KO\n")
cat("  - Late:  Chromatin reorganization/compensation with strengthened boundaries\n")
cat("\n")
cat(sprintf("Output directory: %s\n", output_dir))
cat("\nPlots generated:\n")
cat("  1. enrichment_direction_flip.pdf    - Main finding (stacked bar)\n")
cat("  2. enrichment_direction_counts.pdf  - Raw counts comparison\n")
cat("  3. net_direction_diverging.pdf      - Diverging bar (deviation from 50%)\n")
cat("  4. boundary_type_comparison.pdf     - Type distribution comparison\n")
cat("  5. combined_comparison_summary.pdf  - Multi-panel summary\n")
cat("  6. biological_model_interpretation.pdf - Temporal model figure\n")
cat("  7. timepoint_comparison_stats.tsv   - Statistics table\n")
cat("\n")
