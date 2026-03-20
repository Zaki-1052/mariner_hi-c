#!/usr/bin/env Rscript
# scripts/concordance_pie_chart.R
# Concordance Pie Charts for ABC-RNA DEG Analysis
# Author: Zakir Alibhai
# Date: 2026-03-14
#
# Purpose:
#   Generate pie charts showing concordant vs discordant proportions
#   and the 4 concordance categories for 957 dysregulated genes.
#
#   Paper panel: Figure 4B
#
# Output:
#   - Pie 1: Concordant vs Discordant (binary)
#   - Pie 2: 4 categories (concordant-up, concordant-down,
#            discordant ABC-up/RNA-down, discordant ABC-down/RNA-up)
#   - Combined dual-panel figure
#   - Summary statistics text file
#
# Usage:
#   Rscript scripts/concordance_pie_chart.R
#
# =============================================================================
# CONFIGURATION
# =============================================================================

cat("\n")
cat("================================================\n")
cat("Concordance Pie Charts (Figure 4B)\n")
cat("================================================\n\n")

INPUT_FILE <- "abc/results/discordant_gene_characteristics.tsv"
OUTPUT_DIR <- "abc/results/figures/concordance_pie"

if (!file.exists(INPUT_FILE)) {
  stop(sprintf("Input file not found: %s", INPUT_FILE))
}

# =============================================================================
# LOAD PACKAGES
# =============================================================================

cat("Loading packages...\n")
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

source("scripts/utils/multi_format_output.R")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# LOAD AND CLASSIFY DATA
# =============================================================================

cat("Loading data...\n")
genes <- read.delim(INPUT_FILE, stringsAsFactors = FALSE)

required_cols <- c("concordance", "max_delta_abc", "log2FC")
if (!all(required_cols %in% names(genes))) {
  stop(sprintf("Missing columns: %s",
               paste(setdiff(required_cols, names(genes)),
                     collapse = ", ")))
}

cat(sprintf("  Loaded %d dysregulated genes\n", nrow(genes)))

# 4-way classification
genes$category <- case_when(
  genes$concordance == "Concordant" & genes$log2FC > 0 ~
    "Gained enhancer +\nupregulated gene",
  genes$concordance == "Concordant" & genes$log2FC < 0 ~
    "Lost enhancer +\ndownregulated gene",
  genes$concordance == "Discordant" & genes$max_delta_abc > 0 ~
    "Gained enhancer +\ndownregulated gene",
  genes$concordance == "Discordant" & genes$max_delta_abc < 0 ~
    "Lost enhancer +\nupregulated gene",
  TRUE ~ "Other"
)

# Validate no "Other"
n_other <- sum(genes$category == "Other")
if (n_other > 0) {
  warning(sprintf("%d genes fell into 'Other' category", n_other))
}

# Tally binary
binary_df <- genes %>%
  count(concordance) %>%
  mutate(pct = n / sum(n) * 100,
         label = sprintf("%s\n%d (%.1f%%)", concordance, n, pct))

# Tally 4-category
cat4_df <- genes %>%
  count(category) %>%
  mutate(pct = n / sum(n) * 100,
         label = sprintf("%s\n%d (%.1f%%)", category, n, pct))

# Print summary
cat("\n  Binary concordance:\n")
for (i in seq_len(nrow(binary_df))) {
  cat(sprintf("    %s: %d (%.1f%%)\n",
              binary_df$concordance[i],
              binary_df$n[i],
              binary_df$pct[i]))
}
cat("\n  4-category breakdown:\n")
for (i in seq_len(nrow(cat4_df))) {
  cat(sprintf("    %s: %d (%.1f%%)\n",
              gsub("\n", " ", cat4_df$category[i]),
              cat4_df$n[i],
              cat4_df$pct[i]))
}

# =============================================================================
# PIE 1: CONCORDANT VS DISCORDANT
# =============================================================================

cat("\nGenerating pie charts...\n")

binary_colors <- c(
  "Concordant" = "#4DAF4A",
  "Discordant" = "#E41A1C"
)

binary_df$bar_label <- sprintf("%d\n(%.1f%%)", binary_df$n, binary_df$pct)

p_binary <- ggplot(binary_df,
                   aes(x = concordance, y = n, fill = concordance)) +
  geom_col(width = 0.6, color = "white", linewidth = 0.5) +
  geom_text(aes(label = bar_label),
            vjust = -0.3, size = 5, fontface = "bold") +
  scale_fill_manual(values = binary_colors, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Dysregulated Genes",
       subtitle = sprintf(
         "%d genes (|log2FC|>0.5, padj<0.05, |dABC|>0.01)",
         sum(binary_df$n)),
       x = NULL, y = "Number of Genes") +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(
      face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(
      size = 11, hjust = 0.5, color = "grey40"),
    legend.position = "none",
    panel.grid.major.x = element_blank()
  )

save_multiformat_ggplot(
  p_binary,
  file.path(OUTPUT_DIR, "concordance_pie_binary"),
  width = 6, height = 6
)

# =============================================================================
# PIE 2: 4 CATEGORIES
# =============================================================================

# Order categories: concordant first, then discordant
cat_order <- c(
  "Gained enhancer +\nupregulated gene",
  "Lost enhancer +\ndownregulated gene",
  "Gained enhancer +\ndownregulated gene",
  "Lost enhancer +\nupregulated gene"
)
cat4_df$category <- factor(cat4_df$category, levels = cat_order)

cat4_colors <- c(
  "Gained enhancer +\nupregulated gene"   = "#D6604D",
  "Lost enhancer +\ndownregulated gene"    = "#4393C3",
  "Gained enhancer +\ndownregulated gene"  = "#FDB863",
  "Lost enhancer +\nupregulated gene"      = "#B2ABD2"
)

p_4cat <- ggplot(cat4_df,
                 aes(x = "", y = n, fill = category)) +
  geom_bar(stat = "identity", width = 1,
           color = "white", linewidth = 1) +
  coord_polar("y", start = 0) +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5),
            size = 3.5, fontface = "bold") +
  scale_fill_manual(values = cat4_colors,
                    name = NULL) +
  labs(title = "Concordance Categories",
       subtitle = "Direction of ABC and RNA changes") +
  theme_void(base_size = 14) +
  theme(
    plot.title = element_text(
      face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(
      size = 12, hjust = 0.5, color = "grey40"),
    legend.position = "none"
  )

save_multiformat_ggplot(
  p_4cat,
  file.path(OUTPUT_DIR, "concordance_pie_4cat"),
  width = 7, height = 6
)

# =============================================================================
# COMBINED DUAL PIE
# =============================================================================

p_combined <- p_binary + p_4cat +
  plot_layout(widths = c(1, 1.3)) +
  plot_annotation(
    title = "ABC-RNA Concordance Analysis",
    subtitle = sprintf(
      "%d Dysregulated Genes (|log2FC|>0.5, padj<0.05, |dABC|>0.01)",
      sum(binary_df$n)),
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(
        face = "bold", size = 18, hjust = 0.5),
      plot.subtitle = element_text(
        size = 13, hjust = 0.5, color = "grey40")
    )
  )

save_multiformat_ggplot(
  p_combined,
  file.path(OUTPUT_DIR, "concordance_pie_combined"),
  width = 13, height = 6
)

# =============================================================================
# SUMMARY
# =============================================================================

summary_lines <- c(
  "ABC-RNA Concordance Pie Chart Summary",
  "=====================================",
  "",
  sprintf("Analysis Date: %s", Sys.Date()),
  sprintf("Input: %s", INPUT_FILE),
  sprintf("Total dysregulated genes: %d", nrow(genes)),
  "",
  "Binary Concordance:",
  sprintf("  Concordant: %d (%.1f%%)",
          binary_df$n[binary_df$concordance == "Concordant"],
          binary_df$pct[binary_df$concordance == "Concordant"]),
  sprintf("  Discordant: %d (%.1f%%)",
          binary_df$n[binary_df$concordance == "Discordant"],
          binary_df$pct[binary_df$concordance == "Discordant"]),
  "",
  "4-Category Breakdown:"
)

for (i in seq_len(nrow(cat4_df))) {
  summary_lines <- c(summary_lines,
    sprintf("  %s: %d (%.1f%%)",
            gsub("\n", " ", as.character(cat4_df$category[i])),
            cat4_df$n[i], cat4_df$pct[i]))
}

summary_file <- file.path(OUTPUT_DIR,
                          "concordance_pie_summary.txt")
writeLines(summary_lines, summary_file)
cat(sprintf("  Saved: %s\n", summary_file))

# =============================================================================
# DONE
# =============================================================================

cat("\n================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("================================================\n\n")

cat(sprintf("Output: %s\n", OUTPUT_DIR))
cat("Files:\n")
cat("  - concordance_pie_binary/{pdf,svg,jpg}\n")
cat("  - concordance_pie_4cat/{pdf,svg,jpg}\n")
cat("  - concordance_pie_combined/{pdf,svg,jpg}\n")
cat("  - concordance_pie_summary.txt\n\n")
