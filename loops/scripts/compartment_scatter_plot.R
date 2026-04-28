# scripts/compartment_scatter_plot.R
# PC1 scatter plot: BAP1-KO mutant vs WT control
#
# Jesse Dixon's suggestion (04/10/26 meeting, slide 29):
#   "Do PC1 mutant vs. PC1 WT scatter, not volcano — a change of +/-1 in PC1
#    doesn't necessarily mean a flip."
#
# This shows quantitative weakening vs full compartment switches directly,
# which a volcano plot cannot capture (it loses absolute PC1 identity).
#
# Usage:
#   Rscript scripts/compartment_scatter_plot.R
#   Rscript scripts/compartment_scatter_plot.R path/to/compartment_annotated.tsv
#   Rscript scripts/compartment_scatter_plot.R --output output/custom_dir/
#   Rscript scripts/compartment_scatter_plot.R --n-labels 20

# =============================================================================
# SETUP
# =============================================================================

library(ggplot2)
library(ggrepel)

source("scripts/utils/multi_format_output.R")

# =============================================================================
# PARSE CLI ARGUMENTS
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)

input_file <- "data/tsvs/figure_1_tads_boundaries_compartments/1D_compartment_all_annotated.tsv"
output_dir <- "output/compartment_analysis"
n_labels <- 15

i <- 1
while (i <= length(args)) {
  if (args[i] == "--output" && i < length(args)) {
    output_dir <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--n-labels" && i < length(args)) {
    n_labels <- as.integer(args[i + 1])
    i <- i + 2
  } else if (!startsWith(args[i], "--")) {
    input_file <- args[i]
    i <- i + 1
  } else {
    i <- i + 1
  }
}

cat("================================================\n")
cat("Compartment PC1 Scatter Plot: Mutant vs WT\n")
cat("================================================\n")
cat(sprintf("  Input:    %s\n", input_file))
cat(sprintf("  Output:   %s\n", output_dir))
cat(sprintf("  Labels:   top %d genes\n", n_labels))

# =============================================================================
# LOAD DATA
# =============================================================================

if (!file.exists(input_file)) {
  stop(sprintf("Input file not found: %s", input_file))
}

df <- read.delim(input_file, header = TRUE, stringsAsFactors = FALSE)
cat(sprintf("\nLoaded %d genomic bins (25kb)\n", nrow(df)))

# Validate required columns
required_cols <- c("ctrl_avg_PC1", "mut_avg_PC1", "Difference",
                   "adj_pvalue", "direction", "Gene_Name")
missing <- setdiff(required_cols, colnames(df))
if (length(missing) > 0) {
  stop(sprintf("Missing columns: %s", paste(missing, collapse = ", ")))
}

# =============================================================================
# CLASSIFY POINTS
# =============================================================================
# Sign convention: Difference = mut_avg_PC1 - ctrl_avg_PC1
#   Positive Difference -> mutant PC1 higher -> shift toward A -> B_to_A -> firebrick3
#   Negative Difference -> mutant PC1 lower  -> shift toward B -> A_to_B -> steelblue

fdr_thresh <- 0.05
diff_thresh <- 0.30

df$category <- "NS"
df$category[df$adj_pvalue < fdr_thresh & df$Difference > diff_thresh] <- "B_to_A"
df$category[df$adj_pvalue < fdr_thresh & df$Difference < -diff_thresh] <- "A_to_B"
df$category <- factor(df$category, levels = c("NS", "A_to_B", "B_to_A"))

n_total <- nrow(df)
n_b2a <- sum(df$category == "B_to_A")
n_a2b <- sum(df$category == "A_to_B")

cat(sprintf("\nSignificant bins (FDR < %.2f, |diff| > %.2f):\n", fdr_thresh, diff_thresh))
cat(sprintf("  B->A (more active in mutant):   %d\n", n_b2a))
cat(sprintf("  A->B (more inactive in mutant): %d\n", n_a2b))

# Count actual quadrant switches (cross the zero line)
a_to_b_switch <- sum(df$ctrl_avg_PC1 > 0 & df$mut_avg_PC1 < 0)
b_to_a_switch <- sum(df$ctrl_avg_PC1 < 0 & df$mut_avg_PC1 > 0)

cat(sprintf("\nTrue compartment switches (cross zero):\n"))
cat(sprintf("  A -> B (ctrl>0, mut<0): %d (%.1f%% of total)\n",
            a_to_b_switch, 100 * a_to_b_switch / n_total))
cat(sprintf("  B -> A (ctrl<0, mut>0): %d (%.1f%% of total)\n",
            b_to_a_switch, 100 * b_to_a_switch / n_total))

# Pearson and Spearman correlation
r_pearson <- cor(df$ctrl_avg_PC1, df$mut_avg_PC1, method = "pearson")
r_spearman <- cor(df$ctrl_avg_PC1, df$mut_avg_PC1, method = "spearman")
cat(sprintf("\nCorrelation (ctrl vs mut PC1):\n"))
cat(sprintf("  Pearson r  = %.4f\n", r_pearson))
cat(sprintf("  Spearman p = %.4f\n", r_spearman))

# =============================================================================
# SELECT GENE LABELS
# =============================================================================

df$label <- ""
sig_df <- df[df$category != "NS", ]
if (nrow(sig_df) > 0 && n_labels > 0) {
  top_genes <- head(sig_df[order(sig_df$adj_pvalue), ], n_labels)
  df$label[rownames(df) %in% rownames(top_genes)] <- df$Gene_Name[rownames(df) %in% rownames(top_genes)]
}

# =============================================================================
# BUILD SCATTER PLOT
# =============================================================================

cat("\nBuilding scatter plot...\n")

# Split data for layering (non-sig background, sig foreground)
df_ns <- df[df$category == "NS", ]
df_sig <- df[df$category != "NS", ]

# Axis range for equal scaling
pc1_range <- range(c(df$ctrl_avg_PC1, df$mut_avg_PC1), na.rm = TRUE)
axis_pad <- 0.15
axis_lim <- c(pc1_range[1] - axis_pad, pc1_range[2] + axis_pad)

p <- ggplot() +
  # Non-significant points (background)
  geom_point(data = df_ns,
             aes(x = ctrl_avg_PC1, y = mut_avg_PC1),
             color = "grey80", alpha = 0.12, size = 0.3) +
  # Significant points (foreground)
  geom_point(data = df_sig,
             aes(x = ctrl_avg_PC1, y = mut_avg_PC1, color = category),
             alpha = 0.5, size = 0.6) +
  # Identity line (y = x: no change)
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey40", linewidth = 0.5) +
  # Quadrant lines at 0 (A/B compartment boundary)
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey60", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey60", linewidth = 0.4) +
  # Gene labels
  geom_text_repel(data = df[df$label != "", ],
                  aes(x = ctrl_avg_PC1, y = mut_avg_PC1, label = label),
                  size = 2.8, max.overlaps = 20,
                  segment.color = "grey30", segment.size = 0.3,
                  min.segment.length = 0.2,
                  box.padding = 0.4, point.padding = 0.2) +
  # Color scale
  scale_color_manual(
    values = c("A_to_B" = "steelblue", "B_to_A" = "firebrick3"),
    labels = c("A_to_B" = sprintf("More Inactive A->B (%d)", n_a2b),
               "B_to_A" = sprintf("More Active B->A (%d)", n_b2a))
  ) +
  # Equal aspect ratio

  coord_equal(xlim = axis_lim, ylim = axis_lim) +
  # Quadrant labels
  annotate("text", x = axis_lim[2] - 0.1, y = axis_lim[2] - 0.1,
           label = "A | A", color = "grey50", size = 3, fontface = "italic",
           hjust = 1, vjust = 1) +
  annotate("text", x = axis_lim[1] + 0.1, y = axis_lim[1] + 0.1,
           label = "B | B", color = "grey50", size = 3, fontface = "italic",
           hjust = 0, vjust = 0) +
  annotate("text", x = axis_lim[2] - 0.1, y = axis_lim[1] + 0.1,
           label = sprintf("A->B\n(%d)", a_to_b_switch), color = "steelblue",
           size = 2.5, fontface = "bold", hjust = 1, vjust = 0) +
  annotate("text", x = axis_lim[1] + 0.1, y = axis_lim[2] - 0.1,
           label = sprintf("B->A\n(%d)", b_to_a_switch), color = "firebrick3",
           size = 2.5, fontface = "bold", hjust = 0, vjust = 1) +
  # Correlation annotation
  annotate("text", x = axis_lim[1] + 0.15, y = axis_lim[2] - 0.4,
           label = sprintf("r = %.4f", r_pearson),
           color = "black", size = 3.5, hjust = 0) +
  # Labels
  labs(
    x = "WT PC1",
    y = "BAP1-KO PC1",
    title = "Compartment PC1: BAP1-KO vs WT",
    subtitle = sprintf("%s bins (25kb) | colored: FDR < %.2f, |diff| > %.2f | top %d genes labeled",
                       format(n_total, big.mark = ","), fdr_thresh, diff_thresh, n_labels)
  ) +
  # Theme (matching project conventions)
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0,
                              margin = margin(b = 2)),
    plot.subtitle = element_text(size = 9, hjust = 0,
                                 margin = margin(b = 10)),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.8),
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 10)),
    axis.text = element_text(size = 10),
    legend.position = "top",
    legend.justification = "center",
    legend.box = "horizontal",
    legend.text = element_text(size = 10),
    legend.title = element_blank(),
    plot.margin = margin(15, 15, 15, 15)
  )

# =============================================================================
# SAVE OUTPUTS
# =============================================================================

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

cat("\nSaving plots...\n")
output_base <- file.path(output_dir, "compartment_pc1_scatter")
save_multiformat_ggplot(p, output_base, width = 10, height = 10)

# Save summary statistics
summary_file <- file.path(output_dir, "compartment_pc1_scatter_summary.txt")
summary_lines <- c(
  "Compartment PC1 Scatter Plot Summary",
  "=====================================",
  sprintf("Input: %s", input_file),
  sprintf("Date: %s", Sys.Date()),
  "",
  sprintf("Total bins: %s", format(n_total, big.mark = ",")),
  sprintf("Significance: FDR < %.2f, |diff| > %.2f", fdr_thresh, diff_thresh),
  sprintf("  B->A (more active in mutant):   %d", n_b2a),
  sprintf("  A->B (more inactive in mutant): %d", n_a2b),
  "",
  "True compartment switches (cross zero):",
  sprintf("  A -> B (ctrl>0, mut<0): %d (%.2f%%)", a_to_b_switch, 100 * a_to_b_switch / n_total),
  sprintf("  B -> A (ctrl<0, mut>0): %d (%.2f%%)", b_to_a_switch, 100 * b_to_a_switch / n_total),
  sprintf("  Total switches: %d (%.2f%%)", a_to_b_switch + b_to_a_switch,
          100 * (a_to_b_switch + b_to_a_switch) / n_total),
  "",
  "Correlation (ctrl vs mut PC1):",
  sprintf("  Pearson r  = %.4f", r_pearson),
  sprintf("  Spearman rho = %.4f", r_spearman),
  "",
  "Interpretation:",
  "  High correlation confirms Jesse's observation: changes are mostly",
  "  quantitative weakening along the diagonal, not full A<->B switches.",
  sprintf("  Only %.1f%% of bins actually cross the zero boundary.",
          100 * (a_to_b_switch + b_to_a_switch) / n_total)
)
writeLines(summary_lines, summary_file)
cat(sprintf("  Saved: %s\n", summary_file))

cat("\nDone.\n")
