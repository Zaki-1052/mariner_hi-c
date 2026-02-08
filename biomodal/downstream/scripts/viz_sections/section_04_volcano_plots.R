# biomodal/downstream/scripts/viz_sections/section_04_volcano_plots.R
# Section 4: Volcano Plots
# Standalone script - sources shared config for all dependencies and data

# Run from downstream/ directory
source("scripts/viz_sections/_shared_config.R")

# =============================================================================
# SECTION 4: VOLCANO PLOTS
# =============================================================================

cat("================================================================================\n")
cat("SECTION 4: VOLCANO PLOTS\n")
cat("================================================================================\n\n")

# Function to create volcano plot
create_volcano <- function(df, title, subtitle, meth_type = "5mC") {
  # Identify genes to label
  df$label <- ""
  df$label[df$gene %in% KEY_GENES & df$significant] <- df$gene[df$gene %in% KEY_GENES & df$significant]

  # Color based on significance and direction
  df$color_group <- case_when(
    !df$significant ~ "Not Significant",
    df$mod_difference > 0 ~ "Hypermethylated",
    TRUE ~ "Hypomethylated"
  )

  # Cap -log10(q) for visualization
  df$neg_log10_q_capped <- pmin(df$neg_log10_q, 300)

  p <- ggplot(df, aes(x = mod_difference * 100, y = neg_log10_q_capped)) +
    geom_point(aes(color = color_group), alpha = 0.6, size = 1.5) +
    geom_hline(yintercept = -log10(Q_THRESHOLD), linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = 0, linetype = "solid", color = "grey40") +
    geom_text_repel(
      aes(label = label),
      size = 3.5,
      max.overlaps = 20,
      box.padding = 0.5,
      segment.color = "grey50"
    ) +
    scale_color_manual(
      values = c("Hypermethylated" = "#D7191C",
                 "Hypomethylated" = "#2C7BB6",
                 "Not Significant" = "grey70"),
      name = "Direction"
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = paste(meth_type, "Difference (Mutant - Control, %)"),
      y = expression(-log[10](q-value))
    ) +
    theme_biomodal() +
    theme(legend.position = "bottom")

  return(p)
}

cat("Creating 5mC volcano plot...\n")
p_volcano_mc <- create_volcano(
  mc_dmr,
  "Gene Body 5mC Differential Methylation",
  sprintf("BAP1-KO vs Control | %d significant (q < %.2f)", sum(mc_dmr$significant), Q_THRESHOLD),
  "5mC"
)

save_multiformat_ggplot(p_volcano_mc, file.path(OUTPUT_DIR, "04a_volcano_mc"),
                        width = 10, height = 8)

cat("Creating 5hmC volcano plot...\n")
p_volcano_hmc <- create_volcano(
  hmc_dmr,
  "Gene Body 5hmC Differential Methylation",
  sprintf("BAP1-KO vs Control | %d significant (q < %.2f)", sum(hmc_dmr$significant), Q_THRESHOLD),
  "5hmC"
)

save_multiformat_ggplot(p_volcano_hmc, file.path(OUTPUT_DIR, "04b_volcano_hmc"),
                        width = 10, height = 8)

# Combined volcano panel
cat("Creating combined volcano panel...\n")
p_volcano_combined <- p_volcano_mc + p_volcano_hmc +
  plot_annotation(
    title = "Differential Methylation: 5mC vs 5hmC",
    subtitle = "Key genes labeled: Syt1, Zbtb20, Trpm3, Cntnap2",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
                  plot.subtitle = element_text(hjust = 0.5))
  )

save_multiformat_ggplot(p_volcano_combined, file.path(OUTPUT_DIR, "04_volcano_plots"),
                        width = 16, height = 8)

cat("Section 4 complete.\n\n")
