# biomodal/downstream/scripts/viz_sections/section_03_dmr_statistics.R
# Section 3: DMR Statistics by Genomic Region
# Standalone script - sources shared config for all dependencies and data

# Run from downstream/ directory
source("scripts/viz_sections/_shared_config.R")

# =============================================================================
# SECTION 3: DMR STATISTICS BY GENOMIC REGION
# =============================================================================

cat("================================================================================\n")
cat("SECTION 3: DMR STATISTICS BY GENOMIC REGION\n")
cat("================================================================================\n\n")

# Create region comparison data frame
region_stats <- data.frame(
  Region = character(),
  Total = integer(),
  Significant = integer(),
  Percentage = numeric(),
  stringsAsFactors = FALSE
)

for (region in names(region_dmrs)) {
  if (!is.null(region_dmrs[[region]])) {
    region_stats <- rbind(region_stats, data.frame(
      Region = region,
      Total = nrow(region_dmrs[[region]]),
      Significant = sum(region_dmrs[[region]]$significant),
      Percentage = 100 * sum(region_dmrs[[region]]$significant) / nrow(region_dmrs[[region]])
    ))
  }
}

# Order by significant count
region_stats <- region_stats %>%
  mutate(Region = factor(Region, levels = Region[order(Significant, decreasing = TRUE)]))

cat("Creating region comparison bar chart...\n")

p_region <- ggplot(region_stats, aes(x = Region, y = Significant, fill = Region)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", Significant, Percentage)),
            vjust = -0.3, size = 3.5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Significant mC DMRs by Genomic Region (CG Context)",
    subtitle = sprintf("q-value < %.2f | Gene bodies dominate differential methylation", Q_THRESHOLD),
    x = "", y = "Number of Significant DMRs"
  ) +
  theme_biomodal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

save_multiformat_ggplot(p_region, file.path(OUTPUT_DIR, "03a_dmr_by_region"),
                        width = 10, height = 7)

# mC vs hmC comparison for gene bodies
cat("Creating mC vs hmC comparison...\n")

mc_hmc_compare <- data.frame(
  Type = c("5mC", "5hmC"),
  Total = c(nrow(mc_dmr), nrow(hmc_dmr)),
  Significant = c(sum(mc_dmr$significant), sum(hmc_dmr$significant)),
  Percentage = c(100 * sum(mc_dmr$significant) / nrow(mc_dmr),
                 100 * sum(hmc_dmr$significant) / nrow(hmc_dmr))
)

p_mc_hmc <- ggplot(mc_hmc_compare, aes(x = Type, y = Significant, fill = Type)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", Significant, Percentage)),
            vjust = -0.3, size = 4) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  scale_fill_manual(values = COLORS$methylation) +
  labs(
    title = "Gene Body DMRs: 5mC vs 5hmC",
    subtitle = "5hmC shows more differential genes than 5mC",
    x = "", y = "Number of Significant DMRs"
  ) +
  theme_biomodal() +
  theme(legend.position = "none")

# Context comparison (CG vs CHG vs CHH)
cat("Creating methylation context comparison...\n")

context_data <- data.frame(
  Context = c("CG (CpG)", "CHG", "CHH"),
  Baseline_Methylation = c(72, 0.65, 0.87),  # From analysis_summary
  Significant_DMRs = c(sum(mc_dmr$significant), 0, 0),  # CHG/CHH have 0 significant
  Label = c("Primary Signal", "No Signal", "No Signal")
)

p_context_meth <- ggplot(context_data, aes(x = Context, y = Baseline_Methylation, fill = Context)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%", Baseline_Methylation)), vjust = -0.3, size = 4) +
  scale_fill_manual(values = c("CG (CpG)" = "#E41A1C", "CHG" = "#999999", "CHH" = "#CCCCCC")) +
  labs(
    title = "Baseline Methylation by Context",
    subtitle = "CpG methylation is ~100x more abundant than non-CpG",
    x = "Methylation Context", y = "Methylation Level (%)"
  ) +
  theme_biomodal() +
  theme(legend.position = "none")

p_context_dmr <- ggplot(context_data, aes(x = Context, y = Significant_DMRs, fill = Context)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = Label), vjust = -0.3, size = 4) +
  scale_fill_manual(values = c("CG (CpG)" = "#E41A1C", "CHG" = "#999999", "CHH" = "#CCCCCC")) +
  labs(
    title = "Significant DMRs by Context",
    subtitle = "No significant changes in non-CpG methylation",
    x = "Methylation Context", y = "Significant Gene Body DMRs"
  ) +
  theme_biomodal() +
  theme(legend.position = "none")

# Combine into panel
p_region_combined <- (p_region | p_mc_hmc) / (p_context_meth | p_context_dmr) +
  plot_annotation(
    title = "DMR Statistics Summary",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
  )

save_multiformat_ggplot(p_region_combined, file.path(OUTPUT_DIR, "03_dmr_region_statistics"),
                        width = 14, height = 12)

# -------------------------------------------------------------------------
# FIGURE 3b: Direction comparison - mC vs hmC asymmetry
# -------------------------------------------------------------------------
cat("Creating Figure 3b: Methylation direction comparison...\n")

# Calculate direction statistics for mC
mc_sig_stats <- mc_dmr %>%
  dplyr::filter(significant) %>%
  dplyr::summarise(
    Hypermethylated = sum(direction == "Hypermethylated"),
    Hypomethylated = sum(direction == "Hypomethylated"),
    Total = n()
  ) %>%
 dplyr::mutate(
    Hyper_pct = 100 * Hypermethylated / Total,
    Hypo_pct = 100 * Hypomethylated / Total
  )

# Calculate direction statistics for hmC
hmc_sig_stats <- hmc_dmr %>%
  dplyr::filter(significant) %>%
  dplyr::summarise(
    Increased = sum(direction == "Hypermethylated"),
    Decreased = sum(direction == "Hypomethylated"),
    Total = n()
  ) %>%
  dplyr::mutate(
    Inc_pct = 100 * Increased / Total,
    Dec_pct = 100 * Decreased / Total
  )

# Create data frame for plotting
direction_data <- data.frame(
  Modification = c("5mC", "5mC", "5hmC", "5hmC"),
  Direction = c("Increased", "Decreased", "Increased", "Decreased"),
  Count = c(mc_sig_stats$Hypermethylated, mc_sig_stats$Hypomethylated,
            hmc_sig_stats$Increased, hmc_sig_stats$Decreased),
  Percentage = c(mc_sig_stats$Hyper_pct, mc_sig_stats$Hypo_pct,
                 hmc_sig_stats$Inc_pct, hmc_sig_stats$Dec_pct)
)

direction_data$Modification <- factor(direction_data$Modification, levels = c("5mC", "5hmC"))
direction_data$Direction <- factor(direction_data$Direction, levels = c("Increased", "Decreased"))

# Create the bar chart
p_direction <- ggplot(direction_data, aes(x = Modification, y = Percentage, fill = Direction)) +
 geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.5, width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", Percentage, Count)),
            position = position_dodge(width = 0.7), vjust = -0.3, size = 4, fontface = "bold") +
  scale_fill_manual(values = c("Increased" = "#D7191C", "Decreased" = "#2C7BB6"),
                    name = "Direction in\nBAP1-KO Mutant") +
  scale_y_continuous(limits = c(0, 105), expand = c(0, 0)) +
  labs(
    title = "Asymmetric Direction of Methylation Changes",
    subtitle = sprintf("5mC: %d significant genes | 5hmC: %d significant genes",
                       mc_sig_stats$Total, hmc_sig_stats$Total),
    x = "Methylation Type",
    y = "Percentage of Significant DMRs (%)"
  ) +
  theme_biomodal(base_size = 14) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16)
  ) +
  # Add annotation for key finding
  annotate("text", x = 1.5, y = 98,
           label = "mC increases while hmC decreases\n= Blocked TET-mediated demethylation",
           size = 4, fontface = "italic", color = "grey30")

save_multiformat_ggplot(p_direction, file.path(OUTPUT_DIR, "03b_direction_comparison"),
                        width = 10, height = 8)

cat("Section 3 complete.\n\n")
