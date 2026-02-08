# biomodal/downstream/scripts/viz_sections/section_07_effect_size.R
# Section 7: Effect Size Distributions
# Standalone script - sources shared config for all dependencies and data

# Run from downstream/ directory
source("scripts/viz_sections/_shared_config.R")

# =============================================================================
# SECTION 7: EFFECT SIZE DISTRIBUTIONS
# =============================================================================

cat("================================================================================\n")
cat("SECTION 7: EFFECT SIZE DISTRIBUTIONS\n")
cat("================================================================================\n\n")

cat("Creating effect size distribution plots...\n")

# Histogram of mC differences
p_hist_mc <- ggplot(mc_dmr, aes(x = mod_difference * 100, fill = significant)) +
  geom_histogram(bins = 50, alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_fill_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "grey70"),
                    labels = c("TRUE" = "Significant", "FALSE" = "Not Significant"),
                    name = "") +
  labs(
    title = "5mC Effect Size Distribution",
    subtitle = sprintf("Mean change: %+.2f%% (significant genes)",
                       mean(mc_dmr$mod_difference[mc_dmr$significant]) * 100),
    x = "5mC Change (Mutant - Control, %)", y = "Count"
  ) +
  theme_biomodal()

# Histogram of hmC differences
p_hist_hmc <- ggplot(hmc_dmr, aes(x = mod_difference * 100, fill = significant)) +
  geom_histogram(bins = 50, alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_fill_manual(values = c("TRUE" = "#377EB8", "FALSE" = "grey70"),
                    labels = c("TRUE" = "Significant", "FALSE" = "Not Significant"),
                    name = "") +
  labs(
    title = "5hmC Effect Size Distribution",
    subtitle = sprintf("Mean change: %+.2f%% (significant genes)",
                       mean(hmc_dmr$mod_difference[hmc_dmr$significant]) * 100),
    x = "5hmC Change (Mutant - Control, %)", y = "Count"
  ) +
  theme_biomodal()

# Violin plot comparison
violin_data <- data.frame(
  Type = c(rep("5mC", sum(mc_dmr$significant)), rep("5hmC", sum(hmc_dmr$significant))),
  Change = c(mc_dmr$mod_difference[mc_dmr$significant] * 100,
             hmc_dmr$mod_difference[hmc_dmr$significant] * 100)
)

p_violin <- ggplot(violin_data, aes(x = Type, y = Change, fill = Type)) +
  geom_violin(alpha = 0.7, draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_fill_manual(values = COLORS$methylation) +
  labs(
    title = "Effect Size Comparison: 5mC vs 5hmC",
    subtitle = "5mC shifts positive (hypermethylation) | 5hmC shifts negative (hypomethylation)",
    x = "", y = "Change (Mutant - Control, %)"
  ) +
  theme_biomodal() +
  theme(legend.position = "none")

p_effect_combined <- (p_hist_mc | p_hist_hmc) / p_violin +
  plot_annotation(
    title = "Effect Size Distributions",
    subtitle = "Significant genes show opposite directions: 5mC\u2191 / 5hmC\u2193",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
                  plot.subtitle = element_text(hjust = 0.5))
  )

save_multiformat_ggplot(p_effect_combined, file.path(OUTPUT_DIR, "07_effect_size_distributions"),
                        width = 14, height = 12)

cat("Section 7 complete.\n\n")
