# biomodal/downstream/scripts/viz_sections/section_05_coordinated_changes.R
# Section 5: Coordinated Changes Analysis (Key Insight)
# Standalone script - sources shared config for all dependencies and data

# Run from downstream/ directory
source("scripts/viz_sections/_shared_config.R")

# =============================================================================
# SECTION 5: COORDINATED CHANGES ANALYSIS
# =============================================================================

cat("================================================================================\n")
cat("SECTION 5: COORDINATED CHANGES ANALYSIS (KEY INSIGHT)\n")
cat("================================================================================\n\n")

# Merge mC and hmC data by gene
cat("Identifying genes with coordinated mC/hmC changes...\n")

mc_sig <- mc_dmr %>%
  dplyr::filter(significant) %>%
  dplyr::select(gene, mc_diff = mod_difference, mc_q = dmr_qvalue,
                mc_ctrl = mean_mod_group1, mc_mut = mean_mod_group2)

hmc_sig <- hmc_dmr %>%
  dplyr::filter(significant) %>%
  dplyr::select(gene, hmc_diff = mod_difference, hmc_q = dmr_qvalue,
                hmc_ctrl = mean_mod_group1, hmc_mut = mean_mod_group2)

coordinated <- inner_join(mc_sig, hmc_sig, by = "gene")

# Identify genes with opposite directions (mC up, hmC down)
coordinated <- coordinated %>%
  mutate(
    coordinated_pattern = (mc_diff > 0 & hmc_diff < 0),
    combined_effect = abs(mc_diff) + abs(hmc_diff)
  ) %>%
  arrange(desc(combined_effect))

cat(sprintf("  %d genes significant in both mC and hmC\n", nrow(coordinated)))
cat(sprintf("  %d genes show coordinated mC\u2191/hmC\u2193 pattern (%.1f%%)\n",
            sum(coordinated$coordinated_pattern),
            100 * sum(coordinated$coordinated_pattern) / nrow(coordinated)))

# Save coordinated changes table
write.table(coordinated, file.path(TABLES_DIR, "coordinated_changes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved coordinated_changes.tsv\n")

# Scatter plot: mC vs hmC change
cat("Creating mC vs hmC scatter plot...\n")

coordinated$label <- ""
coordinated$label[coordinated$gene %in% KEY_GENES] <- coordinated$gene[coordinated$gene %in% KEY_GENES]

p_scatter <- ggplot(coordinated, aes(x = mc_diff * 100, y = hmc_diff * 100)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(aes(color = coordinated_pattern), alpha = 0.6, size = 2) +
  geom_text_repel(
    aes(label = label),
    size = 4,
    max.overlaps = 15,
    box.padding = 0.5,
    fontface = "bold"
  ) +
  scale_color_manual(
    values = c("TRUE" = "#8E44AD", "FALSE" = "grey60"),
    labels = c("TRUE" = "mC\u2191 / hmC\u2193", "FALSE" = "Other"),
    name = "Pattern"
  ) +
  labs(
    title = "Coordinated 5mC and 5hmC Changes",
    subtitle = sprintf("%d genes significant in both | %.0f%% show mC\u2191/hmC\u2193 pattern",
                       nrow(coordinated), 100 * mean(coordinated$coordinated_pattern)),
    x = "5mC Change (%)",
    y = "5hmC Change (%)"
  ) +
  theme_biomodal() +
  annotate("text", x = 15, y = -12, label = "COORDINATED\nmC\u2191 / hmC\u2193",
           color = "#8E44AD", fontface = "bold", size = 4) +
  annotate("rect", xmin = 0, xmax = 25, ymin = -20, ymax = 0,
           alpha = 0.1, fill = "#8E44AD")

save_multiformat_ggplot(p_scatter, file.path(OUTPUT_DIR, "05a_mc_hmc_scatter"),
                        width = 10, height = 9)

# Top coordinated genes bar chart
cat("Creating top coordinated genes bar chart...\n")

top_coordinated <- coordinated %>%
  filter(coordinated_pattern) %>%
  head(20) %>%
  mutate(gene = factor(gene, levels = rev(gene)))

top_coord_long <- top_coordinated %>%
  dplyr::select(gene, mc_diff, hmc_diff) %>%
  tidyr::pivot_longer(cols = c(mc_diff, hmc_diff), names_to = "type", values_to = "change") %>%
  dplyr::mutate(
    type = ifelse(type == "mc_diff", "5mC", "5hmC"),
    change_pct = change * 100
  )

p_top_coord <- ggplot(top_coord_long, aes(x = gene, y = change_pct, fill = type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(yintercept = 0, color = "black") +
  scale_fill_manual(values = COLORS$methylation, name = "Modification") +
  coord_flip() +
  labs(
    title = "Top 20 Genes with Coordinated mC\u2191/hmC\u2193 Pattern",
    subtitle = "Consistent with impaired TET-mediated demethylation",
    x = "", y = "Change (Mutant - Control, %)"
  ) +
  theme_biomodal()

save_multiformat_ggplot(p_top_coord, file.path(OUTPUT_DIR, "05b_top_coordinated_genes"),
                        width = 10, height = 10)

# Detailed locus plots for key genes
cat("Creating detailed locus plots...\n")

create_locus_plot <- function(gene_name, mc_data, hmc_data) {
  mc_gene <- mc_data %>% filter(gene == gene_name)
  hmc_gene <- hmc_data %>% filter(gene == gene_name)

  if (nrow(mc_gene) == 0 || nrow(hmc_gene) == 0) {
    return(NULL)
  }

  # Create data for plotting
  plot_data <- data.frame(
    Condition = rep(c("Control", "Mutant"), 2),
    Type = c("5mC", "5mC", "5hmC", "5hmC"),
    Methylation = c(mc_gene$mean_mod_group1 * 100, mc_gene$mean_mod_group2 * 100,
                    hmc_gene$mean_mod_group1 * 100, hmc_gene$mean_mod_group2 * 100)
  )

  # Calculate changes
  mc_change <- (mc_gene$mean_mod_group2 - mc_gene$mean_mod_group1) * 100
  hmc_change <- (hmc_gene$mean_mod_group2 - hmc_gene$mean_mod_group1) * 100

  p <- ggplot(plot_data, aes(x = Type, y = Methylation, fill = Condition)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    scale_fill_manual(values = COLORS$condition) +
    labs(
      title = sprintf("%s Gene Body Methylation", gene_name),
      subtitle = sprintf("mC: %+.1f%% | hmC: %+.1f%% | Coordinated pattern", mc_change, hmc_change),
      x = "", y = "Methylation Level (%)"
    ) +
    theme_biomodal() +
    # Add significance annotations
    annotate("text", x = 1, y = max(plot_data$Methylation) * 1.1,
             label = sprintf("q = %.1e", mc_gene$dmr_qvalue), size = 3) +
    annotate("text", x = 2, y = max(plot_data$Methylation) * 1.1,
             label = sprintf("q = %.1e", hmc_gene$dmr_qvalue), size = 3)

  return(p)
}

# Create Syt1 plot (highest effect)
p_syt1 <- create_locus_plot("Syt1", mc_dmr, hmc_dmr)
if (!is.null(p_syt1)) {
  save_multiformat_ggplot(p_syt1, file.path(OUTPUT_DIR, "05c_syt1_detail"),
                          width = 8, height = 6)
}

# Create Zbtb20 plot
p_zbtb20 <- create_locus_plot("Zbtb20", mc_dmr, hmc_dmr)
if (!is.null(p_zbtb20)) {
  save_multiformat_ggplot(p_zbtb20, file.path(OUTPUT_DIR, "05d_zbtb20_detail"),
                          width = 8, height = 6)
}

# Combined key genes panel
key_plots <- list()
for (gene in c("Syt1", "Zbtb20", "Trpm3", "Cntnap2")) {
  p <- create_locus_plot(gene, mc_dmr, hmc_dmr)
  if (!is.null(p)) {
    key_plots[[gene]] <- p
  }
}

if (length(key_plots) == 4) {
  p_key_genes <- (key_plots[[1]] | key_plots[[2]]) / (key_plots[[3]] | key_plots[[4]]) +
    plot_annotation(
      title = "Key Genes: Coordinated Methylation Changes",
      subtitle = "5mC increase + 5hmC decrease pattern in BAP1-KO",
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
                    plot.subtitle = element_text(hjust = 0.5))
    )

  save_multiformat_ggplot(p_key_genes, file.path(OUTPUT_DIR, "05_coordinated_changes"),
                          width = 14, height = 12)
}

cat("Section 5 complete.\n\n")
