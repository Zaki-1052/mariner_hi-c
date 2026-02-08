# biomodal/downstream/scripts/viz_sections/section_06_top_genes.R
# Section 6: Top Differentially Methylated Genes
# Standalone script - sources shared config for all dependencies and data

# Run from downstream/ directory
source("scripts/viz_sections/_shared_config.R")

# =============================================================================
# SECTION 6: TOP DIFFERENTIALLY METHYLATED GENES
# =============================================================================

cat("================================================================================\n")
cat("SECTION 6: TOP DIFFERENTIALLY METHYLATED GENES\n")
cat("================================================================================\n\n")

# Top 20 mC DMRs
cat("Creating top mC DMRs bar chart...\n")

top_mc <- mc_dmr %>%
  dplyr::filter(significant) %>%
  dplyr::arrange(dmr_qvalue) %>%
  dplyr::distinct(gene, .keep_all = TRUE) %>%
  head(20) %>%
  dplyr::mutate(
    gene = factor(gene, levels = rev(unique(gene))),
    direction = ifelse(mod_difference > 0, "Hypermethylated", "Hypomethylated")
  )

p_top_mc <- ggplot(top_mc, aes(x = gene, y = mod_difference * 100, fill = direction)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_hline(yintercept = 0, color = "black") +
  scale_fill_manual(values = COLORS$direction, name = "Direction") +
  coord_flip() +
  labs(
    title = "Top 20 Gene Body 5mC DMRs",
    subtitle = "Ranked by q-value | BAP1-KO vs Control",
    x = "", y = "5mC Change (%)"
  ) +
  theme_biomodal()

# Top 20 hmC DMRs
cat("Creating top hmC DMRs bar chart...\n")

top_hmc <- hmc_dmr %>%
  dplyr::filter(significant) %>%
  dplyr::arrange(dmr_qvalue) %>%
  dplyr::distinct(gene, .keep_all = TRUE) %>%
  head(20) %>%
  dplyr::mutate(
    gene = factor(gene, levels = rev(unique(gene))),
    direction = ifelse(mod_difference > 0, "Hypermethylated", "Hypomethylated")
  )

p_top_hmc <- ggplot(top_hmc, aes(x = gene, y = mod_difference * 100, fill = direction)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_hline(yintercept = 0, color = "black") +
  scale_fill_manual(values = COLORS$direction, name = "Direction") +
  coord_flip() +
  labs(
    title = "Top 20 Gene Body 5hmC DMRs",
    subtitle = "Ranked by q-value | BAP1-KO vs Control",
    x = "", y = "5hmC Change (%)"
  ) +
  theme_biomodal()

p_top_combined <- p_top_mc | p_top_hmc
save_multiformat_ggplot(p_top_combined, file.path(OUTPUT_DIR, "06a_top_dmrs"),
                        width = 14, height = 10)

# Save top gene lists
write.table(top_mc, file.path(TABLES_DIR, "top_mc_dmrs.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(top_hmc, file.path(TABLES_DIR, "top_hmc_dmrs.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved top_mc_dmrs.tsv and top_hmc_dmrs.tsv\n")

# Venn diagram of overlap
cat("Creating overlap Venn diagram...\n")

mc_genes <- mc_dmr$gene[mc_dmr$significant]
hmc_genes <- hmc_dmr$gene[hmc_dmr$significant]

venn_list <- list(
  "5mC Significant" = mc_genes,
  "5hmC Significant" = hmc_genes
)

p_venn <- ggVennDiagram(venn_list, label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "#4DBBD5") +
  labs(
    title = "Overlap of Significant Genes",
    subtitle = sprintf("5mC: %d | 5hmC: %d | Both: %d",
                       length(mc_genes), length(hmc_genes), length(intersect(mc_genes, hmc_genes)))
  ) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))

save_multiformat_ggplot(p_venn, file.path(OUTPUT_DIR, "06b_venn_overlap"),
                        width = 8, height = 7)

p_top_all <- (p_top_mc | p_top_hmc) / p_venn +
  plot_layout(heights = c(2, 1)) +
  plot_annotation(
    title = "Top Differentially Methylated Genes",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
  )

save_multiformat_ggplot(p_top_all, file.path(OUTPUT_DIR, "06_top_genes"),
                        width = 14, height = 14)

cat("Section 6 complete.\n\n")
