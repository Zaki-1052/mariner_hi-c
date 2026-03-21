# biomodal/downstream/scripts/viz_sections/section_32_chg_exploratory_analysis.R
# Section 32: CHG Context Exploratory Analysis
# Characterize 70 significant CHG mC DMRs and 7 hmC DMRs from run-3 (deep-seq)
# Standalone script - sources shared config for all dependencies and data

# Run from downstream/ directory
source("scripts/viz_sections/_shared_config.R")

# =============================================================================
# SECTION 32: CHG CONTEXT EXPLORATORY ANALYSIS
# =============================================================================

cat("================================================================================\n")
cat("SECTION 32: CHG CONTEXT EXPLORATORY ANALYSIS\n")
cat("================================================================================\n\n")

# -----------------------------------------------------------------------------
# SETUP: Load CHG data
# -----------------------------------------------------------------------------

CHG_DATA_PATHS <- list(
  mc_dmr = file.path(BASE_DIR, "modality/outputs/run-3/outputs_CHG/Results/gencode.vM25.mouse.genes.annotation/DMR_20260221_195318/DMR_mc_control__mutant_20260221_195318.bed"),
  hmc_dmr = file.path(BASE_DIR, "modality/outputs/run-3/outputs_CHG/Results/gencode.vM25.mouse.genes.annotation/DMR_20260221_195318/DMR_hmc_control__mutant_20260221_195318.bed")
)

stopifnot("CHG mC DMR file not found" = file.exists(CHG_DATA_PATHS$mc_dmr))
stopifnot("CHG hmC DMR file not found" = file.exists(CHG_DATA_PATHS$hmc_dmr))

cat("Loading CHG gene body mC DMRs...\n")
chg_mc_dmr <- load_dmr_bed(CHG_DATA_PATHS$mc_dmr)
cat(sprintf("  Loaded %d genes, %d significant (q < %.2f)\n",
            nrow(chg_mc_dmr), sum(chg_mc_dmr$significant), Q_THRESHOLD))

cat("Loading CHG gene body hmC DMRs...\n")
chg_hmc_dmr <- load_dmr_bed(CHG_DATA_PATHS$hmc_dmr)
cat(sprintf("  Loaded %d genes, %d significant (q < %.2f)\n",
            nrow(chg_hmc_dmr), sum(chg_hmc_dmr$significant), Q_THRESHOLD))

# Extract significant subsets
chg_mc_sig <- chg_mc_dmr %>%
  dplyr::filter(significant) %>%
  dplyr::distinct(gene, .keep_all = TRUE)

chg_hmc_sig <- chg_hmc_dmr %>%
  dplyr::filter(significant) %>%
  dplyr::distinct(gene, .keep_all = TRUE)

cat(sprintf("\nCHG mC significant: %d genes (%d hyper, %d hypo)\n",
            nrow(chg_mc_sig),
            sum(chg_mc_sig$mod_difference > 0),
            sum(chg_mc_sig$mod_difference <= 0)))
cat(sprintf("CHG hmC significant: %d genes\n", nrow(chg_hmc_sig)))

# CHG-specific colors
CHG_COLORS <- list(
  context = c("CG" = "#E41A1C", "CHG" = "#FF7F00"),
  concordance = c("Concordant" = "#4DAF4A", "Discordant" = "#E41A1C")
)

# =============================================================================
# FIGURE 32a/b: CHG VOLCANO PLOTS
# =============================================================================

cat("\n--- Figure 32a/b: CHG Volcano Plots ---\n")

create_chg_volcano <- function(df, sig_df, title, subtitle, meth_type = "CHG mC",
                               n_labels = 15) {
  # Label top genes by q-value
  top_genes <- sig_df %>%
    dplyr::arrange(dmr_qvalue) %>%
    head(n_labels) %>%
    pull(gene)

  df$label <- ""
  df$label[df$gene %in% top_genes & df$significant] <- df$gene[df$gene %in% top_genes & df$significant]

  df$color_group <- case_when(
    !df$significant ~ "Not Significant",
    df$mod_difference > 0 ~ "Hypermethylated",
    TRUE ~ "Hypomethylated"
  )

  df$neg_log10_q_capped <- pmin(df$neg_log10_q, 300)

  # Use permille scale for CHG (effect sizes are ~0.001)
  p <- ggplot(df, aes(x = mod_difference * 1000, y = neg_log10_q_capped)) +
    geom_point(aes(color = color_group), alpha = 0.6, size = 1.5) +
    geom_hline(yintercept = -log10(Q_THRESHOLD), linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = 0, linetype = "solid", color = "grey40") +
    geom_text_repel(
      aes(label = label),
      size = 3.5,
      max.overlaps = 25,
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
      x = paste(meth_type, "Difference (permille)"),
      y = expression(-log[10](q-value))
    ) +
    theme_biomodal() +
    theme(legend.position = "bottom")

  return(p)
}

n_hyper_mc <- sum(chg_mc_sig$mod_difference > 0)
n_hypo_mc <- sum(chg_mc_sig$mod_difference <= 0)

p_chg_volcano_mc <- create_chg_volcano(
  chg_mc_dmr, chg_mc_sig,
  "Gene Body CHG 5mC Differential Methylation",
  sprintf("BAP1-KO vs Control | %d significant (q < %.2f) | %d hyper, %d hypo",
          nrow(chg_mc_sig), Q_THRESHOLD, n_hyper_mc, n_hypo_mc),
  "CHG mC", n_labels = 15
)

save_multiformat_ggplot(p_chg_volcano_mc, file.path(OUTPUT_DIR, "32a_chg_volcano_mc"),
                        width = 10, height = 8)

p_chg_volcano_hmc <- create_chg_volcano(
  chg_hmc_dmr, chg_hmc_sig,
  "Gene Body CHG 5hmC Differential Methylation",
  sprintf("BAP1-KO vs Control | %d significant (q < %.2f)",
          nrow(chg_hmc_sig), Q_THRESHOLD),
  "CHG hmC", n_labels = 7
)

save_multiformat_ggplot(p_chg_volcano_hmc, file.path(OUTPUT_DIR, "32b_chg_volcano_hmc"),
                        width = 10, height = 8)

cat("  Saved 32a_chg_volcano_mc and 32b_chg_volcano_hmc\n")

# =============================================================================
# FIGURE 32c: DIRECTION BREAKDOWN
# =============================================================================

cat("\n--- Figure 32c: Direction Breakdown ---\n")

direction_df <- data.frame(
  context = c("CHG mC", "CHG mC", "CHG hmC", "CHG hmC"),
  direction = c("Hypermethylated", "Hypomethylated", "Hypermethylated", "Hypomethylated"),
  count = c(
    sum(chg_mc_sig$mod_difference > 0),
    sum(chg_mc_sig$mod_difference <= 0),
    sum(chg_hmc_sig$mod_difference > 0),
    sum(chg_hmc_sig$mod_difference <= 0)
  ),
  stringsAsFactors = FALSE
)

direction_df <- direction_df %>%
  dplyr::group_by(context) %>%
  dplyr::mutate(
    total = sum(count),
    pct = count / total * 100,
    label = sprintf("%d (%.1f%%)", count, pct)
  ) %>%
  dplyr::ungroup()

p_direction <- ggplot(direction_df, aes(x = context, y = count, fill = direction)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_text(aes(label = label), position = position_dodge(width = 0.7),
            vjust = -0.3, size = 4) +
  scale_fill_manual(values = COLORS$direction, name = "Direction") +
  labs(
    title = "CHG Significant DMR Direction",
    subtitle = "Gene body DMRs | BAP1-KO vs Control",
    x = "", y = "Number of Genes"
  ) +
  theme_biomodal() +
  theme(legend.position = "bottom")

save_multiformat_ggplot(p_direction, file.path(OUTPUT_DIR, "32c_chg_direction_breakdown"),
                        width = 8, height = 6)
cat("  Saved 32c_chg_direction_breakdown\n")

# =============================================================================
# FIGURE 32d/e: CG vs CHG VENN DIAGRAMS
# =============================================================================

cat("\n--- Figure 32d/e: CG vs CHG Venn Diagrams ---\n")

# CG significant gene sets (loaded by shared config)
cg_mc_genes <- mc_dmr$gene[mc_dmr$significant]
cg_hmc_genes <- hmc_dmr$gene[hmc_dmr$significant]
chg_mc_genes <- chg_mc_sig$gene

# 2-set Venn: CG mC vs CHG mC
venn_2set <- list(
  "CG mC Significant" = cg_mc_genes,
  "CHG mC Significant" = chg_mc_genes
)

n_overlap <- length(intersect(cg_mc_genes, chg_mc_genes))
n_chg_unique <- length(setdiff(chg_mc_genes, cg_mc_genes))

p_venn_2 <- ggVennDiagram(venn_2set, label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "#FF7F00") +
  labs(
    title = "CHG vs CG: Significant mC Gene Overlap",
    subtitle = sprintf("%d overlap | %d CHG-unique | %d CG-only",
                       n_overlap, n_chg_unique,
                       length(setdiff(cg_mc_genes, chg_mc_genes)))
  ) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))

save_multiformat_ggplot(p_venn_2, file.path(OUTPUT_DIR, "32d_chg_cg_venn_mc"),
                        width = 8, height = 7)

# 3-set Venn: CG mC + CG hmC + CHG mC
venn_3set <- list(
  "CG mC" = cg_mc_genes,
  "CG hmC" = cg_hmc_genes,
  "CHG mC" = chg_mc_genes
)

p_venn_3 <- ggVennDiagram(venn_3set, label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "#4DBBD5") +
  labs(
    title = "Three-Way Overlap: CG mC, CG hmC, and CHG mC",
    subtitle = "Significant genes across methylation contexts"
  ) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))

save_multiformat_ggplot(p_venn_3, file.path(OUTPUT_DIR, "32e_chg_cg_venn_three_way"),
                        width = 9, height = 8)

cat("  Saved 32d_chg_cg_venn_mc and 32e_chg_cg_venn_three_way\n")

# =============================================================================
# FIGURE 32f/g: DIRECTION CONCORDANCE (mC + hmC)
# =============================================================================

cat("\n--- Figure 32f/g: Direction Concordance ---\n")

# --- mC concordance ---
cg_mc_sig <- mc_dmr %>%
  dplyr::filter(significant) %>%
  dplyr::distinct(gene, .keep_all = TRUE) %>%
  dplyr::select(gene,
                cg_mc_difference = mod_difference,
                cg_mc_qvalue = dmr_qvalue,
                cg_mc_direction = direction)

concordance_mc <- chg_mc_sig %>%
  dplyr::select(gene,
                chg_mc_difference = mod_difference,
                chg_mc_qvalue = dmr_qvalue,
                chg_mc_direction = direction) %>%
  dplyr::inner_join(cg_mc_sig, by = "gene") %>%
  dplyr::mutate(
    concordant = sign(chg_mc_difference) == sign(cg_mc_difference),
    concordance = ifelse(concordant, "Concordant", "Discordant")
  )

n_conc <- sum(concordance_mc$concordant)
n_disc <- sum(!concordance_mc$concordant)

cat(sprintf("  mC concordance: %d concordant, %d discordant (%.1f%% discordant)\n",
            n_conc, n_disc, 100 * n_disc / nrow(concordance_mc)))

# Select genes to label (top by combined absolute effect)
concordance_mc <- concordance_mc %>%
  dplyr::mutate(combined_effect = abs(chg_mc_difference) + abs(cg_mc_difference))

label_genes <- concordance_mc %>%
  dplyr::arrange(desc(combined_effect)) %>%
  head(20) %>%
  pull(gene)

concordance_mc$label <- ifelse(concordance_mc$gene %in% label_genes,
                                concordance_mc$gene, "")

p_concordance_scatter <- ggplot(concordance_mc,
                                 aes(x = cg_mc_difference * 100,
                                     y = chg_mc_difference * 1000)) +
  geom_hline(yintercept = 0, color = "grey40", linetype = "solid") +
  geom_vline(xintercept = 0, color = "grey40", linetype = "solid") +
  geom_point(aes(color = concordance), size = 3, alpha = 0.8) +
  geom_text_repel(
    aes(label = label),
    size = 3, max.overlaps = 25, box.padding = 0.4, segment.color = "grey50"
  ) +
  scale_color_manual(values = CHG_COLORS$concordance, name = "") +
  annotate("text", x = Inf, y = Inf,
           label = sprintf("%d/%d discordant (%.1f%%)",
                           n_disc, nrow(concordance_mc),
                           100 * n_disc / nrow(concordance_mc)),
           hjust = 1.1, vjust = 1.5, size = 4.5, fontface = "bold", color = "#E41A1C") +
  labs(
    title = "CHG vs CG mC Direction Concordance",
    subtitle = sprintf("%d overlapping genes | CHG and CG methylation change direction", nrow(concordance_mc)),
    x = "CG mC Difference (%)",
    y = "CHG mC Difference (permille)"
  ) +
  theme_biomodal() +
  theme(legend.position = "bottom")

save_multiformat_ggplot(p_concordance_scatter,
                        file.path(OUTPUT_DIR, "32f_chg_cg_concordance_scatter"),
                        width = 10, height = 9)

# --- hmC concordance ---
cg_hmc_sig_df <- hmc_dmr %>%
  dplyr::filter(significant) %>%
  dplyr::distinct(gene, .keep_all = TRUE) %>%
  dplyr::select(gene,
                cg_hmc_difference = mod_difference,
                cg_hmc_qvalue = dmr_qvalue,
                cg_hmc_direction = direction)

concordance_hmc <- chg_hmc_sig %>%
  dplyr::select(gene,
                chg_hmc_difference = mod_difference,
                chg_hmc_qvalue = dmr_qvalue,
                chg_hmc_direction = direction) %>%
  dplyr::inner_join(cg_hmc_sig_df, by = "gene")

if (nrow(concordance_hmc) > 0) {
  concordance_hmc <- concordance_hmc %>%
    dplyr::mutate(
      concordant = sign(chg_hmc_difference) == sign(cg_hmc_difference),
      concordance = ifelse(concordant, "Concordant", "Discordant")
    )

  n_conc_hmc <- sum(concordance_hmc$concordant)
  n_disc_hmc <- sum(!concordance_hmc$concordant)
  cat(sprintf("  hmC concordance: %d concordant, %d discordant out of %d overlap\n",
              n_conc_hmc, n_disc_hmc, nrow(concordance_hmc)))
}

# Summary bar chart (mC + hmC)
conc_summary <- data.frame(
  context = c("mC", "mC", "hmC", "hmC"),
  concordance = c("Concordant", "Discordant", "Concordant", "Discordant"),
  count = c(n_conc, n_disc,
            if (nrow(concordance_hmc) > 0) n_conc_hmc else 0,
            if (nrow(concordance_hmc) > 0) n_disc_hmc else 0),
  stringsAsFactors = FALSE
)

conc_summary <- conc_summary %>%
  dplyr::group_by(context) %>%
  dplyr::mutate(
    total = sum(count),
    pct = ifelse(total > 0, count / total * 100, 0),
    label = sprintf("%d (%.0f%%)", count, pct)
  ) %>%
  dplyr::ungroup()

p_concordance_bar <- ggplot(conc_summary, aes(x = context, y = count, fill = concordance)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_text(aes(label = label), position = position_dodge(width = 0.7),
            vjust = -0.3, size = 4) +
  scale_fill_manual(values = CHG_COLORS$concordance, name = "") +
  labs(
    title = "CHG vs CG Direction Concordance Summary",
    subtitle = "Do CHG and CG methylation change in the same direction?",
    x = "Methylation Type", y = "Number of Genes"
  ) +
  theme_biomodal() +
  theme(legend.position = "bottom")

save_multiformat_ggplot(p_concordance_bar,
                        file.path(OUTPUT_DIR, "32g_chg_cg_concordance_bar"),
                        width = 8, height = 6)

cat("  Saved 32f_chg_cg_concordance_scatter and 32g_chg_cg_concordance_bar\n")

# =============================================================================
# FIGURE 32h: ALL 70 CHG DMR GENES (BAR CHART)
# =============================================================================

cat("\n--- Figure 32h: All CHG mC DMR Genes ---\n")

# Annotate with CG overlap status
chg_mc_all <- chg_mc_sig %>%
  dplyr::mutate(
    in_cg_mc = gene %in% cg_mc_genes,
    cg_status = ifelse(in_cg_mc, "Also CG sig", "CHG-unique"),
    gene = factor(gene, levels = rev(gene[order(dmr_qvalue)]))
  )

p_all_genes <- ggplot(chg_mc_all, aes(x = gene, y = mod_difference * 1000, fill = direction)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_hline(yintercept = 0, color = "black") +
  geom_point(aes(y = mod_difference * 1000, shape = cg_status),
             size = 2, color = "black", show.legend = TRUE,
             position = position_nudge(y = ifelse(chg_mc_all$mod_difference > 0, 0.1, -0.1))) +
  scale_fill_manual(values = COLORS$direction, name = "Direction") +
  scale_shape_manual(values = c("Also CG sig" = 16, "CHG-unique" = 4),
                     name = "CG Overlap") +
  coord_flip() +
  labs(
    title = sprintf("All %d Significant CHG mC DMR Genes", nrow(chg_mc_all)),
    subtitle = "Ranked by q-value | BAP1-KO vs Control",
    x = "", y = "CHG mC Difference (permille)"
  ) +
  theme_biomodal() +
  theme(axis.text.y = element_text(size = 8))

save_multiformat_ggplot(p_all_genes, file.path(OUTPUT_DIR, "32h_chg_top_dmr_genes"),
                        width = 10, height = 18)

cat("  Saved 32h_chg_top_dmr_genes\n")

# =============================================================================
# FIGURE 32i: CHG-UNIQUE GENES
# =============================================================================

cat("\n--- Figure 32i: CHG-Unique Genes ---\n")

chg_unique <- chg_mc_sig %>%
  dplyr::filter(!gene %in% cg_mc_genes) %>%
  dplyr::mutate(
    in_cg_hmc = gene %in% cg_hmc_genes,
    cg_hmc_status = ifelse(in_cg_hmc, "In CG hmC sig", "Not in CG"),
    gene = factor(gene, levels = rev(gene[order(dmr_qvalue)]))
  )

cat(sprintf("  %d CHG-unique genes (%d also in CG hmC)\n",
            nrow(chg_unique), sum(chg_unique$in_cg_hmc)))

p_unique <- ggplot(chg_unique, aes(x = gene, y = mod_difference * 1000)) +
  geom_segment(aes(x = gene, xend = gene, y = 0, yend = mod_difference * 1000,
                   color = direction), linewidth = 1.2) +
  geom_point(aes(color = direction, shape = cg_hmc_status), size = 4) +
  geom_hline(yintercept = 0, color = "grey40") +
  scale_color_manual(values = COLORS$direction, name = "Direction") +
  scale_shape_manual(values = c("In CG hmC sig" = 17, "Not in CG" = 16),
                     name = "CG hmC Status") +
  coord_flip() +
  labs(
    title = sprintf("%d CHG-Unique Genes (Not in CG mC Significant)", nrow(chg_unique)),
    subtitle = "Genes with significant CHG mC changes only",
    x = "", y = "CHG mC Difference (permille)"
  ) +
  theme_biomodal()

save_multiformat_ggplot(p_unique, file.path(OUTPUT_DIR, "32i_chg_unique_genes"),
                        width = 10, height = 7)

cat("  Saved 32i_chg_unique_genes\n")

# =============================================================================
# FIGURE 32j/k: CHROMOSOME DISTRIBUTION
# =============================================================================

cat("\n--- Figure 32j/k: Chromosome Distribution ---\n")

# Count significant CHG mC DMRs per chromosome
chr_counts <- chg_mc_sig %>%
  dplyr::count(chr, name = "sig_count")

# Expected counts based on total genes per chromosome
chr_total <- chg_mc_dmr %>%
  dplyr::count(chr, name = "total_genes")

chr_dist <- chr_total %>%
  dplyr::left_join(chr_counts, by = "chr") %>%
  dplyr::mutate(
    sig_count = replace(sig_count, is.na(sig_count), 0),
    expected = total_genes / sum(total_genes) * nrow(chg_mc_sig),
    highlight = ifelse(chr %in% c("chr7", "chr8"), "Enriched", "Other"),
    chr_label = gsub("chr", "", chr),
    chr_num = suppressWarnings(as.numeric(chr_label)),
    chr_order = ifelse(is.na(chr_num), 100 + match(chr_label, c("X", "Y", "M")), chr_num)
  ) %>%
  dplyr::arrange(chr_order) %>%
  dplyr::mutate(chr = factor(chr, levels = chr))

# Fisher's exact test per chromosome
chr_dist$fisher_p <- NA_real_
for (i in seq_len(nrow(chr_dist))) {
  sig_on_chr <- chr_dist$sig_count[i]
  sig_off_chr <- nrow(chg_mc_sig) - sig_on_chr
  total_on_chr <- chr_dist$total_genes[i]
  total_off_chr <- sum(chr_dist$total_genes) - total_on_chr
  nonsig_on_chr <- total_on_chr - sig_on_chr
  nonsig_off_chr <- total_off_chr - sig_off_chr

  mat <- matrix(c(sig_on_chr, nonsig_on_chr, sig_off_chr, nonsig_off_chr),
                nrow = 2)
  ft <- fisher.test(mat)
  chr_dist$fisher_p[i] <- ft$p.value
}

chr_dist$fisher_label <- ifelse(chr_dist$fisher_p < 0.001,
                                 sprintf("p=%.1e", chr_dist$fisher_p),
                                 ifelse(chr_dist$fisher_p < 0.05,
                                        sprintf("p=%.3f", chr_dist$fisher_p), ""))

p_chr_dist <- ggplot(chr_dist, aes(x = chr, y = sig_count)) +
  geom_bar(aes(fill = highlight), stat = "identity", width = 0.7) +
  geom_point(aes(y = expected), color = "black", shape = 4, size = 3) +
  geom_text(aes(label = fisher_label), vjust = -0.5, size = 3, color = "#E41A1C",
            fontface = "bold") +
  scale_fill_manual(values = c("Enriched" = "#FF7F00", "Other" = "grey70"),
                    guide = "none") +
  labs(
    title = "Chromosome Distribution of Significant CHG mC DMRs",
    subtitle = sprintf("Orange = enriched chromosomes | X = expected count (genome-wide rate) | %d/%d on chr7+chr8",
                       sum(chr_dist$sig_count[chr_dist$highlight == "Enriched"]),
                       nrow(chg_mc_sig)),
    x = "Chromosome", y = "Significant DMR Count"
  ) +
  theme_biomodal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_multiformat_ggplot(p_chr_dist, file.path(OUTPUT_DIR, "32j_chg_chromosome_distribution"),
                        width = 12, height = 7)

# Comparison with CG mC chromosome distribution
cg_chr_counts <- mc_dmr %>%
  dplyr::filter(significant) %>%
  dplyr::count(chr, name = "cg_sig_count")

chr_compare <- chr_dist %>%
  dplyr::left_join(cg_chr_counts, by = "chr") %>%
  dplyr::mutate(cg_sig_count = replace(cg_sig_count, is.na(cg_sig_count), 0)) %>%
  dplyr::mutate(
    chg_pct = sig_count / sum(sig_count) * 100,
    cg_pct = cg_sig_count / sum(cg_sig_count) * 100
  ) %>%
  tidyr::pivot_longer(cols = c(chg_pct, cg_pct),
                      names_to = "context", values_to = "pct") %>%
  dplyr::mutate(context = ifelse(context == "chg_pct", "CHG mC", "CG mC"))

p_chr_compare <- ggplot(chr_compare, aes(x = chr, y = pct, fill = context)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(values = CHG_COLORS$context, name = "Context") +
  labs(
    title = "Chromosome Distribution: CHG vs CG Significant mC DMRs",
    subtitle = "Normalized to percentage | CHG shows strong chr7+chr8 enrichment",
    x = "Chromosome", y = "% of Significant DMRs"
  ) +
  theme_biomodal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

save_multiformat_ggplot(p_chr_compare, file.path(OUTPUT_DIR, "32k_chg_vs_cg_chr_comparison"),
                        width = 12, height = 7)

cat("  Saved 32j_chg_chromosome_distribution and 32k_chg_vs_cg_chr_comparison\n")

# =============================================================================
# FIGURE 32l/m: GO/KEGG ENRICHMENT
# =============================================================================

cat("\n--- Figure 32l/m: GO/KEGG Enrichment ---\n")

all_chg_genes <- unique(chg_mc_sig$gene)
cat(sprintf("Running enrichment on %d CHG mC significant genes...\n", length(all_chg_genes)))

tryCatch({
  gene_ids <- bitr(all_chg_genes, fromType = "SYMBOL", toType = "ENTREZID",
                   OrgDb = org.Mm.eg.db)

  if (nrow(gene_ids) > 10) {
    cat(sprintf("  Converted %d / %d genes to Entrez IDs\n", nrow(gene_ids), length(all_chg_genes)))

    # GO Biological Process (relaxed threshold for exploratory analysis)
    cat("Running GO Biological Process enrichment...\n")
    go_bp <- enrichGO(gene = gene_ids$ENTREZID,
                      OrgDb = org.Mm.eg.db,
                      ont = "BP",
                      pAdjustMethod = "BH",
                      qvalueCutoff = 0.1,
                      readable = TRUE)

    if (!is.null(go_bp) && nrow(go_bp@result[go_bp@result$qvalue < 0.1, ]) > 0) {
      p_go_bp <- dotplot(go_bp, showCategory = 15) +
        labs(title = "GO Biological Process: CHG mC Significant Genes",
             subtitle = sprintf("%d genes | Exploratory (q < 0.1)", length(all_chg_genes))) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              plot.subtitle = element_text(hjust = 0.5))

      save_multiformat_ggplot(p_go_bp, file.path(OUTPUT_DIR, "32l_chg_enrichment_go"),
                              width = 10, height = 10)

      write.table(go_bp@result, file.path(TABLES_DIR, "chg_enrichment_go_bp.tsv"),
                  sep = "\t", quote = FALSE, row.names = FALSE)
      cat("  Saved GO BP results\n")
    } else {
      cat("  No significant GO BP terms at q < 0.1\n")
    }

    # KEGG Pathway (relaxed threshold)
    cat("Running KEGG pathway enrichment...\n")
    kegg <- enrichKEGG(gene = gene_ids$ENTREZID,
                       organism = "mmu",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 1,
                       qvalueCutoff = 0.2)

    kegg_sig <- kegg@result[kegg@result$qvalue < 0.2, ]
    if (!is.null(kegg) && nrow(kegg_sig) > 0) {
      # Build manual dotplot since dotplot() uses internal pvalueCutoff for rendering
      kegg_plot_df <- kegg_sig %>%
        head(15) %>%
        dplyr::mutate(
          GeneRatio_num = sapply(GeneRatio, function(x) {
            parts <- as.numeric(strsplit(x, "/")[[1]])
            parts[1] / parts[2]
          }),
          Description = factor(Description, levels = rev(Description))
        )

      p_kegg <- ggplot(kegg_plot_df, aes(x = GeneRatio_num, y = Description)) +
        geom_point(aes(size = Count, color = qvalue)) +
        scale_color_gradient(low = "#E41A1C", high = "#377EB8", name = "q-value") +
        scale_size_continuous(range = c(3, 8), name = "Count") +
        labs(
          title = "KEGG Pathway: CHG mC Significant Genes",
          subtitle = sprintf("%d genes | Exploratory (q < 0.2) | %d pathways",
                             length(all_chg_genes), nrow(kegg_sig)),
          x = "Gene Ratio", y = ""
        ) +
        theme_biomodal()

      save_multiformat_ggplot(p_kegg, file.path(OUTPUT_DIR, "32m_chg_enrichment_kegg"),
                              width = 10, height = 10)

      write.table(kegg@result[kegg@result$qvalue < 0.5, ],
                  file.path(TABLES_DIR, "chg_enrichment_kegg.tsv"),
                  sep = "\t", quote = FALSE, row.names = FALSE)
      cat("  Saved KEGG results\n")
    } else {
      cat("  No significant KEGG pathways at q < 0.2\n")
    }
  } else {
    cat("  Not enough genes converted for enrichment analysis\n")
  }
}, error = function(e) {
  cat(sprintf("  Enrichment analysis error: %s\n", e$message))
})

# =============================================================================
# TABLE EXPORTS
# =============================================================================

cat("\n--- Table Exports ---\n")

# Master table: all 70 CHG mC significant genes with CG cross-reference
cg_mc_for_join <- mc_dmr %>%
  dplyr::distinct(gene, .keep_all = TRUE) %>%
  dplyr::select(gene,
                cg_mc_difference = mod_difference,
                cg_mc_qvalue = dmr_qvalue,
                cg_mc_significant = significant,
                cg_mc_direction = direction)

cg_hmc_for_join <- hmc_dmr %>%
  dplyr::distinct(gene, .keep_all = TRUE) %>%
  dplyr::select(gene,
                cg_hmc_difference = mod_difference,
                cg_hmc_qvalue = dmr_qvalue,
                cg_hmc_significant = significant,
                cg_hmc_direction = direction)

master_table <- chg_mc_sig %>%
  dplyr::select(gene, chr, start, end, num_contexts, mean_coverage,
                chg_mc_group1 = mean_mod_group1, chg_mc_group2 = mean_mod_group2,
                chg_mc_difference = mod_difference, chg_mc_qvalue = dmr_qvalue,
                chg_mc_direction = direction) %>%
  dplyr::left_join(cg_mc_for_join, by = "gene") %>%
  dplyr::left_join(cg_hmc_for_join, by = "gene") %>%
  dplyr::mutate(
    concordant_with_cg_mc = dplyr::case_when(
      is.na(cg_mc_significant) | !cg_mc_significant ~ NA,
      sign(chg_mc_difference) == sign(cg_mc_difference) ~ TRUE,
      TRUE ~ FALSE
    )
  ) %>%
  dplyr::arrange(chg_mc_qvalue)

write.table(master_table, file.path(TABLES_DIR, "chg_mc_significant_dmrs.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved chg_mc_significant_dmrs.tsv (%d genes)\n", nrow(master_table)))

# CHG hmC significant genes
if (nrow(chg_hmc_sig) > 0) {
  hmc_table <- chg_hmc_sig %>%
    dplyr::select(gene, chr, start, end, num_contexts, mean_coverage,
                  chg_hmc_group1 = mean_mod_group1, chg_hmc_group2 = mean_mod_group2,
                  chg_hmc_difference = mod_difference, chg_hmc_qvalue = dmr_qvalue,
                  chg_hmc_direction = direction)

  write.table(hmc_table, file.path(TABLES_DIR, "chg_hmc_significant_dmrs.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  Saved chg_hmc_significant_dmrs.tsv (%d genes)\n", nrow(hmc_table)))
}

# CHG-unique gene list
unique_table <- chg_unique %>%
  dplyr::select(gene, chr, start, end,
                chg_mc_difference = mod_difference, chg_mc_qvalue = dmr_qvalue,
                chg_mc_direction = direction, in_cg_hmc)

write.table(unique_table, file.path(TABLES_DIR, "chg_unique_genes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved chg_unique_genes.tsv (%d genes)\n", nrow(unique_table)))

# Print summary
cat("\n")
cat("================================================================================\n")
cat("CHG EXPLORATORY ANALYSIS SUMMARY\n")
cat("================================================================================\n")
cat(sprintf("CHG mC significant DMRs:  %d (q < %.2f)\n", nrow(chg_mc_sig), Q_THRESHOLD))
cat(sprintf("  Hypermethylated:        %d (%.1f%%)\n", n_hyper_mc, 100 * n_hyper_mc / nrow(chg_mc_sig)))
cat(sprintf("  Hypomethylated:         %d (%.1f%%)\n", n_hypo_mc, 100 * n_hypo_mc / nrow(chg_mc_sig)))
cat(sprintf("CHG hmC significant DMRs: %d\n", nrow(chg_hmc_sig)))
cat(sprintf("Overlap with CG mC:       %d / %d (%.1f%%)\n",
            n_overlap, nrow(chg_mc_sig), 100 * n_overlap / nrow(chg_mc_sig)))
cat(sprintf("CHG-unique genes:         %d\n", n_chg_unique))
cat(sprintf("Direction discordance:    %d / %d (%.1f%%)\n",
            n_disc, nrow(concordance_mc), 100 * n_disc / nrow(concordance_mc)))
cat(sprintf("Chr7+Chr8 concentration:  %d / %d (%.1f%%)\n",
            sum(chr_dist$sig_count[chr_dist$highlight == "Enriched"]),
            nrow(chg_mc_sig),
            100 * sum(chr_dist$sig_count[chr_dist$highlight == "Enriched"]) / nrow(chg_mc_sig)))
cat("================================================================================\n")

cat("Section 32 complete.\n\n")
