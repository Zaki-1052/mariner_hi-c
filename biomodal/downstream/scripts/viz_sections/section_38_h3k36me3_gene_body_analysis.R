# biomodal/downstream/scripts/viz_sections/section_38_h3k36me3_gene_body_analysis.R
# H3K36me3 Gene Body Mark Analysis (SETD2/DNMT3B axis)
#
# H3K36me3 is deposited by SETD2 exclusively in actively transcribed gene bodies,
# coupled to RNA Pol II elongation. It recruits DNMT3B for gene-body methylation.
# This section tests whether H3K36me3 changes track with mC/hmC changes at gene
# bodies, addressing the DNMT3B recruitment hypothesis.
#
# Data: Pre-annotated DiffBind results (39,635 peaks with gene annotations)
# Figures: 38a-38h (8 total)
#
# Usage:
#   cd downstream/
#   Rscript scripts/viz_sections/section_38_h3k36me3_gene_body_analysis.R

source("scripts/viz_sections/_shared_config.R")
suppressPackageStartupMessages(library(ChIPseeker))

SECTION_DIR <- "38_h3k36me3_gene_body"
FDR_THRESHOLD <- 0.05

dir.create(file.path(OUTPUT_DIR, SECTION_DIR), recursive = TRUE, showWarnings = FALSE)

cat("\n")
cat("================================================================================\n")
cat("SECTION 38: H3K36me3 Gene Body Analysis (SETD2/DNMT3B Axis)\n")
cat("================================================================================\n\n")

# =============================================================================
# DATA LOADING
# =============================================================================

cat("Loading H3K36me3 DiffBind data...\n")
me3_all <- load_diffbind_flex(H3K36ME3_FILES$master_annotated, "H3K36me3 (all peaks)")

# Load demethylation ratio table for gene-level methylation data
RATIO_TABLE <- file.path(TABLES_DIR, "demethylation_ratio_all_genes.tsv")
stopifnot(file.exists(RATIO_TABLE))
ratio_data <- read.table(RATIO_TABLE, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat(sprintf("  Ratio table: %d genes loaded\n", nrow(ratio_data)))

# Load differential BED files for coordinate-level overlap
cat("Loading H3K36me3 differential BEDs...\n")
me3_up_gr <- load_chip_peaks(H3K36ME3_FILES$up, "H3K36me3 Up")
me3_down_gr <- load_chip_peaks(H3K36ME3_FILES$down, "H3K36me3 Down")

# =============================================================================
# GENE-LEVEL AGGREGATION (nearest-to-TSS approach)
# =============================================================================

cat("\nAggregating H3K36me3 peaks to gene level...\n")

me3_gene <- me3_all %>%
  dplyr::filter(!is.na(SYMBOL) & SYMBOL != "") %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::summarise(
    me3_fold = Fold[which.min(abs(distanceToTSS))],
    me3_fdr  = FDR[which.min(abs(distanceToTSS))],
    me3_conc = Conc[which.min(abs(distanceToTSS))],
    me3_n_peaks = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::rename(gene = SYMBOL)

cat(sprintf("  %d unique genes with H3K36me3 peaks\n", nrow(me3_gene)))

# Merge with methylation data
merged <- ratio_data %>%
  left_join(me3_gene, by = "gene")

n_with_me3 <- sum(!is.na(merged$me3_fold))
cat(sprintf("  %d genes with both methylation and me3 data (%.1f%%)\n",
            n_with_me3, 100 * n_with_me3 / nrow(merged)))

# =============================================================================
# FIGURE 38a: H3K36me3 DiffBind Overview (Volcano + Annotation Pie)
# =============================================================================

cat("\n--- FIGURE 38a: H3K36me3 DiffBind Overview ---\n\n")

me3_all$sig <- me3_all$FDR < FDR_THRESHOLD
me3_all$neg_log10_fdr <- -log10(pmax(me3_all$FDR, 1e-300))
me3_all$neg_log10_fdr[is.infinite(me3_all$neg_log10_fdr)] <- 300

me3_all$color_group <- case_when(
  me3_all$sig & me3_all$Fold > 0 ~ "Up (Gained)",
  me3_all$sig & me3_all$Fold < 0 ~ "Down (Lost)",
  TRUE ~ "Not Significant"
)

n_up <- sum(me3_all$color_group == "Up (Gained)")
n_down <- sum(me3_all$color_group == "Down (Lost)")
n_ns <- sum(me3_all$color_group == "Not Significant")

# Subsample non-significant for rendering
set.seed(42)
ns_subset <- me3_all %>% dplyr::filter(color_group == "Not Significant") %>% dplyr::sample_frac(0.2)
sig_peaks <- me3_all %>% dplyr::filter(color_group != "Not Significant")
plot_data_volcano <- bind_rows(ns_subset, sig_peaks)

p_volcano <- ggplot(plot_data_volcano, aes(x = Fold, y = neg_log10_fdr, color = color_group)) +
  geom_point(alpha = 0.4, size = 0.6) +
  scale_color_manual(
    values = c("Up (Gained)" = "#D95F02", "Down (Lost)" = "#1B9E77", "Not Significant" = "grey70"),
    name = "Direction"
  ) +
  geom_hline(yintercept = -log10(FDR_THRESHOLD), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  labs(
    title = "H3K36me3 DiffBind: All Peaks",
    subtitle = sprintf("%s up, %s down, %s NS (FDR<%.2f)",
                       format(n_up, big.mark = ","), format(n_down, big.mark = ","),
                       format(n_ns, big.mark = ","), FDR_THRESHOLD),
    x = "Log2 Fold Change (Mutant / Control)",
    y = "-log10(FDR)"
  ) +
  theme_biomodal()

# Annotation pie for significant peaks
sig_anno <- me3_all %>%
  dplyr::filter(sig & !is.na(annotation)) %>%
  mutate(
    anno_simple = case_when(
      grepl("Promoter", annotation) ~ "Promoter",
      grepl("Exon", annotation)     ~ "Exon",
      grepl("3' UTR", annotation)   ~ "3' UTR",
      grepl("5' UTR", annotation)   ~ "5' UTR",
      grepl("Intron", annotation)   ~ "Intron",
      grepl("Downstream", annotation) ~ "Downstream",
      grepl("Intergenic", annotation) ~ "Distal Intergenic",
      TRUE ~ "Other"
    )
  )

anno_counts <- sig_anno %>%
  count(anno_simple) %>%
  arrange(desc(n)) %>%
  mutate(pct = 100 * n / sum(n),
         label = sprintf("%s\n%.1f%%", anno_simple, pct))

anno_colors <- c("Promoter" = "#E41A1C", "Exon" = "#377EB8", "Intron" = "#4DAF4A",
                 "3' UTR" = "#984EA3", "5' UTR" = "#FF7F00", "Distal Intergenic" = "#A65628",
                 "Downstream" = "#F781BF", "Other" = "#999999")

p_anno <- ggplot(anno_counts, aes(x = "", y = n, fill = anno_simple)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y") +
  scale_fill_manual(values = anno_colors, name = "Annotation") +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 2.5, lineheight = 0.85) +
  labs(title = "Significant H3K36me3 Peak Annotations",
       subtitle = sprintf("n = %d differential peaks", sum(me3_all$sig))) +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10))

p_38a <- p_volcano + p_anno + plot_layout(widths = c(2, 1))

save_multiformat_ggplot(p_38a, file.path(OUTPUT_DIR, SECTION_DIR, "38a_h3k36me3_volcano_annotation"),
                        width = 14, height = 7)

# =============================================================================
# FIGURE 38b: O/E Heatmap — me3 Direction x mC Direction
# =============================================================================

cat("\n--- FIGURE 38b: me3 x mC Direction O/E Heatmap ---\n\n")

# Gene-level direction assignment
oe_mc_df <- merged %>%
  dplyr::filter(mc_sig & !is.na(me3_fold) & me3_fdr < FDR_THRESHOLD) %>%
  mutate(
    mc_direction = ifelse(mc_diff > 0, "mC Up", "mC Down"),
    me3_direction = ifelse(me3_fold > 0, "me3 Gained", "me3 Lost")
  )

cat(sprintf("  Genes with both sig mC DMR and sig me3: %d\n", nrow(oe_mc_df)))

build_oe_heatmap <- function(direction_a, direction_b, levels_a, levels_b,
                             title, x_label, y_label, predicted_pairs = list()) {
  quad_table <- table(factor(direction_a, levels = levels_a),
                      factor(direction_b, levels = levels_b))
  fisher_result <- fisher.test(quad_table)
  total_n <- sum(quad_table)
  row_totals <- rowSums(quad_table)
  col_totals <- colSums(quad_table)
  expected <- outer(row_totals, col_totals) / total_n
  enrichment <- quad_table / expected

  heatmap_df <- as.data.frame(as.table(quad_table))
  colnames(heatmap_df) <- c("row_dir", "col_dir", "Observed")
  expected_df <- as.data.frame(as.table(round(expected, 1)))
  colnames(expected_df) <- c("row_dir", "col_dir", "Expected")
  enrichment_df <- as.data.frame(as.table(round(enrichment, 2)))
  colnames(enrichment_df) <- c("row_dir", "col_dir", "Enrichment")

  heatmap_data <- heatmap_df %>%
    left_join(expected_df, by = c("row_dir", "col_dir")) %>%
    left_join(enrichment_df, by = c("row_dir", "col_dir")) %>%
    mutate(
      label = sprintf("Obs: %d\nExp: %.0f\nO/E: %.2f", Observed, Expected, Enrichment),
      row_dir = factor(row_dir, levels = levels_a),
      col_dir = factor(col_dir, levels = levels_b),
      is_predicted = FALSE
    )

  for (pair in predicted_pairs) {
    heatmap_data$is_predicted[heatmap_data$row_dir == pair[1] &
                              heatmap_data$col_dir == pair[2]] <- TRUE
  }

  sub_text <- sprintf("Fisher OR = %.2f, p = %.2e (n = %d genes)",
                      fisher_result$estimate, fisher_result$p.value, total_n)
  if (length(predicted_pairs) > 0) {
    sub_text <- paste0(sub_text, " | Black borders = predicted")
  }

  p <- ggplot(heatmap_data, aes(x = col_dir, y = row_dir, fill = Enrichment)) +
    geom_tile(color = "white", linewidth = 1.5) +
    geom_text(aes(label = label), size = 4, lineheight = 1.2)

  if (any(heatmap_data$is_predicted)) {
    p <- p + geom_tile(data = heatmap_data %>% dplyr::filter(is_predicted),
                       aes(x = col_dir, y = row_dir),
                       fill = NA, color = "black", linewidth = 1.5, linetype = "solid")
  }

  p <- p +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 1, name = "O/E",
                         limits = c(min(heatmap_data$Enrichment) * 0.9,
                                    max(heatmap_data$Enrichment) * 1.1)) +
    labs(title = title, subtitle = sub_text, x = x_label, y = y_label) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey40"),
      axis.text = element_text(size = 11, face = "bold"),
      axis.title = element_text(size = 11, face = "bold"),
      panel.grid = element_blank(),
      legend.position = "right"
    )

  list(plot = p, data = heatmap_data, fisher = fisher_result, table = quad_table)
}

if (nrow(oe_mc_df) >= 4) {
  # Predicted: if DNMT3B pathway, me3 Gained + mC Up enriched
  result_38b <- build_oe_heatmap(
    oe_mc_df$mc_direction, oe_mc_df$me3_direction,
    c("mC Up", "mC Down"), c("me3 Gained", "me3 Lost"),
    title = "H3K36me3 Direction x 5mC Direction (Gene Level)",
    x_label = "H3K36me3 Direction", y_label = "5mC DMR Direction",
    predicted_pairs = list(c("mC Up", "me3 Gained"))
  )

  save_multiformat_ggplot(result_38b$plot,
    file.path(OUTPUT_DIR, SECTION_DIR, "38b_me3_x_mc_oe_heatmap"),
    width = 8, height = 6)

  cat(sprintf("  Fisher OR = %.2f, p = %.2e\n",
              result_38b$fisher$estimate, result_38b$fisher$p.value))
} else {
  cat("  WARNING: Too few genes for mC x me3 O/E analysis\n")
}

# =============================================================================
# FIGURE 38c: O/E Heatmap — me3 Direction x hmC Direction
# =============================================================================

cat("\n--- FIGURE 38c: me3 x hmC Direction O/E Heatmap ---\n\n")

oe_hmc_df <- merged %>%
  dplyr::filter(hmc_sig & !is.na(me3_fold) & me3_fdr < FDR_THRESHOLD) %>%
  mutate(
    hmc_direction = ifelse(hmc_diff > 0, "hmC Up", "hmC Down"),
    me3_direction = ifelse(me3_fold > 0, "me3 Gained", "me3 Lost")
  )

cat(sprintf("  Genes with both sig hmC DMR and sig me3: %d\n", nrow(oe_hmc_df)))

if (nrow(oe_hmc_df) >= 4) {
  result_38c <- build_oe_heatmap(
    oe_hmc_df$hmc_direction, oe_hmc_df$me3_direction,
    c("hmC Up", "hmC Down"), c("me3 Gained", "me3 Lost"),
    title = "H3K36me3 Direction x 5hmC Direction (Gene Level)",
    x_label = "H3K36me3 Direction", y_label = "5hmC DMR Direction",
    predicted_pairs = list(c("hmC Down", "me3 Lost"))
  )

  save_multiformat_ggplot(result_38c$plot,
    file.path(OUTPUT_DIR, SECTION_DIR, "38c_me3_x_hmc_oe_heatmap"),
    width = 8, height = 6)

  cat(sprintf("  Fisher OR = %.2f, p = %.2e\n",
              result_38c$fisher$estimate, result_38c$fisher$p.value))
} else {
  cat("  WARNING: Too few genes for hmC x me3 O/E analysis\n")
}

# =============================================================================
# FIGURE 38d: Scatter — me3 Fold vs mC Diff + me3 Fold vs hmC Diff
# =============================================================================

cat("\n--- FIGURE 38d: me3 Fold vs Methylation Change Scatters ---\n\n")

scatter_df <- merged %>%
  dplyr::filter(!is.na(me3_fold)) %>%
  mutate(
    both_sig = case_when(
      mc_sig & me3_fdr < FDR_THRESHOLD ~ "Both Significant",
      mc_sig | me3_fdr < FDR_THRESHOLD ~ "One Significant",
      TRUE ~ "Neither"
    )
  )

# mC scatter
cor_mc <- cor.test(scatter_df$me3_fold, scatter_df$mc_diff, method = "spearman")
scatter_df$is_key <- scatter_df$gene %in% KEY_GENES

p_scatter_mc <- ggplot(scatter_df, aes(x = me3_fold, y = mc_diff)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  geom_point(aes(color = both_sig), alpha = 0.3, size = 0.8) +
  geom_smooth(method = "lm", color = "#D95F02", fill = "#D95F02", alpha = 0.15, linewidth = 0.8) +
  scale_color_manual(values = c("Both Significant" = "#D95F02", "One Significant" = "#999999",
                                "Neither" = "grey80"), name = "Significance") +
  {if (any(scatter_df$is_key))
    ggrepel::geom_text_repel(data = scatter_df %>% dplyr::filter(is_key),
                             aes(label = gene), size = 2.5, fontface = "italic",
                             max.overlaps = 15, segment.size = 0.2)
  } +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.3, size = 3.5,
           label = sprintf("rho = %+.3f\np = %.2e\nn = %d",
                           cor_mc$estimate, cor_mc$p.value, nrow(scatter_df))) +
  labs(title = "H3K36me3 Fold vs 5mC Change", x = "H3K36me3 Log2FC", y = "5mC Difference") +
  theme_biomodal() +
  theme(legend.position = "none")

# hmC scatter
cor_hmc <- cor.test(scatter_df$me3_fold, scatter_df$hmc_diff, method = "spearman")

p_scatter_hmc <- ggplot(scatter_df, aes(x = me3_fold, y = hmc_diff)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  geom_point(aes(color = both_sig), alpha = 0.3, size = 0.8) +
  geom_smooth(method = "lm", color = "#1B9E77", fill = "#1B9E77", alpha = 0.15, linewidth = 0.8) +
  scale_color_manual(values = c("Both Significant" = "#1B9E77", "One Significant" = "#999999",
                                "Neither" = "grey80"), name = "Significance") +
  {if (any(scatter_df$is_key))
    ggrepel::geom_text_repel(data = scatter_df %>% dplyr::filter(is_key),
                             aes(label = gene), size = 2.5, fontface = "italic",
                             max.overlaps = 15, segment.size = 0.2)
  } +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.3, size = 3.5,
           label = sprintf("rho = %+.3f\np = %.2e\nn = %d",
                           cor_hmc$estimate, cor_hmc$p.value, nrow(scatter_df))) +
  labs(title = "H3K36me3 Fold vs 5hmC Change", x = "H3K36me3 Log2FC", y = "5hmC Difference") +
  theme_biomodal() +
  theme(legend.position = "none")

p_38d <- p_scatter_mc + p_scatter_hmc +
  plot_annotation(title = "Gene-Level H3K36me3 Fold vs Methylation Changes",
                  subtitle = "Spearman correlation | Key genes labeled",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
                                plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40")))

save_multiformat_ggplot(p_38d, file.path(OUTPUT_DIR, SECTION_DIR, "38d_me3_vs_methylation_scatter"),
                        width = 14, height = 7)

cat(sprintf("  mC correlation:  rho = %+.3f, p = %.2e\n", cor_mc$estimate, cor_mc$p.value))
cat(sprintf("  hmC correlation: rho = %+.3f, p = %.2e\n", cor_hmc$estimate, cor_hmc$p.value))

# =============================================================================
# FIGURE 38e: Violin — me3 Fold at Coordinated mC-up/hmC-down Genes
# =============================================================================

cat("\n--- FIGURE 38e: me3 Fold at Coordinated Genes ---\n\n")

violin_df <- merged %>%
  dplyr::filter(!is.na(me3_fold)) %>%
  mutate(
    gene_group = case_when(
      mc_sig & hmc_sig & mc_diff > 0 & hmc_diff < 0 ~ "Coordinated\n(mC\u2191/hmC\u2193)",
      mc_sig | hmc_sig ~ "Other DMR\nGenes",
      TRUE ~ "Non-DMR\nGenes"
    ),
    gene_group = factor(gene_group, levels = c("Coordinated\n(mC\u2191/hmC\u2193)",
                                                "Other DMR\nGenes", "Non-DMR\nGenes"))
  )

group_counts <- violin_df %>% count(gene_group)
cat("  Group sizes:\n")
for (i in seq_len(nrow(group_counts))) {
  cat(sprintf("    %s: %d genes\n", gsub("\n", " ", group_counts$gene_group[i]), group_counts$n[i]))
}

# Wilcoxon: coordinated vs non-DMR
coord_vals <- violin_df$me3_fold[violin_df$gene_group == "Coordinated\n(mC\u2191/hmC\u2193)"]
nondmr_vals <- violin_df$me3_fold[violin_df$gene_group == "Non-DMR\nGenes"]
wilcox_result <- wilcox.test(coord_vals, nondmr_vals)

p_38e <- ggplot(violin_df, aes(x = gene_group, y = me3_fold, fill = gene_group)) +
  geom_violin(alpha = 0.5, scale = "width") +
  geom_boxplot(width = 0.2, outlier.shape = NA, fill = "white", alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("#D95F02", "#999999", "grey80")) +
  annotate("text", x = 2, y = max(violin_df$me3_fold, na.rm = TRUE) * 0.95,
           label = sprintf("Coord vs Non-DMR\nWilcoxon p = %.2e", wilcox_result$p.value),
           size = 3.5, fontface = "italic", lineheight = 0.9) +
  labs(
    title = "H3K36me3 Fold Change at Coordinated vs Other Genes",
    subtitle = "Coordinated = mC up + hmC down (predicted DNMT3B targets if me3 gained)",
    x = "", y = "H3K36me3 Log2 Fold Change"
  ) +
  theme_biomodal() +
  theme(legend.position = "none")

save_multiformat_ggplot(p_38e, file.path(OUTPUT_DIR, SECTION_DIR, "38e_me3_fold_coordinated_violin"),
                        width = 8, height = 7)

cat(sprintf("  Wilcoxon p = %.2e (coord vs non-DMR)\n", wilcox_result$p.value))
cat(sprintf("  Median me3 fold — Coordinated: %+.3f, Non-DMR: %+.3f\n",
            median(coord_vals, na.rm = TRUE), median(nondmr_vals, na.rm = TRUE)))

# =============================================================================
# FIGURE 38f: Bar — DMR Overlap with me3 Differential BEDs
# =============================================================================

cat("\n--- FIGURE 38f: DMR Overlap with me3 Up/Down BEDs ---\n\n")

mc_sig_gr <- dmr_to_granges(mc_dmr %>% dplyr::filter(significant))
mc_hyper_gr <- mc_sig_gr[mc_sig_gr$mod_difference > 0]
mc_hypo_gr <- mc_sig_gr[mc_sig_gr$mod_difference <= 0]

overlap_hyper_up <- sum(countOverlaps(mc_hyper_gr, me3_up_gr) > 0)
overlap_hyper_down <- sum(countOverlaps(mc_hyper_gr, me3_down_gr) > 0)
overlap_hypo_up <- sum(countOverlaps(mc_hypo_gr, me3_up_gr) > 0)
overlap_hypo_down <- sum(countOverlaps(mc_hypo_gr, me3_down_gr) > 0)

overlap_df <- data.frame(
  dmr_type = rep(c("mC Hyper", "mC Hypo"), each = 2),
  me3_type = rep(c("me3 Gained", "me3 Lost"), 2),
  count = c(overlap_hyper_up, overlap_hyper_down, overlap_hypo_up, overlap_hypo_down)
)

# Add totals for percentage
overlap_df <- overlap_df %>%
  group_by(dmr_type) %>%
  mutate(total = c(length(mc_hyper_gr), length(mc_hyper_gr),
                   length(mc_hypo_gr), length(mc_hypo_gr))[row_number()],
         pct = 100 * count / total) %>%
  ungroup()

# Fix total assignment
overlap_df$total <- ifelse(overlap_df$dmr_type == "mC Hyper", length(mc_hyper_gr), length(mc_hypo_gr))
overlap_df$pct <- 100 * overlap_df$count / overlap_df$total

# Fisher's test on the 2x2
fisher_overlap <- fisher.test(matrix(overlap_df$count, nrow = 2, byrow = TRUE))

p_38f <- ggplot(overlap_df, aes(x = dmr_type, y = pct, fill = me3_type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6,
           color = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", pct, count)),
            position = position_dodge(width = 0.7), vjust = -0.3,
            size = 3, lineheight = 0.85) +
  scale_fill_manual(values = c("me3 Gained" = "#D95F02", "me3 Lost" = "#1B9E77"),
                    name = "H3K36me3 Direction") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(
    title = "mC DMR Overlap with H3K36me3 Differential Peaks",
    subtitle = sprintf("Fisher OR = %.2f, p = %.2e | Coordinate-level overlap",
                       fisher_overlap$estimate, fisher_overlap$p.value),
    x = "mC DMR Direction", y = "% of DMRs Overlapping"
  ) +
  theme_biomodal()

save_multiformat_ggplot(p_38f, file.path(OUTPUT_DIR, SECTION_DIR, "38f_dmr_me3_bed_overlap"),
                        width = 8, height = 7)

cat(sprintf("  Fisher OR = %.2f, p = %.2e\n", fisher_overlap$estimate, fisher_overlap$p.value))

# =============================================================================
# FIGURE 38g: Chromosome Distribution of Significant me3 Peaks
# =============================================================================

cat("\n--- FIGURE 38g: Chromosome Distribution ---\n\n")

chr_dist <- me3_all %>%
  dplyr::filter(sig) %>%
  mutate(
    direction = ifelse(Fold > 0, "Gained", "Lost"),
    chr_clean = gsub("^chr", "", Chr)
  )

# Order chromosomes naturally
chr_order <- c(as.character(1:19), "X", "Y")
chr_dist$chr_clean <- factor(chr_dist$chr_clean, levels = chr_order)
chr_dist <- chr_dist %>% dplyr::filter(!is.na(chr_clean))

chr_summary <- chr_dist %>%
  count(chr_clean, direction) %>%
  mutate(direction = factor(direction, levels = c("Gained", "Lost")))

# Highlight chr13 (histone gene cluster)
chr_summary$highlight <- chr_summary$chr_clean == "13"

p_38g <- ggplot(chr_summary, aes(x = chr_clean, y = n, fill = direction)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7, color = "black", linewidth = 0.2) +
  # Highlight chr13
  annotate("rect", xmin = which(chr_order == "13") - 0.5,
           xmax = which(chr_order == "13") + 0.5,
           ymin = -Inf, ymax = Inf, alpha = 0.15, fill = "red") +
  annotate("text", x = which(chr_order == "13"), y = Inf, vjust = -0.5,
           label = "Histone cluster", size = 3, fontface = "italic", color = "red3") +
  scale_fill_manual(values = c("Gained" = "#D95F02", "Lost" = "#1B9E77"), name = "Direction") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Chromosome Distribution of Significant H3K36me3 Peaks",
    subtitle = sprintf("%d total differential peaks (FDR < %.2f)", sum(me3_all$sig), FDR_THRESHOLD),
    x = "Chromosome", y = "Number of Differential Peaks"
  ) +
  theme_biomodal() +
  theme(axis.text.x = element_text(size = 9))

save_multiformat_ggplot(p_38g, file.path(OUTPUT_DIR, SECTION_DIR, "38g_chromosome_distribution"),
                        width = 12, height = 6)

# chr13 stats
n_chr13 <- sum(chr_dist$chr_clean == "13", na.rm = TRUE)
n_total_sig <- nrow(chr_dist)
cat(sprintf("  chr13: %d / %d sig peaks (%.1f%%)\n", n_chr13, n_total_sig, 100 * n_chr13 / n_total_sig))

# =============================================================================
# FIGURE 38h: me3 Fold by Chromatin State
# =============================================================================

cat("\n--- FIGURE 38h: me3 Fold by Chromatin State ---\n\n")

state_df <- merged %>%
  dplyr::filter(!is.na(me3_fold) & !is.na(chromatin_state)) %>%
  mutate(chromatin_state = factor(chromatin_state, levels = CHROMATIN_STATE_ORDER))

state_counts <- state_df %>% count(chromatin_state)
kruskal_result <- kruskal.test(me3_fold ~ chromatin_state, data = state_df)

p_38h <- ggplot(state_df, aes(x = chromatin_state, y = me3_fold, fill = chromatin_state)) +
  geom_violin(alpha = 0.5, scale = "width") +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = CHROMATIN_STATE_COLORS) +
  # Annotate sample sizes
  geom_text(data = state_counts, aes(x = chromatin_state, y = -Inf, label = sprintf("n=%d", n)),
            inherit.aes = FALSE, vjust = -0.5, size = 2.8, color = "grey40") +
  labs(
    title = "H3K36me3 Fold Change by Chromatin State",
    subtitle = sprintf("Kruskal-Wallis p = %.2e (across %d chromatin states)",
                       kruskal_result$p.value, length(unique(state_df$chromatin_state))),
    x = "", y = "H3K36me3 Log2 Fold Change"
  ) +
  theme_biomodal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1, size = 10))

save_multiformat_ggplot(p_38h, file.path(OUTPUT_DIR, SECTION_DIR, "38h_me3_fold_by_chromatin_state"),
                        width = 10, height = 7)

cat(sprintf("  Kruskal-Wallis p = %.2e\n", kruskal_result$p.value))

# =============================================================================
# EXPORT TABLES
# =============================================================================

cat("\n--- Exporting Tables ---\n\n")

# Gene-level summary
export_gene <- me3_gene %>%
  left_join(ratio_data %>% dplyr::select(gene, mc_diff, mc_q, mc_sig, hmc_diff, hmc_q, hmc_sig,
                                          chromatin_state, dmr_status),
            by = "gene")
write.table(export_gene, file.path(TABLES_DIR, "h3k36me3_gene_level_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved h3k36me3_gene_level_summary.tsv (%d rows)\n", nrow(export_gene)))

# Direction concordance summary
if (exists("result_38b")) {
  conc_summary <- data.frame(
    comparison = c("me3 x mC", "me3 x hmC"),
    fisher_or = c(result_38b$fisher$estimate,
                  if (exists("result_38c")) result_38c$fisher$estimate else NA),
    fisher_p = c(result_38b$fisher$p.value,
                 if (exists("result_38c")) result_38c$fisher$p.value else NA),
    n_genes = c(nrow(oe_mc_df),
                if (exists("oe_hmc_df")) nrow(oe_hmc_df) else NA),
    spearman_rho = c(cor_mc$estimate, cor_hmc$estimate),
    spearman_p = c(cor_mc$p.value, cor_hmc$p.value)
  )
  write.table(conc_summary, file.path(TABLES_DIR, "h3k36me3_direction_concordance.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  Saved h3k36me3_direction_concordance.tsv (%d rows)\n", nrow(conc_summary)))
}

cat("\n================================================================================\n")
cat("SECTION 38 COMPLETE\n")
cat("================================================================================\n")
