# biomodal/downstream/scripts/viz_sections/section_20_coordinated_rnaseq.R
# Section 20: Coordinated Gene RNA-seq Expression Integration
# Standalone script - sources shared config for all dependencies and data
#
# Investigates gene expression outcomes for coordinated mC up/hmC down genes.
# Methylation gain does NOT necessarily mean repression -- we check RNA-seq
# to quantify what fraction are repressed, activated, or unchanged.
#
# Figures:
#   20a: Stacked bar -- expression direction breakdown (Tier 1 & 2)
#   20b: Scatter -- combined methylation effect vs log2FC
#   20c: Violin -- log2FC distributions (coordinated vs other DMR genes)
#   20d: Heatmap -- 2x2 enrichment (mC direction x expression direction)
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_20_coordinated_rnaseq.R

# Run from downstream/ directory
source("scripts/viz_sections/_shared_config.R")

suppressPackageStartupMessages(library(readxl))

# =============================================================================
# SECTION 20 CONFIGURATION
# =============================================================================

RNA_SEQ_FILE <- "/Users/zakiralibhai/Documents/GitHub/mariner_hi-c/tads/adult_timepoint_rna-seq-BAP1_WT_KO_v2_Results.xlsx"

# Two-tier expression classification thresholds
TIER1_PADJ <- 0.05
TIER2_PADJ <- 0.05
TIER2_LFC  <- 0.3

# Expression direction colors
EXPR_COLORS <- c("Up" = "#D7191C", "Down" = "#2C7BB6", "Unchanged" = "#999999")

# =============================================================================
# VALIDATE INPUTS
# =============================================================================

cat("================================================================================\n")
cat("SECTION 20: COORDINATED GENE RNA-seq EXPRESSION INTEGRATION\n")
cat("================================================================================\n\n")

stopifnot(file.exists(RNA_SEQ_FILE))
cat(sprintf("  RNA-seq file: %s\n", RNA_SEQ_FILE))

# =============================================================================
# RECOMPUTE COORDINATED CHANGES (independent from Section 5)
# =============================================================================

cat("\nRecomputing coordinated mC/hmC changes...\n")

mc_sig <- mc_dmr %>%
  dplyr::filter(significant) %>%
  dplyr::select(gene, mc_diff = mod_difference, mc_q = dmr_qvalue)

hmc_sig <- hmc_dmr %>%
  dplyr::filter(significant) %>%
  dplyr::select(gene, hmc_diff = mod_difference, hmc_q = dmr_qvalue)

coordinated <- inner_join(mc_sig, hmc_sig, by = "gene") %>%
  mutate(
    quadrant = case_when(
      mc_diff < 0 & hmc_diff < 0 ~ "Q1: mC dn/hmC dn",
      mc_diff < 0 & hmc_diff > 0 ~ "Q2: mC dn/hmC up",
      mc_diff > 0 & hmc_diff > 0 ~ "Q3: mC up/hmC up",
      mc_diff > 0 & hmc_diff < 0 ~ "Q4: mC up/hmC dn"
    ),
    coordinated_pattern = (mc_diff > 0 & hmc_diff < 0),
    combined_effect = abs(mc_diff) + abs(hmc_diff)
  ) %>%
  arrange(desc(combined_effect))

n_coord_total <- nrow(coordinated)
n_q4 <- sum(coordinated$coordinated_pattern)
cat(sprintf("  %d genes significant in both mC and hmC\n", n_coord_total))
cat(sprintf("  %d genes show coordinated mC up/hmC dn pattern (%.1f%%)\n",
            n_q4, 100 * n_q4 / n_coord_total))

# =============================================================================
# LOAD RNA-seq DATA (all genes, not just DEGs)
# =============================================================================

cat("\nLoading RNA-seq data...\n")

rna_all <- read_excel(RNA_SEQ_FILE) %>%
  dplyr::select(gene = ensembl_gene_id, log2FC = log2FoldChange, padj, baseMean) %>%
  dplyr::filter(!is.na(gene) & gene != "")

cat(sprintf("  Loaded %d genes from RNA-seq\n", nrow(rna_all)))
cat(sprintf("  Significant (padj < 0.05): %d\n", sum(rna_all$padj < 0.05, na.rm = TRUE)))
cat(sprintf("  DEG (padj < 0.05 & |log2FC| > 0.3): %d\n",
            sum(rna_all$padj < 0.05 & abs(rna_all$log2FC) > 0.3, na.rm = TRUE)))

# =============================================================================
# JOIN AND CLASSIFY EXPRESSION OUTCOMES
# =============================================================================

cat("\nJoining coordinated genes with RNA-seq...\n")

# Q4 coordinated genes joined with RNA-seq
q4_genes <- coordinated %>% dplyr::filter(coordinated_pattern)
q4_rna <- q4_genes %>% left_join(rna_all, by = "gene")
n_q4_matched <- sum(!is.na(q4_rna$log2FC))
cat(sprintf("  Q4 genes with RNA-seq data: %d / %d (%.1f%%)\n",
            n_q4_matched, nrow(q4_rna), 100 * n_q4_matched / nrow(q4_rna)))

# Other significant mC DMR genes (not in Q4)
other_mc_genes <- mc_dmr %>%
  dplyr::filter(significant & !(gene %in% q4_genes$gene)) %>%
  dplyr::select(gene, mc_diff = mod_difference) %>%
  distinct(gene, .keep_all = TRUE) %>%
  left_join(rna_all, by = "gene")

cat(sprintf("  Other mC DMR genes with RNA-seq data: %d\n",
            sum(!is.na(other_mc_genes$log2FC))))

# All coordinated genes (all quadrants) joined with RNA-seq for output
all_coord_rna <- coordinated %>%
  left_join(rna_all, by = "gene")

# =============================================================================
# HELPER: Build expression breakdown for a gene set
# =============================================================================

build_expr_breakdown <- function(df, group_label) {
  df_valid <- df %>% dplyr::filter(!is.na(padj) & !is.na(log2FC))
  n <- nrow(df_valid)
  if (n == 0) return(NULL)

  # Tier 1: padj only
  t1_up <- sum(df_valid$padj < TIER1_PADJ & df_valid$log2FC > 0)
  t1_dn <- sum(df_valid$padj < TIER1_PADJ & df_valid$log2FC < 0)
  t1_ns <- n - t1_up - t1_dn

  # Tier 2: padj + |log2FC|
  t2_up <- sum(df_valid$padj < TIER2_PADJ & df_valid$log2FC > TIER2_LFC)
  t2_dn <- sum(df_valid$padj < TIER2_PADJ & df_valid$log2FC < -TIER2_LFC)
  t2_ns <- n - t2_up - t2_dn

  data.frame(
    group = group_label,
    tier = rep(c("Tier 1: padj < 0.05",
                 "Tier 2: padj < 0.05 & |log2FC| > 0.3"), each = 3),
    direction = rep(c("Up", "Down", "Unchanged"), 2),
    count = c(t1_up, t1_dn, t1_ns, t2_up, t2_dn, t2_ns),
    total = n,
    stringsAsFactors = FALSE
  ) %>%
    mutate(pct = 100 * count / total)
}

# =============================================================================
# FIGURE 20a: Stacked Bar Chart (Expression Direction Breakdown)
# =============================================================================

cat("\nCreating Figure 20a: Expression direction breakdown...\n")

bar_q4    <- build_expr_breakdown(q4_rna, "Coordinated\n(mC up/hmC dn)")
bar_other <- build_expr_breakdown(other_mc_genes, "Other mC\nDMR Genes")
bar_all   <- build_expr_breakdown(rna_all, "All Genes")

bar_df <- bind_rows(bar_q4, bar_other, bar_all) %>%
  mutate(
    group = factor(group, levels = c("Coordinated\n(mC up/hmC dn)",
                                     "Other mC\nDMR Genes", "All Genes")),
    direction = factor(direction, levels = c("Up", "Down", "Unchanged"))
  )

# Fisher's exact tests: coordinated vs other (Up vs Down among significant)
fisher_tier1 <- fisher.test(matrix(c(
  bar_q4$count[bar_q4$tier == "Tier 1: padj < 0.05" & bar_q4$direction == "Up"],
  bar_q4$count[bar_q4$tier == "Tier 1: padj < 0.05" & bar_q4$direction == "Down"],
  bar_other$count[bar_other$tier == "Tier 1: padj < 0.05" & bar_other$direction == "Up"],
  bar_other$count[bar_other$tier == "Tier 1: padj < 0.05" & bar_other$direction == "Down"]
), nrow = 2, byrow = TRUE,
dimnames = list(c("Coordinated", "Other"), c("Up", "Down"))))

fisher_tier2 <- fisher.test(matrix(c(
  bar_q4$count[bar_q4$tier == "Tier 2: padj < 0.05 & |log2FC| > 0.3" & bar_q4$direction == "Up"],
  bar_q4$count[bar_q4$tier == "Tier 2: padj < 0.05 & |log2FC| > 0.3" & bar_q4$direction == "Down"],
  bar_other$count[bar_other$tier == "Tier 2: padj < 0.05 & |log2FC| > 0.3" & bar_other$direction == "Up"],
  bar_other$count[bar_other$tier == "Tier 2: padj < 0.05 & |log2FC| > 0.3" & bar_other$direction == "Down"]
), nrow = 2, byrow = TRUE,
dimnames = list(c("Coordinated", "Other"), c("Up", "Down"))))

cat(sprintf("  Fisher's (Tier 1, coord vs other): OR = %.2f, p = %.2e\n",
            fisher_tier1$estimate, fisher_tier1$p.value))
cat(sprintf("  Fisher's (Tier 2, coord vs other): OR = %.2f, p = %.2e\n",
            fisher_tier2$estimate, fisher_tier2$p.value))

# Build annotation for facets
fisher_labels <- data.frame(
  tier = c("Tier 1: padj < 0.05", "Tier 2: padj < 0.05 & |log2FC| > 0.3"),
  label = c(
    sprintf("Fisher's OR = %.2f, p = %.2e\n(Coordinated vs Other, Up/Down ratio)",
            fisher_tier1$estimate, fisher_tier1$p.value),
    sprintf("Fisher's OR = %.2f, p = %.2e\n(Coordinated vs Other, Up/Down ratio)",
            fisher_tier2$estimate, fisher_tier2$p.value)
  ),
  stringsAsFactors = FALSE
)

p_20a <- ggplot(bar_df, aes(x = group, y = pct, fill = direction)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7,
           color = "black", linewidth = 0.2) +
  geom_text(aes(label = ifelse(pct > 3, sprintf("%.1f%%\n(n=%d)", pct, count), "")),
            position = position_stack(vjust = 0.5), size = 3) +
  facet_wrap(~tier, ncol = 2) +
  scale_fill_manual(values = EXPR_COLORS, name = "Expression Direction") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 105)) +
  geom_text(data = fisher_labels,
            aes(x = 1.5, y = 103, label = label),
            inherit.aes = FALSE, size = 2.8, color = "grey30",
            vjust = 1, hjust = 0.5) +
  labs(
    title = "Gene Expression Outcomes for Coordinated mC up/hmC dn Genes",
    subtitle = sprintf("n = %d coordinated | %d other mC DMR | %d all genes",
                       nrow(q4_rna %>% dplyr::filter(!is.na(log2FC))),
                       nrow(other_mc_genes %>% dplyr::filter(!is.na(log2FC))),
                       nrow(rna_all %>% dplyr::filter(!is.na(padj)))),
    x = "", y = "Percentage of Genes (%)"
  ) +
  theme_biomodal() +
  theme(legend.position = "top",
        axis.text.x = element_text(size = 10))

save_multiformat_ggplot(p_20a, file.path(OUTPUT_DIR, "20a_coordinated_expression_breakdown"),
                        width = 13, height = 8)

# =============================================================================
# FIGURE 20b: Scatter -- Combined Methylation Effect vs log2FC
# =============================================================================

cat("\nCreating Figure 20b: Combined methylation effect vs log2FC...\n")

scatter_data <- q4_rna %>%
  dplyr::filter(!is.na(log2FC)) %>%
  mutate(
    combined_pct = combined_effect * 100,
    sig_category = case_when(
      padj < 0.05 & log2FC > TIER2_LFC ~ "Up (DEG)",
      padj < 0.05 & log2FC < -TIER2_LFC ~ "Down (DEG)",
      padj < 0.05 ~ "Significant (small effect)",
      TRUE ~ "Not Significant"
    ),
    label = ifelse(gene %in% KEY_GENES, gene, "")
  )

# Spearman correlation
cor_test <- cor.test(scatter_data$combined_pct, scatter_data$log2FC,
                     method = "spearman")
cat(sprintf("  Spearman rho = %.3f, p = %.2e\n", cor_test$estimate, cor_test$p.value))

scatter_colors <- c("Up (DEG)" = "#D7191C", "Down (DEG)" = "#2C7BB6",
                     "Significant (small effect)" = "#FDB863",
                     "Not Significant" = "grey70")

p_20b <- ggplot(scatter_data, aes(x = combined_pct, y = log2FC)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_point(aes(color = sig_category), alpha = 0.5, size = 1.8) +
  geom_smooth(method = "lm", color = "black", linewidth = 0.7, se = TRUE, alpha = 0.15) +
  geom_text_repel(aes(label = label), size = 3.5, max.overlaps = 15,
                  fontface = "italic", color = "grey20",
                  segment.color = "grey60", segment.size = 0.3) +
  scale_color_manual(values = scatter_colors, name = "Expression Status") +
  labs(
    title = "Combined Methylation Effect vs Gene Expression",
    subtitle = sprintf("Coordinated mC up/hmC dn genes (n = %d) | Spearman rho = %.3f, p = %.2e",
                       nrow(scatter_data), cor_test$estimate, cor_test$p.value),
    x = "Combined Methylation Effect (|mC| + |hmC| change, %)",
    y = "RNA-seq log2FC (Mutant / Control)"
  ) +
  theme_biomodal() +
  theme(legend.position = "right")

save_multiformat_ggplot(p_20b, file.path(OUTPUT_DIR, "20b_methylation_vs_expression_scatter"),
                        width = 11, height = 9)

# =============================================================================
# FIGURE 20c: Violin -- log2FC Distributions (Coordinated vs Other)
# =============================================================================

cat("\nCreating Figure 20c: log2FC distribution comparison...\n")

# Build comparison data: Q4 coordinated vs all other significant mC DMR genes
violin_data <- bind_rows(
  q4_rna %>%
    dplyr::filter(!is.na(log2FC)) %>%
    dplyr::select(gene, log2FC) %>%
    mutate(category = "Coordinated\n(mC up/hmC dn)"),
  other_mc_genes %>%
    dplyr::filter(!is.na(log2FC)) %>%
    dplyr::select(gene, log2FC) %>%
    mutate(category = "Other Significant\nmC DMR Genes")
) %>%
  mutate(category = factor(category, levels = c("Coordinated\n(mC up/hmC dn)",
                                                 "Other Significant\nmC DMR Genes")))

wilcox_lfc <- wilcox.test(
  log2FC ~ category, data = violin_data
)
cat(sprintf("  Wilcoxon log2FC: p = %.2e\n", wilcox_lfc$p.value))
cat(sprintf("  Median log2FC coordinated: %.4f\n",
            median(violin_data$log2FC[violin_data$category == "Coordinated\n(mC up/hmC dn)"])))
cat(sprintf("  Median log2FC other: %.4f\n",
            median(violin_data$log2FC[violin_data$category == "Other Significant\nmC DMR Genes"])))

p_20c <- ggplot(violin_data, aes(x = category, y = log2FC, fill = category)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_violin(alpha = 0.6, scale = "width") +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.7) +
  scale_fill_manual(values = c("Coordinated\n(mC up/hmC dn)" = "#8E44AD",
                                "Other Significant\nmC DMR Genes" = "grey70")) +
  annotate("text",
           x = 1.5,
           y = max(violin_data$log2FC, na.rm = TRUE) * 0.95,
           label = sprintf("Wilcoxon p = %.2e", wilcox_lfc$p.value),
           size = 4, fontface = "italic") +
  labs(
    title = "RNA-seq log2FC: Coordinated vs Other DMR Genes",
    subtitle = sprintf("Does mC up + hmC dn predict expression repression? (n = %d vs %d)",
                       sum(violin_data$category == "Coordinated\n(mC up/hmC dn)"),
                       sum(violin_data$category == "Other Significant\nmC DMR Genes")),
    x = "", y = "RNA-seq log2FC (Mutant / Control)"
  ) +
  theme_biomodal() +
  theme(legend.position = "none")

save_multiformat_ggplot(p_20c, file.path(OUTPUT_DIR, "20c_log2fc_violin_comparison"),
                        width = 8, height = 8)

# =============================================================================
# FIGURE 20d: 2x2 Enrichment Heatmap (mC Direction x Expression Direction)
# =============================================================================

cat("\nCreating Figure 20d: mC direction x expression direction heatmap...\n")

# Use ALL significant mC DMR genes with RNA-seq data (not just coordinated)
mc_with_rna <- mc_dmr %>%
  dplyr::filter(significant) %>%
  dplyr::select(gene, mc_diff = mod_difference) %>%
  distinct(gene, .keep_all = TRUE) %>%
  inner_join(rna_all, by = "gene") %>%
  dplyr::filter(padj < TIER1_PADJ)  # Only genes significant in BOTH mC and expression

mc_with_rna <- mc_with_rna %>%
  mutate(
    mc_direction = ifelse(mc_diff > 0, "mC Up", "mC Down"),
    expr_direction = ifelse(log2FC > 0, "Expr Up", "Expr Down")
  )

cat(sprintf("  Genes significant in both mC and expression: %d\n", nrow(mc_with_rna)))

# Contingency table
heatmap_tab <- table(mc_with_rna$mc_direction, mc_with_rna$expr_direction)
cat("  Contingency table:\n")
print(heatmap_tab)

fisher_heatmap <- fisher.test(heatmap_tab)
cat(sprintf("  Fisher's exact test: OR = %.2f, p = %.2e\n",
            fisher_heatmap$estimate, fisher_heatmap$p.value))

# Compute expected counts and enrichment
total_n <- sum(heatmap_tab)
row_totals <- rowSums(heatmap_tab)
col_totals <- colSums(heatmap_tab)
expected_mat <- outer(row_totals, col_totals) / total_n
enrichment_mat <- heatmap_tab / expected_mat

heatmap_data <- as.data.frame(as.table(heatmap_tab))
colnames(heatmap_data) <- c("mc_direction", "expr_direction", "Observed")

expected_flat <- as.data.frame(as.table(round(expected_mat, 1)))
colnames(expected_flat) <- c("mc_direction", "expr_direction", "Expected")

enrichment_flat <- as.data.frame(as.table(round(enrichment_mat, 2)))
colnames(enrichment_flat) <- c("mc_direction", "expr_direction", "Enrichment")

heatmap_data <- heatmap_data %>%
  left_join(expected_flat, by = c("mc_direction", "expr_direction")) %>%
  left_join(enrichment_flat, by = c("mc_direction", "expr_direction")) %>%
  mutate(
    label = sprintf("Obs: %d\nExp: %.0f\nOR: %.2f", Observed, Expected, Enrichment),
    mc_direction = factor(mc_direction, levels = c("mC Up", "mC Down")),
    expr_direction = factor(expr_direction, levels = c("Expr Up", "Expr Down")),
    is_predicted = (mc_direction == "mC Up" & expr_direction == "Expr Down") |
                   (mc_direction == "mC Down" & expr_direction == "Expr Up")
  )

p_20d <- ggplot(heatmap_data, aes(x = expr_direction, y = mc_direction, fill = Enrichment)) +
  geom_tile(color = "white", linewidth = 1.5) +
  geom_text(aes(label = label), size = 4.5, lineheight = 1.2) +
  geom_tile(data = heatmap_data %>% dplyr::filter(is_predicted),
            aes(x = expr_direction, y = mc_direction),
            fill = NA, color = "black", linewidth = 1.5, linetype = "solid") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 1, name = "Enrichment\n(Obs/Exp)",
                       limits = c(min(heatmap_data$Enrichment) * 0.9,
                                  max(heatmap_data$Enrichment) * 1.1)) +
  labs(
    title = "5mC Direction x Expression Direction",
    subtitle = sprintf("Fisher's OR = %.2f, p = %.2e | Black borders = predicted quadrants\n(among %d genes significant in both mC and expression)",
                       fisher_heatmap$estimate, fisher_heatmap$p.value, nrow(mc_with_rna)),
    x = "Expression Direction (RNA-seq)",
    y = "5mC DMR Direction"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"),
    axis.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    panel.grid = element_blank(),
    legend.position = "right"
  )

save_multiformat_ggplot(p_20d, file.path(OUTPUT_DIR, "20d_mc_expression_heatmap"),
                        width = 9, height = 7)

# =============================================================================
# OUTPUT TABLE
# =============================================================================

cat("\nSaving output tables...\n")

# Per-gene table for all coordinated genes (all quadrants) with RNA-seq
export_df <- all_coord_rna %>%
  mutate(
    expr_tier1 = case_when(
      is.na(padj) | is.na(log2FC) ~ "No RNA-seq Data",
      padj < TIER1_PADJ & log2FC > 0 ~ "Up",
      padj < TIER1_PADJ & log2FC < 0 ~ "Down",
      TRUE ~ "Not Significant"
    ),
    expr_tier2 = case_when(
      is.na(padj) | is.na(log2FC) ~ "No RNA-seq Data",
      padj < TIER2_PADJ & log2FC > TIER2_LFC ~ "Up",
      padj < TIER2_PADJ & log2FC < -TIER2_LFC ~ "Down",
      TRUE ~ "Not DEG"
    )
  ) %>%
  dplyr::select(gene, mc_diff, hmc_diff, combined_effect, quadrant,
                log2FC, padj, baseMean, expr_tier1, expr_tier2) %>%
  arrange(desc(combined_effect))

write.table(export_df, file.path(TABLES_DIR, "coordinated_rnaseq_expression.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved: coordinated_rnaseq_expression.tsv (%d rows)\n", nrow(export_df)))

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("SECTION 20 SUMMARY\n")
cat("================================================================================\n")

# Q4 expression breakdown
q4_valid <- q4_rna %>% dplyr::filter(!is.na(padj) & !is.na(log2FC))
cat(sprintf("\nCoordinated mC up/hmC dn genes (Q4): %d total, %d with RNA-seq\n",
            nrow(q4_rna), nrow(q4_valid)))

cat("\n--- Tier 1 (padj < 0.05) ---\n")
t1_up <- sum(q4_valid$padj < TIER1_PADJ & q4_valid$log2FC > 0)
t1_dn <- sum(q4_valid$padj < TIER1_PADJ & q4_valid$log2FC < 0)
t1_ns <- nrow(q4_valid) - t1_up - t1_dn
cat(sprintf("  Up:    %d (%.1f%%)\n", t1_up, 100 * t1_up / nrow(q4_valid)))
cat(sprintf("  Down:  %d (%.1f%%)\n", t1_dn, 100 * t1_dn / nrow(q4_valid)))
cat(sprintf("  NS:    %d (%.1f%%)\n", t1_ns, 100 * t1_ns / nrow(q4_valid)))
cat(sprintf("  Up/Down ratio: %.2f\n", t1_up / max(t1_dn, 1)))

cat("\n--- Tier 2 (padj < 0.05 & |log2FC| > 0.3) ---\n")
t2_up <- sum(q4_valid$padj < TIER2_PADJ & q4_valid$log2FC > TIER2_LFC)
t2_dn <- sum(q4_valid$padj < TIER2_PADJ & q4_valid$log2FC < -TIER2_LFC)
t2_ns <- nrow(q4_valid) - t2_up - t2_dn
cat(sprintf("  Up:    %d (%.1f%%)\n", t2_up, 100 * t2_up / nrow(q4_valid)))
cat(sprintf("  Down:  %d (%.1f%%)\n", t2_dn, 100 * t2_dn / nrow(q4_valid)))
cat(sprintf("  NS:    %d (%.1f%%)\n", t2_ns, 100 * t2_ns / nrow(q4_valid)))
cat(sprintf("  Up/Down ratio: %.2f\n", t2_up / max(t2_dn, 1)))

cat("\n--- Statistical Tests ---\n")
cat(sprintf("  Fisher's Tier 1 (coord vs other, Up/Down): OR = %.2f, p = %.2e\n",
            fisher_tier1$estimate, fisher_tier1$p.value))
cat(sprintf("  Fisher's Tier 2 (coord vs other, Up/Down): OR = %.2f, p = %.2e\n",
            fisher_tier2$estimate, fisher_tier2$p.value))
cat(sprintf("  Wilcoxon log2FC (coord vs other): p = %.2e\n", wilcox_lfc$p.value))
cat(sprintf("  Spearman combined effect vs log2FC: rho = %.3f, p = %.2e\n",
            cor_test$estimate, cor_test$p.value))
cat(sprintf("  Fisher's mC x expression heatmap: OR = %.2f, p = %.2e\n",
            fisher_heatmap$estimate, fisher_heatmap$p.value))

cat("\n--- Biological Interpretation ---\n")
if (t1_dn > t1_up) {
  cat("  mC up/hmC dn genes are biased toward REPRESSION (more Down than Up)\n")
  cat("  Consistent with: methylation gain leading to transcriptional silencing\n")
} else if (t1_up > t1_dn) {
  cat("  mC up/hmC dn genes are biased toward ACTIVATION (more Up than Down)\n")
  cat("  Suggests: methylation gain at gene bodies may mark active transcription\n")
} else {
  cat("  No clear directional bias in expression of coordinated genes\n")
}

cat("\nSection 20 complete.\n\n")
