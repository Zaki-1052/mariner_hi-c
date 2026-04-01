# biomodal/downstream/scripts/viz_sections/section_40_h3k36me2_me3_combined.R
# H3K36me2/me3 Combined Analysis and Conversion Dynamics
#
# H3K36me2 and me3 are interconverted: SETD2 converts me2 to me3 in actively
# transcribed gene bodies. This section analyzes both marks together to understand
# the H3K36 methylation landscape in BAP1-KO and its relationship to DNA methylation.
#
# Data: Gene-level summaries from sections 38 and 39 (or re-derived inline)
# Figures: 40a-40f (6 total)
#
# Usage:
#   cd downstream/
#   Rscript scripts/viz_sections/section_40_h3k36me2_me3_combined.R

source("scripts/viz_sections/_shared_config.R")
suppressPackageStartupMessages({
  library(ChIPseeker)
  library(ggVennDiagram)
})

SECTION_DIR <- "40_h3k36me2_me3_combined"
FDR_THRESHOLD <- 0.05

dir.create(file.path(OUTPUT_DIR, SECTION_DIR), recursive = TRUE, showWarnings = FALSE)

cat("\n")
cat("================================================================================\n")
cat("SECTION 40: H3K36me2/me3 Combined Analysis\n")
cat("================================================================================\n\n")

# =============================================================================
# DATA LOADING — Gene-level tables from sections 38/39 or re-derive
# =============================================================================

RATIO_TABLE <- file.path(TABLES_DIR, "demethylation_ratio_all_genes.tsv")
ME3_TABLE <- file.path(TABLES_DIR, "h3k36me3_gene_level_summary.tsv")
ME2_TABLE <- file.path(TABLES_DIR, "h3k36me2_gene_level_summary.tsv")

stopifnot(file.exists(RATIO_TABLE))
ratio_data <- read.table(RATIO_TABLE, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat(sprintf("  Ratio table: %d genes\n", nrow(ratio_data)))

# Load me3 gene-level (from section 38 or re-derive)
if (file.exists(ME3_TABLE)) {
  me3_gene <- read.table(ME3_TABLE, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  cat(sprintf("  me3 gene table loaded: %d genes\n", nrow(me3_gene)))
} else {
  cat("  me3 gene table not found — re-deriving from master file...\n")
  me3_all <- load_diffbind_flex(H3K36ME3_FILES$master_annotated, "H3K36me3")
  me3_gene <- me3_all %>%
    dplyr::filter(!is.na(SYMBOL) & SYMBOL != "") %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::summarise(
      me3_fold = Fold[which.min(abs(distanceToTSS))],
      me3_fdr  = FDR[which.min(abs(distanceToTSS))],
      me3_n_peaks = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::rename(gene = SYMBOL)
  cat(sprintf("  Derived %d genes\n", nrow(me3_gene)))
}

# Load me2 gene-level (from section 39 or re-derive)
if (file.exists(ME2_TABLE)) {
  me2_gene <- read.table(ME2_TABLE, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  cat(sprintf("  me2 gene table loaded: %d genes\n", nrow(me2_gene)))
} else {
  cat("  me2 gene table not found — re-deriving from master file...\n")
  me2_all <- load_diffbind_flex(H3K36ME2_FILES$master, "H3K36me2")
  me2_sig <- me2_all %>% dplyr::filter(FDR < FDR_THRESHOLD)

  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  me2_sig_gr <- GRanges(seqnames = me2_sig$Chr,
                         ranges = IRanges(start = me2_sig$Start, end = me2_sig$End),
                         Fold = me2_sig$Fold, FDR = me2_sig$FDR)
  me2_anno <- suppressMessages(annotatePeak(me2_sig_gr, TxDb = txdb, annoDb = "org.Mm.eg.db", level = "gene"))
  me2_anno_df <- as.data.frame(me2_anno)

  me2_gene <- me2_anno_df %>%
    dplyr::filter(!is.na(SYMBOL) & SYMBOL != "") %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::summarise(
      me2_fold = Fold[which.min(abs(distanceToTSS))],
      me2_fdr  = FDR[which.min(abs(distanceToTSS))],
      me2_n_peaks = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::rename(gene = SYMBOL)
  cat(sprintf("  Derived %d genes\n", nrow(me2_gene)))
}

# Load existing DiffBind marks for expanded correlation
DIFFBIND_TABLES <- list(
  atac   = file.path(TABLES_DIR, "diffbind_gene_level_all_marks.tsv"),
  multi  = file.path(TABLES_DIR, "diffbind_gene_level_all_marks.tsv")
)

# Build combined profile
combined <- ratio_data %>%
  left_join(me3_gene %>% dplyr::select(gene, me3_fold, me3_fdr, me3_n_peaks), by = "gene") %>%
  left_join(me2_gene %>% dplyr::select(gene, me2_fold, me2_fdr, me2_n_peaks), by = "gene")

# Load existing multi-mark table if available for correlation expansion
multi_mark_table <- file.path(TABLES_DIR, "diffbind_gene_level_all_marks.tsv")
if (file.exists(multi_mark_table)) {
  multi_marks <- read.table(multi_mark_table, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  combined <- combined %>%
    left_join(multi_marks %>% dplyr::select(gene, any_of(c("atac_fold", "k27ac_fold",
                                                            "k27me3_fold", "k119ub_fold"))),
              by = "gene")
  cat(sprintf("  Multi-mark data joined (%d marks)\n",
              sum(c("atac_fold", "k27ac_fold", "k27me3_fold", "k119ub_fold") %in% colnames(combined))))
}

n_both <- sum(!is.na(combined$me2_fold) & !is.na(combined$me3_fold))
cat(sprintf("  Genes with both me2 and me3: %d\n", n_both))
cat(sprintf("  Total combined profile: %d genes\n", nrow(combined)))

# =============================================================================
# FIGURE 40a: Expanded Cross-Mark Correlation Heatmap
# =============================================================================

cat("\n--- FIGURE 40a: Expanded Correlation Heatmap ---\n\n")

# Select columns for correlation
cor_cols <- c("mc_diff", "hmc_diff", "delta_ratio", "me3_fold", "me2_fold")
cor_labels <- c("5mC diff", "5hmC diff", "Delta ratio", "H3K36me3", "H3K36me2")

# Add other marks if available
if ("k119ub_fold" %in% colnames(combined)) {
  cor_cols <- c(cor_cols, "k119ub_fold", "k27me3_fold", "k27ac_fold", "atac_fold")
  cor_labels <- c(cor_labels, "H2AK119ub", "H3K27me3", "H3K27ac", "ATAC-seq")
}

cor_data <- combined[, cor_cols, drop = FALSE]
colnames(cor_data) <- cor_labels

# Pairwise Spearman correlation
n_vars <- ncol(cor_data)
cor_mat <- matrix(NA, n_vars, n_vars, dimnames = list(cor_labels, cor_labels))
p_mat <- matrix(NA, n_vars, n_vars, dimnames = list(cor_labels, cor_labels))

for (i in seq_len(n_vars)) {
  for (j in seq_len(n_vars)) {
    if (i == j) {
      cor_mat[i, j] <- 1
      p_mat[i, j] <- 0
    } else {
      complete <- complete.cases(cor_data[, c(i, j)])
      if (sum(complete) >= 10) {
        ct <- cor.test(cor_data[complete, i], cor_data[complete, j], method = "spearman")
        cor_mat[i, j] <- ct$estimate
        p_mat[i, j] <- ct$p.value
      }
    }
  }
}

sig_stars <- function(p) {
  ifelse(is.na(p), "",
  ifelse(p < 0.001, "***",
  ifelse(p < 0.01, "**",
  ifelse(p < 0.05, "*", ""))))
}

display_mat <- matrix(
  sprintf("%.2f%s", cor_mat, sig_stars(p_mat)),
  nrow = n_vars, dimnames = list(cor_labels, cor_labels)
)
display_mat[is.na(cor_mat)] <- ""

hm_colors <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
hm_breaks <- seq(-1, 1, length.out = 101)

hm_call <- bquote(pheatmap(
  .(cor_mat),
  display_numbers = .(display_mat),
  color = .(hm_colors),
  breaks = .(hm_breaks),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_method = "complete",
  fontsize = 11,
  fontsize_number = 9,
  main = "Cross-Mark Spearman Correlations (Including H3K36me2/me3)",
  border_color = "grey80",
  na_col = "grey90",
  number_color = "black"
))

save_multiformat_pheatmap(
  hm_call,
  file.path(OUTPUT_DIR, SECTION_DIR, "40a_expanded_correlation_heatmap"),
  width = 10, height = 9
)

cat("  Correlation heatmap saved\n")

# =============================================================================
# FIGURE 40b: me2 vs me3 Fold-Change Scatter
# =============================================================================

cat("\n--- FIGURE 40b: me2 vs me3 Fold-Change Scatter ---\n\n")

dual_df <- combined %>%
  dplyr::filter(!is.na(me2_fold) & !is.na(me3_fold)) %>%
  mutate(
    mc_group = case_when(
      mc_sig & mc_diff > 0  ~ "mC Hyper",
      mc_sig & mc_diff <= 0 ~ "mC Hypo",
      TRUE ~ "Non-DMR"
    )
  )

cor_me2_me3 <- cor.test(dual_df$me2_fold, dual_df$me3_fold, method = "spearman")

# Quadrant counts
q1 <- sum(dual_df$me2_fold > 0 & dual_df$me3_fold > 0)
q2 <- sum(dual_df$me2_fold <= 0 & dual_df$me3_fold > 0)
q3 <- sum(dual_df$me2_fold <= 0 & dual_df$me3_fold <= 0)
q4 <- sum(dual_df$me2_fold > 0 & dual_df$me3_fold <= 0)

p_40b <- ggplot(dual_df, aes(x = me2_fold, y = me3_fold, color = mc_group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  geom_point(alpha = 0.35, size = 1) +
  geom_smooth(data = dual_df, aes(x = me2_fold, y = me3_fold),
              inherit.aes = FALSE, method = "lm",
              color = "black", fill = "grey80", alpha = 0.2, linewidth = 0.8) +
  scale_color_manual(values = c("mC Hyper" = "#D7191C", "mC Hypo" = "#2C7BB6", "Non-DMR" = "grey70"),
                     name = "mC Status") +
  # Quadrant annotations
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.3, size = 3,
           label = sprintf("Q1: Both gain\nn=%d", q1), color = "grey30") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.3, size = 3,
           label = sprintf("Q2: me3 gain/me2 loss\n(conversion shift)\nn=%d", q2), color = "grey30") +
  annotate("text", x = -Inf, y = -Inf, hjust = -0.1, vjust = -0.3, size = 3,
           label = sprintf("Q3: Both loss\nn=%d", q3), color = "grey30") +
  annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.3, size = 3,
           label = sprintf("Q4: me2 gain/me3 loss\n(conversion block)\nn=%d", q4), color = "grey30") +
  labs(
    title = "H3K36me2 vs H3K36me3 Fold Change (Gene Level)",
    subtitle = sprintf("Spearman rho = %+.3f, p = %.2e | n = %d genes",
                       cor_me2_me3$estimate, cor_me2_me3$p.value, nrow(dual_df)),
    x = "H3K36me2 Log2FC", y = "H3K36me3 Log2FC"
  ) +
  theme_biomodal()

save_multiformat_ggplot(p_40b, file.path(OUTPUT_DIR, SECTION_DIR, "40b_me2_vs_me3_scatter"),
                        width = 10, height = 9)

cat(sprintf("  me2 vs me3 rho = %+.3f, p = %.2e\n", cor_me2_me3$estimate, cor_me2_me3$p.value))
cat(sprintf("  Quadrants: Q1=%d, Q2=%d, Q3=%d, Q4=%d\n", q1, q2, q3, q4))

# =============================================================================
# FIGURE 40c: me2-me3 Ratio Delta by DMR Status
# =============================================================================

cat("\n--- FIGURE 40c: me2-me3 Ratio Delta by DMR Status ---\n\n")

ratio_delta_df <- dual_df %>%
  mutate(
    me2_me3_delta = me2_fold - me3_fold,
    dmr_group = case_when(
      mc_sig & hmc_sig & mc_diff > 0 & hmc_diff < 0 ~ "Coordinated\n(mC\u2191/hmC\u2193)",
      mc_sig & mc_diff > 0  ~ "mC Hyper\nOnly",
      mc_sig & mc_diff <= 0 ~ "mC Hypo\nOnly",
      TRUE ~ "Non-DMR"
    ),
    dmr_group = factor(dmr_group, levels = c("Coordinated\n(mC\u2191/hmC\u2193)",
                                              "mC Hyper\nOnly", "mC Hypo\nOnly", "Non-DMR"))
  )

kw_result <- kruskal.test(me2_me3_delta ~ dmr_group, data = ratio_delta_df)

group_medians <- ratio_delta_df %>%
  group_by(dmr_group) %>%
  summarise(median_delta = median(me2_me3_delta, na.rm = TRUE), n = n(), .groups = "drop")

p_40c <- ggplot(ratio_delta_df, aes(x = dmr_group, y = me2_me3_delta, fill = dmr_group)) +
  geom_violin(alpha = 0.5, scale = "width") +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("#D95F02", "#E6AB02", "#66A61E", "grey80")) +
  geom_text(data = group_medians,
            aes(x = dmr_group, y = -Inf, label = sprintf("n=%d\nmed=%+.2f", n, median_delta)),
            inherit.aes = FALSE, vjust = -0.3, size = 2.8, lineheight = 0.85, color = "grey30") +
  labs(
    title = "H3K36me2 - H3K36me3 Fold Difference by DMR Status",
    subtitle = sprintf("Kruskal-Wallis p = %.2e | Positive = me2 shift, Negative = me3 shift",
                       kw_result$p.value),
    x = "", y = "me2 Fold - me3 Fold (Delta)"
  ) +
  theme_biomodal() +
  theme(legend.position = "none")

save_multiformat_ggplot(p_40c, file.path(OUTPUT_DIR, SECTION_DIR, "40c_me2_me3_ratio_delta"),
                        width = 10, height = 7)

# =============================================================================
# FIGURE 40d: Venn — Genes with sig me2, sig me3, sig mC Changes
# =============================================================================

cat("\n--- FIGURE 40d: Three-Way Venn ---\n\n")

genes_me2_sig <- me2_gene$gene[!is.na(me2_gene$me2_fdr) & me2_gene$me2_fdr < FDR_THRESHOLD]
genes_me3_sig <- me3_gene$gene[!is.na(me3_gene$me3_fdr) & me3_gene$me3_fdr < FDR_THRESHOLD]
genes_mc_sig  <- ratio_data$gene[ratio_data$mc_sig]

venn_list <- list(
  "H3K36me2\nSig" = genes_me2_sig,
  "H3K36me3\nSig" = genes_me3_sig,
  "mC DMR" = genes_mc_sig
)

p_40d <- ggVennDiagram(venn_list, label = "both", label_alpha = 0,
                       edge_size = 0.5, set_size = 4) +
  scale_fill_gradient(low = "white", high = "#D95F02", name = "Count") +
  labs(title = "Gene Overlap: H3K36me2, H3K36me3, and mC DMR Significance") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))

save_multiformat_ggplot(p_40d, file.path(OUTPUT_DIR, SECTION_DIR, "40d_three_way_venn"),
                        width = 9, height = 8)

# Overlap stats
triple <- length(intersect(intersect(genes_me2_sig, genes_me3_sig), genes_mc_sig))
cat(sprintf("  Triple overlap (me2 + me3 + mC): %d genes\n", triple))
cat(sprintf("  me2 sig: %d | me3 sig: %d | mC DMR: %d\n",
            length(genes_me2_sig), length(genes_me3_sig), length(genes_mc_sig)))

# =============================================================================
# FIGURE 40e: Direction Flow for Triple-Significant Genes
# =============================================================================

cat("\n--- FIGURE 40e: Direction Flow (Grouped Bar) ---\n\n")

triple_genes <- intersect(intersect(genes_me2_sig, genes_me3_sig), genes_mc_sig)

if (length(triple_genes) >= 5) {
  flow_df <- combined %>%
    dplyr::filter(gene %in% triple_genes) %>%
    mutate(
      me2_dir = ifelse(me2_fold > 0, "me2\u2191", "me2\u2193"),
      me3_dir = ifelse(me3_fold > 0, "me3\u2191", "me3\u2193"),
      mc_dir  = ifelse(mc_diff > 0, "mC\u2191", "mC\u2193"),
      pattern = paste(me2_dir, me3_dir, mc_dir, sep = " + ")
    )

  pattern_counts <- flow_df %>%
    count(pattern, sort = TRUE) %>%
    mutate(pct = 100 * n / sum(n))

  # Show top patterns
  cat("  Direction patterns (triple-significant genes):\n")
  for (i in seq_len(min(8, nrow(pattern_counts)))) {
    cat(sprintf("    %s: %d genes (%.1f%%)\n",
                pattern_counts$pattern[i], pattern_counts$n[i], pattern_counts$pct[i]))
  }

  p_40e <- ggplot(pattern_counts %>% slice_head(n = 10),
                  aes(x = reorder(pattern, n), y = n, fill = n)) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%d (%.0f%%)", n, pct)), hjust = -0.1, size = 3.5) +
    scale_fill_gradient(low = "#FEE0D2", high = "#D95F02", guide = "none") +
    coord_flip() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
    labs(
      title = "Direction Patterns for Triple-Significant Genes",
      subtitle = sprintf("Genes significant in me2, me3, AND mC (n = %d)", length(triple_genes)),
      x = "Direction Combination", y = "Number of Genes"
    ) +
    theme_biomodal()

  save_multiformat_ggplot(p_40e, file.path(OUTPUT_DIR, SECTION_DIR, "40e_direction_flow"),
                          width = 10, height = 7)
} else {
  cat(sprintf("  Only %d triple-significant genes — skipping direction flow plot\n",
              length(triple_genes)))
}

# =============================================================================
# FIGURE 40f: GO Enrichment — me2-only vs me3-only vs Shared
# =============================================================================

cat("\n--- FIGURE 40f: GO Enrichment Comparison ---\n\n")

me2_only_genes <- setdiff(genes_me2_sig, genes_me3_sig)
me3_only_genes <- setdiff(genes_me3_sig, genes_me2_sig)
shared_genes <- intersect(genes_me2_sig, genes_me3_sig)

cat(sprintf("  me2-only: %d | me3-only: %d | shared: %d\n",
            length(me2_only_genes), length(me3_only_genes), length(shared_genes)))

# Convert to Entrez IDs for enrichment
convert_to_entrez <- function(symbols) {
  mapped <- suppressMessages(bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID",
                                   OrgDb = org.Mm.eg.db))
  mapped$ENTREZID
}

run_go <- function(genes, label) {
  entrez <- convert_to_entrez(genes)
  if (length(entrez) < 5) return(NULL)
  res <- suppressMessages(enrichGO(
    gene = entrez,
    OrgDb = org.Mm.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05,
    readable = TRUE
  ))
  if (is.null(res) || nrow(as.data.frame(res)) == 0) return(NULL)
  res_df <- as.data.frame(res) %>%
    slice_head(n = 10) %>%
    mutate(group = label)
  res_df
}

go_me2 <- run_go(me2_only_genes, "H3K36me2 Only")
go_me3 <- run_go(me3_only_genes, "H3K36me3 Only")
go_shared <- run_go(shared_genes, "Shared")

go_combined <- bind_rows(go_me2, go_me3, go_shared)

if (nrow(go_combined) > 0) {
  go_combined$Description <- factor(go_combined$Description,
                                     levels = rev(unique(go_combined$Description)))

  p_40f <- ggplot(go_combined, aes(x = -log10(p.adjust), y = Description,
                                    size = Count, color = group)) +
    geom_point(alpha = 0.8) +
    scale_color_manual(values = c("H3K36me2 Only" = "#E6AB02", "H3K36me3 Only" = "#D95F02",
                                   "Shared" = "#984EA3"),
                       name = "Gene Set") +
    scale_size_continuous(range = c(2, 6), name = "Gene Count") +
    facet_wrap(~ group, scales = "free_y", ncol = 1) +
    labs(
      title = "GO Biological Process: me2-Only vs me3-Only vs Shared Differential Genes",
      subtitle = "Top 10 terms per group (BH-adjusted q < 0.05)",
      x = "-log10(Adjusted p-value)", y = ""
    ) +
    theme_biomodal() +
    theme(strip.text = element_text(size = 11, face = "bold"))

  save_multiformat_ggplot(p_40f, file.path(OUTPUT_DIR, SECTION_DIR, "40f_go_comparison"),
                          width = 12, height = 14)
} else {
  cat("  No significant GO terms — skipping plot\n")
}

# =============================================================================
# EXPORT TABLES
# =============================================================================

cat("\n--- Exporting Tables ---\n\n")

export_combined <- combined %>%
  dplyr::select(gene, mc_diff, mc_sig, hmc_diff, hmc_sig, delta_ratio, dmr_status,
                chromatin_state, me3_fold, me3_fdr, me3_n_peaks,
                me2_fold, me2_fdr, me2_n_peaks,
                any_of(c("k119ub_fold", "k27me3_fold", "k27ac_fold", "atac_fold",
                         "k119ub_fdr", "k27me3_fdr", "k27ac_fdr", "atac_fdr")))

write.table(export_combined, file.path(TABLES_DIR, "h3k36_combined_gene_profile.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved h3k36_combined_gene_profile.tsv (%d rows)\n", nrow(export_combined)))

# Correlation matrix export
cor_export <- as.data.frame(cor_mat) %>%
  tibble::rownames_to_column("variable")
write.table(cor_export, file.path(TABLES_DIR, "h3k36_expanded_correlations.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved h3k36_expanded_correlations.tsv\n")

cat("\n================================================================================\n")
cat("SECTION 40 COMPLETE\n")
cat("================================================================================\n")
