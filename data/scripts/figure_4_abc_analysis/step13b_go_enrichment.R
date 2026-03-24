# abc/scripts/step13b_go_enrichment.R
#
# GO/KEGG pathway enrichment comparison: concordant vs discordant gene sets.
# Extends step13 by asking whether concordant and discordant genes are enriched
# in different biological pathways, using clusterProfiler's compareCluster.

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(dplyr)
})

cat("=== Step 13b: GO/KEGG Enrichment (Concordant vs Discordant) ===\n")
cat(sprintf("  Started: %s\n\n", Sys.time()))

# =============================================================================
# CONFIGURATION
# =============================================================================

DISCORDANT_FILE <- "data/tsvs/figure_4_abc_analysis/4B_discordant_gene_characteristics.tsv"  # Original: results/discordant_gene_characteristics.tsv
OUTPUT_DIR      <- "data/plots/figure_4_abc_analysis"  # Original: results/figures/discordant_analysis
TSV_DIR         <- "data/tsvs/figure_4_abc_analysis"  # Original: results/figures/discordant_analysis (TSV outputs)
UTIL_PATH       <- "data/scripts/_shared/multi_format_output.R"  # Original: ../scripts/utils/multi_format_output.R

stopifnot(file.exists(DISCORDANT_FILE))
stopifnot(file.exists(UTIL_PATH))
source(UTIL_PATH)

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TSV_DIR, recursive = TRUE, showWarnings = FALSE)

# Shared aesthetics
theme_pub <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", color = NA),
    plot.title = element_text(size = 12, face = "bold"),
    legend.position = "bottom"
  )

CONCORDANCE_COLORS <- c(Concordant = "#4DAF4A", Discordant = "#E41A1C")

save_plot <- function(p, name, w = 7, h = 6) {
  save_multiformat_ggplot(p, base_path = file.path(OUTPUT_DIR, name),
                          width = w, height = h, dpi = 300,
                          verbose = TRUE, use_subfolders = TRUE)
}

# =============================================================================
# LOAD DATA
# =============================================================================

cat("Loading data...\n")

genes <- read.table(DISCORDANT_FILE, sep = "\t", header = TRUE,
                    stringsAsFactors = FALSE)
cat(sprintf("  Total dysregulated genes: %d\n", nrow(genes)))
cat(sprintf("  Concordant: %d | Discordant: %d\n",
            sum(genes$concordance == "Concordant"),
            sum(genes$concordance == "Discordant")))

# =============================================================================
# GENE ID CONVERSION (SYMBOL → ENTREZID)
# =============================================================================

cat("\nConverting gene symbols to Entrez IDs...\n")

conc_genes <- genes$TargetGene[genes$concordance == "Concordant"]
disc_genes <- genes$TargetGene[genes$concordance == "Discordant"]

conc_ids <- bitr(conc_genes, fromType = "SYMBOL", toType = "ENTREZID",
                 OrgDb = org.Mm.eg.db)
disc_ids <- bitr(disc_genes, fromType = "SYMBOL", toType = "ENTREZID",
                 OrgDb = org.Mm.eg.db)

conc_rate <- nrow(conc_ids) / length(conc_genes)
disc_rate <- nrow(disc_ids) / length(disc_genes)

cat(sprintf("  Concordant: %d / %d mapped (%.1f%%)\n",
            nrow(conc_ids), length(conc_genes), 100 * conc_rate))
cat(sprintf("  Discordant: %d / %d mapped (%.1f%%)\n",
            nrow(disc_ids), length(disc_genes), 100 * disc_rate))

stopifnot(conc_rate >= 0.50)
stopifnot(disc_rate >= 0.50)

# =============================================================================
# compareCluster — GO BP
# =============================================================================

cat("\nRunning compareCluster (GO BP)...\n")

gene_list <- list(
  Concordant = unique(conc_ids$ENTREZID),
  Discordant = unique(disc_ids$ENTREZID)
)

cc_bp <- compareCluster(
  geneClusters = gene_list,
  fun = "enrichGO",
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.10,
  readable = TRUE
)

n_bp_terms <- if (!is.null(cc_bp) && nrow(as.data.frame(cc_bp)) > 0) {
  nrow(as.data.frame(cc_bp))
} else {
  0
}
cat(sprintf("  GO BP enriched terms: %d\n", n_bp_terms))

if (n_bp_terms > 0) {
  # Save results table
  bp_df <- as.data.frame(cc_bp)
  write.table(bp_df, file.path(TSV_DIR, "4B_discordant_go_bp.tsv"),  # Original: OUTPUT_DIR/go_bp_enrichment_results.tsv
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  Concordant-specific terms: %d\n",
              sum(bp_df$Cluster == "Concordant")))
  cat(sprintf("  Discordant-specific terms: %d\n",
              sum(bp_df$Cluster == "Discordant")))

  # Dotplot
  p_bp <- dotplot(cc_bp, showCategory = 10,
                  title = "GO Biological Process: Concordant vs Discordant") +
    scale_color_continuous(low = "#E41A1C", high = "#377EB8", name = "p.adjust") +
    theme_pub +
    theme(axis.text.y = element_text(size = 9))

  save_plot(p_bp, "4B_go_bp_compareCluster", w = 10, h = 8)  # Original: 09_go_bp_compareCluster
} else {
  cat("  No significant GO BP terms — skipping plot.\n")
}

# =============================================================================
# compareCluster — KEGG
# =============================================================================

cat("\nRunning compareCluster (KEGG)...\n")

cc_kegg <- compareCluster(
  geneClusters = gene_list,
  fun = "enrichKEGG",
  organism = "mmu",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.10
)

n_kegg_terms <- if (!is.null(cc_kegg) && nrow(as.data.frame(cc_kegg)) > 0) {
  nrow(as.data.frame(cc_kegg))
} else {
  0
}
cat(sprintf("  KEGG enriched terms: %d\n", n_kegg_terms))

if (n_kegg_terms > 0) {
  # Convert Entrez → readable gene names in KEGG results
  cc_kegg_read <- setReadable(cc_kegg, OrgDb = org.Mm.eg.db,
                              keyType = "ENTREZID")
  kegg_df <- as.data.frame(cc_kegg_read)
  write.table(kegg_df, file.path(TSV_DIR, "4B_discordant_kegg.tsv"),  # Original: OUTPUT_DIR/kegg_enrichment_results.tsv
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  Concordant-specific pathways: %d\n",
              sum(kegg_df$Cluster == "Concordant")))
  cat(sprintf("  Discordant-specific pathways: %d\n",
              sum(kegg_df$Cluster == "Discordant")))

  p_kegg <- dotplot(cc_kegg_read, showCategory = 10,
                    title = "KEGG Pathways: Concordant vs Discordant") +
    scale_color_continuous(low = "#E41A1C", high = "#377EB8", name = "p.adjust") +
    theme_pub +
    theme(axis.text.y = element_text(size = 9))

  save_plot(p_kegg, "4B_kegg_compareCluster", w = 10, h = 8)  # Original: 10_kegg_compareCluster
} else {
  cat("  No significant KEGG pathways — skipping plot.\n")
}

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n=== Summary ===\n")
cat(sprintf("  Gene lists: Concordant = %d, Discordant = %d (Entrez IDs)\n",
            length(gene_list$Concordant), length(gene_list$Discordant)))
cat(sprintf("  GO BP terms: %d | KEGG pathways: %d\n", n_bp_terms, n_kegg_terms))
if (n_bp_terms > 0) {
  bp_df_summary <- as.data.frame(cc_bp)
  top_conc <- head(bp_df_summary[bp_df_summary$Cluster == "Concordant", "Description"], 3)
  top_disc <- head(bp_df_summary[bp_df_summary$Cluster == "Discordant", "Description"], 3)
  if (length(top_conc) > 0) {
    cat("  Top concordant GO BP:\n")
    for (t in top_conc) cat(sprintf("    - %s\n", t))
  }
  if (length(top_disc) > 0) {
    cat("  Top discordant GO BP:\n")
    for (t in top_disc) cat(sprintf("    - %s\n", t))
  }
}

cat(sprintf("\nCompleted: %s\n", Sys.time()))
cat("=== Step 13b complete ===\n")
