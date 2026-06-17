# biomodal/downstream/scripts/viz_sections/section_61h_mecp2_up_quadrant_enrichment.R
# Section 61h: GO enrichment of MeCP2-Up + K119ub-Up quadrant genes
#
# Extracts genes from the top-right quadrant of the 59a scatter
# (K119ub log2fc > 0 AND MeCP2 significantly up) and runs GO BP enrichment
# to test whether this consistent subset is neuronally enriched.
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_61h_mecp2_up_quadrant_enrichment.R

source("scripts/viz_sections/_shared_config.R")

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
})

cat("================================================================================\n")
cat("SECTION 61h: MeCP2-Up + K119ub-Up QUADRANT GO ENRICHMENT\n")
cat("================================================================================\n\n")

SEC61H_DIR <- file.path(OUTPUT_DIR, "61h_mecp2_up_quadrant_enrichment")
dir.create(SEC61H_DIR, recursive = TRUE, showWarnings = FALSE)

QUADRANT_PATH <- file.path(TABLES_DIR, "59_quadrant_master.tsv")
if (!file.exists(QUADRANT_PATH)) stop("Run section 59 first: ", QUADRANT_PATH)

qm <- read.table(QUADRANT_PATH, header = TRUE, sep = "\t",
                  stringsAsFactors = FALSE, quote = "")
cat(sprintf("  Loaded %d genes from quadrant master\n", nrow(qm)))

# --- Filter to top-right quadrant ---

qm$mecp2_up <- !is.na(qm$mecp2_nearest_fdr) &
  qm$mecp2_nearest_fdr < Q_THRESHOLD & qm$mecp2_mean_fold > 0

qm$k119ub_up <- !is.na(qm$gb_log2fc) & qm$gb_log2fc > 0 &
  qm$gb_signal_class == "quantifiable"

quadrant_genes <- qm$gene[qm$mecp2_up & qm$k119ub_up]
mecp2_up_all   <- qm$gene[qm$mecp2_up]
background     <- qm$gene[!is.na(qm$gb_log2fc) & qm$gb_signal_class == "quantifiable"]

cat(sprintf("\n  MeCP2-Up total:               %d genes\n", length(mecp2_up_all)))
cat(sprintf("  K119ub quantifiable:          %d genes\n", length(background)))
cat(sprintf("  Top-right quadrant (Up+Up):   %d genes\n", length(quadrant_genes)))
cat(sprintf("  Top-left quadrant (Up+Down):  %d genes\n",
            sum(qm$mecp2_up & !is.na(qm$gb_log2fc) & qm$gb_log2fc <= 0 &
                qm$gb_signal_class == "quantifiable")))

# Print the gene list
cat(sprintf("\n  Quadrant genes: %s\n", paste(sort(quadrant_genes), collapse = ", ")))

# Save gene list
gene_df <- qm[qm$gene %in% quadrant_genes,
               c("gene", "mecp2_mean_fold", "mecp2_nearest_fdr",
                 "gb_log2fc", "mc_diff", "hmc_diff", "chromatin_state")]
gene_df <- gene_df[order(-gene_df$mecp2_mean_fold), ]

write.table(gene_df, file.path(TABLES_DIR, "61h_mecp2_up_k119ub_up_genes.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("  Saved gene table: 61h_mecp2_up_k119ub_up_genes.tsv (%d genes)\n", nrow(gene_df)))

# --- Convert to Entrez IDs ---

cat("\n--- Running GO BP enrichment ---\n")

entrez_map <- bitr(quadrant_genes, fromType = "SYMBOL", toType = "ENTREZID",
                   OrgDb = org.Mm.eg.db)
cat(sprintf("  Mapped %d / %d symbols to Entrez IDs\n",
            nrow(entrez_map), length(quadrant_genes)))

bg_entrez <- bitr(background, fromType = "SYMBOL", toType = "ENTREZID",
                  OrgDb = org.Mm.eg.db)
cat(sprintf("  Background: %d Entrez IDs\n", nrow(bg_entrez)))

# --- enrichGO ---

ego <- enrichGO(
  gene         = entrez_map$ENTREZID,
  universe     = bg_entrez$ENTREZID,
  OrgDb        = org.Mm.eg.db,
  ont          = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

n_sig <- sum(ego@result$p.adjust < 0.05)
cat(sprintf("  Significant GO BP terms (q < 0.05): %d\n", n_sig))

if (n_sig > 0) {
  top_terms <- head(ego@result[ego@result$p.adjust < 0.05, ], 20)
  cat("\n  Top GO BP terms:\n")
  for (i in seq_len(nrow(top_terms))) {
    r <- top_terms[i, ]
    cat(sprintf("    %-55s q=%.2e  %s\n",
                substr(r$Description, 1, 55), r$p.adjust, r$GeneRatio))
  }

  # Flag neuronal terms
  neuronal_mask <- grepl("synap|neuron|axon|dendrit|nervous|brain|cerebel",
                         top_terms$Description, ignore.case = TRUE)
  n_neuronal <- sum(neuronal_mask)
  cat(sprintf("\n  Of top 20 terms: %d are neuronal/synaptic (%.0f%%)\n",
              n_neuronal, 100 * n_neuronal / nrow(top_terms)))
}

# Save full enrichment table
write.table(ego@result, file.path(TABLES_DIR, "61h_quadrant_go_bp.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("  Saved: 61h_quadrant_go_bp.tsv (%d terms)\n", nrow(ego@result)))

# --- Dot plot ---

if (n_sig >= 3) {
  p_dot <- dotplot(ego, showCategory = min(20, n_sig)) +
    labs(title = "GO BP enrichment: MeCP2-Up + K119ub-Up quadrant genes",
         subtitle = sprintf("n = %d genes | Top-right quadrant from section 59a",
                            length(quadrant_genes))) +
    theme_biomodal()

  save_multiformat_ggplot(p_dot,
    base_path = file.path(SEC61H_DIR, "61h_go_dotplot"),
    width = 11, height = 9, dpi = 300, verbose = TRUE, use_subfolders = TRUE)
}

# --- Also try KEGG ---

cat("\n--- Running KEGG enrichment ---\n")

ekegg <- enrichKEGG(
  gene     = entrez_map$ENTREZID,
  universe = bg_entrez$ENTREZID,
  organism = "mmu",
  pvalueCutoff = 0.05
)

n_kegg <- sum(ekegg@result$p.adjust < 0.05)
cat(sprintf("  Significant KEGG pathways (q < 0.05): %d\n", n_kegg))

if (n_kegg > 0) {
  top_kegg <- head(ekegg@result[ekegg@result$p.adjust < 0.05, ], 10)
  cat("\n  Top KEGG pathways:\n")
  for (i in seq_len(nrow(top_kegg))) {
    r <- top_kegg[i, ]
    cat(sprintf("    %-50s q=%.2e  %s\n",
                substr(r$Description, 1, 50), r$p.adjust, r$GeneRatio))
  }
}

write.table(ekegg@result, file.path(TABLES_DIR, "61h_quadrant_kegg.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("  Saved: 61h_quadrant_kegg.tsv\n"))

# =============================================================================
# RELAXED ANALYSIS: quadrant position only (no MeCP2 FDR cutoff)
# =============================================================================

cat("\n================================================================================\n")
cat("RELAXED ANALYSIS: MeCP2 fold > 0 + K119ub log2fc > 0 (no FDR filter)\n")
cat("================================================================================\n\n")

qm$mecp2_fold_up <- !is.na(qm$mecp2_mean_fold) & qm$mecp2_mean_fold > 0

relaxed_genes <- qm$gene[qm$mecp2_fold_up & qm$k119ub_up]
cat(sprintf("  Relaxed top-right quadrant: %d genes (vs %d strict)\n",
            length(relaxed_genes), length(quadrant_genes)))

# Check key genes
for (g in c("Syt1", "Zbtb20", "Trpm3", "Ampd3")) {
  in_set <- g %in% relaxed_genes
  row <- qm[qm$gene == g, ]
  if (nrow(row) > 0) {
    cat(sprintf("  %s: fold=%.3f  FDR=%.3f  K119ub=%.3f  in_relaxed=%s\n",
                g, row$mecp2_mean_fold, row$mecp2_nearest_fdr, row$gb_log2fc,
                ifelse(in_set, "YES", "no")))
  }
}

# Save relaxed gene list
relaxed_df <- qm[qm$gene %in% relaxed_genes,
                  c("gene", "mecp2_mean_fold", "mecp2_nearest_fdr",
                    "gb_log2fc", "mc_diff", "hmc_diff", "chromatin_state")]
relaxed_df <- relaxed_df[order(-relaxed_df$mecp2_mean_fold), ]

write.table(relaxed_df, file.path(TABLES_DIR, "61h_relaxed_quadrant_genes.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("  Saved: 61h_relaxed_quadrant_genes.tsv (%d genes)\n", nrow(relaxed_df)))

# --- Relaxed GO enrichment ---

cat("\n--- Relaxed GO BP enrichment ---\n")

relax_entrez <- bitr(relaxed_genes, fromType = "SYMBOL", toType = "ENTREZID",
                     OrgDb = org.Mm.eg.db)
cat(sprintf("  Mapped %d / %d symbols to Entrez IDs\n",
            nrow(relax_entrez), length(relaxed_genes)))

ego_relax <- enrichGO(
  gene         = relax_entrez$ENTREZID,
  universe     = bg_entrez$ENTREZID,
  OrgDb        = org.Mm.eg.db,
  ont          = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

n_sig_relax <- sum(ego_relax@result$p.adjust < 0.05)
cat(sprintf("  Significant GO BP terms (q < 0.05): %d\n", n_sig_relax))

if (n_sig_relax > 0) {
  top_relax <- head(ego_relax@result[ego_relax@result$p.adjust < 0.05, ], 25)
  cat("\n  Top GO BP terms (relaxed):\n")
  for (i in seq_len(nrow(top_relax))) {
    r <- top_relax[i, ]
    cat(sprintf("    %-55s q=%.2e  %s\n",
                substr(r$Description, 1, 55), r$p.adjust, r$GeneRatio))
  }

  neuronal_mask_r <- grepl("synap|neuron|axon|dendrit|nervous|brain|cerebel",
                           top_relax$Description, ignore.case = TRUE)
  n_neuronal_r <- sum(neuronal_mask_r)
  cat(sprintf("\n  Of top 25 terms: %d are neuronal/synaptic (%.0f%%)\n",
              n_neuronal_r, 100 * n_neuronal_r / nrow(top_relax)))
}

write.table(ego_relax@result, file.path(TABLES_DIR, "61h_relaxed_go_bp.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("  Saved: 61h_relaxed_go_bp.tsv (%d terms)\n", nrow(ego_relax@result)))

# Dot plot for relaxed
if (n_sig_relax >= 3) {
  p_dot_relax <- dotplot(ego_relax, showCategory = min(25, n_sig_relax)) +
    labs(title = "GO BP: Relaxed top-right quadrant (MeCP2 fold>0 + K119ub>0)",
         subtitle = sprintf("n = %d genes | No FDR filter on MeCP2",
                            length(relaxed_genes))) +
    theme_biomodal()

  save_multiformat_ggplot(p_dot_relax,
    base_path = file.path(SEC61H_DIR, "61h_relaxed_go_dotplot"),
    width = 12, height = 10, dpi = 300, verbose = TRUE, use_subfolders = TRUE)
}

# --- Relaxed KEGG ---

cat("\n--- Relaxed KEGG enrichment ---\n")

ekegg_relax <- enrichKEGG(
  gene     = relax_entrez$ENTREZID,
  universe = bg_entrez$ENTREZID,
  organism = "mmu",
  pvalueCutoff = 0.05
)

n_kegg_relax <- sum(ekegg_relax@result$p.adjust < 0.05)
cat(sprintf("  Significant KEGG pathways (q < 0.05): %d\n", n_kegg_relax))

if (n_kegg_relax > 0) {
  top_kegg_r <- head(ekegg_relax@result[ekegg_relax@result$p.adjust < 0.05, ], 10)
  cat("\n  Top KEGG pathways (relaxed):\n")
  for (i in seq_len(nrow(top_kegg_r))) {
    r <- top_kegg_r[i, ]
    cat(sprintf("    %-50s q=%.2e  %s\n",
                substr(r$Description, 1, 50), r$p.adjust, r$GeneRatio))
  }
}

write.table(ekegg_relax@result, file.path(TABLES_DIR, "61h_relaxed_kegg.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("  Saved: 61h_relaxed_kegg.tsv\n"))

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n================================================================================\n")
cat("SECTION 61h COMPLETE\n")
cat("================================================================================\n")
cat(sprintf("\n  STRICT (MeCP2 FDR<0.05 + K119ub>0):  %d genes, %d GO terms\n",
            length(quadrant_genes), n_sig))
cat(sprintf("  RELAXED (MeCP2 fold>0 + K119ub>0):   %d genes, %d GO terms\n",
            length(relaxed_genes), n_sig_relax))
cat("================================================================================\n")
