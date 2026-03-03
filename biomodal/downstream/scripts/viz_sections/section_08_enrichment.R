# biomodal/downstream/scripts/viz_sections/section_08_enrichment.R
# Section 8: GO/KEGG Enrichment Analysis
# Standalone script - sources shared config for all dependencies and data

# Run from downstream/ directory
source("scripts/viz_sections/_shared_config.R")

# =============================================================================
# SECTION 8: GO/KEGG ENRICHMENT ANALYSIS
# =============================================================================

cat("================================================================================\n")
cat("SECTION 8: GO/KEGG ENRICHMENT ANALYSIS\n")
cat("================================================================================\n\n")

# Get hypermethylated genes (mC up in mutant)
hyper_genes <- mc_dmr %>%
  filter(significant & mod_difference > 0) %>%
  pull(gene) %>%
  unique()

cat(sprintf("Running enrichment on %d hypermethylated genes...\n", length(hyper_genes)))

# Convert gene symbols to Entrez IDs
tryCatch({
  gene_ids <- bitr(hyper_genes, fromType = "SYMBOL", toType = "ENTREZID",
                   OrgDb = org.Mm.eg.db)

  if (nrow(gene_ids) > 10) {
    cat(sprintf("  Converted %d genes to Entrez IDs\n", nrow(gene_ids)))

    # GO Biological Process
    cat("Running GO Biological Process enrichment...\n")
    go_bp <- enrichGO(gene = gene_ids$ENTREZID,
                      OrgDb = org.Mm.eg.db,
                      ont = "BP",
                      pAdjustMethod = "BH",
                      qvalueCutoff = 0.05,
                      readable = TRUE)

    if (!is.null(go_bp) && nrow(go_bp@result) > 0) {
      p_go_bp <- dotplot(go_bp, showCategory = 15) +
        labs(title = "GO Biological Process Enrichment",
             subtitle = "Hypermethylated genes in BAP1-KO") +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))

      save_multiformat_ggplot(p_go_bp, file.path(OUTPUT_DIR, "08a_enrichment_go_bp"),
                              width = 10, height = 10)

      write.table(go_bp@result, file.path(TABLES_DIR, "enrichment_go_bp.tsv"),
                  sep = "\t", quote = FALSE, row.names = FALSE)
      cat("  Saved GO BP results\n")
    }

    # GO Cellular Component
    cat("Running GO Cellular Component enrichment...\n")
    go_cc <- enrichGO(gene = gene_ids$ENTREZID,
                      OrgDb = org.Mm.eg.db,
                      ont = "CC",
                      pAdjustMethod = "BH",
                      qvalueCutoff = 0.05,
                      readable = TRUE)

    if (!is.null(go_cc) && nrow(go_cc@result) > 0) {
      p_go_cc <- dotplot(go_cc, showCategory = 15) +
        labs(title = "GO Cellular Component Enrichment",
             subtitle = "Hypermethylated genes in BAP1-KO") +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))

      save_multiformat_ggplot(p_go_cc, file.path(OUTPUT_DIR, "08b_enrichment_go_cc"),
                              width = 10, height = 10)

      write.table(go_cc@result, file.path(TABLES_DIR, "enrichment_go_cc.tsv"),
                  sep = "\t", quote = FALSE, row.names = FALSE)
      cat("  Saved GO CC results\n")
    }

    # GO Molecular Function
    cat("Running GO Molecular Function enrichment...\n")
    go_mf <- enrichGO(gene = gene_ids$ENTREZID,
                      OrgDb = org.Mm.eg.db,
                      ont = "MF",
                      pAdjustMethod = "BH",
                      qvalueCutoff = 0.05,
                      readable = TRUE)

    if (!is.null(go_mf) && nrow(go_mf@result) > 0) {
      p_go_mf <- dotplot(go_mf, showCategory = 15) +
        labs(title = "GO Molecular Function Enrichment",
             subtitle = "Hypermethylated genes in BAP1-KO") +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))

      save_multiformat_ggplot(p_go_mf, file.path(OUTPUT_DIR, "08c_enrichment_go_mf"),
                              width = 10, height = 10)

      write.table(go_mf@result, file.path(TABLES_DIR, "enrichment_go_mf.tsv"),
                  sep = "\t", quote = FALSE, row.names = FALSE)
      cat("  Saved GO MF results\n")
    }

    # KEGG Pathway
    cat("Running KEGG pathway enrichment...\n")
    kegg <- enrichKEGG(gene = gene_ids$ENTREZID,
                       organism = "mmu",
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.1)

    if (!is.null(kegg) && nrow(kegg@result) > 0) {
      p_kegg <- dotplot(kegg, showCategory = 15) +
        labs(title = "KEGG Pathway Enrichment",
             subtitle = "Hypermethylated genes in BAP1-KO") +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))

      save_multiformat_ggplot(p_kegg, file.path(OUTPUT_DIR, "08d_enrichment_kegg"),
                              width = 10, height = 10)

      write.table(kegg@result, file.path(TABLES_DIR, "enrichment_kegg.tsv"),
                  sep = "\t", quote = FALSE, row.names = FALSE)
      cat("  Saved KEGG results\n")
    }
  } else {
    cat("  Not enough genes converted for enrichment analysis\n")
  }
}, error = function(e) {
  cat(sprintf("  Enrichment analysis error: %s\n", e$message))
})

# =============================================================================
# FIGURES 08e-08h: DELTA-RATIO DECILE ENRICHMENT
# =============================================================================
# Genes ranked by delta_ratio (TET efficiency change) from Section 22.
# D10 = most negative (most TET-impaired), D1 = most positive (least impaired).
# Compare top vs bottom decile pathways with compareCluster + standalone GO.

cat("\n")
cat("================================================================================\n")
cat("DELTA-RATIO DECILE ENRICHMENT (Figures 08e-08h)\n")
cat("================================================================================\n\n")

delta_ratio_path <- file.path(TABLES_DIR, "demethylation_ratio_all_genes.tsv")
stopifnot("demethylation_ratio_all_genes.tsv not found — run Section 22 first" =
            file.exists(delta_ratio_path))

delta_ratio_df <- read.table(delta_ratio_path, header = TRUE, sep = "\t",
                             stringsAsFactors = FALSE)
cat(sprintf("Loaded %d genes with delta_ratio values\n", nrow(delta_ratio_df)))

# Compute decile boundaries
decile_breaks <- quantile(delta_ratio_df$delta_ratio, probs = seq(0, 1, 0.1),
                          na.rm = TRUE)
delta_ratio_df$decile <- cut(delta_ratio_df$delta_ratio,
                             breaks = decile_breaks,
                             labels = paste0("D", 1:10),
                             include.lowest = TRUE)

# D10 = most negative delta_ratio (most TET-impaired)
# D1 = most positive delta_ratio (least impaired / TET-enhanced)
top_decile_genes <- delta_ratio_df$gene[delta_ratio_df$decile == "D10"]
bottom_decile_genes <- delta_ratio_df$gene[delta_ratio_df$decile == "D1"]

cat(sprintf("Top decile (D10, most TET-impaired): %d genes\n", length(top_decile_genes)))
cat(sprintf("  delta_ratio range: [%.3f, %.3f]\n",
            min(delta_ratio_df$delta_ratio[delta_ratio_df$decile == "D10"]),
            max(delta_ratio_df$delta_ratio[delta_ratio_df$decile == "D10"])))
cat(sprintf("Bottom decile (D1, least impaired): %d genes\n", length(bottom_decile_genes)))
cat(sprintf("  delta_ratio range: [%.3f, %.3f]\n",
            min(delta_ratio_df$delta_ratio[delta_ratio_df$decile == "D1"]),
            max(delta_ratio_df$delta_ratio[delta_ratio_df$decile == "D1"])))

# Convert gene symbols to Entrez IDs for each decile
tryCatch({
  top_ids <- bitr(top_decile_genes, fromType = "SYMBOL", toType = "ENTREZID",
                  OrgDb = org.Mm.eg.db)
  bottom_ids <- bitr(bottom_decile_genes, fromType = "SYMBOL", toType = "ENTREZID",
                     OrgDb = org.Mm.eg.db)

  cat(sprintf("  D10 converted: %d / %d genes\n", nrow(top_ids), length(top_decile_genes)))
  cat(sprintf("  D1 converted: %d / %d genes\n", nrow(bottom_ids), length(bottom_decile_genes)))

  # -------------------------------------------------------------------------
  # FIGURE 08e: compareCluster GO BP — Top vs Bottom Decile
  # -------------------------------------------------------------------------
  if (nrow(top_ids) > 10 && nrow(bottom_ids) > 10) {
    cat("\nCreating Figure 08e: compareCluster GO BP...\n")

    gene_clusters <- list(
      "D10 (Most TET-impaired)" = top_ids$ENTREZID,
      "D1 (Least impaired)" = bottom_ids$ENTREZID
    )

    cc_go_bp <- compareCluster(gene_clusters, fun = "enrichGO",
                               OrgDb = org.Mm.eg.db,
                               ont = "BP",
                               pAdjustMethod = "BH",
                               qvalueCutoff = 0.05,
                               readable = TRUE)

    if (!is.null(cc_go_bp) && nrow(cc_go_bp@compareClusterResult) > 0) {
      p_08e <- dotplot(cc_go_bp, showCategory = 10) +
        labs(title = "GO Biological Process: Delta-Ratio Decile Comparison",
             subtitle = "D10 (most TET-impaired) vs D1 (least impaired)") +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"))

      save_multiformat_ggplot(p_08e, file.path(OUTPUT_DIR, "08e_enrichment_delta_ratio_compare_go_bp"),
                              width = 12, height = 10)

      write.table(cc_go_bp@compareClusterResult,
                  file.path(TABLES_DIR, "enrichment_delta_ratio_compare_go_bp.tsv"),
                  sep = "\t", quote = FALSE, row.names = FALSE)
      cat("  Saved GO BP compareCluster results\n")
    } else {
      cat("  No significant GO BP terms in compareCluster\n")
    }

    # -----------------------------------------------------------------------
    # FIGURE 08f: compareCluster KEGG — Top vs Bottom Decile
    # -----------------------------------------------------------------------
    cat("Creating Figure 08f: compareCluster KEGG...\n")

    cc_kegg <- compareCluster(gene_clusters, fun = "enrichKEGG",
                              organism = "mmu",
                              pAdjustMethod = "BH",
                              qvalueCutoff = 0.1)

    if (!is.null(cc_kegg) && nrow(cc_kegg@compareClusterResult) > 0) {
      p_08f <- dotplot(cc_kegg, showCategory = 10) +
        labs(title = "KEGG Pathway: Delta-Ratio Decile Comparison",
             subtitle = "D10 (most TET-impaired) vs D1 (least impaired)") +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"))

      save_multiformat_ggplot(p_08f, file.path(OUTPUT_DIR, "08f_enrichment_delta_ratio_compare_kegg"),
                              width = 12, height = 10)

      write.table(cc_kegg@compareClusterResult,
                  file.path(TABLES_DIR, "enrichment_delta_ratio_compare_kegg.tsv"),
                  sep = "\t", quote = FALSE, row.names = FALSE)
      cat("  Saved KEGG compareCluster results\n")
    } else {
      cat("  No significant KEGG pathways in compareCluster\n")
    }

    # -----------------------------------------------------------------------
    # FIGURE 08g: Standalone GO BP — Top Decile (D10)
    # -----------------------------------------------------------------------
    cat("Creating Figure 08g: GO BP for top decile (D10)...\n")

    go_bp_top <- enrichGO(gene = top_ids$ENTREZID,
                          OrgDb = org.Mm.eg.db,
                          ont = "BP",
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.05,
                          readable = TRUE)

    if (!is.null(go_bp_top) && nrow(go_bp_top@result) > 0) {
      p_08g <- dotplot(go_bp_top, showCategory = 15) +
        labs(title = "GO Biological Process: Top Decile (D10)",
             subtitle = "Most TET-impaired genes (most negative delta-ratio)") +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"))

      save_multiformat_ggplot(p_08g, file.path(OUTPUT_DIR, "08g_enrichment_delta_ratio_top_decile_go_bp"),
                              width = 10, height = 10)

      write.table(go_bp_top@result,
                  file.path(TABLES_DIR, "enrichment_delta_ratio_top_decile_go_bp.tsv"),
                  sep = "\t", quote = FALSE, row.names = FALSE)
      cat("  Saved D10 GO BP results\n")
    } else {
      cat("  No significant GO BP terms for D10\n")
    }

    # -----------------------------------------------------------------------
    # FIGURE 08h: Standalone GO BP — Bottom Decile (D1)
    # -----------------------------------------------------------------------
    cat("Creating Figure 08h: GO BP for bottom decile (D1)...\n")

    go_bp_bottom <- enrichGO(gene = bottom_ids$ENTREZID,
                             OrgDb = org.Mm.eg.db,
                             ont = "BP",
                             pAdjustMethod = "BH",
                             qvalueCutoff = 0.05,
                             readable = TRUE)

    if (!is.null(go_bp_bottom) && nrow(go_bp_bottom@result) > 0) {
      p_08h <- dotplot(go_bp_bottom, showCategory = 15) +
        labs(title = "GO Biological Process: Bottom Decile (D1)",
             subtitle = "Least TET-impaired genes (most positive delta-ratio)") +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"))

      save_multiformat_ggplot(p_08h, file.path(OUTPUT_DIR, "08h_enrichment_delta_ratio_bottom_decile_go_bp"),
                              width = 10, height = 10)

      write.table(go_bp_bottom@result,
                  file.path(TABLES_DIR, "enrichment_delta_ratio_bottom_decile_go_bp.tsv"),
                  sep = "\t", quote = FALSE, row.names = FALSE)
      cat("  Saved D1 GO BP results\n")
    } else {
      cat("  No significant GO BP terms for D1\n")
    }

  } else {
    cat("  Not enough genes converted in one or both deciles for enrichment analysis\n")
  }
}, error = function(e) {
  cat(sprintf("  Delta-ratio decile enrichment error: %s\n", e$message))
})

cat("Section 8 complete.\n\n")
