# biomodal/downstream/scripts/viz_sections/section_08_enrichment.R
# Section 8: GO/KEGG Enrichment Analysis
# Standalone script - sources shared config for all dependencies and data

# Source shared config (handles both Rscript and source())
local({
  script_dir <- NULL
  args <- commandArgs(trailingOnly = FALSE)
  f <- grep("--file=", args, value = TRUE)
  if (length(f) > 0) script_dir <- dirname(normalizePath(sub("--file=", "", f)))
  if (is.null(script_dir)) for (i in seq_len(sys.nframe())) {
    fi <- sys.frame(i)$ofile
    if (!is.null(fi)) { script_dir <- dirname(normalizePath(fi)); break }
  }
  if (is.null(script_dir)) script_dir <- "/Users/zakiralibhai/Documents/GitHub/mariner_hi-c/biomodal/downstream/scripts/viz_sections"
  source(file.path(script_dir, "_shared_config.R"))
})

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

cat("Section 8 complete.\n\n")
