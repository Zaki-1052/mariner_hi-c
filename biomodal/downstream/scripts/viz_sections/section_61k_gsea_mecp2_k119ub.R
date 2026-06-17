# biomodal/downstream/scripts/viz_sections/section_61k_gsea_mecp2_k119ub.R
# Section 61k: GSEA on MeCP2 and K119ub ranked gene lists
#
# Ranks ALL genes by continuous metrics (no threshold), then tests whether
# neuronal/synaptic gene sets are enriched at the extremes.
# GSEA handles the weighting problem: genes with stronger MeCP2 effects
# naturally rank higher and contribute more.
#
# Rankings:
#   1. mecp2_mean_fold — MeCP2 binding change
#   2. K119ub gb_log2fc — ubiquitination change
#   3. Combined: mecp2 * k119ub interaction
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_61k_gsea_mecp2_k119ub.R

source("scripts/viz_sections/_shared_config.R")

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
  library(patchwork)
})

cat("================================================================================\n")
cat("SECTION 61k: GSEA — MeCP2 AND K119ub RANKED GENE LISTS\n")
cat("  No thresholds. All genes ranked. Proper weighting.\n")
cat("================================================================================\n\n")

SEC61K_DIR <- file.path(OUTPUT_DIR, "61k_gsea_mecp2_k119ub")
dir.create(SEC61K_DIR, recursive = TRUE, showWarnings = FALSE)

save_plot <- function(p, name, w = 12, h = 10) {
  save_multiformat_ggplot(p,
    base_path = file.path(SEC61K_DIR, name),
    width = w, height = h, dpi = 300,
    verbose = TRUE, use_subfolders = TRUE
  )
}

# --- Load data ---

qm <- read.table(file.path(TABLES_DIR, "59_quadrant_master.tsv"),
                  header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")
cat(sprintf("  Loaded %d genes from quadrant master\n", nrow(qm)))

# --- Build ranked lists ---

cat("\n--- Building ranked gene lists ---\n")

build_ranked_list <- function(genes, values, label) {
  valid <- !is.na(genes) & !is.na(values) & genes != ""
  g <- genes[valid]
  v <- values[valid]

  entrez <- bitr(g, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

  merged <- data.frame(SYMBOL = g, value = v, stringsAsFactors = FALSE)
  merged <- merge(merged, entrez, by = "SYMBOL")

  # Deduplicate: keep max absolute value per Entrez ID
  merged <- merged[order(-abs(merged$value)), ]
  merged <- merged[!duplicated(merged$ENTREZID), ]

  ranked <- setNames(merged$value, merged$ENTREZID)
  ranked <- sort(ranked, decreasing = TRUE)

  cat(sprintf("  %s: %d genes ranked (range: %.3f to %.3f)\n",
              label, length(ranked), max(ranked), min(ranked)))
  ranked
}

rank_mecp2 <- build_ranked_list(qm$gene, qm$mecp2_mean_fold, "MeCP2 fold")

k119_mask <- qm$gb_signal_class == "quantifiable"
rank_k119ub <- build_ranked_list(qm$gene[k119_mask], qm$gb_log2fc[k119_mask], "K119ub log2fc")

# --- Run GSEA ---

run_gsea <- function(ranked_list, label) {
  cat(sprintf("\n--- GSEA: %s ---\n", label))

  gsea_result <- gseGO(
    geneList     = ranked_list,
    OrgDb        = org.Mm.eg.db,
    ont          = "BP",
    minGSSize    = 15,
    maxGSSize    = 500,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    verbose      = FALSE,
    seed         = TRUE
  )

  n_sig <- sum(gsea_result@result$p.adjust < 0.05)
  n_pos <- sum(gsea_result@result$p.adjust < 0.05 & gsea_result@result$NES > 0)
  n_neg <- sum(gsea_result@result$p.adjust < 0.05 & gsea_result@result$NES < 0)

  cat(sprintf("  Significant terms (q < 0.05): %d (positive NES: %d, negative: %d)\n",
              n_sig, n_pos, n_neg))

  if (n_sig > 0) {
    sig_results <- gsea_result@result[gsea_result@result$p.adjust < 0.05, ]
    sig_results <- sig_results[order(-abs(sig_results$NES)), ]
    top <- head(sig_results, 25)

    cat("\n  Top terms by |NES|:\n")
    for (i in seq_len(nrow(top))) {
      r <- top[i, ]
      dir <- ifelse(r$NES > 0, "+", "-")
      cat(sprintf("    [%s] %-50s NES=%+.2f  q=%.2e  size=%d\n",
                  dir, substr(r$Description, 1, 50), r$NES, r$p.adjust, r$setSize))
    }

    neuronal_mask <- grepl("synap|neuron|axon|dendrit|nervous|brain|cerebel",
                           sig_results$Description, ignore.case = TRUE)
    n_neuronal <- sum(neuronal_mask)
    cat(sprintf("\n  Total neuronal terms (q<0.05): %d / %d significant\n",
                n_neuronal, n_sig))

    if (n_neuronal > 0) {
      neur_terms <- sig_results[neuronal_mask, ]
      neur_terms <- neur_terms[order(-abs(neur_terms$NES)), ]
      cat("  Neuronal terms:\n")
      for (i in seq_len(min(15, nrow(neur_terms)))) {
        r <- neur_terms[i, ]
        dir <- ifelse(r$NES > 0, "+", "-")
        cat(sprintf("    [%s] %-50s NES=%+.2f  q=%.2e\n",
                    dir, substr(r$Description, 1, 50), r$NES, r$p.adjust))
      }
    }
  }

  list(result = gsea_result, n_sig = n_sig)
}

gsea_mecp2  <- run_gsea(rank_mecp2, "MeCP2 fold change")
gsea_k119ub <- run_gsea(rank_k119ub, "K119ub log2FC")

# =============================================================================
# PLOTS
# =============================================================================

cat("\n--- Saving plots ---\n")

# Dotplots
if (gsea_mecp2$n_sig >= 3) {
  p_mecp2 <- dotplot(gsea_mecp2$result, showCategory = 20, split = ".sign") +
    facet_grid(. ~ .sign) +
    labs(title = "GSEA: Genes ranked by MeCP2 fold change",
         subtitle = sprintf("%d significant terms", gsea_mecp2$n_sig)) +
    theme_biomodal(base_size = 10)
  save_plot(p_mecp2, "61k_gsea_mecp2_dotplot", w = 14, h = 10)
}

if (gsea_k119ub$n_sig >= 3) {
  p_k119 <- dotplot(gsea_k119ub$result, showCategory = 20, split = ".sign") +
    facet_grid(. ~ .sign) +
    labs(title = "GSEA: Genes ranked by K119ub log2FC",
         subtitle = sprintf("%d significant terms", gsea_k119ub$n_sig)) +
    theme_biomodal(base_size = 10)
  save_plot(p_k119, "61k_gsea_k119ub_dotplot", w = 14, h = 10)
}

# GSEA running score plots for top neuronal terms
for (gsea_name in c("mecp2", "k119ub")) {
  gsea_obj <- if (gsea_name == "mecp2") gsea_mecp2$result else gsea_k119ub$result

  sig_terms <- gsea_obj@result[gsea_obj@result$p.adjust < 0.05, ]
  neur_terms <- sig_terms[grepl("synap|neuron|axon", sig_terms$Description,
                                ignore.case = TRUE), ]

  if (nrow(neur_terms) > 0) {
    top_neur_id <- neur_terms$ID[which.max(abs(neur_terms$NES))]
    top_neur_desc <- neur_terms$Description[which.max(abs(neur_terms$NES))]

    p_running <- gseaplot2(gsea_obj, geneSetID = top_neur_id, title = top_neur_desc,
                           pvalue_table = TRUE)
    save_plot(p_running, sprintf("61k_%s_running_%s", gsea_name,
              gsub("[^a-zA-Z0-9]", "_", substr(top_neur_desc, 1, 40))), w = 10, h = 7)
  }
}

# =============================================================================
# SAVE TABLES
# =============================================================================

cat("\n--- Saving tables ---\n")

write.table(gsea_mecp2$result@result,
            file.path(TABLES_DIR, "61k_gsea_mecp2_go_bp.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("  Saved: 61k_gsea_mecp2_go_bp.tsv\n")

write.table(gsea_k119ub$result@result,
            file.path(TABLES_DIR, "61k_gsea_k119ub_go_bp.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("  Saved: 61k_gsea_k119ub_go_bp.tsv\n")

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n================================================================================\n")
cat("SECTION 61k SUMMARY: GSEA\n")
cat("================================================================================\n\n")
cat(sprintf("  MeCP2 ranking:  %d genes, %d significant GO terms\n",
            length(rank_mecp2), gsea_mecp2$n_sig))
cat(sprintf("  K119ub ranking: %d genes, %d significant GO terms\n",
            length(rank_k119ub), gsea_k119ub$n_sig))
cat("\n  Key question: are neuronal gene sets enriched among genes with\n")
cat("  the HIGHEST MeCP2 and K119ub increases? GSEA answers this with\n")
cat("  proper weighting — no arbitrary thresholds.\n")
cat("\nSection 61k complete.\n")
