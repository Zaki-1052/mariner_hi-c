# biomodal/downstream/scripts/viz_sections/section_61jk_composites.R
# Composites for sections 61j (peak-level ORA) and 61k (GSEA).
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_61jk_composites.R

source("scripts/viz_sections/_shared_config.R")

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
  library(patchwork)
})

cat("================================================================================\n")
cat("SECTION 61j/k COMPOSITES\n")
cat("================================================================================\n\n")

save_plot <- function(p, dir, name, w = 18, h = 14) {
  save_multiformat_ggplot(p,
    base_path = file.path(dir, name),
    width = w, height = h, dpi = 300,
    verbose = TRUE, use_subfolders = TRUE
  )
}

# =============================================================================
# 61j COMPOSITE: Peak-level ORA (5 analyses)
# =============================================================================

cat("--- 61j Composite: Peak-level ORA ---\n\n")

SEC61J_DIR <- file.path(OUTPUT_DIR, "61j_peak_level_enrichment")

j_files <- list(
  "MeCP2-Up, custom bg"       = file.path(TABLES_DIR, "61j_mecp2up_custom_go_bp.tsv"),
  "MeCP2-Up, genome-wide"     = file.path(TABLES_DIR, "61j_mecp2up_genomewide_go_bp.tsv"),
  "Up+K119ub, custom bg"      = file.path(TABLES_DIR, "61j_mecp2up_k119up_custom_go_bp.tsv"),
  "Up+K119ub, genome-wide"    = file.path(TABLES_DIR, "61j_mecp2up_k119up_genomewide_go_bp.tsv"),
  "MeCP2-Down, genome-wide"   = file.path(TABLES_DIR, "61j_mecp2down_genomewide_go_bp.tsv")
)

j_panels <- list()
for (label in names(j_files)) {
  path <- j_files[[label]]
  if (!file.exists(path)) next

  df <- read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")
  n_sig <- sum(df$p.adjust < 0.05, na.rm = TRUE)
  if (n_sig < 3) {
    cat(sprintf("  %s: %d sig terms (skipping)\n", label, n_sig))
    next
  }

  top <- head(df[df$p.adjust < 0.05, ], 12)
  top$Description <- substr(top$Description, 1, 45)
  top <- top[order(top$p.adjust), ]
  top$Description <- factor(top$Description, levels = rev(top$Description))

  is_neuronal <- grepl("synap|neuron|axon|dendrit|nervous", top$Description, ignore.case = TRUE)

  p <- ggplot(top, aes(x = -log10(p.adjust), y = Description)) +
    geom_col(aes(fill = is_neuronal), alpha = 0.8, show.legend = FALSE) +
    scale_fill_manual(values = c("TRUE" = "#D73027", "FALSE" = "#4393C3")) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
    labs(title = label, x = "-log10(q)", y = NULL) +
    theme_biomodal(base_size = 9) +
    theme(plot.title = element_text(size = 10, face = "bold"))

  j_panels[[label]] <- p
  cat(sprintf("  %s: %d sig terms, %d neuronal\n", label, n_sig, sum(is_neuronal)))
}

if (length(j_panels) >= 2) {
  p_61j <- wrap_plots(j_panels, ncol = 2) +
    plot_annotation(
      title = "Section 61j: Peak-Level GO Enrichment (ORA)",
      subtitle = "7,686 MeCP2-Up peaks -> 2,107 genes | Red = neuronal terms",
      theme = theme(plot.title = element_text(size = 16, face = "bold"))
    )
  save_plot(p_61j, SEC61J_DIR, "61j_composite", w = 18, h = 14)
}

# =============================================================================
# 61k COMPOSITE: GSEA (MeCP2 + K119ub rankings)
# =============================================================================

cat("\n--- 61k Composite: GSEA ---\n\n")

SEC61K_DIR <- file.path(OUTPUT_DIR, "61k_gsea_mecp2_k119ub")

k_files <- list(
  mecp2  = file.path(TABLES_DIR, "61k_gsea_mecp2_go_bp.tsv"),
  k119ub = file.path(TABLES_DIR, "61k_gsea_k119ub_go_bp.tsv")
)

k_panels <- list()
for (ranking in names(k_files)) {
  path <- k_files[[ranking]]
  if (!file.exists(path)) next

  df <- read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")
  sig <- df[!is.na(df$p.adjust) & df$p.adjust < 0.05, ]

  if (nrow(sig) < 3) {
    cat(sprintf("  %s: %d sig terms (skipping)\n", ranking, nrow(sig)))
    next
  }

  # Top positive and negative NES
  sig <- sig[!is.na(sig$NES), ]
  pos_df <- sig[sig$NES > 0, ]
  pos_df <- pos_df[order(-pos_df$NES), ]
  neg_df <- sig[sig$NES < 0, ]
  neg_df <- neg_df[order(neg_df$NES), ]
  top <- rbind(head(pos_df, 8), head(neg_df, 8))
  top$Description <- substr(top$Description, 1, 45)
  top$Description <- make.unique(top$Description, sep = " ")
  top <- top[order(top$NES), ]
  top$Description <- factor(top$Description, levels = top$Description)

  top$is_neuronal <- grepl("synap|neuron|axon|dendrit|nervous|presynap|postsynap",
                           top$Description, ignore.case = TRUE)

  title_label <- ifelse(ranking == "mecp2", "Ranked by MeCP2 fold change",
                        "Ranked by K119ub log2FC")

  p <- ggplot(top, aes(x = NES, y = Description, fill = is_neuronal)) +
    geom_col(alpha = 0.8, show.legend = FALSE) +
    scale_fill_manual(values = c("TRUE" = "#D73027", "FALSE" = "#4393C3")) +
    geom_vline(xintercept = 0, color = "grey30") +
    labs(title = title_label,
         subtitle = sprintf("%d sig terms (pos: %d, neg: %d)",
                            nrow(sig), sum(sig$NES > 0), sum(sig$NES < 0)),
         x = "Normalized Enrichment Score (NES)", y = NULL) +
    theme_biomodal(base_size = 9) +
    theme(plot.title = element_text(size = 11, face = "bold"))

  k_panels[[ranking]] <- p
  cat(sprintf("  %s: %d sig terms, %d neuronal in top\n",
              ranking, nrow(sig), sum(is_neuronal)))
}

if (length(k_panels) == 2) {
  p_61k <- (k_panels[["mecp2"]] | k_panels[["k119ub"]]) +
    plot_annotation(
      title = "Section 61k: GSEA — All Genes Ranked (No Thresholds)",
      subtitle = "Red = neuronal/synaptic terms | Left = MeCP2 ranking | Right = K119ub ranking",
      theme = theme(plot.title = element_text(size = 16, face = "bold"))
    )
  save_plot(p_61k, SEC61K_DIR, "61k_composite", w = 18, h = 10)
}

cat("\nComposites complete.\n")
