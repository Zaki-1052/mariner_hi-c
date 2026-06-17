# biomodal/downstream/scripts/viz_sections/section_61j_peak_level_enrichment.R
# Section 61j: Peak-level GO enrichment for MeCP2-Up peaks
#
# Annotates 7,686 MeCP2-Up DiffBind peaks to nearest genes, then runs GO
# enrichment on the unique gene set. Compares: all MeCP2-Up peaks vs
# MeCP2-Up + K119ub-Up peaks.
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_61j_peak_level_enrichment.R

source("scripts/viz_sections/_shared_config.R")

suppressPackageStartupMessages({
  library(ChIPseeker)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
  library(patchwork)
})

cat("================================================================================\n")
cat("SECTION 61j: PEAK-LEVEL MeCP2-Up GO ENRICHMENT\n")
cat("  7,686 peaks -> nearest genes -> GO BP\n")
cat("================================================================================\n\n")

SEC61J_DIR <- file.path(OUTPUT_DIR, "61j_peak_level_enrichment")
dir.create(SEC61J_DIR, recursive = TRUE, showWarnings = FALSE)

SIGNAL_PATH <- file.path(TABLES_DIR, "62_mecp2_peak_chromatin_signal.tsv")
if (!file.exists(SIGNAL_PATH)) stop("Run section 62 first: ", SIGNAL_PATH)

save_plot <- function(p, name, w = 12, h = 10) {
  save_multiformat_ggplot(p,
    base_path = file.path(SEC61J_DIR, name),
    width = w, height = h, dpi = 300,
    verbose = TRUE, use_subfolders = TRUE
  )
}

# --- Load peaks ---

peaks <- read.table(SIGNAL_PATH, header = TRUE, sep = "\t",
                    stringsAsFactors = FALSE, quote = "")
cat(sprintf("  Loaded %d peaks from signal cache\n", nrow(peaks)))

peaks$mecp2_up   <- peaks$FDR < 0.05 & peaks$Fold > 0
peaks$mecp2_down <- peaks$FDR < 0.05 & peaks$Fold < 0
peaks$k119ub_up  <- !is.na(peaks$k119ub_log2fc) & peaks$k119ub_log2fc > 0

n_up   <- sum(peaks$mecp2_up)
n_down <- sum(peaks$mecp2_down)
n_up_k119 <- sum(peaks$mecp2_up & peaks$k119ub_up)

cat(sprintf("  MeCP2-Up peaks:              %d\n", n_up))
cat(sprintf("  MeCP2-Down peaks:            %d\n", n_down))
cat(sprintf("  MeCP2-Up + K119ub-Up peaks:  %d\n", n_up_k119))

# --- Annotate peaks to nearest gene ---

cat("\n--- Annotating peaks to nearest gene ---\n")

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene

annotate_peaks_to_genes <- function(peak_df, label) {
  gr <- GRanges(seqnames = peak_df$Chr,
                ranges = IRanges(start = peak_df$Start, end = peak_df$End))
  ann <- annotatePeak(gr, TxDb = txdb, level = "gene",
                      annoDb = "org.Mm.eg.db", verbose = FALSE)
  ann_df <- as.data.frame(ann)
  genes <- unique(na.omit(ann_df$SYMBOL))
  cat(sprintf("  %s: %d peaks -> %d unique genes\n", label, nrow(peak_df), length(genes)))
  list(genes = genes, annotation = ann_df)
}

up_result      <- annotate_peaks_to_genes(peaks[peaks$mecp2_up, ], "MeCP2-Up")
up_k119_result <- annotate_peaks_to_genes(peaks[peaks$mecp2_up & peaks$k119ub_up, ], "MeCP2-Up+K119ub-Up")
down_result    <- annotate_peaks_to_genes(peaks[peaks$mecp2_down, ], "MeCP2-Down")
all_result     <- annotate_peaks_to_genes(peaks, "All peaks (background)")

up_genes      <- up_result$genes
up_k119_genes <- up_k119_result$genes
down_genes    <- down_result$genes
bg_genes      <- all_result$genes

# --- Convert to Entrez ---

cat("\n--- Converting to Entrez IDs ---\n")

up_entrez <- bitr(up_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
up_k119_entrez <- bitr(up_k119_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
down_entrez <- bitr(down_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
bg_entrez <- bitr(bg_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

cat(sprintf("  MeCP2-Up:             %d Entrez IDs\n", nrow(up_entrez)))
cat(sprintf("  MeCP2-Up+K119ub-Up:   %d Entrez IDs\n", nrow(up_k119_entrez)))
cat(sprintf("  MeCP2-Down:           %d Entrez IDs\n", nrow(down_entrez)))
cat(sprintf("  Background (all):     %d Entrez IDs\n", nrow(bg_entrez)))

# =============================================================================
# GO ENRICHMENT: 4 analyses
# =============================================================================

run_enrichment <- function(gene_entrez, bg_entrez_ids, label, use_bg = TRUE) {
  cat(sprintf("\n--- %s ---\n", label))

  ego <- enrichGO(
    gene         = gene_entrez,
    universe     = if (use_bg) bg_entrez_ids else NULL,
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
    top <- head(ego@result[ego@result$p.adjust < 0.05, ], 25)
    cat("\n  Top terms:\n")
    for (i in seq_len(nrow(top))) {
      r <- top[i, ]
      cat(sprintf("    %-55s q=%.2e  %s\n",
                  substr(r$Description, 1, 55), r$p.adjust, r$GeneRatio))
    }

    neuronal_mask <- grepl("synap|neuron|axon|dendrit|nervous|brain|cerebel",
                           top$Description, ignore.case = TRUE)
    cat(sprintf("\n  Neuronal in top 25: %d (%.0f%%)\n",
                sum(neuronal_mask), 100 * sum(neuronal_mask) / nrow(top)))
  }

  list(ego = ego, n_sig = n_sig)
}

# 1. MeCP2-Up peaks, custom background
r1 <- run_enrichment(up_entrez$ENTREZID, bg_entrez$ENTREZID,
                     sprintf("MeCP2-Up peaks (%d genes), custom bg", length(up_genes)))

# 2. MeCP2-Up peaks, genome-wide background
r2 <- run_enrichment(up_entrez$ENTREZID, NULL,
                     sprintf("MeCP2-Up peaks (%d genes), genome-wide bg", length(up_genes)),
                     use_bg = FALSE)

# 3. MeCP2-Up + K119ub-Up peaks, custom background
r3 <- run_enrichment(up_k119_entrez$ENTREZID, bg_entrez$ENTREZID,
                     sprintf("MeCP2-Up+K119ub-Up peaks (%d genes), custom bg", length(up_k119_genes)))

# 4. MeCP2-Up + K119ub-Up peaks, genome-wide background
r4 <- run_enrichment(up_k119_entrez$ENTREZID, NULL,
                     sprintf("MeCP2-Up+K119ub-Up peaks (%d genes), genome-wide bg", length(up_k119_genes)),
                     use_bg = FALSE)

# 5. MeCP2-Down for comparison (genome-wide)
r5 <- run_enrichment(down_entrez$ENTREZID, NULL,
                     sprintf("MeCP2-Down peaks (%d genes), genome-wide bg", length(down_genes)),
                     use_bg = FALSE)

# =============================================================================
# SAVE TABLES + PLOTS
# =============================================================================

cat("\n--- Saving results ---\n")

save_results <- list(
  "61j_mecp2up_custom_go_bp"      = r1,
  "61j_mecp2up_genomewide_go_bp"  = r2,
  "61j_mecp2up_k119up_custom_go_bp" = r3,
  "61j_mecp2up_k119up_genomewide_go_bp" = r4,
  "61j_mecp2down_genomewide_go_bp" = r5
)

for (nm in names(save_results)) {
  write.table(save_results[[nm]]$ego@result,
              file.path(TABLES_DIR, paste0(nm, ".tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat(sprintf("  Saved: %s.tsv\n", nm))
}

# Gene lists
write.table(data.frame(gene = up_genes), file.path(TABLES_DIR, "61j_mecp2up_peak_genes.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(gene = up_k119_genes), file.path(TABLES_DIR, "61j_mecp2up_k119up_peak_genes.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("  Saved gene lists: %d MeCP2-Up genes, %d MeCP2-Up+K119ub-Up genes\n",
            length(up_genes), length(up_k119_genes)))

# Dot plots for significant analyses
plot_list <- list()
plot_labels <- c(
  "MeCP2-Up, custom bg", "MeCP2-Up, genome-wide",
  "Up+K119ub, custom bg", "Up+K119ub, genome-wide",
  "MeCP2-Down, genome-wide"
)
plot_files <- c("61j_mecp2up_custom_dotplot", "61j_mecp2up_genomewide_dotplot",
                "61j_up_k119up_custom_dotplot", "61j_up_k119up_genomewide_dotplot",
                "61j_mecp2down_genomewide_dotplot")

for (i in seq_along(save_results)) {
  r <- save_results[[i]]
  if (r$n_sig >= 3) {
    p <- dotplot(r$ego, showCategory = min(20, r$n_sig)) +
      labs(title = sprintf("GO BP: %s", plot_labels[i])) +
      theme_biomodal()
    save_plot(p, plot_files[i])
    plot_list[[plot_labels[i]]] <- p
  }
}

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n================================================================================\n")
cat("SECTION 61j SUMMARY: PEAK-LEVEL ENRICHMENT\n")
cat("================================================================================\n\n")

summary_df <- data.frame(
  analysis = plot_labels,
  n_genes = c(length(up_genes), length(up_genes),
              length(up_k119_genes), length(up_k119_genes),
              length(down_genes)),
  n_peaks = c(n_up, n_up, n_up_k119, n_up_k119, n_down),
  sig_go_terms = sapply(save_results, function(r) r$n_sig),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(summary_df))) {
  r <- summary_df[i, ]
  cat(sprintf("  %-35s %5d peaks -> %4d genes -> %4d GO terms\n",
              r$analysis, r$n_peaks, r$n_genes, r$sig_go_terms))
}

cat("\n  Plots saved to:", SEC61J_DIR, "\n")
cat("  Tables saved to:", TABLES_DIR, "\n")
cat("\nSection 61j complete.\n")
