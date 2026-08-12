# cluster/scripts/12_go_kegg_enrichment.R
#
# Phase 9: Per-cluster GO/KEGG pathway enrichment for Hi-C loop anchor genes.
# Tests whether clust5 (Polycomb gain) and clust6 (anchor disruption) gene
# sets are enriched for distinct biological pathways -- developmental/Polycomb
# targets vs synaptic/neural-function terms.
#
# Input:  outputs/bap1_late/figures/annotation/clust{1..6}_annotation.txt
# Output: outputs/bap1_late/figures/go_enrichment/{TSVs + dotplots}
#
# Run via: bash scripts/run_go_kegg.sh
# Or directly: /usr/local/bin/Rscript scripts/12_go_kegg_enrichment.R

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
  library(ggplot2)
})

cat("=== Phase 9: Per-Cluster GO/KEGG Enrichment ===\n")
cat(sprintf("  Started: %s\n", Sys.time()))
cat(sprintf("  clusterProfiler: %s\n", packageVersion("clusterProfiler")))
cat(sprintf("  org.Mm.eg.db: %s\n", packageVersion("org.Mm.eg.db")))
cat("\n")

# =============================================================================
# CONFIGURATION
# =============================================================================

ANNOT_DIR  <- "outputs/bap1_late/figures/annotation"
OUTPUT_DIR <- "outputs/bap1_late/figures/go_enrichment"
UTIL_PATH  <- "scripts/utils/multi_format_output.R"
CLUSTERS   <- paste0("clust", 1:6)

stopifnot(file.exists(ANNOT_DIR))
stopifnot(file.exists(UTIL_PATH))
source(UTIL_PATH)

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

theme_pub <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", color = NA),
    plot.title = element_text(size = 12, face = "bold"),
    legend.position = "bottom"
  )

save_plot <- function(p, name, w = 10, h = 8) {
  save_multiformat_ggplot(p, base_path = file.path(OUTPUT_DIR, name),
                          width = w, height = h, dpi = 300,
                          verbose = TRUE, use_subfolders = TRUE)
}

has_results <- function(res) {
  !is.null(res) && nrow(as.data.frame(res)) > 0
}

# =============================================================================
# LOAD GENE LISTS
# =============================================================================

cat("Loading per-cluster gene lists...\n")

load_cluster_genes <- function(cluster_name) {
  fpath <- file.path(ANNOT_DIR, paste0(cluster_name, "_annotation.txt"))
  stopifnot(file.exists(fpath))
  df <- read.table(fpath, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
                   quote = "", comment.char = "", fill = TRUE)
  ids <- unique(c(as.character(df$chr1_gene_id), as.character(df$chr2_gene_id)))
  ids <- ids[!is.na(ids) & nchar(trimws(ids)) > 0]
  cat(sprintf("  %s: %d unique Entrez IDs\n", cluster_name, length(ids)))
  ids
}

gene_lists <- lapply(setNames(CLUSTERS, CLUSTERS), load_cluster_genes)
bg_entrez  <- unique(unlist(gene_lists))
cat(sprintf("  Background universe: %d unique Entrez IDs (all anchor-proximal genes)\n\n",
            length(bg_entrez)))

for (cl in CLUSTERS) {
  if (length(gene_lists[[cl]]) < 10) {
    cat(sprintf("  WARNING: %s has only %d genes -- enrichment may be sparse\n",
                cl, length(gene_lists[[cl]])))
  }
}

# =============================================================================
# SECTION 1: compareCluster -- GO BP (all 6 clusters)
# =============================================================================

cat("Running compareCluster GO BP (all 6 clusters)...\n")

cc_bp <- tryCatch(
  compareCluster(
    geneClusters = gene_lists,
    fun          = "enrichGO",
    OrgDb        = org.Mm.eg.db,
    ont          = "BP",
    keyType      = "ENTREZID",
    universe     = bg_entrez,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.10,
    readable     = TRUE
  ),
  error = function(e) { cat("  ERROR:", e$message, "\n"); NULL }
)

if (has_results(cc_bp)) {
  bp_df <- as.data.frame(cc_bp)
  write.table(bp_df, file.path(OUTPUT_DIR, "go_bp_all_clusters.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  Total GO BP terms: %d\n", nrow(bp_df)))
  for (cl in CLUSTERS) {
    n <- sum(bp_df$Cluster == cl)
    if (n > 0) cat(sprintf("    %s: %d terms\n", cl, n))
  }

  p_bp <- dotplot(cc_bp, showCategory = 15,
                  title = "GO Biological Process: Per-Cluster Enrichment") +
    scale_color_continuous(low = "#E41A1C", high = "#377EB8", name = "p.adjust") +
    theme_pub +
    theme(axis.text.y = element_text(size = 9))
  save_plot(p_bp, "go_bp_all_clusters", w = 14, h = 10)
} else {
  cat("  No significant GO BP terms across all clusters.\n")
}
cat("\n")

# =============================================================================
# SECTION 2: compareCluster -- KEGG (all 6 clusters)
# =============================================================================

cat("Running compareCluster KEGG (all 6 clusters)...\n")
cat("  NOTE: KEGG requires internet (REST API) -- will skip if offline\n")

cc_kegg <- tryCatch(
  compareCluster(
    geneClusters = gene_lists,
    fun          = "enrichKEGG",
    organism     = "mmu",
    universe     = bg_entrez,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.10
  ),
  error = function(e) { cat("  ERROR:", e$message, "\n"); NULL }
)

if (has_results(cc_kegg)) {
  cc_kegg_read <- setReadable(cc_kegg, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
  kegg_df <- as.data.frame(cc_kegg_read)
  write.table(kegg_df, file.path(OUTPUT_DIR, "kegg_all_clusters.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  Total KEGG pathways: %d\n", nrow(kegg_df)))
  for (cl in CLUSTERS) {
    n <- sum(kegg_df$Cluster == cl)
    if (n > 0) cat(sprintf("    %s: %d pathways\n", cl, n))
  }

  p_kegg <- dotplot(cc_kegg_read, showCategory = 15,
                    title = "KEGG Pathways: Per-Cluster Enrichment") +
    scale_color_continuous(low = "#E41A1C", high = "#377EB8", name = "p.adjust") +
    theme_pub +
    theme(axis.text.y = element_text(size = 9))
  save_plot(p_kegg, "kegg_all_clusters", w = 14, h = 10)
} else {
  cat("  No significant KEGG pathways across all clusters.\n")
}
cat("\n")

# =============================================================================
# SECTION 3: Focused GO BP -- clust5 (strong gain, 97% up, Polycomb)
# =============================================================================

cat("Running focused enrichment for clust5 (strong gain, 97% up)...\n")
cat(sprintf("  Gene list size: %d\n", length(gene_lists[["clust5"]])))

ego5 <- tryCatch(
  enrichGO(
    gene          = gene_lists[["clust5"]],
    universe      = bg_entrez,
    OrgDb         = org.Mm.eg.db,
    ont           = "BP",
    keyType       = "ENTREZID",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.10,
    readable      = TRUE
  ),
  error = function(e) { cat("  GO BP error:", e$message, "\n"); NULL }
)

if (has_results(ego5)) {
  ego5_df <- as.data.frame(ego5)
  write.table(ego5_df, file.path(OUTPUT_DIR, "go_bp_clust5.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  GO BP: %d terms\n", nrow(ego5_df)))
  p5_go <- dotplot(ego5, showCategory = 20,
                   title = "GO BP: clust5 (strong gain, 97% up, Polycomb)") +
    theme_pub
  save_plot(p5_go, "go_bp_clust5", w = 10, h = 8)
} else {
  cat("  clust5: 0 significant GO BP terms (237 genes, sparse -- expected)\n")
}

# =============================================================================
# SECTION 4: Focused KEGG -- clust5
# =============================================================================

ekegg5 <- tryCatch(
  enrichKEGG(
    gene          = gene_lists[["clust5"]],
    universe      = bg_entrez,
    organism      = "mmu",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.10
  ),
  error = function(e) { cat("  KEGG error:", e$message, "\n"); NULL }
)

if (has_results(ekegg5)) {
  ekegg5_read <- setReadable(ekegg5, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
  ekegg5_df <- as.data.frame(ekegg5_read)
  write.table(ekegg5_df, file.path(OUTPUT_DIR, "kegg_clust5.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  KEGG: %d pathways\n", nrow(ekegg5_df)))
  p5_kegg <- dotplot(ekegg5_read, showCategory = 20,
                     title = "KEGG: clust5 (strong gain, 97% up, Polycomb)") +
    theme_pub
  save_plot(p5_kegg, "kegg_clust5", w = 10, h = 8)
} else {
  cat("  clust5: 0 significant KEGG pathways\n")
}
cat("\n")

# =============================================================================
# SECTION 5: Focused GO BP -- clust6 (strong loss, 78% down)
# =============================================================================

cat("Running focused enrichment for clust6 (strong loss, 78% down)...\n")
cat(sprintf("  Gene list size: %d\n", length(gene_lists[["clust6"]])))

ego6 <- tryCatch(
  enrichGO(
    gene          = gene_lists[["clust6"]],
    universe      = bg_entrez,
    OrgDb         = org.Mm.eg.db,
    ont           = "BP",
    keyType       = "ENTREZID",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.10,
    readable      = TRUE
  ),
  error = function(e) { cat("  GO BP error:", e$message, "\n"); NULL }
)

if (has_results(ego6)) {
  ego6_df <- as.data.frame(ego6)
  write.table(ego6_df, file.path(OUTPUT_DIR, "go_bp_clust6.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  GO BP: %d terms\n", nrow(ego6_df)))
  p6_go <- dotplot(ego6, showCategory = 20,
                   title = "GO BP: clust6 (strong loss, 78% down)") +
    theme_pub
  save_plot(p6_go, "go_bp_clust6", w = 10, h = 8)
} else {
  cat("  clust6: 0 significant GO BP terms\n")
}

# =============================================================================
# SECTION 6: Focused KEGG -- clust6
# =============================================================================

ekegg6 <- tryCatch(
  enrichKEGG(
    gene          = gene_lists[["clust6"]],
    universe      = bg_entrez,
    organism      = "mmu",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.10
  ),
  error = function(e) { cat("  KEGG error:", e$message, "\n"); NULL }
)

if (has_results(ekegg6)) {
  ekegg6_read <- setReadable(ekegg6, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
  ekegg6_df <- as.data.frame(ekegg6_read)
  write.table(ekegg6_df, file.path(OUTPUT_DIR, "kegg_clust6.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  KEGG: %d pathways\n", nrow(ekegg6_df)))
  p6_kegg <- dotplot(ekegg6_read, showCategory = 20,
                     title = "KEGG: clust6 (strong loss, 78% down)") +
    theme_pub
  save_plot(p6_kegg, "kegg_clust6", w = 10, h = 8)
} else {
  cat("  clust6: 0 significant KEGG pathways\n")
}
cat("\n")

# =============================================================================
# SUMMARY
# =============================================================================

cat("=== Summary ===\n\n")

cat("Gene list sizes:\n")
for (cl in CLUSTERS) {
  cat(sprintf("  %s: %d genes\n", cl, length(gene_lists[[cl]])))
}
cat(sprintf("  Background: %d genes\n\n", length(bg_entrez)))

cat("Enrichment results:\n")

summarize_result <- function(label, res) {
  if (has_results(res)) {
    cat(sprintf("  %s: %d terms/pathways\n", label, nrow(as.data.frame(res))))
  } else {
    cat(sprintf("  %s: 0 (no significant terms)\n", label))
  }
}

summarize_result("compareCluster GO BP", cc_bp)
summarize_result("compareCluster KEGG", cc_kegg)
summarize_result("clust5 GO BP", ego5)
summarize_result("clust5 KEGG", ekegg5)
summarize_result("clust6 GO BP", ego6)
summarize_result("clust6 KEGG", ekegg6)
cat("\n")

print_top_terms <- function(label, res, n = 3) {
  if (!has_results(res)) return(invisible(NULL))
  df <- as.data.frame(res)
  top <- head(df$Description, n)
  if (length(top) > 0) {
    cat(sprintf("Top %d %s:\n", min(n, length(top)), label))
    for (term in top) cat(sprintf("  - %s\n", term))
    cat("\n")
  }
}

print_top_terms("clust5 GO BP terms", ego5)
print_top_terms("clust5 KEGG pathways", ekegg5)
print_top_terms("clust6 GO BP terms", ego6)
print_top_terms("clust6 KEGG pathways", ekegg6)

cat(sprintf("Completed: %s\n", Sys.time()))
cat("=== Phase 9 GO/KEGG complete ===\n")
