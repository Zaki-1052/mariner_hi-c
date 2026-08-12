#!/usr/bin/env Rscript
# stripes/stripenn/scripts/deg_stripe_violin.R
#
# T1.3 - RNA-seq x stripe direction integration. Direct analog of
# scripts/deg_loop_anchor_violin.R, adapted for 1D stripes.
#
# Two gene pools per stripe:
#   - anchor-proximal: nearest_gene when distance_to_tss <= 2000 bp
#                      (matches the 7-category anchor classification threshold)
#   - body-resident:   body_genes column (comma-separated symbols)
#
# Usage:
#   Rscript scripts/deg_stripe_violin.R --timepoint late --resolution 5000
#
# Inputs:
#   - Stripes: stripes/stripenn/outputs/{tp}/visualizations/{tp}_annotated_stripes.tsv
#   - RNA-seq (late):  data/upstream/rna_seq/adult_rnaseq_results.xlsx
#     RNA-seq (early): data/upstream/rna_seq/young_rnaseq_results.xlsx
#     Columns: ensembl_gene_id (actually gene SYMBOLS), baseMean, log2FoldChange,
#              lfcSE, stat, pvalue, padj, and per-sample counts
#
# Outputs (outputs/stripe_integration/{label}/rna_integration/):
#   - stripe_gene_expression.tsv       Joined per-gene table
#   - wilcoxon_tests.tsv               Pairwise Wilcoxon by direction, faceted
#   - spearman_summary.tsv             Spearman rho: stripe logFC vs gene log2FC
#   - deg_stripe_violin                Violins (anchor vs body, by direction)
#   - stripe_rna_concordance           Scatter: stripe logFC vs body-gene log2FC
#   - summary.txt

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readxl)
  library(ggpubr)
})

if (dir.exists("/expanse/lustre/projects/csd940/zalibhai")) {
  CODE_DIR  <- "/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/stripes/stripenn"
  DATA_DIR  <- "/expanse/lustre/projects/csd940/zalibhai/stripes/stripenn"
  REPO_ROOT <- "/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c"
} else {
  CODE_DIR  <- normalizePath(file.path(getwd()))
  DATA_DIR  <- CODE_DIR
  REPO_ROOT <- normalizePath(file.path(CODE_DIR, "..", ".."))
}

source(file.path(REPO_ROOT, "scripts/utils/multi_format_output.R"))

args <- commandArgs(trailingOnly = TRUE)
parse_arg <- function(flag, default = NULL) {
  i <- which(args == flag); if (length(i) == 0) return(default)
  if (i == length(args)) stop(sprintf("Missing value for %s", flag))
  args[i + 1]
}
TP_LABEL   <- parse_arg("--timepoint", "late")
RESOLUTION <- as.integer(parse_arg("--resolution", "5000"))
TP_MAP <- list(late = "250402", early = "250831")
if (!TP_LABEL %in% names(TP_MAP)) stop("--timepoint must be 'late' or 'early'")
TP_ID <- TP_MAP[[TP_LABEL]]

STRIPE_FILE <- file.path(CODE_DIR, sprintf("outputs/%s/visualizations/%s_annotated_stripes.tsv",
                      TP_ID, TP_ID))
RNA_FILE <- if (TP_LABEL == "late") {
  file.path(REPO_ROOT, "data/upstream/rna_seq/adult_rnaseq_results.xlsx")
} else {
  file.path(REPO_ROOT, "data/upstream/rna_seq/young_rnaseq_results.xlsx")
}
OUTPUT_DIR <- file.path(REPO_ROOT, sprintf("outputs/stripe_integration/%s/rna_integration", TP_LABEL))
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Parameters (match scripts/deg_loop_anchor_violin.R)
PROMOTER_DIST_THRESHOLD <- 2000
DEG_PADJ_THRESHOLD <- 0.05
DEG_LFC_THRESHOLD  <- 0.3

cat("==============================================================\n")
cat("T1.3 - RNA-seq x Stripe Direction Violin Analysis\n")
cat("==============================================================\n")
cat(sprintf("Timepoint : %s (%s)\n", TP_LABEL, TP_ID))
cat(sprintf("Stripes   : %s\n", STRIPE_FILE))
cat(sprintf("RNA-seq   : %s\n", RNA_FILE))
cat(sprintf("Output    : %s\n\n", OUTPUT_DIR))

stopifnot(file.exists(STRIPE_FILE), file.exists(RNA_FILE))

# -------- LOAD --------
stripes <- read.table(STRIPE_FILE, header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, quote = "", comment.char = "",
                      na.strings = c("NA", ""))
rna <- as.data.frame(read_excel(RNA_FILE))
if (!"ensembl_gene_id" %in% colnames(rna)) {
  stop("Expected column 'ensembl_gene_id' (actually gene symbols) not found")
}
# Column is named ensembl_gene_id but contains gene SYMBOLS
# (comment confirmed in scripts/deg_loop_anchor_violin.R line 162)
rna <- rna %>%
  transmute(symbol = ensembl_gene_id,
            log2FoldChange = as.numeric(log2FoldChange),
            padj           = as.numeric(padj),
            baseMean       = as.numeric(baseMean)) %>%
  filter(!is.na(symbol))
cat(sprintf("RNA-seq: %d genes loaded\n", nrow(rna)))
cat(sprintf("Stripes: %d loaded\n", nrow(stripes)))

# -------- BUILD GENE -> STRIPE ASSIGNMENTS --------

# Anchor-proximal: single nearest gene when distance_to_tss <= 2000
anchor_assign <- stripes %>%
  transmute(stripe_id, direction, direction_confidence, significant_FDR05,
            stripe_logFC = logFC, stripe_FDR = FDR,
            anchor_type, symbol = nearest_gene,
            dist_to_tss = as.numeric(distance_to_tss)) %>%
  filter(!is.na(symbol), symbol != "",
         !is.na(dist_to_tss), dist_to_tss <= PROMOTER_DIST_THRESHOLD) %>%
  mutate(gene_pool = "anchor")

# Body-resident: explode comma-separated body_genes
body_assign <- stripes %>%
  filter(!is.na(body_genes), body_genes != "") %>%
  select(stripe_id, direction, direction_confidence, significant_FDR05,
         stripe_logFC = logFC, stripe_FDR = FDR, anchor_type, body_genes) %>%
  mutate(symbol = strsplit(body_genes, ",", fixed = TRUE)) %>%
  tidyr::unnest(symbol) %>%
  mutate(symbol = trimws(symbol), dist_to_tss = NA_real_) %>%
  filter(symbol != "") %>%
  select(-body_genes) %>%
  mutate(gene_pool = "body")

gene_assign <- bind_rows(anchor_assign, body_assign)
cat(sprintf("Gene-stripe pairs: %d anchor-proximal, %d body-resident\n",
            sum(gene_assign$gene_pool == "anchor"),
            sum(gene_assign$gene_pool == "body")))

# Join RNA-seq
joined <- gene_assign %>%
  left_join(rna, by = "symbol") %>%
  filter(!is.na(log2FoldChange))
cat(sprintf("After RNA-seq join: %d pairs with expression data\n", nrow(joined)))

joined$is_deg <- !is.na(joined$padj) &
  joined$padj < DEG_PADJ_THRESHOLD &
  abs(joined$log2FoldChange) > DEG_LFC_THRESHOLD

write.table(joined,
            file.path(OUTPUT_DIR, "stripe_gene_expression.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

# -------- WILCOXON TESTS (pairwise by direction, per pool) --------
target_dirs <- c("gained", "lost", "unchanged")
wx_tests <- list()
for (pool in c("anchor", "body")) {
  for (d1 in target_dirs) {
    for (d2 in target_dirs) {
      if (d1 >= d2) next
      x <- joined %>%
        filter(gene_pool == pool, direction == d1) %>%
        pull(log2FoldChange)
      y <- joined %>%
        filter(gene_pool == pool, direction == d2) %>%
        pull(log2FoldChange)
      x <- x[!is.na(x)]; y <- y[!is.na(y)]
      if (length(x) < 5 || length(y) < 5) next
      w <- suppressWarnings(wilcox.test(x, y))
      wx_tests[[length(wx_tests) + 1L]] <- data.frame(
        gene_pool = pool,
        group1 = d1, group2 = d2,
        n1 = length(x), n2 = length(y),
        median1 = median(x), median2 = median(y),
        W = as.numeric(w$statistic),
        p_value = w$p.value,
        stringsAsFactors = FALSE
      )
    }
  }
}
if (length(wx_tests) > 0) {
  wx <- do.call(rbind, wx_tests)
  wx$p_adj_BH <- p.adjust(wx$p_value, method = "BH")
  write.table(wx,
              file.path(OUTPUT_DIR, "wilcoxon_tests.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
} else {
  wx <- NULL
}

# -------- SPEARMAN: stripe logFC vs gene log2FC --------
# Restrict to high-confidence significant stripes (cleaner signal)
sp_lines <- list()
for (pool in c("anchor", "body")) {
  d <- joined %>%
    filter(gene_pool == pool,
           significant_FDR05 == TRUE,
           direction_confidence == "high",
           direction %in% c("gained", "lost"))
  if (nrow(d) < 20) next
  r <- suppressWarnings(cor.test(d$stripe_logFC, d$log2FoldChange,
                                 method = "spearman", exact = FALSE))
  sp_lines[[length(sp_lines) + 1L]] <- data.frame(
    gene_pool = pool,
    n = nrow(d),
    spearman_rho = as.numeric(r$estimate),
    p_value = r$p.value,
    stringsAsFactors = FALSE
  )
}
if (length(sp_lines) > 0) {
  sp <- do.call(rbind, sp_lines)
  sp$p_adj_BH <- p.adjust(sp$p_value, method = "BH")
  write.table(sp, file.path(OUTPUT_DIR, "spearman_summary.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
} else {
  sp <- NULL
}

# -------- VIOLIN PLOT --------
plot_df <- joined %>%
  filter(direction %in% target_dirs)

if (nrow(plot_df) >= 30) {
  p <- ggplot(plot_df,
              aes(x = direction, y = log2FoldChange, fill = direction)) +
    geom_violin(trim = FALSE, alpha = 0.7, scale = "width") +
    geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    facet_wrap(~ gene_pool, labeller = labeller(gene_pool = c(
      anchor = "Anchor-proximal (<=2kb to TSS)",
      body   = "Body-resident"))) +
    scale_fill_manual(values = c(gained = "#D73027", lost = "#4575B4",
                                 unchanged = "#999999")) +
    labs(title = sprintf("Gene expression log2FC by stripe direction (%s)", TP_LABEL),
         x = "Stripe direction",
         y = "RNA-seq log2(mut/ctrl)") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "none",
          strip.background = element_rect(fill = "grey90"))

  # Add pairwise significance if we have wilcox results
  if (!is.null(wx)) {
    # Use ggpubr for concise pairwise bracket layer
    p <- p + stat_compare_means(
      comparisons = list(c("gained", "lost"),
                         c("gained", "unchanged"),
                         c("lost", "unchanged")),
      method = "wilcox.test", label = "p.format", size = 3)
  }
  save_multiformat_ggplot(p,
                          file.path(OUTPUT_DIR, "deg_stripe_violin"),
                          width = 10, height = 6, use_subfolders = TRUE)
}

# -------- CONCORDANCE SCATTER --------
scat_df <- joined %>%
  filter(gene_pool == "body",
         significant_FDR05 == TRUE,
         direction_confidence == "high",
         direction %in% c("gained", "lost"))
if (nrow(scat_df) >= 20) {
  rho_body <- sp[sp$gene_pool == "body", "spearman_rho"]
  p_body   <- sp[sp$gene_pool == "body", "p_value"]
  rho_txt <- if (length(rho_body) == 1)
    sprintf("Spearman rho = %.3f (p = %.2g, n = %d)",
            rho_body, p_body, nrow(scat_df)) else ""

  p2 <- ggplot(scat_df,
               aes(x = stripe_logFC, y = log2FoldChange, colour = direction)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_smooth(method = "lm", se = TRUE, colour = "black",
                linetype = "dashed", linewidth = 0.6) +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "grey50") +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "grey50") +
    scale_colour_manual(values = c(gained = "#D73027", lost = "#4575B4")) +
    labs(title = sprintf("Stripe logFC vs body-gene RNA log2FC (%s)", TP_LABEL),
         subtitle = rho_txt,
         x = "Stripe logFC (mut/ctrl)",
         y = "Body-gene RNA log2(mut/ctrl)",
         colour = "Stripe direction") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5))
  save_multiformat_ggplot(p2,
                          file.path(OUTPUT_DIR, "stripe_rna_concordance"),
                          width = 8, height = 6, use_subfolders = TRUE)
}

# -------- SUMMARY --------
summary_lines <- c(
  sprintf("Stripe x RNA-seq integration (%s, res=%d)", TP_LABEL, RESOLUTION),
  sprintf("Generated: %s", Sys.time()),
  sprintf("Stripes:       %d", nrow(stripes)),
  sprintf("RNA-seq genes: %d", nrow(rna)),
  sprintf("Joined pairs:  %d (anchor=%d, body=%d)",
          nrow(joined),
          sum(joined$gene_pool == "anchor"),
          sum(joined$gene_pool == "body")),
  sprintf("DEG pairs (padj<%.2f, |lfc|>%.2f): %d",
          DEG_PADJ_THRESHOLD, DEG_LFC_THRESHOLD, sum(joined$is_deg)),
  ""
)
if (!is.null(wx)) {
  summary_lines <- c(summary_lines, "Wilcoxon pairwise (BH-adjusted):", "",
                     capture.output(print(wx, row.names = FALSE)))
}
if (!is.null(sp)) {
  summary_lines <- c(summary_lines, "", "Spearman concordance:", "",
                     capture.output(print(sp, row.names = FALSE)))
}
writeLines(summary_lines, file.path(OUTPUT_DIR, "summary.txt"))
cat("\nDone.\n")
