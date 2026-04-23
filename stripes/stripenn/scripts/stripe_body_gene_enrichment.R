#!/usr/bin/env Rscript
# stripes/stripenn/scripts/stripe_body_gene_enrichment.R
#
# T1.1 - GO/KEGG enrichment on stripe BODY genes (genes resident inside
# the stripe span pos3-pos4), as a complement to the Stage 7 anchor-nearest-gene
# enrichment in stripes/stripenn/scripts/stripenn_visualizations.R Section 10.
#
# Anchor-gene enrichment captures CTCF-neighborhood regulators. Body-gene
# enrichment captures loci under the stripe-contact "sweep" - a mechanistically
# distinct gene pool.
#
# Usage:
#   Rscript scripts/stripe_body_gene_enrichment.R --timepoint late --resolution 5000
#
# Inputs:
#   stripes/stripenn/outputs/{tp}/visualizations/{tp}_annotated_stripes.tsv
#     - body_genes (comma-separated symbols), direction, direction_confidence,
#       significant_FDR05
#
# Outputs (under outputs/stripe_integration/{label}/body_gene_enrichment/):
#   - stripe_body_genes.tsv          Per-gene-per-category table
#   - go_{bp,cc,mf}.tsv / kegg.tsv   Enrichment results
#   - go_{bp,cc,mf}_dotplot/         Multi-format dotplots (pdf, svg, jpg)
#   - kegg_dotplot/
#   - summary.txt                    Run summary

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(clusterProfiler)
  library(enrichplot)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
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

# -------- CLI --------
args <- commandArgs(trailingOnly = TRUE)
parse_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  if (i == length(args)) stop(sprintf("Missing value for %s", flag))
  args[i + 1]
}
TP_LABEL   <- parse_arg("--timepoint", "late")
RESOLUTION <- as.integer(parse_arg("--resolution", "5000"))

TP_MAP <- list(late = "250402", early = "250831")
if (!TP_LABEL %in% names(TP_MAP)) stop("--timepoint must be 'late' or 'early'")
TP_ID <- TP_MAP[[TP_LABEL]]

ANNOT_FILE <- file.path(CODE_DIR, sprintf("outputs/%s/visualizations/%s_annotated_stripes.tsv",
                     TP_ID, TP_ID))
OUTPUT_DIR <- file.path(REPO_ROOT, sprintf("outputs/stripe_integration/%s/body_gene_enrichment", TP_LABEL))
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("==============================================================\n")
cat("T1.1 - Stripe Body-Gene GO/KEGG Enrichment\n")
cat("==============================================================\n")
cat(sprintf("Timepoint : %s (%s)\n", TP_LABEL, TP_ID))
cat(sprintf("Resolution: %d bp\n", RESOLUTION))
cat(sprintf("Input     : %s\n", ANNOT_FILE))
cat(sprintf("Output    : %s\n\n", OUTPUT_DIR))

stopifnot(file.exists(ANNOT_FILE))

# -------- LOAD --------
ann <- read.table(ANNOT_FILE, header = TRUE, sep = "\t",
                  stringsAsFactors = FALSE, quote = "", comment.char = "",
                  na.strings = c("NA", ""))
cat(sprintf("Loaded %d stripes with %d columns\n", nrow(ann), ncol(ann)))

required <- c("body_genes", "direction", "direction_confidence", "significant_FDR05")
missing <- setdiff(required, colnames(ann))
if (length(missing) > 0) stop("Missing columns: ", paste(missing, collapse = ", "))

# body_genes is comma-separated symbols (per stripenn_visualizations.R line 725:
# paste(syms, collapse = ","))
explode_body <- function(df) {
  df <- df[!is.na(df$body_genes) & df$body_genes != "", , drop = FALSE]
  if (nrow(df) == 0) return(data.frame(stripe_id = character(), symbol = character()))
  out <- df %>%
    dplyr::select(stripe_id, direction, direction_confidence, significant_FDR05, body_genes) %>%
    mutate(symbol = strsplit(body_genes, ",", fixed = TRUE)) %>%
    tidyr::unnest(symbol) %>%
    mutate(symbol = trimws(symbol)) %>%
    filter(symbol != "")
  out
}

body_long <- explode_body(ann)
cat(sprintf("Exploded to %d (stripe, body-gene) pairs (%d unique symbols)\n",
            nrow(body_long), length(unique(body_long$symbol))))

# -------- SYMBOL -> ENTREZ --------
uniq_syms <- unique(body_long$symbol)
sym_map <- tryCatch(
  AnnotationDbi::select(org.Mm.eg.db, keys = uniq_syms,
                        columns = "ENTREZID", keytype = "SYMBOL"),
  error = function(e) {
    cat("SYMBOL lookup failed: ", e$message, "\n")
    data.frame(SYMBOL = character(), ENTREZID = character())
  })
# first non-NA ENTREZID per symbol
sym_map <- sym_map %>%
  filter(!is.na(ENTREZID)) %>%
  distinct(SYMBOL, .keep_all = TRUE)
cat(sprintf("Mapped %d/%d symbols to Entrez IDs\n",
            nrow(sym_map), length(uniq_syms)))

body_long$entrez <- sym_map$ENTREZID[match(body_long$symbol, sym_map$SYMBOL)]
body_long <- body_long[!is.na(body_long$entrez), , drop = FALSE]

# -------- CATEGORY CLUSTERS --------
# Match the "high-confidence" convention used by generate_bedpe.py:
# direction_confidence == "high" AND significant_FDR05 == TRUE for direction in (lost, gained).
# For strengthened/weakened just require significance (populations are tiny).
high_conf_gained <- body_long %>%
  filter(direction == "gained", direction_confidence == "high", significant_FDR05 == TRUE) %>%
  distinct(entrez) %>% pull(entrez)
high_conf_lost <- body_long %>%
  filter(direction == "lost", direction_confidence == "high", significant_FDR05 == TRUE) %>%
  distinct(entrez) %>% pull(entrez)
strengthened <- body_long %>%
  filter(direction == "strengthened", significant_FDR05 == TRUE) %>%
  distinct(entrez) %>% pull(entrez)
weakened <- body_long %>%
  filter(direction == "weakened", significant_FDR05 == TRUE) %>%
  distinct(entrez) %>% pull(entrez)

# Universe = all body genes across the whole union set (timepoint-matched)
universe <- unique(body_long$entrez)

gene_clusters <- list(
  Gained_high = high_conf_gained,
  Lost_high   = high_conf_lost,
  Strengthened = strengthened,
  Weakened    = weakened
)
# Drop small clusters (n < 10) per the plan
gene_clusters <- gene_clusters[sapply(gene_clusters, length) >= 10]
cat("\nGene clusters (Entrez counts):\n")
for (nm in names(gene_clusters)) {
  cat(sprintf("  %-14s %5d\n", nm, length(gene_clusters[[nm]])))
}
cat(sprintf("  %-14s %5d\n", "UNIVERSE", length(universe)))

if (length(gene_clusters) == 0) {
  cat("\nNo cluster passed n>=10; exiting without enrichment.\n")
  quit(save = "no", status = 0)
}

# -------- PER-GENE-PER-CATEGORY TABLE --------
body_export <- body_long %>%
  transmute(
    stripe_id, direction, direction_confidence, significant_FDR05,
    symbol, entrez_id = entrez
  )
write.table(body_export,
            file.path(OUTPUT_DIR, "stripe_body_genes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# -------- GO ENRICHMENT --------
go_config <- list(
  BP = list(label = "Biological Process", show = 20, w = 12, h = 10),
  CC = list(label = "Cellular Component", show = 15, w = 10, h = 8),
  MF = list(label = "Molecular Function", show = 15, w = 10, h = 8)
)

summary_lines <- c(
  sprintf("Stripe body-gene enrichment summary (%s, res=%d)", TP_LABEL, RESOLUTION),
  sprintf("Generated: %s", Sys.time()),
  sprintf("Total stripes: %d", nrow(ann)),
  sprintf("Universe (body-gene Entrez IDs): %d", length(universe)),
  ""
)
for (nm in names(gene_clusters)) {
  summary_lines <- c(summary_lines, sprintf("  %s cluster: n = %d", nm, length(gene_clusters[[nm]])))
}
summary_lines <- c(summary_lines, "")

for (ont in names(go_config)) {
  cfg <- go_config[[ont]]
  cat(sprintf("\nGO %s (%s)...\n", ont, cfg$label))
  go_res <- tryCatch(
    compareCluster(geneCluster = gene_clusters,
                   fun = "enrichGO", OrgDb = org.Mm.eg.db, ont = ont,
                   keyType = "ENTREZID", universe = universe,
                   pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05),
    error = function(e) { cat(sprintf("  compareCluster error: %s\n", e$message)); NULL })

  if (!is.null(go_res) && !is.null(go_res@compareClusterResult) &&
      nrow(go_res@compareClusterResult) > 0) {
    write.table(go_res@compareClusterResult,
                file.path(OUTPUT_DIR, sprintf("go_%s.tsv", tolower(ont))),
                sep = "\t", quote = FALSE, row.names = FALSE)
    p <- dotplot(go_res, showCategory = cfg$show) +
      labs(title = sprintf("GO %s: body-gene enrichment (%s)", cfg$label, TP_LABEL)) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    save_multiformat_ggplot(p,
                            file.path(OUTPUT_DIR, sprintf("go_%s_dotplot", tolower(ont))),
                            width = cfg$w, height = cfg$h, use_subfolders = TRUE)
    msg <- sprintf("GO %s: %d significant terms", ont, nrow(go_res@compareClusterResult))
    cat("  ", msg, "\n")
    summary_lines <- c(summary_lines, msg)
  } else {
    msg <- sprintf("GO %s: no significant terms", ont)
    cat("  ", msg, "\n")
    summary_lines <- c(summary_lines, msg)
  }
}

# -------- KEGG ENRICHMENT --------
cat("\nKEGG pathways...\n")
kegg_res <- tryCatch(
  compareCluster(geneCluster = gene_clusters, fun = "enrichKEGG",
                 organism = "mmu", universe = universe,
                 pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05),
  error = function(e) { cat(sprintf("  compareCluster error: %s\n", e$message)); NULL })

if (!is.null(kegg_res) && !is.null(kegg_res@compareClusterResult) &&
    nrow(kegg_res@compareClusterResult) > 0) {
  write.table(kegg_res@compareClusterResult,
              file.path(OUTPUT_DIR, "kegg.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  p <- dotplot(kegg_res, showCategory = 20) +
    labs(title = sprintf("KEGG: body-gene enrichment (%s)", TP_LABEL)) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  save_multiformat_ggplot(p,
                          file.path(OUTPUT_DIR, "kegg_dotplot"),
                          width = 12, height = 10, use_subfolders = TRUE)
  msg <- sprintf("KEGG: %d significant pathways", nrow(kegg_res@compareClusterResult))
  cat("  ", msg, "\n")
  summary_lines <- c(summary_lines, msg)
} else {
  msg <- "KEGG: no significant pathways"
  cat("  ", msg, "\n")
  summary_lines <- c(summary_lines, msg)
}

writeLines(summary_lines, file.path(OUTPUT_DIR, "summary.txt"))
cat("\nDone.\n")
