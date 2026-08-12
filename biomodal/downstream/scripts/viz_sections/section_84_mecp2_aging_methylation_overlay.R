# scripts/viz_sections/section_84_mecp2_aging_methylation_overlay.R
# Section 84: MeCP2 Aging x Methylation Overlap (Fisher test)
#
# Do the 1,654 mutant-unique aging MeCP2 genes co-localize with
# coordinated mC-up / hmC-down gene bodies?
#
# Re-derives mutant-unique aging-UP genes from MeCP2 aging DiffBind files,
# then tests Fisher enrichment against the coordinated DMR gene set over
# the 20,915-gene universe.
#
# Output: plots/visualizations/84_mecp2_aging_methylation_overlay/
# Table:  plots/visualizations/tables/84_aging_methylation_overlap.tsv
# =============================================================================

source("scripts/viz_sections/_shared_config.R")
source("scripts/viz_sections/_figure_config.R")

suppressPackageStartupMessages({
  library(ggplot2)
})

SEC84_DIR <- file.path(OUTPUT_DIR, "84_mecp2_aging_methylation_overlay")
dir.create(SEC84_DIR, recursive = TRUE, showWarnings = FALSE)

save_plot <- function(p, name, w = 10, h = 7) {
  save_multiformat_ggplot(p,
    base_path = file.path(SEC84_DIR, name),
    width = w, height = h, dpi = 300,
    verbose = TRUE, use_subfolders = TRUE
  )
}

cat("\n=== Section 84: MeCP2 Aging x Methylation Overlap ===\n\n")

# Load an MeCP2 aging DiffBind file and return age-increased gene symbols.
# Replicates section_77's loader: negates Fold to log2(adult/young), filters to
# FDR < threshold & Fold > 0, annotates each peak to its nearest gene.
aging_up_genes <- function(filepath, label) {
  if (!file.exists(filepath)) {
    stop("aging file not found: ", filepath)
  }
  df <- read.table(filepath, header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE, quote = "", fill = TRUE)
  if (!all(c("Fold", "FDR") %in% colnames(df))) {
    stop("aging file missing Fold/FDR columns: ", filepath)
  }
  df$Fold <- -df$Fold
  up <- df[df$FDR < Q_THRESHOLD & df$Fold > 0, , drop = FALSE]
  cat(sprintf("  %s: %d age-increased peaks (FDR < %.2f, adult > young)\n",
              label, nrow(up), Q_THRESHOLD))

  chr_col   <- if ("seqnames" %in% colnames(up)) "seqnames" else "chr"
  start_col <- if ("start" %in% colnames(up)) "start" else "Start"
  end_col   <- if ("end" %in% colnames(up)) "end" else "End"

  gr <- GenomicRanges::GRanges(
    seqnames = up[[chr_col]],
    ranges   = IRanges::IRanges(start = up[[start_col]], end = up[[end_col]])
  )
  anno <- ChIPseeker::annotatePeak(
    gr, TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene,
    annoDb = "org.Mm.eg.db", level = "gene", verbose = FALSE
  )
  syms <- as.data.frame(anno)$SYMBOL
  genes <- unique(syms[!is.na(syms) & syms != ""])
  cat(sprintf("  %s: %d unique age-increased genes\n", label, length(genes)))
  genes
}

panel_84a <- tryCatch({
  ctrl_aging_genes <- aging_up_genes(MECP2_FILES$ctrl_aging, "Control aging")
  mut_aging_genes  <- aging_up_genes(MECP2_FILES$mut_aging,  "Mutant aging")

  mut_unique_aging <- setdiff(mut_aging_genes, ctrl_aging_genes)
  cat(sprintf("  Mutant-unique aging-UP genes = %d\n", length(mut_unique_aging)))

  coord_path <- file.path(TABLES_DIR, "mecp2_coordinated_genes.tsv")
  stopifnot("mecp2_coordinated_genes.tsv not found" = file.exists(coord_path))
  coord_tbl <- read.table(coord_path, header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE, quote = "", check.names = FALSE,
                          comment.char = "")
  stopifnot("gene" %in% colnames(coord_tbl))
  coordinated_genes <- unique(coord_tbl$gene[!is.na(coord_tbl$gene) & coord_tbl$gene != ""])

  univ_path <- file.path(TABLES_DIR, "demethylation_ratio_all_genes.tsv")
  stopifnot("demethylation_ratio_all_genes.tsv not found" = file.exists(univ_path))
  universe_tbl <- read.table(univ_path, header = TRUE, sep = "\t",
                             stringsAsFactors = FALSE, quote = "", check.names = FALSE,
                             comment.char = "")
  stopifnot("gene" %in% colnames(universe_tbl))
  universe <- unique(universe_tbl$gene[!is.na(universe_tbl$gene) & universe_tbl$gene != ""])
  cat(sprintf("  Universe = %d genes; coordinated = %d\n",
              length(universe), length(coordinated_genes)))

  mut_unique_in_univ <- intersect(mut_unique_aging, universe)
  coord_in_univ      <- intersect(coordinated_genes, universe)

  is_aging <- universe %in% mut_unique_in_univ
  is_coord <- universe %in% coord_in_univ

  a <- sum(is_aging  &  is_coord)
  b <- sum(is_aging  & !is_coord)
  c <- sum(!is_aging &  is_coord)
  d <- sum(!is_aging & !is_coord)

  fmat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  ft <- fisher.test(fmat)
  or_val  <- unname(ft$estimate)
  p_val   <- ft$p.value
  ci_val  <- ft$conf.int

  n_univ   <- length(universe)
  expected <- (sum(is_aging) * sum(is_coord)) / n_univ
  observed <- a
  cat(sprintf("  Fisher 2x2: a=%d b=%d c=%d d=%d; OR=%.3f [%.3f, %.3f] p=%.3e\n",
              a, b, c, d, or_val, ci_val[1], ci_val[2], p_val))
  cat(sprintf("  Observed overlap = %d vs expected = %.1f\n", observed, expected))

  overlap_out <- data.frame(
    test = "mut_unique_aging_UP x coordinated_mCup_hmCdown",
    universe_n = n_univ,
    mut_unique_aging_n = sum(is_aging),
    coordinated_n = sum(is_coord),
    overlap_observed = observed,
    overlap_expected = round(expected, 2),
    a = a, b = b, c = c, d = d,
    odds_ratio = or_val,
    ci_lower = ci_val[1],
    ci_upper = ci_val[2],
    p_value = p_val,
    stringsAsFactors = FALSE
  )
  out_path <- file.path(TABLES_DIR, "84_aging_methylation_overlap.tsv")
  write.table(overlap_out, out_path, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("  Saved:", out_path, "\n")

  obs_exp_df <- data.frame(
    kind  = factor(c("Observed", "Expected"), levels = c("Expected", "Observed")),
    count = c(observed, expected)
  )

  ggplot(obs_exp_df, aes(x = kind, y = count, fill = kind)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = ifelse(kind == "Observed",
                                 as.character(round(count)),
                                 sprintf("%.1f", count))),
              vjust = -0.3, size = 2.8) +
    scale_fill_manual(values = c("Expected" = "grey70",
                                 "Observed" = COLORS$condition[["Mutant"]]),
                      guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
    annotate("text", x = 1.5, y = max(obs_exp_df$count) * 1.12,
             label = sprintf("Fisher OR = %.2f\np = %.1e", or_val, p_val),
             size = 2.8, lineheight = 0.95) +
    labs(
      title = "Mut-unique aging genes are coordinated DMRs",
      subtitle = sprintf("%s mut-unique aging genes x %s coordinated genes",
                         format(sum(is_aging), big.mark = ","),
                         format(sum(is_coord), big.mark = ",")),
      x = NULL, y = "Coordinated-DMR genes"
    ) +
    theme_pub()
}, error = function(e) {
  warning("Section 84 skipped: ", conditionMessage(e))
  ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = paste0("Section 84 unavailable:\n", conditionMessage(e)),
             size = 2.8, lineheight = 0.95) +
    theme_void() +
    labs(title = "Aging x methylation overlap")
})

save_plot(panel_84a, "84a_aging_methylation_overlap", w = 7, h = 6)
save_plot(panel_84a, "84_composite", w = 7, h = 6)

cat("\nSection 84 complete.\n")
