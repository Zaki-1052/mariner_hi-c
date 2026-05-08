# loops/scripts/utils/se_utils.R
#
# Shared utility functions for the Superenhancer Hub Analysis pipeline.
# Sourced by se_deg_proximity.R, se_stripe_overlap.R, se_hub_classification.R.
#
# Provides: SE BED loading, DEG loading with invariant gene matching,
#           gene coordinate lookup, bedtools overlap, K27ac classification.

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(GenomicRanges)
  library(IRanges)
  library(GenomeInfoDb)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(org.Mm.eg.db)
})

BEDTOOLS_PATH <- "/opt/homebrew/bin/bedtools"
VALID_CHROMS <- paste0("chr", c(1:19, "X", "Y"))

# =============================================================================
# SE BED Loading
# =============================================================================

load_se_bed <- function(se_file, label) {
  if (!file.exists(se_file)) stop(sprintf("SE file not found: %s", se_file))

  bed <- read_tsv(se_file, col_names = c("chr", "start", "end", "se_name"),
                  col_types = "ciic", show_col_types = FALSE)
  bed <- bed %>% dplyr::filter(chr %in% VALID_CHROMS)

  gr <- GRanges(
    seqnames = bed$chr,
    ranges = IRanges(start = bed$start, end = bed$end),
    se_name = bed$se_name,
    se_source = label
  )
  cat(sprintf("  Loaded %d SEs from %s (%s)\n", length(gr), label, basename(se_file)))
  return(gr)
}

load_all_ses <- function(p60_file, encode_file) {
  p60_gr <- load_se_bed(p60_file, "P60")
  encode_gr <- load_se_bed(encode_file, "ENCODE")
  combined <- c(p60_gr, encode_gr)
  cat(sprintf("  Combined SE set: %d total (%d P60 + %d ENCODE)\n",
              length(combined), length(p60_gr), length(encode_gr)))
  return(combined)
}

# =============================================================================
# RNA-seq / DEG Loading
# =============================================================================

load_rna_for_se <- function(rna_file, padj_thresh = 0.05, lfc_thresh = 0.3) {
  if (!file.exists(rna_file)) stop(sprintf("RNA-seq file not found: %s", rna_file))

  rna_df <- read_excel(rna_file)

  if (!"ensembl_gene_id" %in% colnames(rna_df)) {
    stop("Expected column 'ensembl_gene_id' not found in RNA-seq file")
  }

  all_genes <- rna_df %>%
    dplyr::filter(!is.na(ensembl_gene_id), ensembl_gene_id != "") %>%
    dplyr::select(gene_symbol = ensembl_gene_id, log2FoldChange, padj, baseMean) %>%
    dplyr::distinct(gene_symbol, .keep_all = TRUE)

  degs <- all_genes %>%
    dplyr::filter(!is.na(padj) & padj < padj_thresh) %>%
    dplyr::filter(abs(log2FoldChange) > lfc_thresh) %>%
    dplyr::mutate(
      deg_direction = ifelse(log2FoldChange > 0, "DEG_up", "DEG_down")
    )

  cat(sprintf("  Loaded RNA-seq: %d total genes, %d DEGs (padj < %g, |lfc| > %g)\n",
              nrow(all_genes), nrow(degs), padj_thresh, lfc_thresh))
  cat(sprintf("    - DEG_up: %d, DEG_down: %d\n",
              sum(degs$deg_direction == "DEG_up"),
              sum(degs$deg_direction == "DEG_down")))

  return(list(degs = degs, all_genes = all_genes))
}

# =============================================================================
# Gene Coordinate Lookup
# =============================================================================

get_gene_coordinates <- function(gene_symbols) {
  symbol_to_entrez <- AnnotationDbi::select(
    org.Mm.eg.db,
    keys = unique(gene_symbols),
    columns = "ENTREZID",
    keytype = "SYMBOL"
  ) %>%
    dplyr::filter(!is.na(ENTREZID)) %>%
    dplyr::distinct(SYMBOL, .keep_all = TRUE)

  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  genes_gr <- genes(txdb)

  matched_entrez <- intersect(symbol_to_entrez$ENTREZID, names(genes_gr))
  matched_gr <- genes_gr[matched_entrez]

  tss_pos <- ifelse(as.character(strand(matched_gr)) == "-",
                    end(matched_gr), start(matched_gr))

  gene_coords <- tibble(
    entrez_id = names(matched_gr),
    chr = as.character(seqnames(matched_gr)),
    tss = tss_pos,
    strand = as.character(strand(matched_gr))
  ) %>%
    dplyr::inner_join(symbol_to_entrez, by = c("entrez_id" = "ENTREZID")) %>%
    dplyr::rename(gene_symbol = SYMBOL) %>%
    dplyr::filter(chr %in% VALID_CHROMS)

  n_input <- length(unique(gene_symbols))
  n_mapped <- nrow(gene_coords)
  pct_mapped <- round(100 * n_mapped / n_input, 1)
  cat(sprintf("  Mapped %d / %d gene symbols to coordinates (%s%%)\n",
              n_mapped, n_input, pct_mapped))

  if (pct_mapped < 80) {
    stop(sprintf("Only %s%% of gene symbols mapped — check gene symbol format", pct_mapped))
  }

  return(gene_coords)
}

# =============================================================================
# Invariant Gene Selection
# =============================================================================

select_invariant_genes <- function(degs, all_genes, gene_coords,
                                   n_per_deg = 3, seed = 42) {
  deg_symbols <- unique(degs$gene_symbol)

  all_with_coords <- all_genes %>%
    dplyr::inner_join(gene_coords %>% dplyr::select(gene_symbol, chr),
                      by = "gene_symbol")

  invariant_pool <- all_with_coords %>%
    dplyr::filter(
      !gene_symbol %in% deg_symbols,
      !is.na(padj),
      padj > 0.5,
      abs(log2FoldChange) < 0.1
    )

  degs_with_coords <- degs %>%
    dplyr::inner_join(gene_coords %>% dplyr::select(gene_symbol, chr),
                      by = "gene_symbol")

  breaks <- quantile(all_with_coords$baseMean, probs = seq(0, 1, 0.2), na.rm = TRUE)
  breaks[1] <- 0
  degs_with_coords <- degs_with_coords %>%
    dplyr::mutate(expr_bin = cut(baseMean, breaks = breaks, include.lowest = TRUE, labels = FALSE))
  invariant_pool <- invariant_pool %>%
    dplyr::mutate(expr_bin = cut(baseMean, breaks = breaks, include.lowest = TRUE, labels = FALSE))

  set.seed(seed)
  matched_invariants <- list()

  for (i in seq_len(nrow(degs_with_coords))) {
    deg_row <- degs_with_coords[i, ]
    candidates <- invariant_pool %>%
      dplyr::filter(chr == deg_row$chr, expr_bin == deg_row$expr_bin) %>%
      dplyr::filter(!gene_symbol %in% unlist(matched_invariants))

    if (nrow(candidates) < n_per_deg) {
      candidates <- invariant_pool %>%
        dplyr::filter(expr_bin == deg_row$expr_bin) %>%
        dplyr::filter(!gene_symbol %in% unlist(matched_invariants))
    }

    n_select <- min(n_per_deg, nrow(candidates))
    if (n_select > 0) {
      selected <- candidates %>% dplyr::slice_sample(n = n_select)
      matched_invariants[[i]] <- selected$gene_symbol
    }
  }

  invariant_symbols <- unique(unlist(matched_invariants))
  invariant_df <- all_with_coords %>%
    dplyr::filter(gene_symbol %in% invariant_symbols) %>%
    dplyr::mutate(deg_direction = "invariant")

  cat(sprintf("  Selected %d invariant genes matched to %d DEGs (%.1f:1 ratio)\n",
              nrow(invariant_df), nrow(degs_with_coords),
              nrow(invariant_df) / max(1, nrow(degs_with_coords))))

  return(invariant_df)
}

# =============================================================================
# Bedtools Overlap
# =============================================================================

bedtools_intersect_gr <- function(query_gr, subject_gr,
                                  query_name = "query", subject_name = "subject",
                                  min_overlap = 1, extra_args = "") {
  if (!file.exists(BEDTOOLS_PATH)) {
    stop(sprintf("bedtools not found at %s", BEDTOOLS_PATH))
  }

  query_bed <- tempfile(fileext = ".bed")
  subject_bed <- tempfile(fileext = ".bed")
  output_file <- tempfile(fileext = ".tsv")
  on.exit(file.remove(c(query_bed, subject_bed, output_file)), add = TRUE)

  query_df <- tibble(
    chr = as.character(seqnames(query_gr)),
    start = start(query_gr),
    end = end(query_gr),
    name = seq_along(query_gr)
  )
  subject_df <- tibble(
    chr = as.character(seqnames(subject_gr)),
    start = start(subject_gr),
    end = end(subject_gr),
    name = seq_along(subject_gr)
  )

  write_tsv(query_df, query_bed, col_names = FALSE)
  write_tsv(subject_df, subject_bed, col_names = FALSE)

  cmd <- sprintf("%s intersect -wa -wb %s -a %s -b %s > %s",
                 BEDTOOLS_PATH, extra_args, query_bed, subject_bed, output_file)
  exit_code <- system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)

  if (exit_code != 0) stop(sprintf("bedtools intersect failed (exit %d)", exit_code))

  result_lines <- readLines(output_file)
  if (length(result_lines) == 0) {
    cat(sprintf("  No overlaps found between %s and %s\n", query_name, subject_name))
    return(tibble(query_idx = integer(), subject_idx = integer(),
                  query_chr = character(), query_start = integer(), query_end = integer(),
                  subject_chr = character(), subject_start = integer(), subject_end = integer()))
  }

  result <- read_tsv(output_file, col_names = FALSE, show_col_types = FALSE,
                     col_types = "ciicciic")
  colnames(result) <- c(
    paste0(query_name, "_chr"), paste0(query_name, "_start"),
    paste0(query_name, "_end"), "query_idx",
    paste0(subject_name, "_chr"), paste0(subject_name, "_start"),
    paste0(subject_name, "_end"), "subject_idx"
  )
  result <- result %>%
    dplyr::mutate(query_idx = as.integer(query_idx),
                  subject_idx = as.integer(subject_idx))

  cat(sprintf("  Found %d overlaps between %s (%d) and %s (%d)\n",
              nrow(result), query_name, length(query_gr),
              subject_name, length(subject_gr)))

  return(result)
}

# =============================================================================
# K27ac DiffBind Loading and Classification
# =============================================================================

load_k27ac_diffbind <- function(diffbind_file) {
  if (!file.exists(diffbind_file)) {
    stop(sprintf("DiffBind file not found: %s", diffbind_file))
  }

  db <- read_tsv(diffbind_file, show_col_types = FALSE)

  required_cols <- c("Summit_Chr", "Summit_Start", "Summit_End", "Fold", "FDR")
  missing <- setdiff(required_cols, colnames(db))
  if (length(missing) > 0) {
    stop(sprintf("Missing DiffBind columns: %s", paste(missing, collapse = ", ")))
  }

  db <- db %>% dplyr::filter(Summit_Chr %in% VALID_CHROMS)

  gr <- GRanges(
    seqnames = db$Summit_Chr,
    ranges = IRanges(start = db$Summit_Start, end = db$Summit_End),
    fold = db$Fold,
    fdr = db$FDR
  )

  cat(sprintf("  Loaded %d K27ac DiffBind peaks\n", length(gr)))
  return(gr)
}

classify_k27ac_change <- function(peaks_gr, diffbind_gr,
                                  lfc_thresh = 0.5, fdr_thresh = 0.05) {
  hits <- findOverlaps(peaks_gr, diffbind_gr, ignore.strand = TRUE)

  if (length(hits) == 0) {
    return(tibble(
      peak_idx = seq_along(peaks_gr),
      k27ac_class = rep("no_diffbind_peak", length(peaks_gr))
    ))
  }

  hit_df <- tibble(
    peak_idx = queryHits(hits),
    fold = diffbind_gr$fold[subjectHits(hits)],
    fdr = diffbind_gr$fdr[subjectHits(hits)]
  )

  best_hit <- hit_df %>%
    dplyr::group_by(peak_idx) %>%
    dplyr::slice_min(fdr, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      k27ac_class = dplyr::case_when(
        fdr < fdr_thresh & fold > lfc_thresh  ~ "gained_k27ac",
        fdr < fdr_thresh & fold < -lfc_thresh ~ "lost_k27ac",
        TRUE                                   ~ "stable_k27ac"
      )
    )

  result <- tibble(peak_idx = seq_along(peaks_gr)) %>%
    dplyr::left_join(best_hit %>% dplyr::select(peak_idx, k27ac_class),
                     by = "peak_idx") %>%
    dplyr::mutate(k27ac_class = ifelse(is.na(k27ac_class), "no_diffbind_peak", k27ac_class))

  cat(sprintf("  K27ac classification: %s\n",
              paste(names(table(result$k27ac_class)),
                    table(result$k27ac_class), sep = "=", collapse = ", ")))

  return(result)
}

# =============================================================================
# Enhancer BED Loading
# =============================================================================

load_enhancer_bed <- function(bed_file, label, n_cols = NULL) {
  if (!file.exists(bed_file)) stop(sprintf("Enhancer file not found: %s", bed_file))

  first_line <- readLines(bed_file, n = 1)
  detected_cols <- length(strsplit(first_line, "\t")[[1]])

  if (detected_cols >= 6) {
    bed <- read_tsv(bed_file, col_names = FALSE, show_col_types = FALSE) %>%
      dplyr::select(chr = X1, start = X2, end = X3)
  } else {
    bed <- read_tsv(bed_file,
                    col_names = c("chr", "start", "end")[1:min(3, detected_cols)],
                    show_col_types = FALSE)
  }

  bed <- bed %>% dplyr::filter(chr %in% VALID_CHROMS)

  gr <- GRanges(
    seqnames = bed$chr,
    ranges = IRanges(start = bed$start, end = bed$end)
  )

  cat(sprintf("  Loaded %d enhancer peaks from %s\n", length(gr), label))
  return(gr)
}

# =============================================================================
# Shared Color Schemes
# =============================================================================

DIRECTION_COLORS <- c("lost" = "#4575b4", "gained" = "#e0a730")
DEG_CLASS_COLORS <- c("DEG_up" = "#d73027", "DEG_down" = "#4575b4", "invariant" = "#999999")
K27AC_CLASS_COLORS <- c(
  "gained_k27ac" = "#e0a730",
  "lost_k27ac" = "#4575b4",
  "stable_k27ac" = "#999999",
  "no_diffbind_peak" = "#cccccc"
)
