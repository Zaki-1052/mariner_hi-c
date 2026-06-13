#!/usr/bin/env Rscript
# loops/scripts/build_master_loop_table.R
#
# Build the MASTER ANNOTATED LOOP LIST (Sheet 1) for the paper.
#
# Spine: all FDR<0.05 significant differential loops from the merged-all result set
# (NOT the |logFC|>0.3-filtered "characterized" set). logFC is retained so the
# stricter cut can be applied downstream on demand.
#
# For each significant loop this assembles:
#   - the differential statistics carried from merged_all_results.tsv (incl. logFC)
#   - span_length        : interior gap between anchors  = start2 - end1
#   - loop_footprint      : full extent end2 - start1 (used only for the clust6 split)
#   - cluster / cluster_label : k-means cluster (clust6 -> clust6_short / clust6_long @ 800kb)
#   - loop_class          : structural / CRE / mixed / unclassified (CTCF + E/P rule)
#   - SE_anchor           : TRUE if an anchor overlaps a P60 super-enhancer
#   - shared_anchor       : TRUE if an anchor sits at a loop-switching locus (re-derived on
#                           this significant set: a locus with >=1 lost AND >=1 gained anchor)
#   - DEG_status / DEG_genes : TRUE if an anchor overlaps the gene body of a padj<0.05 DEG
#   - "previous information": per-anchor ChIP overlaps + 8-state chromatin type + loop_type
#                           + nearest gene (re-derived for ALL significant loops)
#
# Design: runs from the REPO ROOT with explicit paths (loops/peaks does not exist, so the
# CWD-relative paths in annotate_loops_extended.R are not reused; its two PURE classification
# functions are reproduced verbatim below). Everything tunable lives in CONFIG / CLI flags.
#
# Usage (from repo root, with the loops/ Bioconductor R, e.g. mariner_env):
#   Rscript loops/scripts/build_master_loop_table.R
#   Rscript loops/scripts/build_master_loop_table.R --fdr 0.05 --deg-padj 0.05 --deg-lfc 0
#
# Outputs (CONFIG$out_dir):
#   significant_loops.tsv        7,981 spine rows (ids + stats) — also the Sheet 2 input
#   master_annotated_loops.tsv   the full master table

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IRanges)
  library(rtracklayer)
  library(readxl)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
})

# =============================================================================
# CONFIG — every path / threshold is here so the table is easy to re-target
# =============================================================================

CONFIG <- list(
  spine_tsv     = "loops/outputs/250402-late_outputs/merged_loops/merged_all_results.tsv",
  clusters_txt  = "cluster/outputs/bap1_late/cluster3/k-6/data/combined-clusters.txt",
  peak_files = list(
    h3k27ac    = "peaks/beds/H3K27acCerebellumLate2.bed",
    h3k27me3   = "peaks/beds/H3K27me3CerebellumLate1.bed",
    h3k4me1    = "peaks/beds/H3K4me1CerebellumLate1.bed",
    h3k4me3    = "peaks/beds/H3K4me3CerebellumLate2.bed",
    bivalent   = "peaks/beds/Bivalent_Cerebellum_Late.bed",
    ctcf       = "peaks/CTCF.bed",
    ctcf_motif = "peaks/ctcf_motifs_mm10.bed"
  ),
  promoter_bed  = "cluster/data/mm10_knownGene_pp.bed",                 # TSS+-750bp, for loop_class E/P
  se_bed        = "peaks/Superenhancers_P60.bed",                       # adult/P60 super-enhancers
  deg_xlsx      = "tads/adult_timepoint_rna-seq-BAP1_WT_KO_v2_Results.xlsx",
  gencode_genes = "biomodal/downstream/modality/mm10/gencode.vM25.mouse.genes.annotation.bed.gz",
  out_dir       = "loops/output/master_loop_table/late",

  fdr_threshold    = 0.05,     # spine = significant (FDR < this)
  deg_padj         = 0.05,     # DEG significance
  deg_lfc          = 0,        # |log2FC| cut for DEGs (0 = padj only, per PI)
  clust6_threshold = 800000,   # clust6_short < thr <= clust6_long  (on loop_footprint = end2-start1)
  shared_maxgap    = 10000,    # loop-switching locus tolerance
  tss_threshold    = 2000,     # promoter distance for chromatin-state classification
  std_chroms       = paste0("chr", c(1:19, "X"))
)

# =============================================================================
# CLI OVERRIDES  (--key value); unknown flags abort loudly
# =============================================================================

parse_cli <- function(cfg) {
  args <- commandArgs(trailingOnly = TRUE)
  flag_map <- c(
    "--spine"            = "spine_tsv",
    "--clusters"         = "clusters_txt",
    "--se-bed"           = "se_bed",
    "--deg-xlsx"         = "deg_xlsx",
    "--gencode"          = "gencode_genes",
    "--out-dir"          = "out_dir",
    "--fdr"              = "fdr_threshold",
    "--deg-padj"         = "deg_padj",
    "--deg-lfc"          = "deg_lfc",
    "--clust6-threshold" = "clust6_threshold",
    "--shared-maxgap"    = "shared_maxgap",
    "--tss-threshold"    = "tss_threshold"
  )
  numeric_keys <- c("fdr_threshold", "deg_padj", "deg_lfc",
                    "clust6_threshold", "shared_maxgap", "tss_threshold")
  i <- 1
  while (i <= length(args)) {
    key <- flag_map[[args[i]]]
    if (is.null(key) || is.na(key)) stop(sprintf("Unknown argument: %s", args[i]))
    if (i == length(args)) stop(sprintf("Missing value for %s", args[i]))
    val <- args[i + 1]
    cfg[[key]] <- if (key %in% numeric_keys) as.numeric(val) else val
    i <- i + 2
  }
  cfg
}

# =============================================================================
# SMALL HELPERS
# =============================================================================

# Integer-safe stringification so coordinate join keys never hit scientific notation.
fmt_int   <- function(x) sprintf("%.0f", as.numeric(x))
key6      <- function(c1, s1, e1, c2, s2, e2)
  paste(c1, fmt_int(s1), fmt_int(e1), c2, fmt_int(s2), fmt_int(e2), sep = "|")
overlaps_any <- function(query_gr, subject_gr) countOverlaps(query_gr, subject_gr) > 0L

require_file <- function(path, what) {
  if (!file.exists(path)) stop(sprintf("%s not found: %s", what, path))
  invisible(path)
}

# BED -> GRanges using only the first three columns (matches load_chip_peaks()).
load_bed_granges <- function(path, label) {
  require_file(path, label)
  df <- read.table(path, sep = "\t", header = FALSE, quote = "",
                   comment.char = "", stringsAsFactors = FALSE)
  gr <- GRanges(seqnames = df$V1, ranges = IRanges(start = df$V2, end = df$V3))
  cat(sprintf("  %-22s %7d regions  (%s)\n", paste0(label, ":"), length(gr), basename(path)))
  gr
}

# -----------------------------------------------------------------------------
# Chromatin-state classification — reproduced VERBATIM (pure functions, no I/O)
# from loops/scripts/annotate_loops_extended.R so the master uses identical logic
# without that script's CWD-relative source()/peak paths.
# -----------------------------------------------------------------------------

ANCHOR_TYPE_ORDER <- c("Active_Promoter", "Repressed_Promoter", "Bivalent_Promoter",
                       "Polycomb", "Active_Enhancer", "Poised_Enhancer",
                       "CTCF_Site", "Other")

classify_anchor_type_extended <- function(h3k27ac_overlap, h3k27me3_overlap,
                                          h3k4me1_overlap, h3k4me3_overlap,
                                          bivalent_overlap,
                                          ctcf_overlap, ctcf_motif_overlap,
                                          distance_to_tss,
                                          tss_threshold = 2000,
                                          use_motif_for_ctcf = TRUE) {
  n <- length(h3k27ac_overlap)
  anchor_type <- rep("Other", n)

  is_active_promoter <- h3k4me3_overlap & !h3k27me3_overlap &
                        !is.na(distance_to_tss) & distance_to_tss <= tss_threshold
  anchor_type[is_active_promoter] <- "Active_Promoter"

  is_repressed_promoter <- !is_active_promoter &
                           h3k27me3_overlap & !h3k27ac_overlap &
                           !is.na(distance_to_tss) & distance_to_tss <= tss_threshold
  anchor_type[is_repressed_promoter] <- "Repressed_Promoter"

  is_bivalent <- !is_active_promoter & !is_repressed_promoter & bivalent_overlap
  anchor_type[is_bivalent] <- "Bivalent_Promoter"

  is_polycomb <- !is_active_promoter & !is_repressed_promoter & !is_bivalent &
                 h3k27me3_overlap &
                 (is.na(distance_to_tss) | distance_to_tss > tss_threshold)
  anchor_type[is_polycomb] <- "Polycomb"

  is_active_enhancer <- !is_active_promoter & !is_repressed_promoter &
                        !is_bivalent & !is_polycomb &
                        h3k27ac_overlap &
                        (is.na(distance_to_tss) | distance_to_tss > tss_threshold)
  anchor_type[is_active_enhancer] <- "Active_Enhancer"

  is_poised_enhancer <- !is_active_promoter & !is_repressed_promoter &
                        !is_bivalent & !is_polycomb & !is_active_enhancer &
                        h3k4me1_overlap & !h3k27ac_overlap & !h3k27me3_overlap &
                        (is.na(distance_to_tss) | distance_to_tss > tss_threshold)
  anchor_type[is_poised_enhancer] <- "Poised_Enhancer"

  if (use_motif_for_ctcf) {
    is_ctcf_site <- !is_active_promoter & !is_repressed_promoter &
                    !is_bivalent & !is_polycomb & !is_active_enhancer &
                    !is_poised_enhancer & ctcf_motif_overlap
  } else {
    is_ctcf_site <- !is_active_promoter & !is_repressed_promoter &
                    !is_bivalent & !is_polycomb & !is_active_enhancer &
                    !is_poised_enhancer & ctcf_overlap
  }
  anchor_type[is_ctcf_site] <- "CTCF_Site"

  anchor_type
}

classify_loop_type_extended <- function(anchor1_type, anchor2_type) {
  n <- length(anchor1_type)
  loop_type <- rep(NA_character_, n)
  for (i in seq_len(n)) {
    t1 <- anchor1_type[i]; t2 <- anchor2_type[i]
    if (t1 == t2) {
      loop_type[i] <- sprintf("%s-%s", t1, t2)
    } else {
      pos1 <- match(t1, ANCHOR_TYPE_ORDER); pos2 <- match(t2, ANCHOR_TYPE_ORDER)
      loop_type[i] <- if (pos1 < pos2) sprintf("%s-%s", t1, t2) else sprintf("%s-%s", t2, t1)
    }
  }
  loop_type
}

# =============================================================================
# STEP FUNCTIONS
# =============================================================================

# Load merged_all and keep FDR<fdr (== the `significant` flag). Order is preserved
# so Sheet 1 and Sheet 2 align row-for-row.
load_significant_spine <- function(path, fdr) {
  require_file(path, "Spine (merged_all_results.tsv)")
  df <- read.delim(path, header = TRUE, sep = "\t", quote = "",
                   comment.char = "", stringsAsFactors = FALSE)
  need <- c("loop_id", "chr1", "start1", "end1", "chr2", "start2", "end2",
            "logFC", "FDR", "significant", "direction")
  miss <- setdiff(need, colnames(df))
  if (length(miss) > 0) stop(sprintf("Spine missing columns: %s", paste(miss, collapse = ", ")))

  is_sig <- as.character(df$significant) == "TRUE" & !is.na(df$FDR) & df$FDR < fdr
  sig <- df[is_sig, , drop = FALSE]
  for (cc in c("start1", "end1", "start2", "end2")) sig[[cc]] <- as.integer(sig[[cc]])
  cat(sprintf("  Spine: %d total -> %d significant (FDR < %g)\n", nrow(df), nrow(sig), fdr))
  cat(sprintf("    up_in_mutant=%d  down_in_mutant=%d\n",
              sum(sig$direction == "up_in_mutant"), sum(sig$direction == "down_in_mutant")))
  sig
}

# Per-anchor ChIP overlaps + TSS distance + 8-state anchor/loop types (previous information).
annotate_chromatin_state <- function(a1, a2, peak_files, tss_threshold) {
  cat("  Loading ChIP peak beds...\n")
  pk <- lapply(names(peak_files), function(k) load_bed_granges(peak_files[[k]], k))
  names(pk) <- names(peak_files)

  flags <- function(anchor) data.frame(
    H3K27ac           = overlaps_any(anchor, pk$h3k27ac),
    H3K27me3          = overlaps_any(anchor, pk$h3k27me3),
    H3K4me1           = overlaps_any(anchor, pk$h3k4me1),
    H3K4me3           = overlaps_any(anchor, pk$h3k4me3),
    Bivalent_Promoter = overlaps_any(anchor, pk$bivalent),
    CTCF              = overlaps_any(anchor, pk$ctcf),
    CTCF_motif        = overlaps_any(anchor, pk$ctcf_motif)
  )
  f1 <- flags(a1); f2 <- flags(a2)

  cat("  Computing TSS distances (TxDb mm10 knownGene)...\n")
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  tss  <- resize(suppressWarnings(genes(txdb)), width = 1, fix = "start")
  tss_dist <- function(anchor) {
    d <- suppressWarnings(distanceToNearest(anchor, tss))
    out <- rep(NA_real_, length(anchor))
    if (length(d) > 0) out[queryHits(d)] <- mcols(d)$distance
    out
  }
  d1 <- tss_dist(a1); d2 <- tss_dist(a2)

  t1 <- classify_anchor_type_extended(f1$H3K27ac, f1$H3K27me3, f1$H3K4me1, f1$H3K4me3,
                                      f1$Bivalent_Promoter, f1$CTCF, f1$CTCF_motif,
                                      d1, tss_threshold, use_motif_for_ctcf = TRUE)
  t2 <- classify_anchor_type_extended(f2$H3K27ac, f2$H3K27me3, f2$H3K4me1, f2$H3K4me3,
                                      f2$Bivalent_Promoter, f2$CTCF, f2$CTCF_motif,
                                      d2, tss_threshold, use_motif_for_ctcf = TRUE)
  lt <- classify_loop_type_extended(t1, t2)

  colnames(f1) <- paste0("anchor1_", colnames(f1), "_overlap")
  colnames(f2) <- paste0("anchor2_", colnames(f2), "_overlap")
  cbind(
    f1, f2,
    data.frame(
      anchor1_distance_to_tss = d1, anchor2_distance_to_tss = d2,
      anchor1_type = t1, anchor2_type = t2, loop_type = lt,
      stringsAsFactors = FALSE
    )
  )
}

# loop_class per cluster/scripts/05_grouped_analyses.py 4.2 (CTCF-only structural rule).
# Reuses the already-computed CTCF ChIP flags; promoter from mm10_knownGene_pp.bed.
compute_loop_class <- function(a1, a2, a1_ctcf, a2_ctcf, a1_h3k27ac, a2_h3k27ac, promoter_bed) {
  prom <- load_bed_granges(promoter_bed, "promoter(pp)")
  a1_EorP <- a1_h3k27ac | overlaps_any(a1, prom)
  a2_EorP <- a2_h3k27ac | overlaps_any(a2, prom)
  ctcf_sum <- as.integer(a1_ctcf) + as.integer(a2_ctcf)
  eorp_sum <- as.integer(a1_EorP) + as.integer(a2_EorP)

  cls <- rep("unclassified", length(a1))
  cls[ctcf_sum == 2 & eorp_sum <  2] <- "structural"
  cls[ctcf_sum <  2 & eorp_sum == 2] <- "CRE"
  cls[ctcf_sum == 2 & eorp_sum == 2] <- "mixed"
  list(loop_class = cls, anchor1_is_EorP = a1_EorP, anchor2_is_EorP = a2_EorP)
}

# shared_anchor re-derived on THIS significant set: a maxgap-merged locus that holds at least
# one lost-loop anchor AND one gained-loop anchor; flag loops touching such loci.
compute_shared_anchor <- function(loops_df, a1, a2, maxgap) {
  all_anchors <- c(a1, a2)
  dir2 <- rep(loops_df$direction, 2L)
  lost   <- all_anchors[dir2 == "down_in_mutant"]
  gained <- all_anchors[dir2 == "up_in_mutant"]

  loci <- reduce(all_anchors, min.gapwidth = maxgap)
  n_lost   <- countOverlaps(loci, lost)
  n_gained <- countOverlaps(loci, gained)
  shared_loci <- loci[n_lost >= 1 & n_gained >= 1]
  cat(sprintf("  Shared (loop-switching) loci: %d (from %d lost + %d gained anchors)\n",
              length(shared_loci), length(lost), length(gained)))

  overlaps_any(a1, shared_loci) | overlaps_any(a2, shared_loci)
}

# DEG gene bodies (gencode vM25) for padj<thr DEGs from the adult RNA-seq XLSX.
load_deg_and_genes <- function(deg_xlsx, gencode_path, padj_thr, lfc_thr) {
  require_file(deg_xlsx, "DEG XLSX"); require_file(gencode_path, "gencode genes bed")

  rna <- as.data.frame(read_excel(deg_xlsx))
  if (!"ensembl_gene_id" %in% colnames(rna))           # column holds SYMBOLS (per load_rnaseq_degs)
    stop("DEG XLSX missing 'ensembl_gene_id' column")
  for (cc in c("log2FoldChange", "padj"))
    if (!cc %in% colnames(rna)) stop(sprintf("DEG XLSX missing '%s'", cc))

  keep <- !is.na(rna$padj) & rna$padj < padj_thr
  if (lfc_thr > 0)                                   # 0 = padj-only (PI choice); NA-safe when applied
    keep <- keep & !is.na(rna$log2FoldChange) & abs(rna$log2FoldChange) > lfc_thr
  deg_symbols <- unique(rna$ensembl_gene_id[keep])
  deg_symbols <- deg_symbols[!is.na(deg_symbols)]
  cat(sprintf("  DEGs: %d (padj < %g, |log2FC| > %g)\n", length(deg_symbols), padj_thr, lfc_thr))

  g <- read.table(gzfile(gencode_path), header = TRUE, sep = "\t", quote = "",
                  comment.char = "", stringsAsFactors = FALSE)
  g <- g[g$Annotation == "Gene", , drop = FALSE]
  chrom <- ifelse(grepl("^chr", g$Chromosome), g$Chromosome, paste0("chr", g$Chromosome))
  gene_gr_all <- GRanges(seqnames = chrom,
                         ranges = IRanges(start = as.integer(g$Start), end = as.integer(g$End)),
                         Name = g$Name)

  matched <- intersect(deg_symbols, mcols(gene_gr_all)$Name)
  cat(sprintf("  DEG symbol -> gencode coordinate match: %d / %d (%.1f%%)\n",
              length(matched), length(deg_symbols),
              100 * length(matched) / max(1, length(deg_symbols))))
  deg_gene_gr <- gene_gr_all[mcols(gene_gr_all)$Name %in% deg_symbols]

  list(deg_gene_gr = deg_gene_gr, gene_gr_all = gene_gr_all)
}

# DEG_status (either anchor overlaps a DEG gene body) + the contributing gene symbols.
compute_deg_status <- function(a1, a2, deg_gene_gr) {
  names_for <- function(anchor) {
    h <- findOverlaps(anchor, deg_gene_gr)
    nm <- rep(NA_character_, length(anchor))
    if (length(h) > 0) {
      by_q <- split(mcols(deg_gene_gr)$Name[subjectHits(h)], queryHits(h))
      nm[as.integer(names(by_q))] <- vapply(by_q, function(v) paste(sort(unique(v)), collapse = ","), "")
    }
    nm
  }
  n1 <- names_for(a1); n2 <- names_for(a2)
  status <- !is.na(n1) | !is.na(n2)
  genes <- mapply(function(x, y) {
    v <- unique(c(strsplit(ifelse(is.na(x), "", x), ",")[[1]],
                  strsplit(ifelse(is.na(y), "", y), ",")[[1]]))
    v <- v[v != ""]
    if (length(v) == 0) "" else paste(sort(v), collapse = ",")
  }, n1, n2, USE.NAMES = FALSE)
  list(DEG_status = status, DEG_genes = genes)
}

# Nearest gene name + distance (gencode gene bodies).
nearest_gene <- function(anchor, gene_gr_all) {
  d <- suppressWarnings(distanceToNearest(anchor, gene_gr_all))
  nm <- rep(NA_character_, length(anchor)); di <- rep(NA_real_, length(anchor))
  if (length(d) > 0) {
    nm[queryHits(d)] <- mcols(gene_gr_all)$Name[subjectHits(d)]
    di[queryHits(d)] <- mcols(d)$distance
  }
  list(name = nm, dist = di)
}

# k-means cluster join (on the 6 coordinate columns) + clust6 short/long split.
join_clusters <- function(loops_df, clusters_txt, clust6_threshold) {
  require_file(clusters_txt, "combined-clusters.txt")
  cl <- read.delim(clusters_txt, header = TRUE, sep = "\t", quote = "",
                   comment.char = "", stringsAsFactors = FALSE)
  cmap <- setNames(cl$GROUP, key6(cl$chr1, cl$x1, cl$x2, cl$chr2, cl$y1, cl$y2))
  cluster <- unname(cmap[key6(loops_df$chr1, loops_df$start1, loops_df$end1,
                              loops_df$chr2, loops_df$start2, loops_df$end2)])

  footprint <- loops_df$end2 - loops_df$start1
  label <- cluster
  is6 <- !is.na(cluster) & cluster == "clust6"
  label[is6] <- ifelse(footprint[is6] < clust6_threshold, "clust6_short", "clust6_long")

  cat(sprintf("  Cluster join: %d matched, %d NA\n",
              sum(!is.na(cluster)), sum(is.na(cluster))))
  list(cluster = cluster, cluster_label = label, loop_footprint = footprint)
}

# =============================================================================
# MAIN
# =============================================================================

main <- function() {
  cfg <- parse_cli(CONFIG)
  cat("================================================================================\n")
  cat("BUILD MASTER ANNOTATED LOOP LIST\n")
  cat("================================================================================\n")
  dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)

  # --- spine ---------------------------------------------------------------
  cat("\n[1] Significant spine\n")
  loops <- load_significant_spine(cfg$spine_tsv, cfg$fdr_threshold)
  if (any(loops$chr1 != loops$chr2))
    stop("Trans-chromosomal loop(s) in significant set — span_length undefined; investigate.")

  a1 <- GRanges(loops$chr1, IRanges(loops$start1, loops$end1))
  a2 <- GRanges(loops$chr2, IRanges(loops$start2, loops$end2))

  # --- geometry ------------------------------------------------------------
  cat("\n[2] Geometry\n")
  span_length <- loops$start2 - loops$end1     # interior gap (PI-specified span length)
  cat(sprintf("  span_length (start2-end1): median=%.0f  min=%.0f  max=%.0f bp\n",
              median(span_length), min(span_length), max(span_length)))

  # --- previous information: chromatin state -------------------------------
  cat("\n[3] Chromatin-state annotation (per-anchor ChIP + 8-state types)\n")
  state <- annotate_chromatin_state(a1, a2, cfg$peak_files, cfg$tss_threshold)

  # --- loop_class ----------------------------------------------------------
  cat("\n[4] Loop class (structural / CRE / mixed / unclassified)\n")
  klass <- compute_loop_class(a1, a2,
                              state$anchor1_CTCF_overlap, state$anchor2_CTCF_overlap,
                              state$anchor1_H3K27ac_overlap, state$anchor2_H3K27ac_overlap,
                              cfg$promoter_bed)
  print(table(klass$loop_class))

  # --- SE_anchor -----------------------------------------------------------
  cat("\n[5] SE_anchor\n")
  se_gr <- load_bed_granges(cfg$se_bed, "super-enhancers")
  SE_anchor <- overlaps_any(a1, se_gr) | overlaps_any(a2, se_gr)
  cat(sprintf("  SE_anchor TRUE: %d (%.1f%%)\n", sum(SE_anchor), 100 * mean(SE_anchor)))

  # --- shared_anchor -------------------------------------------------------
  cat("\n[6] shared_anchor (loop-switching loci, re-derived on this set)\n")
  shared_anchor <- compute_shared_anchor(loops, a1, a2, cfg$shared_maxgap)
  cat(sprintf("  shared_anchor TRUE: %d (%.1f%%)\n", sum(shared_anchor), 100 * mean(shared_anchor)))

  # --- DEG_status + nearest gene ------------------------------------------
  cat("\n[7] DEG_status + nearest gene (gencode vM25)\n")
  degset <- load_deg_and_genes(cfg$deg_xlsx, cfg$gencode_genes, cfg$deg_padj, cfg$deg_lfc)
  deg <- compute_deg_status(a1, a2, degset$deg_gene_gr)
  cat(sprintf("  DEG_status TRUE: %d (%.1f%%)\n", sum(deg$DEG_status), 100 * mean(deg$DEG_status)))
  ng1 <- nearest_gene(a1, degset$gene_gr_all); ng2 <- nearest_gene(a2, degset$gene_gr_all)

  # --- clusters ------------------------------------------------------------
  cat("\n[8] k-means cluster + clust6 short/long\n")
  clust <- join_clusters(loops, cfg$clusters_txt, cfg$clust6_threshold)

  # --- assemble ------------------------------------------------------------
  cat("\n[9] Assembling master table\n")
  carry <- intersect(c(
    "loop_id", "chr1", "start1", "end1", "chr2", "start2", "end2", "coord_string",
    "logFC", "logCPM", "F", "PValue", "FDR", "significant", "exploratory", "category",
    "direction", "resolution_kb", "n_resolutions_detected", "resolutions_list",
    "is_multi_resolution", "FDR_5kb", "FDR_10kb", "FDR_25kb",
    "logFC_5kb", "logFC_10kb", "logFC_25kb", "kept_from_resolution", "n_overlaps"
  ), colnames(loops))
  master <- loops[, carry, drop = FALSE]

  master$span_length     <- span_length
  master$loop_footprint  <- clust$loop_footprint
  master$cluster         <- clust$cluster
  master$cluster_label   <- clust$cluster_label
  master$loop_class      <- klass$loop_class
  master$SE_anchor       <- SE_anchor
  master$shared_anchor   <- shared_anchor
  master$DEG_status      <- deg$DEG_status
  master$DEG_genes       <- deg$DEG_genes
  master <- cbind(master, state)
  master$anchor1_is_EorP <- klass$anchor1_is_EorP
  master$anchor2_is_EorP <- klass$anchor2_is_EorP
  master$anchor1_nearest_gene      <- ng1$name
  master$anchor1_nearest_gene_dist <- ng1$dist
  master$anchor2_nearest_gene      <- ng2$name
  master$anchor2_nearest_gene_dist <- ng2$dist

  # --- write ---------------------------------------------------------------
  spine_out  <- file.path(cfg$out_dir, "significant_loops.tsv")
  master_out <- file.path(cfg$out_dir, "master_annotated_loops.tsv")
  write.table(loops[, carry, drop = FALSE], spine_out, sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(master, master_out, sep = "\t", quote = FALSE, row.names = FALSE)

  cat("\n--------------------------------------------------------------------------------\n")
  cat(sprintf("  Wrote %s  (%d rows)\n", spine_out, nrow(loops)))
  cat(sprintf("  Wrote %s  (%d rows x %d cols)\n", master_out, nrow(master), ncol(master)))
  cat("  cluster_label:\n"); print(table(master$cluster_label, useNA = "ifany"))
  cat("Done.\n")
}

if (!interactive()) main()
