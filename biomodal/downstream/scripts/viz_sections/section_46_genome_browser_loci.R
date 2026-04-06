# biomodal/downstream/scripts/viz_sections/section_46_genome_browser_loci.R
# Multi-omic genome browser views at key gene loci using Gviz
# Renders stacked tracks: 5mC/5hmC (condition-averaged BigWigs), difference,
# RNA-seq coverage, ChIP/epigenomic peak annotations, CpG density, Hi-C loop arcs
# Run from downstream/ directory

source("scripts/viz_sections/_shared_config.R")

cat("================================================================================\n")
cat("SECTION 46: MULTI-OMIC GENOME BROWSER LOCUS VIEWS\n")
cat("================================================================================\n\n")

# Additional packages for genome browser rendering
suppressPackageStartupMessages({
  library(Gviz)
  library(GenomicInteractions)
  library(AnnotationDbi)
  library(readxl)
  library(grid)
})

# GenomicInteractions lacks a seqinfo method, which Gviz requires for
# InteractionTrack. Define it by delegating to the underlying regions GRanges.
if (!hasMethod("seqinfo", "GenomicInteractions")) {
  setMethod("seqinfo", "GenomicInteractions", function(x) seqinfo(regions(x)))
}

# =============================================================================
# CONFIGURATION
# =============================================================================

SECTION_OUTPUT_DIR <- file.path(BASE_DIR, "plots/visualizations/46_genome_browser")
dir.create(SECTION_OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

CACHE_DIR <- file.path(SECTION_OUTPUT_DIR, ".cache")
dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)

# --- Methylation BigWig paths ------------------------------------------------
MC_BW_DIR  <- "/Users/zakiralibhai/Documents/BIO_LAB/methylation-tracks/bigwigs/mc"
HMC_BW_DIR <- "/Users/zakiralibhai/Documents/BIO_LAB/methylation-tracks/bigwigs/hmc"

MC_BW_CTRL <- file.path(MC_BW_DIR, c(
  "evoC-Bap1-ctrl-F.genome.mc_bedgraph.bw",
  "evoC-Bap1-ctrl-M.genome.mc_bedgraph.bw",
  "evoC-Bap1-ctrl-F-B2.genome.mc_bedgraph.bw",
  "evoC-Bap1-ctrl-M-B2.genome.mc_bedgraph.bw"))

MC_BW_MUT <- file.path(MC_BW_DIR, c(
  "evoC-Bap1-mut-F.genome.mc_bedgraph.bw",
  "evoC-Bap1-mut-M.genome.mc_bedgraph.bw",
  "evoC-Bap1-mut-F-B2.genome.mc_bedgraph.bw",
  "evoC-Bap1-mut-M-B2.genome.mc_bedgraph.bw"))

HMC_BW_CTRL <- file.path(HMC_BW_DIR, c(
  "evoC-Bap1-ctrl-F.genome.hmc_bedgraph.bw",
  "evoC-Bap1-ctrl-M.genome.hmc_bedgraph.bw",
  "evoC-Bap1-ctrl-F-B2.genome.hmc_bedgraph.bw",
  "evoC-Bap1-ctrl-M-B2.genome.hmc_bedgraph.bw"))

HMC_BW_MUT <- file.path(HMC_BW_DIR, c(
  "evoC-Bap1-mut-F.genome.hmc_bedgraph.bw",
  "evoC-Bap1-mut-M.genome.hmc_bedgraph.bw",
  "evoC-Bap1-mut-F-B2.genome.hmc_bedgraph.bw",
  "evoC-Bap1-mut-M-B2.genome.hmc_bedgraph.bw"))

# --- RNA-seq BigWig paths ----------------------------------------------------
RNASEQ_BW_CTRL <- file.path(BASE_DIR, "peaks/RNActrl.bw")
RNASEQ_BW_MUT  <- file.path(BASE_DIR, "peaks/RNAmut.bw")

# --- ChIP/ATAC/MeCP2 BigWig paths ---------------------------------------------
# Each is a named list with ctrl/mut paths. A continuous coverage DataTrack
# replaces the BED AnnotationTrack for that mark.
HISTONE_BW_DIR <- "/Users/zakiralibhai/Documents/BIO_LAB/methylation-tracks/histone_mods"
H2AK119UB_BW <- list(ctrl = file.path(HISTONE_BW_DIR, "H2AK119ubCtrl.bw"),
                      mut  = file.path(HISTONE_BW_DIR, "H2AK119ubMut.bw"))
H3K27ME3_BW  <- list(ctrl = file.path(HISTONE_BW_DIR, "H3K27me3Ctrl.bw"),
                      mut  = file.path(HISTONE_BW_DIR, "H3K27me3Mut.bw"))
H3K4ME3_BW   <- list(ctrl = file.path(HISTONE_BW_DIR, "H3K4me3Ctrl.bw"),
                      mut  = file.path(HISTONE_BW_DIR, "H3K4me3Mut.bw"))
H3K27AC_BW   <- list(ctrl = file.path(HISTONE_BW_DIR, "H3K27acCtrl.bw"),
                      mut  = file.path(HISTONE_BW_DIR, "H3K27acMut.bw"))
ATAC_BW      <- list(ctrl = file.path(HISTONE_BW_DIR, "ATACctrl.bw"),
                      mut  = file.path(HISTONE_BW_DIR, "ATACmut.bw"))
H3K27ME1_BW  <- list(ctrl = file.path(HISTONE_BW_DIR, "H3K27me1Ctrl.bw"),
                      mut  = file.path(HISTONE_BW_DIR, "H3K27me1Mut.bw"))
MECP2_BW     <- list(ctrl = file.path(HISTONE_BW_DIR, "MeCP2Ctrl.bw"),
                      mut  = file.path(HISTONE_BW_DIR, "MeCP2Mut.bw"))

# --- View parameters ---------------------------------------------------------
EXTEND_BP <- 50000   # 50kb flanking on each side

# --- Track colors (extending shared COLORS palette) --------------------------
TRACK_COL <- list(
  ctrl       = "#2166AC",
  mut        = "#B2182B",
  diff_hyper = "#D7191C",
  diff_hypo  = "#2C7BB6",
  k119ub_ctrl = "#BDBDBD",
  k119ub_mut  = "#756BB1",
  k27ac_ctrl  = "#1F78B4",
  k27ac_mut   = "#FF7F00",
  k27me3      = "#4DAF4A",
  ctcf        = "#1B9E77",
  atac_up     = "#E6AB02",
  atac_down   = "#66A61E",
  mecp2_up    = "#D95F02",
  mecp2_down  = "#7570B3",
  cpg_island  = "#B2182B",
  cpg_shore   = "#EF8A62",
  cpg_shelf   = "#FDDBC7",
  hic_lost    = "#B2182B",
  hic_gained  = "#2166AC",
  gene        = "#08306B",
  rnaseq      = "#984EA3"
)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' Look up gene coordinates from TxDb and extend by flanking region
get_gene_region <- function(gene_symbol, extend_bp = EXTEND_BP) {
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  entrez <- AnnotationDbi::mapIds(
    org.Mm.eg.db, keys = gene_symbol,
    column = "ENTREZID", keytype = "SYMBOL")

  stopifnot(!is.na(entrez))

  all_genes <- suppressMessages(genes(txdb))
  gene_gr <- all_genes[names(all_genes) == entrez]
  stopifnot(length(gene_gr) == 1)

  # Extend by flanking
  extended <- gene_gr
  start(extended) <- max(1, start(gene_gr) - extend_bp)
  end(extended)   <- end(gene_gr) + extend_bp

  mcols(extended)$symbol <- gene_symbol
  mcols(extended)$gene_start <- start(gene_gr)
  mcols(extended)$gene_end   <- end(gene_gr)
  return(extended)
}


#' Average methylation BigWig signal across samples within a genomic region.
#' Uses viewSums(weighted)/viewSums(count) for mean-at-CpG-sites (NOT viewMeans).
average_bigwigs_in_region <- function(bw_paths, region_gr, bin_size = NULL) {
  chr <- as.character(seqnames(region_gr)[1])
  region_start <- start(region_gr)
  region_end   <- end(region_gr)
  region_width <- width(region_gr)

  if (is.null(bin_size)) {
    bin_size <- max(50, round(region_width / 1000))
  }

  bin_starts <- seq(region_start, region_end - bin_size, by = bin_size)
  bins_gr <- GRanges(seqnames = chr,
                     ranges = IRanges(start = bin_starts, width = bin_size))

  valid_paths <- bw_paths[file.exists(bw_paths)]
  if (length(valid_paths) == 0) {
    stop("No BigWig files found. Paths checked:\n",
         paste(" ", bw_paths, collapse = "\n"))
  }

  score_matrix <- sapply(valid_paths, function(path) {
    bw <- rtracklayer::import.bw(path, which = region_gr)
    if (length(bw) == 0) return(rep(NA_real_, length(bins_gr)))

    weighted_cov <- coverage(bw, weight = "score")
    count_cov    <- coverage(bw)
    if (!chr %in% names(weighted_cov)) return(rep(NA_real_, length(bins_gr)))

    v_weighted <- Views(weighted_cov[[chr]], ranges(bins_gr))
    v_count    <- Views(count_cov[[chr]], ranges(bins_gr))

    sum_scores <- as.numeric(viewSums(v_weighted))
    num_cpgs   <- as.numeric(viewSums(v_count))

    ifelse(num_cpgs > 0, sum_scores / num_cpgs, NA_real_)
  })

  if (is.matrix(score_matrix)) {
    bins_gr$score <- rowMeans(score_matrix, na.rm = TRUE)
  } else {
    bins_gr$score <- score_matrix
  }
  bins_gr$score[is.na(bins_gr$score)] <- 0
  return(bins_gr)
}


#' Import a single RNA-seq BigWig and bin using viewMeans.
#' Unlike methylation BigWigs (sparse CpG-site scores), RNA-seq BigWigs have
#' continuous coverage — viewMeans (signal/bin_width) is the correct approach.
import_rnaseq_bigwig <- function(bw_path, region_gr, bin_size = NULL) {
  chr <- as.character(seqnames(region_gr)[1])
  region_width <- width(region_gr)
  if (is.null(bin_size)) bin_size <- max(50, round(region_width / 1000))

  bin_starts <- seq(start(region_gr), end(region_gr) - bin_size, by = bin_size)
  bins_gr <- GRanges(seqnames = chr,
                     ranges = IRanges(start = bin_starts, width = bin_size))

  bw <- rtracklayer::import.bw(bw_path, which = region_gr)
  if (length(bw) == 0) { bins_gr$score <- 0; return(bins_gr) }

  cov <- coverage(bw, weight = "score")
  if (!chr %in% names(cov)) { bins_gr$score <- 0; return(bins_gr) }

  bins_gr$score <- as.numeric(viewMeans(Views(cov[[chr]], ranges(bins_gr))))
  bins_gr$score[is.na(bins_gr$score)] <- 0
  return(bins_gr)
}


#' Load a BED file as GRanges, adding chr prefix if needed.
#' Filters out non-numeric coordinate rows (header rows, NAs).
load_peak_as_granges <- function(path, name = "peak") {
  if (!file.exists(path)) {
    warning("Peak file not found: ", path)
    return(GRanges())
  }
  df <- read.table(path, header = FALSE, sep = "\t",
                   stringsAsFactors = FALSE, fill = TRUE, comment.char = "#")

  # Filter non-numeric coordinate rows (header rows, NAs)
  df$V2 <- suppressWarnings(as.numeric(df$V2))
  df$V3 <- suppressWarnings(as.numeric(df$V3))
  df <- df[!is.na(df$V1) & !is.na(df$V2) & !is.na(df$V3), ]
  if (nrow(df) == 0) return(GRanges())

  gr <- GRanges(seqnames = df$V1,
                ranges = IRanges(start = df$V2, end = df$V3))

  if (length(gr) > 0 && !any(grepl("^chr", seqnames(gr)))) {
    seqlevels(gr) <- paste0("chr", seqlevels(gr))
  }
  mcols(gr)$name <- name
  return(gr)
}


#' Load CpG annotation BED (has header, no chr prefix)
load_cpg_annotation <- function(path) {
  if (grepl("\\.gz$", path)) {
    df <- read.table(gzfile(path), header = TRUE, sep = "\t",
                     stringsAsFactors = FALSE)
  } else {
    df <- read.table(path, header = TRUE, sep = "\t",
                     stringsAsFactors = FALSE)
  }
  df$Chromosome <- paste0("chr", df$Chromosome)
  GRanges(seqnames = df$Chromosome,
          ranges = IRanges(start = df$Start, end = df$End),
          annotation = df$Annotation)
}


#' Build Hi-C loop InteractionTracks for a region, split by direction
build_hic_tracks <- function(loops_df, region_gr) {
  chr <- as.character(seqnames(region_gr)[1])
  rstart <- start(region_gr)
  rend   <- end(region_gr)

  in_region <- loops_df %>%
    dplyr::filter(
      (chr1 == chr & ((start1 >= rstart & start1 <= rend) |
                      (end1 >= rstart & end1 <= rend))) |
      (chr2 == chr & ((start2 >= rstart & start2 <= rend) |
                      (end2 >= rstart & end2 <= rend)))
    )

  if (nrow(in_region) == 0) return(NULL)

  tracks <- list()

  lost <- in_region %>% dplyr::filter(direction == "down_in_mutant")
  if (nrow(lost) > 0) {
    a1 <- GRanges(lost$chr1, IRanges(lost$start1, lost$end1))
    a2 <- GRanges(lost$chr2, IRanges(lost$start2, lost$end2))
    genome(a1) <- "mm10"
    genome(a2) <- "mm10"
    gi_lost <- GenomicInteractions(a1, a2)
    gi_lost$name <- lost$loop_id
    tracks$lost <- InteractionTrack(gi_lost, name = "Lost Loops",
                                    chromosome = chr)
    displayPars(tracks$lost) <- list(
      col.interactions = TRACK_COL$hic_lost,
      col.anchors.fill = TRACK_COL$hic_lost,
      col.anchors.line = TRACK_COL$hic_lost,
      interaction.dimension = "height",
      interaction.measure = "counts",
      plot.anchors = FALSE,
      plot.trans = FALSE,
      plot.outside = TRUE,
      col.outside = adjustcolor(TRACK_COL$hic_lost, alpha.f = 0.3),
      anchor.height = 0.1
    )
  }

  gained <- in_region %>% dplyr::filter(direction == "up_in_mutant")
  if (nrow(gained) > 0) {
    a1 <- GRanges(gained$chr1, IRanges(gained$start1, gained$end1))
    a2 <- GRanges(gained$chr2, IRanges(gained$start2, gained$end2))
    genome(a1) <- "mm10"
    genome(a2) <- "mm10"
    gi_gained <- GenomicInteractions(a1, a2)
    gi_gained$name <- gained$loop_id
    tracks$gained <- InteractionTrack(gi_gained, name = "Gained Loops",
                                      chromosome = chr)
    displayPars(tracks$gained) <- list(
      col.interactions = TRACK_COL$hic_gained,
      col.anchors.fill = TRACK_COL$hic_gained,
      col.anchors.line = TRACK_COL$hic_gained,
      interaction.dimension = "height",
      interaction.measure = "counts",
      plot.anchors = FALSE,
      plot.trans = FALSE,
      plot.outside = TRUE,
      col.outside = adjustcolor(TRACK_COL$hic_gained, alpha.f = 0.3),
      anchor.height = 0.1
    )
  }

  return(tracks)
}


#' Create a difference DataTrack with baseline at zero.
#' NOTE: Gviz histogram DataTracks don't support per-bar coloring via groups
#' (collapseTrack triggers seqnames NAs). Using single color with baseline.
make_diff_track <- function(diff_gr, name, color, ylim,
                            bg_title, cex_title = 0.7) {
  DataTrack(
    range = diff_gr, genome = "mm10", type = "histogram",
    name = name,
    fill.histogram = color, col.histogram = color,
    ylim = ylim, baseline = 0, col.baseline = "black", lwd.baseline = 0.5,
    background.title = bg_title, col.title = "white", cex.title = cex_title
  )
}


# =============================================================================
# CACHING: Precompute and cache BigWig averages per gene
# =============================================================================

#' Compute all BigWig averages for a gene region (methylation + RNA-seq)
precompute_bigwig_averages <- function(region_gr) {
  mc_ctrl  <- average_bigwigs_in_region(MC_BW_CTRL, region_gr)
  mc_mut   <- average_bigwigs_in_region(MC_BW_MUT, region_gr)
  hmc_ctrl <- average_bigwigs_in_region(HMC_BW_CTRL, region_gr)
  hmc_mut  <- average_bigwigs_in_region(HMC_BW_MUT, region_gr)

  # Compute differences before scaling
  mc_diff <- mc_ctrl
  mc_diff$score <- mc_mut$score - mc_ctrl$score
  hmc_diff <- hmc_ctrl
  hmc_diff$score <- hmc_mut$score - hmc_ctrl$score

  # Scale to percentage (BigWig values are 0-1 fractions)
  mc_ctrl$score  <- mc_ctrl$score * 100
  mc_mut$score   <- mc_mut$score * 100
  mc_diff$score  <- mc_diff$score * 100
  hmc_ctrl$score <- hmc_ctrl$score * 100
  hmc_mut$score  <- hmc_mut$score * 100
  hmc_diff$score <- hmc_diff$score * 100

  result <- list(
    mc_ctrl = mc_ctrl, mc_mut = mc_mut, mc_diff = mc_diff,
    hmc_ctrl = hmc_ctrl, hmc_mut = hmc_mut, hmc_diff = hmc_diff
  )

  # RNA-seq BigWigs (single file per condition, different import method)
  if (!is.null(RNASEQ_BW_CTRL) && file.exists(RNASEQ_BW_CTRL)) {
    result$rnaseq_ctrl <- import_rnaseq_bigwig(RNASEQ_BW_CTRL, region_gr)
    result$rnaseq_mut  <- import_rnaseq_bigwig(RNASEQ_BW_MUT, region_gr)
  }

  # ChIP/ATAC BigWigs (ctrl/mut pairs, same import method as RNA-seq)
  chip_bw_configs <- list(
    h2ak119ub = H2AK119UB_BW, h3k27me3 = H3K27ME3_BW,
    h3k4me3 = H3K4ME3_BW, h3k27ac = H3K27AC_BW,
    atac = ATAC_BW, h3k27me1 = H3K27ME1_BW,
    mecp2 = MECP2_BW
  )
  for (mark_name in names(chip_bw_configs)) {
    bw_config <- chip_bw_configs[[mark_name]]
    if (!is.null(bw_config) && file.exists(bw_config$ctrl) && file.exists(bw_config$mut)) {
      result[[paste0(mark_name, "_ctrl")]] <- import_rnaseq_bigwig(bw_config$ctrl, region_gr)
      result[[paste0(mark_name, "_mut")]]  <- import_rnaseq_bigwig(bw_config$mut, region_gr)
    }
  }

  return(result)
}

#' Load cached BigWig averages or compute and cache them
get_cached_bigwig_averages <- function(gene_symbol, region_gr) {
  cache_file <- file.path(CACHE_DIR, sprintf("%s_bw_cache.rds", gene_symbol))

  # Check if cache exists and is fresh (newer than all source BigWig files)
  if (file.exists(cache_file)) {
    cache_time <- file.mtime(cache_file)
    all_bw_paths <- c(MC_BW_CTRL, MC_BW_MUT, HMC_BW_CTRL, HMC_BW_MUT)
    if (!is.null(RNASEQ_BW_CTRL)) all_bw_paths <- c(all_bw_paths, RNASEQ_BW_CTRL, RNASEQ_BW_MUT)
    for (bw_cfg in list(H2AK119UB_BW, H3K27ME3_BW, H3K4ME3_BW, H3K27AC_BW, ATAC_BW, H3K27ME1_BW, MECP2_BW)) {
      if (!is.null(bw_cfg)) all_bw_paths <- c(all_bw_paths, bw_cfg$ctrl, bw_cfg$mut)
    }
    existing_bw <- all_bw_paths[file.exists(all_bw_paths)]
    newest_bw <- max(file.mtime(existing_bw))

    if (cache_time > newest_bw) {
      cat("    Loading cached BigWig averages...\n")
      return(readRDS(cache_file))
    }
  }

  # Compute fresh
  cat("    Averaging 5mC BigWigs...\n")
  cat("    Averaging 5hmC BigWigs...\n")
  if (!is.null(RNASEQ_BW_CTRL)) cat("    Averaging RNA-seq BigWigs...\n")
  result <- precompute_bigwig_averages(region_gr)
  saveRDS(result, cache_file)
  return(result)
}


# =============================================================================
# MAIN PLOTTING FUNCTION
# =============================================================================

#' Generate a multi-omic genome browser figure for a gene locus
#'
#' @param gene_symbol Gene name (e.g., "Syt1")
#' @param variant "full" (~15 tracks) or "compact" (~8 tracks)
#' @param loops_df Hi-C loop data.frame
#' @param rnaseq_df RNA-seq DESeq2 results data.frame
#' @param cpg_islands_gr CpG island GRanges
#' @param cpg_shores_gr CpG shore GRanges
#' @param cpg_shelves_gr CpG shelf GRanges
#' @param peak_data Pre-loaded peak GRanges list
#' @param precomputed_bw Pre-computed BigWig averages list
#' @param return_tracks If TRUE, return track list instead of rendering
plot_locus_browser <- function(gene_symbol, variant = "full",
                               loops_df, rnaseq_df,
                               cpg_islands_gr, cpg_shores_gr,
                               cpg_shelves_gr,
                               peak_data, precomputed_bw,
                               return_tracks = FALSE) {

  cat(sprintf("  Generating %s browser view for %s...\n", variant, gene_symbol))

  # --- 1. Get gene region ---
  region_gr <- get_gene_region(gene_symbol, extend_bp = EXTEND_BP)

  chr   <- as.character(seqnames(region_gr)[1])
  view_start <- start(region_gr)
  view_end   <- end(region_gr)
  gene_start <- mcols(region_gr)$gene_start
  gene_end   <- mcols(region_gr)$gene_end

  cat(sprintf("    Region: %s:%s-%s (%s kb)\n",
              chr, format(view_start, big.mark = ","),
              format(view_end, big.mark = ","),
              round(width(region_gr) / 1000)))

  # --- 2. RNA-seq annotation ---
  rnaseq_label <- ""
  if (!is.null(rnaseq_df)) {
    gene_expr <- rnaseq_df %>%
      dplyr::filter(ensembl_gene_id == gene_symbol)
    if (nrow(gene_expr) > 0) {
      lfc <- round(gene_expr$log2FoldChange[1], 2)
      padj <- signif(gene_expr$padj[1], 3)
      rnaseq_label <- sprintf("  [RNA-seq: log2FC=%s, padj=%s]", lfc, padj)
    }
  }

  # --- 3. Use precomputed BigWig averages ---
  mc_ctrl_gr  <- precomputed_bw$mc_ctrl
  mc_mut_gr   <- precomputed_bw$mc_mut
  mc_diff_gr  <- precomputed_bw$mc_diff
  hmc_ctrl_gr <- precomputed_bw$hmc_ctrl
  hmc_mut_gr  <- precomputed_bw$hmc_mut
  hmc_diff_gr <- precomputed_bw$hmc_diff

  # --- 4. Build Gviz tracks ---
  cat("    [DEBUG] Starting track construction...\n")
  track_list <- list()
  sizes <- c()

  # Track 0: Genome axis
  track_list$axis <- GenomeAxisTrack(name = "Position")
  sizes <- c(sizes, 0.5)

  # Track 1: Gene model (collapsed meta-transcripts)
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  track_list$gene <- GeneRegionTrack(
    txdb, genome = "mm10", chromosome = chr,
    name = paste0(gene_symbol, rnaseq_label),
    transcriptAnnotation = "symbol",
    stacking = "dense",
    col = TRACK_COL$gene,
    fill = TRACK_COL$gene,
    fontcolor.group = TRACK_COL$gene,
    background.title = "#08306B",
    col.title = "white",
    cex.title = 0.8
  )
  sizes <- c(sizes, 1.2)

  # --- Methylation block ---
  mc_ylim <- c(0, max(100, max(c(mc_ctrl_gr$score, mc_mut_gr$score), na.rm = TRUE)))
  hmc_ylim <- c(0, max(50, max(c(hmc_ctrl_gr$score, hmc_mut_gr$score), na.rm = TRUE)))
  diff_max <- max(abs(c(mc_diff_gr$score, hmc_diff_gr$score)), na.rm = TRUE)
  diff_ylim <- c(-diff_max, diff_max)

  # 5mC Control
  track_list$mc_ctrl <- DataTrack(
    range = mc_ctrl_gr, genome = "mm10", type = "histogram",
    name = "5mC%\nControl",
    fill.histogram = TRACK_COL$ctrl,
    col.histogram = TRACK_COL$ctrl,
    ylim = mc_ylim,
    background.title = TRACK_COL$ctrl,
    col.title = "white", cex.title = 0.7
  )
  sizes <- c(sizes, 1.5)

  # 5mC Mutant
  track_list$mc_mut <- DataTrack(
    range = mc_mut_gr, genome = "mm10", type = "histogram",
    name = "5mC%\nMutant",
    fill.histogram = TRACK_COL$mut,
    col.histogram = TRACK_COL$mut,
    ylim = mc_ylim,
    background.title = TRACK_COL$mut,
    col.title = "white", cex.title = 0.7
  )
  sizes <- c(sizes, 1.5)

  cat("    [DEBUG] mc ctrl/mut tracks OK\n")
  # 5mC Difference (Mut - Ctrl)
  track_list$mc_diff <- make_diff_track(
    mc_diff_gr, name = "5mC%\nDifference",
    color = TRACK_COL$diff_hyper,
    ylim = diff_ylim, bg_title = "#8C510A"
  )
  sizes <- c(sizes, 1.5)

  # 5hmC Control
  track_list$hmc_ctrl <- DataTrack(
    range = hmc_ctrl_gr, genome = "mm10", type = "histogram",
    name = "5hmC%\nControl",
    fill.histogram = TRACK_COL$ctrl,
    col.histogram = TRACK_COL$ctrl,
    ylim = hmc_ylim,
    background.title = TRACK_COL$ctrl,
    col.title = "white", cex.title = 0.7
  )
  sizes <- c(sizes, 1.5)

  if (variant == "full") {
    # 5hmC Mutant (full version only)
    track_list$hmc_mut <- DataTrack(
      range = hmc_mut_gr, genome = "mm10", type = "histogram",
      name = "5hmC%\nMutant",
      fill.histogram = TRACK_COL$mut,
      col.histogram = TRACK_COL$mut,
      ylim = hmc_ylim,
      background.title = TRACK_COL$mut,
      col.title = "white", cex.title = 0.7
    )
    sizes <- c(sizes, 1.5)

    # 5hmC Difference (Mut - Ctrl)
    track_list$hmc_diff <- make_diff_track(
      hmc_diff_gr, name = "5hmC%\nDifference",
      color = TRACK_COL$diff_hypo,
      ylim = diff_ylim, bg_title = "#01665E"
    )
    sizes <- c(sizes, 1.5)
  }

  cat("    [DEBUG] methylation tracks OK\n")
  # --- RNA-seq BigWig tracks (both variants) ---
  if (!is.null(RNASEQ_BW_CTRL) && !is.null(precomputed_bw$rnaseq_ctrl)) {
    rnaseq_ctrl_gr <- precomputed_bw$rnaseq_ctrl
    rnaseq_mut_gr  <- precomputed_bw$rnaseq_mut
    rnaseq_ylim <- c(0, max(c(rnaseq_ctrl_gr$score, rnaseq_mut_gr$score), na.rm = TRUE) * 1.05)
    if (rnaseq_ylim[2] == 0) rnaseq_ylim[2] <- 1

    track_list$rnaseq_ctrl <- DataTrack(
      range = rnaseq_ctrl_gr, genome = "mm10", type = "histogram",
      name = "RNA-seq\nControl",
      fill.histogram = TRACK_COL$rnaseq,
      col.histogram = TRACK_COL$rnaseq,
      ylim = rnaseq_ylim,
      background.title = TRACK_COL$ctrl,
      col.title = "white", cex.title = 0.7
    )
    sizes <- c(sizes, 1.2)

    track_list$rnaseq_mut <- DataTrack(
      range = rnaseq_mut_gr, genome = "mm10", type = "histogram",
      name = "RNA-seq\nMutant",
      fill.histogram = TRACK_COL$rnaseq,
      col.histogram = TRACK_COL$rnaseq,
      ylim = rnaseq_ylim,
      background.title = TRACK_COL$mut,
      col.title = "white", cex.title = 0.7
    )
    sizes <- c(sizes, 1.2)
  }

  cat("    [DEBUG] RNA-seq tracks OK\n")
  # --- Chromatin / epigenomic tracks ---
  # For each mark: use BigWig DataTrack if available, otherwise BED AnnotationTrack.
  # BigWig tracks show continuous coverage (like PI's IGV view); BED tracks show
  # peak calls as colored bars.

  # Helper: add paired ctrl/mut BigWig DataTracks for a chromatin mark
  add_chip_bw_tracks <- function(mark_name, mark_key, bg_color) {
    ctrl_key <- paste0(mark_key, "_ctrl")
    mut_key  <- paste0(mark_key, "_mut")
    if (!is.null(precomputed_bw[[ctrl_key]])) {
      ctrl_gr <- precomputed_bw[[ctrl_key]]
      mut_gr  <- precomputed_bw[[mut_key]]
      ylim <- c(0, max(c(ctrl_gr$score, mut_gr$score), na.rm = TRUE) * 1.05)
      if (ylim[2] == 0) ylim[2] <- 1

      track_list[[paste0(mark_key, "_ctrl")]] <<- DataTrack(
        range = ctrl_gr, genome = "mm10", type = "histogram",
        name = paste0(mark_name, "\nControl"),
        fill.histogram = TRACK_COL$ctrl, col.histogram = TRACK_COL$ctrl,
        ylim = ylim, background.title = bg_color,
        col.title = "white", cex.title = 0.7
      )
      sizes <<- c(sizes, 1.0)

      track_list[[paste0(mark_key, "_mut")]] <<- DataTrack(
        range = mut_gr, genome = "mm10", type = "histogram",
        name = paste0(mark_name, "\nMutant"),
        fill.histogram = TRACK_COL$mut, col.histogram = TRACK_COL$mut,
        ylim = ylim, background.title = bg_color,
        col.title = "white", cex.title = 0.7
      )
      sizes <<- c(sizes, 1.0)
      return(TRUE)
    }
    return(FALSE)
  }

  # H2AK119ub (always shown -- central to BAP1 mechanism)
  if (!add_chip_bw_tracks("H2AK119ub", "h2ak119ub", "#756BB1")) {
    k119_combined <- c(peak_data$k119ub_ctrl, peak_data$k119ub_mut)
    if (length(k119_combined) > 0) {
      k119_combined$feature <- c(
        rep("Control", length(peak_data$k119ub_ctrl)),
        rep("Mutant", length(peak_data$k119ub_mut))
      )
      track_list$k119ub <- AnnotationTrack(
        k119_combined, genome = "mm10", chromosome = chr,
        name = "H2AK119ub",
        feature = k119_combined$feature,
        Control = TRACK_COL$k119ub_ctrl,
        Mutant  = TRACK_COL$k119ub_mut,
        stacking = "dense",
        background.title = "#756BB1",
        col.title = "white", cex.title = 0.7
      )
      sizes <- c(sizes, 0.6)
    }
  }

  if (variant == "full") {
    # H3K27me3
    if (!add_chip_bw_tracks("H3K27me3", "h3k27me3", TRACK_COL$k27me3)) {
      if (length(peak_data$k27me3) > 0) {
        track_list$k27me3 <- AnnotationTrack(
          peak_data$k27me3, genome = "mm10", chromosome = chr,
          name = "H3K27me3",
          fill = TRACK_COL$k27me3, col = TRACK_COL$k27me3,
          stacking = "dense",
          background.title = TRACK_COL$k27me3,
          col.title = "white", cex.title = 0.7
        )
        sizes <- c(sizes, 0.4)
      }
    }

    # H3K4me3
    add_chip_bw_tracks("H3K4me3", "h3k4me3", "#E41A1C")

    # H3K27ac
    if (!add_chip_bw_tracks("H3K27ac", "h3k27ac", "#FF7F00")) {
      k27ac_combined <- c(peak_data$k27ac_ctrl, peak_data$k27ac_mut)
      if (length(k27ac_combined) > 0) {
        k27ac_combined$feature <- c(
          rep("Control", length(peak_data$k27ac_ctrl)),
          rep("Mutant", length(peak_data$k27ac_mut))
        )
        track_list$k27ac <- AnnotationTrack(
          k27ac_combined, genome = "mm10", chromosome = chr,
          name = "H3K27ac",
          feature = k27ac_combined$feature,
          Control = TRACK_COL$k27ac_ctrl,
          Mutant  = TRACK_COL$k27ac_mut,
          stacking = "dense",
          background.title = "#FF7F00",
          col.title = "white", cex.title = 0.7
        )
        sizes <- c(sizes, 0.6)
      }
    }

    # ATAC-seq
    if (!add_chip_bw_tracks("ATAC-seq", "atac", TRACK_COL$atac_up)) {
      atac_combined <- c(peak_data$atac_up, peak_data$atac_down)
      if (length(atac_combined) > 0) {
        atac_combined$feature <- c(
          rep("Up", length(peak_data$atac_up)),
          rep("Down", length(peak_data$atac_down))
        )
        track_list$atac <- AnnotationTrack(
          atac_combined, genome = "mm10", chromosome = chr,
          name = "ATAC-seq",
          feature = atac_combined$feature,
          Up   = TRACK_COL$atac_up,
          Down = TRACK_COL$atac_down,
          stacking = "dense",
          background.title = TRACK_COL$atac_up,
          col.title = "white", cex.title = 0.7
        )
        sizes <- c(sizes, 0.5)
      }
    }

    # H3K27me1
    add_chip_bw_tracks("H3K27me1", "h3k27me1", "#66A61E")

    # MeCP2 (BigWig signal if available, otherwise BED differential peaks)
    if (!add_chip_bw_tracks("MeCP2", "mecp2", TRACK_COL$mecp2_up)) {
      mecp2_combined <- c(peak_data$mecp2_up, peak_data$mecp2_down)
      if (length(mecp2_combined) > 0) {
        mecp2_combined$feature <- c(
          rep("Up", length(peak_data$mecp2_up)),
          rep("Down", length(peak_data$mecp2_down))
        )
        track_list$mecp2 <- AnnotationTrack(
          mecp2_combined, genome = "mm10", chromosome = chr,
          name = "MeCP2",
          feature = mecp2_combined$feature,
          Up   = TRACK_COL$mecp2_up,
          Down = TRACK_COL$mecp2_down,
          stacking = "dense",
          background.title = TRACK_COL$mecp2_up,
          col.title = "white", cex.title = 0.7
        )
        sizes <- c(sizes, 0.5)
      }
    }
  }

  cat("    [DEBUG] ChIP/epigenomic tracks OK\n")
  # CTCF peaks (both variants)
  if (length(peak_data$ctcf) > 0) {
    track_list$ctcf <- AnnotationTrack(
      peak_data$ctcf, genome = "mm10", chromosome = chr,
      name = "CTCF",
      fill = TRACK_COL$ctcf, col = TRACK_COL$ctcf,
      stacking = "dense",
      background.title = TRACK_COL$ctcf,
      col.title = "white", cex.title = 0.7
    )
    sizes <- c(sizes, 0.4)
  }

  cat("    [DEBUG] CTCF track OK\n")
  # --- CpG annotation track ---
  region_for_overlap <- GRanges(chr, IRanges(view_start, view_end))
  cpg_in_view <- subsetByOverlaps(cpg_islands_gr, region_for_overlap)
  shore_in_view <- subsetByOverlaps(cpg_shores_gr, region_for_overlap)
  shelf_in_view <- subsetByOverlaps(cpg_shelves_gr, region_for_overlap)

  cpg_all <- c(cpg_in_view, shore_in_view, shelf_in_view)
  if (length(cpg_all) > 0) {
    cpg_all$feature <- c(
      rep("Island", length(cpg_in_view)),
      rep("Shore", length(shore_in_view)),
      rep("Shelf", length(shelf_in_view))
    )
    track_list$cpg <- AnnotationTrack(
      cpg_all, genome = "mm10", chromosome = chr,
      name = "CpG",
      feature = cpg_all$feature,
      Island = TRACK_COL$cpg_island,
      Shore  = TRACK_COL$cpg_shore,
      Shelf  = TRACK_COL$cpg_shelf,
      stacking = "dense",
      background.title = TRACK_COL$cpg_island,
      col.title = "white", cex.title = 0.7
    )
    sizes <- c(sizes, 0.4)
  }

  cat("    [DEBUG] CpG annotation track OK\n")
  # --- Hi-C loop arcs ---
  hic_tracks <- build_hic_tracks(loops_df, region_gr)
  if (!is.null(hic_tracks)) {
    if (!is.null(hic_tracks$lost)) {
      track_list$hic_lost <- hic_tracks$lost
      sizes <- c(sizes, 1.5)
    }
    if (!is.null(hic_tracks$gained)) {
      track_list$hic_gained <- hic_tracks$gained
      sizes <- c(sizes, 1.5)
    }
  }

  cat("    [DEBUG] Hi-C tracks OK\n")
  # --- Return tracks for composite figure assembly ---
  if (return_tracks) {
    return(list(
      tracks = track_list, sizes = sizes,
      chr = chr, from = view_start, to = view_end,
      gene_symbol = gene_symbol, variant = variant
    ))
  }

  # --- 5. Render the plot ---
  suffix <- ifelse(variant == "compact", "_compact", "")
  out_name <- sprintf("%s_locus%s", gene_symbol, suffix)
  out_dir  <- file.path(SECTION_OUTPUT_DIR, out_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  file_prefix <- file.path(out_dir, out_name)

  fig_height <- sum(sizes) * 0.7 + 1
  fig_width  <- 14

  render_plot <- function() {
    plotTracks(
      track_list,
      chromosome = chr,
      from = view_start,
      to = view_end,
      sizes = sizes,
      main = sprintf("%s Locus -- BAP1-KO Multi-Omic View (%s)", gene_symbol, variant),
      cex.main = 1.0,
      title.width = 1.6,
      col.axis = "black",
      fontsize = 12
    )
  }

  # Render once to off-screen device, capture grob for replay
  pdf(NULL, width = fig_width, height = fig_height)
  tryCatch({
    render_plot()
    plot_grob <- grid.grab()
  }, finally = dev.off())

  # Replay captured grob to each output device (fast — no re-rendering)
  pdf(paste0(file_prefix, ".pdf"), width = fig_width, height = fig_height)
  tryCatch({ grid.newpage(); grid.draw(plot_grob) }, finally = dev.off())

  svglite::svglite(paste0(file_prefix, ".svg"), width = fig_width, height = fig_height)
  tryCatch({ grid.newpage(); grid.draw(plot_grob) }, finally = dev.off())

  jpeg(paste0(file_prefix, ".jpg"),
       width = fig_width * 300, height = fig_height * 300,
       res = 300, quality = 95)
  tryCatch({ grid.newpage(); grid.draw(plot_grob) }, finally = dev.off())

  cat(sprintf("    Saved: %s/{pdf,svg,jpg}\n", out_name))
  return(invisible(file_prefix))
}


# =============================================================================
# DATA LOADING
# =============================================================================

cat("Loading Hi-C loop data...\n")
loops_df <- read.table(
  file.path(REPO_ROOT, "peaks/loop_annotation_extended/late/extended_characterized_loops.tsv"),
  header = TRUE, sep = "\t", stringsAsFactors = FALSE
)
cat(sprintf("  %d differential loops loaded\n", nrow(loops_df)))

cat("Loading RNA-seq DESeq2 results...\n")
rnaseq_path <- file.path(REPO_ROOT, "tads/adult_timepoint_rna-seq-BAP1_WT_KO_v2_Results.xlsx")
rnaseq_df <- NULL
if (file.exists(rnaseq_path)) {
  rnaseq_df <- readxl::read_excel(rnaseq_path)
  cat(sprintf("  %d genes loaded\n", nrow(rnaseq_df)))
} else {
  warning("RNA-seq results file not found: ", rnaseq_path)
}

cat("Loading CpG annotations...\n")
cpg_base <- file.path(BASE_DIR, "modality/mm10")
cpg_islands_gr <- load_cpg_annotation(file.path(cpg_base, "gencode.vM25.mouse.cpg_islands.annotation.bed"))
cpg_shores_gr  <- load_cpg_annotation(file.path(cpg_base, "gencode.vM25.mouse.cpg_shores.annotation.bed.gz"))
cpg_shelves_gr <- load_cpg_annotation(file.path(cpg_base, "gencode.vM25.mouse.cpg_shelves.annotation.bed.gz"))
cat(sprintf("  CpG: %d islands, %d shores, %d shelves\n",
            length(cpg_islands_gr), length(cpg_shores_gr), length(cpg_shelves_gr)))

# Verify methylation BigWig files exist
mc_bw_ok  <- all(file.exists(c(MC_BW_CTRL, MC_BW_MUT)))
hmc_bw_ok <- all(file.exists(c(HMC_BW_CTRL, HMC_BW_MUT)))
if (!mc_bw_ok)  stop("Missing 5mC BigWig files. Check MC_BW_DIR path.")
if (!hmc_bw_ok) stop("Missing 5hmC BigWig files. Check HMC_BW_DIR path.")
cat("  All 16 methylation BigWig files verified\n")

# Verify RNA-seq BigWig files
if (!is.null(RNASEQ_BW_CTRL) && file.exists(RNASEQ_BW_CTRL) && file.exists(RNASEQ_BW_MUT)) {
  cat("  RNA-seq BigWig files verified\n")
} else {
  warning("RNA-seq BigWig files not found, disabling RNA-seq tracks")
  RNASEQ_BW_CTRL <- NULL
  RNASEQ_BW_MUT  <- NULL
}

# Pre-load all peak files once (instead of per-gene)
cat("Loading peak files...\n")
peak_data <- list(
  k119ub_ctrl = load_peak_as_granges(K119UB_FILES$ctrl, "K119ub_Ctrl"),
  k119ub_mut  = load_peak_as_granges(K119UB_FILES$mut, "K119ub_Mut"),
  k27me3      = load_peak_as_granges(CHIP_PEAK_FILES$h3k27me3, "H3K27me3"),
  k27ac_ctrl  = load_peak_as_granges(H3K27AC_FILES$ctrl, "K27ac_Ctrl"),
  k27ac_mut   = load_peak_as_granges(H3K27AC_FILES$mut, "K27ac_Mut"),
  atac_up     = load_peak_as_granges(ATAC_FILES$up, "ATAC_Up"),
  atac_down   = load_peak_as_granges(ATAC_FILES$down, "ATAC_Down"),
  mecp2_up    = load_peak_as_granges(MECP2_FILES$up, "MeCP2_Up"),
  mecp2_down  = load_peak_as_granges(MECP2_FILES$down, "MeCP2_Down"),
  ctcf        = load_peak_as_granges(CHIP_PEAK_FILES$ctcf, "CTCF")
)
cat("  Peak files loaded\n")

# =============================================================================
# GENERATE BROWSER VIEWS
# =============================================================================

# CLI argument: optional gene name to render only that gene (for HPC parallelism)
# Usage: Rscript section_46_genome_browser_loci.R [gene_name]
# No argument = all KEY_GENES + composite figure
args <- commandArgs(trailingOnly = TRUE)
target_genes <- KEY_GENES
run_composite <- TRUE
if (length(args) >= 1 && nchar(args[1]) > 0) {
  stopifnot(args[1] %in% KEY_GENES)
  target_genes <- args[1]
  run_composite <- FALSE
  cat(sprintf("Single-gene mode: %s\n", args[1]))
}

cat("\n--- Generating genome browser views ---\n\n")

for (gene in target_genes) {
  cat(sprintf("Processing %s...\n", gene))

  # Get gene region and precompute BigWig averages ONCE per gene
  region_gr <- get_gene_region(gene, extend_bp = EXTEND_BP)
  precomputed <- get_cached_bigwig_averages(gene, region_gr)

  # Full version
  plot_locus_browser(
    gene, variant = "full",
    loops_df = loops_df, rnaseq_df = rnaseq_df,
    cpg_islands_gr = cpg_islands_gr,
    cpg_shores_gr = cpg_shores_gr,
    cpg_shelves_gr = cpg_shelves_gr,
    peak_data = peak_data,
    precomputed_bw = precomputed
  )

  # Compact version (reuses same precomputed data)
  plot_locus_browser(
    gene, variant = "compact",
    loops_df = loops_df, rnaseq_df = rnaseq_df,
    cpg_islands_gr = cpg_islands_gr,
    cpg_shores_gr = cpg_shores_gr,
    cpg_shelves_gr = cpg_shelves_gr,
    peak_data = peak_data,
    precomputed_bw = precomputed
  )

  cat("\n")
}

# =============================================================================
# COMPOSITE MULTI-PANEL FIGURE
# =============================================================================

if (run_composite) {
cat("\n--- Generating composite figure ---\n\n")

# ---- Panel A: Syt1 compact browser view as grob ----------------------------
cat("  Building Panel A (Syt1 browser)...\n")
syt1_region <- get_gene_region("Syt1", extend_bp = EXTEND_BP)
syt1_bw <- get_cached_bigwig_averages("Syt1", syt1_region)

track_info <- plot_locus_browser(
  "Syt1", variant = "compact",
  loops_df = loops_df, rnaseq_df = rnaseq_df,
  cpg_islands_gr = cpg_islands_gr,
  cpg_shores_gr = cpg_shores_gr,
  cpg_shelves_gr = cpg_shelves_gr,
  peak_data = peak_data,
  precomputed_bw = syt1_bw,
  return_tracks = TRUE
)

# Capture Gviz plotTracks output as a grid grob
panel_a_grob <- NULL
pdf(NULL)
tryCatch({
  grid.newpage()
  plotTracks(
    track_info$tracks,
    chromosome = track_info$chr,
    from = track_info$from,
    to = track_info$to,
    sizes = track_info$sizes,
    main = "Syt1 Locus -- BAP1-KO Multi-Omic View",
    cex.main = 1.0,
    title.width = 1.6,
    col.axis = "black",
    fontsize = 12
  )
  panel_a_grob <- grid.grab()
}, finally = dev.off())

# ---- Panel B: Coordinated mC/hmC scatter -----------------------------------
cat("  Building Panel B (coordinated scatter)...\n")
coord_path <- file.path(TABLES_DIR, "coordinated_changes.tsv")
stopifnot(file.exists(coord_path))
coord_df <- read.table(coord_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

panel_b <- ggplot(coord_df, aes(x = mc_diff * 100, y = hmc_diff * 100)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(aes(color = as.character(coordinated_pattern)), alpha = 0.3, size = 1) +
  scale_color_manual(values = c("TRUE" = "#D7191C", "FALSE" = "grey60"), guide = "none") +
  geom_point(data = . %>% dplyr::filter(gene %in% KEY_GENES),
             size = 2.5, shape = 21, fill = "gold", color = "black", stroke = 0.8) +
  ggrepel::geom_text_repel(
    data = . %>% dplyr::filter(gene %in% KEY_GENES),
    aes(label = gene), fontface = "italic", size = 3,
    max.overlaps = 20, seed = 42
  ) +
  labs(x = "5mC change (%)", y = "5hmC change (%)",
       title = "Coordinated mC/hmC Changes") +
  theme_biomodal()

# ---- Panel C: Cross-gene summary heatmap -----------------------------------
cat("  Building Panel C (cross-gene heatmap)...\n")

# Merge multi-omic metrics for KEY_GENES
summary_df <- data.frame(gene = KEY_GENES, stringsAsFactors = FALSE) %>%
  dplyr::left_join(
    mc_dmr %>%
      dplyr::filter(gene %in% KEY_GENES) %>%
      dplyr::group_by(gene) %>%
      dplyr::slice_min(dmr_qvalue, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::select(gene, mC_diff = mod_difference),
    by = "gene"
  ) %>%
  dplyr::left_join(
    hmc_dmr %>%
      dplyr::filter(gene %in% KEY_GENES) %>%
      dplyr::group_by(gene) %>%
      dplyr::slice_min(dmr_qvalue, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::select(gene, hmC_diff = mod_difference),
    by = "gene"
  )

# Add RNA-seq log2FC
if (!is.null(rnaseq_df)) {
  summary_df <- summary_df %>%
    dplyr::left_join(
      rnaseq_df %>%
        dplyr::select(gene = ensembl_gene_id, log2FC = log2FoldChange),
      by = "gene"
    )
} else {
  summary_df$log2FC <- NA_real_
}

# Add Hi-C loop counts per gene locus
loop_counts <- sapply(KEY_GENES, function(g) {
  rg <- get_gene_region(g, extend_bp = EXTEND_BP)
  ch <- as.character(seqnames(rg)[1])
  rs <- start(rg); re <- end(rg)
  in_region <- loops_df %>%
    dplyr::filter(
      (chr1 == ch & start1 >= rs & end1 <= re) |
      (chr2 == ch & start2 >= rs & end2 <= re)
    )
  lost <- sum(in_region$direction == "down_in_mutant")
  gained <- sum(in_region$direction == "up_in_mutant")
  c(lost = lost, gained = gained)
})
summary_df$loops_lost   <- loop_counts["lost", ]
summary_df$loops_gained <- loop_counts["gained", ]

# Scale numeric columns to z-scores for heatmap
heat_cols <- c("mC_diff", "hmC_diff", "log2FC", "loops_lost", "loops_gained")
heat_df <- summary_df %>%
  dplyr::mutate(dplyr::across(
    dplyr::all_of(heat_cols),
    ~ as.numeric(scale(.x))
  ))

# Pivot to long format
heat_long <- heat_df %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(heat_cols),
    names_to = "metric",
    values_to = "z_score"
  ) %>%
  dplyr::mutate(
    metric = factor(metric,
      levels = heat_cols,
      labels = c("5mC", "5hmC", "RNA log2FC", "Loops Lost", "Loops Gained")
    ),
    gene = factor(gene, levels = rev(KEY_GENES))
  )

panel_c <- ggplot(heat_long, aes(x = metric, y = gene, fill = z_score)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, na.value = "grey90",
                       name = "Z-score") +
  labs(x = NULL, y = NULL, title = "Multi-Omic Summary") +
  theme_biomodal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(face = "italic", size = 10),
    panel.grid = element_blank()
  )

# ---- Assembly --------------------------------------------------------------
cat("  Assembling composite figure...\n")

composite_assembled <- FALSE

if (!is.null(panel_a_grob)) {
  panel_a_wrapped <- patchwork::wrap_elements(panel_a_grob)
  composite <- panel_a_wrapped / (panel_b | panel_c) +
    patchwork::plot_layout(heights = c(2, 1)) +
    patchwork::plot_annotation(tag_levels = "A")
  save_multiformat_ggplot(
    composite,
    file.path(SECTION_OUTPUT_DIR, "composite_syt1_panel"),
    width = 14, height = 16
  )
  composite_assembled <- TRUE
  cat("  Saved: composite_syt1_panel/{pdf,svg,jpg}\n")
}

# Fallback: save panels separately for manual Illustrator assembly
if (!composite_assembled) {
  cat("  Saving panels separately for manual assembly...\n")
  save_multiformat_ggplot(
    panel_b, file.path(SECTION_OUTPUT_DIR, "composite_panel_b"),
    width = 7, height = 6
  )
  save_multiformat_ggplot(
    panel_c, file.path(SECTION_OUTPUT_DIR, "composite_panel_c"),
    width = 7, height = 6
  )
  cat("  Saved: composite_panel_b/{pdf,svg,jpg}, composite_panel_c/{pdf,svg,jpg}\n")
  cat("  Panel A (Syt1 browser): use Syt1_locus_compact/ output for assembly\n")
}
} # end if (run_composite)

cat("\n================================================================================\n")
cat("SECTION 46 COMPLETE\n")
cat(sprintf("Output directory: %s\n", SECTION_OUTPUT_DIR))
cat("================================================================================\n")
