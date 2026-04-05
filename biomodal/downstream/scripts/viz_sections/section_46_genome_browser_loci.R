# biomodal/downstream/scripts/viz_sections/section_46_genome_browser_loci.R
# Multi-omic genome browser views at key gene loci using Gviz
# Renders stacked tracks: 5mC/5hmC (condition-averaged BigWigs), difference,
# ChIP/epigenomic peak annotations, CpG density, Hi-C loop arcs
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

# --- RNA-seq BigWig paths (set paths when available) -------------------------
# To enable RNA-seq coverage tracks, set these to character vectors of BigWig
# file paths (one per sample). When NULL, a gene-level text annotation with
# log2FC/padj from DESeq2 results is shown instead.
RNASEQ_BW_CTRL <- NULL
RNASEQ_BW_MUT  <- NULL

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
  gene        = "#08306B"
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


#' Average BigWig signal across samples within a genomic region
#' Returns a binned GRanges with mean score per bin
average_bigwigs_in_region <- function(bw_paths, region_gr, bin_size = NULL) {
  chr <- as.character(seqnames(region_gr)[1])
  region_start <- start(region_gr)
  region_end   <- end(region_gr)
  region_width <- width(region_gr)

  # Auto-select bin size: ~1000 bins across the view

  if (is.null(bin_size)) {
    bin_size <- max(50, round(region_width / 1000))
  }

  # Create bins across the region
  bin_starts <- seq(region_start, region_end - bin_size, by = bin_size)
  bins_gr <- GRanges(seqnames = chr,
                     ranges = IRanges(start = bin_starts, width = bin_size))

  # Compute mean methylation per bin for each sample
  # Key: BigWig entries are at individual CpG sites (1-2bp wide), so we must

  # compute mean(scores at CpGs) NOT mean(signal over bin width)
  valid_paths <- bw_paths[file.exists(bw_paths)]
  if (length(valid_paths) == 0) {
    stop("No BigWig files found. Paths checked:\n",
         paste(" ", bw_paths, collapse = "\n"))
  }

  score_matrix <- sapply(valid_paths, function(path) {
    bw <- rtracklayer::import.bw(path, which = region_gr)
    if (length(bw) == 0) return(rep(NA_real_, length(bins_gr)))

    # Use weighted coverage / unweighted coverage = mean score at covered sites
    weighted_cov <- coverage(bw, weight = "score")
    count_cov    <- coverage(bw)
    if (!chr %in% names(weighted_cov)) return(rep(NA_real_, length(bins_gr)))

    v_weighted <- Views(weighted_cov[[chr]], ranges(bins_gr))
    v_count    <- Views(count_cov[[chr]], ranges(bins_gr))

    sum_scores <- as.numeric(viewSums(v_weighted))
    num_cpgs   <- as.numeric(viewSums(v_count))

    ifelse(num_cpgs > 0, sum_scores / num_cpgs, NA_real_)
  })

  # Average across samples
  if (is.matrix(score_matrix)) {
    bins_gr$score <- rowMeans(score_matrix, na.rm = TRUE)
  } else {
    bins_gr$score <- score_matrix
  }
  # Replace any remaining NAs with 0 for display
  bins_gr$score[is.na(bins_gr$score)] <- 0
  return(bins_gr)
}


#' Load a BED file as GRanges, adding chr prefix if needed
load_peak_as_granges <- function(path, name = "peak") {
  if (!file.exists(path)) {
    warning("Peak file not found: ", path)
    return(GRanges())
  }
  df <- read.table(path, header = FALSE, sep = "\t",
                   stringsAsFactors = FALSE, fill = TRUE)
  gr <- GRanges(seqnames = df$V1,
                ranges = IRanges(start = df$V2, end = df$V3))

  # Ensure chr prefix
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
  # Add chr prefix
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

  # Filter: at least one anchor overlaps the view region
  in_region <- loops_df %>%
    dplyr::filter(
      (chr1 == chr & ((start1 >= rstart & start1 <= rend) |
                      (end1 >= rstart & end1 <= rend))) |
      (chr2 == chr & ((start2 >= rstart & start2 <= rend) |
                      (end2 >= rstart & end2 <= rend)))
    )

  if (nrow(in_region) == 0) return(NULL)

  tracks <- list()

  # Lost loops (down_in_mutant)
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

  # Gained loops (up_in_mutant)
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
plot_locus_browser <- function(gene_symbol, variant = "full",
                               loops_df, rnaseq_df,
                               cpg_islands_gr, cpg_shores_gr,
                               cpg_shelves_gr) {

  cat(sprintf("  Generating %s browser view for %s...\n", variant, gene_symbol))

  # --- 1. Get gene region ---
  region_gr <- tryCatch(
    get_gene_region(gene_symbol, extend_bp = EXTEND_BP),
    error = function(e) {
      warning(sprintf("Could not find gene %s in TxDb: %s", gene_symbol, e$message))
      return(NULL)
    }
  )
  if (is.null(region_gr)) return(invisible(NULL))

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

  # --- 3. Average methylation BigWigs ---
  cat("    Averaging 5mC BigWigs...\n")
  mc_ctrl_gr  <- average_bigwigs_in_region(MC_BW_CTRL, region_gr)
  mc_mut_gr   <- average_bigwigs_in_region(MC_BW_MUT, region_gr)

  cat("    Averaging 5hmC BigWigs...\n")
  hmc_ctrl_gr <- average_bigwigs_in_region(HMC_BW_CTRL, region_gr)
  hmc_mut_gr  <- average_bigwigs_in_region(HMC_BW_MUT, region_gr)

  # Compute difference tracks (mutant - control)
  mc_diff_gr <- mc_ctrl_gr
  mc_diff_gr$score <- mc_mut_gr$score - mc_ctrl_gr$score

  hmc_diff_gr <- hmc_ctrl_gr
  hmc_diff_gr$score <- hmc_mut_gr$score - hmc_ctrl_gr$score

  # Scale to percentage (BigWig values are 0-1 fractions)
  mc_ctrl_gr$score  <- mc_ctrl_gr$score * 100
  mc_mut_gr$score   <- mc_mut_gr$score * 100
  mc_diff_gr$score  <- mc_diff_gr$score * 100
  hmc_ctrl_gr$score <- hmc_ctrl_gr$score * 100
  hmc_mut_gr$score  <- hmc_mut_gr$score * 100
  hmc_diff_gr$score <- hmc_diff_gr$score * 100

  # --- 4. Build Gviz tracks ---
  track_list <- list()
  sizes <- c()

  # Track 0: Genome axis
  track_list$axis <- GenomeAxisTrack(name = "Position")
  sizes <- c(sizes, 0.5)

  # Track 1: Gene model (collapse to single meta-transcript per gene)
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  track_list$gene <- GeneRegionTrack(
    txdb, genome = "mm10", chromosome = chr,
    name = paste0(gene_symbol, rnaseq_label),
    transcriptAnnotation = "symbol",
    collapseTranscripts = "meta",
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

  # 5mC Difference
  track_list$mc_diff <- DataTrack(
    range = mc_diff_gr, genome = "mm10", type = "histogram",
    name = "5mC%\nDifference",
    fill.histogram = TRACK_COL$diff_hyper,
    col.histogram = TRACK_COL$diff_hyper,
    ylim = diff_ylim,
    baseline = 0, col.baseline = "black", lwd.baseline = 0.5,
    background.title = "#8C510A",
    col.title = "white", cex.title = 0.7
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

    # 5hmC Difference (full version only)
    track_list$hmc_diff <- DataTrack(
      range = hmc_diff_gr, genome = "mm10", type = "histogram",
      name = "5hmC%\nDifference",
      fill.histogram = TRACK_COL$diff_hypo,
      col.histogram = TRACK_COL$diff_hypo,
      ylim = diff_ylim,
      baseline = 0, col.baseline = "black", lwd.baseline = 0.5,
      background.title = "#01665E",
      col.title = "white", cex.title = 0.7
    )
    sizes <- c(sizes, 1.5)
  }

  # --- RNA-seq BigWig tracks (modular: only if paths are set) ---
  if (!is.null(RNASEQ_BW_CTRL)) {
    cat("    Averaging RNA-seq BigWigs...\n")
    rnaseq_ctrl_gr <- average_bigwigs_in_region(RNASEQ_BW_CTRL, region_gr)
    rnaseq_mut_gr  <- average_bigwigs_in_region(RNASEQ_BW_MUT, region_gr)

    track_list$rnaseq_ctrl <- DataTrack(
      range = rnaseq_ctrl_gr, genome = "mm10", type = "histogram",
      name = "RNA-seq\nControl",
      fill.histogram = TRACK_COL$ctrl,
      col.histogram = TRACK_COL$ctrl,
      background.title = TRACK_COL$ctrl,
      col.title = "white", cex.title = 0.7
    )
    sizes <- c(sizes, 1.2)

    track_list$rnaseq_mut <- DataTrack(
      range = rnaseq_mut_gr, genome = "mm10", type = "histogram",
      name = "RNA-seq\nMutant",
      fill.histogram = TRACK_COL$mut,
      col.histogram = TRACK_COL$mut,
      background.title = TRACK_COL$mut,
      col.title = "white", cex.title = 0.7
    )
    sizes <- c(sizes, 1.2)
  }

  # --- Chromatin / epigenomic peak tracks ---

  # H2AK119ub (always shown -- central to BAP1 mechanism)
  k119ub_ctrl_gr <- load_peak_as_granges(K119UB_FILES$ctrl, "K119ub_Ctrl")
  k119ub_mut_gr  <- load_peak_as_granges(K119UB_FILES$mut, "K119ub_Mut")
  k119_combined  <- c(k119ub_ctrl_gr, k119ub_mut_gr)
  if (length(k119_combined) > 0) {
    k119_combined$feature <- c(
      rep("Control", length(k119ub_ctrl_gr)),
      rep("Mutant", length(k119ub_mut_gr))
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

  if (variant == "full") {
    # H3K27me3 peaks
    k27me3_gr <- load_peak_as_granges(CHIP_PEAK_FILES$h3k27me3, "H3K27me3")
    if (length(k27me3_gr) > 0) {
      track_list$k27me3 <- AnnotationTrack(
        k27me3_gr, genome = "mm10", chromosome = chr,
        name = "H3K27me3",
        fill = TRACK_COL$k27me3, col = TRACK_COL$k27me3,
        stacking = "dense",
        background.title = TRACK_COL$k27me3,
        col.title = "white", cex.title = 0.7
      )
      sizes <- c(sizes, 0.4)
    }

    # H3K27ac condition-specific
    k27ac_ctrl_gr <- load_peak_as_granges(H3K27AC_FILES$ctrl, "K27ac_Ctrl")
    k27ac_mut_gr  <- load_peak_as_granges(H3K27AC_FILES$mut, "K27ac_Mut")
    k27ac_combined <- c(k27ac_ctrl_gr, k27ac_mut_gr)
    if (length(k27ac_combined) > 0) {
      k27ac_combined$feature <- c(
        rep("Control", length(k27ac_ctrl_gr)),
        rep("Mutant", length(k27ac_mut_gr))
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

    # ATAC-seq differential
    atac_up_gr   <- load_peak_as_granges(ATAC_FILES$up, "ATAC_Up")
    atac_down_gr <- load_peak_as_granges(ATAC_FILES$down, "ATAC_Down")
    atac_combined <- c(atac_up_gr, atac_down_gr)
    if (length(atac_combined) > 0) {
      atac_combined$feature <- c(
        rep("Up", length(atac_up_gr)),
        rep("Down", length(atac_down_gr))
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

    # MeCP2 differential
    mecp2_up_gr   <- load_peak_as_granges(MECP2_FILES$up, "MeCP2_Up")
    mecp2_down_gr <- load_peak_as_granges(MECP2_FILES$down, "MeCP2_Down")
    mecp2_combined <- c(mecp2_up_gr, mecp2_down_gr)
    if (length(mecp2_combined) > 0) {
      mecp2_combined$feature <- c(
        rep("Up", length(mecp2_up_gr)),
        rep("Down", length(mecp2_down_gr))
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

  # CTCF peaks (both variants)
  ctcf_gr <- load_peak_as_granges(CHIP_PEAK_FILES$ctcf, "CTCF")
  if (length(ctcf_gr) > 0) {
    track_list$ctcf <- AnnotationTrack(
      ctcf_gr, genome = "mm10", chromosome = chr,
      name = "CTCF",
      fill = TRACK_COL$ctcf, col = TRACK_COL$ctcf,
      stacking = "dense",
      background.title = TRACK_COL$ctcf,
      col.title = "white", cex.title = 0.7
    )
    sizes <- c(sizes, 0.4)
  }

  # --- CpG annotation track ---
  # Subset to view region for efficiency
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

  # --- 5. Render the plot ---
  suffix <- ifelse(variant == "compact", "_compact", "")
  out_name <- sprintf("%s_locus%s", gene_symbol, suffix)
  out_dir  <- file.path(SECTION_OUTPUT_DIR, out_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  file_prefix <- file.path(out_dir, out_name)

  fig_height <- sum(sizes) * 0.7 + 1  # scale factor + margins
  fig_width  <- 14

  # Closure captures all variables from this scope, avoiding eval() issues
  render_plot <- function() {
    plotTracks(
      track_list,
      chromosome = chr,
      from = view_start,
      to = view_end,
      sizes = sizes,
      main = sprintf("%s Locus \u2014 BAP1-KO Multi-Omic View (%s)", gene_symbol, variant),
      cex.main = 1.0,
      title.width = 1.2,
      col.axis = "black",
      fontsize = 10
    )
  }

  # PDF
  pdf(paste0(file_prefix, ".pdf"), width = fig_width, height = fig_height)
  tryCatch(render_plot(), finally = dev.off())

  # SVG
  svglite::svglite(paste0(file_prefix, ".svg"), width = fig_width, height = fig_height)
  tryCatch(render_plot(), finally = dev.off())

  # JPEG
  jpeg(paste0(file_prefix, ".jpg"),
       width = fig_width * 300, height = fig_height * 300,
       res = 300, quality = 95)
  tryCatch(render_plot(), finally = dev.off())

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

# =============================================================================
# GENERATE BROWSER VIEWS
# =============================================================================

cat("\n--- Generating genome browser views ---\n\n")

for (gene in KEY_GENES) {
  cat(sprintf("Processing %s...\n", gene))

  tryCatch({
    # Full version
    plot_locus_browser(
      gene, variant = "full",
      loops_df = loops_df, rnaseq_df = rnaseq_df,
      cpg_islands_gr = cpg_islands_gr,
      cpg_shores_gr = cpg_shores_gr,
      cpg_shelves_gr = cpg_shelves_gr
    )

    # Compact version
    plot_locus_browser(
      gene, variant = "compact",
      loops_df = loops_df, rnaseq_df = rnaseq_df,
      cpg_islands_gr = cpg_islands_gr,
      cpg_shores_gr = cpg_shores_gr,
      cpg_shelves_gr = cpg_shelves_gr
    )
  }, error = function(e) {
    cat(sprintf("  ERROR for %s: %s\n", gene, e$message))
  })

  cat("\n")
}

cat("================================================================================\n")
cat("SECTION 46 COMPLETE\n")
cat(sprintf("Output directory: %s\n", SECTION_OUTPUT_DIR))
cat("================================================================================\n")
