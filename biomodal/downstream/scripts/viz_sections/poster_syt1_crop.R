# biomodal/downstream/scripts/viz_sections/poster_syt1_crop.R
# Poster-cropped Syt1 locus: gene model + H2AK119ub (ctrl/mut) + Hi-C loop arcs
# Generates two color variants: original deep shades + loop-script-aligned brighter palette
# Run from downstream/ directory

source("scripts/viz_sections/_shared_config.R")

cat("================================================================================\n")
cat("POSTER: Syt1 LOCUS CROP (gene + K119ub + loops)\n")
cat("================================================================================\n\n")

suppressPackageStartupMessages({
  library(Gviz)
  library(GenomicInteractions)
  library(grid)
})

if (!hasMethod("seqinfo", "GenomicInteractions")) {
  setMethod("seqinfo", "GenomicInteractions", function(x) seqinfo(regions(x)))
}

# =============================================================================
# CONFIGURATION
# =============================================================================

REPO_ROOT   <- normalizePath(file.path(BASE_DIR, "../.."))
OUT_DIR_V1 <- file.path(REPO_ROOT, "poster/figures/syt1_locus_original_colors")
OUT_DIR_V2 <- file.path(REPO_ROOT, "poster/figures/syt1_locus_aligned_colors")
dir.create(OUT_DIR_V1, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR_V2, recursive = TRUE, showWarnings = FALSE)

EXTEND_BP <- 50000

HISTONE_BW_DIR <- "/Users/zakiralibhai/sdsc/bigwigs"
H2AK119UB_BW <- list(
  ctrl = file.path(HISTONE_BW_DIR, "H2AK119ubCtrl.bw"),
  mut  = file.path(HISTONE_BW_DIR, "H2AK119ubMut.bw")
)

COL_PALETTES <- list(
  original = list(
    ctrl       = "#2166AC",
    mut        = "#B2182B",
    hic_lost   = "#B2182B",
    hic_gained = "#2166AC",
    gene       = "#08306B"
  ),
  aligned = list(
    ctrl       = "#4575b4",
    mut        = "#d73027",
    hic_lost   = "#d73027",
    hic_gained = "#4575b4",
    gene       = "#08306B"
  )
)

COL <- COL_PALETTES$original

# =============================================================================
# HELPERS (adapted from section_46)
# =============================================================================

get_gene_region <- function(gene_symbol, extend_bp = EXTEND_BP) {
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  entrez <- AnnotationDbi::mapIds(
    org.Mm.eg.db, keys = gene_symbol,
    column = "ENTREZID", keytype = "SYMBOL")
  stopifnot(!is.na(entrez))

  all_genes <- suppressMessages(genes(txdb))
  gene_gr <- all_genes[names(all_genes) == entrez]
  stopifnot(length(gene_gr) == 1)

  extended <- gene_gr
  start(extended) <- max(1, start(gene_gr) - extend_bp)
  end(extended)   <- end(gene_gr) + extend_bp
  mcols(extended)$symbol <- gene_symbol
  return(extended)
}

import_bigwig <- function(bw_path, region_gr, bin_size = NULL) {
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
    genome(a1) <- "mm10"; genome(a2) <- "mm10"
    gi <- GenomicInteractions(a1, a2)
    gi$name <- lost$loop_id
    tracks$lost <- InteractionTrack(gi, name = "Lost Loops", chromosome = chr)
    displayPars(tracks$lost) <- list(
      col.interactions = COL$hic_lost,
      col.anchors.fill = COL$hic_lost,
      col.anchors.line = COL$hic_lost,
      interaction.dimension = "height",
      interaction.measure = "counts",
      plot.anchors = FALSE, plot.trans = FALSE,
      plot.outside = TRUE,
      col.outside = adjustcolor(COL$hic_lost, alpha.f = 0.7),
      anchor.height = 0.1,
      cex.title = 1.2
    )
  }

  gained <- in_region %>% dplyr::filter(direction == "up_in_mutant")
  if (nrow(gained) > 0) {
    a1 <- GRanges(gained$chr1, IRanges(gained$start1, gained$end1))
    a2 <- GRanges(gained$chr2, IRanges(gained$start2, gained$end2))
    genome(a1) <- "mm10"; genome(a2) <- "mm10"
    gi <- GenomicInteractions(a1, a2)
    gi$name <- gained$loop_id
    tracks$gained <- InteractionTrack(gi, name = "Gained Loops", chromosome = chr)
    displayPars(tracks$gained) <- list(
      col.interactions = COL$hic_gained,
      col.anchors.fill = COL$hic_gained,
      col.anchors.line = COL$hic_gained,
      interaction.dimension = "height",
      interaction.measure = "counts",
      plot.anchors = FALSE, plot.trans = FALSE,
      plot.outside = TRUE,
      col.outside = adjustcolor(COL$hic_gained, alpha.f = 0.7),
      anchor.height = 0.1,
      cex.title = 1.2
    )
  }

  return(tracks)
}

# =============================================================================
# DATA LOADING
# =============================================================================

cat("Loading Hi-C loops...\n")
loops_df <- read.table(LOOP_FILES$late, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE)
cat(sprintf("  %d differential loops\n", nrow(loops_df)))

stopifnot(file.exists(H2AK119UB_BW$ctrl), file.exists(H2AK119UB_BW$mut))
cat("  H2AK119ub BigWig files verified\n")

# =============================================================================
# BUILD TRACKS
# =============================================================================

cat("\nBuilding Syt1 poster tracks...\n")
region_gr <- get_gene_region("Syt1", extend_bp = EXTEND_BP)
chr        <- as.character(seqnames(region_gr)[1])
view_start <- start(region_gr)
view_end   <- end(region_gr)
cat(sprintf("  Region: %s:%s-%s (%s kb)\n",
            chr, format(view_start, big.mark = ","),
            format(view_end, big.mark = ","),
            round(width(region_gr) / 1000)))

cat("  Importing H2AK119ub BigWigs...\n")
k119_ctrl_gr <- import_bigwig(H2AK119UB_BW$ctrl, region_gr)
k119_mut_gr  <- import_bigwig(H2AK119UB_BW$mut, region_gr)
k119_ylim <- c(0, max(c(k119_ctrl_gr$score, k119_mut_gr$score), na.rm = TRUE) * 1.05)
if (k119_ylim[2] == 0) k119_ylim[2] <- 1

track_list <- list()
sizes <- c()

track_list$axis <- GenomeAxisTrack(name = "", scale = 0.2, labelPos = "below")
sizes <- c(sizes, 0.5)

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
track_list$gene <- GeneRegionTrack(
  txdb, genome = "mm10", chromosome = chr,
  name = "Syt1",
  transcriptAnnotation = "symbol",
  stacking = "dense",
  col = COL$gene, fill = COL$gene,
  fontcolor.group = COL$gene,
  background.title = COL$gene,
  col.title = "white", cex.title = 1.2
)
sizes <- c(sizes, 1.2)

track_list$k119ub_ctrl <- DataTrack(
  range = k119_ctrl_gr, genome = "mm10", type = "histogram",
  name = "H2AK119ub\nControl",
  fill.histogram = COL$ctrl, col.histogram = COL$ctrl,
  ylim = k119_ylim,
  background.title = COL$ctrl,
  col.title = "white", cex.title = 1.2
)
sizes <- c(sizes, 1.0)

track_list$k119ub_mut <- DataTrack(
  range = k119_mut_gr, genome = "mm10", type = "histogram",
  name = "H2AK119ub\nMutant",
  fill.histogram = COL$mut, col.histogram = COL$mut,
  ylim = k119_ylim,
  background.title = COL$mut,
  col.title = "white", cex.title = 1.2
)
sizes <- c(sizes, 1.0)

cat("  Building Hi-C loop arcs...\n")
hic_tracks <- build_hic_tracks(loops_df, region_gr)
if (!is.null(hic_tracks)) {
  if (!is.null(hic_tracks$lost)) {
    track_list$hic_lost <- hic_tracks$lost
    sizes <- c(sizes, 1.5)
  }
  if (!is.null(hic_tracks$gained)) {
    track_list$hic_gained <- hic_tracks$gained
    sizes <- c(sizes, 3.0)
  }
}

cat(sprintf("  %d tracks built\n", length(track_list)))

# =============================================================================
# RENDER — generate both color variants
# =============================================================================

fig_width  <- 12
fig_height <- sum(sizes) * 1.0 + 2.0

render_plot <- function() {
  plotTracks(
    track_list,
    chromosome = chr,
    from = view_start,
    to = view_end,
    sizes = sizes,
    main = "Syt1 Locus",
    cex.main = 1.4,
    title.width = 1.6,
    col.axis = "black",
    fontsize = 18
  )
}

save_variant <- function(out_dir, file_stem) {
  file_prefix <- file.path(out_dir, file_stem)

  pdf(NULL, width = fig_width, height = fig_height)
  tryCatch({
    render_plot()
    plot_grob <- grid.grab()
  }, finally = dev.off())

  pdf(paste0(file_prefix, ".pdf"), width = fig_width, height = fig_height)
  tryCatch({ grid.newpage(); grid.draw(plot_grob) }, finally = dev.off())

  svglite::svglite(paste0(file_prefix, ".svg"), width = fig_width, height = fig_height)
  tryCatch({ grid.newpage(); grid.draw(plot_grob) }, finally = dev.off())

  jpeg(paste0(file_prefix, ".jpg"),
       width = fig_width * 300, height = fig_height * 300,
       res = 300, quality = 95)
  tryCatch({ grid.newpage(); grid.draw(plot_grob) }, finally = dev.off())

  cat(sprintf("  Saved: %s/{pdf,svg,jpg}\n", out_dir))
}

recolor_tracks <- function(palette) {
  displayPars(track_list$k119ub_ctrl) <- list(
    fill.histogram = palette$ctrl, col.histogram = palette$ctrl,
    background.title = palette$ctrl
  )
  displayPars(track_list$k119ub_mut) <- list(
    fill.histogram = palette$mut, col.histogram = palette$mut,
    background.title = palette$mut
  )
  if (!is.null(track_list$hic_lost)) {
    displayPars(track_list$hic_lost) <- list(
      col.interactions = palette$hic_lost,
      col.anchors.fill = palette$hic_lost,
      col.anchors.line = palette$hic_lost,
      col.outside = adjustcolor(palette$hic_lost, alpha.f = 0.7)
    )
  }
  if (!is.null(track_list$hic_gained)) {
    displayPars(track_list$hic_gained) <- list(
      col.interactions = palette$hic_gained,
      col.anchors.fill = palette$hic_gained,
      col.anchors.line = palette$hic_gained,
      col.outside = adjustcolor(palette$hic_gained, alpha.f = 0.7)
    )
  }
  track_list <<- track_list
}

cat("\nRendering variant 1 (original colors)...\n")
recolor_tracks(COL_PALETTES$original)
save_variant(OUT_DIR_V1, "Syt1_locus_poster_v2")

cat("Rendering variant 2 (aligned colors)...\n")
recolor_tracks(COL_PALETTES$aligned)
save_variant(OUT_DIR_V2, "Syt1_locus_poster_v2_aligned")

cat("\nDone — compare both variants to choose colors.\n")
