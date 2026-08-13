# poster/pub_browser/pub_browser.R
# Publication-quality genome browser locus views using karyoploteR
#
# Renders multi-track locus figures directly to SVG/PDF/PNG/JPEG in the style of
# manually-polished Illustrator figures: filled area curves with same color per
# mark, rotated colored mark labels, gene models with strand arrows, italic gene
# symbols, scale-bar coordinate marker, and optional gray highlight rectangles.
#
# Usage (CLI):
#   Rscript pub_browser.R \
#     --gene Syt1 --genome mm10 \
#     --mark 'H2AK119ub:#2171B5:ctrl.bw:mut.bw' \
#     --mark 'H3K27me3:#CB181D:ctrl.bw:mut.bw' \
#     --mark 'H3K27ac:#41AB5D:ctrl.bw:mut.bw' \
#     --labels 'control,Bap1-KO' \
#     --output out/syt1_browser
#
# Usage (YAML config):
#   Rscript pub_browser.R --config locus_config.yaml
#
# Mark spec format: 'name:color:ctrl_path:mut_path'
#   - name:  display name (rotated label text)
#   - color: hex color, used for area fill, label, and scale annotation
#   - ctrl:  path to control BigWig (or comma-sep list for replicates)
#   - mut:   path to mutant BigWig (or comma-sep list for replicates)
#
# Each --mark flag adds one mark group (ctrl on top, mut below, shared y-axis).
# Marks render in the order specified. Use YAML --config for complex setups,
# multi-replicate averaging, or batch processing.

# =============================================================================
# LIBRARIES
# =============================================================================

suppressPackageStartupMessages({
  library(karyoploteR)
  library(rtracklayer)
  library(GenomicRanges)
  library(GenomicFeatures)
  library(dplyr)
  library(svglite)
})

load_optional <- function(pkg) {
  suppressPackageStartupMessages(suppressWarnings(
    requireNamespace(pkg, quietly = TRUE)
  ))
}

# =============================================================================
# CLI PARSING
# =============================================================================

parse_int <- function(val, flag) {
  out <- suppressWarnings(as.integer(val))
  if (is.na(out)) stop(sprintf("%s requires an integer, got '%s'", flag, val))
  out
}
parse_num <- function(val, flag) {
  out <- suppressWarnings(as.numeric(val))
  if (is.na(out)) stop(sprintf("%s requires a number, got '%s'", flag, val))
  out
}

parse_cli_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  cfg <- list(
    gene             = NULL,
    region           = NULL,
    genome           = "mm10",
    marks            = list(),
    config_file      = NULL,
    labels           = NULL,
    label_list       = list(),
    output           = "out/locus",
    highlights       = list(),
    highlight_labels = list(),
    hic_loops        = NULL,
    extend_bp        = 50000L,
    scale_bar_kb     = NULL,
    fig_width        = 7.0,
    track_height     = 0.40,
    txdb_pkg         = "TxDb.Mmusculus.UCSC.mm10.knownGene",
    orgdb_pkg        = "org.Mm.eg.db",
    gene_symbol      = NULL,
    focal_genes      = list(),
    show_all_genes   = FALSE,
    bin_size         = NULL,
    genotype_italic  = FALSE,
    help             = FALSE
  )

  i <- 1L
  while (i <= length(args)) {
    arg <- args[i]
    val <- if (i < length(args)) args[i + 1L] else NULL

    switch(arg,
      "--help" = , "-h" = { cfg$help <- TRUE; i <- i + 1L },
      "--gene"            = { cfg$gene <- val; i <- i + 2L },
      "--region"          = { cfg$region <- val; i <- i + 2L },
      "--genome"          = { cfg$genome <- val; i <- i + 2L },
      "--mark"            = { cfg$marks[[length(cfg$marks) + 1L]] <- val; i <- i + 2L },
      "--config"          = { cfg$config_file <- val; i <- i + 2L },
      "--labels"          = { cfg$labels <- val; i <- i + 2L },
      "--label"           = { cfg$label_list[[length(cfg$label_list) + 1L]] <- val; i <- i + 2L },
      "--output"          = { cfg$output <- val; i <- i + 2L },
      "--highlight"       = { cfg$highlights[[length(cfg$highlights) + 1L]] <- val; i <- i + 2L },
      "--highlight-label" = { cfg$highlight_labels[[length(cfg$highlight_labels) + 1L]] <- val; i <- i + 2L },
      "--hic-loops"       = { cfg$hic_loops <- val; i <- i + 2L },
      "--extend"          = { cfg$extend_bp <- parse_int(val, "--extend"); i <- i + 2L },
      "--scale-bar"       = { cfg$scale_bar_kb <- parse_num(val, "--scale-bar"); i <- i + 2L },
      "--width"           = { cfg$fig_width <- parse_num(val, "--width"); i <- i + 2L },
      "--track-height"    = { cfg$track_height <- parse_num(val, "--track-height"); i <- i + 2L },
      "--txdb"            = { cfg$txdb_pkg <- val; i <- i + 2L },
      "--orgdb"           = { cfg$orgdb_pkg <- val; i <- i + 2L },
      "--gene-symbol"     = { cfg$gene_symbol <- val; i <- i + 2L },
      "--focal-gene"      = { cfg$focal_genes[[length(cfg$focal_genes) + 1L]] <- val; i <- i + 2L },
      "--show-all-genes"  = { cfg$show_all_genes <- TRUE; i <- i + 1L },
      "--bin-size"        = { cfg$bin_size <- parse_int(val, "--bin-size"); i <- i + 2L },
      "--genotype-italic" = { cfg$genotype_italic <- TRUE; i <- i + 1L },
      {
        stop(sprintf("Unknown argument: %s\nUse --help for usage.", arg))
      }
    )
  }

  cfg
}

print_help <- function() {
  cat("
pub_browser.R -- Publication-quality genome browser locus tool (karyoploteR)

REQUIRED (one of):
  --gene <SYMBOL>            Gene symbol (resolved via TxDb)
  --region <chr:start-end>   Explicit genomic region
  --config <path.yaml>       YAML config file (overrides CLI defaults)

MARKS (at least one required):
  --mark 'name:color:ctrl.bw:mut.bw'
                             Display name, hex color, BigWig paths.
                             Repeatable. Multiple replicates allowed via
                             comma-separated paths in each slot (averaged).

CORE OPTIONS:
  --genome <build>           Genome assembly (default: mm10)
  --label <text>             Condition label (repeatable; pass exactly twice).
  --labels 'ctrl,mut'        Legacy shortcut: comma-sep condition labels.
                             Cannot contain commas in either label -- use
                             repeatable --label for genotypes with commas.
                             (default: 'control,mutant')
  --output <path-prefix>     Output path without extension
                             (default: out/locus)
  --extend <bp>              Flanking bp on each side of gene
                             (default: 50000)
  --scale-bar <kb>           Scale bar length in kb (auto if not set)
  --txdb <pkg>               TxDb package name
                             (default: TxDb.Mmusculus.UCSC.mm10.knownGene)
  --orgdb <pkg>              OrgDb package for symbol->entrez mapping
                             (default: org.Mm.eg.db)
  --bin-size <bp>            BigWig binning (auto: region_width/1000)

VISUAL OPTIONS:
  --highlight <chr:start-end>      Gray highlight box (repeatable)
  --highlight-label <text>         Label for highlight (repeatable, paired)
  --show-all-genes                 Show all overlapping genes (not just the
                                   focal gene) in the gene model panel
  --gene-symbol <SYMBOL>           Override displayed gene symbol
                                   (useful with --region)
  --genotype-italic                Render second condition label as italic

LAYOUT OPTIONS:
  --width <inches>           Figure width (default: 7.0)
  --track-height <inches>    Per-track height (default: 0.40)

OPTIONAL OVERLAYS:
  --hic-loops <path.tsv>     Hi-C loops TSV with columns:
                             chr1,start1,end1,chr2,start2,end2,direction
                             where direction in {up_in_mutant, down_in_mutant}

EXAMPLES:

  # Minimal three-mark view of Syt1
  Rscript pub_browser.R \\
    --gene Syt1 \\
    --mark 'H2AK119ub:#2171B5:H2AK119ubCtrl.bw:H2AK119ubMut.bw' \\
    --mark 'H3K27ac:#41AB5D:H3K27acCtrl.bw:H3K27acMut.bw' \\
    --mark 'H3K27me3:#CB181D:H3K27me3Ctrl.bw:H3K27me3Mut.bw' \\
    --labels 'control,Bap1-KO' \\
    --output out/syt1_pub

  # With highlight and italic genotype label
  Rscript pub_browser.R \\
    --gene Anxa2r1 \\
    --mark 'H2AK119ub:#2171B5:ctrl.bw:mut.bw' \\
    --label 'control' \\
    --label 'Bap1-KO' \\
    --highlight 'chr3:9505000-9510000' \\
    --highlight-label 'differential region' \\
    --output out/anxa2r1_pub

  # Region-based view with a focal gene resolved from symbol
  Rscript pub_browser.R \\
    --region 'chr3:9450000-9520000' \\
    --gene-symbol Anxa2r1 \\
    --mark 'H3K27me3:#CB181D:ctrl.bw:mut.bw' \\
    --output out/anxa2r1_extended

  # YAML config (recommended for complex setups)
  Rscript pub_browser.R --config configs/syt1_locus.yaml
")
}

parse_mark_spec <- function(spec) {
  parts <- strsplit(spec, ":", fixed = TRUE)[[1]]
  if (length(parts) != 4L) {
    stop(sprintf(
      "Bad --mark spec: '%s'\nExpected 'name:color:ctrl.bw:mut.bw'", spec))
  }
  list(
    name    = parts[1],
    color   = parts[2],
    ctrl    = strsplit(parts[3], ",", fixed = TRUE)[[1]],
    mut     = strsplit(parts[4], ",", fixed = TRUE)[[1]],
    average = TRUE,
    sparse  = FALSE,
    diff    = FALSE
  )
}

parse_region_string <- function(s) {
  m <- regmatches(s, regexec("^([^:]+):([0-9,]+)-([0-9,]+)$", s))[[1]]
  if (length(m) != 4L) {
    stop(sprintf("Bad region string: '%s'\nExpected 'chr:start-end'", s))
  }
  list(
    chr   = m[2],
    start = as.integer(gsub(",", "", m[3], fixed = TRUE)),
    end   = as.integer(gsub(",", "", m[4], fixed = TRUE))
  )
}

# =============================================================================
# CONFIG RESOLUTION (CLI + YAML merge)
# =============================================================================

resolve_config <- function(cli) {
  cfg <- cli

  if (!is.null(cli$config_file)) {
    if (!load_optional("yaml")) {
      stop("Package 'yaml' required for --config. Install with install.packages('yaml').")
    }
    yml <- yaml::read_yaml(cli$config_file)

    cfg$genome          <- if (cli$genome != "mm10") cli$genome else (yml$genome %||% cfg$genome)
    cfg$gene            <- cli$gene %||% yml$gene
    cfg$region          <- cli$region %||% yml$region
    cfg$gene_symbol     <- cli$gene_symbol %||% yml$gene_symbol
    cfg$extend_bp       <- if (cli$extend_bp != 50000L) cli$extend_bp else (yml$extend_bp %||% cfg$extend_bp)
    cfg$scale_bar_kb    <- cli$scale_bar_kb %||% yml$scale_bar_kb
    cfg$output          <- if (cli$output != "out/locus") cli$output else (yml$output %||% cfg$output)
    cfg$hic_loops       <- cli$hic_loops %||% yml$hic_loops
    cfg$fig_width       <- if (cli$fig_width != 7.0) cli$fig_width else (yml$width %||% cfg$fig_width)
    cfg$track_height    <- if (cli$track_height != 0.40) cli$track_height else (yml$track_height %||% cfg$track_height)
    cfg$txdb_pkg        <- if (cli$txdb_pkg != "TxDb.Mmusculus.UCSC.mm10.knownGene") cli$txdb_pkg else (yml$txdb %||% cfg$txdb_pkg)
    cfg$orgdb_pkg       <- if (cli$orgdb_pkg != "org.Mm.eg.db") cli$orgdb_pkg else (yml$orgdb %||% cfg$orgdb_pkg)
    cfg$show_all_genes  <- cli$show_all_genes || isTRUE(yml$show_all_genes)
    cfg$genotype_italic <- cli$genotype_italic || isTRUE(yml$genotype_italic)
    cfg$bin_size        <- cli$bin_size %||% yml$bin_size

    if (length(cli$focal_genes) == 0L && !is.null(yml$focal_genes)) {
      cfg$focal_genes <- as.list(as.character(yml$focal_genes))
    }

    if (length(cli$label_list) == 0L && !is.null(yml$labels)) {
      cfg$label_list <- as.list(as.character(yml$labels))
    }

    if (length(cli$marks) == 0L && !is.null(yml$marks)) {
      cfg$marks <- lapply(yml$marks, function(m) {
        list(
          name    = m$name,
          color   = m$color,
          ctrl    = as.character(m$ctrl),
          mut     = as.character(m$mut),
          average = m$average %||% TRUE,
          sparse  = isTRUE(m$sparse),
          diff    = isTRUE(m$diff)
        )
      })
    } else if (length(cli$marks) > 0L) {
      cfg$marks <- lapply(cli$marks, parse_mark_spec)
    }

    if (length(cli$highlights) == 0L && !is.null(yml$highlights)) {
      cfg$highlights <- lapply(yml$highlights, function(h) {
        if (is.list(h)) parse_region_string(h$region) else parse_region_string(h)
      })
      cfg$highlight_labels <- lapply(yml$highlights, function(h) {
        if (is.list(h)) h$label %||% "" else ""
      })
    } else if (length(cli$highlights) > 0L) {
      cfg$highlights <- lapply(cli$highlights, parse_region_string)
      while (length(cfg$highlight_labels) < length(cfg$highlights)) {
        cfg$highlight_labels[[length(cfg$highlight_labels) + 1L]] <- ""
      }
    }
  } else {
    cfg$marks <- lapply(cli$marks, parse_mark_spec)
    cfg$highlights <- lapply(cli$highlights, parse_region_string)
    while (length(cfg$highlight_labels) < length(cfg$highlights)) {
      cfg$highlight_labels[[length(cfg$highlight_labels) + 1L]] <- ""
    }
  }

  if (length(cfg$label_list) >= 1L) {
    cfg$label_vec <- unlist(cfg$label_list, use.names = FALSE)
  } else if (!is.null(cfg$labels)) {
    cfg$label_vec <- strsplit(cfg$labels, ",", fixed = TRUE)[[1]]
  } else {
    cfg$label_vec <- c("control", "mutant")
  }
  if (length(cfg$label_vec) != 2L) {
    stop(sprintf(
      "Need exactly two condition labels. Got %d: %s\nUse repeatable --label flags for labels containing commas:\n  --label 'control' --label 'Bap1-KO'",
      length(cfg$label_vec),
      paste(shQuote(cfg$label_vec), collapse = ", ")
    ))
  }

  cfg
}

`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1L && is.na(a))) b else a

nice_ceiling <- function(x) {
  if (x <= 0) return(1)
  x_padded <- x * 1.30
  mag <- 10^floor(log10(x_padded))
  norm <- x_padded / mag
  nice <- if (norm <= 1)      1
          else if (norm <= 2) 2
          else if (norm <= 5) 5
          else                10
  nice * mag
}

validate_config <- function(cfg) {
  if (is.null(cfg$gene) && is.null(cfg$region)) {
    stop("Must specify either --gene or --region (or both via --config).")
  }
  if (length(cfg$marks) == 0L) {
    stop("At least one --mark required (or via YAML config).")
  }

  for (mk in cfg$marks) {
    all_bw <- c(mk$ctrl, mk$mut)
    missing <- all_bw[!file.exists(all_bw)]
    if (length(missing) > 0L) {
      stop(sprintf("Missing BigWig files for mark '%s':\n  %s",
                   mk$name, paste(missing, collapse = "\n  ")))
    }
    if (!grepl("^#[0-9A-Fa-f]{6}$", mk$color)) {
      warning(sprintf("Mark '%s' color '%s' is not a 6-digit hex code", mk$name, mk$color))
    }
  }

  if (!is.null(cfg$hic_loops) && !file.exists(cfg$hic_loops)) {
    stop(sprintf("Hi-C loops file not found: %s", cfg$hic_loops))
  }

  invisible(TRUE)
}

# =============================================================================
# DATA LAYER -- region resolution + BigWig binning
# =============================================================================

resolve_gene_symbol <- function(symbol, cfg) {
  if (!load_optional(cfg$txdb_pkg)) {
    stop(sprintf(
      "TxDb package not installed: %s\nInstall with:\n  BiocManager::install('%s')",
      cfg$txdb_pkg, cfg$txdb_pkg))
  }
  if (!load_optional(cfg$orgdb_pkg)) {
    stop(sprintf(
      "OrgDb package not installed: %s\nInstall with:\n  BiocManager::install('%s')",
      cfg$orgdb_pkg, cfg$orgdb_pkg))
  }
  txdb  <- get(cfg$txdb_pkg, envir = asNamespace(cfg$txdb_pkg))
  orgdb <- get(cfg$orgdb_pkg, envir = asNamespace(cfg$orgdb_pkg))

  entrez <- suppressMessages(
    AnnotationDbi::mapIds(orgdb, keys = symbol,
                          column = "ENTREZID", keytype = "SYMBOL")
  )
  if (is.na(entrez)) stop(sprintf("Gene symbol not found in %s: %s", cfg$orgdb_pkg, symbol))

  all_genes <- suppressMessages(GenomicFeatures::genes(txdb))
  gene_gr <- all_genes[names(all_genes) == entrez]
  if (length(gene_gr) != 1L) {
    stop(sprintf("Could not locate gene %s (entrez %s) in %s",
                 symbol, entrez, cfg$txdb_pkg))
  }

  list(
    entrez = entrez,
    chr    = as.character(seqnames(gene_gr)[1]),
    start  = start(gene_gr),
    end    = end(gene_gr),
    strand = as.character(strand(gene_gr))
  )
}

resolve_region <- function(cfg) {
  if (!is.null(cfg$region)) {
    r <- parse_region_string(cfg$region)
    gr <- GRanges(seqnames = r$chr, ranges = IRanges(r$start, r$end))
    focal_symbol <- cfg$gene_symbol %||% cfg$gene
    mcols(gr)$gene_symbol <- focal_symbol %||% ""

    if (!is.null(focal_symbol) && nzchar(focal_symbol)) {
      info <- tryCatch(resolve_gene_symbol(focal_symbol, cfg),
                       error = function(e) {
                         warning(sprintf(
                           "Could not resolve --gene-symbol '%s': %s\nGene model will show all overlapping genes.",
                           focal_symbol, conditionMessage(e)), call. = FALSE)
                         NULL
                       })
      if (!is.null(info)) {
        mcols(gr)$entrez      <- info$entrez
        mcols(gr)$gene_start  <- info$start
        mcols(gr)$gene_end    <- info$end
        mcols(gr)$gene_strand <- info$strand
        if (info$chr != r$chr) {
          warning(sprintf(
            "Focal gene %s is on %s but --region is on %s. Showing overlapping genes only.",
            focal_symbol, info$chr, r$chr), call. = FALSE)
          mcols(gr)$entrez <- NA_character_
        } else if (info$end < r$start || info$start > r$end) {
          warning(sprintf(
            "Focal gene %s (%s:%d-%d) does not overlap --region %s:%d-%d. Rendering anyway.",
            focal_symbol, info$chr, info$start, info$end,
            r$chr, r$start, r$end), call. = FALSE)
        }
      }
    }
    return(gr)
  }

  info <- resolve_gene_symbol(cfg$gene, cfg)

  ext <- GRanges(seqnames = info$chr,
                 ranges = IRanges(max(1L, info$start - cfg$extend_bp),
                                  info$end + cfg$extend_bp))
  mcols(ext)$gene_symbol  <- cfg$gene
  mcols(ext)$gene_start   <- info$start
  mcols(ext)$gene_end     <- info$end
  mcols(ext)$gene_strand  <- info$strand
  mcols(ext)$entrez       <- info$entrez

  ext
}

import_bigwig_binned <- function(bw_path, region_gr, bin_size = NULL,
                                 sparse = FALSE) {
  chr   <- as.character(seqnames(region_gr)[1])
  rstart <- start(region_gr)
  rend   <- end(region_gr)
  rwidth <- width(region_gr)

  if (is.null(bin_size)) bin_size <- max(50L, as.integer(round(rwidth / 1000)))

  bin_starts <- seq.int(rstart, rend - bin_size, by = bin_size)
  bins_gr <- GRanges(seqnames = chr,
                     ranges = IRanges(start = bin_starts, width = bin_size))

  bw <- rtracklayer::import.bw(bw_path, which = region_gr)
  if (length(bw) == 0L) { bins_gr$score <- 0; return(bins_gr) }

  cov <- coverage(bw, weight = "score")
  if (!chr %in% names(cov)) { bins_gr$score <- 0; return(bins_gr) }

  if (sparse) {
    score_sums <- as.numeric(viewSums(Views(cov[[chr]], ranges(bins_gr))))
    occ <- coverage(bw)
    cpg_counts <- as.numeric(viewSums(Views(occ[[chr]], ranges(bins_gr))))
    bins_gr$score <- ifelse(cpg_counts > 0L, score_sums / cpg_counts, 0)
  } else {
    bins_gr$score <- as.numeric(viewMeans(Views(cov[[chr]], ranges(bins_gr))))
    bins_gr$score[is.na(bins_gr$score)] <- 0
  }
  bins_gr
}

average_bigwigs <- function(bw_paths, region_gr, bin_size = NULL,
                            sparse = FALSE) {
  if (length(bw_paths) == 1L) {
    return(import_bigwig_binned(bw_paths[1], region_gr, bin_size,
                                sparse = sparse))
  }

  bins_list <- lapply(bw_paths, import_bigwig_binned,
                      region_gr = region_gr, bin_size = bin_size,
                      sparse = sparse)
  scores <- do.call(cbind, lapply(bins_list, function(g) g$score))

  result <- bins_list[[1]]
  result$score <- rowMeans(scores, na.rm = TRUE)
  result$score[is.na(result$score)] <- 0
  result
}

bigwig_to_df <- function(gr) {
  data.frame(
    pos   = (start(gr) + end(gr)) / 2,
    score = as.numeric(gr$score)
  )
}

# =============================================================================
# LAYOUT COMPUTATION -- proportional r0/r1 positions for karyoploteR tracks
# =============================================================================

compute_track_layout <- function(mark_data_list, has_hic, has_highlight_labels) {
  signal_w    <- 1.0
  diff_w      <- 0.80
  pair_gap_w  <- 0.12
  mark_gap_w  <- 0.45
  gene_w      <- 0.80
  gene_gap_w  <- 0.30
  scalebar_w  <- 0.30
  hic_w       <- 0.80
  hic_gap_w   <- 0.10
  header_w    <- if (has_highlight_labels) 0.14 else 0

  n_marks <- length(mark_data_list)
  n_diff <- sum(sapply(mark_data_list, function(m) isTRUE(m$diff)))

  total <- scalebar_w + header_w +
           n_marks * (2 * signal_w + pair_gap_w) +
           n_diff * (diff_w + pair_gap_w) +
           max(0, n_marks - 1) * mark_gap_w +
           gene_gap_w + gene_w +
           if (has_hic) (hic_gap_w + hic_w) else 0

  positions <- list()
  cursor <- total

  # Scale bar (topmost)
  r1_sb <- cursor / total
  cursor <- cursor - scalebar_w
  r0_sb <- cursor / total
  scalebar_pos <- c(r0_sb, r1_sb)

  # Highlight header labels
  header_pos <- NULL
  if (has_highlight_labels) {
    r1_hdr <- cursor / total
    cursor <- cursor - header_w
    r0_hdr <- cursor / total
    header_pos <- c(r0_hdr, r1_hdr)
  }

  # Signal tracks (ctrl + mut + optional diff per mark)
  for (i in seq_along(mark_data_list)) {
    md <- mark_data_list[[i]]

    r1_ctrl <- cursor / total
    cursor <- cursor - signal_w
    r0_ctrl <- cursor / total

    cursor <- cursor - pair_gap_w

    r1_mut <- cursor / total
    cursor <- cursor - signal_w
    r0_mut <- cursor / total

    diff_pos <- NULL
    if (isTRUE(md$diff)) {
      cursor <- cursor - pair_gap_w
      r1_diff <- cursor / total
      cursor <- cursor - diff_w
      r0_diff <- cursor / total
      diff_pos <- c(r0_diff, r1_diff)
    }

    positions[[i]] <- list(
      ctrl    = c(r0_ctrl, r1_ctrl),
      mut     = c(r0_mut, r1_mut),
      diff    = diff_pos,
      mark_r0 = if (!is.null(diff_pos)) diff_pos[1] else r0_mut,
      mark_r1 = r1_ctrl
    )

    if (i < n_marks) cursor <- cursor - mark_gap_w
  }

  # Gene model
  cursor <- cursor - gene_gap_w
  r1_gene <- cursor / total
  cursor <- cursor - gene_w
  r0_gene <- max(0, cursor / total)

  # Hi-C arcs
  hic_pos <- NULL
  if (has_hic) {
    cursor <- cursor - hic_gap_w
    r1_hic <- cursor / total
    cursor <- cursor - hic_w
    r0_hic <- max(0, cursor / total)
    hic_pos <- c(r0_hic, r1_hic)
  }

  list(
    tracks   = positions,
    gene     = c(r0_gene, r1_gene),
    scalebar = scalebar_pos,
    header   = header_pos,
    hic      = hic_pos
  )
}

# =============================================================================
# GENE MODEL DATA -- resolve target genes and their exons for rendering
# =============================================================================

resolve_gene_model_data <- function(view_start, view_end, chr, region_gr, cfg) {
  txdb <- get(cfg$txdb_pkg, envir = asNamespace(cfg$txdb_pkg))

  view_gr <- GRanges(seqnames = chr, ranges = IRanges(view_start, view_end))
  all_genes <- suppressMessages(GenomicFeatures::genes(txdb))
  overlap_genes <- subsetByOverlaps(all_genes, view_gr)

  focal_entrez <- if ("entrez" %in% names(mcols(region_gr))) {
    e <- mcols(region_gr)$entrez[1]
    if (is.na(e) || !nzchar(as.character(e))) NULL else as.character(e)
  } else NULL

  focal_symbols <- unlist(cfg$focal_genes, use.names = FALSE)
  if (length(focal_symbols) > 0L) {
    orgdb <- get(cfg$orgdb_pkg, envir = asNamespace(cfg$orgdb_pkg))
    focal_entrez_ids <- suppressMessages(
      AnnotationDbi::mapIds(orgdb, keys = focal_symbols,
                            column = "ENTREZID", keytype = "SYMBOL",
                            multiVals = "first"))
    focal_entrez_ids <- focal_entrez_ids[!is.na(focal_entrez_ids)]
    target_genes <- overlap_genes[names(overlap_genes) %in% focal_entrez_ids]
    if (length(target_genes) == 0L) target_genes <- overlap_genes
  } else if (cfg$show_all_genes) {
    target_genes <- overlap_genes
  } else if (!is.null(focal_entrez)) {
    target_genes <- overlap_genes[names(overlap_genes) == focal_entrez]
    if (length(target_genes) == 0L) target_genes <- overlap_genes
  } else if (!is.null(cfg$gene)) {
    orgdb <- get(cfg$orgdb_pkg, envir = asNamespace(cfg$orgdb_pkg))
    entrez <- suppressMessages(
      AnnotationDbi::mapIds(orgdb, keys = cfg$gene,
                            column = "ENTREZID", keytype = "SYMBOL"))
    target_genes <- overlap_genes[names(overlap_genes) == entrez]
    if (length(target_genes) == 0L) target_genes <- overlap_genes
  } else {
    target_genes <- overlap_genes
  }

  orgdb <- if (load_optional(cfg$orgdb_pkg)) {
    get(cfg$orgdb_pkg, envir = asNamespace(cfg$orgdb_pkg))
  } else NULL

  gene_rows <- lapply(seq_along(target_genes), function(i) {
    g <- target_genes[i]
    entrez_id <- names(g)
    symbol <- if (!is.null(orgdb)) {
      tryCatch(
        AnnotationDbi::mapIds(orgdb, keys = entrez_id,
                              column = "SYMBOL", keytype = "ENTREZID"),
        error = function(e) entrez_id
      )
    } else entrez_id
    if (is.na(symbol)) symbol <- entrez_id

    exons_by_gene <- suppressMessages(
      GenomicFeatures::exonsBy(txdb, by = "gene")
    )
    if (entrez_id %in% names(exons_by_gene)) {
      ex <- GenomicRanges::reduce(exons_by_gene[[entrez_id]])
    } else {
      ex <- GRanges()
    }

    list(
      symbol     = as.character(symbol),
      gene_start = start(g),
      gene_end   = end(g),
      strand     = as.character(strand(g)),
      exons      = ex
    )
  })

  # Assign y-rows (greedy packing to avoid gene overlaps)
  if (length(gene_rows) <= 1L || !cfg$show_all_genes) {
    for (i in seq_along(gene_rows)) gene_rows[[i]]$y <- 0
    n_rows <- 1L
  } else {
    sorted_idx <- order(sapply(gene_rows, function(g) g$gene_start))
    row_ends <- integer(0)
    y_assignments <- integer(length(gene_rows))
    pad <- (view_end - view_start) * 0.01
    for (idx in sorted_idx) {
      g <- gene_rows[[idx]]
      slot <- which(row_ends < g$gene_start - pad)
      if (length(slot) == 0L) {
        row_ends <- c(row_ends, g$gene_end)
        y_assignments[idx] <- length(row_ends) - 1L
      } else {
        y_assignments[idx] <- slot[1] - 1L
        row_ends[slot[1]] <- g$gene_end
      }
    }
    for (i in seq_along(gene_rows)) gene_rows[[i]]$y <- y_assignments[i]
    n_rows <- max(y_assignments) + 1L
  }

  list(genes = gene_rows, n_rows = n_rows)
}

# =============================================================================
# RENDERING -- karyoploteR-based figure construction
# =============================================================================

render_locus <- function(region_gr, mark_data_list, layout, gene_data,
                         loops_df, cfg) {
  chr        <- as.character(seqnames(region_gr)[1])
  view_start <- start(region_gr)
  view_end   <- end(region_gr)

  ctrl_label <- cfg$label_vec[1]
  mut_label  <- cfg$label_vec[2]

  # Margins: left=label space, top=room for scale bar, right+bottom=small
  old_par <- par(mar = c(1.5, 5.5, 2.5, 1.5))
  on.exit(par(old_par))

  kp <- plotKaryotype(
    genome      = cfg$genome,
    chromosomes = chr,
    zoom        = region_gr,
    plot.type   = 4,
    labels.plotter   = NULL,
    ideogram.plotter = NULL,
    main        = ""
  )

  # -- Highlight rectangles (full-height, behind everything) --
  if (length(cfg$highlights) > 0L) {
    for (h in cfg$highlights) {
      if (h$chr != chr) next
      kpRect(kp,
             chr = chr, x0 = h$start, x1 = h$end,
             y0 = 0, y1 = 1, r0 = 0, r1 = 1,
             col = adjustcolor("grey85", alpha.f = 0.55),
             border = NA)
    }
  }

  label_x <- view_start + (view_end - view_start) * 0.015

  # -- Highlight header labels --
  if (!is.null(layout$header)) {
    for (i in seq_along(cfg$highlights)) {
      if (i > length(cfg$highlight_labels)) next
      lbl <- cfg$highlight_labels[[i]]
      if (!nzchar(lbl)) next
      h <- cfg$highlights[[i]]
      mid_x <- (h$start + h$end) / 2
      kpText(kp,
             chr = chr, x = mid_x, y = 0.5,
             r0 = layout$header[1], r1 = layout$header[2],
             labels = lbl, cex = 0.7, col = "black")
    }
  }

  # -- Scale bar (top-right) --
  region_kb <- (view_end - view_start) / 1000
  scale_kb <- cfg$scale_bar_kb
  if (is.null(scale_kb)) {
    candidates <- c(1, 2, 5, 10, 20, 50, 100, 200, 500, 1000)
    target <- region_kb * 0.08
    scale_kb <- candidates[which.min(abs(candidates - target))]
  }
  scale_bp  <- scale_kb * 1000
  bar_end   <- view_end - (view_end - view_start) * 0.02
  bar_start <- bar_end - scale_bp
  bar_label <- if (scale_kb >= 1) sprintf("%g kb", scale_kb) else sprintf("%g bp", scale_bp)

  sb <- layout$scalebar
  kpSegments(kp,
             chr = chr, x0 = bar_start, x1 = bar_end,
             y0 = 0.70, y1 = 0.70,
             r0 = sb[1], r1 = sb[2],
             lwd = 2.0, col = "black")
  kpText(kp,
         chr = chr, x = (bar_start + bar_end) / 2, y = 0.30,
         r0 = sb[1], r1 = sb[2],
         labels = bar_label, cex = 0.75, col = "black")

  # -- Signal tracks (per mark: ctrl + mut + optional diff) --
  for (i in seq_along(mark_data_list)) {
    md  <- mark_data_list[[i]]
    pos <- layout$tracks[[i]]

    ctrl_df <- bigwig_to_df(md$ctrl)
    mut_df  <- bigwig_to_df(md$mut)
    ylim    <- md$ylim
    color   <- md$color

    chr_vec <- rep(chr, nrow(ctrl_df))

    # Cap scores at ylim
    ctrl_df$score <- pmin(ctrl_df$score, ylim)
    mut_df$score  <- pmin(mut_df$score, ylim)

    # Ctrl track
    kpArea(kp,
           chr = chr_vec, x = ctrl_df$pos, y = ctrl_df$score,
           ymin = 0, ymax = ylim, base.y = 0,
           r0 = pos$ctrl[1], r1 = pos$ctrl[2],
           col = color,
           border = color, lwd = 0.3)
    # Baseline
    kpAbline(kp, h = 0,
             r0 = pos$ctrl[1], r1 = pos$ctrl[2],
             ymin = 0, ymax = ylim,
             col = color, lwd = 0.6)

    # Scale annotation on ctrl track
    scale_text <- sprintf("0-%g", ylim)
    kpText(kp,
           chr = chr, x = label_x, y = ylim * 0.94,
           r0 = pos$ctrl[1], r1 = pos$ctrl[2],
           ymin = 0, ymax = ylim,
           labels = scale_text, cex = 0.70, col = color,
           pos = 4, offset = 0)

    # Condition label on ctrl track
    kpText(kp,
           chr = chr, x = label_x, y = ylim * 0.50,
           r0 = pos$ctrl[1], r1 = pos$ctrl[2],
           ymin = 0, ymax = ylim,
           labels = ctrl_label, cex = 0.60, col = "black",
           pos = 4, offset = 0)

    # Mut track
    kpArea(kp,
           chr = chr_vec, x = mut_df$pos, y = mut_df$score,
           ymin = 0, ymax = ylim, base.y = 0,
           r0 = pos$mut[1], r1 = pos$mut[2],
           col = color,
           border = color, lwd = 0.3)
    kpAbline(kp, h = 0,
             r0 = pos$mut[1], r1 = pos$mut[2],
             ymin = 0, ymax = ylim,
             col = color, lwd = 0.6)

    # Condition label on mut track
    kpText(kp,
           chr = chr, x = label_x, y = ylim * 0.70,
           r0 = pos$mut[1], r1 = pos$mut[2],
           ymin = 0, ymax = ylim,
           labels = mut_label, cex = 0.60, col = "black",
           pos = 4, offset = 0,
           font = if (cfg$genotype_italic) 3 else 1)

    # Optional diff track
    if (isTRUE(md$diff) && !is.null(pos$diff)) {
      diff_vec <- mut_df$score - ctrl_df$score
      if (isTRUE(md$sparse)) diff_vec <- diff_vec * 100

      abs_max <- max(abs(diff_vec), na.rm = TRUE)
      if (abs_max <= 0 || !is.finite(abs_max)) abs_max <- 1
      diff_ylim <- nice_ceiling(abs_max)

      loss_color <- if (load_optional("colorspace")) {
        colorspace::desaturate(colorspace::lighten(color, 0.3), 0.4)
      } else {
        adjustcolor(color, alpha.f = 0.4)
      }

      # Gains (mut > ctrl): positive, above zero line
      gain_y <- pmax(diff_vec, 0)
      kpArea(kp,
             chr = chr_vec, x = ctrl_df$pos, y = gain_y,
             ymin = -diff_ylim, ymax = diff_ylim, base.y = 0,
             r0 = pos$diff[1], r1 = pos$diff[2],
             col = adjustcolor(color, alpha.f = 0.90),
             border = color, lwd = 0.3)

      # Losses (mut < ctrl): negative, below zero line
      loss_y <- pmin(diff_vec, 0)
      kpArea(kp,
             chr = chr_vec, x = ctrl_df$pos, y = loss_y,
             ymin = -diff_ylim, ymax = diff_ylim, base.y = 0,
             r0 = pos$diff[1], r1 = pos$diff[2],
             col = adjustcolor(loss_color, alpha.f = 0.90),
             border = loss_color, lwd = 0.3)

      # Zero line
      kpAbline(kp, h = 0,
               r0 = pos$diff[1], r1 = pos$diff[2],
               ymin = -diff_ylim, ymax = diff_ylim,
               col = "black", lwd = 0.5)

      # Diff scale label
      unit_label <- if (isTRUE(md$sparse)) "%" else "Δ"
      diff_scale_text <- sprintf("±%g%s", diff_ylim, unit_label)
      kpText(kp,
             chr = chr, x = label_x, y = diff_ylim * 0.90,
             r0 = pos$diff[1], r1 = pos$diff[2],
             ymin = -diff_ylim, ymax = diff_ylim,
             labels = diff_scale_text, cex = 0.55, col = color,
             pos = 4, offset = 0)

      kpText(kp,
             chr = chr, x = label_x, y = diff_ylim * 0.20,
             r0 = pos$diff[1], r1 = pos$diff[2],
             ymin = -diff_ylim, ymax = diff_ylim,
             labels = "difference", cex = 0.50, col = "grey40",
             pos = 4, offset = 0)
    }

    # Rotated mark label (left side)
    kpAddLabels(kp,
                labels       = md$name,
                r0           = pos$mark_r0,
                r1           = pos$mark_r1,
                col          = color,
                cex          = 0.85,
                srt          = 90,
                pos          = 2,
                label.margin = 0.04,
                font         = 2)
  }

  # -- Gene model --
  gene_r0 <- layout$gene[1]
  gene_r1 <- layout$gene[2]
  n_gene_rows <- gene_data$n_rows

  for (g in gene_data$genes) {
    # Compute per-row y positioning within the gene panel
    row_height <- 1 / max(1, n_gene_rows)
    y_center   <- 1 - (g$y + 0.5) * row_height

    # Intron backbone
    kpSegments(kp,
               chr = chr,
               x0 = max(g$gene_start, view_start),
               x1 = min(g$gene_end, view_end),
               y0 = y_center, y1 = y_center,
               r0 = gene_r0, r1 = gene_r1,
               lwd = 0.3, col = "black")

    # Filled exon rectangles
    if (length(g$exons) > 0L) {
      exon_height <- row_height * 0.35
      kpRect(kp,
             chr  = as.character(seqnames(g$exons)),
             x0   = start(g$exons),
             x1   = end(g$exons),
             y0   = y_center - exon_height,
             y1   = y_center + exon_height,
             r0   = gene_r0, r1 = gene_r1,
             col  = "black", border = NA)

      # Strand arrowheads on introns
      if (length(g$exons) >= 2L) {
        intron_starts <- end(g$exons)[-length(g$exons)]
        intron_ends   <- start(g$exons)[-1L]
        midpoints <- (intron_starts + intron_ends) / 2
        intron_widths <- intron_ends - intron_starts
        min_arrow_span <- (view_end - view_start) * 0.005
        keep <- intron_widths > min_arrow_span
        if (any(keep)) {
          arrow_len <- (view_end - view_start) * 0.003
          mids <- midpoints[keep]
          if (g$strand == "+") {
            kpSegments(kp,
                       chr = rep(chr, length(mids)),
                       x0 = mids - arrow_len, x1 = mids + arrow_len,
                       y0 = y_center, y1 = y_center,
                       r0 = gene_r0, r1 = gene_r1,
                       lwd = 0.3, col = "black")
            kpArrows(kp,
                     chr = rep(chr, length(mids)),
                     x0 = mids - arrow_len, x1 = mids + arrow_len,
                     y0 = y_center, y1 = y_center,
                     r0 = gene_r0, r1 = gene_r1,
                     length = 0.03, col = "black", lwd = 0.3)
          } else {
            kpSegments(kp,
                       chr = rep(chr, length(mids)),
                       x0 = mids + arrow_len, x1 = mids - arrow_len,
                       y0 = y_center, y1 = y_center,
                       r0 = gene_r0, r1 = gene_r1,
                       lwd = 0.3, col = "black")
            kpArrows(kp,
                     chr = rep(chr, length(mids)),
                     x0 = mids + arrow_len, x1 = mids - arrow_len,
                     y0 = y_center, y1 = y_center,
                     r0 = gene_r0, r1 = gene_r1,
                     length = 0.03, col = "black", lwd = 0.3)
          }
        }
      }
    }

    # Italic gene symbol above gene body (keep within gene panel bounds)
    sym_x <- (max(g$gene_start, view_start) + min(g$gene_end, view_end)) / 2
    sym_y <- y_center + row_height * 0.28
    kpText(kp,
           chr = chr, x = sym_x, y = min(sym_y, 0.85),
           r0 = gene_r0, r1 = gene_r1,
           labels = g$symbol, cex = 0.9, col = "black",
           font = 3)
  }

  # -- Hi-C loop arcs --
  if (!is.null(loops_df) && !is.null(layout$hic)) {
    in_view <- loops_df %>%
      dplyr::filter(
        (chr1 == chr & ((start1 >= view_start & start1 <= view_end) |
                        (end1 >= view_start & end1 <= view_end))) |
        (chr2 == chr & ((start2 >= view_start & start2 <= view_end) |
                        (end2 >= view_start & end2 <= view_end)))
      )

    if (nrow(in_view) > 0L) {
      lost_color   <- "#CB181D"
      gained_color <- "#2171B5"

      anchor1 <- GRanges(seqnames = in_view$chr1,
                         ranges = IRanges(start = in_view$start1,
                                          end = in_view$end1))
      anchor2 <- GRanges(seqnames = in_view$chr2,
                         ranges = IRanges(start = in_view$start2,
                                          end = in_view$end2))
      loop_colors <- ifelse(in_view$direction == "up_in_mutant",
                            gained_color, lost_color)

      # Scale arc heights by span relative to max span
      span_max <- max(in_view$start2 - in_view$start1)
      arc_heights <- (abs(in_view$start2 - in_view$start1) / span_max) * 0.9

      kpPlotLinks(kp,
                  data = anchor1, data2 = anchor2,
                  r0 = layout$hic[1], r1 = layout$hic[2],
                  col = adjustcolor(loop_colors, alpha.f = 0.60),
                  border = loop_colors,
                  arch.height = arc_heights,
                  lwd = 0.8)
    }
  }

  invisible(kp)
}

# =============================================================================
# OUTPUT -- save in PDF, SVG, PNG, JPEG formats
# =============================================================================

save_multiformat_kp <- function(render_fn, base_path, width, height,
                                dpi = 300) {
  out_dir <- dirname(base_path)
  if (!dir.exists(out_dir) && nzchar(out_dir)) {
    tryCatch(
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE),
      error = function(e) stop(sprintf(
        "Could not create output directory '%s': %s",
        out_dir, conditionMessage(e)))
    )
  }
  stem <- basename(base_path)
  parent <- dirname(base_path)
  fig_dir <- file.path(parent, stem)
  if (!dir.exists(fig_dir)) {
    tryCatch(
      dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE),
      error = function(e) stop(sprintf(
        "Could not create figure subfolder '%s': %s",
        fig_dir, conditionMessage(e)))
    )
  }
  prefix <- file.path(fig_dir, stem)

  # PDF
  pdf(paste0(prefix, ".pdf"), width = width, height = height)
  tryCatch(render_fn(), finally = dev.off())

  # SVG (via svglite for Illustrator-friendly output)
  svglite::svglite(paste0(prefix, ".svg"), width = width, height = height)
  tryCatch(render_fn(), finally = dev.off())

  # PNG
  png(paste0(prefix, ".png"),
      width = width * dpi, height = height * dpi, res = dpi)
  tryCatch(render_fn(), finally = dev.off())

  # JPEG
  jpeg(paste0(prefix, ".jpg"),
       width = width * dpi, height = height * dpi,
       res = dpi, quality = 95)
  tryCatch(render_fn(), finally = dev.off())

  cat(sprintf("  Saved: %s/{pdf,svg,png,jpg}\n", fig_dir))
}

# =============================================================================
# MAIN
# =============================================================================

main <- function() {
  raw_args <- commandArgs(trailingOnly = TRUE)
  if (length(raw_args) == 0L) {
    print_help()
    quit(save = "no", status = 0L)
  }
  cli <- parse_cli_args(raw_args)

  if (cli$help) {
    print_help()
    quit(save = "no", status = 0L)
  }

  cat("================================================================\n")
  cat("pub_browser.R -- publication-quality genome browser (karyoploteR)\n")
  cat("================================================================\n\n")

  cfg <- resolve_config(cli)
  validate_config(cfg)

  # ---- 1. Resolve view region ----
  cat("Resolving region...\n")
  region_gr <- resolve_region(cfg)
  chr        <- as.character(seqnames(region_gr)[1])
  view_start <- start(region_gr)
  view_end   <- end(region_gr)

  cat(sprintf("  %s:%s-%s (%.1f kb)\n",
              chr,
              format(view_start, big.mark = ","),
              format(view_end, big.mark = ","),
              (view_end - view_start) / 1000))

  focal_symbols <- unlist(cfg$focal_genes, use.names = FALSE)
  if (length(focal_symbols) > 0L) {
    cat(sprintf("  gene model: focal genes %s\n",
                paste(focal_symbols, collapse = ", ")))
  } else if (cfg$show_all_genes) {
    cat("  gene model: all overlapping genes\n")
  } else if ("entrez" %in% names(mcols(region_gr))) {
    fe <- mcols(region_gr)$entrez[1]
    if (!is.na(fe) && nzchar(as.character(fe))) {
      focal_sym <- mcols(region_gr)$gene_symbol[1]
      cat(sprintf("  gene model: focal gene %s (entrez %s)\n",
                  focal_sym, fe))
    } else {
      cat("  gene model: all overlapping genes (no focal gene resolved)\n")
    }
  } else {
    cat("  gene model: all overlapping genes\n")
  }

  for (h in cfg$highlights) {
    if (h$chr != chr) {
      warning(sprintf("Highlight chr '%s' != view chr '%s' -- ignoring",
                      h$chr, chr))
    }
  }

  # ---- 2. Load BigWig data ----
  cat(sprintf("Loading %d mark(s)...\n", length(cfg$marks)))
  mark_data <- list()
  for (mk in cfg$marks) {
    cat(sprintf("  %s (%d ctrl, %d mut)...\n",
                mk$name, length(mk$ctrl), length(mk$mut)))
    ctrl_gr <- average_bigwigs(mk$ctrl, region_gr, cfg$bin_size,
                               sparse = mk$sparse)
    mut_gr  <- average_bigwigs(mk$mut, region_gr, cfg$bin_size,
                               sparse = mk$sparse)
    raw_max <- max(c(ctrl_gr$score, mut_gr$score), na.rm = TRUE)
    if (raw_max <= 0 || !is.finite(raw_max)) raw_max <- 1
    auto_ylim <- nice_ceiling(raw_max)
    mark_data[[length(mark_data) + 1L]] <- list(
      name   = mk$name,
      ctrl   = ctrl_gr,
      mut    = mut_gr,
      color  = mk$color,
      ylim   = if (!is.null(mk$ylim)) mk$ylim else auto_ylim,
      diff   = mk$diff,
      sparse = mk$sparse
    )
  }

  # ---- 3. Resolve gene model data ----
  cat("Resolving gene model...\n")
  gene_data <- resolve_gene_model_data(view_start, view_end, chr,
                                       region_gr, cfg)

  # ---- 4. Load optional Hi-C loops ----
  loops_df <- NULL
  if (!is.null(cfg$hic_loops)) {
    cat("Loading Hi-C loops...\n")
    loops_df <- read.table(cfg$hic_loops, header = TRUE, sep = "\t",
                           stringsAsFactors = FALSE)
  }

  # ---- 5. Compute layout ----
  has_hic <- !is.null(loops_df)
  has_highlight_labels <- any(nzchar(unlist(cfg$highlight_labels)))
  layout <- compute_track_layout(mark_data, has_hic, has_highlight_labels)

  # ---- 6. Compute figure dimensions ----
  n_signal_pairs <- length(mark_data) * 2L
  n_diff_tracks  <- sum(sapply(mark_data, function(m) isTRUE(m$diff)))
  n_marks <- length(mark_data)

  total_height <- n_signal_pairs * cfg$track_height +
                  n_diff_tracks * cfg$track_height * 0.60 +
                  0.40 +                                  # gene model
                  0.35 +                                  # scale bar
                  max(0, n_marks - 1) * 0.25 +            # inter-mark gaps
                  n_marks * 0.08 +                        # intra-mark gaps
                  0.20 +                                  # gene-track gap
                  (if (has_hic) 1.0 else 0) +             # Hi-C panel
                  (if (has_highlight_labels) 0.18 else 0) +
                  0.4                                     # margins

  # ---- 7. Render and save ----
  cat("Rendering figure...\n")
  do_render <- function() {
    render_locus(region_gr, mark_data, layout, gene_data, loops_df, cfg)
  }

  cat("Saving outputs...\n")
  save_multiformat_kp(do_render, cfg$output, cfg$fig_width, total_height)

  cat(sprintf("\nDone. Figure dimensions: %.1f x %.1f inches\n",
              cfg$fig_width, total_height))
}

if (sys.nframe() == 0L) {
  main()
}
