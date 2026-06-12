#!/usr/bin/env Rscript
# scripts/pub_browser.R
# Publication-quality genome browser locus views (lab-style aesthetic)
#
# Renders multi-track locus figures directly to SVG/PDF/JPEG in the style of
# manually-polished Illustrator figures: filled area curves with same color
# per mark, rotated colored mark labels, slim black gene models with arrowed
# strand direction, italic gene symbols, scale-bar coordinate marker, and
# optional gray highlight rectangles. No post-processing required.
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
  library(ggplot2)
  library(patchwork)
  library(cowplot)
  library(grid)
  library(rtracklayer)
  library(GenomicRanges)
  library(GenomicFeatures)
  library(dplyr)
  library(svglite)
})

# Optional packages loaded on demand
load_optional <- function(pkg) {
  suppressPackageStartupMessages(suppressWarnings(
    requireNamespace(pkg, quietly = TRUE)
  ))
}

# =============================================================================
# CLI PARSING
# =============================================================================

# Safe numeric parsing: stop with a useful message instead of returning NA
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

# Manual parser since we need repeatable --mark / --highlight / --label flags
parse_cli_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  cfg <- list(
    gene             = NULL,
    region           = NULL,
    genome           = "mm10",
    marks            = list(),
    config_file      = NULL,
    labels           = NULL,        # legacy 'ctrl,mut' shortcut
    label_list       = list(),      # preferred: repeatable --label flag
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
pub_browser.R -- Publication-quality genome browser locus tool

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
                             Supports markdown for italic/superscript via
                             ggtext, e.g. '*Bap1*^f/f^,Math1-cre'.
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
                                   (e.g., 'Bap1^f/f,Math1-cre')

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

  # With highlight and richtext genotype label (italic + superscript)
  Rscript pub_browser.R \\
    --gene Anxa2r1 \\
    --mark 'H2AK119ub:#2171B5:ctrl.bw:mut.bw' \\
    --label 'control' \\
    --label '*Bap1*^f/f^,Math1-cre' \\
    --highlight 'chr3:9505000-9510000' \\
    --highlight-label 'differential region' \\
    --output out/anxa2r1_pub

  # Region-based view with a focal gene resolved from symbol
  # (PI workflow: render with --gene first, then extend the view with --region)
  Rscript pub_browser.R \\
    --region 'chr3:9450000-9520000' \\
    --gene-symbol Anxa2r1 \\
    --mark 'H3K27me3:#CB181D:ctrl.bw:mut.bw' \\
    --output out/anxa2r1_extended

  # Region-based with all overlapping genes
  Rscript pub_browser.R \\
    --region 'chr11:108400000-109100000' \\
    --gene-symbol Syt1 \\
    --show-all-genes \\
    --mark 'H3K27me1:#2171B5:rep1.bw:rep2.bw' \\
    --labels 'replicate 1,replicate 2' \\
    --output out/syt1_region

  # YAML config (recommended for complex setups)
  Rscript pub_browser.R --config configs/syt1_locus.yaml

YAML CONFIG FORMAT:

  genome: mm10
  gene: Syt1
  extend_bp: 50000
  labels: ['control', 'Bap1-KO']
  output: out/syt1_browser
  scale_bar_kb: 50

  marks:
    - name: H2AK119ub
      color: '#2171B5'
      ctrl: /path/to/H2AK119ubCtrl.bw
      mut:  /path/to/H2AK119ubMut.bw
    - name: '5mC'
      color: '#2166AC'
      ctrl: [rep1_ctrl.bw, rep2_ctrl.bw, rep3_ctrl.bw, rep4_ctrl.bw]
      mut:  [rep1_mut.bw,  rep2_mut.bw,  rep3_mut.bw,  rep4_mut.bw]
      average: true

  highlights:
    - region: 'chr11:108700000-108800000'
      label:  'differential region'

  hic_loops: /path/to/loops.tsv
")
}

# Parse 'name:color:ctrl:mut' mark spec into structured list
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

# Parse 'chr:start-end' string into list(chr, start, end)
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
  # Start with CLI values
  cfg <- cli

  # Merge YAML if provided (CLI overrides YAML for scalars; lists are replaced)
  if (!is.null(cli$config_file)) {
    if (!load_optional("yaml")) {
      stop("Package 'yaml' required for --config. Install with install.packages('yaml').")
    }
    yml <- yaml::read_yaml(cli$config_file)

    # Scalars: CLI > YAML > default
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

    # Focal genes: CLI --focal-gene list, then YAML focal_genes array
    if (length(cli$focal_genes) == 0L && !is.null(yml$focal_genes)) {
      cfg$focal_genes <- as.list(as.character(yml$focal_genes))
    }

    # Labels: prefer CLI --label list, then YAML labels array, then CLI --labels (legacy)
    if (length(cli$label_list) == 0L && !is.null(yml$labels)) {
      cfg$label_list <- as.list(as.character(yml$labels))
    }

    # Marks: YAML overrides CLI marks if YAML has them
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

    # Highlights from YAML
    if (length(cli$highlights) == 0L && !is.null(yml$highlights)) {
      cfg$highlights <- lapply(yml$highlights, function(h) {
        if (is.list(h)) parse_region_string(h$region) else parse_region_string(h)
      })
      cfg$highlight_labels <- lapply(yml$highlights, function(h) {
        if (is.list(h)) h$label %||% "" else ""
      })
    } else if (length(cli$highlights) > 0L) {
      cfg$highlights <- lapply(cli$highlights, parse_region_string)
      # Pad highlight_labels to match length
      while (length(cfg$highlight_labels) < length(cfg$highlights)) {
        cfg$highlight_labels[[length(cfg$highlight_labels) + 1L]] <- ""
      }
    }
  } else {
    # No YAML — parse mark specs from CLI
    cfg$marks <- lapply(cli$marks, parse_mark_spec)
    cfg$highlights <- lapply(cli$highlights, parse_region_string)
    while (length(cfg$highlight_labels) < length(cfg$highlights)) {
      cfg$highlight_labels[[length(cfg$highlight_labels) + 1L]] <- ""
    }
  }

  # Resolve final label vector: prefer --label list, then --labels comma-split, else default
  if (length(cfg$label_list) >= 1L) {
    cfg$label_vec <- unlist(cfg$label_list, use.names = FALSE)
  } else if (!is.null(cfg$labels)) {
    cfg$label_vec <- strsplit(cfg$labels, ",", fixed = TRUE)[[1]]
  } else {
    cfg$label_vec <- c("control", "mutant")
  }
  if (length(cfg$label_vec) != 2L) {
    stop(sprintf(
      "Need exactly two condition labels. Got %d: %s\nUse repeatable --label flags for labels containing commas:\n  --label 'control' --label '*Bap1*^f/f^,Math1-cre'",
      length(cfg$label_vec),
      paste(shQuote(cfg$label_vec), collapse = ", ")
    ))
  }

  cfg
}

`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1L && is.na(a))) b else a

# Round y-axis max up to a "nice" number (e.g., 18.4 -> 50, 113.3 -> 200, 52.6 -> 100).
# Choices restricted to {1, 2, 5, 10} * 10^n to match PI's typical track scales
# (which use 0-10, 0-50, 0-100, 0-200, etc.). ~30% headroom for label space.
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

# Resolve a gene symbol to its (entrez, start, end, strand) tuple. Stops with
# an actionable message if the TxDb/OrgDb packages are missing.
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

# Resolve gene symbol to extended GRanges, or build a GRanges from --region.
# When --region is paired with --gene-symbol, also resolve the symbol so the
# gene model panel knows the focal gene (correct strand/exons).
resolve_region <- function(cfg) {
  if (!is.null(cfg$region)) {
    r <- parse_region_string(cfg$region)
    gr <- GRanges(seqnames = r$chr, ranges = IRanges(r$start, r$end))
    focal_symbol <- cfg$gene_symbol %||% cfg$gene
    mcols(gr)$gene_symbol <- focal_symbol %||% ""

    # Resolve focal gene info when a symbol is provided alongside --region
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

# Bin a single BigWig over the view region using viewMeans (correct for
# continuous coverage tracks like RNA-seq, ChIP-seq, ATAC-seq, etc.)
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

# Average multiple BigWig files (replicates) into a single binned GRanges
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

# Convert binned GRanges to data.frame for ggplot (use bin midpoints)
bigwig_to_df <- function(gr) {
  data.frame(
    pos   = (start(gr) + end(gr)) / 2,
    score = as.numeric(gr$score)
  )
}

# =============================================================================
# TRACK BUILDERS -- each returns one ggplot object
# =============================================================================

# Detect if a label string contains richtext markers (*italic*, ^super^, _sub_, HTML tags).
# Returns TRUE only if ggtext is also available; warns once when fallback used.
.RICHTEXT_WARNED <- new.env(parent = emptyenv())
.RICHTEXT_WARNED$shown <- FALSE
is_richtext_label <- function(s) {
  if (is.null(s) || !nzchar(s)) return(FALSE)
  has_markup <- grepl("\\*.+?\\*|\\^.+?\\^|<[a-z/]", s)
  if (!has_markup) return(FALSE)
  if (!load_optional("ggtext")) {
    if (!.RICHTEXT_WARNED$shown) {
      warning(
        "Label contains markdown markers (*, ^, <tags>) but 'ggtext' is not installed. ",
        "Falling back to plain text. Install with install.packages('ggtext') for italic/superscript support.",
        call. = FALSE
      )
      .RICHTEXT_WARNED$shown <- TRUE
    }
    return(FALSE)
  }
  TRUE
}

# Convert ergonomic markdown markers to the HTML that ggtext/gridtext understand.
# - *text* stays as-is (markdown italic, handled natively)
# - ^text^ -> <sup>text</sup>  (caret superscript, not native markdown)
# - ~text~ -> <sub>text</sub>  (tilde subscript, not native markdown)
preprocess_richtext <- function(s) {
  if (is.null(s) || !nzchar(s)) return(s)
  # Superscript: ^...^  (lazy match so '^f/f^,Math1-cre' captures just 'f/f')
  s <- gsub("\\^([^\\^]+?)\\^", "<sup>\\1</sup>", s, perl = TRUE)
  # Subscript: ~...~ (avoid matching markdown strikethrough ~~text~~)
  s <- gsub("(?<!~)~([^~]+?)~(?!~)", "<sub>\\1</sub>", s, perl = TRUE)
  s
}

# Add a condition/scale label as plain text or richtext depending on label content
add_label_layer <- function(p, label, x, y, color = "black",
                            size = 2.2, hjust = 0, vjust = 1,
                            italic = FALSE) {
  if (is.null(label) || !nzchar(label)) return(p)
  if (is_richtext_label(label)) {
    df <- data.frame(x = x, y = y, label = preprocess_richtext(label))
    p + ggtext::geom_richtext(
      data = df, inherit.aes = FALSE,
      aes(x = x, y = y, label = label),
      hjust = hjust, vjust = vjust,
      size = size, color = color,
      fill = NA, label.color = NA,
      label.padding = unit(0, "pt"),
      label.margin = unit(0, "pt")
    )
  } else {
    p + annotate("text",
                 x = x, y = y, label = label,
                 hjust = hjust, vjust = vjust,
                 size = size, color = color,
                 fontface = if (italic) "italic" else "plain")
  }
}

# Single condition signal panel: filled area curve over the view region
build_signal_panel <- function(signal_df, color, ylim,
                               view_start, view_end,
                               condition_label = NULL,
                               scale_label = NULL,
                               highlights = list(),
                               italic_label = FALSE) {
  # Cap signal at ylim for visual clarity (rare overshoots stay at ceiling)
  signal_df$score <- pmin(signal_df$score, ylim)

  p <- ggplot(signal_df, aes(x = pos, y = score))

  # Highlight rectangles first (render behind signal)
  if (length(highlights) > 0L) {
    for (h in highlights) {
      p <- p + annotate("rect",
                        xmin = h$start, xmax = h$end,
                        ymin = -Inf, ymax = Inf,
                        fill = "grey85", alpha = 0.55, color = NA)
    }
  }

  # Filled area curve with thin matching outline -- the PI aesthetic
  p <- p +
    geom_area(fill = color, color = color,
              alpha = 0.90, linewidth = 0.15) +
    # Thin baseline at y=0 in mark color (subtle)
    geom_hline(yintercept = 0, color = color, linewidth = 0.3) +
    scale_y_continuous(limits = c(0, ylim),
                       expand = expansion(mult = c(0, 0.02))) +
    scale_x_continuous(limits = c(view_start, view_end),
                       expand = c(0, 0)) +
    theme_void() +
    theme(
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin      = margin(0, 2, 0, 0, "pt")
    )

  # Scale annotation in mark color (ctrl track only).
  # Positioned at the very top of the panel; nice_ceiling rounding provides
  # ~25%+ headroom above data so this overlays empty space, not peaks.
  if (!is.null(scale_label) && nzchar(scale_label)) {
    p <- add_label_layer(
      p, scale_label,
      x = view_start + (view_end - view_start) * 0.05,
      y = ylim * 0.94,
      color = color, size = 2.2,
      hjust = 0, vjust = 1
    )
  }

  # Condition label: small black text just below the scale annotation
  if (!is.null(condition_label) && nzchar(condition_label)) {
    text_x <- view_start + (view_end - view_start) * 0.015
    # Place lower on mut track (no scale label above) to mirror PI style
    y_pos <- if (!is.null(scale_label)) ylim * 0.72 else ylim * 0.85
    p <- add_label_layer(
      p, condition_label,
      x = text_x, y = y_pos,
      color = "black", size = 2.2,
      hjust = 0, vjust = 1,
      italic = italic_label
    )
  }

  p
}

# Percent difference panel: (mut - ctrl) as diverging area chart.
# Positive values (gain in mut) filled in mark color above zero line;
# negative values (loss in mut) filled in a desaturated shade below.
build_diff_panel <- function(ctrl_df, mut_df, color,
                             view_start, view_end,
                             highlights = list(),
                             is_fraction = FALSE) {
  diff_df <- data.frame(
    pos  = ctrl_df$pos,
    diff = mut_df$score - ctrl_df$score
  )
  if (is_fraction) diff_df$diff <- diff_df$diff * 100

  diff_df$gain <- pmax(diff_df$diff, 0)
  diff_df$loss <- pmin(diff_df$diff, 0)

  abs_max <- max(abs(diff_df$diff), na.rm = TRUE)
  if (abs_max <= 0 || !is.finite(abs_max)) abs_max <- 1
  ylim <- nice_ceiling(abs_max)

  loss_color <- colorspace::desaturate(colorspace::lighten(color, 0.3), 0.4)

  p <- ggplot(diff_df, aes(x = pos))

  if (length(highlights) > 0L) {
    for (h in highlights) {
      p <- p + annotate("rect",
                        xmin = h$start, xmax = h$end,
                        ymin = -Inf, ymax = Inf,
                        fill = "grey85", alpha = 0.55, color = NA)
    }
  }

  p <- p +
    geom_area(aes(y = gain), fill = color, alpha = 0.90, linewidth = 0.15) +
    geom_area(aes(y = loss), fill = loss_color, alpha = 0.90, linewidth = 0.15) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.25) +
    scale_y_continuous(limits = c(-ylim, ylim),
                       expand = expansion(mult = c(0.02, 0.02))) +
    scale_x_continuous(limits = c(view_start, view_end),
                       expand = c(0, 0)) +
    theme_void() +
    theme(
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin      = margin(0, 2, 0, 0, "pt")
    )

  unit_label <- if (is_fraction) "%" else "Δ"
  scale_text <- sprintf("±%g%s", ylim, unit_label)
  p <- p + annotate("text",
                     x = view_start + (view_end - view_start) * 0.05,
                     y = ylim * 0.85,
                     label = scale_text,
                     hjust = 0, vjust = 1,
                     size = 2.0, color = color)

  p <- p + annotate("text",
                     x = view_start + (view_end - view_start) * 0.015,
                     y = ylim * 0.55,
                     label = "mut-ctrl",
                     hjust = 0, vjust = 1,
                     size = 1.8, color = "grey40")

  p
}

# Gene model: thin black backbone, filled exon rectangles, strand arrowheads,
# italic gene symbol(s). Supports single-gene or multi-gene display.
build_gene_model <- function(view_start, view_end, chr, region_gr,
                             cfg, highlights = list()) {
  txdb <- get(cfg$txdb_pkg, envir = asNamespace(cfg$txdb_pkg))

  # Query genes overlapping the view region
  view_gr <- GRanges(seqnames = chr, ranges = IRanges(view_start, view_end))
  all_genes <- suppressMessages(GenomicFeatures::genes(txdb))
  overlap_genes <- subsetByOverlaps(all_genes, view_gr)

  # Pull a resolved focal entrez from the region's mcols if available (set by
  # resolve_region when --gene or --region + --gene-symbol was used).
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

  # Build per-gene rows
  orgdb <- if (load_optional(cfg$orgdb_pkg)) {
    get(cfg$orgdb_pkg, envir = asNamespace(cfg$orgdb_pkg))
  } else NULL

  gene_rows <- lapply(seq_along(target_genes), function(i) {
    g <- target_genes[i]
    entrez <- names(g)
    symbol <- if (!is.null(orgdb)) {
      tryCatch(
        AnnotationDbi::mapIds(orgdb, keys = entrez,
                              column = "SYMBOL", keytype = "ENTREZID"),
        error = function(e) entrez
      )
    } else entrez
    if (is.na(symbol)) symbol <- entrez

    # Get exons for this gene, reduced (union of isoforms)
    exons_by_gene <- suppressMessages(
      GenomicFeatures::exonsBy(txdb, by = "gene")
    )
    if (entrez %in% names(exons_by_gene)) {
      ex <- GenomicRanges::reduce(exons_by_gene[[entrez]])
    } else {
      ex <- GRanges()
    }

    list(
      symbol = as.character(symbol),
      gene_start = start(g),
      gene_end = end(g),
      strand = as.character(strand(g)),
      exons = ex
    )
  })

  # Assign each gene to a y-row to avoid overlap (simple layered layout)
  if (length(gene_rows) <= 1L || !cfg$show_all_genes) {
    for (i in seq_along(gene_rows)) gene_rows[[i]]$y <- 0
    n_rows <- 1L
  } else {
    # Greedy assignment: each gene goes on the lowest available row
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

  # Build data frames for ggplot layers
  backbone_df <- do.call(rbind, lapply(gene_rows, function(g) {
    data.frame(x = g$gene_start, xend = g$gene_end, y = g$y)
  }))

  exon_df_list <- lapply(gene_rows, function(g) {
    if (length(g$exons) == 0L) return(NULL)
    data.frame(
      xmin = start(g$exons),
      xmax = end(g$exons),
      ymin = g$y - 0.22,
      ymax = g$y + 0.22
    )
  })
  exon_df <- do.call(rbind, Filter(Negate(is.null), exon_df_list))

  # Symbol label positions (italic, near gene center, slightly above gene line)
  symbol_df <- do.call(rbind, lapply(gene_rows, function(g) {
    data.frame(
      x = (max(g$gene_start, view_start) + min(g$gene_end, view_end)) / 2,
      y = g$y + 0.85,
      label = g$symbol
    )
  }))

  # Strand arrowheads on introns (between exons), direction depends on strand
  # Only add for focal genes (when not show_all_genes), or for genes wider than threshold
  arrow_df_list <- lapply(gene_rows, function(g) {
    if (length(g$exons) < 2L) return(NULL)
    intron_starts <- end(g$exons)[-length(g$exons)]
    intron_ends   <- start(g$exons)[-1L]
    # Place 1-2 arrows per intron (at intron midpoint)
    midpoints <- (intron_starts + intron_ends) / 2
    intron_widths <- intron_ends - intron_starts
    # Skip introns too narrow to host an arrow
    keep <- intron_widths > (view_end - view_start) * 0.005
    if (!any(keep)) return(NULL)
    arrow_len <- (view_end - view_start) * 0.003
    if (g$strand == "+") {
      data.frame(
        x = midpoints[keep] - arrow_len,
        xend = midpoints[keep] + arrow_len,
        y = g$y, yend = g$y
      )
    } else {
      data.frame(
        x = midpoints[keep] + arrow_len,
        xend = midpoints[keep] - arrow_len,
        y = g$y, yend = g$y
      )
    }
  })
  arrow_df <- do.call(rbind, Filter(Negate(is.null), arrow_df_list))

  # Compose plot
  p <- ggplot()

  # Highlights first (behind everything)
  if (length(highlights) > 0L) {
    for (h in highlights) {
      p <- p + annotate("rect",
                        xmin = h$start, xmax = h$end,
                        ymin = -Inf, ymax = Inf,
                        fill = "grey85", alpha = 0.55, color = NA)
    }
  }

  if (!is.null(backbone_df) && nrow(backbone_df) > 0L) {
    p <- p +
      geom_segment(data = backbone_df,
                   aes(x = x, xend = xend, y = y, yend = y),
                   color = "black", linewidth = 0.15) +
      geom_rect(data = exon_df,
                aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                fill = "black", color = NA)

    if (!is.null(arrow_df) && nrow(arrow_df) > 0L) {
      p <- p + geom_segment(
        data = arrow_df,
        aes(x = x, xend = xend, y = y, yend = yend),
        arrow = arrow(length = unit(2.0, "pt"), type = "open", ends = "last"),
        color = "black", linewidth = 0.15
      )
    }

    p <- p + geom_text(
      data = symbol_df,
      aes(x = x, y = y, label = label),
      fontface = "italic", size = 2.7, color = "black"
    )
  }

  p <- p +
    scale_x_continuous(limits = c(view_start, view_end), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-0.7, n_rows + 0.3),
                       expand = expansion(mult = c(0.05, 0.05))) +
    theme_void() +
    theme(
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin      = margin(2, 2, 0, 0, "pt")
    )

  list(plot = p, n_rows = n_rows)
}

# Scale bar at top-right: short black segment with "X kb" label below
build_scale_bar <- function(view_start, view_end, scale_kb = NULL) {
  region_kb <- (view_end - view_start) / 1000

  # Auto-choose a nice scale bar length if not specified
  if (is.null(scale_kb)) {
    candidates <- c(1, 2, 5, 10, 20, 50, 100, 200, 500, 1000)
    target <- region_kb * 0.08
    scale_kb <- candidates[which.min(abs(candidates - target))]
  }

  scale_bp <- scale_kb * 1000
  # Place bar near right edge with small margin
  bar_end <- view_end - (view_end - view_start) * 0.02
  bar_start <- bar_end - scale_bp

  label <- if (scale_kb >= 1) {
    sprintf("%g kb", scale_kb)
  } else {
    sprintf("%g bp", scale_bp)
  }

  ggplot() +
    geom_segment(
      aes(x = bar_start, xend = bar_end, y = 0.35, yend = 0.35),
      linewidth = 1.2, color = "black"
    ) +
    annotate("text",
             x = (bar_start + bar_end) / 2, y = 0.78,
             label = label, size = 2.6, color = "black",
             hjust = 0.5, vjust = 0.5) +
    scale_x_continuous(limits = c(view_start, view_end), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme_void() +
    theme(
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin      = margin(0, 2, 0, 0, "pt")
    )
}

# Highlight header labels: optional text above highlight regions (e.g.,
# "differential region", "active enhancer"). Placed in a thin spacer panel
# above the first signal track.
build_highlight_headers <- function(view_start, view_end, highlights, labels) {
  if (length(highlights) == 0L) return(NULL)
  has_labels <- any(nzchar(unlist(labels)))
  if (!has_labels) return(NULL)

  label_df <- do.call(rbind, lapply(seq_along(highlights), function(i) {
    if (i > length(labels) || !nzchar(labels[[i]])) return(NULL)
    data.frame(
      x = (highlights[[i]]$start + highlights[[i]]$end) / 2,
      y = 0.5,
      label = labels[[i]]
    )
  }))

  if (is.null(label_df) || nrow(label_df) == 0L) return(NULL)

  ggplot(label_df, aes(x = x, y = y, label = label)) +
    geom_text(size = 2.6, color = "black", fontface = "plain",
              hjust = 0.5, vjust = 0.5) +
    scale_x_continuous(limits = c(view_start, view_end), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme_void() +
    theme(plot.background = element_rect(fill = "white", color = NA),
          plot.margin = margin(0, 2, 0, 0, "pt"))
}

# Hi-C arc panel: quadratic Bezier paths colored by direction
build_hic_arcs <- function(loops_df, view_start, view_end, chr,
                           lost_color = "#CB181D", gained_color = "#2171B5") {
  # Filter loops with at least one anchor in view
  in_view <- loops_df %>%
    dplyr::filter(
      (chr1 == chr & ((start1 >= view_start & start1 <= view_end) |
                      (end1 >= view_start & end1 <= view_end))) |
      (chr2 == chr & ((start2 >= view_start & start2 <= view_end) |
                      (end2 >= view_start & end2 <= view_end)))
    )

  if (nrow(in_view) == 0L) return(NULL)

  # Build Bezier arcs
  span_max <- max(in_view$start2 - in_view$start1)
  arc_df <- do.call(rbind, lapply(seq_len(nrow(in_view)), function(i) {
    x0 <- (in_view$start1[i] + in_view$end1[i]) / 2
    x1 <- (in_view$start2[i] + in_view$end2[i]) / 2
    # Height proportional to span
    h <- (abs(x1 - x0) / span_max)
    t_seq <- seq(0, 1, length.out = 60)
    x_arc <- (1 - t_seq)^2 * x0 + 2 * (1 - t_seq) * t_seq * ((x0 + x1) / 2) + t_seq^2 * x1
    y_arc <- 2 * (1 - t_seq) * t_seq * h
    data.frame(
      x = x_arc, y = y_arc,
      loop_id = i,
      direction = in_view$direction[i]
    )
  }))

  ggplot(arc_df, aes(x = x, y = y, group = loop_id, color = direction)) +
    geom_path(linewidth = 0.4, alpha = 0.75) +
    scale_color_manual(values = c(
      "down_in_mutant" = lost_color,
      "up_in_mutant"   = gained_color
    )) +
    scale_x_continuous(limits = c(view_start, view_end), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1.1), expand = c(0, 0)) +
    theme_void() +
    theme(
      legend.position = "none",
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin      = margin(0, 2, 0, 0, "pt")
    )
}

# =============================================================================
# LEFT-COLUMN LABELS -- rotated colored mark name + scale annotation
# =============================================================================

# Single rotated mark name panel, height-matched to its ctrl+mut pair
build_mark_label_panel <- function(mark_name, color) {
  cowplot::ggdraw() +
    cowplot::draw_label(
      mark_name,
      x = 0.6, y = 0.5,
      angle = 90,
      color = color,
      fontface = "bold",
      size = 10
    )
}

# Empty spacer panel for left column (matches non-track rows)
build_empty_label_panel <- function() {
  cowplot::ggdraw() +
    theme(plot.background = element_rect(fill = "white", color = NA))
}

# =============================================================================
# FIGURE ASSEMBLY
# =============================================================================

# Layout constants (overridable via cfg)
make_layout <- function(cfg) {
  list(
    track_h    = cfg$track_height,
    pair_gap   = 0.02,
    mark_gap   = 0.08,
    gene_h     = 0.22,  # thinner gene panel (~half a track) to match PI's style
    scalebar_h = 0.22,
    hic_h      = 1.00,
    header_h   = 0.18,
    fig_width  = cfg$fig_width,
    label_frac = 0.10
  )
}

# Assemble all panels into final figure
assemble_figure <- function(track_panels, mark_configs,
                            gene_result, scalebar_panel,
                            highlight_header_panel = NULL,
                            hic_panel = NULL,
                            cfg) {
  L <- make_layout(cfg)

  # ---- Right column: vertical stack of all data panels ----
  right_panels  <- list()
  right_heights <- numeric(0)

  # Top: scale bar
  right_panels[[length(right_panels) + 1L]] <- scalebar_panel
  right_heights <- c(right_heights, L$scalebar_h)

  # Optional highlight header row
  if (!is.null(highlight_header_panel)) {
    right_panels[[length(right_panels) + 1L]] <- highlight_header_panel
    right_heights <- c(right_heights, L$header_h)
  }

  # Mark groups: ctrl, mut, [optional diff], then mark_gap spacer
  mark_names <- names(track_panels)
  for (i in seq_along(mark_names)) {
    mn <- mark_names[i]
    right_panels[[length(right_panels) + 1L]] <- track_panels[[mn]]$ctrl
    right_heights <- c(right_heights, L$track_h)
    right_panels[[length(right_panels) + 1L]] <- track_panels[[mn]]$mut
    right_heights <- c(right_heights, L$track_h)
    if (!is.null(track_panels[[mn]]$diff)) {
      right_panels[[length(right_panels) + 1L]] <- track_panels[[mn]]$diff
      right_heights <- c(right_heights, L$track_h)
    }
    if (i < length(mark_names)) {
      right_panels[[length(right_panels) + 1L]] <- patchwork::plot_spacer()
      right_heights <- c(right_heights, L$mark_gap)
    }
  }

  # Gene model (height scales with number of gene rows)
  gene_h_eff <- L$gene_h * max(1, gene_result$n_rows)
  right_panels[[length(right_panels) + 1L]] <- gene_result$plot
  right_heights <- c(right_heights, gene_h_eff)

  # Optional Hi-C arcs at the bottom
  if (!is.null(hic_panel)) {
    right_panels[[length(right_panels) + 1L]] <- hic_panel
    right_heights <- c(right_heights, L$hic_h)
  }

  right_col <- patchwork::wrap_plots(
    right_panels, ncol = 1, heights = right_heights
  )

  # ---- Left column: rotated mark labels + scale annotations ----
  left_panels  <- list()
  left_heights <- numeric(0)

  # Spacer for scale bar row
  left_panels[[length(left_panels) + 1L]] <- build_empty_label_panel()
  left_heights <- c(left_heights, L$scalebar_h)

  # Spacer for highlight header row (if present)
  if (!is.null(highlight_header_panel)) {
    left_panels[[length(left_panels) + 1L]] <- build_empty_label_panel()
    left_heights <- c(left_heights, L$header_h)
  }

  # One rotated label panel per mark group (spans ctrl+mut+diff)
  for (i in seq_along(mark_names)) {
    mn <- mark_names[i]
    cfg_mark <- mark_configs[[mn]]
    n_rows <- if (cfg_mark$has_diff) 3L else 2L
    label_panel <- build_mark_label_panel(
      mark_name = cfg_mark$display_name,
      color     = cfg_mark$color
    )
    left_panels[[length(left_panels) + 1L]] <- label_panel
    left_heights <- c(left_heights, L$track_h * n_rows)

    if (i < length(mark_names)) {
      left_panels[[length(left_panels) + 1L]] <- build_empty_label_panel()
      left_heights <- c(left_heights, L$mark_gap)
    }
  }

  # Spacer for gene model
  left_panels[[length(left_panels) + 1L]] <- build_empty_label_panel()
  left_heights <- c(left_heights, gene_h_eff)

  # Spacer for Hi-C
  if (!is.null(hic_panel)) {
    left_panels[[length(left_panels) + 1L]] <- build_empty_label_panel()
    left_heights <- c(left_heights, L$hic_h)
  }

  left_col <- patchwork::wrap_plots(
    left_panels, ncol = 1, heights = left_heights
  )

  # ---- Combine left + right ----
  # Wrap each multi-panel column so patchwork treats them as single units
  final <- cowplot::plot_grid(
    left_col, right_col,
    ncol = 2,
    rel_widths = c(L$label_frac, 1 - L$label_frac),
    align = "h", axis = "tb"
  )

  # Total figure height (sum of all row heights + small margin)
  total_height <- sum(right_heights) + 0.3

  list(
    figure       = final,
    width        = L$fig_width,
    height       = total_height
  )
}

# =============================================================================
# OUTPUT -- save in PDF, SVG, JPEG formats
# =============================================================================

save_multiformat <- function(plot_obj, base_path, width, height, dpi = 300) {
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

  ggplot2::ggsave(paste0(prefix, ".pdf"), plot_obj, width = width, height = height)
  ggplot2::ggsave(paste0(prefix, ".svg"), plot_obj, width = width, height = height,
                  device = svglite::svglite)
  ggplot2::ggsave(paste0(prefix, ".png"), plot_obj, width = width, height = height,
                  dpi = dpi, device = "png")
  ggplot2::ggsave(paste0(prefix, ".jpg"), plot_obj, width = width, height = height,
                  dpi = dpi, device = "jpeg")

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
  cat("pub_browser.R -- publication-quality genome browser\n")
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

  # Echo gene model mode so the user can see what will be shown
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

  # Parse highlight regions (validate they fall within view)
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
    mark_data[[mk$name]] <- list(
      ctrl   = ctrl_gr,
      mut    = mut_gr,
      color  = mk$color,
      ylim   = if (!is.null(mk$ylim)) mk$ylim else auto_ylim,
      diff   = mk$diff,
      sparse = mk$sparse
    )
  }

  # ---- 3. Build track panels ----
  cat("Building track panels...\n")
  ctrl_label <- cfg$label_vec[1]
  mut_label  <- cfg$label_vec[2]

  track_panels <- list()
  mark_configs <- list()

  for (mn in names(mark_data)) {
    md <- mark_data[[mn]]
    ctrl_df <- bigwig_to_df(md$ctrl)
    mut_df  <- bigwig_to_df(md$mut)

    # Scale annotation shown above the ctrl track in mark color (PI style)
    scale_label_text <- sprintf("0-%g", md$ylim)

    ctrl_panel <- build_signal_panel(
      ctrl_df, md$color, md$ylim, view_start, view_end,
      condition_label = ctrl_label,
      scale_label     = scale_label_text,
      highlights      = cfg$highlights,
      italic_label    = FALSE
    )

    mut_panel <- build_signal_panel(
      mut_df, md$color, md$ylim, view_start, view_end,
      condition_label = mut_label,
      scale_label     = NULL,
      highlights      = cfg$highlights,
      italic_label    = cfg$genotype_italic
    )

    diff_panel <- NULL
    if (isTRUE(md$diff)) {
      diff_panel <- build_diff_panel(
        ctrl_df, mut_df, md$color,
        view_start, view_end,
        highlights    = cfg$highlights,
        is_fraction   = isTRUE(md$sparse)
      )
    }

    track_panels[[mn]] <- list(ctrl = ctrl_panel, mut = mut_panel,
                               diff = diff_panel)
    mark_configs[[mn]] <- list(
      display_name = mn,
      color        = md$color,
      ylim         = md$ylim,
      has_diff     = !is.null(diff_panel)
    )
  }

  # ---- 4. Build gene model ----
  cat("Building gene model...\n")
  gene_result <- build_gene_model(view_start, view_end, chr, region_gr,
                                  cfg, cfg$highlights)

  # ---- 5. Build scale bar ----
  scalebar_panel <- build_scale_bar(view_start, view_end, cfg$scale_bar_kb)

  # ---- 6. Optional highlight header labels ----
  highlight_header_panel <- build_highlight_headers(
    view_start, view_end, cfg$highlights, cfg$highlight_labels
  )

  # ---- 7. Optional Hi-C arcs ----
  hic_panel <- NULL
  if (!is.null(cfg$hic_loops)) {
    cat("Loading Hi-C loops...\n")
    loops_df <- read.table(cfg$hic_loops, header = TRUE, sep = "\t",
                           stringsAsFactors = FALSE)
    hic_panel <- build_hic_arcs(loops_df, view_start, view_end, chr)
  }

  # ---- 8. Assemble ----
  cat("Assembling figure...\n")
  result <- assemble_figure(
    track_panels, mark_configs,
    gene_result, scalebar_panel,
    highlight_header_panel = highlight_header_panel,
    hic_panel = hic_panel,
    cfg = cfg
  )

  # ---- 9. Save ----
  cat("Saving outputs...\n")
  save_multiformat(result$figure, cfg$output, result$width, result$height)

  cat(sprintf("\nDone. Figure dimensions: %.1f x %.1f inches\n",
              result$width, result$height))
}

# Run if invoked as script
if (sys.nframe() == 0L) {
  main()
}
