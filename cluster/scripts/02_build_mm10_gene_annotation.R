#!/usr/bin/env Rscript
# cluster/scripts/02_build_mm10_gene_annotation.R
# Build mm10 promoter BED (TSS +/- 750bp, 7-col gencode-style format) for
# use by Popay's bedpe_analysis.bedtools_annotation() and the loop
# classifier in Phase 4. Output schema matches
# cluster/clustering_example_data/gencode.v25.annotation_pp.bed:
#   chr  start  end  gene_id  score(=0)  strand  gene_name
# Width per row: 1500bp (end - start). BED is 0-based half-open.

suppressPackageStartupMessages({
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(org.Mm.eg.db)
  library(GenomicFeatures)
  library(GenomicRanges)
  library(GenomeInfoDb)
})

# Resolve repo root from this script's location (Rscript-friendly).
arg_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(arg_file) >= 1) {
  script_path <- normalizePath(sub("^--file=", "", arg_file[1]),
                               winslash = "/", mustWork = TRUE)
} else {
  stop("Could not determine script path; run via Rscript, not source().")
}
script_dir <- dirname(script_path)
repo_root  <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
out_path   <- file.path(repo_root, "cluster", "data", "mm10_knownGene_pp.bed")

cat("[1/5] Loading TxDb and pulling gene bodies...\n")
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
genes_gr <- genes(txdb, single.strand.genes.only = TRUE)
cat("      genes:", length(genes_gr), "\n")

cat("[2/5] Building promoters (TSS +/- 750bp, strand-aware)...\n")
prom_gr <- promoters(genes_gr, upstream = 750, downstream = 750)

# Drop entries that would fall off chromosome ends (preserving the strict
# 1500bp width Popay's downstream tooling assumes).
chrom_lengths <- seqlengths(prom_gr)
keep_in_bounds <- start(prom_gr) >= 1 &
                  end(prom_gr) <= chrom_lengths[as.character(seqnames(prom_gr))]
n_drop_oob <- sum(!keep_in_bounds)
prom_gr <- prom_gr[keep_in_bounds]
cat("      promoters:", length(prom_gr),
    "(dropped", n_drop_oob, "for falling outside chromosome bounds)\n")

cat("[3/5] Filtering to standard chromosomes (chr1..19, chrX, chrY, chrM)...\n")
std_chroms <- paste0("chr", c(1:19, "X", "Y", "M"))
prom_gr <- keepSeqlevels(prom_gr, std_chroms, pruning.mode = "coarse")
cat("      retained:", length(prom_gr), "\n")

cat("[4/5] Mapping Entrez IDs to gene symbols (org.Mm.eg.db)...\n")
entrez_ids <- names(prom_gr)
symbols <- suppressMessages(mapIds(
  org.Mm.eg.db, keys = entrez_ids, keytype = "ENTREZID",
  column = "SYMBOL", multiVals = "first"
))
symbols[is.na(symbols)] <- ""
n_no_symbol <- sum(symbols == "")
cat("      genes without symbol mapping (gene_name left blank):",
    n_no_symbol, "of", length(symbols), "\n")

cat("[5/5] Writing 7-col BED:", out_path, "\n")
bed_df <- data.frame(
  chr       = as.character(seqnames(prom_gr)),
  start     = start(prom_gr) - 1L,
  end       = as.integer(end(prom_gr)),
  gene_id   = entrez_ids,
  score     = 0L,
  strand    = as.character(strand(prom_gr)),
  gene_name = as.character(symbols),
  stringsAsFactors = FALSE
)

# Sanity: every row must be 1500bp wide before sort.
bad_width <- which((bed_df$end - bed_df$start) != 1500L)
if (length(bad_width) > 0) {
  stop(sprintf("FATAL: %d rows do not have width 1500 (e.g. row %d width %d)",
               length(bad_width), bad_width[1],
               bed_df$end[bad_width[1]] - bed_df$start[bad_width[1]]))
}

# Sort by chrom (numeric for chr1..19, then X/Y/M), then start.
chrom_order <- factor(bed_df$chr, levels = std_chroms)
bed_df <- bed_df[order(chrom_order, bed_df$start), ]

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
write.table(
  bed_df, file = out_path,
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE,
  eol = "\n"
)

cat("\n=== Summary ===\n")
cat("Output:", out_path, "\n")
cat("Rows:", nrow(bed_df), "\n")
cat("Unique chromosomes:", paste(sort(unique(bed_df$chr)), collapse = ", "), "\n")
strand_tab <- table(bed_df$strand)
cat("Strand counts:",
    paste(names(strand_tab), strand_tab, sep = "=", collapse = ", "), "\n")
cat("Width invariant (end - start == 1500):",
    all((bed_df$end - bed_df$start) == 1500L), "\n")
cat("Rows with empty gene_name:", sum(bed_df$gene_name == ""),
    "(expected: matches mapping-failure count above)\n")
cat("\nPhase 1.2 complete.\n")
