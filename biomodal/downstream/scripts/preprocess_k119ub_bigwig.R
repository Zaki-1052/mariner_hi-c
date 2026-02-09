# biomodal/downstream/scripts/preprocess_k119ub_bigwig.R
# HPC preprocessing: Extract K119ub continuous signal from merged bigwigs per gene
#
# Computes mean K119ub signal over gene bodies and promoters for ctrl/mut,
# applies median-ratio normalization, and outputs a portable TSV.
#
# Usage (on Expanse HPC):
#   Rscript scripts/preprocess_k119ub_bigwig.R \
#     --ctrl /expanse/lustre/projects/csd940/zalibhai/loop-class/H2AK119ub_ctrl_merged.bw \
#     --mut /expanse/lustre/projects/csd940/zalibhai/loop-class/H2AK119ub_mut_merged.bw \
#     --output data/k119ub_gene_signal.tsv

# =============================================================================
# PARSE CLI ARGUMENTS
# =============================================================================

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  defaults <- list(
    ctrl = "/expanse/lustre/projects/csd940/zalibhai/loop-class/H2AK119ub_ctrl_merged.bw",
    mut  = "/expanse/lustre/projects/csd940/zalibhai/loop-class/H2AK119ub_mut_merged.bw",
    output = "data/k119ub_gene_signal.tsv"
  )

  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--ctrl" && i < length(args)) {
      defaults$ctrl <- args[i + 1]; i <- i + 2
    } else if (args[i] == "--mut" && i < length(args)) {
      defaults$mut <- args[i + 1]; i <- i + 2
    } else if (args[i] == "--output" && i < length(args)) {
      defaults$output <- args[i + 1]; i <- i + 2
    } else {
      stop(sprintf("Unknown argument: %s\nUsage: Rscript preprocess_k119ub_bigwig.R --ctrl <bw> --mut <bw> --output <tsv>", args[i]))
    }
  }
  defaults
}

params <- parse_args()

cat("================================================================================\n")
cat("K119UB BIGWIG SIGNAL PREPROCESSING\n")
cat("================================================================================\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Ctrl bigwig:", params$ctrl, "\n")
cat("Mut bigwig: ", params$mut, "\n")
cat("Output:     ", params$output, "\n\n")

# =============================================================================
# LOAD PACKAGES
# =============================================================================

cat("Loading packages...\n")
suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicFeatures)
  library(GenomicRanges)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
})
cat("Packages loaded.\n\n")

# =============================================================================
# VALIDATE INPUTS
# =============================================================================

if (!file.exists(params$ctrl)) stop("Ctrl bigwig not found: ", params$ctrl)
if (!file.exists(params$mut))  stop("Mut bigwig not found: ", params$mut)

output_dir <- dirname(params$output)
if (output_dir != "" && output_dir != ".") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

# =============================================================================
# EXTRACT GENE AND PROMOTER COORDINATES
# =============================================================================

cat("Extracting gene body coordinates from TxDb...\n")
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# Gene bodies: standard chromosomes, >= 500bp
gene_gr <- genes(txdb)
std_chroms <- paste0("chr", c(1:19, "X", "Y"))
gene_gr <- keepSeqlevels(gene_gr, std_chroms, pruning.mode = "coarse")
gene_gr <- gene_gr[width(gene_gr) >= 500]
cat(sprintf("  Gene bodies: %d genes (standard chroms, >= 500bp)\n", length(gene_gr)))

# Promoters: +/- 2kb around TSS
cat("Extracting promoter coordinates (TSS +/- 2kb)...\n")
prom_gr <- promoters(txdb, upstream = 2000, downstream = 2000)
prom_gr <- keepSeqlevels(prom_gr, std_chroms, pruning.mode = "coarse")

# Deduplicate promoters to one per gene (take the first transcript per gene)
prom_gene_ids <- AnnotationDbi::select(txdb,
  keys = names(prom_gr),
  keytype = "TXNAME",
  columns = "GENEID"
)
prom_gr$gene_id <- prom_gene_ids$GENEID[match(names(prom_gr), prom_gene_ids$TXNAME)]
prom_gr <- prom_gr[!is.na(prom_gr$gene_id)]

# Keep only one promoter per gene (first occurrence)
prom_gr <- prom_gr[!duplicated(prom_gr$gene_id)]
names(prom_gr) <- prom_gr$gene_id

# Restrict to genes that passed the gene body filter
shared_genes <- intersect(names(gene_gr), names(prom_gr))
gene_gr <- gene_gr[shared_genes]
prom_gr <- prom_gr[shared_genes]
cat(sprintf("  Promoters: %d (one per gene, matched to gene bodies)\n", length(prom_gr)))

# Map ENTREZID to SYMBOL
cat("Mapping ENTREZ IDs to gene symbols...\n")
symbols <- mapIds(org.Mm.eg.db,
  keys = names(gene_gr),
  keytype = "ENTREZID",
  column = "SYMBOL",
  multiVals = "first"
)
cat(sprintf("  Mapped: %d / %d genes have symbols\n",
            sum(!is.na(symbols)), length(symbols)))

# =============================================================================
# COMPUTE MEAN SIGNAL PER REGION
# =============================================================================

compute_mean_signal <- function(bw_path, regions, label) {
  cat(sprintf("  Computing mean signal from %s (%d regions)...\n", label, length(regions)))

  # Import bigwig only for queried regions (memory efficient)
  bw <- import.bw(bw_path, which = regions)

  # Build coverage weighted by score
  cov <- coverage(bw, weight = "score")

  # Extract views for each region and compute means
  v <- Views(cov, split(ranges(regions), seqnames(regions)))
  means <- viewMeans(v)

  # Flatten from RleViewsList to numeric vector aligned to regions
  result <- numeric(length(regions))
  chr_names <- as.character(seqnames(regions))
  for (chr in names(means)) {
    idx <- which(chr_names == chr)
    result[idx] <- as.numeric(means[[chr]])
  }

  cat(sprintf("    Range: [%.4f, %.4f], Median: %.4f, Non-zero: %d / %d\n",
              min(result), max(result), median(result),
              sum(result > 0), length(result)))
  result
}

cat("\n--- Gene body signal extraction ---\n")
gb_ctrl_raw <- compute_mean_signal(params$ctrl, gene_gr, "ctrl gene bodies")
gb_mut_raw  <- compute_mean_signal(params$mut,  gene_gr, "mut gene bodies")

cat("\n--- Promoter signal extraction ---\n")
pr_ctrl_raw <- compute_mean_signal(params$ctrl, prom_gr, "ctrl promoters")
pr_mut_raw  <- compute_mean_signal(params$mut,  prom_gr, "mut promoters")

# =============================================================================
# BETWEEN-CONDITION NORMALIZATION (median ratio)
# =============================================================================

cat("\n--- Between-condition normalization ---\n")

# Use gene body signals for normalization (more stable than promoters)
ctrl_nonzero <- gb_ctrl_raw[gb_ctrl_raw > 0]
mut_nonzero  <- gb_mut_raw[gb_mut_raw > 0]

median_ctrl <- median(ctrl_nonzero)
median_mut  <- median(mut_nonzero)
scale_factor <- median_ctrl / median_mut

cat(sprintf("  Ctrl median (non-zero gene bodies): %.6f\n", median_ctrl))
cat(sprintf("  Mut median (non-zero gene bodies):  %.6f\n", median_mut))
cat(sprintf("  Scale factor (ctrl/mut): %.6f\n", scale_factor))
cat(sprintf("  Interpretation: mut signals multiplied by %.4f to match ctrl scale\n", scale_factor))

# Apply scaling to mut signals
gb_mut_norm <- gb_mut_raw * scale_factor
pr_mut_norm <- pr_mut_raw * scale_factor

# =============================================================================
# COMPUTE LOG2FC WITH ADAPTIVE PSEUDOCOUNT
# =============================================================================

cat("\n--- Computing log2 fold changes ---\n")

# Adaptive pseudocount: 5th percentile of non-zero signals, floor at 0.01
all_nonzero <- c(gb_ctrl_raw[gb_ctrl_raw > 0], gb_mut_norm[gb_mut_norm > 0])
pseudocount <- max(quantile(all_nonzero, 0.05), 0.01)
cat(sprintf("  Adaptive pseudocount: %.6f (5th pctl of non-zero signals)\n", pseudocount))

compute_log2fc <- function(ctrl, mut, pc) {
  log2fc <- rep(NA_real_, length(ctrl))
  signal_class <- rep("no_signal", length(ctrl))

  # Quantifiable: both above pseudocount
  both_above <- (ctrl > pc) & (mut > pc)
  log2fc[both_above] <- log2((mut[both_above] + pc) / (ctrl[both_above] + pc))
  signal_class[both_above] <- "quantifiable"

  # One condition: one above pseudocount
  one_above <- xor(ctrl > pc, mut > pc)
  log2fc[one_above] <- log2((mut[one_above] + pc) / (ctrl[one_above] + pc))
  signal_class[one_above] <- "one_condition"

  # No signal: both <= pseudocount → NA
  list(log2fc = log2fc, signal_class = signal_class)
}

gb_result <- compute_log2fc(gb_ctrl_raw, gb_mut_norm, pseudocount)
pr_result <- compute_log2fc(pr_ctrl_raw, pr_mut_norm, pseudocount)

# Report class distribution
for (region_label in c("Gene body", "Promoter")) {
  classes <- if (region_label == "Gene body") gb_result$signal_class else pr_result$signal_class
  cat(sprintf("  %s signal classes:\n", region_label))
  for (cls in c("quantifiable", "one_condition", "no_signal")) {
    n <- sum(classes == cls)
    cat(sprintf("    %-15s: %d (%.1f%%)\n", cls, n, 100 * n / length(classes)))
  }
}

# =============================================================================
# BUILD OUTPUT TABLE
# =============================================================================

cat("\n--- Building output table ---\n")

gene_coords <- as.data.frame(gene_gr)

output_df <- data.frame(
  entrez_id      = names(gene_gr),
  symbol         = as.character(symbols[names(gene_gr)]),
  chr            = gene_coords$seqnames,
  start          = gene_coords$start,
  end            = gene_coords$end,
  width          = gene_coords$width,
  gb_ctrl_signal = round(gb_ctrl_raw, 6),
  gb_mut_signal  = round(gb_mut_norm, 6),
  gb_log2fc      = round(gb_result$log2fc, 6),
  gb_signal_class = gb_result$signal_class,
  pr_ctrl_signal = round(pr_ctrl_raw, 6),
  pr_mut_signal  = round(pr_mut_norm, 6),
  pr_log2fc      = round(pr_result$log2fc, 6),
  pr_signal_class = pr_result$signal_class,
  stringsAsFactors = FALSE
)

# Remove genes without symbols
n_before <- nrow(output_df)
output_df <- output_df[!is.na(output_df$symbol), ]
cat(sprintf("  Removed %d genes without symbols (%d remaining)\n",
            n_before - nrow(output_df), nrow(output_df)))

# =============================================================================
# WRITE OUTPUT
# =============================================================================

cat(sprintf("\nWriting output to: %s\n", params$output))
write.table(output_df, params$output, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  %d genes written\n", nrow(output_df)))

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n================================================================================\n")
cat("PREPROCESSING SUMMARY\n")
cat("================================================================================\n")
cat(sprintf("  Total genes:           %d\n", nrow(output_df)))
cat(sprintf("  Pseudocount:           %.6f\n", pseudocount))
cat(sprintf("  Normalization factor:  %.6f (applied to mut)\n", scale_factor))
cat(sprintf("  Raw median ctrl:       %.6f\n", median_ctrl))
cat(sprintf("  Raw median mut:        %.6f\n", median_mut))

# Gene body log2FC distribution (quantifiable only)
gb_q <- output_df$gb_log2fc[output_df$gb_signal_class == "quantifiable"]
if (length(gb_q) > 0) {
  cat(sprintf("\n  Gene body log2FC (quantifiable, n=%d):\n", length(gb_q)))
  cat(sprintf("    Median: %+.4f\n", median(gb_q)))
  cat(sprintf("    Mean:   %+.4f\n", mean(gb_q)))
  cat(sprintf("    SD:     %.4f\n", sd(gb_q)))
  cat(sprintf("    Range:  [%+.4f, %+.4f]\n", min(gb_q), max(gb_q)))
  cat(sprintf("    Positive (mut > ctrl): %d (%.1f%%)\n",
              sum(gb_q > 0), 100 * sum(gb_q > 0) / length(gb_q)))
}

# Promoter log2FC distribution (quantifiable only)
pr_q <- output_df$pr_log2fc[output_df$pr_signal_class == "quantifiable"]
if (length(pr_q) > 0) {
  cat(sprintf("\n  Promoter log2FC (quantifiable, n=%d):\n", length(pr_q)))
  cat(sprintf("    Median: %+.4f\n", median(pr_q)))
  cat(sprintf("    Mean:   %+.4f\n", mean(pr_q)))
  cat(sprintf("    SD:     %.4f\n", sd(pr_q)))
  cat(sprintf("    Range:  [%+.4f, %+.4f]\n", min(pr_q), max(pr_q)))
}

cat("\nPreprocessing complete.\n")
